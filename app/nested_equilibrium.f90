program nested_equilibrium
    ! Nested 3D equilibrium from a boundary, the toroidal flux and the
    ! profiles iota(s), p(s) in the normalized toroidal flux s, solved by the
    ! native map Newton (docs/nested-equilibrium-input.md). With free_boundary
    ! the boundary is the initial guess and the plasma sits in the vacuum field
    ! of coil groups tabulated in an mgrid file, scaled by extcur
    ! (milca_ideal_free_boundary). Writes the solved state's global quantities
    ! for cross-code comparison.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh, hermite_profile_evaluate, &
        ideal_seed_3d_angle_density
    use milca_ideal_map_leaves, only: ideal_point_gradient_jets, &
        ideal_point_field_gradient_jets
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce, ideal_map_prolong, ideal_map_prolong_correction, &
        ideal_map_freeze, ideal_map_sensitivity, IDEAL_MAP_OK
    use milca_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use milca_map_3d_force, only: map_3d_force, map_3d_force_field
    use milca_flux_map_force, only: flux_map_force_t, flux_map_force_print
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use milca_grad_shafranov, only: gs_state_t
    use milca_curved_triangle, only: triangle_locate
    use milca_result_io, only: run_t, result_file_t
    use fortnum_quadrature, only: gauss_legendre
    use milca_map_3d_loop, only: map_3d_loop_t, map_3d_loop_setup, map_3d_loop_integral
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use milca_ideal_seed_3d_leaves, only: seed_hermite_profile_jvp
    use milca_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t, curve_chart_t
    use milca_deck_chart, only: deck_axis_chart, deck_boundary_to_chart
    use milca_ideal_free_boundary, only: step_cap_fraction, fb_prof_seconds, fb_prof_calls, ideal_fb_t, ideal_fb_init, ideal_fb_solve, &
        ideal_fb_report_t, ideal_fb_drive_traces, ideal_fb_set_traces, &
        ideal_fb_gradient, ideal_fb_sensitivity, IDEAL_FB_OK, IDEAL_FB_FOLDED, &
        ideal_fb_init_bem, ideal_fb_bem_settle, ideal_fb_fix_vertical, ideal_fb_set_coil_split, ideal_fb_enable_angle, &
        ideal_fb_solve_shell_map, ideal_fb_solve_vacuum, ideal_fb_solve_multipliers, &
        ideal_fb_hessian, ideal_fb_step
    use milca_ideal_fb_bem, only: fb_coil_flux_t
    use milca_ideal_fb_circuits, only: ideal_fb_circuits_t, ideal_fb_circuits_build
    use milca_ideal_fb_mgrid, only: mgrid_basis_t, mgrid_basis_init, mgrid_current_nodes
    use milca_ideal_fb_nested, only: fb_nested_geometry, fb_nested_mesh, &
        fb_nested_seed, fb_nested_transfer, fb_nested_settle, filament_basis_t, &
        fb_nested_plasma_mesh, fb_nested_coil_flux
    use milca_mgrid_file, only: mgrid_t, mgrid_read
    use fortfem_gs_single_layer_shape_2d, only: sl_near_order => near_order, &
        sl_far_order => far_order, sl_radial_n => radial_n
    implicit none
    ! namelist storage: the toroidal modes a boundary and map need follow
    ! from the spectrum and fourier_tolerance, memory grows like map_ntor^2
    integer, parameter :: mmax = 64, nmax = 256, kmax = 32, knot_max = 1024
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi
    integer :: nfp, mpol, ntor, boundary_nodes, degree, max_iterations, angle_mpol
    integer :: unit, ios, status, solve_status, step, steps, nodes, order, n
    integer :: ladder_nodes(8), ladder_degree(8)
    integer :: map_ntor, ntor_level
    ! iota of the current matching: polynomial through iota_nodes
    ! Chebyshev-Gauss surfaces, stored as a Hermite profile on store_knots
    integer :: iota_nodes, iota_level
    integer, parameter :: store_knots = 129, iota_nodes_max = 32
    integer, parameter :: mesh_nodes(7) = [16, 24, 32, 48, 64, 96, 128]
    logical :: automatic_ntor, automatic_mesh, automatic_iota
    logical :: mesh_confirmed, modes_confirmed, iota_confirmed, converged
    logical :: iota_vanishes = .false.
    ! the solved map of the level folds at the volume-integral points
    logical :: level_folded = .false.
    real(dp) :: fourier_tolerance, fourier_error, mesh_tolerance, mesh_error
    real(dp) :: iota_tolerance, iota_error
    real(dp) :: observed(4), previous_observed(4)
    real(dp), allocatable :: boundary_input(:, :, :)
    ! free boundary: seeds from the solved boundary (axisymmetric decks)
    integer, parameter :: max_reseeds = 10, reseed_modes = 16, coarse_nodes = 16
    real(dp), parameter :: reseed_tolerance = 1e-5_dp
    integer :: reseed_step
    real(dp) :: reseed_move, reseed_change, reseed_seen(4), reseed_last(4)
    ! the free-boundary residual target of a level: loose for the seeds of the
    ! reseed loop, whose boundaries only start the next seed
    real(dp) :: level_tolerance = 0
    real(dp), parameter :: seed_tolerance = 1e-5_dp
    type(ideal_map_t) :: previous_map
    real(dp), allocatable :: previous_seed(:, :, :)
    real(dp) :: phiedge, tolerance
    real(dp) :: rbc(-nmax:nmax, 0:mmax), rbs(-nmax:nmax, 0:mmax)
    real(dp) :: zbc(-nmax:nmax, 0:mmax), zbs(-nmax:nmax, 0:mmax)
    real(dp) :: iota(0:kmax), pressure(0:kmax), current(0:kmax)
    character(16) :: iota_type, pressure_type, current_type
    integer :: n_knots
    real(dp) :: profile_knots(knot_max), iota_values(knot_max), iota_slopes(knot_max)
    real(dp) :: pressure_values(knot_max), pressure_slopes(knot_max)
    real(dp) :: current_values(knot_max), current_slopes(knot_max)
    logical :: prescribe_current, stationary
    integer :: max_current_iterations
    real(dp) :: current_tolerance, current_error, solve_tolerance
    integer :: current_iterations
    real(dp), allocatable :: boundary(:, :, :), seed_x(:, :, :), coarse_seed(:, :, :)
    real(dp) :: residual, reference, volume, area, magnetic, kinetic, raxis
    real(dp) :: iota_axis, iota_edge, p_axis, unused, ctor
    type(run_t) :: run
    type(result_file_t) :: result
    type(gs_state_t) :: mesh
    type(ideal_map_t) :: map, coarse
    type(toroidal_seed_t) :: seed
    type(toroidal_profile_t) :: profile
    type(ideal_map_solve_report_t) :: report
    integer :: newton_total = 0 ! Newton steps of all fixed-boundary solves
    character(256) :: message
    ! Free boundary: coil groups of mgrid_file scaled by extcur, mgrid
    ! interpolation degree in (R, Z), and the clearance of the coupling
    ! surface Gamma_c around the guess in units of the minor radius.
    integer, parameter :: extcur_max = 512
    logical :: free_boundary
    ! the free-boundary solver is active (off for the fixed-boundary interior
    ! of the seed, which matches the current on the guess)
    logical :: fb_solving = .false.
    character(512) :: mgrid_file
    real(dp) :: extcur(extcur_max), coupling_margin, margin_used
    integer :: mgrid_degree
    type(ideal_fb_t) :: fb
    type(ideal_fb_report_t) :: fb_report
    type(mgrid_basis_t) :: coil_basis
    type(filament_basis_t) :: filament
    type(ideal_fb_circuits_t) :: circuits
    type(cylindrical_chart_t) :: cylinder
    ! Chart of the map: 'cylindrical' (R, Z, phi) or 'deck', the curve chart
    ! along the boundary's m = 0 curve (milca_deck_chart), whose sections
    ! stay centred on a helical boundary.
    character(16) :: chart
    type(curve_chart_t) :: deck_frame
    class(ideal_chart_t), allocatable :: map_chart
    ! poloidal modes of the boundary in the deck chart's geometric angle
    integer, parameter :: deck_mpol = 32
    real(dp), allocatable :: filament_drive(:, :), filament_load(:, :)
    real(dp) :: fb_centre(2), fb_axes(2), circulation, coil_closure
    real(dp), allocatable :: fb_curve(:), fit_curve(:)
    ! Vacuum of the free boundary: 'bem' (sheet current on Gamma_p, no vacuum
    ! mesh; axisymmetric decks), 'shell' (shell to Gamma_c) or 'auto' (bem for
    ! axisymmetric decks, otherwise shell).
    character(8) :: vacuum
    logical :: use_bem = .false., split_coils = .false.
    type(fb_coil_flux_t) :: coil_flux
    real(dp) :: bem_f_vac, bem_f_spread
    namelist /nested_input/ nfp, mpol, ntor, rbc, rbs, zbc, zbs, phiedge, &
        iota, pressure, boundary_nodes, degree, tolerance, max_iterations, angle_mpol, &
        ladder_nodes, ladder_degree, stationary, prescribe_current, current, max_current_iterations, &
        current_tolerance, iota_type, pressure_type, &
        current_type, n_knots, profile_knots, iota_values, iota_slopes, pressure_values, &
        pressure_slopes, current_values, current_slopes, map_ntor, fourier_tolerance, mesh_tolerance, &
        iota_nodes, iota_tolerance, free_boundary, mgrid_file, extcur, mgrid_degree, &
        coupling_margin, vacuum, chart

    nfp = 1; mpol = 1; ntor = 0; phiedge = 1
    rbc = 0; rbs = 0; zbc = 0; zbs = 0; iota = 0; pressure = 0
    boundary_nodes = -1; degree = 3; tolerance = 1e-8_dp; max_iterations = 200
    angle_mpol = -1; ladder_nodes = 0; ladder_degree = 0
    stationary = .false.
    prescribe_current = .false.; current = 0; max_current_iterations = 30
    current_tolerance = 1e-5_dp
    map_ntor = -1; fourier_tolerance = 1e-4_dp
    mesh_tolerance = 1e-4_dp
    iota_nodes = -1; iota_tolerance = 1e-4_dp
    free_boundary = .false.; mgrid_file = ''; extcur = 0; mgrid_degree = 5
    coupling_margin = 0.5_dp; vacuum = 'auto'; chart = 'cylindrical'
    iota_type = 'power_series'; pressure_type = 'power_series'; current_type = 'power_series'
    n_knots = 0; profile_knots = 0; iota_values = 0; iota_slopes = 0
    pressure_values = 0; pressure_slopes = 0; current_values = 0; current_slopes = 0
    call run%start('nested_equilibrium')
    call run%open_config(unit)
    ios = 0; message = ''
    if (run%has_config()) read (unit, nml=nested_input, iostat=ios, iomsg=message)
    call run%close_config(unit, 'nested_input', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: nested_equilibrium [--config FILE]'
    if (nfp < 1 .or. mpol < 1 .or. mpol > mmax .or. ntor < 0 .or. ntor > nmax) &
        error stop 'nested_equilibrium: mpol, ntor or nfp out of range'
    if (.not. phiedge > 0 .or. .not. ieee_is_finite(phiedge)) &
        error stop 'nested_equilibrium: phiedge must be positive'
    if (.not. all([kind_ok(iota_type), kind_ok(pressure_type), kind_ok(current_type)])) &
        error stop 'nested_equilibrium: profile types are power_series or hermite'
    if (any([iota_type, pressure_type, current_type] == 'hermite') .and. &
        (n_knots < 2 .or. n_knots > knot_max)) &
        error stop 'nested_equilibrium: hermite profiles need 2 to 1024 knots'
    if (map_ntor > nmax .or. (map_ntor >= 0 .and. map_ntor < ntor)) &
        error stop 'nested_equilibrium: map_ntor below the boundary ntor or above 256'
    if (.not. fourier_tolerance > 0) error stop 'nested_equilibrium: fourier_tolerance'
    if (.not. mesh_tolerance > 0) error stop 'nested_equilibrium: mesh_tolerance'
    if (.not. iota_tolerance > 0) error stop 'nested_equilibrium: iota_tolerance'
    if (iota_nodes > iota_nodes_max .or. (iota_nodes >= 0 .and. iota_nodes < 2)) &
        error stop 'nested_equilibrium: iota_nodes from 2 to 32, or -1'
    if (stationary .and. (free_boundary .or. prescribe_current)) &
        error stop 'nested_equilibrium: stationary requires a fixed boundary and prescribed iota'
    if (free_boundary .and. len_trim(mgrid_file) == 0) &
        error stop 'nested_equilibrium: free_boundary needs mgrid_file'
    if (free_boundary .and. .not. coupling_margin > 0) &
        error stop 'nested_equilibrium: coupling_margin must be positive'
    if (vacuum /= 'auto' .and. vacuum /= 'shell' .and. vacuum /= 'bem' .and. &
        vacuum /= 'split') error stop 'nested_equilibrium: vacuum is auto, bem, shell or split'
    split_coils = vacuum == 'split'
    if (free_boundary .and. split_coils .and. (ntor /= 0 .or. map_ntor > 0)) &
        error stop 'nested_equilibrium: vacuum = split needs an axisymmetric deck'
    if (chart /= 'cylindrical' .and. chart /= 'deck') &
        error stop 'nested_equilibrium: chart is cylindrical or deck'
    if (chart == 'deck' .and. (free_boundary .or. prescribe_current)) &
        error stop 'nested_equilibrium: chart = deck needs a fixed boundary and prescribed iota'
    use_bem = vacuum == 'bem' .or. (vacuum == 'auto' .and. ntor == 0 .and. map_ntor <= 0)
    if (free_boundary .and. use_bem .and. (ntor /= 0 .or. map_ntor > 0)) &
        error stop 'nested_equilibrium: vacuum = bem needs an axisymmetric deck'
    call run%open_resolved(unit)
    write (unit, nml=nested_input)
    call run%close_resolved(unit)

    ! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp phi)
    allocate (boundary_input(4, -ntor:ntor, 0:mpol))
    boundary_input(1, :, :) = rbc(-ntor:ntor, 0:mpol)
    boundary_input(2, :, :) = rbs(-ntor:ntor, 0:mpol)
    boundary_input(3, :, :) = zbc(-ntor:ntor, 0:mpol)
    boundary_input(4, :, :) = zbs(-ntor:ntor, 0:mpol)
    if (chart == 'deck') then
        call deck_axis_chart(boundary_input, ntor, nfp, deck_frame, status)
        if (status /= 0) call fail('the chart of the boundary''s m = 0 curve is not admitted')
        allocate (map_chart, source=deck_frame)
    else
        allocate (map_chart, source=cylinder)
    end if
    profile%iota = iota(0:last(iota))
    ! the solver uses mu0 = 1 pressure units
    profile%pressure = mu0*pressure(0:last(pressure))
    if (iota_type == 'hermite') call profile%set_hermite_iota(profile_knots(:n_knots), &
        iota_values(:n_knots), iota_slopes(:n_knots))
    if (pressure_type == 'hermite') call profile%set_hermite_pressure( &
        profile_knots(:n_knots), mu0*pressure_values(:n_knots), mu0*pressure_slopes(:n_knots))
    call profile%evaluate(0.0_dp, iota_axis, p_axis, unused)
    call profile%evaluate(1.0_dp, iota_edge, unused, unused)
    if (prescribe_current .and. all(boundary_input(:, -ntor:-1, :) == 0) .and. &
        all(boundary_input(:, 1:ntor, :) == 0)) call admit_axisymmetric_current()

    ! Resolution under error control: the mesh (boundary_nodes < 0), the
    ! interior toroidal modes (map_ntor < 0; a 3D boundary's interior needs more
    ! harmonics than the boundary) and the nodes of the matched iota
    ! (iota_nodes < 0) are raised until the observables change by at most
    ! their tolerance between successive levels. The active ladders alternate
    ! until each is confirmed by one refinement or is at its finest level;
    ! axisymmetric data keep ntor = 0.
    automatic_mesh = boundary_nodes < 0
    automatic_ntor = map_ntor < 0 .and. ntor > 0
    ! a change between levels estimates the truncation error only when each
    ! level's current is matched well below it
    if (automatic_ntor) current_tolerance = min(current_tolerance, fourier_tolerance/10)
    if (automatic_mesh) current_tolerance = min(current_tolerance, mesh_tolerance/10)
    automatic_iota = prescribe_current .and. iota_nodes < 0
    if (automatic_iota) current_tolerance = min(current_tolerance, iota_tolerance/10)
    iota_level = merge(6, iota_nodes, automatic_iota)
    if (automatic_mesh) boundary_nodes = mesh_nodes(1)
    ntor_level = merge(ntor, map_ntor, map_ntor < 0)
    mesh_error = 0; fourier_error = 0; iota_error = 0
    call set_toroidal_level(ntor_level)
    if (free_boundary) then
        call set_bem_orders(boundary_nodes)
        call free_boundary_setup()
        ! the frozen field carries the seed's labelling: seed again from the
        ! solved boundary until it stops moving (axisymmetric decks), on a
        ! coarse mesh first, then on the requested one from that boundary
        ! the first solves only start the seeds: loose residual target
        if (ntor == 0) level_tolerance = seed_tolerance
        if (ntor == 0 .and. boundary_nodes > coarse_nodes) then
            call fb_level(coarse_nodes, .true.)
            call reseed_loop(coarse_nodes, max_reseeds)
            call reseed_from_solution(reseed_move)
            call fb_level(boundary_nodes, .true.)
            call reseed_loop(boundary_nodes, 2)
        else
            call fb_level(boundary_nodes, .true.)
            if (ntor == 0) call reseed_loop(boundary_nodes, max_reseeds)
        end if
    else
        call mesh_ladder()
    end if
    call finish_level()
    call refine_folded()
    do
        ! an iota without a resolved sign ends the resolution ladders
        if (iota_vanishes) exit
        mesh_confirmed = .true.; modes_confirmed = .true.; iota_confirmed = .true.
        if (automatic_mesh) call raise_mesh(mesh_confirmed)
        if (iota_vanishes) exit
        if (automatic_ntor) call raise_modes(modes_confirmed)
        if (iota_vanishes) exit
        if (automatic_iota) call raise_iota(iota_confirmed)
        if (iota_vanishes) exit
        if (count([automatic_mesh, automatic_ntor, automatic_iota]) < 2) exit
        if (mesh_confirmed .and. modes_confirmed .and. iota_confirmed) exit
    end do
    call profile%evaluate(0.0_dp, iota_axis, unused, unused)
    call profile%evaluate(1.0_dp, iota_edge, unused, unused)
    if (free_boundary) then
        residual = fb_report%residual
    else
        residual = free_gradient()/max(tiny(1.0_dp), reference)
    end if
    call integrals(volume, area, magnetic, kinetic, status)
    if (status /= 0) call fail('a volume integral met a folded map')
    call axis_radius(raxis, status)
    if (status /= 0) call fail('the magnetic axis was not found')
    ctor = edge_current(status)
    if (status /= 0) call fail('the edge current readback met a folded map')

    call result%create(run)
    call result%attribute('physical_model', merge('free-boundary ', 'fixed-boundary', &
        free_boundary)//' nested ideal MHD, prescribed iota')
    converged = solve_status == 0 .and. residual <= tolerance .and. current_error <= current_tolerance &
        .and. .not. iota_vanishes &
        .and. (.not. automatic_ntor .or. fourier_error <= fourier_tolerance) &
        .and. (.not. automatic_mesh .or. mesh_error <= mesh_tolerance) &
        .and. (.not. automatic_iota .or. iota_error <= iota_tolerance)
    call result%attribute('converged', merge('true ', 'false', converged))
    call result%attribute('chart', trim(chart))
    call result%write('volume_p', volume, 'm^3', 'plasma volume')
    call result%write('aminor_p', sqrt(area/pi), 'm', 'sqrt of mean cross-section area / pi')
    call result%write('rmajor_p', volume/(2*pi*area), 'm', 'volume / (2 pi mean area)')
    call result%write('aspect', volume/(2*pi*area)/sqrt(area/pi), '1', 'rmajor_p / aminor_p')
    call result%write('raxis_cc', raxis, 'm', 'toroidal mean of the axis major radius')
    call result%write('iotaf_axis', iota_axis, '1', 'rotational transform on axis')
    call result%write('iotaf_edge', iota_edge, '1', 'rotational transform at the edge')
    call result%write('pressure_axis', p_axis/mu0, 'Pa', 'pressure on axis')
    call result%write('wb', magnetic/(4*pi**2), 'T^2 m^3', &
        'VMEC normalisation int B^2/2 dV / (2 pi)^2 (mu0 = 1)')
    call result%write('betatotal', kinetic/magnetic, '1', 'int p dV / int B^2/2 dV')
    call result%write('ctor', ctor, 'A', 'enclosed toroidal current at the edge')
    if (free_boundary) then
        call result%attribute('vacuum', merge('bem  ', 'shell', use_bem))
        call result%write('coupling_margin', margin_used, '1', &
            'clearance of the coupling surface around the guess / minor radius')
        call result%write('coil_closure', coil_closure, 'T m', &
            'largest curl-free mismatch of the coil traces on the coupling surface')
        call result%write('circulation', circulation, 'T m', &
            'poloidal circulation of the vacuum field around the plasma')
        call result%write('vacuum_f', fb%f_vac, 'T m', 'R B_phi of the coil field')
    end if
    if (prescribe_current) then
        call result%attribute('current_status', merge('iota vanishes', 'matched      ', &
            iota_vanishes))
        if (iota_vanishes) call result%attribute('failure', 'the iota matched to the' &
            //' prescribed current comes within current_tolerance of zero, which the' &
            //' frozen map does not admit')
        call result%write('current_error', current_error, '1', &
            'largest iota change of the last Newton step on the prescribed current')
        call result%write('current_iterations', real(current_iterations, dp), '1', &
            'iota updates to match the prescribed current')
    end if
    call result%write('map_ntor', real(ntor_level, dp), '1', 'interior toroidal modes')
    call result%write('boundary_nodes', real(boundary_nodes, dp), '1', &
        'boundary nodes of the reference disk mesh')
    if (automatic_mesh) call result%write('mesh_error', mesh_error, '1', &
        'largest change of wb (relative), raxis/aminor_p and iota on axis and edge' &
        //' from the previous mesh')
    if (prescribe_current) call result%write('iota_nodes', real(iota_level, dp), '1', &
        'Chebyshev-Gauss surfaces of the polynomial iota matched to the current')
    if (automatic_iota) call result%write('iota_error', iota_error, '1', &
        'largest change of wb (relative), raxis/aminor_p and iota on axis and edge' &
        //' from iota_nodes - 2')
    if (automatic_ntor) call result%write('fourier_error', fourier_error, '1', &
        'largest change of wb (relative), raxis/aminor_p and iota on axis and edge' &
        //' from map_ntor - 2')
    call result%write('force_residual', residual, '1', 'relative free reduced-gradient l2')
    call result%write('tolerance', tolerance, '1', 'requested force_residual')
    if (free_boundary) report%iterations = fb_report%iterations
    call result%write('iterations', real(report%iterations, dp), '1', 'Newton iterations')
    call result%write('energy_seconds', report%energy_s, 's', 'final solve: energy and gradient')
    call result%write('hessian_seconds', report%hessian_s, 's', 'final solve: Hessian assembly')
    call result%write('factor_seconds', report%factor_s, 's', 'final solve: factorisation and solves')
    call result%write('factorizations', real(report%factorizations, dp), '1', &
        'final solve: factorisations')
    call result%write('lu_fallbacks', real(report%lu_fallbacks, dp), '1', &
        'final solve: Cholesky failures refactored by pivoted LU')
    call result%write('analyses', real(report%analyses, dp), '1', 'final solve: symbolic analyses')
    call result%close()
    if (.not. free_boundary) call benchmark_export()

    block
        integer :: venv
        call get_environment_variable('MILCA_FB_PROFILE', status=venv)
        if (venv == 0) print '(a,13(i6,f8.3))', 'profile calls/s: map energy, toroidal, ' &
            //'bem terms, map hessian, toroidal hessian, bem hessian, factor, refactor, ' &
            //'newton factor, slide system, exact blocks, map steps, newton solve', &
            (fb_prof_calls(venv), fb_prof_seconds(venv), venv=1, 13)
    end block
contains

    subroutine fail(reason)
        ! A state that cannot be solved or read back is a result, not a crash:
        ! write converged = false with the reason and stop.
        character(*), intent(in) :: reason
        print '(2a)', 'nested_equilibrium: ', reason
        call result%create(run)
        call result%attribute('physical_model', 'fixed-boundary nested ideal MHD')
        call result%attribute('converged', 'false')
        call result%attribute('failure', reason)
        call result%close()
        stop
    end subroutine

    subroutine benchmark_export()
        ! Diagnostic for cross-code benchmarks (benchmarks/cemv-vs-vmec): with
        ! MILCA_NESTED_EXPORT=DIR the sizes, the force balance J X B - grad p
        ! on the measure of the flux map (milca_map_3d_force) are printed and
        ! axis.txt, surface_s0.25.txt, surface_s0.50.txt, summary.json are
        ! written to DIR in the layout of benchmarks/label-qa/compare.py (64
        ! cylindrical planes of one period, 128 reference angles).
        integer, parameter :: nout = 64, nray = 128
        character(512) :: dir
        type(flux_map_force_t) :: force
        real(dp) :: centre(2), y(2), xi(2), point(3), values(EQ_VALUES), level, orientation
        real(dp) :: rays(2, nray)
        integer :: cells(nray)
        integer :: length, u, j, q, l, local, cell
        call get_environment_variable('MILCA_NESTED_EXPORT', dir, length)
        if (length == 0) return
        if (chart /= 'cylindrical') then
            print '(a)', 'benchmark: the export samples cylindrical planes, chart = cylindrical only'
            return
        end if
        print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,f5.1,a,i0,a,i0)', 'benchmark: degree ', degree, &
            '  boundary_nodes ', boundary_nodes, '  map_ntor ', ntor_level, '  cells ', &
            map%cells, '  coefficients ', map%ndof, '  unknowns ', count(.not. map%fixed), &
            '  sigma', map%sigma, '  Newton (first ladder) ', newton_total, &
            '  Newton (final solve) ', report%iterations
        orientation = benchmark_field_orientation()
        call map_3d_force(map, seed, cemv_profiles, orientation, 10, force, local)
        if (local /= 0) call fail('benchmark: invalid map sample in force diagnostic')
        call flux_map_force_print(force, 'Lagrangian map')
        print '(a,5f10.3)', 'benchmark seconds: energy, Hessian, linear, symbolic, numeric ', &
            report%energy_s, report%hessian_s, report%factor_s, report%analysis_s, report%numeric_s
        print '(a,5i8)', 'benchmark counts: Hessians, factorizations, analyses, Krylov solves and steps ', &
            report%hessians, report%factorizations, report%analyses, report%krylov_solves, report%krylov_iterations
        call execute_command_line('mkdir -p '//trim(dir))
        centre = seed%centre
        open (newunit=u, file=trim(dir)//'/axis.txt', status='replace')
        write (u, '(a)') '# phi R Z'
        call locate(centre, cell, xi)
        if (cell == 0) error stop 'benchmark export: axis'
        do j = 0, nout - 1
            call map_3d_sample(map, seed, cell, xi, 2*pi*j/(nfp*nout), point, values, local)
            write (u, '(3es24.16)') point(3), point(1), point(2)
        end do
        close (u)
        do l = 1, 2
            level = 0.25_dp*l
            open (newunit=u, file=trim(dir)//merge('/surface_s0.25.txt', '/surface_s0.50.txt', &
                l == 1), status='replace')
            write (u, '(a)') '# phi angle R Z |B|'
            do q = 1, nray
                y = centre + sqrt(level)*seed%axes*[cos(2*pi*(q - 1)/nray), sin(2*pi*(q - 1)/nray)]
                call locate(y, cells(q), rays(:, q))
                if (cells(q) == 0) error stop 'benchmark export: surface point'
            end do
            do j = 0, nout - 1
                do q = 1, nray
                    call map_3d_sample(map, seed, cells(q), rays(:, q), 2*pi*j/(nfp*nout), &
                        point, values, local)
                    if (local /= 0) cycle
                    write (u, '(5es24.16)') point(3), 2*pi*(q - 1)/nray, point(1), point(2), &
                        norm2(values([EQ_BR, EQ_BZ, EQ_BPHI]))
                end do
            end do
            close (u)
        end do
        open (newunit=u, file=trim(dir)//'/summary.json', status='replace')
        ! the energy of the flux map and of the references: int (B^2/2 - mu0 p) dV
        write (u, '(a,es24.16,a,i0,a,i0,3a)') '{"energy_T2m3": ', magnetic - kinetic, &
            ', "unknowns": ', count(.not. map%fixed), ', "newton_steps": ', &
            max(newton_total, report%iterations), ', "converged": ', trim(merge('true ', 'false', converged)), '}'
        close (u)
    end subroutine

    real(dp) function benchmark_field_orientation() result(orientation)
        ! Match the frozen-in field independently of force balance.
        real(dp), parameter :: xi(2) = [0.3_dp, 0.2_dp], signs(2) = [1.0_dp, -1.0_dp]
        real(dp) :: b(3), j(3), f(3), point(3), values(EQ_VALUES), reference(3)
        real(dp) :: error(2), polarity(2), phi
        integer :: e, k, o, local

        error = 0
        do e = 1, map%cells, max(1, map%cells/50)
            do k = 0, 2
                phi = 2*pi*k/(nfp*3.7_dp)
                call map_3d_sample(map, seed, e, xi, phi, point, values, local)
                if (local /= 0) call fail('benchmark: cannot sample the frozen-in field')
                reference = values([EQ_BR, EQ_BPHI, EQ_BZ])
                do o = 1, 2
                    call map_3d_force_field(map, seed, cemv_profiles, signs(o), e, xi, phi, b, j, f, local)
                    if (local /= 0) call fail('benchmark: force diagnostic requires a cylindrical chart')
                    if (e == 1 .and. k == 0) polarity(o) = sign(1.0_dp, b(2)*reference(2))
                    error(o) = max(error(o), norm2(polarity(o)*b - reference)/norm2(reference))
                end do
            end do
        end do
        print '(a,2es11.3)', 'benchmark field errors for the two orientations ', error
        o = minloc(error, dim=1)
        if (.not. error(o) < 1e-10_dp) call fail('benchmark: inconsistent magnetic field')
        orientation = signs(o)
    end function benchmark_field_orientation

    subroutine perturb_start()
        ! Diagnostic for convergence basins (benchmarks/cemv-vs-vmec):
        ! MILCA_NESTED_START_SHIFT=d moves the interior of the start by
        ! d (1 - s) in R (the axis by d metres, the wall stays);
        ! MILCA_NESTED_START_ANGLE=a starts the boundary field-line angle with
        ! eta = a sin(theta) at the wall (the m = 1, n = 0 mode).
        character(64) :: value
        real(dp) :: shift, angle, s
        integer :: length, a, i
        call get_environment_variable('MILCA_NESTED_START_SHIFT', value, length)
        if (length > 0) then
            read (value, *) shift
            do a = 1, map%nodes
                if (map%boundary(a)) cycle
                s = sum(((map%reference(:, a) - seed%centre)/seed%axes)**2)
                map%x(1, 1, a) = map%x(1, 1, a) + shift*(1 - s)
            end do
        end if
        call get_environment_variable('MILCA_NESTED_START_ANGLE', value, length)
        if (length > 0 .and. map%nangle > 0) then
            read (value, *) angle
            do i = 1, size(map%angle_m)
                if (map%angle_m(i) == 1 .and. map%angle_n(i) == 0) map%angle(2*i) = angle
            end do
        end if
    end subroutine

    subroutine locate(y, cell, xi)
        ! cell and local coordinates of a reference point (0: not found)
        real(dp), intent(in) :: y(2)
        integer, intent(out) :: cell
        real(dp), intent(out) :: xi(2)
        integer :: e, local
        cell = 0
        do e = 1, map%cells
            call triangle_locate(map%reference(:, map%conn(:, e)), y, xi, local)
            if (local /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1e-10_dp) cycle
            cell = e
            return
        end do
    end subroutine

    subroutine cemv_profiles(s, tau, gp, p)
        ! This force diagnostic needs only the first derivatives in s.
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: iota_value, pressure_value, integral, iota_slope, pressure_slope

        call profile%evaluate(s, iota_value, pressure_value, integral, iota_slope, pressure_slope)
        tau = [phiedge/pi, 0.0_dp, 0.0_dp]
        gp = phiedge/(2*pi)*[iota_value, iota_slope, 0.0_dp]
        p = [pressure_value, pressure_slope, 0.0_dp]
    end subroutine

    subroutine mesh_ladder()
        ! Resolution ladder (docs/equilibrium-harness.md): each level starts from
        ! its own seed plus the correction of the solved coarser level and stops
        ! at the same absolute gradient relative to its uncorrected seed.
        steps = count(ladder_nodes > 0)
        do step = 1, steps + 1
            nodes = boundary_nodes; order = degree
            if (step <= steps) then
                nodes = ladder_nodes(step); order = ladder_degree(step)
            end if
            call ideal_seed_disk_mesh(mesh, nodes, order, status)
            if (status /= 0) error stop 'nested_equilibrium: disk mesh'
            call init_seed(status)
            if (status /= 0) call fail('the seed is not admitted (folded or invalid input)')
            reference = free_gradient()
            seed_x = map%x
            call perturb_start()
            if (step > 1) then
                call continue_fixed_level(coarse, coarse_seed)
            end if
            call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
                reference_gradient=reference, stationary=stationary, &
                model_agreement_damping=map%ntor > 0 .and. .not. prescribe_current)
            if (stationary .and. step <= steps .and. solve_status /= IDEAL_MAP_OK) &
                call fail('the coarse stationary solve did not converge')
            newton_total = newton_total + report%iterations
            coarse = map; coarse_seed = seed_x
        end do
    end subroutine

    subroutine continue_fixed_level(start_map, start_seed)
        type(ideal_map_t), intent(in) :: start_map
        real(dp), intent(in) :: start_seed(:, :, :)
        integer :: local

        if (stationary) then
            ! A branch start needs both the solved geometry and the field angle.
            call ideal_map_prolong(map, start_map, local)
        else
            call ideal_map_prolong_correction(map, start_map, start_seed, local)
        end if
        if (local /= IDEAL_MAP_OK) call fail('the resolution level transfer failed')
    end subroutine

    subroutine free_boundary_setup()
        ! coil groups of the mgrid file and the coupling surface around the guess
        type(mgrid_t) :: grid
        real(dp), allocatable :: coils_inside(:, :)
        integer :: local
        if (degree < 2 .or. degree > 3) &
            call fail('the free boundary takes degree 2 or 3')
        call mgrid_read(trim(mgrid_file), grid, local, message)
        if (local /= 0) call fail('mgrid: '//trim(message))
        if (grid%nextcur > extcur_max) &
            call fail('mgrid: more coil groups than extcur holds')
        if (grid%kp > 1 .and. mod(nfp, grid%nfp) /= 0) &
            call fail('mgrid: the field periods of the grid do not divide nfp')
        call mgrid_basis_init(coil_basis, grid%b, grid%nfp, grid%rmin, grid%rmax, &
            grid%zmin, grid%zmax, grid%unit_scale(), mgrid_degree, local)
        if (local /= 0) call fail('mgrid: invalid grid or mgrid_degree (1, 3 or 5)')
        print '(a,3i5,a,i0,a,i0)', 'mgrid ir jz kp', grid%ir, grid%jz, grid%kp, &
            ', nfp ', grid%nfp, ', coil groups ', grid%nextcur
        call mgrid_current_nodes(coil_basis, extcur(:grid%nextcur), coils_inside)
        call fb_nested_geometry(boundary_input, nfp, coupling_margin, [grid%rmin, &
            grid%rmax, grid%zmin, grid%zmax], coils_inside, fb_centre, fb_axes, &
            fb_curve, margin_used, local, message)
        if (local /= 0) call fail(trim(message))
        fit_curve = fb_curve
        if (split_coils) then
            ! coils split off: Gamma_c only has to enclose the plasma current
            block
                real(dp) :: centre(2), axes(2), margin_wide
                real(dp), allocatable :: none(:, :)
                allocate (none(2, 0))
                call fb_nested_geometry(boundary_input, nfp, coupling_margin, [grid%rmin, &
                    grid%rmax, grid%zmin, grid%zmax], none, centre, axes, fb_curve, &
                    margin_wide, local, message)
                if (local /= 0) call fail(trim(message))
                print '(a,f6.3)', 'coils split off: coupling margin', margin_wide
            end block
        end if
        if (use_bem .or. split_coils .or. split_env()) then
            ! the fit holds on the coil-free region inside fit_curve
            call fb_nested_coil_flux(coil_basis, extcur(:grid%nextcur), fb_centre, &
                fit_curve, coil_flux, bem_f_vac, bem_f_spread, local)
            if (local /= 0) call fail('the coil flux fit inside Gamma_c failed')
            print '(a,es10.3,a,es12.4,a,es9.2)', 'coil flux fit misfit ', &
                coil_flux%misfit, ', F_vac (T m) ', bem_f_vac, ', R B_phi spread ', &
                bem_f_spread
        end if
        print '(a,2f9.4,a,2f8.4,a,f6.3)', 'reference disc centre', fb_centre, &
            ' half-widths', fb_axes, ', coupling margin', margin_used
        filament%source = fb_centre
        seed%profile = profile; seed%flux = phiedge
        seed%centre = fb_centre; seed%axes = fb_axes
        if (.not. profile%valid() .or. .not. single_signed()) &
            call fail('the frozen field needs an iota profile without zeros')
        circulation = 0
        if (prescribe_current) circulation = mu0*current_target(1.0_dp)
    end subroutine

    subroutine fb_level(nodes, first)
        ! Free boundary at the mesh of nodes boundary nodes and the current
        ! toroidal modes: coupling-surface operators and traces, the start
        ! from the guess (first) or from the solved level, and the solve.
        integer, intent(in) :: nodes
        logical, intent(in) :: first
        type(ideal_fb_t) :: next
        integer :: local
        seed%profile = profile
        if (use_bem) then
            call set_bem_orders(nodes)
            call bem_level(nodes, next)
        else
            call shell_level(nodes, next)
        end if
        if (.not. next%f_vac > 0) call fail('the coil toroidal field is not positive;' &
            //' reverse extcur and the current for a negative toroidal flux')
        if (first) then
            call fb_nested_seed(next, boundary, fb_centre, fb_axes, local)
            if (local /= IDEAL_FB_OK) &
                call fail('the free-boundary seed is not admitted')
        else
            call fb_nested_transfer(next, fb, local)
            if (local /= IDEAL_FB_OK) call fail('the resolution level transfer failed')
            if (next%na0 > 0 .and. fb%na0 > 0) next%plasma%angle = fb%plasma%angle
        end if
        if (first .and. state_cached(next)) then
            call install_traces(next)
            call fb_nested_settle(next, local)
        else if (first) then
            if (prescribe_current) call seed_current(next)
            call interior_equilibrium(next)
            if (.not. prescribe_current) then
                ! prescribed iota: start from the circulation of the seeded plasma
                map = next%plasma
                circulation = mu0*edge_current(local)
                if (local /= 0) call fail('the edge current readback met a folded map')
                call install_traces(next)
            end if
            call fb_nested_settle(next, local)
            if (local /= IDEAL_FB_OK) &
                call fail('the free-boundary seed is not admitted')
            call state_store(next)
            if (split_env() .and. .not. use_bem) call split_check(next)
        end if
        fb = next
        call fb_solve(max(tolerance, level_tolerance))
        block
            character(len=8) :: env
            integer :: n_env, st_env
            call get_environment_variable('MILCA_FB_FIRST_ONLY', env, n_env, st_env)
            if (first .and. st_env == 0 .and. n_env > 0) stop
        end block
    end subroutine

    logical function state_cached(next) result(found)
        type(ideal_fb_t), intent(inout) :: next
        character(len=512) :: path
        integer :: n_env, st_env, unit, io
        found = .false.
        call get_environment_variable('MILCA_FB_STATE', path, n_env, st_env)
        if (st_env /= 0 .or. n_env == 0) return
        inquire (file=trim(path), exist=found)
        if (.not. found) return
        open (newunit=unit, file=trim(path), form='unformatted', access='stream', &
            status='old', action='read', iostat=io)
        read (unit, iostat=io) next%plasma%x, next%plasma%b, next%plasma%p, circulation
        close (unit)
        found = io == 0
        if (found) print '(a)', 'state cache: loaded '//trim(path)
    end function

    subroutine state_store(next)
        type(ideal_fb_t), intent(in) :: next
        character(len=512) :: path
        integer :: n_env, st_env, unit
        logical :: there
        call get_environment_variable('MILCA_FB_STATE', path, n_env, st_env)
        if (st_env /= 0 .or. n_env == 0) return
        inquire (file=trim(path), exist=there)
        if (there) return
        open (newunit=unit, file=trim(path), form='unformatted', access='stream', &
            status='new', action='write')
        write (unit) next%plasma%x, next%plasma%b, next%plasma%p, circulation
        close (unit)
        print '(a)', 'state cache: stored '//trim(path)
    end subroutine

    logical function gauge_env()
        ! position gauge with the field angle; the transversal gauge of the
        ! seed interior is a testbed option only (it loses the edge current)
        integer :: n_env, st_env
        character(len=8) :: env
        call get_environment_variable('MILCA_TRANSVERSAL', env, n_env, st_env)
        gauge_env = .not. (st_env == 0 .and. n_env > 0)
    end function

    logical function split_env()
        integer :: n_env, st_env
        character(len=8) :: env
        call get_environment_variable('MILCA_SPLIT_CHECK', env, n_env, st_env)
        split_env = st_env == 0 .and. n_env > 0
    end function

    subroutine split_check(target)
        ! Reduced gradient (shell map, Phi~, multipliers exact at the plasma
        ! map) with the coil traces on Gamma_c against the coil split.
        type(ideal_fb_t), intent(in) :: target
        type(ideal_fb_t) :: one, two
        real(dp), allocatable :: g1(:), g2(:), d(:, :), l(:, :)
        real(dp) :: lag1, lag2
        integer :: local, np, k
        allocate (g1(target%n), g2(target%n))
        one = target
        two = target
        d = circulation*filament_drive; l = circulation*filament_load
        call ideal_fb_set_traces(two, d, l, target%f_vac, local)
        call ideal_fb_set_coil_split(two, coil_flux, local)
        if (local /= IDEAL_FB_OK) call fail('coil split set-up failed')
        do k = 1, 2
            call ideal_fb_solve_shell_map(one, local); call ideal_fb_solve_vacuum(one, local)
            call ideal_fb_solve_multipliers(one, local)
            call ideal_fb_solve_shell_map(two, local); call ideal_fb_solve_vacuum(two, local)
            call ideal_fb_solve_multipliers(two, local)
        end do
        call ideal_fb_gradient(one, g1, local, lag1)
        call ideal_fb_gradient(two, g2, local, lag2)
        np = target%np
        print '(a,4es12.4)', 'reduced plasma |g1| |g2| |g1-g2| cos:', norm2(g1(:np)), &
            norm2(g2(:np)), norm2(g1(:np) - g2(:np)), &
            dot_product(g1(:np), g2(:np))/(norm2(g1(:np))*norm2(g2(:np)))
        print '(a,2es12.4)', 'reduced rest |g1| |g2|:', norm2(g1(np+1:)), norm2(g2(np+1:))
        print '(a,2es22.14)', 'reduced L1 L2:', lag1, lag2
        call hvp_check(one, 'traces')
        call hvp_check(two, 'split ')
        call reduced_fd(one, 'traces')
        call reduced_fd(two, 'split ')
        stop
    end subroutine

    subroutine hvp_check(state, label)
        ! H v against the central difference of the gradient along v, per block.
        type(ideal_fb_t), intent(in) :: state
        character(*), intent(in) :: label
        type(ideal_fb_t) :: work
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: val(:), v(:), hv(:), gp(:), gm(:), fd(:)
        real(dp) :: h
        integer :: i, p, local, np
        np = state%np
        work = state
        call ideal_fb_hessian(work, rp, ci, val)
        allocate (v(state%n), hv(state%n), gp(state%n), gm(state%n))
        v = [(sin(0.37_dp*i), i = 1, state%n)]
        v(np + 1:) = 0
        hv = 0
        do i = 1, state%n
            do p = rp(i), rp(i + 1) - 1
                hv(i) = hv(i) + val(p)*v(ci(p))
            end do
        end do
        h = 1.0e-6_dp
        call ideal_fb_step(work, v, h); call ideal_fb_gradient(work, gp, local)
        call ideal_fb_step(work, v, -2*h); call ideal_fb_gradient(work, gm, local)
        fd = (gp - gm)/(2*h)
        print '(a,a,2es12.4)', label, ' HVP plasma-direction: rel err plasma rows, rest rows', &
            norm2(fd(:np) - hv(:np))/norm2(fd(:np)), norm2(fd(np+1:) - hv(np+1:))/max(norm2(fd(np+1:)), 1e-300_dp)
    end subroutine

    subroutine reduced_fd(state, label)
        ! dE_red/dv by central differences (blocks re-solved) against g.v.
        type(ideal_fb_t), intent(in) :: state
        character(*), intent(in) :: label
        type(ideal_fb_t) :: work
        real(dp), allocatable :: v(:), g(:)
        real(dp) :: ep, em, e0, h
        integer :: i, local, np
        np = state%np
        allocate (v(state%n), g(state%n))
        v = 0
        v(:np) = [(sin(0.71_dp*i), i = 1, np)]
        work = state
        call ideal_fb_gradient(work, g, local, e0)
        do i = 4, 7
            h = 10.0_dp**(-i)
            work = state; call ideal_fb_step(work, v, h); call settle(work)
            call ideal_fb_gradient(work, g, local, ep)
            work = state; call ideal_fb_step(work, v, -h); call settle(work)
            call ideal_fb_gradient(work, g, local, em)
            work = state
            call ideal_fb_gradient(work, g, local, e0)
            print '(a,a,es9.1,2es16.8)', label, ' reduced dE/dv h, FD, g.v:', h, (ep - em)/(2*h), &
                dot_product(g(:np), v(:np))
        end do
    end subroutine

    subroutine settle(w)
        type(ideal_fb_t), intent(inout) :: w
        integer :: k, lo
        do k = 1, 2
            call ideal_fb_solve_shell_map(w, lo); call ideal_fb_solve_vacuum(w, lo)
            call ideal_fb_solve_multipliers(w, lo)
        end do
    end subroutine

    logical function talk_env()
        integer :: n_env, st_env
        character(len=8) :: env
        call get_environment_variable('MILCA_FB_TALK', env, n_env, st_env)
        talk_env = st_env == 0 .and. n_env > 0
    end function

    subroutine shell_level(nodes, next)
        ! shell to Gamma_c: mesh, per-n exterior and the coil traces
        integer, intent(in) :: nodes
        type(ideal_fb_t), intent(out) :: next
        real(dp), allocatable :: chain(:, :), d(:, :, :), l(:, :, :)
        real(dp) :: closure(1), f(1), spread(1)
        integer :: local
        call fb_nested_mesh(mesh, chain, fb_centre, fb_axes, fb_curve, nodes, degree, &
            local)
        if (local /= 0) call fail('the free-boundary mesh failed')
        seed%profile = profile
        call ideal_fb_init(next, mesh, chain, cylinder, seed, nfp, ntor_level, 0.0_dp, &
            circulation, local, nphi=max(1, 4*ntor_level + 1))
        if (local /= IDEAL_FB_OK) call fail('the free-boundary set-up failed')
        if (split_coils) then
            call ideal_fb_set_coil_split(next, coil_flux, local)
            if (local /= IDEAL_FB_OK) call fail('the coil split set-up failed')
            if (.not. any(abs(boundary(2:3, :, :)) > 0)) call ideal_fb_fix_vertical(next)
        else
            call ideal_fb_circuits_build(next, coil_basis, circuits, local)
            if (local /= IDEAL_FB_OK) &
                call fail('the coil field grid does not cover the coupling surface')
        end if
        allocate (d(next%nt, next%nmode, 1), l(next%nt, next%nmode, 1))
        call ideal_fb_drive_traces(next, filament, [1.0_dp], d, l, closure, f, spread, &
            local)
        if (local /= IDEAL_FB_OK) call fail('the plasma-current traces failed')
        filament_drive = d(:, :, 1); filament_load = l(:, :, 1)
        ! bound of the closure of the summed coil field: a coil inside Gamma_c
        ! (a toroidal current it encloses) shows as an n = 0 closure
        coil_closure = 0
        if (.not. split_coils) coil_closure = sum(abs(extcur(:circuits%ncircuit))*circuits%closure)
        call install_traces(next)
        print '(a,es12.4,a,es10.2,a,i0,a,i0)', 'F_vac (T m)', next%f_vac, &
            ', coil closure', coil_closure, ', unknowns ', next%n, ', modes ', &
            ntor_level
    end subroutine

    subroutine bem_level(nodes, next)
        ! sheet current on Gamma_p: the plasma mesh alone; an up-down
        ! symmetric deck (no rbs, zbc) holds the vertical position like
        ! VMEC with lasym = F
        integer, intent(in) :: nodes
        type(ideal_fb_t), intent(out) :: next
        integer :: local
        call fb_nested_plasma_mesh(mesh, fb_centre, fb_axes, nodes, degree, local)
        if (local /= 0) call fail('the plasma mesh failed')
        call ideal_fb_init_bem(next, mesh, cylinder, seed, nfp, bem_f_vac, &
            circulation, coil_flux, degree - 1, local, nphi=1)
        if (local /= IDEAL_FB_OK) call fail('the free-boundary set-up failed')
        if (.not. any(abs(boundary(2:3, :, :)) > 0)) call ideal_fb_fix_vertical(next)
        if (angle_mpol > 0 .and. ntor == 0) then
            call ideal_fb_enable_angle(next, angle_mpol, phiedge, fb_centre, fb_axes, local)
            if (local /= IDEAL_FB_OK) call fail('the field angle of the free boundary failed')
        end if
        print '(a,es12.4,a,i0)', 'F_vac (T m)', next%f_vac, ', unknowns ', next%n
    end subroutine

    subroutine interior_equilibrium(target)
        ! The seed's interior in fixed-boundary equilibrium (position gauge,
        ! Gamma_p held), with a prescribed current matched there by the
        ! fixed-boundary Newton on iota: the free-boundary Newton then starts
        ! with the force imbalance on Gamma_p only and with the plasma's edge
        ! current equal to the vacuum circulation, instead of moving the whole
        ! plasma to relax the Winslow interior and a current sheet.
        type(ideal_fb_t), intent(inout) :: target
        integer :: local
        map = target%plasma
        ! transversal gauge with omega: the relabelling along the surfaces is a
        ! gauge instead of a soft material slide (18 instead of 43 free-boundary
        ! iterations at 48 boundary nodes; both gauges converge to the same
        ! volume under refinement)
        map%position_gauge = gauge_env()
        block
            character(len=16) :: env
            integer :: n_env, st_env
            call get_environment_variable('MILCA_GAUGE_WEIGHT', env, n_env, st_env)
            if (st_env == 0 .and. n_env > 0) read (env, *) map%gauge_weight
        end block
        seed%profile = profile
        call ideal_map_freeze(map, seed, local)
        if (local /= 0) call fail('the seed folded when frozen for its interior')
        reference = free_gradient()
        call ideal_map_solve(map, tolerance, max_iterations, report, local, &
            reference_gradient=reference)
        print '(a,i0,a,i0)', 'fixed-boundary interior of the seed: status ', local, &
            ', iterations ', report%iterations
        fb_solving = .false.
        if (prescribe_current) then
            solve_tolerance = tolerance
            call match_current(current_error)
        end if
        fb_solving = .true.
        target%plasma%x = map%x
        seed%profile = profile
        call ideal_map_freeze(target%plasma, seed, local)
        if (local /= 0) call fail('the seed folded when frozen with the matched iota')
    end subroutine

    subroutine seed_current(target)
        ! Iota of the seed matched to the prescribed current on the seed's
        ! geometry, where I(s_j) = A_j c_j + C_j is exact: the first solve then
        ! starts with the plasma current equal to the vacuum circulation.
        type(ideal_fb_t), intent(inout) :: target
        type(map_3d_loop_t), allocatable :: loops(:)
        real(dp), allocatable :: s(:), c(:), now(:), shifted(:)
        real(dp) :: unused1, unused2
        integer :: j, local, n
        n = iota_level
        s = iota_surfaces(n)
        allocate (loops(n), c(n), now(n), shifted(n))
        map = target%plasma
        do j = 1, n
            call map_3d_loop_setup(map, s(j), loops(j), local, elliptic=free_boundary)
            if (local /= 0) call fail('a current loop met a curved cell edge')
            call profile%evaluate(s(j), c(j), unused1, unused2)
        end do
        call set_iota(c); seed%profile = profile
        call measure(loops, now)
        call set_iota(c + 1); seed%profile = profile
        call measure(loops, shifted)
        do j = 1, n
            shifted(j) = c(j) + (current_target(s(j)) - now(j))/(shifted(j) - now(j))
        end do
        call set_iota(shifted)
        ! a seed whose current needs iota of both signs keeps its iota
        if (.not. single_signed()) call set_iota(c)
        seed%profile = profile
        call ideal_map_freeze(target%plasma, seed, local)
        if (local /= 0) call fail('the seed folded when frozen with the matched iota')
    end subroutine

    subroutine install_traces(target)
        ! vacuum traces: the coil groups at extcur plus the plasma current
        type(ideal_fb_t), intent(inout) :: target
        real(dp), allocatable :: d(:, :), l(:, :)
        real(dp) :: f
        integer :: k, local
        if (use_bem) then
            target%circulation = circulation
            call ideal_fb_bem_settle(target, local)
            if (local /= IDEAL_FB_OK) call fail('the sheet current is not finite')
            return
        end if
        d = circulation*filament_drive; l = circulation*filament_load; f = 0
        if (split_coils) f = bem_f_vac
        if (split_coils) then
            target%circulation = circulation
            call ideal_fb_set_traces(target, d, l, f, local)
            if (local /= IDEAL_FB_OK) call fail('the vacuum traces are not finite')
            return
        end if
        do k = 1, circuits%ncircuit
            d = d + extcur(k)*circuits%drive(:, :, k)
            l = l + extcur(k)*circuits%load(:, :, k)
            f = f + extcur(k)*circuits%f_vac(k)
        end do
        target%circulation = circulation
        call ideal_fb_set_traces(target, d, l, f, local)
        if (local /= IDEAL_FB_OK) call fail('the vacuum traces are not finite')
    end subroutine

    subroutine fb_solve(target)
        ! Newton on the free-boundary residual to the relative residual target.
        ! A solve that stalls has driven its Levenberg damping up; it is
        ! restarted from its last state with the damping reset while each
        ! restart lowers the residual by at least a quarter.
        real(dp), intent(in) :: target
        real(dp) :: last
        integer :: restart
        last = huge(1.0_dp)
        do restart = 0, 20
            call ideal_fb_solve(fb, 0.1_dp*target, max_iterations, fb_report, &
                solve_status, verbose=talk_env(), residual_tolerance=target)
            if (solve_status == IDEAL_FB_OK .or. solve_status == IDEAL_FB_FOLDED) exit
            if (.not. fb_report%residual < 0.75_dp*last) exit
            last = fb_report%residual
        end do
        map = fb%plasma
        print '(a,i0,a,i0,a,es10.3)', 'free-boundary solve: status ', solve_status, &
            ', iterations ', fb_report%iterations, ', residual ', fb_report%residual
    end subroutine

    subroutine match_circulation(error)
        ! Prescribed iota: the vacuum circulation c equals mu0 times the edge
        ! current of the solved plasma, so that no current sheet remains on
        ! Gamma_p. error: the relative change of c in the last update.
        ! The edge current responds to c through the plasma position with a
        ! gain near -1, so the fixed point c = mu0 I(c) is found by secants.
        real(dp), intent(out) :: error
        real(dp) :: defect, previous_c, previous_defect, next_c
        integer :: k, local
        error = 0
        previous_c = 0; previous_defect = 0
        do k = 1, max_current_iterations
            defect = mu0*edge_current(local) - circulation
            if (local /= 0) call fail('the edge current readback met a folded map')
            error = abs(defect)/max(abs(circulation + defect), tiny(1.0_dp))
            print '(a,i3,a,es12.4,a,es10.3)', 'circulation iteration ', k, ': ', &
                circulation, ' T m, relative defect ', error
            current_iterations = k
            if (error <= current_tolerance) exit
            next_c = circulation + defect
            if (k > 1 .and. abs(defect - previous_defect) > 0) next_c = circulation &
                - defect*(circulation - previous_c)/(defect - previous_defect)
            previous_c = circulation; previous_defect = defect
            circulation = next_c
            call install_traces(fb)
            call fb_solve(tolerance)
        end do
    end subroutine

    subroutine fb_current_jacobian(c, loops, measured, jacobian)
        ! current_jacobian with the free-boundary response: dg/dc_k on all
        ! unknowns (plasma rows only) and dx_k = -K^-1 dg/dc_k from the
        ! slide-gauged KKT matrix K of the solved state (ideal_fb_sensitivity).
        real(dp), intent(in) :: c(:), measured(:)
        type(map_3d_loop_t), intent(in) :: loops(:)
        real(dp), intent(out) :: jacobian(:, :)
        real(dp), allocatable :: rhs(:, :), dx(:, :), plus(:), minus(:), x0(:, :, :)
        real(dp), allocatable :: step(:, :, :)
        real(dp) :: ip(size(c)), im(size(c)), e(size(c)), h, length
        integer :: k, local, n
        n = size(c)
        allocate (rhs(fb%n, n), dx(fb%n, n), plus(fb%n), minus(fb%n))
        call set_iota(c + 1)
        seed%profile = profile
        call measure(loops, ip)
        jacobian = 0
        do k = 1, n
            jacobian(k, k) = ip(k) - measured(k)
        end do
        do k = 1, n
            e = 0; e(k) = 1
            call fb_frozen_gradient(c + e, plus)
            call fb_frozen_gradient(c - e, minus)
            rhs(:, k) = -(plus - minus)/2
        end do
        call fb_frozen_gradient(c, plus)
        call ideal_fb_sensitivity(fb, rhs, dx, local)
        if (local /= IDEAL_FB_OK) &
            call fail('the free-boundary sensitivity solve failed')
        map = fb%plasma
        x0 = map%x
        length = maxval(maxval(map%reference, dim=2) - minval(map%reference, dim=2))
        allocate (step(2, fb%nmode, fb%plasma%nodes))
        do k = 1, n
            step = reshape(dx(:fb%np, k), shape(step))
            h = 1e-6_dp*length/max(maxval(abs(step)), tiny(1.0_dp))
            map%x(1:2, :, :) = x0(1:2, :, :) + h*step
            call measure(loops, ip)
            map%x(1:2, :, :) = x0(1:2, :, :) - h*step
            call measure(loops, im)
            jacobian(:, k) = jacobian(:, k) + (ip - im)/(2*h)
        end do
        map%x = x0
    end subroutine

    subroutine fb_frozen_gradient(c, gradient)
        ! free-boundary residual with the plasma frozen at the iota values c
        real(dp), intent(in) :: c(:)
        real(dp), intent(out) :: gradient(:)
        integer :: local
        call set_iota(c)
        seed%profile = profile
        call ideal_map_freeze(fb%plasma, seed, local)
        if (local /= 0) call fail('the map folded when frozen in the current Jacobian')
        call ideal_fb_gradient(fb, gradient, local)
        if (local /= IDEAL_FB_OK) &
            call fail('the free-boundary residual met a folded map')
    end subroutine


    logical function single_signed()
        ! the seed's admission rule: iota without zeros on [0, 1]
        real(dp) :: value, first, unused1, unused2
        integer :: a
        call profile%evaluate(0.0_dp, first, unused1, unused2)
        single_signed = abs(first) > tiny(1.0_dp)
        do a = 1, 128
            call profile%evaluate(a/128.0_dp, value, unused1, unused2)
            single_signed = single_signed .and. value*first > 0
        end do
    end function

    real(dp) function min_abs_iota()
        ! smallest |iota| on the grid of the seed's iota admission
        real(dp) :: value, unused1, unused2
        integer :: a
        min_abs_iota = huge(1.0_dp)
        do a = 0, 128
            call profile%evaluate(a/128.0_dp, value, unused1, unused2)
            min_abs_iota = min(min_abs_iota, abs(value))
        end do
    end function

    subroutine raise_mesh(confirmed)
        ! mesh levels at the current toroidal modes until mesh_error is below
        ! mesh_tolerance; confirmed when the first refinement already is
        logical, intent(out) :: confirmed
        integer :: k, steps
        real(dp) :: last
        confirmed = .false.
        steps = 0; last = huge(1.0_dp)
        do k = 1, size(mesh_nodes)
            if (mesh_nodes(k) <= boundary_nodes) cycle
            call next_level(mesh_nodes(k), ntor_level)
            call refine_folded()
            mesh_error = change()
            steps = steps + 1
            print '(a,i4,a,es10.3)', 'boundary nodes ', boundary_nodes, ': mesh error ', &
                mesh_error
            if (iota_vanishes .or. mesh_error <= mesh_tolerance) exit
            ! a change that does not shrink is limited by another resolution
            if (mesh_error >= last) exit
            last = mesh_error
        end do
        ! a ladder at its finest level has nothing left to confirm
        confirmed = steps == 0 .or. (steps == 1 .and. mesh_error <= mesh_tolerance)
    end subroutine

    subroutine raise_modes(confirmed)
        ! toroidal levels + 2 at the current mesh until fourier_error is below
        ! fourier_tolerance; confirmed when the first one already is
        logical, intent(out) :: confirmed
        integer :: steps
        real(dp) :: last
        confirmed = .false.
        steps = 0; last = huge(1.0_dp)
        do while (ntor_level + 2 <= nmax)
            call next_level(boundary_nodes, ntor_level + 2)
            call refine_folded()
            fourier_error = change()
            steps = steps + 1
            print '(a,i3,a,es10.3)', 'toroidal modes ', ntor_level, ': fourier error ', &
                fourier_error
            if (iota_vanishes .or. fourier_error <= fourier_tolerance) exit
            ! a change that does not shrink is limited by another resolution
            if (fourier_error >= last) exit
            last = fourier_error
        end do
        ! a ladder at its finest level has nothing left to confirm
        confirmed = steps == 0 .or. (steps == 1 .and. fourier_error <= fourier_tolerance)
    end subroutine

    subroutine raise_iota(confirmed)
        ! iota polynomials of 2 more nodes on the solved map until iota_error is
        ! below iota_tolerance; confirmed when the first one already is
        logical, intent(out) :: confirmed
        integer :: steps
        real(dp) :: last
        confirmed = .false.
        steps = 0; last = huge(1.0_dp)
        do while (iota_level + 2 <= iota_nodes_max)
            previous_observed = observed
            iota_level = iota_level + 2
            call finish_level()
            call refine_folded()
            iota_error = change()
            steps = steps + 1
            print '(a,i3,a,es10.3)', 'iota nodes ', iota_level, ': iota error ', iota_error
            if (iota_vanishes .or. iota_error <= iota_tolerance) exit
            ! a change that does not shrink is limited by another resolution
            if (iota_error >= last) exit
            last = iota_error
        end do
        ! a ladder at its finest level has nothing left to confirm
        confirmed = steps == 0 .or. (steps == 1 .and. iota_error <= iota_tolerance)
    end subroutine

    subroutine next_level(nodes_next, modes_next)
        ! seed at the new resolution plus the correction of the solved level
        integer, intent(in) :: nodes_next, modes_next
        logical :: continued
        ! a folded level is neither a reference for the change nor a start
        continued = .not. level_folded
        if (continued) then
            previous_map = map; previous_observed = observed
            if (allocated(seed_x)) previous_seed = seed_x
        end if
        boundary_nodes = nodes_next; ntor_level = modes_next
        call set_toroidal_level(ntor_level)
        if (free_boundary) then
            call fb_level(boundary_nodes, .false.)
            call finish_level()
            return
        end if
        call ideal_seed_disk_mesh(mesh, boundary_nodes, degree, status)
        if (status /= 0) error stop 'nested_equilibrium: disk mesh'
        call init_seed(status)
        if (status /= 0) call fail('the seed of a resolution level is not admitted')
        reference = free_gradient()
        seed_x = map%x
        if (continued) then
            call continue_fixed_level(previous_map, previous_seed)
        end if
        call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
            reference_gradient=reference, stationary=stationary, &
            model_agreement_damping=map%ntor > 0 .and. .not. prescribe_current)
        call finish_level()
        ! The energy sees the map only at its quadrature points, and the
        ! continued start can reach a stationary point that folds between
        ! them (ideal_map_prolong_correction); the level is solved again from
        ! its injective Winslow seed.
        if (continued .and. level_folded) then
            print '(a,i4,a,i3,a)', 'boundary nodes ', boundary_nodes, ', toroidal modes ', &
                ntor_level, ': the continued map folds, solving from the seed'
            call init_seed(status)
            if (status /= 0) call fail('the seed of a resolution level is not admitted')
            reference = free_gradient()
            seed_x = map%x
            call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
                reference_gradient=reference, stationary=stationary, &
                model_agreement_damping=map%ntor > 0 .and. .not. prescribe_current)
            call finish_level()
        end if
    end subroutine

    subroutine refine_folded()
        ! A level whose map folds at the volume-integral points does not
        ! resolve the map. A fold on a plane of the energy's toroidal
        ! quadrature lies between in-plane points and asks for a finer mesh;
        ! a fold only between those planes asks for more toroidal modes,
        ! whose 4 ntor + 1 planes follow. Each such level starts from its seed.
        logical :: in_plane
        do while (level_folded)
            in_plane = folds_on_energy_planes()
            print '(a,i4,a,i3,a,a)', 'boundary nodes ', boundary_nodes, ', toroidal modes ', &
                ntor_level, ': the solved map folds ', &
                merge('on a plane of the energy quadrature        ', &
                'between the planes of the energy quadrature', in_plane)
            if (in_plane) then
                if (.not. automatic_mesh .or. boundary_nodes >= mesh_nodes(size(mesh_nodes))) &
                    call fail('the solved map folds between the in-plane quadrature points' &
                    //' at the finest admitted mesh')
                call next_level(minval(mesh_nodes, mask=mesh_nodes > boundary_nodes), ntor_level)
            else
                if (.not. automatic_ntor .or. ntor_level + 2 > nmax) &
                    call fail('the solved map folds between the planes of the energy' &
                    //' quadrature at the most admitted toroidal modes')
                call next_level(boundary_nodes, ntor_level + 2)
            end if
        end do
    end subroutine

    logical function folds_on_energy_planes()
        ! a fold at the volume-integral points of a plane of the energy's
        ! toroidal quadrature
        integer, parameter :: nq = 8
        real(dp) :: x(nq), w(nq), xi(2), jac, magnetic, pressure, radius, along
        integer :: e, i, j, k, local
        call gauss_legendre(nq, x, w)
        x = (x + 1)/2
        folds_on_energy_planes = .true.
        do e = 1, map%cells
            do k = 1, map%nphi
                do j = 1, nq
                    do i = 1, nq
                        xi = [x(i)*(1 - x(j)), x(j)]
                        call chart_point(e, xi, map%phi(k), [0.0_dp, 0.0_dp], jac, magnetic, &
                            pressure, radius, along, local)
                        if (local /= 0) return
                    end do
                end do
            end do
        end do
        folds_on_energy_planes = .false.
    end function

    subroutine finish_level()
        ! the prescribed current and the observables of the solved level
        current_error = 0; current_iterations = 0; solve_tolerance = tolerance
        if (prescribe_current) then
            call match_current(current_error)
        else if (free_boundary) then
            call match_circulation(current_error)
        end if
        call observe(observed)
    end subroutine

    real(dp) function change()
        ! wb relative, raxis/aminor_p and iota on axis and edge
        change = max(abs(observed(1)/previous_observed(1) - 1), &
            maxval(abs(observed(2:) - previous_observed(2:))))
    end function

    subroutine set_bem_orders(nodes)
        !! Quadrature of the plasma-boundary single layer by mesh: at 16
        !! boundary nodes near/far Gauss orders 4/3 and 8 graded radial points
        !! leave the exact Solov'ev free boundary unchanged against 8/6/20
        !! (3.50e-4 against 3.48e-4 in volume) at 60 % of the time; finer
        !! meshes take the higher orders so the quadrature stays below the
        !! plasma discretisation.
        integer, intent(in) :: nodes
        character(16) :: env
        integer :: venv
        if (nodes > 0 .and. nodes <= 24) then
            sl_near_order = 4; sl_far_order = 3; sl_radial_n = 8
        else if (nodes > 0 .and. nodes <= 48) then
            sl_near_order = 6; sl_far_order = 4; sl_radial_n = 12
        else
            sl_near_order = 8; sl_far_order = 6; sl_radial_n = 20
        end if
        call get_environment_variable('MILCA_BEM_ORDERS', env, status=venv)
        if (venv == 0) read (env, *) sl_near_order, sl_far_order, sl_radial_n
        call get_environment_variable('MILCA_STEP_CAP', env, status=venv)
        if (venv == 0) read (env, *) step_cap_fraction
    end subroutine set_bem_orders

    subroutine reseed_loop(nodes, steps)
        !! Seeds again from the solved boundary at this mesh until wb and the
        !! axis change by less than reseed_tolerance between seeds (at most
        !! steps seeds).
        integer, intent(in) :: nodes, steps
        call observe(reseed_seen)
        level_tolerance = seed_tolerance
        do reseed_step = 1, steps
            call reseed_from_solution(reseed_move)
            if (reseed_move < reseed_tolerance) exit
            call fb_level(nodes, .true.)
            reseed_last = reseed_seen
            call observe(reseed_seen)
            ! wb relative, the axis relative to the minor radius
            reseed_change = max(abs(reseed_seen(1)/reseed_last(1) - 1), &
                                abs(reseed_seen(2) - reseed_last(2)))
            print '(a,i2,a,es10.3,a,es10.3)', 'reseed', reseed_step, &
                ': boundary moved', reseed_move, ', wb and axis changed', reseed_change
            if (reseed_change < reseed_tolerance) exit
        end do
        level_tolerance = 0
        call fb_solve(tolerance)
    end subroutine reseed_loop

    subroutine reseed_from_solution(move)
        !! The solved plasma boundary (its boundary nodes) as the new seed:
        !! least-squares Fourier series of R and Z in the geometric angle about
        !! the nodes' centroid (cosines of R, sines of Z; all four when the
        !! deck is up-down asymmetric), at most reseed_modes and half the
        !! nodes. move: the largest change of a node's seed position against
        !! the previous seed, relative to the minor radius.
        real(dp), intent(out) :: move
        real(dp), allocatable :: x(:, :), th(:), a(:, :), b(:, :), sv(:), work(:), old(:, :)
        real(dp) :: c(2), radius, value
        integer :: nb, i, m, mm, ncol, rank, info, lwork, a_node
        logical :: asym
        interface
            subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
                import :: dp
                integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                real(dp), intent(out) :: s(*), work(*)
                real(dp), intent(in) :: rcond
                integer, intent(out) :: rank, info
            end subroutine dgelss
        end interface
        nb = count(fb%plasma%boundary)
        allocate (x(2, nb), th(nb))
        i = 0
        do a_node = 1, fb%plasma%nodes
            if (.not. fb%plasma%boundary(a_node)) cycle
            i = i + 1
            x(:, i) = fb%plasma%x(1:2, 1, a_node)
        end do
        c = sum(x, 2)/nb
        th = atan2(x(2, :) - c(2), x(1, :) - c(1))
        radius = sqrt(maxval(sum((x - spread(c, 2, nb))**2, 1)))
        asym = any(abs(boundary_input(2:3, :, :)) > 0)
        mm = min(reseed_modes, nb/2 - 1)
        ncol = mm + 1
        if (asym) ncol = 2*mm + 1
        allocate (a(nb, ncol), b(max(nb, ncol), 2), sv(ncol))
        do i = 1, nb
            a(i, 1:mm + 1) = cos([(m*th(i), m=0, mm)])
            if (asym) a(i, mm + 2:) = sin([(m*th(i), m=1, mm)])
        end do
        b(:nb, 1) = x(1, :)
        b(:nb, 2) = x(2, :)
        ! R in cosines (and sines), Z in sines (and cosines): two fits
        allocate (old(2, nb))
        do i = 1, nb
            old(:, i) = seed_point(th(i))
        end do
        lwork = -1
        allocate (work(1))
        call dgelss(nb, ncol, 1, a, nb, b(:, 1), max(nb, ncol), sv, -1.0_dp, rank, work, &
                    lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))
        block
            real(dp), allocatable :: ar(:, :), az(:, :), bz(:)
            ar = a
            call dgelss(nb, ncol, 1, ar, nb, b(:, 1), max(nb, ncol), sv, -1.0_dp, rank, &
                        work, lwork, info)
            ! Z: sines first (m = 1..mm), then cosines when asymmetric
            allocate (az(nb, ncol))
            do i = 1, nb
                az(i, 1) = 0
                az(i, 1:mm) = sin([(m*th(i), m=1, mm)])
                if (asym) az(i, mm + 1:) = cos([(m*th(i), m=0, mm)])
            end do
            bz = b(:, 2)
            call dgelss(nb, merge(ncol, mm, asym), 1, az, nb, bz, max(nb, ncol), sv, &
                        -1.0_dp, rank, work, lwork, info)
            deallocate (boundary_input)
            mpol = mm
            allocate (boundary_input(4, -ntor:ntor, 0:mpol), source=0.0_dp)
            boundary_input(1, 0, 0:mm) = b(1:mm + 1, 1)
            boundary_input(4, 0, 1:mm) = bz(1:mm)
            if (asym) then
                boundary_input(2, 0, 1:mm) = b(mm + 2:2*mm + 1, 1)
                boundary_input(3, 0, 0:mm) = bz(mm + 1:2*mm + 1)
            end if
        end block
        call set_toroidal_level(ntor_level)
        move = 0
        do i = 1, nb
            value = norm2(seed_point(th(i)) - old(:, i))
            move = max(move, value)
        end do
        move = move/radius
    end subroutine reseed_from_solution

    function seed_point(theta) result(p)
        !! (R, Z) of the present seed boundary at the angle theta (n = 0).
        real(dp), intent(in) :: theta
        real(dp) :: p(2)
        integer :: m
        p = 0
        do m = 0, ubound(boundary_input, 3)
            p(1) = p(1) + boundary_input(1, 0, m)*cos(m*theta) + boundary_input(2, 0, m)*sin(m*theta)
            p(2) = p(2) + boundary_input(3, 0, m)*cos(m*theta) + boundary_input(4, 0, m)*sin(m*theta)
        end do
    end function seed_point

    subroutine set_toroidal_level(modes)
        ! boundary padded with zero coefficients up to the map's modes
        integer, intent(in) :: modes
        integer :: local
        if (allocated(boundary)) deallocate (boundary)
        if (chart == 'deck') then
            ! the frame turns with the curve: the sections need the map's modes
            call deck_boundary_to_chart(boundary_input, ntor, nfp, deck_frame, deck_mpol, &
                modes, boundary, local)
            if (local /= 0) call fail('the boundary does not convert to the deck chart')
            return
        end if
        allocate (boundary(4, -modes:modes, 0:mpol), source=0.0_dp)
        boundary(:, -ntor:ntor, :) = boundary_input
    end subroutine

    subroutine observe(values)
        ! wb (relative), raxis relative to the minor radius, iota on axis and edge
        real(dp), intent(out) :: values(4)
        real(dp) :: v, a, w, k, r, i0, i1, unused1, unused2
        integer :: local
        call integrals(v, a, w, k, local)
        level_folded = local /= 0
        values = 0
        if (level_folded) return
        call axis_radius(r, local)
        if (local /= 0) call fail('the magnetic axis was not found')
        call profile%evaluate(0.0_dp, i0, unused1, unused2)
        call profile%evaluate(1.0_dp, i1, unused1, unused2)
        values = [w, r/sqrt(a/pi), i0, i1]
    end subroutine

    subroutine match_current(error)
        ! Prescribed enclosed toroidal current I(s), matched pointwise at the
        ! iota_level Chebyshev-Gauss surfaces s_j, which exclude the axis.
        ! Iota is the polynomial through its values c_j there, each level
        ! starting from the previous iota at the s_j. On a fixed map B on the
        ! surface s depends on iota(s) only, so I(s_j) = A_j c_j + C_j; with
        ! the map in equilibrium Newton solves F(c) = I(s; x(c), c) - I_target
        ! = 0 with dF/dc = diag(A) + (dI/dx) dx/dc, dx/dc = -H^-1 dg/dc
        ! (current_jacobian). error: the iota increment of the last step.
        real(dp), intent(out) :: error
        type(map_3d_loop_t), allocatable :: loops(:)
        real(dp), allocatable :: s(:), c(:), target(:), measured(:), step(:, :)
        real(dp), allocatable :: jacobian(:, :)
        real(dp), allocatable :: c0(:)
        real(dp) :: unused1, unused2, tau
        type(toroidal_profile_t) :: previous
        type(ideal_map_t) :: map0
        type(ideal_fb_t) :: fb0
        integer :: j, k, local, n, info, sign_attempt
        integer, allocatable :: pivot(:)
        n = iota_level
        s = iota_surfaces(n)
        allocate (loops(n), c(n), target(n), measured(n), step(n, 1), jacobian(n, n), &
            pivot(n))
        do j = 1, n
            call map_3d_loop_setup(map, s(j), loops(j), local, elliptic=free_boundary)
            if (local /= 0) call fail('a current loop met a curved cell edge')
            target(j) = current_target(s(j))
            call profile%evaluate(s(j), c(j), unused1, unused2)
        end do
        ! the map in equilibrium with the represented iota
        call set_iota(c)
        call resolve(local)
        if (local /= 0) call fail('the map folded when frozen with the represented iota')
        c0 = c; map0 = map
        if (fb_solving) fb0 = fb
        k = 0
        sign_attempt = 1
        do
            call measure(loops, measured)
            if (fb_solving) then
                call fb_current_jacobian(c, loops, measured, jacobian)
            else
                call current_jacobian(c, loops, measured, jacobian)
            end if
            step(:, 1) = target - measured
            call dgesv(n, 1, jacobian, n, pivot, step, n, info)
            if (info /= 0) call fail('the current Jacobian is singular')
            error = maxval(abs(step(:, 1)))
            print '(a,i3,a,es10.3)', 'current iteration ', k, ': iota increment ', error
            if (error <= current_tolerance .or. k == max_current_iterations) then
                ! a loosely solved map is solved fully and the step taken again
                if (solve_tolerance <= tolerance .or. k == max_current_iterations) exit
                solve_tolerance = tolerance
                if (fb_solving) then
                    call fb_solve(solve_tolerance)
                else
                    call ideal_map_solve(map, solve_tolerance, max_iterations, report, &
                        solve_status, reference_gradient=reference)
                end if
                cycle
            end if
            k = k + 1
            current_iterations = k
            ! The frozen map admits iota of one sign only (docs/ideal-seed-3d.md).
            ! A step may reverse that sign on every surface; one that would
            ! make iota change sign is halved until iota keeps one sign, which
            ! holds for a short enough step from the admitted c.
            previous = profile
            tau = 1
            do
                call set_iota(c + tau*step(:, 1))
                if (single_signed() .or. tau < epsilon(1.0_dp)) exit
                tau = tau/2
            end do
            if (tau < 1) print '(a,es10.3)', 'current iteration: step scaled by ', tau
            c = c + tau*step(:, 1)
            ! Iota within current_tolerance of zero has no resolved sign. The
            ! admitted iota form two components, iota > 0 and iota < 0, and
            ! the damped steps stay in one, so the other is searched from the
            ! mirrored start before the match is given up.
            if (.not. single_signed() .or. min_abs_iota() <= current_tolerance) then
                profile = previous
                print '(a,es10.3,a)', 'current iteration: the matched iota comes within ', &
                    current_tolerance, ' of zero'
                if (sign_attempt == 2) then
                    iota_vanishes = .true.
                    print '(a)', 'current iteration: on both signs, which the frozen map' &
                        //' does not admit'
                    exit
                end if
                sign_attempt = 2
                print '(a)', 'current iteration: starting from the iota of the other sign'
                c = -c0; map = map0
                if (fb_solving) fb = fb0
                call set_iota(c)
                solve_tolerance = tolerance
                call resolve(local)
                if (local /= 0) call fail('the map folded when frozen with the mirrored iota')
                cycle
            end if
            ! intermediate maps need only the accuracy of the iota increment
            solve_tolerance = max(tolerance, min(1e-4_dp, 1e-3_dp*error))
            call resolve(local)
            if (local /= 0) call fail('the map folded when frozen with the updated iota')
        end do
    end subroutine

    subroutine current_jacobian(c, loops, measured, jacobian)
        ! dI(s_j)/dc_k of the map in equilibrium. Fixed map: I(s_j) is affine in
        ! c_j alone, so one shift of every c by 1 gives the exact diagonal A.
        ! The energy gradient is quadratic in c, so dg/dc_k from the unit steps
        ! +-e_k is exact; the map's response dx_k = -H^-1 dg/dc_k
        ! (ideal_map_sensitivity) enters through a central difference of the
        ! loop integrals along dx_k, whose O(h^2) error changes only the rate of
        ! convergence, not the matched current.
        real(dp), intent(in) :: c(:), measured(:)
        type(map_3d_loop_t), intent(in) :: loops(:)
        real(dp), intent(out) :: jacobian(:, :)
        real(dp), allocatable :: dg(:, :), dx(:, :), plus(:), minus(:), x0(:, :, :)
        real(dp), allocatable :: angle0(:)
        real(dp) :: ip(size(c)), im(size(c)), e(size(c)), h, length
        integer :: k, local, n, ncoef
        n = size(c)
        ncoef = size(map%x)
        allocate (dg(map%ndof, n), dx(map%ndof, n), plus(map%ndof), minus(map%ndof))
        call set_iota(c + 1)
        seed%profile = profile
        call measure(loops, ip)
        jacobian = 0
        do k = 1, n
            jacobian(k, k) = ip(k) - measured(k)
        end do
        do k = 1, n
            e = 0; e(k) = 1
            call frozen_gradient(c + e, plus)
            call frozen_gradient(c - e, minus)
            dg(:, k) = (plus - minus)/2
        end do
        ! the map frozen with c again for its Hessian and the loop integrals
        call frozen_gradient(c, plus)
        call ideal_map_sensitivity(map, dg, dx, local)
        if (local /= 0) call fail('the map sensitivity solve failed')
        x0 = map%x
        if (map%nangle > 0) angle0 = map%angle
        length = maxval(maxval(map%reference, dim=2) - minval(map%reference, dim=2))
        do k = 1, n
            h = 1e-6_dp*length/max(maxval(abs(dx(:ncoef, k))), tiny(1.0_dp))
            call displace(x0, angle0, h, dx(:, k))
            call measure(loops, ip)
            call displace(x0, angle0, -h, dx(:, k))
            call measure(loops, im)
            jacobian(:, k) = jacobian(:, k) + (ip - im)/(2*h)
        end do
        map%x = x0
        if (map%nangle > 0) map%angle = angle0
    end subroutine

    subroutine frozen_gradient(c, gradient)
        ! energy gradient of the fixed map frozen with the iota values c
        real(dp), intent(in) :: c(:)
        real(dp), intent(out) :: gradient(:)
        real(dp) :: energy
        integer :: local
        call set_iota(c)
        seed%profile = profile
        call ideal_map_freeze(map, seed, local)
        if (local /= 0) call fail('the map folded when frozen in the current Jacobian')
        call ideal_map_energy(map, map%x, energy, local, gradient)
        if (local /= 0) call fail('the energy met a folded map in the current Jacobian')
    end subroutine

    subroutine displace(x0, angle0, t, direction)
        ! map coefficients and field-line angle moved by t along direction
        real(dp), intent(in) :: x0(:, :, :), t, direction(:)
        real(dp), allocatable, intent(in) :: angle0(:)
        integer :: ncoef
        ncoef = size(x0)
        map%x = x0 + t*reshape(direction(:ncoef), shape(x0))
        if (map%nangle > 0) map%angle = angle0 + t*direction(ncoef + 1:)
    end subroutine

    subroutine measure(loops, values)
        ! enclosed toroidal current (A) through the loops
        type(map_3d_loop_t), intent(in) :: loops(:)
        real(dp), intent(out) :: values(:)
        integer :: j, local
        do j = 1, size(loops)
            values(j) = map_3d_loop_integral(map, seed, loops(j), local)/mu0
            if (local /= 0) call fail('the current readback met a folded map')
        end do
    end subroutine

    subroutine set_iota(c)
        ! the polynomial through the values c at the Chebyshev-Gauss surfaces,
        ! stored by its values and slopes on store_knots uniform knots
        real(dp), intent(in) :: c(:)
        real(dp) :: t(store_knots), value(store_knots), slope(store_knots)
        real(dp) :: w(size(c)), dw(size(c))
        integer :: k
        do k = 1, store_knots
            t(k) = real(k - 1, dp)/(store_knots - 1)
            call lagrange_weights(size(c), t(k), iota_surfaces(size(c)), w)
            call lagrange_deriv_weights(size(c), t(k), iota_surfaces(size(c)), dw)
            value(k) = dot_product(w, c); slope(k) = dot_product(dw, c)
        end do
        call profile%set_hermite_iota(t, value, slope)
    end subroutine

    pure function iota_surfaces(n) result(x)
        ! Chebyshev-Gauss surfaces, which exclude the axis
        integer, intent(in) :: n
        real(dp) :: x(n)
        integer :: k
        x = [((1 - cos((2*k - 1)*pi/(2*n)))/2, k=1, n)]
    end function

    real(dp) function current_target(t)
        ! prescribed enclosed toroidal current (A) at normalised toroidal flux t
        real(dp), intent(in) :: t
        real(dp) :: unused
        integer :: k
        if (current_type == 'hermite') then
            call hermite_profile_evaluate(profile_knots(:n_knots), current_values(:n_knots), &
                current_slopes(:n_knots), t, current_target, unused)
        else
            current_target = sum(current*t**[(k, k=0, kmax)])
        end if
    end function

    subroutine admit_axisymmetric_current()
        ! On an axisymmetric boundary Ampere's law on the surface s gives
        ! mu0 I(s) = psi'(s) K(s) with K(s) = oint |grad s|/R dl > 0 and K ~ s
        ! near the axis, so iota = psi'/Phi' has the zeros of I(s)/s: a zero of
        ! I on (0, 1] or I'(0) = 0 makes iota vanish, which the frozen map does
        ! not admit. Sampled on the grid of the seed's iota admission.
        real(dp) :: slope, value
        integer :: a
        character(16) :: at
        slope = current_slope_axis()
        if (.not. abs(slope) > 0) call fail('on an axisymmetric boundary the prescribed' &
            //' current has dI/ds = 0 on the axis, so iota vanishes there, which the' &
            //' frozen map does not admit')
        do a = 1, 128
            value = current_target(a/128.0_dp)
            if (value*slope > 0) cycle
            write (at, '(f8.5)') a/128.0_dp
            call fail('on an axisymmetric boundary the prescribed current changes sign' &
                //' or vanishes by s = '//trim(adjustl(at))//', so iota does, which the' &
                //' frozen map does not admit')
        end do
    end subroutine

    real(dp) function current_slope_axis()
        ! dI/ds at s = 0 of the prescribed current
        real(dp) :: h, unused
        if (current_type == 'hermite') then
            h = profile_knots(2) - profile_knots(1)
            call seed_hermite_profile_jvp(-profile_knots(1)/h, h, current_values(1), &
                current_slopes(1), current_values(2), current_slopes(2), 1.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, current_slope_axis, unused)
            current_slope_axis = current_slope_axis/h
        else
            current_slope_axis = current(1)
        end if
    end function

    subroutine init_seed(status)
        integer, intent(out) :: status
        character(8) :: lagged
        integer :: length

        ! every plane from the Winslow map of its cross-section
        if (angle_mpol >= 0) then
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, angle_mpol=angle_mpol, winslow=.true., chart=map_chart)
        else
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, winslow=.true., chart=map_chart)
        end if
        call get_environment_variable('MILCA_NESTED_BENCHMARK_LAGGED', lagged, length)
        if (status == 0 .and. length > 0) map%lagged_hessian = lagged == '1'
    end subroutine

    subroutine resolve(local)
        ! freeze the new iota on the solved map and solve again
        integer, intent(out) :: local
        seed%profile = profile
        if (fb_solving) then
            call ideal_map_freeze(fb%plasma, seed, local)
            if (local /= 0) return
            call fb_solve(solve_tolerance)
            return
        end if
        call ideal_map_freeze(map, seed, local)
        if (local /= 0) return
        ! the same absolute target as the first solve, relative to its seed
        call ideal_map_solve(map, solve_tolerance, max_iterations, report, solve_status, &
            reference_gradient=reference)
    end subroutine

    real(dp) function free_gradient()
        ! l2 norm of the free reduced energy gradient of the current map
        real(dp), allocatable :: gradient(:), reduced(:)
        real(dp) :: energy
        integer :: local
        allocate (gradient(map%ndof), reduced(map%nred))
        call ideal_map_energy(map, map%x, energy, local, gradient)
        if (local /= 0) call fail('the energy met a folded map')
        call ideal_map_reduce(map, gradient, reduced)
        where (map%fixed) reduced = 0
        free_gradient = norm2(reduced)
    end function

    logical function kind_ok(kind)
        character(*), intent(in) :: kind
        kind_ok = kind == 'power_series' .or. kind == 'hermite'
    end function

    integer function last(c)
        ! highest nonzero coefficient, at least the constant
        real(dp), intent(in) :: c(0:)
        do last = ubound(c, 1), 1, -1
            if (c(last) /= 0) return
        end do
    end function

    subroutine chart_point(e, xi, phi, tangent, jac, magnetic, pressure, radius, along, status)
        ! The solved map at master point xi of cell e and angle phi-hat, in any
        ! chart: jac = dV/d(xi, phi-hat), magnetic = B^2/2 jac, the pressure,
        ! the major radius and along = B . dX/dxi tangent. The energy density
        ! at p = 0 is |DX b|^2/(2 mu0 J) and its b-gradient DX^T B/mu0.
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2), phi, tangent(2)
        real(dp), intent(out) :: jac, magnetic, pressure, radius, along
        integer, intent(out) :: status
        real(dp) :: q(3), d(3, 3), ref(2), g(2, 2), chi, cg(2), b3, dcg(2), db3, b(3), &
                    pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), gv(12), gb(3)
        status = 1
        call map_3d_coordinates(map, e, xi, phi, q, d, ref, g)
        call seed%evaluate(ref, chi, cg, b3, pressure)
        call ideal_seed_3d_angle_density(map, seed, ref, phi, dcg, db3)
        cg = matmul(transpose(g), cg + dcg)
        b = [-cg(2), cg(1), (b3 + db3)*(g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))]
        call map%chart%jet(q, pos, d1, d2, d3)
        call ideal_point_gradient_jets(d1, d2, d3, d, b, 0.0_dp, map%mu0, map%sigma, &
            magnetic, jac, gv)
        call ideal_point_field_gradient_jets(d1, d2, d3, d, b, 0.0_dp, map%mu0, map%sigma, gb)
        radius = hypot(pos(1), pos(2))
        along = map%mu0*dot_product(gb(1:2), tangent)
        if (.not. (jac > 0 .and. radius > 0)) return
        if (.not. all(ieee_is_finite([magnetic, along]))) return
        status = 0
    end subroutine

    subroutine integrals(volume, area, magnetic, kinetic, status)
        ! Gauss-Legendre on the collapsed reference triangle and the trapezoidal
        ! rule in phi, exact for the periodic Fourier dependence of the map.
        real(dp), intent(out) :: volume, area, magnetic, kinetic
        integer, intent(out) :: status
        integer, parameter :: nq = 8
        real(dp) :: x(nq), w(nq), xi(2), weight, phi, jac, sums(4), density, pressure, &
                    radius, along
        integer :: e, i, j, k, nphi, failed, local
        call gauss_legendre(nq, x, w)
        x = (x + 1)/2; w = w/2
        nphi = max(16, 4*map%ntor + 4)
        sums = 0; failed = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(e,i,j,k,xi,weight,phi,jac,density,pressure,radius,along,local) &
        !$omp reduction(+:sums,failed)
        do e = 1, map%cells
            do k = 0, nphi - 1
                phi = 2*pi*k/(nfp*nphi)
                do j = 1, nq
                    do i = 1, nq
                        xi = [x(i)*(1 - x(j)), x(j)]
                        weight = w(i)*w(j)*(1 - x(j))
                        call chart_point(e, xi, phi, [0.0_dp, 0.0_dp], jac, density, &
                            pressure, radius, along, local)
                        failed = failed + local
                        if (local /= 0) cycle
                        sums = sums + weight*[jac, jac/radius, density, jac*pressure]
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = failed
        ! mean cross-section area = int dV/R / (2 pi)
        volume = 2*pi*sums(1)/nphi; area = sums(2)/nphi
        magnetic = 2*pi*sums(3)/nphi; kinetic = 2*pi*sums(4)/nphi
    end subroutine

    real(dp) function edge_current(status)
        ! (1/mu0) oint B.dl along the discrete wall at phi = 0, the enclosed
        ! toroidal current; Gauss-Legendre on every wall edge
        integer, intent(out) :: status
        integer, parameter :: ng = 8
        real(dp) :: centre(2), radius(2), xi(2)
        real(dp) :: corner(2, 3)
        real(dp) :: x(ng), w(ng), q(3), d(3, 3), ref(2), g(2, 2), tangent(2), along
        real(dp) :: sense, jac, density, pressure, major
        integer :: k, e, local, i, a
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        radius = (maxval(map%reference(1, :)) - minval(map%reference(1, :)))/2
        ! the free boundary's reference disc is an ellipse
        if (free_boundary) radius = (maxval(map%reference, dim=2) &
            - minval(map%reference, dim=2))/2
        edge_current = 0; status = 1
        corner = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        call gauss_legendre(ng, x, w)
        x = (x + 1)/2; w = w/2
        do e = 1, map%cells
            do i = 1, 3
                a = mod(i, 3) + 1
                if (.not. (map%boundary(map%conn(i, e)) .and. map%boundary(map%conn(a, e)))) cycle
                ! a wall edge, not a chord between two wall nodes
                call map_3d_coordinates(map, e, (corner(:, i) + corner(:, a))/2, 0.0_dp, q, d, ref, g)
                if (abs(norm2((ref - centre)/radius) - 1) > 1e-3_dp) cycle
                tangent = corner(:, a) - corner(:, i)
                sense = sign(1.0_dp, (ref(1) - centre(1))*dot_product(g(2, :), tangent) &
                    - (ref(2) - centre(2))*dot_product(g(1, :), tangent))
                do k = 1, ng
                    xi = corner(:, i) + x(k)*tangent
                    call chart_point(e, xi, 0.0_dp, tangent, jac, density, pressure, major, &
                        along, local)
                    if (local /= 0) return
                    edge_current = edge_current + sense*w(k)*along
                end do
            end do
        end do
        edge_current = edge_current/mu0
        status = 0
    end function

    subroutine axis_radius(raxis, status)
        ! image of the centre of the reference disc (the label's O-point),
        ! averaged over the cylindrical angle
        real(dp), intent(out) :: raxis
        integer, intent(out) :: status
        real(dp) :: centre(2), xi(2), q(3), d(3, 3), ref(2), g(2, 2), weight, dphi
        real(dp) :: pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), tangent(3)
        integer :: e, k, nphi, local
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        status = 1
        do e = 1, map%cells
            call triangle_locate(map%reference(:, map%conn(:, e)), centre, xi, local)
            if (local /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1e-10_dp) cycle
            nphi = max(16, 4*map%ntor + 4)
            ! mean over the cylindrical angle, which the map displaces:
            ! dphi_cyl/dphi-hat = (X1 dX2 - X2 dX1)/R^2 along the axis
            raxis = 0; weight = 0
            do k = 0, nphi - 1
                call map_3d_coordinates(map, e, xi, 2*pi*k/(nfp*nphi), q, d, ref, g)
                call map%chart%jet(q, pos, d1, d2, d3)
                tangent = matmul(d1, d(:, 3))
                dphi = (pos(1)*tangent(2) - pos(2)*tangent(1))/(pos(1)**2 + pos(2)**2)
                raxis = raxis + hypot(pos(1), pos(2))*dphi; weight = weight + dphi
            end do
            raxis = raxis/weight
            status = 0
            return
        end do
    end subroutine
end program nested_equilibrium
