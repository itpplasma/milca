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
        ideal_seed_3d_init, ideal_seed_disk_mesh, hermite_profile_evaluate
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce, ideal_map_prolong_correction, ideal_map_freeze, &
        ideal_map_sensitivity
    use milca_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use milca_grad_shafranov, only: gs_state_t
    use milca_curved_triangle, only: triangle_locate
    use milca_result_io, only: run_t, result_file_t
    use fortnum_quadrature, only: gauss_legendre
    use milca_map_3d_loop, only: map_3d_loop_t, map_3d_loop_setup, map_3d_loop_integral
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use milca_ideal_seed_3d_leaves, only: seed_hermite_profile_jvp
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, ideal_fb_solve, &
        ideal_fb_report_t, ideal_fb_drive_traces, ideal_fb_set_traces, &
        ideal_fb_gradient, ideal_fb_sensitivity, IDEAL_FB_OK, IDEAL_FB_FOLDED, &
        ideal_fb_boundary_polygon
    use milca_ideal_fb_circuits, only: ideal_fb_circuits_t, ideal_fb_circuits_build
    use milca_ideal_fb_mgrid, only: mgrid_basis_t, mgrid_basis_init, mgrid_current_nodes
    use milca_ideal_fb_nested, only: fb_nested_offset_distance, fb_nested_layer_mesh, &
        fb_nested_section, &
        fb_nested_geometry, fb_nested_mesh, &
        fb_nested_seed, fb_nested_transfer, fb_nested_settle, filament_basis_t
    use milca_mgrid_file, only: mgrid_t, mgrid_read
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
    logical :: mesh_confirmed, modes_confirmed, iota_confirmed
    logical :: iota_vanishes = .false.
    ! the solved map of the level folds at the volume-integral points
    logical :: level_folded = .false.
    real(dp) :: fourier_tolerance, fourier_error, mesh_tolerance, mesh_error
    real(dp) :: iota_tolerance, iota_error
    real(dp) :: observed(4), previous_observed(4)
    real(dp), allocatable :: boundary_input(:, :, :)
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
    logical :: prescribe_current
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
    type(cylindrical_chart_t) :: chart
    real(dp), allocatable :: filament_drive(:, :), filament_load(:, :)
    real(dp) :: fb_centre(2), fb_axes(2), circulation, coil_closure
    real(dp), allocatable :: fb_curve(:)
    ! shell_layers > 0: the shell is layered along segments from Gamma_p to
    ! Gamma_c (fb_nested_layer_mesh), solved by the inertia-controlled
    ! Newton; 0: the harmonic-extension shell to the curve of
    ! fb_nested_geometry.
    integer :: shell_layers
    real(dp) :: shell_distance
    real(dp), allocatable :: shell_positions(:, :), shell_t(:)
    integer, allocatable :: shell_lines(:, :)
    logical, allocatable :: shell_cells(:)
    ! Gamma_c: the convex hull of shell_section (the guess's section, or a
    ! solved Gamma_p when Gamma_c is re-placed) offset by shell_distance.
    real(dp), allocatable :: shell_section(:, :), shell_obstacles(:, :)
    real(dp), allocatable :: shell_section_reference(:, :)
    real(dp) :: shell_box(4)
    integer :: replacements = 0
    namelist /nested_input/ nfp, mpol, ntor, rbc, rbs, zbc, zbs, phiedge, &
        iota, pressure, boundary_nodes, degree, tolerance, max_iterations, angle_mpol, &
        ladder_nodes, ladder_degree, prescribe_current, current, max_current_iterations, &
        current_tolerance, iota_type, pressure_type, &
        current_type, n_knots, profile_knots, iota_values, iota_slopes, pressure_values, &
        pressure_slopes, current_values, current_slopes, map_ntor, fourier_tolerance, mesh_tolerance, &
        iota_nodes, iota_tolerance, free_boundary, mgrid_file, extcur, mgrid_degree, &
        coupling_margin, shell_layers

    nfp = 1; mpol = 1; ntor = 0; phiedge = 1
    rbc = 0; rbs = 0; zbc = 0; zbs = 0; iota = 0; pressure = 0
    boundary_nodes = -1; degree = 3; tolerance = 1e-8_dp; max_iterations = 200
    angle_mpol = -1; ladder_nodes = 0; ladder_degree = 0
    prescribe_current = .false.; current = 0; max_current_iterations = 30
    current_tolerance = 1e-5_dp
    map_ntor = -1; fourier_tolerance = 1e-4_dp
    mesh_tolerance = 1e-4_dp
    iota_nodes = -1; iota_tolerance = 1e-4_dp
    free_boundary = .false.; mgrid_file = ''; extcur = 0; mgrid_degree = 5
    coupling_margin = 0.5_dp; shell_layers = 2
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
    if (free_boundary .and. len_trim(mgrid_file) == 0) &
        error stop 'nested_equilibrium: free_boundary needs mgrid_file'
    if (free_boundary .and. .not. coupling_margin > 0) &
        error stop 'nested_equilibrium: coupling_margin must be positive'
    call run%open_resolved(unit)
    write (unit, nml=nested_input)
    call run%close_resolved(unit)

    ! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp phi)
    allocate (boundary_input(4, -ntor:ntor, 0:mpol))
    boundary_input(1, :, :) = rbc(-ntor:ntor, 0:mpol)
    boundary_input(2, :, :) = rbs(-ntor:ntor, 0:mpol)
    boundary_input(3, :, :) = zbc(-ntor:ntor, 0:mpol)
    boundary_input(4, :, :) = zbs(-ntor:ntor, 0:mpol)
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
        call free_boundary_setup()
        call fb_level(boundary_nodes, .true.)
    else
        call mesh_ladder()
    end if
    call finish_level()
    call require_unfolded()
    do
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
    call result%attribute('converged', merge('true ', 'false', &
        solve_status == 0 .and. residual <= tolerance .and. current_error <= current_tolerance &
        .and. .not. iota_vanishes &
        .and. (.not. automatic_ntor .or. fourier_error <= fourier_tolerance) &
        .and. (.not. automatic_mesh .or. mesh_error <= mesh_tolerance) &
        .and. (.not. automatic_iota .or. iota_error <= iota_tolerance)))
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
            if (step > 1) then
                call ideal_map_prolong_correction(map, coarse, coarse_seed, status)
                if (status /= 0) call fail('the resolution ladder correction failed')
            end if
            call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
                reference_gradient=reference)
            coarse = map; coarse_seed = seed_x
        end do
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
        if (shell_layers > 0) then
            call fb_nested_section(boundary_input, fb_centre, fb_axes, shell_section)
            shell_box = [grid%rmin, grid%rmax, grid%zmin, grid%zmax]
            shell_obstacles = coils_inside
            call fb_nested_offset_distance(shell_section, shell_box, shell_obstacles, &
                coupling_margin*sqrt(fb_axes(1)*fb_axes(2)), shell_distance, local)
            if (local /= 0) call fail('no room for the coupling surface between the' &
                //' guess and the coils')
            margin_used = shell_distance/sqrt(fb_axes(1)*fb_axes(2))
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

    recursive subroutine fb_level(nodes, first)
        ! Free boundary at the mesh of nodes boundary nodes and the current
        ! toroidal modes: coupling-surface operators and traces, the start
        ! from the guess (first) or from the solved level, and the solve.
        integer, intent(in) :: nodes
        logical, intent(in) :: first
        type(ideal_fb_t) :: next
        real(dp), allocatable :: chain(:, :), d(:, :, :), l(:, :, :)
        real(dp) :: closure(1), f(1), spread(1)
        integer :: local
        seed%profile = profile
        if (shell_layers > 0) then
            call fb_nested_layer_mesh(boundary, fb_centre, fb_axes, nodes, degree, &
                shell_layers, shell_distance, mesh, chain, shell_positions, &
                shell_lines, shell_t, local, plasma_cells=shell_cells, &
                section=shell_section, section_reference=shell_section_reference)
            if (local /= 0) call fail('the free-boundary mesh failed')
            call ideal_fb_init(next, mesh, chain, chart, seed, nfp, ntor_level, &
                0.0_dp, circulation, local, nphi=max(1, 4*ntor_level + 1), &
                positions=shell_positions, lines=shell_lines, line_t=shell_t, &
                plasma_cells=shell_cells)
        else
            call fb_nested_mesh(mesh, chain, fb_centre, fb_axes, fb_curve, nodes, &
                degree, local)
            if (local /= 0) call fail('the free-boundary mesh failed')
            call ideal_fb_init(next, mesh, chain, chart, seed, nfp, ntor_level, &
                0.0_dp, circulation, local, nphi=max(1, 4*ntor_level + 1))
        end if
        if (local /= IDEAL_FB_OK) call fail('the free-boundary set-up failed')
        call ideal_fb_circuits_build(next, coil_basis, circuits, local)
        if (local /= IDEAL_FB_OK) &
            call fail('the coil field grid does not cover the coupling surface')
        allocate (d(next%nt, next%nmode, 1), l(next%nt, next%nmode, 1))
        call ideal_fb_drive_traces(next, filament, [1.0_dp], d, l, closure, f, spread, &
            local)
        if (local /= IDEAL_FB_OK) call fail('the plasma-current traces failed')
        filament_drive = d(:, :, 1); filament_load = l(:, :, 1)
        ! bound of the closure of the summed coil field: a coil inside Gamma_c
        ! (a toroidal current it encloses) shows as an n = 0 closure
        coil_closure = sum(abs(extcur(:circuits%ncircuit))*circuits%closure)
        call install_traces(next)
        print '(a,es12.4,a,es10.2,a,i0,a,i0)', 'F_vac (T m)', next%f_vac, &
            ', coil closure', coil_closure, ', unknowns ', next%n, ', modes ', &
            ntor_level
        if (.not. next%f_vac > 0) call fail('the coil toroidal field is not positive;' &
            //' reverse extcur and the current for a negative toroidal flux')
        if (first) then
            call fb_nested_seed(next, boundary, fb_centre, fb_axes, local)
            if (local /= IDEAL_FB_OK) &
                call fail('the free-boundary seed is not admitted')
        else
            call fb_nested_transfer(next, fb, local)
            if (local /= IDEAL_FB_OK) call fail('the resolution level transfer failed')
        end if
        if (first) then
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
        end if
        fb = next
        ! The layered shell is solved with the inertia-controlled Newton
        ! (the stable state, k = 0) and its fold barrier.
        fb%inertia_control = shell_layers > 0
        call fb_solve(tolerance)
        ! A solve that converges only with the fold barrier active (or folds)
        ! has Gamma_c too close: Gamma_c is re-placed around the solved
        ! boundary and the level is solved again from the transferred map.
        if (solve_status /= IDEAL_FB_OK .and. shell_layers > 0 .and. &
            replacements < 3) then
            replacements = replacements + 1
            call ideal_fb_boundary_polygon(fb, 16, shell_section, local, &
                reference=shell_section_reference)
            if (local /= 0) call fail('the solved boundary is not a closed chain')
            call fb_nested_offset_distance(shell_section, shell_box, shell_obstacles, &
                coupling_margin*sqrt(fb_axes(1)*fb_axes(2)), shell_distance, local)
            if (local /= 0) call fail('no room for the coupling surface around the' &
                //' solved boundary')
            print '(a,f6.3)', 'Gamma_c re-placed around the solved boundary at d/a ', &
                shell_distance/sqrt(fb_axes(1)*fb_axes(2))
            call fb_level(nodes, .false.)
            replacements = replacements - 1
        end if
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
        map%position_gauge = .true.
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
        d = circulation*filament_drive; l = circulation*filament_load; f = 0
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
                solve_status, verbose=.false., residual_tolerance=target)
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
            call require_unfolded()
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
        previous_map = map; previous_observed = observed
        if (allocated(seed_x)) previous_seed = seed_x
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
        call ideal_map_prolong_correction(map, previous_map, previous_seed, status)
        if (status /= 0) call fail('the resolution level correction failed')
        call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
            reference_gradient=reference)
        call finish_level()
        ! The energy sees the map only at its quadrature points, and the
        ! continued start can reach a stationary point that folds between
        ! them (ideal_map_prolong_correction); the level is solved again from
        ! its injective Winslow seed.
        if (level_folded) then
            print '(a,i4,a,i3,a)', 'boundary nodes ', boundary_nodes, ', toroidal modes ', &
                ntor_level, ': the continued map folds, solving from the seed'
            call init_seed(status)
            if (status /= 0) call fail('the seed of a resolution level is not admitted')
            reference = free_gradient()
            seed_x = map%x
            call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
                reference_gradient=reference)
            call finish_level()
        end if
        call require_unfolded()
    end subroutine

    subroutine require_unfolded()
        if (level_folded) call fail('a volume integral met a folded map')
    end subroutine

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

    subroutine set_toroidal_level(modes)
        ! boundary padded with zero coefficients up to the map's modes
        integer, intent(in) :: modes
        if (allocated(boundary)) deallocate (boundary)
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
        type(ideal_fb_t) :: fb0, fb_good
        real(dp), allocatable :: c_good(:)
        integer :: j, k, local, n, info, sign_attempt, halvings
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
        if (fb_solving .and. solve_status /= IDEAL_FB_OK) &
            call fail('the free-boundary solve did not converge before current matching')
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
            c_good = c
            if (fb_solving) fb_good = fb
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
            ! An unconverged free-boundary solve has no current to match: the
            ! iota step is halved from the last converged state until the
            ! solve converges.
            halvings = 0
            do while (fb_solving .and. solve_status /= IDEAL_FB_OK)
                halvings = halvings + 1
                if (halvings > 6) call fail('the free-boundary solve does not converge' &
                    //' along the current-matching step')
                tau = tau/2
                print '(a,es10.3)', 'current iteration: unconverged solve, step scaled by ', &
                    tau
                c = c_good + tau*step(:, 1)
                fb = fb_good
                call set_iota(c)
                call resolve(local)
                if (local /= 0) call fail('the map folded when frozen with the updated iota')
            end do
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
        ! every plane from the Winslow map of its cross-section
        if (angle_mpol >= 0) then
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, angle_mpol=angle_mpol, winslow=.true.)
        else
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, winslow=.true.)
        end if
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

    subroutine integrals(volume, area, magnetic, kinetic, status)
        ! Gauss-Legendre on the collapsed reference triangle and the trapezoidal
        ! rule in phi, exact for the periodic Fourier dependence of the map.
        real(dp), intent(out) :: volume, area, magnetic, kinetic
        integer, intent(out) :: status
        integer, parameter :: nq = 8
        real(dp) :: x(nq), w(nq), xi(2), weight, q(3), d(3, 3), ref(2), g(2, 2)
        real(dp) :: point(3), values(EQ_VALUES), phi, jac, sums(4)
        integer :: e, i, j, k, nphi, failed, local
        call gauss_legendre(nq, x, w)
        x = (x + 1)/2; w = w/2
        nphi = max(16, 4*map%ntor + 4)
        sums = 0; failed = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(e,i,j,k,xi,weight,phi,q,d,ref,g,point,values,jac,local) &
        !$omp reduction(+:sums,failed)
        do e = 1, map%cells
            do k = 0, nphi - 1
                phi = 2*pi*k/(nfp*nphi)
                do j = 1, nq
                    do i = 1, nq
                        xi = [x(i)*(1 - x(j)), x(j)]
                        weight = w(i)*w(j)*(1 - x(j))
                        call map_3d_coordinates(map, e, xi, phi, q, d, ref, g)
                        ! dV = R |det dq/d(xi, phi)|; the map also displaces phi
                        jac = abs(d(1, 1)*(d(2, 2)*d(3, 3) - d(2, 3)*d(3, 2)) &
                            - d(1, 2)*(d(2, 1)*d(3, 3) - d(2, 3)*d(3, 1)) &
                            + d(1, 3)*(d(2, 1)*d(3, 2) - d(2, 2)*d(3, 1)))
                        call map_3d_sample(map, seed, e, xi, phi, point, values, local)
                        failed = failed + local
                        sums = sums + weight*[q(1)*jac, jac, &
                            q(1)*jac*sum(values([EQ_BR, EQ_BZ, EQ_BPHI])**2)/2, &
                            q(1)*jac*values(EQ_P)]
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
        real(dp) :: centre(2), radius(2), xi(2), point(3), values(EQ_VALUES)
        real(dp) :: corner(2, 3)
        real(dp) :: x(ng), w(ng), q(3), d(3, 3), ref(2), g(2, 2), tangent(2), along(3)
        real(dp) :: sense
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
                    call map_3d_coordinates(map, e, xi, 0.0_dp, q, d, ref, g)
                    call map_3d_sample(map, seed, e, xi, 0.0_dp, point, values, local)
                    if (local /= 0) return
                    along = matmul(d(:, 1:2), tangent); along(3) = along(3)*q(1)
                    edge_current = edge_current + sense*w(k)* &
                        dot_product(values([EQ_BR, EQ_BZ, EQ_BPHI]), along)
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
        real(dp) :: centre(2), xi(2), q(3), d(3, 3), ref(2), g(2, 2), weight
        integer :: e, k, nphi, local
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        status = 1
        do e = 1, map%cells
            call triangle_locate(map%reference(:, map%conn(:, e)), centre, xi, local)
            if (local /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1e-10_dp) cycle
            nphi = max(16, 4*map%ntor + 4)
            ! mean over the cylindrical angle, which the map displaces: dphi_cyl = d(3, 3) dphi
            raxis = 0; weight = 0
            do k = 0, nphi - 1
                call map_3d_coordinates(map, e, xi, 2*pi*k/(nfp*nphi), q, d, ref, g)
                raxis = raxis + q(1)*d(3, 3); weight = weight + d(3, 3)
            end do
            raxis = raxis/weight
            status = 0
            return
        end do
    end subroutine
end program nested_equilibrium
