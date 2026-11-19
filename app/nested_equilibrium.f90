program nested_equilibrium
    ! Fixed-boundary nested 3D equilibrium from a boundary, the toroidal flux and
    ! the profiles iota(s), p(s) in the normalized toroidal flux s, solved by the
    ! native map Newton (docs/nested-equilibrium-input.md). Writes the solved
    ! state's global quantities for cross-code comparison.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh, hermite_profile_evaluate
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce, ideal_map_prolong_correction, ideal_map_freeze
    use milca_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use milca_grad_shafranov, only: gs_state_t
    use milca_curved_triangle, only: triangle_locate
    use milca_result_io, only: run_t, result_file_t
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    implicit none
    integer, parameter :: mmax = 64, nmax = 32, kmax = 32, knot_max = 1024
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi
    integer :: nfp, mpol, ntor, boundary_nodes, degree, max_iterations, angle_mpol
    integer :: unit, ios, status, solve_status, step, steps, nodes, order, n
    integer :: ladder_nodes(8), ladder_degree(8)
    integer :: map_ntor, ntor_level, level, current_nodes
    logical :: automatic_ntor, iota_crosses_zero = .false.
    real(dp) :: fourier_tolerance, fourier_error, observed(4), previous_observed(4)
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
    namelist /nested_input/ nfp, mpol, ntor, rbc, rbs, zbc, zbs, phiedge, &
        iota, pressure, boundary_nodes, degree, tolerance, max_iterations, angle_mpol, &
        ladder_nodes, ladder_degree, prescribe_current, current, max_current_iterations, &
        current_tolerance, iota_type, pressure_type, &
        current_type, n_knots, profile_knots, iota_values, iota_slopes, pressure_values, &
        pressure_slopes, current_values, current_slopes, map_ntor, fourier_tolerance, current_nodes

    nfp = 1; mpol = 1; ntor = 0; phiedge = 1
    rbc = 0; rbs = 0; zbc = 0; zbs = 0; iota = 0; pressure = 0
    boundary_nodes = 16; degree = 3; tolerance = 1e-8_dp; max_iterations = 200
    angle_mpol = -1; ladder_nodes = 0; ladder_degree = 0
    prescribe_current = .false.; current = 0; max_current_iterations = 30
    current_tolerance = 1e-5_dp
    map_ntor = -1; fourier_tolerance = 1e-4_dp; current_nodes = 16
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
        error stop 'nested_equilibrium: map_ntor below the boundary ntor or above 32'
    if (.not. fourier_tolerance > 0) error stop 'nested_equilibrium: fourier_tolerance'
    if (current_nodes < 2 .or. current_nodes > 64) &
        error stop 'nested_equilibrium: current_nodes from 2 to 64'
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

    ! Interior toroidal modes: the solution of a 3D boundary needs more
    ! toroidal harmonics than the boundary has. With map_ntor < 0 they are
    ! chosen by the change of the observables between levels ntor_level - 2
    ! and ntor_level (fourier_error); axisymmetric data keep ntor = 0.
    automatic_ntor = map_ntor < 0 .and. ntor > 0
    ! the change between levels estimates the truncation error only when each
    ! level's current is matched well below it
    if (automatic_ntor) current_tolerance = min(current_tolerance, fourier_tolerance/10)
    ntor_level = merge(ntor, map_ntor, map_ntor < 0)
    fourier_error = 0
    level = 0
    do
        call set_toroidal_level(ntor_level)
        if (level == 0) then
            call mesh_ladder()
        else
            ! seed with the new modes plus the correction of the solved level
            call ideal_seed_disk_mesh(mesh, boundary_nodes, degree, status)
            if (status /= 0) error stop 'nested_equilibrium: disk mesh'
            call init_seed(status)
            if (status /= 0) call fail('the seed of a toroidal level is not admitted')
            reference = free_gradient()
            seed_x = map%x
            call ideal_map_prolong_correction(map, previous_map, previous_seed, status)
            if (status /= 0) call fail('the toroidal level correction failed')
            call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
                reference_gradient=reference)
        end if
        current_error = 0; current_iterations = 0; solve_tolerance = tolerance
        if (prescribe_current) call match_current(current_error)
        call observe(observed)
        if (level > 0) fourier_error = max(abs(observed(1)/previous_observed(1) - 1), &
            maxval(abs(observed(2:) - previous_observed(2:))))
        print '(a,i3,a,es10.3)', 'toroidal modes ', ntor_level, ': fourier error ', &
            fourier_error
        if (iota_crosses_zero .or. .not. automatic_ntor) exit
        if (level > 0 .and. fourier_error <= fourier_tolerance) exit
        if (ntor_level + 2 > nmax) exit
        previous_map = map; previous_seed = seed_x; previous_observed = observed
        ntor_level = ntor_level + 2; level = level + 1
    end do
    call profile%evaluate(0.0_dp, iota_axis, unused, unused)
    call profile%evaluate(1.0_dp, iota_edge, unused, unused)
    residual = free_gradient()/max(tiny(1.0_dp), reference)
    call integrals(volume, area, magnetic, kinetic, status)
    if (status /= 0) call fail('a volume integral met a folded map')
    call axis_radius(raxis, status)
    if (status /= 0) call fail('the magnetic axis was not found')
    ctor = enclosed_current(1.0_dp, status)
    if (status /= 0) call fail('the edge current readback met a folded map')

    call result%create(run)
    call result%attribute('physical_model', 'fixed-boundary nested ideal MHD, prescribed iota')
    call result%attribute('converged', merge('true ', 'false', &
        solve_status == 0 .and. residual <= tolerance .and. current_error <= current_tolerance &
        .and. .not. iota_crosses_zero &
        .and. (.not. automatic_ntor .or. fourier_error <= fourier_tolerance)))
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
    if (prescribe_current) then
        call result%attribute('current_status', merge('iota crosses zero', 'matched          ', &
            iota_crosses_zero))
        call result%write('current_error', current_error, '1', &
            'largest enclosed-current error relative to the largest prescribed current')
        call result%write('current_iterations', real(current_iterations, dp), '1', &
            'iota updates to match the prescribed current')
    end if
    call result%write('map_ntor', real(ntor_level, dp), '1', 'interior toroidal modes')
    if (automatic_ntor) call result%write('fourier_error', fourier_error, '1', &
        'largest change of wb (relative), raxis/aminor_p and iota on axis and edge' &
        //' from map_ntor - 2')
    call result%write('force_residual', residual, '1', 'relative free reduced-gradient l2')
    call result%write('tolerance', tolerance, '1', 'requested force_residual')
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
        if (local /= 0) call fail('a volume integral met a folded map')
        call axis_radius(r, local)
        if (local /= 0) call fail('the magnetic axis was not found')
        call profile%evaluate(0.0_dp, i0, unused1, unused2)
        call profile%evaluate(1.0_dp, i1, unused1, unused2)
        values = [w, r/sqrt(a/pi), i0, i1]
    end subroutine

    subroutine match_current(error)
        ! Prescribed enclosed toroidal current I(s). For a fixed map I on each
        ! surface is affine in iota there; its slope is read on the solved map
        ! with the frozen form's iota shifted by delta, so every update is the
        ! exact Newton step for that map at the current_nodes Chebyshev-Gauss
        ! surfaces, interpolated (refit) and re-solved.
        real(dp), intent(out) :: error
        real(dp), parameter :: delta = 1e-3_dp
        real(dp) :: s(current_nodes), target(current_nodes), measured(current_nodes)
        real(dp) :: shifted(current_nodes), slope(current_nodes), iota_s(current_nodes)
        real(dp) :: scale, unused1, unused2
        type(toroidal_profile_t) :: previous
        integer :: j, k, local, ns
        ns = current_nodes
        do j = 1, ns
            s(j) = (1 - cos((2*j - 1)*pi/(2*ns)))/2
            if (current_type == 'hermite') then
                call hermite_profile_evaluate(profile_knots(:n_knots), &
                    current_values(:n_knots), current_slopes(:n_knots), s(j), target(j), &
                    unused1)
            else
                target(j) = sum(current*s(j)**[(k, k=0, kmax)])
            end if
        end do
        scale = maxval(abs(target))
        do k = 0, max_current_iterations
            call measure(s, measured)
            if (k == 0) then
                ! the current per unit iota of the starting map scales a
                ! current-free target, so that the error is an iota error there
                call seed%profile%shift_iota(delta)
                call measure(s, shifted)
                call seed%profile%shift_iota(-delta)
                scale = max(scale, maxval(abs(shifted - measured))/delta, tiny(1.0_dp))
            end if
            error = maxval(abs(measured - target))/scale
            print '(a,i3,a,es10.3)', 'current iteration ', k, ': relative error ', error
            if (error <= current_tolerance .or. k == max_current_iterations) then
                ! the current of a loosely solved map is measured again on the
                ! fully solved one
                if (solve_tolerance <= tolerance .or. k == max_current_iterations) exit
                solve_tolerance = tolerance
                call ideal_map_solve(map, solve_tolerance, max_iterations, report, &
                    solve_status, reference_gradient=reference)
                call measure(s, measured)
                error = maxval(abs(measured - target))/scale
                print '(a,i3,a,es10.3)', 'current iteration ', k, ': relative error ', error
                if (error <= current_tolerance) exit
            end if
            current_iterations = k + 1
            call seed%profile%shift_iota(delta)
            call measure(s, shifted)
            call seed%profile%shift_iota(-delta)
            slope = (shifted - measured)/delta
            do j = 1, ns
                call profile%evaluate(s(j), iota_s(j), unused1, unused2)
            end do
            ! Newton step on the solved map, where the current is affine in iota
            ! on each surface. The frozen map admits iota of one sign only
            ! (docs/ideal-seed-3d.md); the step may reverse that sign.
            iota_s = iota_s + (target - measured)/slope
            previous = profile
            call refit(s, iota_s)
            if (.not. single_signed()) then
                profile = previous
                iota_crosses_zero = .true.
                print '(a)', 'current iteration: the matched iota crosses zero, which the' &
                    //' frozen map does not admit'
                exit
            end if
            ! intermediate maps need only the accuracy of the current mismatch
            solve_tolerance = max(tolerance, min(1e-4_dp, 1e-3_dp*error))
            call resolve(local)
            if (local /= 0) call fail('the map folded when frozen with the updated iota')
        end do
    end subroutine

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

    subroutine measure(s, values)
        ! enclosed toroidal current on the surfaces s
        real(dp), intent(in) :: s(:)
        real(dp), intent(out) :: values(:)
        integer :: i, local
        do i = 1, size(s)
            values(i) = enclosed_current(s(i), local)
            if (local /= 0) call fail('the current readback met a folded map')
        end do
    end subroutine

    subroutine refit(s, values)
        ! iota as the polynomial interpolant of the matched values at the
        ! Chebyshev-Gauss surfaces s, which exclude the axis where the matched
        ! iota is only a limit; stored as values and slopes on 65 Hermite knots
        real(dp), intent(in) :: s(:), values(:)
        integer, parameter :: knots = 65
        real(dp) :: t(knots), value(knots), slope(knots), w(size(s)), dw(size(s))
        integer :: k
        do k = 1, knots
            t(k) = real(k - 1, dp)/(knots - 1)
            call lagrange_weights(size(s), t(k), s, w)
            call lagrange_deriv_weights(size(s), t(k), s, dw)
            value(k) = dot_product(w, values); slope(k) = dot_product(dw, values)
        end do
        call profile%set_hermite_iota(t, value, slope)
    end subroutine

    subroutine resolve(local)
        ! freeze the new iota on the solved map and solve again
        integer, intent(out) :: local
        seed%profile = profile
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

    real(dp) function enclosed_current(s, status)
        ! (1/mu0) oint B.dl around a poloidal loop on the flux surface s at phi = 0;
        ! every such loop encloses the same current. Inside, the loop is the image
        ! of the reference circle; at s = 1 it runs along the discrete wall.
        real(dp), intent(in) :: s
        integer, intent(out) :: status
        integer, parameter :: nt = 512, ng = 8
        real(dp) :: centre(2), radius, y(2), xi(2), trial(2), point(3, 0:nt), b(3, 0:nt)
        real(dp) :: values(EQ_VALUES), dl(3), corner(2, 3), x(ng), w(ng), q(3), d(3, 3)
        real(dp) :: ref(2), g(2, 2), tangent(2), along(3), sense
        integer :: k, e, local, owner, i, a
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        radius = sqrt(s)*(maxval(map%reference(1, :)) - minval(map%reference(1, :)))/2
        enclosed_current = 0; status = 1
        if (s >= 1) then
            corner = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
            call gauss_legendre(ng, x, w)
            x = (x + 1)/2; w = w/2
            do e = 1, map%cells
                do i = 1, 3
                    a = mod(i, 3) + 1
                    if (.not. (map%boundary(map%conn(i, e)) .and. map%boundary(map%conn(a, e)))) cycle
                    ! a wall edge, not a chord between two wall nodes
                    call map_3d_coordinates(map, e, (corner(:, i) + corner(:, a))/2, 0.0_dp, q, d, ref, g)
                    if (abs(norm2(ref - centre)/radius - 1) > 1e-3_dp) cycle
                    tangent = corner(:, a) - corner(:, i)
                    sense = sign(1.0_dp, (ref(1) - centre(1))*dot_product(g(2, :), tangent) &
                        - (ref(2) - centre(2))*dot_product(g(1, :), tangent))
                    do k = 1, ng
                        xi = corner(:, i) + x(k)*tangent
                        call map_3d_coordinates(map, e, xi, 0.0_dp, q, d, ref, g)
                        call map_3d_sample(map, seed, e, xi, 0.0_dp, point(:, 0), values, local)
                        if (local /= 0) return
                        along = matmul(d(:, 1:2), tangent); along(3) = along(3)*q(1)
                        enclosed_current = enclosed_current + sense*w(k)* &
                            dot_product(values([EQ_BR, EQ_BZ, EQ_BPHI]), along)
                    end do
                end do
            end do
        else
            do k = 0, nt
                y = centre + radius*[cos(2*pi*k/nt), sin(2*pi*k/nt)]
                owner = 0
                do e = 1, map%cells
                    call triangle_locate(map%reference(:, map%conn(:, e)), y, trial, local)
                    if (local /= 0) cycle
                    if (minval([trial, 1 - sum(trial)]) < -1e-10_dp) cycle
                    owner = e; xi = trial; exit
                end do
                if (owner == 0) return
                call map_3d_sample(map, seed, owner, xi, 0.0_dp, point(:, k), values, local)
                if (local /= 0) return
                b(:, k) = values([EQ_BR, EQ_BZ, EQ_BPHI])
            end do
            do k = 0, nt - 1
                ! (R, Z, phi) increments; the toroidal length is R dphi
                dl = point(:, k + 1) - point(:, k)
                dl(3) = dl(3)*(point(1, k + 1) + point(1, k))/2
                enclosed_current = enclosed_current + dot_product((b(:, k + 1) + b(:, k))/2, dl)
            end do
        end if
        enclosed_current = enclosed_current/mu0
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
