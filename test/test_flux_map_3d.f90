program test_flux_map_3d
    !! Gates F2a and F2b of docs/flux-map.md: the Cartesian flux map in 3D
    !! (milca_flux_map_3d), unknown u = (x, y) per toroidal plane, no axis
    !! unknowns, the meshes the harmonic pullbacks of a disc O-grid
    !! (flux_map_sections).
    !!
    !! solovev (F2a): the Solov'ev case of gate F1 through the 3D path (3
    !! identical planes, start u = G). Passed when Newton converges and B
    !! agrees with the F1 numbers within 2 %.
    !!
    !! landreman (F2b): the exact sheared Landreman family of gates E4 and
    !! A1-3D (eps 0.6, 13 planes over the torus). The wall curves are the
    !! only geometric input; the first level starts from the normalised polar
    !! map of each section, the next from the previous level (stationary).
    !! Passed when Newton converges, B against the exact field is
    !! within 20 % of the fixed-label fit of E4a on the exact-axis mesh
    !! (4.98e-2, 1.161e-2 at 2, 4 cells) and converges at order p; the axis
    !! (zero of u) is reported (it lies along an unstable direction).
    !!
    !! Arguments: case, levels, planes, degree, toroidal quadrature points per
    !! plane (2), extra section quadrature points (4), start of the levels after
    !! the first: 'cold' (the polar start; default for solovev), 'nested'
    !! (the previous level, minimisation) or 'stationary' (the previous level,
    !! Newton on stationarity; default for landreman, whose exact state is a
    !! saddle); cells of the first level (2); verbose. With the
    !! environment variable MILCA_FM3_HESSIAN=prefix every landreman level
    !! writes its Hessian at the solution (fm3_dump_hessian)
    !! (benchmarks/flux-map-spike/hessian_spectrum.py). With
    !! MILCA_FM3_START_PLANES=m a landreman run first solves the first level
    !! on m planes from the polar start and starts its first level from that
    !! state (fm3_prolong, trigonometric in t) with the chosen start mode.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use milca_nested_oracles, only: NESTED_SHEARED
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_solve
    use flux_map_sections
    use solovev_label_fields, only: solovev_profiles, exact_b, boundary_samples, r0, alpha, &
        psi_b, chart
    use landreman_sections, only: section_oracle, section_exact_b, section_axis, &
        section_boundary_radius, section_fit_profiles, section_profiles
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: degree = 2, levels = 0, nt = 0, failures = 0, level, tpoints = 2, extra = 4, first = 2
    logical :: cold = .true., stationary = .false., mode_given = .false.
    type(fm3_case_t), target :: previous
    integer :: cells(3) = [2, 4, 8], start_planes = 0, planes
    logical :: verbose = .false.
    character(16) :: which = 'solovev', arg
    real(dp) :: err(3), axis_err(3)
    ! F1 (pointwise B_phi, R B_phi = 3) and the fixed-label fit of E4a
    real(dp), parameter :: f1_b(3) = [2.42e-2_dp, 7.11e-3_dp, 1.93e-3_dp], &
                           e4a_fit(2) = [4.98e-2_dp, 1.161e-2_dp]

    if (command_argument_count() >= 1) call get_command_argument(1, which)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) degree
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) tpoints
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read (arg, *) extra
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        cold = trim(arg) == 'cold'
        stationary = trim(arg) == 'stationary'
        mode_given = .true.
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, arg)
        read (arg, *) first
        cells = first*[1, 2, 4]
    end if
    verbose = command_argument_count() >= 9
    clebsch_nfp = 1
    select case (trim(which))
    case ('solovev')
        if (levels == 0) levels = 3
        if (nt == 0) nt = 3
        call curve_chart_cylindrical(chart, 1)
        do level = 1, levels
            call run_solovev(cells(level), err(level), axis_err(level))
            if (degree == 2) call check(abs(err(level)/f1_b(level) - 1) < 0.02_dp, &
                                        'B error equals gate F1 within 2 %')
        end do
    case ('landreman')
        if (levels == 0) levels = 2
        if (nt == 0) nt = 13
        ! the exact state is a saddle (pressure-driven, docs/flux-map.md):
        ! the levels after the first seek the stationary point next to the
        ! previous one
        if (.not. mode_given) then
            cold = .false.
            stationary = .true.
        end if
        section_oracle%family = NESTED_SHEARED
        section_oracle%eps = 0.6_dp
        section_oracle%saxis = 2.4_dp
        section_oracle%delta = 0.2_dp
        call section_fit_profiles()
        call get_environment_variable('MILCA_FM3_START_PLANES', arg, level)
        if (level > 0) read (arg, *) start_planes
        if (start_planes > 0) then
            planes = nt
            nt = start_planes
            call run_landreman(cells(1), degree, err(1), axis_err(1))
            nt = planes
        end if
        if (degree > 2) call run_landreman(cells(1), 2, err(1), axis_err(1))
        do level = 1, levels
            call run_landreman(cells(level), degree, err(level), axis_err(level))
            if (degree == 2 .and. level <= 2) call check(err(level) < 1.2_dp*e4a_fit(level), &
                'B within 20 % of the fixed-label fit of E4a')
        end do
    case default
        error stop 'unknown case'
    end select
    do level = 2, levels
        print '(a,2f7.2)', 'rates of B, axis', log(err(level - 1)/err(level))/log(2.0_dp), &
            log(axis_err(level - 1)/axis_err(level))/log(2.0_dp)
        call check(log(err(level - 1)/err(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'B converges at order p')
        if (trim(which) == 'solovev') &
            call check(axis_err(level) < axis_err(level - 1), 'the axis converges')
    end do
    if (failures > 0) error stop 'flux map 3D gate failed'
    print '(2a)', 'PASS flux map 3D gate ', trim(which)

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            print '(2a)', 'PASS ', what
        else
            print '(2a)', 'FAIL ', what
            failures = failures + 1
        end if
    end subroutine check

    subroutine solovev_fm(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: f(3)

        call solovev_profiles(s, f, gp, p)
        tau = -2*f
        p = pressure_scale()*p
    end subroutine solovev_fm

    subroutine solovev_b(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)

        b = exact_b([x(1), 0.0_dp, x(3)])
    end subroutine solovev_b

    subroutine landreman_fm(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: f(3)

        call section_profiles(s, f, gp, p)
        tau = -2*f
        p = pressure_scale()*p
    end subroutine landreman_fm

    real(dp) function pressure_scale()
        !! diagnostic: MILCA_FM3_PSCALE multiplies the pressure (the state is
        !! then no longer the exact one)
        character(32) :: v
        integer :: length
        real(dp), save :: scale = -huge(1.0_dp)
        if (scale == -huge(1.0_dp)) then
            scale = 1
            call get_environment_variable('MILCA_FM3_PSCALE', v, length)
            if (length > 0) read (v, *) scale
        end if
        pressure_scale = scale
    end function pressure_scale

    subroutine run_solovev(n, e, a_err)
        integer, intent(in) :: n
        real(dp), intent(out) :: e, a_err
        type(fm3_case_t) :: c
        real(dp), allocatable :: one(:, :), samples(:, :, :)
        real(dp) :: w, ferr(2), axis(2), t0, t1
        integer :: st, k, clock0, clock1, rate

        call system_clock(clock0, rate)
        call boundary_samples(one)
        allocate (samples(2, size(one, 2), nt))
        do k = 1, nt
            samples(:, :, k) = one
        end do
        ! the mesh of gate F1: default core, centre between the wall radii
        call fm3_build(c, samples, n, degree, 4, st)
        if (st /= 0) error stop 'meshes'
        call fm3_start(c)
        call fm3_gauge_symmetric(c, st)
        if (st /= 0) error stop 'gauge rows'
        c%fm%verbose = verbose
        call flux_map3d_solve(c%fm, c%derham, c%quad, solovev_fm, c%base, 1e-9_dp, 60, st)
        call check(st == 0, 'Newton converges')
        print '(a,i0,a,i0,a,l1)', '  Hessian inertia: negative ', c%fm%negative_eigenvalues, &
            ', zero ', c%fm%zero_eigenvalues, '  reliable ', c%fm%inertia_reliable
        ! the Solovev states are saddles of W too: negative eigenvalues of the
        ! dense spectrum of the dumped Hessian (hessian_spectrum.py): 2, 2, 4
        ! at 2, 4, 8 cells
        call check(c%fm%inertia_reliable .and. c%fm%zero_eigenvalues == 0 .and. &
                   (n /= 2 .or. c%fm%negative_eigenvalues == 2) .and. &
                   (n /= 4 .or. c%fm%negative_eigenvalues == 2) .and. &
                   (n /= 8 .or. c%fm%negative_eigenvalues == 4), &
                   'Hessian inertia equals the dense spectrum of the dumped Hessian')
        call flux_map3d_system(c%fm, c%derham, c%quad, solovev_fm, c%base, c%fm%z, w, &
                               exact_field=solovev_b, field_error=ferr)
        e = sqrt(ferr(1)/ferr(2))
        a_err = 0
        do k = 1, nt
            call fm3_axis(c, (k - 1)/real(nt, dp), solovev_fm, axis, st)
            if (st /= 0) error stop 'axis'
            a_err = max(a_err, norm2(axis - [r0, 0.0_dp]))
        end do
        call fm3_dump_hessian(c, solovev_fm)
        call system_clock(clock1)
        t0 = 0
        t1 = real(clock1 - clock0, dp)/rate
        print '(a,i2,a,i0,a,i3,a,es10.3,a,es10.3,a,es16.9,a,f7.2,a)', 'cells', n, '  unknowns ', &
            2*c%n0*nt, '  Newton', c%fm%steps, '  B error', e, '  axis', a_err, '  W', w, &
            '  (', t1 - t0, ' s)'
    end subroutine run_solovev

    subroutine landreman_case(n, deg, np, c, a, minor)
        !! Meshes of np planes from the wall curves (rays about the exact
        !! axis a, which enters only as their origin) and the polar start.
        integer, intent(in) :: n, deg, np
        type(fm3_case_t), intent(out) :: c
        real(dp), allocatable, intent(out) :: a(:, :)
        real(dp), intent(out) :: minor
        integer, parameter :: nsamples = 128
        real(dp) :: samples(2, nsamples, np), quality(np), tk, plane_phi, sym
        integer :: st, k, i

        allocate (a(2, np))
        do k = 1, np
            plane_phi = 2*pi*(k - 1)/np
            a(:, k) = section_axis(plane_phi)
            do i = 1, nsamples
                tk = 2*pi*(i - 1)/nsamples
                samples(:, i, k) = a(:, k) + section_boundary_radius(plane_phi, a(:, k), tk)* &
                                   [cos(tk), sin(tk)]
            end do
        end do
        minor = 0.5_dp*(maxval(samples(1, :, 1)) - minval(samples(1, :, 1)))
        call fm3_build(c, samples, n, deg, extra, st)
        if (st /= 0) error stop 'meshes'
        c%fm%toroidal_points = tpoints
        call fm3_mesh_quality(c, quality)
        call fm3_start_polar(c)
        ! the gauge rows need plane 1 up-down symmetric
        sym = maxval(abs(samples(2, 2:nsamples/2, 1) + samples(2, nsamples:nsamples/2 + 2:-1, 1)))
        if (sym < 1e-8_dp*minor .and. abs(a(2, 1)) < 1e-10_dp) then
            call fm3_gauge_symmetric(c, st)
            if (st /= 0) error stop 'gauge rows'
        end if
        print '(a,i0,a,es9.2,a,f6.3,a,f7.3,a,i0)', '  planes ', np, '  plane 1 asymmetry', sym, &
            '  worst Jacobian ratio', minval(quality), '  mesh centre off the axis / minor radius', &
            maxval(norm2(c%centre - a, 1))/minor, '  gauge rows ', size(c%fm%gauge_target)
        c%fm%verbose = verbose
    end subroutine landreman_case

    subroutine run_landreman(n, deg, e, a_err)
        !! The first level from the polar start (fm3_start_polar; the start
        !! u = G goes to the state of lower energy). Default (the gate): a
        !! further level seeks the stationary point next to the previous
        !! solution; with 'cold' it starts from the polar start, with
        !! 'nested' it starts from the previous solution (fm3_prolong)
        !! and minimises; with 'stationary' it seeks the stationary point next
        !! to it. At 8 cells, at degree 3 or on 21 planes the minimisation
        !! leaves the exact state for one of lower energy (docs/flux-map.md).
        integer, intent(in) :: n, deg
        real(dp), intent(out) :: e, a_err
        type(fm3_case_t) :: c
        real(dp), allocatable :: a(:, :), flux(:)
        real(dp) :: w, ferr(2), axis(2), t1, minor, bmin, bmax
        integer :: st, k, clock0, clock1, rate

        call system_clock(clock0, rate)
        call landreman_case(n, deg, nt, c, a, minor)
        if (previous%n0 > 0 .and. .not. cold) then
            call fm3_prolong(previous, c, st)
            if (st /= 0) error stop 'prolongation'
            c%fm%stationary = stationary
        end if
        call flux_map3d_solve(c%fm, c%derham, c%quad, landreman_fm, c%base, 1e-9_dp, 60, st)
        call check(st == 0, 'Newton converges')
        call check(c%fm%negative_eigenvalues >= 0 .and. c%fm%inertia_reliable, &
                   'Hessian inertia of the converged state is reliable')
        ! F2b: the exact state is a saddle, six negative eigenvalues (dense
        ! spectrum of the dumped Hessian, hessian_spectrum.py)
        if (n == 4 .and. deg == 2 .and. nt == 13 .and. stationary) &
            call check(c%fm%negative_eigenvalues == 6, &
                       'Hessian of the 4-cell, 13-plane state has 6 negative eigenvalues')
        call flux_map3d_system(c%fm, c%derham, c%quad, landreman_fm, c%base, c%fm%z, w, &
                               exact_field=section_exact_b, field_error=ferr)
        e = sqrt(ferr(1)/ferr(2))
        allocate (flux(nt))
        call fm3_flux(c, landreman_fm, flux, bmin, bmax)
        print '(a,2es13.5,a,2es11.3)', '  toroidal flux through the planes, min max', &
            minval(flux), maxval(flux), '  B_phi range', bmin, bmax
        a_err = 0
        do k = 1, nt
            call fm3_axis(c, (k - 1)/real(nt, dp), landreman_fm, axis, st)
            if (st /= 0) error stop 'axis'
            a_err = max(a_err, norm2(axis - a(:, k)))
            if (verbose) print '(a,i3,2es13.5)', '  axis minus exact axis, plane', k, axis - a(:, k)
        end do
        call system_clock(clock1)
        t1 = real(clock1 - clock0, dp)/rate
        print '(a,i2,a,i0,a,i3,a,es10.3,a,es10.3,a,es16.9,a,f7.2,a)', 'cells', n, '  unknowns ', &
            2*c%n0*nt, '  Newton', c%fm%steps, '  B error', e, '  axis', a_err, '  W', w, &
            '  (', t1, ' s)'
        print '(a,i0,a,i0,a,l1)', '  Hessian inertia: negative eigenvalues ', &
            c%fm%negative_eigenvalues, ', zero ', c%fm%zero_eigenvalues, '  reliable ', &
            c%fm%inertia_reliable
        call fm3_dump_hessian(c, landreman_fm)
        previous = c
    end subroutine run_landreman

end program test_flux_map_3d
