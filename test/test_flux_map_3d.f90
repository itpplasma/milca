program test_flux_map_3d
    !! Gates F2a and F2b of docs/flux-map.md: the Cartesian flux map in 3D
    !! (milca_flux_map_3d), unknown u = (x, y) per toroidal plane, no axis
    !! unknowns, the meshes the harmonic pullbacks of a disc O-grid
    !! (flux_map_sections).
    !!
    !! solovev (F2a): the Solov'ev case of gate F1 through the 3D path (3
    !! identical planes, start u = G). Compare the discrete energy and
    !! gradient against the axisymmetric operator on the same mesh, then
    !! check the field error against the exact solution.
    !!
    !! landreman (F2b): the exact sheared Landreman family of gates E4 and
    !! A1-3D (eps 0.6, 13 planes over the torus, shifted in phi so that no
    !! plane is a symmetry plane). The wall curves are the only geometric
    !! input; the first level starts from the normalised polar map of each
    !! section, the next from the previous level (stationary). The gauge is
    !! the moment constraint of milca_flux_map_3d; no symmetry row is used
    !! in either case.
    !! Passed when Newton converges, B against the exact field is
    !! within 20 % of the fixed-label fit of E4a on the exact-axis mesh
    !! (4.98e-2, 1.161e-2 at 2, 4 cells) and converges at order p; the axis
    !! (zero of u) is reported (it lies along an unstable direction).
    !!
    !! Both cases report, per level, the force balance J X B - grad p with
    !! J = curl B_h from the second derivatives of u_h (milca_flux_map_force),
    !! the jump of B^2/2 at the C0 block interfaces, the sign of grad phi .
    !! (grad x X grad y) at all quadrature points and the winding of the wall
    !! angle; one sign and winding one are checked, and at degree 2 that the
    !! force balance converges.
    !!
    !! Arguments: case, levels, planes, degree, toroidal quadrature points per
    !! plane (2), extra section quadrature points (4), start of the levels after
    !! the first: 'cold' (the polar start; default for solovev), 'nested'
    !! (the previous level, minimisation) or 'stationary' (the previous level,
    !! Newton on stationarity); cells of the first level (2); verbose. With the
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
    use milca_flux_map, only: flux_map_t, flux_map_init, flux_map_system
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_solve
    use flux_map_sections
    use milca_flux_map_force, only: flux_map_force_t, flux_map3d_force, flux_map_force_print, &
        flux_map_force_norm, FORCE_F, FORCE_GRAD_B2, FORCE_GRAD_P
    use solovev_label_fields, only: solovev_profiles, exact_b, boundary_samples, r0, alpha, &
        psi_b, chart, flux_scale
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
    real(dp) :: err(3), axis_err(3), force_err(3)
    ! landreman: the planes lie at phi = 2 pi (k - 1 + plane_shift)/planes, so
    ! that none of them is a symmetry plane of the family
    ! (MILCA_FM3_PLANE_SHIFT overrides it; 0: plane 1 up-down symmetric)
    real(dp) :: plane_shift = 0.37_dp
    ! F1 resolved-grid field accuracy and the fixed-label fit of E4a.
    real(dp), parameter :: f1_b(2:3) = [7.11e-3_dp, 1.93e-3_dp], &
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
        ! diagnostic: MILCA_FM3_GUIDE multiplies the toroidal field (R B_phi =
        ! 3 times it; the state stays exact)
        call get_environment_variable('MILCA_FM3_GUIDE', arg, level)
        if (level > 0) read (arg, *) flux_scale
        do level = 1, levels
            call run_solovev(cells(level), err(level), axis_err(level))
            if (degree == 2 .and. flux_scale == 1 .and. level > 1) &
                call check(err(level) <= 1.02_dp*f1_b(level), 'B error meets the resolved F1 accuracy')
        end do
    case ('landreman')
        if (levels == 0) levels = 2
        if (nt == 0) nt = 13
        section_oracle%family = NESTED_SHEARED
        section_oracle%eps = 0.6_dp
        section_oracle%saxis = 2.4_dp
        section_oracle%delta = 0.2_dp
        call section_fit_profiles()
        call get_environment_variable('MILCA_FM3_PLANE_SHIFT', arg, level)
        if (level > 0) read (arg, *) plane_shift
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
        print '(a,3f7.2)', 'rates of B, axis, force balance', &
            log(err(level - 1)/err(level))/log(2.0_dp), &
            log(axis_err(level - 1)/axis_err(level))/log(2.0_dp), &
            log(force_err(level - 1)/force_err(level))/log(2.0_dp)
        call check(log(err(level - 1)/err(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'B converges at order p')
        if (trim(which) == 'solovev') &
            call check(axis_err(level) < axis_err(level - 1), 'the axis converges')
        ! J = curl B_h has the order p - 1; measured at degree 2: Solovev
        ! 0.68, 0.79 (0.93, 0.98 on 16, 32 cells), Landreman 1.05
        if (degree == 2 .and. level <= 2) &
            call check(log(force_err(level - 1)/force_err(level))/log(2.0_dp) > 0.5_dp, &
                       'the force balance |J X B - grad p|/|grad p| converges')
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
        b(2) = flux_scale*b(2)
    end subroutine solovev_b

    subroutine landreman_fm(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: f(3)

        call section_profiles(s, f, gp, p)
        tau = -2*f
        p = pressure_scale()*p
    end subroutine landreman_fm

    subroutine landreman_b(x, b)
        !! The exact field at the toroidal angle of the shifted planes.
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)

        call section_exact_b([x(1), x(2) + 2*pi*plane_shift/nt, x(3)], b)
    end subroutine landreman_b

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
        call fm3_build(c, samples, n, degree, 4, st)
        if (st /= 0) error stop 'meshes'
        call fm3_start(c)
        if (n == cells(1)) call check_axisymmetric_reduction(c)
        c%fm%verbose = verbose
        call flux_map3d_solve(c%fm, c%derham, c%quad, solovev_fm, c%base, 1e-9_dp, 60, st)
        call check(st == 0, 'Newton converges')
        print '(a,i0,a,i0,a,l1)', '  Hessian inertia: negative ', c%fm%negative_eigenvalues, &
            ', zero ', c%fm%zero_eigenvalues, '  reliable ', c%fm%inertia_reliable
        ! the Solovev states are saddles of W too: negative eigenvalues of the
        ! dense spectrum of the dumped Hessian (hessian_spectrum.py): 2, 2, 4
        ! at 2, 4, 8 cells
        if (degree == 2) &
            call check(c%fm%inertia_reliable .and. c%fm%zero_eigenvalues == 0 .and. &
                   (n /= 2 .or. c%fm%negative_eigenvalues == 2) .and. &
                   (n /= 4 .or. c%fm%negative_eigenvalues == 2) .and. &
                   (n /= 8 .or. c%fm%negative_eigenvalues == 4), &
                   'Hessian inertia equals the dense spectrum of the dumped Hessian')
        call flux_map3d_system(c%fm, c%derham, c%quad, solovev_fm, c%base, c%fm%z, w, &
                               exact_field=solovev_b, field_error=ferr)
        e = sqrt(ferr(1)/ferr(2))
        call report_force(c, solovev_fm)
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

    subroutine check_axisymmetric_reduction(c)
        type(fm3_case_t), intent(in) :: c
        type(flux_map_t) :: poloidal
        real(dp), allocatable :: grad_2d(:), grad_3d(:), reduced(:)
        real(dp) :: energy_2d, energy_3d
        integer :: i, n, nt, ntot

        n = c%n0
        nt = c%nt
        ntot = n*nt
        call flux_map_init(poloidal, n, pack([(i, i=1, n)], c%fm%is_wall))
        poloidal%z = [c%fm%z(:n), c%fm%z(ntot + 1:ntot + n)]
        call flux_map_system(poloidal, c%derham, c%quad, solovev_fm, poloidal%z, energy_2d, grad_2d)
        call flux_map3d_system(c%fm, c%derham, c%quad, solovev_fm, c%base, c%fm%z, energy_3d, grad_3d)
        reduced = [sum(reshape(grad_3d(:ntot), [n, nt]), dim=2), &
                   sum(reshape(grad_3d(ntot + 1:), [n, nt]), dim=2)]
        print '(a,2es11.3)', '  axisymmetric reduction: energy and gradient differences ', &
            abs(energy_3d - energy_2d), norm2(reduced - grad_2d)
        call check(abs(energy_3d - energy_2d) < 1e-12_dp*max(abs(energy_2d), 1.0_dp), &
                   '3D and axisymmetric energies on the same mesh agree')
        call check(norm2(reduced - grad_2d) < 1e-11_dp*max(norm2(grad_2d), 1.0_dp), &
                   '3D gradient restricted to axisymmetry matches the axisymmetric gradient')
    end subroutine check_axisymmetric_reduction

    subroutine report_force(c, profiles)
        !! Force balance J X B - grad p of the state, J = curl B from the
        !! second derivatives of u_h (milca_flux_map_force), the sign of the
        !! Jacobian and the winding of the wall angle.
        type(fm3_case_t), intent(in) :: c
        procedure(solovev_fm) :: profiles
        type(flux_map_force_t) :: force
        real(dp) :: winding(c%nt), least_step(c%nt), jump
        integer :: st

        call flux_map3d_force(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, 10, force, st)
        if (st /= 0) error stop 'force balance'
        call flux_map_force_print(force, 'CEMV')
        force_err(level) = flux_map_force_norm(force, FORCE_F, 0)/ &
                           flux_map_force_norm(force, FORCE_GRAD_P, 0)
        call fm3_wall_winding(c, profiles, winding, least_step, st)
        if (st /= 0) error stop 'wall winding'
        call fm3_interface_jump(c, profiles, jump, st)
        if (st /= 0) error stop 'interface jump'
        print '(a,2f10.6,a,es10.2,a,es10.2)', '  winding of the wall angle per plane, min max', &
            minval(winding), maxval(winding), '  least step of the wall angle / mean step', &
            minval(least_step), '  jump of B^2/2 at the block interfaces (relative rms)', jump
        call check(force%jacobian_min*force%jacobian_max > 0, &
                   'grad phi . (grad x X grad y) has one sign at all quadrature points')
        call check(all(abs(winding - 1) < 1e-8_dp), 'the wall angle winds once')
    end subroutine report_force

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
            plane_phi = 2*pi*(k - 1 + plane_shift)/np
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
        ! no plane is a symmetry plane (plane_shift), and no symmetry is used
        sym = maxval(abs(samples(2, 2:nsamples/2, 1) + samples(2, nsamples:nsamples/2 + 2:-1, 1)))
        print '(a,i0,a,es9.2,a,f6.3,a,f7.3,a,l1)', '  planes ', np, '  plane 1 asymmetry', sym, &
            '  worst Jacobian ratio', minval(quality), '  mesh centre off the axis / minor radius', &
            maxval(norm2(c%centre - a, 1))/minor, '  gauge constraint ', c%fm%gauge_constraint
        c%fm%verbose = verbose
    end subroutine landreman_case

    subroutine run_landreman(n, deg, e, a_err)
        !! The first level from the polar start (fm3_start_polar; the start
        !! u = G goes to the state of lower energy). The registered gate
        !! starts each level independently from the polar start. With
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
        if (n == 4 .and. deg == 2 .and. nt == 13 .and. (stationary .or. .not. mode_given)) &
            call check(c%fm%negative_eigenvalues == 6, &
                       'Hessian of the 4-cell, 13-plane state has 6 negative eigenvalues')
        call flux_map3d_system(c%fm, c%derham, c%quad, landreman_fm, c%base, c%fm%z, w, &
                               exact_field=landreman_b, field_error=ferr)
        e = sqrt(ferr(1)/ferr(2))
        call report_force(c, landreman_fm)
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
