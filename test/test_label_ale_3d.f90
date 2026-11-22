program test_label_ale_3d
    !! Gate A1-3D of docs/equilibrium-roadmap.md: the 3D label lane on the
    !! ALE mesh (milca_label_ale_3d) on the exact sheared Landreman family of
    !! gate E4 (eps 0.6, two field periods, magnetic shear). Unlike E4 no
    !! exact information enters the solve except the profiles: the planes are
    !! meshed with their cores on an axis guess offset helically from the
    !! exact axis, the labels start from the geometric guess (|x - X| /
    !! rho_wall)^2, the winding frames come from the discrete labels, and the
    !! X columns of Newton are exact (shape derivatives).
    !!
    !! Passed when the solve converges in at most 4 passes, B is within 20 %
    !! of the fixed-label fit on the mesh of the exact axis (E4a, the
    !! interpolation level) and converges at order p, and the axis converges.
    !!
    !! Arguments: levels (default 2), planes (13), axis offset in units of
    !! the minor radius (0.05).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, spline_derham_project
    use milca_nested_oracles, only: NESTED_SHEARED
    use milca_spline_ale, only: spline_set_geometry
    use milca_label_ale_3d
    use landreman_sections, only: section_oracle, section_label, section_exact_b, &
        section_axis, section_boundary_radius, section_fit_profiles, section_profiles
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: degree = 2, nsamples = 128
    integer :: nt = 13, levels = 2, level, failures
    integer :: cells(3) = [2, 4, 8]
    real(dp) :: offset = 0.05_dp, plane_phi
    integer :: max_passes = 4, max_steps = 16
    real(dp), allocatable :: err(:), err_fixed(:), axis_err(:)
    type(label3d_t), target :: probe
    character(16) :: arg

    failures = 0
    section_oracle%family = NESTED_SHEARED
    section_oracle%eps = 0.6_dp
    section_oracle%saxis = 2.4_dp
    section_oracle%delta = 0.2_dp
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) offset
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(4, arg)
        read (arg, *) max_passes
        call get_command_argument(5, arg)
        read (arg, *) max_steps
    end if
    call section_fit_profiles()
    allocate (err(levels), err_fixed(levels), axis_err(levels))
    do level = 1, levels
        call run(cells(level), err(level), err_fixed(level), axis_err(level))
        call check(err(level) < 1.2_dp*err_fixed(level), &
                   'B within 20 % of the fixed-label fit on the exact-axis mesh')
    end do
    do level = 2, levels
        print '(a,2f7.2)', 'rates of B, axis', log(err(level - 1)/err(level))/log(2.0_dp), &
            log(axis_err(level - 1)/axis_err(level))/log(2.0_dp)
        call check(log(err(level - 1)/err(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'B converges at order p')
        call check(axis_err(level) < axis_err(level - 1), 'the axis converges')
    end do
    if (failures > 0) error stop 'gate A1-3D failed'
    print '(a)', 'PASS gate A1-3D'

contains

    subroutine run(n, e, e_fixed, a_err)
        integer, intent(in) :: n
        real(dp), intent(out) :: e, e_fixed, a_err
        type(label3d_t) :: s, ref
        real(dp) :: samples(2, nsamples, nt), a(2, nt), x0(2, nt), minor, t0, t1, tk
        integer :: k, st, i

        call cpu_time(t0)
        do k = 1, nt
            plane_phi = 2*pi*(k - 1)/nt
            a(:, k) = section_axis(plane_phi)
            do i = 1, nsamples
                tk = 2*pi*(i - 1)/nsamples
                samples(:, i, k) = a(:, k) + section_boundary_radius(plane_phi, a(:, k), tk)* &
                                   [cos(tk), sin(tk)]
            end do
        end do
        minor = 0.5_dp*(maxval(samples(1, :, 1)) - minval(samples(1, :, 1)))
        do k = 1, nt
            tk = (k - 1)/real(nt, dp)
            x0(:, k) = a(:, k) + offset*minor*[cos(2*pi*tk), sin(2*pi*tk)]
        end do

        ! the interpolation level: exact labels on the mesh of the exact axis
        call label3d_init(ref, samples, a, n, st)
        if (st /= 0) error stop 'reference mesh'
        call exact_labels(ref)
        call label3d_frames(ref)
        call label3d_solve_g(ref, section_profiles, st)
        if (st /= 0) error stop 'reference g'
        e_fixed = label3d_field_error(ref, section_profiles, section_exact_b)

        ! the solve from the axis guess and geometric labels
        call label3d_init(s, samples, x0, n, st)
        if (st /= 0) error stop 'mesh'
        call label3d_initial_labels(s)
        call label3d_frames(s)
        call label3d_solve_g(s, section_profiles, st)
        if (st /= 0) error stop 'initial g'
        s%verbose = .true.
        call label3d_solve(s, section_profiles, 1e-9_dp, 0.05_dp, max_passes, max_steps, st)
        call cpu_time(t1)
        e = label3d_field_error(s, section_profiles, section_exact_b)
        a_err = maxval(norm2(s%x - a, 1))
        print '(a,i2,a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,es10.3,a,f7.2,a)', 'cells', n, &
            '  status ', st, '  passes ', s%passes, '  Newton ', s%newton_total, &
            '  B error', e, '  fixed-label fit', e_fixed, '  axis', a_err, '  (', t1 - t0, ' s)'
        call check(st == 0, 'ALE loop converges')
        call check(s%passes <= 4, 'at most 4 passes')
    end subroutine run

    subroutine exact_labels(s)
        !! chi = Pi0 of the exact label on every plane.
        type(label3d_t), intent(inout) :: s
        type(spline_derham_t) :: d
        integer :: k, n0, st
        n0 = size(s%w)
        if (allocated(s%chi)) deallocate (s%chi)
        allocate (s%chi(n0*nt))
        do k = 1, nt
            plane_phi = 2*pi*(k - 1)/nt
            call spline_set_geometry(s%space, s%base(:, :, k))
            call spline_derham_build(s%space, d, st, 0)
            if (st /= 0) error stop 'plane derham'
            probe%space = s%space
            call spline_derham_project(d, 0, label_cb, s%chi((k - 1)*n0 + 1:k*n0), st, 12)
            if (st /= 0) error stop 'label projection'
        end do
        call spline_set_geometry(s%space, s%base(:, :, 1))
        s%chi(s%wall) = 1
    end subroutine exact_labels

    subroutine label_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st
        call spline_geometry_at(probe%space, block, xi, jet, st)
        values(1) = section_label(jet(1, 1), plane_phi, jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

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

end program test_label_ale_3d
