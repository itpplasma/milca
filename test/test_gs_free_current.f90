program test_gs_free_current
    !! Contour integrals of the current-profile mode (milca_gs_free_current)
    !! on psi = psi0 - c ((R - R0)^2 + Z^2), which P2 elements hold exactly:
    !! on the circle of radius r, a1 = oint R dl/|grad psi| = pi R0/c and
    !! a2 = oint dl/(R |grad psi|) = pi/(c sqrt(R0^2 - r^2)).
    !! The profile update on the same field under a constant law (p' = P,
    !! FF' = G) with the exact p(s), I(s) and Phi_edge of that law as target
    !! must return the law and the boundary level.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_grad_shafranov, only: gs_state_t, initialize_gs
    use milca_gs_free_boundary, only: gs_fb_surface_integrals, GS_FB_OK
    use milca_gs_profile_law, only: gs_profile_law_t
    use milca_gs_free_current, only: gs_current_target_t, gs_current_update, &
        gs_current_profiles, GS_CUR_OK
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 1.7_dp, c = 2.0_dp, psi0 = 1.0_dp
    type(gs_state_t) :: state
    real(dp) :: contour(2, 64), radii(4), levels(4), a1(4), a2(4), e1, e2
    integer :: k, status

    do k = 1, size(contour, 2)
        contour(:, k) = [r0 + 0.6_dp*cos(2*pi*(k - 1)/size(contour, 2)), &
                         0.6_dp*sin(2*pi*(k - 1)/size(contour, 2))]
    end do
    call initialize_gs(state, contour, 0.004_dp, status, 2)
    if (status /= 0) error stop 'mesh'
    state%psi = psi0 - c*((state%points(1, :) - r0)**2 + state%points(2, :)**2)
    radii = [0.05_dp, 0.2_dp, 0.4_dp, 0.55_dp]
    levels = psi0 - c*radii**2
    call gs_fb_surface_integrals(state, state%psi, [r0, 0.0_dp], levels, 1.0_dp, &
                                 256, a1, a2, status)
    if (status /= GS_FB_OK) error stop 'surface integrals failed'
    e1 = maxval(abs(a1 - pi*r0/c))/(pi*r0/c)
    e2 = maxval(abs(a2 - pi/(c*sqrt(r0**2 - radii**2)))/(pi/(c*sqrt(r0**2 - radii**2))))
    print '(a,2es10.2)', 'relative errors a1 a2', e1, e2
    if (e1 > 1.0e-7_dp .or. e2 > 1.0e-7_dp) error stop 'contour integrals'
    call update_fixed_point()
    print '(a)', 'gs_free_current contour integrals and update pass'
contains
    subroutine update_fixed_point()
        real(dp), parameter :: mu0 = 4.0e-7_dp*pi, p_coef = 2.0e4_dp, g_coef = -0.5_dp, &
            f_edge = 3.0_dp, p_edge = 500.0_dp, r_x = 0.58_dp, r_b = 0.45_dp
        integer, parameter :: nk = 21, nq = 4001
        type(gs_profile_law_t) :: law, next
        type(gs_current_target_t) :: target, made
        real(dp) :: ub, ux, cut, next_cut, phi_deep, uq(nq), fa2(nq), di(nq), phi(nq), &
            current(nq), f, a1u, a2u, err_p, err_f, err_c
        integer :: j, st

        ub = c*r_b**2
        ux = c*r_x**2
        cut = ub/ux
        allocate (law%knots(nk), law%pp(4, nk - 1), law%ff(4, nk - 1))
        law%knots = [(real(j, dp)/(nk - 1), j=0, nk - 1)]
        law%pp = 0; law%pp(1, :) = p_coef
        law%ff = 0; law%ff(1, :) = g_coef
        law%f_edge = f_edge; law%p_edge = p_edge; law%normalized = .true.
        law%axis = ub
        ! Exact target by the trapezoidal rule on nq points in u (h^2 = 4e-9 of u_b^2)
        do j = 1, nq
            uq(j) = ub*real(j - 1, dp)/(nq - 1)
            f = sqrt(f_edge**2 + 2*g_coef*(ub - uq(j)))
            a1u = pi*r0/c
            a2u = pi/(c*sqrt(r0**2 - uq(j)/c))
            fa2(j) = f*a2u
            di(j) = a1u*p_coef + a2u*g_coef/mu0
        end do
        phi(1) = 0; current(1) = 0
        do j = 2, nq
            phi(j) = phi(j - 1) + 0.5_dp*(uq(j) - uq(j - 1))*(fa2(j) + fa2(j - 1))
            current(j) = current(j - 1) + 0.5_dp*(uq(j) - uq(j - 1))*(di(j) + di(j - 1))
        end do
        target%phi_edge = phi(nq)
        target%s = phi/phi(nq)
        target%p = p_edge + p_coef*(ub - uq)
        target%dp = -p_coef/(fa2/phi(nq))
        target%current = current
        target%dcurrent = di/(fa2/phi(nq))
        call gs_current_update(state, state%psi, [r0, 0.0_dp], psi0, psi0 - ub, cut, &
            1.0_dp, mu0, law, target, 256, 0.995_dp, next, next_cut, phi_deep, st)
        if (st /= GS_CUR_OK) error stop 'update failed'
        err_p = maxval(abs(next%pp(1, :) - p_coef))/p_coef
        err_f = maxval(abs(next%ff(1, :) - g_coef))/abs(g_coef)
        err_c = abs(next_cut - cut)/cut
        print '(a,3es10.2)', 'update errors pp ff cut', err_p, err_f, err_c
        if (max(err_p, err_f, err_c) > 1.0e-6_dp) error stop 'update is not a fixed point'
        ! The inverse: profiles of the field under the law against the exact ones
        call gs_current_profiles(state, state%psi, [r0, 0.0_dp], psi0, psi0 - ub, 1.0_dp, &
            mu0, law, 200, 256, made, st)
        if (st /= GS_CUR_OK) error stop 'profiles failed'
        err_c = abs(made%phi_edge - target%phi_edge)/target%phi_edge
        err_p = 0; err_f = 0
        do j = 1, size(made%s)
            err_p = max(err_p, abs(made%p(j) - interp(target%s, target%p, made%s(j))) &
                /target%p(1))
            err_f = max(err_f, abs(made%current(j) - interp(target%s, target%current, &
                made%s(j)))/abs(target%current(nq)))
        end do
        print '(a,3es10.2)', 'profile errors Phi_edge p I', err_c, err_p, err_f
        if (max(err_p, err_f, err_c) > 1.0e-6_dp) error stop 'profiles of the field'
    end subroutine

    real(dp) function interp(x, y, at) result(v)
        real(dp), intent(in) :: x(:), y(:), at
        integer :: i
        i = min(max(count(x < at), 1), size(x) - 1)
        v = y(i) + (y(i + 1) - y(i))*(at - x(i))/(x(i + 1) - x(i))
    end function
end program test_gs_free_current
