module milca_gs_spline_current
    !! Current-profile update for the free-boundary Grad-Shafranov solve on
    !! splines (VMEC's p(s), I(s), Phi_edge; s = Phi/Phi_edge). Same equations
    !! as milca_gs_free_current for the triangle solver:
    !!     dPhi/du = F a2,   dI/du = a1 p' + a2 FF'/mu0,
    !! with u = sigma (psi_a - psi) and the ray contour integrals a1, a2 of
    !! gs_spline_level_integrals. The boundary level u_b solves
    !! Phi(u_b) = Phi_edge, so it is a result here, not a cut fraction.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_space, only: spline_space_t
    use milca_gs_profile_law, only: gs_profile_law_t
    use milca_gs_spline, only: gs_spline_level_integrals, GS_SPLINE_OK
    use milca_gs_free_current, only: gs_current_target_t, hermite_law, interpolate4, &
        cubic_integral, flux_function
    use milca_grad_shafranov_leaves, only: gs_current_profile_update, gs_linear_extrapolation
    implicit none
    private
    public :: gs_spline_current_update
    integer, parameter, public :: GS_SPLINE_CUR_OK = 0, GS_SPLINE_CUR_FAILED = 1, &
        GS_SPLINE_CUR_FLUX_UNREACHABLE = 2
    integer, parameter :: table = 48

contains

    subroutine gs_spline_current_update(space, coef, origin, psi_axis, psi_bound, sigma, &
                                        mu0, law, target, n_rays, next, psi_bound_next, &
                                        status)
        !! The next law (on the knots of law) and boundary flux from the flux
        !! coef with axis value psi_axis and boundary psi_bound, solved under
        !! law. The table of Phi(u) reaches 0.5 % beyond the current boundary
        !! (2 or 5 % if Phi_edge lies further out), inside the vacuum annulus.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :), origin(2), psi_axis, psi_bound, sigma, mu0
        type(gs_profile_law_t), intent(in) :: law
        type(gs_current_target_t), intent(in) :: target
        integer, intent(in) :: n_rays
        type(gs_profile_law_t), intent(out) :: next
        real(dp), intent(out) :: psi_bound_next
        integer, intent(out) :: status
        type(gs_profile_law_t) :: current
        real(dp) :: t(table), u(table), a1(table), a2(table), f(table), phi(table)
        real(dp), allocatable :: uk(:), k1(:), k2(:), fk(:), sk(:), pp(:), ff(:)
        real(dp) :: ub, ub_old, top, s_tor, dpds, dids
        real(dp), parameter :: tops(3) = [1.005_dp, 1.02_dp, 1.05_dp]
        integer :: j, nk, st, attempt

        status = GS_SPLINE_CUR_FAILED
        psi_bound_next = psi_bound
        nk = size(law%knots)
        allocate (uk(nk), k1(nk), k2(nk), fk(nk), sk(nk), pp(nk), ff(nk))
        current = law
        current%axis = psi_axis - psi_bound
        ub_old = sigma*(psi_axis - psi_bound)
        ! The table reaches just beyond the current boundary: Gamma_p may sit
        ! within a few per cent of a separatrix, beyond which the levels open.
        do attempt = 1, 3
            top = tops(attempt)
            do j = 1, table
                t(j) = top*(1 - (1 - real(j, dp)/table)**2)
            end do
            u = t*ub_old
            call gs_spline_level_integrals(space, coef, origin, psi_axis - sigma*u, sigma, &
                                           n_rays, a1, a2, st)
            if (st /= GS_SPLINE_OK) return
            do j = 1, table
                f(j) = flux_function(current, u(j), ub_old)
            end do
            phi = 0
            do j = 1, table
                phi(j) = merge(0.0_dp, phi(max(j - 1, 1)), j == 1) + &
                    cubic_integral(u, f*a2, merge(0.0_dp, u(max(j - 1, 1)), j == 1), u(j))
            end do
            if (phi(table) >= target%phi_edge) exit
        end do
        if (phi(table) < target%phi_edge) then
            status = GS_SPLINE_CUR_FLUX_UNREACHABLE
            return
        end if
        ub = interpolate4(phi, u, target%phi_edge, count(phi < target%phi_edge))
        uk = law%knots*ub
        uk(1) = 1.0e-2_dp*uk(2)
        call gs_spline_level_integrals(space, coef, origin, psi_axis - sigma*uk, sigma, &
                                       n_rays, k1, k2, st)
        if (st /= GS_SPLINE_OK) return
        do j = 1, nk
            fk(j) = flux_function(current, uk(j), ub_old)
            sk(j) = interpolate4(u, phi, uk(j), count(u < uk(j)))/target%phi_edge
        end do
        sk(nk) = 1
        do j = 1, nk
            s_tor = min(max(sk(j), 0.0_dp), 1.0_dp)
            dpds = interpolate4(target%s, target%dp, s_tor, count(target%s < s_tor))
            dids = interpolate4(target%s, target%dcurrent, s_tor, count(target%s < s_tor))
            call gs_current_profile_update(sigma, mu0, fk(j), k1(j), k2(j), target%phi_edge, &
                                           dpds, dids, pp(j), ff(j))
        end do
        call gs_linear_extrapolation(0.0_dp, uk(1), uk(2), pp(1), pp(2), pp(1))
        call gs_linear_extrapolation(0.0_dp, uk(1), uk(2), ff(1), ff(2), ff(1))
        next = law
        call hermite_law(law%knots, pp, next%pp)
        call hermite_law(law%knots, ff, next%ff)
        next%p_edge = target%p(size(target%p))
        psi_bound_next = psi_axis - sigma*ub
        status = GS_SPLINE_CUR_OK
    end subroutine gs_spline_current_update
end module milca_gs_spline_current
