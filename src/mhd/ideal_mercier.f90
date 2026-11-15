module milca_ideal_mercier
    ! Local interchange criteria on one surface of a flux-coordinate
    ! equilibrium (docs/stability.md, M3): the finite-radius Mercier criterion
    ! D_Merc = D_shear + D_curr + D_well + D_geod > 0 of Bauer et al. (1984)
    ! as written by Landreman and Jorge, J. Plasma Phys. 86, 905860510 (2020),
    ! eqs. (4.16)-(4.20), its equivalent form (4.1), and the magnetic well
    ! d2V/dpsi2 (stable for d2V/dpsi2 < 0). psi is the toroidal flux per
    ! radian of milca_ideal_flux, mu0 J = curl B and mu0 p enters as mu p.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_flux, only: flux_equilibrium_t, flux_point_t
    use milca_ideal_flux_leaves, only: flux_jets_array, flux_metric_array
    implicit none
    private
    public :: mercier_t, flux_mercier

    !> Surface quantities, derivatives with respect to psi.
    type :: mercier_t
        real(dp) :: iota = 0, shear = 0, dpdpsi = 0, vp = 0, vpp = 0
        real(dp) :: d_shear = 0, d_curr = 0, d_well = 0, d_geod = 0
        !> Sum of the four terms and the value of form (4.1).
        real(dp) :: total = 0, total_alt = 0
    end type mercier_t
contains
    subroutine flux_mercier(equilibrium, rho, nfp, ntheta, nzeta, mu, r, status)
        ! Trapezoidal surface integrals over one field period, times nfp.
        class(flux_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: rho, mu
        integer, intent(in) :: nfp, ntheta, nzeta
        type(mercier_t), intent(out) :: r
        integer, intent(out) :: status
        real(dp), parameter :: pi = acos(-1.0_dp), twopi = 2*pi
        type(flux_point_t) :: e
        real(dp) :: jac, b(3), g(3, 3), q(3), t(12, 12), jacr, grho(3), gzeta(3)
        real(dp) :: mj(3), b2, jb, gpsi2, w, w0, theta, zeta, psi1, psi2, iota1
        real(dp) :: a0, a1, s_b2, s_1b2, s_jb, s_jb2, s_j2, itor, gcov, ip, sg, &
            spsi, mp, xib, xi2
        integer :: i, j
        status = 1
        if (rho <= 0 .or. rho >= 1 .or. nfp < 1 .or. ntheta < 1 .or. nzeta < 1) return
        a0 = 0; a1 = 0; s_b2 = 0; s_1b2 = 0; s_jb = 0; s_jb2 = 0; s_j2 = 0
        itor = 0; gcov = 0
        w0 = nfp*(twopi/ntheta)*(twopi/(nfp*nzeta))
        do j = 1, nzeta
            do i = 1, ntheta
                theta = twopi*(i - 1)/ntheta
                zeta = twopi*(j - 1)/(nfp*nzeta)
                call equilibrium%evaluate(rho, theta, zeta, e)
                call flux_jets_array(e%xu, e%xuu, e%psi1, e%psi2, e%iota, &
                    e%iota1, e%pp, e%lam, jac, b, g, q, t)
                call flux_metric_array(e%xu, e%xuu, jacr, grho, gzeta)
                mj = [g(3, 2) - g(2, 3), g(1, 3) - g(3, 1), g(2, 1) - g(1, 2)]
                b2 = dot_product(b, b); jb = dot_product(mj, b)
                gpsi2 = e%psi1**2*dot_product(grho, grho)
                ! dS/|grad psi| = |sqrt(g)| dtheta dzeta / |psi'|.
                w = w0*abs(jac)/abs(e%psi1)
                a0 = a0 + w0*abs(jac)
                a1 = a1 + w0*sign(1.0_dp, jac)*jacr
                s_b2 = s_b2 + w*b2/gpsi2
                s_1b2 = s_1b2 + w/b2
                s_jb = s_jb + w*jb/gpsi2
                s_jb2 = s_jb2 + w*jb**2/(b2*gpsi2)
                s_j2 = s_j2 + w*dot_product(mj, mj)/gpsi2
                ! mu0 dI_tor/drho through x_rho x x_theta = sqrt(g) grad zeta.
                itor = itor + w0*dot_product(mj, jac*gzeta)
                gcov = gcov + w0*dot_product(b, e%xu(:, 3))
            end do
        end do
        psi1 = e%psi1; psi2 = e%psi2; iota1 = e%iota1
        r%iota = e%iota
        r%shear = iota1/psi1
        r%dpdpsi = e%pp/psi1
        r%vp = a0/psi1
        r%vpp = (a1/psi1 - a0*psi2/psi1**2)/psi1
        ! I(psi) = mu0 I_tor/(2 pi), the angle mean of the covariant B_theta.
        ip = itor/(twopi*twopi)/psi1
        sg = sign(1.0_dp, gcov); spsi = sign(1.0_dp, psi1)
        mp = mu*r%dpdpsi
        r%d_shear = r%shear**2/(16*pi**2)
        r%d_curr = -sg*r%shear/twopi**4*(s_jb - ip*s_b2)
        r%d_well = mp/twopi**6*(spsi*r%vpp - mp*s_1b2)*s_b2
        r%d_geod = (s_jb**2 - s_b2*s_jb2)/twopi**6
        r%total = r%d_shear + r%d_curr + r%d_well + r%d_geod
        ! Form (4.1) with Xi = mu0 J - I'(psi) B.
        xib = s_jb - ip*s_b2
        xi2 = s_j2 - 2*ip*s_jb + ip**2*s_b2
        r%total_alt = ((sg*2*pi**2*r%shear - xib)**2 &
            + (spsi*mp*r%vpp - xi2)*s_b2)/twopi**6
        if (.not. all(ieee_is_finite([r%total, r%total_alt, r%vpp]))) return
        status = 0
    end subroutine flux_mercier
end module milca_ideal_mercier
