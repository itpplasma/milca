module ideal_flux_charts
    ! Analytic test equilibria for milca_ideal_flux.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_flux, only: flux_equilibrium_t, flux_point_t
    implicit none
    private
    public :: cylinder_t, relabeled_t, rippled_t
    !> Straight periodic cylinder x = (a rho cos t, a rho sin t, R0 zeta).
    type, extends(flux_equilibrium_t) :: cylinder_t
        real(dp) :: a = 1, r0 = 1
        integer :: profile = 1
    contains
        procedure :: evaluate => cylinder
    end type cylinder_t
    !> The same cylinder in the angle theta with phi = theta + k rho^2 sin theta
    !> the geometric angle: lambda = k rho^2 sin theta.
    type, extends(cylinder_t) :: relabeled_t
        real(dp) :: k = 0.5_dp
    contains
        procedure :: evaluate => relabeled
    end type relabeled_t
    !> A chart with two field periods, R = 3 + rho cos t + 0.15 rho cos(t - 2z),
    !> Z = rho sin t - 0.15 rho sin(t - 2z), phi = z; not an equilibrium.
    type, extends(flux_equilibrium_t) :: rippled_t
    contains
        procedure :: evaluate => rippled
    end type rippled_t
contains
    subroutine cylinder(self, rho, theta, zeta, point)
        class(cylinder_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp) :: r, bz, bzp, p, pp, c, s, cyl_a, cyl_r0
        cyl_a = self%a; cyl_r0 = self%r0
        r = cyl_a*rho; c = cos(theta + 0*zeta); s = sin(theta)
        if (self%profile == 1) then
            bz = sqrt(0.25_dp + r*r); bzp = r/bz
            p = 1.5_dp*(cyl_a**2 - r*r); pp = -3*r
        else
            bz = 1; bzp = 0; p = 0.3_dp; pp = 0
        end if
        point%xu(:, 1) = [cyl_a*c, cyl_a*s, 0.0_dp]
        point%xu(:, 2) = [-r*s, r*c, 0.0_dp]
        point%xu(:, 3) = [0.0_dp, 0.0_dp, cyl_r0]
        point%xuu = 0
        point%xuu(:, 1, 2) = [-cyl_a*s, cyl_a*c, 0.0_dp]
        point%xuu(:, 2, 1) = point%xuu(:, 1, 2)
        point%xuu(:, 2, 2) = [-r*c, -r*s, 0.0_dp]
        ! psi' = a r Bz, iota = R0 Btheta/(r Bz) in rho; zeta = z/R0.
        point%psi1 = cyl_a*r*bz
        point%psi2 = cyl_a**2*(bz + r*bzp)
        if (self%profile == 1) then
            point%iota = cyl_r0/bz
            point%iota1 = -cyl_a*cyl_r0*bzp/bz**2
        else
            point%iota = 0; point%iota1 = 0
        end if
        point%lam = 0
        point%p = p; point%pp = cyl_a*pp; point%density = 1
    end subroutine cylinder

    subroutine relabeled(self, rho, theta, zeta, point)
        class(relabeled_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp) :: phi, fr, ft, frr, frt, ftt, xr(3), xf(3), xrf(3), xff(3)
        phi = theta + self%k*rho**2*sin(theta)
        fr = 2*self%k*rho*sin(theta); ft = 1 + self%k*rho**2*cos(theta)
        frr = 2*self%k*sin(theta); frt = 2*self%k*rho*cos(theta)
        ftt = -self%k*rho**2*sin(theta)
        call cylinder(self, rho, phi, zeta, point)
        xr = point%xu(:, 1); xf = point%xu(:, 2)
        xrf = point%xuu(:, 1, 2); xff = point%xuu(:, 2, 2)
        point%xu(:, 1) = xr + xf*fr
        point%xu(:, 2) = xf*ft
        point%xuu(:, 1, 1) = 2*xrf*fr + xff*fr**2 + xf*frr
        point%xuu(:, 1, 2) = xrf*ft + xff*fr*ft + xf*frt
        point%xuu(:, 2, 1) = point%xuu(:, 1, 2)
        point%xuu(:, 2, 2) = xff*ft**2 + xf*ftt
        point%lam = 0
        point%lam(1) = self%k*rho**2*cos(theta)
        point%lam(3) = 2*self%k*rho*cos(theta)
        point%lam(4) = -self%k*rho**2*sin(theta)
    end subroutine relabeled

    subroutine rippled(self, rho, theta, zeta, point)
        class(rippled_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp), parameter :: r0 = 3, d = 0.15_dp
        real(dp) :: c, s, c2, s2, big, rr(3), zz(3), rz(3, 3), zzz(3, 3), e(3), &
            ep(3)
        integer :: i, j
        c = cos(theta); s = sin(theta)
        c2 = cos(theta - 2*zeta); s2 = sin(theta - 2*zeta)
        rr = [c + d*c2, -rho*s - d*rho*s2, 2*d*rho*s2]
        zz = [s - d*s2, rho*c - d*rho*c2, 2*d*rho*c2]
        rz = 0; zzz = 0
        rz(1, 2) = -s - d*s2; rz(1, 3) = 2*d*s2
        rz(2, 2) = -rho*c - d*rho*c2; rz(2, 3) = 2*d*rho*c2; rz(3, 3) = -4*d*rho*c2
        zzz(1, 2) = c - d*c2; zzz(1, 3) = 2*d*c2
        zzz(2, 2) = -rho*s + d*rho*s2; zzz(2, 3) = -2*d*rho*s2; zzz(3, 3) = 4*d*rho*s2
        do j = 1, 3
            do i = 1, j - 1
                rz(j, i) = rz(i, j); zzz(j, i) = zzz(i, j)
            end do
        end do
        big = r0 + rho*c + d*rho*c2
        e = [cos(zeta), sin(zeta), 0.0_dp]; ep = [-sin(zeta), cos(zeta), 0.0_dp]
        do i = 1, 3
            point%xu(:, i) = rr(i)*e + [0.0_dp, 0.0_dp, zz(i)]
        end do
        point%xu(:, 3) = point%xu(:, 3) + big*ep
        do j = 1, 3
            do i = 1, 3
                point%xuu(:, i, j) = rz(i, j)*e + [0.0_dp, 0.0_dp, zzz(i, j)]
                if (i == 3) point%xuu(:, i, j) = point%xuu(:, i, j) + rr(j)*ep
                if (j == 3) point%xuu(:, i, j) = point%xuu(:, i, j) + rr(i)*ep
                if (i == 3 .and. j == 3) point%xuu(:, i, j) = point%xuu(:, i, j) &
                    - big*e
            end do
        end do
        point%psi1 = 1.3_dp*rho; point%psi2 = 1.3_dp
        point%iota = 0.4_dp + 0.3_dp*rho**2; point%iota1 = 0.6_dp*rho
        point%lam = 0
        point%p = 0.1_dp*(1 - rho**2); point%pp = -0.2_dp*rho; point%density = 1
    end subroutine rippled
end module ideal_flux_charts
