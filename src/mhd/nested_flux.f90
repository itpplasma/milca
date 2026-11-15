module milca_nested_flux
    ! The sheared exact equilibria of milca_nested_oracles (Landreman, Issan
    ! et al.) as flux-coordinate equilibria for milca_ideal_flux and the
    ! local criteria (docs/stability.md). Chart (rho, theta, zeta) is the
    ! oracle's surface map with rho scaled to 1 on the boundary psi = delta and
    ! zeta its toroidal parameter t; B, grad B and p are the exact fields. The
    ! toroidal flux per radian psi(rho), iota and the stream-function jets
    ! follow from sqrt(g) B^zeta and sqrt(g) B^theta (flux_stream): psi' is
    ! the theta-mean of sqrt(g) B^zeta at zeta = 0 and psi' iota the
    ! zeta-mean of sqrt(g) B^theta at theta = 0, by periodic trapezoidal
    ! quadrature with nquad points. mu0 = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_flux, only: flux_equilibrium_t, flux_point_t
    use milca_ideal_flux_leaves, only: flux_stream_array
    use milca_nested_oracles, only: nested_oracle_t, NESTED_SHEARED
    use milca_nested_oracle_leaves, only: nested_sheared_chart
    implicit none
    private
    public :: nested_flux_t

    type, extends(flux_equilibrium_t) :: nested_flux_t
        type(nested_oracle_t) :: oracle
        integer :: nquad = 64
    contains
        procedure :: evaluate => nested_evaluate
        procedure :: valid => nested_valid
    end type nested_flux_t
contains
    logical function nested_valid(self)
        class(nested_flux_t), intent(in) :: self
        nested_valid = self%oracle%family == NESTED_SHEARED .and. &
            self%oracle%valid() .and. self%nquad >= 8
    end function nested_valid

    subroutine nested_evaluate(self, rho, theta, zeta, point)
        class(nested_flux_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp), parameter :: twopi = 2*acos(-1.0_dp)
        real(dp) :: b(3), g(3, 3), gp(3), p, sb(4), lam(7), sums(4), psi1, psi2
        real(dp) :: iota, iota1, angle, bk(3), gk(3, 3), gpk(3), pk, xu(3, 3), &
            xuu(3, 3, 3)
        integer :: k
        sums = 0
        do k = 1, self%nquad
            angle = twopi*(k - 1)/self%nquad
            call jets(rho, angle, 0.0_dp, xu, xuu, bk, gk, gpk, pk)
            call flux_stream_array(xu, xuu, bk, gk, 1.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, lam, sb)
            sums([2, 4]) = sums([2, 4]) + sb([2, 4])
            call jets(rho, 0.0_dp, angle, xu, xuu, bk, gk, gpk, pk)
            call flux_stream_array(xu, xuu, bk, gk, 1.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, lam, sb)
            sums([1, 3]) = sums([1, 3]) + sb([1, 3])
        end do
        sums = sums/self%nquad
        psi1 = sums(2); psi2 = sums(4)
        iota = sums(1)/psi1
        iota1 = (sums(3) - iota*psi2)/psi1
        call jets(rho, theta, zeta, point%xu, point%xuu, b, g, gp, p)
        call flux_stream_array(point%xu, point%xuu, b, g, psi1, psi2, iota, &
            iota1, point%lam, sb)
        point%psi1 = psi1; point%psi2 = psi2
        point%iota = iota; point%iota1 = iota1
        point%p = p; point%pp = dot_product(gp, point%xu(:, 1))
        point%density = 1
    contains
        subroutine jets(r, th, t, xu_, xuu_, bv, gv, gpv, pv)
            ! Chart jets in the scaled radius and the exact field there.
            real(dp), intent(in) :: r, th, t
            real(dp), intent(out) :: xu_(3, 3), xuu_(3, 3, 3), bv(3), gv(3, 3), &
                gpv(3), pv
            real(dp) :: scale, v(20)
            integer :: ok
            scale = sqrt(2*self%oracle%delta)
            associate(o => self%oracle)
                call nested_sheared_chart(scale*r, th, t, o%eps, o%tau, o%saxis, &
                    o%lam, xu_(1, 1), xu_(2, 1), xu_(3, 1), xu_(1, 2), xu_(2, 2), &
                    xu_(3, 2), xu_(1, 3), xu_(2, 3), xu_(3, 3), xuu_(1, 1, 1), &
                    xuu_(2, 1, 1), xuu_(3, 1, 1), xuu_(1, 1, 2), xuu_(2, 1, 2), &
                    xuu_(3, 1, 2), xuu_(1, 2, 2), xuu_(2, 2, 2), xuu_(3, 2, 2), &
                    xuu_(1, 1, 3), xuu_(2, 1, 3), xuu_(3, 1, 3), xuu_(1, 2, 3), &
                    xuu_(2, 2, 3), xuu_(3, 2, 3), xuu_(1, 3, 3), xuu_(2, 3, 3), &
                    xuu_(3, 3, 3))
            end associate
            xuu_(:, 2, 1) = xuu_(:, 1, 2); xuu_(:, 3, 1) = xuu_(:, 1, 3)
            xuu_(:, 3, 2) = xuu_(:, 2, 3)
            xu_(:, 1) = scale*xu_(:, 1)
            xuu_(:, 1, :) = scale*xuu_(:, 1, :)
            xuu_(:, :, 1) = scale*xuu_(:, :, 1)
            call self%oracle%point(position(scale*r, th, t), v, ok)
            pv = v(2); bv = v(3:5); gpv = v(9:11); gv = reshape(v(12:20), [3, 3])
        end subroutine jets

        function position(r, th, t) result(x)
            real(dp), intent(in) :: r, th, t
            real(dp) :: x(3), jac
            call self%oracle%surface(r, th, t, x, jac)
        end function position
    end subroutine nested_evaluate
end module milca_nested_flux
