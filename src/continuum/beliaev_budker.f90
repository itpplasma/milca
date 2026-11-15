! Relativistic test-particle collision coefficients (bb_test of
! docs/collisions.md): with the Beliaev-Budker kernel U of
! milca_beliaev_budker_leaves and a background f',
!   C(f) = Gamma d/du_i (D^ij d f/du_j - A^i f),
!   D^ij = (1/2) int U^ij f' d^3u',  A^i = (1/2) int U^ij d f'/du'_j d^3u',
! per Gamma n for a Juttner background of temperature theta = kT/(m c^2),
! u = p/(m c). Where the Juttner weight at u is not negligible, the integral
! is centred on u, so that the 1/|u - u'| singularity of the kernel is
! integrated in spherical coordinates around it; otherwise it is centred on
! the background, which is then smooth and far from the singular point.
module milca_beliaev_budker
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_beliaev_budker_leaves, only: bb_kernel
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: juttner_bb_coefficients

    real(dp), parameter :: pi = acos(-1.0_dp)
    !> Juttner tail cut: exp(-(gamma - 1)/theta) below exp(-TAIL).
    real(dp), parameter :: TAIL = 60

contains

    !> D(3,3) and A(3) per Gamma n at test momentum u for a Juttner
    !> background of temperature theta, with n quadrature points per
    !> dimension.
    subroutine juttner_bb_coefficients(u, theta, n, d, a)
        real(dp), intent(in) :: u(3), theta
        integer, intent(in) :: n
        real(dp), intent(out) :: d(3, 3), a(3)
        real(dp) :: norm, gu
        gu = sqrt(1 + dot_product(u, u))
        norm = juttner_norm(theta, n)
        if ((gu - 1)/theta < TAIL) then
            call centred_on_test(u, theta, n, d, a)
        else
            call centred_on_background(u, theta, n, d, a)
        end if
        d = d/norm
        a = a/norm
    end subroutine juttner_bb_coefficients

    ! int exp(-(gamma - 1)/theta) d^3u, kinetic energy theta t^2
    real(dp) function juttner_norm(theta, n) result(norm)
        real(dp), intent(in) :: theta
        integer, intent(in) :: n
        real(dp) :: x(4*n), w(4*n), t, g, p
        integer :: i
        call gauss_legendre(4*n, x, w)
        norm = 0
        do i = 1, 4*n
            t = sqrt(TAIL)*(x(i) + 1)/2
            g = 1 + theta*t**2
            p = sqrt(g**2 - 1)
            ! d^3u = 4 pi p^2 dp, p dp = g dg, dg = 2 theta t dt
            norm = norm + w(i)*sqrt(TAIL)/2*4*pi*p*g*2*theta*t*exp(-t**2)
        end do
    end function juttner_norm

    ! the integrand at background momentum v: U f' and U grad' f'
    subroutine accumulate(u, v, theta, weight, d, a)
        real(dp), intent(in) :: u(3), v(3), theta, weight
        real(dp), intent(inout) :: d(3, 3), a(3)
        real(dp) :: k(6), m(3, 3), gv, f, grad(3)
        gv = sqrt(1 + dot_product(v, v))
        f = exp(-(gv - 1)/theta)
        if (f == 0) return
        call bb_kernel(u(1), u(2), u(3), v(1), v(2), v(3), k(1), k(2), k(3), k(4), k(5), k(6))
        m = reshape([k(1), k(2), k(4), k(2), k(3), k(5), k(4), k(5), k(6)], [3, 3])
        grad = -v/(gv*theta)*f
        d = d + weight*m*f/2
        a = a + weight*matmul(m, grad)/2
    end subroutine accumulate

    subroutine centred_on_test(u, theta, n, d, a)
        real(dp), intent(in) :: u(3), theta
        integer, intent(in) :: n
        real(dp), intent(out) :: d(3, 3), a(3)
        real(dp) :: x(2*n), w(2*n), xa(n), wa(n), rho, rmax, c, s, phi, dir(3)
        integer :: i, j, l
        call gauss_legendre(2*n, x, w)
        call gauss_legendre(n, xa, wa)
        rmax = norm2(u) + sqrt((1 + TAIL*theta)**2 - 1)
        d = 0
        a = 0
        do i = 1, 2*n
            rho = rmax*(x(i) + 1)/2
            do j = 1, n
                c = xa(j)
                s = sqrt(1 - c**2)
                do l = 1, n
                    phi = pi*(xa(l) + 1)
                    dir = [s*cos(phi), s*sin(phi), c]
                    call accumulate(u, u + rho*dir, theta, &
                        w(i)*rmax/2*rho**2*wa(j)*wa(l)*pi, d, a)
                end do
            end do
        end do
    end subroutine centred_on_test

    subroutine centred_on_background(u, theta, n, d, a)
        real(dp), intent(in) :: u(3), theta
        integer, intent(in) :: n
        real(dp), intent(out) :: d(3, 3), a(3)
        real(dp) :: x(2*n), w(2*n), xa(n), wa(n), t, g, p, c, s, phi, jac
        integer :: i, j, l
        call gauss_legendre(2*n, x, w)
        call gauss_legendre(n, xa, wa)
        d = 0
        a = 0
        do i = 1, 2*n
            t = sqrt(TAIL)*(x(i) + 1)/2
            g = 1 + theta*t**2
            p = sqrt(g**2 - 1)
            ! p^2 dp = p g dg = p g 2 theta t dt
            jac = sqrt(TAIL)/2*p*g*2*theta*t
            do j = 1, n
                c = xa(j)
                s = sqrt(1 - c**2)
                do l = 1, n
                    phi = pi*(xa(l) + 1)
                    call accumulate(u, p*[s*cos(phi), s*sin(phi), c], theta, &
                        w(i)*jac*wa(j)*wa(l)*pi, d, a)
                end do
            end do
        end do
    end subroutine centred_on_background

end module milca_beliaev_budker
