! Photon emission of an ultrarelativistic electron in a strong constant field
! in the locally constant field approximation (Baier and Katkov 1968; Ritus
! 1985), spin-averaged. With the quantum parameter chi = gamma |F_perp|/F_crit,
! F_crit = m^2 c^3/(e hbar), and the photon energy fraction delta, the
! emission rate is
!   dW/d delta = alpha m c^2/(hbar gamma) S(chi, delta),
!   S = [(1 - delta + 1/(1 - delta)) K_{2/3}(y) - int_y^inf K_{1/3}]/(sqrt(3) pi),
!   y = 2 delta/(3 chi (1 - delta)).
! The emitted power over the classical (Landau-Lifshitz) power
! (2/3) alpha chi^2 m^2 c^4/hbar is the Gaunt factor
!   g(chi) = 3/(2 chi^2) int_0^1 delta S d delta,
! derived limits g -> 1 - (55 sqrt(3)/16) chi as chi -> 0, and the
! classical synchrotron spectrum int_y^inf K_{5/3} at delta << 1.
! Pair creation by a photon of energy fraction delta to the electron, with the
! photon's chi and gamma = photon energy/(m c^2), the crossed rate
!   dW/d delta = alpha m c^2/(hbar gamma) P(chi, delta),
!   P = [(delta/(1 - delta) + (1 - delta)/delta) K_{2/3}(y) + int_y^inf K_{1/3}]/(sqrt(3) pi),
!   y = 2/(3 chi delta (1 - delta)),
! whose total at chi << 1 is, by the saddle at delta = 1/2,
! (3 sqrt(3)/(16 sqrt(2))) chi exp(-8/(3 chi)) (Erber 1966).
module milca_strong_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: bessel_k_frac, bessel_k_tail, quantum_synchrotron_rate, quantum_gaunt
    public :: pair_creation_spectrum, pair_creation_rate

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: nk = 160

contains

    !> exp(x) K_nu(x) for real nu, from int_0^inf exp(-x (cosh t - 1)) cosh(nu t) dt.
    real(dp) function bessel_k_frac(nu, x) result(k)
        real(dp), intent(in) :: nu, x
        k = k_integral(nu, x, .false.)
    end function bessel_k_frac

    !> exp(x) int_x^inf K_nu(s) ds = exp(x) int_0^inf exp(-x cosh t) cosh(nu t)/cosh t dt.
    real(dp) function bessel_k_tail(nu, x) result(k)
        real(dp), intent(in) :: nu, x
        k = k_integral(nu, x, .true.)
    end function bessel_k_tail

    real(dp) function k_integral(nu, x, tail) result(k)
        real(dp), intent(in) :: nu, x
        logical, intent(in) :: tail
        real(dp) :: t(nk), w(nk), tmax, tt, f
        integer :: i
        call gauss_legendre(nk, t, w)
        tmax = acosh(1 + 60/x)
        k = 0
        do i = 1, nk
            ! t = tmax u^2 resolves the region t << 1, where small x puts its weight
            tt = tmax*((t(i) + 1)/2)**2
            f = exp(-x*(cosh(tt) - 1))*cosh(nu*tt)
            if (tail) f = f/cosh(tt)
            k = k + w(i)*tmax*(t(i) + 1)/2*f
        end do
    end function k_integral

    !> S(chi, delta) of the emission rate.
    real(dp) function quantum_synchrotron_rate(chi, delta) result(s)
        real(dp), intent(in) :: chi, delta
        real(dp) :: y
        y = 2*delta/(3*chi*(1 - delta))
        s = ((1 - delta + 1/(1 - delta))*bessel_k_frac(2.0_dp/3, y) - &
            bessel_k_tail(1.0_dp/3, y))*exp(-y)/(sqrt(3.0_dp)*pi)
    end function quantum_synchrotron_rate

    !> Emitted power over the classical power, with n points in u, where
    !> y = 60 u^3 and delta = 3 chi y/(2 + 3 chi y).
    real(dp) function quantum_gaunt(chi, n) result(g)
        real(dp), intent(in) :: chi
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), u, y, delta, jac
        integer :: i
        call gauss_legendre(n, x, w)
        g = 0
        do i = 1, n
            u = (x(i) + 1)/2
            y = 60*u**3
            delta = 3*chi*y/(2 + 3*chi*y)
            jac = 0.5_dp*180*u**2*6*chi/(2 + 3*chi*y)**2
            g = g + w(i)*jac*delta*quantum_synchrotron_rate(chi, delta)
        end do
        g = 1.5_dp*g/chi**2
    end function quantum_gaunt

    !> P(chi, delta) of the pair-creation rate.
    real(dp) function pair_creation_spectrum(chi, delta) result(s)
        real(dp), intent(in) :: chi, delta
        real(dp) :: y
        y = 2/(3*chi*delta*(1 - delta))
        s = ((delta/(1 - delta) + (1 - delta)/delta)*bessel_k_frac(2.0_dp/3, y) + &
            bessel_k_tail(1.0_dp/3, y))*exp(-y)/(sqrt(3.0_dp)*pi)
    end function pair_creation_spectrum

    !> int_0^1 P d delta times exp(8/(3 chi)), with n points on each half,
    !> delta = (1 - cos(pi u/2))/2 clustering at the ends.
    real(dp) function pair_creation_rate(chi, n) result(r)
        real(dp), intent(in) :: chi
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), u, delta, y
        integer :: i
        call gauss_legendre(n, x, w)
        r = 0
        do i = 1, n
            u = (x(i) + 1)/2
            delta = (1 - cos(pi*u/2))/2
            y = 2/(3*chi*delta*(1 - delta))
            ! twice the half 0 < delta < 1/2; exp(-y + 8/(3 chi)) stays finite
            r = r + w(i)*pi/8*sin(pi*u/2)*2* &
                ((delta/(1 - delta) + (1 - delta)/delta)*bessel_k_frac(2.0_dp/3, y) + &
                bessel_k_tail(1.0_dp/3, y))*exp(-y + 8/(3*chi))/(sqrt(3.0_dp)*pi)
        end do
    end function pair_creation_rate

end module milca_strong_field
