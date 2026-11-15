! Compton scattering of a photon of energy eps (in m c^2) on a Juttner
! electron gas at theta = kT/(m c^2), with the Klein-Nishina cross section in
! the electron rest frame,
!   d sigma/d Omega_r = (r_e^2/2) (e'/e)^2 (e'/e + e/e' - sin^2 Theta),
!   e' = e/(1 + e (1 - cos Theta)),  e = eps gamma (1 - beta cos alpha),
! and the scattered energy back in the lab, gamma e' (1 + beta cos alpha'),
! alpha' the rest-frame angle between scattered photon and electron velocity.
! The scattering rate per photon and the rate of its energy change, per n_e
! r_e^2 c, average (1 - beta cos alpha) over the electrons. Derived limits:
! the mean fractional gain per scattering is (4/3) <p^2> = 4 theta K3/K2 as
! eps -> 0 (Thomson), -eps for cold electrons at eps << 1 (recoil), and 4
! theta - eps when both are small (Kompaneets).
module milca_thermal_compton
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    private

    public :: compton_rates

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    !> Scattering rate and energy-change rate of a photon of energy eps,
    !> n points per integration variable.
    subroutine compton_rates(eps, theta, n, rate, heating)
        real(dp), intent(in) :: eps, theta
        integer, intent(in) :: n
        real(dp), intent(out) :: rate, heating
        real(dp) :: x(n), w(n), umax, t, u, gam, p, beta, fw, norm, ca, flux, er, car, sar
        real(dp) :: ct, st, ep, cphi, cap, de, dsig, wout
        integer :: i, ia, it, ip
        call gauss_legendre(n, x, w)
        umax = 60
        norm = 0
        rate = 0
        heating = 0
        do i = 1, n
            ! kinetic energy theta u, u = t^2 smooth at rest
            t = sqrt(umax)*(x(i) + 1)/2
            u = t**2
            gam = 1 + theta*u
            p = sqrt(gam**2 - 1)
            beta = p/gam
            ! p^2 dp exp(-(gam - 1)/theta) = p gam theta exp(-u) du
            fw = w(i)*sqrt(umax)/2*2*t*p*gam*theta*exp(-u)
            norm = norm + fw
            do ia = 1, n
                ca = x(ia)
                flux = 1 - beta*ca
                er = eps*gam*flux
                ! the photon direction in the rest frame, to the electron velocity
                car = (ca - beta)/flux
                sar = sqrt(max(0.0_dp, 1 - car**2))
                do it = 1, n
                    ct = x(it)
                    st = sqrt(max(0.0_dp, 1 - ct**2))
                    ep = er/(1 + er*(1 - ct))
                    dsig = 0.5_dp*(ep/er)**2*(ep/er + er/ep - st**2)
                    wout = 0
                    do ip = 1, n
                        cphi = cos(pi*(x(ip) + 1))
                        cap = car*ct + sar*st*cphi
                        de = gam*ep*(1 + beta*cap) - eps
                        wout = wout + w(ip)*pi*de
                    end do
                    rate = rate + fw*w(ia)/2*w(it)*flux*2*pi*dsig
                    heating = heating + fw*w(ia)/2*w(it)*flux*dsig*wout
                end do
            end do
        end do
        rate = rate/norm
        heating = heating/norm
    end subroutine compton_rates

end module milca_thermal_compton
