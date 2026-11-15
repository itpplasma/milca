! Escape of cyclotron and synchrotron radiation from a uniform slab of
! thickness L with the magnetic field in its plane and specularly reflecting
! walls of reflectivity R, which keep the angle to the field. Absorption
! follows from the thermal emissivity by Kirchhoff's law against the
! Rayleigh-Jeans intensity omega^2 kT/(8 pi^3 c^2) (both polarizations,
! vacuum propagation), so along a ray at angle psi to the wall normal
!   tau = pi Lambda j(w, mu)/(theta w^2 |cos psi|),  Lambda = omega_pe^2 L/(omega_c c),
! with j of milca_cyclotron_emission. A ray bouncing between the walls leaves
! the intensity B (1 - e^-tau)/(1 - R e^-tau) at the wall, and the power
! escaping per volume, over the optically thin emission, is weighted by
!   E(tau) = (1 - R)(1 - e^-tau)/(tau (1 - R e^-tau)):
!   P/(n_e K omega_c^2) = int dw int dOmega j(w, mu) E(tau),
! cos psi = sin theta cos phi about the field. Derived limits: the thin
! emission as tau -> 0 for any R, and the blackbody flux of both walls,
! (theta/(pi Lambda)) 2 pi int w^2 dw, where tau >> 1.
module milca_synchrotron_transport
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_cyclotron_emission, only: thermal_emissivity
    implicit none
    private

    public :: slab_synchrotron_power, escape_factor

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    pure real(dp) function escape_factor(tau, refl) result(e)
        real(dp), intent(in) :: tau, refl
        real(dp) :: em, frac
        em = exp(-tau)
        ! (1 - e^-tau)/tau without cancellation at small tau
        if (tau < 1e-4_dp) then
            frac = 1 - tau/2 + tau**2/6
        else
            frac = (1 - em)/tau
        end if
        e = (1 - refl)*frac/(1 - refl*em)
    end function escape_factor

    !> P/(n_e K omega_c^2) escaping between w = wmin and wmax, harmonics up to
    !> nmax; nmu points in the angle to the field, nphi in the azimuth. At
    !> angle theta, harmonic n is centred near w = n and ends at n/sin theta
    !> with a square-root edge, so the w intervals break at both, each split
    !> into nsub panels of 8 points clustered towards its upper end.
    real(dp) function slab_synchrotron_power(theta, lambda, refl, wmin, wmax, nmax, nsub, &
            nmu, nphi) result(p)
        real(dp), intent(in) :: theta, lambda, refl, wmin, wmax
        integer, intent(in) :: nmax, nsub, nmu, nphi
        real(dp) :: xw(8), ww(8), xm(nmu), wm(nmu), xp(nphi), wp(nphi)
        real(dp) :: ang, mu, st, lo, hi, a, b, v, w, jac, j, az, tau, e, pmu
        real(dp) :: brk(2*nmax + 2*ceiling(wmax) + 2)
        integer :: l, n, ip, i, k, nb, ib
        call gauss_legendre(8, xw, ww)
        call gauss_legendre(nmu, xm, wm)
        call gauss_legendre(nphi, xp, wp)
        p = 0
        do l = 1, nmu
            ! the emissivity is even in mu: theta in (0, pi/2), twice
            ang = pi/4*(xm(l) + 1)
            mu = cos(ang)
            st = sin(ang)
            pmu = 0
            nb = 0
            call add(wmin)
            call add(wmax)
            do n = 1, ceiling(wmax)
                call add(real(n, dp))
            end do
            do n = 1, nmax
                call add(n/st)
            end do
            call sort(brk(1:nb))
            do ib = 1, nb - 1
                lo = brk(ib)
                hi = brk(ib + 1)
                if (hi <= lo) cycle
                do ip = 1, nsub
                    ! v in (0, 1) from hi down to lo, w = hi - (hi - lo) v^2
                    a = real(ip - 1, dp)/nsub
                    b = real(ip, dp)/nsub
                    do i = 1, 8
                        v = a + (b - a)*(xw(i) + 1)/2
                        w = hi - (hi - lo)*v**2
                        jac = (b - a)/2*2*(hi - lo)*v
                        j = thermal_emissivity(theta, w, mu, nmax, 32)
                        if (j <= 0) cycle
                        e = 0
                        do k = 1, nphi
                            ! the azimuth in (0, pi/2), four times
                            az = pi/4*(xp(k) + 1)
                            tau = pi*lambda*j/(theta*w**2*st*cos(az))
                            e = e + wp(k)*pi/4*4*escape_factor(tau, refl)
                        end do
                        pmu = pmu + ww(i)*jac*j*e
                    end do
                end do
            end do
            p = p + wm(l)*pi/4*2*st*pmu
        end do
    contains
        subroutine add(x)
            real(dp), intent(in) :: x
            if (x < wmin .or. x > wmax) return
            nb = nb + 1
            brk(nb) = x
        end subroutine add
    end function slab_synchrotron_power

    pure subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: t
        integer :: i, k
        do i = 2, size(x)
            t = x(i)
            k = i - 1
            do while (k >= 1)
                if (x(k) <= t) exit
                x(k + 1) = x(k)
                k = k - 1
            end do
            x(k + 1) = t
        end do
    end subroutine sort

end module milca_synchrotron_transport
