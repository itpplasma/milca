! Cyclotron and synchrotron emission of electrons in a uniform magnetic field,
! harmonic by harmonic (Schott 1912; Bekefi, Radiation Processes in Plasmas,
! 1966), with vacuum propagation and both polarizations summed. Units: the
! nonrelativistic cyclotron frequency omega_c = e B/m, momenta in m c, powers
! in K omega_c^2, K = e^2/(8 pi^2 epsilon_0 c). An electron with momentum
! (p_par, p_perp) emits harmonic n at angle theta to B (mu = cos theta) with
! frequency w_n = n/(gamma - p_par mu) (in omega_c) the power per unit solid
! angle and per unit particle time
!   dP_n/dOmega = (w_n gamma)^2 [((mu - beta_par)^2/(1 - mu^2)) J_n^2(x)
!                 + beta_perp^2 J_n'^2(x)]/(gamma^2 (1 - beta_par mu)),
! x = w_n p_perp sqrt(1 - mu^2); summed over n and directions it is the Larmor
! power (4 pi/3) p_perp^2. The thermal emissivity j(w, mu) of a Juttner
! distribution resolves the resonance delta(w - w_n) in p_perp; with s = n/w,
!   j = 2 pi sum_n n int dp_par f(p) gamma^2 F_n/s,  gamma = s + mu p_par,
! over the p_par where p_perp^2 = gamma^2 - 1 - p_par^2 >= 0, per electron and
! per unit w and solid angle, F_n the bracket above.
module milca_cyclotron_emission
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    private

    public :: harmonic_power, thermal_emissivity

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    !> dP_n/dOmega of one electron and its frequency w_n.
    subroutine harmonic_power(n, p_par, p_perp, mu, power, w)
        integer, intent(in) :: n
        real(dp), intent(in) :: p_par, p_perp, mu
        real(dp), intent(out) :: power, w
        real(dp) :: gam, bpar, bperp, x, d
        gam = sqrt(1 + p_par**2 + p_perp**2)
        bpar = p_par/gam
        bperp = p_perp/gam
        d = 1 - bpar*mu
        w = n/(gam*d)
        x = w*p_perp*sqrt(1 - mu**2)
        power = (w*gam)**2*strength(n, x, bpar, bperp, mu)/(gam**2*d)
    end subroutine harmonic_power

    pure real(dp) function strength(n, x, bpar, bperp, mu) result(f)
        integer, intent(in) :: n
        real(dp), intent(in) :: x, bpar, bperp, mu
        real(dp) :: jb(3)
        ! J_(n-1), J_n, J_(n+1) from one recurrence
        jb = bessel_jn(n - 1, n + 1, x)
        f = (mu - bpar)**2/(1 - mu**2)*jb(2)**2 + bperp**2*((jb(1) - jb(3))/2)**2
    end function strength

    !> j(w, mu) per electron of a Juttner distribution at theta = kT/(m c^2),
    !> harmonics nmin..nmax (nmin = 1 by default), nq points along each
    !> resonance. Only harmonics that resonate contribute: s = n/w lies
    !> between sin(theta_B) and gamma_max (1 + |mu|), gamma_max = 1 + 60 theta.
    !> On a resonance gamma = s + mu p_par is linear in p_par, so the Juttner
    !> weight decays from the low-energy end over theta/|mu|; the quadrature
    !> starts there and stops after 60 e-foldings. Where more than nexact + 2
    !> nsamp harmonics resonate (overlapping lines of a hot plasma), the
    !> summand is smooth in n beyond the first nexact, which are summed
    !> exactly, and the rest is a trapezoid sum over nsamp of them.
    real(dp) function thermal_emissivity(theta, w, mu, nmax, nq, nmin) result(j)
        real(dp), intent(in) :: theta, w, mu
        integer, intent(in) :: nmax, nq
        integer, intent(in), optional :: nmin
        integer, parameter :: nexact = 16, nsamp = 96
        real(dp) :: x(nq), wq(nq), s2, norm, smax, hp, hn, first, tsum
        integer :: n, n0, nlo, nhi, nn, k
        call gauss_legendre(nq, x, wq)
        n0 = 1
        if (present(nmin)) n0 = nmin
        s2 = 1 - mu**2
        norm = 4*pi*theta*scaled_bessel_k(2, 1/theta)
        smax = (1 + 60*theta)*(1 + abs(mu))
        nlo = max(n0, floor(w*sqrt(s2)) + 1)
        nhi = min(nmax, floor(w*smax))
        j = 0
        if (nhi < nlo) return
        if (nhi - nlo + 1 <= nexact + 2*nsamp) then
            do n = nlo, nhi
                j = j + harmonic(n)
            end do
        else
            do n = nlo, nlo + nexact - 1
                j = j + harmonic(n)
            end do
            ! sum_(n=a..b) f(n) = trapezoid integral + (f(a) + f(b))/2, on
            ! samples n_k = a + (b - a)(k/nsamp)^2 clustered at the lower end,
            ! where the summand rises like a square root
            n = nlo + nexact
            hp = harmonic(n)
            first = hp
            tsum = 0
            do k = 1, nsamp
                nn = nlo + nexact + nint((nhi - nlo - nexact)*(real(k, dp)/nsamp)**2)
                if (nn <= n) cycle
                hn = harmonic(nn)
                tsum = tsum + (nn - n)*(hp + hn)/2
                n = nn
                hp = hn
            end do
            j = j + tsum + (first + hp)/2
        end if
        j = j*2*pi
    contains
        real(dp) function harmonic(n) result(h)
            integer, intent(in) :: n
            real(dp) :: s, q, pa, pb, plo, dir, len, ppar, pperp2, gam
            integer :: i
            h = 0
            s = n/w
            if (s**2 <= s2) return
            q = sqrt(s**2 - s2)
            pa = (s*mu - q)/s2
            pb = (s*mu + q)/s2
            ! gamma grows with mu p_par: start at the end with the smaller gamma
            plo = pa
            dir = 1
            if (mu < 0) then
                plo = pb
                dir = -1
            end if
            len = pb - pa
            if (abs(mu)*len > 60*theta) len = 60*theta/abs(mu)
            do i = 1, nq
                ppar = plo + dir*len*(x(i) + 1)/2
                gam = s + mu*ppar
                pperp2 = max(gam**2 - 1 - ppar**2, 0.0_dp)
                h = h + wq(i)*len/2*exp(-(gam - 1)/theta)/norm*gam**2/s* &
                    strength(n, w*sqrt(pperp2*s2), ppar/gam, sqrt(pperp2)/gam, mu)*n
            end do
        end function harmonic
    end function thermal_emissivity

end module milca_cyclotron_emission
