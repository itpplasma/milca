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
        f = (mu - bpar)**2/(1 - mu**2)*bessel_jn(n, x)**2 + &
            bperp**2*((bessel_jn(n - 1, x) - bessel_jn(n + 1, x))/2)**2
    end function strength

    !> j(w, mu) per electron of a Juttner distribution at theta = kT/(m c^2),
    !> harmonics nmin..nmax (nmin = 1 by default), nq points along the
    !> resonance. On it gamma = s + mu p_par is linear in p_par, so the
    !> Juttner weight decays from the low-energy end over theta/|mu|; the
    !> quadrature starts there and stops after 60 e-foldings.
    real(dp) function thermal_emissivity(theta, w, mu, nmax, nq, nmin) result(j)
        real(dp), intent(in) :: theta, w, mu
        integer, intent(in) :: nmax, nq
        integer, intent(in), optional :: nmin
        real(dp) :: x(nq), wq(nq), s, s2, q, pa, pb, plo, dir, len, ppar, pperp2, gam, norm
        integer :: n, i, n0
        call gauss_legendre(nq, x, wq)
        n0 = 1
        if (present(nmin)) n0 = nmin
        s2 = 1 - mu**2
        norm = 4*pi*theta*scaled_bessel_k(2, 1/theta)
        j = 0
        do n = n0, nmax
            s = n/w
            if (s**2 <= s2) cycle
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
                j = j + wq(i)*len/2*exp(-(gam - 1)/theta)/norm*gam**2/s* &
                    strength(n, w*sqrt(pperp2*s2), ppar/gam, sqrt(pperp2)/gam, mu)*n
            end do
        end do
        j = j*2*pi
    end function thermal_emissivity

end module milca_cyclotron_emission
