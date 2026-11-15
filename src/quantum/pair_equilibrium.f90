! Electron-positron pairs in chemical equilibrium with photons (mu_+ + mu_- = 0,
! photons at mu_gamma = 0) at temperature theta = kT/(m c^2), with the net
! charge density n_- - n_+ = n_charge fixed by the ions. Fermi-Dirac
! occupation, spin degeneracy 2; densities in lambda_C^-3, lambda_C = hbar/(m c),
! chemical potentials including the rest energy, in m c^2. The nondegenerate
! limit is the Juttner density theta K2(1/theta)/pi^2 per species and the
! law of mass action n_+ n_- = n_eq^2; the ultrarelativistic limit is
! 3 zeta(3) theta^3/(2 pi^2).
module milca_pair_equilibrium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: fermion_density, pair_equilibrium

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: nq = 160

contains

    !> Density of spin-1/2 fermions of unit mass with chemical potential mu:
    !> (1/pi^2) int p^2 dp/(exp((E - mu)/theta) + 1), E = sqrt(1 + p^2);
    !> the antiparticles of the same plasma have -mu.
    real(dp) function fermion_density(theta, mu) result(n)
        real(dp), intent(in) :: theta, mu
        real(dp) :: x(nq), w(nq), pf, pmax, emax
        call gauss_legendre(nq, x, w)
        ! split at the Fermi momentum, where the occupation steps
        pf = 0
        if (mu > 1) pf = sqrt(mu**2 - 1)
        emax = max(mu, 1.0_dp) + 60*theta
        pmax = sqrt(emax**2 - 1)
        n = (piece(0.0_dp, pf) + piece(pf, pmax))/pi**2
    contains
        real(dp) function piece(a, b) result(s)
            real(dp), intent(in) :: a, b
            real(dp) :: p, arg
            integer :: i
            s = 0
            if (b <= a) return
            do i = 1, nq
                p = a + (b - a)*(x(i) + 1)/2
                arg = (sqrt(1 + p**2) - mu)/theta
                s = s + w(i)*(b - a)/2*p**2*occupation(arg)
            end do
        end function piece
    end function fermion_density

    !> Positron and electron densities at temperature theta for the net
    !> charge density n_charge >= 0 (bisection on the electron chemical
    !> potential to rounding).
    subroutine pair_equilibrium(theta, n_charge, n_plus, n_minus)
        real(dp), intent(in) :: theta, n_charge
        real(dp), intent(out) :: n_plus, n_minus
        real(dp) :: lo, hi, mu
        integer :: it
        lo = 0
        hi = 1
        do while (excess(hi) < n_charge)
            hi = 2*hi
        end do
        do it = 1, 200
            mu = (lo + hi)/2
            if (excess(mu) < n_charge) then
                lo = mu
            else
                hi = mu
            end if
            if (hi - lo <= 4*epsilon(1.0_dp)*hi) exit
        end do
        mu = (lo + hi)/2
        n_minus = fermion_density(theta, mu)
        n_plus = fermion_density(theta, -mu)
    contains
        real(dp) function excess(m)
            real(dp), intent(in) :: m
            excess = fermion_density(theta, m) - fermion_density(theta, -m)
        end function excess
    end subroutine pair_equilibrium

    pure real(dp) function occupation(arg)
        real(dp), intent(in) :: arg
        if (arg > 700) then
            occupation = 0
        else if (arg > 0) then
            occupation = exp(-arg)/(1 + exp(-arg))
        else
            occupation = 1/(1 + exp(arg))
        end if
    end function occupation

end module milca_pair_equilibrium
