! Radiative recombination into the ground state of a hydrogenic two-body
! system of reduced mass mu (hydrogen, positronium, muonium), from the exact
! 1s photoionization cross section of Stobbe (Ann. Phys. 399, 661, 1930),
!   sigma_PI = (2^9 pi^2/3) alpha a^2 (1 + eps)^-4 exp(-4 nu arccot nu)/(1 - exp(-2 pi nu)),
! eps = E/I the relative kinetic energy over the binding energy
! I = mu c^2 alpha^2/2, nu = eps^(-1/2), a = hbar/(mu c alpha), and the Milne
! relation sigma_rec = 2 (k_gamma/k_e)^2 sigma_PI, where the 2 is
! g_bound g_photon/(g_free pair) for hydrogen (2 2/2) and positronium (4 2/4).
! The Maxwellian rate over the relative energy at tau = kT/I is
!   <sigma v> = a^2 alpha^4 c R(tau),
!   R = (2^9 pi^2/6) (2/sqrt(pi)) tau^(-3/2) int exp(-eps/tau) (1 + eps)^-2 g(eps) d eps,
! g = exp(-4 nu arccot nu)/(1 - exp(-2 pi nu)). Since a^2 scales as mu^-2 and
! I as mu, positronium (mu = m/2) recombines as 4 alpha_H(T I_H/I_Ps).
module milca_hydrogenic_recombination
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: ground_state_recombination

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    !> R(tau) of the ground-state recombination rate, n quadrature points.
    real(dp) function ground_state_recombination(tau, n) result(r)
        real(dp), intent(in) :: tau
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), eps, emax, root, g
        integer :: i
        call gauss_legendre(n, x, w)
        emax = 60*tau
        r = 0
        do i = 1, n
            eps = emax*(x(i) + 1)/2
            root = sqrt(eps)
            ! nu arccot nu = atan(sqrt(eps))/sqrt(eps), finite at threshold
            g = exp(-4*atan(root)/root)/(1 - exp(-2*pi/root))
            r = r + w(i)*emax/2*exp(-eps/tau)*g/(1 + eps)**2
        end do
        r = r*2.0_dp**9*pi**2/6*2/sqrt(pi)*tau**(-1.5_dp)
    end function ground_state_recombination

end module milca_hydrogenic_recombination
