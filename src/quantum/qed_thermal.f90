! Thermal rates of two-body QED processes for Juttner distributions:
! <sigma v_Moller> of electron-positron annihilation by the invariant
! reduction of Gondolo and Gelmini (Nucl. Phys. B 360, 145, 1991),
!   <sigma v> = int sigma(s) (s - 4) sqrt(s) K1(sqrt(s)/theta) ds
!               / (8 theta K2(1/theta)^2),   s >= 4,
! units m = c = 1, rate coefficients in r_e^2 c. The modified Bessel
! functions enter exponentially scaled, so low temperatures do not underflow.
module milca_qed_thermal
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_leaves, only: qed_dirac_annihilation_total
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: scaled_bessel_k, annihilation_rate, pair_energy_rate, kernel_s

    abstract interface
        !> A function of the invariant s of a colliding pair.
        real(dp) function kernel_s(s)
            import :: dp
            real(dp), intent(in) :: s
        end function kernel_s
    end interface

contains

    !> exp(x) K_nu(x) from K_nu(x) = int_0^inf exp(-x cosh t) cosh(nu t) dt,
    !> with t up to where x (cosh t - 1) reaches 60.
    real(dp) function scaled_bessel_k(nu, x) result(k)
        integer, intent(in) :: nu
        real(dp), intent(in) :: x
        integer, parameter :: n = 96
        real(dp) :: t(n), w(n), tmax, tt
        integer :: i
        call gauss_legendre(n, t, w)
        tmax = acosh(1 + 60/x)
        k = 0
        do i = 1, n
            tt = tmax*(t(i) + 1)/2
            k = k + w(i)*tmax/2*exp(-x*(cosh(tt) - 1))*cosh(nu*tt)
        end do
    end function scaled_bessel_k

    !> <sigma v_Moller> of e- e+ -> 2 gamma at temperature theta = kT/(m c^2).
    real(dp) function annihilation_rate(theta, n) result(rate)
        real(dp), intent(in) :: theta
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), y, root_s, s, sigma, gam, k2, smax, jac
        integer :: i
        call gauss_legendre(n, x, w)
        k2 = scaled_bessel_k(2, 1/theta)
        ! sqrt(s) = 2 + y^2 theta: smooth at threshold, the Boltzmann factor
        ! exp(-(sqrt(s) - 2)/theta) = exp(-y^2) cut at y^2 = 60
        smax = sqrt(60.0_dp)
        rate = 0
        do i = 1, n
            y = smax*(x(i) + 1)/2
            root_s = 2 + y**2*theta
            s = root_s**2
            jac = smax/2*2*y*theta*2*root_s
            ! the positron Lorentz factor in the electron rest frame
            gam = s/2 - 1
            if (gam <= 1) cycle
            call qed_dirac_annihilation_total(gam, sigma)
            rate = rate + w(i)*jac*sigma*(s - 4)*root_s*scaled_bessel_k(1, root_s/theta)* &
                exp(-y**2)
        end do
        rate = rate/(8*theta*k2**2)
    end function annihilation_rate

    !> <v_Moller gamma_cm w(s)> over two Juttner distributions of equal mass:
    !> the lab energy of the pair is gamma_cm sqrt(s), so the reduction of
    !> Gondolo and Gelmini carries K2(sqrt(s)/theta) in place of K1. With w
    !> the centre-of-mass energy radiated per collision (a cross section
    !> times energy, symmetric emission) this is the lab power per pair
    !> density product.
    real(dp) function pair_energy_rate(theta, n, w) result(rate)
        real(dp), intent(in) :: theta
        integer, intent(in) :: n
        procedure(kernel_s) :: w
        real(dp) :: x(n), wq(n), y, root_s, s, k2, smax, jac
        integer :: i
        call gauss_legendre(n, x, wq)
        k2 = scaled_bessel_k(2, 1/theta)
        smax = sqrt(60.0_dp)
        rate = 0
        do i = 1, n
            y = smax*(x(i) + 1)/2
            root_s = 2 + y**2*theta
            s = root_s**2
            jac = smax/2*2*y*theta*2*root_s
            rate = rate + wq(i)*jac*w(s)*(s - 4)*root_s*scaled_bessel_k(2, root_s/theta)* &
                exp(-y**2)
        end do
        rate = rate/(8*theta*k2**2)
    end function pair_energy_rate

end module milca_qed_thermal
