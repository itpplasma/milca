program test_positronium_formation
    ! Ground-state radiative recombination of hydrogenic systems: hydrogen at
    ! 10^4 K against the tabulated case-A minus case-B coefficient
    ! 1.58e-13 cm^3/s (Osterbrock and Ferland, Astrophysics of Gaseous
    ! Nebulae, 2006, table 2.1); the threshold limit
    ! R -> (2^9 pi^2/6) e^-4 (2/sqrt(pi)) tau^(-1/2) at low temperature;
    ! positronium as 4 alpha_H at twice the temperature.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_hydrogenic_recombination, only: ground_state_recombination
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), alpha = 7.2973525693e-3_dp
    real(dp), parameter :: a0_cm = 5.29177210903e-9_dp, c_cm = 2.99792458e10_dp
    ! hydrogen binding energy with the proton's reduced mass, and k_B, in eV
    real(dp), parameter :: i_h = 13.598434_dp, kb = 8.617333262e-5_dp
    real(dp), parameter :: mu_h = 1/(1 + 1/1836.15267343_dp)
    real(dp) :: rate, tau, lim

    tau = kb*1e4_dp/i_h
    rate = (a0_cm/mu_h)**2*alpha**4*c_cm*ground_state_recombination(tau, 96)
    print '(a, es12.4)', 'hydrogen ground-state recombination at 1e4 K (cm^3/s):', rate
    call require(abs(rate/1.58e-13_dp - 1) < 1e-2_dp, 'hydrogen at 10^4 K')

    tau = 1e-5_dp
    lim = 2.0_dp**9*pi**2/6*exp(-4.0_dp)*2/sqrt(pi)/sqrt(tau)
    print '(a, es10.2)', 'threshold limit defect at tau = 1e-5:', ground_state_recombination(tau, 96)/lim - 1
    call require(abs(ground_state_recombination(tau, 96)/lim - 1) < 1e-3_dp, 'threshold limit')

    ! positronium: a = 2 a0, I = 6.8 eV
    tau = kb*5e3_dp/(i_h/mu_h/2)
    rate = (2*a0_cm)**2*alpha**4*c_cm*ground_state_recombination(tau, 96)
    print '(a, es12.4)', 'positronium formation at 5000 K (cm^3/s):', rate
    call require(abs(rate/(4*1.58e-13_dp) - 1) < 1.5e-2_dp, 'positronium scaling')
    print '(a)', 'ground-state recombination matches hydrogen and its limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_positronium_formation
