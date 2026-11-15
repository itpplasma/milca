program test_qed_elwert
    ! The Coulomb-corrected electron-nucleus spectrum: it reduces to Born as
    ! the nuclear charge vanishes; at the tip, where Born vanishes like p',
    ! it approaches the derived nonrelativistic value
    ! d sigma/d omega = (32/3) Z^2 alpha r_e^2 (1 - exp(-2 pi nu))/(p^2 omega),
    ! nu = Z alpha/beta. Units m = c = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_bremsstrahlung, only: nucleus_photon_spectrum
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), alpha = 7.2973525693e-3_dp
    real(dp) :: e, p, omega, born, coul, tip, nu

    e = 1.2_dp
    omega = 0.1_dp
    born = nucleus_photon_spectrum(e, omega, 24)
    coul = nucleus_photon_spectrum(e, omega, 24, 1e-8_dp)
    call require(abs(coul/born - 1) < 1e-8_dp, 'Born limit for vanishing charge')

    p = 0.02_dp
    e = sqrt(1 + p**2)
    omega = (e - 1)*(1 - 1e-6_dp)
    nu = 5*alpha*e/p
    coul = nucleus_photon_spectrum(e, omega, 24, 5.0_dp)
    tip = 32.0_dp/3*(1 - exp(-2*pi*nu))/(p**2*omega)
    print '(a, f6.3, 2es14.6)', 'boron, nu, tip spectrum and derived tip:', nu, coul, tip
    call require(abs(coul/tip - 1) < 2e-3_dp, 'nonrelativistic tip of the Coulomb-corrected spectrum')
    print '(a)', 'Coulomb-corrected bremsstrahlung matches its limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_qed_elwert
