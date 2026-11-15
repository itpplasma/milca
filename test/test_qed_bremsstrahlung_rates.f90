program test_qed_bremsstrahlung_rates
    ! Electron-nucleus bremsstrahlung rates against their limits: the photon
    ! spectrum approaches the ultrarelativistic Bethe-Heitler formula as the
    ! energy grows; the energy radiated per path is v (16/3) in the
    ! nonrelativistic Born limit (the integral of s ln((1+s)/(1-s)) over
    ! (0,1) is 1); the thermal power is (16/3) sqrt(8 theta/pi), the Born
    ! Gaunt factor 2 sqrt(3)/pi, with a defect linear in theta.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_bremsstrahlung
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), energies(2) = [20.0_dp, 100.0_dp]
    real(dp) :: e, omega, ep, ref, ratio(2), loss, p, d(2), theta(2)
    integer :: i

    ! ultrarelativistic Bethe-Heitler, unscreened
    do i = 1, 2
        e = energies(i)
        omega = 0.5_dp*(e - 1)
        ep = e - omega
        ref = 4*(e**2 + ep**2 - 2*e*ep/3)*(log(2*e*ep/omega) - 0.5_dp)/(e**2*omega)
        ratio(i) = nucleus_photon_spectrum(e, omega, 48)/ref
    end do
    print '(a, 2f9.5)', 'spectrum over ultrarelativistic Bethe-Heitler at E = 20, 100:', ratio
    call require(abs(ratio(2) - 1) < 2e-3_dp .and. abs(ratio(2) - 1) < abs(ratio(1) - 1), &
        'spectrum approaches the ultrarelativistic limit')

    p = 0.02_dp
    e = sqrt(1 + p**2)
    loss = nucleus_energy_loss(e, 24, 16)
    call require(abs(loss/(16.0_dp/3) - 1) < 1e-3_dp, 'nonrelativistic energy loss per path')

    theta = [1e-3_dp, 4e-3_dp]
    do i = 1, 2
        d(i) = nucleus_thermal_power(theta(i), 20, 12, 16, 40.0_dp)/ &
            (16.0_dp/3*sqrt(8*theta(i)/pi)) - 1
    end do
    print '(a, 2es10.2)', 'thermal power defect at theta = 1e-3, 4e-3:', d
    call require(abs(d(1)) < 1e-2_dp, 'Born thermal power, Gaunt factor 2 sqrt(3)/pi')
    call require(abs(d(2)/d(1) - 4) < 0.5_dp, 'thermal power defect linear in theta')
    print '(a)', 'bremsstrahlung rates match their limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_qed_bremsstrahlung_rates
