! p-11B fusion power against bremsstrahlung at equal electron and ion
! temperature, per n_e^2 in W m^3, for boron fractions x = n_B/n_p; CSV on
! standard output. Fusion from the S factor of Wang et al. (2026), electron-ion
! bremsstrahlung from the Born amplitudes with the Coulomb correction of
! Elwert for each nuclear charge, electron-electron bremsstrahlung from the
! eight-diagram amplitudes.
! Synchrotron losses are not included.
!   pb11_power_balance T_keV ...
program pb11_power_balance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_fusion_reactivity, only: maxwellian_reactivity, pb11_s_factor, pb11_breaks, &
        PB11_GAMOW_SQRT_KEV, PB11_MRC2_KEV
    use milca_qed_bremsstrahlung, only: nucleus_thermal_power, electron_energy_radiated
    use milca_qed_thermal, only: pair_energy_rate
    implicit none
    ! CODATA 2018
    real(dp), parameter :: alpha = 7.2973525693e-3_dp, r_e = 2.8179403262e-15_dp
    real(dp), parameter :: mc2_j = 8.1871057769e-14_dp, c = 299792458.0_dp
    real(dp), parameter :: mc2_kev = 510.99895_dp, mev_j = 1.602176634e-13_dp
    ! p + 11B -> 3 alpha releases 8.68 MeV
    real(dp), parameter :: q_fusion = 8.68_dp*mev_j
    real(dp), parameter :: fractions(4) = [0.1_dp, 0.15_dp, 0.2_dp, 0.3_dp]
    real(dp) :: unit_power, t, theta, sigv, p_ep, p_eb, p_ee, np, nb, fus, brem
    integer :: i, j, nargs
    character(64) :: arg

    unit_power = alpha*r_e**2*mc2_j*c
    nargs = command_argument_count()
    print '(a)', 'T_keV,x,sigv_m3s,P_fus,P_ei,P_ee,ratio'
    do i = 1, nargs
        call get_command_argument(i, arg)
        read (arg, *) t
        theta = t/mc2_kev
        sigv = 1e-6_dp*maxwellian_reactivity(t, PB11_GAMOW_SQRT_KEV, PB11_MRC2_KEV, &
            pb11_s_factor, 10000.0_dp, pb11_breaks, 0.05_dp)
        ! per n_e n_i Z^2 and per n_e^2 (the 1/2 counts each electron pair once)
        p_ep = nucleus_thermal_power(theta, 24, 16, 12, 40.0_dp, 1.0_dp)*unit_power
        p_eb = nucleus_thermal_power(theta, 24, 16, 12, 40.0_dp, 5.0_dp)*unit_power
        p_ee = 0.5_dp*pair_energy_rate(theta, 16, radiated)*unit_power
        do j = 1, size(fractions)
            ! quasi-neutrality n_e = n_p + 5 n_B with n_e = 1
            np = 1/(1 + 5*fractions(j))
            nb = fractions(j)*np
            fus = np*nb*sigv*q_fusion
            brem = np*p_ep + 25*nb*p_eb + p_ee
            print '(f0.1, ",", f0.2, 5(",", es12.5))', t, fractions(j), sigv, fus, &
                np*p_ep + 25*nb*p_eb, p_ee, fus/brem
        end do
        flush (6)
    end do

contains

    real(dp) function radiated(s)
        real(dp), intent(in) :: s
        radiated = electron_energy_radiated(s, 10)
    end function radiated

end program pb11_power_balance
