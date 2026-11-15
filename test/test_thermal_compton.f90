program test_thermal_compton
    ! Thermal Compton scattering against derived limits: the Thomson cross
    ! section 8 pi/3 r_e^2 and the gain (4/3) <p^2> = 4 theta K3/K2 per
    ! scattering as eps -> 0; the recoil -eps on cold electrons; the
    ! Kompaneets energy exchange 4 theta - eps when both are small.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_thermal_compton, only: compton_rates
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: rate, heat, theta, ref, d1, d2

    theta = 0.3_dp
    call compton_rates(1e-7_dp, theta, 32, rate, heat)
    ref = 4*theta*scaled_bessel_k(3, 1/theta)/scaled_bessel_k(2, 1/theta)
    print '(a, 2es11.3)', 'Thomson rate and gain per scattering, relative defects:', &
        rate/(8*pi/3) - 1, heat/(rate*1e-7_dp)/ref - 1
    call require(abs(rate/(8*pi/3) - 1) < 1e-5_dp, 'Thomson cross section')
    call require(abs(heat/(rate*1e-7_dp)/ref - 1) < 1e-5_dp, 'Thomson gain (4/3) <p^2>')

    call compton_rates(1e-3_dp, 1e-9_dp, 32, rate, heat)
    d1 = heat/(rate*1e-3_dp)/(-1e-3_dp) - 1
    call compton_rates(5e-4_dp, 1e-9_dp, 32, rate, heat)
    d2 = heat/(rate*5e-4_dp)/(-5e-4_dp) - 1
    print '(a, 2es11.3)', 'cold recoil over -eps - 1 at eps = 1e-3, 5e-4:', d1, d2
    call require(abs(d1) < 1e-2_dp .and. abs(d2/d1 - 0.5_dp) < 0.05_dp, 'recoil limit')

    call compton_rates(2e-3_dp, 1e-3_dp, 32, rate, heat)
    d1 = heat/(rate*2e-3_dp)/(4e-3_dp - 2e-3_dp) - 1
    print '(a, es11.3)', 'Kompaneets 4 theta - eps, relative defect:', d1
    call require(abs(d1) < 2e-2_dp, 'Kompaneets limit')
    print '(a)', 'thermal Compton scattering matches its derived limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_thermal_compton
