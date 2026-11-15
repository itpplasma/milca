program test_fusion_reactivity
    ! Maxwellian reactivities from S factors: the D-T S factor of Bosch and
    ! Hale (Table IV) integrated over a Maxwellian against their independent
    ! reactivity parameterization (Table VII); the three pieces of the p-11B
    ! S factor of Wang et al. (2026) join continuously; the p-11B reactivity
    ! is converged in the quadrature, the 148 keV resonance included.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_fusion_reactivity
    use milca_fusion_reactivity_leaves, only: bosch_hale_s, bosch_hale_reactivity, &
        pb11_s_low, pb11_s_mid, pb11_s_high
    implicit none
    ! Bosch and Hale (1992), D + T -> 4He + n
    real(dp), parameter :: dt_bg = 34.3827_dp, dt_mrc2 = 1124656.0_dp
    real(dp), parameter :: dt_a(5) = [6.927e4_dp, 7.454e8_dp, 2.050e6_dp, 5.2002e4_dp, 0.0_dp]
    real(dp), parameter :: dt_b(4) = [63.8_dp, -0.995_dp, 6.981e-5_dp, 1.728e-4_dp]
    real(dp), parameter :: dt_c(7) = [1.17302e-9_dp, 1.51361e-2_dp, 7.51886e-2_dp, &
        4.60643e-3_dp, 1.35e-2_dp, -1.0675e-4_dp, 1.366e-5_dp]
    real(dp), parameter :: temps(5) = [2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp]
    real(dp), parameter :: pb_temps(4) = [50.0_dp, 150.0_dp, 300.0_dp, 600.0_dp]
    real(dp) :: none(0), integrated, fit, lo, hi, fine, coarse
    integer :: i

    do i = 1, size(temps)
        integrated = maxwellian_reactivity(temps(i), dt_bg, dt_mrc2, dt_s, 550.0_dp, none, 0.05_dp)
        call bosch_hale_reactivity(temps(i), dt_c(1), dt_c(2), dt_c(3), dt_c(4), dt_c(5), &
            dt_c(6), dt_c(7), dt_bg, dt_mrc2, fit)
        print '(a, f6.1, 2es14.6, f9.5)', 'D-T T, integrated, Table VII, ratio:', temps(i), &
            integrated, fit, integrated/fit
        ! Tables IV and VII are separate fits to the same R-matrix results,
        ! each good to a fraction of a percent
        call require(abs(integrated/fit - 1) < 1e-2_dp, "D-T reactivity against Bosch-Hale fit")
    end do

    call pb11_s_low(400.0_dp, lo)
    call pb11_s_mid(400.0_dp, hi)
    ! the table gives four significant digits: 330.247 against 330.2 MeV b
    call require(abs(lo/hi - 1) < 5e-4_dp, 'p-11B S factor continuous at 400 keV')
    call pb11_s_mid(700.0_dp, lo)
    call pb11_s_high(700.0_dp, hi)
    call require(abs(lo/hi - 1) < 1e-3_dp, 'p-11B S factor continuous at 700 keV')

    do i = 1, size(pb_temps)
        coarse = maxwellian_reactivity(pb_temps(i), PB11_GAMOW_SQRT_KEV, PB11_MRC2_KEV, &
            pb11_s_factor, 10000.0_dp, pb11_breaks, 0.1_dp)
        fine = maxwellian_reactivity(pb_temps(i), PB11_GAMOW_SQRT_KEV, PB11_MRC2_KEV, &
            pb11_s_factor, 10000.0_dp, pb11_breaks, 0.025_dp)
        print '(a, f6.1, es14.6)', 'p-11B T, <sigma v> (cm^3/s):', pb_temps(i), fine
        call require(abs(coarse/fine - 1) < 1e-8_dp, 'p-11B reactivity converged')
    end do
    print '(a)', 'fusion reactivities match their oracles'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    ! keV b from the Pade form in keV mb
    real(dp) function dt_s(e)
        real(dp), intent(in) :: e
        call bosch_hale_s(e, dt_a(1), dt_a(2), dt_a(3), dt_a(4), dt_a(5), dt_b(1), dt_b(2), &
            dt_b(3), dt_b(4), dt_s)
        dt_s = dt_s*1e-3_dp
    end function dt_s

end program test_fusion_reactivity
