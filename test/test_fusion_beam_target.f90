program test_fusion_beam_target
    ! Beam-target reactivity (N4) against derived limits: on a cold target it
    ! is sigma(E_cm) v_beam; averaged over a beam Maxwellian of the target's
    ! temperature it is the Maxwellian reactivity of N3 (with the reduced mass),
    ! for D-T from the S factor of Bosch and Hale.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_fusion_reactivity, only: beam_target_reactivity, maxwellian_reactivity
    use milca_fusion_reactivity_leaves, only: bosch_hale_s
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), c_cm = 2.99792458e10_dp, barn = 1e-24_dp
    real(dp), parameter :: dt_bg = 34.3827_dp, md = 1875612.9_dp, mt = 2808921.1_dp
    real(dp), parameter :: dt_a(5) = [6.927e4_dp, 7.454e8_dp, 2.050e6_dp, 5.2002e4_dp, 0.0_dp]
    real(dp), parameter :: dt_b(4) = [63.8_dp, -0.995_dp, 6.981e-5_dp, 1.728e-4_dp]
    real(dp) :: none(0), mu, eb, cold, ref, t, avg, maxw, x(64), w(64), e, ecm
    integer :: i

    mu = md*mt/(md + mt)
    ! a 100 keV deuteron beam on a 1 eV triton target
    eb = 100
    cold = beam_target_reactivity(eb, md, mt, 1e-3_dp, dt_bg, dt_s, 200)
    ecm = eb*mt/(md + mt)
    ref = dt_s(ecm)/ecm*exp(-dt_bg/sqrt(ecm))*barn*sqrt(2*eb/md)*c_cm
    print '(a, 2es14.6)', 'cold target and sigma(E_cm) v_b:', cold, ref
    call require(abs(cold/ref - 1) < 1e-4_dp, 'cold-target limit')

    ! a beam Maxwellian of temperature t: the beam energy distribution
    ! 2 sqrt(E/pi) T^(-3/2) exp(-E/T), E = T y^2
    t = 10
    call gauss_legendre(64, x, w)
    avg = 0
    do i = 1, 64
        e = t*(5*(x(i) + 1)/2)**2
        avg = avg + w(i)*5.0_dp/2*2*(5*(x(i) + 1)/2)*t* &
            2*sqrt(e/pi)*t**(-1.5_dp)*exp(-e/t)*beam_target_reactivity(e, md, mt, t, dt_bg, dt_s, 200)
    end do
    maxw = maxwellian_reactivity(t, dt_bg, mu, dt_s, 550.0_dp, none, 0.05_dp)
    print '(a, 2es14.6)', 'beam Maxwellian average and Maxwellian reactivity:', avg, maxw
    call require(abs(avg/maxw - 1) < 1e-4_dp, 'Maxwellian average reproduces N3')
    print '(a)', 'beam-target reactivity matches its limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    real(dp) function dt_s(e)
        real(dp), intent(in) :: e
        call bosch_hale_s(e, dt_a(1), dt_a(2), dt_a(3), dt_a(4), dt_a(5), dt_b(1), dt_b(2), &
            dt_b(3), dt_b(4), dt_s)
        dt_s = dt_s*1e-3_dp
    end function dt_s

end program test_fusion_beam_target
