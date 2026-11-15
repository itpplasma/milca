!> Layer response Delta of the cylinder kinetic layer
!> (milca_kinetic_layer_delta): the response to the real field of the
!> harmonic (-m, -n) is the complex conjugate of that of (m, n); Delta does
!> not depend on the window once the core holds the layer; and a
!> collisional layer has the dissipative sign of the resistive layer (L5)
!> at the same Doppler frequency; and the resonant radius of a linear
!> |q| profile is found exactly.
program test_kinetic_layer_delta
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_delta, only: kl_delta, kl_resonance
    use milca_kinetic_layer_cylinder, only: KL_OK
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_window
    use milca_kinetic_layer_cylinder, only: kl_background_t
    use milca_mhd_layer_n, only: layer_params_t, layer_delta
    implicit none
    real(dp), parameter :: rm = 55
    integer :: failures = 0
    type(kl_profiles_t) :: p
    complex(dp) :: d1, d2, i1, i2
    integer :: st

    call profiles(p, 6, 2)
    call kl_delta(p, rm, 1.5_dp, 3.0_dp, 256, 96, 1, d1, i1, st)
    call check('solve', st == KL_OK)
    if (failures > 0) error stop 'FAIL'
    print '(a,2es14.6)', 'Delta [1/cm] ', d1
    call check('nonzero dissipative part', abs(aimag(d1)) > 1e-3_dp*abs(d1))

    call profiles(p, -6, -2)
    call kl_delta(p, rm, 1.5_dp, 3.0_dp, 256, 96, 1, d2, i2, st)
    call close_to('reality: Delta(-m,-n) = conj Delta(m,n)', d2, conjg(d1), 1e-10_dp)

    call profiles(p, 6, 2)
    call kl_delta(p, rm, 2.0_dp, 4.0_dp, 384, 144, 1, d2, i2, st)
    call close_to('window independence', d2, d1, 1e-2_dp)

    call resonance()
    call doppler_sign(0.3_dp)
    call doppler_sign(-0.3_dp)

    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    !> Deuterium plasma with gradients, a radial electric field and a
    !> resonance q = m/n at r = rm.
    subroutine profiles(p, m, n)
        type(kl_profiles_t), intent(out) :: p
        integer, intent(in) :: m, n
        integer, parameter :: np = 401
        integer :: i
        allocate (p%r(np), p%ne(np), p%te(np), p%ti(np), p%q(np), p%er(np))
        do i = 1, np
            p%r(i) = 40 + 30*(i - 1)/real(np - 1, dp)
        end do
        p%ne = 3e13_dp*exp(-(p%r - 45)/30)
        p%te = 800*exp(-(p%r - 45)/25)
        p%ti = 900*exp(-(p%r - 45)/20)
        p%q = -3*(1 + 0.01_dp*(p%r - rm))
        p%er = -0.3_dp + 0.01_dp*(p%r - rm)
        p%btor = -18000; p%r0 = 165; p%m = m; p%n = n; p%zi = 1; p%ai = 2
    end subroutine profiles

    !> q = -3 (1 + 0.01 (r - rm)): |q| = 3.03 at rm + 1, slope 0.03.
    subroutine resonance()
        type(kl_profiles_t) :: f
        real(dp) :: r, dq
        integer :: s1
        call profiles(f, 6, 2)
        call kl_resonance(f, 3.03_dp, r, dq, s1)
        call check('resonance found', s1 == KL_OK)
        call check('resonant radius', abs(r - (rm + 1)) < 1e-9_dp)
        call check('shear', abs(dq - 0.03_dp) < 1e-9_dp)
    end subroutine resonance

    !> Flat collisional plasma (no diamagnetic frequencies) in the radial
    !> field er: the static field is seen at the Doppler frequency omega_E,
    !> and Im Delta has the sign of the resistive-inertial layer with
    !> gamma = -i (0 - omega_E) (L5 convention gamma = -i (omega - n Omega)).
    subroutine doppler_sign(er)
        real(dp), intent(in) :: er
        type(kl_profiles_t) :: f
        type(kl_background_t) :: bg
        type(layer_params_t) :: par
        complex(dp) :: dk, dr, ik
        integer :: s1, s2
        call profiles(f, 6, 2)
        f%ne = 5e13_dp; f%te = 30; f%ti = 30; f%er = er
        call kl_window(f, rm, 1.5_dp, 3.0_dp, 256, bg, s1)
        call kl_delta(f, rm, 1.5_dp, 3.0_dp, 256, 96, 1, dk, ik, s2)
        call check('collisional solve', s1 == 0 .and. s2 == KL_OK)
        par = layer_params_t(cmplx(0.0_dp, bg%omega_e(129), dp), 1.0_dp, 1.0e-7_dp, &
            0.0_dp, 1.0_dp)
        dr = layer_delta(par, s1)
        print '(a,es10.2,a,2es12.4,a,2es12.4)', 'omega_E ', bg%omega_e(129), &
            '  Delta kinetic ', dk, '  resistive ', dr
        call check('dissipative sign of the collisional layer', &
            sign(1.0_dp, aimag(dk)) == sign(1.0_dp, aimag(dr)))
    end subroutine doppler_sign

    subroutine close_to(name, got, want, tol)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: tol
        if (abs(got - want) <= tol*abs(want)) return
        failures = failures + 1
        print '(a,a,4es16.7)', 'FAIL ', name, got, want
    end subroutine close_to

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_kinetic_layer_delta
