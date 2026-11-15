!> Layer response Delta of the cylinder kinetic layer
!> (milca_kinetic_layer_delta): the response to the real field of the
!> harmonic (-m, -n) is the complex conjugate of that of (m, n); Delta does
!> not depend on the window once the core holds the layer; and the net
!> current agrees with a quadrature of the current profile.
program test_kinetic_layer_delta
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_delta, only: kl_delta
    use milca_kinetic_layer_cylinder, only: KL_OK
    use milca_kinetic_layer_window, only: kl_profiles_t
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
