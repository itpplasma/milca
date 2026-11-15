program test_synchrotron_transport
    ! Slab escape of synchrotron radiation against its derived limits: a
    ! transparent slab emits the Juttner-averaged Larmor power
    ! (8 pi/3) theta K3/K2 whatever the wall reflectivity; where it is opaque,
    ! the walls carry the blackbody flux (theta/(pi Lambda)) 2 pi int w^2 dw.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_synchrotron_transport
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: theta, ref, p0, p9, pt, bb

    theta = 0.02_dp
    ref = 8*pi/3*theta*scaled_bessel_k(3, 1/theta)/scaled_bessel_k(2, 1/theta)
    p0 = slab_synchrotron_power(theta, 1e-12_dp, 0.0_dp, 0.0_dp, 16.0_dp, 16, 2, 24, 8)
    p9 = slab_synchrotron_power(theta, 1e-12_dp, 0.9_dp, 0.0_dp, 16.0_dp, 16, 2, 24, 8)
    print '(a, 2es11.3)', 'thin slab over Larmor power - 1, R = 0 and 0.9:', p0/ref - 1, p9/ref - 1
    call require(abs(p0/ref - 1) < 1e-5_dp .and. abs(p9/ref - 1) < 1e-5_dp, 'transparent slab')

    ! the first harmonic band of a 25 keV plasma at large opacity
    theta = 0.05_dp
    pt = slab_synchrotron_power(theta, 1e12_dp, 0.0_dp, 0.95_dp, 1.0_dp, 4, 2, 32, 16)
    bb = theta/(pi*1e12_dp)*2*pi*(1.0_dp - 0.95_dp**3)/3
    print '(a, es11.3)', 'opaque band over blackbody flux - 1:', pt/bb - 1
    call require(abs(pt/bb - 1) < 1e-3_dp, 'opaque slab radiates as a blackbody')
    print '(a)', 'synchrotron escape matches its derived limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_synchrotron_transport
