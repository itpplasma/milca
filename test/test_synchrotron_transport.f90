program test_synchrotron_transport
    ! Slab escape of synchrotron radiation: where the slab is opaque the
    ! walls carry the blackbody flux (theta/(pi Lambda)) 2 pi int w^2 dw; the
    ! escaping power as the Larmor power minus the absorbed part equals the
    ! direct integral of the emission times the escape factor; on a hot
    ! plasma the strided sum over overlapping harmonics matches the sum of
    ! every harmonic.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_synchrotron_transport
    use milca_cyclotron_emission, only: thermal_emissivity
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: theta, pt, bb, direct, deficit, strided, exact, s0, smax, w, mu, larmor
    integer :: n

    theta = 0.05_dp
    pt = slab_band_power(theta, 1e12_dp, 0.0_dp, 0.95_dp, 1.0_dp, 4, 2, 32, 16)
    bb = theta/(pi*1e12_dp)*2*pi*(1.0_dp - 0.95_dp**3)/3
    print '(a, es11.3)', 'opaque band over blackbody flux - 1:', pt/bb - 1
    call require(abs(pt/bb - 1) < 1e-3_dp, 'opaque slab radiates as a blackbody')

    theta = 0.02_dp
    direct = slab_band_power(theta, 10.0_dp, 0.9_dp, 0.0_dp, 16.0_dp, 16, 2, 24, 8)
    deficit = slab_synchrotron_power(theta, 10.0_dp, 0.9_dp, 24, 8, 1e-12_dp)
    ! both carry the emission quadrature error, about 1e-5 of the Larmor power
    larmor = 8*pi/3*theta*scaled_bessel_k(3, 1/theta)/scaled_bessel_k(2, 1/theta)
    print '(a, es11.3)', '(Larmor minus absorbed - direct escape)/Larmor:', (deficit - direct)/larmor
    call require(abs(deficit - direct) < 2e-5_dp*larmor, 'absorbed-part form of the escaping power')

    theta = 0.4_dp
    w = 40
    mu = 0.3_dp
    strided = thermal_emissivity(theta, w, mu, 1000000, 32)
    s0 = sqrt(1 - mu**2)
    smax = (1 + 60*theta)*(1 + mu)
    exact = 0
    do n = floor(w*s0) + 1, floor(w*smax)
        exact = exact + thermal_emissivity(theta, w, mu, n, 32, nmin=n)
    end do
    print '(a, es11.3)', 'strided over exact harmonic sum - 1:', strided/exact - 1
    call require(abs(strided/exact - 1) < 1e-4_dp, 'strided harmonic sum')
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
