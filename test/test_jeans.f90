program test_jeans
    ! The kinetic Jeans instability of a self-gravitating Maxwellian (G3 of
    ! docs/gravity.md), the plasma response with the coupling reversed: the
    ! root makes the gravitational dielectric 1 - (1 + zeta Z(zeta))/k^2 vanish
    ! with the series Z of milca_vlasov_response; derived limits gamma = 1 -
    ! (3/2) k^2 for small k and gamma = sqrt(2) (1 - k^2)/sqrt(pi) near the
    ! Jeans length k = 1; the initial-value response grows at gamma. Units
    ! omega_J, v_th and the Jeans length.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vlasov_response, only: plasma_dispersion, jeans_root, vlasov_density_response
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: gamma, k, h, slope
    complex(dp) :: zeta, z, dz, eps
    real(dp), allocatable :: n(:)
    integer :: m
    logical :: ok

    k = 0.5_dp
    call jeans_root(k, gamma, ok)
    call require(ok .and. gamma > 0, 'Jeans root at k = 0.5')
    zeta = cmplx(0.0_dp, gamma/(sqrt(2.0_dp)*k), dp)
    call plasma_dispersion(zeta, z, dz)
    eps = 1 - (1 + zeta*z)/k**2
    print '(a, f12.9, es10.2)', 'gamma(k = 0.5) and |epsilon|:', gamma, abs(eps)
    call require(abs(eps) < 1e-12_dp, 'gravitational dielectric vanishes at the root')

    k = 1e-2_dp
    call jeans_root(k, gamma, ok)
    call require(abs(gamma - (1 - 1.5_dp*k**2)) < 10*k**4, 'cold Jeans limit')
    k = 0.999_dp
    call jeans_root(k, gamma, ok)
    call require(abs(gamma/(sqrt(2/pi)*(1 - k**2)) - 1) < 1e-2_dp, 'marginal Jeans length')
    call jeans_root(1.0_dp, gamma, ok)
    call require(.not. ok, 'no growth at and beyond the Jeans length')

    ! initial value: the late-time logarithmic slope of the density
    k = 0.5_dp
    call jeans_root(k, gamma, ok)
    h = 0.01_dp
    m = nint(30/gamma/h)
    allocate (n(m))
    call vlasov_density_response(k, h, n, coupling=-1.0_dp)
    slope = (log(n(m)) - log(n(m - 100)))/(100*h)
    print '(a, 2f12.8)', 'late-time growth and Jeans root:', slope, gamma
    call require(abs(slope/gamma - 1) < 1e-4_dp, 'initial-value response grows at the root')
    print '(a)', 'Jeans instability matches its derived limits'

contains

    subroutine require(cond, what)
        logical, intent(in) :: cond
        character(*), intent(in) :: what
        if (.not. cond) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_jeans
