program test_pair_equilibrium
    ! Thermal electron-positron pairs against derived limits: the Juttner
    ! density theta K2(1/theta)/pi^2 when nondegenerate; the ultrarelativistic
    ! density 3 zeta(3) theta^3/(2 pi^2); with ions the law of mass action
    ! n_+ n_- = n_eq^2 and charge neutrality. The gamma gamma -> e+ e- rate of a
    ! Wien photon gas equals the annihilation rate n_eq^2 <sigma v> of the
    ! Juttner pairs (detailed balance between the Breit-Wheeler and Dirac
    ! cross sections).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_pair_equilibrium, only: fermion_density, pair_equilibrium
    use milca_qed_thermal, only: scaled_bessel_k, annihilation_rate, photon_pair_rate
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), zeta3 = 1.2020569031595942854_dp
    real(dp) :: theta, neq, np, nm, nc, ref, prod, ann, ths(3)
    integer :: i

    theta = 0.05_dp
    neq = theta*scaled_bessel_k(2, 1/theta)*exp(-1/theta)/pi**2
    call pair_equilibrium(theta, 0.0_dp, np, nm)
    print '(a, 2es14.6)', 'theta = 0.05 pairs and Juttner density:', np, neq
    call require(abs(np/neq - 1) < 1e-8_dp .and. abs(nm/np - 1) < 1e-8_dp, 'nondegenerate limit')

    theta = 1000
    ref = 3*zeta3*theta**3/(2*pi**2)
    call pair_equilibrium(theta, 0.0_dp, np, nm)
    print '(a, f12.9)', 'theta = 1000 pairs over 3 zeta(3) theta^3/(2 pi^2):', np/ref
    call require(abs(np/ref - 1) < 1e-5_dp, 'ultrarelativistic limit')

    theta = 0.05_dp
    neq = theta*scaled_bessel_k(2, 1/theta)*exp(-1/theta)/pi**2
    nc = 3*neq
    call pair_equilibrium(theta, nc, np, nm)
    print '(a, 2es10.2)', 'mass action and neutrality defects:', np*nm/neq**2 - 1, (nm - np)/nc - 1
    call require(abs(np*nm/neq**2 - 1) < 1e-7_dp .and. abs((nm - np)/nc - 1) < 1e-10_dp, &
        'law of mass action with ions')

    theta = 2
    call pair_equilibrium(theta, 0.5_dp, np, nm)
    call require(abs((nm - np)/0.5_dp - 1) < 1e-10_dp, 'neutrality, relativistic and degenerate')

    ths = [0.1_dp, 1.0_dp, 10.0_dp]
    do i = 1, 3
        theta = ths(i)
        neq = theta*scaled_bessel_k(2, 1/theta)*exp(-1/theta)/pi**2
        prod = photon_pair_rate(theta, 200)
        ann = neq**2*annihilation_rate(theta, 200)
        print '(a, f6.2, 2es16.8)', 'theta, gamma gamma -> e+e- and annihilation rates:', theta, &
            prod, ann
        call require(abs(prod/ann - 1) < 1e-8_dp, 'detailed balance of pair creation')
    end do
    print '(a)', 'pair equilibrium matches its derived limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_pair_equilibrium
