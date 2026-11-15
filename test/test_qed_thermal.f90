program test_qed_thermal
    ! Thermal pair annihilation: the scaled Bessel functions against their
    ! large-argument expansion; the one-dimensional invariant rate against a
    ! direct quadrature over both Juttner momenta and the angle between them;
    ! the nonrelativistic limit pi r_e^2 c of Dirac's sigma v, approached
    ! as 1 - 9 theta^2 (no first-order term).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_thermal
    use milca_qed_leaves, only: qed_dirac_annihilation_total
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: x, k, ref, rate, direct, d(2)

    ! e^x K_nu(x) = sqrt(pi/(2x)) (1 + (mu-1)/(8x) + (mu-1)(mu-9)/(2 (8x)^2)
    !     + (mu-1)(mu-9)(mu-25)/(6 (8x)^3) + ...), mu = 4 nu^2
    x = 400
    k = scaled_bessel_k(2, x)
    ref = sqrt(pi/(2*x))*(1 + 15/(8*x) + 15*7/(2*(8*x)**2) + 15*7*(-9)/(6*(8*x)**3))
    call require(abs(k/ref - 1) < 1e-9_dp, 'scaled K2 at large argument')
    rate = annihilation_rate(0.3_dp, 48)
    direct = direct_rate(0.3_dp, 40)
    print '(a, 2es16.8)', 'annihilation rate at theta = 0.3, invariant and direct:', rate, direct
    call require(abs(rate/direct - 1) < 1e-8_dp, 'invariant reduction equals direct quadrature')
    d(1) = annihilation_rate(1e-4_dp, 48)/pi - 1
    d(2) = annihilation_rate(4e-4_dp, 48)/pi - 1
    print '(a, 2es10.2)', 'defect against pi r_e^2 c at theta = 1e-4, 4e-4:', d
    ! sigma v/(pi r_e^2 c) = 1 - (3/5) eps^2 + ..., eps = gamma_rel - 1, and
    ! <eps^2> = 15 theta^2 for two Maxwellians: the defect is -9 theta^2
    call require(abs(d(1)/(-9e-8_dp) - 1) < 0.05_dp .and. abs(d(2)/(-9*16e-8_dp) - 1) < 0.05_dp, &
        'nonrelativistic limit with defect -9 theta^2')
    print '(a)', 'thermal annihilation rate matches its oracles'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    ! <sigma v_Moller> over two Juttner distributions: kinetic energies in
    ! t^2 theta, the relative angle by Gauss-Legendre
    real(dp) function direct_rate(theta, n) result(rate)
        real(dp), intent(in) :: theta
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), tmax, t1, t2, g1, g2, p1, p2, w1, w2, c, dot, sigma, vm, norm
        integer :: i, j, l
        call gauss_legendre(n, x, w)
        tmax = sqrt(60.0_dp)
        rate = 0
        norm = 0
        do i = 1, n
            t1 = tmax*(x(i) + 1)/2
            g1 = 1 + theta*t1**2
            p1 = sqrt(g1**2 - 1)
            w1 = w(i)*tmax/2*2*t1*theta*p1*g1*exp(-t1**2)
            do j = 1, n
                t2 = tmax*(x(j) + 1)/2
                g2 = 1 + theta*t2**2
                p2 = sqrt(g2**2 - 1)
                w2 = w(j)*tmax/2*2*t2*theta*p2*g2*exp(-t2**2)
                norm = norm + w1*w2*2
                do l = 1, n
                    c = x(l)
                    dot = g1*g2 - p1*p2*c
                    vm = sqrt(dot**2 - 1)/(g1*g2)
                    call qed_dirac_annihilation_total(dot, sigma)
                    rate = rate + w1*w2*w(l)*sigma*vm
                end do
            end do
        end do
        rate = rate/norm
    end function direct_rate

end program test_qed_thermal
