program test_qed_ee_power
    ! Electron-electron bremsstrahlung power of a Juttner plasma: the
    ! energy-weighted invariant reduction (K2 in place of K1) against a
    ! direct quadrature over both momenta for a test kernel; the power from
    ! the eight-diagram amplitudes against the fit of Svensson (Astrophys. J.
    ! 258, 335, 1982), n_e^2 r_e^2 alpha m c^3 (20/(9 sqrt(pi))) (44 - 3 pi^2)
    ! theta^(3/2) (1 + 1.1 theta + theta^2 - 1.25 theta^(5/2)), within the
    ! few percent of that fit.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_thermal, only: pair_energy_rate
    use milca_qed_bremsstrahlung, only: electron_energy_radiated
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: reduced, direct, theta, mine, fit

    reduced = pair_energy_rate(0.3_dp, 48, test_kernel)
    direct = direct_energy_rate(0.3_dp, 40)
    print '(a, 2es16.8)', 'energy-weighted rate, reduced and direct:', reduced, direct
    call require(abs(reduced/direct - 1) < 1e-8_dp, 'K2 reduction equals direct quadrature')

    ! the nonrelativistic coefficient (20/(9 sqrt(pi)))(44 - 3 pi^2) is
    ! analytic: agreement to the first relativistic correction at low theta
    theta = 0.003_dp
    mine = 0.5_dp*pair_energy_rate(theta, 24, radiated)
    print '(a, f8.4)', 'electron-electron power over Svensson at theta = 0.003:', &
        mine/svensson(theta)
    call require(abs(mine/svensson(theta) - 1) < 5e-3_dp, &
        'nonrelativistic electron-electron power coefficient')
    ! the relativistic correction of the fit is approximate
    theta = 0.1_dp
    ! the 1/2 counts each pair of identical electrons once
    mine = 0.5_dp*pair_energy_rate(theta, 24, radiated)
    fit = svensson(theta)
    print '(a, 2es14.6, f8.4)', 'electron-electron power and Svensson fit at theta = 0.1:', &
        mine, fit, mine/fit
    call require(abs(mine/fit - 1) < 0.05_dp, 'electron-electron power within the fit accuracy')
    print '(a)', 'electron-electron bremsstrahlung power matches its oracles'

contains

    real(dp) function svensson(theta)
        real(dp), intent(in) :: theta
        svensson = 20/(9*sqrt(pi))*(44 - 3*pi**2)*theta**1.5_dp* &
            (1 + 1.1_dp*theta + theta**2 - 1.25_dp*theta**2.5_dp)
    end function svensson

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    real(dp) function test_kernel(s)
        real(dp), intent(in) :: s
        test_kernel = s - 4
    end function test_kernel

    real(dp) function radiated(s)
        real(dp), intent(in) :: s
        radiated = electron_energy_radiated(s, 12)
    end function radiated

    ! <v_Moller gamma_cm w(s)>: kinetic energies in t^2 theta, the relative
    ! angle by Gauss-Legendre; gamma_cm = (E1 + E2)/sqrt(s)
    real(dp) function direct_energy_rate(theta, n) result(rate)
        real(dp), intent(in) :: theta
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), tmax, t1, t2, g1, g2, p1, p2, w1, w2, dot, s, vm, norm
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
                    dot = g1*g2 - p1*p2*x(l)
                    s = 2 + 2*dot
                    vm = sqrt(dot**2 - 1)/(g1*g2)
                    rate = rate + w1*w2*w(l)*test_kernel(s)*vm*(g1 + g2)/sqrt(s)
                end do
            end do
        end do
        rate = rate/norm
    end function direct_energy_rate

end program test_qed_ee_power
