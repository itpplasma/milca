program test_cyclotron_emission
    ! Harmonic emission against the Larmor power: summed over harmonics and
    ! directions it is (4 pi/3) p_perp^2 for any parallel momentum; the
    ! thermal emissivity integrated over frequency and directions is the
    ! Juttner average (4 pi/3) <p_perp^2> = (8 pi/3) theta K3(1/theta)/K2(1/theta).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_cyclotron_emission
    use milca_qed_thermal, only: scaled_bessel_k
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nmu = 200
    real(dp) :: xm(nmu), wm(nmu), cases(2, 3), tot, pw, w, theta, ref
    integer :: k, n, i

    call gauss_legendre(nmu, xm, wm)
    cases = reshape([0.6_dp, 0.0_dp, 0.5_dp, 0.4_dp, 0.3_dp, -0.9_dp], [2, 3])
    do k = 1, 3
        tot = 0
        do n = 1, 120
            do i = 1, nmu
                call harmonic_power(n, cases(2, k), cases(1, k), xm(i), pw, w)
                tot = tot + wm(i)*2*pi*pw
            end do
        end do
        print '(a, 2f6.2, es12.3)', 'p_perp, p_par, harmonic sum/Larmor - 1:', cases(:, k), &
            tot/(4*pi/3*cases(1, k)**2) - 1
        call require(abs(tot/(4*pi/3*cases(1, k)**2) - 1) < 1e-10_dp, 'Larmor power')
    end do

    theta = 0.02_dp
    tot = total(theta)
    ref = 8*pi/3*theta*scaled_bessel_k(3, 1/theta)/scaled_bessel_k(2, 1/theta)
    print '(a, es12.3)', 'thermal emission over Juttner Larmor power - 1:', tot/ref - 1
    call require(abs(tot/ref - 1) < 1e-6_dp, 'thermal emissivity integrates to the Larmor power')
    print '(a)', 'cyclotron emission matches the Larmor power'

contains

    real(dp) function total(th)
        real(dp), intent(in) :: th
        integer, parameter :: nw = 64, nh = 30, na = 48
        real(dp) :: xw(nw), ww(nw), xa(na), wa(na), mu, s0, smax, s, u, phi, jac
        integer :: h, i, l
        call gauss_legendre(nw, xw, ww)
        call gauss_legendre(na, xa, wa)
        total = 0
        ! each harmonic separately: w = h/s over its resonant range of s; the
        ! emissivity is even in mu, and mu = cos(phi) on 0 < phi < pi/2 is smooth
        do h = 1, nh
            do l = 1, na
                phi = pi/4*(xa(l) + 1)
                mu = cos(phi)
                jac = 2*pi/4*sin(phi)
                s0 = sin(phi)
                smax = (1 + 60*th)*(1 + mu)
                do i = 1, nw
                    u = (xw(i) + 1)/2
                    s = s0 + (smax - s0)*u**2
                    total = total + wa(l)*jac*2*pi*ww(i)*(smax - s0)*u*h/s**2* &
                        harmonic_emissivity(th, h, h/s, mu)
                end do
            end do
        end do
    end function total

    ! the emissivity of harmonic h alone
    real(dp) function harmonic_emissivity(th, h, w, mu)
        real(dp), intent(in) :: th, w, mu
        integer, intent(in) :: h
        harmonic_emissivity = thermal_emissivity(th, w, mu, h, 48, nmin=h)
    end function harmonic_emissivity

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_cyclotron_emission
