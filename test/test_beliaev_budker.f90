program test_beliaev_budker
    ! Relativistic test-particle coefficients against derived properties:
    ! the kernel is symmetric and annihilates u/gamma - u'/gamma' (so a
    ! Juttner distribution is a steady state); detailed balance A = -D u/(gamma
    ! theta) on a Juttner background; the nonrelativistic limit is the Rosenbluth
    ! diffusion D_par = g''/2, D_perp = g'/(2 v) with the isotropic potential
    ! g(v) = 4 pi int f(v') v'^2 (v_> + v_<^2/(3 v_>)) dv'; a fast electron on a
    ! cold background approaches the generated cold-background coefficients.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_beliaev_budker
    use milca_beliaev_budker_leaves, only: bb_kernel, bb_cold
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: u(3), v(3), k(6), m(3, 3), null(3), d(3, 3), a(3), theta, gu, e(3)
    real(dp) :: dpar, dperp, gpp, gp, speed, c(9)

    u = [0.3_dp, -0.7_dp, 1.2_dp]
    v = [-0.4_dp, 0.5_dp, 0.2_dp]
    m = kernel(u, v)
    null = matmul(m, u/sqrt(1 + dot_product(u, u)) - v/sqrt(1 + dot_product(v, v)))
    call require(maxval(abs(null)) < 1e-13_dp*maxval(abs(m)), 'kernel null vector')
    call require(maxval(abs(kernel(v, u) - m)) < 1e-13_dp*maxval(abs(m)), &
        'kernel symmetric in u and u-prime')

    theta = 0.1_dp
    u = [0.3_dp, 0.2_dp, 0.6_dp]
    call juttner_bb_coefficients(u, theta, 32, d, a)
    gu = sqrt(1 + dot_product(u, u))
    print '(a, es10.2)', 'detailed balance defect:', norm2(a + matmul(d, u)/(gu*theta))/norm2(a)
    call require(norm2(a + matmul(d, u)/(gu*theta)) < 1e-6_dp*norm2(a), &
        'detailed balance on a Juttner background')

    theta = 1e-4_dp
    speed = 1.5_dp*sqrt(theta)
    u = [0.0_dp, 0.0_dp, speed]
    call juttner_bb_coefficients(u, theta, 32, d, a)
    call rosenbluth(speed, theta, gp, gpp)
    dpar = d(3, 3)
    dperp = d(1, 1)
    print '(a, 2f10.6)', 'nonrelativistic D_par, D_perp over Rosenbluth:', dpar/(gpp/2), &
        dperp/(gp/(2*speed))
    call require(abs(dpar/(gpp/2) - 1) < 1e-3_dp .and. abs(dperp/(gp/(2*speed)) - 1) < 1e-3_dp, &
        'nonrelativistic Rosenbluth limit')

    u = [0.0_dp, 0.0_dp, 2.0_dp]
    call juttner_bb_coefficients(u, theta, 32, d, a)
    call bb_cold(u(1), u(2), u(3), c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9))
    e = [c(7), c(8), c(9)]
    print '(a, 2f10.6)', 'fast electron, D_perp and A_par over cold background:', &
        d(1, 1)/c(1), a(3)/e(3)
    call require(abs(d(1, 1)/c(1) - 1) < 1e-3_dp .and. abs(a(3)/e(3) - 1) < 1e-3_dp, &
        'cold-background limit')
    print '(a)', 'Beliaev-Budker coefficients match their derived limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function kernel(p, q) result(mat)
        real(dp), intent(in) :: p(3), q(3)
        real(dp) :: mat(3, 3), kk(6)
        call bb_kernel(p(1), p(2), p(3), q(1), q(2), q(3), kk(1), kk(2), kk(3), kk(4), &
            kk(5), kk(6))
        mat = reshape([kk(1), kk(2), kk(4), kk(2), kk(3), kk(5), kk(4), kk(5), kk(6)], [3, 3])
    end function kernel

    ! g'(v) and g''(v) for a normalized Maxwellian of temperature theta
    ! (velocity variance theta per component), from the angle average
    ! <|v - v'|> = v_> + v_<^2/(3 v_>)
    subroutine rosenbluth(vel, th, g1, g2)
        real(dp), intent(in) :: vel, th
        real(dp), intent(out) :: g1, g2
        integer, parameter :: n = 200
        real(dp) :: x(n), w(n), vp, f, vmax
        integer :: i
        call gauss_legendre(n, x, w)
        g1 = 0
        g2 = 0
        do i = 1, n
            ! inner part 0 < v' < v
            vp = vel*(x(i) + 1)/2
            f = exp(-vp**2/(2*th))/(2*pi*th)**1.5_dp
            g1 = g1 + w(i)*vel/2*4*pi*f*vp**2*(1 - vp**2/(3*vel**2))
            g2 = g2 + w(i)*vel/2*4*pi*f*vp**2*(2*vp**2/(3*vel**3))
        end do
        vmax = vel + 12*sqrt(th)
        do i = 1, n
            ! outer part v < v' < vmax
            vp = vel + (vmax - vel)*(x(i) + 1)/2
            f = exp(-vp**2/(2*th))/(2*pi*th)**1.5_dp
            g1 = g1 + w(i)*(vmax - vel)/2*4*pi*f*vp**2*(2*vel/(3*vp))
            g2 = g2 + w(i)*(vmax - vel)/2*4*pi*f*vp*(2.0_dp/3)
        end do
    end subroutine rosenbluth

end program test_beliaev_budker
