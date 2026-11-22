program test_ideal_fb_bem_analytic
    !! The vacuum energy of the axisymmetric free-boundary BEM
    !! (theory/derivations/free-boundary-plasma-boundary-bem.md) against
    !! analytic values, no coils: the minimum over sheet currents tau with
    !! oint tau = c of 2 pi tau'V tau/2 is pi c^2/(s'V^-1 s), and for a thin
    !! circular loop of minor radius a about R0 the exterior field energy is
    !! c^2 R0 (ln(8 R0/a) - 2)/2 (O((a/R0)^2 ln) corrections), so
    !! s'V^-1 s -> 2 pi/(R0 (ln(8 R0/a) - 2)).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_fb_bem, only: fb_bem_t, fb_bem_setup, fb_bem_terms, fb_coil_flux_value
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 4.0_dp
    integer :: failures, ia, ib
    real(dp) :: a, rel
    integer, parameter :: nbs(3) = [16, 32, 64]
    real(dp), parameter :: as(2) = [1e-2_dp, 1e-3_dp]

    failures = 0
    do ia = 1, size(as)
        a = as(ia)*r0
        do ib = 1, size(nbs)
            rel = loop_error(nbs(ib), a)
            print '(a,es9.2,a,i4,a,es10.2)', 'a/R0', as(ia), '  boundary nodes', nbs(ib), &
                '  relative error of s''V^-1 s', rel
        end do
        ! thin-loop formula error O((a/R0)^2 ln(8 R0/a))
        if (rel > 10*as(ia)**2*log(8/as(ia)) + 1e-6_dp) then
            print '(a)', 'FAIL thin-loop energy'
            failures = failures + 1
        end if
    end do
    do ib = 1, size(nbs)
        rel = force_error(nbs(ib), 1e-3_dp*r0)
        print '(a,i4,a,es10.2)', 'loop in a coil field, boundary nodes', nbs(ib), &
            '  relative error of dL/dR0', rel
    end do
    if (rel > 1e-3_dp) then
        print '(a)', 'FAIL radial force of a thin loop in a coil field'
        failures = failures + 1
    end if
    if (failures > 0) error stop 'BEM vacuum energy against the thin loop'
    print '(a)', 'PASS BEM vacuum energy of a thin loop'

contains

    real(dp) function loop_error(nb, a) result(rel)
        integer, intent(in) :: nb
        real(dp), intent(in) :: a
        type(fb_bem_t) :: bem
        real(dp), allocatable :: x(:, :), gx(:), gt(:), htt(:, :), htp(:), v(:, :), w(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: energy, gp, th, exact, sv
        integer :: i, st, nt, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        allocate (x(2, nb))
        do i = 1, nb
            th = 2*pi*(i - 1)/nb
            x(:, i) = [r0 + a*cos(th), a*sin(th)]
        end do
        bem%coil%source = reshape([20.0_dp, 0.0_dp], [2, 1])
        bem%coil%coef = [0.0_dp]
        call fb_bem_setup(bem, [(i, i=1, nb)], x, 2, 1, st)
        if (st /= 0) error stop 'bem setup'
        nt = bem%nt
        allocate (gx(2*nb), gt(nt), htt(nt, nt), htp(nt))
        call fb_bem_terms(bem, x(:, bem%node), 1.0_dp, 1.0_dp, energy, gx, gt, gp, htt=htt, &
                          htp=htp)
        v = -htt/(2*pi)
        w = htp
        allocate (ipiv(nt))
        call dgesv(nt, 1, v, nt, ipiv, w, nt, info)
        if (info /= 0) error stop 'dgesv'
        sv = dot_product(htp, w)
        exact = 2*pi/(r0*(log(8*r0/a) - 2))
        rel = abs(sv - exact)/exact
    end function loop_error

    real(dp) function force_error(nb, a) result(rel)
        !! Sum of dL/dR over the chain nodes at the stationary (tau, psi_hat)
        !! against -de/dR0, e = c^2 R0 (ln(8 R0/a) - 2)/2 + 2 pi c psi_c(R0),
        !! one ring coil at (6, 1.5) carrying coef 0.3 in units of G, c = 1.
        integer, intent(in) :: nb
        real(dp), intent(in) :: a
        type(fb_bem_t) :: bem
        real(dp), allocatable :: x(:, :), gx(:), gt(:), htt(:, :), htp(:), m(:, :), b(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: energy, gp, th, psi_c, grad(2), hess(3), dedr, c
        integer :: i, st, nt, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        c = 1
        allocate (x(2, nb))
        do i = 1, nb
            th = 2*pi*(i - 1)/nb
            x(:, i) = [r0 + a*cos(th), a*sin(th)]
        end do
        bem%coil%source = reshape([6.0_dp, 1.5_dp], [2, 1])
        bem%coil%coef = [0.3_dp]
        call fb_bem_setup(bem, [(i, i=1, nb)], x, 2, 1, st)
        if (st /= 0) error stop 'bem setup'
        nt = bem%nt
        allocate (gx(2*nb), gt(nt), htt(nt, nt), htp(nt))
        bem%tau = 0
        bem%psi_hat = 0
        call fb_bem_terms(bem, x(:, bem%node), 1.0_dp, c, energy, gx, gt, gp, htt=htt, htp=htp)
        ! stationary point: [-2 pi V, s; s', 0] [tau; psi_hat] = [2 pi g; c], 2 pi g = -gt
        allocate (m(nt + 1, nt + 1), b(nt + 1), ipiv(nt + 1))
        m = 0
        m(:nt, :nt) = htt
        m(:nt, nt + 1) = htp
        m(nt + 1, :nt) = htp
        b(:nt) = -gt
        b(nt + 1) = c
        call dgesv(nt + 1, 1, m, nt + 1, ipiv, b, nt + 1, info)
        if (info /= 0) error stop 'dgesv'
        bem%tau = b(:nt)
        bem%psi_hat = b(nt + 1)
        call fb_bem_terms(bem, x(:, bem%node), 1.0_dp, c, energy, gx, gt, gp)
        call fb_coil_flux_value(bem%coil, [r0, 0.0_dp], psi_c, grad, hess)
        dedr = c**2*(log(8*r0/a) - 1)/2 + 2*pi*c*grad(1)
        rel = abs(sum(gx(1::2)) + dedr)/abs(dedr)
        print '(a,2es14.6)', '  sum dL/dR, -de/dR0', sum(gx(1::2)), -dedr
    end function force_error

end program test_ideal_fb_bem_analytic
