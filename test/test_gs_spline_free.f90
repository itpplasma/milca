program test_gs_spline_free
    !! Free-boundary coupling of the spline Grad-Shafranov energy
    !! (gs_spline_vacuum_setup, gs_spline_free_linear) on a circular wall:
    !! - coils only: a filament outside the wall, psi must converge to its
    !!   flux psi_c = I G(x, x_c) inside at the optimal rate;
    !! - coils and plasma: the plasma flux w = S(t) G(x, x0) (C5 smoothstep of
    !!   the Green function between radii 0.15 and 0.6 about x0, the field of
    !!   FortFEM's coupling test) with its source f = -div(grad w/R); then
    !!   psi = w + psi_c, L2 error at the optimal rate p + 1 = 4;
    !! - the nonlinear free-boundary Newton (gs_spline_newton_free) with
    !!   constant sources p' = A, FF' = B must reproduce the linear solve with
    !!   f = mu0 R A + B/R to round-off.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_green_hessian, gs_green_gradient
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use milca_gs_spline, only: gs_spline_kernel_t, gs_spline_vacuum_t, &
        gs_spline_vacuum_setup, gs_spline_free_linear, GS_SPLINE_OK, &
        gs_spline_law_kernel_t, gs_spline_newton_free
    use fortfem_spline_eval, only: spline_locate
    implicit none
    real(dp), parameter :: x0(2) = [3.1_dp, 0.05_dp], a = 0.15_dp, b = 0.6_dp, &
        coil(2) = [4.3_dp, 0.9_dp], current = 2.0_dp
    real(dp) :: err(3), vac_err, rate
    integer :: k

    vac_err = solve(8, .false.)
    rate = log(vac_err/solve(16, .false.))/log(2.0_dp)
    print '(a,es10.3,a,f6.2)', 'coils only: L2 error at 8 cells', vac_err, ', rate', rate
    if (rate < 3.5_dp) error stop 'coil flux not reproduced at the optimal rate'
    do k = 1, 3
        err(k) = solve(2**(k + 1), .true.)
        print '(a,i3,a,es10.3)', 'coils and plasma: cells', 2**(k + 1), '  L2 error', err(k)
    end do
    rate = log(err(2)/err(3))/log(2.0_dp)
    print '(a,f6.2)', 'rate', rate
    if (rate < 3.5_dp) error stop 'free-boundary coupling below the optimal rate'
    call newton_against_linear()
    print '(a)', 'spline free-boundary coupling pass'
contains
    real(dp) function solve(n, plasma) result(error)
        integer, intent(in) :: n
        logical, intent(in) :: plasma
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad, check
        type(gs_spline_kernel_t) :: kernel
        type(gs_spline_vacuum_t) :: vac
        real(dp), allocatable :: coef(:, :), jet(:, :), psi_c(:), dn(:)
        real(dp) :: w, exact, g, gg(4)
        integer :: q, st

        call layout_ogrid([3.1_dp, 0.0_dp], 0.8_dp, 0.5_dp, layout, st)
        options%degree = 3
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        call spline_quadrature_build(space, 4, quad, st)
        kernel%order = 1
        kernel%r = quad%s(1, :)
        allocate (kernel%f(quad%point_count))
        kernel%f = 0
        if (plasma) then
            do q = 1, quad%point_count
                kernel%f(q) = source(quad%s(:, q))
            end do
        end if
        call gs_spline_vacuum_setup(space, LAYOUT_WALL, 2, 1.0_dp, vac, st)
        if (st /= GS_SPLINE_OK) error stop 'vacuum setup'
        allocate (psi_c(size(vac%nodes, 2)))
        do q = 1, size(psi_c)
            psi_c(q) = current*gs_green(vac%nodes(1, q), vac%nodes(2, q), coil(1), coil(2))
        end do
        allocate (dn(vac%side%point_count))
        do q = 1, size(dn)
            call gs_green_gradient(vac%side%s(1, q), vac%side%s(2, q), coil(1), coil(2), g, gg)
            dn(q) = current*dot_product(gg(1:2), vac%normal(:, q))
        end do
        allocate (coef(1, space%dof_count))
        call gs_spline_free_linear(space, quad, kernel, vac, psi_c, dn, coef, st)
        if (st /= GS_SPLINE_OK) error stop 'free solve'
        call spline_quadrature_build(space, 6, check, st)
        allocate (jet(1, check%point_count))
        call spline_evaluate(space, check, coef, 0, jet, st)
        error = 0
        do q = 1, check%point_count
            exact = current*gs_green(check%s(1, q), check%s(2, q), coil(1), coil(2))
            if (plasma) then
                call field(check%s(:, q), w)
                exact = exact + w
            end if
            error = error + check%weight(q)*(jet(1, q) - exact)**2
        end do
        error = sqrt(error)
    end function

    subroutine newton_against_linear()
        real(dp), parameter :: mu0 = 4.0e-7_dp*acos(-1.0_dp), a_law = 3.0e5_dp, &
            b_law = 0.4_dp
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(gs_spline_kernel_t) :: plain
        type(gs_spline_law_kernel_t) :: law
        type(gs_spline_vacuum_t) :: vac1, vac0
        real(dp), allocatable :: c1(:, :), c2(:, :), psi_c(:), dn(:)
        real(dp) :: g, gg(4), axis(2), residual, txi(2, 1)
        integer :: q, st, iterations, tblock(1)
        call layout_ogrid([3.1_dp, 0.0_dp], 0.8_dp, 0.5_dp, layout, st)
        options%degree = 3
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = 6
        call spline_space_build(layout, options, space, st)
        call spline_quadrature_build(space, 4, quad, st)
        call gs_spline_vacuum_setup(space, LAYOUT_WALL, 2, 1.0_dp, vac1, st)
        call gs_spline_vacuum_setup(space, LAYOUT_WALL, 2, mu0, vac0, st)
        allocate (psi_c(size(vac1%nodes, 2)), dn(vac1%side%point_count))
        do q = 1, size(psi_c)
            psi_c(q) = current*gs_green(vac1%nodes(1, q), vac1%nodes(2, q), coil(1), coil(2))
        end do
        do q = 1, size(dn)
            call gs_green_gradient(vac1%side%s(1, q), vac1%side%s(2, q), coil(1), coil(2), g, gg)
            dn(q) = current*dot_product(gg(1:2), vac1%normal(:, q))
        end do
        ! a weak coil (1/40 of the coupling test) so that the plasma has an
        ! interior flux maximum, which the free-boundary Newton needs
        psi_c = psi_c/40
        dn = dn/40
        ! linear reference in mu0 = 1 units: the coil terms scale with mu0
        plain%order = 1
        plain%r = quad%s(1, :)
        plain%f = mu0*quad%s(1, :)*a_law + b_law/quad%s(1, :)
        allocate (c1(1, space%dof_count), c2(1, space%dof_count))
        call gs_spline_free_linear(space, quad, plain, vac1, psi_c, dn, c1, st)
        if (st /= GS_SPLINE_OK) error stop 'linear reference'
        law%order = 1
        law%mu0 = mu0
        law%r = quad%s(1, :)
        law%plasma = spread(.true., 1, quad%point_count)
        allocate (law%law%knots(2), law%law%pp(4, 1), law%law%ff(4, 1))
        law%law%knots = [-1.0e3_dp, 1.0e3_dp]
        law%law%pp = 0; law%law%pp(1, 1) = a_law
        law%law%ff = 0; law%law%ff(1, 1) = b_law
        law%law%f_edge = 3.0_dp
        ! the whole wall interior is plasma here: an edge pressure keeps
        ! p(psi - psi_b) >= 0 where psi < psi_b (the gradient is unchanged)
        law%law%p_edge = 1.0e6_dp
        law%law%normalized = .false.
        call spline_locate(space, 1, [3.4_dp, 0.0_dp], tblock(1), txi(:, 1), st)
        ! a start with an extremum: half the linear solution
        c2 = 0.5_dp*c1
        call gs_spline_newton_free(space, quad, law, vac0, psi_c, dn, 1.0_dp, tblock, txi, &
                                   c2, 1.0e-13_dp, 10, iterations, residual, axis, st)
        if (st /= GS_SPLINE_OK) error stop 'free-boundary Newton'
        print '(a,es10.3,a,i0)', 'Newton against linear: max coefficient difference', &
            maxval(abs(c2 - c1))/maxval(abs(c1)), ', iterations ', iterations
        if (maxval(abs(c2 - c1)) > 1.0e-9_dp*maxval(abs(c1))) &
            error stop 'free-boundary Newton differs from the linear solve'
    end subroutine

    subroutine field_all(x, w, grad, hess)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: w, grad(2), hess(2, 2)
        real(dp) :: g, g4(4), h4(4, 4), d(2), r, t, st, dst, d2st, dt(2)
        integer :: i, j
        d = x - x0
        r = norm2(d)
        t = (r - a)/(b - a)
        if (t <= 0.0_dp) then
            w = 0; grad = 0; hess = 0
            return
        end if
        call gs_green_hessian(x(1), x(2), x0(1), x0(2), g, g4, h4)
        if (t >= 1.0_dp) then
            w = g; grad = g4(1:2); hess = h4(1:2, 1:2)
            return
        end if
        st = t**6*(462.0_dp + t*(-1980.0_dp + t*(3465.0_dp + t*(-3080.0_dp + &
            t*(1386.0_dp - 252.0_dp*t)))))
        dst = 2772.0_dp*(t*(1.0_dp - t))**5
        d2st = 13860.0_dp*(t*(1.0_dp - t))**4*(1.0_dp - 2.0_dp*t)
        dt = d/(r*(b - a))
        w = st*g
        grad = dst*dt*g + st*g4(1:2)
        do j = 1, 2
            do i = 1, 2
                hess(i, j) = d2st*dt(i)*dt(j)*g + dst*g*(merge(1.0_dp, 0.0_dp, i == j) &
                    - d(i)*d(j)/r**2)/(r*(b - a)) + dst*(dt(i)*g4(j) + dt(j)*g4(i)) &
                    + st*h4(i, j)
            end do
        end do
    end subroutine

    subroutine field(x, w)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: w
        real(dp) :: grad(2), hess(2, 2)
        call field_all(x, w, grad, hess)
    end subroutine

    real(dp) function source(x) result(f)
        real(dp), intent(in) :: x(2)
        real(dp) :: w, grad(2), hess(2, 2)
        call field_all(x, w, grad, hess)
        f = -(hess(1, 1) + hess(2, 2))/x(1) + grad(1)/x(1)**2
    end function
end program test_gs_spline_free
