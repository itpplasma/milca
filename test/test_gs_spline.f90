program test_gs_spline
    !! Grad-Shafranov energy on the spline O-grid (milca_gs_spline) against the
    !! exact solution psi = R^4/8 - 0.05 (R^4 - 4 R^2 Z^2) + 0.2 R^2 with
    !! Delta* psi = R^2 (f = -R) on a D-shaped wall (kappa 1.8, delta 0.4),
    !! Dirichlet data on the wall: L2 error at p = 3 under refinement must
    !! fall at the optimal rate p + 1 = 4. Then the law kernel with
    !! constant p' = A, FF' = B and Newton: psi = -(mu0 A/8)(R^2 - R0^2)^2 -
    !! (B/2) Z^2 + c0 on its own zero level set, same rate. Last, the
    !! normalised law p', FF' ~ (1 - s) in s = 1 - psi/psi_a with the
    !! self-consistent axis flux: Newton with the axis term converges in at
    !! most 6 iterations and R_axis converges at rate >= 3.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use milca_gs_spline, only: gs_spline_kernel_t, gs_spline_dirichlet_solve, GS_SPLINE_OK, &
        gs_spline_law_kernel_t, gs_spline_newton, gs_spline_newton_axis, gs_spline_extremum
    use milca_gs_profile_law, only: gs_profile_law_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp*pi, r0 = 1.7_dp, &
        a_law = 2.0e5_dp, b_law = 0.6_dp, c0 = 0.04_dp
    real(dp) :: err(3), rate
    integer :: k

    do k = 1, 3
        err(k) = solve(2**(k + 1))
        print '(a,i3,a,es10.3)', 'cells per side', 2**(k + 1), '  L2 error', err(k)
    end do
    rate = log(err(2)/err(3))/log(2.0_dp)
    print '(a,f6.2)', 'rate', rate
    if (rate < 3.5_dp) error stop 'spline Grad-Shafranov below the optimal rate'
    if (err(3) > 1.0e-5_dp) error stop 'spline Grad-Shafranov error too large'
    do k = 1, 3
        err(k) = solve_law(2**(k + 1))
        print '(a,i3,a,es10.3)', 'law kernel: cells per side', 2**(k + 1), '  L2 error', err(k)
    end do
    rate = log(err(2)/err(3))/log(2.0_dp)
    print '(a,f6.2)', 'rate', rate
    if (rate < 3.5_dp) error stop 'law kernel below the optimal rate'
    block
        real(dp) :: raxis(3), psia
        integer :: its(3)
        do k = 1, 3
            call solve_axis(2**(k + 2), raxis(k), psia, its(k))
            print '(a,i3,a,f14.10,a,i0)', 'normalised law: cells', 2**(k + 2), '  R_axis', &
                raxis(k), '  Newton iterations ', its(k)
        end do
        if (maxval(its) > 6) error stop 'axis Newton not quadratic'
        rate = log(abs(raxis(1) - raxis(2))/abs(raxis(2) - raxis(3)))/log(2.0_dp)
        print '(a,f6.2)', 'R_axis rate', rate
        if (rate < 3.0_dp) error stop 'R_axis below rate 3'
    end block
    print '(a)', 'spline Grad-Shafranov Solov''ev pass'
contains
    pure function exact(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi
        psi = r**4/8 - 0.05_dp*(r**4 - 4*r**2*z**2) + 0.2_dp*r**2
    end function

    real(dp) function solve(n) result(error)
        integer, intent(in) :: n
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad, check
        type(gs_spline_kernel_t) :: kernel
        real(dp) :: samples(2, 256), t
        real(dp), allocatable :: coef(:, :), jet(:, :)
        integer :: j, st, q

        call layout_ogrid([1.7_dp, 0.0_dp], 0.6_dp, 0.5_dp, layout, st)
        if (st /= 0) error stop 'layout'
        do j = 1, size(samples, 2)
            t = 2*pi*(j - 1)/size(samples, 2)
            samples(:, j) = [1.7_dp + 0.6_dp*cos(t + 0.4_dp*sin(t)), 1.08_dp*sin(t)]
        end do
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'wall curve'
        options%degree = 3
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_quadrature_build(space, 4, quad, st)
        if (st /= 0) error stop 'quadrature'
        kernel%order = 1
        kernel%r = quad%s(1, :)
        kernel%f = -quad%s(1, :)
        allocate (coef(1, space%dof_count))
        call gs_spline_dirichlet_solve(space, quad, kernel, LAYOUT_WALL, exact, coef, st)
        if (st /= GS_SPLINE_OK) error stop 'solve'
        call spline_quadrature_build(space, 6, check, st)
        allocate (jet(1, check%point_count))
        call spline_evaluate(space, check, coef, 0, jet, st)
        if (st /= 0) error stop 'evaluate'
        error = 0
        do q = 1, check%point_count
            error = error + check%weight(q)*(jet(1, q) - exact(check%s(1, q), check%s(2, q)))**2
        end do
        error = sqrt(error)
    end function

    pure function solovev(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi
        psi = -(mu0*a_law/8)*(r**2 - r0**2)**2 - 0.5_dp*b_law*z**2 + c0
    end function

    real(dp) function solve_law(n) result(error)
        integer, intent(in) :: n
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad, check
        type(gs_spline_law_kernel_t) :: kernel
        real(dp) :: samples(2, 256), t, lo, hi, mid
        real(dp), allocatable :: coef(:, :), jet(:, :)
        integer :: j, st, q, it, iterations
        real(dp) :: residual

        ! wall: the zero level of the exact solution along rays from (R0, 0)
        do j = 1, size(samples, 2)
            t = 2*pi*(j - 1)/size(samples, 2)
            lo = 0; hi = 1.5_dp
            do it = 1, 80
                mid = 0.5_dp*(lo + hi)
                if (solovev(r0 + mid*cos(t), mid*sin(t)) > 0) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            samples(:, j) = [r0 + lo*cos(t), lo*sin(t)]
        end do
        call layout_ogrid([r0, 0.0_dp], 0.4_dp, 0.5_dp, layout, st)
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'law wall curve'
        options%degree = 3
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        call spline_quadrature_build(space, 4, quad, st)
        kernel%order = 1
        kernel%r = quad%s(1, :)
        kernel%mu0 = mu0
        allocate (kernel%law%knots(2), kernel%law%pp(4, 1), kernel%law%ff(4, 1))
        kernel%law%knots = [-1.0_dp, 1.0_dp]
        kernel%law%pp = 0; kernel%law%pp(1, 1) = a_law
        kernel%law%ff = 0; kernel%law%ff(1, 1) = b_law
        kernel%law%f_edge = 3.0_dp
        kernel%law%normalized = .false.
        allocate (coef(1, space%dof_count))
        coef = 0
        call gs_spline_newton(space, quad, kernel, LAYOUT_WALL, coef, 1.0e-12_dp, 10, &
                              iterations, residual, st)
        if (st /= GS_SPLINE_OK) error stop 'law Newton'
        call spline_quadrature_build(space, 6, check, st)
        allocate (jet(1, check%point_count))
        call spline_evaluate(space, check, coef, 0, jet, st)
        error = 0
        do q = 1, check%point_count
            error = error + check%weight(q)*(jet(1, q) - solovev(check%s(1, q), check%s(2, q)))**2
        end do
        error = sqrt(error)
    end function

    subroutine solve_axis(n, raxis, psia, iterations)
        integer, intent(in) :: n
        real(dp), intent(out) :: raxis, psia
        integer, intent(out) :: iterations
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(gs_spline_law_kernel_t) :: kernel
        type(gs_spline_kernel_t) :: seed
        real(dp) :: samples(2, 256), t, axis(2), residual, value, check(2)
        real(dp), allocatable :: coef(:, :)
        integer :: j, st

        call layout_ogrid([1.7_dp, 0.0_dp], 0.6_dp, 0.5_dp, layout, st)
        do j = 1, size(samples, 2)
            t = 2*pi*(j - 1)/size(samples, 2)
            samples(:, j) = [1.7_dp + 0.6_dp*cos(t + 0.4_dp*sin(t)), 1.08_dp*sin(t)]
        end do
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        options%degree = 3
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        call spline_quadrature_build(space, 4, quad, st)
        kernel%order = 1
        kernel%mu0 = mu0
        kernel%r = quad%s(1, :)
        allocate (kernel%law%knots(2), kernel%law%pp(4, 1), kernel%law%ff(4, 1))
        kernel%law%knots = [0.0_dp, 1.0_dp]
        kernel%law%pp = 0; kernel%law%pp(1:2, 1) = [-4.0e4_dp, 4.0e4_dp]
        kernel%law%ff = 0; kernel%law%ff(1:2, 1) = [-1.0_dp, 1.0_dp]
        kernel%law%f_edge = 3.0_dp
        kernel%law%axis = -0.3_dp
        kernel%law%normalized = .true.
        seed%order = 1
        seed%r = quad%s(1, :)
        seed%f = spread(1.0_dp, 1, quad%point_count)
        allocate (coef(1, space%dof_count))
        call gs_spline_dirichlet_solve(space, quad, seed, LAYOUT_WALL, zero, coef, st)
        call gs_spline_extremum(space, quad, coef, 1.0_dp, axis, value, st)
        coef = coef*kernel%law%axis/value
        call gs_spline_newton_axis(space, quad, kernel, LAYOUT_WALL, -1.0_dp, coef, &
                                   1.0e-12_dp, 20, iterations, residual, axis, st)
        if (st /= GS_SPLINE_OK) error stop 'axis Newton failed'
        call gs_spline_extremum(space, quad, coef, -1.0_dp, check, psia, st)
        if (abs(psia - kernel%law%axis) > 1.0e-12_dp*abs(psia)) &
            error stop 'law axis flux is not the extremum of psi'
        raxis = axis(1)
    end subroutine

    pure function zero(r, z) result(value)
        real(dp), intent(in) :: r, z
        real(dp) :: value
        value = 0*(r + z)
    end function
end program test_gs_spline
