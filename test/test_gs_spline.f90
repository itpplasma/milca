program test_gs_spline
    !! Grad-Shafranov energy on the spline O-grid (milca_gs_spline) against the
    !! exact solution psi = R^4/8 - 0.05 (R^4 - 4 R^2 Z^2) + 0.2 R^2 with
    !! Delta* psi = R^2 (f = -R) on a D-shaped wall (kappa 1.8, delta 0.4),
    !! Dirichlet data on the wall: L2 error at p = 3 under refinement must
    !! fall at the optimal rate p + 1 = 4.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use milca_gs_spline, only: gs_spline_kernel_t, gs_spline_dirichlet_solve, GS_SPLINE_OK
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
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
end program test_gs_spline
