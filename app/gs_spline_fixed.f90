program gs_spline_fixed_main
    !! Fixed-boundary Grad-Shafranov on the spline O-grid (milca_gs_spline):
    !! psi = 0 on the wall curve, normalised profile law in psi_N with the
    !! self-consistent axis flux (Newton at fixed law%axis, then law%axis
    !! from the extremum, repeated). Prints axis, psi_axis, volume, plasma
    !! current and wall time.
    !!
    !! Namelist &gs_spline_fixed: boundary_file (count, then R Z rows,
    !! counter-clockwise from the outboard midplane), profile_file
    !! (normalised law, as gs_phase1), degree, cells (per block side),
    !! core_fraction, mu0, tolerance.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use milca_gs_profile_law, only: read_gs_profile_law, gs_law_sources
    use milca_gs_spline, only: gs_spline_law_kernel_t, gs_spline_newton_axis, &
        gs_spline_extremum, gs_spline_kernel_t, gs_spline_dirichlet_solve, GS_SPLINE_OK
    implicit none
    character(512) :: boundary_file = 'boundary.dat', profile_file = 'profile.law'
    character(512) :: config = 'input.nml'
    integer :: degree = 3, cells = 8
    real(dp) :: core_fraction = 0.5_dp, mu0 = 4.0e-7_dp*acos(-1.0_dp), tolerance = 1.0e-11_dp
    namelist /gs_spline_fixed/ boundary_file, profile_file, degree, cells, core_fraction, &
        mu0, tolerance
    type(block_layout_t) :: layout
    type(spline_space_options_t) :: options
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad
    type(gs_spline_law_kernel_t) :: kernel
    type(gs_spline_kernel_t) :: seed_kernel
    real(dp), allocatable :: samples(:, :), coef(:, :), jet(:, :)
    real(dp) :: axis(2), psi_axis, sigma, residual, volume, current, pp, ff, centre(2), &
        radius
    integer :: u, n, k, st, iterations, q
    integer(8) :: t0, t1, rate

    if (command_argument_count() >= 2) call get_command_argument(2, config)
    open (newunit=u, file=trim(config), status='old', action='read')
    read (u, nml=gs_spline_fixed)
    close (u)
    call system_clock(t0, rate)
    open (newunit=u, file=trim(boundary_file), status='old', action='read')
    read (u, *) n
    allocate (samples(2, n))
    do k = 1, n
        read (u, *) samples(:, k)
    end do
    close (u)
    call read_gs_profile_law(trim(profile_file), kernel%law, st)
    if (st /= 0 .or. .not. kernel%law%normalized) error stop 'invalid normalised profile_file'
    centre = [0.5_dp*(minval(samples(1, :)) + maxval(samples(1, :))), &
              0.5_dp*(minval(samples(2, :)) + maxval(samples(2, :)))]
    radius = 0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :)))
    call layout_ogrid(centre, radius, core_fraction, layout, st)
    if (st /= 0) error stop 'layout'
    call layout_set_curve(layout, LAYOUT_WALL, samples, st)
    if (st /= 0) error stop 'wall curve'
    options%degree = degree
    options%regularity = SPLINE_C0
    allocate (options%knots(2, layout%block_count))
    options%knots(:, :)%element_count = cells
    call spline_space_build(layout, options, space, st)
    if (st /= 0) error stop 'space'
    call spline_quadrature_build(space, degree + 1, quad, st)
    if (st /= 0) error stop 'quadrature'
    kernel%order = 1
    kernel%mu0 = mu0
    kernel%r = quad%s(1, :)
    allocate (coef(1, space%dof_count), jet(1, quad%point_count))
    ! seed: the solution for a uniform source, scaled to the axis flux of the file
    sigma = sign(1.0_dp, kernel%law%axis)
    seed_kernel%order = 1
    seed_kernel%r = quad%s(1, :)
    seed_kernel%f = spread(1.0_dp, 1, quad%point_count)
    call gs_spline_dirichlet_solve(space, quad, seed_kernel, LAYOUT_WALL, zero, coef, st)
    if (st /= GS_SPLINE_OK) error stop 'seed'
    call gs_spline_extremum(space, quad, coef, 1.0_dp, axis, psi_axis, st)
    if (st /= GS_SPLINE_OK) error stop 'seed axis'
    coef = coef*kernel%law%axis/psi_axis
    call gs_spline_newton_axis(space, quad, kernel, LAYOUT_WALL, sigma, coef, tolerance, &
                               40, iterations, residual, axis, st)
    write (output_unit, '(a,i0,a,i0,a,es10.3)') 'Newton status ', st, ' iterations ', &
        iterations, ' residual ', residual
    if (st /= GS_SPLINE_OK) error stop 'Newton failed'
    psi_axis = kernel%law%axis
    call spline_evaluate(space, quad, coef, 0, jet, st)
    volume = 0
    current = 0
    do q = 1, quad%point_count
        volume = volume + 2*acos(-1.0_dp)*quad%s(1, q)*quad%weight(q)
        call gs_law_sources(kernel%law, jet(1, q), pp, ff)
        current = current + (quad%s(1, q)*pp + ff/(mu0*quad%s(1, q)))*quad%weight(q)
    end do
    call system_clock(t1)
    write (output_unit, '(a,2f13.8,a,es18.10,a,f12.7,a,f14.3,a,i0,a,f8.3,a)') 'axis', axis, &
        ' psi_axis', psi_axis, ' volume', volume, ' current', current, ' dofs ', &
        space%dof_count, ' wall ', real(t1 - t0, dp)/rate, ' s'
contains
    pure function zero(r, z) result(value)
        real(dp), intent(in) :: r, z
        real(dp) :: value
        value = 0*(r + z)
    end function
end program gs_spline_fixed_main
