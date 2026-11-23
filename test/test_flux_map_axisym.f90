program test_flux_map_axisym
    !! Gate F1 of docs/flux-map.md: the axisymmetric Solov'ev equilibrium from
    !! the Cartesian flux map u = (x, y): section -> unit disc, s = x^2 + y^2
    !! = psi_N, with the profiles of the exact equilibrium (as gate E1):
    !!
    !!   psi(s) = psi_b s,  p(s) = source psi_b (1 - s),  tau(s) = Phi'(s)/pi.
    !!
    !! The O-grid is centred 15 % of the minor radius from the magnetic axis
    !! in R and, by default, 10 % of it above the midplane, so that the mesh
    !! has no symmetry; it stays fixed. The start is the polar map about that
    !! centre, u = (x - c)/rho_wall(angle), so no axis information enters:
    !! the axis is the zero of u_h after the solve. The gauge theta* ->
    !! theta* + const is removed by the moment constraint orthogonal to it
    !! (milca_gauge_moments); there are no symmetry rows.
    !! Part 1 (pointwise B_phi): order p at R B_phi = 3. Part 2 (B_phi
    !! projected per cell onto its mean): from the cold start to R B_phi =
    !! 300 (B_phi/B_pol about 150) by the continuation driver
    !! (fortfem_parameter_continuation, R B_phi = 3 * 100^s), which chooses
    !! its stages; the error must not grow with the guide field (R B_phi = 30
    !! is solved from the state at 300).
    !! Arguments (measurements): projection degree (default 0; negative:
    !! pointwise), quadrature points per direction beyond the degree,
    !! spline degree, offset of the mesh centre from the midplane in per cent
    !! of the minor radius (10; 0: the symmetric mesh), 'relaxed' for
    !! relaxed steps of the minimisation ('energy': the default), verbose.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, spline_locate, &
        spline_point_matrix, spline_boundary_nodes
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_boundary_dofs, spline_derham_project
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_spline_ale, only: spline_geometry_dofs
    use fortfem_parameter_continuation, only: continuation_options_t, &
        continuation_report_t, continuation_run
    use milca_flux_map, only: flux_map_t, flux_map_init, flux_map_system, &
        flux_map_solve
    use solovev_label_fields
    implicit none
    integer, parameter :: levels = 3, cells(levels) = [2, 4, 8]
    integer :: degree = 2
    real(dp), parameter :: scales(3) = [1.0_dp, 10.0_dp, 100.0_dp]
    ! interpolation error of B on the axis-centred meshes of gate E1
    real(dp), parameter :: interp(3) = [2.196e-2_dp, 4.367e-3_dp, 9.017e-4_dp]
    integer :: projection = 0, extra_points = 4, level, k, failures, offset = 10
    logical :: verbose = .false., relaxed = .false.
    real(dp) :: centre(2), rmin, rmax
    ! the case of the current level (solve) and the errors of its stages
    type(spline_quadrature_t) :: quad
    type(spline_derham_t) :: derham
    type(flux_map_t) :: fm
    real(dp) :: stage_error(6, size(scales)), flux0
    integer :: stage_cells
    real(dp) :: b_err(levels), s_err(levels), axis_err(levels), pol(levels, size(scales)), &
                e(6, size(scales)), ref_pol
    logical :: ok
    character(16) :: arg

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) projection
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) extra_points
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) degree
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) offset
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        relaxed = trim(arg) == 'relaxed'
    end if
    if (command_argument_count() >= 6) verbose = .true.
    call curve_chart_cylindrical(chart, 1)
    rmin = sqrt(r0**2 - sqrt(psi_b/alpha))
    rmax = sqrt(r0**2 + sqrt(psi_b/alpha))
    centre = [0.5_dp*(rmin + rmax), 0.01_dp*offset*0.5_dp*(rmax - rmin)]
    print '(a,f7.4,a,f7.4,a,2f7.3)', 'mesh centre R ', centre(1), '  magnetic axis R ', r0, &
        '  offset of the mesh centre / minor radius (R, Z) ', &
        (r0 - centre(1))/(0.5_dp*(rmax - rmin)), centre(2)/(0.5_dp*(rmax - rmin))

    print '(a)', 'pointwise B_phi, R B_phi = 3'
    do level = 1, levels
        call solve(cells(level), -1, scales(1:1), e, ok, ref_pol)
        b_err(level) = e(1, 1)
        s_err(level) = e(3, 1)
        axis_err(level) = e(6, 1)
        print '(a,es9.2)', '  B_pol of the projected exact label on this mesh', ref_pol
        call check(ok, 'Newton converges')
        call check(e(2, 1) < 2.5_dp*ref_pol, &
                   'B_pol within a factor 2.5 of the projected exact label')
    end do
    do level = 2, levels
        print '(a,f6.2,a,f6.2,a,f6.2)', 'rate of B', rate(b_err, level), '  of s', &
            rate(s_err, level), '  of the axis', rate(axis_err, level)
        call check(rate(b_err, level) > degree - 0.3_dp, 'B converges at order p')
        call check(axis_err(level) < axis_err(level - 1), 'the axis converges')
    end do

    print '(a,i0,a)', 'B_phi projected onto degree ', projection, &
        ' per cell, cold start at R B_phi = 300 by the continuation driver'
    do level = 1, levels
        call solve(cells(level), projection, scales, e, ok, ref_pol)
        pol(level, :) = e(2, :)
        call check(ok, 'cold start at R B_phi = 300: the continuation converges, no hand-set stage')
        call check(abs(pol(level, 3)/pol(level, 2) - 1) < 0.05_dp, &
                   'B_pol error saturates: R B_phi = 300 within 5 % of 30')
        call check(pol(level, 3) < 1.25_dp*pol(level, 1), &
                   'B_pol error at R B_phi = 300 within 25 % of 3')
    end do
    do level = 2, levels
        print '(a,3f6.2)', 'rate of B_pol at the three guide fields', &
            (log(pol(level - 1, k)/pol(level, k))/log(2.0_dp), k=1, size(scales))
    end do
    call check(rate(pol(:, size(scales)), levels) > degree - 0.3_dp, &
               'B_pol converges at order p at the strongest guide field')
    flux_scale = 1
    if (failures > 0) error stop 'flux map gate F1 failed'

contains

    real(dp) function rate(err, level)
        real(dp), intent(in) :: err(:)
        integer, intent(in) :: level

        rate = log(err(level - 1)/err(level))/log(2.0_dp)
    end function rate

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            print '(2a)', 'PASS ', what
        else
            print '(2a)', 'FAIL ', what
            failures = failures + 1
        end if
    end subroutine check

    subroutine profiles(s, tau, gp, p)
        !! tau = Phi'/pi = -2 F' in the convention of solovev_profiles.
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: f(3)

        call solovev_profiles(s, f, gp, p)
        tau = -2*f
    end subroutine profiles

    real(dp) function wall_radius(theta) result(rho)
        !! Distance from the mesh centre to the wall psi = psi_b along a ray.
        real(dp), intent(in) :: theta
        real(dp) :: lo, hi
        integer :: k

        lo = 0
        hi = 2
        do k = 1, 60
            rho = 0.5_dp*(lo + hi)
            if (psi(centre(1) + rho*cos(theta), centre(2) + rho*sin(theta)) < psi_b) then
                lo = rho
            else
                hi = rho
            end if
        end do
    end function wall_radius

    subroutine solve(n, proj, guide, e, converged, ref_pol)
        !! One guide field (guide of size one): the solve from the polar
        !! start. Several: from the polar start to the last of them by the
        !! continuation driver in s, R B_phi = 3 guide(last)^s, then the
        !! guide fields between from the state at the last. Per guide field e
        !! = relative L2 errors of B and B_pol, the max error of s, the
        !! relative L2 errors of B_phi pointwise and projected, the axis
        !! error. ref_pol: B_pol of the projected exact label psi_N on this
        !! mesh.
        integer, intent(in) :: n, proj
        real(dp), intent(in) :: guide(:)
        real(dp), intent(out) :: e(:, :), ref_pol
        logical, intent(out) :: converged
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(continuation_options_t) :: path
        type(continuation_report_t) :: report
        real(dp), allocatable :: samples(:, :), ctrl(:, :), bt(:), z(:), sl(:), gs(:, :), bb(:)
        integer, allocatable :: wall(:)
        real(dp) :: d(2), theta, energy, x(3), j3(3, 3), w, num(5), den(3)
        real(dp), allocatable :: chi(:), v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        integer :: st, n0, i, p, stage, cell, q

        call boundary_samples(samples)
        call layout_ogrid(centre, 0.5_dp*(rmax - rmin), 0.5_dp, layout, st)
        if (st /= 0) error stop 'layout'
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'wall curve'
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_quadrature_build(space, degree + extra_points, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
        n0 = derham%dimension(0)
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall)
        call flux_map_init(fm, n0, wall)
        fm%projection = proj
        fm%verbose = verbose
        fm%relaxed_steps = relaxed
        stage_cells = n
        ! start: the polar map about the mesh centre, from the control points
        call spline_geometry_dofs(space, ctrl, st)
        if (st /= 0) error stop 'control points'
        do i = 1, n0
            d = ctrl(:, i) - centre
            theta = atan2(d(2), d(1))
            if (fm%is_wall(i)) then
                fm%z(i) = theta
            else
                fm%z(i) = d(1)/wall_radius(theta)
                fm%z(n0 + i) = d(2)/wall_radius(theta)
            end if
        end do
        allocate (bt(quad%point_count), sl(quad%point_count), gs(2, quad%point_count), &
                  bb(quad%point_count))
        flux_scale = guide(1)
        call flux_map_system(fm, derham, quad, profiles, fm%z, energy, label=sl, grad_label=gs, &
                             bphi=bt, bphi_projected=bb)
        flux0 = sum(quad%weight*bt)/(toroidal_flux(1.0_dp)*flux_scale)
        if (size(guide) == 1) then
            call flux_map_solve(fm, derham, quad, profiles, 1e-9_dp*guide(1), 150, st)
            converged = st == 0
            print '(a,i0,a,i0,a,l1)', '  Hessian inertia on the gauge constraint: negative ', &
                fm%negative_eigenvalues, ', zero ', fm%zero_eigenvalues, '  reliable ', &
                fm%inertia_reliable
            call errors(1, fm%steps)
        else
            z = fm%z
            call continuation_run(guide_stage, z, path, report, st)
            fm%z = z
            converged = st == 0
            print '(a,i0,a,i0,a,i0,a,f6.3)', '  continuation: stages ', report%stages, &
                ', failed ', report%failed_stages, ', Newton steps ', report%total_steps, &
                ', reached s = ', report%reached
            call errors(size(guide), report%total_steps)
            do stage = size(guide) - 1, 2, -1
                flux_scale = guide(stage)
                call flux_map_solve(fm, derham, quad, profiles, 1e-9_dp*guide(stage), 60, st)
                converged = converged .and. st == 0
                call errors(stage, fm%steps)
            end do
        end if
        e(:, :size(guide)) = stage_error(:, :size(guide))
        ! reference: B_pol of the projected exact label on this mesh
        allocate (chi(n0))
        call spline_derham_project(derham, 0, label, chi, st, 12)
        if (st /= 0) error stop 'label projection'
        num = 0
        den = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                w = quad%weight(p)/x(1)
                num(1) = num(1) + w*sum((matmul(gn(:, :, q), chi(d0)) - grad_psi(x(1), x(3))/psi_b)**2)
                den(1) = den(1) + w*sum((grad_psi(x(1), x(3))/psi_b)**2)
            end do
        end do
        ref_pol = sqrt(num(1)/den(1))
    end subroutine solve

    subroutine guide_stage(s, x, steps, status)
        !! The stage R B_phi = 3 scales(last)^s of the continuation, from x.
        real(dp), intent(in) :: s
        real(dp), intent(inout) :: x(:)
        integer, intent(out) :: steps, status

        flux_scale = scales(size(scales))**s
        fm%z = x
        call flux_map_solve(fm, derham, quad, profiles, 1e-9_dp*flux_scale, 150, status)
        steps = fm%steps
        if (status /= 0) return
        x = fm%z
        if (s == 0) call errors(1, steps)
    end subroutine guide_stage

    subroutine errors(stage, steps)
        !! Errors of the state fm%z at the guide field flux_scale into
        !! stage_error(:, stage).
        integer, intent(in) :: stage, steps
        real(dp), allocatable :: s(:), gs(:, :), bt(:), bb(:), val(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, x(3), j3(3, 3), b_ex(3), w, bpol(2), num(5), den(3), flux1, &
                    axis(2), uv(2), jac(2, 2), det, xi(2)
        integer :: st, n0, k, block, it, p

        n0 = fm%n0
        allocate (s(quad%point_count), gs(2, quad%point_count), bt(quad%point_count), &
                  bb(quad%point_count))
        call flux_map_system(fm, derham, quad, profiles, fm%z, energy, label=s, &
                             grad_label=gs, bphi=bt, bphi_projected=bb)
        flux1 = sum(quad%weight*bt)/(toroidal_flux(1.0_dp)*flux_scale)
        num = 0
        den = 0
        stage_error(3, stage) = 0
        do p = 1, quad%point_count
            call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
            b_ex = exact_b(x)
            b_ex(2) = flux_scale*b_ex(2)
            w = quad%weight(p)*x(1)
            bpol = psi_b*[-gs(2, p), gs(1, p)]/x(1)
            num(1) = num(1) + w*((bpol(1) - b_ex(1))**2 + (bpol(2) - b_ex(3))**2)
            den(1) = den(1) + w*(b_ex(1)**2 + b_ex(3)**2)
            num(2) = num(2) + w*(bt(p) - b_ex(2))**2
            num(3) = num(3) + w*(bb(p) - b_ex(2))**2
            den(2) = den(2) + w*b_ex(2)**2
            stage_error(3, stage) = max(stage_error(3, stage), &
                                        abs(s(p) - psi(x(1), x(3))/psi_b))
        end do
        stage_error(2, stage) = sqrt(num(1)/den(1))
        stage_error(4, stage) = sqrt(num(2)/den(2))
        stage_error(5, stage) = sqrt(num(3)/den(2))
        ! B with the pointwise B_phi (the projection enters the energy only)
        stage_error(1, stage) = sqrt((num(1) + num(2))/(den(1) + den(2)))
        ! magnetic axis: the zero of u_h, Newton from the mesh centre
        axis = centre
        block = 1
        do it = 1, 30
            call spline_locate(space, block, axis, block, xi, st)
            if (st /= 0) error stop 'axis location'
            do k = 0, 2
                call spline_point_matrix(space, [block], reshape(xi, [2, 1]), k, rowptr, &
                                         colidx, val, st)
                if (any(fm%is_wall(colidx) .and. abs(val) > 1e-14_dp)) &
                    error stop 'axis at the wall'
                if (k == 0) then
                    uv = [dot_product(val, fm%z(colidx)), dot_product(val, fm%z(n0 + colidx))]
                else
                    jac(:, k) = [dot_product(val, fm%z(colidx)), &
                                 dot_product(val, fm%z(n0 + colidx))]
                end if
            end do
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            axis = axis - [jac(2, 2)*uv(1) - jac(1, 2)*uv(2), &
                           jac(1, 1)*uv(2) - jac(2, 1)*uv(1)]/det
            if (norm2(uv) < 1e-13_dp) exit
        end do
        stage_error(6, stage) = norm2(axis - [r0, 0.0_dp])
        print '(a,i2,a,f6.0,a,i3,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,es8.1)', &
            'cells', stage_cells, ' R B_phi', f0*flux_scale, ' steps', steps, ' B', &
            stage_error(1, stage), ' B_pol', stage_error(2, stage), ' s', stage_error(3, stage), &
            ' B_phi', stage_error(4, stage), ' proj', stage_error(5, stage), ' axis', &
            stage_error(6, stage), ' flux change', abs(flux1 - flux0)
    end subroutine errors

end program test_flux_map_axisym
