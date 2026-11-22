program test_label_angle_forms
    !! Angle forms of the label lane (docs/eulerian-label-plan.md, section
    !! 4b): E1 (Solov'ev, Newton on chi, g, X, mu) with the guide field
    !! scaled, f0 = R B_phi (any constant is an exact Solov'ev state), and the
    !! winding form either geometric (grad theta_X) or adapted: grad theta_X +
    !! grad psi_h with psi_h single valued in V0, renewed between Newton passes
    !! from the previous state (psi <- psi + g/F'(chi), g <- 0, g re-solved),
    !! so that g stays small and B_phi needs no cancellation between g and the
    !! form. The class of the form (one poloidal turn about X) is topological;
    !! its geometry is free.
    !!
    !! Arguments: f0 (3), mode (geometric | adapted), levels (3), passes (4).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0, &
        spline_space_update_geometry
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, spline_locate, &
        spline_point_matrix, &
        spline_quadrature_point_singular
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_evaluate
    use fortsparse_csr, only: csr_t, csr_matvec
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_axisym_system, clebsch_node_gradients
    use milca_gs_spline, only: gs_spline_extremum
    use solovev_label_fields
    implicit none
    integer, parameter :: levels = 3
    integer :: degree = 2
    real(dp) :: amplitude = 0.05_dp, ridge = 1e-12_dp
    logical :: adapted = .false., use_frame = .false.
    real(dp) :: frame_m(2, 2) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
    integer :: passes = 4
    integer :: max_steps = 12, extra_points = 8, base_cells = 0
    type(block_layout_t) :: base_layout
    integer :: level, cells(levels) = [2, 4, 8], failures, steps
    real(dp) :: error(levels), chi_error(levels), rate
    ! interpolation error of B on these meshes (G2a, projection of the exact A)
    real(dp), parameter :: interp(3) = [2.196e-2_dp, 4.367e-3_dp, 9.017e-4_dp]
    real(dp), allocatable :: history(:)
    character(16) :: arg
    integer :: k, nlev

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) f0
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        adapted = trim(arg) == 'adapted'
    end if
    nlev = levels
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) nlev
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) passes
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) extra_points
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        use_frame = trim(arg) == 'frame'
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) degree
    end if
    print '(a,f8.2,a,l1,a,i0)', 'f0', f0, '  adapted angle', adapted, '  extra points ', extra_points
    call curve_chart_cylindrical(chart, 1)
    do level = 1, nlev
        call solve(cells(level), error(level), chi_error(level), steps, history)
        print '(a,i3,a,es10.3,a,es10.3,a,i0)', 'cells', cells(level), ' L2 error of B', &
            error(level), ' max |chi - Pi0 psi_N|', chi_error(level), ' Newton steps ', steps
        print '(a,10es9.1)', '  |grad W| history (last 10)', history(max(1, size(history) - 9):)
        call check(steps <= 10*passes, 'Newton converges')
        call check(error(level) < 1.1_dp*interp(level), &
                   'B within 10 % of the interpolation error')
    end do
    do level = 2, nlev
        rate = log(error(level - 1)/error(level))/log(2.0_dp)
        print '(a,f6.2,a,f6.2)', 'rate of B', rate, '  rate of chi', &
            log(chi_error(level - 1)/chi_error(level))/log(2.0_dp)
        call check(rate > degree - 0.3_dp, 'L2 error of B converges at order p')
        call check(log(chi_error(level - 1)/chi_error(level))/log(2.0_dp) > degree, &
                   'labels converge faster than order p')
    end do
    if (failures > 0) error stop 'Clebsch Newton gate E1 failed'

contains

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

    subroutine perturbed_label(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), js(2, 2), pn

        call section_jet(block, xi, s, js)
        pn = psi(s(1), s(2))/psi_b
        values(1) = pn + amplitude*pn*(1 - pn)*(s(1) - r0 + s(2))
        associate (unused => t)
        end associate
    end subroutine perturbed_label

    subroutine build_mesh(n, x, samples, quad, aquad, derham)
        !! ALE: the four core corners of the O-grid follow the axis, so X stays
        !! the central vertex of the core block (n even) and the polar and
        !! graded cells are the same for every X. The wall corners and the
        !! wall curve stay fixed, the ruled blocks blend between them.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2), samples(:, :)
        type(spline_quadrature_t), intent(out) :: quad, aquad
        type(spline_derham_t), intent(out) :: derham
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        integer :: st

        if (base_cells /= n) then
            call layout_ogrid([r0, 0.0_dp], 0.5_dp*(maxval(samples(1, :)) - &
                              minval(samples(1, :))), 0.5_dp, base_layout, st)
            if (st /= 0) error stop 'layout'
            call layout_set_curve(base_layout, LAYOUT_WALL, samples, st)
            if (st /= 0) error stop 'wall curve'
            options%degree = degree
            options%regularity = SPLINE_C0
            allocate (options%knots(2, base_layout%block_count))
            options%knots(:, :)%element_count = n
            call spline_space_build(base_layout, options, space, st)
            if (st /= 0) error stop 'space'
            base_cells = n
        end if
        layout = base_layout
        layout%vertices(:, 1:4) = base_layout%vertices(:, 1:4) + spread(x - [r0, 0.0_dp], 2, 4)
        call spline_space_update_geometry(space, layout, st)
        if (st /= 0) error stop 'ALE geometry'
        call spline_quadrature_build(space, degree + extra_points, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
        call spline_quadrature_point_singular(space, quad, x, aquad, st)
        if (st /= 0) error stop 'axis quadrature'
    end subroutine build_mesh

    subroutine system_call(derham, aquad, axis, chi, g, energy, grad, hessian, shift)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: aquad
        real(dp), intent(in) :: axis(2), chi(:), g(:), shift(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: grad(:)
        type(csr_t), intent(out) :: hessian
        if (use_frame) then
            call clebsch_axisym_system(derham, aquad, solovev_profiles, axis, chi, g, energy, &
                                       grad, hessian, shift, frame_m)
        else
            call clebsch_axisym_system(derham, aquad, solovev_profiles, axis, chi, g, energy, &
                                       grad, hessian, shift)
        end if
    end subroutine system_call

    subroutine renew_frame(derham, axis, chi)
        !! frame_m = H^(1/2) of chi_h at X (second derivative rows of the
        !! cell holding X), identity where H is not positive definite.
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: axis(2), chi(:)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        real(dp) :: xi(2), h(3), hm(2, 2), tr, dt, sq
        integer :: blk, st, c
        if (.not. use_frame) return
        call spline_locate(derham%space, 1, axis, blk, xi, st)
        if (st /= 0) error stop 'frame: axis location'
        do c = 3, 5
            call spline_point_matrix(derham%space, [blk], reshape(xi, [2, 1]), c, rp, ci, v, st)
            h(c - 2) = dot_product(v, chi(ci))
        end do
        hm = reshape([h(1), h(2), h(2), h(3)], [2, 2])
        tr = hm(1, 1) + hm(2, 2)
        dt = hm(1, 1)*hm(2, 2) - hm(1, 2)**2
        if (dt > 0 .and. tr > 0) then
            sq = sqrt(dt)
            frame_m = (hm + sq*reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))/sqrt(tr + 2*sq)
        end if
        print '(a,4f9.4)', '  frame', frame_m
    end subroutine renew_frame

    subroutine solve_g(derham, aquad, axis, chi, g, shift)
        !! g minimising W at fixed labels (quadratic and convex in g).
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: aquad
        real(dp), intent(in) :: axis(2), chi(:), shift(:)
        real(dp), intent(inout) :: g(:)
        type(csr_t) :: hessian
        real(dp), allocatable :: grad(:), delta(:)
        real(dp) :: energy
        integer :: n0, i, st
        n0 = size(chi)
        call system_call(derham, aquad, axis, chi, g, energy, grad, &
                                   hessian, shift)
        call add_g_ridge(hessian, n0)
        allocate (delta(2*n0))
        call eu_constrained_solve(hessian, -grad, [(i, i=1, n0)], 0*chi, delta, st)
        if (st /= 0) error stop 'g solve'
        g = g + delta(n0 + 1:)
    end subroutine solve_g

    subroutine ale_residual(n, x, samples, chi, g, mu, wall, r, shift)
        !! Newton residual on the mesh centred at x.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2), samples(:, :), chi(:), g(:), mu, shift(:)
        integer, intent(in) :: wall(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian
        real(dp), allocatable :: nrow(:), ngrad(:, :), grad(:)
        real(dp) :: energy

        call build_mesh(n, x, samples, quad, aquad, derham)
        call axis_jets(derham, x, nrow, ngrad)
        call system_call(derham, aquad, x, chi, g, energy, grad, hessian, &
                                   shift)
        r = [grad + mu*[nrow, 0*nrow], matmul(chi, ngrad), dot_product(nrow, chi)]
        r(wall) = 0
    end subroutine ale_residual

    subroutine solve(n, err, chi_err, steps, history)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, chi_err
        integer, intent(out) :: steps
        real(dp), allocatable, intent(out) :: history(:)
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian, jac_t
        real(dp), allocatable :: samples(:, :), chi(:), g(:), chi_exact(:), grad(:), &
                                 delta(:), nrow(:), ngrad(:, :), r(:), rp(:), rm(:), dgdx(:, :), &
                                 rt(:), shift(:)
        integer, allocatable :: wall(:)
        real(dp) :: energy, axis(2), mu, alpha, r0n, fp(3), gp3(3), pp(3)
        real(dp), parameter :: hx = 1e-6_dp, unit2(2, 2) = reshape([1.0_dp, 0.0_dp, &
                                                                    0.0_dp, 1.0_dp], [2, 2])
        integer :: st, n0, it, i, k, pass, nsc, isc
        logical :: converged

        call boundary_samples(samples)
        axis = [r0, 0.0_dp]
        call build_mesh(n, axis, samples, quad, aquad, derham)
        n0 = derham%dimension(0)
        allocate (chi(n0), chi_exact(n0), g(n0), shift(n0))
        call spline_derham_project(derham, 0, perturbed_label, chi, st, 12)
        g = 0
        shift = 0
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall)
        chi(wall) = 1
        call renew_frame(derham, axis, chi)
        g = 0
        call solve_g(derham, aquad, axis, chi, g, shift)
        mu = 0
        allocate (history(0), dgdx(2*n0 + 3, 2))
        steps = 0
        ! continuation in the guide field: R B_phi = scale f0 from 3 upward
        nsc = max(1, ceiling(log(max(f0/3, 1.0_dp))/log(3.0_dp)) + 1)
        do isc = 1, nsc
            flux_scale = min(1.0_dp, 3**(isc - 1)*3/f0)
            if (nsc == 1) flux_scale = 1
            print '(a,f9.2)', '  guide field R B_phi', flux_scale*f0
            g = 0
            call solve_g(derham, aquad, axis, chi, g, shift)
            do pass = 1, passes
                if (adapted) then
                    ! the angle renewed from the state: theta* = theta - g/F'(chi) in this sign convention
                    do i = 1, n0
                        call solovev_profiles(chi(i), fp, gp3, pp)
                        shift(i) = shift(i) - g(i)/fp(1)
                    end do
                    g = 0
                    call build_mesh(n, axis, samples, quad, aquad, derham)
                    call renew_frame(derham, axis, chi)
                    call solve_g(derham, aquad, axis, chi, g, shift)
                    print '(a,i2,a,es10.2,a,es10.2)', '  pass', pass, '  max |g| after renewal', &
                        maxval(abs(g)), '  max |psi|', maxval(abs(shift))
                end if
                converged = .false.
                r0n = -1
                do it = 1, max_steps
                    call build_mesh(n, axis, samples, quad, aquad, derham)
                    call axis_jets(derham, axis, nrow, ngrad)
                    call system_call(derham, aquad, axis, chi, g, energy, &
                                               grad, hessian, shift)
                    r = [grad + mu*[nrow, 0*nrow], matmul(chi, ngrad), dot_product(nrow, chi)]
                    r(wall) = 0
                    history = [history, norm2(r)]
                    if (r0n < 0) r0n = norm2(r)
                    if (it > 1 .and. norm2(r) < 1e-8_dp*max(r0n, 1.0_dp)) then
                        converged = .true.
                        exit
                    end if
                    do k = 1, 2
                        call ale_residual(n, axis + hx*unit2(:, k), samples, chi, g, mu, wall, rp, shift)
                        call ale_residual(n, axis - hx*unit2(:, k), samples, chi, g, mu, wall, rm, shift)
                        dgdx(:, k) = (rp - rm)/(2*hx)
                    end do
                    call add_g_ridge(hessian, n0)
                    call ale_system(hessian, nrow, ngrad, dgdx, jac_t)
                    allocate (delta(2*n0 + 3))
                    call eu_constrained_solve(jac_t, -r, wall, 0*chi(wall), delta, st)
                    if (st /= 0) error stop 'Newton solve'
                    alpha = 1
                    do k = 1, 12
                        call ale_residual(n, axis + alpha*delta(2*n0 + 1:2*n0 + 2), samples, &
                                          chi + alpha*delta(1:n0), g + alpha*delta(n0 + 1:2*n0), &
                                          mu + alpha*delta(2*n0 + 3), wall, rt, shift)
                        if (norm2(rt) < (1 - alpha/4)*norm2(r)) exit
                        alpha = alpha/2
                    end do
                    print '(a,i3,a,es10.2,a,es9.2,a,es9.2)', '    newton', it, '  |r|', norm2(r), &
                        '  alpha', alpha, '  |dchi|', maxval(abs(delta(1:n0)))
                    chi = chi + alpha*delta(1:n0)
                    g = g + alpha*delta(n0 + 1:2*n0)
                    axis = axis + alpha*delta(2*n0 + 1:2*n0 + 2)
                    mu = mu + alpha*delta(2*n0 + 3)
                    deallocate (delta)
                    steps = steps + 1
                end do
                if (converged .and. (.not. adapted .or. maxval(abs(g)) < 1e-3_dp*maxval(abs(shift)) &
                                     .or. pass == passes)) exit
                if (converged .and. .not. adapted) exit
            end do
        end do
        call build_mesh(n, axis, samples, quad, aquad, derham)
        call spline_derham_project(derham, 0, label, chi_exact, st, 12)
        err = field_error(derham, aquad, chi, g, axis, shift)
        chi_err = maxval(abs(chi - chi_exact))
        print '(a,2es11.3,a,l1)', '  axis - exact axis', axis - [r0, 0.0_dp], '  converged ', &
            converged

    end subroutine solve

    subroutine ale_system(hessian, nrow, ngrad, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, g, X_R, X_Z, mu): the (chi, g)
        !! block is the Hessian, the X columns are dgdx (all rows, with the
        !! mesh motion), the mu column is [N(X); 0; 0].
        type(csr_t), intent(in) :: hessian
        real(dp), intent(in) :: nrow(:), ngrad(:, :), dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        type(triplet_list) :: t
        integer :: n0, nf, i, k

        n0 = size(nrow)
        nf = 2*n0
        call t%init(size(hessian%val) + 4*size(dgdx, 1) + 2*n0 + 4)
        call t%add_csr(hessian, 0, 0)
        do k = 1, 2
            do i = 1, nf + 3
                ! J(i, X_k) = dgdx(i, k); in the transpose: row X_k, column i
                if (dgdx(i, k) /= 0) call t%add(nf + k, i, dgdx(i, k))
            end do
            call t%add(nf + k, nf + k, 0.0_dp)
            do i = 1, n0
                ! J(X_k, chi_i) = dN_i/dX_k (the X rows of the residual)
                call t%add(i, nf + k, ngrad(i, k))
            end do
        end do
        do i = 1, n0
            ! J(chi_i, mu) = N_i(X) and J(mu, chi_i) = N_i(X)
            call t%add(nf + 3, i, nrow(i))
            call t%add(i, nf + 3, nrow(i))
        end do
        call t%add(nf + 3, nf + 3, 0.0_dp)
        call t%to_csr(nf + 3, nf + 3, jac_t)
    end subroutine ale_system

    subroutine axis_jets(derham, axis, nrow, ngrad)
        !! Node function values at X and their X derivatives (central
        !! differences).
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: axis(2)
        real(dp), allocatable, intent(out) :: nrow(:), ngrad(:, :)
        real(dp), parameter :: h = 1e-5_dp
        integer :: k
        real(dp) :: e(2)

        nrow = axis_values(derham, axis)
        allocate (ngrad(size(nrow), 2))
        do k = 1, 2
            e = 0
            e(k) = h
            ngrad(:, k) = (axis_values(derham, axis + e) - axis_values(derham, axis - e))/(2*h)
        end do
    end subroutine axis_jets

    function axis_values(derham, axis) result(row)
        !! Values of the V0 node functions at the axis.
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: axis(2)
        real(dp), allocatable :: row(:), unit(:)
        real(dp) :: xi(2), v(1)
        integer :: block, st, i

        call spline_locate(space, 1, axis, block, xi, st)
        if (st /= 0) error stop 'axis location'
        allocate (row(derham%dimension(0)), unit(derham%dimension(0)))
        unit = 0
        do i = 1, size(row)
            unit(i) = 1
            call spline_derham_evaluate(derham, 0, unit, block, xi, 0.0_dp, v)
            row(i) = v(1)
            unit(i) = 0
        end do
    end function axis_values

    subroutine add_g_ridge(hessian, n0)
        !! g -> g + k(chi) leaves B unchanged: relative ridge on the g block.
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: n0
        real(dp) :: shift
        integer :: i, p

        shift = ridge*maxval(abs(hessian%val))
        do i = n0 + 1, 2*n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine add_g_ridge

    real(dp) function field_error(derham, quad, chi, g, axis, shift) result(rel)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi(:), g(:), axis(2), shift(:)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        real(dp) :: uf(2)
        real(dp) :: s(2), js(2, 2), x(3), j3(3, 3), c, gc(2), gg(2), th(2), f(3), gp(3), &
                    pr(3), bphi, bpol(2), b_ex(3), w, num, den
        integer :: cell, q, p

        num = 0
        den = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                gg = matmul(gn(:, :, q), g(d0))
                if (use_frame) then
                    uf = matmul(frame_m, s - axis)
                    th = matmul(transpose(frame_m), [-uf(2), uf(1)])/(uf(1)**2 + uf(2)**2)
                else
                    th = [-(s(2) - axis(2)), s(1) - axis(1)]/((s(1) - axis(1))**2 + &
                                                                (s(2) - axis(2))**2)
                end if
                th = th + matmul(gn(:, :, q), shift(d0))
                call solovev_profiles(c, f, gp, pr)
                bphi = gg(2)*gc(1) - gg(1)*gc(2) - f(1)*(gc(1)*th(2) - gc(2)*th(1))
                bpol = gp(1)*[-gc(2), gc(1)]/s(1)
                b_ex = exact_b(x)
                w = 2*pi*s(1)*quad%weight(p)
                num = num + w*((bphi - b_ex(2))**2 + (bpol(1) - b_ex(1))**2 + &
                               (bpol(2) - b_ex(3))**2)
                den = den + w*sum(b_ex**2)
            end do
        end do
        rel = sqrt(num/den)
    end function field_error

end program test_label_angle_forms
