program test_eulerian_clebsch_newton
    !! Gate E1 of docs/eulerian-label-plan.md: the axisymmetric Solov'ev
    !! equilibrium from the label lane, with the labels chi (chi = 1 on the
    !! wall), g in V0 and the magnetic axis X as unknowns and the profiles of
    !! the exact equilibrium for chi = psi_N:
    !!
    !!   G(chi) = psi_b chi,  p(chi) = source psi_b (1 - chi),
    !!   F(chi) = -Phi(chi)/(2 pi),  Phi(chi) = f0 int_{S(chi)} dA/R,
    !!
    !! since int_S B_phi dA = -2 pi F(chi) (the g term integrates to zero over
    !! a level set, Green's theorem on F' grad chi x grad theta). Newton on
    !!
    !!   grad W + mu [N(X); 0] = 0,  grad chi_h(X) = 0,  chi_h(X) = 0,
    !!
    !! W = int |B|^2/2 - int p: X is the label minimum, and chi_h(X) = 0 removes
    !! the line current 2 pi F(chi(X)) delta_X that no quadrature sees. The
    !! mesh follows X (ALE: the core corners of the O-grid move with X, so X is
    !! the central vertex of the core block), and FortFEM's polar rule
    !! integrates the 1/rho terms of the winding form in the four cells at X.
    !! From labels perturbed by 5 %, B and chi are compared with the exact
    !! equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0, &
        spline_space_update_geometry
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, spline_locate, &
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
    integer, parameter :: degree = 2, levels = 3, table_size = 1100
    real(dp), parameter :: table_step = 1e-3_dp
    real(dp) :: table(0:table_size), amplitude = 0.05_dp, ridge = 1e-12_dp
    integer :: max_steps = 12, extra_points = 8, base_cells = 0
    type(block_layout_t) :: base_layout
    integer :: level, cells(levels) = [2, 4, 8], failures, steps
    real(dp) :: error(levels), chi_error(levels), rate
    ! interpolation error of B on these meshes (G2a, projection of the exact A)
    real(dp), parameter :: interp(3) = [2.196e-2_dp, 4.367e-3_dp, 9.017e-4_dp]
    real(dp), allocatable :: history(:)
    character(16) :: arg
    integer :: k

    failures = 0
    ! arguments for measurements: perturbation amplitude, Newton steps, extra
    ! Gauss points per direction (the polar rule needs at least 8 in total)
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) amplitude
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) max_steps
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) extra_points
    end if
    call curve_chart_cylindrical(chart, 1)
    do k = 0, table_size
        table(k) = -toroidal_flux(k*table_step)/(2*pi)
    end do
    do level = 1, levels
        call solve(cells(level), error(level), chi_error(level), steps, history)
        print '(a,i3,a,es10.3,a,es10.3,a,i0)', 'cells', cells(level), ' L2 error of B', &
            error(level), ' max |chi - Pi0 psi_N|', chi_error(level), ' Newton steps ', steps
        print '(a,10es9.1)', '  |grad W| history (last 10)', history(max(1, size(history) - 9):)
        call check(steps <= 10, 'Newton converges in at most 10 steps')
        call check(error(level) < 1.1_dp*interp(level), &
                   'B within 10 % of the interpolation error')
    end do
    do level = 2, levels
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

    pure real(dp) function toroidal_flux(chi) result(phi)
        !! f0 int_{S(chi)} dA/R by rays from the axis, 8-point Gauss in rho.
        real(dp), intent(in) :: chi
        real(dp), parameter :: x8(8) = [-0.9602898564975363_dp, -0.7966664774136267_dp, &
                                         -0.5255324099163290_dp, -0.1834346424956498_dp, &
                                         0.1834346424956498_dp, 0.5255324099163290_dp, &
                                         0.7966664774136267_dp, 0.9602898564975363_dp]
        real(dp), parameter :: w8(8) = [0.1012285362903763_dp, 0.2223810344533745_dp, &
                                         0.3137066458778873_dp, 0.3626837833783620_dp, &
                                         0.3626837833783620_dp, 0.3137066458778873_dp, &
                                         0.2223810344533745_dp, 0.1012285362903763_dp]
        integer, parameter :: rays = 128
        real(dp) :: theta, rho1, rho
        integer :: i, j

        phi = 0
        if (chi <= 0) return
        do i = 1, rays
            theta = 2*pi*(i - 1)/rays
            rho1 = ray_radius(chi*psi_b, theta)
            do j = 1, 8
                rho = rho1*(x8(j) + 1)/2
                phi = phi + w8(j)*rho1/2*rho/(r0 + rho*cos(theta))
            end do
        end do
        phi = f0*phi*2*pi/rays
    end function toroidal_flux

    real(dp) function cubic(i0, x)
        integer, intent(in) :: i0
        real(dp), intent(in) :: x
        real(dp) :: l
        integer :: i, j

        cubic = 0
        do i = 0, 3
            l = 1
            do j = 0, 3
                if (j /= i) l = l*(x/table_step - i0 - j)/(i - j)
            end do
            cubic = cubic + l*table(i0 + i)
        end do
    end function cubic

    subroutine profiles(chi, f, g, p)
        !! F from the table (cubic pieces, F''' exact per piece, linear below
        !! the axis label 0), G and p of the Solov'ev equilibrium.
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: f(3), g(3), p(3)
        real(dp), parameter :: h = 1e-4_dp
        integer :: i0

        g = [psi_b, 0.0_dp, 0.0_dp]
        p = [source*psi_b*(1 - chi), -source*psi_b, 0.0_dp]
        if (chi <= 0) then
            f = [(-11*table(0) + 18*table(1) - 9*table(2) + 2*table(3))/(6*table_step), &
                 0.0_dp, 0.0_dp]
            return
        end if
        i0 = min(max(int(chi/table_step) - 1, 0), table_size - 3)
        f(1) = (cubic(i0, chi + h) - cubic(i0, chi - h))/(2*h)
        f(2) = (cubic(i0, chi + h) - 2*cubic(i0, chi) + cubic(i0, chi - h))/h**2
        f(3) = (table(i0 + 3) - 3*table(i0 + 2) + 3*table(i0 + 1) - table(i0))/table_step**3
    end subroutine profiles

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

    subroutine ale_residual(n, x, samples, chi, g, mu, wall, r)
        !! Newton residual on the mesh centred at x.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2), samples(:, :), chi(:), g(:), mu
        integer, intent(in) :: wall(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian
        real(dp), allocatable :: nrow(:), ngrad(:, :), grad(:)
        real(dp) :: energy

        call build_mesh(n, x, samples, quad, aquad, derham)
        call axis_jets(derham, x, nrow, ngrad)
        call clebsch_axisym_system(derham, aquad, profiles, x, chi, g, energy, grad, hessian)
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
                                 rt(:)
        integer, allocatable :: wall(:)
        real(dp) :: energy, axis(2), mu, alpha
        real(dp), parameter :: hx = 1e-6_dp, unit2(2, 2) = reshape([1.0_dp, 0.0_dp, &
                                                                    0.0_dp, 1.0_dp], [2, 2])
        integer :: st, n0, it, i, k

        call boundary_samples(samples)
        axis = [r0, 0.0_dp]
        call build_mesh(n, axis, samples, quad, aquad, derham)
        n0 = derham%dimension(0)
        allocate (chi(n0), chi_exact(n0), g(n0))
        call spline_derham_project(derham, 0, perturbed_label, chi, st, 12)
        g = 0
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall)
        chi(wall) = 1
        ! initial g: W is quadratic and convex in g at fixed labels
        call clebsch_axisym_system(derham, aquad, profiles, axis, chi, g, energy, grad, hessian)
        call add_g_ridge(hessian, n0)
        allocate (delta(2*n0))
        call eu_constrained_solve(hessian, -grad, [(i, i=1, n0)], 0*chi, delta, st)
        if (st /= 0) error stop 'g solve'
        g = g + delta(n0 + 1:)
        deallocate (delta)
        ! Newton on (chi, g, X, mu): grad W + mu [N(X); 0] = 0, grad chi_h(X) = 0,
        ! chi_h(X) = 0. X is the magnetic axis and the centre of the mesh.
        mu = 0
        allocate (history(0), dgdx(2*n0 + 3, 2))
        steps = 0
        do it = 1, max_steps
            call build_mesh(n, axis, samples, quad, aquad, derham)
            call axis_jets(derham, axis, nrow, ngrad)
            call clebsch_axisym_system(derham, aquad, profiles, axis, chi, g, energy, grad, &
                                       hessian)
            r = [grad + mu*[nrow, 0*nrow], matmul(chi, ngrad), dot_product(nrow, chi)]
            r(wall) = 0
            history = [history, norm2(r)]
            if (it > 1) then
                ! the X columns are central differences: 1e-8 is their floor
                if (history(it) < 1e-8_dp*history(1)) exit
            end if
            ! X columns by central differences, with the mesh moving with X
            do k = 1, 2
                call ale_residual(n, axis + hx*unit2(:, k), samples, chi, g, mu, wall, rp)
                call ale_residual(n, axis - hx*unit2(:, k), samples, chi, g, mu, wall, rm)
                dgdx(:, k) = (rp - rm)/(2*hx)
            end do
            call add_g_ridge(hessian, n0)
            call ale_system(hessian, nrow, ngrad, dgdx, jac_t)
            allocate (delta(2*n0 + 3))
            call eu_constrained_solve(jac_t, -r, wall, 0*chi(wall), delta, st)
            if (st /= 0) error stop 'Newton solve'
            alpha = 1
            do k = 1, 8
                call ale_residual(n, axis + alpha*delta(2*n0 + 1:2*n0 + 2), samples, &
                                  chi + alpha*delta(1:n0), g + alpha*delta(n0 + 1:2*n0), &
                                  mu + alpha*delta(2*n0 + 3), wall, rt)
                if (norm2(rt) < (1 - alpha/4)*history(it)) exit
                alpha = alpha/2
            end do
            chi = chi + alpha*delta(1:n0)
            g = g + alpha*delta(n0 + 1:2*n0)
            axis = axis + alpha*delta(2*n0 + 1:2*n0 + 2)
            mu = mu + alpha*delta(2*n0 + 3)
            deallocate (delta)
            steps = steps + 1
        end do
        call build_mesh(n, axis, samples, quad, aquad, derham)
        call spline_derham_project(derham, 0, label, chi_exact, st, 12)
        err = field_error(derham, aquad, chi, g, axis)
        chi_err = maxval(abs(chi - chi_exact))
        print '(a,2es11.3)', '  axis - exact axis', axis - [r0, 0.0_dp]
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

    real(dp) function field_error(derham, quad, chi, g, axis) result(rel)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi(:), g(:), axis(2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
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
                th = [-(s(2) - axis(2)), s(1) - axis(1)]/((s(1) - axis(1))**2 + &
                                                            (s(2) - axis(2))**2)
                call profiles(c, f, gp, pr)
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

end program test_eulerian_clebsch_newton
