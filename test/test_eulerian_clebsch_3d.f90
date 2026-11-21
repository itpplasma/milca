program test_eulerian_clebsch_3d
    !! Gates E3 and E3b of docs/eulerian-label-plan.md: gate E1 through the
    !! 3D label lane. chi and g are node values on nt = 2M + 1 toroidal
    !! planes of the section x trigonometric-nodal complex (M = 1), the
    !! energy comes from the 3D FortSym leaf clebsch_energy_3d in physical
    !! (R, phi, Z) components. The axis is a curve X(t) with a point X_k per
    !! plane, and the section core moves with it (ALE), so chi_k(X_k) = 0
    !! and grad chi_k(X_k) = 0 are linear in chi. From labels and axes that
    !! differ per plane, Newton must reach E1 with an axisymmetric solution.
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
    use milca_eulerian_clebsch, only: clebsch_node_gradients, clebsch_3d_system
    use milca_gs_spline, only: gs_spline_extremum
    use solovev_label_fields
    implicit none
    integer, parameter :: degree = 2, levels = 2, nt = 3
    real(dp) :: amplitude = 0.05_dp, ridge = 1e-12_dp, label_amplitude, axis_offset = 2e-3_dp, wall_shift = 0
    integer :: max_steps = 12, extra_points = 8, base_cells = 0, component = 1
    type(block_layout_t) :: base_layout
    real(dp), allocatable :: rv_base(:), zv_base(:)
    integer :: level, cells(levels) = [2, 4], failures, steps, top = levels
    real(dp) :: error(levels), chi_error(levels), rate
    ! interpolation error of B on these meshes (G2a, projection of the exact A)
    real(dp), parameter :: interp(2) = [2.196e-2_dp, 4.367e-3_dp]
    real(dp), allocatable :: history(:)
    character(16) :: arg
    real(dp) :: shift_arg
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
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) axis_offset
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) top
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read (arg, *) wall_shift
    end if
    call curve_chart_cylindrical(chart, 1)
    ! the identities hold for the unshifted wall; the shift is exploration only
    shift_arg = wall_shift
    wall_shift = 0
    call geometry_identity(2, samples_once())
    wall_shift = shift_arg
    do level = 1, top
        call solve(cells(level), error(level), chi_error(level), steps, history)
        print '(a,i3,a,es10.3,a,es10.3,a,i0)', 'cells', cells(level), ' L2 error of B', &
            error(level), ' max |chi - Pi0 psi_N|', chi_error(level), ' Newton steps ', steps
        print '(a,10es9.1)', '  |grad W| history (last 10)', history(max(1, size(history) - 9):)
        call check(steps <= 10, 'Newton converges in at most 10 steps')
        call check(error(level) < 1.1_dp*interp(level), &
                   'B within 10 % of the interpolation error')
    end do
    do level = 2, top
        rate = log(error(level - 1)/error(level))/log(2.0_dp)
        print '(a,f6.2,a,f6.2)', 'rate of B', rate, '  rate of chi', &
            log(chi_error(level - 1)/chi_error(level))/log(2.0_dp)
        call check(rate > degree - 0.3_dp, 'L2 error of B converges at order p')
        call check(log(chi_error(level - 1)/chi_error(level))/log(2.0_dp) > degree, &
                   'labels converge faster than order p')
    end do
    if (failures > 0) error stop 'Clebsch 3D gate E3b failed'

contains

    function samples_once() result(samples)
        real(dp), allocatable :: samples(:, :)

        call boundary_samples(samples)
    end function samples_once

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
        values(1) = pn + label_amplitude*pn*(1 - pn)*(s(1) - r0 + s(2))
        associate (unused => t)
        end associate
    end subroutine perturbed_label

    subroutine coordinate(block, xi, t, values)
        !! R or Z of the current section geometry, in V0 by isoparametry.
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), js(2, 2)

        call section_jet(block, xi, s, js)
        values(1) = s(component)
        associate (unused => t)
        end associate
    end subroutine coordinate

    subroutine zero_profiles(chi, f, g, p)
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: f(3), g(3), p(3)

        f = 0*chi
        g = 0
        p = 0
    end subroutine zero_profiles

    subroutine build_mesh(n, x, samples, quad, aquad, derham)
        !! The four core corners of the O-grid sit at x, so x is the central
        !! vertex of the core block (n even). The wall corners and the wall
        !! curve stay fixed, the ruled blocks blend between them.
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

    subroutine motion_weights(n, samples, axis0, w, rv, zv)
        !! V0 node values w of the section displacement per unit axis shift
        !! (1 on the core block, 0 on the wall), and the node values rv, zv
        !! of R and Z on the mesh centred at axis0. The displacement is the
        !! same for shifts in R and in Z.
        integer, intent(in) :: n
        real(dp), intent(in) :: samples(:, :), axis0(2)
        real(dp), allocatable, intent(out) :: w(:), rv(:), zv(:)
        real(dp), parameter :: shift = 1e-2_dp
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        real(dp), allocatable :: rs(:), zs(:)
        integer :: st

        call build_mesh(n, axis0 + [shift, 0.0_dp], samples, quad, aquad, derham)
        allocate (rs(derham%dimension(0)), zs(derham%dimension(0)), &
                  rv(derham%dimension(0)), zv(derham%dimension(0)))
        component = 1
        call spline_derham_project(derham, 0, coordinate, rs, st, 12)
        call build_mesh(n, axis0 + [0.0_dp, shift], samples, quad, aquad, derham)
        component = 2
        call spline_derham_project(derham, 0, coordinate, zs, st, 12)
        call build_mesh(n, axis0, samples, quad, aquad, derham)
        call spline_derham_project(derham, 0, coordinate, zv, st, 12)
        component = 1
        call spline_derham_project(derham, 0, coordinate, rv, st, 12)
        w = (rs - rv)/shift
        rv_base = rv
        zv_base = zv
        call check(maxval(abs((zs - zv)/shift - w)) < 1e-8_dp, &
                   'the mesh moves rigidly with the axis in R and in Z')
    end subroutine motion_weights

    function wall_offset(k) result(dk)
        !! Helical translation of the whole section on plane k (exploration
        !! of a 3D equilibrium; zero in the gate).
        integer, intent(in) :: k
        real(dp) :: dk(2)

        dk = wall_shift*[cos(2*pi*(k - 1)/nt), sin(2*pi*(k - 1)/nt)]
    end function wall_offset

    function core_nodes(axis0, w, x) result(nodes)
        !! (R, Z) node values per plane of the mesh whose core follows X_k.
        real(dp), intent(in) :: axis0(2), w(:), x(:, :)
        real(dp) :: nodes(2, size(w), size(x, 2))
        integer :: k

        real(dp) :: dk(2)

        do k = 1, size(x, 2)
            dk = wall_offset(k)
            nodes(1, :, k) = rv_base + dk(1) + w*(x(1, k) - axis0(1) - dk(1))
            nodes(2, :, k) = zv_base + dk(2) + w*(x(2, k) - axis0(2) - dk(2))
        end do
    end function core_nodes

    subroutine geometry_identity(n, samples)
        !! g = R and chi = Z on the moving mesh are exact in V0 x trig, so
        !! B = grad g x grad chi = -e_phi and W = int R/2 dR dZ dphi for every
        !! axis curve: the moved Jacobian, the mesh velocity in d/dphi and
        !! the volume weight are exact.
        integer, intent(in) :: n
        real(dp), intent(in) :: samples(:, :)
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian
        real(dp), allocatable :: w(:), rv(:), zv(:), chi(:), g(:), grad(:)
        real(dp) :: axis0(2), plane_axis(2, nt), e_base, e_moved, tk
        real(dp), allocatable :: nodes(:, :, :)
        integer :: k

        axis0 = [r0, 0.0_dp]
        call motion_weights(n, samples, axis0, w, rv, zv)
        call build_mesh(n, axis0, samples, quad, aquad, derham)
        allocate (nodes(2, size(w), nt))
        do k = 1, nt
            tk = (k - 1)/real(nt, dp)
            plane_axis(:, k) = axis0 + 0.1_dp*[cos(2*pi*tk), sin(2*pi*tk)]
        end do
        chi = [(zv, k=1, nt)]
        g = [(rv, k=1, nt)]
        call clebsch_3d_system(derham, aquad, nt, zero_profiles, axis0, chi, g, e_base, &
                               grad, hessian)
        chi = [(zv + w*(plane_axis(2, k) - axis0(2)), k=1, nt)]
        g = [(rv + w*(plane_axis(1, k) - axis0(1)), k=1, nt)]
        call clebsch_3d_system(derham, aquad, nt, zero_profiles, axis0, chi, g, e_moved, &
                               grad, hessian, plane_axis, core_nodes(axis0, w, plane_axis))
        print '(a,es10.2)', 'moving mesh: relative change of int R/2 dV', &
            e_moved/e_base - 1
        call check(abs(e_moved/e_base - 1) < 1e-12_dp, &
                   'energy of an exact field is independent of the axis curve')
        ! the whole section, wall included, translated by 0.1 (cos phi, sin
        ! phi): int R dA per plane is A R_c(phi), whose phi average is A r0
        do k = 1, nt
            nodes(1, :, k) = rv + plane_axis(1, k) - axis0(1)
            nodes(2, :, k) = zv + plane_axis(2, k) - axis0(2)
        end do
        chi = [(nodes(2, :, k), k=1, nt)]
        g = [(nodes(1, :, k), k=1, nt)]
        call clebsch_3d_system(derham, aquad, nt, zero_profiles, axis0, chi, g, e_moved, &
                               grad, hessian, plane_axis, nodes)
        print '(a,es10.2)', 'helical section: relative change of int R/2 dV', &
            e_moved/e_base - 1
        call check(abs(e_moved/e_base - 1) < 1e-12_dp, &
                   'energy of an exact field in a section that moves with phi')
    end subroutine geometry_identity

    subroutine residual(derham, quad, axis0, w, nrow, ngrad, wall, chi, g, x, mu, r, hessian)
        !! Newton residual in (chi, g, X_k, mu_k): stationarity with the
        !! multipliers of chi_k(X_k) = 0, then grad chi_k(X_k) = 0 and
        !! chi_k(X_k) = 0 per plane. X_k is the central vertex of the moving
        !! core, so the rows nrow and ngrad of the base mesh hold for all X.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: axis0(2), w(:), nrow(:), ngrad(:, :), chi(:), g(:), &
                                x(:, :), mu(:)
        integer, intent(in) :: wall(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), allocatable :: grad(:)
        real(dp) :: energy
        integer :: n0, k, a, nf

        n0 = size(nrow)
        nf = 2*n0*nt
        call clebsch_3d_system(derham, quad, nt, solovev_profiles, axis0, chi, g, energy, grad, &
                               hessian, x, core_nodes(axis0, w, x))
        allocate (r(nf + 3*nt))
        r(1:nf) = grad
        do k = 1, nt
            associate (chik => chi((k - 1)*n0 + 1:k*n0))
                r((k - 1)*n0 + 1:k*n0) = r((k - 1)*n0 + 1:k*n0) + mu(k)*nrow
                do a = 1, 2
                    r(nf + 2*(k - 1) + a) = dot_product(ngrad(:, a), chik)
                end do
                r(nf + 2*nt + k) = dot_product(nrow, chik)
            end associate
        end do
        r(wall) = 0
    end subroutine residual

    subroutine solve(n, err, chi_err, steps, history)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, chi_err
        integer, intent(out) :: steps
        real(dp), allocatable, intent(out) :: history(:)
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian, jac_t
        real(dp), allocatable :: samples(:, :), chi2(:), chi(:), g(:), chi_exact(:), grad(:), &
                                 delta(:), nrow(:), ngrad(:, :), r(:), rp(:), rm(:), dgdx(:, :), &
                                 rt(:), chi_avg(:), g_avg(:), w(:), rv(:), zv(:)
        integer, allocatable :: wall2(:), wall(:)
        real(dp) :: energy, axis0(2), x(2, nt), xp(2, nt), mu(nt), alpha, spread_planes, &
                    axis_mean(2), tk
        real(dp), parameter :: hx = 1e-6_dp
        integer :: st, n0, ntot, nf, it, i, k, j

        call boundary_samples(samples)
        axis0 = [r0, 0.0_dp]
        call motion_weights(n, samples, axis0, w, rv, zv)
        call build_mesh(n, axis0, samples, quad, aquad, derham)
        call axis_jets(derham, axis0, nrow, ngrad)
        n0 = derham%dimension(0)
        ntot = n0*nt
        nf = 2*ntot
        allocate (chi2(n0), chi_exact(n0), chi(ntot), g(ntot))
        ! a non-axisymmetric start: labels and axes differ per plane
        do k = 1, nt
            tk = (k - 1)/real(nt, dp)
            label_amplitude = amplitude*(1 + 0.5_dp*cos(2*pi*tk))
            call spline_derham_project(derham, 0, perturbed_label, chi2, st, 12)
            chi((k - 1)*n0 + 1:k*n0) = chi2
            x(:, k) = axis0 + wall_offset(k) + axis_offset*[cos(2*pi*tk), sin(2*pi*tk)]
        end do
        g = 0
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall2)
        wall = [((k - 1)*n0 + wall2, k=1, nt)]
        chi(wall) = 1
        ! initial g: W is quadratic and convex in g at fixed labels
        call clebsch_3d_system(derham, aquad, nt, solovev_profiles, axis0, chi, g, energy, grad, &
                               hessian, x, core_nodes(axis0, w, x))
        call add_g_ridge(hessian, ntot)
        allocate (delta(nf))
        call eu_constrained_solve(hessian, -grad, [(i, i=1, ntot)], 0*chi, delta, st)
        if (st /= 0) error stop 'g solve'
        g = g + delta(ntot + 1:)
        deallocate (delta)
        mu = 0
        allocate (history(0), dgdx(nf + 3*nt, 2*nt))
        steps = 0
        do it = 1, max_steps
            call residual(derham, aquad, axis0, w, nrow, ngrad, wall, chi, g, x, mu, r, hessian)
            history = [history, norm2(r)]
            if (it > 1) then
                ! the X columns are central differences: 1e-8 is their floor
                if (history(it) < 1e-8_dp*history(1)) exit
            end if
            do j = 1, 2*nt
                xp = x
                xp(mod(j - 1, 2) + 1, (j + 1)/2) = x(mod(j - 1, 2) + 1, (j + 1)/2) + hx
                call residual(derham, aquad, axis0, w, nrow, ngrad, wall, chi, g, xp, mu, rp)
                xp(mod(j - 1, 2) + 1, (j + 1)/2) = x(mod(j - 1, 2) + 1, (j + 1)/2) - hx
                call residual(derham, aquad, axis0, w, nrow, ngrad, wall, chi, g, xp, mu, rm)
                dgdx(:, j) = (rp - rm)/(2*hx)
            end do
            call add_g_ridge(hessian, ntot)
            call ale_system(hessian, n0, nrow, ngrad, dgdx, jac_t)
            allocate (delta(nf + 3*nt))
            call eu_constrained_solve(jac_t, -r, wall, 0*chi(wall), delta, st)
            if (st /= 0) error stop 'Newton solve'
            ! non-monotone: from labels that differ per plane, |r| rises once
            ! before the quadratic phase, and a monotone search stalls
            alpha = 1
            do k = 1, 8
                call residual(derham, aquad, axis0, w, nrow, ngrad, wall, &
                              chi + alpha*delta(1:ntot), g + alpha*delta(ntot + 1:nf), &
                              x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt]), &
                              mu + alpha*delta(nf + 2*nt + 1:), rt)
                if (norm2(rt) < 10*maxval(history)) exit
                alpha = alpha/2
            end do
            chi = chi + alpha*delta(1:ntot)
            g = g + alpha*delta(ntot + 1:nf)
            x = x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt])
            mu = mu + alpha*delta(nf + 2*nt + 1:)
            deallocate (delta)
            steps = steps + 1
        end do
        if (wall_shift /= 0) then
            do k = 1, nt
                print '(a,i2,a,2f12.7,a,es12.4)', '  plane', k, ' X - wall offset', &
                    x(:, k) - wall_offset(k), '  chi spread', &
                    maxval(abs(chi((k - 1)*n0 + 1:k*n0) - chi(1:n0)))
            end do
            print '(a,es24.15)', '  energy', energy_at(derham, aquad, axis0, w, chi, g, x)
            err = 0
            chi_err = 0
            return
        end if
        axis_mean = sum(x, 2)/nt
        ! all X_k equal: the moved mesh is the mesh centred at their mean
        call build_mesh(n, axis_mean, samples, quad, aquad, derham)
        call spline_derham_project(derham, 0, label, chi_exact, st, 12)
        chi_avg = sum(reshape(chi, [n0, nt]), 2)/nt
        g_avg = sum(reshape(g, [n0, nt]), 2)/nt
        spread_planes = maxval(abs(reshape(chi, [n0, nt]) - spread(chi_avg, 2, nt)))
        err = field_error(derham, aquad, chi_avg, g_avg, axis_mean)
        chi_err = maxval(abs(chi_avg - chi_exact))
        print '(a,2es11.3,a,es10.2,a,es10.2)', '  axis - exact axis', axis_mean - axis0, &
            '  max plane deviation of X', maxval(abs(x - spread(axis_mean, 2, nt))), &
            ' of chi', spread_planes
        call check(spread_planes < 1e-8_dp, 'labels become axisymmetric')
        call check(maxval(abs(x - spread(axis_mean, 2, nt))) < 1e-8_dp, &
                   'the axis curve becomes a circle')
    end subroutine solve

    real(dp) function energy_at(derham, quad, axis0, w, chi, g, x) result(energy)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: axis0(2), w(:), chi(:), g(:), x(:, :)
        real(dp), allocatable :: grad(:)

        call clebsch_3d_system(derham, quad, nt, solovev_profiles, axis0, chi, g, energy, &
                               grad, plane_axis=x, plane_nodes=core_nodes(axis0, w, x))
    end function energy_at

    subroutine ale_system(hessian, n0, nrow, ngrad, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, g, X_k, mu_k): the (chi, g)
        !! block is the Hessian, the X columns are dgdx (all rows, with the
        !! mesh motion), the mu_k column is N(X) on plane k.
        type(csr_t), intent(in) :: hessian
        integer, intent(in) :: n0
        real(dp), intent(in) :: nrow(:), ngrad(:, :), dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        type(triplet_list) :: t
        integer :: nf, nx, i, j, k, a

        nf = 2*n0*nt
        nx = 2*nt
        call t%init(size(hessian%val) + nx*size(dgdx, 1) + 6*n0*nt + 4*nt)
        call t%add_csr(hessian, 0, 0)
        do j = 1, nx
            do i = 1, nf + 3*nt
                ! J(i, X_j) = dgdx(i, j); in the transpose: row X_j, column i
                if (dgdx(i, j) /= 0) call t%add(nf + j, i, dgdx(i, j))
            end do
            call t%add(nf + j, nf + j, 0.0_dp)
        end do
        do k = 1, nt
            do i = 1, n0
                do a = 1, 2
                    ! J(X_k rows, chi_i on plane k) = d N_i/dX at X
                    call t%add((k - 1)*n0 + i, nf + 2*(k - 1) + a, ngrad(i, a))
                end do
                call t%add(nf + nx + k, (k - 1)*n0 + i, nrow(i))
                call t%add((k - 1)*n0 + i, nf + nx + k, nrow(i))
            end do
            call t%add(nf + nx + k, nf + nx + k, 0.0_dp)
        end do
        call t%to_csr(nf + 3*nt, nf + 3*nt, jac_t)
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

end program test_eulerian_clebsch_3d
