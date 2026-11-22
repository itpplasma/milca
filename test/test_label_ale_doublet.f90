program test_label_ale_doublet
    !! Gate A2a of docs/equilibrium-roadmap.md: two magnetic axes on an ALE
    !! mesh. The exact up-down symmetric Solov'ev doublet (support module
    !! solovev_doublet: psi = psi_O + s chi, Delta* psi = R^2, f = f0) on
    !! FortFEM's separatrix-aligned doublet layout (18 blocks, X a vertex,
    !! cores 1 and 6 about the O-points). milca_label_ale_multi solves the
    !! Grad-Shafranov form of the label lane with chi (chi = 2 on the wall),
    !! both O-points and their multipliers as unknowns; each core follows
    !! its O-point, the X-point vertex follows the discrete saddle between
    !! passes. Start: labels perturbed by 5 %, O-points offset by 10 % of the
    !! lobe half widths in R and Z, in different directions, the X vertex
    !! by half of that (the symmetry is broken).
    !!
    !! Passed when Newton converges in at most 8 steps per pass and 6 passes,
    !! the poloidal field is within 20 % of the projected exact label on the
    !! same final mesh, converges at order p in the lobes (the outer ring to
    !! the wall is pre-asymptotic on these meshes: 1.1 from 4 to 8 cells,
    !! 2.5 from 8 to 16, the same for the projected exact label), both
    !! O-points converge to the exact ones, and the X vertex ends at the
    !! fitted saddle of chi_h and converges to the exact X at order p.
    !!
    !! Arguments: levels (default 3: 2, 4, 8 cells per block side), offset
    !! (default 0.1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_SEPARATRIX, LAYOUT_WALL
    use fortfem_spline_separatrix, only: doublet_curves_t, layout_doublet
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_boundary_dofs, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_project
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_label_ale_multi, only: label_multi_t, label_multi_init, label_multi_solve, &
        label_multi_saddle
    use milca_spline_ale, only: spline_geometry_dofs
    use solovev_doublet, only: solovev_doublet_t
    implicit none
    integer, parameter :: degree = 2, max_levels = 4
    real(dp), parameter :: pi = acos(-1.0_dp), amplitude = 0.05_dp
    integer :: cells(max_levels) = [2, 4, 8, 16], levels = 3, level, failures
    real(dp) :: offset = 0.1_dp, half_width(2)
    real(dp) :: err(max_levels), err_proj(max_levels), o_err(max_levels), x_err(max_levels), &
                lobe_err(max_levels), xv_err(max_levels)
    type(solovev_doublet_t) :: d
    type(block_layout_t) :: layout
    type(spline_space_t), pointer :: space_of_label => null()
    character(16) :: arg

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) offset
    end if
    call d%init()
    print '(a,2f10.6,a,2f10.6,a,2f10.6,a,es11.4)', 'X', d%xx, '  O1', d%xo(:, 1), '  O2', &
        d%xo(:, 2), '  s', d%psi_x - d%psi_o
    call build_layout()
    do level = 1, levels
        call run(cells(level), err(level), err_proj(level), o_err(level), x_err(level), &
                 lobe_err(level), xv_err(level))
    end do
    do level = 2, levels
        print '(a,i3,a,5f7.2)', 'to', cells(level), &
            ' rates of B_pol, in the lobes, projection, O, X vertex', &
            rate(err), rate(lobe_err), rate(err_proj), rate(o_err), rate(xv_err)
        call check(rate(xv_err) > degree - 0.3_dp, 'the X vertex converges at order p')
        call check(rate(lobe_err) > degree - 0.3_dp, 'B_pol in the lobes converges at order p')
        call check(o_err(level) < o_err(level - 1)/2, 'the O-points converge')
    end do
    if (failures > 0) error stop 'gate A2a failed'
    print '(a)', 'PASS gate A2a'

contains

    real(dp) function rate(e)
        real(dp), intent(in) :: e(:)
        rate = log(e(level - 1)/e(level))/log(2.0_dp)
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

    subroutine build_layout()
        type(doublet_curves_t) :: curves
        integer :: st, k
        curves%xpoint = d%xx
        curves%opoint = d%xo
        do k = 1, 2
            call d%lobe_samples(k, 2e-3_dp, curves%lobe(k)%points)
        end do
        call d%wall_samples(curves%wall)
        call layout_doublet(curves, 0.4_dp, layout, st)
        if (st /= 0) error stop 'doublet layout'
        half_width = 0.5_dp*(maxval(curves%lobe(1)%points, 2) - minval(curves%lobe(1)%points, 2))
        print '(a,2f9.5)', 'lobe half widths in R, Z', half_width
    end subroutine build_layout

    subroutine gs_solovev(chi, p, ff)
        !! p(chi) = -s chi (Delta* psi = R^2 = -R^2 p'(psi)), f = f0.
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: p(3), ff(3)
        p = [-(d%psi_x - d%psi_o)*chi, -(d%psi_x - d%psi_o), 0.0_dp]
        ff = 0
    end subroutine gs_solovev

    subroutine run(n, err, err_proj, o_err, x_err, lobe_err, xv_err)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, err_proj, o_err, x_err, lobe_err, xv_err
        type(spline_space_options_t) :: options
        type(spline_space_t), target :: space
        type(spline_derham_t) :: derham
        type(label_multi_t), target :: s
        real(dp), allocatable :: x(:, :), chi(:)
        integer, allocatable :: sep(:)
        real(dp) :: axes(2, 2), dir(2, 2), t0, t1, sc, e2(2)
        integer :: st, k, xdof, i, st2, solve_status
        real(dp) :: xs(2)

        call cpu_time(t0)
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_boundary_dofs(space, LAYOUT_SEPARATRIX, sep)
        call spline_geometry_dofs(space, x, st)
        xdof = minloc([(norm2(x(:, i) - d%xx), i=1, size(x, 2))], 1)
        if (norm2(x(:, xdof) - d%xx) > 1e-10_dp) error stop 'X is not a control point'
        ! offsets in units of the lobe half widths in R and Z
        dir(:, 1) = [1.0_dp, 0.5_dp]
        dir(:, 2) = [-0.5_dp, 1.0_dp]
        do k = 1, 2
            axes(:, k) = d%xo(:, k) + offset*half_width*dir(:, k)
        end do
        sc = d%psi_x - d%psi_o
        ! the X vertex starts off the separatrix crossing as well
        call label_multi_init(s, space, [1, 6], axes, pack(sep, sep /= xdof), xdof, sc, &
                              degree + 3, st, d%xx + offset*half_width*[0.5_dp, 0.25_dp])
        if (st /= 0) then
            print *, 'label_multi_init status', st
            error stop 'init'
        end if
        ! start: the exact label projected on the moved mesh, perturbed
        space_of_label => s%derham%space
        allocate (chi(s%derham%dimension(0)))
        call spline_derham_project(s%derham, 0, perturbed_cb, chi, st, 12)
        if (st /= 0) error stop 'label projection'
        chi(s%wall) = 2
        s%chi = chi
        call label_multi_solve(s, gs_solovev, 1e-10_dp, 1e-3_dp, 8, 8, solve_status)
        call cpu_time(t1)
        print '(a,i3,a,i0,a,i0,a,i0)', 'cells', n, ': status ', solve_status, ', passes ', s%passes, &
            ', Newton steps ', s%newton_total
        print '(a,12es9.1)', '  |r| history', s%newton_history(1:min(12, size(s%newton_history)))
        o_err = max(norm2(s%axes(:, 1) - d%xo(:, 1)), norm2(s%axes(:, 2) - d%xo(:, 2)))
        ! the X vertex against the saddle of the final chi_h
        call label_multi_saddle(s, s%xpoint, xs, st2)
        x_err = norm2(s%xpoint - xs)
        if (st2 /= 0) x_err = huge(1.0_dp)
        print '(a,2es11.3)', '  X vertex - exact X', s%xpoint - d%xx
        xv_err = norm2(s%xpoint - d%xx)
        e2 = bpol_error(s%derham, s%quad, s%chi, sc)
        err = e2(1)
        lobe_err = e2(2)
        ! the projected exact label on the same final mesh
        space_of_label => s%derham%space
        derham = s%derham
        call spline_derham_project(derham, 0, label_cb, chi, st, 12)
        if (st /= 0) error stop 'exact projection'
        e2 = bpol_error(s%derham, s%quad, chi, sc)
        err_proj = e2(1)
        print '(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,f6.2,a)', '  B_pol error', err, &
            ' (projection ', err_proj, '), in the lobes', lobe_err, ', O error', o_err, &
            ', |X - saddle|', x_err, ', ', t1 - t0, ' s'
        call check(solve_status == 0, 'the ALE loop converges')
        call check(s%passes <= 6, 'at most 6 passes')
        call check(s%newton_total <= 8*s%passes, 'at most 8 Newton steps per pass')
        call check(err < 1.2_dp*err_proj, 'B_pol within 20 % of the projected exact label')
        call check(x_err < 0.05_dp*s%cell_size, 'the X vertex is at the saddle')
    end subroutine run

    function bpol_error(derham, quad, chi, sc) result(rel)
        !! Relative L2 errors (weight 2 pi R) of B_pol = s grad chi x grad phi:
        !! over the section, and over the lobes (blocks 1-10) relative to the
        !! whole section.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi(:), sc
        real(dp) :: rel(2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        real(dp) :: gc(2), ge(2), r, w, num(2), den, e2
        integer :: cell, q, p
        num = 0
        den = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                r = quad%s(1, p)
                w = 2*pi*r*quad%weight(p)
                gc = matmul(gn(:, :, q), chi(d0))
                ge = d%label_gradient(quad%s(1, p), quad%s(2, p))
                e2 = w*sc**2*sum((gc - ge)**2)/r**2
                num(1) = num(1) + e2
                if (quad%block(p) <= 10) num(2) = num(2) + e2
                den = den + w*sc**2*sum(ge**2)/r**2
            end do
        end do
        rel = sqrt(num/den)
    end function bpol_error

    subroutine label_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st
        call spline_geometry_at(space_of_label, block, xi, jet, st)
        values(1) = d%label(jet(1, 1), jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

    subroutine perturbed_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6), c
        integer :: st
        call spline_geometry_at(space_of_label, block, xi, jet, st)
        c = d%label(jet(1, 1), jet(2, 1))
        values(1) = c + amplitude*c*(2 - c)*(jet(1, 1) - 1 + jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine perturbed_cb

end program test_label_ale_doublet
