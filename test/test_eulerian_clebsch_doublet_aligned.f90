program test_eulerian_clebsch_doublet_aligned
    !! Gate E2b of docs/eulerian-label-plan.md: gate E2a on a mesh aligned with
    !! the separatrix. The exact Solov'ev doublet (test/support/solovev_doublet)
    !! is meshed by FortFEM's layout_doublet (18 blocks, X a vertex of valence
    !! 8, the separatrix lobes as block sides), with knots graded algebraically
    !! towards the separatrix and X (power grading). The label chi lies in the
    !! continuous V0 of that mesh, g in V0 cut along the separatrix, so g may
    !! jump there by the region-dependent function of the exact solution. With
    !! the label held fixed, the toroidal field f0/R is fitted in L2 by
    !!
    !!   B_phi = (grad g x grad chi)_phi + sum_j F_j'(chi) (grad chi x h_j)_phi
    !!
    !! with the exact profiles per region (lobe j: only h_j, F_j = -Phi_lobe/
    !! (2 pi); outside: both, each with half of the toroidal flux); a point's
    !! region is the region of its block, not a value of chi. Three labels
    !! separate the error of a misplaced separatrix from that of g: the exact
    !! chi, its projection chi_h, and chi_h with its separatrix nodes set to 1
    !! (exactly 1 on the mesh separatrix). The L2 error of B_phi is measured
    !! away from the separatrix (lobes chi < 0.3, outside 1.7 < chi < 1.95),
    !! in a band |chi - 1| < 0.2 outside a disc of radius 0.04 about X, and in
    !! that disc. A spurious current sheet on the separatrix would be a jump of
    !! R B_phi (= f0 exactly), which no L2 measure of fixed width sees as an
    !! order loss; the layer measure is the RMS of R B_phi/f0 - 1 in the layer
    !! 0 < |chi - 1| < 1/(2n) outside the disc, relative to the layer itself.
    !! It stays O(1) under refinement with a sheet and falls without one. (The
    !! one-sided traces themselves are no measure: each Clebsch term is
    !! log-singular on the separatrix, only their sum is finite.) The mesh
    !! separatrix is the degree-p isoparametric image of the traced one; its
    !! misplacement max |chi - 1| on it is reported, and with the exact label
    !! it puts points of the band on the wrong side of the singular profiles.
    !!
    !! Grading: FortFEM's power grading (breakpoints (i/n)^beta along a whole
    !! block chain) leaves the cells away from the separatrix of size about
    !! beta/n, so the smooth regions lose their order before the singular ones
    !! gain it. The graded ends are therefore restricted to a layer: a
    !! fraction theta of the chain's cells at x_i = tau (i/m)^beta, m = theta n,
    !! the rest uniform, with tau = theta/(beta (1 - theta) + theta) so that the
    !! last layer cell and the uniform cells have the same size. theta = 0 is
    !! the plain power grading.
    !!
    !! Arguments: levels (default 3: 4, 8, 16 cells per block side), grading
    !! power (default 4), layer fraction theta (default 0.5), degree
    !! (default 2), radius of the disc about X (default 0.04).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_SEPARATRIX
    use fortfem_spline_separatrix, only: doublet_curves_t, layout_doublet
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_grade_knots, knot_spec_t, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate, &
        spline_side_quadrature_t, spline_interface_quadrature
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_evaluate, spline_derham_boundary_dofs
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use solovev_label_fields, only: triplet_list
    use solovev_doublet, only: solovev_doublet_t, doublet_f0
    implicit none
    integer, parameter :: exact_label = 1, projected_label = 2, constrained_label = 3
    character(11), parameter :: variant_name(3) = ['exact      ', 'projected  ', 'constrained']
    ! measures: total, lobes, outside, band, disc about X, separatrix layer
    integer, parameter :: nmeasure = 6
    character(8), parameter :: measure_name(nmeasure) = ['total   ', 'lobes   ', 'outside ', &
                                                         'band    ', 'X disc  ', 'layer   ']
    real(dp), parameter :: pi = acos(-1.0_dp), band_width = 0.2_dp
    integer :: cells(4) = [4, 8, 16, 32], levels = 3, degree = 2, level, variant, failures, m
    real(dp) :: power = 4, theta = 0.5_dp, ridge = 1e-12_dp, disc_radius = 0.04_dp
    real(dp), allocatable :: err(:, :, :), misplaced(:)
    real(dp) :: last_misplacement
    real(dp) :: rate(nmeasure)
    type(solovev_doublet_t) :: d
    type(block_layout_t) :: layout
    ! the continuous space of the label, for the projection callback
    type(spline_space_t), pointer :: space_of_label => null()
    character(16) :: arg

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) power
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) theta
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) degree
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) disc_radius
    end if
    call d%init()
    print '(a,2f12.8,a,2f12.8,a,es12.4)', 'X', d%xx, '  O1', d%xo(:, 1), &
        '  depth', d%psi_x - d%psi_o
    call build_layout()
    print '(a,i2,a,f5.1,a,f5.2)', 'degree', degree, '  grading power', power, '  layer fraction', theta
    allocate (err(levels, nmeasure, 3), misplaced(levels))
    do level = 1, levels
        do variant = 1, 3
            call fit(cells(level), variant, err(level, :, variant))
            if (variant == exact_label) misplaced(level) = last_misplacement
            print '(i4,1x,a,6es11.3)', cells(level), variant_name(variant), err(level, :, variant)
        end do
    end do
    print '(a,6(1x,a8))', 'rates             ', measure_name
    do variant = 1, 3
        do level = 2, levels
            print '(a,i3,1x,a,6f9.2)', 'to', cells(level), variant_name(variant), &
                log(err(level - 1, :, variant)/err(level, :, variant))/log(2.0_dp)
        end do
    end do
    ! projected and constrained labels: order p in every region, the X disc
    ! included, and a layer error that falls (no current sheet)
    do variant = projected_label, constrained_label
        rate = log(err(levels - 1, :, variant)/err(levels, :, variant))/log(2.0_dp)
        do m = 1, 5
            call check(rate(m) > degree - 0.5_dp, trim(variant_name(variant))//': '// &
                       trim(measure_name(m))//' converges at order p')
        end do
        call check(rate(6) > 0.5_dp, trim(variant_name(variant))// &
                   ': the separatrix layer error falls (no current sheet)')
    end do
    ! exact label: order p away from the separatrix; near it the mesh
    ! separatrix, misplaced at order p + 1, puts points on the wrong side
    rate = log(err(levels - 1, :, exact_label)/err(levels, :, exact_label))/log(2.0_dp)
    call check(rate(2) > degree - 0.5_dp .and. rate(3) > degree - 0.5_dp, &
               'exact: lobes and outside converge at order p')
    call check(log(misplaced(levels - 1)/misplaced(levels))/log(2.0_dp) > degree + 0.5_dp, &
               'the mesh separatrix converges at order p + 1')
    if (failures > 0) error stop 'aligned doublet gate E2b failed'

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

    subroutine build_layout()
        !! The separatrix-aligned doublet layout from the traced lobes and the
        !! wall chi = 2; checks that blocks 1-5 lie in lobe 1 (above X), 6-10
        !! in lobe 2 and 11-18 outside.
        type(doublet_curves_t) :: curves
        real(dp) :: centre(2)
        integer :: st, k, j, region

        curves%xpoint = d%xx
        curves%opoint = d%xo
        do k = 1, 2
            call d%lobe_samples(k, 2e-3_dp, curves%lobe(k)%points)
        end do
        call d%wall_samples(curves%wall)
        call layout_doublet(curves, 0.4_dp, layout, st)
        if (st /= 0) error stop 'doublet layout'
        ! lobe blocks: the corner mean lies in their lobe; ring blocks wrap
        ! around a lobe, so check instead that they reach the wall chi = 2
        do k = 1, layout%block_count
            centre = sum(layout%vertices(:, layout%corners(:, k)), 2)/4
            if (k <= 10) then
                region = 3
                if (d%label(centre(1), centre(2)) < 1) region = merge(1, 2, centre(2) > d%xx(2))
                if (region /= region_of(k)) error stop 'doublet layout: lobe block order'
            else
                if (maxval([(d%label(layout%vertices(1, layout%corners(j, k)), &
                                     layout%vertices(2, layout%corners(j, k))), j=1, 4)]) &
                    < 1.9_dp) error stop 'doublet layout: ring block order'
            end if
        end do
        print '(a,i3,a,i3)', 'blocks', layout%block_count, '  separatrix sides', &
            count(layout%side_tag == LAYOUT_SEPARATRIX)
    end subroutine build_layout

    subroutine layer_grading(knots, n)
        !! Replace the power-graded ends that spline_grade_knots set (a first
        !! or last cell below half the uniform size) by graded layers.
        type(knot_spec_t), intent(inout) :: knots(:, :)
        integer, intent(in) :: n
        logical :: at0, at1
        integer :: k, dd
        real(dp) :: mirrored(n/2 + 1)

        do k = 1, size(knots, 2)
            do dd = 1, 2
                if (.not. allocated(knots(dd, k)%breakpoints)) cycle
                associate (b => knots(dd, k)%breakpoints)
                    at0 = b(2) - b(1) < 0.5_dp/n
                    at1 = b(n + 1) - b(n) < 0.5_dp/n
                end associate
                if (at0 .and. at1) then
                    mirrored = 1 - 0.5_dp*layered_reverse(n/2)
                    knots(dd, k)%breakpoints = [0.5_dp*layered(n/2), mirrored(2:)]
                else if (at0) then
                    knots(dd, k)%breakpoints = layered(n)
                else if (at1) then
                    knots(dd, k)%breakpoints = 1 - [layered_reverse(n)]
                end if
            end do
        end do
    end subroutine layer_grading

    function layered(n) result(x)
        !! n + 1 breakpoints on [0, 1], graded towards 0 in a layer.
        integer, intent(in) :: n
        real(dp) :: x(n + 1), tau
        integer :: m, i

        m = max(1, nint(theta*n))
        if (m >= n) then
            x = [((real(i, dp)/n)**power, i=0, n)]
            return
        end if
        tau = real(m, dp)/n/(power*(1 - real(m, dp)/n) + real(m, dp)/n)
        do i = 0, m
            x(i + 1) = tau*(real(i, dp)/m)**power
        end do
        do i = m + 1, n
            x(i + 1) = tau + (1 - tau)*real(i - m, dp)/(n - m)
        end do
    end function layered

    function layered_reverse(n) result(x)
        !! The breakpoints of layered(n) from the far end: 1 - x(n + 1 - i),
        !! so 1 - layered_reverse(n) runs from the layer at 1 back to 0
        !! after reversal; returned in decreasing order of distance to 1.
        integer, intent(in) :: n
        real(dp) :: x(n + 1), y(n + 1)
        integer :: i

        y = layered(n)
        x = [(y(n + 2 - i), i=1, n + 1)]
    end function layered_reverse

    pure integer function region_of(block)
        integer, intent(in) :: block

        region_of = 3
        if (block <= 10) region_of = (block - 1)/5 + 1
    end function region_of

    subroutine build_space(n, cut, space, derham)
        integer, intent(in) :: n
        logical, intent(in) :: cut
        type(spline_space_t), intent(out) :: space
        type(spline_derham_t), intent(out) :: derham
        type(spline_space_options_t) :: options
        integer :: st

        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_grade_knots(layout, LAYOUT_SEPARATRIX, 1.0_dp, options%knots, st, power=power)
        if (st /= 0) error stop 'grading'
        if (theta > 0) call layer_grading(options%knots, n)
        if (cut) options%cut_tags = [LAYOUT_SEPARATRIX]
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
    end subroutine build_space

    subroutine fit(n, variant, measures)
        integer, intent(in) :: n, variant
        real(dp), intent(out) :: measures(nmeasure)
        type(spline_space_t), target :: space_c, space_g
        type(spline_derham_t) :: derham_c, derham_g
        type(spline_quadrature_t) :: quad, q1, aquad
        type(triplet_list) :: tp
        type(csr_t) :: mp
        real(dp), allocatable :: chi_h(:), v0(:, :, :), gn(:, :, :), vc(:, :, :), gc(:, :, :), &
                                 rhs(:), g(:), row(:)
        integer, allocatable :: d0(:), dc(:), none(:), sep(:)
        real(dp) :: s(2), jet(2, 6), w, gx(2), c, b_phi, target, oh(2, 2), msq(2, 2, 2), &
                    acc(6), den, layer_den, e2
        integer :: st, n0, cell, q, p, i, j, pass, region

        call build_space(n, .false., space_c, derham_c)
        call build_space(n, .true., space_g, derham_g)
        n0 = derham_g%dimension(0)
        allocate (chi_h(derham_c%dimension(0)), rhs(n0), g(n0), none(0))
        if (variant /= exact_label) then
            space_of_label => space_c
            call spline_derham_project(derham_c, 0, label_cb, chi_h, st, 12)
            if (st /= 0) error stop 'label projection'
            if (variant == constrained_label) then
                call spline_derham_boundary_dofs(derham_c, LAYOUT_SEPARATRIX, 0, sep)
                if (size(sep) == 0) error stop 'no separatrix nodes'
                chi_h(sep) = 1
            end if
        end if
        ! winding centres where grad chi = 0, in the frame H^(1/2) of chi there
        do j = 1, 2
            if (variant == exact_label) then
                oh(:, j) = d%xo(:, j)
            else
                oh(:, j) = minimum_h(derham_c, space_c, chi_h, d%xo(:, j))
            end if
            msq(:, :, j) = sqrt_hessian(derham_c, space_c, chi_h, oh(:, j), variant)
        end do
        call spline_quadrature_build(space_g, degree + 8, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_quadrature_point_singular(space_g, quad, oh(:, 1), q1, st)
        if (st /= 0) error stop 'polar O1'
        call spline_quadrature_point_singular(space_g, q1, oh(:, 2), aquad, st)
        if (st /= 0) error stop 'polar O2'
        rhs = 0
        g = 0
        acc = 0
        den = 0
        layer_den = 0
        do pass = 1, 2
            if (pass == 1) call tp%init(aquad%point_count*(degree + 1)**4)
            do cell = 1, size(aquad%point_offset)
                call clebsch_node_gradients(derham_g, aquad, cell, v0, d0, gn)
                if (variant /= exact_label) &
                    call clebsch_node_gradients(derham_c, aquad, cell, vc, dc, gc)
                if (allocated(row)) deallocate (row)
                allocate (row(size(d0)))
                do q = 1, aquad%points_per_cell
                    p = aquad%point_offset(cell) + q - 1
                    call spline_geometry_at(space_g, aquad%block(p), aquad%xi(:, p), jet, st)
                    s = jet(:, 1)
                    w = 2*pi*s(1)*aquad%weight(p)
                    region = region_of(aquad%block(p))
                    if (variant == exact_label) then
                        c = d%label(s(1), s(2))
                        gx = d%label_gradient(s(1), s(2))
                    else
                        c = dot_product(vc(1, :, q), chi_h(dc))
                        gx = matmul(gc(:, :, q), chi_h(dc))
                    end if
                    target = doublet_f0/s(1)
                    b_phi = winding_part(region, c, gx, s, oh, msq)
                    row = gn(2, :, q)*gx(1) - gn(1, :, q)*gx(2)
                    if (pass == 1) then
                        rhs(d0) = rhs(d0) + w*row*(target - b_phi)
                        do i = 1, size(d0)
                            do j = 1, size(d0)
                                call tp%add(d0(i), d0(j), w*row(i)*row(j))
                            end do
                        end do
                    else
                        b_phi = b_phi + dot_product(row, g(d0))
                        e2 = w*(b_phi - target)**2
                        c = d%label(s(1), s(2))
                        acc(1) = acc(1) + e2
                        if (c < 0.3_dp) acc(2) = acc(2) + e2
                        if (c > 1.7_dp .and. c < 1.95_dp) acc(3) = acc(3) + e2
                        if (norm2(s - d%xx) < disc_radius) then
                            acc(5) = acc(5) + e2
                        else if (abs(c - 1) < band_width) then
                            acc(4) = acc(4) + e2
                            if (abs(c - 1) < 0.5_dp/n) then
                                acc(6) = acc(6) + e2
                                layer_den = layer_den + w*target**2
                            end if
                        end if
                        den = den + w*target**2
                    end if
                end do
            end do
            if (pass == 1) then
                do i = 1, n0
                    call tp%add(i, i, 0.0_dp)
                end do
                call tp%to_csr(n0, n0, mp)
                call scaled_solve(mp, n0, rhs, none, g, st)
                if (st /= 0) error stop 'g fit'
            end if
        end do
        measures(1:5) = sqrt(acc(1:5)/den)
        measures(6) = sqrt(acc(6)/layer_den)
        if (variant == exact_label) then
            last_misplacement = misplacement(space_g)
            print '(a,i4,a,es10.3)', 'cells', n, '  separatrix misplacement max |chi - 1|', &
                last_misplacement
        end if
    end subroutine fit

    real(dp) function winding_part(region, c, gx, s, oh, msq) result(b_phi)
        !! sum_j F_j'(chi) (grad chi x h_j)_phi with the profiles of the region,
        !! evaluated at |chi - 1| also where the label of a point lies on the
        !! other side of the separatrix than its block.
        integer, intent(in) :: region
        real(dp), intent(in) :: c, gx(2), s(2), oh(2, 2), msq(2, 2, 2)
        real(dp) :: fpp(2), ex(2), th(2)
        integer :: j

        fpp = 0
        if (region <= 2) then
            fpp(region) = -d%lobe_slope(c)/(2*pi)
        else
            fpp = -d%outer_slope(c)/(4*pi)
        end if
        b_phi = 0
        do j = 1, 2
            if (fpp(j) == 0) cycle
            ! theta_j = atan2 of u = M_j (s - O_j), grad theta = M_j (-u2, u1)/|u|^2
            ex = matmul(msq(:, :, j), s - oh(:, j))
            th = matmul(msq(:, :, j), [-ex(2), ex(1)])/sum(ex**2)
            b_phi = b_phi - fpp(j)*(gx(1)*th(2) - gx(2)*th(1))
        end do
    end function winding_part

    real(dp) function misplacement(space) result(worst)
        !! max |chi - 1| of the exact label at Gauss points of the mesh
        !! separatrix.
        type(spline_space_t), intent(in) :: space
        type(spline_side_quadrature_t) :: side
        integer :: q, st

        call spline_interface_quadrature(space, LAYOUT_SEPARATRIX, degree + 2, side, st)
        if (st /= 0) error stop 'separatrix quadrature'
        worst = 0
        do q = 1, side%point_count
            worst = max(worst, abs(d%label(side%s(1, q), side%s(2, q)) - 1))
        end do
    end function misplacement

    function minimum_h(derham, space, chi_h, x0) result(x)
        !! Newton on grad chi_h = 0 by central differences of the spline.
        type(spline_derham_t), intent(in) :: derham
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: chi_h(:), x0(2)
        real(dp) :: x(2), gr(2), hm(2, 2), det
        integer :: it

        x = x0
        do it = 1, 20
            gr = fd_gradient(derham, space, chi_h, x)
            hm = fd_hessian(derham, space, chi_h, x)
            det = hm(1, 1)*hm(2, 2) - hm(1, 2)**2
            x = x - [hm(2, 2)*gr(1) - hm(1, 2)*gr(2), -hm(2, 1)*gr(1) + hm(1, 1)*gr(2)]/det
        end do
    end function minimum_h

    function sqrt_hessian(derham, space, chi_h, x, variant) result(m)
        !! H^(1/2) of the label at x: (H + sqrt(det H) I)/sqrt(tr H + 2 sqrt(det H)).
        type(spline_derham_t), intent(in) :: derham
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: chi_h(:), x(2)
        integer, intent(in) :: variant
        real(dp) :: m(2, 2), hm(2, 2), sd
        real(dp), parameter :: h = 1e-4_dp

        if (variant == exact_label) then
            hm(:, 1) = (d%label_gradient(x(1) + h, x(2)) - d%label_gradient(x(1) - h, x(2)))/(2*h)
            hm(:, 2) = (d%label_gradient(x(1), x(2) + h) - d%label_gradient(x(1), x(2) - h))/(2*h)
            hm(1, 2) = (hm(1, 2) + hm(2, 1))/2
            hm(2, 1) = hm(1, 2)
        else
            hm = fd_hessian(derham, space, chi_h, x)
        end if
        sd = sqrt(hm(1, 1)*hm(2, 2) - hm(1, 2)**2)
        m = hm
        m(1, 1) = m(1, 1) + sd
        m(2, 2) = m(2, 2) + sd
        m = m/sqrt(hm(1, 1) + hm(2, 2) + 2*sd)
    end function sqrt_hessian

    function fd_gradient(derham, space, chi_h, x) result(gr)
        type(spline_derham_t), intent(in) :: derham
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: chi_h(:), x(2)
        real(dp) :: gr(2)
        real(dp), parameter :: h = 1e-4_dp

        gr = [value_h(derham, space, chi_h, x + [h, 0.0_dp]) - &
              value_h(derham, space, chi_h, x - [h, 0.0_dp]), &
              value_h(derham, space, chi_h, x + [0.0_dp, h]) - &
              value_h(derham, space, chi_h, x - [0.0_dp, h])]/(2*h)
    end function fd_gradient

    function fd_hessian(derham, space, chi_h, x) result(hm)
        type(spline_derham_t), intent(in) :: derham
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: chi_h(:), x(2)
        real(dp) :: hm(2, 2), v0
        real(dp), parameter :: h = 1e-4_dp

        v0 = value_h(derham, space, chi_h, x)
        hm(1, 1) = (value_h(derham, space, chi_h, x + [h, 0.0_dp]) - 2*v0 + &
                    value_h(derham, space, chi_h, x - [h, 0.0_dp]))/h**2
        hm(2, 2) = (value_h(derham, space, chi_h, x + [0.0_dp, h]) - 2*v0 + &
                    value_h(derham, space, chi_h, x - [0.0_dp, h]))/h**2
        hm(1, 2) = (value_h(derham, space, chi_h, x + [h, h]) - &
                    value_h(derham, space, chi_h, x + [h, -h]) - &
                    value_h(derham, space, chi_h, x + [-h, h]) + &
                    value_h(derham, space, chi_h, x - [h, h]))/(4*h**2)
        hm(2, 1) = hm(1, 2)
    end function fd_hessian

    real(dp) function value_h(derham, space, chi_h, y)
        type(spline_derham_t), intent(in) :: derham
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: chi_h(:), y(2)
        real(dp) :: xi(2), v(1)
        integer :: block, st

        call spline_locate(space, 1, y, block, xi, st)
        if (st /= 0) error stop 'locate'
        call spline_derham_evaluate(derham, 0, chi_h, block, xi, 0.0_dp, v)
        value_h = v(1)
    end function value_h

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

    subroutine scaled_solve(m, n, rhs, none, x, status)
        !! Solve the normal equations with symmetric Jacobi scaling: the graded
        !! cells spread the diagonal over many orders of magnitude, and a ridge
        !! relative to the largest entry swamps the rows of the smallest cells
        !! (the error then grows under refinement). With D = diag(M), solve
        !! (D^-1/2 M D^-1/2 + ridge I) y = D^-1/2 rhs and x = D^-1/2 y, refined
        !! against M; the ridge removes g -> g + k(chi) by a uniformly small
        !! relative perturbation.
        type(csr_t), intent(inout) :: m
        integer, intent(in) :: n
        real(dp), intent(in) :: rhs(:)
        integer, intent(in) :: none(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        real(dp) :: scale(n), r(n), dx(n)
        type(csr_t) :: m0
        integer :: i, p, it

        m0 = m
        scale = 1
        do i = 1, n
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                if (m%col(p) == i .and. m%val(p) > 0) scale(i) = 1/sqrt(m%val(p))
            end do
        end do
        do i = 1, n
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                m%val(p) = scale(i)*m%val(p)*scale(m%col(p))
                if (m%col(p) == i) m%val(p) = m%val(p) + ridge
            end do
        end do
        call eu_constrained_solve(m, scale*rhs, none, [real(dp) ::], x, status)
        x = scale*x
        ! iterative refinement against the unscaled matrix: at strong grading
        ! the condition number reaches 1e14 (power 8, 32 cells: residual 2e-2
        ! after the factor, 1e-5 after two corrections)
        do it = 1, 5
            r = rhs - matvec(m0, x)
            if (norm2(r) <= 1e-10_dp*norm2(rhs)) exit
            call eu_constrained_solve(m, scale*r, none, [real(dp) ::], dx, status)
            x = x + scale*dx
        end do
    end subroutine scaled_solve

    function matvec(m, x) result(y)
        type(csr_t), intent(in) :: m
        real(dp), intent(in) :: x(:)
        real(dp) :: y(m%nrow)
        integer :: i, p

        y = 0
        do i = 1, m%nrow
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                y(i) = y(i) + m%val(p)*x(m%col(p))
            end do
        end do
    end function matvec

end program test_eulerian_clebsch_doublet_aligned
