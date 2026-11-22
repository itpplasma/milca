program test_label_e5_unaligned
    !! Gate E5u of docs/eulerian-label-plan.md: can the direct Clebsch label
    !! lane carry a separatrix (island chain, doublet) on a mesh NOT aligned
    !! with it, with only the O- and X-points as tracked points? Case: the
    !! Solov'ev doublet of E2a (two O-points, one X-point, one separatrix;
    !! the axisymmetric proxy of a two-island chain: in helical symmetry the
    !! label representation is the same with a different metric). Labels
    !! fixed to chi_h = Pi0 chi on a plain O-grid centred at X (unaligned),
    !! exact profiles per region, the toroidal field f0/R fitted in L2 by
    !!
    !!   mode 1: g in continuous V0 (E2a),
    !!   mode 2: g broken across the level set chi_h = 1 (one copy of V0 per
    !!           region, region of each quadrature point from chi_h and the
    !!           side of X: Eulerian jump capture, the mesh unaligned),
    !!   mode 3: mode 2 with the log-singular flux slope factored out,
    !!           g = F_r'(chi_h) zeta, zeta broken: B_phi = F_r' (grad zeta x
    !!           grad chi_h + grad chi_h x h_r), zeta = theta - theta* bounded.
    !!
    !! Errors of B_phi: total, lobes (chi < 0.3), outside (1.7 < chi < 1.95),
    !! band (|chi - 1| < 0.2). The interpolation level of a smooth field on
    !! the same mesh is order p (lobes and outside show it).
    !!
    !! Arguments: mode (default 3), levels (3: 4, 8, 16; 4 adds 32), Gauss
    !! points per direction (degree + 8).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_cell_basis, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use solovev_label_fields, only: triplet_list
    use solovev_doublet, only: solovev_doublet_t, doublet_f0
    implicit none
    integer, parameter :: degree = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: levels = 3, level, cells(4) = [4, 8, 16, 32], failures, mode = 3, &
               points = degree + 8
    real(dp) :: err(4, 4), rate(4)
    real(dp) :: ridge = 1e-12_dp
    logical :: elliptic = .true.
    character(16) :: arg
    type(spline_space_t) :: space
    type(solovev_doublet_t) :: d

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) mode
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) points
    end if
    call d%init()
    print '(a,i0,a,i0)', 'mode ', mode, '  Gauss points per direction ', points
    print '(a)', 'cells  L2 error of B_phi: total, lobes, outside, band |chi - 1| < 0.2'
    do level = 1, levels
        call fit(cells(level), err(level, 1), err(level, 2), err(level, 3), err(level, 4))
        print '(i5,4es12.3)', cells(level), err(level, :)
    end do
    do level = 2, levels
        print '(a,i3,4f8.2)', 'rates to', cells(level), &
            log(err(level - 1, :)/err(level, :))/log(2.0_dp)
    end do
    rate = log(err(levels - 1, :)/err(levels, :))/log(2.0_dp)
    call check(rate(2) > degree - 0.3_dp, 'lobes converge at order p')
    call check(rate(3) > degree - 0.3_dp, 'outside the separatrix converges at order p')
    ! the verdict measurement: on the unaligned mesh the separatrix band
    ! stays below order 1 for every mode (log layer normal to the separatrix)
    call check(rate(4) < 1, 'the separatrix band stays below order 1 on the unaligned mesh')
    if (failures > 0) error stop 'gate E5u failed'

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

    subroutine fit(n, rel_phi, rel_lobe, rel_out, rel_band)
        integer, intent(in) :: n
        real(dp), intent(out) :: rel_phi, rel_lobe, rel_out, rel_band
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad, q1, aquad
        type(spline_derham_t) :: derham
        type(triplet_list) :: tp
        type(csr_t) :: mp
        real(dp), allocatable :: samples(:, :), chi_h(:), v0(:, :, :), gn(:, :, :), rhs(:), &
                                 g(:), row(:)
        integer, allocatable :: d0(:), none(:)
        real(dp) :: s(2), jet(2, 6), w, gx(2), c, fpp(2), b_phi, num, den, &
                    target, ex(2), far, outer, oh(2, 2), msq(2, 2, 2), th(2), band, scale
        integer :: st, n0, cell, q, p, i, j, pass, region, nr
        integer, allocatable :: dr(:)

        call d%wall_samples(samples)
        call layout_ogrid(d%xx, 0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :))), 0.5_dp, &
                          layout, st)
        if (st /= 0) error stop 'layout'
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'wall curve'
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_quadrature_build(space, points, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
        n0 = derham%dimension(0)
        ! unknowns: one V0 copy (mode 1) or one per region (lobe 1, lobe 2,
        ! outside)
        nr = merge(1, 3, mode == 1)
        allocate (chi_h(n0), rhs(nr*n0), g(nr*n0), none(0))
        call spline_derham_project(derham, 0, label_cb, chi_h, st, 12)
        if (st /= 0) error stop 'label projection'
        ! winding centres at the minima of chi_h: there grad chi_h = 0, so
        ! grad chi_h x h_j stays bounded
        do j = 1, 2
            oh(:, j) = minimum_h(derham, chi_h, d%xo(:, j))
            msq(:, :, j) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
            if (elliptic) msq(:, :, j) = sqrt_hessian(derham, chi_h, oh(:, j))
        end do
        call spline_quadrature_point_singular(space, quad, oh(:, 1), q1, st)
        if (st /= 0) error stop 'polar O1'
        call spline_quadrature_point_singular(space, q1, oh(:, 2), aquad, st)
        if (st /= 0) error stop 'polar O2'
        rhs = 0
        g = 0
        num = 0
        den = 0
        far = 0
        outer = 0
        band = 0
        do pass = 1, 2
            if (pass == 1) call tp%init(aquad%point_count*81)
            do cell = 1, size(aquad%point_offset)
                call clebsch_node_gradients(derham, aquad, cell, v0, d0, gn)
                if (allocated(row)) deallocate (row)
                allocate (row(size(d0)))
                do q = 1, aquad%points_per_cell
                    p = aquad%point_offset(cell) + q - 1
                    call spline_geometry_at(space, aquad%block(p), aquad%xi(:, p), jet, st)
                    s = jet(:, 1)
                    w = 2*pi*s(1)*aquad%weight(p)
                    c = dot_product(v0(1, :, q), chi_h(d0))
                    gx = matmul(gn(:, :, q), chi_h(d0))
                    ! F_j' per region, F = -Phi/(2 pi)
                    if (c < 1) then
                        fpp = 0
                        j = merge(1, 2, s(2) > d%xx(2))
                        fpp(j) = -d%lobe_slope(c)/(2*pi)
                        region = j
                        scale = fpp(j)
                    else
                        fpp = -d%outer_slope(c)/(4*pi)
                        region = 3
                        scale = 2*fpp(1)
                    end if
                    if (mode == 1) region = 1
                    dr = (region - 1)*n0 + d0
                    target = doublet_f0/s(1)
                    ! the h_j part of B_phi: -F_j' (chi_R theta_j,Z - chi_Z theta_j,R)
                    b_phi = 0
                    do j = 1, 2
                        ! theta_j = atan2 of u = M_j (s - O_j), grad theta = M_j (-u2, u1)/|u|^2
                        ex = matmul(msq(:, :, j), s - oh(:, j))
                        th = matmul(msq(:, :, j), [-ex(2), ex(1)])/sum(ex**2)
                        b_phi = b_phi - fpp(j)*(gx(1)*th(2) - gx(2)*th(1))
                    end do
                    row = gn(2, :, q)*gx(1) - gn(1, :, q)*gx(2)
                    ! mode 3: g = F_r'(chi_h) zeta, grad g x grad chi_h = F_r' grad zeta x grad chi_h
                    if (mode == 3) row = scale*row
                    if (pass == 1) then
                        rhs(dr) = rhs(dr) + w*row*(target - b_phi)
                        do i = 1, size(d0)
                            do j = 1, size(d0)
                                call tp%add(dr(i), dr(j), w*row(i)*row(j))
                            end do
                        end do
                    else
                        b_phi = b_phi + dot_product(row, g(dr))
                        num = num + w*(b_phi - target)**2
                        if (c < 0.3_dp) far = far + w*(b_phi - target)**2
                        if (c > 1.7_dp .and. c < 1.95_dp) outer = outer + w*(b_phi - target)**2
                        if (abs(c - 1) < 0.2_dp) band = band + w*(b_phi - target)**2
                        den = den + w*target**2
                    end if
                end do
            end do
            if (pass == 1) then
                do i = 1, nr*n0
                    call tp%add(i, i, 0.0_dp)
                end do
                call tp%to_csr(nr*n0, nr*n0, mp)
                call jacobi_ridge(mp, nr*n0, ridge)
                call eu_constrained_solve(mp, rhs, none, [real(dp) ::], g, st)
                if (st /= 0) error stop 'g fit'
            end if
        end do
        rel_phi = sqrt(num/den)
        rel_lobe = sqrt(far/den)
        rel_out = sqrt(outer/den)
        rel_band = sqrt(band/den)
    end subroutine fit

    subroutine jacobi_ridge(m, n, eps)
        !! Ridge eps times the largest diagonal on every diagonal entry, and
        !! that diagonal on empty rows (V0 copies unused by their region).
        type(csr_t), intent(inout) :: m
        integer, intent(in) :: n
        real(dp), intent(in) :: eps
        real(dp) :: dmax
        integer :: i, p
        dmax = 0
        do i = 1, n
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                if (m%col(p) == i) dmax = max(dmax, m%val(p))
            end do
        end do
        do i = 1, n
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                if (m%col(p) == i) then
                    if (m%val(p) <= 0) then
                        m%val(p) = dmax
                    else
                        m%val(p) = m%val(p) + eps*dmax
                    end if
                end if
            end do
        end do
    end subroutine jacobi_ridge

    function minimum_h(derham, chi_h, x0) result(x)
        !! Newton on grad chi_h = 0 by central differences of the spline.
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: chi_h(:), x0(2)
        real(dp) :: x(2), g(2), hm(2, 2), det
        real(dp), parameter :: h = 1e-4_dp
        integer :: it

        x = x0
        do it = 1, 20
            g = [value_h(derham, chi_h, x + [h, 0.0_dp]) - value_h(derham, chi_h, x - [h, 0.0_dp]), &
                 value_h(derham, chi_h, x + [0.0_dp, h]) - value_h(derham, chi_h, x - [0.0_dp, h])]/(2*h)
            hm(1, 1) = (value_h(derham, chi_h, x + [h, 0.0_dp]) - 2*value_h(derham, chi_h, x) + value_h(derham, chi_h, x - [h, 0.0_dp]))/h**2
            hm(2, 2) = (value_h(derham, chi_h, x + [0.0_dp, h]) - 2*value_h(derham, chi_h, x) + value_h(derham, chi_h, x - [0.0_dp, h]))/h**2
            hm(1, 2) = (value_h(derham, chi_h, x + [h, h]) - value_h(derham, chi_h, x + [h, -h]) - value_h(derham, chi_h, x + [-h, h]) &
                        + value_h(derham, chi_h, x - [h, h]))/(4*h**2)
            hm(2, 1) = hm(1, 2)
            det = hm(1, 1)*hm(2, 2) - hm(1, 2)**2
            x = x - [hm(2, 2)*g(1) - hm(1, 2)*g(2), -hm(2, 1)*g(1) + hm(1, 1)*g(2)]/det
        end do

    end function minimum_h

    function sqrt_hessian(derham, chi_h, x) result(m)
        !! H^(1/2) of chi_h at x: (H + sqrt(det H) I)/sqrt(tr H + 2 sqrt(det H)).
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: chi_h(:), x(2)
        real(dp) :: m(2, 2), hm(2, 2), sd
        real(dp), parameter :: h = 1e-4_dp

        hm(1, 1) = (value_h(derham, chi_h, x + [h, 0.0_dp]) - 2*value_h(derham, chi_h, x) &
                    + value_h(derham, chi_h, x - [h, 0.0_dp]))/h**2
        hm(2, 2) = (value_h(derham, chi_h, x + [0.0_dp, h]) - 2*value_h(derham, chi_h, x) &
                    + value_h(derham, chi_h, x - [0.0_dp, h]))/h**2
        hm(1, 2) = (value_h(derham, chi_h, x + [h, h]) - value_h(derham, chi_h, x + [h, -h]) &
                    - value_h(derham, chi_h, x + [-h, h]) + value_h(derham, chi_h, x - [h, h])) &
                   /(4*h**2)
        hm(2, 1) = hm(1, 2)
        sd = sqrt(hm(1, 1)*hm(2, 2) - hm(1, 2)**2)
        m = hm
        m(1, 1) = m(1, 1) + sd
        m(2, 2) = m(2, 2) + sd
        m = m/sqrt(hm(1, 1) + hm(2, 2) + 2*sd)
    end function sqrt_hessian

    real(dp) function value_h(derham, chi_h, y)
        type(spline_derham_t), intent(in) :: derham
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

        call spline_geometry_at(space, block, xi, jet, st)
        values(1) = d%label(jet(1, 1), jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

    subroutine add_ridge(m, n, eps)
        type(csr_t), intent(inout) :: m
        integer, intent(in) :: n
        real(dp), intent(in) :: eps
        integer :: i, p

        do i = 1, n
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                if (m%col(p) == i) m%val(p) = m%val(p) + eps
            end do
        end do
    end subroutine add_ridge

end program test_label_e5_unaligned
