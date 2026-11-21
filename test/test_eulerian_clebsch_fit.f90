program test_eulerian_clebsch_fit
    !! Gate E0 of docs/eulerian-label-plan.md: approximation property of the
    !! direct Clebsch representation with fixed discrete labels chi_h = Pi0
    !! chi. The exact Solov'ev field is fitted in L2 by
    !!
    !!   B = grad g x grad chi_h + F'(chi_h) grad chi_h x h_theta
    !!       + G'(chi_h) grad chi_h x grad phi,
    !!
    !! g in V0, F and G cubic B-splines in chi with F(chi_axis) = 0, h_theta
    !! the winding form about the exact axis. In axisymmetry (e_R x e_Z =
    !! -e_phi) this is
    !!
    !!   B_phi = g_Z chi_R - g_R chi_Z - F' (chi_R theta_Z - chi_Z theta_R),
    !!   (B_R, B_Z) = G' (-chi_Z, chi_R)/R.
    !!
    !! The multiplier formulation locks for the nonlinear label (13 % at
    !! every mesh, docs/eulerian-3d.md). Here the fit error must converge at
    !! order p for both labels.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_locate
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_cell_basis, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use solovev_label_fields
    implicit none
    integer, parameter :: degree = 2, levels = 4
    integer :: level, cells(levels) = [2, 4, 8, 16], failures
    real(dp) :: error(levels), rate
    character(16) :: arg
    ! relative ridge: removes g -> g + k(chi_h) and unsupported profile splines
    real(dp) :: ridge = 1e-12_dp

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        nonlinear_label = arg == 'nonlinear'
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) ridge
    end if
    call curve_chart_cylindrical(chart, 1)
    do level = 1, levels
        error(level) = fit_error(cells(level))
        print '(a,i3,a,es10.3)', 'cells', cells(level), ' L2 error of the Clebsch fit', &
            error(level)
    end do
    do level = 2, levels
        rate = log(error(level - 1)/error(level))/log(2.0_dp)
        print '(a,f6.2)', 'rate', rate
        call check(rate > degree - 0.3_dp, 'Clebsch fit converges at order p')
    end do
    call check(error(levels) < 1e-3_dp, 'fit error below 1e-3 at 16 cells')
    if (failures > 0) error stop 'Clebsch fit gate E0 failed'

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

    pure subroutine cubic_bspline(u, b, db)
        !! Cardinal cubic B-spline on [0, 4] and its derivative.
        real(dp), intent(in) :: u
        real(dp), intent(out) :: b, db

        b = 0
        db = 0
        if (u < 0 .or. u >= 4) return
        if (u < 1) then
            b = u**3/6
            db = u**2/2
        else if (u < 2) then
            b = (-3*u**3 + 12*u**2 - 12*u + 4)/6
            db = (-9*u**2 + 24*u - 12)/6
        else if (u < 3) then
            b = (3*u**3 - 24*u**2 + 60*u - 44)/6
            db = (9*u**2 - 48*u + 60)/6
        else
            b = (4 - u)**3/6
            db = -(4 - u)**2/2
        end if
    end subroutine cubic_bspline

    real(dp) function fit_error(n) result(rel)
        integer, intent(in) :: n
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad
        type(spline_derham_t) :: derham
        type(triplet_list) :: tp, tq
        type(csr_t) :: mp, mq
        real(dp), allocatable :: samples(:, :), chi_h(:), v0(:, :, :), v1(:, :, :), &
                                 gn(:, :), rhs_p(:), rhs_q(:), xp(:), xq(:), bp(:), bq(:, :)
        integer, allocatable :: d0(:), d1(:), s0(:), s1(:), cols(:), none(:)
        integer :: iq(4)
        real(dp) :: s(2), js(2, 2), jinv(2, 2), w, gx(2), chi_q, hx(2), x(3), j3(3, 3), &
                    b_ex(3), lo, hi, step, chi_axis, u, bk, dbk, det, v(1), num, den, &
                    bphi, bpol(2)
        integer :: st, n0, nk, np, cell, q, p, i, j, k, e, kk, block, blk, nloc

        call boundary_samples(samples)
        call layout_ogrid([0.5_dp*(minval(samples(1, :)) + maxval(samples(1, :))), 0.0_dp], &
                          0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :))), 0.5_dp, &
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
        call spline_quadrature_build(space, degree + 2, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
        n0 = derham%dimension(0)
        allocate (chi_h(n0))
        call spline_derham_project(derham, 0, label, chi_h, st, 12)
        if (st /= 0) error stop 'label projection'
        call spline_locate(space, 1, [r0, 0.0_dp], block, s, st)
        if (st /= 0) error stop 'axis'
        call spline_derham_evaluate(derham, 0, chi_h, block, s, 0.0_dp, v)
        chi_axis = v(1)
        ! profile splines: nk uniform intervals on [lo, hi], functions k = 1..nk+3
        lo = min(chi_axis, minval(chi_h)) - 0.02_dp
        hi = max(1.0_dp, maxval(chi_h)) + 0.02_dp
        nk = 8 + 4*n
        step = (hi - lo)/nk
        np = nk + 3
        ! phi problem: unknowns g (n0) and F (np); poloidal problem: G (np)
        call tp%init(quad%point_count*(9 + 4)**2)
        call tq%init(quad%point_count*16)
        allocate (rhs_p(n0 + np), rhs_q(np), gn(2, 9), bp(9 + 4), bq(2, 4), cols(9 + 4))
        rhs_p = 0
        rhs_q = 0
        num = 0
        den = 0
        do cell = 1, size(quad%point_offset)
            call spline_derham_cell_basis(derham, quad, cell, 0, v0, d0, s0)
            call spline_derham_cell_basis(derham, quad, cell, 1, v1, d1, s1)
            nloc = size(d0)
            if (size(gn, 2) /= nloc) then
                deallocate (gn, bp, cols)
                allocate (gn(2, nloc), bp(nloc + 4), cols(nloc + 4))
            end if
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                det = js(1, 1)*js(2, 2) - js(1, 2)*js(2, 1)
                jinv = reshape([js(2, 2), -js(2, 1), -js(1, 2), js(1, 1)], [2, 2])/det
                w = 2*pi*s(1)*quad%weight(p)
                ! physical gradients of the node functions from the incidence
                gn = 0
                do j = 1, size(d1)
                    e = d1(j)
                    do kk = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                        do i = 1, nloc
                            if (d0(i) == derham%grad%colidx(kk)) gn(:, i) = gn(:, i) + &
                                derham%grad%val(kk)*s1(j)*v1(:, j, q)
                        end do
                    end do
                end do
                gn = matmul(transpose(jinv), gn)
                chi_q = dot_product(v0(1, :, q), chi_h(d0))
                gx = matmul(gn, chi_h(d0))
                hx = winding(s(1), s(2))
                b_ex = exact_b(x)
                ! B_phi row: g part, then the four F splines at chi_q
                bp(1:nloc) = gn(2, :)*gx(1) - gn(1, :)*gx(2)
                cols(1:nloc) = d0
                k = floor((chi_q - lo)/step)
                do i = 1, 4
                    blk = k + i
                    u = (chi_q - lo)/step - (blk - 4)
                    call cubic_bspline(u, bk, dbk)
                    ! F_blk = (chi - chi_axis) B_blk, F' = B_blk + (chi - chi_axis) B_blk'
                    bp(nloc + i) = -(bk + (chi_q - chi_axis)*dbk/step) &
                                   *(gx(1)*hx(2) - gx(2)*hx(1))
                    cols(nloc + i) = n0 + blk
                    bq(:, i) = dbk/step*[-gx(2), gx(1)]/s(1)
                    iq(i) = blk
                end do
                do i = 1, nloc + 4
                    rhs_p(cols(i)) = rhs_p(cols(i)) + w*bp(i)*b_ex(2)
                    do j = 1, nloc + 4
                        call tp%add(cols(i), cols(j), w*bp(i)*bp(j))
                    end do
                end do
                do i = 1, 4
                    rhs_q(iq(i)) = rhs_q(iq(i)) + w*dot_product(bq(:, i), [b_ex(1), b_ex(3)])
                    do j = 1, 4
                        call tq%add(iq(i), iq(j), w*dot_product(bq(:, i), bq(:, j)))
                    end do
                end do
                den = den + w*sum(b_ex**2)
            end do
        end do
        ! g -> g + k(chi_h) leaves B unchanged, and profile splines outside
        ! the range of chi_h have no support: small ridge on all unknowns.
        ! G is defined up to a constant: fix its first coefficient
        do i = 1, n0 + np
            call tp%add(i, i, 0.0_dp)
        end do
        do i = 1, np
            call tq%add(i, i, 0.0_dp)
        end do
        call tp%to_csr(n0 + np, n0 + np, mp)
        call add_ridge(mp, n0 + np, ridge*maxval(abs(mp%val)))
        call tq%to_csr(np, np, mq)
        call add_ridge(mq, np, ridge*maxval(abs(mq%val)))
        allocate (xp(n0 + np), xq(np), none(0))
        call eu_constrained_solve(mp, rhs_p, none, [real(dp) ::], xp, st)
        if (st /= 0) error stop 'phi fit'
        call eu_constrained_solve(mq, rhs_q, [1], [0.0_dp], xq, st)
        if (st /= 0) error stop 'poloidal fit'
        ! error, second pass over the quadrature points
        do cell = 1, size(quad%point_offset)
            call spline_derham_cell_basis(derham, quad, cell, 0, v0, d0, s0)
            call spline_derham_cell_basis(derham, quad, cell, 1, v1, d1, s1)
            nloc = size(d0)
            if (size(gn, 2) /= nloc) then
                deallocate (gn)
                allocate (gn(2, nloc))
            end if
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                det = js(1, 1)*js(2, 2) - js(1, 2)*js(2, 1)
                jinv = reshape([js(2, 2), -js(2, 1), -js(1, 2), js(1, 1)], [2, 2])/det
                w = 2*pi*s(1)*quad%weight(p)
                gn = 0
                do j = 1, size(d1)
                    e = d1(j)
                    do kk = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                        do i = 1, nloc
                            if (d0(i) == derham%grad%colidx(kk)) gn(:, i) = gn(:, i) + &
                                derham%grad%val(kk)*s1(j)*v1(:, j, q)
                        end do
                    end do
                end do
                gn = matmul(transpose(jinv), gn)
                chi_q = dot_product(v0(1, :, q), chi_h(d0))
                gx = matmul(gn, chi_h(d0))
                hx = winding(s(1), s(2))
                b_ex = exact_b(x)
                bphi = dot_product(gn(2, :)*gx(1) - gn(1, :)*gx(2), xp(d0))
                bpol = 0
                k = floor((chi_q - lo)/step)
                do i = 1, 4
                    blk = k + i
                    u = (chi_q - lo)/step - (blk - 4)
                    call cubic_bspline(u, bk, dbk)
                    bphi = bphi - xp(n0 + blk)*(bk + (chi_q - chi_axis)*dbk/step) &
                           *(gx(1)*hx(2) - gx(2)*hx(1))
                    bpol = bpol + xq(blk)*dbk/step*[-gx(2), gx(1)]/s(1)
                end do
                num = num + w*((bphi - b_ex(2))**2 + (bpol(1) - b_ex(1))**2 + &
                               (bpol(2) - b_ex(3))**2)
            end do
        end do
        rel = sqrt(num/den)
    end function fit_error

    subroutine add_ridge(m, n, eps)
        !! m(i, i) += eps for i <= n (the diagonal entries exist).
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

end program test_eulerian_clebsch_fit
