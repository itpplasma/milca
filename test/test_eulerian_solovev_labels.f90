program test_eulerian_solovev_labels
    !! Gate G2a of docs/eulerian-3d.md: axisymmetric Solov'ev through the
    !! A + chi formulation with the labels fixed to chi = psi/psi_b. Unknowns
    !! A (V1, n x A of the exact field on the wall) and the single-valued
    !! current potential lambda (V0); the toroidal current enters through
    !! C(chi) times the winding form about the axis, C = -I(chi)/(2 pi):
    !!
    !!   int curl dA . B + int lambda curl dA . grad chi = int curl dA . C h_theta
    !!   int nu B . grad chi = 0.
    !!
    !! lambda and nu are splines of degree p - 1, which keeps the pair
    !! inf-sup stable and free of locking for a label that is nonlinear in
    !! psi (argument `nonlinear`). Checks: B within 25 % of the interpolation
    !! error and converging at order p, and the equation not imposed,
    !! B . (grad lambda + C' h_theta) = -dp/dchi, holding weakly with a
    !! residual that falls with the mesh.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_tree, spline_derham_evaluate, &
        spline_derham_cell_basis
    use fortsparse_csr, only: csr_t, csr_matvec, csr_transpose, csr_multiply, csr_budget_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_vacuum, only: eu_curl_curl, eu_constrained_solve, &
        eu_to_csr
    use solovev_label_fields
    implicit none
    integer, parameter :: degree = 2, levels = 4
    ! -eps I on the multiplier block, relative to max|K|: needed only for
    ! an equal-degree multiplier (drop = 0), and it pollutes B on fine meshes
    real(dp) :: epsilon = 0
    ! degree of the multiplier space (lambda, nu): degree - drop. drop = 1
    ! is inf-sup stable here and does not lock for nonlinear labels.
    integer :: drop = 1
    integer :: level, cells(levels) = [2, 4, 8, 16], failures
    real(dp) :: error(levels), residual(levels), constraint, rate, error_interp
    character(16) :: env

    failures = 0
    ! arguments: label (linear | nonlinear), then eps and drop for measurements
    if (command_argument_count() >= 1) then
        call get_command_argument(1, env)
        nonlinear_label = env == 'nonlinear'
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, env)
        read (env, *) drop
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, env)
        read (env, *) epsilon
    end if
    call curve_chart_cylindrical(chart, 1)
    do level = 1, levels
        call solve(cells(level), error(level), error_interp, residual(level), constraint)
        print '(a,i3,a,2es10.3,a,es10.3,a,es10.3)', 'cells', cells(level), &
            ' L2 error of B (Galerkin, interpolant)', error(level), error_interp, &
            ' residual of (4)', residual(level), ' |D a_I|/(|D||a_I|)', constraint
        call check(error(level) < 1.25_dp*error_interp, &
                   'B within 25 % of the interpolation error')
    end do
    do level = 2, levels
        rate = log(error(level - 1)/error(level))/log(2.0_dp)
        print '(a,f6.2,a,f6.2)', 'rate of B', rate, '  rate of (4)', &
            log(residual(level - 1)/residual(level))/log(2.0_dp)
        call check(rate > degree - 0.3_dp, 'L2 error of B converges at order p')
        call check(residual(level) < residual(level - 1)/2, &
                   'residual of the pressure equation falls with the mesh')
    end do
    if (failures > 0) error stop 'Solov''ev labels gate G2a failed'

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

    subroutine solve(n, err, err_i, res, cons)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, err_i, res, cons
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad
        type(spline_derham_t) :: derham, derham_m
        type(spline_space_t) :: space_m
        type(csr_t) :: k, d2t, c, ct, dt, d, saddle
        type(csr_budget_t) :: budget
        type(fortsparse_status_t) :: info
        integer, allocatable :: nodes(:), edges(:), tree(:), fixed(:)
        real(dp), allocatable :: samples(:, :), a_i(:), x(:), rhs(:), load2(:), lam(:), &
                                 a_h(:), da(:), chi_h(:), gchi(:)
        logical, allocatable :: tree_mask(:)
        integer :: st, n0, n1, i

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
        options%degree = degree - drop
        call spline_space_build(layout, options, space_m, st)
        if (st /= 0) error stop 'multiplier space'
        call spline_derham_build(space_m, derham_m, st, 0)
        if (st /= 0) error stop 'multiplier derham'
        n0 = derham_m%dimension(0)
        n1 = derham%dimension(1)
        allocate (a_i(n1))
        call spline_derham_project(derham, 1, potential, a_i, st, 12)
        if (st /= 0) error stop 'projection'
        allocate (chi_h(derham%dimension(0)))
        call spline_derham_project(derham, 0, label, chi_h, st, 12)
        if (st /= 0) error stop 'label projection'
        gchi = apply(derham%grad%n, derham%grad%rowptr, derham%grad%colidx, &
                     derham%grad%val, chi_h)
        call eu_curl_curl(derham, quad, metric2, k, st, 2)
        if (st /= 0) error stop 'curl-curl'
        call assemble_pairings(derham, derham_m, quad, chi_h, gchi, d2t, load2)
        call eu_to_csr(derham%curl, n1, c)
        call csr_transpose(c, ct, budget, info)
        call csr_multiply(ct, d2t, dt, budget, info)
        if (info%code /= 0) error stop 'D^T'
        call csr_transpose(dt, d, budget, info)
        allocate (da(n0))
        call csr_matvec(d, a_i, da)
        cons = norm2(da)/(norm2(d%val)*norm2(a_i))
        call saddle_matrix(k, dt, d, epsilon*maxval(abs(k%val)), saddle)
        allocate (rhs(n1 + n0), x(n1 + n0))
        rhs = 0
        call csr_matvec(ct, load2, rhs(1:n1))
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, nodes)
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 1, edges)
        call spline_derham_tree(derham, tree_mask, roots=nodes)
        tree = pack([(i, i=1, size(tree_mask))], tree_mask)
        ! constants are the exact kernel of D^T on the free A (chi_h = 1 on the
        ! wall): pin one multiplier node
        fixed = [edges, tree, n1 + 1]
        call eu_constrained_solve(saddle, rhs, fixed, [a_i(edges), 0*a_i(tree), 0.0_dp], &
                                  x, st)
        if (st /= 0) print *, 'saddle status', st, ' sizes', n1, n0, size(fixed)
        if (st /= 0) error stop 'saddle solve'
        a_h = x(1:n1)
        lam = x(n1 + 1:)
        err = l2_error(derham, quad, curl(derham, a_h))
        err_i = l2_error(derham, quad, curl(derham, a_i))
        res = pressure_residual(derham, derham_m, quad, chi_h, gchi, curl(derham, a_h), lam)
    end subroutine solve

    subroutine saddle_matrix(k, dt, d, eps, s)
        !! [[K, D^T], [D, -eps I]] as CSR from triplets (symmetric, so CSC =
        !! CSR). eps removes the flux-function kernel of D^T.
        type(csr_t), intent(in) :: k, dt, d
        real(dp), intent(in) :: eps
        type(csr_t), intent(out) :: s
        type(csc_t) :: m
        type(fortsparse_status_t) :: info
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: n1, n0, nz, i

        n1 = k%nrow
        n0 = d%nrow
        nz = size(k%val) + 2*size(d%val) + d%nrow
        allocate (rows(nz), cols(nz), vals(nz))
        nz = 0
        rows(1:n0) = [(n1 + i, i=1, n0)]
        cols(1:n0) = rows(1:n0)
        vals(1:n0) = -eps
        nz = n0
        call append(k, 0, 0, rows, cols, vals, nz)
        call append(dt, 0, n1, rows, cols, vals, nz)
        call append(d, n1, 0, rows, cols, vals, nz)
        call csc_from_triplet(n1 + n0, n1 + n0, rows, cols, vals, m, info)
        if (info%code /= 0) error stop 'saddle triplets'
        s%nrow = n1 + n0
        s%ncol = n1 + n0
        s%row_ptr = m%col_ptr
        s%col = m%row_idx
        s%val = m%val
    end subroutine saddle_matrix

    subroutine append(a, row0, col0, rows, cols, vals, nz)
        type(csr_t), intent(in) :: a
        integer, intent(in) :: row0, col0
        integer, intent(inout) :: rows(:), cols(:), nz
        real(dp), intent(inout) :: vals(:)
        integer :: i, p

        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                nz = nz + 1
                rows(nz) = row0 + i
                cols(nz) = col0 + a%col(p)
                vals(nz) = a%val(p)
            end do
        end do
    end subroutine append

    subroutine assemble_pairings(derham, derham_m, quad, chi_h, gchi, d2t, load2)
        !! Axisymmetric section pairings with the E x T' block of V2, whose
        !! dof on a section edge with 1-form basis (v1, v2) has the 2-form
        !! proxy (b1, b2) = (-v2, v1). With B = J b/det J:
        !!   B . grad chi dV = b . grad_xi chi dxi dt,
        !!   B . w dV = b . J^T w dxi dt,
        !! so D2^T(f, nu) = int N_nu b_f . grad_xi chi dxi and
        !! load2(f) = int b_f . J^T (C h_theta) dxi. nu from derham_m.
        type(spline_derham_t), intent(in) :: derham, derham_m
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi_h(:), gchi(:)
        type(csr_t), intent(out) :: d2t
        real(dp), allocatable, intent(out) :: load2(:)
        real(dp), allocatable :: v0(:, :, :), v1(:, :, :)
        integer, allocatable :: dofs0(:), dofs1(:), s0(:), s1(:), rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(csc_t) :: m
        type(fortsparse_status_t) :: info
        real(dp) :: s(2), js(2, 2), w, g(2), u(2), chi(1), gp(3), bf, pn, dpn
        integer :: cell, q, p, i, j, nf, nz

        nf = derham%face_count
        allocate (load2(derham%dimension(2)))
        load2 = 0
        nz = 0
        nz = quad%point_count*2*derham%degree*(derham%degree + 1)*(derham%degree + 1)**2
        allocate (rows(nz), cols(nz), vals(nz))
        nz = 0
        do cell = 1, size(quad%point_offset)
            call spline_derham_cell_basis(derham_m, quad, cell, 0, v0, dofs0, s0)
            call spline_derham_cell_basis(derham, quad, cell, 1, v1, dofs1, s1)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                w = quad%weight(p)/abs(js(1, 1)*js(2, 2) - js(1, 2)*js(2, 1))
                call spline_derham_evaluate(derham, 0, chi_h, quad%block(p), quad%xi(:, p), &
                                            0.0_dp, chi)
                call spline_derham_evaluate(derham, 1, gchi, quad%block(p), quad%xi(:, p), &
                                            0.0_dp, gp)
                g = gp(1:2)
                call label_inverse(chi(1), pn, dpn)
                u = matmul(transpose(js), enclosed_current(pn)/(2*pi)*winding(s(1), s(2)))
                do j = 1, size(dofs1)
                    bf = s1(j)*(-v1(2, j, q)*g(1) + v1(1, j, q)*g(2))
                    load2(nf + dofs1(j)) = load2(nf + dofs1(j)) + w*s1(j)* &
                                           (-v1(2, j, q)*u(1) + v1(1, j, q)*u(2))
                    do i = 1, size(dofs0)
                        nz = nz + 1
                        rows(nz) = nf + dofs1(j)
                        cols(nz) = dofs0(i)
                        vals(nz) = w*v0(1, i, q)*bf
                    end do
                end do
            end do
        end do
        call csc_from_triplet(derham_m%dimension(0), derham%dimension(2), cols(:nz), &
                              rows(:nz), vals(:nz), m, info)
        if (info%code /= 0) error stop 'pairing triplets'
        d2t%nrow = derham%dimension(2)
        d2t%ncol = derham_m%dimension(0)
        d2t%row_ptr = m%col_ptr
        d2t%col = m%row_idx
        d2t%val = m%val
    end subroutine assemble_pairings

    real(dp) function pressure_residual(derham, derham_m, quad, chi_h, gchi, b2, lam) &
            result(rel)
        !! r(nu) = int nu (B . grad lambda + C'(chi) B . h_theta + dp/dchi) dV
        !! over V0, relative to the same with the pressure term alone.
        type(spline_derham_t), intent(in) :: derham, derham_m
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi_h(:), gchi(:), b2(:), lam(:)
        real(dp), allocatable :: v0(:, :, :), r(:), rp(:), glam(:)
        integer, allocatable :: dofs0(:), s0(:), wall(:)
        real(dp) :: s(2), js(2, 2), w, proxy(3), gl(3), chi(1), dc, dpdchi, f, x(3), &
                    j3(3, 3), pn, dpn
        real(dp), parameter :: h = 1e-5_dp
        integer :: cell, q, p, i, e, kk

        allocate (r(derham%dimension(0)), rp(derham%dimension(0)), &
                  glam(derham_m%dimension(1)))
        do e = 1, derham_m%grad%n
            glam(e) = 0
            do kk = derham_m%grad%rowptr(e), derham_m%grad%rowptr(e + 1) - 1
                glam(e) = glam(e) + derham_m%grad%val(kk)*lam(derham_m%grad%colidx(kk))
            end do
        end do
        r = 0
        rp = 0
        do cell = 1, size(quad%point_offset)
            call spline_derham_cell_basis(derham, quad, cell, 0, v0, dofs0, s0)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                w = quad%weight(p)/abs(js(1, 1)*js(2, 2) - js(1, 2)*js(2, 1))
                call spline_derham_evaluate(derham, 0, chi_h, quad%block(p), quad%xi(:, p), &
                                            0.0_dp, chi)
                call label_inverse(chi(1), pn, dpn)
                dc = (enclosed_current(pn + h) - enclosed_current(pn - h))/(4*pi*h)*dpn
                dpdchi = -source*psi_b*dpn
                call spline_derham_evaluate(derham, 2, b2, quad%block(p), quad%xi(:, p), &
                                            0.0_dp, proxy)
                call spline_derham_evaluate(derham_m, 1, glam, quad%block(p), quad%xi(:, p), &
                                            0.0_dp, gl)
                f = dot_product(proxy, gl) + dc*dot_product(proxy(1:2), &
                                                            matmul(transpose(js), winding(s(1), s(2))))
                do i = 1, size(dofs0)
                    r(dofs0(i)) = r(dofs0(i)) + w*v0(1, i, q)*(f + dpdchi*abs(det3(j3)))
                    rp(dofs0(i)) = rp(dofs0(i)) + w*v0(1, i, q)*dpdchi*abs(det3(j3))
                end do
            end do
        end do
        ! (4) comes from the chi variation, which vanishes on the wall
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall)
        r(wall) = 0
        rp(wall) = 0
        rel = norm2(r)/norm2(rp)
    end function pressure_residual

end program test_eulerian_solovev_labels
