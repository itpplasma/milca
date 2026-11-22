module milca_spline_ale
    !! ALE mesh of a C0 spline section (docs/equilibrium-roadmap.md): the mesh
    !! moves by its own equation, never as a physical unknown.
    !!
    !! - spline_winslow: the block maps minimise the Winslow energy
    !!   int |dX/dxi|_F^2 / det(dX/dxi) dxi (Dirichlet energy of the inverse
    !!   map, generated FortSym leaves of gen_winslow_seed) over the control
    !!   points, with fixed control points (the wall) and pinned ones moved to
    !!   targets (a mesh centre attracted to a magnetic axis). Newton with the
    !!   exact Hessian; steps keep det > 0 at every quadrature point.
    !! - spline_transfer_v0: a V0 field of one mesh interpolated on another
    !!   mesh of the same section (Greville interpolation of the old field at
    !!   the physical points of the new mesh).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_space, only: spline_space_t, spline_cell_broken
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_cell_xi_basis, spline_geometry_at, spline_locate, spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_project, &
        spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_winslow_seed_leaves, only: winslow_density, winslow_gradient, winslow_hessian
    implicit none
    private
    public :: spline_geometry_dofs, spline_set_geometry, spline_point_dof, spline_winslow, &
              spline_transfer_v0, spline_min_jacobian, spline_blend_weights, spline_cell_at, &
              spline_block_dofs

    ! state of spline_transfer_v0's projection callback
    type(spline_derham_t) :: old_derham
    type(spline_space_t) :: new_space
    real(dp), allocatable :: old_coef(:)
    integer :: hint = 1, transfer_status = 0

contains

    subroutine spline_geometry_dofs(space, x, status)
        !! Control points per global function, x(2, dof_count), of a C0 space
        !! (every broken function is one global function with weight 1).
        type(spline_space_t), intent(in) :: space
        real(dp), allocatable, intent(out) :: x(:, :)
        integer, intent(out) :: status
        integer :: b, p
        status = 1
        allocate (x(2, space%dof_count), source=0.0_dp)
        do b = 1, space%broken_count
            do p = space%broken_ptr(b), space%broken_ptr(b + 1) - 1
                if (abs(space%broken_weight(p) - 1) > 1e-12_dp) return
                x(:, space%broken_dof(p)) = space%geometry(:, b)
            end do
        end do
        status = 0
    end subroutine spline_geometry_dofs

    subroutine spline_set_geometry(space, x)
        type(spline_space_t), intent(inout) :: space
        real(dp), intent(in) :: x(:, :)
        integer :: b, p
        do b = 1, space%broken_count
            space%geometry(:, b) = 0
            do p = space%broken_ptr(b), space%broken_ptr(b + 1) - 1
                space%geometry(:, b) = space%geometry(:, b) + &
                                       space%broken_weight(p)*x(:, space%broken_dof(p))
            end do
        end do
    end subroutine spline_set_geometry

    integer function spline_point_dof(space, block, xi) result(dof)
        !! The global function whose control point is the mesh vertex at
        !! (block, xi) (a cell corner, where only that function is nonzero);
        !! 0 when (block, xi) is not a vertex.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        real(dp) :: jet(2, 6), x(2, 1)
        real(dp), allocatable :: g(:, :)
        integer :: st, j
        dof = 0
        call spline_geometry_at(space, block, xi, jet, st)
        if (st /= 0) return
        call spline_geometry_dofs(space, g, st)
        if (st /= 0) return
        x(:, 1) = jet(:, 1)
        j = minloc(sum((g - spread(x(:, 1), 2, size(g, 2)))**2, 1), 1)
        if (norm2(g(:, j) - x(:, 1)) < 1e-10_dp*(1 + norm2(x(:, 1)))) dof = j
    end function spline_point_dof

    subroutine spline_winslow(space, fixed, pin_block, pin_xi, target, tolerance, max_steps, &
                              steps, status)
        !! Winslow mesh of the section: control points `fixed` (global
        !! functions, the wall) stay, the map point at (pin_block, pin_xi) is
        !! moved to `target` (two linear constraints with multipliers), the
        !! other control points minimise the Winslow energy. Status 0
        !! converged, 1 no valid step, 2 not converged, 3 linear solve,
        !! 4 non-C0 space.
        type(spline_space_t), intent(inout) :: space
        integer, intent(in) :: fixed(:), pin_block, max_steps
        real(dp), intent(in) :: pin_xi(2), target(2), tolerance
        integer, intent(out) :: steps, status
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: x(:, :), xt(:, :), delta(:), grad(:), xj(:, :, :, :), &
                                 wref(:), cval(:), r(:)
        integer, allocatable :: dofs(:, :), clamp(:), rowptr(:), ccol(:)
        type(csr_t) :: kkt, hs
        real(dp) :: shift, dmax
        real(dp) :: energy, trial, alpha, scale, mindet, lam(2), viol(2)
        integer :: n, nloc, nq, cell, q, a, st, k
        logical :: moving, debug
        character(8) :: env

        call get_environment_variable('MILCA_WINSLOW_DEBUG', env)
        debug = len_trim(env) > 0
        call spline_geometry_dofs(space, x, st)
        status = 4
        if (st /= 0) return
        n = space%dof_count
        nloc = (space%degree + 1)**2
        if (pin_block > 0) then
            call spline_point_matrix(space, [pin_block], reshape(pin_xi, [2, 1]), 0, rowptr, &
                                     ccol, cval, st)
            status = 1
            if (st /= 0) return
        else
            ! no pin: a zero constraint row, its multipliers clamped
            ccol = [1]
            cval = [0.0_dp]
        end if
        call spline_quadrature_build(space, space%degree + 2, quad, st)
        if (st /= 0) return
        nq = quad%points_per_cell
        allocate (xj(nloc, 6, nq, space%cell_count), dofs(nloc, space%cell_count), &
                  wref(quad%point_count))
        do cell = 1, space%cell_count
            call spline_cell_xi_basis(space, quad, cell, xj(:, :, :, cell))
            call spline_cell_broken(space, cell, dofs(:, cell))
            do a = 1, nloc
                dofs(a, cell) = space%broken_dof(space%broken_ptr(dofs(a, cell)))
            end do
            do q = 1, nq
                k = quad%point_offset(cell) + q - 1
                wref(k) = quad%weight(k)/abs(cell_det(x, cell, q))
            end do
        end do
        scale = maxval(abs(x))
        clamp = [fixed, n + fixed]
        if (pin_block <= 0) clamp = [clamp, 2*n + 1, 2*n + 2]
        lam = 0
        steps = 0
        status = 2
        energy = 0
        do while (steps <= max_steps)
            call assemble(x, energy, grad, kkt, mindet)
            viol = [dot_product(cval, x(1, ccol)), dot_product(cval, x(2, ccol))] - target
            if (pin_block <= 0) viol = 0
            moving = norm2(viol) > tolerance*scale
            r = [grad, -viol]
            r(1:2*n) = r(1:2*n) + [lam(1)*scatter(cval), lam(2)*scatter(cval)]
            r(1:2*n) = -r(1:2*n)
            r(clamp) = 0
            ! Newton, Levenberg-shifted until the step descends (the Winslow
            ! energy is polyconvex, not convex: near folds H is indefinite)
            shift = 0
            dmax = maxval(abs(kkt%val))
            allocate (delta(2*n + 2))
            do k = 1, 20
                hs = kkt
                call shift_diagonal(hs, 2*n, shift)
                call eu_constrained_solve(hs, r, clamp, 0*real(clamp, dp), delta, st)
                if (st /= 0) then
                    status = 3
                    return
                end if
                if (moving .or. dot_product(grad, delta(1:2*n)) < 0) exit
                shift = max(10*shift, 1e-6_dp*dmax)
            end do
            if (.not. moving .and. maxval(abs(delta(1:2*n))) < tolerance*scale) then
                status = 0
                deallocate (delta)
                exit
            end if
            if (steps == max_steps) exit
            if (debug) then
                block
                    real(dp) :: ep, em, md2, dd(2, n)
                    dd = transpose(reshape(delta(1:2*n), [n, 2]))
                    call energy_only(x + 1e-6_dp*dd, ep, md2)
                    call energy_only(x - 1e-6_dp*dd, em, md2)
                    print '(a,2es14.6)', '  dE: fd, grad.d', (ep - em)/2e-6_dp, &
                        dot_product(grad, delta(1:2*n))
                end block
            end if
            alpha = 1
            do k = 1, 30
                xt = transpose(reshape(delta(1:2*n), [n, 2]))*alpha + x
                call energy_only(xt, trial, mindet)
                if (mindet > 0 .and. (moving .or. trial <= energy + 1e-14_dp*abs(energy))) exit
                alpha = alpha/2
            end do
            if (.not. mindet > 0) then
                status = 1
                return
            end if
            if (debug) print '(a,i3,a,es14.6,a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '  winslow', steps, &
                ' E', energy, ' viol', norm2(viol), ' alpha', alpha, ' det', mindet, &
                ' |dx|', maxval(abs(delta(1:2*n)))
            x = xt
            lam = lam + alpha*delta(2*n + 1:2*n + 2)
            steps = steps + 1
            deallocate (delta)
        end do
        call spline_set_geometry(space, x)

    contains

        subroutine shift_diagonal(h, nd, sh)
            type(csr_t), intent(inout) :: h
            integer, intent(in) :: nd
            real(dp), intent(in) :: sh
            integer :: i, p
            do i = 1, nd
                do p = h%row_ptr(i), h%row_ptr(i + 1) - 1
                    if (h%col(p) == i) h%val(p) = h%val(p) + sh
                end do
            end do
        end subroutine shift_diagonal

        function scatter(c) result(v)
            real(dp), intent(in) :: c(:)
            real(dp) :: v(n)
            v = 0
            v(ccol) = c
        end function scatter

        real(dp) function cell_det(xc, c, qq) result(d)
            real(dp), intent(in) :: xc(:, :)
            integer, intent(in) :: c, qq
            real(dp) :: am(2, 2)
            am(:, 1) = matmul(xc(:, dofs(:, c)), xj(:, 2, qq, c))
            am(:, 2) = matmul(xc(:, dofs(:, c)), xj(:, 3, qq, c))
            d = am(1, 1)*am(2, 2) - am(1, 2)*am(2, 1)
        end function cell_det

        subroutine energy_only(xc, e, md)
            real(dp), intent(in) :: xc(:, :)
            real(dp), intent(out) :: e, md
            real(dp) :: am(2, 2), w, det
            integer :: c, qq
            e = 0
            md = huge(1.0_dp)
            do c = 1, space%cell_count
                do qq = 1, nq
                    am(:, 1) = matmul(xc(:, dofs(:, c)), xj(:, 2, qq, c))
                    am(:, 2) = matmul(xc(:, dofs(:, c)), xj(:, 3, qq, c))
                    det = am(1, 1)*am(2, 2) - am(1, 2)*am(2, 1)
                    md = min(md, det)
                    if (det <= 0) cycle
                    call winslow_density(am(1, 1), am(2, 1), am(1, 2), am(2, 2), 0.0_dp, w, det)
                    e = e + wref(quad%point_offset(c) + qq - 1)*w
                end do
            end do
        end subroutine energy_only

        subroutine assemble(xc, e, gr, hm, md)
            !! Winslow gradient and the KKT matrix [H C^T; C 0] (symmetric).
            real(dp), intent(in) :: xc(:, :)
            real(dp), intent(out) :: e, md
            real(dp), allocatable, intent(out) :: gr(:)
            type(csr_t), intent(out) :: hm
            real(dp) :: am(2, 2), w, det, g4(4), h(10), h4(4, 4), du(4, 2*nloc), wq
            integer, allocatable :: rows(:), cols(:), idx(:)
            real(dp), allocatable :: vals(:)
            type(csc_t) :: m
            type(fortsparse_status_t) :: info
            integer :: c, qq, i, j, nz, kk, cap
            allocate (gr(2*n), source=0.0_dp)
            cap = space%cell_count*(2*nloc)**2 + 4*size(ccol) + 2
            allocate (rows(cap), cols(cap), vals(cap), idx(2*nloc))
            e = 0
            md = huge(1.0_dp)
            nz = 0
            do c = 1, space%cell_count
                idx = [dofs(:, c), n + dofs(:, c)]
                kk = nz
                do i = 1, 2*nloc
                    do j = 1, 2*nloc
                        kk = kk + 1
                        rows(kk) = idx(i)
                        cols(kk) = idx(j)
                        vals(kk) = 0
                    end do
                end do
                do qq = 1, nq
                    wq = wref(quad%point_offset(c) + qq - 1)
                    am(:, 1) = matmul(xc(:, dofs(:, c)), xj(:, 2, qq, c))
                    am(:, 2) = matmul(xc(:, dofs(:, c)), xj(:, 3, qq, c))
                    call winslow_density(am(1, 1), am(2, 1), am(1, 2), am(2, 2), 0.0_dp, w, det)
                    md = min(md, det)
                    e = e + wq*w
                    call winslow_gradient(am(1, 1), am(2, 1), am(1, 2), am(2, 2), 0.0_dp, &
                                          g4(1), g4(2), g4(3), g4(4))
                    call winslow_hessian(am(1, 1), am(2, 1), am(1, 2), am(2, 2), 0.0_dp, &
                                         h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10))
                    h4 = reshape([h(1), h(2), h(4), h(7), h(2), h(3), h(5), h(8), &
                                  h(4), h(5), h(6), h(9), h(7), h(8), h(9), h(10)], [4, 4])
                    ! (a11, a21, a12, a22) = (x_1, y_1, x_2, y_2)
                    du = 0
                    du(1, 1:nloc) = xj(:, 2, qq, c)
                    du(3, 1:nloc) = xj(:, 3, qq, c)
                    du(2, nloc + 1:) = xj(:, 2, qq, c)
                    du(4, nloc + 1:) = xj(:, 3, qq, c)
                    do i = 1, 2*nloc
                        gr(idx(i)) = gr(idx(i)) + wq*dot_product(g4, du(:, i))
                    end do
                    kk = nz
                    do i = 1, 2*nloc
                        do j = 1, 2*nloc
                            kk = kk + 1
                            vals(kk) = vals(kk) + wq*dot_product(du(:, i), matmul(h4, du(:, j)))
                        end do
                    end do
                end do
                nz = kk
            end do
            do i = 1, size(ccol)
                rows(nz + 1:nz + 4) = [ccol(i), 2*n + 1, n + ccol(i), 2*n + 2]
                cols(nz + 1:nz + 4) = [2*n + 1, ccol(i), 2*n + 2, n + ccol(i)]
                vals(nz + 1:nz + 4) = cval(i)
                nz = nz + 4
            end do
            rows(nz + 1:nz + 2) = [2*n + 1, 2*n + 2]
            cols(nz + 1:nz + 2) = [2*n + 1, 2*n + 2]
            vals(nz + 1:nz + 2) = 0
            nz = nz + 2
            call csc_from_triplet(2*n + 2, 2*n + 2, cols(:nz), rows(:nz), vals(:nz), m, info)
            if (info%code /= 0) error stop 'Winslow KKT triplets'
            hm%nrow = 2*n + 2
            hm%ncol = 2*n + 2
            hm%row_ptr = m%col_ptr
            hm%col = m%row_idx
            hm%val = m%val
        end subroutine assemble

    end subroutine spline_winslow

    subroutine spline_blend_weights(space, inner_cells, wall, w, status)
        !! Weights per control point for a translation that follows an
        !! interior point: w = 1 on the control points of inner_cells, 0 on
        !! `wall`, harmonic in the block parameters between (reference
        !! Laplace, int grad_xi w . grad_xi v dxi = 0).
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: inner_cells(:), wall(:)
        real(dp), allocatable, intent(out) :: w(:)
        integer, intent(out) :: status
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: xj(:, :, :), vals(:), rhs(:), val(:)
        integer, allocatable :: br(:), rows(:), cols(:), inner(:), fixed(:)
        type(csc_t) :: m
        type(csr_t) :: k
        type(fortsparse_status_t) :: info
        integer :: n, nloc, cell, q, a, b, nz, i
        logical, allocatable :: is_inner(:)
        real(dp) :: t
        n = space%dof_count
        nloc = (space%degree + 1)**2
        call spline_quadrature_build(space, space%degree + 1, quad, status)
        if (status /= 0) return
        allocate (xj(nloc, 6, quad%points_per_cell), br(nloc), is_inner(n))
        is_inner = .false.
        allocate (rows(space%cell_count*nloc**2), cols(space%cell_count*nloc**2), &
                  vals(space%cell_count*nloc**2))
        nz = 0
        do cell = 1, space%cell_count
            call spline_cell_xi_basis(space, quad, cell, xj)
            call spline_cell_broken(space, cell, br)
            do a = 1, nloc
                br(a) = space%broken_dof(space%broken_ptr(br(a)))
            end do
            if (any(inner_cells == cell)) is_inner(br) = .true.
            do a = 1, nloc
                do b = 1, nloc
                    t = 0
                    do q = 1, quad%points_per_cell
                        t = t + quad%t_weight(1 + mod(q - 1, quad%points_per_direction))* &
                            quad%t_weight(1 + (q - 1)/quad%points_per_direction)* &
                            (xj(a, 2, q)*xj(b, 2, q) + xj(a, 3, q)*xj(b, 3, q))
                    end do
                    nz = nz + 1
                    rows(nz) = br(a)
                    cols(nz) = br(b)
                    vals(nz) = t
                end do
            end do
        end do
        call csc_from_triplet(n, n, cols(:nz), rows(:nz), vals(:nz), m, info)
        k%nrow = n
        k%ncol = n
        k%row_ptr = m%col_ptr
        k%col = m%row_idx
        k%val = m%val
        inner = pack([(i, i=1, n)], is_inner)
        fixed = [inner, wall]
        val = [(1.0_dp, i=1, size(inner)), (0.0_dp, i=1, size(wall))]
        allocate (w(n), rhs(n), source=0.0_dp)
        call eu_constrained_solve(k, rhs, fixed, val, w, status)
    end subroutine spline_blend_weights

    subroutine spline_block_dofs(space, block, cells, dofs)
        !! The cells of `block` and the global functions supported on them.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block
        integer, allocatable, intent(out) :: cells(:), dofs(:)
        integer, allocatable :: br(:)
        logical, allocatable :: on(:)
        integer :: c, a, i
        allocate (br((space%degree + 1)**2), on(space%dof_count))
        on = .false.
        cells = pack([(c, c=1, space%cell_count)], space%cell_block == block)
        do i = 1, size(cells)
            call spline_cell_broken(space, cells(i), br)
            do a = 1, size(br)
                on(space%broken_dof(space%broken_ptr(br(a)))) = .true.
            end do
        end do
        dofs = pack([(i, i=1, space%dof_count)], on)
    end subroutine spline_block_dofs

    subroutine spline_cell_at(space, block, xi, cell, centre)
        !! Cell of the point (block, xi) and the parameters of its centre.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        integer, intent(out) :: cell
        real(dp), intent(out) :: centre(2)
        integer :: e(2), c
        associate (b1 => space%blocks(block)%breaks1, b2 => space%blocks(block)%breaks2)
            e(1) = max(1, min(size(b1) - 1, count(b1(2:size(b1) - 1) <= xi(1)) + 1))
            e(2) = max(1, min(size(b2) - 1, count(b2(2:size(b2) - 1) <= xi(2)) + 1))
            centre = [0.5_dp*(b1(e(1)) + b1(e(1) + 1)), 0.5_dp*(b2(e(2)) + b2(e(2) + 1))]
        end associate
        cell = 0
        do c = 1, space%cell_count
            if (space%cell_block(c) == block .and. all(space%cell_element(:, c) == e)) then
                cell = c
                return
            end if
        end do
    end subroutine spline_cell_at

    real(dp) function spline_min_jacobian(space) result(md)
        !! Smallest det(dX/dxi) over the Gauss points (degree + 2 per direction).
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: xj(:, :, :), x(:, :)
        integer, allocatable :: br(:)
        real(dp) :: am(2, 2)
        integer :: cell, q, nloc, st
        nloc = (space%degree + 1)**2
        call spline_quadrature_build(space, space%degree + 2, quad, st)
        allocate (xj(nloc, 6, quad%points_per_cell), br(nloc))
        md = huge(1.0_dp)
        do cell = 1, space%cell_count
            call spline_cell_xi_basis(space, quad, cell, xj)
            call spline_cell_broken(space, cell, br)
            x = space%geometry(:, br)
            do q = 1, quad%points_per_cell
                am(:, 1) = matmul(x, xj(:, 2, q))
                am(:, 2) = matmul(x, xj(:, 3, q))
                md = min(md, am(1, 1)*am(2, 2) - am(1, 2)*am(2, 1))
            end do
        end do
    end function spline_min_jacobian

    subroutine spline_transfer_v0(from, coef_from, to, coef_to, status)
        !! coef_to: V0 interpolant on `to` of the V0 field coef_from of `from`
        !! (same physical section). Status 1 when a point of `to` lies
        !! outside `from`.
        type(spline_derham_t), intent(in) :: from, to
        real(dp), intent(in) :: coef_from(:)
        real(dp), intent(out) :: coef_to(:)
        integer, intent(out) :: status
        old_derham = from
        new_space = to%space
        old_coef = coef_from
        hint = 1
        transfer_status = 0
        call spline_derham_project(to, 0, old_value, coef_to, status, to%degree + 3)
        if (status == 0) status = transfer_status
    end subroutine spline_transfer_v0

    subroutine old_value(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6), oxi(2)
        integer :: ob, st
        call spline_geometry_at(new_space, block, xi, jet, st)
        call spline_locate(old_derham%space, hint, jet(:, 1), ob, oxi, st)
        if (st /= 0) then
            transfer_status = 1
            values = 0
            return
        end if
        hint = ob
        call spline_derham_evaluate(old_derham, 0, old_coef, ob, oxi, t, values)
    end subroutine old_value

end module milca_spline_ale
