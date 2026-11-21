module milca_gs_spline
    !! Grad-Shafranov energy on FortFEM's multipatch spline spaces
    !! (docs/free-boundary-ale.md, section 0): the Eulerian unknown psi on
    !! block-structured splines whose blocks follow the plasma boundary and
    !! the separatrix, so that sources never cut a cell.
    !!
    !! Weak form of Delta* psi = -R f in the (R, Z) plane, Delta* psi =
    !! R div(grad psi/R): the energy
    !!
    !!     E(psi) = int [ |grad psi|^2/(2 R) - f psi ] dR dZ
    !!
    !! is stationary at the solution. Here f(R, Z) is a fixed source per
    !! quadrature point; nonlinear plasma sources come with the free-boundary
    !! solver.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_space, only: spline_space_t, spline_boundary_dofs
    use fortfem_spline_eval, only: spline_quadrature_t, spline_point_matrix, &
        spline_boundary_nodes, spline_evaluate, spline_evaluate_at, spline_locate
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_sources, gs_law_primitives
    use fortfem_spline_energy, only: pointwise_energy_t, spline_energy, &
        spline_energy_hessian_csr
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free
    implicit none
    private
    public :: gs_spline_kernel_t, gs_spline_dirichlet_solve, boundary_value_t
    public :: gs_spline_law_kernel_t, gs_spline_newton, gs_spline_extremum, &
        gs_spline_newton_axis
    integer, parameter, public :: GS_SPLINE_OK = 0, GS_SPLINE_FAILED = 1

    type, extends(pointwise_energy_t) :: gs_spline_kernel_t
        !> R and the source f at every quadrature point
        real(dp), allocatable :: r(:), f(:)
    contains
        procedure :: evaluate => gs_kernel_evaluate
    end type

    type, extends(pointwise_energy_t) :: gs_spline_law_kernel_t
        !! E = int [ |grad psi|^2/(2 mu0 R) - R p(psi) - F(psi)^2/(2 mu0 R) ]
        !! dR dZ with p, F^2/2 the primitives of the law (psi measured from
        !! the boundary): its Euler-Lagrange equation is
        !! Delta* psi = -mu0 R^2 p' - F F'.
        type(gs_profile_law_t) :: law
        real(dp) :: mu0 = 4.0e-7_dp*acos(-1.0_dp)
        real(dp), allocatable :: r(:)
    contains
        procedure :: evaluate => law_kernel_evaluate
    end type

    type, extends(pointwise_energy_t) :: axis_derivative_kernel_t
        !! Its energy gradient is u_i = d g_i/d psi_a = -int (R dp'/dpsi_a +
        !! dFF'/dpsi_a/(mu0 R)) v_i: the dependence of the law sources on the
        !! axis flux through psi_N.
        type(gs_profile_law_t) :: law
        real(dp) :: mu0 = 4.0e-7_dp*acos(-1.0_dp)
        real(dp), allocatable :: r(:)
    contains
        procedure :: evaluate => axis_kernel_evaluate
    end type

    abstract interface
        pure function boundary_value_t(r, z) result(value)
            import :: dp
            real(dp), intent(in) :: r, z
            real(dp) :: value
        end function
    end interface

contains

    subroutine gs_kernel_evaluate(self, cell, first_point, npoint, jet, weight, need, &
                                  w, gradient, hessian, status)
        class(gs_spline_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r, f
        integer :: q
        status = 0
        do q = 1, npoint
            r = self%r(first_point + q - 1)
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2)/r - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [-f, jet(2, q)/r, jet(3, q)/r]
            if (need >= 2) hessian(:, q) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1/r, &
                                            0.0_dp, 0.0_dp, 0.0_dp, 1/r]
        end do
    end subroutine gs_kernel_evaluate

    subroutine law_kernel_evaluate(self, cell, first_point, npoint, jet, weight, need, &
                                   w, gradient, hessian, status)
        class(gs_spline_law_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r, psi, pp, ff, dpp, dff, pressure, f, a
        integer :: q, st
        status = 0
        do q = 1, npoint
            r = self%r(first_point + q - 1)
            psi = jet(1, q)
            a = 1/(self%mu0*r)
            call gs_law_primitives(self%law, psi, pressure, f, st)
            if (st /= 0) then
                status = 1
                return
            end if
            w(q) = 0.5_dp*a*(jet(2, q)**2 + jet(3, q)**2) - r*pressure - 0.5_dp*a*f*f
            if (need >= 1) then
                call gs_law_sources(self%law, psi, pp, ff, dpp, dff)
                gradient(:, q) = [-(r*pp + a*ff), a*jet(2, q), a*jet(3, q)]
                if (need >= 2) hessian(:, q) = [-(r*dpp + a*dff), 0.0_dp, 0.0_dp, &
                                                0.0_dp, a, 0.0_dp, 0.0_dp, 0.0_dp, a]
            end if
        end do
    end subroutine law_kernel_evaluate

    subroutine axis_kernel_evaluate(self, cell, first_point, npoint, jet, weight, need, &
                                    w, gradient, hessian, status)
        class(axis_derivative_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r, pp, ff, dpp, dff, app, aff
        integer :: q
        status = 0
        do q = 1, npoint
            r = self%r(first_point + q - 1)
            w(q) = 0
            if (need >= 1) then
                call gs_law_sources(self%law, jet(1, q), pp, ff, dpp, dff, app, aff)
                gradient(:, q) = [-(r*app + aff/(self%mu0*r)), 0.0_dp, 0.0_dp]
            end if
            if (need >= 2) hessian(:, q) = 0
        end do
    end subroutine axis_kernel_evaluate

    subroutine gs_spline_newton_axis(space, quad, kernel, tag, sigma, coef, tolerance, &
                                     max_iterations, iterations, residual, axis, status)
        !! Fixed-boundary Newton with the normalised law in psi_N = psi/psi_a
        !! and the self-consistent axis flux psi_a = sigma-extremum of psi
        !! (psi = 0 on the tagged sides). The Jacobian is H + u v^T with
        !! u = d g/d psi_a and v the basis values at the axis (d psi_a =
        !! v . dc there, grad psi = 0), solved by Sherman-Morrison on one
        !! factorisation. kernel%law%axis follows psi_a.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(gs_spline_law_kernel_t), intent(inout) :: kernel
        integer, intent(in) :: tag, max_iterations
        real(dp), intent(in) :: sigma, tolerance
        real(dp), intent(inout) :: coef(:, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out) :: residual, axis(2)
        type(axis_derivative_kernel_t) :: dkernel
        type(sparse_direct_factor_t) :: factor
        type(csr_matrix_t) :: h
        integer, allocatable :: bdofs(:), interior(:), map(:), rowptr(:), colidx(:)
        real(dp), allocatable :: g(:, :), u(:, :), x(:), y(:), v(:), trial(:, :), val(:)
        real(dp) :: energy, norm0, norm, damping, trial_norm, psi_a, xi(2), vx, vy
        integer :: ndof, k, st, back, block

        status = GS_SPLINE_FAILED
        ndof = size(coef, 2)
        call spline_boundary_dofs(space, tag, bdofs)
        allocate (map(ndof))
        map = 0
        map(bdofs) = 1
        interior = pack([(k, k=1, ndof)], map == 0)
        map = 0
        map(interior) = [(k, k=1, size(interior))]
        coef(1, bdofs) = 0
        dkernel%order = 1
        dkernel%mu0 = kernel%mu0
        dkernel%r = kernel%r
        allocate (g(1, ndof), u(1, ndof), v(size(interior)))
        call residual_at(coef, norm, st)
        if (st /= 0) return
        norm0 = max(norm, tiny(1.0_dp))
        do iterations = 0, max_iterations
            residual = norm/norm0
            if (residual <= tolerance) then
                status = GS_SPLINE_OK
                return
            end if
            if (iterations == max_iterations) exit
            call spline_energy_hessian_csr(space, quad, kernel, coef, h, st)
            if (st /= 0) return
            call factor_submatrix(h, interior, factor, st)
            if (st /= 0) return
            dkernel%law = kernel%law
            call spline_energy(space, quad, dkernel, coef, energy, st, u)
            if (st /= 0) return
            allocate (x(size(interior)), y(size(interior)))
            call sparse_direct_solve_factored(factor, -g(1, interior), x, st)
            if (st == 0) call sparse_direct_solve_factored(factor, u(1, interior), y, st)
            call sparse_direct_free(factor)
            if (st /= 0) return
            ! v: the basis values at the axis, on the interior unknowns
            call spline_locate(space, 1, axis, block, xi, st)
            if (st /= 0) return
            call spline_point_matrix(space, [block], reshape(xi, [2, 1]), 0, rowptr, &
                                     colidx, val, st)
            if (st /= 0) return
            v = 0
            do k = rowptr(1), rowptr(2) - 1
                if (map(colidx(k)) > 0) v(map(colidx(k))) = val(k)
            end do
            vx = dot_product(v, x)
            vy = dot_product(v, y)
            x = x - y*vx/(1 + vy)
            damping = 1
            do back = 0, 10
                trial = coef
                trial(1, interior) = trial(1, interior) + damping*x
                call residual_at(trial, trial_norm, st)
                if (st == 0 .and. trial_norm < norm) exit
                damping = 0.5_dp*damping
            end do
            deallocate (x, y)
            if (.not. trial_norm < norm) return
            coef = trial
            norm = trial_norm
        end do
    contains
        subroutine residual_at(c, value, st)
            !! Axis of c, the law at its flux, and the interior gradient norm.
            real(dp), intent(in) :: c(:, :)
            real(dp), intent(out) :: value
            integer, intent(out) :: st
            call gs_spline_extremum(space, quad, c, sigma, axis, psi_a, st)
            if (st /= 0) return
            kernel%law%axis = psi_a
            call spline_energy(space, quad, kernel, c, energy, st, g)
            if (st /= 0) return
            value = norm2(g(1, interior))
        end subroutine residual_at
    end subroutine gs_spline_newton_axis

    subroutine factor_submatrix(h, rows, factor, status)
        type(csr_matrix_t), intent(in) :: h
        integer, intent(in) :: rows(:)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(out) :: status
        integer, allocatable :: map(:), colptr(:), rowind(:)
        real(dp), allocatable :: values(:)
        integer :: ni, i, k, j, nnz
        ni = size(rows)
        allocate (map(h%n))
        map = 0
        map(rows) = [(k, k=1, ni)]
        allocate (colptr(ni + 1))
        colptr(1) = 1
        nnz = 0
        do i = 1, ni
            do k = h%rowptr(rows(i)), h%rowptr(rows(i) + 1) - 1
                if (map(h%colidx(k)) > 0) nnz = nnz + 1
            end do
            colptr(i + 1) = nnz + 1
        end do
        allocate (rowind(nnz), values(nnz))
        j = 0
        do i = 1, ni
            do k = h%rowptr(rows(i)), h%rowptr(rows(i) + 1) - 1
                if (map(h%colidx(k)) == 0) cycle
                j = j + 1
                rowind(j) = map(h%colidx(k))
                values(j) = h%val(k)
            end do
        end do
        call sparse_direct_factor_csc(factor, ni, colptr, rowind, values, status)
    end subroutine factor_submatrix

    subroutine gs_spline_newton(space, quad, kernel, tag, coef, tolerance, max_iterations, &
                                iterations, residual, status)
        !! Newton for the stationary point of the kernel energy with psi = 0
        !! on the sides with the tag: H_II dc = -g_I, steps halved until the
        !! gradient norm falls (at most ten times). Converged when
        !! |g_I| <= tolerance |g_I(first iterate)|.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        integer, intent(in) :: tag, max_iterations
        real(dp), intent(inout) :: coef(:, :)
        real(dp), intent(in) :: tolerance
        integer, intent(out) :: iterations, status
        real(dp), intent(out) :: residual
        integer, allocatable :: bdofs(:), interior(:), map(:)
        real(dp), allocatable :: g(:, :), step(:), trial(:, :)
        real(dp) :: energy, norm0, norm, damping, trial_norm
        integer :: ndof, k, st, back

        status = GS_SPLINE_FAILED
        ndof = size(coef, 2)
        call spline_boundary_dofs(space, tag, bdofs)
        allocate (map(ndof))
        map = 0
        map(bdofs) = 1
        interior = pack([(k, k=1, ndof)], map == 0)
        coef(1, bdofs) = 0
        allocate (g(1, ndof))
        call spline_energy(space, quad, kernel, coef, energy, st, g)
        if (st /= 0) return
        norm0 = max(norm2(g(1, interior)), tiny(1.0_dp))
        norm = norm0
        do iterations = 0, max_iterations
            residual = norm/norm0
            if (residual <= tolerance) then
                status = GS_SPLINE_OK
                return
            end if
            if (iterations == max_iterations) exit
            call interior_newton_step(space, quad, kernel, coef, g, interior, step, st)
            if (st /= 0) return
            damping = 1
            do back = 0, 10
                trial = coef
                trial(1, interior) = trial(1, interior) + damping*step
                call spline_energy(space, quad, kernel, trial, energy, st, g)
                trial_norm = huge(1.0_dp)
                if (st == 0) trial_norm = norm2(g(1, interior))
                if (trial_norm < norm) exit
                damping = 0.5_dp*damping
            end do
            if (.not. trial_norm < norm) return
            coef = trial
            norm = trial_norm
        end do
    end subroutine gs_spline_newton

    subroutine interior_newton_step(space, quad, kernel, coef, g, interior, step, status)
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :), g(:, :)
        integer, intent(in) :: interior(:)
        real(dp), allocatable, intent(out) :: step(:)
        integer, intent(out) :: status
        type(csr_matrix_t) :: h
        call spline_energy_hessian_csr(space, quad, kernel, coef, h, status)
        if (status /= 0) return
        call solve_submatrix(h, interior, -g(1, interior), step, status)
    end subroutine interior_newton_step

    subroutine solve_submatrix(h, rows, rhs, solution, status)
        !! Solve H(rows, rows) x = rhs for symmetric H in CSR (its CSR rows are
        !! the CSC columns) with the sparse direct solver.
        type(csr_matrix_t), intent(in) :: h
        integer, intent(in) :: rows(:)
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        type(sparse_direct_factor_t) :: factor
        integer, allocatable :: map(:), colptr(:), rowind(:)
        real(dp), allocatable :: values(:)
        integer :: ni, i, k, j, nnz
        ni = size(rows)
        allocate (map(h%n))
        map = 0
        map(rows) = [(k, k=1, ni)]
        allocate (colptr(ni + 1))
        colptr(1) = 1
        nnz = 0
        do i = 1, ni
            do k = h%rowptr(rows(i)), h%rowptr(rows(i) + 1) - 1
                if (map(h%colidx(k)) > 0) nnz = nnz + 1
            end do
            colptr(i + 1) = nnz + 1
        end do
        allocate (rowind(nnz), values(nnz))
        j = 0
        do i = 1, ni
            do k = h%rowptr(rows(i)), h%rowptr(rows(i) + 1) - 1
                if (map(h%colidx(k)) == 0) cycle
                j = j + 1
                rowind(j) = map(h%colidx(k))
                values(j) = h%val(k)
            end do
        end do
        allocate (solution(ni))
        call sparse_direct_factor_csc(factor, ni, colptr, rowind, values, status)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, rhs, solution, status)
        call sparse_direct_free(factor)
    end subroutine solve_submatrix

    subroutine gs_spline_extremum(space, quad, coef, sigma, point, value, status)
        !! Location and value of the extremum of sigma*psi (sigma = +1 a
        !! maximum): the best quadrature point, then Newton on grad psi = 0
        !! with the exact spline Hessian.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: coef(:, :), sigma
        real(dp), intent(out) :: point(2), value
        integer, intent(out) :: status
        real(dp), allocatable :: jet(:, :)
        real(dp) :: xi(2), j6(6), h(2, 2), dx(2), det
        integer :: q, block, it, st
        status = GS_SPLINE_FAILED
        allocate (jet(1, quad%point_count))
        call spline_evaluate(space, quad, coef, 0, jet, st)
        if (st /= 0) return
        q = maxloc(sigma*jet(1, :), dim=1)
        point = quad%s(:, q)
        block = quad%block(q)
        do it = 1, 30
            call spline_locate(space, block, point, block, xi, st)
            if (st /= 0) return
            call spline_evaluate_at(space, coef, block, xi, 2, j6, st)
            if (st /= 0) return
            h = reshape([j6(4), j6(5), j6(5), j6(6)], [2, 2])
            det = h(1, 1)*h(2, 2) - h(1, 2)**2
            if (abs(det) <= tiny(1.0_dp)) return
            dx = -[h(2, 2)*j6(2) - h(1, 2)*j6(3), -h(1, 2)*j6(2) + h(1, 1)*j6(3)]/det
            point = point + dx
            if (norm2(dx) <= 1.0e-13_dp*max(1.0_dp, norm2(point))) exit
        end do
        call spline_locate(space, block, point, block, xi, st)
        if (st /= 0) return
        call spline_evaluate_at(space, coef, block, xi, 0, j6(1:1), st)
        if (st /= 0) return
        value = j6(1)
        status = GS_SPLINE_OK
    end subroutine gs_spline_extremum

    subroutine gs_spline_dirichlet_solve(space, quad, kernel, tag, boundary_value, &
                                         coef, status)
        !! Minimiser of the energy with psi = boundary_value on the sides with
        !! the tag: the boundary coefficients fit the data in least squares
        !! at the boundary nodes (four panels per cell side), the interior
        !! ones solve H_II c_I = -g_I with the sparse direct solver.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(gs_spline_kernel_t), intent(in) :: kernel
        integer, intent(in) :: tag
        procedure(boundary_value_t) :: boundary_value
        real(dp), intent(out) :: coef(:, :)
        integer, intent(out) :: status
        integer, allocatable :: bdofs(:), blk(:), rowptr(:), colidx(:), map(:), &
            colptr(:), rowind(:)
        real(dp), allocatable :: xi(:, :), s(:, :), val(:), a(:, :), b(:), x(:), &
            g(:, :), rhs(:), sol(:), values(:)
        type(csr_matrix_t) :: h
        type(sparse_direct_factor_t) :: factor
        real(dp) :: energy
        integer :: ndof, nb, npt, i, k, j, ni, st, info, nnz
        integer, allocatable :: interior(:)

        status = GS_SPLINE_FAILED
        ndof = size(coef, 2)
        coef = 0
        call spline_boundary_dofs(space, tag, bdofs)
        nb = size(bdofs)
        call spline_boundary_nodes(space, tag, 4, blk, xi, s, st)
        if (st /= 0) return
        call spline_point_matrix(space, blk, xi, 0, rowptr, colidx, val, st)
        if (st /= 0) return
        npt = size(blk)
        allocate (map(ndof))
        map = 0
        map(bdofs) = [(k, k=1, nb)]
        allocate (a(npt, nb), b(npt))
        a = 0
        do i = 1, npt
            do k = rowptr(i), rowptr(i + 1) - 1
                if (map(colidx(k)) > 0) a(i, map(colidx(k))) = val(k)
            end do
            b(i) = boundary_value(s(1, i), s(2, i))
        end do
        call least_squares(a, b, x, info)
        if (info /= 0) return
        coef(1, bdofs) = x
        ! interior: H_II c_I = -g_I at c with the boundary values
        allocate (g(1, ndof))
        call spline_energy(space, quad, kernel, coef, energy, st, g)
        if (st /= 0) return
        call spline_energy_hessian_csr(space, quad, kernel, coef, h, st)
        if (st /= 0) return
        interior = pack([(k, k=1, ndof)], map == 0)
        ni = size(interior)
        map = 0
        map(interior) = [(k, k=1, ni)]
        ! H is symmetric, so the CSR rows of H_II are its CSC columns
        allocate (colptr(ni + 1))
        colptr(1) = 1
        nnz = 0
        do i = 1, ni
            do k = h%rowptr(interior(i)), h%rowptr(interior(i) + 1) - 1
                if (map(h%colidx(k)) > 0) nnz = nnz + 1
            end do
            colptr(i + 1) = nnz + 1
        end do
        allocate (rowind(nnz), values(nnz))
        j = 0
        do i = 1, ni
            do k = h%rowptr(interior(i)), h%rowptr(interior(i) + 1) - 1
                if (map(h%colidx(k)) == 0) cycle
                j = j + 1
                rowind(j) = map(h%colidx(k))
                values(j) = h%val(k)
            end do
        end do
        call sparse_direct_factor_csc(factor, ni, colptr, rowind, values, st)
        if (st /= 0) return
        rhs = -g(1, interior)
        allocate (sol(ni))
        call sparse_direct_solve_factored(factor, rhs, sol, st)
        call sparse_direct_free(factor)
        if (st /= 0) return
        coef(1, interior) = coef(1, interior) + sol
        status = GS_SPLINE_OK
    end subroutine gs_spline_dirichlet_solve

    subroutine least_squares(a, b, x, info)
        !! min |a x - b| by the normal equations with Cholesky (a has full
        !! column rank and is well conditioned for B-spline traces).
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: n(:, :)
        external :: dposv
        n = matmul(transpose(a), a)
        x = matmul(transpose(a), b)
        call dposv('U', size(n, 1), 1, n, size(n, 1), x, size(x), info)
    end subroutine least_squares
end module milca_gs_spline
