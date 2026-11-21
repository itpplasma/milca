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
        spline_boundary_nodes, spline_side_quadrature_t, spline_side_quadrature, spline_evaluate, spline_evaluate_at, spline_locate
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_sources, gs_law_primitives
    use milca_grad_shafranov_leaves, only: gs_spline_energy, gs_source, gs_l2_projection, &
        gs_ray_coarea_weights
    use fortfem_spline_energy, only: pointwise_energy_t, spline_energy, &
        spline_energy_hessian_csr
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
        gs_bem_steklov_symmetric
    implicit none
    private
    public :: gs_spline_kernel_t, gs_spline_dirichlet_solve, boundary_value_t
    public :: gs_spline_law_kernel_t, gs_spline_newton, gs_spline_extremum, &
        gs_spline_newton_axis
    public :: gs_spline_vacuum_t, gs_spline_vacuum_setup, gs_spline_free_linear
    public :: gs_spline_newton_free, gs_spline_project, gs_spline_level_integrals, &
        ray_value
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
        !> free boundary: sources only where plasma(point), at psi - psi_b
        logical, allocatable :: plasma(:)
        real(dp) :: psi_b = 0
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
        logical, allocatable :: plasma(:)
        real(dp) :: psi_b = 0
        !> .true.: the gradient is D_i = int (R dp'/dpsi + dFF'/dpsi/(mu0 R)) v_i
        logical :: flux_derivative = .false.
    contains
        procedure :: evaluate => axis_kernel_evaluate
    end type

    type, extends(pointwise_energy_t) :: projection_kernel_t
        !! (1/2) (u - v)^2 with v sampled at every quadrature point
        real(dp), allocatable :: v(:)
    contains
        procedure :: evaluate => projection_kernel_evaluate
    end type

    type :: gs_spline_vacuum_t
        !! Exterior of the wall: psi = psi_c + psi_p, psi_c the coil flux
        !! (Delta* psi_c = 0 inside), psi_p decaying outside. The wall terms
        !! of the energy are
        !!     (1/2) (T psi - psi_c)^T S (T psi - psi_c)/mu0
        !!         - oint (d_n psi_c/R) psi dl/mu0
        !! with the trace T at the quadratic panel nodes of the wall and the
        !! symmetric discrete Steklov-Poincare matrix S of FortFEM (operator
        !! -div(grad/R)); the second term cancels the coil part of the
        !! domain flux, so psi = psi_c is stationary without plasma.
        real(dp), allocatable :: nodes(:, :) !< (2, nb) wall nodes
        integer, allocatable :: dofs(:) !< columns of T (trace dofs)
        real(dp), allocatable :: trace(:, :) !< T on dofs, (nb, ntrace)
        real(dp), allocatable :: steklov(:, :) !< S, (nb, nb)
        real(dp), allocatable :: block(:, :) !< T^T S T/mu0 on dofs
        type(spline_side_quadrature_t) :: side !< Gauss points on the wall
        real(dp), allocatable :: normal(:, :) !< outward unit normals there
        integer, allocatable :: side_rowptr(:), side_colidx(:)
        real(dp), allocatable :: side_val(:) !< basis values there (CSR)
        real(dp) :: mu0 = 1.0_dp
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
        !! w = |grad psi|^2/(2 R) - f psi: the generated energy leaf with
        !! mu0 = 1, pressure f psi/R and pprime f/R.
        class(gs_spline_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r, f, g_psi, g_r, g_z, h_psi, h_grad
        integer :: q
        associate (unused_cell => cell, unused_weight => weight)
        end associate
        status = 0
        do q = 1, npoint
            r = self%r(first_point + q - 1)
            f = self%f(first_point + q - 1)
            call gs_spline_energy(r, 1.0_dp, jet(2, q), jet(3, q), f*jet(1, q)/r, 0.0_dp, &
                                  f/r, 0.0_dp, 0.0_dp, 0.0_dp, w(q), g_psi, g_r, g_z, &
                                  h_psi, h_grad)
            if (need >= 1) gradient(:, q) = [g_psi, g_r, g_z]
            if (need >= 2) hessian(:, q) = [h_psi, 0.0_dp, 0.0_dp, 0.0_dp, h_grad, 0.0_dp, &
                                            0.0_dp, 0.0_dp, h_grad]
        end do
    end subroutine gs_kernel_evaluate

    subroutine law_kernel_evaluate(self, cell, first_point, npoint, jet, weight, need, &
                                   w, gradient, hessian, status)
        !! Generated energy leaf (gs_spline_energy): the law at psi - psi_b on
        !! plasma points, no sources on vacuum points.
        class(gs_spline_law_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r, psi, pp, ff, dpp, dff, pressure, f, g_psi, g_r, g_z, h_psi, h_grad
        integer :: q, st
        logical :: plasma
        associate (unused_cell => cell, unused_weight => weight)
        end associate
        status = 0
        do q = 1, npoint
            r = self%r(first_point + q - 1)
            psi = jet(1, q) - self%psi_b
            plasma = .true.
            if (allocated(self%plasma)) plasma = self%plasma(first_point + q - 1)
            pressure = 0; f = 0; pp = 0; ff = 0; dpp = 0; dff = 0
            if (plasma) then
                call gs_law_primitives(self%law, psi, pressure, f, st)
                if (st /= 0) then
                    status = 1
                    return
                end if
                if (need >= 1) call gs_law_sources(self%law, psi, pp, ff, dpp, dff)
            end if
            call gs_spline_energy(r, self%mu0, jet(2, q), jet(3, q), pressure, f, pp, ff, &
                                  dpp, dff, w(q), g_psi, g_r, g_z, h_psi, h_grad)
            if (need >= 1) gradient(:, q) = [g_psi, g_r, g_z]
            if (need >= 2) hessian(:, q) = [h_psi, 0.0_dp, 0.0_dp, 0.0_dp, h_grad, 0.0_dp, &
                                            0.0_dp, 0.0_dp, h_grad]
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
                gradient(:, q) = 0
                if (allocated(self%plasma)) then
                    if (.not. self%plasma(first_point + q - 1)) cycle
                end if
                call gs_law_sources(self%law, jet(1, q) - self%psi_b, pp, ff, dpp, dff, &
                                    app, aff)
                ! the generated source leaf is mu0 (R pprime + ffprime/(mu0 R))
                if (self%flux_derivative) then
                    call gs_source(r, self%mu0, dpp, dff, gradient(1, q))
                    gradient(1, q) = gradient(1, q)/self%mu0
                else
                    call gs_source(r, self%mu0, app, aff, gradient(1, q))
                    gradient(1, q) = -gradient(1, q)/self%mu0
                end if
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

    subroutine gs_spline_vacuum_setup(space, tag, panels_per_cell, mu0, vac, status)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: tag, panels_per_cell
        real(dp), intent(in) :: mu0
        type(gs_spline_vacuum_t), intent(out) :: vac
        integer, intent(out) :: status
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        integer, allocatable :: blk(:), rowptr(:), colidx(:), map(:)
        real(dp), allocatable :: xi(:, :), val(:)
        integer :: i, k, nb
        logical, allocatable :: used(:)
        status = GS_SPLINE_FAILED
        vac%mu0 = mu0
        call spline_boundary_nodes(space, tag, panels_per_cell, blk, xi, vac%nodes, status)
        if (status /= 0) return
        call spline_point_matrix(space, blk, xi, 0, rowptr, colidx, val, status)
        if (status /= 0) return
        nb = size(blk)
        allocate (used(space%dof_count))
        used = .false.
        used(colidx) = .true.
        vac%dofs = pack([(k, k=1, space%dof_count)], used)
        allocate (map(space%dof_count))
        map = 0
        map(vac%dofs) = [(k, k=1, size(vac%dofs))]
        allocate (vac%trace(nb, size(vac%dofs)))
        vac%trace = 0
        do i = 1, nb
            do k = rowptr(i), rowptr(i + 1) - 1
                vac%trace(i, map(colidx(k))) = val(k)
            end do
        end do
        call gs_bem_curve_from_nodes(vac%nodes, curve, status)
        if (status /= 0) return
        call gs_bem_assemble(curve, 1, ops, status)
        if (status /= 0) return
        call gs_bem_steklov_symmetric(ops, vac%steklov, status)
        if (status /= 0) return
        vac%block = matmul(transpose(vac%trace), matmul(vac%steklov, vac%trace))/mu0
        call spline_side_quadrature(space, tag, space_degree_points(), vac%side, status)
        if (status /= 0) return
        ! counter-clockwise tangent (t_R, t_Z): the outward normal is (t_Z, -t_R)
        vac%normal = reshape([vac%side%tangent(2, :), -vac%side%tangent(1, :)], &
                             [vac%side%point_count, 2])
        vac%normal = transpose(vac%normal)
        call spline_point_matrix(space, vac%side%block, vac%side%xi, 0, vac%side_rowptr, &
                                 vac%side_colidx, vac%side_val, status)
        if (status /= 0) return
        status = GS_SPLINE_OK
    contains
        integer function space_degree_points() result(np)
            np = 8
        end function
    end subroutine gs_spline_vacuum_setup

    subroutine gs_spline_free_linear(space, quad, kernel, vac, psi_c, dn_psi_c, coef, &
                                     status)
        !! Minimiser of the domain energy of a quadratic kernel plus the wall
        !! terms with the coil flux psi_c at the wall nodes and its outward
        !! normal derivative dn_psi_c at vac%side: one Newton step from zero,
        !! (H + T^T S T/mu0) c = -g(0) + T^T S psi_c/mu0
        !!                       + oint (d_n psi_c/R) v dl/mu0.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        type(gs_spline_vacuum_t), intent(in) :: vac
        real(dp), intent(in) :: psi_c(:), dn_psi_c(:)
        real(dp), intent(out) :: coef(:, :)
        integer, intent(out) :: status
        type(csr_matrix_t) :: h
        type(sparse_direct_factor_t) :: factor
        real(dp), allocatable :: g(:, :), rhs(:), sol(:)
        real(dp) :: energy
        integer :: n, st, q, k
        status = GS_SPLINE_FAILED
        n = size(coef, 2)
        coef = 0
        allocate (g(1, n))
        call spline_energy(space, quad, kernel, coef, energy, st, g)
        if (st /= 0) return
        call spline_energy_hessian_csr(space, quad, kernel, coef, h, st)
        if (st /= 0) return
        rhs = -g(1, :)
        rhs(vac%dofs) = rhs(vac%dofs) + matmul(transpose(vac%trace), &
                                               matmul(vac%steklov, psi_c))/vac%mu0
        do q = 1, vac%side%point_count
            do k = vac%side_rowptr(q), vac%side_rowptr(q + 1) - 1
                rhs(vac%side_colidx(k)) = rhs(vac%side_colidx(k)) + vac%side%weight(q) &
                    *dn_psi_c(q)/vac%side%s(1, q)*vac%side_val(k)/vac%mu0
            end do
        end do
        call factor_with_block(h, vac%dofs, vac%block, factor, st)
        if (st /= 0) return
        allocate (sol(n))
        call sparse_direct_solve_factored(factor, rhs, sol, st)
        call sparse_direct_free(factor)
        if (st /= 0) return
        coef(1, :) = sol
        status = GS_SPLINE_OK
    end subroutine gs_spline_free_linear

    subroutine factor_with_block(h, dofs, block, factor, status)
        !! Factor H + B with B dense on the index set dofs (symmetric, so the
        !! CSR rows are the CSC columns).
        type(csr_matrix_t), intent(in) :: h
        integer, intent(in) :: dofs(:)
        real(dp), intent(in) :: block(:, :)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(out) :: status
        integer, allocatable :: colptr(:), rowind(:), local(:), mark(:)
        real(dp), allocatable :: values(:), work(:)
        integer :: n, i, k, j, nnz, c
        n = h%n
        allocate (local(n), mark(n), work(n))
        local = 0
        local(dofs) = [(k, k=1, size(dofs))]
        mark = 0
        allocate (colptr(n + 1))
        colptr(1) = 1
        ! count: union of the CSR row and, for a trace row, all trace dofs
        nnz = 0
        do i = 1, n
            c = 0
            do k = h%rowptr(i), h%rowptr(i + 1) - 1
                mark(h%colidx(k)) = i
                c = c + 1
            end do
            if (local(i) > 0) then
                do k = 1, size(dofs)
                    if (mark(dofs(k)) /= i) c = c + 1
                end do
            end if
            nnz = nnz + c
            colptr(i + 1) = nnz + 1
        end do
        allocate (rowind(nnz), values(nnz))
        mark = 0
        work = 0
        j = 0
        do i = 1, n
            do k = h%rowptr(i), h%rowptr(i + 1) - 1
                work(h%colidx(k)) = work(h%colidx(k)) + h%val(k)
                mark(h%colidx(k)) = i
            end do
            if (local(i) > 0) then
                do k = 1, size(dofs)
                    work(dofs(k)) = work(dofs(k)) + block(local(i), k)
                    mark(dofs(k)) = i
                end do
            end if
            do k = h%rowptr(i), h%rowptr(i + 1) - 1
                if (mark(h%colidx(k)) /= i) cycle
                j = j + 1
                rowind(j) = h%colidx(k)
                values(j) = work(h%colidx(k))
                work(h%colidx(k)) = 0
                mark(h%colidx(k)) = 0
            end do
            if (local(i) > 0) then
                do k = 1, size(dofs)
                    if (mark(dofs(k)) /= i) cycle
                    j = j + 1
                    rowind(j) = dofs(k)
                    values(j) = work(dofs(k))
                    work(dofs(k)) = 0
                    mark(dofs(k)) = 0
                end do
            end if
            colptr(i + 1) = j + 1
            call sort_segment(rowind(colptr(i):j), values(colptr(i):j))
        end do
        call sparse_direct_factor_csc(factor, n, colptr, rowind(:j), values(:j), status)
    contains
        pure subroutine sort_segment(idx, v)
            integer, intent(inout) :: idx(:)
            real(dp), intent(inout) :: v(:)
            integer :: a, b, ti
            real(dp) :: tv
            do a = 2, size(idx)
                ti = idx(a); tv = v(a); b = a - 1
                do while (b >= 1)
                    if (idx(b) <= ti) exit
                    idx(b + 1) = idx(b); v(b + 1) = v(b); b = b - 1
                end do
                idx(b + 1) = ti; v(b + 1) = tv
            end do
        end subroutine sort_segment
    end subroutine factor_with_block

    subroutine projection_kernel_evaluate(self, cell, first_point, npoint, jet, weight, &
                                          need, w, gradient, hessian, status)
        class(projection_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: g, h
        integer :: q
        associate (unused_cell => cell, unused_weight => weight)
        end associate
        status = 0
        do q = 1, npoint
            call gs_l2_projection(jet(1, q), self%v(first_point + q - 1), w(q), g, h)
            if (need >= 1) gradient(1, q) = g
            if (need >= 2) hessian(1, q) = h
        end do
    end subroutine projection_kernel_evaluate

    subroutine gs_spline_level_integrals(space, coef, origin, levels, sigma, n_rays, a1, a2, &
                                         status)
        !! a1 = oint R dl/|grad psi| and a2 = oint dl/(R |grad psi|) on the
        !! levels psi = levels(j), star shaped about origin and ordered from
        !! the inside out: on the ray of angle theta, dl/|grad psi| =
        !! rho dtheta/|d psi/d rho| (generated weights), trapezoidal rule on
        !! n_rays rays (spectrally accurate). Each root by Newton from the
        !! previous level's root and slope, bracketed from below by that root;
        !! an outward march in 1 cm steps only where Newton leaves the bracket.
        !! The rays are independent (OpenMP).
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :), origin(2), levels(:), sigma
        integer, intent(in) :: n_rays
        real(dp), intent(out) :: a1(:), a2(:)
        integer, intent(out) :: status
        real(dp) :: theta, e(2), lo, hi, rho, f, df, w1, w2, r_prev, s_prev
        integer :: k, j, it, block, st
        logical :: bracketed, failed
        a1 = 0
        a2 = 0
        failed = .false.
        !$omp parallel do default(shared) schedule(dynamic, 4) &
        !$omp private(k, j, it, theta, e, lo, hi, rho, f, df, w1, w2, r_prev, s_prev, &
        !$omp& block, st, bracketed) reduction(+:a1, a2) reduction(.or.:failed)
        do k = 1, n_rays
            if (failed) cycle
            theta = 2*acos(-1.0_dp)*(k - 1)/n_rays
            e = [cos(theta), sin(theta)]
            block = 1
            r_prev = 0
            s_prev = 0
            do j = 1, size(levels)
                lo = r_prev
                hi = huge(1.0_dp)
                bracketed = .false.
                if (s_prev < 0) then
                    rho = r_prev + sigma*(levels(j) - levels(max(j - 1, 1)))/s_prev
                    if (j == 1 .or. rho <= lo) rho = lo + 0.01_dp
                else
                    rho = lo + 0.01_dp
                end if
                st = 1
                do it = 1, 80
                    call ray_value(space, coef, origin, e, levels(j), sigma, rho, block, f, df, &
                                   st)
                    if (st /= 0) then
                        if (rho - lo <= 1.0e-12_dp) exit
                        rho = 0.5_dp*(lo + min(rho, hi))
                        cycle
                    end if
                    if (f > 0) then
                        lo = rho
                    else
                        hi = rho
                        bracketed = .true.
                    end if
                    if (abs(f) <= 1.0e-15_dp*max(1.0_dp, abs(levels(j))) .or. &
                        (bracketed .and. hi - lo <= 1.0e-14_dp)) exit
                    if (df < 0 .and. rho - f/df > lo .and. rho - f/df < hi) then
                        rho = rho - f/df
                    else if (bracketed) then
                        rho = 0.5_dp*(lo + hi)
                    else
                        rho = lo + 0.01_dp
                    end if
                end do
                if (st == 0) call ray_value(space, coef, origin, e, levels(j), sigma, rho, &
                                            block, f, df, st)
                if (st /= 0) then
                    failed = .true.
                    exit
                end if
                call gs_ray_coarea_weights(origin(1) + rho*e(1), rho, abs(df), w1, w2)
                a1(j) = a1(j) + w1
                a2(j) = a2(j) + w2
                r_prev = rho
                s_prev = df
            end do
        end do
        !$omp end parallel do
        status = merge(GS_SPLINE_FAILED, GS_SPLINE_OK, failed)
        a1 = a1*2*acos(-1.0_dp)/n_rays
        a2 = a2*2*acos(-1.0_dp)/n_rays
    end subroutine gs_spline_level_integrals

    subroutine ray_value(space, coef, origin, e, level, sigma, radius, block, f, df, status)
        !! f = sigma (psi - level) and df = d f/d rho at origin + radius e;
        !! block is the hint in and the located block out.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :), origin(2), e(2), level, sigma, radius
        integer, intent(inout) :: block
        real(dp), intent(out) :: f, df
        integer, intent(out) :: status
        real(dp) :: xi(2), j3(3)
        integer :: found
        f = 0
        df = 0
        call spline_locate(space, block, origin + radius*e, found, xi, status)
        if (status /= 0) return
        block = found
        call spline_evaluate_at(space, coef, block, xi, 1, j3, status)
        if (status /= 0) return
        f = sigma*(j3(1) - level)
        df = sigma*dot_product(j3(2:3), e)
    end subroutine ray_value

    subroutine gs_spline_project(space, quad, values, coef, status)
        !! L2 projection of values sampled at the quadrature points.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: coef(:, :)
        integer, intent(out) :: status
        type(projection_kernel_t) :: kernel
        type(csr_matrix_t) :: h
        real(dp), allocatable :: g(:, :), sol(:)
        real(dp) :: energy
        integer :: k
        kernel%order = 0
        kernel%v = values
        coef = 0
        allocate (g(1, size(coef, 2)))
        call spline_energy(space, quad, kernel, coef, energy, status, g)
        if (status /= 0) return
        call spline_energy_hessian_csr(space, quad, kernel, coef, h, status)
        if (status /= 0) return
        call solve_submatrix(h, [(k, k=1, size(coef, 2))], -g(1, :), sol, status)
        if (status /= 0) return
        coef(1, :) = sol
    end subroutine gs_spline_project

    subroutine gs_spline_newton_free(space, quad, kernel, vac, psi_c, dn_psi_c, sigma, &
                                     target_block, target_xi, coef, tolerance, &
                                     max_iterations, iterations, residual, axis, status)
        !! Free-boundary Newton: the law kernel on its plasma points with
        !! psi_b = mean of psi at the target points and psi_a the sigma-
        !! extremum over the plasma points, plus the wall terms of vac. The
        !! Jacobian H + T^T S T/mu0 + u_a v_a^T + u_b v_b^T is solved by
        !! Woodbury on one factorisation (u_b = D - u_a with D the flux
        !! derivative of the sources, v_a and v_b the basis values at the axis
        !! and the mean over the targets).
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(gs_spline_law_kernel_t), intent(inout) :: kernel
        type(gs_spline_vacuum_t), intent(in) :: vac
        real(dp), intent(in) :: psi_c(:), dn_psi_c(:), sigma, tolerance
        integer, intent(in) :: target_block(:), max_iterations
        real(dp), intent(in) :: target_xi(:, :)
        real(dp), intent(inout) :: coef(:, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out) :: residual, axis(2)
        type(axis_derivative_kernel_t) :: dkernel
        type(sparse_direct_factor_t) :: factor
        type(csr_matrix_t) :: h
        integer, allocatable :: rowptr(:), colidx(:), trow(:), tcol(:)
        real(dp), allocatable :: g(:, :), ua(:, :), dd(:, :), load(:), x(:), ya(:), yb(:), &
            va(:), vb(:), trial(:, :), val(:), tval(:)
        real(dp) :: energy, norm0, norm, damping, trial_norm, psi_a, psi_b, xi(2), m(2, 2), &
            rhs2(2), det, z(2)
        integer :: n, k, st, back, block, q

        status = GS_SPLINE_FAILED
        n = size(coef, 2)
        dkernel%order = 1
        dkernel%mu0 = kernel%mu0
        dkernel%r = kernel%r
        if (allocated(kernel%plasma)) dkernel%plasma = kernel%plasma
        ! the coil load oint (d_n psi_c/R) v dl/mu0
        allocate (load(n))
        load = 0
        do q = 1, vac%side%point_count
            do k = vac%side_rowptr(q), vac%side_rowptr(q + 1) - 1
                load(vac%side_colidx(k)) = load(vac%side_colidx(k)) + vac%side%weight(q) &
                    *dn_psi_c(q)/vac%side%s(1, q)*vac%side_val(k)/vac%mu0
            end do
        end do
        ! v_b: mean of the basis rows at the targets; without targets psi_b is
        ! the fixed kernel%psi_b (v_b = 0)
        allocate (vb(n))
        vb = 0
        if (size(target_block) > 0) then
            call spline_point_matrix(space, target_block, target_xi, 0, trow, tcol, tval, st)
            if (st /= 0) return
        else
            allocate (trow(1), tcol(0), tval(0))
            trow = 1
        end if
        do q = 1, size(target_block)
            do k = trow(q), trow(q + 1) - 1
                vb(tcol(k)) = vb(tcol(k)) + tval(k)/size(target_block)
            end do
        end do
        allocate (g(1, n), ua(1, n), dd(1, n), va(n), x(n), ya(n), yb(n))
        call residual_at(coef, norm, st)
        if (st /= 0) return
        ! scale: the larger of the first residual and the coil drive at the
        ! wall, so that a warm start is not judged against its own small
        ! residual
        norm0 = max(norm, norm2(matmul(transpose(vac%trace), matmul(vac%steklov, psi_c))) &
                    /vac%mu0 + norm2(load), tiny(1.0_dp))
        do iterations = 0, max_iterations
            residual = norm/norm0
            if (residual <= tolerance) then
                status = GS_SPLINE_OK
                return
            end if
            if (iterations == max_iterations) exit
            call spline_energy_hessian_csr(space, quad, kernel, coef, h, st)
            if (st /= 0) return
            call factor_with_block(h, vac%dofs, vac%block, factor, st)
            if (st /= 0) return
            dkernel%law = kernel%law
            dkernel%psi_b = kernel%psi_b
            dkernel%flux_derivative = .false.
            call spline_energy(space, quad, dkernel, coef, energy, st, ua)
            if (st /= 0) return
            dkernel%flux_derivative = .true.
            call spline_energy(space, quad, dkernel, coef, energy, st, dd)
            if (st /= 0) return
            call sparse_direct_solve_factored(factor, -g(1, :), x, st)
            if (st == 0) call sparse_direct_solve_factored(factor, ua(1, :), ya, st)
            if (st == 0) call sparse_direct_solve_factored(factor, dd(1, :) - ua(1, :), yb, st)
            call sparse_direct_free(factor)
            if (st /= 0) return
            call spline_locate(space, 1, axis, block, xi, st)
            if (st /= 0) return
            call spline_point_matrix(space, [block], reshape(xi, [2, 1]), 0, rowptr, &
                                     colidx, val, st)
            if (st /= 0) return
            va = 0
            va(colidx(rowptr(1):rowptr(2) - 1)) = val(rowptr(1):rowptr(2) - 1)
            ! Woodbury: (A + Y-terms) dc = -g with U = [u_a, u_b], V = [v_a, v_b]
            m(1, :) = [1 + dot_product(va, ya), dot_product(va, yb)]
            m(2, :) = [dot_product(vb, ya), 1 + dot_product(vb, yb)]
            rhs2 = [dot_product(va, x), dot_product(vb, x)]
            det = m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)
            if (abs(det) <= tiny(1.0_dp)) return
            z = [m(2, 2)*rhs2(1) - m(1, 2)*rhs2(2), -m(2, 1)*rhs2(1) + m(1, 1)*rhs2(2)]/det
            x = x - ya*z(1) - yb*z(2)
            damping = 1
            do back = 0, 10
                trial = coef
                trial(1, :) = trial(1, :) + damping*x
                call residual_at(trial, trial_norm, st)
                if (st == 0 .and. trial_norm < norm) exit
                damping = 0.5_dp*damping
            end do
            if (.not. trial_norm < norm) return
            coef = trial
            norm = trial_norm
        end do
    contains
        subroutine residual_at(c, value, st)
            !! Axis and boundary flux of c, the law there, and the norm of the
            !! full gradient (domain + wall).
            real(dp), intent(in) :: c(:, :)
            real(dp), intent(out) :: value
            integer, intent(out) :: st
            if (allocated(kernel%plasma)) then
                call gs_spline_extremum(space, quad, c, sigma, axis, psi_a, st, kernel%plasma)
            else
                call gs_spline_extremum(space, quad, c, sigma, axis, psi_a, st)
            end if
            if (st /= 0) return
            if (size(target_block) > 0) then
                psi_b = dot_product(vb, c(1, :))
            else
                psi_b = kernel%psi_b
            end if
            kernel%psi_b = psi_b
            kernel%law%axis = psi_a - psi_b
            call spline_energy(space, quad, kernel, c, energy, st, g)
            if (st /= 0) return
            g(1, vac%dofs) = g(1, vac%dofs) + matmul(transpose(vac%trace), &
                matmul(vac%steklov, matmul(vac%trace, c(1, vac%dofs)) - psi_c))/vac%mu0
            g(1, :) = g(1, :) - load
            value = norm2(g(1, :))
        end subroutine residual_at
    end subroutine gs_spline_newton_free

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

    subroutine gs_spline_extremum(space, quad, coef, sigma, point, value, status, mask)
        !! Location and value of the extremum of sigma*psi (sigma = +1 a
        !! maximum): the best quadrature point, then Newton on grad psi = 0
        !! with the exact spline Hessian.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: coef(:, :), sigma
        real(dp), intent(out) :: point(2), value
        integer, intent(out) :: status
        logical, intent(in), optional :: mask(:)
        real(dp), allocatable :: jet(:, :)
        real(dp) :: xi(2), j6(6), h(2, 2), dx(2), det
        integer :: q, block, it, st
        status = GS_SPLINE_FAILED
        allocate (jet(1, quad%point_count))
        call spline_evaluate(space, quad, coef, 0, jet, st)
        if (st /= 0) return
        if (present(mask)) then
            q = maxloc(sigma*jet(1, :), dim=1, mask=mask)
        else
            q = maxloc(sigma*jet(1, :), dim=1)
        end if
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
