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
        spline_boundary_nodes
    use fortfem_spline_energy, only: pointwise_energy_t, spline_energy, &
        spline_energy_hessian_csr
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free
    implicit none
    private
    public :: gs_spline_kernel_t, gs_spline_dirichlet_solve, boundary_value_t
    integer, parameter, public :: GS_SPLINE_OK = 0, GS_SPLINE_FAILED = 1

    type, extends(pointwise_energy_t) :: gs_spline_kernel_t
        !> R and the source f at every quadrature point
        real(dp), allocatable :: r(:), f(:)
    contains
        procedure :: evaluate => gs_kernel_evaluate
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
