module milca_eulerian_vacuum
    !! Vector-potential vacuum fields on FortFEM's spline de Rham complex
    !! (docs/eulerian-3d.md, gate G1): curl-curl operator K = C^T M2 C and
    !! the direct solve of K a = rhs with fixed (Dirichlet and gauge) dofs.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_eval, only: spline_quadrature_t
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_mass, derham_metric
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free, sparse_direct_refactor_csc, &
        sparse_direct_factor_symmetric_csc, sparse_direct_inertia, SPARSE_DIRECT_LOWER
    use fortsparse_csr, only: csr_t, csr_transpose, csr_multiply, csr_budget_t, csr_matvec
    use fortsparse_status, only: fortsparse_status_t
    use fortsparse_pcg, only: pcg_params_t, pcg_handle_t, pcg_setup_ams, pcg_solve, &
        pcg_destroy, PCG_AMS
    implicit none
    private

    public :: eu_curl_curl, eu_constrained_solve, eu_to_csr, eu_ams_solve
    public :: eu_constrained_t, eu_constrained_factor, eu_constrained_apply, &
              eu_constrained_free, eu_constrained_refactor, eu_constrained_inertia

    type :: eu_constrained_t
        !! A factorised free block of K (fixed rows and columns removed), for
        !! several right-hand sides with homogeneous fixed values: a Newton
        !! step and its simplified corrections (affine-covariant damping).
        type(sparse_direct_factor_t) :: factor
        integer, allocatable :: free(:)
        integer :: n = 0
        !! K given as its upper triangle by rows, factorised by LDL^T
        logical :: symmetric = .false.
    end type eu_constrained_t

contains

    subroutine eu_to_csr(a, ncol, b)
        type(csr_matrix_t), intent(in) :: a
        integer, intent(in) :: ncol
        type(csr_t), intent(out) :: b

        b%nrow = a%n
        b%ncol = ncol
        b%row_ptr = a%rowptr(1:a%n + 1)
        b%col = a%colidx(1:a%rowptr(a%n + 1) - 1)
        b%val = a%val(1:a%rowptr(a%n + 1) - 1)
    end subroutine eu_to_csr

    subroutine eu_curl_curl(derham, quad, metric, k, status, points)
        !! K = C^T M2 C with the 2-form mass weight `metric` (G/|J| for
        !! |curl A|^2 dV), so a^T K a = int |B|^2 dV.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(derham_metric) :: metric
        type(csr_t), intent(out) :: k
        integer, intent(out) :: status
        integer, intent(in), optional :: points
        type(csr_matrix_t) :: m2
        type(csr_t) :: c, ct, m, mc
        type(csr_budget_t) :: budget
        type(fortsparse_status_t) :: info

        call spline_derham_mass(derham, quad, 2, metric, m2, status, points)
        if (status /= 0) return
        call eu_to_csr(derham%curl, derham%dimension(1), c)
        call eu_to_csr(m2, derham%dimension(2), m)
        call csr_multiply(m, c, mc, budget, info)
        if (info%code == 0) call csr_transpose(c, ct, budget, info)
        if (info%code == 0) call csr_multiply(ct, mc, k, budget, info)
        status = info%code
    end subroutine eu_curl_curl

    subroutine eu_constrained_solve(k, rhs, fixed, values, a, status)
        !! Solves K a = rhs on the free rows with a(fixed) = values, by a
        !! sparse direct factorisation of the free block (K symmetric, so its
        !! CSR rows are the CSC columns).
        type(csr_t), intent(in) :: k
        real(dp), intent(in) :: rhs(:)
        integer, intent(in) :: fixed(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: a(:)
        integer, intent(out) :: status
        type(sparse_direct_factor_t) :: factor
        integer, allocatable :: local(:), free(:), colptr(:), rowind(:)
        real(dp), allocatable :: val(:), b(:), x(:)
        integer :: n, nf, i, j, p, nz

        n = k%nrow
        allocate (local(n))
        local = 1
        local(fixed) = 0
        free = pack([(i, i=1, n)], local == 1)
        nf = size(free)
        local = 0
        local(free) = [(i, i=1, nf)]
        a = 0
        a(fixed) = values
        allocate (colptr(nf + 1), rowind(size(k%val)), val(size(k%val)), b(nf))
        colptr(1) = 1
        nz = 0
        do j = 1, nf
            i = free(j)
            b(j) = rhs(i)
            do p = k%row_ptr(i), k%row_ptr(i + 1) - 1
                if (local(k%col(p)) > 0) then
                    nz = nz + 1
                    rowind(nz) = local(k%col(p))
                    val(nz) = k%val(p)
                else
                    b(j) = b(j) - k%val(p)*a(k%col(p))
                end if
            end do
            colptr(j + 1) = nz + 1
        end do
        call sparse_direct_factor_csc(factor, nf, colptr, rowind(:nz), val(:nz), status)
        if (status /= 0) return
        allocate (x(nf))
        call sparse_direct_solve_factored(factor, b, x, status)
        call sparse_direct_free(factor)
        if (status == 0) a(free) = x
    end subroutine eu_constrained_solve

    subroutine eu_constrained_factor(k, fixed, c, status, symmetric)
        !! Factorises the free block of K (K stored by rows of K^T, as in
        !! eu_constrained_solve). With symmetric, K is symmetric and k holds
        !! only its upper triangle (columns >= row, diagonal included); the
        !! free block is factorised by LDL^T, which also gives the inertia
        !! (eu_constrained_inertia).
        type(csr_t), intent(in) :: k
        integer, intent(in) :: fixed(:)
        type(eu_constrained_t), intent(out) :: c
        integer, intent(out) :: status
        logical, intent(in), optional :: symmetric
        integer, allocatable :: local(:), colptr(:), rowind(:)
        real(dp), allocatable :: val(:)
        integer :: nf, i, j, p, nz

        c%n = k%nrow
        c%symmetric = .false.
        if (present(symmetric)) c%symmetric = symmetric
        allocate (local(c%n))
        local = 1
        local(fixed) = 0
        c%free = pack([(i, i=1, c%n)], local == 1)
        nf = size(c%free)
        local = 0
        local(c%free) = [(i, i=1, nf)]
        allocate (colptr(nf + 1), rowind(size(k%val)), val(size(k%val)))
        colptr(1) = 1
        nz = 0
        do j = 1, nf
            i = c%free(j)
            do p = k%row_ptr(i), k%row_ptr(i + 1) - 1
                if (local(k%col(p)) > 0) then
                    nz = nz + 1
                    rowind(nz) = local(k%col(p))
                    val(nz) = k%val(p)
                end if
            end do
            colptr(j + 1) = nz + 1
        end do
        if (c%symmetric) then
            ! the rows of the upper triangle are the columns of the lower one
            call sparse_direct_factor_symmetric_csc(c%factor, nf, colptr, rowind(:nz), &
                val(:nz), status, SPARSE_DIRECT_LOWER)
        else
            call sparse_direct_factor_csc(c%factor, nf, colptr, rowind(:nz), val(:nz), status)
        end if
    end subroutine eu_constrained_factor

    subroutine eu_constrained_refactor(k, fixed, c, status, symmetric)
        !! As eu_constrained_factor for a K whose sparsity pattern and fixed
        !! slots equal those of the last factorisation of c (Newton steps,
        !! shifted Hessians): FortSparse keeps the fill-reducing ordering and
        !! the symbolic analysis and redoes only the numerical factorisation.
        !! A first call (c not factorised), another size, or another
        !! symmetric flag factorises afresh. symmetric as in
        !! eu_constrained_factor (k the upper triangle).
        type(csr_t), intent(in) :: k
        integer, intent(in) :: fixed(:)
        type(eu_constrained_t), intent(inout) :: c
        integer, intent(out) :: status
        logical, intent(in), optional :: symmetric
        integer, allocatable :: local(:), colptr(:), rowind(:)
        real(dp), allocatable :: val(:)
        integer :: nf, i, j, p, nz
        logical :: sym

        sym = .false.
        if (present(symmetric)) sym = symmetric
        if (c%n /= k%nrow .or. .not. allocated(c%free) .or. (c%symmetric .neqv. sym)) then
            call eu_constrained_free(c)
            call eu_constrained_factor(k, fixed, c, status, sym)
            return
        end if
        nf = size(c%free)
        allocate (local(c%n))
        local = 0
        local(c%free) = [(i, i=1, nf)]
        allocate (colptr(nf + 1), rowind(size(k%val)), val(size(k%val)))
        colptr(1) = 1
        nz = 0
        do j = 1, nf
            i = c%free(j)
            do p = k%row_ptr(i), k%row_ptr(i + 1) - 1
                if (local(k%col(p)) > 0) then
                    nz = nz + 1
                    rowind(nz) = local(k%col(p))
                    val(nz) = k%val(p)
                end if
            end do
            colptr(j + 1) = nz + 1
        end do
        call sparse_direct_refactor_csc(c%factor, nf, colptr, rowind(:nz), val(:nz), status)
    end subroutine eu_constrained_refactor

    subroutine eu_constrained_apply(c, rhs, a, status)
        !! a = K_ff^(-1) rhs_f on the free rows, 0 on the fixed rows.
        type(eu_constrained_t), intent(inout) :: c
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: a(:)
        integer, intent(out) :: status
        real(dp), allocatable :: x(:)

        allocate (x(size(c%free)))
        call sparse_direct_solve_factored(c%factor, rhs(c%free), x, status)
        a = 0
        if (status == 0) a(c%free) = x
    end subroutine eu_constrained_apply

    subroutine eu_constrained_inertia(c, negative, zero, positive, reliable, status)
        !! Eigenvalue counts of the free block of the last symmetric
        !! factorisation of c (Sylvester inertia of its LDL^T); status /= 0
        !! for an LU factor.
        type(eu_constrained_t), intent(inout) :: c
        integer, intent(out) :: negative, zero, positive
        logical, intent(out) :: reliable
        integer, intent(out) :: status

        call sparse_direct_inertia(c%factor, negative, zero, positive, status, reliable)
    end subroutine eu_constrained_inertia

    subroutine eu_constrained_free(c)
        type(eu_constrained_t), intent(inout) :: c

        call sparse_direct_free(c%factor)
        if (allocated(c%free)) deallocate (c%free)
        c%n = 0
        c%symmetric = .false.
    end subroutine eu_constrained_free

    subroutine eu_ams_solve(k, grad, interpolation, n0, fixed, values, fixed_nodes, a, &
                            rtol, iterations, status, mass)
        !! Ungauged solve of K a = 0 with a(fixed) = values by PCG with the AMS
        !! preconditioner on the free edges. Gradients of functions that
        !! vanish on fixed_nodes are the kernel of the free block, and the
        !! right-hand side -K_fd a_d is orthogonal to them, so PCG converges in
        !! range(K). grad is n1 x n0 and interpolation n1 x 3 n0. With the
        !! 1-form mass matrix, AMS solves its nodal components with the vector
        !! Laplacian G^T M1 G (Hiptmair-Xu) instead of blocks of Pi^T K Pi.
        type(csr_t), intent(in) :: k
        type(csr_matrix_t), intent(in) :: grad, interpolation
        integer, intent(in) :: n0, fixed(:), fixed_nodes(:)
        real(dp), intent(in) :: values(:), rtol
        real(dp), intent(out) :: a(:)
        integer, intent(out) :: iterations, status
        type(csr_matrix_t), intent(in), optional :: mass
        type(csr_t) :: kf, gf, pf, g, p, mf, m
        type(pcg_params_t) :: params
        type(pcg_handle_t) :: handle
        type(fortsparse_status_t) :: info
        integer, allocatable :: row_map(:), node_map(:), pi_map(:)
        real(dp), allocatable :: rhs(:), x(:), kd(:)
        integer :: n1, i, d, nf, ni

        n1 = k%nrow
        allocate (row_map(n1), node_map(n0))
        row_map = 1
        row_map(fixed) = 0
        nf = 0
        do i = 1, n1
            if (row_map(i) > 0) then
                nf = nf + 1
                row_map(i) = nf
            end if
        end do
        node_map = 1
        node_map(fixed_nodes) = 0
        ni = 0
        do i = 1, n0
            if (node_map(i) > 0) then
                ni = ni + 1
                node_map(i) = ni
            end if
        end do
        allocate (pi_map(3*n0))
        do d = 1, 3
            pi_map((d - 1)*n0 + 1:d*n0) = merge(node_map + (d - 1)*ni, 0, node_map > 0)
        end do
        call submatrix(k, row_map, row_map, nf, kf)
        call eu_to_csr(grad, n0, g)
        call submatrix(g, row_map, node_map, ni, gf)
        call eu_to_csr(interpolation, 3*n0, p)
        call submatrix(p, row_map, pi_map, 3*ni, pf)
        a = 0
        a(fixed) = values
        allocate (kd(n1), rhs(nf), x(nf))
        call csr_matvec(k, a, kd)
        rhs = -pack(kd, row_map > 0)
        params%preconditioner = PCG_AMS
        params%rtol = rtol
        params%ams%smoothing_steps = 2
        if (present(mass)) then
            call eu_to_csr(mass, n1, m)
            call submatrix(m, row_map, row_map, nf, mf)
            call pcg_setup_ams(kf, gf, pf, params, handle, info, mass=mf)
        else
            call pcg_setup_ams(kf, gf, pf, params, handle, info)
        end if
        status = info%code
        if (status /= 0) return
        x = 0
        call pcg_solve(handle, params, rhs, x, info, iterations)
        status = info%code
        call pcg_destroy(handle)
        do i = 1, n1
            if (row_map(i) > 0) a(i) = x(row_map(i))
        end do
    end subroutine eu_ams_solve

    subroutine submatrix(a, row_map, col_map, ncol, b)
        !! Rows and columns with a positive map entry, renumbered by it.
        type(csr_t), intent(in) :: a
        integer, intent(in) :: row_map(:), col_map(:), ncol
        type(csr_t), intent(out) :: b
        integer :: i, p, nz, nr

        nr = count(row_map > 0)
        b%nrow = nr
        b%ncol = ncol
        allocate (b%row_ptr(nr + 1), b%col(size(a%val)), b%val(size(a%val)))
        b%row_ptr(1) = 1
        nz = 0
        do i = 1, a%nrow
            if (row_map(i) == 0) cycle
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (col_map(a%col(p)) == 0) cycle
                nz = nz + 1
                b%col(nz) = col_map(a%col(p))
                b%val(nz) = a%val(p)
            end do
            b%row_ptr(row_map(i) + 1) = nz + 1
        end do
        b%col = b%col(:nz)
        b%val = b%val(:nz)
    end subroutine submatrix

end module milca_eulerian_vacuum
