module milca_eulerian_vacuum
    !! Vector-potential vacuum fields on FortFEM's spline de Rham complex
    !! (docs/eulerian-3d.md, gate G1): curl-curl operator K = C^T M2 C,
    !! the tree gauge relative to a Dirichlet boundary and the direct solve
    !! of K a = 0 with n x A given on the boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_eval, only: spline_quadrature_t
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_mass, derham_metric
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free
    use fortsparse_csr, only: csr_t, csr_transpose, csr_multiply, csr_budget_t
    use fortsparse_status, only: fortsparse_status_t
    implicit none
    private

    public :: eu_curl_curl, eu_boundary_tree, eu_constrained_solve, eu_to_csr

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

    subroutine eu_boundary_tree(derham, boundary_nodes, tree_edges)
        !! Gauge edges for A with n x A fixed: a spanning forest of the grad
        !! graph grown by breadth-first search from all boundary nodes at
        !! once, so the boundary acts as one root. Its edges (one per
        !! interior node) carry exactly the gauge A -> A + grad f, f = 0 on
        !! the boundary.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: boundary_nodes(:)
        integer, allocatable, intent(out) :: tree_edges(:)
        type(csr_t) :: g, gt
        type(csr_budget_t) :: budget
        type(fortsparse_status_t) :: info
        logical, allocatable :: seen(:)
        integer, allocatable :: queue(:), parent(:)
        integer :: head, tail, v, e, w, p, q, nv

        call eu_to_csr(derham%grad, derham%dimension(0), g)
        call csr_transpose(g, gt, budget, info)
        nv = derham%dimension(0)
        allocate (seen(nv), queue(nv), parent(nv))
        seen = .false.
        parent = 0
        seen(boundary_nodes) = .true.
        tail = size(boundary_nodes)
        queue(1:tail) = boundary_nodes
        head = 0
        do while (head < tail)
            head = head + 1
            v = queue(head)
            do p = gt%row_ptr(v), gt%row_ptr(v + 1) - 1
                e = gt%col(p)
                do q = g%row_ptr(e), g%row_ptr(e + 1) - 1
                    w = g%col(q)
                    if (seen(w) .or. g%val(q) == 0) cycle
                    seen(w) = .true.
                    parent(w) = e
                    tail = tail + 1
                    queue(tail) = w
                end do
            end do
        end do
        tree_edges = pack(parent, parent > 0)
    end subroutine eu_boundary_tree

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

end module milca_eulerian_vacuum
