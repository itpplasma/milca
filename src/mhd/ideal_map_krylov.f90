module milca_ideal_map_krylov
    ! Deflated preconditioned conjugate gradients with recycled Ritz vectors
    ! for A x = b, A symmetric. Given an orthonormal basis W of directions
    ! that the preconditioner M leaves soft (small eigenvalues of M^-1 A),
    ! the iteration solves exactly on span(W) through E = W^T A W and runs
    ! conjugate gradients on the A-orthogonal complement (Saad, Yeung, Erhel
    ! and Guyomarc'h, SIAM J. Sci. Comput. 21, 1909 (2000)): x0 = W E^-1 W^T b,
    ! p = z - W E^-1 (A W)^T z. From the Lanczos tridiagonal of the run (CG
    ! coefficients alpha, beta) it returns the Ritz vectors of the smallest
    ! Ritz values of the deflated, preconditioned operator, which the next
    ! solve (the next Newton step) adds to W. A direction of non-positive
    ! curvature p^T A p <= 0 ends the solve with negative = .true.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: deflated_pcg, krylov_operator, krylov_preconditioner

    abstract interface
        subroutine krylov_operator(v, w)
            import :: dp
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: w(:)
        end subroutine
        subroutine krylov_preconditioner(v, w)
            import :: dp
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: w(:)
        end subroutine
    end interface
contains

    subroutine deflated_pcg(apply_a, apply_m, b, tolerance, max_iterations, w, x, &
            iterations, negative, nritz, ritz)
        ! x: solution to |b - A x| <= tolerance |b|; w(n, k): deflation basis
        ! (k may be 0, orthonormal columns); ritz(n, nritz): Ritz vectors of the
        ! nritz smallest Ritz values of the run (fewer if the run is shorter).
        procedure(krylov_operator) :: apply_a
        procedure(krylov_preconditioner) :: apply_m
        real(dp), intent(in) :: b(:), tolerance, w(:, :)
        integer, intent(in) :: max_iterations, nritz
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: iterations
        logical, intent(out) :: negative
        real(dp), allocatable, intent(out) :: ritz(:, :)
        real(dp), allocatable :: aw(:, :), e(:, :), r(:), z(:), p(:), q(:), t(:)
        real(dp), allocatable :: lanczos(:, :), alphas(:), betas(:)
        real(dp) :: rz, rz_new, alpha, beta, b0, pq
        integer :: n, k, j, info
        integer, allocatable :: pivot(:)
        n = size(b); k = size(w, 2)
        negative = .false.
        x = 0
        allocate (aw(n, k), e(k, k), pivot(k), t(k), r(n), z(n), p(n), q(n))
        allocate (lanczos(n, max_iterations + 1), alphas(max_iterations), &
            betas(max_iterations))
        do j = 1, k
            call apply_a(w(:, j), aw(:, j))
        end do
        if (k > 0) then
            e = matmul(transpose(w), aw)
            call dgetrf(k, k, e, k, pivot, info)
            if (info /= 0) then
                ! a singular coarse matrix: drop the deflation
                k = 0
            end if
        end if
        if (k > 0) then
            t = matmul(transpose(w(:, :k)), b)
            call dgetrs('N', k, 1, e, k, pivot, t, k, info)
            x = matmul(w(:, :k), t)
            r = b - matmul(aw(:, :k), t)
        else
            r = b
        end if
        b0 = norm2(b)
        iterations = 0
        if (.not. b0 > 0) then
            allocate (ritz(n, 0))
            return
        end if
        call apply_m(r, z)
        p = z
        call deflate(p, z)
        rz = dot_product(r, z)
        lanczos(:, 1) = z/sqrt(max(rz, tiny(1.0_dp)))
        do iterations = 1, max_iterations
            call apply_a(p, q)
            pq = dot_product(p, q)
            if (.not. pq > 0) then
                negative = .true.
                exit
            end if
            alpha = rz/pq
            x = x + alpha*p
            r = r - alpha*q
            alphas(iterations) = alpha
            if (norm2(r) <= tolerance*b0) exit
            call apply_m(r, z)
            rz_new = dot_product(r, z)
            beta = rz_new/rz
            betas(iterations) = beta
            p = z + beta*p
            call deflate(p, z)
            rz = rz_new
            lanczos(:, iterations + 1) = z/sqrt(max(rz, tiny(1.0_dp)))
        end do
        iterations = min(iterations, max_iterations)
        call ritz_vectors()
    contains
        subroutine deflate(v, zz)
            ! v <- v - W E^-1 (A W)^T zz
            real(dp), intent(inout) :: v(:)
            real(dp), intent(in) :: zz(:)
            real(dp) :: c(k)
            integer :: inf
            if (k == 0) return
            c = matmul(transpose(aw(:, :k)), zz)
            call dgetrs('N', k, 1, e, k, pivot, c, k, inf)
            v = v - matmul(w(:, :k), c)
        end subroutine

        subroutine ritz_vectors()
            ! Lanczos tridiagonal of the run: T(i, i) = 1/alpha_i +
            ! beta_(i-1)/alpha_(i-1), T(i, i+1) = -sqrt(beta_i)/alpha_i
            real(dp), allocatable :: tt(:, :), ev(:), work(:)
            integer :: m, i, nr, lwork, inf
            m = iterations
            if (negative) m = iterations - 1
            m = max(0, min(m, size(alphas)))
            nr = min(nritz, m)
            allocate (ritz(n, nr))
            if (nr == 0) return
            allocate (tt(m, m), ev(m))
            tt = 0
            do i = 1, m
                tt(i, i) = 1/alphas(i)
                if (i > 1) tt(i, i) = tt(i, i) + betas(i - 1)/alphas(i - 1)
                if (i < m) then
                    tt(i, i + 1) = -sqrt(betas(i))/alphas(i)
                    tt(i + 1, i) = tt(i, i + 1)
                end if
            end do
            lwork = max(1, 3*m)
            allocate (work(lwork))
            call dsyev('V', 'U', m, tt, m, ev, work, lwork, inf)
            if (inf /= 0) then
                deallocate (ritz)
                allocate (ritz(n, 0))
                return
            end if
            ! eigenvalues ascending: the first nr are the softest
            ritz = matmul(lanczos(:, :m), tt(:, :nr))
        end subroutine
    end subroutine
end module
