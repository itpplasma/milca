module milca_ideal_map_mode_pc
    ! Overlapping additive Schwarz preconditioner of the map Hessian in the
    ! toroidal number (position gauge). The Hessian couples the Fourier modes
    ! m, m' through sum_j t_m t'_m' d2w/dv2 over the toroidal planes. Each
    ! subdomain p = 0 .. ntor - 1 holds the cos/sin modes of n = p and p + 1
    ! with their exact couplings, a principal submatrix of H and so positive
    ! definite where H is; the preconditioner is
    !   M^-1 = sum_p W R_p^T (R_p H R_p + S_p)^-1 R_p W,
    ! W = diag(1/sqrt(multiplicity of n)), symmetric for conjugate gradients.
    ! It keeps the first toroidal harmonic of the plane variation, which
    ! carries the near-null relabelling directions that a plane average
    ! loses, in O(nodes nmode) memory. Each subdomain matrix is a 2D mesh
    ! matrix over (q1, q2) of its nodes and modes, factored by Cholesky with
    ! a Levenberg-Marquardt shift S_p raised until positive definite;
    ! field-angle coefficients take the diagonal of their block. For an
    ! axisymmetric state (no coupling between different n) it is H^-1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_t, ideal_map_rdof, ideal_map_angle_rdof, &
        ideal_map_angle_diagonal, ideal_map_hessian_store_t, IDEAL_MAP_OK, &
        IDEAL_MAP_BAD_INPUT, IDEAL_MAP_SOLVE_FAILED
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, csr_assemble_blocks
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_SPD, &
        FORTSPARSE_ORDER_AMD
    implicit none
    private
    public :: ideal_map_mode_pc_t, ideal_map_mode_pc_build, ideal_map_mode_pc_apply, &
        ideal_map_mode_pc_destroy

    type :: subdomain_t
        ! reduced dofs of the subdomain and their weights 1/sqrt(multiplicity)
        integer, allocatable :: rdof(:)
        real(dp), allocatable :: weight(:)
        type(csr_matrix_t) :: matrix
        ! the Levenberg-Marquardt shift added to each subdomain dof
        real(dp), allocatable :: shift(:)
        type(sparse_solver_t) :: solver
        logical :: factored = .false.
    end type

    type :: ideal_map_mode_pc_t
        type(subdomain_t), allocatable :: part(:)
        integer, allocatable :: angle_rdof(:)
        real(dp), allocatable :: angle_diag(:)
    end type
contains

    subroutine ideal_map_mode_pc_build(map, store, mu, pc, status)
        ! Subdomain matrices from the stored point Hessians with the diagonal
        ! shift mu |diag| (raised by 10 until positive definite).
        type(ideal_map_t), intent(in) :: map
        type(ideal_map_hessian_store_t), intent(in) :: store
        real(dp), intent(in) :: mu
        type(ideal_map_mode_pc_t), intent(inout) :: pc
        integer, intent(out) :: status
        integer :: npart, p, k
        status = IDEAL_MAP_BAD_INPUT
        if (.not. map%position_gauge .or. .not. store%ready) return
        call ideal_map_mode_pc_destroy(pc)
        npart = max(1, map%ntor)
        allocate (pc%part(0:npart - 1))
        do p = 0, npart - 1
            call build_part(map, store, p, mu, pc%part(p), status)
            if (status /= IDEAL_MAP_OK) return
        end do
        if (map%nangle > 0) then
            allocate (pc%angle_rdof(map%nangle), pc%angle_diag(map%nangle))
            call ideal_map_angle_diagonal(map, map%x, pc%angle_diag)
            do k = 1, map%nangle
                pc%angle_rdof(k) = ideal_map_angle_rdof(map, k)
            end do
            pc%angle_diag = (1 + mu)*abs(pc%angle_diag) &
                + epsilon(1.0_dp)*maxval(abs(pc%angle_diag))
            where (.not. pc%angle_diag > 0) pc%angle_diag = 1
        end if
        status = IDEAL_MAP_OK
    end subroutine

    pure subroutine part_modes(ntor, p, mode)
        ! Fourier modes of the toroidal numbers p and p + 1 (p = 0 alone when
        ! ntor = 0); mode 1 is n = 0, modes 2 n and 2 n + 1 the cos and sin of n
        integer, intent(in) :: ntor, p
        integer, allocatable, intent(out) :: mode(:)
        if (ntor == 0) then
            mode = [1]
        else if (p == 0) then
            mode = [1, 2, 3]
        else
            mode = [2*p, 2*p + 1, 2*p + 2, 2*p + 3]
        end if
    end subroutine

    subroutine build_part(map, store, p, mu, part, status)
        ! Subdomain of the toroidal numbers p and p + 1 (p alone for ntor = 0).
        type(ideal_map_t), intent(in) :: map
        type(ideal_map_hessian_store_t), intent(in) :: store
        integer, intent(in) :: p
        real(dp), intent(in) :: mu
        type(subdomain_t), intent(inout) :: part
        integer, intent(out) :: status
        real(dp), allocatable :: ke(:, :, :)
        integer, allocatable :: conn(:, :), mode(:), local(:)
        logical, allocatable :: fixed(:)
        type(fortsparse_status_t) :: sstatus
        real(dp) :: scale, sh
        integer :: nm, ne, row, k, tries, e, node, i, c, n
        status = IDEAL_MAP_SOLVE_FAILED
        call part_modes(map%ntor, p, mode)
        nm = size(mode)
        ne = 2*nm*map%nloc
        allocate (part%rdof(2*nm*map%nodes), part%weight(2*nm*map%nodes), &
            fixed(2*nm*map%nodes), local(map%nred))
        local = 0
        do node = 1, map%nodes
            do i = 1, nm
                n = mode(i)/2
                do c = 1, 2
                    k = c + 2*(i - 1) + 2*nm*(node - 1)
                    part%rdof(k) = ideal_map_rdof(map, c, mode(i), node)
                    local(part%rdof(k)) = k
                    fixed(k) = map%fixed(part%rdof(k))
                    part%weight(k) = 1
                    if (n > 0 .and. n < map%ntor) part%weight(k) = 1/sqrt(2.0_dp)
                end do
            end do
        end do
        allocate (conn(ne, map%cells), ke(ne, ne, map%cells))
        !$omp parallel do default(shared) private(e) schedule(dynamic)
        do e = 1, map%cells
            call mode_element(map, store, e, mode, conn(:, e), ke(:, :, e))
            conn(:, e) = local(conn(:, e))
        end do
        !$omp end parallel do
        call csr_build_pattern(size(part%rdof), conn, .false., part%matrix, &
            with_slots=.true.)
        call csr_assemble_blocks(part%matrix, conn, ke)
        deallocate (ke, conn)
        do row = 1, part%matrix%n
            do k = part%matrix%rowptr(row), part%matrix%rowptr(row + 1) - 1
                if (fixed(row) .or. fixed(part%matrix%colidx(k))) part%matrix%val(k) = 0
                if (fixed(row) .and. part%matrix%colidx(k) == row) part%matrix%val(k) = 1
            end do
        end do
        allocate (part%shift(part%matrix%n))
        scale = 0
        do row = 1, part%matrix%n
            if (fixed(row)) cycle
            do k = part%matrix%rowptr(row), part%matrix%rowptr(row + 1) - 1
                if (part%matrix%colidx(k) == row) scale = max(scale, abs(part%matrix%val(k)))
            end do
        end do
        sh = mu
        do tries = 1, 40
            call shifted_factor(sh, sstatus)
            if (status_ok(sstatus)) then
                status = IDEAL_MAP_OK
                return
            end if
            sh = max(10*sh, 1.0e-12_dp)
        end do
    contains
        subroutine shifted_factor(t, st)
            ! factor of the matrix plus t (|diag| + 1e-12 scale) on the free diagonal
            real(dp), intent(in) :: t
            type(fortsparse_status_t), intent(out) :: st
            real(dp), allocatable :: val(:)
            integer :: rr, kk
            val = part%matrix%val
            part%shift = 0
            do rr = 1, part%matrix%n
                if (fixed(rr)) cycle
                do kk = part%matrix%rowptr(rr), part%matrix%rowptr(rr + 1) - 1
                    if (part%matrix%colidx(kk) /= rr) cycle
                    part%shift(rr) = t*(abs(val(kk)) + 1.0e-12_dp*scale)
                    val(kk) = val(kk) + part%shift(rr)
                end do
            end do
            if (part%factored) call sparse_destroy(part%solver)
            part%factored = .false.
            part%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            part%solver%options%symmetry = FORTSPARSE_SYM_SPD
            part%solver%options%ordering = FORTSPARSE_ORDER_AMD
            call sparse_factor(part%solver, part%matrix%n, part%matrix%n, size(val), &
                part%matrix%rowptr, part%matrix%colidx, val, st)
            part%factored = status_ok(st)
        end subroutine
    end subroutine

    subroutine mode_element(map, store, e, mode, conn, ke)
        ! Element of cell e for the given Fourier modes: their exact Hessian
        ! couplings sum_j,g w phi^T d2w/dv2 phi, exactly symmetric; conn holds
        ! the reduced dofs (c, mode, node).
        type(ideal_map_t), intent(in) :: map
        type(ideal_map_hessian_store_t), intent(in) :: store
        integer, intent(in) :: e, mode(:)
        integer, intent(out) :: conn(:)
        real(dp), intent(out) :: ke(:, :)
        real(dp), allocatable :: phi(:, :), hphi(:, :)
        real(dp) :: hv(12, 12), weight, dphi
        integer :: nm, ne, a, i, c, l, g, j, k, pi, pj
        nm = size(mode)
        ne = 2*nm*map%nloc
        do a = 1, map%nloc
            do i = 1, nm
                do c = 1, 2
                    conn(c + 2*(i - 1) + 2*nm*(a - 1)) = &
                        ideal_map_rdof(map, c, mode(i), map%conn(a, e))
                end do
            end do
        end do
        ke = 0
        allocate (phi(12, ne), hphi(12, ne))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        do j = 1, map%nphi
            do g = 1, map%nq
                weight = map%weight(g)*dphi
                do k = 1, 78
                    call packed(k, pi, pj)
                    hv(pi, pj) = store%hp(k, g, j, e); hv(pj, pi) = hv(pi, pj)
                end do
                phi = 0
                do a = 1, map%nloc
                    do i = 1, nm
                        do c = 1, 2
                            l = c + 2*(i - 1) + 2*nm*(a - 1)
                            phi(c, l) = map%shape(a, g)*map%trig(mode(i), j)
                            phi(3 + c, l) = map%dshape(1, a, g)*map%trig(mode(i), j)
                            phi(6 + c, l) = map%dshape(2, a, g)*map%trig(mode(i), j)
                            phi(9 + c, l) = map%shape(a, g)*map%dtrig(mode(i), j)
                        end do
                    end do
                end do
                hphi = weight*matmul(hv, phi)
                ke = ke + matmul(transpose(phi), hphi)
            end do
        end do
        ke = (ke + transpose(ke))/2
    contains
        pure subroutine packed(n, ri, rj)
            ! (row, column) of entry n of the hp packing (gen/gen_ideal_map.f90)
            integer, intent(in) :: n
            integer, intent(out) :: ri, rj
            integer :: tt
            if (n <= 45) then
                rj = 1
                do while (rj*(rj + 1)/2 < n)
                    rj = rj + 1
                end do
                ri = n - rj*(rj - 1)/2
            else if (n <= 72) then
                tt = n - 46
                rj = 10 + tt/9
                ri = 1 + mod(tt, 9)
            else
                tt = n - 73
                rj = 10
                do while ((rj - 9)*(rj - 8)/2 < tt + 1)
                    rj = rj + 1
                end do
                ri = 10 + tt - (rj - 10)*(rj - 9)/2
            end if
        end subroutine
    end subroutine

    subroutine ideal_map_mode_pc_apply(map, pc, r, z, status)
        ! z = M^-1 r in the reduced dofs; fixed dofs pass through.
        type(ideal_map_t), intent(in) :: map
        type(ideal_map_mode_pc_t), intent(inout) :: pc
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        type(fortsparse_status_t) :: sstatus
        real(dp), allocatable :: rb(:), zb(:)
        integer :: p
        status = IDEAL_MAP_SOLVE_FAILED
        z = 0
        do p = lbound(pc%part, 1), ubound(pc%part, 1)
            rb = pc%part(p)%weight*r(pc%part(p)%rdof)
            allocate (zb(size(rb)))
            call sparse_solve(pc%part(p)%solver, rb, zb, sstatus)
            if (.not. status_ok(sstatus)) return
            z(pc%part(p)%rdof) = z(pc%part(p)%rdof) + pc%part(p)%weight*zb
            deallocate (zb)
        end do
        if (allocated(pc%angle_rdof)) z(pc%angle_rdof) = r(pc%angle_rdof)/pc%angle_diag
        where (map%fixed) z = r
        status = IDEAL_MAP_OK
    end subroutine

    subroutine ideal_map_mode_pc_destroy(pc)
        type(ideal_map_mode_pc_t), intent(inout) :: pc
        integer :: p
        if (allocated(pc%part)) then
            do p = lbound(pc%part, 1), ubound(pc%part, 1)
                if (pc%part(p)%factored) call sparse_destroy(pc%part(p)%solver)
            end do
            deallocate (pc%part)
        end if
        if (allocated(pc%angle_rdof)) deallocate (pc%angle_rdof)
        if (allocated(pc%angle_diag)) deallocate (pc%angle_diag)
    end subroutine
end module
