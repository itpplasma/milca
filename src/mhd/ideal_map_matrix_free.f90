module milca_ideal_map_matrix_free
    ! Matrix-free damped Newton for the map energy (position gauge), the
    ! counterpart of ideal_map_solve without the assembled, mode-coupled
    ! Hessian. Per iteration the point Hessians are stored
    ! (ideal_map_store_hessians), and each trial solves
    ! (H + mu diag|H|) s = -g inexactly (forcing min(0.5, sqrt(|g|/|g0|))) by
    ! deflated conjugate gradients with the overlapping Schwarz preconditioner
    ! in n (milca_ideal_map_mode_pc) and the Ritz vectors of the previous
    ! solve. Non-positive curvature in CG raises mu like a failed Cholesky in
    ! the direct path. The globalisation is that of ideal_map_solve: the
    ! displacement cap of 5 % of the reference extent, acceptance by a tenth of
    ! the predicted reduction -(g.s + s.H s/2) outside the roundoff band of W
    ! and by the trapezoidal energy change inside it, and the damping update of
    ! Nielsen (1999). Memory: the point Hessians, O(cells nq nphi), and the
    ! Schwarz factors, O(nodes nmode).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve_report_t, ideal_map_energy, &
        ideal_map_solve, &
        ideal_map_reduce, ideal_map_expand, ideal_map_hvp, ideal_map_hessian_store_t, &
        ideal_map_store_hessians, ideal_map_rdof, ideal_map_angle_rdof, &
        ideal_map_angle_diagonal, ideal_map_dof, IDEAL_MAP_OK, IDEAL_MAP_BAD_INPUT, &
        IDEAL_MAP_FOLDED, IDEAL_MAP_NOT_CONVERGED, IDEAL_MAP_SOLVE_FAILED
    use milca_ideal_map_mode_pc, only: ideal_map_mode_pc_t, ideal_map_mode_pc_build, &
        ideal_map_mode_pc_apply, ideal_map_mode_pc_destroy, ideal_map_mode_pc_refactor
    use milca_ideal_map_krylov, only: deflated_pcg
    implicit none
    private
    public :: ideal_map_solve_matrix_free, ideal_map_hessian_diagonal
    public :: ideal_map_hessian_nonzeros, ideal_map_solve_auto
    public :: IDEAL_MAP_DIRECT_BYTES_PER_NONZERO
    ! peak resident bytes of ideal_map_solve per nonzero of the reduced Hessian:
    ! slope of 70, 81, 124, 163 MB at 0.42, 0.67, 1.12, 1.53 million nonzeros
    ! (Landreman P3, 8 to 24 nodes) is 84 B
    real(dp), parameter :: IDEAL_MAP_DIRECT_BYTES_PER_NONZERO = 90
    ! recycled Ritz directions kept across Newton steps
    integer, parameter :: max_recycled = 40
contains

    subroutine ideal_map_solve_matrix_free(map, tolerance, max_iterations, report, status, &
            reference_gradient, verbose)
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        type(ideal_map_solve_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: reference_gradient
        logical, intent(in), optional :: verbose
        type(ideal_map_hessian_store_t) :: store
        type(ideal_map_mode_pc_t) :: pc
        real(dp), allocatable :: grad(:), gred(:), diag(:), shift(:), step(:), full(:), &
            trial(:, :, :), hs(:), rhs(:), w(:, :), ritz(:, :), trial_angle(:)
        real(dp) :: energy, trial_energy, g0, gnorm, ginit, mu, scale, alpha, cap, eta, pred
        real(dp) :: nu, ratio
        real(dp) :: t0
        integer :: it, tries, estatus, nred, its, nc
        logical :: accepted, negative, talk
        status = IDEAL_MAP_BAD_INPUT
        talk = .false.
        if (present(verbose)) talk = verbose
        if (.not. map%position_gauge) return
        nred = map%nred
        nc = 3*map%nmode*map%nodes
        allocate (grad(map%ndof), gred(nred), diag(nred), shift(nred), step(nred), &
            full(map%ndof), hs(nred), rhs(nred), w(nred, 0), trial_angle(map%nangle))
        if (map%nangle > 0) trial_angle = map%angle
        call ideal_map_energy(map, map%x, energy, estatus, grad)
        if (estatus /= IDEAL_MAP_OK) then
            status = IDEAL_MAP_FOLDED; return
        end if
        call ideal_map_reduce(map, grad, gred)
        where (map%fixed) gred = 0
        gnorm = norm2(gred); ginit = gnorm; g0 = gnorm
        report%initial_gradient = gnorm
        if (present(reference_gradient)) g0 = reference_gradient
        mu = max(map%initial_damping, 1.0e-14_dp)
        cap = 0.05_dp*maxval(maxval(map%reference, dim=2) - minval(map%reference, dim=2))
        status = IDEAL_MAP_NOT_CONVERGED
        report%energy = energy; report%gradient = gnorm
        do it = 1, max_iterations
            if (.not. gnorm > tolerance*g0) then
                status = IDEAL_MAP_OK; exit
            end if
            report%iterations = it
            t0 = wall()
            call ideal_map_store_hessians(map, map%x, store)
            call ideal_map_hessian_diagonal(map, store, diag)
            report%hessian_s = report%hessian_s + wall() - t0
            report%hessians = report%hessians + 1
            scale = maxval(abs(diag), mask=.not. map%fixed)
            eta = min(0.5_dp, sqrt(gnorm/max(ginit, tiny(1.0_dp))))
            accepted = .false.
            nu = 2; ratio = 1
            do tries = 1, 30
                report%trials = report%trials + 1
                t0 = wall()
                shift = merge(mu*(abs(diag) + 1.0e-12_dp*scale), 0.0_dp, .not. map%fixed)
                ! one assembly of the subdomain matrices per Newton iteration
                if (tries == 1) then
                    call ideal_map_mode_pc_build(map, store, mu, pc, estatus)
                else
                    call ideal_map_mode_pc_refactor(pc, mu, estatus)
                end if
                if (estatus /= IDEAL_MAP_OK) then
                    report%factor_s = report%factor_s + wall() - t0
                    mu = mu*nu; nu = 2*nu; cycle
                end if
                report%factorizations = report%factorizations + max(1, map%ntor)
                rhs = -gred
                call deflated_pcg(operator, preconditioner, rhs, eta, 500, w, step, its, &
                    negative, max_recycled, ritz)
                report%krylov_solves = report%krylov_solves + 1
                report%krylov_iterations = report%krylov_iterations + its
                report%factor_s = report%factor_s + wall() - t0
                if (negative) then
                    report%indefinite_trials = report%indefinite_trials + 1
                    mu = mu*nu; nu = 2*nu; cycle
                end if
                where (map%fixed) step = 0
                call ideal_map_expand(map, step, full)
                alpha = 1
                if (maxval(abs(full(:nc))) > cap) alpha = cap/maxval(abs(full(:nc)))
                if (map%nangle > 0) alpha = min(alpha, 1/max(tiny(1.0_dp), &
                    maxval(abs(full(nc + 1:)))))
                full = alpha*full; step = alpha*step
                trial = map%x + reshape(full(:nc), shape(map%x))
                if (map%nangle > 0) trial_angle = map%angle + full(nc + 1:)
                t0 = wall()
                call ideal_map_energy(map, trial, trial_energy, estatus, angle=trial_angle)
                report%energy_s = report%energy_s + wall() - t0
                if (estatus == IDEAL_MAP_OK .and. trial_energy < energy &
                    - 1024*epsilon(1.0_dp)*abs(energy)) then
                    call ideal_map_hvp(map, map%x, step, hs, estatus, store)
                    pred = -dot_product(gred, step) - dot_product(step, hs)/2
                    ratio = (energy - trial_energy)/pred
                    if (ratio >= 0.1_dp) then
                        accepted = .true.; exit
                    end if
                end if
                if (estatus == IDEAL_MAP_OK .and. abs(trial_energy - energy) &
                    <= 1024*epsilon(1.0_dp)*abs(energy)) then
                    call ideal_map_energy(map, trial, trial_energy, estatus, grad, &
                        angle=trial_angle)
                    call ideal_map_reduce(map, grad, rhs)
                    where (map%fixed) rhs = 0
                    if (estatus == IDEAL_MAP_OK .and. dot_product(gred + rhs, step) <= 0) then
                        ratio = 1; accepted = .true.; exit
                    end if
                end if
                mu = mu*nu; nu = 2*nu
            end do
            if (.not. accepted) exit
            report%step = maxval(abs(trial - map%x))
            map%x = trial
            if (map%nangle > 0) map%angle = trial_angle
            call ideal_map_energy(map, map%x, energy, estatus, grad)
            call ideal_map_reduce(map, grad, gred)
            where (map%fixed) gred = 0
            gnorm = norm2(gred)
            call recycle(ritz)
            report%energy = energy; report%gradient = gnorm/g0; report%damping = mu
            if (talk) print '(a,i4,es14.6,es11.3,es10.2,i4,i5,es10.2)', &
                'matrix-free newton', it, energy, gnorm/g0, mu, tries, its, alpha
            if (.not. alpha < 1) mu = max(mu*max(1.0_dp/3, 1 - (2*ratio - 1)**3), 1.0e-14_dp)
            if (report%step <= 1.0e-13_dp*maxval(abs(map%reference))) then
                status = IDEAL_MAP_OK; exit
            end if
        end do
        if (.not. gnorm > tolerance*g0) status = IDEAL_MAP_OK
        call ideal_map_mode_pc_destroy(pc)
    contains
        subroutine operator(v, y)
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: y(:)
            integer :: st
            call ideal_map_hvp(map, map%x, v, y, st, store)
            y = y + shift*v
            where (map%fixed) y = v
        end subroutine

        subroutine preconditioner(v, y)
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: y(:)
            integer :: st
            call ideal_map_mode_pc_apply(map, pc, v, y, st)
        end subroutine

        subroutine recycle(new)
            ! W <- orthonormal basis of [new, W], at most max_recycled columns
            real(dp), intent(in) :: new(:, :)
            real(dp), allocatable :: v(:, :)
            integer :: i, j, k
            v = reshape([new, w], [nred, size(new, 2) + size(w, 2)])
            where (spread(map%fixed, 2, size(v, 2))) v = 0
            k = 0
            do j = 1, size(v, 2)
                if (k == max_recycled) exit
                do i = 1, k
                    v(:, j) = v(:, j) - dot_product(v(:, i), v(:, j))*v(:, i)
                end do
                if (.not. norm2(v(:, j)) > 1.0e-8_dp) cycle
                k = k + 1
                v(:, k) = v(:, j)/norm2(v(:, j))
            end do
            w = v(:, :k)
        end subroutine

        real(dp) function wall()
            integer(int64) :: count, rate
            call system_clock(count, rate)
            wall = real(count, dp)/rate
        end function
    end subroutine

    subroutine ideal_map_hessian_diagonal(map, store, diag)
        ! Diagonal of the reduced Hessian (position gauge) from the stored
        ! point Hessians: sum w phi^T d2w/dv2 phi for the test function of every
        ! (component, mode, node), with the angle block's diagonal.
        type(ideal_map_t), intent(in) :: map
        type(ideal_map_hessian_store_t), intent(in) :: store
        real(dp), intent(out) :: diag(:)
        real(dp), allocatable :: part(:)
        real(dp) :: f(0:3), weight, dphi, h(4, 4)
        integer :: e, g, j, a, m, c, r, rp, row
        integer :: idx(4, 3)
        diag = 0
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        ! v index of test function r = 0..3 and component c: 3 r + c
        do r = 0, 3
            do c = 1, 3
                idx(r + 1, c) = 3*r + c
            end do
        end do
        !$omp parallel default(shared) private(e, g, j, a, m, c, r, rp, row, f, weight, &
        !$omp h, part)
        allocate (part(size(diag)))
        part = 0
        !$omp do schedule(dynamic)
        do e = 1, map%cells
            do j = 1, map%nphi
                do g = 1, map%nq
                    weight = map%weight(g)*dphi
                    do c = 1, 2
                        do r = 1, 4
                            do rp = 1, 4
                                h(r, rp) = unpacked(store%hp(:, g, j, e), idx(r, c), &
                                    idx(rp, c))
                            end do
                        end do
                        do a = 1, map%nloc
                            do m = 1, map%nmode
                                f(0) = map%shape(a, g)*map%trig(m, j)
                                f(1) = map%dshape(1, a, g)*map%trig(m, j)
                                f(2) = map%dshape(2, a, g)*map%trig(m, j)
                                f(3) = map%shape(a, g)*map%dtrig(m, j)
                                row = ideal_map_rdof(map, c, m, map%conn(a, e))
                                part(row) = part(row) + weight*dot_product(f, matmul(h, f))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp critical
        diag = diag + part
        !$omp end critical
        deallocate (part)
        !$omp end parallel
        if (map%nangle > 0) then
            block
                real(dp) :: angle_diag(map%nangle)
                integer :: k
                call ideal_map_angle_diagonal(map, map%x, angle_diag)
                do k = 1, map%nangle
                    diag(ideal_map_angle_rdof(map, k)) = angle_diag(k)
                end do
            end block
        end if
    contains
        pure real(dp) function unpacked(hp, i, k)
            ! entry (i, k) of the symmetric 12 x 12 matrix packed as hp
            real(dp), intent(in) :: hp(78)
            integer, intent(in) :: i, k
            integer :: lo, hi
            lo = min(i, k); hi = max(i, k)
            if (hi <= 9) then
                unpacked = hp(lo + hi*(hi - 1)/2)
            else if (lo <= 9) then
                unpacked = hp(45 + lo + 9*(hi - 10))
            else
                unpacked = hp(72 + (lo - 9) + (hi - 10)*(hi - 9)/2)
            end if
        end function
    end subroutine

    integer(int64) function ideal_map_hessian_nonzeros(map) result(nnz)
        ! Nonzeros of the reduced Hessian that ideal_map_solve assembles: the
        ! (2 nmode)^2 block of every pair of nodes sharing a cell, the angle
        ! columns and rows (dense against all dofs) and the angle block.
        type(ideal_map_t), intent(in) :: map
        integer, allocatable :: first(:), next(:), cell_of(:), mark(:)
        integer(int64) :: pairs
        integer :: a, e, i, item
        allocate (first(map%nodes), next(map%cells*map%nloc), cell_of(map%cells*map%nloc), &
            mark(map%nodes))
        first = 0; item = 0
        do e = 1, map%cells
            do i = 1, map%nloc
                item = item + 1
                cell_of(item) = e
                next(item) = first(map%conn(i, e))
                first(map%conn(i, e)) = item
            end do
        end do
        mark = 0; pairs = 0
        do a = 1, map%nodes
            item = first(a)
            do while (item > 0)
                e = cell_of(item); item = next(item)
                do i = 1, map%nloc
                    if (mark(map%conn(i, e)) == a) cycle
                    mark(map%conn(i, e)) = a
                    pairs = pairs + 1
                end do
            end do
        end do
        nnz = pairs*int(2*map%nmode, int64)**2 &
            + 2*int(map%nangle, int64)*(map%nred - map%nangle) + int(map%nangle, int64)**2
    end function

    subroutine ideal_map_solve_auto(map, tolerance, max_iterations, report, status, &
            memory_budget, matrix_free, reference_gradient, verbose)
        ! ideal_map_solve when its estimated peak memory,
        ! IDEAL_MAP_DIRECT_BYTES_PER_NONZERO times the Hessian nonzeros, fits
        ! memory_budget (bytes), ideal_map_solve_matrix_free otherwise. The
        ! matrix-free path needs the position gauge, without it the direct path
        ! runs. matrix_free tells which path ran.
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: tolerance, memory_budget
        integer, intent(in) :: max_iterations
        type(ideal_map_solve_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(out) :: matrix_free
        real(dp), intent(in), optional :: reference_gradient
        logical, intent(in), optional :: verbose
        matrix_free = map%position_gauge .and. IDEAL_MAP_DIRECT_BYTES_PER_NONZERO &
            *real(ideal_map_hessian_nonzeros(map), dp) > memory_budget
        if (matrix_free) then
            call ideal_map_solve_matrix_free(map, tolerance, max_iterations, report, status, &
                reference_gradient, verbose)
        else
            call ideal_map_solve(map, tolerance, max_iterations, report, status, verbose, &
                reference_gradient=reference_gradient)
        end if
    end subroutine
end module
