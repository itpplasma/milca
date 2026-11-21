module milca_inertia_newton
    !! Inertia-controlled Newton steps on a symmetric (KKT) system, written
    !! against FortSparse only so that it can move into FortSparse unchanged.
    !! inertia_newton_factor factors K + delta diag(shift) with LDL^T after a
    !! positive diagonal congruence (units), reads the inertia and raises
    !! delta until it equals the target (n_pos, n_neg, 0), with the IPOPT
    !! update (Waechter and Biegler, Math. Program. 106 (2006) 25, sec. 3.1):
    !! delta = 0 first, then max(delta_min, delta_last/3), or delta_first when the
    !! last step needed none, then delta_factor times per rejection. shift
    !! selects the rows that may be regularised (primal blocks) with their
    !! weight; rows with shift 0 (multipliers, a maximised block) keep their
    !! exact values. inertia_newton_accept is the ratio test of a merit:
    !! accepted when the actual decrease is at least eta_accept of the
    !! predicted one, and inside the roundoff band of the merit when the
    !! stationarity norm decreases. inertia_watchdog_t is the watchdog of
    !! Chamberlain, Powell, Lemarechal and Pedersen (Math. Program. Stud. 16
    !! (1982) 1; Nocedal and Wright, Numerical Optimization, sec. 15.6): a full
    !! step the ratio test rejects is taken anyway when the merit has to cross
    !! a small barrier, and the iteration returns to the checkpoint unless the
    !! merit falls below it within max_steps further steps.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse, only: sparse_solver_t, sparse_refactor, sparse_solve, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_AMD, sparse_inertia, &
        FORTSPARSE_INERTIA_RELIABLE, FORTSPARSE_INERTIA_UNCERTAIN, csc_normalize_symmetric_raw
    implicit none
    private
    public :: inertia_newton_t, inertia_newton_factor, inertia_newton_solve, &
        inertia_newton_accept, inertia_watchdog_t, inertia_watchdog_start, &
        inertia_watchdog_check
    integer, parameter, public :: INERTIA_NEWTON_OK = 0, INERTIA_NEWTON_FAILED = 1
    integer, parameter, public :: WATCHDOG_IDLE = 0, WATCHDOG_COMMIT = 1, &
        WATCHDOG_RESTORE = 2
    type :: inertia_newton_t
        real(dp) :: delta_min = 1.0e-20_dp, delta_first = 1.0e-4_dp
        real(dp) :: delta_factor = 8, delta_max = 1.0e20_dp
        real(dp) :: eta_accept = 0.1_dp, roundoff_band = 1024
        ! last delta used, for the first trial of the next factorisation
        real(dp) :: delta_last = 0
        ! congruence scaling of the last factor (solve reads it)
        real(dp), allocatable :: scale(:)
        ! (positive, negative, zero) of the last factor
        integer :: inertia(3) = -1
    end type
    type :: inertia_watchdog_t
        !> further steps allowed to fall below the checkpoint, steps without a
        !! new watch after a return, and the decrease required (fraction of the
        !! predicted decrease of the tentative step)
        integer :: max_steps = 3, cooldown = 5
        real(dp) :: eta = 0.25_dp
        !> > 0: steps since the tentative step; < 0: cooldown left
        integer :: steps = 0
        real(dp) :: merit = 0, predicted = 0
        real(dp), allocatable :: state(:), gradient(:)
    end type
contains
    subroutine inertia_newton_factor(opt, rowptr, colidx, val, shift, n_pos, n_neg, &
                                     solver, delta, status, minimum, any_inertia)
        !! K in full-storage CSR (symmetric, so also its CSC). On success solver
        !! holds the factor of D (K + delta diag(shift)) D. minimum: smallest
        !! delta to try (a step rejected by the ratio test passes the larger
        !! delta it needs).
        type(inertia_newton_t), intent(inout) :: opt
        integer, intent(in) :: rowptr(:), colidx(:), n_pos, n_neg
        real(dp), intent(in) :: val(:), shift(:)
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), intent(out) :: delta
        integer, intent(out) :: status
        real(dp), intent(in), optional :: minimum
        !> factor once at delta = 0 whatever its inertia (a Newton step near
        !! the solution, accepted by the caller on the stationarity norm)
        logical, intent(in), optional :: any_inertia
        real(dp), allocatable :: scaled(:)
        type(fortsparse_status_t) :: sst
        integer(int64) :: signs(3)
        real(dp) :: correction
        integer :: n, i, p, reliability
        logical :: factored
        status = INERTIA_NEWTON_FAILED
        n = size(rowptr) - 1
        if (size(shift) /= n .or. n_pos + n_neg /= n) return
        if (allocated(opt%scale)) deallocate (opt%scale)
        allocate (opt%scale(n), scaled(size(val)))
        do i = 1, n
            opt%scale(i) = maxval(abs(val(rowptr(i):rowptr(i + 1) - 1)))
            if (abs(shift(i)) > 0) opt%scale(i) = max(opt%scale(i), shift(i))
        end do
        if (any(.not. opt%scale > 0)) return
        opt%scale = 1/sqrt(opt%scale)
        ! delta = 0 first (the Newton step whenever the inertia allows it),
        ! then from max(delta_min, delta_last/3)
        delta = 0
        if (present(minimum)) delta = minimum
        factored = .false.
        do
            do i = 1, n
                do p = rowptr(i), rowptr(i + 1) - 1
                    scaled(p) = val(p)
                    if (colidx(p) == i) scaled(p) = scaled(p) + delta*shift(i)
                    scaled(p) = scaled(p)*opt%scale(i)*opt%scale(colidx(p))
                end do
            end do
            ! LDL^T reads an exactly symmetric matrix: assembly roundoff is
            ! averaged out (rejected above a relative 1e-8).
            call csc_normalize_symmetric_raw(n, rowptr, colidx, scaled, 1.0e-8_dp, &
                correction, sst)
            if (.not. status_ok(sst)) then
                print '(2a)', 'inertia_newton_factor: ', trim(sst%msg)
                return
            end if
            ! multifrontal LDL^T; the same pattern over shifts and Newton steps
            ! reuses its ordering and symbolic analysis
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            call sparse_refactor(solver, n, n, size(scaled), rowptr, colidx, scaled, sst)
            factored = status_ok(sst)
            opt%inertia = -1
            if (factored) then
                call sparse_inertia(solver, signs(1), signs(2), signs(3), reliability, &
                                    sst)
                ! An uncertain count (a pivot within roundoff of zero) is
                ! accepted only when it equals the target exactly.
                if (status_ok(sst) .and. (reliability == FORTSPARSE_INERTIA_RELIABLE &
                    .or. reliability == FORTSPARSE_INERTIA_UNCERTAIN)) &
                    opt%inertia = int(signs([3, 1, 2]))
            end if
            if (opt%inertia(1) == n_pos .and. opt%inertia(2) == n_neg .and. &
                opt%inertia(3) == 0) exit
            if (present(any_inertia)) then
                if (any_inertia .and. factored .and. opt%inertia(3) == 0) exit
            end if
            if (.not. delta > 0 .and. opt%delta_last > 0) then
                delta = max(opt%delta_min, opt%delta_last/3)
            else if (.not. delta > 0) then
                delta = opt%delta_first
            else
                delta = opt%delta_factor*delta
            end if
            if (delta > opt%delta_max) return
        end do
        opt%delta_last = delta
        status = INERTIA_NEWTON_OK
    end subroutine

    subroutine inertia_newton_solve(opt, solver, b, x, status)
        !! x = (K + delta diag(shift))^-1 b with the factor of the last
        !! inertia_newton_factor.
        type(inertia_newton_t), intent(in) :: opt
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        type(fortsparse_status_t) :: sst
        status = INERTIA_NEWTON_FAILED
        call sparse_solve(solver, opt%scale*b, x, sst)
        if (.not. status_ok(sst)) return
        x = opt%scale*x
        if (all(ieee_is_finite(x))) status = INERTIA_NEWTON_OK
    end subroutine

    logical function inertia_newton_accept(opt, merit, merit_trial, predicted, &
                                           gnorm, gnorm_trial) result(accept)
        !! Ratio test: (merit - merit_trial) >= eta_accept predicted, predicted
        !! > 0. When the predicted change is inside roundoff_band ulps of the
        !! merit, the merit cannot resolve it and the stationarity norm decides.
        type(inertia_newton_t), intent(in) :: opt
        real(dp), intent(in) :: merit, merit_trial, predicted, gnorm, gnorm_trial
        real(dp) :: band
        accept = .false.
        if (.not. ieee_is_finite(merit_trial)) return
        band = opt%roundoff_band*epsilon(1.0_dp)*max(abs(merit), abs(merit_trial))
        if (predicted <= band) then
            accept = gnorm_trial < gnorm
            return
        end if
        if (.not. predicted > 0) return
        accept = merit - merit_trial >= opt%eta_accept*predicted
    end function

    logical function inertia_watchdog_start(wd, merit, predicted, state, gradient) &
            result(start)
        !! .true.: take the rejected full step tentatively; the checkpoint
        !! (merit, state, gradient) is kept for a return.
        type(inertia_watchdog_t), intent(inout) :: wd
        real(dp), intent(in) :: merit, predicted, state(:), gradient(:)
        start = wd%steps == 0 .and. predicted > 0
        if (.not. start) return
        wd%steps = 1
        wd%merit = merit; wd%predicted = predicted
        wd%state = state; wd%gradient = gradient
    end function

    integer function inertia_watchdog_check(wd, merit) result(action)
        !! After each accepted step that follows the tentative one: COMMIT
        !! once the merit is below the checkpoint by eta predicted, RESTORE
        !! (to wd%state, wd%gradient, wd%merit) after max_steps without.
        type(inertia_watchdog_t), intent(inout) :: wd
        real(dp), intent(in) :: merit
        action = WATCHDOG_IDLE
        if (wd%steps < 0) wd%steps = wd%steps + 1
        if (wd%steps <= 0) return
        if (merit <= wd%merit - wd%eta*wd%predicted) then
            wd%steps = 0
            action = WATCHDOG_COMMIT
        else if (wd%steps >= wd%max_steps) then
            wd%steps = -wd%cooldown
            action = WATCHDOG_RESTORE
        else
            wd%steps = wd%steps + 1
        end if
    end function
end module milca_inertia_newton
