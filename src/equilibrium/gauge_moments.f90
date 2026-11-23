module milca_gauge_moments
    !! Gauge by constraints (docs/flux-map.md, "Gauge constraint, one solver
    !! loop and continuation"): moment constraints orthogonal to a continuous
    !! symmetry.
    !!
    !! An energy W(z) that is invariant under a continuous group has, at a
    !! state z, the tangents g_j(z) of the group orbit (the generators) as
    !! directions of no restoring force: at a stationary point they are null
    !! vectors of the Hessian, and Newton creeps along them unless they are
    !! removed. They are removed here by one linear constraint per generator,
    !!
    !!     <g_j(z_ref), z - z_ref>_M = 0,   j = 1 .. m,                    (1)
    !!
    !! with M a positive weight (a mass matrix): the state may differ from
    !! the reference z_ref only by something M-orthogonal to the orbit through
    !! z_ref. (1) selects one point of every orbit near z_ref as long as the
    !! Gram matrix <g_i(z_ref), g_j(z)>_M is regular (the constraints are
    !! transverse to the orbits); it costs no energy and changes no physical
    !! quantity, and its multipliers vanish at a stationary point when the
    !! symmetry is exact for the discrete energy. Unlike a penalty or pinned
    !! rows it needs no weight, no symmetric mesh and no symmetric state.
    !!
    !! The caller supplies the moments, i.e. the rows z -> <g_j(z_ref), z>_M
    !! as dense vectors, and their Gram matrix with the generators; this
    !! module drops generators that are dependent at the resolution of the
    !! moments (pivoted Cholesky of the Gram matrix scaled to unit diagonal)
    !! and returns the rows of (1) for the Newton solver
    !! (fortfem_constrained_newton). With the rows rebuilt at every iterate
    !! (z_ref the iterate) every step is M-orthogonal to the orbit through the
    !! current state. User: the rotation gauge of the flux map (milca_flux_map,
    !! milca_flux_map_3d); test_gauge_moments checks the module on the rigid
    !! motions of a point cloud.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_constrained_newton, only: newton_constraints_t, newton_constraints_init, &
        newton_constraints_add
    implicit none
    private

    public :: gauge_moment_constraints

contains

    subroutine gauge_moment_constraints(moments, gram, z_ref, constraints, kept, tolerance)
        !! moments(:, j): the row of generator j; gram(i, j) = <g_i, g_j>_M at
        !! z_ref (symmetric, positive semidefinite). constraints: on entry
        !! empty, or the rows of an earlier call with the same number of
        !! generators, whose pattern is then kept (rows that follow the
        !! state); on return the rows of
        !! (1), one per generator on all slots, each scaled to largest entry
        !! one; kept(j) tells which are kept, the others are rows of zeros
        !! (inactive in the solver). A generator is dropped when its pivot in the
        !! Cholesky factorisation of the scaled Gram matrix is below tolerance
        !! (default 1e-8): it is, to that accuracy, a combination of kept ones.
        real(dp), intent(in) :: moments(:, :), gram(:, :), z_ref(:)
        type(newton_constraints_t), intent(inout) :: constraints
        logical, intent(out) :: kept(:)
        real(dp), intent(in), optional :: tolerance
        real(dp) :: a(size(gram, 1), size(gram, 1)), d(size(gram, 1)), tol, big
        integer, allocatable :: cols(:)
        integer :: m, j, k, i, n

        m = size(gram, 1)
        n = size(moments, 1)
        tol = 1.0e-8_dp
        if (present(tolerance)) tol = tolerance
        if (.not. allocated(constraints%target)) call newton_constraints_init(constraints)
        do j = 1, m
            d(j) = sqrt(max(gram(j, j), 0.0_dp))
        end do
        kept = d > 0
        do j = 1, m
            do i = 1, m
                a(i, j) = 0
                if (kept(i) .and. kept(j)) a(i, j) = gram(i, j)/(d(i)*d(j))
            end do
        end do
        ! Cholesky with diagonal pivoting; the remaining diagonal is the
        ! squared distance of a generator from the span of the kept ones
        d = [(a(j, j), j=1, m)]
        kept = .false.
        do k = 1, m
            j = maxloc(d, 1)
            if (.not. d(j) > tol) exit
            kept(j) = .true.
            a(:, j) = a(:, j)/sqrt(d(j))
            do i = 1, m
                if (kept(i)) cycle
                a(:, i) = a(:, i) - a(i, j)*a(:, j)
                d(i) = a(i, i)
            end do
            d(j) = -1
        end do
        if (size(constraints%target) == m) then
            ! rows rebuilt at another state keep their pattern
            do j = 1, m
                associate (c => constraints%col(constraints%ptr(j):constraints%ptr(j + 1) - 1), &
                           v => constraints%val(constraints%ptr(j):constraints%ptr(j + 1) - 1))
                    v = 0
                    constraints%target(j) = 0
                    if (.not. kept(j)) cycle
                    big = maxval(abs(moments(:, j)))
                    v = moments(c, j)/big
                    constraints%target(j) = dot_product(v, z_ref(c))
                end associate
            end do
            return
        end if
        ! the pattern: the slots on which some moment does not vanish; a
        ! dropped generator is a row of zeros
        call newton_constraints_init(constraints)
        cols = pack([(i, i=1, n)], any(moments /= 0, 2))
        do j = 1, m
            if (kept(j)) then
                big = maxval(abs(moments(:, j)))
                call newton_constraints_add(constraints, cols, moments(cols, j)/big, &
                                            dot_product(moments(cols, j), z_ref(cols))/big)
            else
                call newton_constraints_add(constraints, cols, 0*z_ref(cols), 0.0_dp)
            end if
        end do
    end subroutine gauge_moment_constraints

end module milca_gauge_moments
