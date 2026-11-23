program test_gauge_moments
    !! milca_gauge_moments: moment constraints orthogonal to a continuous
    !! symmetry, on an energy whose symmetry group is known exactly.
    !!
    !! Four points in space with W = sum_(i<j) (|x_i - x_j|^2 - d_ij^2)^2/4,
    !! d_ij the edge lengths of a tetrahedron: W is invariant under the six
    !! rigid motions, the minimum is the orbit of the tetrahedron, and the
    !! Hessian there has six zero eigenvalues. Checked:
    !!
    !! 1. without constraints the factorisation of the Hessian at the
    !!    minimum is singular (six zero pivots, or not reliable);
    !! 2. with the six moment constraints, following the state, Newton
    !!    converges from a deformed start, the matrix is regular and positive
    !!    on the constraints, the multipliers vanish (the symmetry is exact)
    !!    and the edge lengths are reached;
    !! 3. the same with the constraints fixed at the start: the solution
    !!    satisfies them, i.e. the centroid and the first angular moments of
    !!    the start are kept;
    !! 4. a seventh generator that is a combination of two others is dropped
    !!    (a row of zeros, inactive in the solver) and the solve is unchanged;
    !! 5. for collinear points the rotation about their line is no generator:
    !!    five of six are kept;
    !! 6. rows rebuilt at another state keep their pattern.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsparse_csr, only: csr_t
    use fortfem_constrained_newton, only: newton_constraints_t, newton_options_t, &
        newton_report_t, newton_constraints_init, constrained_newton, NEWTON_CONVERGED
    use milca_gauge_moments, only: gauge_moment_constraints
    implicit none
    integer, parameter :: np = 4, n = 3*np
    real(dp), parameter :: tetra(3, np) = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.1_dp, &
        0.0_dp, 0.3_dp, 0.9_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.8_dp], [3, np])
    integer :: generators = 6, failures = 0
    real(dp) :: z0(n), z(n), moments(n, 7), gram(7, 7)
    type(newton_constraints_t) :: c, none
    type(newton_options_t) :: options
    type(newton_report_t) :: report
    logical :: kept(7)
    integer :: status, fixed(0), i
    integer, allocatable :: pattern(:)

    ! a deformed, moved and turned copy of the tetrahedron
    do i = 1, np
        z0(3*i - 2:3*i) = [0.9_dp*tetra(1, i) - 0.3_dp*tetra(2, i) + 0.2_dp, &
                           0.4_dp*tetra(1, i) + 1.1_dp*tetra(2, i) - 0.1_dp, &
                           1.2_dp*tetra(3, i) + 0.05_dp*i]
    end do
    call newton_constraints_init(none)

    ! 1
    z = reshape(tetra, [n])
    call constrained_newton(energy, z, fixed, none, options, report, status)
    call check(status == NEWTON_CONVERGED .and. &
               (report%zero_eigenvalues == 6 .or. .not. report%inertia_reliable), &
               'the Hessian at the minimum is singular along the six rigid motions')

    ! 2
    z = z0
    call rows(z, c)
    call constrained_newton(energy, z, fixed, c, options, report, status, rows)
    call check(status == NEWTON_CONVERGED, 'Newton converges with the moment constraints')
    call check(report%negative_eigenvalues == 0 .and. report%zero_eigenvalues == 0 .and. &
               report%inertia_reliable, 'regular and positive on the constraints')
    call check(report%multiplier_norm < 1e-9_dp, 'the multipliers of an exact symmetry vanish')
    call check(edge_defect(z) < 1e-9_dp, 'the edge lengths of the tetrahedron are reached')
    print '(a,i0,a,es9.2,a,es9.2)', 'following constraints: steps ', report%steps, &
        '  multiplier', report%multiplier_norm, '  edge defect', edge_defect(z)

    ! 3
    z = z0
    call newton_constraints_init(c)
    call rows(z0, c)
    call constrained_newton(energy, z, fixed, c, options, report, status)
    call check(status == NEWTON_CONVERGED .and. edge_defect(z) < 1e-9_dp, &
               'Newton converges with the constraints fixed at the start')
    call generator_moments(z0, moments(:, :6), gram(:6, :6))
    call check(maxval(abs(matmul(transpose(moments(:, :6)), z - z0))) < 1e-10_dp, &
               'the solution differs from the start by something orthogonal to the orbit')

    ! 4
    generators = 7
    z = z0
    call newton_constraints_init(c)
    call rows(z, c)
    call generator_moments(z, moments, gram)
    call gauge_moment_constraints(moments, gram, z, c, kept)
    call check(count(kept) == 6, 'a generator that depends on two others is dropped')
    pattern = c%col
    call constrained_newton(energy, z, fixed, c, options, report, status, rows)
    call check(status == NEWTON_CONVERGED .and. report%zero_eigenvalues == 0 .and. &
               report%negative_eigenvalues == 0 .and. edge_defect(z) < 1e-9_dp, &
               'the solve with the dropped generator as an inactive row')
    ! 6
    call check(size(c%col) == size(pattern) .and. all(c%col == pattern) .and. &
               size(c%target) == 7, 'rebuilt rows keep their pattern')
    generators = 6

    ! 5
    do i = 1, np
        z(3*i - 2:3*i) = 0.3_dp*i*[1.0_dp, 2.0_dp, -1.0_dp] + [0.1_dp, 0.0_dp, 0.2_dp]
    end do
    call generator_moments(z, moments(:, :6), gram(:6, :6))
    call newton_constraints_init(c)
    call gauge_moment_constraints(moments(:, :6), gram(:6, :6), z, c, kept(:6))
    call check(count(kept(:6)) == 5, 'collinear points: the rotation about their line is dropped')

    if (failures > 0) error stop 'gauge moments failed'
    print '(a)', 'PASS gauge moments'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            print '(2a)', 'PASS ', what
        else
            print '(2a)', 'FAIL ', what
            failures = failures + 1
        end if
    end subroutine check

    real(dp) function edge_defect(x) result(d)
        real(dp), intent(in) :: x(:)
        integer :: a, b

        d = 0
        do a = 1, np
            do b = a + 1, np
                d = max(d, abs(norm2(x(3*a - 2:3*a) - x(3*b - 2:3*b)) - &
                               norm2(tetra(:, a) - tetra(:, b))))
            end do
        end do
    end function edge_defect

    subroutine energy(x, w, gradient, hessian)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: w
        real(dp), intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp) :: g(n), h(n, n), d(3), e, blk(3, 3)
        integer :: a, b, p, q, k

        w = 0
        g = 0
        h = 0
        do a = 1, np
            do b = a + 1, np
                d = x(3*a - 2:3*a) - x(3*b - 2:3*b)
                e = sum(d**2) - sum((tetra(:, a) - tetra(:, b))**2)
                w = w + 0.25_dp*e**2
                g(3*a - 2:3*a) = g(3*a - 2:3*a) + e*d
                g(3*b - 2:3*b) = g(3*b - 2:3*b) - e*d
                do p = 1, 3
                    do q = 1, 3
                        blk(p, q) = 2*d(p)*d(q)
                    end do
                    blk(p, p) = blk(p, p) + e
                end do
                h(3*a - 2:3*a, 3*a - 2:3*a) = h(3*a - 2:3*a, 3*a - 2:3*a) + blk
                h(3*b - 2:3*b, 3*b - 2:3*b) = h(3*b - 2:3*b, 3*b - 2:3*b) + blk
                h(3*a - 2:3*a, 3*b - 2:3*b) = h(3*a - 2:3*a, 3*b - 2:3*b) - blk
                h(3*b - 2:3*b, 3*a - 2:3*a) = h(3*b - 2:3*b, 3*a - 2:3*a) - blk
            end do
        end do
        if (present(gradient)) gradient = g
        if (.not. present(hessian)) return
        hessian%nrow = n
        hessian%ncol = n
        allocate (hessian%row_ptr(n + 1), hessian%col(n*(n + 1)/2), hessian%val(n*(n + 1)/2))
        k = 0
        do p = 1, n
            hessian%row_ptr(p) = k + 1
            do q = p, n
                k = k + 1
                hessian%col(k) = q
                hessian%val(k) = h(p, q)
            end do
        end do
        hessian%row_ptr(n + 1) = k + 1
    end subroutine energy

    subroutine generator_moments(x, m, g)
        !! The generators of the rigid motions at x: three translations, three
        !! rotations about the axes and, as a seventh, the sum of the first
        !! translation and the third rotation. With the identity as weight
        !! the moments are the generators themselves.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: m(:, :), g(:, :)
        integer :: a

        m = 0
        do a = 1, np
            m(3*a - 2, 1) = 1
            m(3*a - 1, 2) = 1
            m(3*a, 3) = 1
            m(3*a - 2:3*a, 4) = [0.0_dp, -x(3*a), x(3*a - 1)]
            m(3*a - 2:3*a, 5) = [x(3*a), 0.0_dp, -x(3*a - 2)]
            m(3*a - 2:3*a, 6) = [-x(3*a - 1), x(3*a - 2), 0.0_dp]
        end do
        if (size(m, 2) == 7) m(:, 7) = m(:, 1) + m(:, 6)
        g = matmul(transpose(m), m)
    end subroutine generator_moments

    subroutine rows(x, constraints)
        !! The constraints at the state x (the update of the Newton solver).
        real(dp), intent(in) :: x(:)
        type(newton_constraints_t), intent(inout) :: constraints
        real(dp) :: m(n, generators), g(generators, generators)
        logical :: keep(generators)

        call generator_moments(x, m, g)
        call gauge_moment_constraints(m, g, x, constraints, keep)
    end subroutine rows

end program test_gauge_moments
