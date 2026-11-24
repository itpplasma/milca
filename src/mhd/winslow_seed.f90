module milca_winslow_seed
    ! Seed of the nested map by the Winslow map of each toroidal plane: the
    ! plane coordinates X = (q1, q2) of the interior nodes minimise
    ! E = int |A|_F^2 / det A dy, A = dX/dy on the reference disk, which is the
    ! Dirichlet energy of the inverse map y(X). Its minimiser is injective for
    ! any simply connected cross-section, so it seeds boundaries whose harmonic
    ! extension folds. A folded start is untangled with the regularised
    ! determinant (det + sqrt(det^2 + eps^2))/2 and eps -> 0 (Escobar et al.,
    ! Comput. Methods Appl. Mech. Engrg. 192, 2775 (2003)). The Winslow map
    ! exists exactly when the section is bounded by a Jordan curve (Rado-
    ! Kneser-Choquet applied to the harmonic inverse), so a plane is admitted
    ! only if its boundary polygon is simple. Boundary nodes keep their
    ! positions; the plane solutions are projected onto the map's toroidal
    ! modes by least squares.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_map, only: ideal_map_t
    use milca_triangle_basis, only: triangle_basis_jet
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use milca_winslow_seed_leaves, only: winslow_density, winslow_gradient, &
        winslow_hessian, winslow_folded_density, winslow_folded_gradient, &
        winslow_folded_hessian
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_GENERAL, FORTSPARSE_ORDER_AMD
    implicit none
    private
    public :: winslow_seed, winslow_plane_energy, WINSLOW_OK, WINSLOW_FAILED, &
        WINSLOW_NOT_SIMPLE

    ! WINSLOW_NOT_SIMPLE: a section boundary is not a simple closed polygon,
    ! so no injective map exists
    integer, parameter :: WINSLOW_OK = 0, WINSLOW_FAILED = 1, WINSLOW_NOT_SIMPLE = 2

    type :: plane_t
        ! reference gradients dny(2, nloc, nq, cells) and weights wq(nq, cells)
        real(dp), allocatable :: dny(:, :, :, :), wq(:, :)
        integer :: nq = 0
        ! unknown index of (component, node), 0 for fixed boundary nodes
        integer, allocatable :: dof(:, :)
        integer :: n = 0
        integer, allocatable :: rowptr(:), colidx(:)
        ! boundary nodes in the order of their reference angle
        integer, allocatable :: ring(:)
        ! reference disk area
        real(dp) :: area = 0
    end type
contains

    subroutine winslow_seed(map, status, tolerance)
        ! Replace map%x(1:2, :, :) by the plane-wise Winslow maps; status
        ! WINSLOW_FAILED when a plane does not reach a map with positive
        ! determinants and a gradient below tolerance times the energy per
        ! length of the section (default 1e-10).
        type(ideal_map_t), intent(inout) :: map
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tolerance
        type(plane_t) :: plane
        real(dp), allocatable :: planes(:, :, :), normal(:, :), rhs(:, :)
        real(dp) :: tol, orientation, area
        integer :: j, m, a, info, plane_status(map%nphi)
        integer, allocatable :: pivot(:)
        status = WINSLOW_FAILED
        tol = 1e-10_dp
        if (present(tolerance)) tol = tolerance
        call plane_setup(map, plane)
        allocate (planes(2, map%nodes, map%nphi))
        !$omp parallel do default(shared) private(a, area, orientation) schedule(static) if(map%nphi > 1)
        do j = 1, map%nphi
            associate (x => planes(:, :, j))
                do a = 1, map%nodes
                    x(:, a) = matmul(map%x(1:2, :, a), map%trig(:, j))
                end do
                area = ring_area(plane, x)
                if (.not. abs(area) > 0 .or. .not. ring_simple(plane, x)) then
                    plane_status(j) = WINSLOW_NOT_SIMPLE
                    cycle
                end if
                ! a clockwise boundary has det A < 0: solve the reflected plane
                orientation = sign(1.0_dp, area)
                x(2, :) = orientation*x(2, :)
                call plane_solve(map, plane, x, tol, plane_status(j))
                x(2, :) = orientation*x(2, :)
            end associate
        end do
        !$omp end parallel do
        do j = 1, map%nphi
            if (plane_status(j) == WINSLOW_OK) cycle
            status = plane_status(j)
            return
        end do
        ! least squares on the plane values: (T T^t) c = T X
        allocate (normal(map%nmode, map%nmode), rhs(map%nmode, 2*map%nodes), &
            pivot(map%nmode))
        normal = matmul(map%trig, transpose(map%trig))
        do a = 1, map%nodes
            do m = 1, 2
                rhs(:, 2*(a - 1) + m) = matmul(map%trig, planes(m, a, :))
            end do
        end do
        call dgesv(map%nmode, 2*map%nodes, normal, map%nmode, pivot, rhs, map%nmode, info)
        if (info /= 0) return
        do a = 1, map%nodes
            if (map%boundary(a)) cycle
            do m = 1, 2
                map%x(m, :, a) = rhs(:, 2*(a - 1) + m)
            end do
        end do
        status = WINSLOW_OK
    end subroutine

    subroutine winslow_plane_energy(map, x, eps, energy, min_det)
        ! E and the smallest determinant of A of plane coordinates x(2, nodes)
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :), eps
        real(dp), intent(out) :: energy, min_det
        type(plane_t) :: plane
        call plane_setup(map, plane)
        call plane_energy(map, plane, x, eps, energy, min_det)
    end subroutine

    subroutine plane_setup(map, plane)
        ! The energy grows like 1/det A, so its quadrature weights must be
        ! positive: a negative weight (the centroid of the map's symmetric
        ! degree-7 rule) rewards a vanishing determinant there, the energy is
        ! unbounded below and the untangling diverges. The plane therefore
        ! takes a collapsed Gauss-Legendre rule of the map's degree, 2 p + 1.
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(out) :: plane
        real(dp) :: jr(2, 2), inv(2, 2), det, shape(map%nloc), dshape(2, map%nloc)
        real(dp), allocatable :: qxi(:), qeta(:), qw(:)
        integer :: e, g, a, b, k, l, n, row, nnz, p, info
        integer, allocatable :: count(:), mark(:), cols(:, :)
        p = nint((sqrt(8.0_dp*map%nloc + 1) - 3)/2)
        call triangle_duffy_quadrature(max(7, 2*p + 1), qxi, qeta, qw, info)
        if (info /= 0) error stop 'winslow_seed: plane quadrature'
        plane%nq = size(qw)
        allocate (plane%dny(2, map%nloc, plane%nq, map%cells), plane%wq(plane%nq, map%cells))
        do g = 1, plane%nq
            call triangle_basis_jet([qxi(g), qeta(g)], shape, dshape)
            do e = 1, map%cells
                ! jr(i, r) = dy_i/dxi_r
                jr = matmul(map%reference(:, map%conn(:, e)), transpose(dshape))
                det = jr(1, 1)*jr(2, 2) - jr(1, 2)*jr(2, 1)
                inv = reshape([jr(2, 2), -jr(2, 1), -jr(1, 2), jr(1, 1)], [2, 2])/det
                ! dN/dy = jr^-t dN/dxi
                plane%dny(:, :, g, e) = matmul(transpose(inv), dshape)
                plane%wq(g, e) = qw(g)*abs(det)
            end do
        end do
        plane%area = sum(plane%wq)
        call setup_ring(map, plane)
        allocate (plane%dof(2, map%nodes))
        plane%dof = 0
        n = 0
        do a = 1, map%nodes
            if (map%boundary(a)) cycle
            do k = 1, 2
                n = n + 1
                plane%dof(k, a) = n
            end do
        end do
        plane%n = n
        ! CSR pattern of the node pairs sharing a cell, both components
        allocate (count(n), mark(n), cols(n, 2*map%nloc*12))
        count = 0; mark = 0
        do e = 1, map%cells
            do a = 1, map%nloc
                do k = 1, 2
                    row = plane%dof(k, map%conn(a, e))
                    if (row == 0) cycle
                    do b = 1, map%nloc
                        do l = 1, 2
                            call add(row, plane%dof(l, map%conn(b, e)))
                        end do
                    end do
                end do
            end do
        end do
        allocate (plane%rowptr(n + 1))
        plane%rowptr(1) = 1
        do row = 1, n
            plane%rowptr(row + 1) = plane%rowptr(row) + count(row)
        end do
        nnz = plane%rowptr(n + 1) - 1
        allocate (plane%colidx(nnz))
        do row = 1, n
            call sort(cols(row, :count(row)))
            plane%colidx(plane%rowptr(row):plane%rowptr(row + 1) - 1) = cols(row, :count(row))
        end do
    contains
        subroutine add(r, c)
            integer, intent(in) :: r, c
            if (c == 0) return
            if (any(cols(r, :count(r)) == c)) return
            if (count(r) == size(cols, 2)) error stop 'winslow_seed: node degree'
            count(r) = count(r) + 1
            cols(r, count(r)) = c
        end subroutine
    end subroutine

    pure subroutine sort(v)
        integer, intent(inout) :: v(:)
        integer :: i, j, t
        do i = 2, size(v)
            t = v(i); j = i - 1
            do while (j >= 1)
                if (v(j) <= t) exit
                v(j + 1) = v(j); j = j - 1
            end do
            v(j + 1) = t
        end do
    end subroutine

    subroutine cell_jet(map, plane, x, e, g, amat)
        ! A(k, i) = dX_k/dy_i at quadrature point g of cell e
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: e, g
        real(dp), intent(out) :: amat(2, 2)
        amat = matmul(x(:, map%conn(:, e)), transpose(plane%dny(:, :, g, e)))
    end subroutine

    subroutine plane_energy(map, plane, x, eps, energy, min_det)
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: x(:, :), eps
        real(dp), intent(out) :: energy, min_det
        real(dp) :: amat(2, 2), w, det
        integer :: e, g
        energy = 0; min_det = huge(1.0_dp)
        do e = 1, map%cells
            do g = 1, plane%nq
                call cell_jet(map, plane, x, e, g, amat)
                if (amat(1, 1)*amat(2, 2) - amat(1, 2)*amat(2, 1) < 0) then
                    call winslow_folded_density(amat(1, 1), amat(2, 1), amat(1, 2), &
                        amat(2, 2), eps, w, det)
                else
                    call winslow_density(amat(1, 1), amat(2, 1), amat(1, 2), amat(2, 2), &
                        eps, w, det)
                end if
                min_det = min(min_det, det)
                energy = energy + plane%wq(g, e)*w
            end do
        end do
        if (eps <= 0 .and. min_det <= 0) energy = huge(1.0_dp)
    end subroutine

    subroutine plane_assemble(map, plane, x, eps, gradient, values)
        ! gradient(n) and CSR values of the Hessian with respect to the unknowns
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: x(:, :), eps
        real(dp), intent(out) :: gradient(:), values(:)
        real(dp) :: amat(2, 2), gp(4), hp(10), h(4, 4), wgt
        integer :: e, g, a, b, k, l, i, j, p, q, row, col, idx
        integer, parameter :: packed(4, 4) = reshape([1, 2, 4, 7, 2, 3, 5, 8, 4, 5, 6, 9, &
            7, 8, 9, 10], [4, 4])
        gradient = 0; values = 0
        do e = 1, map%cells
            do g = 1, plane%nq
                call cell_jet(map, plane, x, e, g, amat)
                if (amat(1, 1)*amat(2, 2) - amat(1, 2)*amat(2, 1) < 0) then
                    call winslow_folded_gradient(amat(1, 1), amat(2, 1), amat(1, 2), &
                        amat(2, 2), eps, gp(1), gp(2), gp(3), gp(4))
                    call winslow_folded_hessian(amat(1, 1), amat(2, 1), amat(1, 2), &
                        amat(2, 2), eps, hp(1), hp(2), hp(3), hp(4), hp(5), hp(6), hp(7), &
                        hp(8), hp(9), hp(10))
                else
                    call winslow_gradient(amat(1, 1), amat(2, 1), amat(1, 2), amat(2, 2), &
                        eps, gp(1), gp(2), gp(3), gp(4))
                    call winslow_hessian(amat(1, 1), amat(2, 1), amat(1, 2), amat(2, 2), &
                        eps, hp(1), hp(2), hp(3), hp(4), hp(5), hp(6), hp(7), hp(8), hp(9), &
                        hp(10))
                end if
                do q = 1, 4
                    do p = 1, 4
                        h(p, q) = hp(packed(p, q))
                    end do
                end do
                wgt = plane%wq(g, e)
                do a = 1, map%nloc
                    do k = 1, 2
                        row = plane%dof(k, map%conn(a, e))
                        if (row == 0) cycle
                        ! p = k + 2 (i - 1) indexes A(k, i)
                        gradient(row) = gradient(row) + wgt*(gp(k)*plane%dny(1, a, g, e) &
                            + gp(k + 2)*plane%dny(2, a, g, e))
                        do b = 1, map%nloc
                            do l = 1, 2
                                col = plane%dof(l, map%conn(b, e))
                                if (col == 0) cycle
                                idx = locate(row, col)
                                do i = 1, 2
                                    do j = 1, 2
                                        values(idx) = values(idx) + wgt*h(k + 2*(i - 1), &
                                            l + 2*(j - 1))*plane%dny(i, a, g, e) &
                                            *plane%dny(j, b, g, e)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    contains
        integer function locate(r, c)
            integer, intent(in) :: r, c
            integer :: lo, hi, mid
            lo = plane%rowptr(r); hi = plane%rowptr(r + 1) - 1
            do while (lo < hi)
                mid = (lo + hi)/2
                if (plane%colidx(mid) < c) then
                    lo = mid + 1
                else
                    hi = mid
                end if
            end do
            locate = lo
        end function
    end subroutine

    subroutine plane_solve(map, plane, x, tol, status)
        ! Newton on E_eps with eps continuation to 0, Levenberg-Marquardt
        ! damping and an energy-decrease line search
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(in) :: plane
        real(dp), intent(inout) :: x(:, :)
        real(dp), intent(in) :: tol
        integer, intent(out) :: status
        integer, parameter :: max_newton = 200
        real(dp), parameter :: max_mu = 1e8_dp
        real(dp), allocatable :: gradient(:), values(:), step(:), trial(:, :), diag(:)
        real(dp) :: eps, energy, min_det, trial_energy, trial_det, mu, scale, g0, t
        integer :: it, k, a, row, info
        logical :: accepted, roundoff
        real(dp), allocatable :: trial_gradient(:), scratch(:)
        status = WINSLOW_FAILED
        allocate (gradient(plane%n), trial_gradient(plane%n), values(size(plane%colidx)), step(plane%n), &
            trial(size(x, 1), size(x, 2)), diag(plane%n), scratch(size(plane%colidx)))
        call plane_energy(map, plane, x, 0.0_dp, energy, min_det)
        ! the determinant scale of the plane: its area per reference area,
        ! positive for an admitted boundary
        scale = ring_area(plane, x)/plane%area
        eps = 0
        if (min_det <= 0) eps = scale
        do
            ! every continuation level starts from the same damping
            mu = 1e-6_dp
            call plane_energy(map, plane, x, eps, energy, min_det)
            if (.not. ieee_is_finite(energy)) return
            ! gradient scale: energy per length of the section
            g0 = energy/sqrt(scale*plane%area)
            do it = 1, max_newton
                call plane_assemble(map, plane, x, eps, gradient, values)
                if (.not. all(ieee_is_finite(gradient))) return
                if (.not. all(ieee_is_finite(values))) return
                if (norm2(gradient) <= tol*g0) exit
                do row = 1, plane%n
                    diag(row) = values(diagonal(row))
                end do
                accepted = .false.
                ! beyond max_mu the step is a vanishing gradient step
                do while (mu <= max_mu)
                    call damped_solve(plane, values, diag, mu, -gradient, step, info)
                    if (info == 0) then
                        trial = x
                        do a = 1, size(x, 2)
                            do k = 1, 2
                                row = plane%dof(k, a)
                                if (row > 0) trial(k, a) = x(k, a) + step(row)
                            end do
                        end do
                        call plane_energy(map, plane, trial, eps, trial_energy, trial_det)
                        ! below the roundoff of the energy sum a step is judged by
                        ! the gradient norm
                        roundoff = trial_energy < huge(1.0_dp) .and. &
                            abs(trial_energy - energy) <= 1e3_dp*epsilon(1.0_dp)*abs(energy)
                        if (roundoff) then
                            call plane_assemble(map, plane, trial, eps, trial_gradient, scratch)
                            roundoff = norm2(trial_gradient) < norm2(gradient)
                        end if
                        if (trial_energy < energy .or. roundoff) then
                            x = trial; energy = trial_energy; min_det = trial_det
                            mu = max(mu/10, 1e-12_dp)
                            accepted = .true.
                            exit
                        end if
                    end if
                    mu = mu*10
                end do
                if (.not. accepted) exit
            end do
            if (eps <= 0) exit
            ! continue towards eps = 0 once the determinants are positive
            if (min_det > 0) then
                t = eps
                eps = 0
                call plane_energy(map, plane, x, eps, energy, min_det)
                if (energy < huge(1.0_dp)) cycle
                eps = t
            end if
            eps = eps/10
            if (eps < 1e-14_dp*scale) return
        end do
        if (min_det <= 0) return
        if (.not. all(ieee_is_finite(x))) return
        call plane_assemble(map, plane, x, 0.0_dp, gradient, values)
        if (norm2(gradient) > tol*g0) return
        status = WINSLOW_OK
    contains
        integer function diagonal(r)
            integer, intent(in) :: r
            integer :: i
            do i = plane%rowptr(r), plane%rowptr(r + 1) - 1
                if (plane%colidx(i) == r) then
                    diagonal = i; return
                end if
            end do
            diagonal = 0
        end function
    end subroutine

    subroutine setup_ring(map, plane)
        ! the reference disk is convex, so its boundary nodes are ordered by
        ! their angle about the boundary centroid
        type(ideal_map_t), intent(in) :: map
        type(plane_t), intent(inout) :: plane
        real(dp), allocatable :: angle(:)
        real(dp) :: centre(2), t
        integer :: a, i, j, n, k
        n = count(map%boundary)
        allocate (plane%ring(n), angle(n))
        plane%ring = pack([(a, a=1, map%nodes)], map%boundary)
        do i = 1, 2
            centre(i) = sum(map%reference(i, plane%ring))/n
        end do
        do i = 1, n
            a = plane%ring(i)
            angle(i) = atan2(map%reference(2, a) - centre(2), map%reference(1, a) - centre(1))
        end do
        do i = 2, n
            t = angle(i); k = plane%ring(i); j = i - 1
            do while (j >= 1)
                if (angle(j) <= t) exit
                angle(j + 1) = angle(j); plane%ring(j + 1) = plane%ring(j); j = j - 1
            end do
            angle(j + 1) = t; plane%ring(j + 1) = k
        end do
    end subroutine

    real(dp) function ring_area(plane, x) result(area)
        ! signed area enclosed by the boundary polygon (shoelace)
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: x(:, :)
        integer :: i, a, b
        area = 0
        do i = 1, size(plane%ring)
            a = plane%ring(i); b = plane%ring(modulo(i, size(plane%ring)) + 1)
            area = area + x(1, a)*x(2, b) - x(1, b)*x(2, a)
        end do
        area = area/2
    end function

    logical function ring_simple(plane, x) result(simple)
        ! no two non-adjacent edges of the boundary polygon meet; a sampled
        ! test of the Jordan property, not a certificate for the curve
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: x(:, :)
        real(dp) :: p(2), q(2), r(2), u(2)
        integer :: i, k, n
        n = size(plane%ring)
        simple = .false.
        do i = 1, n
            p = x(:, plane%ring(i)); q = x(:, plane%ring(modulo(i, n) + 1))
            do k = i + 2, n
                if (i == 1 .and. k == n) cycle
                r = x(:, plane%ring(k)); u = x(:, plane%ring(modulo(k, n) + 1))
                if (segments_meet(p, q, r, u)) return
            end do
        end do
        simple = .true.
    end function

    pure logical function segments_meet(p, q, r, u) result(meet)
        ! closed segments pq and ru share a point
        real(dp), intent(in) :: p(2), q(2), r(2), u(2)
        real(dp) :: o(4)
        integer :: i
        o = [turn(p, q, r), turn(p, q, u), turn(r, u, p), turn(r, u, q)]
        if (o(1)*o(2) > 0 .or. o(3)*o(4) > 0) then
            meet = .false.
        else if (maxval(abs(o)) > 0) then
            meet = .true.
        else
            ! collinear: the projections overlap on both axes
            meet = .true.
            do i = 1, 2
                meet = meet .and. max(p(i), q(i)) >= min(r(i), u(i)) &
                    .and. max(r(i), u(i)) >= min(p(i), q(i))
            end do
        end if
    contains
        pure real(dp) function turn(a, b, c)
            real(dp), intent(in) :: a(2), b(2), c(2)
            turn = (b(1) - a(1))*(c(2) - a(2)) - (b(2) - a(2))*(c(1) - a(1))
        end function
    end function

    subroutine damped_solve(plane, values, diag, mu, rhs, step, status)
        ! (H + mu max|diag| I) step = rhs by sparse LU
        type(plane_t), intent(in) :: plane
        real(dp), intent(in) :: values(:), diag(:), mu, rhs(:)
        real(dp), intent(out) :: step(:)
        integer, intent(out) :: status
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sstatus
        real(dp), allocatable :: shifted(:)
        real(dp) :: shift
        integer :: r, i
        status = 1
        shifted = values
        shift = mu*maxval(abs(diag))
        do r = 1, plane%n
            do i = plane%rowptr(r), plane%rowptr(r + 1) - 1
                if (plane%colidx(i) == r) shifted(i) = shifted(i) + shift
            end do
        end do
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        call sparse_factor(solver, plane%n, plane%n, size(shifted), plane%rowptr, &
            plane%colidx, shifted, sstatus)
        if (.not. status_ok(sstatus)) return
        call sparse_solve(solver, rhs, step, sstatus)
        call sparse_destroy(solver)
        if (.not. status_ok(sstatus)) return
        if (.not. all(ieee_is_finite(step))) return
        status = 0
    end subroutine
end module milca_winslow_seed
