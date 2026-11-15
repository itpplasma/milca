!> QUMOND, the quasi-linear formulation of MOND (Milgrom 2010), for an
!> axisymmetric mass distribution. Two linear Poisson problems with an
!> algebraic step between them:
!>   lap Phi_N = 4 pi G rho,          g_N = -grad Phi_N,
!>   gt = nu(|g_N|/a0) g_N,           lap Phi = -div gt,     g = -grad Phi.
!> The interpolating function is a declared model choice from the family
!>   nu_n(y) = (1/2 + sqrt(1/4 + y^-n))^(1/n),
!> n = 1 "simple", n = 2 "standard" (Famaey and McGaugh 2012).
!> Discretization: P1 finite elements on the meridional half plane (R, z)
!> with the volume weight 2 pi R, on a polar mesh with geometric radial and
!> uniform polar spacing. The equator theta = pi/2 is a mesh line, so a
!> midplane surface density enters as an edge load. Step 1 takes the monopole
!> potential -G M/r_max on the outer sphere. Step 3 is the natural problem
!> n.grad Phi = -n.gt with Phi = 0 at the origin: g is the R-weighted L2
!> projection of gt onto discrete gradients, so the curl of gt is removed and
!> ||g||^2 + ||g - gt||^2 = ||gt||^2 holds to round-off.
module milca_qumond
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    use fortsparse, only: csc_t, csc_from_triplet, sparse_solver_t, sparse_factor, &
                          sparse_solve, sparse_destroy, fortsparse_status_t, status_ok, &
                          FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_SPD
    implicit none
    private
    public :: qumond_t, qumond_field_t, mass_density
    public :: qumond_setup, qumond_solve, qumond_free, qumond_nu, qumond_node, qumond_cell

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    !> Collapsed Gauss-Legendre points per direction on a triangle.
    integer, parameter :: nquad = 4

    abstract interface
        !> Volume density rho(R, z), or a midplane surface density Sigma(R, 0).
        function mass_density(r_cyl, z) result(rho)
            import :: dp
            real(dp), intent(in) :: r_cyl, z
            real(dp) :: rho
        end function mass_density
    end interface

    !> Polar P1 mesh and the two factored Laplacians. Node 1 is the origin,
    !> node(i, j) sits at radius r(i), polar angle j pi/nth.
    type :: qumond_t
        integer :: nr = 0, nth = 0, nnode = 0, ntri = 0
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: x(:, :)
        integer, allocatable :: tri(:, :)
        !> Basis gradients (2, 3, ntri), volume weights 2 pi R_c |T|, centroids.
        real(dp), allocatable :: grad(:, :, :), weight(:), centroid(:, :)
        !> Row sums of the stiffness over the outer-sphere nodes.
        real(dp), allocatable :: outer_coupling(:)
        type(sparse_solver_t) :: dirichlet, neumann
        logical :: ok = .false.
    end type

    !> Nodal potentials and element accelerations: Newtonian g_n, algebraic
    !> gt = nu g_n, and QUMOND g = -grad Phi. Mass is the discrete total.
    type :: qumond_field_t
        real(dp) :: mass = 0
        real(dp), allocatable :: phi_n(:), phi(:)
        real(dp), allocatable :: g_n(:, :), g_alg(:, :), g(:, :)
        logical :: ok = .false.
    end type

contains

    !> nu_n(y); the y < 1 branch avoids overflow of y^-n.
    elemental function qumond_nu(n, y) result(nu)
        integer, intent(in) :: n
        real(dp), intent(in) :: y
        real(dp) :: nu, yn
        if (y >= 1) then
            nu = (0.5_dp + sqrt(0.25_dp + y**(-n)))**(1.0_dp/n)
        else
            yn = y**n
            nu = (0.5_dp*yn + sqrt(0.25_dp*yn**2 + yn))**(1.0_dp/n)/y
        end if
    end function qumond_nu

    pure integer function qumond_node(q, i, j) result(k)
        type(qumond_t), intent(in) :: q
        integer, intent(in) :: i, j
        if (i == 0) then
            k = 1
        else
            k = 1 + (i - 1)*(q%nth + 1) + j + 1
        end if
    end function qumond_node

    !> Triangle in the cell between rings i and i + 1 and polar interval j;
    !> outer selects the one whose edge lies on ring i + 1. Ring 0 is the
    !> origin, whose fan has one triangle per interval.
    pure integer function qumond_cell(q, i, j, outer) result(t)
        type(qumond_t), intent(in) :: q
        integer, intent(in) :: i, j
        logical, intent(in) :: outer
        if (i == 0) then
            t = j
        else
            t = q%nth + 2*((i - 1)*q%nth + j - 1) + merge(1, 2, outer)
        end if
    end function qumond_cell

    subroutine qumond_setup(q, r_min, r_max, nr, nth)
        type(qumond_t), intent(inout) :: q
        real(dp), intent(in) :: r_min, r_max
        integer, intent(in) :: nr, nth
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        real(dp) :: theta, p(2, 3), area2, kab
        integer :: i, j, t, a, b, n, na, nb
        logical :: outer_a, outer_b

        call qumond_free(q)
        if (nr < 2 .or. nth < 4 .or. mod(nth, 2) /= 0) return
        if (.not. (r_min > 0 .and. r_max > r_min)) return
        q%nr = nr
        q%nth = nth
        q%nnode = 1 + nr*(nth + 1)
        q%ntri = nth + 2*nth*(nr - 1)
        allocate (q%r(nr), q%x(2, q%nnode), q%tri(3, q%ntri))
        allocate (q%grad(2, 3, q%ntri), q%weight(q%ntri), q%centroid(2, q%ntri))
        do i = 1, nr
            q%r(i) = r_min*(r_max/r_min)**(real(i - 1, dp)/(nr - 1))
        end do
        q%x(:, 1) = 0
        do i = 1, nr
            do j = 0, nth
                theta = pi*j/nth
                q%x(:, qumond_node(q, i, j)) = q%r(i)*[sin(theta), cos(theta)]
            end do
            ! axis and equator exactly on their lines
            q%x(1, qumond_node(q, i, 0)) = 0
            q%x(1, qumond_node(q, i, nth)) = 0
            q%x(2, qumond_node(q, i, nth/2)) = 0
        end do
        do j = 1, nth
            q%tri(:, j) = [1, qumond_node(q, 1, j - 1), qumond_node(q, 1, j)]
        end do
        do i = 1, nr - 1
            do j = 1, nth
                q%tri(:, qumond_cell(q, i, j, .true.)) = [qumond_node(q, i, j - 1), &
                    qumond_node(q, i + 1, j - 1), qumond_node(q, i + 1, j)]
                q%tri(:, qumond_cell(q, i, j, .false.)) = [qumond_node(q, i, j - 1), &
                    qumond_node(q, i + 1, j), qumond_node(q, i, j)]
            end do
        end do
        do t = 1, q%ntri
            p = q%x(:, q%tri(:, t))
            area2 = (p(1, 2) - p(1, 1))*(p(2, 3) - p(2, 1)) - (p(1, 3) - p(1, 1))*(p(2, 2) - p(2, 1))
            q%grad(:, 1, t) = [p(2, 2) - p(2, 3), p(1, 3) - p(1, 2)]/area2
            q%grad(:, 2, t) = [p(2, 3) - p(2, 1), p(1, 1) - p(1, 3)]/area2
            q%grad(:, 3, t) = [p(2, 1) - p(2, 2), p(1, 2) - p(1, 1)]/area2
            q%centroid(:, t) = sum(p, dim=2)/3
            ! R is linear on the triangle, so the centroid rule is exact
            q%weight(t) = 2*pi*q%centroid(1, t)*abs(area2)/2
        end do

        ! Step 1 eliminates the outer-sphere nodes, step 3 pins the origin.
        allocate (q%outer_coupling(q%nnode))
        q%outer_coupling = 0
        allocate (rows(9*q%ntri + q%nnode), cols(9*q%ntri + q%nnode), vals(9*q%ntri + q%nnode))
        n = 0
        do t = 1, q%ntri
            do a = 1, 3
                na = q%tri(a, t)
                outer_a = na > q%nnode - (nth + 1)
                do b = 1, 3
                    nb = q%tri(b, t)
                    outer_b = nb > q%nnode - (nth + 1)
                    kab = q%weight(t)*dot_product(q%grad(:, a, t), q%grad(:, b, t))
                    if (outer_b .and. .not. outer_a) q%outer_coupling(na) = q%outer_coupling(na) + kab
                    if (outer_a .or. outer_b) cycle
                    n = n + 1
                    rows(n) = na
                    cols(n) = nb
                    vals(n) = kab
                end do
            end do
        end do
        do na = q%nnode - nth, q%nnode
            n = n + 1
            rows(n) = na
            cols(n) = na
            vals(n) = 1
        end do
        call factor(q%dirichlet, rows(:n), cols(:n), vals(:n), q%nnode, q%ok)
        if (.not. q%ok) return
        n = 0
        do t = 1, q%ntri
            do a = 1, 3
                do b = 1, 3
                    if (q%tri(a, t) == 1 .or. q%tri(b, t) == 1) cycle
                    n = n + 1
                    rows(n) = q%tri(a, t)
                    cols(n) = q%tri(b, t)
                    vals(n) = q%weight(t)*dot_product(q%grad(:, a, t), q%grad(:, b, t))
                end do
            end do
        end do
        n = n + 1
        rows(n) = 1
        cols(n) = 1
        vals(n) = 1
        call factor(q%neumann, rows(:n), cols(:n), vals(:n), q%nnode, q%ok)
    end subroutine qumond_setup

    subroutine factor(solver, rows, cols, vals, n, ok)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: rows(:), cols(:), n
        real(dp), intent(in) :: vals(:)
        logical, intent(out) :: ok
        type(csc_t) :: matrix
        type(fortsparse_status_t) :: status
        call csc_from_triplet(n, n, rows, cols, vals, matrix, status)
        ok = status_ok(status)
        if (.not. ok) return
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%refine = .false.
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(solver, matrix, status)
        ok = status_ok(status)
    end subroutine factor

    !> QUMOND field of the volume density rho(R, z) plus an optional midplane
    !> surface density Sigma(R), for the coupling G, acceleration scale a0 and
    !> interpolating function nu_n.
    subroutine qumond_solve(q, grav, a0, n_nu, density, field, surface)
        type(qumond_t), intent(inout) :: q
        real(dp), intent(in) :: grav, a0
        integer, intent(in) :: n_nu
        procedure(mass_density) :: density
        type(qumond_field_t), intent(out) :: field
        procedure(mass_density), optional :: surface
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: load(:), rhs(:)
        real(dp) :: boundary, gn, y
        integer :: t, a, outer0

        if (.not. q%ok .or. .not. (grav > 0 .and. a0 > 0) .or. n_nu < 1) return
        allocate (load(q%nnode), rhs(q%nnode))
        call mass_load(q, density, load)
        if (present(surface)) call surface_load(q, surface, load)
        field%mass = sum(load)

        ! step 1: Newtonian potential, monopole value on the outer sphere
        outer0 = q%nnode - q%nth
        boundary = -grav*field%mass/q%r(q%nr)
        rhs = -4*pi*grav*load - boundary*q%outer_coupling
        rhs(outer0:) = boundary
        allocate (field%phi_n(q%nnode), field%phi(q%nnode))
        call sparse_solve(q%dirichlet, rhs, field%phi_n, status)
        if (.not. status_ok(status)) return
        allocate (field%g_n(2, q%ntri), field%g_alg(2, q%ntri), field%g(2, q%ntri))
        call accelerations(q, field%phi_n, field%g_n)

        ! step 2: algebraic relation
        do t = 1, q%ntri
            gn = norm2(field%g_n(:, t))
            if (gn > 0) then
                y = gn/a0
                field%g_alg(:, t) = qumond_nu(n_nu, y)*field%g_n(:, t)
            else
                field%g_alg(:, t) = 0
            end if
        end do

        ! step 3: Phi from the divergence of gt, natural boundary condition
        rhs = 0
        do t = 1, q%ntri
            do a = 1, 3
                rhs(q%tri(a, t)) = rhs(q%tri(a, t)) - &
                    q%weight(t)*dot_product(field%g_alg(:, t), q%grad(:, a, t))
            end do
        end do
        rhs(1) = 0
        call sparse_solve(q%neumann, rhs, field%phi, status)
        if (.not. status_ok(status)) return
        call accelerations(q, field%phi, field%g)
        field%ok = all(ieee_is_finite(field%phi)) .and. all(ieee_is_finite(field%g))
    end subroutine qumond_solve

    subroutine accelerations(q, phi, g)
        type(qumond_t), intent(in) :: q
        real(dp), intent(in) :: phi(:)
        real(dp), intent(out) :: g(:, :)
        integer :: t
        do t = 1, q%ntri
            g(:, t) = -matmul(q%grad(:, :, t), phi(q%tri(:, t)))
        end do
    end subroutine accelerations

    !> load(a) = int rho phi_a dV by collapsed Gauss-Legendre on each triangle.
    subroutine mass_load(q, density, load)
        type(qumond_t), intent(in) :: q
        procedure(mass_density) :: density
        real(dp), intent(out) :: load(:)
        real(dp) :: u(nquad), w(nquad), xi, eta, wq, shape(3), p(2, 3), point(2), area2
        integer :: t, k, l
        call gauss_legendre(nquad, u, w)
        u = (u + 1)/2
        w = w/2
        load = 0
        do t = 1, q%ntri
            p = q%x(:, q%tri(:, t))
            area2 = abs((p(1, 2) - p(1, 1))*(p(2, 3) - p(2, 1)) - (p(1, 3) - p(1, 1))*(p(2, 2) - p(2, 1)))
            do k = 1, nquad
                do l = 1, nquad
                    xi = u(k)
                    eta = u(l)*(1 - u(k))
                    wq = w(k)*w(l)*(1 - u(k))*area2
                    shape = [1 - xi - eta, xi, eta]
                    point = matmul(p, shape)
                    load(q%tri(:, t)) = load(q%tri(:, t)) + &
                        wq*2*pi*point(1)*density(point(1), point(2))*shape
                end do
            end do
        end do
    end subroutine mass_load

    !> Edge load of a midplane surface density along theta = pi/2.
    subroutine surface_load(q, surface, load)
        type(qumond_t), intent(in) :: q
        procedure(mass_density) :: surface
        real(dp), intent(inout) :: load(:)
        real(dp) :: u(nquad), w(nquad), r0, r1, rq, f
        integer :: i, k, n0, n1
        call gauss_legendre(nquad, u, w)
        u = (u + 1)/2
        w = w/2
        do i = 0, q%nr - 1
            n0 = qumond_node(q, i, q%nth/2)
            n1 = qumond_node(q, i + 1, q%nth/2)
            r0 = q%x(1, n0)
            r1 = q%x(1, n1)
            do k = 1, nquad
                rq = r0 + u(k)*(r1 - r0)
                f = w(k)*(r1 - r0)*2*pi*rq*surface(rq, 0.0_dp)
                load(n0) = load(n0) + f*(1 - u(k))
                load(n1) = load(n1) + f*u(k)
            end do
        end do
    end subroutine surface_load

    subroutine qumond_free(q)
        type(qumond_t), intent(inout) :: q
        if (q%ok .or. allocated(q%r)) then
            call sparse_destroy(q%dirichlet)
            call sparse_destroy(q%neumann)
        end if
        if (allocated(q%r)) deallocate (q%r, q%x, q%tri, q%grad, q%weight, q%centroid)
        if (allocated(q%outer_coupling)) deallocate (q%outer_coupling)
        q%ok = .false.
    end subroutine qumond_free

end module milca_qumond
