module milca_map_3d_loop
    ! Line integral of B along the image of a reference circle in the plane
    ! phi = 0: the toroidal current enclosed by the flux surface through that
    ! circle (every poloidal loop on the surface encloses the same current).
    ! The circle is cut at its crossings with the straight cell edges, found
    ! exactly from a quadratic, and every arc lies in one cell, where the
    ! integrand is smooth, so Gauss-Legendre on each arc converges
    ! spectrally. The arcs depend only on the reference mesh; they are built
    ! once per mesh and surface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_ideal_map, only: ideal_map_t, ideal_frozen_t
    use milca_curved_triangle, only: triangle_locate
    use milca_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    implicit none
    private
    public :: map_3d_loop_t, map_3d_loop_setup, map_3d_loop_integral

    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Gauss points per arc
    integer, parameter :: arc_points = 8

    type :: map_3d_loop_t
        ! quadrature points: owning cell, local coordinates, weight in theta
        ! and the reference tangent dy/dtheta
        integer, allocatable :: cell(:)
        real(dp), allocatable :: xi(:, :), weight(:), tangent(:, :)
    end type
contains

    subroutine map_3d_loop_setup(map, s, loop, status, elliptic)
        ! the circle of radius sqrt(s) times the reference radius about the
        ! reference centre, 0 < s < 1; with elliptic, the ellipse of sqrt(s)
        ! times the reference half-widths (an elliptic reference disc); status 1
        ! when it crosses a curved interior edge or a point is found in no cell
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: s
        type(map_3d_loop_t), intent(out) :: loop
        integer, intent(out) :: status
        logical, intent(in), optional :: elliptic
        real(dp), allocatable :: angles(:), sorted(:)
        real(dp) :: centre(2), radius(2), x(arc_points), w(arc_points), t, h
        real(dp) :: y(2), xi(2)
        integer :: e, i, k, n, owner, previous
        status = 1
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        radius = sqrt(s)*(maxval(map%reference(1, :)) - minval(map%reference(1, :)))/2
        if (present(elliptic)) then
            if (elliptic) radius(2) = sqrt(s)*(maxval(map%reference(2, :)) &
                - minval(map%reference(2, :)))/2
        end if
        allocate (angles(0))
        do e = 1, map%cells
            do i = 1, 3
                call edge_crossings(map, e, i, centre, radius, angles, status)
                if (status /= 0) return
            end do
        end do
        sorted = unique_angles(angles)
        if (size(sorted) == 0) sorted = [0.0_dp]
        sorted = [sorted, sorted(1) + 2*pi]
        call gauss_legendre(arc_points, x, w)
        n = (size(sorted) - 1)*arc_points
        allocate (loop%cell(n), loop%xi(2, n), loop%weight(n), loop%tangent(2, n))
        k = 0; previous = 1
        status = 1
        do i = 1, size(sorted) - 1
            h = sorted(i + 1) - sorted(i)
            ! the owner of the arc: the cell holding its midpoint
            t = (sorted(i) + sorted(i + 1))/2
            owner = find_cell(map, centre + radius*[cos(t), sin(t)], previous)
            if (owner == 0) return
            previous = owner
            do e = 1, arc_points
                t = sorted(i) + h*(x(e) + 1)/2
                y = centre + radius*[cos(t), sin(t)]
                call triangle_locate(map%reference(:, map%conn(:, owner)), y, xi, status)
                if (status /= 0) return
                k = k + 1
                loop%cell(k) = owner; loop%xi(:, k) = xi
                loop%weight(k) = h/2*w(e)
                loop%tangent(:, k) = radius*[-sin(t), cos(t)]
            end do
        end do
        status = 0
    end subroutine

    subroutine edge_crossings(map, cell, edge, centre, radius, list, status)
        ! angles at which the ellipse of half-widths radius crosses the edge
        ! from vertex edge to the next vertex of the cell, found on the unit
        ! circle of the scaled coordinates; status 1 when that edge is curved
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: centre(2), radius(2)
        real(dp), allocatable, intent(inout) :: list(:)
        integer, intent(out) :: status
        real(dp) :: a(2), b(2), d(2), f(2), qa, qb, qc, disc, root(2), p(2), node(2)
        integer :: j, m, degree
        status = 0
        a = map%reference(:, map%conn(edge, cell))
        b = map%reference(:, map%conn(mod(edge, 3) + 1, cell))
        d = (b - a)/radius; f = (a - centre)/radius
        qa = dot_product(d, d); qb = 2*dot_product(f, d)
        qc = dot_product(f, f) - 1
        disc = qb**2 - 4*qa*qc
        if (disc < 0) return
        root = [(-qb - sqrt(disc))/(2*qa), (-qb + sqrt(disc))/(2*qa)]
        if (all(root < 0 .or. root > 1)) return
        ! A wall edge lies on the reference boundary, the circle of the
        ! reference radius about the centre, which a loop with s < 1 never meets.
        if (map%boundary(map%conn(edge, cell)) .and. &
            map%boundary(map%conn(mod(edge, 3) + 1, cell))) return
        ! Any other crossed edge must be straight: in MILCA's local order (three
        ! vertices, then the degree - 1 nodes of each edge k -> k + 1) its
        ! interior nodes lie on the chord.
        degree = nint((sqrt(8.0_dp*map%nloc + 1) - 3)/2)
        do m = 4 + (edge - 1)*(degree - 1), 3 + edge*(degree - 1)
            node = (map%reference(:, map%conn(m, cell)) - a)/radius
            if (abs(d(1)*node(2) - d(2)*node(1)) > 1e-12_dp*qa) then
                status = 1
                return
            end if
        end do
        do j = 1, 2
            if (root(j) < 0 .or. root(j) > 1) cycle
            p = f + root(j)*d
            list = [list, modulo(atan2(p(2), p(1)), 2*pi)]
        end do
    end subroutine

    integer function find_cell(map, p, guess) result(owner)
        ! the cell holding the reference point p, searched from guess on
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: p(2)
        integer, intent(in) :: guess
        real(dp) :: local(2)
        integer :: c, j, info
        owner = 0
        do j = 0, map%cells - 1
            c = modulo(guess - 1 + j, map%cells) + 1
            call triangle_locate(map%reference(:, map%conn(:, c)), p, local, info)
            if (info /= 0) cycle
            if (minval([local, 1 - sum(local)]) < -1e-12_dp) cycle
            owner = c
            return
        end do
    end function

    pure function unique_angles(v) result(u)
        ! sorted angles with coincident crossings (an edge shared by two cells,
        ! a crossing at a vertex) merged
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: u(:)
        real(dp) :: a(size(v)), t
        integer :: i, j, n
        a = v
        do i = 2, size(a)
            t = a(i); j = i - 1
            do while (j >= 1)
                if (a(j) <= t) exit
                a(j + 1) = a(j); j = j - 1
            end do
            a(j + 1) = t
        end do
        n = 0
        do i = 1, size(a)
            if (n > 0) then
                if (a(i) - a(n) <= 1e-12_dp) cycle
            end if
            n = n + 1; a(n) = a(i)
        end do
        if (n > 1) then
            if (a(1) + 2*pi - a(n) <= 1e-12_dp) n = n - 1
        end if
        u = a(:n)
    end function

    real(dp) function map_3d_loop_integral(map, frozen, loop, status) result(integral)
        ! oint B . dl along the loop's image at phi = 0
        type(ideal_map_t), intent(in) :: map
        class(ideal_frozen_t), intent(in) :: frozen
        type(map_3d_loop_t), intent(in) :: loop
        integer, intent(out) :: status
        real(dp) :: q(3), d(3, 3), ref(2), g(2, 2), point(3), values(EQ_VALUES)
        real(dp) :: dxi(2), along(3), inverse(2, 2), det
        integer :: k
        integral = 0
        do k = 1, size(loop%cell)
            call map_3d_coordinates(map, loop%cell(k), loop%xi(:, k), 0.0_dp, q, d, ref, g)
            call map_3d_sample(map, frozen, loop%cell(k), loop%xi(:, k), 0.0_dp, point, &
                values, status)
            if (status /= 0) return
            ! dxi/dtheta = (dy/dxi)^-1 dy/dtheta
            det = g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1)
            inverse = reshape([g(2, 2), -g(2, 1), -g(1, 2), g(1, 1)], [2, 2])/det
            dxi = matmul(inverse, loop%tangent(:, k))
            ! (R, Z, phi) increments; the toroidal length is R dphi
            along = matmul(d(:, 1:2), dxi); along(3) = along(3)*q(1)
            integral = integral + loop%weight(k)*dot_product(values([EQ_BR, EQ_BZ, EQ_BPHI]), along)
        end do
        status = 0
    end function
end module milca_map_3d_loop
