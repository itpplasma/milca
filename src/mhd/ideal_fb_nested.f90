module milca_ideal_fb_nested
    !! Set-up of the nested-map free boundary from a VMEC-style boundary guess
    !! (docs/free-boundary-nested-map-3d.md, section 10): the reference disc
    !! of the frozen field, the axisymmetric coupling curve Gamma_c, the mesh
    !! of the region inside it, the seed of the plasma map and its transfer
    !! between resolution levels. The physics stays in
    !! milca_ideal_free_boundary; this module only places and seeds it.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_grad_shafranov, only: gs_state_t, initialize_gs
    use milca_ideal_free_boundary, only: ideal_fb_t, ideal_fb_field_basis_t, &
        ideal_fb_solve_shell_map, ideal_fb_solve_vacuum, ideal_fb_solve_multipliers, &
        IDEAL_FB_OK, IDEAL_FB_BAD_INPUT, IDEAL_FB_FOLDED
    use milca_ideal_seed_3d, only: ideal_seed_boundary_coefficients, ideal_seed_disk_mesh
    use milca_sort_keys, only: sort_key_order
    use milca_winslow_seed, only: winslow_seed, WINSLOW_OK
    use milca_curved_triangle, only: triangle_locate
    use milca_triangle_basis, only: triangle_basis_jet
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    implicit none
    private
    public :: fb_nested_geometry, fb_nested_mesh, fb_nested_seed, fb_nested_transfer
    public :: fb_nested_settle, fb_nested_layer_mesh, fb_nested_offset_distance
    public :: filament_basis_t
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Gamma_c is a star-shaped curve r(alpha) about the centre, sampled on
    ! rays directions and smoothed over 2 window + 1 of them.
    integer, parameter :: rays = 512, window = 8

    !> The field of an axisymmetric filament at (R, Z) = source with unit
    !> poloidal circulation: the multivalued part c theta/2 pi of the shell
    !> potential continued outside Gamma_c, so that the decaying exterior is
    !> single valued (ideal_fb_drive_traces with circulation 1).
    type, extends(ideal_fb_field_basis_t) :: filament_basis_t
        real(dp) :: source(2) = 0
    contains
        procedure :: nfield => filament_nfield
        procedure :: sample => filament_sample
    end type

contains

    integer function filament_nfield(self)
        class(filament_basis_t), intent(in) :: self
        filament_nfield = 1
    end function

    subroutine filament_sample(self, rz, nfp, ntor, bcos, bsin, status)
        class(filament_basis_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
        integer, intent(out) :: status
        real(dp) :: g, gg(4)
        integer :: q
        status = IDEAL_FB_BAD_INPUT
        if (nfp < 1 .or. ntor < 0 .or. size(bcos, 3) /= 1) return
        bcos = 0; bsin = 0
        do q = 1, size(rz, 2)
            call gs_green_gradient(rz(1, q), rz(2, q), self%source(1), self%source(2), &
                g, gg)
            bcos(1, 0, 1, q) = gg(2)/rz(1, q)
            bcos(3, 0, 1, q) = -gg(1)/rz(1, q)
        end do
        status = IDEAL_FB_OK
    end subroutine

    subroutine fb_nested_geometry(boundary, nfp, margin, box, obstacles, centre, axes, &
                                  curve, margin_used, status, message)
        !! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp
        !! phi), counter-clockwise in (R, Z). centre and axes: the elliptic
        !! reference disc inscribed in the mean bounding box of the sections,
        !! so that the shell map starts close to the identity.
        !! curve(k): Gamma_c as r(alpha) about the centre at alpha = 2 pi k/rays,
        !! linear in between. It keeps the clearance margin*a from every section of
        !! the guess, a the minor radius of the mean section area, and stays in
        !! the box (rmin, rmax, zmin, zmax) of the coil field and off the
        !! current-carrying points obstacles(2, :) (coils inside the box);
        !! where they are closer, the clearance shrinks to a share of the room.
        !! margin_used: the smallest clearance / a.
        real(dp), intent(in) :: boundary(:, :, 0:)
        integer, intent(in) :: nfp
        real(dp), intent(in) :: margin, box(4), obstacles(:, :)
        real(dp), intent(out) :: centre(2), axes(2), margin_used
        real(dp), allocatable, intent(out) :: curve(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        integer, parameter :: ntheta = 1024
        real(dp), allocatable :: rz(:, :, :)
        real(dp) :: envelope(rays), fit(rays), room(rays), area, centroid(2), sums(3)
        real(dp) :: cross, minor, low(2), high(2), core(rays)
        real(dp) :: limit(rays)
        integer :: nphi, ntor, j, i, k
        status = IDEAL_FB_BAD_INPUT
        message = ''
        ntor = (size(boundary, 2) - 1)/2
        nphi = 4*ntor + 1
        allocate (rz(2, ntheta, nphi))
        sums = 0; low = 0; high = 0
        do j = 1, nphi
            call section(2*pi*(j - 1)/(nfp*nphi), rz(:, :, j))
            area = 0; centroid = 0
            do i = 1, ntheta
                k = mod(i, ntheta) + 1
                cross = rz(1, i, j)*rz(2, k, j) - rz(1, k, j)*rz(2, i, j)
                area = area + cross/2
                centroid = centroid + cross*(rz(:, i, j) + rz(:, k, j))/6
            end do
            if (.not. area > 0) then
                message = 'the boundary is clockwise or degenerate in (R, Z)'
                return
            end if
            sums = sums + [area, centroid]
            low = low + minval(rz(:, :, j), dim=2)/nphi
            high = high + maxval(rz(:, :, j), dim=2)/nphi
        end do
        centre = (low + high)/2
        axes = (high - low)/2
        ! the clearance scales with the minor radius of the mean section area
        minor = sqrt(sums(1)/nphi/pi)
        ! Gamma_c: the guess's sections and the reference ellipse widened by
        ! the clearance d (their Minkowski sum with a disc of radius d, which
        ! keeps d to every section point), capped at a share of the room to the
        ! box (1 % of it kept free) and to the coils inside it, and smoothed by
        ! two moving averages over 2 window + 1 rays. The room is first taken
        ! as its minimum over 4 window + 1 rays, so every average stays below
        ! the room of its own ray. d = margin*a, and the share of the room is
        ! 0.7.
        do k = 1, rays
            room(k) = free_distance([cos(2*pi*(k - 1)/rays), sin(2*pi*(k - 1)/rays)])
        end do
        do k = 1, rays
            limit(k) = minval(room([(modulo(k - 1 + i, rays) + 1, i=-2*window, &
                2*window)]))
        end do
        call outline(core)
        if (any(core >= room)) then
            message = 'a coil or the edge of the coil field grid touches the' &
                //' boundary guess'
            return
        end if
        call widened(margin*minor, envelope)
        envelope = min(envelope, core + 0.7_dp*(limit - core))
        fit = average(average(envelope))
        if (.not. all(fit > core + 0.02_dp*minor)) then
            message = 'the coils or the edge of the coil field grid leave no room for' &
                //' the coupling surface around the boundary guess'
            return
        end if
        margin_used = minval(fit - core)/minor
        curve = fit
        status = IDEAL_FB_OK
    contains
        subroutine section(phi, points)
            real(dp), intent(in) :: phi
            real(dp), intent(out) :: points(:, :)
            real(dp) :: theta, arg
            integer :: m, n, i
            points = 0
            do i = 1, ntheta
                theta = 2*pi*(i - 1)/ntheta
                do m = 0, ubound(boundary, 3)
                    do n = -ntor, ntor
                        arg = m*theta - n*nfp*phi
                        points(1, i) = points(1, i) + boundary(1, n + ntor + 1, m) &
                            *cos(arg) + boundary(2, n + ntor + 1, m)*sin(arg)
                        points(2, i) = points(2, i) + boundary(3, n + ntor + 1, m) &
                            *cos(arg) + boundary(4, n + ntor + 1, m)*sin(arg)
                    end do
                end do
            end do
        end subroutine

        subroutine widened(clearance, r)
            ! r(alpha) of the union of discs of radius clearance about every
            ! section point and about the reference ellipse
            real(dp), intent(in) :: clearance
            real(dp), intent(out) :: r(:)
            real(dp) :: e(2), q(2), along, across
            integer :: k, i, j
            r = 0
            do k = 1, rays
                e = [cos(2*pi*(k - 1)/rays), sin(2*pi*(k - 1)/rays)]
                do j = 1, nphi + 1
                    do i = 1, ntheta
                        if (j <= nphi) then
                            q = rz(:, i, j) - centre
                        else
                            q = axes*[cos(2*pi*(i - 1)/ntheta), &
                                sin(2*pi*(i - 1)/ntheta)]
                        end if
                        along = dot_product(q, e); across = q(2)*e(1) - q(1)*e(2)
                        if (abs(across) < clearance) &
                            r(k) = max(r(k), along + sqrt(clearance**2 - across**2))
                    end do
                end do
            end do
        end subroutine

        subroutine outline(r)
            ! r(alpha) of the farthest crossing of the ray with every section
            ! and with the reference ellipse
            real(dp), intent(out) :: r(:)
            real(dp) :: e(2), a(2), b(2), f(2), det, along, s
            integer :: k, i, j
            r = 0
            do k = 1, rays
                e = [cos(2*pi*(k - 1)/rays), sin(2*pi*(k - 1)/rays)]
                ! the ellipse (x/axes)^2 = 1 along the ray
                r(k) = 1/norm2(e/axes)
                do j = 1, nphi
                    do i = 1, ntheta
                        a = rz(:, i, j) - centre
                        b = rz(:, mod(i, ntheta) + 1, j) - centre
                        f = b - a
                        det = f(1)*e(2) - f(2)*e(1)
                        if (abs(det) <= tiny(1.0_dp)) cycle
                        ! a + s f = along e
                        s = (a(2)*e(1) - a(1)*e(2))/det
                        if (s < 0 .or. s > 1) cycle
                        along = dot_product(a + s*f, e)
                        r(k) = max(r(k), along)
                    end do
                end do
            end do
        end subroutine

        function average(r) result(mean)
            real(dp), intent(in) :: r(:)
            real(dp) :: mean(size(r))
            integer :: k, i
            do k = 1, rays
                mean(k) = sum(r([(modulo(k - 1 + i, rays) + 1, i=-window, window)])) &
                    /(2*window + 1)
            end do
        end function

        real(dp) function free_distance(d) result(t)
            ! along the ray centre + t d to the box or to a disc of radius pad
            ! about a coil point
            real(dp), intent(in) :: d(2)
            real(dp) :: pad(2), q(2), along, across
            integer :: k
            pad = 0.01_dp*[box(2) - box(1), box(4) - box(3)]
            t = huge(1.0_dp)
            do k = 1, size(obstacles, 2)
                q = obstacles(:, k) - centre
                along = dot_product(q, d); across = q(2)*d(1) - q(1)*d(2)
                if (along > 0 .and. abs(across) < maxval(pad)) &
                    t = min(t, along - sqrt(maxval(pad)**2 - across**2))
            end do
            if (d(1) > 0) t = min(t, (box(2) - pad(1) - centre(1))/d(1))
            if (d(1) < 0) t = min(t, (box(1) + pad(1) - centre(1))/d(1))
            if (d(2) > 0) t = min(t, (box(4) - pad(2) - centre(2))/d(2))
            if (d(2) < 0) t = min(t, (box(3) + pad(2) - centre(2))/d(2))
        end function
    end subroutine

    pure real(dp) function curve_radius(curve, alpha) result(r)
        ! periodic linear interpolation of the ray samples
        real(dp), intent(in) :: curve(:), alpha
        real(dp) :: t
        integer :: k
        t = modulo(alpha/(2*pi), 1.0_dp)*size(curve)
        k = min(int(t), size(curve) - 1)
        t = t - k
        r = (1 - t)*curve(k + 1) + t*curve(mod(k + 1, size(curve)) + 1)
    end function

    subroutine fb_nested_mesh(mesh, chain, centre, axes, curve, nodes, degree, status)
        !! Reference mesh of the region inside Gamma_c with the reference
        !! ellipse (nodes points) as the constrained chain Gamma_p. The cells
        !! are at most the chain spacing and twice the narrowest reference
        !! gap between Gamma_p and Gamma_c, so that a thin shell has cells
        !! across it; Gamma_c gets the same spacing.
        type(gs_state_t), intent(out) :: mesh
        real(dp), allocatable, intent(out) :: chain(:, :)
        real(dp), intent(in) :: centre(2), axes(2), curve(:)
        integer, intent(in) :: nodes, degree
        integer, intent(out) :: status
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h, length, alpha, gap
        integer :: j, nc, nchain
        h = 0
        do j = 1, rays
            alpha = 2*pi*(j - 1)/rays
            h = h + norm2(axes*([cos(alpha + 2*pi/rays), sin(alpha + 2*pi/rays)] &
                - [cos(alpha), sin(alpha)]))/nodes
        end do
        length = 0; gap = huge(1.0_dp)
        do j = 1, rays
            alpha = 2*pi*(j - 1)/rays
            length = length + norm2(point(alpha) - point(alpha + 2*pi/rays))
            gap = min(gap, curve_radius(curve, alpha) - 1/norm2([cos(alpha), &
                sin(alpha)]/axes))
        end do
        ! a chain segment longer than the gap would let its curved edge bulge
        ! to Gamma_c: the chain gets at least one node per 2 gaps
        nchain = max(nodes, 4*ceiling(nodes*h/(8*gap)))
        h = min(h*nodes/nchain, 2*gap)
        nc = max(nchain, ceiling(length/h))
        allocate (chain(2, nchain), contour(2, nc))
        do j = 1, nchain
            alpha = 2*pi*(j - 1)/nchain
            chain(:, j) = centre + axes*[cos(alpha), sin(alpha)]
        end do
        do j = 1, nc
            contour(:, j) = point(2*pi*(j - 1)/nc)
        end do
        call initialize_gs(mesh, contour, h*h/2, status, degree, curved=.true., &
            interior_chain=chain, geometry_only=.true.)
    contains
        function point(a) result(p)
            real(dp), intent(in) :: a
            real(dp) :: p(2)
            p = centre + curve_radius(curve, a)*[cos(a), sin(a)]
        end function
    end subroutine

    subroutine fb_nested_seed(fb, boundary, centre, axes, status)
        !! Plasma map from the boundary guess: harmonic extension on the
        !! reference disc, then the Winslow map of every plane, which is
        !! injective for any simple cross-section; shell map, vacuum potential
        !! and multipliers solved at that geometry.
        type(ideal_fb_t), intent(inout) :: fb
        real(dp), intent(in) :: boundary(:, :, 0:), centre(2), axes(2)
        integer, intent(out) :: status
        integer :: a
        status = IDEAL_FB_BAD_INPUT
        if (size(boundary, 2) /= fb%nmode) return
        fb%plasma%x = 0
        do a = 1, fb%plasma%nodes
            call ideal_seed_boundary_coefficients(fb%plasma%reference(:, a), boundary, &
                fb%plasma%x(:, :, a), centre=centre, axes=axes)
        end do
        call winslow_seed(fb%plasma, a)
        status = IDEAL_FB_FOLDED
        if (a /= WINSLOW_OK) return
        call fb_nested_settle(fb, status)
    end subroutine

    subroutine fb_nested_transfer(fb, old, status)
        !! Plasma map of a solved level interpolated to the new reference mesh
        !! (same reference disc), toroidal modes padded with zeros; shell map,
        !! potential and multipliers solved at the new geometry.
        type(ideal_fb_t), intent(inout) :: fb
        type(ideal_fb_t), intent(in) :: old
        integer, intent(out) :: status
        real(dp) :: xi(2), best(2), lowest, quality
        real(dp), allocatable :: value(:), gradient(:, :)
        integer :: a, e, found, local, m, k
        status = IDEAL_FB_BAD_INPUT
        allocate (value(size(old%plasma%conn, 1)), &
            gradient(2, size(old%plasma%conn, 1)))
        m = min(fb%nmode, old%nmode)
        fb%plasma%x = 0
        do a = 1, fb%plasma%nodes
            ! the old cell where the point is most inside (chain points of the
            ! new level lie on the old chain only up to the spline error)
            found = 0; lowest = -huge(1.0_dp)
            do e = 1, old%plasma%cells
                call triangle_locate(old%plasma%reference(:, old%plasma%conn(:, e)), &
                    fb%plasma%reference(:, a), xi, local)
                if (local /= 0) cycle
                quality = minval([xi, 1 - sum(xi)])
                if (quality > lowest) then
                    lowest = quality; found = e; best = xi
                end if
                if (quality >= 0) exit
            end do
            if (found == 0 .or. lowest < -5e-2_dp) return
            call triangle_basis_jet(best, value, gradient)
            do k = 1, m
                fb%plasma%x(:, k, a) = matmul(old%plasma%x(:, k, &
                    old%plasma%conn(:, found)), value)
            end do
        end do
        call fb_nested_settle(fb, status)
    end subroutine

    subroutine fb_nested_settle(fb, status)
        ! Constraint, vacuum and multiplier rows solved at the plasma geometry.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(out) :: status
        integer :: k
        do k = 1, fb%ns
            if (fb%role(k) == 1) fb%xs(:, :, k) = fb%plasma%x(:, :, fb%plasma_node(k))
        end do
        call ideal_fb_solve_shell_map(fb, status)
        if (status /= IDEAL_FB_OK) return
        fb%phi = 0
        call ideal_fb_solve_vacuum(fb, status)
        if (status /= IDEAL_FB_OK) return
        fb%lam = 0
        call ideal_fb_solve_multipliers(fb, status)
        if (status /= IDEAL_FB_OK) return
        if (.not. all(ieee_is_finite(fb%xs))) status = IDEAL_FB_FOLDED
    end subroutine

    subroutine fb_nested_offset_distance(boundary, centre, axes, box, obstacles, cap, &
                                         distance, status)
        !! Distance d of Gamma_c as an outward normal offset of the guess's
        !! mean section: below the reach of its normal lines (0.9 times the
        !! smallest radius of curvature where the section is concave; a convex
        !! section has unbounded reach), half the smallest gap to the coils
        !! inside the box and to the box edge, and cap.
        real(dp), intent(in) :: boundary(:, :, 0:), centre(2), axes(2), box(4)
        real(dp), intent(in) :: obstacles(:, :), cap
        real(dp), intent(out) :: distance
        integer, intent(out) :: status
        integer, parameter :: ns = 1024
        real(dp) :: x(2, ns), t(2), tt(2), kappa, gap, reach
        integer :: j, jm, jp, k
        status = IDEAL_FB_BAD_INPUT
        do j = 1, ns
            x(:, j) = section_point(boundary, centre, axes, 2*pi*(j - 1)/ns)
        end do
        reach = huge(1.0_dp); gap = huge(1.0_dp)
        do j = 1, ns
            jm = modulo(j - 2, ns) + 1; jp = mod(j, ns) + 1
            t = (x(:, jp) - x(:, jm))/2
            tt = x(:, jp) - 2*x(:, j) + x(:, jm)
            ! signed curvature, positive where the counter-clockwise section
            ! turns left (convex)
            kappa = (t(1)*tt(2) - t(2)*tt(1))/norm2(t)**3
            if (kappa < 0) reach = min(reach, 0.9_dp/(-kappa))
            gap = min(gap, x(1, j) - box(1), box(2) - x(1, j), x(2, j) - box(3), &
                box(4) - x(2, j))
            do k = 1, size(obstacles, 2)
                gap = min(gap, norm2(obstacles(:, k) - x(:, j)))
            end do
        end do
        distance = min(reach, gap/2, cap)
        if (.not. distance > 0) return
        status = IDEAL_FB_OK
    end subroutine

    function section_point(boundary, centre, axes, alpha) result(x)
        ! Mean (n = 0) section of the guess at reference angle alpha.
        real(dp), intent(in) :: boundary(:, :, 0:), centre(2), axes(2), alpha
        real(dp) :: x(2), values(3, size(boundary, 2))
        call ideal_seed_boundary_coefficients(centre + axes*[cos(alpha), sin(alpha)], &
            boundary, values, centre=centre, axes=axes)
        x = values(1:2, 1)
    end function

    subroutine fb_nested_layer_mesh(boundary, centre, axes, nodes, degree, layers, &
                                    distance, mesh, chain, positions, lines, line_t, &
                                    status, spines)
        !! Reference mesh of the plasma disc (the reference ellipse, nodes
        !! boundary nodes) and a layered shell around it for
        !! ideal_fb_init(positions, lines, line_t). Every Gamma_p chain node y_p
        !! (the disc's boundary nodes, degree per edge) starts a segment to its
        !! Gamma_c node; the shell nodes lie on these segments at t = k/(degree
        !! layers), and the cells are the two triangles of every (edge, layer)
        !! quadrilateral. Reference: y = y_p + t delta nu (nu the ellipse
        !! normal, delta = min(axes)/4). Physical rest positions: Gamma_c is
        !! the guess's mean section X(y_p) offset by distance along its outward
        !! normal, the interior nodes x = (1 - t) y_p + t X_c (only the
        !! constraint rows read them). chain: all Gamma_p chain nodes, for the
        !! plasma-cell test of ideal_fb_init.
        real(dp), intent(in) :: boundary(:, :, 0:), centre(2), axes(2), distance
        integer, intent(in) :: nodes, degree, layers
        type(gs_state_t), intent(out) :: mesh
        real(dp), allocatable, intent(out) :: chain(:, :), positions(:, :), line_t(:)
        integer, allocatable, intent(out) :: lines(:, :)
        integer, intent(out) :: status
        !> outward normal of the guess's section at every Gamma_p node: the
        !! spine along which ideal_fb_init lets that node move
        real(dp), allocatable, intent(out), optional :: spines(:, :)
        type(gs_state_t) :: disc
        integer, allocatable :: ring(:), lattice(:, :), tri(:, :)
        real(dp), allocatable :: points(:, :), xc(:, :)
        real(dp) :: nu(2), delta, xp(2), tangent(2), alpha, t, area
        real(dp), parameter :: eps = 1.0e-6_dp
        integer :: np, nring, nr, i, r, k, e, l, nloc, ncell, ntri, a, b, next
        integer :: local(2, 10), corner(2, 3), j
        status = IDEAL_FB_BAD_INPUT
        if (layers < 1 .or. (degree /= 2 .and. degree /= 3)) return
        if (.not. distance > 0) return
        call ideal_seed_disk_mesh(disc, nodes, degree, status)
        if (status /= 0) return
        status = IDEAL_FB_BAD_INPUT
        np = size(disc%points, 2); nloc = size(disc%triangles, 1)
        ! the seed's disk mesh is the unit disc about (2, 0)
        do i = 1, np
            disc%points(:, i) = centre + axes*(disc%points(:, i) - [2.0_dp, 0.0_dp])
        end do
        call boundary_ring(disc, degree, ring, status)
        if (status /= 0) return
        status = IDEAL_FB_BAD_INPUT
        nring = size(ring); nr = degree*layers
        ! Gamma_c of every chain node: offset of the mean section along its
        ! outward normal at the same reference angle.
        allocate (xc(2, nring))
        do i = 1, nring
            alpha = atan2((disc%points(2, ring(i)) - centre(2))/axes(2), &
                (disc%points(1, ring(i)) - centre(1))/axes(1))
            xp = section_point(boundary, centre, axes, alpha)
            tangent = section_point(boundary, centre, axes, alpha + eps) &
                - section_point(boundary, centre, axes, alpha - eps)
            xc(:, i) = xp + distance*[tangent(2), -tangent(1)]/norm2(tangent)
        end do
        ! Lattice (chain node, layer node) -> mesh node; layer 0 is Gamma_p.
        allocate (lattice(nring, 0:nr), points(2, np + nring*nr))
        points(:, :np) = disc%points
        lattice(:, 0) = ring
        delta = minval(axes)/4
        k = np
        do r = 1, nr
            do i = 1, nring
                k = k + 1
                lattice(i, r) = k
                nu = (disc%points(:, ring(i)) - centre)/axes**2
                nu = nu/norm2(nu)
                points(:, k) = disc%points(:, ring(i)) + real(r, dp)/nr*delta*nu
            end do
        end do
        allocate (positions(2, k), lines(2, k), line_t(k))
        if (present(spines)) then
            allocate (spines(2, k))
            spines = 0
            do i = 1, nring
                spines(:, ring(i)) = xc(:, i) - section_point(boundary, centre, axes, &
                    atan2((disc%points(2, ring(i)) - centre(2))/axes(2), &
                    (disc%points(1, ring(i)) - centre(1))/axes(1)))
            end do
        end if
        positions = points; lines = 0; line_t = -1
        do r = 1, nr
            t = real(r, dp)/nr
            do i = 1, nring
                lines(:, lattice(i, r)) = [ring(i), lattice(i, nr)]
                line_t(lattice(i, r)) = t
                positions(:, lattice(i, r)) = (1 - t)*disc%points(:, ring(i)) &
                    + t*xc(:, i)
            end do
        end do
        ! Local nodes of the P2/P3 triangle in lattice steps from its
        ! vertices: vertices, degree - 1 nodes per edge (1-2, 2-3, 3-1), then
        ! the P3 centroid.
        ncell = size(disc%triangles, 2) + 2*(nring/degree)*layers
        allocate (tri(nloc, ncell))
        tri(:, :size(disc%triangles, 2)) = disc%triangles
        ntri = size(disc%triangles, 2)
        do e = 1, nring/degree
            do l = 1, layers
                i = (e - 1)*degree; r = (l - 1)*degree
                do j = 1, 2
                    if (j == 1) then
                        corner = reshape([0, 0, degree, 0, degree, degree], [2, 3])
                    else
                        corner = reshape([0, 0, degree, degree, 0, degree], [2, 3])
                    end if
                    ! counter-clockwise in the reference
                    area = cross2(lattice_point(corner(:, 2)) - lattice_point(corner(:, 1)), &
                        lattice_point(corner(:, 3)) - lattice_point(corner(:, 1)))
                    if (area < 0) corner = corner(:, [1, 3, 2])
                    local(:, 1:3) = corner
                    k = 3
                    do a = 1, 3
                        next = mod(a, 3) + 1
                        do b = 1, degree - 1
                            k = k + 1
                            local(:, k) = corner(:, a) + b*(corner(:, next) - corner(:, a)) &
                                /degree
                        end do
                    end do
                    if (degree == 3) local(:, 10) = sum(corner, dim=2)/3
                    ntri = ntri + 1
                    do k = 1, nloc
                        tri(k, ntri) = lattice(modulo(i + local(1, k), nring) + 1, &
                            r + local(2, k))
                    end do
                end do
            end do
        end do
        mesh%points = points
        mesh%triangles = tri
        allocate (mesh%boundary(size(points, 2)))
        mesh%boundary = .false.
        mesh%boundary(lattice(:, nr)) = .true.
        call mesh_neighbors(tri, mesh%neighbors)
        chain = disc%points(:, ring)
        status = IDEAL_FB_OK
    contains
        function lattice_point(c) result(y)
            integer, intent(in) :: c(2)
            real(dp) :: y(2)
            y = points(:, lattice(modulo(i + c(1), nring) + 1, r + c(2)))
        end function
    end subroutine

    pure real(dp) function cross2(a, b)
        real(dp), intent(in) :: a(2), b(2)
        cross2 = a(1)*b(2) - a(2)*b(1)
    end function

    subroutine boundary_ring(mesh, degree, ring, status)
        ! Boundary nodes of a P2/P3 mesh in counter-clockwise chain order,
        ! degree per boundary edge (edge k of a cell runs from vertex k to
        ! k + 1; neighbors(k, e) = 0 on the boundary).
        type(gs_state_t), intent(in) :: mesh
        integer, intent(in) :: degree
        integer, allocatable, intent(out) :: ring(:)
        integer, intent(out) :: status
        integer, allocatable :: first(:)
        integer :: e, k, nb, i, cell, side, start, j
        status = 1
        allocate (first(size(mesh%points, 2)))
        first = 0; nb = 0
        do e = 1, size(mesh%triangles, 2)
            do k = 1, 3
                if (mesh%neighbors(k, e) > 0) cycle
                first(mesh%triangles(k, e)) = 3*(e - 1) + k
                nb = nb + 1
            end do
        end do
        if (nb < 3) return
        allocate (ring(nb*degree))
        start = findloc(first > 0, .true., dim=1)
        i = start
        do j = 1, nb
            cell = (first(i) - 1)/3 + 1; side = first(i) - 3*(cell - 1)
            ring((j - 1)*degree + 1) = mesh%triangles(side, cell)
            do k = 1, degree - 1
                ring((j - 1)*degree + 1 + k) = &
                    mesh%triangles(3 + (side - 1)*(degree - 1) + k, cell)
            end do
            i = mesh%triangles(mod(side, 3) + 1, cell)
            if (first(i) == 0) return
        end do
        if (i /= start) return
        status = 0
    end subroutine

    subroutine mesh_neighbors(triangles, neighbors)
        ! neighbors(k, e): the cell across edge (vertex k, vertex k + 1), 0 on
        ! the boundary.
        integer, intent(in) :: triangles(:, :)
        integer, allocatable, intent(out) :: neighbors(:, :)
        integer(int64), allocatable :: keys(:)
        integer, allocatable :: order(:)
        integer :: e, k, a, b, n, i, c1, c2
        n = size(triangles, 2)
        allocate (keys(3*n), order(3*n), neighbors(3, n))
        neighbors = 0
        do e = 1, n
            do k = 1, 3
                a = triangles(k, e); b = triangles(mod(k, 3) + 1, e)
                keys(3*(e - 1) + k) = int(min(a, b), int64)*2_int64**31 &
                    + max(a, b)
            end do
        end do
        call sort_key_order(keys, order)
        do i = 1, 3*n - 1
            if (keys(order(i)) /= keys(order(i + 1))) cycle
            c1 = order(i); c2 = order(i + 1)
            neighbors(mod(c1 - 1, 3) + 1, (c1 - 1)/3 + 1) = (c2 - 1)/3 + 1
            neighbors(mod(c2 - 1, 3) + 1, (c2 - 1)/3 + 1) = (c1 - 1)/3 + 1
        end do
    end subroutine
end module milca_ideal_fb_nested
