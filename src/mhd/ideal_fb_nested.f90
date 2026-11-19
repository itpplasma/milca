module milca_ideal_fb_nested
    !! Set-up of the nested-map free boundary from a VMEC-style boundary guess
    !! (docs/free-boundary-nested-map-3d.md, section 10): the reference disc
    !! of the frozen field, the axisymmetric coupling curve Gamma_c, the mesh
    !! of the region inside it, the seed of the plasma map and its transfer
    !! between resolution levels. The physics stays in
    !! milca_ideal_free_boundary; this module only places and seeds it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_grad_shafranov, only: gs_state_t, initialize_gs
    use milca_ideal_free_boundary, only: ideal_fb_t, ideal_fb_field_basis_t, &
        ideal_fb_solve_shell_map, ideal_fb_solve_vacuum, ideal_fb_solve_multipliers, &
        IDEAL_FB_OK, IDEAL_FB_BAD_INPUT, IDEAL_FB_FOLDED
    use milca_ideal_seed_3d, only: ideal_seed_boundary_coefficients
    use milca_winslow_seed, only: winslow_seed, WINSLOW_OK
    use milca_curved_triangle, only: triangle_locate
    use milca_triangle_basis, only: triangle_basis_jet
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    implicit none
    private
    public :: fb_nested_geometry, fb_nested_mesh, fb_nested_seed, fb_nested_transfer
    public :: fb_nested_settle
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
                                  curve, outline_mean, margin_used, status, message)
        !! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp
        !! phi), counter-clockwise in (R, Z). centre and axes: the elliptic
        !! reference disc of the frozen field, inscribed in the mean bounding
        !! box of the sections. outline_mean(k): the mean over the sections of
        !! their distance from the centre along the ray alpha = 2 pi (k - 1)/rays,
        !! the reference of Gamma_p in the shell (fb_nested_mesh).
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
        real(dp), allocatable, intent(out) :: curve(:), outline_mean(:)
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
        call outline(core, outline_mean)
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
            ! section point
            real(dp), intent(in) :: clearance
            real(dp), intent(out) :: r(:)
            real(dp) :: e(2), q(2), along, across
            integer :: k, i, j
            r = 0
            do k = 1, rays
                e = [cos(2*pi*(k - 1)/rays), sin(2*pi*(k - 1)/rays)]
                do j = 1, nphi
                    do i = 1, ntheta
                        q = rz(:, i, j) - centre
                        along = dot_product(q, e); across = q(2)*e(1) - q(1)*e(2)
                        if (abs(across) < clearance) &
                            r(k) = max(r(k), along + sqrt(clearance**2 - across**2))
                    end do
                end do
            end do
        end subroutine

        subroutine outline(r, mean)
            ! r(alpha) of the farthest crossing of the ray with every section,
            ! and the mean over the sections of their farthest crossing
            real(dp), intent(out) :: r(:)
            real(dp), allocatable, intent(out) :: mean(:)
            real(dp) :: e(2), a(2), b(2), f(2), det, along, s, here
            integer :: k, i, j
            r = 0
            allocate (mean(rays))
            mean = 0
            do k = 1, rays
                e = [cos(2*pi*(k - 1)/rays), sin(2*pi*(k - 1)/rays)]
                do j = 1, nphi
                    here = 0
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
                        here = max(here, along)
                    end do
                    r(k) = max(r(k), here)
                    mean(k) = mean(k) + here/nphi
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

    subroutine fb_nested_mesh(mesh, chain, plasma_points, centre, axes, curve, &
                              outline_mean, nodes, degree, status)
        !! Mesh of the region inside Gamma_c in the shell reference, where
        !! Gamma_p is the mean section outline_mean of the guess (nodes chain
        !! points). The plasma cells get their own reference plasma_points,
        !! the elliptic disc of the frozen field, through the radial map
        !!     G(y) = centre + (y - centre) rho_s(alpha)/rho_e(alpha)
        !! between the ellipse rho_e and the mean section rho_s at the polar
        !! angle alpha about the centre. The harmonic extension of the shell
        !! map then starts from the guess itself instead of carrying the
        !! displacement from the ellipse to the section across the shell,
        !! which folds a shell thinner than that displacement. The cells are
        !! at most the chain spacing and twice the narrowest gap between
        !! Gamma_p and Gamma_c; Gamma_c gets the same spacing.
        type(gs_state_t), intent(out) :: mesh
        real(dp), allocatable, intent(out) :: chain(:, :), plasma_points(:, :)
        real(dp), intent(in) :: centre(2), axes(2), curve(:), outline_mean(:)
        integer, intent(in) :: nodes, degree
        integer, intent(out) :: status
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h, length, alpha, gap, y(2)
        integer :: j, nc, nchain, a
        h = 0
        do j = 1, rays
            alpha = 2*pi*(j - 1)/rays
            h = h + norm2(on_section(alpha + 2*pi/rays) - on_section(alpha))/nodes
        end do
        length = 0; gap = huge(1.0_dp)
        do j = 1, rays
            alpha = 2*pi*(j - 1)/rays
            length = length + norm2(point(alpha) - point(alpha + 2*pi/rays))
            gap = min(gap, curve_radius(curve, alpha) - section_radius(alpha))
        end do
        ! a chain segment longer than the gap would let its curved edge bulge
        ! to Gamma_c: the chain gets at least one node per 2 gaps
        nchain = max(nodes, 4*ceiling(nodes*h/(8*gap)))
        h = min(h*nodes/nchain, 2*gap)
        nc = max(nchain, ceiling(length/h))
        allocate (chain(2, nchain), contour(2, nc))
        do j = 1, nchain
            ! the image of the ellipse node at parameter 2 pi (j - 1)/nchain
            alpha = 2*pi*(j - 1)/nchain
            y = axes*[cos(alpha), sin(alpha)]
            chain(:, j) = on_section(atan2(y(2), y(1)))
        end do
        do j = 1, nc
            contour(:, j) = point(2*pi*(j - 1)/nc)
        end do
        call initialize_gs(mesh, contour, h*h/2, status, degree, curved=.true., &
            interior_chain=chain, geometry_only=.true.)
        if (status /= 0) return
        allocate (plasma_points, mold=mesh%points)
        do a = 1, size(mesh%points, 2)
            y = mesh%points(:, a) - centre
            alpha = atan2(y(2), y(1))
            plasma_points(:, a) = centre + y/(section_radius(alpha)*norm2([cos(alpha), &
                sin(alpha)]/axes))
        end do
    contains
        function point(a) result(p)
            real(dp), intent(in) :: a
            real(dp) :: p(2)
            p = centre + curve_radius(curve, a)*[cos(a), sin(a)]
        end function

        function on_section(a) result(p)
            real(dp), intent(in) :: a
            real(dp) :: p(2)
            p = centre + section_radius(a)*[cos(a), sin(a)]
        end function

        real(dp) function section_radius(a) result(r)
            ! the trigonometric interpolant of the ray samples, smooth in alpha
            real(dp), intent(in) :: a
            integer :: k, n
            n = size(outline_mean)
            r = 0
            do k = 1, n
                r = r + outline_mean(k)*dirichlet(a - 2*pi*(k - 1)/n, n)
            end do
        end function
    end subroutine

    pure real(dp) function dirichlet(x, n) result(w)
        ! cardinal function of trigonometric interpolation on n (even)
        ! equispaced points, the Nyquist mode split evenly
        real(dp), intent(in) :: x
        integer, intent(in) :: n
        if (abs(sin(x/2)) < 1e-12_dp) then
            w = 1
            return
        end if
        w = sin(n*x/2)*cos(x/2)/(n*sin(x/2))
    end function

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
end module milca_ideal_fb_nested
