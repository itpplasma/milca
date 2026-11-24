module flux_map_sections
    !! Section meshes, start values and diagnostics for the 3D flux map
    !! (milca_flux_map_3d) from the wall curves of the toroidal planes only:
    !! no magnetic axis enters. Each plane is the pullback of one O-grid of
    !! the unit disc by the harmonic map of its section onto the disc
    !! (fortfem_harmonic_disc: valid for any smooth simple closed wall
    !! curve), and the start is that map, u = G.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    use fortsparse_csr, only: csr_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_cell_broken, spline_space_build, SPLINE_C0
    use fortfem_harmonic_disc, only: harmonic_disc_t, harmonic_disc_build, &
        harmonic_disc_inverse, harmonic_disc_pull_back, harmonic_disc_wall_error, &
        spline_jacobian_certificate
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_point_matrix, spline_locate, spline_boundary_nodes, &
        spline_side_quadrature_t, spline_interface_quadrature, &
        spline_cell_tables, spline_cell_xi_jets, spline_point_jet_matrix
    use fortfem_spline_layout, only: LAYOUT_INTERIOR
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_boundary_dofs, spline_derham_project, spline_derham_t_basis
    use milca_eulerian_clebsch, only: clebsch_node_gradients, clebsch_nfp
    use milca_spline_ale, only: spline_geometry_dofs, spline_set_geometry
    use milca_flux_map, only: flux_map_profiles
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_t, flux_map3d_init, &
        flux_map3d_point, flux_map3d_expand_symmetric
    implicit none
    private

    type, public :: fm3_case_t
        type(spline_space_t) :: space ! geometry of plane 1
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad
        type(flux_map3d_t) :: fm
        integer :: n0 = 0, nt = 0, cells = 0, degree = 0
        real(dp), allocatable :: base(:, :, :) ! node coordinates per plane
        real(dp), allocatable :: centre(:, :) ! mesh centres, G^-1(0)
        real(dp), allocatable :: samples(:, :, :)
        real(dp), allocatable :: disc(:, :) ! control points of the disc O-grid
        ! per plane: distance of the mesh wall from the wall curve, relative
        ! error of the section area, smallest/largest Bernstein coefficient of
        ! the block Jacobians and whether all are positive (mesh valid)
        real(dp), allocatable :: wall_distance(:), area_error(:), jacobian_ratio(:)
        logical, allocatable :: certified(:)
        real(dp) :: map_residual = 0 ! boundary error of the harmonic maps
        ! fm3_start_polar: area centroids and sqrt of the second moments
        real(dp), allocatable :: centroid(:, :), frame(:, :, :)
    end type fm3_case_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    public :: fm3_build, fm3_start, fm3_start_polar, fm3_axis, fm3_mesh_quality, &
              fm3_locate, fm3_state, fm3_flux, fm3_prolong, fm3_preimage, &
              fm3_dump_hessian, fm3_wall_winding, fm3_interface_jump

    ! state of the projection callback of fm3_prolong
    type(fm3_case_t), pointer :: old => null()
    real(dp), allocatable :: old_x(:), old_y(:)
    integer :: old_component = 1

contains

    subroutine fm3_build(c, samples, cells, degree, extra_points, status, preserve_area, grading, match_sample_angle, &
                         radial_cells)
        !! Meshes of all planes from the wall curves samples(2, m, nt):
        !! counter-clockwise, equally spaced in a smooth periodic parameter,
        !! the first sample (mapped to the disc angle 0) a smooth function of
        !! the plane. preserve_area (default true): the mesh wall encloses the
        !! area of the wall curve. grading (default 0): strength of the
        !! curvature weight of the boundary correspondence of the harmonic
        !! maps (more wall cells at tips). match_sample_angle preserves the sample angle.
        !! radial_cells sets the radial count of the outer blocks.
        !! Status 10 + s: status s of the harmonic map
        !! of a plane (3: the samples do not resolve the curve), 20: a mesh is
        !! not certified valid.
        type(fm3_case_t), intent(out) :: c
        real(dp), intent(in) :: samples(:, :, :)
        integer, intent(in) :: cells, degree, extra_points
        integer, intent(out) :: status
        logical, intent(in), optional :: preserve_area
        real(dp), intent(in), optional :: grading
        logical, intent(in), optional :: match_sample_angle
        integer, intent(in), optional :: radial_cells
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(harmonic_disc_t) :: map
        real(dp), allocatable :: xg(:, :)
        integer, allocatable :: wall(:)
        complex(dp) :: z
        logical :: keep
        integer :: k

        c%nt = size(samples, 3)
        c%cells = cells
        c%degree = degree
        c%samples = samples
        keep = .true.
        if (present(preserve_area)) keep = preserve_area
        allocate (c%centre(2, c%nt), c%wall_distance(c%nt), c%area_error(c%nt), &
                  c%jacobian_ratio(c%nt), c%certified(c%nt))
        call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, core_fraction(), layout, status)
        if (status /= 0) return
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        if (present(radial_cells)) options%knots(2, 2:)%element_count = radial_cells
        call spline_space_build(layout, options, c%space, status)
        if (status /= 0) return
        call spline_geometry_dofs(c%space, c%disc, status)
        if (status /= 0) return
        allocate (c%base(2, size(c%disc, 2), c%nt))
        c%map_residual = 0
        do k = c%nt, 1, -1
            call harmonic_disc_build(map, samples(:, :, k), status, grading, match_sample_angle)
            if (status /= 0) then
                status = 10 + status
                return
            end if
            c%map_residual = max(c%map_residual, map%residual)
            call harmonic_disc_pull_back(c%space, map, status, preserve_area=keep)
            if (status /= 0) return
            call harmonic_disc_wall_error(c%space, map, c%wall_distance(k), c%area_error(k), &
                                          status)
            if (status /= 0) return
            call spline_jacobian_certificate(c%space, c%certified(k), c%jacobian_ratio(k), status)
            if (status /= 0) return
            call spline_geometry_dofs(c%space, xg, status)
            if (status /= 0) return
            c%base(:, :, k) = xg
            call harmonic_disc_inverse(map, (0.0_dp, 0.0_dp), z, status)
            if (status /= 0) return
            c%centre(:, k) = [real(z, dp), aimag(z)]
        end do
        status = 20
        if (.not. all(c%certified)) return
        call spline_derham_build(c%space, c%derham, status, 0)
        if (status /= 0) return
        call spline_quadrature_build(c%space, degree + extra_points, c%quad, status)
        if (status /= 0) return
        c%n0 = c%derham%dimension(0)
        call spline_derham_boundary_dofs(c%derham, LAYOUT_WALL, 0, wall)
        call flux_map3d_init(c%fm, c%n0, c%nt, wall)
    end subroutine fm3_build

    real(dp) function core_fraction()
        !! Size of the O-grid core relative to the disc radius (0.5; diagnostic
        !! override MILCA_FM3_CORE).
        character(32) :: v
        integer :: length
        core_fraction = 0.5_dp
        call get_environment_variable('MILCA_FM3_CORE', v, length)
        if (length > 0) read (v, *) core_fraction
    end function core_fraction

    subroutine fm3_start(c)
        !! u = G on every plane: the control points of the disc O-grid; a wall
        !! function takes their angle.
        type(fm3_case_t), intent(inout) :: c
        integer :: i, k, ntot

        ntot = c%n0*c%nt
        do k = 1, c%nt
            do i = 1, c%n0
                if (c%fm%is_wall(i)) then
                    c%fm%z((k - 1)*c%n0 + i) = atan2(c%disc(2, i), c%disc(1, i))
                else
                    c%fm%z((k - 1)*c%n0 + i) = c%disc(1, i)
                    c%fm%z(ntot + (k - 1)*c%n0 + i) = c%disc(2, i)
                end if
            end do
        end do
    end subroutine fm3_start

    pure subroutine moments(poly, centroid, root)
        !! Area centroid and the symmetric square root of the second-moment
        !! tensor of a closed polygon.
        real(dp), intent(in) :: poly(:, :)
        real(dp), intent(out) :: centroid(2), root(2, 2)
        real(dp) :: a, cr, p(2), q(2), m(2, 2), tr, det, sd
        integer :: i, j, n

        n = size(poly, 2)
        a = 0
        centroid = 0
        do i = 1, n
            j = mod(i, n) + 1
            cr = poly(1, i)*poly(2, j) - poly(1, j)*poly(2, i)
            a = a + cr/2
            centroid = centroid + (poly(:, i) + poly(:, j))*cr/6
        end do
        centroid = centroid/a
        m = 0
        do i = 1, n
            j = mod(i, n) + 1
            p = poly(:, i) - centroid
            q = poly(:, j) - centroid
            cr = p(1)*q(2) - q(1)*p(2)
            m(1, 1) = m(1, 1) + cr*(p(1)**2 + p(1)*q(1) + q(1)**2)/12
            m(2, 2) = m(2, 2) + cr*(p(2)**2 + p(2)*q(2) + q(2)**2)/12
            m(1, 2) = m(1, 2) + cr*(p(1)*q(2) + 2*p(1)*p(2) + 2*q(1)*q(2) + q(1)*p(2))/24
        end do
        m(2, 1) = m(1, 2)
        m = m/a
        ! sqrt of a symmetric positive 2 x 2 matrix
        tr = m(1, 1) + m(2, 2)
        det = m(1, 1)*m(2, 2) - m(1, 2)**2
        sd = sqrt(det)
        root = m
        root(1, 1) = root(1, 1) + sd
        root(2, 2) = root(2, 2) + sd
        root = root/sqrt(tr + 2*sd)
    end subroutine moments

    pure function normalised(c, k, p) result(v)
        !! frame^-1 (p - centre) of plane k.
        type(fm3_case_t), intent(in) :: c
        integer, intent(in) :: k
        real(dp), intent(in) :: p(2)
        real(dp) :: v(2), d(2), det

        associate (s => c%frame(:, :, k))
            d = p - c%centroid(:, k)
            det = s(1, 1)*s(2, 2) - s(1, 2)*s(2, 1)
            v = [s(2, 2)*d(1) - s(1, 2)*d(2), -s(2, 1)*d(1) + s(1, 1)*d(2)]/det
        end associate
    end function normalised

    pure real(dp) function wall_radius(c, k, angle) result(rho)
        !! Radius of the wall curve of plane k in the normalised frame along
        !! the ray at angle (the curve star-shaped about the centre there):
        !! intersection with the sample polygon.
        type(fm3_case_t), intent(in) :: c
        integer, intent(in) :: k
        real(dp), intent(in) :: angle
        real(dp) :: p(2), q(2), e(2), den, t, u
        integer :: i, j, n

        n = size(c%samples, 2)
        e = [cos(angle), sin(angle)]
        rho = 0
        do i = 1, n
            j = mod(i, n) + 1
            p = normalised(c, k, c%samples(:, i, k))
            q = normalised(c, k, c%samples(:, j, k)) - p
            den = e(1)*q(2) - e(2)*q(1)
            if (abs(den) < 1e-300_dp) cycle
            t = (p(1)*q(2) - p(2)*q(1))/den
            u = (p(1)*e(2) - p(2)*e(1))/den
            if (t > 0 .and. u >= 0 .and. u < 1) rho = max(rho, t)
        end do
    end function wall_radius

    subroutine fm3_start_polar(c)
        !! Alternative start (gate F2b, whose exact state is a saddle that the
        !! start u = G does not reach): u = normalised position / wall radius
        !! at the control points of every plane, about the area centroid
        !! after an affine normalisation by the second moments of the
        !! section (star-shaped sections only); Theta = the normalised polar
        !! angle on the wall, continuous in t.
        type(fm3_case_t), intent(inout) :: c
        real(dp) :: v(2), theta, prev
        integer :: i, k, ntot

        if (.not. allocated(c%centroid)) then
            allocate (c%centroid(2, c%nt), c%frame(2, 2, c%nt))
            do k = 1, c%nt
                call moments(c%samples(:, :, k), c%centroid(:, k), c%frame(:, :, k))
            end do
        end if
        ntot = c%n0*c%nt
        do i = 1, c%n0
            prev = 0
            do k = 1, c%nt
                v = normalised(c, k, c%base(:, i, k))
                theta = atan2(v(2), v(1))
                if (c%fm%is_wall(i)) then
                    if (k > 1) theta = theta + two_pi*nint((prev - theta)/two_pi)
                    prev = theta
                    c%fm%z((k - 1)*c%n0 + i) = theta
                else
                    c%fm%z((k - 1)*c%n0 + i) = v(1)/wall_radius(c, k, theta)
                    c%fm%z(ntot + (k - 1)*c%n0 + i) = v(2)/wall_radius(c, k, theta)
                end if
            end do
        end do
    end subroutine fm3_start_polar

    subroutine fm3_locate(c, x0, block, d0, nv, gn, status)
        !! Section functions at the point x0 of the reference plane 1: their
        !! indices, values and gradients in the geometry of plane 1.
        type(fm3_case_t), intent(in) :: c
        real(dp), intent(in) :: x0(2)
        integer, intent(inout) :: block
        integer, allocatable, intent(out) :: d0(:)
        real(dp), allocatable, intent(out) :: nv(:), gn(:, :)
        integer, intent(out) :: status
        integer, allocatable :: rowptr(:)
        real(dp), allocatable :: jets(:, :)
        real(dp) :: xi(2)
        integer :: located

        call spline_locate(c%space, block, x0, located, xi, status)
        if (status /= 0) return
        block = located
        call spline_point_jet_matrix(c%space, [block], reshape(xi, [2, 1]), 1, rowptr, d0, jets, status)
        if (status /= 0) return
        nv = jets(1, :)
        gn = jets(2:3, :)
    end subroutine fm3_locate

    subroutine fm3_state(c, x0, block, t, profiles, rz, u, b, status, label_jacobian)
        !! Position, u and B at the point with reference-plane coordinates x0
        !! and toroidal parameter t.
        type(fm3_case_t), intent(in) :: c
        real(dp), intent(in) :: x0(2), t
        integer, intent(inout) :: block
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: rz(2), u(8), b(3)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: label_jacobian(2, 2)
        integer, allocatable :: d0(:)
        real(dp), allocatable :: nv(:), gn(:, :)

        call fm3_locate(c, x0, block, d0, nv, gn, status)
        if (status /= 0) return
        call flux_map3d_point(c%fm, c%base, d0, nv, gn, t, profiles, rz, u, b, label_jacobian)
    end subroutine fm3_state

    subroutine fm3_axis(c, t, profiles, axis, status, x0)
        !! Magnetic axis at the toroidal parameter t: the zero of u_h, Newton
        !! in the coordinates of the reference plane from its mesh centre.
        type(fm3_case_t), intent(in) :: c
        real(dp), intent(in) :: t
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: axis(2)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: x0(2)
        integer, allocatable :: d0(:)
        real(dp), allocatable :: nv(:), gn(:, :)
        real(dp) :: x(2), uv(2), jac(2, 2), det, u(8), b(3), lt(c%nt, 2), ht(c%nt, 2), &
                    zx(c%n0), zy(c%n0)
        integer :: it, block, ntot, k

        ntot = c%n0*c%nt
        call spline_derham_t_basis(c%nt, t, lt, ht)
        zx = 0
        zy = 0
        do k = 1, c%nt
            zx = zx + lt(k, 1)*c%fm%z((k - 1)*c%n0 + 1:k*c%n0)
            zy = zy + lt(k, 1)*c%fm%z(ntot + (k - 1)*c%n0 + 1:ntot + k*c%n0)
        end do
        x = c%centre(:, 1)
        block = 1
        do it = 1, 40
            call fm3_locate(c, x, block, d0, nv, gn, status)
            if (status /= 0) return
            status = 4
            if (any(c%fm%is_wall(d0) .and. abs(nv) > 1e-14_dp)) return
            uv = [dot_product(nv, zx(d0)), dot_product(nv, zy(d0))]
            jac(1, :) = matmul(gn, zx(d0))
            jac(2, :) = matmul(gn, zy(d0))
            if (norm2(uv) < 1e-13_dp) exit
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            x = x - [jac(2, 2)*uv(1) - jac(1, 2)*uv(2), jac(1, 1)*uv(2) - jac(2, 1)*uv(1)]/det
        end do
        status = 5
        if (norm2(uv) >= 1e-10_dp) return
        call flux_map3d_point(c%fm, c%base, d0, nv, gn, t, profiles, axis, u, b)
        if (present(x0)) x0 = x
        status = 0
    end subroutine fm3_axis

    subroutine fm3_preimage(c, t, profiles, target, x, rz, b, status, block_hint, inverse_jacobian)
        !! The point of the plane t that u_h maps to target: Newton in the
        !! coordinates x of the reference plane from
        !! the guess x. Returns its position rz and B. The optional hint is
        !! updated to the converged block. The optional inverse maps a label
        !! displacement to a reference-plane displacement.
        type(fm3_case_t), intent(in) :: c
        real(dp), intent(in) :: t, target(2)
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(inout) :: x(2)
        real(dp), intent(out) :: rz(2), b(3)
        integer, intent(out) :: status
        integer, intent(inout), optional :: block_hint
        real(dp), intent(out), optional :: inverse_jacobian(2, 2)
        real(dp) :: u(8), r(2), jac(2, 2), inverse(2, 2), det, dx(2), lim
        integer :: it, block

        block = 1
        if (present(block_hint)) block = block_hint
        lim = 0.05_dp*(maxval(c%samples(1, :, 1)) - minval(c%samples(1, :, 1)))
        do it = 1, 60
            call fm3_state(c, x, block, t, profiles, rz, u, b, status, jac)
            if (status /= 0) return
            r = u(1:2) - target
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            inverse = 0
            if (abs(det) > tiny(det)) &
                inverse = reshape([jac(2, 2), -jac(2, 1), -jac(1, 2), jac(1, 1)], [2, 2])/det
            if (norm2(r) < 1e-11_dp) then
                if (present(block_hint)) block_hint = block
                if (present(inverse_jacobian)) inverse_jacobian = inverse
                return
            end if
            if (abs(det) <= tiny(det)) exit
            dx = -matmul(inverse, r)
            ! steps limited to a tenth of the section size
            if (norm2(dx) > lim) dx = dx*lim/norm2(dx)
            x = x + dx
        end do
        status = 6
    end subroutine fm3_preimage

    subroutine fm3_prolong(c0, c, status)
        !! Start of c from the state of c0 (the same wall curves and block
        !! layout; other cells, degree or planes): trigonometric interpolation
        !! in t, then per plane the V0 projection of u_h of c0 at the same
        !! block parameters. A wall function takes the angle of its projected
        !! coefficients.
        type(fm3_case_t), intent(in), target :: c0
        type(fm3_case_t), intent(inout) :: c
        integer, intent(out) :: status
        real(dp) :: lt(c0%nt, 2), ht(c0%nt, 2), theta, prev(c%n0)
        real(dp), allocatable :: cx(:), cy(:)
        integer :: k, k0, i, n0, m0

        n0 = c%n0
        m0 = c0%n0
        old => c0
        allocate (cx(n0), cy(n0), old_x(m0), old_y(m0))
        prev = 0
        do k = 1, c%nt
            call spline_derham_t_basis(c0%nt, (k - 1)/real(c%nt, dp), lt, ht)
            old_x = 0
            old_y = 0
            do k0 = 1, c0%nt
                old_x = old_x + lt(k0, 1)*c0%fm%z((k0 - 1)*m0 + 1:k0*m0)
                old_y = old_y + lt(k0, 1)*c0%fm%z(m0*c0%nt + (k0 - 1)*m0 + 1:m0*c0%nt + k0*m0)
            end do
            old_component = 1
            call spline_derham_project(c%derham, 0, old_value, cx, status, c%degree + 3)
            if (status /= 0) return
            old_component = 2
            call spline_derham_project(c%derham, 0, old_value, cy, status, c%degree + 3)
            if (status /= 0) return
            do i = 1, n0
                if (c%fm%is_wall(i)) then
                    theta = atan2(cy(i), cx(i))
                    if (k > 1) theta = theta + two_pi*nint((prev(i) - theta)/two_pi)
                    prev(i) = theta
                    c%fm%z((k - 1)*n0 + i) = theta
                else
                    c%fm%z((k - 1)*n0 + i) = cx(i)
                    c%fm%z(n0*c%nt + (k - 1)*n0 + i) = cy(i)
                end if
            end do
        end do
        deallocate (old_x, old_y)
        old => null()
    end subroutine fm3_prolong

    subroutine old_value(block, xi, t, values)
        !! Component old_component of u_h of the old case at (block, xi).
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:)
        real(dp) :: sig, acc, ref, tb, u(2)
        integer :: i, st
        logical :: first

        call spline_point_matrix(old%space, [block], reshape(xi, [2, 1]), 0, rowptr, colidx, &
                                 val, st)
        if (st /= 0) error stop 'fm3_prolong: old point'
        u = 0
        sig = 0
        acc = 0
        first = .true.
        ref = 0
        do i = 1, size(colidx)
            if (old%fm%is_wall(colidx(i))) then
                if (first) ref = old_x(colidx(i))
                first = .false.
                tb = old_x(colidx(i)) + two_pi*nint((ref - old_x(colidx(i)))/two_pi)
                sig = sig + val(i)
                acc = acc + val(i)*tb
            else
                u = u + val(i)*[old_x(colidx(i)), old_y(colidx(i))]
            end if
        end do
        if (sig > 0) u = u + sig*[cos(acc/sig), sin(acc/sig)]
        values(1) = u(old_component)
        associate (unused => t)
        end associate
    end subroutine old_value

    subroutine fm3_wall_winding(c, profiles, winding, least_step, status)
        !! Per plane the winding number of u_h along the wall (the wall angle
        !! Theta_h followed once around, counter-clockwise, over 2 pi) and
        !! the smallest increment of the angle between neighbouring wall
        !! points (4 per cell side) over the mean increment: positive when
        !! the wall angle is monotone.
        type(fm3_case_t), intent(in) :: c
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: winding(:), least_step(:)
        integer, intent(out) :: status
        integer, allocatable :: wblock(:), rowptr(:), colidx(:)
        real(dp), allocatable :: wxi(:, :), ws(:, :), jets(:, :), angle(:)
        real(dp) :: rz(2), u(8), b(3), d
        integer :: k, i, n, first, last

        call spline_boundary_nodes(c%space, LAYOUT_WALL, 2, wblock, wxi, ws, status)
        if (status /= 0) return
        call spline_point_jet_matrix(c%space, wblock, wxi, 1, rowptr, colidx, jets, status)
        if (status /= 0) return
        n = size(wblock)
        allocate (angle(n))
        do k = 1, c%nt
            do i = 1, n
                first = rowptr(i)
                last = rowptr(i + 1) - 1
                call flux_map3d_point(c%fm, c%base, colidx(first:last), jets(1, first:last), &
                    jets(2:3, first:last), (k - 1)/real(c%nt, dp), profiles, rz, u, b)
                angle(i) = atan2(u(2), u(1))
            end do
            winding(k) = 0
            least_step(k) = huge(1.0_dp)
            do i = 1, n
                d = angle(mod(i, n) + 1) - angle(i)
                d = d - two_pi*nint(d/two_pi)
                winding(k) = winding(k) + d
                least_step(k) = min(least_step(k), d*n/two_pi)
            end do
            winding(k) = winding(k)/two_pi
        end do
    end subroutine fm3_wall_winding

    subroutine fm3_interface_jump(c, profiles, jump, status)
        !! Magnetic-pressure jump at C0 section block interfaces, scaled by
        !! RMS magnetic pressure (Gauss points of every cell side, plane 1
        !! weights and magnetic field values from all planes).
        type(fm3_case_t), intent(in) :: c
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: jump
        integer, intent(out) :: status
        type(spline_side_quadrature_t) :: sq
        integer, allocatable :: rowptr(:), colidx(:), opposite_ptr(:), opposite_idx(:)
        real(dp), allocatable :: jets(:, :), opposite_jets(:, :)
        real(dp) :: b(3, 2), acc(2), t, rz(2), u(8)
        integer :: q, k, first, last, opposite_first, opposite_last

        call spline_interface_quadrature(c%space, LAYOUT_INTERIOR, c%degree + 2, sq, status)
        if (status /= 0) return
        call spline_point_jet_matrix(c%space, sq%block, sq%xi, 1, rowptr, colidx, jets, status)
        if (status /= 0) return
        call spline_point_jet_matrix(c%space, sq%neighbour_block, sq%neighbour_xi, 1, &
            opposite_ptr, opposite_idx, opposite_jets, status)
        if (status /= 0) return
        acc = 0
        do q = 1, sq%point_count
            first = rowptr(q)
            last = rowptr(q + 1) - 1
            opposite_first = opposite_ptr(q)
            opposite_last = opposite_ptr(q + 1) - 1
            do k = 1, c%nt
                t = (k - 1)/real(c%nt, dp)
                call flux_map3d_point(c%fm, c%base, colidx(first:last), jets(1, first:last), &
                    jets(2:3, first:last), t, profiles, rz, u, b(:, 1))
                call flux_map3d_point(c%fm, c%base, opposite_idx(opposite_first:opposite_last), &
                    opposite_jets(1, opposite_first:opposite_last), &
                    opposite_jets(2:3, opposite_first:opposite_last), t, profiles, rz, u, b(:, 2))
                acc = acc + sq%weight(q)*[(sum(b(:, 1)**2) - sum(b(:, 2)**2))**2, &
                                          (sum(b(:, 1)**2) + sum(b(:, 2)**2))**2]/4
            end do
        end do
        jump = sqrt(acc(1)/acc(2))
    end subroutine fm3_interface_jump

    subroutine fm3_flux(c, profiles, flux, bmin, bmax)
        !! Toroidal flux through every plane and the range of B_phi over the
        !! quadrature points of the planes.
        type(fm3_case_t), intent(in) :: c
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: flux(:), bmin, bmax
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        real(dp) :: f(2, 2), rz(2), u(8), b(3)
        integer :: cell, q, p, k

        flux = 0
        bmin = huge(1.0_dp)
        bmax = -huge(1.0_dp)
        !$omp parallel do default(private) shared(c) reduction(+:flux) reduction(min:bmin) reduction(max:bmax)
        do cell = 1, size(c%quad%point_offset)
            call clebsch_node_gradients(c%derham, c%quad, cell, v0, d0, gn)
            do q = 1, c%quad%points_per_cell
                p = c%quad%point_offset(cell) + q - 1
                do k = 1, c%nt
                    f = matmul(c%base(:, d0, k), transpose(gn(:, :, q)))
                    call flux_map3d_point(c%fm, c%base, d0, v0(1, :, q), gn(:, :, q), &
                                          (k - 1)/real(c%nt, dp), profiles, rz, u, b)
                    flux(k) = flux(k) + c%quad%weight(p)*abs(f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1))*b(2)
                    bmin = min(bmin, b(2))
                    bmax = max(bmax, b(2))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine fm3_flux

    subroutine fm3_mesh_quality(c, quality)
        !! Per plane min/max of the Jacobian determinant of the section map
        !! over the quadrature points; negative: folded.
        type(fm3_case_t), intent(in) :: c
        real(dp), intent(out) :: quality(:)
        real(dp) :: minimum(c%nt), maximum(c%nt), det, widths(2)
        real(dp), allocatable :: table1(:, :, :), table2(:, :, :), jets(:, :, :)
        integer, allocatable :: geometry_nodes(:)
        integer :: cell, q, k, p, nquad

        p = c%degree
        nquad = c%quad%points_per_direction
        minimum = huge(1.0_dp)
        maximum = -huge(1.0_dp)
        !$omp parallel default(none) shared(c,p,nquad) &
        !$omp private(geometry_nodes,table1,table2,jets,widths,cell,k,q,det) &
        !$omp reduction(min:minimum) reduction(max:maximum)
        allocate (geometry_nodes((p + 1)**2), table1(p + 1, nquad, 0:2), &
                  table2(p + 1, nquad, 0:2), jets(2, 6, nquad**2))
        !$omp do
        do cell = 1, c%space%cell_count
            call spline_cell_broken(c%space, cell, geometry_nodes)
            geometry_nodes = c%space%broken_dof(c%space%broken_ptr(geometry_nodes))
            call spline_cell_tables(c%space, c%quad, cell, table1, table2, widths)
            do k = 1, c%nt
                call spline_cell_xi_jets(table1, table2, c%base(:, geometry_nodes, k), jets)
                do q = 1, nquad**2
                    det = jets(1, 2, q)*jets(2, 3, q) - jets(1, 3, q)*jets(2, 2, q)
                    minimum(k) = min(minimum(k), det)
                    maximum(k) = max(maximum(k), det)
                end do
            end do
        end do
        !$omp end do
        deallocate (geometry_nodes, table1, table2, jets)
        !$omp end parallel
        quality = minimum/maximum
    end subroutine fm3_mesh_quality

    subroutine fm3_dump_hessian(c, profiles)
        !! With the environment variable MILCA_FM3_HESSIAN=prefix: the Hessian
        !! of W at c%fm%z (without the gauge constraints) as CSR, the state,
        !! the gradient, the fixed slots, which section functions are wall
        !! functions, W, the per-plane L2 mass matrix of the section functions
        !! (fm3_section_mass) and the node coordinates, to
        !! prefix_<cells>_<degree>_<planes>.bin (stream; 4-byte integers,
        !! 8-byte reals). Read by benchmarks/flux-map-spike/hessian_spectrum.py.
        type(fm3_case_t), intent(in) :: c
        procedure(flux_map_profiles) :: profiles
        type(csr_t) :: h, upper
        real(dp), allocatable :: g(:)
        real(dp) :: w
        character(512) :: prefix, path
        integer :: unit, length

        call get_environment_variable('MILCA_FM3_HESSIAN', prefix, length)
        if (length == 0) return
        call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, w, g, upper)
        call flux_map3d_expand_symmetric(upper, h)
        write (path, '(a,a,i0,a,i0,a,i0,a)') trim(prefix), '_', c%cells, '_', c%degree, '_', &
            c%nt, '.bin'
        open (newunit=unit, file=trim(path), access='stream', form='unformatted', &
              status='replace')
        write (unit) int(h%nrow, int32), int(size(h%val), int32), int(c%n0, int32), &
            int(c%nt, int32), int(size(c%fm%fixed), int32), 0_int32
        write (unit) int(h%row_ptr, int32), int(h%col, int32), h%val, c%fm%z, g, &
            int(c%fm%fixed, int32), merge(1_int32, 0_int32, c%fm%is_wall), w
        write (unit) fm3_section_mass(c), c%base
        close (unit)
        print '(2a)', '  Hessian written to ', trim(path)
    end subroutine fm3_dump_hessian

    function fm3_section_mass(c) result(m)
        !! Per plane k the mass matrix int phi_i phi_j R dA (2 pi/(nt nfp)) of
        !! the section functions in the geometry of plane k: with it
        !! H v = mu M v has eigenvalues mu that compare between meshes (an L2
        !! metric of the perturbation of u; the wall angle slot taken as
        !! |delta u|).
        type(fm3_case_t), intent(in) :: c
        real(dp) :: m(c%n0, c%n0, c%nt)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        real(dp) :: f(2, 2), r, dv
        integer :: cell, q, p, k

        m = 0
        do cell = 1, size(c%quad%point_offset)
            call clebsch_node_gradients(c%derham, c%quad, cell, v0, d0, gn)
            do q = 1, c%quad%points_per_cell
                p = c%quad%point_offset(cell) + q - 1
                do k = 1, c%nt
                    f = matmul(c%base(:, d0, k), transpose(gn(:, :, q)))
                    r = dot_product(c%base(1, d0, k), v0(1, :, q))
                    dv = c%quad%weight(p)*abs(f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1))*r*two_pi/ &
                         (c%nt*clebsch_nfp)
                    m(d0, d0, k) = m(d0, d0, k) + dv*spread(v0(1, :, q), 2, size(d0))* &
                                   spread(v0(1, :, q), 1, size(d0))
                end do
            end do
        end do
    end function fm3_section_mass

end module flux_map_sections
