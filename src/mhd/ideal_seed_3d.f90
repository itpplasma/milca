module milca_ideal_seed_3d
    ! Fixed-boundary seed only: no exact interior equilibrium is used.
    ! Boundary convention: R,Z = sum C_mn cos(m theta-n*nfp*phi)
    !                                  + S_mn sin(m theta-n*nfp*phi).
    ! The reference disk is (Rref-2,Zref)=(u,v), s=u^2+v^2. Profiles use
    ! normalized TOROIDAL flux, unlike ideal_seed's poloidal-flux profiles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_map, only: ideal_map_t, ideal_frozen_t, ideal_map_init, &
        ideal_map_freeze, ideal_map_energy, ideal_map_rdof, IDEAL_MAP_OK, &
        IDEAL_MAP_BAD_INPUT, IDEAL_MAP_FOLDED, ideal_map_angle_enable, &
        ideal_map_angle_value
    use milca_ideal_map_leaves, only: ideal_reference_form
    use milca_ideal_field_angle_leaves, only: field_angle_form
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_triangle_basis, only: TRIANGLE_MAX_DEGREE
    use milca_grad_shafranov, only: gs_state_t, initialize_gs
    use milca_winslow_seed, only: winslow_seed, WINSLOW_OK, WINSLOW_NOT_SIMPLE
    use milca_ideal_seed_3d_leaves, only: seed_complex_product, seed_boundary_mode, &
        seed_horner, seed_toroidal_form, seed_disk_coordinates, seed_disk_boundary, &
        seed_primitive_coefficient, seed_primitive_finish, seed_hermite_profile, &
        seed_disc_coordinates, seed_disc_pullback
    implicit none
    private
    public :: toroidal_profile_t, toroidal_seed_t, ideal_seed_3d_init
    public :: hermite_profile_evaluate
    public :: ideal_seed_disk_mesh, ideal_seed_boundary_coefficients
    public :: ideal_seed_3d_angle_density
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: toroidal_profile_t
        ! Ascending monomial coefficients in s, pressure in solver units.
        real(dp), allocatable :: iota(:), pressure(:)
        ! Piecewise cubic Hermite profiles on knots 0 = s_1 < ... < s_K = 1
        ! with values and s-slopes; each replaces the monomials when set by
        ! set_hermite_iota or set_hermite_pressure.
        real(dp), allocatable :: iota_knots(:), iota_values(:), iota_slopes(:)
        real(dp), allocatable :: iota_primitive(:)
        real(dp), allocatable :: pressure_knots(:), pressure_values(:), pressure_slopes(:)
    contains
        procedure :: evaluate => profile_evaluate
        procedure :: set_hermite_iota, set_hermite_pressure, shift_iota, valid => profile_valid
    end type
    type, extends(ideal_frozen_t) :: toroidal_seed_t
        type(toroidal_profile_t) :: profile
        real(dp) :: flux = 0
        ! Reference disc: centre and half-widths in q1 and q2 (an ellipse).
        ! The free boundary fits it to the plasma so that its shell map stays
        ! close to the physical one.
        real(dp) :: centre(2) = [2, 0], axes(2) = [1, 1]
    contains
        procedure :: evaluate => frozen_evaluate
        procedure :: flux_label => seed_flux_label
    end type
contains
    pure subroutine profile_evaluate(profile, s, iota, pressure, integral)
        class(toroidal_profile_t), intent(in) :: profile
        real(dp), intent(in) :: s
        real(dp), intent(out) :: iota, pressure, integral
        real(dp) :: next, coefficient
        integer :: k
        iota = 0; pressure = 0; integral = 0
        if (allocated(profile%pressure_values)) then
            call hermite_profile_evaluate(profile%pressure_knots, profile%pressure_values, &
                profile%pressure_slopes, s, pressure, next)
        else
            do k = size(profile%pressure), 1, -1
                call seed_horner(s, profile%pressure(k), pressure, next)
                pressure = next
            end do
        end if
        if (allocated(profile%iota_values)) then
            call hermite_profile_evaluate(profile%iota_knots, profile%iota_values, &
                profile%iota_slopes, s, iota, integral, profile%iota_primitive)
            return
        end if
        do k = size(profile%iota), 1, -1
            call seed_horner(s, profile%iota(k), iota, next)
            iota = next
            call seed_primitive_coefficient(profile%iota(k), real(k, dp), coefficient)
            call seed_horner(s, coefficient, integral, next)
            integral = next
        end do
        call seed_primitive_finish(s, integral, next)
        integral = next
    end subroutine

    pure subroutine hermite_profile_evaluate(knots, values, slopes, s, value, primitive, &
            cumulative)
        ! Value and integral from 0 of the piecewise cubic Hermite profile at s;
        ! cumulative(k) is the integral up to knots(k) when known. Outside
        ! [knots(1), knots(K)] the end pieces continue as cubics.
        real(dp), intent(in) :: knots(:), values(:), slopes(:), s
        real(dp), intent(out) :: value, primitive
        real(dp), intent(in), optional :: cumulative(:)
        real(dp) :: h, piece, unused
        integer :: low, high, mid, k
        low = 1; high = size(knots)
        do while (high - low > 1)
            mid = (low + high)/2
            if (s < knots(mid)) then
                high = mid
            else
                low = mid
            end if
        end do
        h = knots(low + 1) - knots(low)
        call seed_hermite_profile((s - knots(low))/h, h, values(low), slopes(low), &
            values(low + 1), slopes(low + 1), value, piece)
        if (present(cumulative)) then
            primitive = cumulative(low) + piece
            return
        end if
        primitive = piece
        do k = 1, low - 1
            h = knots(k + 1) - knots(k)
            call seed_hermite_profile(1.0_dp, h, values(k), slopes(k), values(k + 1), &
                slopes(k + 1), unused, piece)
            primitive = primitive + piece
        end do
    end subroutine

    subroutine set_hermite_iota(profile, knots, values, slopes)
        class(toroidal_profile_t), intent(inout) :: profile
        real(dp), intent(in) :: knots(:), values(:), slopes(:)
        profile%iota_knots = knots
        profile%iota_values = values; profile%iota_slopes = slopes
        call hermite_iota_primitive(profile)
    end subroutine

    subroutine set_hermite_pressure(profile, knots, values, slopes)
        class(toroidal_profile_t), intent(inout) :: profile
        real(dp), intent(in) :: knots(:), values(:), slopes(:)
        profile%pressure_knots = knots
        profile%pressure_values = values; profile%pressure_slopes = slopes
    end subroutine

    subroutine hermite_iota_primitive(profile)
        class(toroidal_profile_t), intent(inout) :: profile
        real(dp) :: unused
        integer :: k
        if (allocated(profile%iota_primitive)) deallocate (profile%iota_primitive)
        allocate (profile%iota_primitive(size(profile%iota_knots)))
        profile%iota_primitive(1) = 0
        do k = 2, size(profile%iota_knots)
            call hermite_profile_evaluate(profile%iota_knots(:k), profile%iota_values(:k), &
                profile%iota_slopes(:k), profile%iota_knots(k), unused, &
                profile%iota_primitive(k), profile%iota_primitive(:k))
        end do
    end subroutine

    subroutine shift_iota(profile, delta)
        ! iota -> iota + delta in either representation
        class(toroidal_profile_t), intent(inout) :: profile
        real(dp), intent(in) :: delta
        if (allocated(profile%iota_values)) then
            profile%iota_values = profile%iota_values + delta
            call hermite_iota_primitive(profile)
        else
            profile%iota(1) = profile%iota(1) + delta
        end if
    end subroutine

    logical function profile_valid(profile)
        ! one representation per quantity, finite, and Hermite knots strictly
        ! increasing from 0 to 1
        class(toroidal_profile_t), intent(in) :: profile
        profile_valid = .false.
        if (allocated(profile%iota_values)) then
            if (.not. hermite_ok(profile%iota_knots, profile%iota_values, &
                profile%iota_slopes)) return
        else
            if (.not. allocated(profile%iota)) return
            if (size(profile%iota) < 1 .or. .not. all(ieee_is_finite(profile%iota))) return
        end if
        if (allocated(profile%pressure_values)) then
            if (.not. hermite_ok(profile%pressure_knots, profile%pressure_values, &
                profile%pressure_slopes)) return
        else
            if (.not. allocated(profile%pressure)) return
            if (size(profile%pressure) < 1 .or. .not. all(ieee_is_finite(profile%pressure))) &
                return
        end if
        profile_valid = .true.
    contains
        logical function hermite_ok(knots, values, slopes)
            real(dp), allocatable, intent(in) :: knots(:), values(:), slopes(:)
            integer :: n
            hermite_ok = .false.
            if (.not. allocated(knots) .or. .not. allocated(slopes)) return
            n = size(knots)
            if (n < 2 .or. size(values) /= n .or. size(slopes) /= n) return
            if (.not. all(ieee_is_finite(knots))) return
            if (abs(knots(1)) > 0 .or. abs(knots(n) - 1) > 0) return
            if (any(knots(2:) <= knots(:n - 1))) return
            hermite_ok = all(ieee_is_finite(values)) .and. all(ieee_is_finite(slopes))
        end function
    end function

    subroutine frozen_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(toroidal_seed_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: u, v, s, iota, p, integral, edge, qi, pi_unused, gu, gv, bu
        call seed_disc_coordinates(point(1), point(2), field%centre(1), field%centre(2), &
            field%axes(1), field%axes(2), u, v, s)
        call field%profile%evaluate(s, iota, p, integral)
        call field%profile%evaluate(1.0_dp, qi, pi_unused, edge)
        call seed_toroidal_form(u, v, field%flux, pi, integral, edge, iota, p, &
            chi, gu, gv, bu, pressure)
        call seed_disc_pullback(gu, gv, bu, field%axes(1), field%axes(2), &
            chi_gradient(1), chi_gradient(2), b3)
    end subroutine

    subroutine seed_flux_label(field, point, s, status)
        ! the reference disk radius squared is the normalised toroidal flux
        class(toroidal_seed_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: s
        integer, intent(out) :: status
        real(dp) :: u, v
        call seed_disk_coordinates(point(1), point(2), u, v, s)
        status = 0
    end subroutine

    subroutine ideal_seed_disk_mesh(mesh, boundary_nodes, degree, status)
        type(gs_state_t), intent(out) :: mesh
        integer, intent(in) :: boundary_nodes, degree
        integer, intent(out) :: status
        real(dp) :: points(2, boundary_nodes), h
        integer :: j
        status = IDEAL_MAP_BAD_INPUT
        if (degree < 2 .or. degree > TRIANGLE_MAX_DEGREE) return
        ! Fewer, larger cells are the point of the high orders; the periodic
        ! boundary spline needs four nodes.
        if (boundary_nodes < merge(8, 4, degree <= 3)) return
        do j = 1, boundary_nodes
            call disk_boundary(real(j - 1, dp)/boundary_nodes, points(:, j))
        end do
        h = 2*pi/boundary_nodes
        call initialize_gs(mesh, points, h*h/2, status, degree, curve=disk_boundary, &
            geometry_only=.true.)
    end subroutine

    subroutine disk_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call seed_disk_boundary(t, pi, point(1), point(2))
    end subroutine

    subroutine ideal_seed_boundary_coefficients(point, coefficients, values, axis, &
            centre, axes)
        ! Harmonic extension: rho^m exp(i m theta)=(u+i v)^m is regular
        ! at the axis. Mode sums and their storage are orchestration; the
        ! complex recurrence and real Fourier conversion are generated.
        ! An axis (the m = 0 layout of the coefficients) replaces the
        ! boundary's m = 0 curve at the centre: m = 0 + (1 - s)(axis - m = 0).
        real(dp), intent(in) :: point(2), coefficients(:, :, 0:)
        real(dp), intent(out) :: values(:, :)
        real(dp), intent(in), optional :: axis(:, :)
        ! reference disc (default: centre (2, 0), half-widths 1)
        real(dp), intent(in), optional :: centre(2), axes(2)
        real(dp) :: u, v, s, re, im, nr, ni, term(4), c(4)
        integer :: m, n, ntor, mode
        ntor = (size(coefficients, 2) - 1)/2
        values = 0; re = 1; im = 0
        if (present(centre) .and. present(axes)) then
            call seed_disc_coordinates(point(1), point(2), centre(1), centre(2), axes(1), &
                axes(2), u, v, s)
        else
            call seed_disk_coordinates(point(1), point(2), u, v, s)
        end if
        do m = 0, ubound(coefficients, 3)
            do n = -ntor, ntor
                c = coefficients(:, n + ntor + 1, m)
                if (m == 0 .and. present(axis)) c = c + (1 - s)*(axis(:, n + ntor + 1) - c)
                call seed_boundary_mode(re, im, c(1), c(2), c(3), c(4), &
                    term(1), term(2), term(3), term(4))
                mode = 2*abs(n)
                if (n == 0) then
                    values(1:2, 1) = values(1:2, 1) + term([1, 3])
                else
                    values(1:2, mode) = values(1:2, mode) + term([1, 3])
                    values(1:2, mode + 1) = values(1:2, mode + 1) &
                        + sign(1.0_dp, real(n, dp))*term([2, 4])
                end if
            end do
            call seed_complex_product(re, im, u, v, nr, ni)
            re = nr; im = ni
        end do
    end subroutine

    subroutine ideal_seed_3d_init(map, seed, mesh, coefficients, nfp, flux, &
            profile, status, nphi, angle_mpol, angle_ntor, axis, winslow)
        type(ideal_map_t), intent(out) :: map
        type(toroidal_seed_t), intent(out) :: seed
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: coefficients(:, :, 0:), flux
        integer, intent(in) :: nfp
        type(toroidal_profile_t), intent(in) :: profile
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi
        ! Field-line angle of the frozen field (docs/ideal-seed-3d.md): the
        ! boundary may then be given in any poloidal angle. angle_mpol >= 0
        ! enables it with poloidal modes <= angle_mpol and toroidal modes
        ! <= angle_ntor (default: the boundary ntor). Default: off.
        integer, intent(in), optional :: angle_mpol, angle_ntor
        ! Initial magnetic axis in the m = 0 layout of the coefficients
        ! (default: the boundary's m = 0 curve).
        real(dp), intent(in), optional :: axis(:, :)
        ! Interior of every toroidal plane from the Winslow map of its
        ! cross-section (milca_winslow_seed) instead of the harmonic extension.
        logical, intent(in), optional :: winslow
        type(cylindrical_chart_t) :: chart
        real(dp) :: iota, pressure, integral, previous, energy, minimum_radius
        integer :: a, m, j, ntor, angles
        status = IDEAL_MAP_BAD_INPUT
        if (nfp < 1 .or. .not. flux > 0) return
        if (.not. ieee_is_finite(flux)) return
        if (size(coefficients, 1) /= 4) return
        if (present(axis)) then
            if (any(shape(axis) /= [4, size(coefficients, 2)])) return
            if (.not. all(ieee_is_finite(axis))) return
        end if
        if (mod(size(coefficients, 2), 2) /= 1) return
        if (.not. all(ieee_is_finite(coefficients))) return
        if (.not. profile%valid()) return
        previous = 0
        do a = 0, 128
            call profile%evaluate(a/128.0_dp, iota, pressure, integral)
            if (.not. all(ieee_is_finite([iota, pressure, integral]))) return
            if (abs(iota) <= tiny(iota)) return
            if (a > 0) then
                if ((iota > 0) .neqv. (previous > 0)) return
            end if
            previous = iota
        end do
        ntor = (size(coefficients, 2) - 1)/2
        angles = max(1, 4*ntor + 1)
        if (present(nphi)) angles = nphi
        call ideal_map_init(map, mesh, chart, nfp, ntor, status, nphi=angles)
        if (status /= IDEAL_MAP_OK) return
        seed%profile = profile; seed%flux = flux
        ! A fixed physical toroidal coordinate admits both R/Z increments.
        ! Reference-disk radial directions can become tangent to the seeded
        ! surfaces (already for a clockwise circular boundary). The frozen
        ! b3 is positive, so the phi-fixed representative is admissible.
        map%position_gauge = .true.
        call ideal_map_freeze(map, seed, status)
        if (status /= IDEAL_MAP_OK) return
        minimum_radius=huge(1.0_dp)
        !$omp parallel do schedule(static) private(a,m)
        do a = 1, map%nodes
            call ideal_seed_boundary_coefficients(map%reference(:, a), &
                coefficients, map%x(:, :, a), axis)
            if (.not. map%boundary(a)) cycle
            ! A nonaxisymmetric fixed wall is parametrized by phi-hat. A
            ! boundary omega would move its points off that prescribed wall.
            do m = 1, map%nmode
                map%fixed(ideal_map_rdof(map, 2, m, a)) = .true.
            end do
        end do
        !$omp end parallel do
        if (present(winslow)) then
            if (winslow) then
                call winslow_seed(map, m)
                if (m == WINSLOW_NOT_SIMPLE) then
                    status = IDEAL_MAP_BAD_INPUT
                    return
                else if (m /= WINSLOW_OK) then
                    status = IDEAL_MAP_FOLDED
                    return
                end if
            end if
        end if
        do a = 1, map%nodes
            do j = 1, map%nphi
                minimum_radius = min(minimum_radius, dot_product(map%x(1, :, a), map%trig(:, j)))
            end do
        end do
        if (.not. minimum_radius>0) then
            status=IDEAL_MAP_BAD_INPUT
            return
        end if
        ! Sampled positive-Jacobian admission, not a whole-cell certificate.
        call ideal_map_energy(map, map%x, energy, status)
        if (status == IDEAL_MAP_FOLDED) then
            ! The supplied poloidal angle may be clockwise. Orientation is
            ! fixed by the initialized physical map, not the identity disk
            ! used temporarily by ideal_map_init. Mixed signs still fail.
            map%sigma = -map%sigma
            call ideal_map_energy(map, map%x, energy, status)
        end if
        if (status /= IDEAL_MAP_OK .or. .not. present(angle_mpol)) return
        if (angle_mpol < 0) return
        m = ntor
        if (present(angle_ntor)) m = angle_ntor
        call ideal_map_angle_enable(map, angle_mpol, m, status)
        if (status /= IDEAL_MAP_OK) return
        call seed_angle_tables(map, flux)
    end subroutine

    subroutine seed_angle_tables(map, flux)
        ! Reference harmonics (u + i v)^m and the frozen-form increments of a
        ! unit eta_t and eta_p, pulled back to master coordinates like the
        ! frozen density itself (ideal_map_freeze).
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: flux
        real(dp) :: y(2), dg(2, 2), u, v, s, c(3, 2), b(3)
        integer :: e, g, k
        !$omp parallel do schedule(static) private(e, g, k, y, dg, u, v, s, c, b)
        do e = 1, map%cells
            do g = 1, map%nq
                y = matmul(map%reference(:, map%conn(:, e)), map%shape(:, g))
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(map%dshape(:, :, g)))
                call seed_disk_coordinates(y(1), y(2), u, v, s)
                call harmonics(u, v, map%angle_power(:, :, g, e))
                call field_angle_form(u, v, 1.0_dp, 0.0_dp, flux, pi, &
                    c(1, 1), c(2, 1), c(3, 1))
                call field_angle_form(u, v, 0.0_dp, 1.0_dp, flux, pi, &
                    c(1, 2), c(2, 2), c(3, 2))
                do k = 1, 2
                    call ideal_reference_form(dg(1, 1), dg(2, 1), dg(1, 2), dg(2, 2), &
                        c(1, k), c(2, k), c(3, k), b(1), b(2), b(3))
                    map%angle_form(:, k, g, e) = b
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine

    pure subroutine harmonics(u, v, power)
        ! power(:, m) = (Re, Im) (u + i v)^m, m = 0..ubound.
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: power(:, 0:)
        integer :: m
        power(:, 0) = [1.0_dp, 0.0_dp]
        do m = 1, ubound(power, 2)
            call seed_complex_product(power(1, m - 1), power(2, m - 1), u, v, &
                power(1, m), power(2, m))
        end do
    end subroutine

    subroutine ideal_seed_3d_angle_density(map, seed, point, phi, dcg, db3)
        ! Field-angle increment of the frozen data at a reference point and
        ! toroidal angle phi-hat: chi gradient += dcg, b3 += db3 (zero without
        ! field angle). Add it to toroidal_seed_t%evaluate when sampling B.
        type(ideal_map_t), intent(in) :: map
        type(toroidal_seed_t), intent(in) :: seed
        real(dp), intent(in) :: point(2), phi
        real(dp), intent(out) :: dcg(2), db3
        real(dp) :: u, v, s, eta, eta_t, eta_p
        real(dp), allocatable :: power(:, :)
        dcg = 0; db3 = 0
        if (map%nangle == 0) return
        allocate (power(2, 0:map%angle_mpol))
        call seed_disk_coordinates(point(1), point(2), u, v, s)
        call harmonics(u, v, power)
        call ideal_map_angle_value(map, map%angle, power, phi, eta, eta_t, eta_p)
        call field_angle_form(u, v, eta_t, eta_p, seed%flux, pi, dcg(1), dcg(2), db3)
    end subroutine
end module
