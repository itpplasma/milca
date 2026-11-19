module milca_import_gvec
    ! GVEC parameter files (INI, TOML, YAML) to the intermediate form of
    ! milca_equilibrium_import. GVEC maps (X1, X2, zeta) by its hmap; X1, X2
    ! = sum c cos(m theta - n nfp zeta) + s sin(...), theta counterclockwise in
    ! the (X1, X2) plane (GVEC documentation, "Coordinate conventions", and
    ! src/hmap of gvec-group/gvec, MIT licence). Under phi = s zeta with
    ! s = +-1, iota, the toroidal flux and the toroidal current take the
    ! factor s, as GVEC's own VMEC conversion flips them (python/gvec/util).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_import_text, only: kv_store_t, read_gvec_ini, read_toml, read_yaml, lower, mode_key
    use milca_equilibrium_import, only: imported_equilibrium_t, import_profile_t, import_note, &
        import_refuse, import_fail, IMPORT_OK, boundary_fit, power_series_profile, &
        spline_profile, hermite_profile, hermite_value, IN_MMAX, IN_NMAX
    use milca_import_vmec, only: import_wout
    use fortio_hdf5_reader, only: hdf5_file_t
    use fortio_status, only: fortio_status_t
    implicit none
    private
    public :: import_gvec, gvec_surface_t, cartesian_refit
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: gvec_surface_t
        ! q1, q2 of (theta, zeta): cos/sin coefficients of (m theta - n nfp zeta)
        integer :: nfp = 1
        real(dp), allocatable :: c1(:, :), s1(:, :), c2(:, :), s2(:, :)
        integer, allocatable :: mode_m(:), mode_n(:)
        real(dp), allocatable :: mode_c(:, :)
        ! frame: 0 cylindrical, 20 Frenet, 21 G-frame
        integer :: frame = 0
        ! Frenet axis R0, Z0 = sum rc cos(n nfp zeta) + rs sin(n nfp zeta)
        real(dp), allocatable :: rc(:), rs(:), zc(:), zs(:)
        logical :: omnig = .false.
        ! G-frame: one field period of the axis and its N, B vectors in the
        ! co-rotating frame, as cos/sin coefficients of k nfp zeta
        integer :: winding = 1
        real(dp), allocatable :: hc(:, :, :), hs(:, :, :)
    end type

contains

    subroutine import_gvec(path, eq)
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        type(kv_store_t) :: p
        type(gvec_surface_t) :: surface
        character(len=:), allocatable :: dir, ext
        logical :: ok
        integer :: hmap
        real(dp) :: s_zeta
        eq%format = 'GVEC'
        eq%source = path
        dir = directory(path)
        ext = lower(path(index(path, '.', back=.true.) + 1:))
        select case (ext)
        case ('toml')
            call read_toml(path, p, ok)
        case ('yaml', 'yml')
            call read_yaml(path, p, ok)
        case default
            call read_gvec_ini(path, p, ok)
        end select
        if (.not. ok) then
            call import_fail(eq, 'cannot read '//trim(path))
            return
        end if
        if (p%int('whichinitequilibrium', 0) == 1) then
            ! boundary and profiles of a VMEC wout file (GVEC src/vmec/vmec.F90)
            call import_wout(dir//p%text('vmecwoutfile', ''), eq)
            eq%format = 'GVEC (VMEC wout)'
            eq%source = path
            if (p%has('iota_type') .or. p%has('pres_type') .or. p%has('x1_b_cos(0,0)')) &
                call import_note(eq, 'profile or boundary parameters beside the wout file are ignored')
            return
        end if
        if (p%int('whichinitequilibrium', 0) /= 0) then
            call import_refuse(eq, 'GVEC whichInitEquilibrium other than 0 or 1')
            return
        end if
        hmap = p%int('which_hmap', 1)
        surface%nfp = p%int('nfp', 1)
        select case (hmap)
        case (1, 20, 21)
        case (3)
            call import_refuse(eq, 'GVEC hmap 3 is a straight periodic cylinder; '// &
                'nested_equilibrium solves toroidal equilibria')
            return
        case (10)
            call import_refuse(eq, 'GVEC hmap 10 is a knotted torus, not single-valued over the '// &
                'cylindrical angle')
            return
        case default
            call import_refuse(eq, 'GVEC hmap other than 1, 20 or 21')
            return
        end select
        if (p%int('getboundaryfromfile', 0) == 1) then
            call boundary_from_file(dir//p%text('boundary_filename', ''), surface, eq)
        else
            call boundary_from_parameters(p, surface)
        end if
        if (eq%status /= IMPORT_OK) return
        call compact_modes(surface)
        select case (hmap)
        case (1)
            ! x = R cos(zeta), y = -R sin(zeta): phi = -zeta
            s_zeta = -1
            eq%nfp = surface%nfp
            call cylindrical_modes(surface, eq)
        case (20)
            surface%frame = 20
            call frenet_axis(p, surface, eq)
            if (eq%status /= IMPORT_OK) return
            call cartesian_refit(surface, eq, s_zeta)
        case (21)
            surface%frame = 21
            call gframe_axis(dir//p%text('hmap_ncfile', ''), surface, eq)
            if (eq%status /= IMPORT_OK) return
            call cartesian_refit(surface, eq, s_zeta)
        end select
        if (eq%status /= IMPORT_OK) return
        eq%phiedge = s_zeta*p%real('phiedge', 1.0_dp)
        call gvec_profile(p, 'pres', eq%pressure, eq)
        if (p%has('i_tor_type')) then
            call gvec_profile(p, 'i_tor', eq%current, eq)
            eq%current%coefficients = s_zeta*eq%current%coefficients
            if (allocated(eq%current%values)) then
                eq%current%values = s_zeta*eq%current%values
                eq%current%slopes = s_zeta*eq%current%slopes
            end if
            eq%prescribe_current = .true.
            call import_note(eq, 'I_tor prescribed; GVEC''s iota is the start of the current matching')
        end if
        if (p%has('iota_type')) then
            call gvec_profile(p, 'iota', eq%iota, eq)
            if (allocated(eq%iota%coefficients)) eq%iota%coefficients = s_zeta*eq%iota%coefficients
            if (allocated(eq%iota%values)) then
                eq%iota%values = s_zeta*eq%iota%values
                eq%iota%slopes = s_zeta*eq%iota%slopes
            end if
        else if (.not. eq%prescribe_current) then
            call import_fail(eq, 'GVEC input without an iota or I_tor profile')
        end if
    end subroutine

    function directory(path) result(dir)
        character(*), intent(in) :: path
        character(len=:), allocatable :: dir
        integer :: slash
        slash = index(path, '/', back=.true.)
        dir = path(:slash)
    end function

    subroutine boundary_from_parameters(p, surface)
        ! X1_b_*, X2_b_* within X?_mn_max and X?_sin_cos, plus X?pert_b_* when
        ! boundary_perturb: GVEC's boundary is their sum (mhd3d.F90)
        type(kv_store_t), intent(in) :: p
        type(gvec_surface_t), intent(inout) :: surface
        integer :: mn1(2), mn2(2)
        logical :: perturb
        perturb = p%logical('boundary_perturb', .false.)
        call mode_range(p, 'x1', mn1)
        call mode_range(p, 'x2', mn2)
        allocate (surface%c1(-IN_NMAX:IN_NMAX, 0:IN_MMAX), surface%s1(-IN_NMAX:IN_NMAX, 0:IN_MMAX), &
            surface%c2(-IN_NMAX:IN_NMAX, 0:IN_MMAX), surface%s2(-IN_NMAX:IN_NMAX, 0:IN_MMAX))
        surface%c1 = 0; surface%s1 = 0; surface%c2 = 0; surface%s2 = 0
        call add_modes('x1_b_cos', 'x1', mn1, surface%c1, 'cos')
        call add_modes('x1_b_sin', 'x1', mn1, surface%s1, 'sin')
        call add_modes('x2_b_cos', 'x2', mn2, surface%c2, 'cos')
        call add_modes('x2_b_sin', 'x2', mn2, surface%s2, 'sin')
        if (perturb) then
            call add_modes('x1pert_b_cos', 'x1', mn1, surface%c1, 'cos')
            call add_modes('x1pert_b_sin', 'x1', mn1, surface%s1, 'sin')
            call add_modes('x2pert_b_cos', 'x2', mn2, surface%c2, 'cos')
            call add_modes('x2pert_b_sin', 'x2', mn2, surface%s2, 'sin')
        end if
    contains
        subroutine add_modes(name, variable, limit, c, kind)
            character(*), intent(in) :: name, variable, kind
            integer, intent(in) :: limit(2)
            real(dp), intent(inout) :: c(-IN_NMAX:, 0:)
            integer, allocatable :: m(:), n(:)
            real(dp), allocatable :: v(:)
            character(len=:), allocatable :: basis
            integer :: k
            basis = lower(p%text(variable//'_sin_cos', '_sincos_'))
            if (index(basis, kind) == 0) return
            call p%modes(name, m, n, v)
            do k = 1, size(m)
                if (m(k) < 0 .or. m(k) > limit(1) .or. abs(n(k)) > limit(2)) cycle
                ! GVEC's basis holds m = 0 with n >= 0 only
                if (m(k) == 0 .and. n(k) < 0) cycle
                c(n(k), m(k)) = c(n(k), m(k)) + v(k)
            end do
        end subroutine
    end subroutine

    subroutine mode_range(p, variable, limit)
        type(kv_store_t), intent(in) :: p
        character(*), intent(in) :: variable
        integer, intent(out) :: limit(2)
        real(dp), allocatable :: values(:)
        logical :: ok
        limit = [IN_MMAX, IN_NMAX]
        call p%reals(variable//'_mn_max', values, ok)
        if (ok .and. size(values) == 2) limit = min(nint(values), [IN_MMAX, IN_NMAX])
    end subroutine

    subroutine boundary_from_file(path, surface, eq)
        ! X(theta_i, zeta_j), Y on one field period, projected on GVEC's
        ! modes m <= m_max, |n| <= n_max (src/grid/boundary_from_file)
        character(*), intent(in) :: path
        type(gvec_surface_t), intent(inout) :: surface
        type(imported_equilibrium_t), intent(inout) :: eq
        type(hdf5_file_t) :: h
        type(fortio_status_t) :: status
        real(dp), allocatable :: x(:, :), y(:, :), theta(:), zeta(:)
        integer :: ntheta, nzeta, mmax, nmax, lasym, nfp, i
        call h%open(path, status)
        if (.not. status%ok()) then
            call import_fail(eq, 'cannot read the boundary file '//path)
            return
        end if
        call read_integer(h, 'NFP', nfp, status)
        call read_integer(h, 'boundary/ntheta', ntheta, status)
        call read_integer(h, 'boundary/nzeta', nzeta, status)
        call read_integer(h, 'boundary/m_max', mmax, status)
        call read_integer(h, 'boundary/n_max', nmax, status)
        call read_integer(h, 'boundary/lasym', lasym, status)
        call read_grid(h, 'boundary/X(::)', ntheta, nzeta, x, status)
        call read_grid(h, 'boundary/Y(::)', ntheta, nzeta, y, status)
        if (.not. status%ok()) then
            call import_fail(eq, 'incomplete boundary file '//path)
            return
        end if
        surface%nfp = nfp
        call h%read_r64_1('boundary/theta(:)', theta, status)
        if (.not. status%ok()) theta = [(2*pi*(i - 1)/ntheta, i=1, ntheta)]
        call status%clear()
        call h%read_r64_1('boundary/zeta(:)', zeta, status)
        if (.not. status%ok()) zeta = [(2*pi*(i - 1)/(nfp*nzeta), i=1, nzeta)]
        allocate (surface%c1(-IN_NMAX:IN_NMAX, 0:IN_MMAX), surface%s1(-IN_NMAX:IN_NMAX, 0:IN_MMAX), &
            surface%c2(-IN_NMAX:IN_NMAX, 0:IN_MMAX), surface%s2(-IN_NMAX:IN_NMAX, 0:IN_MMAX))
        call project(x, theta, zeta, nfp, mmax, nmax, surface%c1, surface%s1)
        call project(y, theta, zeta, nfp, mmax, nmax, surface%c2, surface%s2)
        if (lasym == 0) then
            surface%s1 = 0; surface%c2 = 0
        end if
        call h%close(status)
    end subroutine

    subroutine read_integer(h, name, value, status)
        ! GVEC writes 64-bit integers
        type(hdf5_file_t), intent(inout) :: h
        character(*), intent(in) :: name
        integer, intent(out) :: value
        type(fortio_status_t), intent(inout) :: status
        real(dp), allocatable :: values(:)
        value = 0
        if (.not. status%ok()) return
        call h%read_numeric_r64_flat(name, values, status)
        if (status%ok() .and. size(values) == 1) value = nint(values(1))
    end subroutine

    subroutine read_grid(h, name, n1, n2, values, status)
        ! a 2-D dataset written as (n1, n2) by GVEC's Python tools
        type(hdf5_file_t), intent(inout) :: h
        character(*), intent(in) :: name
        integer, intent(in) :: n1, n2
        real(dp), allocatable, intent(out) :: values(:, :)
        type(fortio_status_t), intent(inout) :: status
        real(dp), allocatable :: raw(:, :)
        call h%read_r64_2(name, raw, status)
        if (.not. status%ok()) return
        ! the file holds (n1, n2) in C order, which reads back reversed
        if (size(raw, 1) == n2 .and. size(raw, 2) == n1) then
            values = transpose(raw)
        else
            values = raw
        end if
    end subroutine

    subroutine project(f, theta, zeta, nfp, mmax, nmax, c, s)
        ! least-squares projection of f(theta_i, zeta_j) on the GVEC modes,
        ! exact for uniform grids resolving them
        real(dp), intent(in) :: f(:, :), theta(:), zeta(:)
        integer, intent(in) :: nfp, mmax, nmax
        real(dp), intent(out) :: c(-IN_NMAX:, 0:), s(-IN_NMAX:, 0:)
        integer :: m, n, i, j
        real(dp) :: a, w
        c = 0; s = 0
        do m = 0, mmax
            do n = -nmax, nmax
                if (m == 0 .and. n < 0) cycle
                w = 2.0_dp/(size(theta)*size(zeta))
                if (m == 0 .and. n == 0) w = w/2
                do j = 1, size(zeta)
                    do i = 1, size(theta)
                        a = m*theta(i) - n*nfp*zeta(j)
                        c(n, m) = c(n, m) + w*f(i, j)*cos(a)
                        s(n, m) = s(n, m) + w*f(i, j)*sin(a)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine cylindrical_modes(surface, eq)
        ! hmap 1: R = X1, Z = X2, phi = -zeta, so (m, n) -> (m, -n)
        type(gvec_surface_t), intent(in) :: surface
        type(imported_equilibrium_t), intent(inout) :: eq
        integer :: n
        do n = -IN_NMAX, IN_NMAX
            eq%rbc(-n, :) = surface%c1(n, :); eq%rbs(-n, :) = surface%s1(n, :)
            eq%zbc(-n, :) = surface%c2(n, :); eq%zbs(-n, :) = surface%s2(n, :)
        end do
    end subroutine

    subroutine frenet_axis(p, surface, eq)
        type(kv_store_t), intent(in) :: p
        type(gvec_surface_t), intent(inout) :: surface
        type(imported_equilibrium_t), intent(inout) :: eq
        integer :: nmax
        logical :: ok
        nmax = p%int('hmap_n_max', -1)
        if (nmax < 0) then
            call import_fail(eq, 'GVEC hmap 20 without hmap_n_max')
            return
        end if
        surface%nfp = p%int('hmap_nfp', surface%nfp)
        surface%omnig = p%logical('hmap_omnig', .false.)
        call coefficients('hmap_rc', surface%rc)
        call coefficients('hmap_rs', surface%rs)
        call coefficients('hmap_zc', surface%zc)
        call coefficients('hmap_zs', surface%zs)
    contains
        subroutine coefficients(name, c)
            character(*), intent(in) :: name
            real(dp), allocatable, intent(out) :: c(:)
            real(dp), allocatable :: v(:)
            allocate (c(0:nmax))
            c = 0
            call p%reals(name, v, ok)
            if (ok) c(0:min(nmax, size(v) - 1)) = v(:min(nmax, size(v) - 1) + 1)
        end subroutine
    end subroutine

    subroutine gframe_axis(path, surface, eq)
        ! axis/xyz, Nxyz, Bxyz over one full turn; rotated into the frame
        ! co-rotating by winding*zeta they are field-period periodic and are
        ! held by their Fourier modes up to axis/n_max (hmap_axisNB.F90)
        character(*), intent(in) :: path
        type(gvec_surface_t), intent(inout) :: surface
        type(imported_equilibrium_t), intent(inout) :: eq
        type(hdf5_file_t) :: h
        type(fortio_status_t) :: status
        real(dp), allocatable :: zeta(:), v(:, :, :), raw(:, :)
        real(dp) :: dphi, phi_prev, phi, z, total
        character(16), parameter :: names(3) = ['axis/xyz(::)    ', 'axis/Nxyz(::)   ', &
            'axis/Bxyz(::)   ']
        integer :: nfp, nzeta, nmax, k, q, j, nfull
        call h%open(path, status)
        if (.not. status%ok()) then
            call import_fail(eq, 'cannot read the G-frame file '//path)
            return
        end if
        call read_integer(h, 'NFP', nfp, status)
        call read_integer(h, 'axis/nzeta', nzeta, status)
        call read_integer(h, 'axis/n_max', nmax, status)
        call h%read_r64_1('axis/zeta(:)', zeta, status)
        nfull = nfp*nzeta
        allocate (v(3, nfull, 3))
        do q = 1, 3
            call read_grid(h, trim(names(q)), 3, nfull, raw, status)
            if (status%ok()) v(:, :, q) = raw
        end do
        call h%close(status)
        if (.not. status%ok()) then
            call import_fail(eq, 'incomplete G-frame file '//path)
            return
        end if
        if (nfp /= surface%nfp) call import_note(eq, 'nfp taken from the G-frame file')
        surface%nfp = nfp
        ! winding number of the axis about the vertical axis
        total = 0; phi_prev = atan2(v(2, 1, 1), v(1, 1, 1))
        do k = 2, nfull + 1
            j = modulo(k - 1, nfull) + 1
            phi = atan2(v(2, j, 1), v(1, j, 1))
            dphi = modulo(phi - phi_prev + pi, 2*pi) - pi
            total = total + dphi
            phi_prev = phi
        end do
        surface%winding = nint(total/(2*pi))
        if (abs(surface%winding) /= 1) then
            call import_refuse(eq, 'the G-frame axis winds more than once about the vertical axis '// &
                '(a knot); the boundary is not single-valued over the cylindrical angle')
            return
        end if
        ! co-rotating components, sampled on the first field period
        allocate (surface%hc(3, 0:nmax, 3), surface%hs(3, 0:nmax, 3))
        surface%hc = 0; surface%hs = 0
        do q = 1, 3
            do j = 1, nzeta
                z = zeta(1) + (j - 1)*2*pi/nfull
                raw = reshape(rotate(v(:, j, q), -surface%winding*z), [3, 1])
                do k = 0, nmax
                    surface%hc(:, k, q) = surface%hc(:, k, q) + raw(:, 1)*cos(k*nfp*z)*merge(1, 2, k == 0)
                    surface%hs(:, k, q) = surface%hs(:, k, q) + raw(:, 1)*sin(k*nfp*z)*2
                end do
            end do
        end do
        surface%hc = surface%hc/nzeta; surface%hs = surface%hs/nzeta
    end subroutine

    pure function rotate(x, angle) result(y)
        real(dp), intent(in) :: x(3), angle
        real(dp) :: y(3)
        y = [cos(angle)*x(1) - sin(angle)*x(2), sin(angle)*x(1) + cos(angle)*x(2), x(3)]
    end function

    pure subroutine surface_q(surface, theta, zeta, q1, q2)
        type(gvec_surface_t), intent(in) :: surface
        real(dp), intent(in) :: theta, zeta
        real(dp), intent(out) :: q1, q2
        integer :: k
        real(dp) :: a, c, s
        q1 = 0; q2 = 0
        do k = 1, size(surface%mode_m)
            a = surface%mode_m(k)*theta - surface%mode_n(k)*surface%nfp*zeta
            c = cos(a); s = sin(a)
            q1 = q1 + surface%mode_c(1, k)*c + surface%mode_c(2, k)*s
            q2 = q2 + surface%mode_c(3, k)*c + surface%mode_c(4, k)*s
        end do
    end subroutine

    subroutine compact_modes(surface)
        ! the nonzero modes as a list, for fast evaluation
        type(gvec_surface_t), intent(inout) :: surface
        integer :: m, n
        allocate (surface%mode_m(0), surface%mode_n(0), surface%mode_c(4, 0))
        do m = 0, IN_MMAX
            do n = -IN_NMAX, IN_NMAX
                if (surface%c1(n, m) == 0 .and. surface%s1(n, m) == 0 .and. &
                    surface%c2(n, m) == 0 .and. surface%s2(n, m) == 0) cycle
                surface%mode_m = [surface%mode_m, m]
                surface%mode_n = [surface%mode_n, n]
                surface%mode_c = reshape([surface%mode_c, surface%c1(n, m), surface%s1(n, m), &
                    surface%c2(n, m), surface%s2(n, m)], [4, size(surface%mode_m)])
            end do
        end do
    end subroutine

    pure function surface_xyz(surface, theta, zeta) result(x)
        type(gvec_surface_t), intent(in) :: surface
        real(dp), intent(in) :: theta, zeta
        real(dp) :: x(3), q1, q2, x0(3), nv(3), bv(3), d1(3), d2(3), t(3), sigma
        real(dp) :: r, r1, r2, z1, z2, zz, c, s, w
        integer :: n, q
        call surface_q(surface, theta, zeta, q1, q2)
        select case (surface%frame)
        case (20)
            ! X0 = (R0 cos zeta, R0 sin zeta, Z0), Frenet frame of X0
            r = 0; r1 = 0; r2 = 0; zz = 0; z1 = 0; z2 = 0
            do n = 0, size(surface%rc) - 1
                w = n*surface%nfp
                c = cos(w*zeta); s = sin(w*zeta)
                r = r + surface%rc(n)*c + surface%rs(n)*s
                r1 = r1 + w*(-surface%rc(n)*s + surface%rs(n)*c)
                r2 = r2 - w**2*(surface%rc(n)*c + surface%rs(n)*s)
                zz = zz + surface%zc(n)*c + surface%zs(n)*s
                z1 = z1 + w*(-surface%zc(n)*s + surface%zs(n)*c)
                z2 = z2 - w**2*(surface%zc(n)*c + surface%zs(n)*s)
            end do
            c = cos(zeta); s = sin(zeta)
            x0 = [r*c, r*s, zz]
            d1 = [r1*c - r*s, r1*s + r*c, z1]
            d2 = [r2*c - 2*r1*s - r*c, r2*s + 2*r1*c - r*s, z2]
            t = d1/norm2(d1)
            bv = cross(d1, d2)
            bv = bv/norm2(bv)
            nv = cross(bv, t)
            sigma = 1
            if (surface%omnig) sigma = sign(1.0_dp, sin(surface%nfp*zeta))
            x = x0 + sigma*(q1*nv + q2*bv)
        case (21)
            do q = 1, 3
                d1 = surface%hc(:, 0, q)
                do n = 1, size(surface%hc, 2) - 1
                    d1 = d1 + surface%hc(:, n, q)*cos(n*surface%nfp*zeta) &
                        + surface%hs(:, n, q)*sin(n*surface%nfp*zeta)
                end do
                d1 = rotate(d1, surface%winding*zeta)
                if (q == 1) x0 = d1
                if (q == 2) nv = d1
                if (q == 3) bv = d1
            end do
            x = x0 + q1*nv + q2*bv
        case default
            x = [q1*cos(zeta), -q1*sin(zeta), q2]
        end select
    end function

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    subroutine cartesian_refit(surface, eq, s_zeta)
        ! The surface x(theta, zeta) cut by the planes phi = const: on each
        ! theta line zeta(phi) is found by bisection, then R(theta, phi),
        ! Z(theta, phi) are fitted in VMEC's basis. Admitted only if phi
        ! is monotone in zeta on every theta line. The fit error at the
        ! midpoints of the fit grid is the stated geometry error.
        type(gvec_surface_t), intent(in) :: surface
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), intent(out) :: s_zeta
        integer, parameter :: ntable = 512
        integer :: mfit, nfit, nt, np, i, j, m
        real(dp), allocatable :: r(:, :), z(:, :), table(:, :)
        real(dp) :: theta, phi, rr, zz, rf, zf, error, radius
        logical :: monotone
        mfit = 0
        do m = 0, IN_MMAX
            if (any(surface%c1(:, m) /= 0) .or. any(surface%s1(:, m) /= 0) .or. &
                any(surface%c2(:, m) /= 0) .or. any(surface%s2(:, m) /= 0)) mfit = m
        end do
        ! the change of toroidal parameter couples poloidal modes; modes
        ! beyond nested_equilibrium's limits are cut and counted later
        mfit = min(IN_MMAX/2, 2*mfit + 8)
        nfit = IN_NMAX/2
        nt = 2*mfit + 2; np = 2*nfit + 2
        eq%nfp = surface%nfp
        ! phi on a fine zeta table over one full turn, per theta line
        allocate (table(0:ntable, nt))
        monotone = .true.
        do i = 1, nt
            call phi_line(2*pi*(i - 1)/nt, table(:, i))
        end do
        s_zeta = sign(1.0_dp, table(ntable, 1) - table(0, 1))
        do i = 1, nt
            if (any(s_zeta*(table(1:, i) - table(:ntable - 1, i)) <= 0)) monotone = .false.
            if (abs(abs(table(ntable, i) - table(0, i)) - 2*pi) > 1e-6_dp) monotone = .false.
        end do
        if (.not. monotone) then
            call import_refuse(eq, 'the boundary is not single-valued over the cylindrical angle '// &
                '(phi is not monotone along zeta)')
            return
        end if
        allocate (r(nt, np), z(nt, np))
        do j = 1, np
            phi = 2*pi*(j - 1)/(surface%nfp*np)
            do i = 1, nt
                call section(2*pi*(i - 1)/nt, phi, table(:, i), r(i, j), z(i, j))
            end do
        end do
        call boundary_fit(eq, r, z, mfit, nfit)
        ! check off the fit grid
        error = 0
        radius = sqrt(sum(eq%rbc(:, 1:)**2 + eq%zbs(:, 1:)**2 + eq%rbs(:, 1:)**2 + eq%zbc(:, 1:)**2)/2)
        do i = 1, nt, max(1, nt/16)
            theta = 2*pi*(i - 0.5_dp)/nt
            call phi_line(theta, table(:, 1))
            do j = 1, np
                phi = 2*pi*(j - 0.5_dp)/(surface%nfp*np)
                call section(theta, phi, table(:, 1), rr, zz)
                call fitted(theta, phi, rf, zf)
                error = max(error, hypot(rr - rf, zz - zf))
            end do
        end do
        eq%geometry_error = eq%geometry_error + error
        call import_note(eq, 'boundary refitted on cylindrical planes; largest midpoint deviation '// &
            text(error)//' m, relative to the minor radius '//text(error/radius))
        if (error > 1e-3_dp*radius) call import_refuse(eq, 'the boundary is not resolved by '// &
            text(real(nfit, dp))//' toroidal modes over the cylindrical angle')
    contains
        subroutine phi_line(theta, phis)
            real(dp), intent(in) :: theta
            real(dp), intent(out) :: phis(0:)
            real(dp) :: x(3), previous
            integer :: k
            do k = 0, ntable
                x = surface_xyz(surface, theta, 2*pi*k/ntable)
                phis(k) = atan2(x(2), x(1))
                if (k > 0) then
                    previous = phis(k - 1)
                    phis(k) = previous + modulo(phis(k) - previous + pi, 2*pi) - pi
                end if
            end do
        end subroutine
        subroutine section(theta, phi, phis, r, z)
            real(dp), intent(in) :: theta, phi, phis(0:)
            real(dp), intent(out) :: r, z
            real(dp) :: target, a, b, mid, x(3), f, fa, fb
            integer :: k, it, side
            ! the turn of the table holding phi
            target = phis(0) + modulo(s_zeta*(phi - phis(0)), 2*pi)*s_zeta
            k = 0
            do while (k < ntable - 1)
                if (s_zeta*(phis(k + 1) - target) >= 0) exit
                k = k + 1
            end do
            ! Illinois regula falsi on the bracket of the table
            a = 2*pi*k/ntable; b = 2*pi*(k + 1)/ntable
            fa = phis(k) - target; fb = phis(k + 1) - target
            side = 0
            mid = a
            do it = 1, 100
                if (fb == fa) exit
                mid = (a*fb - b*fa)/(fb - fa)
                x = surface_xyz(surface, theta, mid)
                f = atan2(x(2), x(1))
                f = modulo(f - target + pi, 2*pi) - pi
                if (abs(f) < 1e-15_dp .or. b - a < 1e-15_dp) exit
                if (sign(1.0_dp, f) == sign(1.0_dp, fb)) then
                    b = mid; fb = f
                    if (side == 1) fa = fa/2
                    side = 1
                else
                    a = mid; fa = f
                    if (side == -1) fb = fb/2
                    side = -1
                end if
            end do
            x = surface_xyz(surface, theta, mid)
            r = hypot(x(1), x(2)); z = x(3)
        end subroutine
        subroutine fitted(theta, phi, r, z)
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: r, z
            integer :: m, n
            real(dp) :: a
            r = 0; z = 0
            do m = 0, mfit
                do n = -nfit, nfit
                    a = m*theta - n*surface%nfp*phi
                    r = r + eq%rbc(n, m)*cos(a) + eq%rbs(n, m)*sin(a)
                    z = z + eq%zbc(n, m)*cos(a) + eq%zbs(n, m)*sin(a)
                end do
            end do
        end subroutine
    end subroutine

    function text(x) result(out)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: out
        character(32) :: buffer
        write (buffer, '(es10.3)') x
        out = trim(adjustl(buffer))
    end function

    subroutine gvec_profile(p, name, profile, eq)
        ! polynomial, bspline or interpolation (not-a-knot cubic spline) in
        ! s = rho^2, times <name>_scale
        type(kv_store_t), intent(in) :: p
        character(*), intent(in) :: name
        type(import_profile_t), intent(out) :: profile
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), allocatable :: coefs(:), knots(:), x(:), y(:)
        real(dp) :: scale
        logical :: ok, ok2
        character(len=:), allocatable :: kind
        scale = p%real(name//'_scale', 1.0_dp)
        kind = lower(p%text(name//'_type', 'polynomial'))
        select case (kind)
        case ('polynomial')
            call p%reals(name//'_coefs', coefs, ok)
            if (.not. ok) coefs = [0.0_dp]
            profile = power_series_profile(scale*coefs)
        case ('bspline')
            call p%reals(name//'_coefs', coefs, ok)
            call p%reals(name//'_knots', knots, ok2)
            if (.not. (ok .and. ok2)) then
                call import_fail(eq, name//': B-spline needs coefs and knots')
                return
            end if
            call bspline_profile(knots, scale*coefs, profile, eq, name)
        case ('interpolation')
            call p%reals(name//'_rho2', x, ok)
            call p%reals(name//'_vals', y, ok2)
            if (.not. (ok .and. ok2) .or. size(x) /= size(y)) then
                call import_fail(eq, name//': interpolation needs rho2 and vals of one length')
                return
            end if
            if (lower(p%text(name//'_bc_type_axis', 'not_a_knot')) /= 'not_a_knot' .or. &
                lower(p%text(name//'_bc_type_edge', 'not_a_knot')) /= 'not_a_knot') then
                call import_refuse(eq, name//': only not-a-knot spline ends are imported')
                return
            end if
            profile = spline_profile(x, scale*y, ok)
            if (.not. ok .or. abs(x(1)) > 1e-12_dp .or. abs(x(size(x)) - 1) > 1e-12_dp) &
                call import_fail(eq, name//': spline nodes must run from s = 0 to 1, at least 4')
        case default
            call import_refuse(eq, name//': GVEC profile type '//kind)
        end select
    end subroutine

    subroutine bspline_profile(t, c, profile, eq, name)
        ! Degree p = size(t) - size(c) - 1. With p <= 3 and C1 knots the
        ! spline is a C1 piecewise cubic and its values and slopes at the
        ! distinct knots hold it exactly; otherwise 257 uniform knots are
        ! added and the error noted.
        real(dp), intent(in) :: t(:), c(:)
        type(import_profile_t), intent(out) :: profile
        type(imported_equilibrium_t), intent(inout) :: eq
        character(*), intent(in) :: name
        real(dp), allocatable :: knots(:), values(:), slopes(:)
        real(dp) :: error, scale, x, h, v
        integer :: degree, j, k, multiplicity
        logical :: exact
        degree = size(t) - size(c) - 1
        if (degree < 0) then
            call import_fail(eq, name//': B-spline has more coefficients than knots allow')
            return
        end if
        allocate (knots(0))
        exact = degree <= 3
        do j = 1, size(t)
            if (size(knots) > 0) then
                if (t(j) <= knots(size(knots))) cycle
            end if
            multiplicity = count(t == t(j))
            if (t(j) > t(1) .and. t(j) < t(size(t)) .and. degree - multiplicity < 1) exact = .false.
            knots = [knots, t(j)]
        end do
        if (.not. exact) knots = merged(knots, [((k - 1)/256.0_dp, k=1, 257)])
        allocate (values(size(knots)), slopes(size(knots)))
        do j = 1, size(knots)
            call bspline_eval(t, c, degree, knots(j), values(j), slopes(j))
        end do
        profile = hermite_profile(knots, values, slopes)
        if (.not. exact) then
            error = 0; scale = maxval(abs(values))
            do j = 1, size(knots) - 1
                h = knots(j + 1) - knots(j)
                do k = 1, 3
                    x = knots(j) + k*h/4
                    call hermite_value(knots, values, slopes, x, v)
                    error = max(error, abs(v - bspline_at(x)))
                end do
            end do
            call import_note(eq, name//': B-spline of degree above 3 as Hermite, relative error '// &
                text(error/max(scale, tiny(1.0_dp))))
        end if
    contains
        real(dp) function bspline_at(x)
            real(dp), intent(in) :: x
            real(dp) :: d
            call bspline_eval(t, c, degree, x, bspline_at, d)
        end function
    end subroutine

    pure function merged(a, b) result(out)
        real(dp), intent(in) :: a(:), b(:)
        real(dp), allocatable :: out(:)
        real(dp) :: all(size(a) + size(b)), x
        integer :: i, j
        all = [a, b]
        do i = 2, size(all)
            x = all(i); j = i - 1
            do while (j >= 1)
                if (all(j) <= x) exit
                all(j + 1) = all(j); j = j - 1
            end do
            all(j + 1) = x
        end do
        allocate (out(0))
        do i = 1, size(all)
            if (size(out) > 0) then
                if (all(i) - out(size(out)) <= 1e-12_dp) cycle
            end if
            out = [out, all(i)]
        end do
    end function

    pure subroutine bspline_eval(t, c, p, x, value, slope)
        ! Cox-de Boor basis of degree p and p - 1 at x; the last span is
        ! closed so that x = t(end) is inside
        real(dp), intent(in) :: t(:), c(:), x
        integer, intent(in) :: p
        real(dp), intent(out) :: value, slope
        real(dp) :: b(size(t) - 1), bprev(size(t) - 1), left, right
        integer :: i, k, n, span
        n = size(t) - 1
        b = 0
        span = 0
        do i = 1, n
            if (t(i) < t(i + 1) .and. x >= t(i) .and. x < t(i + 1)) span = i
        end do
        if (span == 0) then
            do i = n, 1, -1
                if (t(i) < t(i + 1)) then
                    span = i
                    exit
                end if
            end do
        end if
        b(span) = 1
        bprev = b
        do k = 1, p
            bprev = b
            do i = 1, n - k
                left = 0; right = 0
                if (t(i + k) > t(i)) left = (x - t(i))/(t(i + k) - t(i))*bprev(i)
                if (t(i + k + 1) > t(i + 1)) right = (t(i + k + 1) - x)/(t(i + k + 1) - t(i + 1))*bprev(i + 1)
                b(i) = left + right
            end do
            b(n - k + 1:) = 0
        end do
        value = sum(c*b(:size(c)))
        ! derivative from the degree p - 1 basis held in bprev
        slope = 0
        if (p >= 1) then
            do i = 1, size(c)
                if (t(i + p) > t(i)) slope = slope + c(i)*p/(t(i + p) - t(i))*bprev(i)
                if (i + 1 <= n) then
                    if (t(i + p + 1) > t(i + 1)) slope = slope - c(i)*p/(t(i + p + 1) - t(i + 1))*bprev(i + 1)
                end if
            end do
        end if
    end subroutine
end module milca_import_gvec
