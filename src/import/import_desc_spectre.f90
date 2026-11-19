module milca_import_desc_spectre
    ! DESC text inputs and SPECTRE TOML inputs to the intermediate form of
    ! milca_equilibrium_import.
    !
    ! DESC (desc/input_reader.py, MIT licence): R, Z = sum c G_m(theta)
    ! G_n(nfp phi), G_k(x) = cos(k x) for k >= 0 and sin(|k| x) for k < 0;
    ! profiles are power series in rho = sqrt(s) indexed by l; its VMEC
    ! reader keeps VMEC's angles and iota and stores the current as
    ! signgs I_V, signgs = -1 for a counterclockwise boundary.
    !
    ! SPECTRE (spectre/converters/vmec_helpers.py, MIT licence): R, Z = sum
    ! c cos(m theta - n nfp phi) + s sin(...) keyed "(m, n)"; theta is VMEC's
    ! reversed while iota is kept; the enclosed volume and sheet currents are
    ! -mu0 I_V at the interfaces of the relaxed (MRxMHD) volumes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_import_text, only: kv_store_t, read_lines, read_toml, lower, parse_reals
    use milca_equilibrium_import, only: imported_equilibrium_t, import_profile_t, import_note, &
        import_refuse, import_fail, IMPORT_OK, boundary_fit, boundary_area, power_series_profile, &
        hermite_profile, spline_profile, hermite_value, profile_value, IN_MMAX, IN_NMAX
    implicit none
    private
    public :: import_desc, import_spectre
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi

contains

    subroutine import_desc(path, eq)
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        character(len=:), allocatable :: lines(:), line, mode
        real(dp), allocatable :: rmn(:, :), zmn(:, :), r(:, :), z(:, :)
        real(dp) :: p(0:64), io(0:64), cu(0:64), value, psi, signgs
        logical :: ok, has_iota, has_current, has_labels
        integer :: i, l, m, n, mmax, nmax, nt, np, j, k, nfp, top
        character(8) :: names(8)
        real(dp) :: values(8)
        integer :: nlabels, labels(3)
        eq%source = path
        eq%format = 'DESC'
        call read_lines(path, lines, ok)
        if (.not. ok) then
            call import_fail(eq, 'cannot read '//trim(path))
            return
        end if
        allocate (rmn(-IN_MMAX:IN_MMAX, -IN_NMAX:IN_NMAX), zmn(-IN_MMAX:IN_MMAX, -IN_NMAX:IN_NMAX))
        rmn = 0; zmn = 0; p = 0; io = 0; cu = 0
        nfp = 1; psi = 1; mode = 'lcfs'
        has_iota = .false.; has_current = .false.
        mmax = 0; nmax = 0; top = 0
        do i = 1, size(lines)
            line = lines(i)
            if (index(line, '#') > 0) line = line(:index(line, '#') - 1)
            if (len_trim(line) == 0) cycle
            call row(line, names, values, nlabels, labels, has_labels)
            if (.not. has_labels) then
                ! global "key = value"
                if (index(line, '=') == 0) cycle
                select case (lower(trim(adjustl(line(:index(line, '=') - 1)))))
                case ('nfp')
                    nfp = nint(number(line(index(line, '=') + 1:)))
                case ('psi')
                    psi = number(line(index(line, '=') + 1:))
                case ('bdry_mode')
                    mode = lower(trim(adjustl(line(index(line, '=') + 1:))))
                end select
                cycle
            end if
            if (any(names(:nlabels) == 'r1')) then
                m = labels(2); n = labels(3)
                if (abs(m) > IN_MMAX .or. abs(n) > IN_NMAX) then
                    call import_refuse(eq, 'DESC boundary mode beyond the working limit')
                    return
                end if
                do k = 1, nlabels
                    if (names(k) == 'r1') rmn(m, n) = values(k)
                    if (names(k) == 'z1') zmn(m, n) = values(k)
                end do
                if (rmn(m, n) /= 0 .or. zmn(m, n) /= 0) then
                    mmax = max(mmax, abs(m)); nmax = max(nmax, abs(n))
                end if
            else if (any(names(:nlabels) == 'p') .or. any(names(:nlabels) == 'i') .or. &
                    any(names(:nlabels) == 'c')) then
                l = labels(1)
                if (l < 0 .or. l > 64) then
                    call import_refuse(eq, 'DESC profile power above rho^64')
                    return
                end if
                top = max(top, l)
                do k = 1, nlabels
                    if (names(k) == 'p') p(l) = values(k)
                    if (names(k) == 'i') then
                        io(l) = values(k); has_iota = .true.
                    end if
                    if (names(k) == 'c') then
                        cu(l) = values(k); has_current = .true.
                    end if
                end do
            end if
        end do
        if (mode /= 'lcfs') then
            call import_refuse(eq, 'DESC bdry_mode '//mode//' (lcfs only)')
            return
        end if
        eq%nfp = nfp
        eq%phiedge = psi
        nt = 2*mmax + 2; np = 2*nmax + 2
        allocate (r(nt, np), z(nt, np))
        r = 0; z = 0
        do j = 1, np
            do i = 1, nt
                do m = -mmax, mmax
                    do n = -nmax, nmax
                        if (rmn(m, n) == 0 .and. zmn(m, n) == 0) cycle
                        value = g(m, 2*pi*(i - 1)/nt)*g(n, 2*pi*(j - 1)/np)
                        r(i, j) = r(i, j) + rmn(m, n)*value
                        z(i, j) = z(i, j) + zmn(m, n)*value
                    end do
                end do
            end do
        end do
        call boundary_fit(eq, r, z, mmax, nmax)
        signgs = merge(-1.0_dp, 1.0_dp, boundary_area(eq) > 0)
        eq%pressure = rho_series(eq, p(:top), 'pressure')
        if (has_current .and. .not. has_iota) then
            eq%current = rho_series(eq, signgs*cu(:top), 'current')
            eq%prescribe_current = .true.
        else
            eq%iota = rho_series(eq, io(:top), 'iota')
            if (has_current) call import_note(eq, 'DESC input gives iota and current; iota is kept')
        end if
        if (.not. (has_iota .or. has_current)) then
            eq%current = power_series_profile([0.0_dp])
            eq%prescribe_current = .true.
            call import_note(eq, 'no iota or current given: DESC''s default zero current')
        end if
    contains
        pure real(dp) function g(k, x)
            integer, intent(in) :: k
            real(dp), intent(in) :: x
            if (k >= 0) then
                g = cos(k*x)
            else
                g = sin(-k*x)
            end if
        end function
        real(dp) function number(text)
            character(*), intent(in) :: text
            real(dp), allocatable :: v(:)
            logical :: parsed
            call parse_reals(text, v, parsed)
            number = 0
            if (parsed .and. size(v) > 0) number = v(1)
        end function
    end subroutine

    function rho_series(eq, c, name) result(profile)
        ! sum c_l rho^l: a power series in s when only even l occur,
        ! otherwise Hermite on 257 knots with the error noted
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), intent(in) :: c(0:)
        character(*), intent(in) :: name
        type(import_profile_t) :: profile
        real(dp) :: knots(257), v(257), d(257), x, exact, approx, error, scale
        integer :: ll, jj
        character(32) :: buffer
        if (all(c(1::2) == 0)) then
            profile = power_series_profile(c(0::2))
            return
        end if
        do jj = 1, 257
            knots(jj) = (jj - 1)/256.0_dp
            v(jj) = series(c, knots(jj)); d(jj) = 0
            do ll = 1, ubound(c, 1)
                if (knots(jj) > 0) d(jj) = d(jj) + c(ll)*ll/2.0_dp*knots(jj)**(ll/2.0_dp - 1)
            end do
        end do
        ! one-sided at the axis, where odd powers of rho have no s-slope
        d(1) = (v(2) - v(1))/knots(2)
        profile = hermite_profile(knots, v, d)
        error = 0; scale = maxval(abs(v))
        do jj = 1, 256
            do ll = 1, 3
                x = knots(jj) + ll*(knots(jj + 1) - knots(jj))/4
                exact = series(c, x)
                call hermite_value(knots, v, d, x, approx)
                error = max(error, abs(exact - approx))
            end do
        end do
        write (buffer, '(es10.3)') error/max(scale, tiny(1.0_dp))
        call import_note(eq, name//': odd powers of rho, Hermite on 257 knots, relative error '// &
            trim(adjustl(buffer)))
    end function

    pure real(dp) function series(c, s)
        real(dp), intent(in) :: c(0:), s
        integer :: q
        series = 0
        do q = 0, ubound(c, 1)
            series = series + c(q)*sqrt(s)**q
        end do
    end function

    subroutine row(line, names, values, nlabels, labels, has_labels)
        ! "l: 0  m: -2  n: 0  R1 = 1.0  Z1 = 2.0" -> labels (l, m, n) and
        ! the named values
        character(*), intent(in) :: line
        character(8), intent(out) :: names(8)
        real(dp), intent(out) :: values(8)
        integer, intent(out) :: nlabels, labels(3)
        logical, intent(out) :: has_labels
        character(len=len(line)) :: t
        character(len=64) :: words(64)
        character(len=:), allocatable :: w
        integer :: count, i, start, ios
        t = line
        do i = 1, len(t)
            if (t(i:i) == achar(9)) t(i:i) = ' '
        end do
        count = 0; i = 1
        do while (i <= len_trim(t) .and. count < 64)
            if (t(i:i) == ' ') then
                i = i + 1
                cycle
            end if
            start = i
            do while (i <= len(t))
                if (t(i:i) == ' ') exit
                i = i + 1
            end do
            count = count + 1
            words(count) = t(start:i - 1)
        end do
        labels = 0; nlabels = 0; has_labels = .false.
        i = 1
        do while (i <= count)
            w = trim(words(i))
            block
                if (len(w) >= 2 .and. w(len(w):len(w)) == ':' .and. i < count) then
                    select case (lower(w))
                    case ('l:')
                        read (words(i + 1), *, iostat=ios) labels(1)
                        has_labels = .true.
                    case ('m:')
                        read (words(i + 1), *, iostat=ios) labels(2)
                        has_labels = .true.
                    case ('n:')
                        read (words(i + 1), *, iostat=ios) labels(3)
                        has_labels = .true.
                    end select
                    i = i + 2
                    cycle
                end if
                if (i + 2 <= count .and. trim(words(min(i + 1, count))) == '=' .and. nlabels < 8) then
                    nlabels = nlabels + 1
                    names(nlabels) = lower(w)
                    read (words(i + 2), *, iostat=ios) values(nlabels)
                    if (ios /= 0) values(nlabels) = 0
                    i = i + 3
                    cycle
                end if
            end block
            i = i + 1
        end do
    end subroutine

    subroutine import_spectre(path, eq)
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        type(kv_store_t) :: p
        integer, allocatable :: m(:), n(:)
        real(dp), allocatable :: v(:), tflux(:), pressure(:), iota(:), ivolume(:), isurf(:)
        real(dp), allocatable :: s(:), enclosed(:)
        integer :: nvol, k, f, lconstraint
        logical :: ok
        character(4), parameter :: fields(4) = ['rbc ', 'rbs ', 'zbc ', 'zbs ']
        eq%source = path
        eq%format = 'SPECTRE'
        call read_toml(path, p, ok)
        if (.not. ok) then
            call import_fail(eq, 'cannot read '//trim(path))
            return
        end if
        if (p%int('physics_igeometry', 3) /= 3) then
            call import_refuse(eq, 'SPECTRE igeometry other than 3 is a slab or a periodic cylinder; '// &
                'nested_equilibrium solves toroidal equilibria')
            return
        end if
        if (p%int('physics_lfreebound', 0) /= 0) then
            call import_refuse(eq, 'free-boundary SPECTRE input; nested_equilibrium is fixed-boundary')
            return
        end if
        nvol = p%int('physics_nvol', 1)
        call p%reals('physics_tflux', tflux, ok)
        call p%reals('physics_pressure', pressure, ok)
        if (.not. ok) pressure = [0.0_dp]
        if (any(p%real('physics_pscale', 1.0_dp)*pressure /= 0)) then
            call import_refuse(eq, 'SPECTRE''s pressure is stepped across its volumes; '// &
                'nested_equilibrium''s profiles are continuously differentiable')
            return
        end if
        eq%nfp = p%int('physics_nfp', 1)
        eq%phiedge = p%real('physics_phiedge', 1.0_dp)
        do f = 1, 4
            call p%modes('physics_'//trim(fields(f)), m, n, v)
            do k = 1, size(v)
                if (m(k) < 0 .or. m(k) > IN_MMAX .or. abs(n(k)) > IN_NMAX) cycle
                select case (f)
                case (1)
                    eq%rbc(n(k), m(k)) = v(k)
                case (2)
                    eq%rbs(n(k), m(k)) = v(k)
                case (3)
                    eq%zbc(n(k), m(k)) = v(k)
                case (4)
                    eq%zbs(n(k), m(k)) = v(k)
                end select
            end do
        end do
        if (size(tflux) < nvol .or. nvol < 3) then
            call import_refuse(eq, 'SPECTRE input with fewer than 3 volumes or interfaces')
            return
        end if
        s = [0.0_dp, tflux(:nvol)/tflux(nvol)]
        call p%reals('physics_iota', iota, ok)
        lconstraint = p%int('physics_lconstraint', -1)
        select case (lconstraint)
        case (1)
            if (size(iota) < nvol + 1) then
                call import_fail(eq, 'SPECTRE lconstraint = 1 needs iota on the axis and interfaces')
                return
            end if
            eq%iota = spline_profile(s, -iota(:nvol + 1), ok)
            call import_note(eq, 'MRxMHD with iota held on the interfaces: iota between them is the '// &
                'not-a-knot spline through the interface values (a model choice)')
        case (3)
            call p%reals('physics_ivolume', ivolume, ok)
            call p%reals('physics_isurf', isurf, ok)
            if (size(ivolume) < nvol .or. size(isurf) < nvol) then
                call import_fail(eq, 'SPECTRE lconstraint = 3 needs ivolume and isurf')
                return
            end if
            enclosed = [0.0_dp, (-(ivolume(k) + sum(isurf(:k)))/mu0, k=1, nvol)]
            eq%current = spline_profile(s, enclosed, ok)
            eq%prescribe_current = .true.
            if (size(iota) >= nvol + 1) eq%iota = spline_profile(s, -iota(:nvol + 1), ok)
            call import_note(eq, 'MRxMHD with currents held: the enclosed current between the '// &
                'interfaces is the not-a-knot spline through the interface values (a model choice)')
        case default
            call import_refuse(eq, 'SPECTRE lconstraint other than 1 (iota) or 3 (currents) has no '// &
                'nested ideal-MHD counterpart')
            return
        end select
        if (.not. ok) call import_fail(eq, 'SPECTRE interface profile needs at least 3 volumes')
    end subroutine
end module milca_import_desc_spectre
