module milca_equilibrium_import
    ! A foreign fixed-boundary equilibrium input in one intermediate form and
    ! its conversion to nested_equilibrium's namelist
    ! (docs/equilibrium-import.md).
    !
    ! Intermediate convention (VMEC's): R, Z = sum c cos(m theta - n nfp phi)
    ! + s sin(m theta - n nfp phi) with phi the cylindrical angle,
    ! counterclockwise seen from above, and theta of either orientation; iota
    ! is measured in these angles; the enclosed toroidal current I_V has the
    ! sign of VMEC's CURTOR; phiedge carries its sign. nested_equilibrium
    ! wants theta counterclockwise in (R, Z), phiedge > 0 and its current
    ! I = -sign(phiedge) I_V (docs/nested-equilibrium-input.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init
    implicit none
    private
    public :: import_profile_t, imported_equilibrium_t, IMPORT_OK, IMPORT_UNSUPPORTED, &
        IMPORT_ERROR, IN_MMAX, IN_NMAX, MILCA_NMAX, MILCA_MMAX, boundary_reserve, mode_bounds, import_note, import_refuse, import_fail, &
        boundary_fit, boundary_point, boundary_area, profile_value, power_series_profile, &
        spline_profile, hermite_profile, hermite_value, import_normalize, import_write, &
        orientation_flip
    integer, parameter :: IMPORT_OK = 0, IMPORT_UNSUPPORTED = 3, IMPORT_ERROR = 1
    ! working limits; nested_equilibrium accepts m <= 64, |n| <= 32
    ! input-side limits of the readers that store foreign modes first
    integer, parameter :: IN_MMAX = 128, IN_NMAX = 128
    ! nested_equilibrium's namelist storage
    integer, parameter :: MILCA_MMAX = 64, MILCA_NMAX = 256, MILCA_KMAX = 32, KNOT_MAX = 1024
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: import_profile_t
        ! none, power_series (coefficients of s^k) or hermite
        character(16) :: kind = 'none'
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: knots(:), values(:), slopes(:)
    end type

    type :: imported_equilibrium_t
        integer :: status = IMPORT_OK
        character(1024) :: reason = ''
        character(1024) :: source = '', format = ''
        character(256) :: notes(64) = ''
        integer :: n_notes = 0
        integer :: nfp = 1
        real(dp) :: phiedge = 1
        ! (-n_max:n_max, 0:m_max), sized by boundary_reserve to the modes held
        real(dp), allocatable :: rbc(:, :), rbs(:, :), zbc(:, :), zbs(:, :)
        type(import_profile_t) :: iota, pressure, current
        logical :: prescribe_current = .false.
        ! largest boundary displacement left out by the conversion (m), and
        ! its accepted bound relative to the minor radius
        real(dp) :: geometry_error = 0, geometry_tolerance = 1e-6_dp
    end type

contains

    subroutine import_note(eq, text)
        type(imported_equilibrium_t), intent(inout) :: eq
        character(*), intent(in) :: text
        if (eq%n_notes == size(eq%notes)) return
        eq%n_notes = eq%n_notes + 1
        eq%notes(eq%n_notes) = text
    end subroutine

    subroutine boundary_reserve(eq, mmax, nmax)
        ! boundary arrays covering m <= mmax, |n| <= nmax, keeping held modes
        type(imported_equilibrium_t), intent(inout) :: eq
        integer, intent(in) :: mmax, nmax
        integer :: m0, n0
        m0 = -1; n0 = -1
        if (allocated(eq%rbc)) then
            m0 = ubound(eq%rbc, 2); n0 = ubound(eq%rbc, 1)
            if (m0 >= mmax .and. n0 >= nmax) return
        end if
        call grow(eq%rbc); call grow(eq%rbs); call grow(eq%zbc); call grow(eq%zbs)
    contains
        subroutine grow(c)
            real(dp), allocatable, intent(inout) :: c(:, :)
            real(dp), allocatable :: wider(:, :)
            allocate (wider(-max(nmax, n0):max(nmax, n0), 0:max(mmax, m0)))
            wider = 0
            if (allocated(c)) wider(-n0:n0, 0:m0) = c
            call move_alloc(wider, c)
        end subroutine
    end subroutine

    pure subroutine mode_bounds(eq, mmax, nmax)
        type(imported_equilibrium_t), intent(in) :: eq
        integer, intent(out) :: mmax, nmax
        mmax = -1; nmax = -1
        if (.not. allocated(eq%rbc)) return
        mmax = ubound(eq%rbc, 2); nmax = ubound(eq%rbc, 1)
    end subroutine

    subroutine import_refuse(eq, reason)
        ! a readable input whose physics nested_equilibrium does not model
        type(imported_equilibrium_t), intent(inout) :: eq
        character(*), intent(in) :: reason
        if (eq%status /= IMPORT_OK) return
        eq%status = IMPORT_UNSUPPORTED
        eq%reason = reason
    end subroutine

    subroutine import_fail(eq, reason)
        type(imported_equilibrium_t), intent(inout) :: eq
        character(*), intent(in) :: reason
        if (eq%status == IMPORT_ERROR) return
        eq%status = IMPORT_ERROR
        eq%reason = reason
    end subroutine

    subroutine boundary_fit(eq, r, z, mfit, nfit)
        ! VMEC-basis coefficients of R, Z sampled at theta_i = 2 pi (i-1)/nt,
        ! phi_j = 2 pi (j-1)/(nfp np); exact for data band-limited to
        ! m <= mfit, |n| <= nfit when nt > 2 mfit and np > 2 nfit.
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), intent(in) :: r(:, :), z(:, :)
        integer, intent(in) :: mfit, nfit
        complex(dp), allocatable :: fr(:, :), fz(:, :), er(:, :), ez(:, :)
        integer :: nt, np, i, j, m, n
        complex(dp) :: w
        nt = size(r, 1); np = size(r, 2)
        if (2*mfit >= nt .or. 2*nfit >= np) &
            error stop 'boundary_fit: grid too coarse for the requested modes'
        ! F(m, n) = <f exp(-i (m theta - n psi))>, psi = nfp phi
        allocate (er(nt, -nfit:nfit), ez(nt, -nfit:nfit))
        er = 0; ez = 0
        do n = -nfit, nfit
            do j = 1, np
                w = exp(cmplx(0, n*2*pi*(j - 1)/np, dp))
                er(:, n) = er(:, n) + r(:, j)*w
                ez(:, n) = ez(:, n) + z(:, j)*w
            end do
        end do
        allocate (fr(-nfit:nfit, 0:mfit), fz(-nfit:nfit, 0:mfit))
        fr = 0; fz = 0
        do m = 0, mfit
            do i = 1, nt
                w = exp(cmplx(0, -m*2*pi*(i - 1)/nt, dp))
                fr(:, m) = fr(:, m) + er(i, :)*w
                fz(:, m) = fz(:, m) + ez(i, :)*w
            end do
        end do
        fr = fr/(nt*np); fz = fz/(nt*np)
        if (allocated(eq%rbc)) deallocate (eq%rbc, eq%rbs, eq%zbc, eq%zbs)
        call boundary_reserve(eq, mfit, nfit)
        ! f = c cos(a) + s sin(a) gives F = (c - i s)/2 on each pair, F = c on (0, 0)
        eq%rbc(0, 0) = real(fr(0, 0)); eq%zbc(0, 0) = real(fz(0, 0))
        do m = 0, mfit
            do n = -nfit, nfit
                if (m == 0 .and. n <= 0) cycle
                eq%rbc(n, m) = 2*real(fr(n, m)); eq%rbs(n, m) = -2*aimag(fr(n, m))
                eq%zbc(n, m) = 2*real(fz(n, m)); eq%zbs(n, m) = -2*aimag(fz(n, m))
            end do
        end do
    end subroutine

    pure subroutine boundary_point(eq, theta, phi, r, z)
        type(imported_equilibrium_t), intent(in) :: eq
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: r, z
        integer :: m, n
        real(dp) :: a, c, s
        r = 0; z = 0
        if (.not. allocated(eq%rbc)) return
        do m = 0, ubound(eq%rbc, 2)
            do n = lbound(eq%rbc, 1), ubound(eq%rbc, 1)
                if (eq%rbc(n, m) == 0 .and. eq%rbs(n, m) == 0 .and. eq%zbc(n, m) == 0 &
                    .and. eq%zbs(n, m) == 0) cycle
                a = m*theta - n*eq%nfp*phi
                c = cos(a); s = sin(a)
                r = r + eq%rbc(n, m)*c + eq%rbs(n, m)*s
                z = z + eq%zbc(n, m)*c + eq%zbs(n, m)*s
            end do
        end do
    end subroutine

    real(dp) function boundary_area(eq) result(area)
        ! signed area of the phi = 0 section, positive counterclockwise
        type(imported_equilibrium_t), intent(in) :: eq
        integer, parameter :: nt = 512
        real(dp) :: r(nt), z(nt)
        integer :: i
        do i = 1, nt
            call boundary_point(eq, 2*pi*(i - 1)/nt, 0.0_dp, r(i), z(i))
        end do
        area = 0.5_dp*sum(r*cshift(z, 1) - cshift(r, 1)*z)
    end function

    subroutine orientation_flip(eq)
        ! theta -> -theta: (n, m) -> (-n, m) for m > 0, sine terms change
        ! sign; iota changes sign, the current is physical and keeps it
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), allocatable, dimension(:, :) :: rbc, rbs, zbc, zbs
        integer :: n, nmax
        nmax = ubound(eq%rbc, 1)
        rbc = eq%rbc; rbs = eq%rbs; zbc = eq%zbc; zbs = eq%zbs
        do n = -nmax, nmax
            eq%rbc(n, 1:) = rbc(-n, 1:); eq%rbs(n, 1:) = -rbs(-n, 1:)
            eq%zbc(n, 1:) = zbc(-n, 1:); eq%zbs(n, 1:) = -zbs(-n, 1:)
        end do
        call scale_profile(eq%iota, -1.0_dp)
    end subroutine

    subroutine scale_profile(p, factor)
        type(import_profile_t), intent(inout) :: p
        real(dp), intent(in) :: factor
        if (allocated(p%coefficients)) p%coefficients = factor*p%coefficients
        if (allocated(p%values)) p%values = factor*p%values
        if (allocated(p%slopes)) p%slopes = factor*p%slopes
    end subroutine

    function power_series_profile(coefficients) result(p)
        real(dp), intent(in) :: coefficients(:)
        type(import_profile_t) :: p
        p%kind = 'power_series'
        p%coefficients = coefficients
    end function

    function hermite_profile(knots, values, slopes) result(p)
        real(dp), intent(in) :: knots(:), values(:), slopes(:)
        type(import_profile_t) :: p
        p%kind = 'hermite'
        p%knots = knots; p%values = values; p%slopes = slopes
    end function

    function spline_profile(x, y, ok) result(p)
        ! not-a-knot cubic spline through (x_j, y_j), held exactly by its
        ! nodal values and slopes
        real(dp), intent(in) :: x(:), y(:)
        logical, intent(out) :: ok
        type(import_profile_t) :: p
        type(cubic_spline_t) :: sp
        integer :: status
        call cubic_spline_init(x, reshape(y, [size(y), 1]), sp, status)
        ok = status == 0
        if (.not. ok) return
        p = hermite_profile(x, y, sp%slope(:, 1))
    end function

    pure subroutine hermite_value(knots, values, slopes, x, value, slope)
        real(dp), intent(in) :: knots(:), values(:), slopes(:), x
        real(dp), intent(out) :: value
        real(dp), intent(out), optional :: slope
        integer :: k
        real(dp) :: h, t, h00, h10, h01, h11
        k = 1
        do while (k < size(knots) - 1)
            if (x < knots(k + 1)) exit
            k = k + 1
        end do
        h = knots(k + 1) - knots(k)
        t = (x - knots(k))/h
        h00 = (1 + 2*t)*(1 - t)**2; h10 = t*(1 - t)**2
        h01 = t**2*(3 - 2*t); h11 = t**2*(t - 1)
        value = h00*values(k) + h*h10*slopes(k) + h01*values(k + 1) + h*h11*slopes(k + 1)
        if (present(slope)) slope = (6*t**2 - 6*t)*(values(k) - values(k + 1))/h &
            + (3*t**2 - 4*t + 1)*slopes(k) + (3*t**2 - 2*t)*slopes(k + 1)
    end subroutine

    pure real(dp) function profile_value(p, s) result(value)
        type(import_profile_t), intent(in) :: p
        real(dp), intent(in) :: s
        integer :: k
        value = 0
        select case (p%kind)
        case ('power_series')
            do k = size(p%coefficients), 1, -1
                value = value*s + p%coefficients(k)
            end do
        case ('hermite')
            call hermite_value(p%knots, p%values, p%slopes, s, value)
        end select
    end function

    subroutine import_normalize(eq)
        ! intermediate convention -> nested_equilibrium's
        type(imported_equilibrium_t), intent(inout) :: eq
        integer :: n
        if (eq%status /= IMPORT_OK) return
        if (.not. (abs(eq%phiedge) > 0 .and. ieee_is_finite(eq%phiedge))) then
            call import_fail(eq, 'the toroidal flux is zero or not finite')
            return
        end if
        if (.not. allocated(eq%rbc)) then
            call import_fail(eq, 'the input holds no boundary')
            return
        end if
        if (boundary_area(eq) < 0) then
            call orientation_flip(eq)
            call import_note(eq, 'theta reversed to run counterclockwise in (R, Z); iota changes sign')
        end if
        call scale_profile(eq%current, -sign(1.0_dp, eq%phiedge))
        eq%phiedge = abs(eq%phiedge)
        ! m = 0 terms on n >= 0 only
        do n = 1, ubound(eq%rbc, 1)
            eq%rbc(n, 0) = eq%rbc(n, 0) + eq%rbc(-n, 0); eq%rbs(n, 0) = eq%rbs(n, 0) - eq%rbs(-n, 0)
            eq%zbc(n, 0) = eq%zbc(n, 0) + eq%zbc(-n, 0); eq%zbs(n, 0) = eq%zbs(n, 0) - eq%zbs(-n, 0)
            eq%rbc(-n, 0) = 0; eq%rbs(-n, 0) = 0; eq%zbc(-n, 0) = 0; eq%zbs(-n, 0) = 0
        end do
        eq%rbs(0, 0) = 0; eq%zbs(0, 0) = 0
        call truncate_boundary(eq)
        call fit_profile(eq, eq%iota, 'iota')
        call fit_profile(eq, eq%pressure, 'pressure')
        call fit_profile(eq, eq%current, 'current')
    end subroutine

    subroutine truncate_boundary(eq)
        ! The smallest ntor and mpol whose dropped modes sum to at most half the
        ! geometry tolerance each: the sum of the dropped |c| bounds the
        ! displacement of the boundary, so the spectrum's own decay sets the
        ! resolution. Single modes below 1e-12 of the minor radius go too.
        type(imported_equilibrium_t), intent(inout) :: eq
        real(dp), allocatable :: by_n(:), by_m(:)
        real(dp) :: radius, dropped, tiny, budget
        integer :: m, n, mmax, nmax, ntor, mpol
        call mode_bounds(eq, mmax, nmax)
        radius = sqrt(sum(eq%rbc(:, 1:)**2 + eq%rbs(:, 1:)**2 + eq%zbc(:, 1:)**2 &
            + eq%zbs(:, 1:)**2)/2)
        tiny = 1e-12_dp*radius
        budget = eq%geometry_tolerance*radius/2
        allocate (by_n(0:nmax), by_m(0:mmax))
        by_n = 0; by_m = 0
        do m = 0, mmax
            do n = -nmax, nmax
                associate (c => abs(eq%rbc(n, m)) + abs(eq%rbs(n, m)) + abs(eq%zbc(n, m)) &
                        + abs(eq%zbs(n, m)))
                    by_n(abs(n)) = by_n(abs(n)) + c
                    by_m(m) = by_m(m) + c
                end associate
            end do
        end do
        ntor = nmax
        do while (ntor > 0)
            if (sum(by_n(ntor:)) > budget) exit
            ntor = ntor - 1
        end do
        mpol = mmax
        do while (mpol > 1)
            if (sum(by_m(mpol:)) > budget) exit
            mpol = mpol - 1
        end do
        dropped = 0
        do m = 0, mmax
            do n = -nmax, nmax
                if (m <= mpol .and. abs(n) <= ntor) then
                    call drop(eq%rbc(n, m), tiny); call drop(eq%rbs(n, m), tiny)
                    call drop(eq%zbc(n, m), tiny); call drop(eq%zbs(n, m), tiny)
                else
                    call drop(eq%rbc(n, m), huge(1.0_dp)); call drop(eq%rbs(n, m), huge(1.0_dp))
                    call drop(eq%zbc(n, m), huge(1.0_dp)); call drop(eq%zbs(n, m), huge(1.0_dp))
                end if
            end do
        end do
        eq%geometry_error = eq%geometry_error + dropped
        call import_note(eq, 'boundary held by m <= '//int_text(mpol)//', |n| <= '//int_text(ntor)// &
            ', dropped |c| sum '//real_text(dropped)//' m')
        if (ntor > MILCA_NMAX .or. mpol > MILCA_MMAX) call import_refuse(eq, 'the boundary needs '// &
            int_text(mpol)//' poloidal and '//int_text(ntor)//' toroidal modes, beyond '// &
            'nested_equilibrium''s namelist storage')
    contains
        subroutine drop(c, limit)
            real(dp), intent(inout) :: c
            real(dp), intent(in) :: limit
            if (abs(c) < limit) then
                dropped = dropped + abs(c)
                c = 0
            end if
        end subroutine
    end subroutine

    function int_text(k) result(text)
        integer, intent(in) :: k
        character(len=:), allocatable :: text
        character(16) :: buffer
        write (buffer, '(i0)') k
        text = trim(buffer)
    end function

    subroutine fit_profile(eq, p, name)
        ! power series beyond s^32 become Hermite profiles; the error of that
        ! representation is noted
        type(imported_equilibrium_t), intent(inout) :: eq
        type(import_profile_t), intent(inout) :: p
        character(*), intent(in) :: name
        real(dp), allocatable :: knots(:), values(:), slopes(:)
        real(dp) :: error, scale, h, x
        integer :: k, j
        if (p%kind == 'power_series') then
            k = size(p%coefficients)
            do while (k > 1)
                if (p%coefficients(k) /= 0) exit
                k = k - 1
            end do
            p%coefficients = p%coefficients(:max(k, 1))
            if (size(p%coefficients) <= MILCA_KMAX + 1) return
            allocate (knots(257), values(257), slopes(257))
            do j = 1, 257
                knots(j) = (j - 1)/256.0_dp
                values(j) = profile_value(p, knots(j))
                slopes(j) = 0
                do k = size(p%coefficients), 2, -1
                    slopes(j) = slopes(j)*knots(j) + (k - 1)*p%coefficients(k)
                end do
            end do
            error = 0; scale = maxval(abs(values))
            do j = 1, 256
                h = knots(j + 1) - knots(j)
                do k = 1, 3
                    x = knots(j) + k*h/4
                    error = max(error, abs(profile_value(p, x) - hermite_eval(x)))
                end do
            end do
            p = hermite_profile(knots, values, slopes)
            call import_note(eq, trim(name)//': power series of degree above 32 as Hermite on 257 knots, '// &
                'relative error '//real_text(error/max(scale, tiny(1.0_dp))))
        end if
        if (p%kind == 'hermite') then
            if (size(p%knots) > KNOT_MAX) call import_refuse(eq, trim(name)//' needs more than 1024 Hermite knots')
        end if
    contains
        real(dp) function hermite_eval(s)
            real(dp), intent(in) :: s
            call hermite_value(knots, values, slopes, s, hermite_eval)
        end function
    end subroutine

    function real_text(x) result(text)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: text
        character(32) :: buffer
        write (buffer, '(es10.3)') x
        text = trim(adjustl(buffer))
    end function

    subroutine import_write(eq, path, ok)
        ! nested_input namelist; notes follow the group, where the read stops
        type(imported_equilibrium_t), intent(in) :: eq
        character(*), intent(in) :: path
        logical, intent(out) :: ok
        real(dp), allocatable :: knots(:)
        type(import_profile_t) :: shared(3)
        character(8), parameter :: names(3) = ['iota    ', 'pressure', 'current ']
        integer :: unit, ios, m, n, mpol, ntor, k
        real(dp) :: start
        ok = .false.
        if (eq%status /= IMPORT_OK) return
        mpol = 1; ntor = 0
        do m = 0, min(MILCA_MMAX, ubound(eq%rbc, 2))
            do n = -min(MILCA_NMAX, ubound(eq%rbc, 1)), min(MILCA_NMAX, ubound(eq%rbc, 1))
                if (any([eq%rbc(n, m), eq%rbs(n, m), eq%zbc(n, m), eq%zbs(n, m)] /= 0)) then
                    mpol = max(mpol, m); ntor = max(ntor, abs(n))
                end if
            end do
        end do
        shared = [eq%iota, eq%pressure, eq%current]
        ! one knot set for every Hermite profile: inserting a knot into a
        ! piecewise cubic leaves it unchanged
        allocate (knots(0))
        do k = 1, 3
            if (shared(k)%kind == 'hermite') knots = merge_knots(knots, shared(k)%knots)
        end do
        if (size(knots) > KNOT_MAX) return
        do k = 1, 3
            if (shared(k)%kind == 'hermite') call reknot(shared(k), knots)
        end do
        open (newunit=unit, file=path, status='replace', action='write', iostat=ios)
        if (ios /= 0) return
        write (unit, '(a)') '&nested_input'
        write (unit, '(a,i0,a,i0,a,i0,a,es24.16e3,a)') '  nfp = ', eq%nfp, ', mpol = ', mpol, &
            ', ntor = ', ntor, ', phiedge = ', eq%phiedge, ','
        do m = 0, mpol
            do n = -ntor, ntor
                call coefficient('rbc', eq%rbc(n, m))
                call coefficient('rbs', eq%rbs(n, m))
                call coefficient('zbc', eq%zbc(n, m))
                call coefficient('zbs', eq%zbs(n, m))
            end do
        end do
        if (size(knots) > 0) then
            write (unit, '(a,i0,a)') '  n_knots = ', size(knots), ','
            call list('profile_knots', knots)
        end if
        do k = 1, 3
            select case (shared(k)%kind)
            case ('power_series')
                do m = 1, size(shared(k)%coefficients)
                    if (shared(k)%coefficients(m) /= 0) write (unit, '(a,a,i0,a,es24.16e3,a)') &
                        '  ', trim(names(k))//'(', m - 1, ') = ', shared(k)%coefficients(m), ','
                end do
            case ('hermite')
                write (unit, '(a)') '  '//trim(names(k))//"_type = 'hermite',"
                call list(trim(names(k))//'_values', shared(k)%values)
                call list(trim(names(k))//'_slopes', shared(k)%slopes)
            end select
        end do
        if (eq%prescribe_current) then
            write (unit, '(a)') '  prescribe_current = .true.,'
            if (eq%iota%kind == 'none') then
                ! a start of either sign: the current matching may reverse it
                start = 0.3_dp
                if (eq%current%kind /= 'none') then
                    if (profile_value(eq%current, 1.0_dp) < 0) start = -0.3_dp
                end if
                write (unit, '(a,es24.16e3,a)') '  iota(0) = ', start, ','
            end if
        end if
        ! a foreign boundary angle is not a field-line angle
        write (unit, '(a)') '  angle_mpol = 4, max_iterations = 400,'
        write (unit, '(a)') '/'
        write (unit, '(a)') 'imported from '//trim(eq%format)//' input '//trim(eq%source)
        write (unit, '(a)') 'largest boundary displacement left out: '//real_text(eq%geometry_error)//' m'
        do k = 1, eq%n_notes
            write (unit, '(a)') trim(eq%notes(k))
        end do
        close (unit)
        ok = .true.
    contains
        subroutine coefficient(name, c)
            character(*), intent(in) :: name
            real(dp), intent(in) :: c
            if (c == 0) return
            write (unit, '(a,a,a,i0,a,i0,a,es24.16e3,a)') '  ', name, '(', n, ',', m, ') = ', c, ','
        end subroutine
        subroutine list(name, values)
            character(*), intent(in) :: name
            real(dp), intent(in) :: values(:)
            integer :: j, q
            write (unit, '(a)') '  '//name//' ='
            do j = 1, size(values), 4
                write (unit, '(4x,4(es24.16e3,a))') (values(q), ',', q=j, min(j + 3, size(values)))
            end do
        end subroutine
    end subroutine

    pure function merge_knots(a, b) result(out)
        real(dp), intent(in) :: a(:), b(:)
        real(dp), allocatable :: out(:), all(:)
        integer :: i, j
        real(dp) :: t
        all = [a, b]
        do i = 2, size(all)
            t = all(i); j = i - 1
            do while (j >= 1)
                if (all(j) <= t) exit
                all(j + 1) = all(j); j = j - 1
            end do
            all(j + 1) = t
        end do
        allocate (out(0))
        do i = 1, size(all)
            if (size(out) > 0) then
                if (all(i) - out(size(out)) <= 1e-12_dp) cycle
            end if
            out = [out, all(i)]
        end do
    end function

    subroutine reknot(p, knots)
        type(import_profile_t), intent(inout) :: p
        real(dp), intent(in) :: knots(:)
        real(dp) :: values(size(knots)), slopes(size(knots))
        integer :: j
        do j = 1, size(knots)
            call hermite_value(p%knots, p%values, p%slopes, knots(j), values(j), slopes(j))
        end do
        p = hermite_profile(knots, values, slopes)
    end subroutine
end module milca_equilibrium_import
