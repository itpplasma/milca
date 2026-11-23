module milca_vmec_state
    !! The interior state of a VMEC wout file as a candidate for MILCA's own
    !! models (residual check, start of a re-solve): R, Z on the full grid,
    !! lambda on the half grid, evaluated at any (s, theta, zeta) and inverted
    !! at a point (R, Z, zeta) -> (s, theta). Stellarator-symmetric states
    !! only (lasym refused).
    !!
    !! Conventions of the wout file: s is the normalized toroidal flux on the
    !! full grid s_j = (j - 1)/(ns - 1); angles m theta - n zeta with xn
    !! including nfp and zeta the cylindrical angle; the straight-field-line
    !! angle is theta* = theta + lambda. Radial interpolation keeps the
    !! regularity at the axis: a coefficient of poloidal mode m >= 1 is
    !! rho^m g(s), rho = sqrt(s), g linear in s between grid surfaces
    !! (extrapolated linearly at the axis and beyond the last surface).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_inquire_variable, &
        nf90_inquire_dimension, NF90_NOWRITE, NF90_NOERR
    implicit none
    private

    type, public :: vmec_state_t
        integer :: ns = 0, nfp = 1, mnmax = 0
        real(dp) :: phiedge = 0
        real(dp), allocatable :: xm(:), xn(:)
        integer, allocatable :: m(:)
        !! g of R, Z (full grid) and of lambda (half grid; index j = 1 holds
        !! the extrapolation to s = 0), (mode, surface)
        real(dp), allocatable :: gr(:, :), gz(:, :), gl(:, :)
        real(dp), allocatable :: sf(:), sh(:), iotaf(:), presf(:)
        !! vmec_state_jet: the coefficients over rho^parity, parity = mod(m, 2)
        !! (smooth functions of s), and the second s-derivatives at the grid
        !! points of the cubic splines through them
        real(dp), allocatable :: hr(:, :), hz(:, :), hl(:, :), mr(:, :), mz(:, :), ml(:, :)
        integer, allocatable :: parity(:)
    end type vmec_state_t

    public :: vmec_state_read, vmec_state_set, vmec_state_eval, vmec_state_invert, vmec_state_jet

contains

    subroutine vmec_state_read(path, st, status, reason)
        !! status 0, or 1 with reason.
        character(*), intent(in) :: path
        type(vmec_state_t), intent(out) :: st
        integer, intent(out) :: status
        character(*), intent(out) :: reason
        integer :: ncid, lasym, j, nfp
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), lmns(:, :), phi(:), xm(:), xn(:)
        logical :: ok

        status = 1
        reason = ''
        if (nf90_open(path, NF90_NOWRITE, ncid) /= NF90_NOERR) then
            reason = 'cannot open '//trim(path)
            return
        end if
        ok = scalar('ns', st%ns) .and. scalar('nfp', nfp) .and. scalar('mnmax', st%mnmax) &
            .and. scalar('lasym__logical__', lasym)
        ok = ok .and. vector('xm', xm) .and. vector('xn', xn) .and. &
            vector('iotaf', st%iotaf) .and. vector('presf', st%presf) .and. vector('phi', phi)
        ok = ok .and. matrix('rmnc', rmnc) .and. matrix('zmns', zmns) .and. matrix('lmns', lmns)
        j = nf90_close(ncid)
        if (.not. ok) then
            reason = 'incomplete wout file '//trim(path)
            return
        end if
        if (lasym /= 0) then
            reason = 'non-symmetric (lasym) wout states are not supported'
            return
        end if
        if (st%ns < 4) then
            reason = 'the wout file needs at least 4 surfaces'
            return
        end if
        call vmec_state_set(st, nfp, xm, xn, rmnc, zmns, lmns, phi(st%ns))
        status = 0
    contains
        logical function scalar(name, value)
            character(*), intent(in) :: name
            integer, intent(out) :: value
            integer :: id
            scalar = nf90_inq_varid(ncid, name, id) == NF90_NOERR
            if (scalar) scalar = nf90_get_var(ncid, id, value) == NF90_NOERR
        end function scalar
        logical function vector(name, values)
            character(*), intent(in) :: name
            real(dp), allocatable, intent(out) :: values(:)
            integer :: id, dims(1), length
            vector = nf90_inq_varid(ncid, name, id) == NF90_NOERR
            if (.not. vector) return
            vector = nf90_inquire_variable(ncid, id, dimids=dims) == NF90_NOERR
            if (.not. vector) return
            vector = nf90_inquire_dimension(ncid, dims(1), len=length) == NF90_NOERR
            allocate (values(length))
            if (vector) vector = nf90_get_var(ncid, id, values) == NF90_NOERR
        end function vector
        logical function matrix(name, values)
            ! (mode, surface)
            character(*), intent(in) :: name
            real(dp), allocatable, intent(out) :: values(:, :)
            integer :: id, dims(2), length(2)
            matrix = nf90_inq_varid(ncid, name, id) == NF90_NOERR
            if (.not. matrix) return
            matrix = nf90_inquire_variable(ncid, id, dimids=dims) == NF90_NOERR
            if (.not. matrix) return
            matrix = nf90_inquire_dimension(ncid, dims(1), len=length(1)) == NF90_NOERR .and. &
                nf90_inquire_dimension(ncid, dims(2), len=length(2)) == NF90_NOERR
            if (.not. matrix) return
            allocate (values(length(1), length(2)))
            matrix = nf90_get_var(ncid, id, values) == NF90_NOERR
            if (matrix .and. size(values, 1) /= st%mnmax) values = transpose(values)
        end function matrix
    end subroutine vmec_state_read

    subroutine vmec_state_set(st, nfp, xm, xn, rmnc, zmns, lmns, phiedge)
        !! The state from its Fourier coefficients (mode, surface): rmnc, zmns
        !! on the full grid, lmns on the half grid (column 1 unused), as in a
        !! wout file. iotaf and presf are left to the caller.
        type(vmec_state_t), intent(inout) :: st
        integer, intent(in) :: nfp
        real(dp), intent(in) :: xm(:), xn(:), rmnc(:, :), zmns(:, :), lmns(:, :), phiedge
        integer :: j, k, last

        st%ns = size(rmnc, 2)
        st%mnmax = size(rmnc, 1)
        st%nfp = nfp
        st%xm = xm
        st%xn = xn
        st%phiedge = phiedge
        if (allocated(st%gr)) deallocate (st%gr, st%gz, st%gl, st%mr, st%mz, st%ml, st%hr, &
                                          st%hz, st%hl)
        st%m = nint(st%xm)
        st%sf = [((j - 1.0_dp)/(st%ns - 1), j=1, st%ns)]
        st%sh = [((j - 1.5_dp)/(st%ns - 1), j=1, st%ns)]
        st%sh(1) = 0
        allocate (st%gr(st%mnmax, st%ns), st%gz(st%mnmax, st%ns), st%gl(st%mnmax, st%ns))
        do k = 1, st%mnmax
            do j = 2, st%ns
                st%gr(k, j) = rmnc(k, j)/sqrt(st%sf(j))**st%m(k)
                st%gz(k, j) = zmns(k, j)/sqrt(st%sf(j))**st%m(k)
                st%gl(k, j) = lmns(k, j)/sqrt(st%sh(j))**st%m(k)
            end do
            if (st%m(k) == 0) then
                st%gr(k, 1) = rmnc(k, 1)
                st%gz(k, 1) = zmns(k, 1)
            else
                st%gr(k, 1) = extrapolate(st%sf(2:3), st%gr(k, 2:3), 0.0_dp)
                st%gz(k, 1) = extrapolate(st%sf(2:3), st%gz(k, 2:3), 0.0_dp)
            end if
            st%gl(k, 1) = extrapolate(st%sh(2:3), st%gl(k, 2:3), 0.0_dp)
        end do
        allocate (st%mr(st%mnmax, st%ns), st%mz(st%mnmax, st%ns), st%ml(st%mnmax, st%ns), &
                  st%hr(st%mnmax, st%ns), st%hz(st%mnmax, st%ns), st%hl(st%mnmax, st%ns))
        st%parity = mod(st%m, 2)
        ! A four-surface state has only three positive-radius samples.
        last = min(5, st%ns)
        do k = 1, st%mnmax
            do j = 2, st%ns
                st%hr(k, j) = rmnc(k, j)/sqrt(st%sf(j))**st%parity(k)
                st%hz(k, j) = zmns(k, j)/sqrt(st%sf(j))**st%parity(k)
                st%hl(k, j) = lmns(k, j)/sqrt(st%sh(j))**st%parity(k)
            end do
            if (st%parity(k) == 0) then
                st%hr(k, 1) = rmnc(k, 1)
                st%hz(k, 1) = zmns(k, 1)
            else
                st%hr(k, 1) = end_value(st%sf(2:last), st%hr(k, 2:last))
                st%hz(k, 1) = end_value(st%sf(2:last), st%hz(k, 2:last))
            end if
            st%hl(k, 1) = end_value(st%sh(2:last), st%hl(k, 2:last))
            call spline_moments(st%sf, st%hr(k, :), st%mr(k, :))
            call spline_moments(st%sf, st%hz(k, :), st%mz(k, :))
            call spline_moments(st%sh, st%hl(k, :), st%ml(k, :))
        end do
    end subroutine vmec_state_set

    pure real(dp) function extrapolate(x, y, x0)
        real(dp), intent(in) :: x(2), y(2), x0
        extrapolate = y(1) + (y(2) - y(1))*(x0 - x(1))/(x(2) - x(1))
    end function extrapolate

    pure subroutine radial(grid, g, s, val, slope)
        !! g(s) linear between grid points (extrapolated outside) and its
        !! s-derivative, for every mode.
        real(dp), intent(in) :: grid(:), g(:, :), s
        real(dp), intent(out) :: val(:), slope(:)
        integer :: j, n
        real(dp) :: w

        n = size(grid)
        j = min(max(1, int(s*(n - 1)) + 1), n - 1)
        if (grid(2) - grid(1) < 0.75_dp/(n - 1)) then
            ! half grid: grid(1) = 0, grid(j) = (j - 1.5)/(n - 1)
            j = min(max(1, int(s*(n - 1) + 1.5_dp)), n - 1)
        end if
        slope = (g(:, j + 1) - g(:, j))/(grid(j + 1) - grid(j))
        w = s - grid(j)
        val = g(:, j) + w*slope
    end subroutine radial

    pure real(dp) function end_value(x, y) result(v)
        !! Value at the axis of the interpolant through the supplied samples.
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: w
        integer :: i, j

        v = 0
        do i = 1, size(x)
            w = y(i)
            do j = 1, size(x)
                if (j /= i) w = w*(0 - x(j))/(x(i) - x(j))
            end do
            v = v + w
        end do
    end function end_value

    pure real(dp) function end_second(x, y) result(d2)
        !! Second derivative at x(1) of the cubic through four points.
        real(dp), intent(in) :: x(4), y(4)
        real(dp) :: w
        integer :: i, j, k, l

        d2 = 0
        do i = 1, 4
            ! second derivative of the Lagrange polynomial of node i at x(1)
            w = 0
            do j = 1, 4
                if (j == i) cycle
                do k = 1, 4
                    if (k == i .or. k == j) cycle
                    l = 10 - i - j - k
                    w = w + (x(1) - x(l))
                end do
            end do
            do j = 1, 4
                if (j /= i) w = w/(x(i) - x(j))
            end do
            d2 = d2 + w*y(i)
        end do
    end function end_second

    pure subroutine spline_moments(x, y, m)
        !! Second derivatives m of the C2 cubic spline through (x, y), the end
        !! values from the cubics through the first and the last four points.
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: m(:)
        real(dp) :: d(size(x)), c(size(x)), rhs(size(x)), h0, h1, w
        integer :: n, i

        n = size(x)
        m(1) = end_second(x(1:4), y(1:4))
        m(n) = end_second(x(n:n - 3:-1), y(n:n - 3:-1))
        ! Thomas algorithm for the interior moments
        c(1) = 0
        rhs(1) = m(1)
        d(1) = 1
        do i = 2, n - 1
            h0 = x(i) - x(i - 1)
            h1 = x(i + 1) - x(i)
            w = h0/d(i - 1)
            d(i) = 2*(h0 + h1) - w*c(i - 1)
            c(i) = h1
            rhs(i) = 6*((y(i + 1) - y(i))/h1 - (y(i) - y(i - 1))/h0) - w*rhs(i - 1)
        end do
        do i = n - 1, 2, -1
            m(i) = (rhs(i) - c(i)*m(i + 1))/d(i)
        end do
    end subroutine spline_moments

    pure subroutine radial_spline(grid, g, mom, s, val, d1, d2)
        !! The cubic spline of every mode at s (the end cubics extrapolated)
        !! with its first and second s-derivatives.
        real(dp), intent(in) :: grid(:), g(:, :), mom(:, :), s
        real(dp), intent(out) :: val(:), d1(:), d2(:)
        integer :: j, n
        real(dp) :: h, a, b

        n = size(grid)
        j = min(max(1, int(s*(n - 1)) + 1), n - 1)
        if (grid(2) - grid(1) < 0.75_dp/(n - 1)) j = min(max(1, int(s*(n - 1) + 1.5_dp)), n - 1)
        h = grid(j + 1) - grid(j)
        a = (grid(j + 1) - s)/h
        b = (s - grid(j))/h
        val = a*g(:, j) + b*g(:, j + 1) + ((a**3 - a)*mom(:, j) + (b**3 - b)*mom(:, j + 1))*h**2/6
        d1 = (g(:, j + 1) - g(:, j))/h + (-(3*a**2 - 1)*mom(:, j) + (3*b**2 - 1)*mom(:, j + 1))*h/6
        d2 = a*mom(:, j) + b*mom(:, j + 1)
    end subroutine radial_spline

    subroutine vmec_state_jet(st, rho, theta, zeta, rj, zj, lj)
        !! 2-jets of R, Z and lambda in (rho, theta, zeta), ordered (f, f_1,
        !! f_2, f_3, f_11, f_12, f_13, f_22, f_23, f_33), so that curl B of
        !! the state can be formed. Radially a coefficient is rho^parity h(s),
        !! parity = mod(m, 2), with h the C2 cubic spline in s through the
        !! grid values (the variables VMEC itself discretises; the rho^m
        !! scaling of vmec_state_eval divides the high modes by numbers that
        !! are tiny near the axis, and the second derivatives of its
        !! piecewise linear g do not exist).
        type(vmec_state_t), intent(in) :: st
        real(dp), intent(in) :: rho, theta, zeta
        real(dp), intent(out) :: rj(10), zj(10), lj(10)
        real(dp), dimension(st%mnmax) :: v, d1, d2, pm, p1, p2, ang, c, sn, g0, g1, g2
        real(dp) :: s
        integer :: pass

        s = rho**2
        pm = rho**st%parity
        p1 = st%parity
        p2 = 0
        ang = st%xm*theta - st%xn*zeta
        c = cos(ang)
        sn = sin(ang)
        do pass = 1, 3
            select case (pass)
            case (1)
                call radial_spline(st%sf, st%hr, st%mr, s, v, d1, d2)
            case (2)
                call radial_spline(st%sf, st%hz, st%mz, s, v, d1, d2)
            case (3)
                call radial_spline(st%sh, st%hl, st%ml, s, v, d1, d2)
            end select
            ! G(rho) = rho^parity h(rho^2) and its rho derivatives
            g0 = pm*v
            g1 = p1*v + pm*d1*2*rho
            g2 = p2*v + 4*rho*p1*d1 + pm*(4*s*d2 + 2*d1)
            if (pass == 1) then
                rj = [sum(g0*c), sum(g1*c), -sum(g0*st%xm*sn), sum(g0*st%xn*sn), sum(g2*c), &
                      -sum(g1*st%xm*sn), sum(g1*st%xn*sn), -sum(g0*st%xm**2*c), &
                      sum(g0*st%xm*st%xn*c), -sum(g0*st%xn**2*c)]
            else
                associate (j => merge(1, 2, pass == 2))
                    if (j == 1) then
                        zj = sine_jet()
                    else
                        lj = sine_jet()
                    end if
                end associate
            end if
        end do

    contains

        function sine_jet() result(j)
            real(dp) :: j(10)
            j = [sum(g0*sn), sum(g1*sn), sum(g0*st%xm*c), -sum(g0*st%xn*c), sum(g2*sn), &
                 sum(g1*st%xm*c), -sum(g1*st%xn*c), -sum(g0*st%xm**2*sn), &
                 sum(g0*st%xm*st%xn*sn), -sum(g0*st%xn**2*sn)]
        end function sine_jet
    end subroutine vmec_state_jet

    subroutine vmec_state_eval(st, rho, theta, zeta, r, z, lam, jac)
        !! R, Z, lambda at (rho = sqrt(s), theta, zeta); jac = d(R, Z)/d(rho, theta).
        type(vmec_state_t), intent(in) :: st
        real(dp), intent(in) :: rho, theta, zeta
        real(dp), intent(out) :: r, z, lam, jac(2, 2)
        real(dp) :: s, vr(st%mnmax), sr(st%mnmax), vz(st%mnmax), sz(st%mnmax), vl(st%mnmax), &
                    sl(st%mnmax), pm(st%mnmax), dpm(st%mnmax), ang(st%mnmax), c(st%mnmax), &
                    sn(st%mnmax)

        s = rho**2
        call radial(st%sf, st%gr, s, vr, sr)
        call radial(st%sf, st%gz, s, vz, sz)
        call radial(st%sh, st%gl, s, vl, sl)
        pm = rho**st%m
        where (st%m > 0)
            dpm = st%m*rho**(st%m - 1)
        elsewhere
            dpm = 0
        end where
        ang = st%xm*theta - st%xn*zeta
        c = cos(ang)
        sn = sin(ang)
        r = sum(pm*vr*c)
        z = sum(pm*vz*sn)
        lam = sum(pm*vl*sn)
        ! d/drho of rho^m g(rho^2) = m rho^(m-1) g + rho^m g' 2 rho
        jac(1, 1) = sum((dpm*vr + pm*sr*2*rho)*c)
        jac(2, 1) = sum((dpm*vz + pm*sz*2*rho)*sn)
        jac(1, 2) = -sum(pm*vr*st%xm*sn)
        jac(2, 2) = sum(pm*vz*st%xm*c)
    end subroutine vmec_state_eval

    subroutine vmec_state_invert(st, rz, zeta, rho, theta, status)
        !! (rho, theta) with (R, Z)(rho, theta, zeta) = rz by damped Newton
        !! from the given (rho, theta); rho may exceed 1 slightly (the last
        !! surface extrapolated). status 0, or 1 when Newton fails.
        type(vmec_state_t), intent(in) :: st
        real(dp), intent(in) :: rz(2), zeta
        real(dp), intent(inout) :: rho, theta
        integer, intent(out) :: status
        real(dp) :: r, z, lam, jac(2, 2), f(2), d(2), det, a, fn, rt, tt, gn, scale
        integer :: it, k

        status = 1
        call vmec_state_eval(st, rho, theta, zeta, r, z, lam, jac)
        f = [r, z] - rz
        fn = norm2(f)
        scale = abs(st%gr(1, st%ns))
        do it = 1, 60
            if (fn < 1e-13_dp*scale) then
                status = 0
                return
            end if
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            if (abs(det) < tiny(1.0_dp)) return
            d = -[jac(2, 2)*f(1) - jac(1, 2)*f(2), -jac(2, 1)*f(1) + jac(1, 1)*f(2)]/det
            a = 1
            do k = 1, 30
                rt = rho + a*d(1)
                tt = theta + a*d(2)
                if (rt < 0) then
                    rt = -rt
                    tt = tt + acos(-1.0_dp)
                end if
                call vmec_state_eval(st, rt, tt, zeta, r, z, lam, jac)
                gn = norm2([r, z] - rz)
                if (gn < (1 - 1e-4_dp*a)*fn .or. gn < 1e-13_dp*scale) exit
                a = a/2
            end do
            if (k > 30) return
            rho = rt
            theta = tt
            f = [r, z] - rz
            fn = gn
        end do
        if (fn < 1e-10_dp*scale) status = 0
    end subroutine vmec_state_invert

end module milca_vmec_state
