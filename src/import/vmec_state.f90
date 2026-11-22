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
    end type vmec_state_t

    public :: vmec_state_read, vmec_state_eval, vmec_state_invert

contains

    subroutine vmec_state_read(path, st, status, reason)
        !! status 0, or 1 with reason.
        character(*), intent(in) :: path
        type(vmec_state_t), intent(out) :: st
        integer, intent(out) :: status
        character(*), intent(out) :: reason
        integer :: ncid, lasym, j, k
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), lmns(:, :), phi(:)
        logical :: ok

        status = 1
        reason = ''
        if (nf90_open(path, NF90_NOWRITE, ncid) /= NF90_NOERR) then
            reason = 'cannot open '//trim(path)
            return
        end if
        ok = scalar('ns', st%ns) .and. scalar('nfp', st%nfp) .and. scalar('mnmax', st%mnmax) &
            .and. scalar('lasym__logical__', lasym)
        ok = ok .and. vector('xm', st%xm) .and. vector('xn', st%xn) .and. &
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
        st%phiedge = phi(st%ns)
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
