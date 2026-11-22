module milca_import_vmec
    ! VMEC inputs to the intermediate form of milca_equilibrium_import, which
    ! is VMEC's own convention: the INDATA namelist with power-series
    ! profiles, its &VSEQ sequence files, VMEX's &MIRROR group (refused) and
    ! wout NetCDF files, whose boundary is the last surface and whose iota
    ! and pressure on the uniform s grid are held by the not-a-knot cubic
    ! spline that GVEC builds from them (src/vmec/vmec.F90).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_import_text, only: kv_store_t, read_namelist_group, namelist_groups, lower
    use milca_equilibrium_import, only: imported_equilibrium_t, import_note, import_refuse, &
        import_fail, IMPORT_OK, power_series_profile, spline_profile, boundary_reserve, &
        import_profile_t
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, nf90_inquire_variable, &
        nf90_inquire_dimension, NF90_NOWRITE, NF90_NOERR
    implicit none
    private
    public :: import_vmec, import_wout

contains

    subroutine import_vmec(path, eq)
        ! a namelist file: &INDATA, &VSEQ or &MIRROR
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        character(32), allocatable :: groups(:)
        logical :: ok
        eq%source = path
        eq%format = 'VMEC'
        call namelist_groups(path, groups, ok)
        if (.not. ok) then
            call import_fail(eq, 'cannot read '//trim(path))
        else if (any(groups == 'indata')) then
            call import_indata(path, eq)
        else if (any(groups == 'vseq')) then
            call import_sequence(path, eq)
        else if (any(groups == 'mirror')) then
            eq%format = 'VMEX mirror'
            call import_refuse(eq, 'VMEX &MIRROR describes an open linear mirror; nested_equilibrium '// &
                'solves toroidal equilibria')
        else
            call import_fail(eq, 'no &INDATA, &VSEQ or &MIRROR group in '//trim(path))
        end if
    end subroutine

    subroutine import_sequence(path, eq)
        ! VMEC runs input.<ext> for every EXTENSION in turn; the last one is
        ! the final equilibrium
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        type(kv_store_t) :: store
        type(imported_equilibrium_t), allocatable :: member
        character(len=:), allocatable :: list, dir, name
        integer :: i, start
        logical :: ok
        call read_namelist_group(path, 'vseq', store, ok)
        list = store%text('extension', '')
        dir = path(:index(path, '/', back=.true.))
        name = ''
        i = 1
        do while (i <= len(list))
            if (list(i:i) == '"' .or. list(i:i) == "'") then
                start = i + 1
                i = index(list(start:), list(i:i)) + start - 1
                if (i < start) exit
                name = list(start:i - 1)
            end if
            i = i + 1
        end do
        if (len(name) == 0) then
            call import_fail(eq, 'VSEQ without an EXTENSION list')
            return
        end if
        allocate (member)
        call import_indata(dir//'input.'//name, member)
        member%format = 'VMEC sequence member'
        call import_note(member, 'last member input.'//name//' of the VMEC sequence '//trim(path))
        eq = member
        if (eq%status /= IMPORT_OK) eq%reason = 'sequence member input.'//name//': '//trim(eq%reason)
    end subroutine

    subroutine import_indata(path, eq)
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        type(kv_store_t) :: p
        integer, allocatable :: n(:), m(:)
        real(dp), allocatable :: v(:), am(:), ai(:), ac(:), enclosed(:)
        character(4), parameter :: fields(4) = ['rbc ', 'rbs ', 'zbc ', 'zbs ']
        character(5), parameter :: profiles(3) = ['pmass', 'piota', 'pcurr']
        character(len=:), allocatable :: kind
        integer :: f, k
        logical :: ok
        eq%source = path
        eq%format = 'VMEC'
        call read_namelist_group(path, 'indata', p, ok)
        if (.not. ok) then
            call import_fail(eq, 'no &INDATA group in '//trim(path))
            return
        end if
        ! VMEC defaults LFREEB to true and runs free-boundary whenever a grid is named
        if (p%logical('lfreeb', .true.) .and. len_trim(p%text('mgrid_file', '')) > 0 .and. &
            lower(p%text('mgrid_file', 'none')) /= 'none') then
            ! the boundary is the initial guess in the field of the coil groups
            eq%free_boundary = .true.
            eq%mgrid_file = p%text('mgrid_file', '')
            call read_extcur(eq%extcur)
            if (size(eq%extcur) == 0 .or. all(eq%extcur == 0)) then
                call import_refuse(eq, 'free-boundary VMEC input without EXTCUR')
                return
            end if
        end if
        if (p%logical('lrfp', .false.)) then
            call import_refuse(eq, 'LRFP: iota of the poloidal flux is not imported')
            return
        end if
        do k = 1, 3
            kind = lower(p%text(trim(profiles(k))//'_type', 'power_series'))
            if (kind == 'power_series') cycle
            ! knot splines of pressure and iota: the not-a-knot cubic through
            ! the AUX knots (dense knots make the end condition immaterial)
            if (k <= 2 .and. (kind == 'cubic_spline' .or. kind == 'akima_spline')) cycle
            call import_refuse(eq, 'VMEC profile type '//kind//' of '//profiles(k)// &
                               ' (power_series; cubic or Akima splines for pressure and iota)')
            return
        end do
        if (p%real('bloat', 1.0_dp) /= 1 .or. abs(p%real('spres_ped', 1.0_dp)) < 1) then
            call import_refuse(eq, 'BLOAT or SPRES_PED (power_series without them only)')
            return
        end if
        eq%nfp = p%int('nfp', 1)
        eq%phiedge = p%real('phiedge', 1.0_dp)
        do f = 1, 4
            call p%modes(trim(fields(f)), n, m, v)
            do k = 1, size(v)
                if (m(k) < 0) cycle
                call boundary_reserve(eq, m(k), abs(n(k)))
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
        call boundary_reserve(eq, 1, 0)
        if (.not. p%logical('lasym', .false.)) then
            eq%rbs = 0; eq%zbc = 0
        end if
        call coefficients('am', am)
        call coefficients('ai', ai)
        call coefficients('ac', ac)
        eq%pressure = power_series_profile(p%real('pres_scale', 1.0_dp)*am)
        if (lower(p%text('pmass_type', 'power_series')) /= 'power_series') then
            call aux_spline('am', p%real('pres_scale', 1.0_dp), eq%pressure, ok)
            if (.not. ok) return
        end if
        if (p%int('ncurr', 0) == 1) then
            ! AC is dI/ds; I(s) is scaled to CURTOR at the edge
            enclosed = [0.0_dp, ac/[(real(k, dp), k=1, size(ac))]]
            if (sum(enclosed) /= 0) then
                enclosed = p%real('curtor', 0.0_dp)*enclosed/sum(enclosed)
            else
                enclosed = [0.0_dp]
            end if
            eq%current = power_series_profile(enclosed)
            eq%prescribe_current = .true.
        else
            eq%iota = power_series_profile(ai)
            if (lower(p%text('piota_type', 'power_series')) /= 'power_series') then
                call aux_spline('ai', 1.0_dp, eq%iota, ok)
                if (.not. ok) return
            end if
        end if
    contains
        subroutine aux_spline(name, scale, profile, ok)
            ! NAME_AUX_S, NAME_AUX_F as a spline profile; the knots end where
            ! VMEC's do (the first non-increasing knot after the first)
            character(*), intent(in) :: name
            real(dp), intent(in) :: scale
            type(import_profile_t), intent(out) :: profile
            logical, intent(out) :: ok
            real(dp), allocatable :: ks(:), kf(:)
            logical :: parsed
            integer :: nk
            call p%reals(name//'_aux_s', ks, parsed)
            call p%reals(name//'_aux_f', kf, parsed)
            nk = 1
            do while (nk < min(size(ks), size(kf)))
                if (ks(nk + 1) <= ks(nk)) exit
                nk = nk + 1
            end do
            ok = nk >= 4
            if (.not. ok) then
                call import_refuse(eq, 'spline profile '//name//' with fewer than four knots')
                return
            end if
            profile = spline_profile(ks(:nk), scale*kf(:nk), ok)
            if (.not. ok) call import_refuse(eq, 'spline profile '//name//' not admitted')
        end subroutine

        subroutine read_extcur(c)
            ! EXTCUR as a list, EXTCUR(k) = v entries, or both (indexed wins)
            real(dp), allocatable, intent(out) :: c(:)
            real(dp), allocatable :: listed(:)
            real(dp) :: v
            logical :: parsed, found
            integer :: k, last
            character(16) :: key
            call p%reals('extcur', listed, parsed)
            if (.not. parsed) then
                deallocate (listed)
                allocate (listed(0))
            end if
            allocate (c(512))
            c = 0
            c(:size(listed)) = listed
            last = size(listed)
            do k = 1, 512
                write (key, '(a,i0,a)') 'extcur(', k, ')'
                if (.not. p%has(trim(key))) cycle
                v = p%real(trim(key), 0.0_dp, found)
                if (found) then
                    c(k) = v
                    last = max(last, k)
                end if
            end do
            c = c(:last)
        end subroutine

        subroutine coefficients(name, c)
            character(*), intent(in) :: name
            real(dp), allocatable, intent(out) :: c(:)
            logical :: parsed
            call p%reals(name, c, parsed)
            if (.not. parsed .or. size(c) == 0) c = [0.0_dp]
        end subroutine
    end subroutine

    subroutine import_wout(path, eq)
        character(*), intent(in) :: path
        type(imported_equilibrium_t), intent(inout) :: eq
        integer :: ncid, ns, nfp, mnmax, lasym, j, k, m, n
        real(dp), allocatable :: xm(:), xn(:), rmnc(:, :), zmns(:, :), rmns(:, :), zmnc(:, :)
        real(dp), allocatable :: iotaf(:), presf(:), phi(:), s(:)
        logical :: ok
        eq%source = path
        eq%format = 'VMEC wout'
        if (nf90_open(path, NF90_NOWRITE, ncid) /= NF90_NOERR) then
            call import_fail(eq, 'cannot open the wout file '//trim(path))
            return
        end if
        ok = scalar('ns', ns) .and. scalar('nfp', nfp) .and. scalar('mnmax', mnmax) &
            .and. scalar('lasym__logical__', lasym)
        ok = ok .and. vector('xm', xm) .and. vector('xn', xn) .and. vector('iotaf', iotaf) &
            .and. vector('presf', presf) .and. vector('phi', phi)
        ok = ok .and. matrix('rmnc', rmnc) .and. matrix('zmns', zmns)
        if (ok .and. lasym /= 0) ok = matrix('rmns', rmns) .and. matrix('zmnc', zmnc)
        j = nf90_close(ncid)
        if (.not. ok) then
            call import_fail(eq, 'incomplete wout file '//trim(path))
            return
        end if
        eq%nfp = nfp
        eq%phiedge = phi(ns)
        do k = 1, mnmax
            m = nint(xm(k)); n = nint(xn(k))/nfp
            call boundary_reserve(eq, m, abs(n))
            eq%rbc(n, m) = rmnc(k, ns); eq%zbs(n, m) = zmns(k, ns)
            if (lasym /= 0) then
                eq%rbs(n, m) = rmns(k, ns); eq%zbc(n, m) = zmnc(k, ns)
            end if
        end do
        s = [((j - 1.0_dp)/(ns - 1), j=1, ns)]
        eq%iota = spline_profile(s, iotaf, ok)
        eq%pressure = spline_profile(s, presf, ok)
        if (.not. ok) call import_fail(eq, 'the wout file needs at least 4 surfaces')
        call import_note(eq, 'iota and pressure: not-a-knot cubic splines through the wout full-grid '// &
            'values, as in GVEC''s VMEC initialisation')
    contains
        logical function scalar(name, value)
            character(*), intent(in) :: name
            integer, intent(out) :: value
            integer :: id
            scalar = nf90_inq_varid(ncid, name, id) == NF90_NOERR
            if (scalar) scalar = nf90_get_var(ncid, id, value) == NF90_NOERR
        end function
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
        end function
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
            if (matrix .and. size(values, 1) /= mnmax) values = transpose(values)
        end function
    end subroutine
end module milca_import_vmec
