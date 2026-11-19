!> MAKEGRID's mgrid NetCDF file: the fields of the coil groups on a
!> cylindrical grid over one field period (variables ir, jz, kp, nfp,
!> nextcur, rmin, rmax, zmin, zmax, mgrid_mode, raw_coil_cur and
!> br_NNN, bp_NNN, bz_NNN of shape (phi, zee, rad) per group NNN).
!> mgrid_mode 'S' stores fields per unit EXTCUR, 'R' fields at the currents
!> raw_coil_cur, so a group scales by EXTCUR/raw_coil_cur. A file without
!> mgrid_mode is read as 'S'.
module milca_mgrid_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_create, nf90_close, nf90_inq_varid, &
        nf90_get_var, &
        nf90_put_var, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_strerror, &
        NF90_NOERR, NF90_NOWRITE, NF90_CLOBBER, NF90_INT, &
        NF90_DOUBLE, NF90_CHAR
    implicit none
    private
    public :: mgrid_t, mgrid_read, mgrid_write

    type :: mgrid_t
        integer :: ir = 0, jz = 0, kp = 0, nfp = 1, nextcur = 0
        real(dp) :: rmin = 0, rmax = 0, zmin = 0, zmax = 0
        character :: mode = 'S'
        real(dp), allocatable :: raw_coil_cur(:)
        !> (B_R, B_phi, B_Z) at (R_i, Z_j, phi_l) for group k: b(:, i, j, l, k)
        real(dp), allocatable :: b(:, :, :, :, :)
    contains
        procedure :: unit_scale => mgrid_unit_scale
    end type

contains

    function mgrid_unit_scale(self) result(scale)
        !! Factor per group that turns the stored field into the field per
        !! unit EXTCUR.
        class(mgrid_t), intent(in) :: self
        real(dp) :: scale(self%nextcur)
        scale = 1
        if (self%mode == 'R' .or. self%mode == 'r') scale = 1/self%raw_coil_cur
    end function

    subroutine mgrid_read(path, grid, status, message)
        character(*), intent(in) :: path
        type(mgrid_t), intent(out) :: grid
        integer, intent(out) :: status
        character(*), intent(out) :: message
        real(dp), allocatable :: component(:, :, :)
        character(3) :: tag
        character(2), parameter :: names(3) = ['br', 'bp', 'bz']
        character(1) :: mode(1)
        integer :: ncid, varid, k, c
        message = ''
        status = nf90_open(path, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            message = trim(path)//': '//trim(nf90_strerror(status))
            return
        end if
        call get_int('ir', grid%ir); if (status /= 0) return
        call get_int('jz', grid%jz); if (status /= 0) return
        call get_int('kp', grid%kp); if (status /= 0) return
        call get_int('nfp', grid%nfp); if (status /= 0) return
        call get_int('nextcur', grid%nextcur); if (status /= 0) return
        call get_real('rmin', grid%rmin); if (status /= 0) return
        call get_real('rmax', grid%rmax); if (status /= 0) return
        call get_real('zmin', grid%zmin); if (status /= 0) return
        call get_real('zmax', grid%zmax); if (status /= 0) return
        if (min(grid%ir, grid%jz, grid%kp, grid%nfp, grid%nextcur) < 1 &
            .or. grid%nextcur > 999) then
            call stop_with('grid dimensions out of range')
            return
        end if
        allocate (grid%raw_coil_cur(grid%nextcur))
        grid%raw_coil_cur = 1
        if (nf90_inq_varid(ncid, 'mgrid_mode', varid) == NF90_NOERR) then
            status = nf90_get_var(ncid, varid, mode)
            if (status /= NF90_NOERR) then
                call stop_with('mgrid_mode')
                return
            end if
            grid%mode = mode(1)
        end if
        if (nf90_inq_varid(ncid, 'raw_coil_cur', varid) == NF90_NOERR) then
            status = nf90_get_var(ncid, varid, grid%raw_coil_cur)
            if (status /= NF90_NOERR) then
                call stop_with('raw_coil_cur')
                return
            end if
        end if
        if ((grid%mode == 'R' .or. grid%mode == 'r') .and. &
            any(.not. abs(grid%raw_coil_cur) > 0)) then
            call stop_with('raw grid with a zero raw_coil_cur')
            return
        end if
        allocate (grid%b(3, grid%ir, grid%jz, grid%kp, grid%nextcur), &
            component(grid%ir, grid%jz, grid%kp))
        do k = 1, grid%nextcur
            write (tag, '(i3.3)') k
            do c = 1, 3
                status = nf90_inq_varid(ncid, names(c)//'_'//tag, varid)
                if (status == NF90_NOERR) status = nf90_get_var(ncid, varid, component)
                if (status /= NF90_NOERR) then
                    call stop_with(names(c)//'_'//tag)
                    return
                end if
                grid%b(c, :, :, :, k) = component
            end do
        end do
        status = nf90_close(ncid)
        if (status /= NF90_NOERR) &
            message = trim(path)//': '//trim(nf90_strerror(status))
    contains
        subroutine get_int(name, value)
            character(*), intent(in) :: name
            integer, intent(out) :: value
            status = nf90_inq_varid(ncid, name, varid)
            if (status == NF90_NOERR) status = nf90_get_var(ncid, varid, value)
            if (status /= NF90_NOERR) call stop_with(name)
        end subroutine

        subroutine get_real(name, value)
            character(*), intent(in) :: name
            real(dp), intent(out) :: value
            status = nf90_inq_varid(ncid, name, varid)
            if (status == NF90_NOERR) status = nf90_get_var(ncid, varid, value)
            if (status /= NF90_NOERR) call stop_with(name)
        end subroutine

        subroutine stop_with(what)
            character(*), intent(in) :: what
            integer :: ignored
            message = trim(path)//': cannot read '//what
            if (status == NF90_NOERR) status = 1
            ignored = nf90_close(ncid)
        end subroutine
    end subroutine

    subroutine mgrid_write(path, grid, status)
        !! The same variables as MAKEGRID writes, classic NetCDF.
        character(*), intent(in) :: path
        type(mgrid_t), intent(in) :: grid
        integer, intent(out) :: status
        character(2), parameter :: names(3) = ['br', 'bp', 'bz']
        character(3) :: tag
        integer :: ncid, drad, dzee, dphi, dcur, dmode, k, c, ids(9), varid
        integer, allocatable :: field(:, :)
        allocate (field(3, grid%nextcur))
        status = nf90_create(path, NF90_CLOBBER, ncid)
        if (status /= NF90_NOERR) return
        call ok(nf90_def_dim(ncid, 'rad', grid%ir, drad))
        call ok(nf90_def_dim(ncid, 'zee', grid%jz, dzee))
        call ok(nf90_def_dim(ncid, 'phi', grid%kp, dphi))
        call ok(nf90_def_dim(ncid, 'external_coils', grid%nextcur, dcur))
        call ok(nf90_def_dim(ncid, 'dim_00001', 1, dmode))
        call ok(nf90_def_var(ncid, 'ir', NF90_INT, varid=ids(1)))
        call ok(nf90_def_var(ncid, 'jz', NF90_INT, varid=ids(2)))
        call ok(nf90_def_var(ncid, 'kp', NF90_INT, varid=ids(3)))
        call ok(nf90_def_var(ncid, 'nfp', NF90_INT, varid=ids(4)))
        call ok(nf90_def_var(ncid, 'nextcur', NF90_INT, varid=ids(5)))
        call ok(nf90_def_var(ncid, 'rmin', NF90_DOUBLE, varid=ids(6)))
        call ok(nf90_def_var(ncid, 'zmin', NF90_DOUBLE, varid=ids(7)))
        call ok(nf90_def_var(ncid, 'rmax', NF90_DOUBLE, varid=ids(8)))
        call ok(nf90_def_var(ncid, 'zmax', NF90_DOUBLE, varid=ids(9)))
        call ok(nf90_def_var(ncid, 'mgrid_mode', NF90_CHAR, [dmode], varid))
        call ok(nf90_def_var(ncid, 'raw_coil_cur', NF90_DOUBLE, [dcur], k))
        do k = 1, grid%nextcur
            write (tag, '(i3.3)') k
            do c = 1, 3
                call ok(nf90_def_var(ncid, names(c)//'_'//tag, NF90_DOUBLE, &
                    [drad, dzee, dphi], field(c, k)))
            end do
        end do
        call ok(nf90_enddef(ncid))
        if (status /= NF90_NOERR) return
        call ok(nf90_put_var(ncid, ids(1), grid%ir))
        call ok(nf90_put_var(ncid, ids(2), grid%jz))
        call ok(nf90_put_var(ncid, ids(3), grid%kp))
        call ok(nf90_put_var(ncid, ids(4), grid%nfp))
        call ok(nf90_put_var(ncid, ids(5), grid%nextcur))
        call ok(nf90_put_var(ncid, ids(6), grid%rmin))
        call ok(nf90_put_var(ncid, ids(7), grid%zmin))
        call ok(nf90_put_var(ncid, ids(8), grid%rmax))
        call ok(nf90_put_var(ncid, ids(9), grid%zmax))
        call ok(nf90_put_var(ncid, varid, [grid%mode]))
        if (nf90_inq_varid(ncid, 'raw_coil_cur', k) == NF90_NOERR) &
            call ok(nf90_put_var(ncid, k, grid%raw_coil_cur))
        do k = 1, grid%nextcur
            do c = 1, 3
                call ok(nf90_put_var(ncid, field(c, k), grid%b(c, :, :, :, k)))
            end do
        end do
        call ok(nf90_close(ncid))
    contains
        subroutine ok(code)
            integer, intent(in) :: code
            if (status == NF90_NOERR) status = code
        end subroutine
    end subroutine
end module milca_mgrid_file
