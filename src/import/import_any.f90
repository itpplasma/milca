module milca_import_any
    ! Format detection and the full conversion of one foreign input.
    use milca_import_text, only: read_lines, lower
    use milca_equilibrium_import, only: imported_equilibrium_t, import_fail, import_normalize
    use milca_import_gvec, only: import_gvec
    use milca_import_vmec, only: import_vmec, import_wout
    use milca_import_desc_spectre, only: import_desc, import_spectre
    implicit none
    private
    public :: import_equilibrium, detect_format

contains

    subroutine import_equilibrium(path, format, eq)
        ! format: auto, gvec, desc, spectre, vmec or wout
        character(*), intent(in) :: path, format
        type(imported_equilibrium_t), intent(out) :: eq
        character(len=:), allocatable :: kind
        kind = lower(format)
        if (kind == 'auto') kind = detect_format(path)
        select case (kind)
        case ('gvec')
            call import_gvec(path, eq)
        case ('desc')
            call import_desc(path, eq)
        case ('spectre')
            call import_spectre(path, eq)
        case ('vmec')
            call import_vmec(path, eq)
        case ('wout')
            call import_wout(path, eq)
        case default
            call import_fail(eq, 'unknown input format of '//trim(path))
            return
        end select
        call import_normalize(eq)
    end subroutine

    function detect_format(path) result(kind)
        character(*), intent(in) :: path
        character(len=:), allocatable :: kind, ext, base, line
        character(len=:), allocatable :: lines(:)
        logical :: ok
        integer :: i
        base = path(index(path, '/', back=.true.) + 1:)
        ext = ''
        if (index(base, '.') > 0) ext = lower(base(index(base, '.', back=.true.) + 1:))
        kind = 'unknown'
        if (ext == 'nc') then
            kind = 'wout'
            return
        end if
        if (ext == 'yaml' .or. ext == 'yml' .or. ext == 'ini') then
            kind = 'gvec'
            return
        end if
        call read_lines(path, lines, ok)
        if (.not. ok) return
        do i = 1, size(lines)
            line = lower(trim(adjustl(lines(i))))
            if (len(line) == 0) cycle
            if (line(1:1) == '!' .or. line(1:1) == '#') cycle
            if (ext == 'toml') then
                if (line == '[physics]') then
                    kind = 'spectre'
                    return
                end if
                cycle
            end if
            if (line(1:1) == '&' .or. line(1:1) == '$') then
                kind = 'vmec'
                return
            end if
            if (index(line, 'whichinitequilibrium') > 0) then
                kind = 'gvec'
                return
            end if
            if (line(1:2) == 'l:' .or. line(1:2) == 'm:' .or. line(1:2) == 'n:') then
                kind = 'desc'
                return
            end if
        end do
        if (ext == 'toml') kind = 'gvec'
    end function
end module milca_import_any
