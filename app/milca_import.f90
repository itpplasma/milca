program milca_import
    ! Converts a GVEC, DESC, SPECTRE or VMEC input (or a VMEC wout file) into
    ! nested_equilibrium's namelist (docs/equilibrium-import.md).
    ! Exit status 0 when written, 3 when readable but not modelled by
    ! nested_equilibrium (reason on stderr), 1 for an unreadable input.
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use milca_equilibrium_import, only: imported_equilibrium_t, import_write, IMPORT_OK, &
        IMPORT_UNSUPPORTED
    use milca_import_any, only: import_equilibrium
    implicit none
    type(imported_equilibrium_t), allocatable :: eq
    character(4096) :: arguments(2), word
    character(len=:), allocatable :: format, input, output
    integer :: count, i, k
    logical :: ok
    count = command_argument_count()
    format = 'auto'
    k = 0
    i = 1
    do while (i <= count)
        call get_command_argument(i, word)
        if (trim(word) == '--format' .and. i < count) then
            call get_command_argument(i + 1, word)
            format = trim(word)
            i = i + 2
            cycle
        end if
        k = k + 1
        if (k <= 2) arguments(k) = word
        i = i + 1
    end do
    if (k /= 2) then
        write (error_unit, '(a)') 'usage: milca_import [--format auto|gvec|desc|spectre|vmec|wout] '// &
            'INPUT NAMELIST'
        stop 1, quiet=.true.
    end if
    input = trim(arguments(1)); output = trim(arguments(2))
    allocate (eq)
    call import_equilibrium(input, format, eq)
    do i = 1, eq%n_notes
        write (output_unit, '(a)') trim(eq%notes(i))
    end do
    if (eq%status == IMPORT_UNSUPPORTED) then
        write (error_unit, '(a)') 'unsupported: '//trim(eq%reason)
        stop 3, quiet=.true.
    else if (eq%status /= IMPORT_OK) then
        write (error_unit, '(a)') 'error: '//trim(eq%reason)
        stop 1, quiet=.true.
    end if
    call import_write(eq, output, ok)
    if (.not. ok) then
        write (error_unit, '(a)') 'error: cannot write '//output
        stop 1, quiet=.true.
    end if
end program milca_import
