program test_nested_equilibrium_app
    ! Checks output/nested_equilibrium.nc against the same input solved by a
    ! VMEC code: usage RESULT RAXIS WB [VOLUME [IOTA_EDGE]]. The axis must agree
    ! to 1.5 mm, the magnetic energy to 2e-5, the volume (when given, > 0) to
    ! 1e-4 and the edge iota (when given) to 1e-3. With RAXIS = failure the
    ! result must instead report converged = false and a failure reason, which
    ! contains the word given as WB (when given).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
        nf90_get_att, NF90_NOWRITE, NF90_NOERR, NF90_GLOBAL
    implicit none
    character(256) :: path, arg
    character(8) :: converged
    character(256) :: failure
    real(dp) :: volume, raxis, wb, iota_edge, expected(4)
    logical :: check_iota
    integer :: ncid, k
    call get_command_argument(1, path)
    expected = -1
    check_iota = command_argument_count() >= 5
    call get_command_argument(2, arg)
    if (trim(arg) /= 'failure') then
        do k = 2, min(command_argument_count(), 5)
            call get_command_argument(k, arg); read (arg, *) expected(k - 1)
        end do
    end if
    if (nf90_open(trim(path), NF90_NOWRITE, ncid) /= NF90_NOERR) error stop 'open result'
    if (nf90_get_att(ncid, NF90_GLOBAL, 'converged', converged) /= NF90_NOERR) &
        error stop 'converged attribute'
    call get_command_argument(2, arg)
    if (trim(arg) == 'failure') then
        failure = ''
        if (nf90_get_att(ncid, NF90_GLOBAL, 'failure', failure) /= NF90_NOERR) &
            error stop 'failure attribute'
        if (trim(converged) /= 'false' .or. len_trim(failure) == 0) error stop 'failure not reported'
        ! an optional third argument is a word the reason must contain
        if (command_argument_count() >= 3) then
            call get_command_argument(3, arg)
            if (index(failure, trim(arg)) == 0) error stop 'unexpected failure reason'
        end if
        print '(2a)', 'PASS failure reported: ', trim(failure)
        stop
    end if
    volume = scalar('volume_p'); raxis = scalar('raxis_cc'); wb = scalar('wb')
    iota_edge = scalar('iotaf_edge')
    if (nf90_close(ncid) /= NF90_NOERR) error stop 'close result'
    print '(a,4es16.8)', 'volume, raxis, wb, iota edge: ', volume, raxis, wb, iota_edge
    if (trim(converged) /= 'true') error stop 'not converged'
    if (abs(raxis - expected(1)) > 1.5e-3_dp) error stop 'axis'
    if (abs(wb/expected(2) - 1) > 2e-5_dp) error stop 'magnetic energy'
    if (expected(3) > 0 .and. abs(volume/expected(3) - 1) > 1e-4_dp) error stop 'volume'
    if (check_iota .and. abs(iota_edge - expected(4)) > 1e-3_dp) error stop 'edge iota'
    print '(a)', 'PASS nested equilibrium'
contains
    real(dp) function scalar(name)
        character(*), intent(in) :: name
        integer :: varid
        if (nf90_inq_varid(ncid, name, varid) /= NF90_NOERR) error stop name
        if (nf90_get_var(ncid, varid, scalar) /= NF90_NOERR) error stop name
    end function
end program
