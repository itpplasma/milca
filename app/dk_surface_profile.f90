!> Local 3D transport on the surfaces of a Boozer `.bc` file, for example a
!! helical-core tokamak (milca_dk_surface_device, input described there).
!! Usage: dk_surface_profile INPUT [LOCAL_PREFIX]
!! Prints per surface: s, status, q, B0 [T], dpsi/ds [T m^2], then per species
!! pair (a, b) the transport matrix d(a, b, i, j), i, j = 1..3; a line
!! 'flows' with <B u_par> per species [T m/s] and <j_par B> [T A/m^2]; with a
!! potential, a line 'phi1' with its <B u_par> per species and <j_par B>.
!! With LOCAL_PREFIX, writes LOCAL_PREFIX_<surface>.dat with theta, phi_B,
!! the local bootstrap current and (with a potential) the potential-driven
!! current [A/m^2] on the grid.
program dk_surface_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_dk_surface_device, only: dk_surface_input_t, dk_surface_result_t, &
        dk_surface_read_input, dk_surface_profile_run => dk_surface_profile, DKS_OK
    implicit none
    real(dp), parameter :: twopi = 8*atan(1.0_dp)
    type(boozer_file_t) :: eq
    type(dk_surface_input_t) :: inp
    type(dk_surface_result_t), allocatable :: r(:)
    character(1024) :: path, prefix, name
    integer :: status, i, a, b, j, k, u, nper

    if (command_argument_count() < 1) error stop 'usage: dk_surface_profile INPUT [LOCAL_PREFIX]'
    call get_command_argument(1, path)
    prefix = ''
    if (command_argument_count() > 1) call get_command_argument(2, prefix)
    call dk_surface_read_input(trim(path), inp, status)
    if (status /= DKS_OK) error stop 'cannot read the input'
    call load_boozer_file(trim(inp%eq_path), eq, status)
    if (status /= BOOZER_OK) error stop 'cannot read the Boozer file'
    nper = eq%d%nper
    call dk_surface_profile_run(eq, inp, r, status)
    do i = 1, size(inp%s)
        print '(a,f12.8,i3,3es18.9)', 'surface ', inp%s(i), r(i)%status, r(i)%q, r(i)%b0, &
            r(i)%psi_prime
        if (r(i)%status /= DKS_OK) cycle
        do a = 1, size(inp%z)
            do b = 1, size(inp%z)
                print '(2i3,9es18.9)', a, b, r(i)%d(a, b, :, :)
            end do
        end do
        print '(a,*(es18.9))', 'flows ', r(i)%flow, r(i)%bootstrap
        if (allocated(r(i)%flow_phi)) print '(a,*(es18.9))', 'phi1 ', r(i)%flow_phi, &
            r(i)%current_phi
        if (len_trim(prefix) == 0) cycle
        write (name, '(a,a,i0,a)') trim(prefix), '_', i, '.dat'
        open (newunit=u, file=trim(name), status='replace', action='write')
        write (u, '(a)') '# theta phi_B j_bs[A/m^2] j_phi[A/m^2]'
        do k = 1, inp%n_zeta
            do j = 1, inp%n_theta
                associate (p => j + inp%n_theta*(k - 1))
                    if (allocated(r(i)%j_phi)) then
                        write (u, '(4es18.9)') twopi*(j - 1)/inp%n_theta, &
                            twopi*(k - 1)/(inp%n_zeta*nper), r(i)%j_bs(p), r(i)%j_phi(p)
                    else
                        write (u, '(4es18.9)') twopi*(j - 1)/inp%n_theta, &
                            twopi*(k - 1)/(inp%n_zeta*nper), r(i)%j_bs(p), 0.0_dp
                    end if
                end associate
            end do
        end do
        close (u)
    end do
    if (status /= DKS_OK) then
        write (error_unit, '(a,i0)') 'dk_surface_profile finished with status ', status
        error stop 1
    end if
end program dk_surface_profile
