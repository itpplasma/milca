!> Kinetic layers at several resonant surfaces matched to an outer solution
!> (milca_kinetic_layer_match), Gaussian units, torques in N m.
!> Usage: kl_layer_match INPUT. INPUT, free format:
!>   directory with n.dat, Te.dat, Ti.dat, q.dat, Er.dat (as kl_profiles)
!>   btor[G] R0[cm] n Z_i A_i
!>   core_half_width[cm] transition_width[cm] window_points M
!>     conservation_model(0..3) boltzmann_gyroaverage(0|1)
!>   N, the number of surfaces
!>   N lines: m r_s[cm] Re C, Im C [G cm]
!>   N lines: row j of the outer tearing matrix E, Re E_j1, Im E_j1, ...
!> Prints per surface m, r_s, Delta [1/cm], sigma [statA/G], |psi/r_s| [G]
!> and the torque [N m], then the sum of the torques.
program kl_layer_match
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_profiles_read, KLW_OK
    use milca_kinetic_layer_cylinder, only: KL_OK
    use milca_kinetic_layer_delta, only: kl_delta
    use milca_kinetic_layer_match, only: klm_solve, KLM_OK
    implicit none
    type(kl_profiles_t) :: p
    character(1024) :: path, dir
    complex(dp), allocatable :: c(:), e(:, :), d(:), sg(:), psi(:)
    real(dp), allocatable :: r_s(:), torque(:), row(:), cr(:, :)
    integer, allocatable :: m(:)
    real(dp) :: a, w
    integer :: u, n, nwin, mmax, model, gyro, nr, j, st

    if (command_argument_count() /= 1) error stop 'usage: kl_layer_match INPUT'
    call get_command_argument(1, path)
    open (newunit=u, file=trim(path), status='old', action='read')
    read (u, '(a)') dir
    read (u, *) p%btor, p%r0, n, p%zi, p%ai
    read (u, *) a, w, nwin, mmax, model, gyro
    read (u, *) nr
    allocate (m(nr), r_s(nr), cr(2, nr), c(nr), e(nr, nr), d(nr), sg(nr), psi(nr), &
        torque(nr), row(2*nr))
    do j = 1, nr
        read (u, *) m(j), r_s(j), cr(:, j)
    end do
    c = cmplx(cr(1, :), cr(2, :), dp)
    do j = 1, nr
        read (u, *) row
        e(j, :) = cmplx(row(1::2), row(2::2), dp)
    end do
    close (u)
    call kl_profiles_read(trim(dir), p, st)
    if (st /= KLW_OK) error stop 'profile files missing, short or with different radii'
    p%n = n
    do j = 1, nr
        p%m = m(j)
        call kl_delta(p, r_s(j), a, w, nwin, mmax, model, d(j), sg(j), st, &
            boltzmann_gyroaverage=gyro /= 0)
        if (st /= KL_OK) then
            write (error_unit, '(a,i0,a,i0)') 'layer response failed at m = ', m(j), &
                ' with status ', st
            error stop 1
        end if
    end do
    call klm_solve(e, c, d, sg, r_s, m, n, p%r0, psi, torque, st)
    if (st /= KLM_OK) then
        write (error_unit, '(a,i0)') 'matching failed with status ', st
        error stop 1
    end if
    print '(a)', '# m r_s[cm] Re,Im Delta[1/cm] Re,Im sigma[statA/G] |psi/r_s|[G] T[N m]'
    do j = 1, nr
        print '(i4,f10.4,4es16.7,2es16.7)', m(j), r_s(j), d(j), sg(j), abs(psi(j))/r_s(j), &
            torque(j)
    end do
    print '(a,es16.7)', '# sum T[N m] ', sum(torque)
end program kl_layer_match
