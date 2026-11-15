!> Cylinder kinetic layer response from radial profiles
!> (milca_kinetic_layer_window, milca_kinetic_layer_cylinder), Gaussian units.
!> Usage: kl_profiles INPUT. INPUT, free format:
!>   directory with n.dat, Te.dat, Ti.dat, q.dat, Er.dat (two columns: r [cm] and
!>     n_e [cm^-3], T [eV], q, E_r [statV/cm])
!>   btor[G] R0[cm] m n Z_i A_i
!>   r_res[cm] core_half_width[cm] transition_width[cm] window_points
!>   M conservation_model(0..3) boltzmann_gyroaverage(0|1) Re(B_r) Im(B_r)[G]
!> Prints per window point r, Re and Im of Phi [statV] and of the parallel
!> current density [statA/cm^2].
program kl_profiles
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_window, KLW_OK
    use milca_kinetic_layer_cylinder, only: kl_background_t, kl_kernels, kl_solve, &
        kl_series, kl_current, KL_OK
    implicit none
    type(kl_profiles_t) :: p
    type(kl_background_t) :: bg
    complex(dp), allocatable :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :), phi(:), jm(:), &
        pr(:), jr(:)
    character(1024) :: path, dir
    real(dp) :: rres, a, w, re, im
    integer :: u, nwin, mmax, model, gyro, j, st

    if (command_argument_count() /= 1) error stop 'usage: kl_profiles INPUT'
    call get_command_argument(1, path)
    open (newunit=u, file=trim(path), status='old', action='read')
    read (u, '(a)') dir
    read (u, *) p%btor, p%r0, p%m, p%n, p%zi, p%ai
    read (u, *) rres, a, w, nwin
    read (u, *) mmax, model, gyro, re, im
    close (u)
    call load(trim(dir)//'/n.dat', p%r, p%ne)
    call load(trim(dir)//'/Te.dat', p%r, p%te)
    call load(trim(dir)//'/Ti.dat', p%r, p%ti)
    call load(trim(dir)//'/q.dat', p%r, p%q)
    call load(trim(dir)//'/Er.dat', p%r, p%er)
    call kl_window(p, rres, a, w, nwin, bg, st)
    if (st /= KLW_OK) then
        write (error_unit, '(a,i0)') 'kl_window failed with status ', st
        error stop 1
    end if
    call kl_kernels(bg, mmax, model, kphi, kb, kjphi, kjb, st, &
        boltzmann_gyroaverage=gyro /= 0)
    if (st /= KL_OK) error stop 'kl_kernels failed'
    call kl_solve(kphi, kb, bg%period, mmax, cmplx(re, im, dp), phi, st)
    if (st /= KL_OK) error stop 'singular Poisson system'
    jm = kl_current(kjphi, kjb, phi, cmplx(re, im, dp), mmax)
    pr = kl_series(phi, bg%period, mmax, bg%r)
    jr = kl_series(jm, bg%period, mmax, bg%r)
    do j = 1, size(bg%r)
        print '(f14.8,4es20.11)', bg%r(j), pr(j), jr(j)
    end do

contains

    !> Two columns r, y; profiles share the radii of the first file read.
    subroutine load(f, r, y)
        character(*), intent(in) :: f
        real(dp), allocatable, intent(inout) :: r(:)
        real(dp), allocatable, intent(out) :: y(:)
        real(dp), allocatable :: rr(:)
        real(dp) :: x1, x2
        integer :: v, n, ios, k
        open (newunit=v, file=f, status='old', action='read')
        n = 0
        do
            read (v, *, iostat=ios) x1, x2
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind (v)
        allocate (rr(n), y(n))
        do k = 1, n
            read (v, *) rr(k), y(k)
        end do
        close (v)
        if (.not. allocated(r)) then
            r = rr
        else if (size(r) /= n) then
            error stop 'profile files differ in their radii'
        else if (maxval(abs(r - rr)) > 1e-9_dp*maxval(abs(r))) then
            error stop 'profile files differ in their radii'
        end if
    end subroutine load

end program kl_profiles
