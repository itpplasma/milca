!> Cylinder kinetic layer response on a periodised window
!> (milca_kinetic_layer_cylinder), Gaussian units.
!> Usage: kl_cylinder INPUT. INPUT, free format:
!>   period[cm] M conservation_model(0..3) boltzmann_gyroaverage(0|1) Re(B_r) Im(B_r)[G]
!>   N n_species
!>   per species: Z mass[g]
!>   N rows: r[cm] k_par k_s[1/cm] omega_E[1/s] B0[G] E_r[statV/cm], then per species
!>           n[cm^-3] T[erg] dn/dr dT/dr nu[1/s]
!> Prints per window point r, Re and Im of Phi [statV], of the total parallel
!> current density [statA/cm^2] and of the parallel current of each species.
program kl_cylinder
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_kinetic_layer_cylinder, only: kl_background_t, kl_kernels, kl_solve, kl_series, &
        KL_OK
    implicit none
    type(kl_background_t) :: bg
    complex(dp), allocatable :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :), phi(:), &
        jm(:, :), kp(:, :), kbb(:, :), stored(:, :, :)
    complex(dp) :: br
    real(dp) :: re, im
    character(1024) :: path
    integer :: u, n, ns, mmax, model, gyro, j, a, st

    if (command_argument_count() /= 1) error stop 'usage: kl_cylinder INPUT'
    call get_command_argument(1, path)
    open (newunit=u, file=trim(path), status='old', action='read')
    read (u, *) bg%period, mmax, model, gyro, re, im
    br = cmplx(re, im, dp)
    read (u, *) n, ns
    allocate (bg%species(ns), bg%r(n), bg%kpar(n), bg%ks(n), bg%omega_e(n), bg%b0(n), &
        bg%er(n))
    do a = 1, ns
        read (u, *) bg%species(a)%z, bg%species(a)%mass
        allocate (bg%species(a)%n(n), bg%species(a)%t(n), bg%species(a)%dn(n), &
            bg%species(a)%dt(n), bg%species(a)%nu(n))
    end do
    do j = 1, n
        read (u, *) bg%r(j), bg%kpar(j), bg%ks(j), bg%omega_e(j), bg%b0(j), bg%er(j), &
            (bg%species(a)%n(j), bg%species(a)%t(j), bg%species(a)%dn(j), &
            bg%species(a)%dt(j), bg%species(a)%nu(j), a=1, ns)
    end do
    close (u)

    ! Species by species, so that each species' current is available; the
    ! kernels of the sum are the sums of the kernels.
    allocate (jm(2*mmax + 1, ns))
    do a = 1, ns
        call kl_kernels(bg, mmax, model, kp, kbb, kjphi, kjb, st, &
            boltzmann_gyroaverage=gyro /= 0, only_species=a)
        if (st /= KL_OK) then
            write (error_unit, '(a,i0)') 'kl_kernels failed with status ', st
            error stop 1
        end if
        if (a == 1) then
            kphi = kp
            kb = kbb
        else
            kphi = kphi + kp
            kb = kb + kbb
        end if
        ! Current kernels kept per species; Phi follows below.
        jm(:, a) = br*kjb(:, mmax + 1)
        if (.not. allocated(stored)) allocate (stored(2*mmax + 1, 2*mmax + 1, ns))
        stored(:, :, a) = kjphi
    end do
    call kl_solve(kphi, kb, bg%period, mmax, br, phi, st)
    if (st /= KL_OK) error stop 'singular Poisson system'
    do a = 1, ns
        jm(:, a) = jm(:, a) + matmul(stored(:, :, a), phi)
    end do
    block
        complex(dp) :: p(n), jt(n), js(n, ns)
        p = kl_series(phi, bg%period, mmax, bg%r)
        do a = 1, ns
            js(:, a) = kl_series(jm(:, a), bg%period, mmax, bg%r)
        end do
        jt = sum(js, 2)
        do j = 1, n
            print '(f14.8,*(es20.11))', bg%r(j), p(j), jt(j), (js(j, a), a=1, ns)
        end do
    end block

end program kl_cylinder
