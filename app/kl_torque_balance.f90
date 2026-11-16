!> Torque balance at one resonant surface with the cylinder kinetic layer
!> response (milca_kinetic_layer_torque), Gaussian units, torques in N m.
!> Usage: kl_torque_balance INPUT. INPUT, free format:
!>   directory with n.dat, Te.dat, Ti.dat, q.dat, Er.dat (as kl_profiles)
!>   btor[G] R0[cm] m n Z_i A_i
!>   r_res[cm] core_half_width[cm] transition_width[cm] window_points
!>   M conservation_model(0..3) boltzmann_gyroaverage(0|1)
!>   V_E min, max, step [km/s]
!>   s = sign(m) sign(B_z), vacuum resonant field B_vac [G] at the reference current
!>   viscous stiffness K [N m s/rad], lever h [m], unperturbed V_E [km/s], path points
!> Prints V_E, Re and Im sigma [statA/G], |B_m/B_vac| and the torque [N m] at the
!> reference current, then each fold as (I/I_ref)^2 and the V_E before and after it.
program kl_torque_balance
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_profiles_read, KLW_OK
    use milca_kinetic_layer_cylinder, only: kl_light, KL_OK
    use milca_kinetic_layer_torque, only: klt_sigma, klt_torque, klt_branch, klt_fold_t
    implicit none
    type(kl_profiles_t) :: p
    type(klt_fold_t), allocatable :: f(:)
    complex(dp), allocatable :: sg(:)
    real(dp), allocatable :: ve(:), tau(:), v(:), y(:)
    character(1024) :: path, dir
    real(dp) :: rres, a, w, v0, v1, dv, s, bvac, kv, h, vu
    integer :: u, nwin, mmax, model, gyro, nv, npts, i, st
    logical :: locked

    if (command_argument_count() /= 1) error stop 'usage: kl_torque_balance INPUT'
    call get_command_argument(1, path)
    open (newunit=u, file=trim(path), status='old', action='read')
    read (u, '(a)') dir
    read (u, *) p%btor, p%r0, p%m, p%n, p%zi, p%ai
    read (u, *) rres, a, w, nwin
    read (u, *) mmax, model, gyro
    read (u, *) v0, v1, dv
    read (u, *) s, bvac
    read (u, *) kv, h, vu, npts
    close (u)
    call kl_profiles_read(trim(dir), p, st)
    if (st /= KLW_OK) error stop 'profile files missing, short or with different radii'
    nv = nint((v1 - v0)/dv) + 1
    ve = [(v0 + dv*(i - 1), i = 1, nv)]
    allocate (sg(nv))
    call klt_sigma(p, rres, a, w, nwin, mmax, model, ve, sg, st, &
        boltzmann_gyroaverage=gyro /= 0)
    if (st /= KL_OK) then
        write (error_unit, '(a,i0)') 'layer response failed with status ', st
        error stop 1
    end if
    tau = klt_torque(sg, rres, p%r0, p%m, p%n, s, bvac)
    print '(a)', '# V_E[km/s] Re,Im sigma[statA/G] |B_m/B_vac| tau[N m]'
    do i = 1, nv
        print '(f12.5,2es18.9,f12.7,es16.7)', ve(i), sg(i), &
            1/abs(1 - (0.0_dp, 1.0_dp)*sg(i)/(kl_light*rres*s)), tau(i)
    end do
    call klt_branch(ve, tau, vu, kv, h, npts, v, y, f, locked)
    print '(a,l1,a,f10.5)', '# path locked at a torque zero: ', locked, ', ends at V_E ', v(size(v))
    do i = 1, size(f)
        print '(a,es14.6,a,f10.5,a,f10.5)', '# fold (I/I_ref)^2 ', f(i)%y, ' V_E before ', &
            f(i)%v_before, ' after ', f(i)%v_after
    end do
end program kl_torque_balance
