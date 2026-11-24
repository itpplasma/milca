program test_equilibrium_3d
    ! test_equilibrium_3d TAU MPOL NTOR NODES DEGREE GATE_B GATE_P GATE_LABEL
    !                     [COARSE_MPOL COARSE_NTOR COARSE_NODES COARSE_DEGREE]
    ! Gates are about twice the observed sampled errors (docs/ideal-seed-3d.md).
    ! Without a coarse level, a reversed-iota control must be clearly worse.
    ! With a coarse level, every sampled error must decrease under refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_harness, only: equilibrium_resolution_t, &
        equilibrium_run_report_t, run_equilibrium_case, write_equilibrium_report
    use milca_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use milca_grad_shafranov, only: gs_state_t
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_prolong, ideal_map_energy, ideal_map_reduce, IDEAL_MAP_OK
    use milca_map_3d_equilibrium, only: map_3d_sample
    use milca_map_3d_readback_leaves, only: map_3d_cartesian
    use milca_equilibrium_3d_driver, only: equilibrium_3d_driver_t, &
        nested_integer_driver, run_boundary_profile_equilibrium
    use milca_nested_oracles, only: nested_oracle_t, NESTED_SHEARED
    use milca_nested_seed_input, only: nested_seed_boundary, nested_seed_profile, NESTED_SEED_NFP, NESTED_SEED_OK
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(nested_oracle_t) :: oracle
    type(equilibrium_resolution_t) :: resolution, coarse
    character(512) :: arg
    character(32) :: deck_family
    integer :: mpol, ntor, cmpol, cntor
    real(dp) :: gates(3), fine(3), base(3), reversed(3)
    oracle%eps=0.5_dp;oracle%delta=1.0_dp/64
    call get_command_argument(1,arg)
    if (arg == 'sheared-forward') then
        call set_sheared_oracle()
        call check_sheared_forward()
        stop
    end if
    if (arg == 'deck') then
        if (command_argument_count() < 2 .or. command_argument_count() > 3) &
            error stop 'usage: deck OUTPUT [landreman-sheared]'
        call get_command_argument(2, arg)
        call get_command_argument(3, deck_family)
        select case (trim(deck_family))
        case ('')
        case ('landreman-sheared')
            call set_sheared_oracle()
        case default
            error stop 'unknown manufactured deck'
        end select
        call write_deck(trim(arg))
        stop
    end if
    if(command_argument_count()/=8 .and. command_argument_count()/=12) &
        error stop 'usage: TAU MPOL NTOR NODES DEGREE GATE_B GATE_P GATE_LABEL [coarse]'
    read(arg,*) oracle%tau
    call get_command_argument(2,arg);read(arg,*) mpol
    call get_command_argument(3,arg);read(arg,*) ntor
    call get_command_argument(4,arg);read(arg,*) resolution%boundary_nodes
    call get_command_argument(5,arg);read(arg,*) resolution%degree
    call get_command_argument(6,arg);read(arg,*) gates(1)
    call get_command_argument(7,arg);read(arg,*) gates(2)
    call get_command_argument(8,arg);read(arg,*) gates(3)
    call solve(mpol,ntor,resolution,1.0_dp,.true.,fine)
    if(fine(1)>gates(1)) error stop 'sampled B error above recorded gate'
    if(fine(2)>gates(2)) error stop 'sampled pressure error above recorded gate'
    if(fine(3)>gates(3)) error stop 'sampled surface-label error above recorded gate'
    if(command_argument_count()==12) then
        call get_command_argument(9,arg);read(arg,*) cmpol
        call get_command_argument(10,arg);read(arg,*) cntor
        call get_command_argument(11,arg);read(arg,*) coarse%boundary_nodes
        call get_command_argument(12,arg);read(arg,*) coarse%degree
        call solve(cmpol,cntor,coarse,1.0_dp,.false.,base)
        print '(a,3es12.4)', 'refinement_ratio_B_p_label=',base/fine
        if(any(fine>=base)) error stop 'sampled errors did not decrease under refinement'
    else
        ! The prescribed iota sign is physical: the reversed transform gives
        ! the same pressure and labels but must not reproduce the oracle B.
        call solve(mpol,ntor,resolution,-1.0_dp,.false.,reversed)
        print '(a,es12.4)', 'reversed_iota_B_sample_rms_relative=',reversed(1)
        if(reversed(1)<5*fine(1)) error stop 'B is insensitive to the iota sign'
    end if
contains
    subroutine set_sheared_oracle()
        oracle%family = NESTED_SHEARED
        oracle%eps = 0.6_dp
        oracle%saxis = 2.4_dp
        oracle%delta = 0.2_dp
    end subroutine

    subroutine check_sheared_forward()
        type(ideal_map_t) :: first, final
        type(toroidal_seed_t) :: first_seed, final_seed
        type(toroidal_profile_t) :: profile
        type(ideal_map_solve_report_t) :: report
        type(gs_state_t) :: mesh
        real(dp), allocatable :: wall(:, :, :), gradient(:), reduced(:)
        real(dp) :: flux, fit, energy, reference, initial_error(2), final_error(2)
        integer :: status

        call nested_seed_boundary(oracle, 8, 8, wall, status)
        if (status /= NESTED_SEED_OK) error stop 'sheared boundary'
        call nested_seed_profile(oracle, 12, 1e-7_dp, profile, flux, fit, status)
        if (status /= NESTED_SEED_OK) error stop 'sheared profiles'
        call ideal_seed_disk_mesh(mesh, 12, 3, status)
        if (status /= 0) error stop 'sheared initial mesh'
        call ideal_seed_3d_init(first, first_seed, mesh, wall, NESTED_SEED_NFP, flux, profile, &
            status, angle_mpol=6, winslow=.true.)
        if (status /= IDEAL_MAP_OK) error stop 'sheared initial seed'
        call ideal_map_solve(first, 1e-5_dp, 48, report, status, stationary=.true.)
        if (status /= IDEAL_MAP_OK) error stop 'sheared initial stationary solve'
        call sheared_field_errors(first, first_seed, initial_error)

        call ideal_seed_disk_mesh(mesh, 16, 4, status)
        if (status /= 0) error stop 'sheared refined mesh'
        call ideal_seed_3d_init(final, final_seed, mesh, wall, NESTED_SEED_NFP, flux, profile, &
            status, angle_mpol=6, winslow=.true.)
        if (status /= IDEAL_MAP_OK) error stop 'sheared refined seed'
        allocate (gradient(final%ndof), reduced(final%nred))
        call ideal_map_energy(final, final%x, energy, status, gradient)
        if (status /= IDEAL_MAP_OK) error stop 'sheared refined seed gradient'
        call ideal_map_reduce(final, gradient, reduced)
        where (final%fixed) reduced = 0
        reference = norm2(reduced)
        call ideal_map_prolong(final, first, status)
        if (status /= IDEAL_MAP_OK) error stop 'sheared stationary transfer'
        call ideal_map_solve(final, 1e-5_dp, 48, report, status, &
            reference_gradient=reference, stationary=.true.)
        print '(a,i0,a,i0,a,es14.6)', 'sheared refined status ', status, &
            ' iterations ', report%iterations, ' relative gradient ', report%gradient
        if (status /= IDEAL_MAP_OK) error stop 'sheared refined stationary solve'
        call sheared_field_errors(final, final_seed, final_error)
        print '(a,2es14.6)', 'sheared initial sampled relative B, p errors ', initial_error
        print '(a,2es14.6)', 'sheared refined sampled relative B, p errors ', final_error
        if (any(final_error >= initial_error)) error stop 'sheared sampled errors did not decrease'
        if (any(final_error > 5e-4_dp)) error stop 'sheared sampled error above 0.05 percent'
    end subroutine

    subroutine sheared_field_errors(map, seed, errors)
        type(ideal_map_t), intent(in) :: map
        type(toroidal_seed_t), intent(in) :: seed
        real(dp), intent(out) :: errors(2)
        real(dp), parameter :: xi(2, 2) = reshape([0.30_dp, 0.40_dp, 0.22_dp, 0.13_dp], [2, 2])
        integer, parameter :: planes = 17
        real(dp) :: point(3), value(EQ_VALUES), x(3), b(3), exact(20), sums(4), phi
        integer :: plane, cell, k, status

        sums = 0
        do plane = 1, planes
            phi = 2*acos(-1.0_dp)*(plane - 0.37_dp)/(map%nfp*planes)
            do cell = 1, map%cells
                do k = 1, size(xi, 2)
                    call map_3d_sample(map, seed, cell, xi(:, k), phi, point, value, status)
                    if (status /= 0) error stop 'sheared candidate sample'
                    call map_3d_cartesian(point(1), point(2), point(3), &
                        value(EQ_BR), value(EQ_BZ), value(EQ_BPHI), x(1), x(2), x(3), b(1), b(2), b(3))
                    call oracle%point(x, exact, status)
                    if (status /= 0) error stop 'sheared analytic sample'
                    sums = sums + [sum((b - exact(3:5))**2), sum(exact(3:5)**2), &
                        (value(EQ_P) - exact(2))**2, exact(2)**2]
                end do
            end do
        end do
        errors = sqrt(sums(1:4:2)/sums(2:4:2))
    end subroutine

    subroutine write_deck(path)
        ! Give every code the same boundary and normalized-flux profiles.
        character(*), intent(in) :: path
        integer, parameter :: mmax = 8, nmax = 8
        real(dp), parameter :: mu0 = 4e-7_dp*acos(-1.0_dp)
        real(dp), allocatable :: boundary(:, :, :)
        type(toroidal_profile_t) :: profile
        real(dp) :: flux, err, deck_sign
        integer :: m, n, u, status, profile_degree
        logical :: sheared

        call nested_seed_boundary(oracle, mmax, nmax, boundary, status)
        if (status /= NESTED_SEED_OK) error stop 'manufactured benchmark boundary'
        if (maxval(abs(boundary(2:3, :, :))) > 1e-12_dp) error stop 'benchmark boundary is not symmetric'
        sheared = oracle%family == NESTED_SHEARED
        profile_degree = merge(12, 1, sheared)
        call nested_seed_profile(oracle, profile_degree, merge(1e-7_dp, 1e-8_dp, sheared), &
                                 profile, flux, err, status)
        if (status /= NESTED_SEED_OK) error stop 'manufactured benchmark profiles'
        open(newunit=u, file=path, status='replace')
        write(u, '(a)') '&INDATA', ' LASYM = F, NCURR = 0, GAMMA = 0,'
        write(u, '(3(a,i0),a)') ' NFP = ', NESTED_SEED_NFP, ', MPOL = ', merge(14, mmax + 1, sheared), &
                               ', NTOR = ', merge(10, nmax, sheared), ','
        ! The integer deck reverses the oracle's flux and straight-angle senses.
        deck_sign = merge(1.0_dp, -1.0_dp, sheared)
        write(u, '(a,es24.16,a)') ' PHIEDGE = ', deck_sign*flux, ','
        write(u, '(a,*(es24.16,1x))') ' AI = ', deck_sign*profile%iota
        write(u, '(a,*(es24.16,1x))') ' AM = ', profile%pressure/mu0
        if (sheared) then
            write(u, '(a)') ' NS_ARRAY = 12, 25, 101,', ' FTOL_ARRAY = 1e-11, 1e-12, 1e-13,', &
                           ' NITER_ARRAY = 1200, 2500, 15000,'
        else
            write(u, '(a)') ' NS_ARRAY = 12, 25, 51,', ' FTOL_ARRAY = 1e-10, 1e-12, 1e-13,', &
                           ' NITER_ARRAY = 4000, 6000, 10000,'
        end if
        do m = 0, mmax
            do n = -nmax, nmax
                if (abs(boundary(1, n, m)) > 1e-14_dp) &
                    write(u, '(a,i0,a,i0,a,es24.16,a)') ' RBC(', n, ',', m, ') = ', boundary(1, n, m), ','
                if (abs(boundary(4, n, m)) > 1e-14_dp) &
                    write(u, '(a,i0,a,i0,a,es24.16,a)') ' ZBS(', n, ',', m, ') = ', boundary(4, n, m), ','
            end do
        end do
        write(u, '(a)') '/'
        close(u)
    end subroutine write_deck

    subroutine solve(mpol,ntor,resolution,iota_sign,gate,errors)
        integer, intent(in) :: mpol, ntor
        type(equilibrium_resolution_t), intent(in) :: resolution
        real(dp), intent(in) :: iota_sign
        logical, intent(in) :: gate
        real(dp), intent(out) :: errors(3)
        type(equilibrium_3d_driver_t) :: driver
        type(equilibrium_state_t) :: equilibrium
        type(equilibrium_run_report_t) :: report
        character(64) :: name
        integer :: status
        real(dp) :: started, finished
        name='landreman_3d_nested'
        if(abs(oracle%tau)>tiny(1.0_dp)) name='issan_3d_nonsymmetric'
        started=omp_get_wtime()
        call nested_integer_driver(oracle,mpol,ntor,driver,status)
        if(status/=0) error stop 'integer-family boundary/profile data'
        if(driver%profile%iota(1)<=0) error stop 'integer-family iota must be positive'
        driver%profile%iota=iota_sign*driver%profile%iota
        call run_equilibrium_case(trim(name),'map',resolution,'legacy', &
            run_boundary_profile_equilibrium,driver,report,equilibrium)
        finished=omp_get_wtime()
        errors=[driver%magnetic_error,driver%pressure_error,driver%label_error]
        print '(a,2i3,2i4,f6.1)', 'case mpol ntor nodes degree iota_sign=', &
            mpol,ntor,resolution%boundary_nodes,resolution%degree,iota_sign
        call write_equilibrium_report(6,report)
        print '(a,es24.16)', 'B_sample_rms_relative=',errors(1)
        print '(a,es24.16)', 'p_sample_rms_relative=',errors(2)
        print '(a,es24.16)', 'surface_label_max_relative=',errors(3)
        print '(a,es24.16)', 'workflow_wall_s=',finished-started
        print '(a,es24.16)', 'energy_gradient_wall_s=',driver%solve%energy_s
        print '(a,es24.16)', 'hessian_wall_s=',driver%solve%hessian_s
        print '(a,es24.16)', 'factor_solve_wall_s=',driver%solve%factor_s
        print '(a,i0)', 'cores=',omp_get_max_threads()
        if(.not. gate) return
        print '(a)', '[V] 1728 equal-weight candidate-coordinate samples; no volume norm'
        print '(a)', '[V] force and a posteriori estimator effectivity unavailable'
        if(report%status/=0) error stop 'native D8 3D equilibrium'
        if(report%residual>driver%tolerance) error stop 'Newton tolerance not attained'
        if(report%oracle_error>=report%seed_error) error stop 'seed did not improve'
    end subroutine
end program
