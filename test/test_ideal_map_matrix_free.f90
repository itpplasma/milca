program test_ideal_map_matrix_free
    ! Matrix-free Newton against the direct ideal_map_solve on the Landreman
    ! nested seed (P3, 8 boundary nodes, toroidal modes 2): the Hessian
    ! diagonal from the stored point Hessians equals that of the assembled
    ! matrix to roundoff, the nonzero estimate of the memory-based selection
    ! equals its pattern, and both solves reach the gradient tolerance with
    ! the same energy to 1e-10.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_t, ideal_map_hessian, ideal_map_solve, &
        ideal_map_solve_report_t, ideal_map_hessian_store_t, ideal_map_store_hessians
    use milca_ideal_map_matrix_free, only: ideal_map_solve_matrix_free, &
        ideal_map_hessian_diagonal, ideal_map_hessian_nonzeros
    use milca_ideal_seed_3d, only: ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_equilibrium_3d_driver, only: equilibrium_3d_driver_t, nested_integer_driver
    use milca_nested_oracles, only: nested_oracle_t
    use milca_grad_shafranov, only: gs_state_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    type(nested_oracle_t) :: oracle
    type(equilibrium_3d_driver_t) :: driver
    type(gs_state_t) :: mesh
    type(ideal_map_t) :: direct, free
    type(ideal_map_solve_report_t) :: rd, rf
    type(ideal_map_hessian_store_t) :: store
    type(csr_matrix_t) :: h
    real(dp), allocatable :: dg(:), dm(:)
    real(dp) :: error
    integer :: status, i, k, nodes, modes
    character(16) :: arg, only
    ! optional arguments for measurements: NODES [direct|free [MPOL=NTOR]]
    nodes = 8; only = 'both'; modes = 2
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read (arg, *) nodes
    end if
    if (command_argument_count() >= 2) call get_command_argument(2, only)
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read (arg, *) modes
    end if
    oracle%eps = 0.5_dp; oracle%delta = 1.0_dp/64; oracle%tau = 0
    call nested_integer_driver(oracle, modes, modes, driver, status)
    if (status /= 0) error stop 'driver'
    call ideal_seed_disk_mesh(mesh, nodes, 3, status)
    if (status /= 0) error stop 'mesh'
    call ideal_seed_3d_init(driver%map, driver%seed, mesh, driver%boundary, driver%nfp, &
        driver%toroidal_flux, driver%profile, status)
    if (status /= 0) error stop 'seed'
    if (only /= 'both') then
        call measure()
        stop
    end if
    ! Hessian diagonal
    call ideal_map_hessian(driver%map, driver%map%x, h)
    allocate (dg(driver%map%nred), dm(driver%map%nred))
    do i = 1, driver%map%nred
        do k = h%rowptr(i), h%rowptr(i + 1) - 1
            if (h%colidx(k) == i) dg(i) = h%val(k)
        end do
    end do
    call ideal_map_store_hessians(driver%map, driver%map%x, store)
    call ideal_map_hessian_diagonal(driver%map, store, dm)
    error = maxval(abs(dg - dm))/maxval(abs(dg))
    print '(a,es10.3)', 'Hessian diagonal: max relative difference ', error
    if (.not. error < 1e-12_dp) error stop 'Hessian diagonal differs'
    print '(a,2i10)', 'Hessian nonzeros: assembled, estimated', h%rowptr(size(h%rowptr)) - 1, &
        ideal_map_hessian_nonzeros(driver%map)
    if (ideal_map_hessian_nonzeros(driver%map) /= h%rowptr(size(h%rowptr)) - 1) &
        error stop 'Hessian nonzero estimate differs from the assembled pattern'
    block
        type(ideal_map_t) :: angled
        type(csr_matrix_t) :: ha
        call ideal_seed_3d_init(angled, driver%seed, mesh, driver%boundary, driver%nfp, &
            driver%toroidal_flux, driver%profile, status, angle_mpol=2)
        if (status /= 0) error stop 'seed with field angle'
        call ideal_map_hessian(angled, angled%x, ha)
        print '(a,2i10)', 'with field angle: assembled, estimated', &
            ha%rowptr(size(ha%rowptr)) - 1, ideal_map_hessian_nonzeros(angled)
        if (ideal_map_hessian_nonzeros(angled) /= ha%rowptr(size(ha%rowptr)) - 1) &
            error stop 'Hessian nonzero estimate differs with field angle'
    end block
    direct = driver%map; free = driver%map
    call ideal_map_solve(direct, 1.0e-8_dp, 200, rd, status)
    if (status /= 0) error stop 'direct solve'
    call ideal_map_solve_matrix_free(free, 1.0e-8_dp, 200, rf, status)
    print '(a,i4,es22.14,a,i4,es22.14,a,i6)', 'direct: iterations', rd%iterations, rd%energy, &
        '; matrix-free:', rf%iterations, rf%energy, ' CG iterations', rf%krylov_iterations
    if (status /= 0) error stop 'matrix-free solve'
    error = abs(rf%energy/rd%energy - 1)
    print '(a,es10.3)', 'relative energy difference ', error
    if (.not. error < 1e-10_dp) error stop 'matrix-free energy differs'
    print '(a)', 'PASS matrix-free Newton'
contains
    subroutine measure()
        integer(8) :: c0, c1, rate
        call system_clock(c0, rate)
        if (only == 'direct') then
            call ideal_map_solve(driver%map, 1.0e-8_dp, 200, rd, status)
        else
            call ideal_map_solve_matrix_free(driver%map, 1.0e-8_dp, 200, rd, status)
        end if
        call system_clock(c1)
        print '(a,a,a,i8,a,i10,a,i4,a,es22.14,a,es10.3,a,i6,a,f8.2)', trim(only), &
            ': unknowns', '', driver%map%nred, ' Hessian nonzeros', &
            ideal_map_hessian_nonzeros(driver%map), ' iterations', rd%iterations, ' energy', &
            rd%energy, ' gradient', rd%gradient, ' CG', rd%krylov_iterations, ' seconds', &
            real(c1 - c0)/rate
    end subroutine
end program
