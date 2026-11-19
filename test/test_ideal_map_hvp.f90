program test_ideal_map_hvp
    ! The matrix-free Hessian-vector product ideal_map_hvp equals the
    ! assembled reduced Hessian times the same direction to roundoff, on the
    ! Landreman nested seed (P3, toroidal modes 2), on a perturbed map and with
    ! field-angle coefficients (angle_mpol 2).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_hessian, ideal_map_hvp
    use milca_ideal_seed_3d, only: ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_equilibrium_3d_driver, only: equilibrium_3d_driver_t, nested_integer_driver
    use milca_nested_oracles, only: nested_oracle_t
    use milca_grad_shafranov, only: gs_state_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    type(nested_oracle_t) :: oracle
    type(equilibrium_3d_driver_t) :: driver
    type(gs_state_t) :: mesh
    integer :: status, pass
    oracle%eps = 0.5_dp; oracle%delta = 1.0_dp/64; oracle%tau = 0
    call nested_integer_driver(oracle, 2, 2, driver, status)
    if (status /= 0) error stop 'driver'
    call ideal_seed_disk_mesh(mesh, 8, 3, status)
    if (status /= 0) error stop 'mesh'
    call ideal_seed_3d_init(driver%map, driver%seed, mesh, driver%boundary, driver%nfp, &
        driver%toroidal_flux, driver%profile, status)
    if (status /= 0) error stop 'seed'
    do pass = 1, 2
        if (pass == 2) driver%map%x = driver%map%x*(1 + 1.0e-3_dp*sin(1.0_dp &
            *reshape([(real(status + 1, dp), status=1, size(driver%map%x))], &
            shape(driver%map%x))))
        call compare(pass)
    end do
    call ideal_seed_3d_init(driver%map, driver%seed, mesh, driver%boundary, driver%nfp, &
        driver%toroidal_flux, driver%profile, status, angle_mpol=2)
    if (status /= 0) error stop 'seed with field angle'
    if (driver%map%nangle == 0) error stop 'no field-angle coefficients'
    driver%map%angle = [(0.01_dp*sin(2.3_dp*status), status=1, driver%map%nangle)]
    call compare(3)
    print '(a)', 'PASS matrix-free Hessian-vector product'
contains
    subroutine compare(pass)
        integer, intent(in) :: pass
        type(csr_matrix_t) :: h
        real(dp), allocatable :: d(:), y1(:), y2(:)
        integer :: i, k, n, st
        real(dp) :: error
        n = driver%map%nred
        allocate (d(n), y1(n), y2(n))
        d = [(sin(1.7_dp*i) + 0.3_dp*cos(5.1_dp*i), i=1, n)]
        call ideal_map_hessian(driver%map, driver%map%x, h)
        y1 = 0
        do i = 1, n
            do k = h%rowptr(i), h%rowptr(i + 1) - 1
                y1(i) = y1(i) + h%val(k)*d(h%colidx(k))
            end do
        end do
        call ideal_map_hvp(driver%map, driver%map%x, d, y2, st)
        if (st /= 0) error stop 'hvp status'
        error = maxval(abs(y1 - y2))/maxval(abs(y1))
        print '(a,i2,a,i6,a,es10.3)', 'pass', pass, ' unknowns', n, &
            ' max |H d - hvp| / max |H d| ', error
        if (.not. error < 1e-12_dp) error stop 'hvp differs from the assembled Hessian'
    end subroutine
end program
