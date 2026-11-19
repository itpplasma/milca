program test_ideal_map_mode_pc
    ! Overlapping Schwarz preconditioner in n of the map Hessian. (1) On an axisymmetric
    ! state carried by a map with toroidal modes 2 the couplings between
    ! different n vanish, so the preconditioner is H^-1: M^-1 H d = d on the
    ! free dofs. (2) On the 3D Landreman seed (P3, 8
    ! boundary nodes, toroidal modes 2, solved), conjugate gradients on
    ! H + mu |diag H| for mu = 1e-2, 1e-4 and 1e-6 need fewer iterations with
    ! it than with Jacobi, at most 50 to the inexact-Newton forcing 1e-3, and
    ! no more when deflated by the recycled Ritz vectors of a first solve.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_hessian, ideal_map_hvp, ideal_map_hessian_store_t, &
        ideal_map_store_hessians, ideal_map_solve, ideal_map_solve_report_t
    use milca_ideal_map_mode_pc, only: ideal_map_mode_pc_t, ideal_map_mode_pc_build, &
        ideal_map_mode_pc_apply
    use milca_ideal_map_krylov, only: deflated_pcg
    use milca_ideal_seed_3d, only: ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_equilibrium_3d_driver, only: equilibrium_3d_driver_t, nested_integer_driver
    use milca_nested_oracles, only: nested_oracle_t
    use milca_grad_shafranov, only: gs_state_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    type(nested_oracle_t) :: oracle
    type(equilibrium_3d_driver_t) :: driver
    type(gs_state_t) :: mesh
    type(ideal_map_hessian_store_t) :: store
    type(ideal_map_mode_pc_t) :: pc
    type(ideal_map_solve_report_t) :: report
    real(dp), allocatable :: dg(:)
    real(dp) :: mu
    logical :: blocks_mode
    integer :: status, k, its_band, its_jacobi, its_newton, its_deflated
    real(dp), parameter :: mus(3) = [1.0e-2_dp, 1.0e-4_dp, 1.0e-6_dp]
    integer :: nodes
    character(16) :: arg
    nodes = 8
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read (arg, *) nodes
    end if
    oracle%eps = 0.5_dp; oracle%delta = 1.0_dp/64; oracle%tau = 0
    call nested_integer_driver(oracle, 2, 2, driver, status)
    if (status /= 0) error stop 'driver'
    call ideal_seed_disk_mesh(mesh, nodes, 3, status)
    if (status /= 0) error stop 'mesh'
    call axisymmetric_check()
    call ideal_seed_3d_init(driver%map, driver%seed, mesh, driver%boundary, driver%nfp, &
        driver%toroidal_flux, driver%profile, status)
    if (status /= 0) error stop '3D seed'
    ! the solved equilibrium, where Newton steps are taken
    call ideal_map_solve(driver%map, 1.0e-6_dp, 160, report, status)
    if (status /= 0) error stop '3D solve'
    call diagonal()
    call ideal_map_store_hessians(driver%map, driver%map%x, store)
    do k = 1, size(mus)
        mu = mus(k)
        call ideal_map_mode_pc_build(driver%map, store, mu, pc, status)
        if (status /= 0) error stop 'build 3D'
        its_band = pcg(.true., 1e-8_dp)
        its_jacobi = pcg(.false., 1e-8_dp)
        its_newton = pcg(.true., 1e-3_dp)
        call recycled(its_deflated)
        print '(a,es9.1,a,i5,a,i5,a,i5,a,i5)', '3D equilibrium mu ', mu, &
            ': CG iterations to 1e-8 Schwarz ', its_band, ', Jacobi ', its_jacobi, &
            '; to the Newton forcing 1e-3 ', its_newton, &
            ', deflated by 40 recycled Ritz vectors ', its_deflated
        if (its_band >= its_jacobi) error stop 'Schwarz does not beat Jacobi'
        if (its_deflated > its_newton) error stop 'deflation adds iterations'
        if (nodes == 8 .and. its_deflated > 50) error stop 'more than 50 CG iterations to 1e-3'
    end do
    print '(a)', 'PASS band-in-n preconditioner'
contains
    subroutine axisymmetric_check()
        real(dp), allocatable :: boundary(:, :, :), d(:), y(:), z(:), s(:)
        real(dp) :: error
        integer :: i, n
        boundary = driver%boundary
        boundary(:, -2:-1, :) = 0; boundary(:, 1:2, :) = 0
        call ideal_seed_3d_init(driver%map, driver%seed, mesh, boundary, driver%nfp, &
            driver%toroidal_flux, driver%profile, status)
        if (status /= 0) error stop 'axisymmetric seed'
        call ideal_map_store_hessians(driver%map, driver%map%x, store)
        call ideal_map_mode_pc_build(driver%map, store, 0.0_dp, pc, status)
        if (status /= 0) error stop 'build'
        n = driver%map%nred
        allocate (d(n), y(n), z(n), s(n))
        d = [(sin(1.7_dp*i) + 0.3_dp*cos(5.1_dp*i), i=1, n)]
        where (driver%map%fixed) d = 0
        call ideal_map_hvp(driver%map, driver%map%x, d, y, status, store)
        if (status /= 0) error stop 'hvp'
        s = 0
        do i = lbound(pc%part, 1), ubound(pc%part, 1)
            s(pc%part(i)%rdof) = max(s(pc%part(i)%rdof), pc%part(i)%shift)
        end do
        if (maxval(s) > 0) then
            ! the 8-node seed of the registered test is positive definite
            if (nodes == 8) error stop 'axisymmetric Hessian not positive definite'
            print '(a)', 'axisymmetric: seed Hessian indefinite, comparison skipped'
            return
        end if
        call ideal_map_mode_pc_apply(driver%map, pc, y, z, status)
        if (status /= 0) error stop 'apply'
        error = maxval(abs(z - d), mask=.not. driver%map%fixed)/maxval(abs(d))
        print '(a,es10.3)', 'axisymmetric: max |M^-1 H d - d| / |d| ', error
        if (.not. error < 1e-8_dp) error stop 'the band differs from the Hessian'
    end subroutine

    subroutine diagonal()
        type(csr_matrix_t) :: h
        integer :: i, kk
        call ideal_map_hessian(driver%map, driver%map%x, h)
        allocate (dg(driver%map%nred))
        do i = 1, driver%map%nred
            do kk = h%rowptr(i), h%rowptr(i + 1) - 1
                if (h%colidx(kk) == i) dg(i) = h%val(kk)
            end do
        end do
    end subroutine

    subroutine recycled(its)
        ! a first solve gives 12 Ritz vectors; the second, other right-hand side
        ! is deflated with them, as in consecutive Newton steps
        integer, intent(out) :: its
        real(dp), allocatable :: b(:), x(:), ritz(:, :), w(:, :)
        logical :: negative
        integer :: n, i, j, k
        blocks_mode = .true.
        n = driver%map%nred
        allocate (b(n), x(n), w(n, 0))
        b = [(cos(0.9_dp*i), i=1, n)]
        where (driver%map%fixed) b = 0
        call deflated_pcg(operator, precondition, b, 1e-8_dp, 2000, w, x, its, negative, &
            40, ritz)
        ! orthonormal basis of the Ritz vectors
        w = ritz
        k = 0
        do j = 1, size(w, 2)
            do i = 1, k
                w(:, j) = w(:, j) - dot_product(w(:, i), w(:, j))*w(:, i)
            end do
            if (.not. norm2(w(:, j)) > 1e-10_dp) cycle
            k = k + 1
            w(:, k) = w(:, j)/norm2(w(:, j))
        end do
        b = [(sin(1.3_dp*i) + 0.2_dp, i=1, n)]
        where (driver%map%fixed) b = 0
        call deflated_pcg(operator, precondition, b, 1e-3_dp, 2000, w(:, :k), x, its, &
            negative, 40, ritz)
        if (negative) error stop 'negative curvature in the deflated solve'
    end subroutine

    integer function pcg(blocks, tol) result(its)
        logical, intent(in) :: blocks
        real(dp), intent(in) :: tol
        real(dp), allocatable :: b(:), x(:), r(:), z(:), p(:), q(:)
        real(dp) :: rz, rz_new, alpha, b0
        integer :: n, i
        blocks_mode = blocks
        n = driver%map%nred
        allocate (b(n), x(n), r(n), z(n), p(n), q(n))
        b = [(cos(0.9_dp*i), i=1, n)]
        where (driver%map%fixed) b = 0
        x = 0; r = b; b0 = norm2(b)
        call precondition(r, z)
        p = z; rz = dot_product(r, z)
        do its = 1, 5000
            call operator(p, q)
            alpha = rz/dot_product(p, q)
            x = x + alpha*p; r = r - alpha*q
            if (norm2(r) <= tol*b0) return
            call precondition(r, z)
            rz_new = dot_product(r, z)
            p = z + (rz_new/rz)*p; rz = rz_new
        end do
    end function

    subroutine operator(v, w)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        call ideal_map_hvp(driver%map, driver%map%x, v, w, status, store)
        w = w + mu*abs(dg)*v
        where (driver%map%fixed) w = v
    end subroutine

    subroutine precondition(v, w)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        if (.not. blocks_mode) then
            w = v/((1 + mu)*abs(dg))
            where (driver%map%fixed) w = v
            return
        end if
        call ideal_map_mode_pc_apply(driver%map, pc, v, w, status)
    end subroutine
end program
