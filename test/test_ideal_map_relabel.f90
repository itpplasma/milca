program test_ideal_map_relabel
    ! ideal_map_relabel_space on the Landreman nested seed (P3, toroidal
    ! modes 2, 8 and 16 radial nodes) after a smooth deformation without
    ! stellarator symmetry. The bases are orthonormal and zero on the wall.
    ! Relative to a radial expansion of the surfaces, the energy slope along
    ! the reshuffles theta -> theta + f(rho) falls with h (a continuum null),
    ! while along the theta-dependent family it stays O(0.1-1) (a physical
    ! force). ideal_map_relabel_project removes the span.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_relabel_space, ideal_map_relabel_project, &
        ideal_map_energy, ideal_map_reduce, ideal_map_rdof
    use milca_ideal_seed_3d, only: ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_equilibrium_3d_driver, only: equilibrium_3d_driver_t, nested_integer_driver
    use milca_nested_oracles, only: nested_oracle_t
    use milca_grad_shafranov, only: gs_state_t
    implicit none
    real(dp) :: null8, null16, theta8, theta16
    call measure(8, null8, theta8)
    call measure(16, null16, theta16)
    if (.not. null16 < null8/4) error stop 'reshuffle slope does not fall with h'
    if (.not. null16 < 2e-2_dp) error stop 'reshuffle not flat'
    if (.not. min(theta8, theta16) > 1e-1_dp) error stop 'theta family flat'
    print '(a)', 'PASS relabel space'
contains
    subroutine measure(nodes, null_ratio, theta_ratio)
        ! largest slope ratio of the reshuffles and of the m = 1, k = 0 family
        integer, intent(in) :: nodes
        real(dp), intent(out) :: null_ratio, theta_ratio
        type(nested_oracle_t) :: oracle
        type(equilibrium_3d_driver_t) :: driver
        type(gs_state_t) :: mesh
        real(dp), allocatable :: basis(:, :), d(:), grad(:), gred(:), r(:)
        real(dp) :: energy, slope_radial
        integer :: status, i, n
        oracle%eps = 0.5_dp; oracle%delta = 1.0_dp/64; oracle%tau = 0
        call nested_integer_driver(oracle, 2, 2, driver, status)
        if (status /= 0) error stop 'driver'
        call ideal_seed_disk_mesh(mesh, nodes, 3, status)
        if (status /= 0) error stop 'mesh'
        call ideal_seed_3d_init(driver%map, driver%seed, mesh, driver%boundary, driver%nfp, &
            driver%toroidal_flux, driver%profile, status)
        if (status /= 0) error stop 'seed'
        ! a smooth deformation without stellarator symmetry, which alone would
        ! hide half of the slopes
        block
            real(dp) :: c(2), rad, rho, th, ext
            integer :: a, m
            associate (map => driver%map)
                c = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
                rad = maxval(norm2(map%reference - spread(c, 2, map%nodes), dim=1))
                ext = maxval(map%x(2, 1, :)) - minval(map%x(2, 1, :))
                do a = 1, map%nodes
                    rho = norm2(map%reference(:, a) - c)/rad
                    th = atan2(map%reference(2, a) - c(2), map%reference(1, a) - c(1))
                    do m = 1, min(3, map%nmode)
                        map%x(1, m, a) = map%x(1, m, a) + 0.02_dp*ext*(1 - rho**2)*rho &
                            *sin(th + 0.4_dp*m)
                        map%x(2, m, a) = map%x(2, m, a) + 0.02_dp*ext*(1 - rho**2)*rho &
                            *cos(2*th + 0.9_dp*m)
                    end do
                end do
            end associate
        end block
        associate (map => driver%map)
            n = map%nred
            allocate (d(n), grad(map%ndof), gred(n), r(n))
            call ideal_map_energy(map, map%x, energy, status, grad)
            if (status /= 0) error stop 'energy'
            call ideal_map_reduce(map, grad, gred)
            where (map%fixed) gred = 0
            ! comparison: radial expansion (1 - rho^2) (X - X_axis) of the surfaces
            r = 0
            block
                real(dp) :: centre(2), radius, rho
                integer :: a, m, axis
                centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
                radius = maxval(norm2(map%reference - spread(centre, 2, map%nodes), dim=1))
                axis = minloc(norm2(map%reference - spread(centre, 2, map%nodes), dim=1), 1)
                do a = 1, map%nodes
                    rho = norm2(map%reference(:, a) - centre)/radius
                    do m = 1, map%nmode
                        r(ideal_map_rdof(map, 1, m, a)) = (1 - rho**2) &
                            *(map%x(1, m, a) - map%x(1, m, axis))
                        r(ideal_map_rdof(map, 2, m, a)) = (1 - rho**2) &
                            *(map%x(2, m, a) - map%x(2, m, axis))
                    end do
                end do
            end block
            where (map%fixed) r = 0
            r = r/norm2(r)
            slope_radial = abs(dot_product(r, gred))
            call ideal_map_relabel_space(map, map%x, basis, status)
            if (status /= 0) error stop 'relabel space status'
            call check_basis(basis, map%fixed, 4)
            null_ratio = maxval(abs(matmul(gred, basis)))/slope_radial
            call ideal_map_relabel_space(map, map%x, basis, status, mpol=1, ntor=1)
            if (status /= 0) error stop 'coarse space status'
            call check_basis(basis, map%fixed, 36)
            ! directions 4..9: f = 1 - rho^2, m = 1, orthogonal to the reshuffles
            theta_ratio = maxval(abs(matmul(gred, basis(:, 4:9))))/slope_radial
            print '(a,i3,a,2es10.2)', 'nodes', nodes, &
                '  slope / radial slope: reshuffles, m = 1', null_ratio, theta_ratio
            d = r
            call ideal_map_relabel_project(basis, d)
            if (.not. maxval(abs(matmul(d, basis))) < 1e-14_dp) error stop 'projection'
        end associate
    end subroutine

    subroutine check_basis(basis, fixed, expected)
        real(dp), intent(in) :: basis(:, :)
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: expected
        real(dp) :: orth
        integer :: i, j
        orth = 0
        do j = 1, size(basis, 2)
            orth = max(orth, maxval(abs(matmul(basis(:, j), basis) &
                - merge(1.0_dp, 0.0_dp, [(i == j, i=1, size(basis, 2))]))))
        end do
        print '(a,i4,a,i6,a,es10.3)', 'relabel directions', size(basis, 2), ' of', &
            size(basis, 1), '  max |B^T B - I| ', orth
        if (size(basis, 2) /= expected) error stop 'number of relabel directions'
        if (.not. orth < 1e-12_dp) error stop 'basis not orthonormal'
        if (any(abs(pack(basis, spread(fixed, 2, size(basis, 2)))) > 0)) &
            error stop 'basis moves the wall'
    end subroutine
end program
