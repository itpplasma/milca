program test_winslow_seed
    ! Winslow seed of the nested map: (1) on a circular section the translation
    ! of the reference disk is conformal, hence the exact minimiser, and a
    ! perturbed interior returns to it; (2) a bean section whose harmonic
    ! extension folds is admitted with the Winslow seed, which has positive
    ! determinants in every plane; (3) a doubled boundary is not a simple curve,
    ! has no injective map and is rejected before any solve.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_ideal_map, only: ideal_map_t, IDEAL_MAP_BAD_INPUT
    use milca_grad_shafranov, only: gs_state_t
    use milca_winslow_seed, only: winslow_seed, winslow_plane_energy, WINSLOW_OK, &
        WINSLOW_NOT_SIMPLE
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(toroidal_profile_t) :: profile
    type(toroidal_seed_t) :: seed
    type(ideal_map_t) :: map
    type(gs_state_t) :: mesh
    real(dp) :: coefficients(4, -1:1, 0:2), x(2, 4096), energy, min_det, error
    integer :: status, a, j
    profile%iota = [0.4_dp]; profile%pressure = [0.0_dp]
    call ideal_seed_disk_mesh(mesh, 16, 3, status)
    if (status /= 0) error stop 'disk mesh'

    ! (1) circle R = 3 + cos theta, Z = sin theta: X = y + (1, 0)
    coefficients = 0
    coefficients(1, 0, 0) = 3; coefficients(1, 0, 1) = 1; coefficients(4, 0, 1) = 1
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 1, 1.0_dp, profile, status)
    if (status /= 0) error stop 'circle seed'
    do a = 1, map%nodes
        if (map%boundary(a)) cycle
        map%x(1:2, 1, a) = map%x(1:2, 1, a) + 0.05_dp*[sin(7.0_dp*a), cos(5.0_dp*a)]
    end do
    call winslow_seed(map, status)
    if (status /= WINSLOW_OK) error stop 'circle Winslow'
    error = 0
    do a = 1, map%nodes
        error = max(error, maxval(abs(map%x(1:2, 1, a) - map%reference(:, a) - [1.0_dp, 0.0_dp])))
    end do
    print '(a,es10.3)', 'circle: largest distance to the conformal map ', error
    if (error > 1e-8_dp) error stop 'circle: not the conformal map'

    ! (2) bean section, shifted with phi (nfp = 2): the harmonic extension folds
    coefficients = 0
    coefficients(1, 0, 0) = 3; coefficients(1, 0, 1) = 1; coefficients(4, 0, 1) = 1
    ! X + i Z - 3 = z + 0.4 z^2 + 0.15 conj(z)^2 on the boundary: its harmonic
    ! extension has the Jacobian |1 + 0.8 z|^2 - 0.09 |z|^2 < 0 near z = -1
    coefficients(1, 0, 2) = 0.55_dp; coefficients(4, 0, 2) = 0.25_dp
    coefficients(1, 1, 1) = 0.1_dp; coefficients(4, 1, 1) = 0.1_dp
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 2, 1.0_dp, profile, status)
    if (status == 0) error stop 'bean: the harmonic extension should fold'
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 2, 1.0_dp, profile, status, &
        winslow=.true.)
    if (status /= 0) error stop 'bean: Winslow seed not admitted'
    do j = 1, map%nphi
        do a = 1, map%nodes
            x(:, a) = matmul(map%x(1:2, :, a), map%trig(:, j))
        end do
        call winslow_plane_energy(map, x(:, :map%nodes), 0.0_dp, energy, min_det)
        print '(a,i0,a,es10.3,a,es10.3)', 'bean plane ', j, ': energy ', energy, &
            ', smallest det ', min_det
        if (.not. min_det > 0) error stop 'bean: folded plane'
    end do

    ! (3) R = 3 + 0.5 cos 2 theta, Z = sin theta traverses one curve twice
    coefficients = 0
    coefficients(1, 0, 0) = 3; coefficients(1, 0, 2) = 0.5_dp; coefficients(4, 0, 1) = 1
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 1, 1.0_dp, profile, status, &
        winslow=.true.)
    if (status /= IDEAL_MAP_BAD_INPUT) error stop 'doubled boundary not rejected as input'
    call winslow_seed(map, status)
    if (status /= WINSLOW_NOT_SIMPLE) error stop 'doubled boundary not found self-intersecting'
    print '(a)', 'PASS Winslow seed'
end program
