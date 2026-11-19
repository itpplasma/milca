program test_map_3d_loop
    ! Loop integral of B on the circular torus R = R0 + rho cos theta,
    ! Z = rho sin theta (a = 1), which the P3 map represents exactly. With the
    ! frozen density (0, iota bt, bt) in (s, theta, phi), s = rho^2 and
    ! 2 pi bt = phiedge, the contravariant B^theta = 2 iota bt/R and
    ! oint B . dl = 2 iota bt s oint dtheta/R = 2 iota phiedge s/sqrt(R0^2 - s).
    ! The arcs must cover the circle once.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_ideal_map, only: ideal_map_t
    use milca_grad_shafranov, only: gs_state_t
    use milca_map_3d_loop, only: map_3d_loop_t, map_3d_loop_setup, map_3d_loop_integral
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3, iota = 0.4_dp, flux = 1.5_dp
    real(dp), parameter :: surfaces(4) = [1e-3_dp, 0.3_dp, 0.77_dp, 0.995_dp]
    type(toroidal_profile_t) :: profile
    type(toroidal_seed_t) :: seed
    type(ideal_map_t) :: map
    type(gs_state_t) :: mesh
    type(map_3d_loop_t) :: loop
    real(dp) :: coefficients(4, 0:0, 0:1), value, exact, error
    integer :: status, k
    call ideal_seed_disk_mesh(mesh, 16, 3, status)
    if (status /= 0) error stop 'disk mesh'
    coefficients = 0
    coefficients(1, 0, 0) = r0; coefficients(1, 0, 1) = 1; coefficients(4, 0, 1) = 1
    profile%iota = [iota]; profile%pressure = [0.0_dp]
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 1, flux, profile, status)
    if (status /= 0) error stop 'circular seed'
    do k = 1, size(surfaces)
        call map_3d_loop_setup(map, surfaces(k), loop, status)
        if (status /= 0) error stop 'loop setup'
        error = abs(sum(loop%weight) - 2*pi)
        if (error > 1e-13_dp) error stop 'arcs do not cover the circle once'
        value = map_3d_loop_integral(map, seed, loop, status)
        if (status /= 0) error stop 'loop integral'
        exact = 2*iota*flux*surfaces(k)/sqrt(r0**2 - surfaces(k))
        print '(a,f6.3,a,i4,a,es10.3)', 's ', surfaces(k), ': arc points ', size(loop%cell), &
            ', relative error ', abs(value/exact - 1)
        if (abs(value/exact - 1) > 1e-12_dp) error stop 'loop integral differs from the exact value'
    end do
    print '(a)', 'PASS map_3d loop'
end program
