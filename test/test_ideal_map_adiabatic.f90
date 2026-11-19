program test_ideal_map_adiabatic
    ! Adiabatic pressure p = K(s)/V'(s)^gamma in the map lane:
    ! 1. energy gradient and Hessian (frozen part at P_h plus the low-rank
    !    term G Q G^T) against energy and gradient differences,
    ! 2. a fixed-pressure equilibrium with pressure p* is stationary for the
    !    mass K = p* V'^gamma of that state, with P_h = p* to roundoff,
    ! 3. Newton with the Woodbury step returns to that state from the seed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_map_t, ideal_map_energy, ideal_map_hessian, &
        ideal_map_reduce, ideal_map_expand, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_set_adiabatic, ideal_map_adiabatic_lowrank, &
        ideal_map_adiabatic_profile, IDEAL_MAP_OK
    use milca_ideal_seed_3d, only: toroidal_seed_t, toroidal_profile_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_nested_oracles, only: nested_oracle_t, NESTED_SHEARED
    use milca_nested_seed_input, only: nested_seed_boundary, nested_seed_profile, &
        NESTED_SEED_OK
    use milca_grad_shafranov, only: gs_state_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    real(dp), parameter :: gamma = 5.0_dp/3, p0 = 2.0e-3_dp
    integer, parameter :: degree = 6
    type(nested_oracle_t) :: oracle
    type(toroidal_profile_t) :: profile
    type(ideal_map_t) :: map, solved, probe
    type(toroidal_seed_t) :: seed
    type(ideal_map_solve_report_t) :: report
    real(dp), allocatable :: boundary(:, :, :), x_seed(:, :, :), g(:), gred(:)
    real(dp) :: flux, fit_error, e_fixed, e_ad, e_back, g_fixed, g_ad, worst
    real(dp) :: s(7), v(7), p(7)
    integer :: status, i
    oracle%family = NESTED_SHEARED; oracle%eps = 0; oracle%delta = 0.14_dp
    call nested_seed_profile(oracle, 12, 1e-9_dp, profile, flux, fit_error, status)
    if (status /= NESTED_SEED_OK) error stop 'profile'
    call nested_seed_boundary(oracle, 8, 1, boundary, status)
    if (status /= NESTED_SEED_OK) error stop 'boundary'
    ! linear pressure: P_h reproduces it exactly for degree >= 1
    profile%pressure = [p0, -p0]
    if (allocated(profile%pressure_values)) deallocate (profile%pressure_values)
    call setup(map)
    x_seed = map%x

    ! fixed-pressure equilibrium p*
    call ideal_map_solve(map, 1e-10_dp, 60, report, status)
    if (status /= IDEAL_MAP_OK) error stop 'fixed-pressure solve'
    allocate (g(map%ndof), gred(map%nred))
    call ideal_map_energy(map, map%x, e_fixed, status, g)
    call ideal_map_reduce(map, g, gred)
    g_fixed = norm2(pack(gred, .not. map%fixed))
    solved = map

    ! mass of that state: K = p* V'^gamma with V' of the solved map
    probe = solved
    call ideal_map_set_adiabatic(probe, seed, gamma, degree, unit_mass, status)
    if (status /= IDEAL_MAP_OK) error stop 'probe mass'
    call ideal_map_set_adiabatic(map, seed, gamma, degree, solved_mass, status)
    if (status /= IDEAL_MAP_OK) error stop 'adiabatic mass'
    s = [(0.05_dp + 0.15_dp*i, i=0, 6)]
    call ideal_map_adiabatic_profile(map, map%x, s, v, p, status)
    if (status /= IDEAL_MAP_OK) error stop 'profile query'
    worst = maxval(abs(p - p0*(1 - s)))/p0
    print '(a,es11.3)', 'P_h against p* at the fixed-pressure state:', worst
    if (worst > 1e-12_dp) error stop 'projected pressure'
    call ideal_map_energy(map, map%x, e_ad, status, g)
    call ideal_map_reduce(map, g, gred)
    g_ad = norm2(pack(gred, .not. map%fixed))
    print '(a,2es11.3)', 'reduced gradient fixed / adiabatic:', g_fixed, g_ad
    if (abs(g_ad - g_fixed) > 1e-10_dp*max(g_fixed, 1e-30_dp) + 1e-14_dp) &
        error stop 'adiabatic gradient at the fixed-pressure state'

    ! derivatives at the seed and after two adiabatic Newton steps
    map%x = x_seed
    call directional(map, 1)
    call ideal_map_solve(map, 1e-3_dp, 2, report, status)
    call directional(map, 2)

    ! Newton (Woodbury step) from the seed back to the stationary state
    map%x = x_seed
    call ideal_map_solve(map, 1e-10_dp, 60, report, status)
    print '(a,i0,a,i0)', 'adiabatic Newton status ', status, ', iterations ', &
        report%iterations
    if (status /= IDEAL_MAP_OK) error stop 'adiabatic solve'
    call ideal_map_energy(map, map%x, e_back, status)
    print '(a,es11.3)', 'energy difference to the fixed-pressure state:', &
        (e_back - e_ad)/abs(e_ad)
    if (abs(e_back - e_ad) > 1e-10_dp*abs(e_ad)) error stop 'adiabatic state'
    call ideal_map_adiabatic_profile(map, map%x, s, v, p, status)
    worst = maxval(abs(p - p0*(1 - s)))/p0
    print '(a,es11.3)', 'P_h against p* after the adiabatic solve:', worst
    if (worst > 1e-6_dp) error stop 'adiabatic pressure profile'
    print '(a)', 'PASS ideal map adiabatic'
contains
    subroutine setup(m)
        type(ideal_map_t), intent(out) :: m
        type(gs_state_t) :: mesh
        integer :: st
        call ideal_seed_disk_mesh(mesh, 8, 3, st)
        if (st /= 0) error stop 'mesh'
        call ideal_seed_3d_init(m, seed, mesh, boundary(:, -1:1, :), 2, flux, &
            profile, st)
        if (st /= IDEAL_MAP_OK) error stop 'seed'
    end subroutine

    function unit_mass(t) result(k)
        real(dp), intent(in) :: t
        real(dp) :: k
        k = 1 + 0*t
    end function

    function solved_mass(t) result(k)
        real(dp), intent(in) :: t
        real(dp) :: k, vt(1), pt(1)
        integer :: st
        call ideal_map_adiabatic_profile(probe, probe%x, [t], vt, pt, st)
        if (st /= IDEAL_MAP_OK) error stop 'probe profile'
        k = p0*(1 - t)*vt(1)**gamma
    end function

    subroutine directional(m, pass)
        type(ideal_map_t), intent(inout) :: m
        integer, intent(in) :: pass
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: d(:), full(:), gp(:), gm(:), hd(:), rp(:), rm(:)
        real(dp), allocatable :: lr(:, :), q(:, :), u(:)
        real(dp) :: e0, ep, em, h, slope, curvature
        integer :: st, k, j
        allocate (d(m%nred), full(m%ndof), gp(m%ndof), gm(m%ndof), hd(m%nred), &
            rp(m%nred), rm(m%nred), u(m%nred))
        do k = 1, m%nred
            d(k) = sin(1.3_dp*k + pass)
        end do
        where (m%fixed) d = 0
        call ideal_map_expand(m, d, full)
        call ideal_map_energy(m, m%x, e0, st, g)
        if (st /= IDEAL_MAP_OK) error stop 'energy'
        call ideal_map_reduce(m, g, gred)
        call ideal_map_hessian(m, m%x, matrix)
        hd = 0
        do k = 1, m%nred
            do j = matrix%rowptr(k), matrix%rowptr(k + 1) - 1
                hd(k) = hd(k) + matrix%val(j)*d(matrix%colidx(j))
            end do
        end do
        call ideal_map_adiabatic_lowrank(m, m%x, lr, q, st)
        if (st /= IDEAL_MAP_OK) error stop 'low rank'
        ! G Q G^T d in reduced coordinates
        do k = lbound(lr, 2), ubound(lr, 2)
            call ideal_map_reduce(m, lr(:, k), u)
            lr(:size(u), k) = u
        end do
        hd = hd + matmul(lr(:m%nred, :), matmul(q, matmul(d, lr(:m%nred, :))))
        h = 1e-6_dp
        call shifted(m, full, h, ep, gp)
        call shifted(m, full, -h, em, gm)
        call ideal_map_reduce(m, gp, rp)
        call ideal_map_reduce(m, gm, rm)
        slope = abs((ep - em)/(2*h) - dot_product(gred, d)) &
            /(norm2(pack(gred, .not. m%fixed))*norm2(d))
        curvature = maxval(abs((rp - rm)/(2*h) - hd), mask=.not. m%fixed) &
            /maxval(abs(hd), mask=.not. m%fixed)
        print '(a,i0,2es11.3)', 'adiabatic derivatives (pass, slope, Hessian action): ', &
            pass, slope, curvature
        if (slope > 1e-6_dp) error stop 'adiabatic gradient vs energy difference'
        if (curvature > 1e-6_dp) error stop 'adiabatic Hessian vs gradient difference'
    end subroutine

    subroutine shifted(m, full, h, e, grad)
        type(ideal_map_t), intent(in) :: m
        real(dp), intent(in) :: full(:), h
        real(dp), intent(out) :: e, grad(:)
        integer :: st
        call ideal_map_energy(m, m%x + h*reshape(full, shape(m%x)), e, st, grad)
        if (st /= IDEAL_MAP_OK) error stop 'shifted energy'
    end subroutine
end program test_ideal_map_adiabatic
