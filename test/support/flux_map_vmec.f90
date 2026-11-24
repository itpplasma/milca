module flux_map_vmec
    !! A VMEC state (milca_vmec_state) as a flux-map state of an fm3 case:
    !! per plane the V0 projection of u = rho (cos theta*, orientation *
    !! sin theta*), theta* = theta + lambda, rho = sqrt(s), the point of plane
    !! k at a block parameter inverted to VMEC's (rho, theta). Wall functions
    !! take the angle of their projected coefficients (as fm3_prolong).
    !! Used to evaluate the flux map's energy and gradient at VMEC's solution.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_project
    use milca_vmec_state, only: vmec_state_t, vmec_state_eval, vmec_state_invert, vmec_state_jet
    use milca_flux_map, only: flux_map_profiles
    use milca_flux_map_force, only: flux_map_force_t, flux_map_force_reset, flux_map_force_add, &
        flux_map_jet_to_cylindrical, flux_map_jet_product, flux_map_jet_cos_sin
    use milca_eulerian_clebsch, only: clebsch_nfp, clebsch_exact_field
    use flux_map_sections, only: fm3_case_t
    implicit none
    private

    public :: fm3_from_vmec, fm3_vmec_force

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    ! Allow a small wall interpolation error while excluding distant roots
    ! of the extrapolated Fourier series.
    real(dp), parameter :: projection_rho_max = 1.05_dp
    type(fm3_case_t), pointer :: cur => null()
    type(vmec_state_t), pointer :: vs => null()
    integer :: plane = 1, component = 1, failures = 0
    real(dp) :: orient = 1, guess(2) = [0.5_dp, 0.0_dp]

    public :: fm3_from_vmec_failures

contains

    integer function fm3_from_vmec_failures()
        fm3_from_vmec_failures = failures
    end function fm3_from_vmec_failures

    subroutine fm3_from_vmec(c, state, orientation, status)
        !! Overwrites c%fm%z. orientation +1 or -1: the sign of the second
        !! component of u relative to VMEC's theta*.
        type(fm3_case_t), intent(inout), target :: c
        type(vmec_state_t), intent(in), target :: state
        real(dp), intent(in) :: orientation
        integer, intent(out) :: status
        real(dp), allocatable :: cx(:), cy(:)
        real(dp) :: theta, prev(c%n0)
        integer :: k, i, n0

        cur => c
        vs => state
        orient = orientation
        failures = 0
        n0 = c%n0
        allocate (cx(n0), cy(n0))
        prev = 0
        do k = 1, c%nt
            plane = k
            component = 1
            call spline_derham_project(c%derham, 0, vmec_value, cx, status, c%degree + 3)
            if (status /= 0) exit
            component = 2
            call spline_derham_project(c%derham, 0, vmec_value, cy, status, c%degree + 3)
            if (status /= 0) exit
            do i = 1, n0
                if (c%fm%is_wall(i)) then
                    theta = atan2(cy(i), cx(i))
                    if (k > 1) theta = theta + two_pi*nint((prev(i) - theta)/two_pi)
                    prev(i) = theta
                    c%fm%z((k - 1)*n0 + i) = theta
                else
                    c%fm%z((k - 1)*n0 + i) = cx(i)
                    c%fm%z(n0*c%nt + (k - 1)*n0 + i) = cy(i)
                end if
            end do
        end do
        cur => null()
        vs => null()
        if (failures > 0) status = 1
    end subroutine fm3_from_vmec

    subroutine fm3_vmec_force(state, profiles, orientation, bands, radial, acc, exact_field, field_error)
        !! The force balance (milca_flux_map_force) of the VMEC state itself,
        !! without a mesh: u = rho (cos theta*, orientation sin theta*),
        !! theta* = theta + lambda, as a function of (R, phi, Z) through the
        !! 2-jets of R, Z, lambda (vmec_state_jet, cubic splines in s), on
        !! the midpoints of `radial` equal intervals of s, 4 mpol_max + 8
        !! angles theta and 4 ntor_max + 8 angles zeta of one field period.
        !! The profiles are those of the flux map (the deck), not read from
        !! the state.
        type(vmec_state_t), intent(in) :: state
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: orientation
        integer, intent(in) :: bands, radial
        type(flux_map_force_t), intent(out) :: acc
        procedure(clebsch_exact_field), optional :: exact_field
        real(dp), intent(out), optional :: field_error(2)
        real(dp) :: rho, theta, zeta, rj(10), zj(10), lj(10), tj(10), cj(10), sj(10), pj(10), &
                    xq(10), yq(10), xc(10), yc(10), w, field(3), reference(3)
        integer :: i, j, k, ntheta, nzeta

        ntheta = 4*maxval(state%m) + 8
        nzeta = 4*nint(maxval(abs(state%xn))/state%nfp) + 8
        call flux_map_force_reset(acc, bands)
        if (present(field_error)) field_error = 0
        do i = 1, radial
            rho = sqrt((i - 0.5_dp)/radial)
            pj = 0
            pj(1:2) = [rho, 1.0_dp]
            do k = 1, nzeta
                zeta = two_pi*(k - 1)/(nzeta*state%nfp)
                do j = 1, ntheta
                    theta = two_pi*(j - 1)/ntheta
                    call vmec_state_jet(state, rho, theta, zeta, rj, zj, lj)
                    tj = lj
                    tj(1) = tj(1) + theta
                    tj(3) = tj(3) + 1
                    call flux_map_jet_cos_sin(tj, cj, sj)
                    xq = flux_map_jet_product(pj, cj)
                    yq = orientation*flux_map_jet_product(pj, sj)
                    call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, xq, xc)
                    call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, yq, yc)
                    ! dV = R |d(R, Z)/d(rho, theta)| drho dtheta dzeta, ds = 2 rho drho
                    w = rj(1)*abs(rj(2)*zj(3) - rj(3)*zj(2))/(2*rho*radial)*two_pi/ntheta &
                        *two_pi/nzeta
                    call flux_map_force_add(acc, w, rj(1), xc, yc, profiles, field)
                    if (present(exact_field) .and. present(field_error)) then
                        call exact_field([rj(1), zeta, zj(1)], reference)
                        field_error = field_error + w*[sum((field - reference)**2), sum(reference**2)]
                    end if
                end do
            end do
        end do
    end subroutine fm3_vmec_force

    subroutine vmec_value(block, xi, t, values)
        !! Component `component` of u at (block, xi) of plane `plane`.
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:)
        real(dp) :: rz(2), zeta, rho, theta, r, z, lam, jac(2, 2), u(2), best, d
        integer :: st, i, j

        call spline_point_matrix(cur%space, [block], reshape(xi, [2, 1]), 0, rowptr, colidx, &
                                 val, st)
        if (st /= 0) error stop 'fm3_from_vmec: point'
        rz = matmul(cur%base(:, colidx, plane), val)
        zeta = two_pi*(plane - 1)/(cur%nt*clebsch_nfp)
        rho = guess(1)
        theta = guess(2)
        call vmec_state_invert(vs, rz, zeta, rho, theta, st, projection_rho_max)
        if (st /= 0) then
            ! coarse search, then Newton again
            best = huge(1.0_dp)
            do i = 1, 21
                do j = 0, 63
                    call vmec_state_eval(vs, 0.05_dp*i, two_pi*j/64, zeta, r, z, lam, jac)
                    d = norm2([r, z] - rz)
                    if (d < best) then
                        best = d
                        rho = 0.05_dp*i
                        theta = two_pi*j/64
                    end if
                end do
            end do
            call vmec_state_invert(vs, rz, zeta, rho, theta, st, projection_rho_max)
            if (st /= 0) failures = failures + 1
        end if
        guess = [rho, theta]
        call vmec_state_eval(vs, rho, theta, zeta, r, z, lam, jac)
        u = rho*[cos(theta + lam), orient*sin(theta + lam)]
        values(1) = u(component)
        associate (unused => t)
        end associate
    end subroutine vmec_value

end module flux_map_vmec
