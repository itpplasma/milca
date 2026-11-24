module solovev_inverse_check_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_map, only: ideal_frozen_t
    use solovev_label_fields, only: psi, grad_psi, psi_b, source, f0
    implicit none
    private
    public :: solovev_frozen_t

    type, extends(ideal_frozen_t) :: solovev_frozen_t
    contains
        procedure :: evaluate
    end type
contains
    subroutine evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(solovev_frozen_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure

        chi = psi(point(1), point(2))
        chi_gradient = grad_psi(point(1), point(2))
        b3 = f0/point(1)
        pressure = source*(psi_b - chi)
    end subroutine
end module solovev_inverse_check_field

program test_equilibrium_map_inverse
    ! Compose independently solved forward and inverse maps on Solov'ev flux surfaces.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_grad_shafranov, only: gs_state_t, initialize_gs
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_curved_triangle, only: triangle_locate
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve_report_t, ideal_map_init, &
                             ideal_map_freeze, ideal_map_solve, ideal_map_point
    use milca_map_3d_equilibrium, only: map_3d_sample
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BPHI, EQ_BZ
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_solve
    use flux_map_sections, only: fm3_case_t, fm3_build, fm3_start, fm3_state
    use solovev_label_fields, only: psi, psi_b, r0, ray_radius, boundary_samples, &
                                   solovev_profiles, toroidal_flux, exact_b
    use solovev_inverse_check_field, only: solovev_frozen_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), poloidal_levels(3) = [0.2_dp, 0.5_dp, 0.8_dp]
    integer, parameter :: angles = 64, planes = 3
    integer :: level, levels, forward_steps, inverse_steps
    real(dp) :: errors(2, 3), order(2), field_errors(2)
    character(8) :: arg

    levels = 3
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (levels < 1 .or. levels > 3) error stop 'refinement levels must be between one and three'
    clebsch_nfp = 1
    print '(a)', ' cells   nodes   N_fwd   N_inv  RMS toroidal s  RMS angle (rad)  relative field errors (fwd, inv)'
    do level = 1, levels
        call comparison(2**level, 6*2**level, errors(:, level), field_errors, forward_steps, inverse_steps)
        print '(4i8,4es16.6)', 2**level, 6*2**level, forward_steps, inverse_steps, errors(:, level), field_errors
    end do
    if (levels == 3) then
        order = log(errors(:, 2)/errors(:, 3))/log(2.0_dp)
        print '(a,2f8.3)', 'last refinement orders (flux, angle)', order
        if (.not. all(order >= 1.7_dp)) error stop 'inverse composition refinement below quadratic field order'
        if (.not. all(field_errors <= 0.01_dp)) error stop 'physical fields must agree with the manufactured field'
    end if
contains
    subroutine comparison(cells, nodes, error, field_errors, forward_steps, inverse_steps)
        integer, intent(in) :: cells, nodes
        real(dp), intent(out) :: error(2), field_errors(2)
        integer, intent(out) :: forward_steps, inverse_steps
        type(cylindrical_chart_t) :: cylindrical
        type(solovev_frozen_t) :: frozen
        type(gs_state_t) :: mesh
        type(ideal_map_t) :: forward
        type(ideal_map_solve_report_t) :: report
        type(fm3_case_t) :: inverse
        real(dp), allocatable :: wall(:, :)
        real(dp) :: contour(2, nodes), point(2), xi(2), mapped(3), ref(2), rz(2), u(8), b(3)
        real(dp) :: step, expected_s, total_flux, phase(angles), exact_angle(0:angles), weights(4), q(4)
        real(dp) :: chi, delta, offset, t, flux_error, angle_error
        real(dp) :: values(EQ_VALUES), expected_b(3), field_scale
        integer :: i, j, cell, block, status

        do i = 1, nodes
            call boundary(real(i - 1, dp)/nodes, contour(:, i))
        end do
        step = 2*pi*sqrt(psi_b)/(r0*nodes)
        call initialize_gs(mesh, contour, step**2/2, status, 2, curve=boundary, geometry_only=.true.)
        if (status /= 0) error stop 'forward reference mesh'
        call ideal_map_init(forward, mesh, cylindrical, 1, 0, status, quadrature_order=7)
        if (status /= 0) error stop 'forward map initialisation'
        call ideal_map_freeze(forward, frozen, status)
        if (status /= 0) error stop 'forward frozen field'
        do i = 1, forward%nodes
            ref = forward%reference(:, i)
            if (forward%boundary(i)) cycle
            delta = 0.02_dp*(1 - psi(ref(1), ref(2))/psi_b)
            forward%x(:2, 1, i) = forward%x(:2, 1, i) + delta*forward%direction(:, i)
        end do
        call ideal_map_solve(forward, 1e-8_dp, 80, report, status)
        if (status /= 0) error stop 'forward equilibrium'
        forward_steps = report%iterations

        call boundary_samples(wall)
        call fm3_build(inverse, spread(wall, 3, planes), cells, 2, 4, status)
        if (status /= 0) error stop 'inverse reference mesh'
        call fm3_start(inverse)
        call flux_map3d_solve(inverse%fm, inverse%derham, inverse%quad, profiles, inverse%base, 1e-9_dp, 60, status)
        if (status /= 0) error stop 'inverse equilibrium'
        inverse_steps = inverse%fm%steps

        block = 1
        error = 0
        field_errors = 0
        field_scale = 0
        total_flux = toroidal_flux(1.0_dp)
        step = 2*pi/angles
        call gauss_legendre_ab(4, 0.0_dp, step, q, weights)
        do j = 1, size(poloidal_levels)
            chi = poloidal_levels(j)
            expected_s = toroidal_flux(chi)/total_flux
            ! psi = ((R**2-r0**2)/2)**2 + (R*Z)**2. The straight angle has
            ! density R**(-3) in the polar angle of these two functions.
            exact_angle(0) = 0
            do i = 1, angles
                t = (i - 1)*step
                exact_angle(i) = exact_angle(i - 1) + sum(weights/radius(chi, t + q)**3)
            end do
            exact_angle = exact_angle*(2*pi/exact_angle(angles))
            flux_error = 0
            do i = 1, angles
                t = (i - 1)*step
                point(1) = radius(chi, t)
                point(2) = sqrt(psi_b*chi)*sin(t)/point(1)
                do cell = 1, forward%cells
                    xi = matmul(transpose(mesh%gradient(:, 2:3, cell)), point - &
                                forward%reference(:, forward%conn(1, cell)))
                    if (minval([xi, 1 - sum(xi)]) < -0.15_dp) cycle
                    call triangle_locate(forward%reference(:, forward%conn(:, cell)), point, xi, status)
                    if (status == 0) exit
                end do
                if (cell > forward%cells) error stop 'exact contour point outside forward reference mesh'
                call ideal_map_point(forward, cell, xi, 0.0_dp, mapped, ref)
                if (norm2(ref - point) > 1e-9_dp) error stop 'forward reference lookup'
                call fm3_state(inverse, mapped(:2), block, mapped(3)/(2*pi), profiles, rz, u, b, status)
                if (status /= 0) error stop 'forward image outside inverse mesh'
                if (norm2(rz - mapped(:2)) > 1e-9_dp) error stop 'inverse physical lookup'
                flux_error = flux_error + (toroidal_flux(sum(u(:2)**2))/total_flux - expected_s)**2
                phase(i) = atan2(u(2), u(1)) - exact_angle(i - 1)
                call map_3d_sample(forward, frozen, cell, xi, 0.0_dp, mapped, values, status)
                if (status /= 0) error stop 'forward field evaluation'
                expected_b = exact_b([rz(1), 0.0_dp, rz(2)])
                field_errors(1) = field_errors(1) + sum((values([EQ_BR, EQ_BPHI, EQ_BZ]) - expected_b)**2)
                field_errors(2) = field_errors(2) + sum((b - expected_b)**2)
                field_scale = field_scale + sum(expected_b**2)
            end do
            offset = atan2(sum(sin(phase)), sum(cos(phase)))
            angle_error = sum((modulo(phase - offset + pi, 2*pi) - pi)**2)
            error = error + [flux_error, angle_error]
        end do
        error = sqrt(error/(angles*size(poloidal_levels)))
        field_errors = sqrt(field_errors/field_scale)
    end subroutine comparison

    elemental real(dp) function radius(chi, t)
        real(dp), intent(in) :: chi, t
        radius = sqrt(r0**2 + 2*sqrt(chi*psi_b)*cos(t))
    end function radius

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*pi*t
        point = [r0, 0.0_dp] + ray_radius(psi_b, angle)*[cos(angle), sin(angle)]
    end subroutine boundary

    subroutine profiles(chi, tau, gp, p)
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        call solovev_profiles(chi, tau, gp, p)
        tau = -2*tau
    end subroutine profiles
end program test_equilibrium_map_inverse
