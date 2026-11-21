program test_eulerian_clebsch_landreman
    !! Gate E4a of docs/eulerian-label-plan.md: the 3D label lane on an exact
    !! nested 3D equilibrium (Landreman integer family, docs/nested-oracles.md:
    !! two field periods, psi = rho^2, iota = 2). The section of every
    !! toroidal plane phi_k is meshed by an O-grid centred on the exact axis
    !! X_k and fitted to the exact boundary psi = delta; the planes give the
    !! node coordinates s_k of the section map. chi = psi/delta is projected
    !! per plane and held fixed; g solves J.grad chi = 0 (the g equations of
    !! W). The winding angle is taken in the frame H^(1/2) of chi at X_k. The
    !! L2 error of B against the exact field must converge at order p.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_3d_system
    use milca_nested_oracles, only: nested_oracle_t, NESTED_INTEGER
    use milca_nested_oracle_leaves, only: nested_integer_surface, nested_ray, &
        nested_integer_flux
    implicit none
    integer, parameter :: degree = 2, samples_n = 128
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(nested_oracle_t) :: o
    type(spline_space_t) :: space
    integer :: nt = 13, levels = 3, level, failures, component, plane
    integer :: cells(3) = [2, 4, 8]
    real(dp) :: err(3), rate, plane_phi, sign_f = 1, sign_g = -1
    real(dp) :: ridge = 1e-12_dp
    character(16) :: arg

    failures = 0
    o%family = NESTED_INTEGER
    o%eps = 0.1_dp
    o%tau = 0
    o%delta = 1.0_dp/64
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) o%eps
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) sign_f
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) sign_g
    end if
    print '(a,f6.3,a,i3,a,2f5.1)', 'eps', o%eps, '  planes', nt, '  signs F, G', sign_f, sign_g
    do level = 1, levels
        err(level) = solve_g(cells(level))
        print '(a,i3,a,es10.3)', 'cells', cells(level), '  L2 error of B', err(level)
    end do
    do level = 2, levels
        rate = log(err(level - 1)/err(level))/log(2.0_dp)
        print '(a,f6.2)', 'rate', rate
        call check(rate > degree - 0.3_dp, 'L2 error of B converges at order p')
    end do
    if (failures > 0) error stop 'Landreman gate E4a failed'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            print '(2a)', 'PASS ', what
        else
            print '(2a)', 'FAIL ', what
            failures = failures + 1
        end if
    end subroutine check

    function axis_at(phi) result(x)
        !! (R, Z) of the magnetic axis in the plane phi: the surface map at
        !! rho = 0, its parameter t found where the ray angle matches phi.
        real(dp), intent(in) :: phi
        real(dp) :: x(2), lo, hi, t, c(3), jac, ray, radius
        integer :: k

        lo = phi - pi/2
        hi = phi + pi/2
        do k = 1, 60
            t = (lo + hi)/2
            call nested_integer_surface(0.0_dp, 0.0_dp, t, o%eps, o%tau, c(1), c(2), c(3), jac)
            call nested_ray(c(1), c(2), phi, ray, radius)
            if (ray < 0) then
                lo = t
            else
                hi = t
            end if
        end do
        x = [sqrt(c(1)**2 + c(2)**2), c(3)]
    end function axis_at

    real(dp) function chi_at(r, phi, z) result(chi)
        real(dp), intent(in) :: r, phi, z
        real(dp) :: v(20)
        integer :: st

        call o%point([r*cos(phi), r*sin(phi), z], v, st)
        if (st /= 0) error stop 'oracle point outside the domain'
        chi = v(1)/o%delta
    end function chi_at

    subroutine exact_b(x, b)
        !! Cylindrical (R, phi, Z) components of the exact field.
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: v(20)
        integer :: st

        call o%point([x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)], v, st)
        if (st /= 0) error stop 'exact field outside the domain'
        b = [v(3)*cos(x(2)) + v(4)*sin(x(2)), -v(3)*sin(x(2)) + v(4)*cos(x(2)), v(5)]
    end subroutine exact_b

    subroutine profiles(chi, f, g, p)
        !! F = -sign_f Phi(psi)/(2 pi), G = sign_g iota Phi/(2 pi), psi = delta chi,
        !! Phi from the closed form; p is not used by the g equations. The
        !! oracle's poloidal angle runs opposite to h_theta's, hence sign_g = -1:
        !! the other three sign pairs miss B by 35 % to 200 % at 2 cells.
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: f(3), g(3), p(3)
        real(dp) :: d(3)

        d = flux_derivatives(chi)
        f = -sign_f*d/(2*pi)
        g = sign_g*2*d/(2*pi)
        p = 0
    end subroutine profiles

    function flux_derivatives(chi) result(d)
        !! dPhi/dchi and two more derivatives by central differences.
        real(dp), intent(in) :: chi
        real(dp) :: d(3), fm2, fm1, f0, fp1, fp2
        real(dp), parameter :: h = 1e-3_dp

        fm2 = phi_of(chi - 2*h)
        fm1 = phi_of(chi - h)
        f0 = phi_of(chi)
        fp1 = phi_of(chi + h)
        fp2 = phi_of(chi + 2*h)
        d(1) = (fp1 - fm1)/(2*h)
        d(2) = (fp1 - 2*f0 + fm1)/h**2
        d(3) = (fp2 - 2*fp1 + 2*fm1 - fm2)/(2*h**3)
    end function flux_derivatives

    real(dp) function phi_of(chi)
        real(dp), intent(in) :: chi
        real(dp) :: psi

        psi = o%delta*chi
        call nested_integer_flux(psi, o%eps, o%tau, phi_of)
    end function phi_of

    subroutine coordinate(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(space, block, xi, jet, st)
        values(1) = jet(component, 1)
        associate (unused => t)
        end associate
    end subroutine coordinate

    subroutine label_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(space, block, xi, jet, st)
        values(1) = chi_at(jet(1, 1), plane_phi, jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

    subroutine plane_layout(phi, centre, layout)
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: centre(2)
        type(block_layout_t), intent(out) :: layout
        real(dp) :: samples(2, samples_n)
        integer :: i, st

        centre = axis_at(phi)
        ! rays from the axis: counter-clockwise from the outboard midplane
        do i = 1, samples_n
            samples(:, i) = centre + boundary_radius(phi, centre, 2*pi*(i - 1)/samples_n) &
                            *[cos(2*pi*(i - 1)/samples_n), sin(2*pi*(i - 1)/samples_n)]
        end do
        call layout_ogrid(centre, 0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :))), &
                          0.5_dp, layout, st)
        if (st /= 0) error stop 'layout'
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'boundary curve'
    end subroutine plane_layout

    real(dp) function boundary_radius(phi, centre, angle) result(r)
        !! Distance from the axis to chi = 1 along a ray in the plane phi.
        real(dp), intent(in) :: phi, centre(2), angle
        real(dp) :: lo, hi, e(2)
        integer :: k

        e = [cos(angle), sin(angle)]
        lo = 0
        hi = 0.05_dp
        do while (chi_at(centre(1) + hi*e(1), phi, centre(2) + hi*e(2)) < 1)
            lo = hi
            hi = 1.5_dp*hi
        end do
        do k = 1, 60
            r = (lo + hi)/2
            if (chi_at(centre(1) + r*e(1), phi, centre(2) + r*e(2)) < 1) then
                lo = r
            else
                hi = r
            end if
        end do
    end function boundary_radius

    function sqrt_hessian(phi, x) result(m)
        !! H^(1/2) of the exact chi in the plane phi at x.
        real(dp), intent(in) :: phi, x(2)
        real(dp) :: m(2, 2), hm(2, 2), sd
        real(dp), parameter :: h = 1e-3_dp

        hm(1, 1) = (chi_at(x(1) + h, phi, x(2)) - 2*chi_at(x(1), phi, x(2)) &
                    + chi_at(x(1) - h, phi, x(2)))/h**2
        hm(2, 2) = (chi_at(x(1), phi, x(2) + h) - 2*chi_at(x(1), phi, x(2)) &
                    + chi_at(x(1), phi, x(2) - h))/h**2
        hm(1, 2) = (chi_at(x(1) + h, phi, x(2) + h) - chi_at(x(1) + h, phi, x(2) - h) &
                    - chi_at(x(1) - h, phi, x(2) + h) + chi_at(x(1) - h, phi, x(2) - h))/(4*h**2)
        hm(2, 1) = hm(1, 2)
        sd = sqrt(hm(1, 1)*hm(2, 2) - hm(1, 2)**2)
        m = hm
        m(1, 1) = m(1, 1) + sd
        m(2, 2) = m(2, 2) + sd
        m = m/sqrt(hm(1, 1) + hm(2, 2) + 2*sd)
    end function sqrt_hessian

    real(dp) function solve_g(n) result(rel)
        integer, intent(in) :: n
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad, aquad
        type(spline_derham_t) :: derham
        type(csr_t) :: hessian
        real(dp), allocatable :: nodes(:, :, :), chi(:), g(:), grad(:), delta(:), values(:)
        real(dp) :: axes(2, nt), frames(2, 2, nt), energy, ferr(2), centre(2)
        integer :: st, n0, ntot, k, i
        integer, allocatable :: free(:)

        do plane = nt, 1, -1
            plane_phi = 2*pi*(plane - 1)/nt
            call plane_layout(plane_phi, centre, layout)
            axes(:, plane) = centre
            frames(:, :, plane) = sqrt_hessian(plane_phi, centre)
            if (plane == nt) then
                options%degree = degree
                options%regularity = SPLINE_C0
                allocate (options%knots(2, layout%block_count))
                options%knots(:, :)%element_count = n
                call spline_space_build(layout, options, space, st)
                if (st /= 0) error stop 'space'
            else
                call spline_space_update_geometry(space, layout, st)
                if (st /= 0) error stop 'plane geometry'
            end if
            call spline_derham_build(space, derham, st, 0)
            if (st /= 0) error stop 'derham'
            n0 = derham%dimension(0)
            if (.not. allocated(nodes)) allocate (nodes(2, n0, nt), chi(n0*nt), values(n0))
            do component = 1, 2
                call spline_derham_project(derham, 0, coordinate, values, st, 12)
                nodes(component, :, plane) = values
            end do
            call spline_derham_project(derham, 0, label_cb, values, st, 12)
            chi((plane - 1)*n0 + 1:plane*n0) = values
        end do
        ! plane 1 is the reference section of the quadrature
        call spline_quadrature_build(space, degree + 8, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_quadrature_point_singular(space, quad, axes(:, 1), aquad, st)
        if (st /= 0) error stop 'axis quadrature'
        ntot = n0*nt
        allocate (g(ntot))
        g = 0
        call clebsch_3d_system(derham, aquad, nt, profiles, axes(:, 1), chi, g, energy, grad, &
                               hessian, plane_axis=axes, plane_nodes=nodes, plane_frames=frames)
        call add_g_ridge(hessian, ntot)
        allocate (delta(2*ntot))
        free = [(i, i=1, ntot)]
        call eu_constrained_solve(hessian, -grad, free, 0*chi, delta, st)
        if (st /= 0) error stop 'g solve'
        g = delta(ntot + 1:)
        call clebsch_3d_system(derham, aquad, nt, profiles, axes(:, 1), chi, g, energy, grad, &
                               plane_axis=axes, plane_nodes=nodes, plane_frames=frames, &
                               exact_field=exact_b, field_error=ferr)
        rel = sqrt(ferr(1)/ferr(2))
        associate (unused => k)
        end associate
    end function solve_g

    subroutine add_g_ridge(hessian, n0)
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: n0
        real(dp) :: shift
        integer :: i, p

        shift = ridge*maxval(abs(hessian%val))
        do i = n0 + 1, 2*n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine add_g_ridge

end program test_eulerian_clebsch_landreman
