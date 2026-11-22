program test_label_ale_axisym
    !! Gate A1 of docs/equilibrium-roadmap.md: the axisymmetric Solov'ev
    !! equilibrium from the label lane on an ALE mesh. The O-grid starts with
    !! its core block off the magnetic axis (offset dR, dZ in units of the
    !! minor radius, default 0.1, 0.05); labels start 5 % off the exact ones
    !! (E1's perturbation) and X at the label minimum. milca_label_ale
    !! alternates Newton on (chi, g, X, mu), with the core following X, and
    !! the mesh update with the transfer of chi and g.
    !!
    !! Passed when the mesh loop ends in at most 2 passes, B is within 10 %
    !! of the interpolation error on 2, 4 and 8 cells (the bound of E1, whose
    !! mesh is tied to the exact axis) and converges at order p. The
    !! Winslow mesh update (2 and 4 cells) must converge in 2 passes and at
    !! order p; its B error is reported against the blend.
    !!
    !! Arguments: dR dZ [levels].
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0, &
        spline_space_update_geometry
    use fortfem_spline_eval, only: spline_quadrature_t, spline_locate
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_project, &
        spline_derham_evaluate
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_label_ale, only: label_ale_t, label_ale_init, label_ale_solve, &
        label_ale_initial_g
    use milca_spline_ale, only: spline_min_jacobian
    use fortfem_spline_space, only: spline_space_t
    use solovev_label_fields
    implicit none
    integer, parameter :: degree = 2, levels = 3
    real(dp), parameter :: interp(3) = [2.196e-2_dp, 4.367e-3_dp, 9.017e-4_dp]
    real(dp) :: amplitude = 0.05_dp, offset(2) = [0.1_dp, 0.05_dp], mesh_tol = 0.05_dp
    integer :: cells(levels) = [2, 4, 8], level, failures, nlevels
    real(dp) :: error(levels), werror(levels), minor
    character(32) :: arg

    failures = 0
    nlevels = levels
    if (command_argument_count() >= 2) then
        call get_command_argument(1, arg)
        read (arg, *) offset(1)
        call get_command_argument(2, arg)
        read (arg, *) offset(2)
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) nlevels
    end if
    call curve_chart_cylindrical(chart, 1)
    do level = 1, nlevels
        call run(cells(level), .false., error(level))
    end do
    do level = 2, nlevels
        print '(a,f6.2)', 'rate of B', log(error(level - 1)/error(level))/log(2.0_dp)
        call check(log(error(level - 1)/error(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'L2 error of B converges at order p')
    end do
    ! Winslow between core and wall
    do level = 1, min(2, nlevels)
        call run(cells(level), .true., werror(level))
        print '(a,f6.2)', '  Winslow / blend error of B', werror(level)/error(level)
    end do
    if (nlevels >= 2) then
        print '(a,f6.2)', 'rate of B (Winslow)', log(werror(1)/werror(2))/log(2.0_dp)
        call check(log(werror(1)/werror(2))/log(2.0_dp) > degree - 0.3_dp, &
                   'Winslow mesh: L2 error of B converges at order p')
    end if
    if (failures > 0) error stop 'gate A1 failed'
    print '(a)', 'PASS gate A1'

contains

    subroutine run(n, winslow, err)
        integer, intent(in) :: n
        logical, intent(in) :: winslow
        real(dp), intent(out) :: err
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(label_ale_t) :: s
        real(dp), allocatable :: samples(:, :)
        real(dp) :: shift(2), t0, t1
        integer :: st, last

        call cpu_time(t0)
        call boundary_samples(samples)
        minor = 0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :)))
        call layout_ogrid([r0, 0.0_dp], minor, 0.5_dp, layout, st)
        if (st /= 0) error stop 'layout'
        call layout_set_curve(layout, LAYOUT_WALL, samples, st)
        if (st /= 0) error stop 'wall'
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        ! the starting mesh: core block moved off the axis
        shift = offset*minor
        layout%vertices(:, 1:4) = layout%vertices(:, 1:4) + spread(shift, 2, 4)
        call spline_space_update_geometry(space, layout, st)
        if (st /= 0) error stop 'shifted mesh'
        call label_ale_init(s, space, 1, [0.5_dp, 0.5_dp], 8, st)
        if (st /= 0) error stop 'label_ale_init'
        allocate (s%chi(s%derham%dimension(0)), s%g(s%derham%dimension(0)))
        call spline_derham_project(s%derham, 0, perturbed_label, s%chi, st, 12)
        if (st /= 0) error stop 'labels'
        s%chi(s%wall) = 1
        s%g = 0
        s%axis = label_minimum(s)
        print '(a,i2,a,l1,a,2f9.5,a,2f9.5)', 'cells', n, ' Winslow ', winslow, &
            ': core shifted by', shift, '  label minimum', s%axis
        call label_ale_init_g(s)
        call label_ale_solve(s, space, solovev_profiles, 1e-9_dp, mesh_tol, 8, 12, st, winslow)
        call cpu_time(t1)
        last = size(s%newton_history)
        print '(a,i0,a,i0,a,i0,a,i0)', '  status ', st, ', passes ', s%passes, &
            ', Newton steps ', s%newton_total, ', Winslow steps ', s%winslow_total
        print '(a,14es9.1)', '  |r| history (first 14)', s%newton_history(1:min(last, 14))
        print '(a,2es11.3,a,es10.2)', '  axis - exact axis', s%axis - [r0, 0.0_dp], &
            '  min det of the mesh', spline_min_jacobian(space)
        space = s%derham%space
        err = field_error(s%derham, s%aquad, s%chi, s%g, s%axis)
        print '(a,es10.3,a,es10.3,a,f6.2,a)', '  L2 error of B', err, ' (interpolation ', &
            interp(level), '), ', t1 - t0, ' s'
        call check(st == 0, 'ALE loop converges')
        call check(s%passes <= 2, 'at most 2 mesh passes')
        if (.not. winslow) call check(err < 1.1_dp*interp(level), &
                                      'B within 10 % of the interpolation error')
    end subroutine run

    subroutine label_ale_init_g(s)
        type(label_ale_t), intent(inout) :: s
        integer :: st
        call label_ale_initial_g(s, solovev_profiles, st)
        if (st /= 0) error stop 'initial g'
    end subroutine label_ale_init_g

    function label_minimum(s) result(x)
        !! Minimum of chi_h over the quadrature points, then Newton on
        !! grad chi_h = 0 by central differences.
        type(label_ale_t), intent(in) :: s
        real(dp) :: x(2), best, v, h(2, 2), gr(2), e(2)
        real(dp), parameter :: d = 1e-4_dp
        integer :: p, it, k
        best = huge(1.0_dp)
        do p = 1, s%quad%point_count
            v = chi_at(s, s%quad%s(:, p))
            if (v < best) then
                best = v
                x = s%quad%s(:, p)
            end if
        end do
        do it = 1, 20
            do k = 1, 2
                e = 0
                e(k) = d
                gr(k) = (chi_at(s, x + e) - chi_at(s, x - e))/(2*d)
                h(k, k) = (chi_at(s, x + e) - 2*chi_at(s, x) + chi_at(s, x - e))/d**2
            end do
            h(1, 2) = (chi_at(s, x + [d, d]) - chi_at(s, x + [d, -d]) - chi_at(s, x + [-d, d]) + &
                       chi_at(s, x - [d, d]))/(4*d*d)
            h(2, 1) = h(1, 2)
            e = [h(2, 2)*gr(1) - h(1, 2)*gr(2), h(1, 1)*gr(2) - h(2, 1)*gr(1)]/ &
                (h(1, 1)*h(2, 2) - h(1, 2)**2)
            x = x - e
            if (norm2(e) < 1e-12_dp) exit
        end do
    end function label_minimum

    real(dp) function chi_at(s, x) result(v)
        type(label_ale_t), intent(in) :: s
        real(dp), intent(in) :: x(2)
        real(dp) :: xi(2), val(1)
        integer :: block, st
        call spline_locate(s%derham%space, 1, x, block, xi, st)
        call spline_derham_evaluate(s%derham, 0, s%chi, block, xi, 0.0_dp, val)
        v = val(1)
    end function chi_at

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

    subroutine perturbed_label(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), js(2, 2), pn
        call section_jet(block, xi, s, js)
        pn = psi(s(1), s(2))/psi_b
        values(1) = pn + amplitude*pn*(1 - pn)*(s(1) - r0 + s(2))
        associate (unused => t)
        end associate
    end subroutine perturbed_label

    real(dp) function field_error(derham, quad, chi, g, axis) result(rel)
        !! Relative L2 error of B against the exact Solov'ev field (E1).
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: chi(:), g(:), axis(2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :)
        integer, allocatable :: d0(:)
        real(dp) :: s(2), js(2, 2), x(3), j3(3, 3), c, gc(2), gg(2), th(2), f(3), gp(3), &
                    pr(3), bphi, bpol(2), b_ex(3), w, num, den
        integer :: cell, q, p
        num = 0
        den = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call section_jet(quad%block(p), quad%xi(:, p), s, js)
                call full(quad%block(p), quad%xi(:, p), 0.0_dp, x, j3)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                gg = matmul(gn(:, :, q), g(d0))
                th = [-(s(2) - axis(2)), s(1) - axis(1)]/((s(1) - axis(1))**2 + &
                                                            (s(2) - axis(2))**2)
                call solovev_profiles(c, f, gp, pr)
                bphi = gg(2)*gc(1) - gg(1)*gc(2) - f(1)*(gc(1)*th(2) - gc(2)*th(1))
                bpol = gp(1)*[-gc(2), gc(1)]/s(1)
                b_ex = exact_b(x)
                w = 2*pi*s(1)*quad%weight(p)
                num = num + w*((bphi - b_ex(2))**2 + (bpol(1) - b_ex(1))**2 + &
                               (bpol(2) - b_ex(3))**2)
                den = den + w*sum(b_ex**2)
            end do
        end do
        rel = sqrt(num/den)
    end function field_error

end program test_label_ale_axisym
