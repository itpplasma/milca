program test_eulerian_clebsch_landreman_newton
    !! Gate E4 of docs/eulerian-label-plan.md: a nonaxisymmetric equilibrium
    !! of the 3D label lane by Newton, with the labels chi, g and the magnetic
    !! axis curve (one point X_k per toroidal plane) as unknowns, on the exact
    !! sheared Landreman family (docs/nested-oracles.md, two field periods,
    !! magnetic shear). The integer family is not used: every surface there is
    !! rational and shearless, and its g block is nearly singular (E4a).
    !!
    !! Label chi = psi/delta (1 on the boundary), profiles of the exact state:
    !!
    !!   F(chi) = -Phi(psi)/(2 pi),  G'(chi) = -iota(psi) Phi'(chi)/(2 pi)
    !!   p(chi) = (delta - psi)/lambda^2,   psi = delta chi,
    !!
    !! (the sign pair that E4a selects) with Phi and iota Chebyshev series fitted to the oracle's flux
    !! and rotation on [0, delta], so the residual is smooth in chi (a
    !! piecewise profile costs quadratic convergence, E3b). Each plane is an
    !! O-grid fitted to its exact boundary (test/support/landreman_sections);
    !! its core follows X_k rigidly (ALE), so X_k stays the central vertex,
    !! the polar rule of the reference plane serves every plane and the rows
    !! N(X_k), grad N(X_k) of each plane stay fixed. Newton on
    !!
    !!   grad W + mu_k [N(X_k); 0] = 0,  grad chi_k(X_k) = 0,  chi_k(X_k) = 0,
    !!
    !! (W = int |B|^2/2 - int p) from labels perturbed per plane and a helical
    !! axis offset. The winding angle is taken in the frame H^(1/2) of the
    !! exact label at the exact axis of each plane, held fixed. B is compared
    !! with the exact field and with the fixed-label fit of E4a on the same
    !! mesh (the interpolation level), the labels and the axis with the exact
    !! ones.
    !!
    !! Arguments: levels (default 2: 2, 4 cells per block side), planes
    !! (default 13), axis offset (default 2e-3), label perturbation
    !! (default 0.02), eps (default 0.6), 'jacobian' to check the Newton
    !! matrix against differences of the residual, S (default 2.4), delta
    !! (default 0.2).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_3d_system
    use milca_nested_oracles, only: NESTED_SHEARED
    use solovev_label_fields, only: triplet_list
    use landreman_sections, only: section_oracle, section_label, section_exact_b, &
        section_sqrt_hessian, section_layout, section_fit_profiles, section_profiles
    implicit none
    integer, parameter :: degree = 2
    ! ridge on the g block: of the fixed-label g solve, and of the Newton
    ! matrix only (not the residual), where it damps the discrete near-null
    ! g -> g + k(chi_h) without changing the solution
    real(dp), parameter :: pi = acos(-1.0_dp), ridge = 1e-12_dp, newton_ridge = 1e-8_dp
    integer :: nt = 13, levels = 2, max_steps = 16, level, failures, component
    integer :: cells(3) = [2, 4, 8]
    real(dp) :: axis_offset = 2e-3_dp, amplitude = 0.02_dp, plane_phi
    logical :: jacobian_check = .false.
    real(dp), allocatable :: err(:), err_fixed(:), chi_err(:), axis_err(:)
    integer, allocatable :: steps(:)
    type(spline_space_t) :: space
    character(16) :: arg

    failures = 0
    ! the gate case: iota 4.73 to 5.11 (shear 0.38), the axis moves from
    ! R = 1.456 to 1.646 over a field period
    section_oracle%family = NESTED_SHEARED
    section_oracle%eps = 0.6_dp
    section_oracle%saxis = 2.4_dp
    section_oracle%delta = 0.2_dp
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) axis_offset
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) amplitude
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) section_oracle%eps
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        jacobian_check = arg == 'jacobian'
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) section_oracle%saxis
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, arg)
        read (arg, *) section_oracle%delta
    end if
    if (.not. section_oracle%valid()) error stop 'sheared oracle parameters'
    call section_fit_profiles()
    print '(a,f5.2,a,f5.2,a,f7.4,a,i3,a,es9.2,a,es9.2)', 'sheared family eps', section_oracle%eps, &
        '  S', section_oracle%saxis, '  delta', section_oracle%delta, &
        '  planes', nt, '  axis offset', axis_offset, '  label perturbation', amplitude
    allocate (err(levels), err_fixed(levels), chi_err(levels), axis_err(levels), steps(levels))
    do level = 1, levels
        call solve(cells(level), err(level), err_fixed(level), chi_err(level), &
                   axis_err(level), steps(level))
        print '(a,i3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,i3)', 'cells', cells(level), &
            '  B error', err(level), '  fixed-label fit', err_fixed(level), &
            '  max |chi - Pi0 chi|', chi_err(level), '  axis', axis_err(level), &
            '  Newton steps', steps(level)
        call check(steps(level) <= 10, 'Newton converges in at most 10 steps')
        call check(err(level) < 1.2_dp*err_fixed(level), &
                   'B within 20 % of the fixed-label fit on the same mesh')
    end do
    do level = 2, levels
        print '(a,3f7.2)', 'rates of B, chi, axis', &
            log(err(level - 1)/err(level))/log(2.0_dp), &
            log(chi_err(level - 1)/chi_err(level))/log(2.0_dp), &
            log(axis_err(level - 1)/axis_err(level))/log(2.0_dp)
        call check(log(err(level - 1)/err(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'B converges at order p')
        call check(chi_err(level) < chi_err(level - 1)/2, 'the labels converge')
        call check(axis_err(level) < axis_err(level - 1)/2, 'the axis curve converges')
    end do
    if (failures > 0) error stop 'Landreman Newton gate E4 failed'

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
        values(1) = section_label(jet(1, 1), plane_phi, jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

    subroutine plane_geometry(k, n, derham, core)
        !! The space on plane k (the boundary of that plane), with its core
        !! centred at core when present, else on the exact axis.
        integer, intent(in) :: k, n
        type(spline_derham_t), intent(out) :: derham
        real(dp), intent(in), optional :: core(2)
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        real(dp) :: centre(2)
        integer :: st
        logical :: rebuild

        plane_phi = 2*pi*(k - 1)/nt
        if (present(core)) then
            call section_layout(plane_phi, centre, layout, core)
        else
            call section_layout(plane_phi, centre, layout)
        end if
        rebuild = .true.
        ! (no short-circuit .or. in Fortran: test the allocation first)
        if (allocated(space%blocks)) rebuild = space%blocks(1)%element_count(1) /= n
        if (rebuild) then
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
    end subroutine plane_geometry

    function projected(derham, which) result(v)
        !! V0 node values of R (1), Z (2) or the exact label (3) on the
        !! current geometry.
        type(spline_derham_t), intent(in) :: derham
        integer, intent(in) :: which
        real(dp) :: v(derham%dimension(0))
        integer :: st

        if (which <= 2) then
            component = which
            call spline_derham_project(derham, 0, coordinate, v, st, 12)
        else
            call spline_derham_project(derham, 0, label_cb, v, st, 12)
        end if
        if (st /= 0) error stop 'projection'
    end function projected

    function axis_values(derham, x) result(row)
        !! Values of the V0 node functions at the physical point x.
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: x(2)
        real(dp), allocatable :: row(:), unit(:)
        real(dp) :: xi(2), v(1)
        integer :: block, st, i

        call spline_locate(space, 1, x, block, xi, st)
        if (st /= 0) error stop 'axis location'
        allocate (row(derham%dimension(0)), unit(derham%dimension(0)))
        unit = 0
        do i = 1, size(row)
            unit(i) = 1
            call spline_derham_evaluate(derham, 0, unit, block, xi, 0.0_dp, v)
            row(i) = v(1)
            unit(i) = 0
        end do
    end function axis_values

    function nodes_at(base, w, a, x) result(nodes)
        !! Node coordinates per plane of the mesh whose core follows X_k.
        real(dp), intent(in) :: base(:, :, :), w(:), a(:, :), x(:, :)
        real(dp) :: nodes(2, size(w), nt)
        integer :: k

        do k = 1, nt
            nodes(1, :, k) = base(1, :, k) + w*(x(1, k) - a(1, k))
            nodes(2, :, k) = base(2, :, k) + w*(x(2, k) - a(2, k))
        end do
    end function nodes_at

    subroutine residual(derham, quad, a, base, w, frames, nrow, ngrad, wall, chi, g, x, mu, &
                        r, hessian, ferr)
        !! Newton residual in (chi, g, X_k, mu_k), see the header.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: a(:, :), base(:, :, :), w(:), frames(:, :, :), nrow(:), &
                                ngrad(:, :, :), chi(:), g(:), x(:, :), mu(:)
        integer, intent(in) :: wall(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), intent(out), optional :: ferr(2)
        real(dp), allocatable :: grad(:)
        real(dp) :: energy
        integer :: n0, k, c, nf

        n0 = size(nrow)
        nf = 2*n0*nt
        if (present(ferr)) then
            call clebsch_3d_system(derham, quad, nt, section_profiles, a(:, 1), chi, g, energy, grad, &
                                   hessian, x, nodes_at(base, w, a, x), frames, &
                                   section_exact_b, ferr)
        else
            call clebsch_3d_system(derham, quad, nt, section_profiles, a(:, 1), chi, g, energy, grad, &
                                   hessian, x, nodes_at(base, w, a, x), frames)
        end if
        allocate (r(nf + 3*nt))
        r(1:nf) = grad
        do k = 1, nt
            associate (chik => chi((k - 1)*n0 + 1:k*n0))
                r((k - 1)*n0 + 1:k*n0) = r((k - 1)*n0 + 1:k*n0) + mu(k)*nrow
                do c = 1, 2
                    r(nf + 2*(k - 1) + c) = dot_product(ngrad(:, c, k), chik)
                end do
                r(nf + 2*nt + k) = dot_product(nrow, chik)
            end associate
        end do
        r(wall) = 0
    end subroutine residual

    subroutine solve(n, err, err_fixed, chi_err, axis_err, steps)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, err_fixed, chi_err, axis_err
        integer, intent(out) :: steps
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad, aquad
        type(csr_t) :: hessian, jac_t
        real(dp), allocatable :: base(:, :, :), w(:), chi0(:), chi(:), g(:), grad(:), delta(:), &
                                 nrow(:), ngrad(:, :, :), r(:), rp(:), rm(:), dgdx(:, :), rt(:), &
                                 rs(:), zs(:), history(:), chi_exact(:)
        integer, allocatable :: wall1(:), wall(:)
        real(dp) :: a(2, nt), x(2, nt), xp(2, nt), frames(2, 2, nt), mu(nt), energy, ferr(2), &
                    alpha, tk, shape
        real(dp), parameter :: hx = 1e-6_dp, hn = 1e-5_dp, shift = 1e-3_dp
        integer :: n0, ntot, nf, k, j, it, i, st

        ! base meshes, exact labels, frames and axis rows per plane
        do k = nt, 1, -1
            call plane_geometry(k, n, derham)
            a(:, k) = section_centre(k)
            frames(:, :, k) = section_sqrt_hessian(2*pi*(k - 1)/nt, a(:, k))
            n0 = derham%dimension(0)
            if (.not. allocated(base)) then
                allocate (base(2, n0, nt), chi0(n0*nt), ngrad(n0, 2, nt))
            end if
            base(1, :, k) = projected(derham, 1)
            base(2, :, k) = projected(derham, 2)
            chi0((k - 1)*n0 + 1:k*n0) = projected(derham, 3)
            ngrad(:, 1, k) = (axis_values(derham, a(:, k) + [hn, 0.0_dp]) - &
                              axis_values(derham, a(:, k) - [hn, 0.0_dp]))/(2*hn)
            ngrad(:, 2, k) = (axis_values(derham, a(:, k) + [0.0_dp, hn]) - &
                              axis_values(derham, a(:, k) - [0.0_dp, hn]))/(2*hn)
            if (k == nt) nrow = axis_values(derham, a(:, k))
            if (maxval(abs(axis_values(derham, a(:, k)) - nrow)) > 1e-10_dp) &
                error stop 'the axis is not a fixed vertex of every plane'
        end do
        print '(a,2f8.4,a,2f8.4)', '  exact axis: R from', minval(a(1, :)), maxval(a(1, :)), &
            '  Z from', minval(a(2, :)), maxval(a(2, :))
        ! motion weights from plane 1: core moved by shift in R, then in Z
        call plane_geometry(1, n, derham, a(:, 1) + [shift, 0.0_dp])
        rs = projected(derham, 1)
        call plane_geometry(1, n, derham, a(:, 1) + [0.0_dp, shift])
        zs = projected(derham, 2)
        w = (rs - base(1, :, 1))/shift
        call check(maxval(abs((zs - base(2, :, 1))/shift - w)) < 1e-8_dp, &
                   'the core moves rigidly with the axis in R and in Z')
        ! the reference plane carries the quadrature and the polar rule at X
        call plane_geometry(1, n, derham)
        call spline_quadrature_build(space, degree + 8, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_quadrature_point_singular(space, quad, a(:, 1), aquad, st)
        if (st /= 0) error stop 'axis quadrature'
        ntot = n0*nt
        nf = 2*ntot
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, wall1)
        wall = [((k - 1)*n0 + wall1, k=1, nt)]

        ! fixed-label fit (E4a) on this mesh: the interpolation level of B
        allocate (g(ntot), delta(nf))
        g = 0
        call clebsch_3d_system(derham, aquad, nt, section_profiles, a(:, 1), chi0, g, energy, grad, &
                               hessian, a, base, frames)
        call add_g_ridge(hessian, ntot)
        call eu_constrained_solve(hessian, -grad, [(i, i=1, ntot)], 0*chi0, delta, st)
        if (st /= 0) error stop 'fixed-label g solve'
        g = delta(ntot + 1:)
        call clebsch_3d_system(derham, aquad, nt, section_profiles, a(:, 1), chi0, g, energy, grad, &
                               plane_axis=a, plane_nodes=base, plane_frames=frames, &
                               exact_field=section_exact_b, field_error=ferr)
        err_fixed = sqrt(ferr(1)/ferr(2))

        ! start: labels perturbed per plane (zero at the axis and the wall),
        ! a helical axis offset, g fitted at fixed labels
        allocate (chi(ntot))
        do k = 1, nt
            tk = (k - 1)/real(nt, dp)
            associate (c0 => chi0((k - 1)*n0 + 1:k*n0))
                chi((k - 1)*n0 + 1:k*n0) = c0 + 4*amplitude*(1 + 0.5_dp*cos(2*pi*tk))*c0*(1 - c0)
            end associate
            x(:, k) = a(:, k) + axis_offset*[cos(2*pi*tk), sin(2*pi*tk)]
        end do
        chi(wall) = 1
        g = 0
        call clebsch_3d_system(derham, aquad, nt, section_profiles, a(:, 1), chi, g, energy, grad, &
                               hessian, x, nodes_at(base, w, a, x), frames)
        call add_g_ridge(hessian, ntot)
        call eu_constrained_solve(hessian, -grad, [(i, i=1, ntot)], 0*chi, delta, st)
        if (st /= 0) error stop 'g solve'
        g = delta(ntot + 1:)
        deallocate (delta)
        mu = 0
        allocate (history(0), dgdx(nf + 3*nt, 2*nt))
        steps = 0
        do it = 1, max_steps
            call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, chi, g, x, mu, &
                          r, hessian)
            history = [history, norm2(r)]
            ! the X columns are central differences: 1e-8 is their floor
            if (it > 1) then
                if (history(it) < 1e-8_dp*history(1)) exit
            end if
            do j = 1, 2*nt
                xp = x
                xp(mod(j - 1, 2) + 1, (j + 1)/2) = x(mod(j - 1, 2) + 1, (j + 1)/2) + hx
                call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, chi, g, xp, &
                              mu, rp)
                xp(mod(j - 1, 2) + 1, (j + 1)/2) = x(mod(j - 1, 2) + 1, (j + 1)/2) - hx
                call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, chi, g, xp, &
                              mu, rm)
                dgdx(:, j) = (rp - rm)/(2*hx)
            end do
            call add_g_ridge(hessian, ntot, newton_ridge)
            call ale_system(hessian, n0, nrow, ngrad, dgdx, jac_t)
            if (jacobian_check .and. it == 1) then
                jacobian: block
                    !! J v against central differences of the residual, per
                    !! block of v (chi, g, X, mu), wall labels held fixed and
                    !! wall rows excluded.
                    real(dp), allocatable :: v(:), jv(:), fd(:), rpv(:), rmv(:)
                    real(dp) :: hv
                    integer :: blk, lo, hi, jj, p

                    allocate (v(nf + 3*nt), jv(nf + 3*nt))
                    do blk = 1, 4
                        select case (blk)
                        case (1)
                            lo = 1; hi = ntot
                        case (2)
                            lo = ntot + 1; hi = nf
                        case (3)
                            lo = nf + 1; hi = nf + 2*nt
                        case default
                            lo = nf + 2*nt + 1; hi = nf + 3*nt
                        end select
                        v = 0
                        do jj = lo, hi
                            v(jj) = sin(1.3_dp*jj) + 0.5_dp
                        end do
                        v(wall) = 0
                        hv = merge(1e-6_dp, 1e-5_dp, blk == 3)/maxval(abs(v))
                        jv = 0
                        do jj = 1, jac_t%nrow
                            do p = jac_t%row_ptr(jj), jac_t%row_ptr(jj + 1) - 1
                                jv(jac_t%col(p)) = jv(jac_t%col(p)) + jac_t%val(p)*v(jj)
                            end do
                        end do
                        call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, &
                                      chi + hv*v(1:ntot), g + hv*v(ntot + 1:nf), &
                                      x + hv*reshape(v(nf + 1:nf + 2*nt), [2, nt]), &
                                      mu + hv*v(nf + 2*nt + 1:), rpv)
                        call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, &
                                      chi - hv*v(1:ntot), g - hv*v(ntot + 1:nf), &
                                      x - hv*reshape(v(nf + 1:nf + 2*nt), [2, nt]), &
                                      mu - hv*v(nf + 2*nt + 1:), rmv)
                        fd = (rpv - rmv)/(2*hv)
                        ! wall rows: the residual is zero there and the solve fixes
                        ! the wall labels, so the matrix rows there are not compared
                        jv(wall) = 0
                        print '(a,a6,a,es10.2,a,es10.2,a,3es10.2)', '  Jacobian check, v in ', &
                            trim(merge('chi', merge('g  ', merge('X  ', 'mu ', blk == 3), blk == 2), &
                                       blk == 1)), ': |Jv - FD|/|FD|', norm2(jv - fd)/norm2(fd), &
                            '  |FD|', norm2(fd), '  rows chi, g, rest', &
                            norm2(jv(1:ntot) - fd(1:ntot))/max(norm2(fd(1:ntot)), tiny(1.0_dp)), &
                            norm2(jv(ntot + 1:nf) - fd(ntot + 1:nf))/max(norm2(fd(ntot + 1:nf)), &
                                                                          tiny(1.0_dp)), &
                            norm2(jv(nf + 1:) - fd(nf + 1:))/max(norm2(fd(nf + 1:)), tiny(1.0_dp))
                    end do
                end block jacobian
            end if
            allocate (delta(nf + 3*nt))
            call eu_constrained_solve(jac_t, -r, wall, 0*chi(wall), delta, st)
            if (st /= 0) error stop 'Newton solve'
            ! non-monotone, as in E3b
            alpha = 1
            do k = 1, 8
                call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, &
                              chi + alpha*delta(1:ntot), g + alpha*delta(ntot + 1:nf), &
                              x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt]), &
                              mu + alpha*delta(nf + 2*nt + 1:), rt)
                if (norm2(rt) < 10*maxval(history)) exit
                alpha = alpha/2
            end do
            print '(a,i3,a,4es9.2,a,f6.3,a,3es9.2)', '  step', it, '  |r| chi, g, axis, mu', &
                norm2(r(1:ntot)), norm2(r(ntot + 1:nf)), norm2(r(nf + 1:nf + 2*nt)), &
                norm2(r(nf + 2*nt + 1:)), '  alpha', alpha, '  |d chi|, |d g|, |d X|', &
                maxval(abs(delta(1:ntot))), maxval(abs(delta(ntot + 1:nf))), &
                maxval(abs(delta(nf + 1:nf + 2*nt)))
            chi = chi + alpha*delta(1:ntot)
            g = g + alpha*delta(ntot + 1:nf)
            x = x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt])
            mu = mu + alpha*delta(nf + 2*nt + 1:)
            deallocate (delta)
            steps = steps + 1
        end do
        print '(a,16es9.1)', '  |r| history', history
        call residual(derham, aquad, a, base, w, frames, nrow, ngrad, wall, chi, g, x, mu, r, &
                      ferr=ferr)
        err = sqrt(ferr(1)/ferr(2))
        axis_err = maxval(norm2(x - a, 1))
        ! labels against the exact label projected on the final (moved) mesh
        chi_err = 0
        do k = 1, nt
            call plane_geometry(k, n, derham, x(:, k))
            chi_exact = projected(derham, 3)
            chi_err = max(chi_err, maxval(abs(chi((k - 1)*n0 + 1:k*n0) - chi_exact)))
        end do
        shape = maxval(abs(x(1, :) - a(1, :)))
        print '(a,es10.3,a,es10.3,a,es10.3)', '  max |X - X_exact|', axis_err, '  (R part', &
            shape, ')  multiplier spread', maxval(mu) - minval(mu)
    end subroutine solve

    function section_centre(k) result(c)
        !! The exact axis of plane k (the centre of its O-grid).
        integer, intent(in) :: k
        real(dp) :: c(2)
        type(block_layout_t) :: layout

        call section_layout(2*pi*(k - 1)/nt, c, layout)
    end function section_centre

    subroutine ale_system(hessian, n0, nrow, ngrad, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, g, X_k, mu_k): the (chi, g)
        !! block is the Hessian, the X columns are dgdx (all rows, with the
        !! mesh motion), the axis rows of plane k use its grad N(X_k), the
        !! mu_k column is N(X_k) on plane k.
        type(csr_t), intent(in) :: hessian
        integer, intent(in) :: n0
        real(dp), intent(in) :: nrow(:), ngrad(:, :, :), dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        type(triplet_list) :: t
        integer :: nf, nx, i, j, k, c

        nf = 2*n0*nt
        nx = 2*nt
        call t%init(size(hessian%val) + nx*size(dgdx, 1) + 6*n0*nt + 4*nt)
        call t%add_csr(hessian, 0, 0)
        do j = 1, nx
            do i = 1, nf + 3*nt
                if (dgdx(i, j) /= 0) call t%add(nf + j, i, dgdx(i, j))
            end do
            call t%add(nf + j, nf + j, 0.0_dp)
        end do
        do k = 1, nt
            do i = 1, n0
                do c = 1, 2
                    call t%add((k - 1)*n0 + i, nf + 2*(k - 1) + c, ngrad(i, c, k))
                end do
                call t%add(nf + nx + k, (k - 1)*n0 + i, nrow(i))
                call t%add((k - 1)*n0 + i, nf + nx + k, nrow(i))
            end do
            call t%add(nf + nx + k, nf + nx + k, 0.0_dp)
        end do
        call t%to_csr(nf + 3*nt, nf + 3*nt, jac_t)
    end subroutine ale_system

    subroutine add_g_ridge(hessian, n0, relative)
        !! g -> g + k(chi) leaves B unchanged: relative ridge on the g block.
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: n0
        real(dp), intent(in), optional :: relative
        real(dp) :: shift_
        integer :: i, p

        shift_ = ridge*maxval(abs(hessian%val))
        if (present(relative)) shift_ = relative*maxval(abs(hessian%val))
        do i = n0 + 1, 2*n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift_
            end do
        end do
    end subroutine add_g_ridge

end program test_eulerian_clebsch_landreman_newton
