module milca_gs_free_current
    !! Free-boundary Grad-Shafranov equilibrium with the profiles of a VMEC
    !! run: pressure p(s) and enclosed toroidal current I(s) in the normalised
    !! toroidal flux s = Phi/Phi_edge, and the toroidal flux Phi_edge inside
    !! the plasma boundary (docs/gs-free-boundary.md, "Current profiles").
    !!
    !! The solver's law is p'(psi) and FF'(psi) in the normalised poloidal
    !! flux. With u = sigma (psi_a - psi) (sigma = axis_sign), the contour
    !! integrals a1 = oint R dl/|grad psi|, a2 = oint dl/(R |grad psi|) of the
    !! current solution give
    !!
    !!     dPhi/du = F a2,                                             (1)
    !!     dI/du   = oint j_phi dl/|grad u| = a1 p_psi + a2 FF'_psi/mu0,  (2)
    !!
    !! with p_psi = dp/dpsi = -sigma (dp/ds)(ds/du), ds/du = F a2/Phi_edge
    !! and dI/du = sigma (dI/ds)(ds/du) (the current flows along sigma). The
    !! boundary level u_b solves Phi(u_b) = Phi_edge with (1); the cut
    !! fraction of the X-point mode becomes u_b/u_X. Equation (2) solved for
    !! FF'_psi on the law knots, and p_psi, give the next law (cubic Hermite
    !! pieces, C1): gs_current_update. gs_fb_current_solve alternates Newton
    !! solves and updates (a Picard iteration, as the current modes of
    !! CHEASE and HELENA).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_polynomial, only: lagrange_weights
    use milca_grad_shafranov, only: gs_state_t
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_primitives, gs_law_sources
    use milca_gs_free_boundary, only: gs_fb_system_t, gs_fb_eval_t, gs_fb_newton, &
        gs_fb_history_t, &
        gs_fb_evaluate, gs_fb_surface_integrals, GS_FB_OK, GS_FB_XPOINT
    implicit none
    private
    public :: gs_current_target_t, gs_current_report_t, gs_read_current_target, &
        gs_current_update, gs_fb_current_solve, gs_current_profiles, &
        gs_write_current_target
    integer, parameter, public :: GS_CUR_OK = 0, GS_CUR_BAD_INPUT = 1, &
        GS_CUR_SOLVE_FAILED = 2, GS_CUR_CONTOUR_FAILED = 3, &
        GS_CUR_NOT_CONVERGED = 4, GS_CUR_FLUX_UNREACHABLE = 5
    integer, parameter :: table = 48

    type :: gs_current_target_t
        !> s, p (Pa), dp/ds, I (A, magnitude), dI/ds on increasing s in [0, 1]
        real(dp), allocatable :: s(:), p(:), dp(:), current(:), dcurrent(:)
        real(dp) :: phi_edge = 0 !< Wb, magnitude
    end type

    type :: gs_current_report_t
        integer :: passes = 0
        real(dp) :: change = 0 !< last relative law change
        real(dp) :: cut_fraction = 1
        real(dp) :: phi_deepest = 0 !< Phi at the deepest traced level
        real(dp) :: phi_edge = 0 !< Phi inside the boundary of the last update
        real(dp), allocatable :: history(:)
    end type

contains

    subroutine gs_read_current_target(path, target, status)
        !! Text file: "n phi_edge", then n rows "s p dp/ds I dI/ds".
        character(*), intent(in) :: path
        type(gs_current_target_t), intent(out) :: target
        integer, intent(out) :: status
        integer :: u, n, k, ios
        real(dp) :: row(5)
        status = GS_CUR_BAD_INPUT
        open (newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (u, *, iostat=ios) n, target%phi_edge
        if (ios /= 0 .or. n < 4) then
            close (u)
            return
        end if
        allocate (target%s(n), target%p(n), target%dp(n), target%current(n), &
            target%dcurrent(n))
        do k = 1, n
            read (u, *, iostat=ios) row
            if (ios /= 0) then
                close (u)
                return
            end if
            target%s(k) = row(1); target%p(k) = row(2); target%dp(k) = row(3)
            target%current(k) = row(4); target%dcurrent(k) = row(5)
        end do
        close (u)
        if (abs(target%s(1)) > 1.0e-12_dp .or. abs(target%s(n) - 1) > 1.0e-12_dp) return
        ! Nodes closer than 1e-9 make the cubic interpolation weights blow up
        if (any(target%s(2:) - target%s(:n - 1) < 1.0e-9_dp)) return
        if (.not. target%phi_edge > 0) return
        status = GS_CUR_OK
    end subroutine

    subroutine gs_current_update(state, psi, axis, psi_axis, psi_bound, cut, sigma, &
                                 mu0, law, target, contour_points, deepest, &
                                 next, next_cut, phi_deepest, status, fixed_cut, &
                                 phi_edge)
        !! One profile update from the flux psi with axis value psi_axis and
        !! boundary psi_bound = psi_axis + cut (psi_X - psi_axis), solved under
        !! law: next law on the knots of law (p_edge from the target, F_edge
        !! kept) and the cut fraction that encloses target%phi_edge.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), axis(2), psi_axis, psi_bound, cut, sigma, mu0
        type(gs_profile_law_t), intent(in) :: law
        type(gs_current_target_t), intent(in) :: target
        integer, intent(in) :: contour_points
        real(dp), intent(in) :: deepest
        type(gs_profile_law_t), intent(out) :: next
        real(dp), intent(out) :: next_cut, phi_deepest
        integer, intent(out) :: status
        !> fixed_cut: keep the boundary level; Phi_edge is then the flux it
        !> encloses (returned in phi_edge) and s is normalised by it.
        logical, intent(in), optional :: fixed_cut
        real(dp), intent(out), optional :: phi_edge
        type(gs_current_target_t) :: goal
        logical :: keep
        type(gs_profile_law_t) :: current
        real(dp) :: t(table), u(table), a1(table), a2(table), f(table), phi(table)
        real(dp), allocatable :: uk(:), k1(:), k2(:), fk(:), sk(:), pp(:), ff(:)
        real(dp) :: ux, ub, ub_old, s_tor, dsdu, dpds, dids
        integer :: j, nk, stat, attempt
        real(dp) :: top

        status = GS_CUR_BAD_INPUT
        next_cut = cut
        phi_deepest = 0
        if (.not. (cut > 0 .and. cut <= 1) .or. .not. (deepest > 0 .and. deepest < 1)) return
        nk = size(law%knots)
        allocate (uk(nk), k1(nk), k2(nk), fk(nk), sk(nk), pp(nk), ff(nk))
        current = law
        current%axis = psi_axis - psi_bound
        ub_old = sigma*(psi_axis - psi_bound)
        ux = ub_old/cut
        ! (1) on levels t = u/u_X graded towards t_top, a quarter of the way from the
        ! current boundary to the separatrix (surfaces near the separatrix
        ! can leave the mesh), or to deepest when that does not reach
        ! Phi_edge.
        top = min(deepest, cut + 0.25_dp*(1 - cut))
        if (present(fixed_cut)) then
            if (fixed_cut) top = min(deepest, cut + 0.1_dp*(1 - cut))
        end if
        do attempt = 1, 2
            do j = 1, table
                t(j) = top*(1 - (1 - real(j, dp)/table)**2)
            end do
            u = t*ux
            status = GS_CUR_CONTOUR_FAILED
            call gs_fb_surface_integrals(state, psi, axis, psi_axis - sigma*u, sigma, &
                                         contour_points, a1, a2, stat)
            if (stat /= GS_FB_OK) return
            do j = 1, table
                f(j) = flux_function(current, u(j), ub_old)
            end do
            ! Phi(u) = int_0^u F a2: on each interval the local cubic
            ! interpolant of F a2 at two Gauss points (O(h^4)); [0, u(1)]
            ! extrapolates it.
            phi = 0
            do j = 1, table
                phi(j) = merge(0.0_dp, phi(max(j - 1, 1)), j == 1) + &
                    cubic_integral(u, f*a2, merge(0.0_dp, u(max(j - 1, 1)), j == 1), u(j))
            end do
            phi_deepest = phi(table)
            if (phi(table) >= target%phi_edge .or. top >= deepest) exit
            if (present(fixed_cut)) then
                if (fixed_cut) exit
            end if
            top = deepest
        end do
        keep = .false.
        if (present(fixed_cut)) keep = fixed_cut
        goal = target
        if (keep) then
            ub = ub_old
            goal%phi_edge = interpolate4(u, phi, ub, count(u < ub))
        else
            if (phi(table) < target%phi_edge) then
                status = GS_CUR_FLUX_UNREACHABLE
                return
            end if
            ub = interpolate4(phi, u, target%phi_edge, count(phi < target%phi_edge))
        end if
        if (present(phi_edge)) phi_edge = goal%phi_edge
        ! (2) on the knots of the next law; the axis knot at a small depth
        uk = law%knots*ub
        uk(1) = 1.0e-2_dp*uk(2)
        call gs_fb_surface_integrals(state, psi, axis, psi_axis - sigma*uk, sigma, &
                                     contour_points, k1, k2, stat)
        if (stat /= GS_FB_OK) return
        do j = 1, nk
            fk(j) = flux_function(current, uk(j), ub_old)
            sk(j) = interpolate4(u, phi, uk(j), count(u < uk(j)))/goal%phi_edge
        end do
        sk(nk) = 1
        do j = 1, nk
            s_tor = min(max(sk(j), 0.0_dp), 1.0_dp)
            dsdu = fk(j)*k2(j)/goal%phi_edge
            dpds = interpolate4(goal%s, goal%dp, s_tor, count(goal%s < s_tor))
            dids = interpolate4(goal%s, goal%dcurrent, s_tor, count(goal%s < s_tor))
            pp(j) = -sigma*dpds*dsdu
            ff(j) = mu0*(sigma*dids*dsdu - k1(j)*pp(j))/k2(j)
        end do
        ! The axis knot was evaluated at the small depth uk(1); linear
        ! extrapolation to u = 0 (error O(uk(1) uk(2)))
        pp(1) = pp(1) - (pp(2) - pp(1))*uk(1)/(uk(2) - uk(1))
        ff(1) = ff(1) - (ff(2) - ff(1))*uk(1)/(uk(2) - uk(1))
        next = law
        call hermite_law(law%knots, pp, next%pp)
        call hermite_law(law%knots, ff, next%ff)
        next%p_edge = goal%p(size(goal%p))
        next_cut = ub/ux
        status = GS_CUR_OK
    end subroutine gs_current_update

    subroutine gs_current_profiles(state, psi, axis, psi_axis, psi_bound, sigma, mu0, &
                                   law, levels, contour_points, target, status)
        !! p(s), I(s) and Phi_edge of a solution under law with boundary
        !! psi_bound, on levels + 1 points uniform in the poloidal depth u
        !! (the inverse of gs_current_update): Phi and I by the fourth-order
        !! rule of (1) and (2) on the contour integrals.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), axis(2), psi_axis, psi_bound, sigma, mu0
        type(gs_profile_law_t), intent(in) :: law
        integer, intent(in) :: levels, contour_points
        type(gs_current_target_t), intent(out) :: target
        integer, intent(out) :: status
        type(gs_profile_law_t) :: current
        real(dp), allocatable :: u(:), a1(:), a2(:), f(:), dpdu(:), didu(:), phi(:), cur(:)
        real(dp) :: ub, pressure, pp, ff
        integer :: j, st

        status = GS_CUR_BAD_INPUT
        if (levels < 4) return
        current = law
        current%axis = psi_axis - psi_bound
        ub = sigma*(psi_axis - psi_bound)
        allocate (u(levels + 1), a1(levels + 1), a2(levels + 1), f(levels + 1), &
            dpdu(levels + 1), didu(levels + 1), phi(levels + 1), cur(levels + 1))
        u = [(ub*real(j, dp)/levels, j=0, levels)]
        u(1) = 1.0e-2_dp*u(2)
        status = GS_CUR_CONTOUR_FAILED
        call gs_fb_surface_integrals(state, psi, axis, psi_axis - sigma*u, sigma, &
                                     contour_points, a1, a2, st)
        if (st /= GS_FB_OK) return
        u(1) = 0
        a1(1) = a1(1) - (a1(2) - a1(1))*1.0e-2_dp/0.99_dp
        a2(1) = a2(1) - (a2(2) - a2(1))*1.0e-2_dp/0.99_dp
        allocate (target%s(levels + 1), target%p(levels + 1), target%dp(levels + 1), &
            target%current(levels + 1), target%dcurrent(levels + 1))
        do j = 1, levels + 1
            call gs_law_primitives(current, current%axis*(1 - u(j)/ub), pressure, f(j), st)
            f(j) = abs(f(j))
            target%p(j) = pressure
            call gs_law_sources(current, current%axis*(1 - u(j)/ub), pp, ff)
            dpdu(j) = -sigma*pp ! dp/du
            didu(j) = sigma*(a1(j)*pp + a2(j)*ff/mu0) ! dI/du
        end do
        phi(1) = 0; cur(1) = 0
        do j = 2, levels + 1
            phi(j) = phi(j - 1) + cubic_integral(u, f*a2, u(j - 1), u(j))
            cur(j) = cur(j - 1) + cubic_integral(u, didu, u(j - 1), u(j))
        end do
        target%phi_edge = phi(levels + 1)
        target%s = phi/phi(levels + 1)
        target%current = cur
        target%dp = dpdu/(f*a2/target%phi_edge)
        target%dcurrent = didu/(f*a2/target%phi_edge)
        status = GS_CUR_OK
    end subroutine gs_current_profiles

    subroutine gs_write_current_target(path, target)
        character(*), intent(in) :: path
        type(gs_current_target_t), intent(in) :: target
        integer :: u, k
        open (newunit=u, file=path, status='replace', action='write')
        write (u, '(i0,1x,es24.16)') size(target%s), target%phi_edge
        do k = 1, size(target%s)
            write (u, '(5es24.16)') target%s(k), target%p(k), target%dp(k), &
                target%current(k), target%dcurrent(k)
        end do
        close (u)
    end subroutine

    subroutine gs_fb_current_solve(state, sys, currents, target, status, report, &
                                   max_passes, tolerance, contour_points, &
                                   deepest, relaxation, verbose, fixed_cut)
        !! Outer iteration from a converged free-boundary solution in state%psi
        !! with sys set up in the X-point mode. Converged when the law values
        !! at the knots and the cut fraction change by less than tolerance
        !! (relative); the final Newton solve uses the converged law.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        real(dp), intent(in) :: currents(:)
        type(gs_current_target_t), intent(in) :: target
        integer, intent(out) :: status
        type(gs_current_report_t), intent(out) :: report
        integer, intent(in), optional :: max_passes, contour_points
        real(dp), intent(in), optional :: tolerance, deepest, relaxation
        logical, intent(in), optional :: verbose, fixed_cut
        type(gs_fb_eval_t) :: ev
        type(gs_profile_law_t) :: next
        real(dp) :: sigma, c_max, tol, omega, change, cut, next_cut, phi_deepest
        real(dp), allocatable :: x(:), g(:), f(:), f_old(:), g_old(:), weight(:), &
            df(:, :), dg(:, :), gamma(:)
        integer, parameter :: memory = 6
        integer :: npass, npts, pass, stat, nk, used, retry
        type(gs_fb_history_t) :: history
        real(dp), allocatable :: psi_accepted(:), x_accepted(:)
        logical :: talk

        status = GS_CUR_BAD_INPUT
        if (.not. sys%ready .or. sys%cfg%mode /= GS_FB_XPOINT) return
        if (.not. allocated(target%s)) return
        npass = 40; if (present(max_passes)) npass = max_passes
        tol = 1.0e-7_dp; if (present(tolerance)) tol = tolerance
        npts = 256; if (present(contour_points)) npts = contour_points
        c_max = 0.995_dp; if (present(deepest)) c_max = deepest
        omega = 1.0_dp; if (present(relaxation)) omega = relaxation
        talk = .false.; if (present(verbose)) talk = verbose
        sigma = real(sys%cfg%axis_sign, dp)
        allocate (report%history(npass))
        nk = size(sys%law%knots)
        allocate (df(2*nk + 1, memory), dg(2*nk + 1, memory), gamma(memory))
        used = 0
        do pass = 1, npass
            ! A failed solve halves the last step back towards the accepted
            ! state (at most six times).
            do retry = 0, 6
                call gs_fb_newton(state, sys, currents, stat, history)
                if (stat == GS_FB_OK .or. pass == 1) exit
                if (retry == 6) exit
                state%psi = psi_accepted
                x = 0.5_dp*(pack_state(sys%law, sys%cfg%cut_fraction) + x_accepted)
                call unpack_state(x, sys%law, sys%cfg%cut_fraction)
                if (talk) print '(a,i0)', '  solve failed, step halved ', retry + 1
            end do
            if (stat /= GS_FB_OK) then
                status = GS_CUR_SOLVE_FAILED
                return
            end if
            psi_accepted = state%psi
            x_accepted = pack_state(sys%law, sys%cfg%cut_fraction)
            call gs_fb_evaluate(state, sys, state%psi, currents, ev, stat, &
                jacobian=.false.)
            if (stat /= GS_FB_OK) then
                status = GS_CUR_SOLVE_FAILED
                return
            end if
            cut = sys%cfg%cut_fraction
            call gs_current_update(state, state%psi, ev%axis, ev%psi_axis, &
                ev%psi_bound, cut, sigma, sys%mu0, sys%law, target, npts, c_max, &
                next, next_cut, phi_deepest, stat, fixed_cut, report%phi_edge)
            report%phi_deepest = phi_deepest
            if (stat /= GS_CUR_OK) then
                status = stat
                return
            end if
            ! Anderson acceleration of the fixed point x = G(x) on the scaled
            ! knot values and cut fraction (Walker and Ni, SIAM J. Numer.
            ! Anal. 49 (2011) 1715): the Picard map is unstable at the fixed
            ! point (deviations grow 2-7 fold per pass on DIII-D).
            x = pack_state(sys%law, cut)
            g = pack_state(next, next_cut)
            if (pass == 1) then
                weight = [spread(max(maxval(abs(x(1:nk))), tiny(1.0_dp)), 1, nk), &
                    spread(max(maxval(abs(x(nk + 1:2*nk))), tiny(1.0_dp)), 1, nk), &
                    1.0e-2_dp]
            end if
            x = x/weight
            g = g/weight
            f = g - x
            change = maxval(abs(f))
            if (pass > 1) then
                call push_history(df, f - f_old, dg, g - g_old, used, memory)
            end if
            f_old = f
            g_old = g
            if (used > 0) then
                call least_squares(df(:, :used), f, gamma)
                g = g - matmul(dg(:, :used), gamma(:used))
            end if
            g = (1 - omega)*x + omega*g
            g = g*weight
            call unpack_state(g, next, next_cut)
            report%history(pass) = change
            report%passes = pass
            report%change = change
            report%cut_fraction = next_cut
            if (talk) print '(a,i3,a,es10.3,a,f9.6,a,es13.6,a,f9.5,a,i0)', 'current pass', &
                pass, ' change', change, ' cut', next_cut, ' Ip', ev%plasma_current, &
                ' Phi(deepest)', phi_deepest, ' Newton iterations ', history%iterations
            sys%law%pp = next%pp
            sys%law%ff = next%ff
            sys%law%p_edge = next%p_edge
            sys%cfg%cut_fraction = next_cut
            if (change <= tol) then
                call gs_fb_newton(state, sys, currents, stat)
                status = merge(GS_CUR_OK, GS_CUR_SOLVE_FAILED, stat == GS_FB_OK)
                return
            end if
        end do
        status = GS_CUR_NOT_CONVERGED
    contains
        function pack_state(law, fraction) result(v)
            type(gs_profile_law_t), intent(in) :: law
            real(dp), intent(in) :: fraction
            real(dp) :: v(2*nk + 1)
            v(1:nk - 1) = law%pp(1, :)
            v(nk) = hermite_end(law%pp, law%knots)
            v(nk + 1:2*nk - 1) = law%ff(1, :)
            v(2*nk) = hermite_end(law%ff, law%knots)
            v(2*nk + 1) = fraction
        end function

        subroutine unpack_state(v, law, fraction)
            real(dp), intent(in) :: v(:)
            type(gs_profile_law_t), intent(inout) :: law
            real(dp), intent(out) :: fraction
            call hermite_law(law%knots, v(1:nk), law%pp)
            call hermite_law(law%knots, v(nk + 1:2*nk), law%ff)
            fraction = min(v(2*nk + 1), c_max)
        end subroutine
    end subroutine gs_fb_current_solve

    pure function hermite_end(c, knots) result(value)
        !! Value of the last piece at the last knot.
        real(dp), intent(in) :: c(:, :), knots(:)
        real(dp) :: value, h
        integer :: n
        n = size(knots)
        h = knots(n) - knots(n - 1)
        value = c(1, n - 1) + h*(c(2, n - 1) + h*(c(3, n - 1) + h*c(4, n - 1)))
    end function

    subroutine push_history(df, dfn, dg, dgn, used, memory)
        !! Append a column, dropping the oldest when memory columns are held.
        real(dp), intent(inout) :: df(:, :), dg(:, :)
        real(dp), intent(in) :: dfn(:), dgn(:)
        integer, intent(inout) :: used
        integer, intent(in) :: memory
        if (used == memory) then
            df(:, 1:memory - 1) = df(:, 2:memory)
            dg(:, 1:memory - 1) = dg(:, 2:memory)
            used = memory - 1
        end if
        used = used + 1
        df(:, used) = dfn
        dg(:, used) = dgn
    end subroutine

    subroutine least_squares(a, b, x)
        !! min |a x - b| by modified Gram-Schmidt QR (a has few columns);
        !! columns that are numerically dependent get x = 0.
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: q(size(a, 1), size(a, 2)), r(size(a, 2), size(a, 2)), y(size(a, 2))
        integer :: i, j, m
        m = size(a, 2)
        q = a
        r = 0
        do j = 1, m
            do i = 1, j - 1
                r(i, j) = dot_product(q(:, i), q(:, j))
                q(:, j) = q(:, j) - r(i, j)*q(:, i)
            end do
            r(j, j) = norm2(q(:, j))
            if (r(j, j) > 1.0e-12_dp*max(norm2(a(:, j)), tiny(1.0_dp))) then
                q(:, j) = q(:, j)/r(j, j)
            else
                q(:, j) = 0
                r(j, j) = 0
            end if
        end do
        y = matmul(transpose(q), b)
        x = 0
        do j = m, 1, -1
            if (r(j, j) == 0) cycle
            x(j) = (y(j) - dot_product(r(j, j + 1:m), x(j + 1:m)))/r(j, j)
        end do
    end subroutine

    real(dp) function flux_function(law, level, bound) result(value)
        !! |F| at depth level below the axis: the law inside the boundary at
        !! depth bound, the vacuum F_edge outside. law%axis = psi_a - psi_b.
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: level, bound
        real(dp) :: pressure
        integer :: st
        value = abs(law%f_edge)
        if (level >= bound) return
        call gs_law_primitives(law, law%axis*(1 - level/bound), pressure, value, st)
        value = abs(value)
    end function

    pure subroutine hermite_law(knots, y, c)
        !! C1 cubic Hermite pieces through y at the knots, slopes by the
        !! three-point (non-uniform) difference, one-sided at the ends.
        real(dp), intent(in) :: knots(:), y(:)
        real(dp), intent(inout) :: c(:, :)
        real(dp) :: m(size(y)), h, d, h0, h1
        integer :: i, n
        n = size(knots)
        do i = 2, n - 1
            h0 = knots(i) - knots(i - 1); h1 = knots(i + 1) - knots(i)
            m(i) = (h0*h0*(y(i + 1) - y(i)) + h1*h1*(y(i) - y(i - 1)))/(h0*h1*(h0 + h1))
        end do
        m(1) = 2*(y(2) - y(1))/(knots(2) - knots(1)) - m(2)
        m(n) = 2*(y(n) - y(n - 1))/(knots(n) - knots(n - 1)) - m(n - 1)
        do i = 1, n - 1
            h = knots(i + 1) - knots(i)
            d = (y(i + 1) - y(i))/h
            c(1, i) = y(i)
            c(2, i) = m(i)
            c(3, i) = (3*d - 2*m(i) - m(i + 1))/h
            c(4, i) = (m(i) + m(i + 1) - 2*d)/(h*h)
        end do
    end subroutine

    real(dp) function cubic_integral(x, y, lo, hi) result(value)
        !! Two-point Gauss rule on [lo, hi] of the cubic Lagrange interpolant
        !! of y(x) through the four nodes around the interval: exact for it.
        real(dp), intent(in) :: x(:), y(:), lo, hi
        real(dp) :: g(2), mid
        integer :: k, below
        mid = 0.5_dp*(lo + hi)
        below = count(x < mid)
        g = mid + 0.5_dp*(hi - lo)*[-1, 1]/sqrt(3.0_dp)
        value = 0
        do k = 1, 2
            value = value + 0.5_dp*(hi - lo)*interpolate4(x, y, g(k), below)
        end do
    end function

    real(dp) function interpolate4(x, y, at, below) result(value)
        !! Cubic Lagrange interpolation of y(x) at 'at' on the four nodes
        !! around the interval (below, below + 1).
        real(dp), intent(in) :: x(:), y(:), at
        integer, intent(in) :: below
        real(dp) :: w(4)
        integer :: first
        first = max(1, min(below - 1, size(x) - 3))
        call lagrange_weights(4, at, x(first:first + 3), w)
        value = dot_product(w, y(first:first + 3))
    end function
end module milca_gs_free_current
