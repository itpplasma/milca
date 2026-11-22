program test_label_ale_doublet_clebsch
    !! Gate E2 (separatrix part) of docs/eulerian-label-plan.md: Newton for
    !! the exact Solov'ev doublet in the direct-Clebsch form
    !!
    !!   A = g grad chi + sum_j F_j(chi) h_j + G(chi) grad phi,
    !!   W = int [ |B|^2/2 - p(chi) ] dV,
    !!
    !! on the separatrix-aligned, graded E2b mesh (FortFEM layout_doublet),
    !! chi in continuous V0 (chi = 1 on the mesh separatrix, 2 on the wall),
    !! g in V0 cut along the separatrix, winding forms h_j about the O-points
    !! in the frame H^(1/2) of the label (lobe j: h_j with F_j' = -kappa
    !! Phi_lobe'/(2 pi); outside: both, each with half the flux). The O- and
    !! X-points stay at the exact points (the mesh motion is gate A2a):
    !! chi_h(O_j) = 0 and grad chi_h(O_j) = 0 by multipliers. The profiles
    !! are least-squares fits of the traced flux slope in t^m, t^m log|t|
    !! (t = chi - 1, m <= 4, 5e-5 relative), so the Newton matrix is exact.
    !!
    !! kappa scales the toroidal field (f = kappa f0; any constant f is an
    !! exact Solov'ev state). The doublet has B_phi/B_pol ~ 1e3 kappa. Finding
    !! (recorded in the plan): Newton converges quadratically at kappa = 1e-3,
    !! slowly at 1e-2 and not at 0.1 and 1, with the log-singular profiles and
    !! equally with their log terms removed; the Jacobian matches differences
    !! of the residual to 1e-9 at kappa = 1. The obstruction is the guide
    !! field, not the separatrix: after elimination of g the label Hessian
    !! has O(B_phi^2) negative grid-scale modes in the lobe cores, which
    !! swamp the poloidal stiffness O(B_pol^2). The single-axis gate E1
    !! shows the same trend (f0 = 30: labels converge at order 0.6; f0 = 300:
    !! Newton diverges). The multiplicative label update (Newton in
    !! log|chi_i - 1|) and g of degree 3 do not help.
    !!
    !! Passed (kappa = 1e-3): Newton converges in at most 12 steps on 4 and 8
    !! cells, B_phi and B_pol within 20 % of the fixed-label fit on the same
    !! mesh (E2b level), the rate of B_phi within 0.3 of the fit's rate, and
    !! the Newton matrix agrees with central differences of the residual.
    !!
    !! Arguments: levels (2), update ('add' or 'mult', add), amplitude (0.1),
    !! grading power (4), layer fraction (0.5), 'jacobian' or 'dump' or
    !! 'none', kappa (1e-3), 'smooth' (log terms removed) or 'log', degree
    !! of g (2).
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_SEPARATRIX, LAYOUT_WALL
    use fortfem_spline_separatrix, only: doublet_curves_t, layout_doublet
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_grade_knots, knot_spec_t, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate, &
        spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_spline_ale, only: spline_geometry_dofs
    use solovev_label_fields, only: triplet_list
    use solovev_doublet, only: solovev_doublet_t, doublet_f0
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: degree = 2, nbasis = 10
    integer :: cells(4) = [4, 8, 16, 32], levels = 2, level, failures, max_steps = 30, &
               g_degree = degree
    real(dp) :: power = 4, theta = 0.5_dp, amplitude = 0.1_dp, ridge = 1e-10_dp
    logical :: multiplicative = .false., jacobian_check = .true., dump_matrix = .false.
    real(dp) :: coef(nbasis, 2), sc, kappa = 1e-3_dp, jac_err = 0
    logical :: smooth_profiles = .false.
    real(dp), allocatable :: err(:, :), err_fit(:, :)
    integer, allocatable :: steps(:)
    type(solovev_doublet_t) :: d
    type(block_layout_t) :: layout
    type(spline_space_t), pointer :: space_of_label => null()
    character(16) :: arg
    ! state of one mesh level
    type(spline_space_t), target :: space_c, space_g
    type(spline_derham_t) :: derham_c, derham_g
    type(spline_quadrature_t) :: aquad
    real(dp), allocatable :: rows(:, :), mu(:)
    integer, allocatable :: fixed(:)
    real(dp) :: msq(2, 2, 2)
    integer :: nc, ng, ntot

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        multiplicative = trim(arg) == 'mult'
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) amplitude
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) power
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) theta
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        jacobian_check = trim(arg) == 'jacobian'
        dump_matrix = trim(arg) == 'dump'
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) kappa
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, arg)
        smooth_profiles = trim(arg) == 'smooth'
    end if
    if (command_argument_count() >= 9) then
        call get_command_argument(9, arg)
        read (arg, *) g_degree
    end if
    call d%init()
    sc = d%psi_x - d%psi_o
    call fit_profiles()
    call build_layout()
    print '(a,i2,a,i2)', 'degree of chi', degree, '  of g', g_degree
    print '(a,es9.2,a,es9.2)', 'toroidal field scale kappa', kappa, &
        '  B_phi/B_pol at the wall about', kappa*doublet_f0/(sc*1.0_dp)
    print '(a,a,a,f6.3,a,f4.1,a,f4.2)', 'update ', merge('multiplicative', 'additive      ', &
                                                         multiplicative), '  amplitude', &
        amplitude, '  grading power', power, '  layer fraction', theta
    allocate (err(2, levels), err_fit(2, levels), steps(levels))
    do level = 1, levels
        call run(cells(level), err(:, level), err_fit(:, level), steps(level))
        call check(steps(level) <= 12, 'Newton converges in at most 12 steps')
        call check(err(1, level) < 1.2_dp*err_fit(1, level), &
                   'B_phi within 20 % of the fixed-label fit on the same mesh')
        call check(err(2, level) < 1.2_dp*err_fit(2, level), &
                   'B_pol within 20 % of the projected exact label on the same mesh')
    end do
    if (jacobian_check) call check(jac_err < 1e-6_dp, 'the Newton matrix matches differences of the residual')
    do level = 2, levels
        print '(a,i3,a,2f7.2,a,2f7.2)', 'to', cells(level), ' rates of B_phi, B_pol', &
            log(err(:, level - 1)/err(:, level))/log(2.0_dp), '  (fixed-label fit', &
            log(err_fit(:, level - 1)/err_fit(:, level))/log(2.0_dp), ')'
        call check(log(err(1, level - 1)/err(1, level))/log(2.0_dp) > &
                   log(err_fit(1, level - 1)/err_fit(1, level))/log(2.0_dp) - 0.3_dp, &
                   'B_phi converges at the rate of the fixed-label fit')
    end do
    if (failures > 0) error stop 'Clebsch doublet gate E2 failed'
    print '(a)', 'PASS gate E2 (Clebsch doublet, separatrix)'

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

    ! ---------------------------------------------------------------- profiles

    subroutine basis(t, b, db, d2b, d3b)
        !! t^m (k = 2m + 1) and t^m log|t| (k = 2m + 2) and their t derivatives.
        real(dp), intent(in) :: t
        real(dp), intent(out) :: b(nbasis), db(nbasis), d2b(nbasis), d3b(nbasis)
        real(dp) :: l, tm(-3:4)
        integer :: m, k
        l = log(max(abs(t), 1e-300_dp))
        do m = -3, 4
            tm(m) = sign(max(abs(t), 1e-300_dp), t)**m
        end do
        do k = 1, nbasis
            m = (k - 1)/2
            if (mod(k, 2) == 1) then
                b(k) = tm(m)
                db(k) = m*tm(m - 1)
                d2b(k) = m*(m - 1)*tm(m - 2)
                d3b(k) = m*(m - 1)*(m - 2)*tm(m - 3)
            else
                b(k) = tm(m)*l
                db(k) = m*tm(m - 1)*l + tm(m - 1)
                d2b(k) = m*(m - 1)*tm(m - 2)*l + (2*m - 1)*tm(m - 2)
                d3b(k) = m*(m - 1)*(m - 2)*tm(m - 3)*l + (3*m*m - 6*m + 2)*tm(m - 3)
            end if
        end do
    end subroutine basis

    subroutine fit_profiles()
        !! dPhi/dchi per side (1 lobe, t < 0; 2 outside, t > 0) by least
        !! squares (quadruple precision normal equations, columns scaled).
        real(qp) :: a(nbasis, nbasis), r(nbasis), f
        real(dp) :: b(nbasis), db(nbasis), d2b(nbasis), d3b(nbasis), t, worst, v
        real(dp), allocatable :: ts(:), vs(:)
        integer :: side, k, i, j, n
        do side = 1, 2
            n = 0
            allocate (ts(size(d%tab_s)), vs(size(d%tab_s)))
            do k = 1, size(d%tab_s)
                t = 10**d%tab_s(k)
                if (side == 1 .and. t > 0.99_dp) cycle
                ! the traced table is noisy below 1e-6; the fit carries the
                ! asymptotic law A log|t| + B inwards
                if (t < 1e-6_dp) cycle
                n = n + 1
                ts(n) = merge(-t, t, side == 1)
                vs(n) = merge(d%tab_lobe(k), d%tab_out(k), side == 1)
            end do
            a = 0
            r = 0
            do k = 1, n
                call basis(ts(k), b, db, d2b, d3b)
                do i = 1, nbasis
                    r(i) = r(i) + real(b(i), qp)*vs(k)
                    do j = 1, nbasis
                        a(i, j) = a(i, j) + real(b(i), qp)*b(j)
                    end do
                end do
            end do
            ! Gauss-Jordan in quadruple precision
            do i = 1, nbasis
                j = maxloc(abs(a(i:, i)), 1) + i - 1
                if (j /= i) then
                    a([i, j], :) = a([j, i], :)
                    r([i, j]) = r([j, i])
                end if
                do j = 1, nbasis
                    if (j == i) cycle
                    f = a(j, i)/a(i, i)
                    a(j, :) = a(j, :) - f*a(i, :)
                    r(j) = r(j) - f*r(i)
                end do
            end do
            do i = 1, nbasis
                coef(i, side) = real(r(i)/a(i, i), dp)
            end do
            worst = 0
            do k = 1, n
                call basis(ts(k), b, db, d2b, d3b)
                v = dot_product(coef(:, side), b)
                worst = max(worst, abs(v - vs(k))/abs(vs(k)))
            end do
            if (smooth_profiles) then
                ! diagnostic: the log terms removed (F smooth at the separatrix)
                coef(2:nbasis:2, side) = 0
            end if
            print '(a,i2,a,es9.2,a,2es11.3)', 'profile fit side', side, '  max rel error', worst, &
                '  log coefficient, constant', coef(2, side), coef(1, side)
            deallocate (ts, vs)
        end do
    end subroutine fit_profiles

    subroutine flux_slopes(region, c, fp)
        !! fp(1:3, j) = (F_j', F_j'', F_j''') of the region at chi = c (the
        !! side of the region: lobes t < 0, outside t > 0, clamped to it).
        integer, intent(in) :: region
        real(dp), intent(in) :: c
        real(dp), intent(out) :: fp(3, 2)
        real(dp) :: b(nbasis), db(nbasis), d2b(nbasis), d3b(nbasis), t, v(3)
        integer :: side
        side = merge(1, 2, region <= 2)
        t = c - 1
        if (side == 1) t = min(t, -1e-300_dp)
        if (side == 2) t = max(t, 1e-300_dp)
        call basis(t, b, db, d2b, d3b)
        v = [dot_product(coef(:, side), b), dot_product(coef(:, side), db), &
             dot_product(coef(:, side), d2b)]
        v = kappa*v
        fp = 0
        if (region <= 2) then
            fp(:, region) = -v/(2*pi)
        else
            fp(:, 1) = -v/(4*pi)
            fp(:, 2) = -v/(4*pi)
        end if
    end subroutine flux_slopes

    ! ------------------------------------------------------------------ mesh

    subroutine build_layout()
        type(doublet_curves_t) :: curves
        integer :: st, k
        curves%xpoint = d%xx
        curves%opoint = d%xo
        do k = 1, 2
            call d%lobe_samples(k, 2e-3_dp, curves%lobe(k)%points)
        end do
        call d%wall_samples(curves%wall)
        call layout_doublet(curves, 0.4_dp, layout, st)
        if (st /= 0) error stop 'doublet layout'
    end subroutine build_layout

    pure integer function region_of(block)
        integer, intent(in) :: block
        region_of = 3
        if (block <= 10) region_of = (block - 1)/5 + 1
    end function region_of

    subroutine build_space(n, cut, space, derham, p)
        integer, intent(in) :: n, p
        logical, intent(in) :: cut
        type(spline_space_t), intent(out) :: space
        type(spline_derham_t), intent(out) :: derham
        type(spline_space_options_t) :: options
        integer :: st
        options%degree = p
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_grade_knots(layout, LAYOUT_SEPARATRIX, 1.0_dp, options%knots, st, power=power)
        if (st /= 0) error stop 'grading'
        if (theta > 0) call layer_grading(options%knots, n)
        if (cut) options%cut_tags = [LAYOUT_SEPARATRIX]
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_derham_build(space, derham, st, 0)
        if (st /= 0) error stop 'derham'
    end subroutine build_space

    subroutine layer_grading(knots, n)
        type(knot_spec_t), intent(inout) :: knots(:, :)
        integer, intent(in) :: n
        logical :: at0, at1
        integer :: k, dd
        real(dp) :: mirrored(n/2 + 1)
        do k = 1, size(knots, 2)
            do dd = 1, 2
                if (.not. allocated(knots(dd, k)%breakpoints)) cycle
                associate (b => knots(dd, k)%breakpoints)
                    at0 = b(2) - b(1) < 0.5_dp/n
                    at1 = b(n + 1) - b(n) < 0.5_dp/n
                end associate
                if (at0 .and. at1) then
                    mirrored = 1 - 0.5_dp*layered_reverse(n/2)
                    knots(dd, k)%breakpoints = [0.5_dp*layered(n/2), mirrored(2:)]
                else if (at0) then
                    knots(dd, k)%breakpoints = layered(n)
                else if (at1) then
                    knots(dd, k)%breakpoints = 1 - [layered_reverse(n)]
                end if
            end do
        end do
    end subroutine layer_grading

    function layered(n) result(x)
        integer, intent(in) :: n
        real(dp) :: x(n + 1), tau
        integer :: m, i
        m = max(1, nint(theta*n))
        if (m >= n) then
            x = [((real(i, dp)/n)**power, i=0, n)]
            return
        end if
        tau = real(m, dp)/n/(power*(1 - real(m, dp)/n) + real(m, dp)/n)
        do i = 0, m
            x(i + 1) = tau*(real(i, dp)/m)**power
        end do
        do i = m + 1, n
            x(i + 1) = tau + (1 - tau)*real(i - m, dp)/(n - m)
        end do
    end function layered

    function layered_reverse(n) result(x)
        integer, intent(in) :: n
        real(dp) :: x(n + 1), y(n + 1)
        integer :: i
        y = layered(n)
        x = [(y(n + 2 - i), i=1, n + 1)]
    end function layered_reverse

    ! ---------------------------------------------------------------- solve

    subroutine run(n, e, e_fit, nsteps)
        integer, intent(in) :: n
        real(dp), intent(out) :: e(2), e_fit(2)
        integer, intent(out) :: nsteps
        type(spline_quadrature_t) :: quad, q1
        real(dp), allocatable :: chi(:), g(:), chi0(:), xg(:, :), r(:), rt(:), delta(:), &
                                 history(:), gfit(:)
        integer, allocatable :: wall(:), sep(:)
        real(dp) :: alpha, t0, t1, vmax
        integer :: st, i, it, k, j
        type(csr_t) :: jac

        call cpu_time(t0)
        call build_space(n, .false., space_c, derham_c, degree)
        call build_space(n, .true., space_g, derham_g, g_degree)
        nc = derham_c%dimension(0)
        ng = derham_g%dimension(0)
        ntot = nc + ng + 6
        call spline_derham_boundary_dofs(derham_c, LAYOUT_WALL, 0, wall)
        call spline_derham_boundary_dofs(derham_c, LAYOUT_SEPARATRIX, 0, sep)
        fixed = [wall, sep]
        do j = 1, 2
            msq(:, :, j) = sqrt_hessian(d%xo(:, j))
        end do
        call spline_quadrature_build(space_g, g_degree + 8 + mod(g_degree, 2), quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_quadrature_point_singular(space_g, quad, d%xo(:, 1), q1, st)
        if (st /= 0) error stop 'polar O1'
        call spline_quadrature_point_singular(space_g, q1, d%xo(:, 2), aquad, st)
        if (st /= 0) error stop 'polar O2'
        ! constraint rows: chi, d chi/dR, d chi/dZ at O_1, O_2
        if (allocated(mu)) deallocate (mu)
        if (allocated(rows)) deallocate (rows)
        allocate (rows(nc, 6), source=0.0_dp)
        do j = 1, 2
            block
                integer, allocatable :: rp(:), ci(:)
                real(dp), allocatable :: v(:)
                real(dp) :: xi(2)
                integer :: b, c
                call spline_locate(space_c, 1, d%xo(:, j), b, xi, st)
                if (st /= 0) error stop 'O location'
                do c = 0, 2
                    call spline_point_matrix(space_c, [b], reshape(xi, [2, 1]), c, rp, ci, v, st)
                    if (st /= 0) error stop 'O rows'
                    rows(ci, 3*(j - 1) + c + 1) = v
                end do
            end block
        end do

        ! the projected exact label, separatrix at 1, wall at 2
        space_of_label => space_c
        allocate (chi0(nc), g(ng), mu(6), source=0.0_dp)
        call spline_derham_project(derham_c, 0, label_cb, chi0, st, 12)
        if (st /= 0) error stop 'label projection'
        chi0(sep) = 1
        chi0(wall) = 2

        ! the fixed-label fit of g (the E2b level) and its errors
        call solve_g(chi0, g)
        gfit = g
        e_fit = field_errors(chi0, gfit)

        ! start: labels perturbed on their side of the separatrix, g fitted
        call spline_geometry_dofs(space_c, xg, st)
        chi = chi0
        do i = 1, nc
            if (any(fixed == i)) cycle
            chi(i) = 1 + (chi0(i) - 1)*(1 + amplitude*(xg(1, i) - 1 + xg(2, i)))
        end do
        call solve_g(chi, g)
        mu = 0
        allocate (history(0))
        nsteps = -1
        do it = 0, max_steps
            call system(chi, g, mu, r, jac)
            history = [history, norm2(r)]
            if (it > 0 .and. history(it + 1) < 1e-10_dp*history(1)) then
                nsteps = it
                exit
            end if
            if (it == max_steps) exit
            call scaled_solve(jac, -r, fixed, delta, st)
            if (st /= 0) error stop 'Newton solve'
            if (it == 0 .and. dump_matrix) then
                block
                    integer :: u, ii, pp
                    open (newunit=u, file='newton_matrix.txt', status='replace')
                    write (u, '(3i8)') ntot, nc, ng
                    do ii = 1, ntot
                        do pp = jac%row_ptr(ii), jac%row_ptr(ii + 1) - 1
                            write (u, '(2i8,es26.17)') ii, jac%col(pp), jac%val(pp)
                        end do
                    end do
                    close (u)
                    open (newunit=u, file='newton_aux.txt', status='replace')
                    do ii = 1, nc
                        write (u, '(i8,4es26.17,i3)') ii, xg(:, ii), chi(ii), delta(ii), &
                            merge(1, 0, any(fixed == ii))
                    end do
                    close (u)
                end block
            end if
            if (it == 0 .and. jacobian_check .and. level == 1) then
                block
                    real(dp), allocatable :: v(:), rp(:), rm(:), jv(:)
                    real(dp) :: hv
                    integer :: blk, lo, hi, ii
                    do blk = 1, 3
                        allocate (v(ntot), source=0.0_dp)
                        lo = merge(1, merge(nc + 1, nc + ng + 1, blk == 2), blk == 1)
                        hi = merge(nc, merge(nc + ng, ntot, blk == 2), blk == 1)
                        do ii = lo, hi
                            v(ii) = sin(1.7_dp*ii) + 0.3_dp
                        end do
                        v(fixed) = 0
                        ! scale chi directions to the layer: relative to |chi - 1|
                        if (blk == 1) v(1:nc) = v(1:nc)*abs(chi - 1)
                        hv = 1e-6_dp/maxval(abs(v))
                        call system(chi + hv*v(1:nc), g + hv*v(nc + 1:nc + ng), mu + hv*v(nc + ng + 1:), rp)
                        call system(chi - hv*v(1:nc), g - hv*v(nc + 1:nc + ng), mu - hv*v(nc + ng + 1:), rm)
                        jv = matvec(jac, v)
                        jv(fixed) = 0
                        jac_err = max(jac_err, norm2(jv - (rp - rm)/(2*hv))/norm2((rp - rm)/(2*hv)))
                        print '(a,i2,a,es10.2,a,es10.2,a,3es10.2)', '  Jacobian check block', blk, &
                            ': |Jv - FD|/|FD|', norm2(jv - (rp - rm)/(2*hv))/norm2((rp - rm)/(2*hv)), &
                            '  |FD|', norm2((rp - rm)/(2*hv)), '  rows chi g mu', &
                            norm2(jv(1:nc) - (rp(1:nc) - rm(1:nc))/(2*hv)), &
                            norm2(jv(nc + 1:nc + ng) - (rp(nc + 1:nc + ng) - rm(nc + 1:nc + ng))/(2*hv)), &
                            norm2(jv(nc + ng + 1:) - (rp(nc + ng + 1:) - rm(nc + ng + 1:))/(2*hv))
                        deallocate (v)
                    end do
                end block
            end if
            alpha = 1
            if (multiplicative) then
                ! at most a factor e per step on |chi_i - 1|
                vmax = 0
                do i = 1, nc
                    if (abs(chi(i) - 1) > 0 .and. .not. any(fixed == i)) &
                        vmax = max(vmax, abs(delta(i)/(chi(i) - 1)))
                end do
                alpha = min(1.0_dp, 1/max(vmax, tiny(1.0_dp)))
                if (it <= 1) then
                    block
                        real(dp) :: ratio(nc)
                        integer :: top, ii
                        ratio = 0
                        do ii = 1, nc
                            if (.not. any(fixed == ii)) ratio(ii) = abs(delta(ii)/(chi(ii) - 1))
                        end do
                        do top = 1, 6
                            ii = maxloc(ratio, 1)
                            print '(a,i6,a,2f9.5,a,es10.2,a,es10.2,a,es10.2)', '      node', ii, &
                                '  at', xg(:, ii), '  chi - 1', chi(ii) - 1, '  d chi', delta(ii), &
                                '  ratio', ratio(ii)
                            ratio(ii) = 0
                        end do
                        ii = maxloc(abs(delta(1:nc)), 1)
                        print '(a,2f9.5,a,es10.2,a,es10.2)', '      max |d chi| at', xg(:, ii), &
                            '  chi', chi(ii), '  d chi', delta(ii)
                        print '(a,es10.2,a,es10.2,a,6es10.2)', '      max |d chi|', maxval(abs(delta(1:nc))), &
                            '  max |d g|', maxval(abs(delta(nc + 1:nc + ng))), '  d mu', delta(nc + ng + 1:)
                    end block
                end if
            end if
            do k = 1, 20
                call system(updated(chi, delta(1:nc), alpha), g + alpha*delta(nc + 1:nc + ng), &
                            mu + alpha*delta(nc + ng + 1:), rt)
                if (norm2(rt) < (1 - alpha/4)*history(it + 1)) exit
                alpha = alpha/2
            end do
            print '(a,i3,a,es10.3,a,es9.2,a,es10.3)', '    step', it, '  |r|', history(it + 1), &
                '  alpha', alpha, '  trial |r|', norm2(rt)
            chi = updated(chi, delta(1:nc), alpha)
            g = g + alpha*delta(nc + 1:nc + ng)
            mu = mu + alpha*delta(nc + ng + 1:)
        end do
        call cpu_time(t1)
        e = field_errors(chi, g)
        if (nsteps < 0) nsteps = 99
        print '(a,i3,a,i3,a,2es10.3,a,2es10.3,a,es9.2,a,f7.2,a)', 'cells', n, '  Newton steps', &
            nsteps, '  B_phi, B_pol error', e, '  (fixed-label fit', e_fit, &
            ')  max |chi - chi0|', maxval(abs(chi - chi0)), '  ', t1 - t0, ' s'
        print '(a,16es9.1)', '  |r| history', history(1:min(16, size(history)))

    end subroutine run

    function updated(c, dc, a) result(cn)
        real(dp), intent(in) :: c(:), dc(:), a
        real(dp) :: cn(size(c))
        integer :: ii
        if (.not. multiplicative) then
            cn = c + a*dc
            cn(fixed) = c(fixed)
            return
        end if
        cn = c
        do ii = 1, size(c)
            if (any(fixed == ii)) cycle
            if (abs(c(ii) - 1) > 0) then
                cn(ii) = 1 + (c(ii) - 1)*exp(a*dc(ii)/(c(ii) - 1))
            else
                cn(ii) = c(ii) + a*dc(ii)
            end if
        end do
    end function updated

    subroutine solve_g(c, gg)
        !! g at fixed labels: W is quadratic and convex in g.
        real(dp), intent(in) :: c(:)
        real(dp), intent(inout) :: gg(:)
        real(dp), allocatable :: rr(:), dd(:)
        type(csr_t) :: jj
        integer :: ii, stg
        call system(c, gg, mu, rr, jj)
        call scaled_solve(jj, -rr, [[(ii, ii=1, nc)], [(nc + ng + ii, ii=1, 6)]], dd, stg)
        if (stg /= 0) error stop 'g solve'
        gg = gg + dd(nc + 1:nc + ng)
    end subroutine solve_g

    subroutine system(c, gg, m, rr, jj)
        !! Residual (grad W + constraint terms; constraint rows) and, with
        !! jj, the Newton matrix (symmetric).
        real(dp), intent(in) :: c(:), gg(:), m(:)
        real(dp), allocatable, intent(out) :: rr(:)
        type(csr_t), intent(out), optional :: jj
        real(dp), allocatable :: vc(:, :, :), gc(:, :, :), vg(:, :, :), gq(:, :, :), &
                                 db(:), hl(:, :)
        integer, allocatable :: dc(:), dg(:), loc(:)
        real(dp) :: s(2), jet(2, 6), w, cv, gx(2), ggx(2), fp(3, 2), u(2), bphi, h(2, 2), &
                    ex(2), kh(2), a2, kn(2)
        integer :: cell, q, p, ii, kk, jw, nl, ncl, region, stq
        type(triplet_list) :: tp
        allocate (rr(ntot), source=0.0_dp)
        if (present(jj)) call tp%init(size(aquad%point_offset)*((degree + 1)**2 + (g_degree + 1)**2)**2 + &
                                      12*nc + ntot)
        do cell = 1, size(aquad%point_offset)
            call clebsch_node_gradients(derham_c, aquad, cell, vc, dc, gc)
            call clebsch_node_gradients(derham_g, aquad, cell, vg, dg, gq)
            ncl = size(dc)
            nl = ncl + size(dg)
            loc = [dc, nc + dg]
            if (allocated(db)) deallocate (db, hl)
            allocate (db(nl), hl(nl, nl))
            hl = 0
            do q = 1, aquad%points_per_cell
                p = aquad%point_offset(cell) + q - 1
                call spline_geometry_at(space_g, aquad%block(p), aquad%xi(:, p), jet, stq)
                s = jet(:, 1)
                w = 2*pi*s(1)*aquad%weight(p)
                region = region_of(aquad%block(p))
                cv = dot_product(vc(1, :, q), c(dc))
                gx = matmul(gc(:, :, q), c(dc))
                ggx = matmul(gq(:, :, q), gg(dg))
                call flux_slopes(region, cv, fp)
                do jw = 1, 2
                    ex = matmul(msq(:, :, jw), s - d%xo(:, jw))
                    h(:, jw) = matmul(msq(:, :, jw), [-ex(2), ex(1)])/sum(ex**2)
                end do
                u = ggx - fp(1, 1)*h(:, 1) - fp(1, 2)*h(:, 2)
                bphi = cross(gx, u)
                ! K(grad chi, h_j), weighted by F_j'' and F_j'''
                kh = [cross(gx, h(:, 1)), cross(gx, h(:, 2))]
                a2 = sc**2/s(1)**2
                ! first derivatives of B_phi
                do ii = 1, ncl
                    kn = [cross(gc(:, ii, q), h(:, 1)), cross(gc(:, ii, q), h(:, 2))]
                    db(ii) = cross(gc(:, ii, q), u) - vc(1, ii, q)*(fp(2, 1)*kh(1) + fp(2, 2)*kh(2))
                end do
                do ii = 1, size(dg)
                    db(ncl + ii) = cross(gx, gq(:, ii, q))
                end do
                ! gradient: s^2 grad chi . grad N / R^2 + B dB + s N (p = -s chi)
                rr(dc) = rr(dc) + w*(a2*matmul(gx, gc(:, :, q)) + bphi*db(1:ncl) + &
                                     sc*vc(1, :, q))
                rr(nc + dg) = rr(nc + dg) + w*bphi*db(ncl + 1:)
                if (present(jj)) then
                    do kk = 1, nl
                        hl(:, kk) = hl(:, kk) + w*db*db(kk)
                    end do
                    do kk = 1, ncl
                        kn = [cross(gc(:, kk, q), h(:, 1)), cross(gc(:, kk, q), h(:, 2))]
                        do ii = 1, ncl
                            hl(ii, kk) = hl(ii, kk) + w*(a2*dot_product(gc(:, ii, q), gc(:, kk, q)) &
                                + bphi*( &
                                -vc(1, kk, q)*(fp(2, 1)*cross(gc(:, ii, q), h(:, 1)) + &
                                               fp(2, 2)*cross(gc(:, ii, q), h(:, 2))) &
                                - vc(1, ii, q)*(fp(2, 1)*kn(1) + fp(2, 2)*kn(2)) &
                                - vc(1, ii, q)*vc(1, kk, q)*(fp(3, 1)*kh(1) + fp(3, 2)*kh(2))))
                        end do
                        do ii = 1, size(dg)
                            hl(ncl + ii, kk) = hl(ncl + ii, kk) + &
                                               w*bphi*cross(gc(:, kk, q), gq(:, ii, q))
                            hl(kk, ncl + ii) = hl(kk, ncl + ii) + &
                                               w*bphi*cross(gc(:, kk, q), gq(:, ii, q))
                        end do
                    end do
                end if
            end do
            if (present(jj)) then
                do kk = 1, nl
                    do ii = 1, nl
                        call tp%add(loc(ii), loc(kk), hl(ii, kk))
                    end do
                end do
            end if
        end do
        ! constraints chi_h(O_j) = 0, grad chi_h(O_j) = 0
        do kk = 1, 6
            rr(1:nc) = rr(1:nc) + m(kk)*rows(:, kk)
            rr(nc + ng + kk) = dot_product(rows(:, kk), c)
            if (present(jj)) then
                do ii = 1, nc
                    if (rows(ii, kk) /= 0) then
                        call tp%add(ii, nc + ng + kk, rows(ii, kk))
                        call tp%add(nc + ng + kk, ii, rows(ii, kk))
                    end if
                end do
            end if
        end do
        rr(fixed) = 0
        if (present(jj)) then
            do ii = 1, ntot
                call tp%add(ii, ii, 0.0_dp)
            end do
            call tp%to_csr(ntot, ntot, jj)
        end if
    end subroutine system

    subroutine scaled_solve(m, rhs, fix, x, status)
        !! Symmetric Jacobi scaling of the chi and g rows (the graded cells
        !! spread the diagonal over decades), a ridge on the g block of
        !! the scaled matrix (g -> g + k(chi) and the region constants),
        !! iterative refinement against the unscaled matrix.
        type(csr_t), intent(in) :: m
        real(dp), intent(in) :: rhs(:)
        integer, intent(in) :: fix(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        real(dp) :: scale(ntot), rr(ntot), dx(ntot)
        type(csr_t) :: ms
        integer :: ii, pp, kt
        ms = m
        scale = 1
        do ii = 1, nc + ng
            do pp = m%row_ptr(ii), m%row_ptr(ii + 1) - 1
                if (m%col(pp) == ii .and. abs(m%val(pp)) > 0) scale(ii) = 1/sqrt(abs(m%val(pp)))
            end do
        end do
        do ii = 1, ntot
            do pp = ms%row_ptr(ii), ms%row_ptr(ii + 1) - 1
                ms%val(pp) = scale(ii)*ms%val(pp)*scale(ms%col(pp))
                if (ms%col(pp) == ii .and. ii > nc .and. ii <= nc + ng) &
                    ms%val(pp) = ms%val(pp) + ridge
            end do
        end do
        allocate (x(ntot))
        call eu_constrained_solve(ms, scale*rhs, fix, 0*real(fix, dp), x, status)
        if (status /= 0) return
        x = scale*x
        do kt = 1, 4
            rr = rhs - matvec(m, x)
            rr(fix) = 0
            if (norm2(rr) <= 1e-12_dp*norm2(rhs)) exit
            call eu_constrained_solve(ms, scale*rr, fix, 0*real(fix, dp), dx, status)
            if (status /= 0) return
            x = x + scale*dx
        end do
    end subroutine scaled_solve

    function field_errors(c, gg) result(er)
        !! Relative L2 errors (weight 2 pi R) of B_phi against f0/R and of
        !! B_pol against the exact poloidal field.
        real(dp), intent(in) :: c(:), gg(:)
        real(dp) :: er(2), num(2), den(2), s(2), jet(2, 6), w, cv, gx(2), ggx(2), fp(3, 2), &
                    h(2, 2), ex(2), bphi, ge(2)
        real(dp), allocatable :: vc(:, :, :), gc(:, :, :), vg(:, :, :), gq(:, :, :)
        integer, allocatable :: dc(:), dg(:)
        integer :: cell, q, p, jw, stq
        num = 0
        den = 0
        do cell = 1, size(aquad%point_offset)
            call clebsch_node_gradients(derham_c, aquad, cell, vc, dc, gc)
            call clebsch_node_gradients(derham_g, aquad, cell, vg, dg, gq)
            do q = 1, aquad%points_per_cell
                p = aquad%point_offset(cell) + q - 1
                call spline_geometry_at(space_g, aquad%block(p), aquad%xi(:, p), jet, stq)
                s = jet(:, 1)
                w = 2*pi*s(1)*aquad%weight(p)
                cv = dot_product(vc(1, :, q), c(dc))
                gx = matmul(gc(:, :, q), c(dc))
                ggx = matmul(gq(:, :, q), gg(dg))
                call flux_slopes(region_of(aquad%block(p)), cv, fp)
                do jw = 1, 2
                    ex = matmul(msq(:, :, jw), s - d%xo(:, jw))
                    h(:, jw) = matmul(msq(:, :, jw), [-ex(2), ex(1)])/sum(ex**2)
                end do
                bphi = cross(gx, ggx - fp(1, 1)*h(:, 1) - fp(1, 2)*h(:, 2))
                ge = d%label_gradient(s(1), s(2))
                num(1) = num(1) + w*(bphi - kappa*doublet_f0/s(1))**2
                den(1) = den(1) + w*(kappa*doublet_f0/s(1))**2
                num(2) = num(2) + w*sc**2*sum((gx - ge)**2)/s(1)**2
                den(2) = den(2) + w*sc**2*sum(ge**2)/s(1)**2
            end do
        end do
        er = sqrt(num/den)
    end function field_errors


    pure real(dp) function cross(a, b)
        real(dp), intent(in) :: a(2), b(2)
        cross = a(1)*b(2) - a(2)*b(1)
    end function cross

    function matvec(m, x) result(y)
        type(csr_t), intent(in) :: m
        real(dp), intent(in) :: x(:)
        real(dp) :: y(m%nrow)
        integer :: i, p
        y = 0
        do i = 1, m%nrow
            do p = m%row_ptr(i), m%row_ptr(i + 1) - 1
                y(i) = y(i) + m%val(p)*x(m%col(p))
            end do
        end do
    end function matvec

    function sqrt_hessian(x) result(m)
        !! H^(1/2) of the exact label at x.
        real(dp), intent(in) :: x(2)
        real(dp) :: m(2, 2), hm(2, 2), sd
        real(dp), parameter :: h = 1e-4_dp
        hm(:, 1) = (d%label_gradient(x(1) + h, x(2)) - d%label_gradient(x(1) - h, x(2)))/(2*h)
        hm(:, 2) = (d%label_gradient(x(1), x(2) + h) - d%label_gradient(x(1), x(2) - h))/(2*h)
        hm(1, 2) = (hm(1, 2) + hm(2, 1))/2
        hm(2, 1) = hm(1, 2)
        sd = sqrt(hm(1, 1)*hm(2, 2) - hm(1, 2)**2)
        m = hm
        m(1, 1) = m(1, 1) + sd
        m(2, 2) = m(2, 2) + sd
        m = m/sqrt(hm(1, 1) + hm(2, 2) + 2*sd)
    end function sqrt_hessian

    subroutine label_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st
        call spline_geometry_at(space_of_label, block, xi, jet, st)
        values(1) = d%label(jet(1, 1), jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_cb

end program test_label_ale_doublet_clebsch
