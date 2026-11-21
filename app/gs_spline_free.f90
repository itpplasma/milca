program gs_spline_free_main
    !! Free-boundary Grad-Shafranov on the spline O-grid with a periodic
    !! annulus (milca_gs_spline): the core blocks are the plasma (law sources),
    !! the annulus up to the wall is vacuum, the wall couples to the coils and
    !! the exterior (FortFEM Steklov-Poincare). psi_b is the mean of psi at the
    !! target points. ALE: after each Newton solve the interface Gamma_p moves
    !! onto the contour psi = psi_b, until it moves less than move_tolerance.
    !!
    !! Namelist &gs_spline_free: wall_file, interface_file, target_file (count,
    !! then R Z rows; curves counter-clockwise from the outboard midplane),
    !! profile_file (normalised law), coil_file, currents(:), seed_file (nr nz,
    !! R(1:nr), Z(1:nz), psi(nr, nz) R fastest), degree, cells, annulus_cells,
    !! interface_fraction, max_moves, move_tolerance, mu0, relaxation (of the
    !! Gamma_p move), quadrature_points (0: degree + 1), edge_grading (radial
    !! grading of the core ring blocks towards Gamma_p, FortFEM SPLINE_GRID_*;
    !! default S2), edge_reversed.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use fortnum_polynomial, only: lagrange_weights
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_locate, spline_evaluate_at
    use milca_gs_profile_law, only: read_gs_profile_law, gs_law_sources, gs_law_primitives
    use milca_gs_coils, only: gs_circuit_set_t, gs_read_coil_table, gs_circuit_flux_vectors
    use milca_gs_spline, only: gs_spline_law_kernel_t, gs_spline_vacuum_t, &
        gs_spline_vacuum_setup, gs_spline_newton_free, gs_spline_project, ray_value, &
        GS_SPLINE_OK
    use milca_gs_free_current, only: gs_current_target_t, gs_read_current_target, &
        hermite_law, hermite_end, push_history, least_squares, GS_CUR_OK
    use milca_gs_spline_current, only: gs_spline_current_update, GS_SPLINE_CUR_OK
    use milca_gs_profile_law, only: gs_profile_law_t
    implicit none
    integer, parameter :: max_circuits = 64, contour_samples = 256
    character(512) :: wall_file = 'wall.dat', interface_file = 'interface.dat', &
        target_file = 'target.dat', profile_file = 'profile.law', coil_file = 'coils.txt', &
        seed_file = 'seed.dat', config = 'input.nml', message, current_file = ''
    integer :: current_passes = 60, current_knots = 129
    real(dp) :: current_tolerance = 1.0e-6_dp
    real(dp) :: currents(max_circuits) = 0, interface_fraction = 0.8_dp, &
        move_tolerance = 1.0e-7_dp, mu0 = 4.0e-7_dp*acos(-1.0_dp), tolerance = 1.0e-11_dp, &
        relaxation = 0.5_dp
    integer :: degree = 3, cells = 8, annulus_cells = 8, max_moves = 30, quadrature_points = 0, &
        edge_grading = 2
    logical :: edge_reversed = .false.
    namelist /gs_spline_free/ wall_file, interface_file, target_file, profile_file, &
        coil_file, currents, seed_file, degree, cells, annulus_cells, interface_fraction, &
        max_moves, move_tolerance, mu0, tolerance, relaxation, quadrature_points, &
        edge_grading, edge_reversed, current_file, current_passes, current_tolerance, &
        current_knots
    type(block_layout_t) :: layout
    type(spline_space_options_t) :: options
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad
    type(gs_spline_law_kernel_t) :: kernel
    type(gs_spline_vacuum_t) :: vac
    type(gs_circuit_set_t) :: set
    real(dp), allocatable :: wall(:, :), gamma(:, :), targets(:, :), coef(:, :), &
        psi_c(:), dn(:), txi(:, :), jet(:, :), seed_r(:), seed_z(:), seed_psi(:, :), &
        values(:), contour(:, :)
    integer, allocatable :: tblock(:), ray_block(:)
    real(dp), allocatable :: ray_radius(:)
    real(dp) :: axis(2), residual, centre(2), radius, volume, current, pp, ff, move, psi_b, &
        phi_area
    integer :: u, st, iterations, nc, move_count, q, k
    integer(8) :: t0, t1, rate

    if (command_argument_count() >= 2) call get_command_argument(2, config)
    open (newunit=u, file=trim(config), status='old', action='read')
    read (u, nml=gs_spline_free)
    close (u)
    call system_clock(t0, rate)
    call read_curve(trim(wall_file), wall)
    call read_curve(trim(interface_file), gamma)
    call read_curve(trim(target_file), targets)
    call read_gs_profile_law(trim(profile_file), kernel%law, st)
    if (st /= 0 .or. .not. kernel%law%normalized) error stop 'invalid normalised profile_file'
    call gs_read_coil_table(trim(coil_file), set, st, message)
    if (st /= 0) error stop 'invalid coil table: '//trim(message)
    nc = size(set%incidence, 2)
    call read_seed(trim(seed_file))
    centre = [0.5_dp*(minval(wall(1, :)) + maxval(wall(1, :))), &
              0.5_dp*(minval(wall(2, :)) + maxval(wall(2, :)))]
    radius = 0.5_dp*(maxval(wall(1, :)) - minval(wall(1, :)))
    call layout_ogrid(centre, radius, 0.5_dp, layout, st, &
                      interface_fraction=interface_fraction)
    if (st /= 0) error stop 'layout'
    call layout_set_curve(layout, LAYOUT_WALL, wall, st)
    if (st /= 0) error stop 'wall curve'
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, gamma, st)
    if (st /= 0) error stop 'interface curve'
    options%degree = degree
    options%regularity = SPLINE_C0
    allocate (options%knots(2, layout%block_count))
    options%knots(:, :)%element_count = cells
    options%knots(2, 6)%element_count = annulus_cells
    ! radial grading of the core ring blocks towards Gamma_p, where the
    ! profiles are steep
    options%knots(2, 2:5)%grading = edge_grading
    options%knots(2, 2:5)%reversed = edge_reversed
    call spline_space_build(layout, options, space, st)
    if (st /= 0) error stop 'space'
    kernel%order = 1
    kernel%mu0 = mu0
    call rebuild(.true.)
    ! seed: the reference flux map projected onto the spline space
    allocate (coef(1, space%dof_count))
    allocate (values(quad%point_count))
    do q = 1, quad%point_count
        values(q) = seed_value(quad%s(1, q), quad%s(2, q))
    end do
    call gs_spline_project(space, quad, values, coef, st)
    if (st /= 0) error stop 'seed projection'
    do move_count = 0, max_moves
        call gs_spline_newton_free(space, quad, kernel, vac, psi_c, dn, 1.0_dp, tblock, txi, &
                                   coef, tolerance, 30, iterations, residual, axis, st)
        write (output_unit, '(a,i3,a,i0,a,i0,a,es10.3,a,2f12.8)') 'move', move_count, &
            ' Newton status ', st, ' iterations ', iterations, ' residual ', residual, &
            ' axis', axis
        if (st /= GS_SPLINE_OK) error stop 'free-boundary Newton failed'
        psi_b = kernel%psi_b
        call trace_contour(psi_b, centre, contour, st)
        if (st /= 0) error stop 'contour trace failed'
        move = maxval(norm2(contour - gamma, dim=1))
        write (output_unit, '(a,es10.3,a,f14.10)') '  Gamma_p moves', move, ' psi_b', psi_b
        if (move <= move_tolerance) exit
        ! The staggered map Gamma_p -> contour has an eigenvalue near -1 on
        ! DIII-D (a two-cycle of 1.5 mm without relaxation).
        ! The first move only reparametrises Gamma_p to the ray samples.
        if (move_count == 0) then
            gamma = contour
        else
            gamma = gamma + relaxation*(contour - gamma)
        end if
        call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, gamma, st)
        if (st /= 0) error stop 'interface update'
        call spline_space_update_geometry(space, layout, st)
        if (st /= 0) error stop 'geometry update'
        call rebuild(.false.)
    end do
    if (len_trim(current_file) > 0) call current_mode()
    allocate (jet(1, quad%point_count))
    call spline_evaluate(space, quad, coef, 0, jet, st)
    volume = 0
    current = 0
    phi_area = 0
    do q = 1, quad%point_count
        if (.not. kernel%plasma(q)) cycle
        volume = volume + 2*acos(-1.0_dp)*quad%s(1, q)*quad%weight(q)
        call gs_law_sources(kernel%law, jet(1, q) - kernel%psi_b, pp, ff)
        current = current + (quad%s(1, q)*pp + ff/(mu0*quad%s(1, q)))*quad%weight(q)
        call gs_law_primitives(kernel%law, jet(1, q) - kernel%psi_b, pp, ff, st)
        phi_area = phi_area + abs(ff)/quad%s(1, q)*quad%weight(q)
    end do
    write (output_unit, '(a,f14.9)') 'toroidal flux inside Gamma_p (area integral)', phi_area
    call system_clock(t1)
    write (output_unit, '(a,2f13.8,a,2es18.10,a,f12.7,a,f14.3,a,i0,a,f8.3,a)') 'axis', &
        axis, ' psi_a psi_b', kernel%law%axis + kernel%psi_b, kernel%psi_b, ' volume', &
        volume, ' current', current, ' dofs ', space%dof_count, ' wall ', &
        real(t1 - t0, dp)/rate, ' s'
contains
    subroutine rebuild(with_wall)
        !! Quadrature, plasma mask and targets of the current geometry; with
        !! with_wall also the wall coupling and the coil data. A move of
        !! Gamma_p leaves the wall, the trace of the basis on it and so T, S
        !! and the coil data unchanged.
        logical, intent(in) :: with_wall
        real(dp), allocatable :: g(:, :), gr(:, :), gz(:, :)
        integer :: j
        call spline_quadrature_build(space, merge(quadrature_points, degree + 1, &
            quadrature_points > 0), quad, st)
        if (st /= 0) error stop 'quadrature'
        kernel%r = quad%s(1, :)
        kernel%plasma = quad%block <= 5
        if (with_wall) then
            call gs_spline_vacuum_setup(space, LAYOUT_WALL, 2, mu0, vac, st)
            if (st /= GS_SPLINE_OK) error stop 'vacuum setup'
            allocate (g(size(vac%nodes, 2), nc))
            call gs_circuit_flux_vectors(set, vac%nodes(1, :), vac%nodes(2, :), g)
            psi_c = matmul(g, currents(:nc))
            deallocate (g)
            allocate (g(vac%side%point_count, nc), gr(vac%side%point_count, nc), &
                      gz(vac%side%point_count, nc))
            call gs_circuit_flux_vectors(set, vac%side%s(1, :), vac%side%s(2, :), g, gr, gz)
            dn = matmul(gr, currents(:nc))*vac%normal(1, :) + &
                 matmul(gz, currents(:nc))*vac%normal(2, :)
        end if
        if (allocated(tblock)) deallocate (tblock, txi)
        allocate (tblock(size(targets, 2)), txi(2, size(targets, 2)))
        do j = 1, size(targets, 2)
            call spline_locate(space, 1, targets(:, j), tblock(j), txi(:, j), st)
            if (st /= 0) call spline_locate(space, 6, targets(:, j), tblock(j), txi(:, j), st)
            if (st /= 0) print '(a,i0,a,2f12.6,a,i0)', 'target ', j, ' at', targets(:, j), &
                ' not located, status ', st
            if (st /= 0) error stop 'target outside the mesh'
        end do
    end subroutine rebuild

    subroutine current_mode()
        !! VMEC's constraints p(s), I(s), Phi_edge from the target-mean state.
        !! Per pass one Newton solve, the profile update (law, psi_b) and the
        !! psi_b contour; law, psi_b and Gamma_p move together in one Anderson
        !! iteration (the plain map is unstable, see milca_gs_free_current).
        !! Gamma_p is scaled by 1 cm in the change measure.
        integer, parameter :: memory = 6
        type(gs_current_target_t) :: target
        type(gs_profile_law_t) :: next
        integer, allocatable :: none_block(:)
        real(dp), allocatable :: none_xi(:, :), x(:), g(:), f(:), f_old(:), g_old(:), &
            weight(:), df(:, :), dg(:, :), gamma_aa(:)
        real(dp) :: psi_a, psi_b_next, change
        integer :: pass, nk, used, ng
        real(dp), parameter :: gamma_scale = 1.0e-2_dp
        call gs_read_current_target(trim(current_file), target, st)
        if (st /= GS_CUR_OK) error stop 'invalid current_file'
        call reknot(kernel%law, current_knots)
        allocate (none_block(0), none_xi(2, 0))
        nk = size(kernel%law%knots)
        ng = size(gamma)
        allocate (df(2*nk + 1 + ng, memory), dg(2*nk + 1 + ng, memory), gamma_aa(memory))
        used = 0
        do pass = 1, current_passes
            ! one Newton solve at the current law, psi_b and Gamma_p, then the
            ! profile update and the psi_b contour: all three move together
            ! in one Anderson iteration (law, psi_b and Gamma_p samples)
            call gs_spline_newton_free(space, quad, kernel, vac, psi_c, dn, 1.0_dp, &
                                       none_block, none_xi, coef, tolerance, 30, &
                                       iterations, residual, axis, st)
            if (st /= GS_SPLINE_OK) error stop 'current mode: Newton failed'
            call trace_contour(kernel%psi_b, centre, contour, st)
            if (st /= 0) error stop 'current mode: contour trace failed'
            psi_a = kernel%law%axis + kernel%psi_b
            call gs_spline_current_update(space, coef, centre, psi_a, kernel%psi_b, 1.0_dp, &
                                          mu0, kernel%law, target, contour_samples, next, &
                                          psi_b_next, st)
            if (st /= GS_SPLINE_CUR_OK) error stop 'current mode: update failed'
            x = [pack_state(kernel%law, kernel%psi_b, nk), reshape(gamma, [ng])]
            g = [pack_state(next, psi_b_next, nk), reshape(contour, [ng])]
            if (pass == 1) weight = [spread(maxval(abs(x(1:nk))), 1, nk), &
                spread(maxval(abs(x(nk + 1:2*nk))), 1, nk), abs(psi_a - kernel%psi_b), &
                spread(gamma_scale, 1, ng)]
            x = x/weight
            g = g/weight
            f = g - x
            change = maxval(abs(f))
            move = maxval(norm2(contour - gamma, dim=1))
            if (pass > 1) call push_history(df, f - f_old, dg, g - g_old, used, memory)
            f_old = f
            g_old = g
            if (used > 0) then
                call least_squares(df(:, :used), f, gamma_aa)
                g = g - matmul(dg(:, :used), gamma_aa(:used))
            end if
            g = g*weight
            write (output_unit, '(a,i3,a,i0,a,es10.3,a,es10.3,a,f14.10)') 'current pass', &
                pass, ' Newton ', iterations, ' change', change, ' Gamma_p moves', move, &
                ' psi_b', kernel%psi_b
            if (change <= current_tolerance) exit
            call hermite_law(kernel%law%knots, g(1:nk), kernel%law%pp)
            call hermite_law(kernel%law%knots, g(nk + 1:2*nk), kernel%law%ff)
            kernel%law%p_edge = next%p_edge
            kernel%psi_b = g(2*nk + 1)
            kernel%law%axis = psi_a - kernel%psi_b
            gamma = reshape(g(2*nk + 2:), shape(gamma))
            call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, gamma, st)
            if (st /= 0) error stop 'interface update'
            call spline_space_update_geometry(space, layout, st)
            if (st /= 0) error stop 'geometry update'
            call rebuild(.false.)
        end do
    end subroutine current_mode

    subroutine reknot(law, n)
        !! The law on n knots in psi_N packed towards the axis and the edge,
        !! s_j = (1 - cos(pi j/(n - 1)))/2, from its own values: the profile
        !! updates then sample the target at n levels.
        type(gs_profile_law_t), intent(inout) :: law
        integer, intent(in) :: n
        real(dp), allocatable :: knots(:), pv(:), fv(:)
        real(dp) :: pp, ff
        integer :: j
        knots = [(0.5_dp*(1 - cos(acos(-1.0_dp)*real(j, dp)/(n - 1))), j=0, n - 1)]
        allocate (pv(n), fv(n))
        do j = 1, n
            call gs_law_sources(law, law%axis*(1 - knots(j)), pp, ff)
            pv(j) = pp
            fv(j) = ff
        end do
        deallocate (law%knots, law%pp, law%ff)
        law%knots = knots
        allocate (law%pp(4, n - 1), law%ff(4, n - 1))
        call hermite_law(law%knots, pv, law%pp)
        call hermite_law(law%knots, fv, law%ff)
    end subroutine reknot

    subroutine anderson_move(gamma, contour, step, used)
        !! Anderson acceleration (memory 4) of the fixed point Gamma_p =
        !! contour(Gamma_p), damped by the relaxation; step 0 resets.
        real(dp), intent(inout) :: gamma(:, :)
        real(dp), intent(in) :: contour(:, :)
        integer, intent(in) :: step
        integer, intent(inout) :: used
        integer, parameter :: memory = 4
        real(dp), allocatable, save :: dfm(:, :), dgm(:, :), f_prev(:), g_prev(:)
        real(dp), allocatable :: xv(:), gv(:), fv(:), coeff(:)
        integer :: n
        n = size(gamma)
        xv = reshape(gamma, [n])
        gv = reshape(contour, [n])
        fv = gv - xv
        if (.not. allocated(dfm)) allocate (dfm(n, memory), dgm(n, memory))
        if (size(dfm, 1) /= n) then
            deallocate (dfm, dgm)
            allocate (dfm(n, memory), dgm(n, memory))
        end if
        if (step == 0) used = 0
        if (step > 0 .and. allocated(f_prev)) call push_history(dfm, fv - f_prev, dgm, gv - g_prev, &
                                                                 used, memory)
        f_prev = fv
        g_prev = gv
        if (used > 0) then
            allocate (coeff(memory))
            call least_squares(dfm(:, :used), fv, coeff)
            gv = gv - matmul(dgm(:, :used), coeff(:used))
            xv = xv - matmul(dgm(:, :used) - dfm(:, :used), coeff(:used))
        end if
        gamma = reshape(xv + relaxation*(gv - xv), shape(gamma))
    end subroutine anderson_move

    function pack_state(law, psi_b, nk) result(v)
        !! Law values at the knots (p', then FF') and psi_b in one vector.
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi_b
        integer, intent(in) :: nk
        real(dp) :: v(2*nk + 1)
        v(1:nk - 1) = law%pp(1, :)
        v(nk) = hermite_end(law%pp, law%knots)
        v(nk + 1:2*nk - 1) = law%ff(1, :)
        v(2*nk) = hermite_end(law%ff, law%knots)
        v(2*nk + 1) = psi_b
    end function pack_state

    subroutine trace_contour(level, origin, points, status)
        !! The level psi = level along contour_samples rays from a fixed origin
        !! inside the plasma (fixed, so that a sample keeps its angle between
        !! moves), counter-clockwise from the outboard midplane. Newton along
        !! each ray with d psi/d rho = grad psi . e, safeguarded by a bracket,
        !! started from the radius of the previous call (first call: a march
        !! to the first sign change). The rays are independent (OpenMP).
        real(dp), intent(in) :: level, origin(2)
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, intent(out) :: status
        real(dp) :: theta, e(2), lo, hi, rho, f, df
        integer :: j, it, rst
        logical :: failed
        allocate (points(2, contour_samples))
        if (.not. allocated(ray_radius)) then
            allocate (ray_radius(contour_samples), ray_block(contour_samples))
            ray_radius = -1
            ray_block = 1
        end if
        failed = .false.
        !$omp parallel do default(shared) schedule(dynamic, 4) &
        !$omp private(j, it, theta, e, lo, hi, rho, f, df, rst) reduction(.or.:failed)
        do j = 1, contour_samples
            if (failed) cycle
            theta = 2*acos(-1.0_dp)*(j - 1)/contour_samples
            e = [cos(theta), sin(theta)]
            if (ray_radius(j) <= 0) then
                lo = 0
                hi = 0
                do
                    hi = hi + 0.02_dp
                    call ray_value(space, coef, origin, e, level, 1.0_dp, hi, ray_block(j), &
                                   f, df, rst)
                    if (rst /= 0 .or. f <= 0) exit
                    lo = hi
                end do
                rho = 0.5_dp*(lo + hi)
            else
                rho = ray_radius(j)
                lo = 0.5_dp*rho
                hi = 1.5_dp*rho
            end if
            rst = 0
            do it = 1, 50
                call ray_value(space, coef, origin, e, level, 1.0_dp, rho, ray_block(j), f, df, &
                               rst)
                if (rst /= 0) exit
                if (f > 0) then
                    lo = max(lo, rho)
                else
                    hi = min(hi, rho)
                end if
                if (abs(f) <= 1.0e-14_dp*max(1.0_dp, abs(level))) exit
                if (df < 0) then
                    rho = rho - f/df
                    if (rho <= lo .or. rho >= hi) rho = 0.5_dp*(lo + hi)
                else
                    rho = 0.5_dp*(lo + hi)
                end if
                if (hi - lo <= 1.0e-13_dp) exit
            end do
            if (rst /= 0) then
                failed = .true.
                cycle
            end if
            ray_radius(j) = rho
            points(:, j) = origin + rho*e
        end do
        !$omp end parallel do
        status = merge(1, 0, failed)
    end subroutine trace_contour

    subroutine read_curve(path, array)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: array(:, :)
        integer :: unit, n, j
        open (newunit=unit, file=path, status='old', action='read')
        read (unit, *) n
        allocate (array(2, n))
        do j = 1, n
            read (unit, *) array(:, j)
        end do
        close (unit)
    end subroutine read_curve

    subroutine read_seed(path)
        character(*), intent(in) :: path
        integer :: unit, nr, nz
        open (newunit=unit, file=path, status='old', action='read')
        read (unit, *) nr, nz
        allocate (seed_r(nr), seed_z(nz), seed_psi(nr, nz))
        read (unit, *) seed_r
        read (unit, *) seed_z
        read (unit, *) seed_psi
        close (unit)
    end subroutine read_seed

    real(dp) function seed_value(r, z) result(value)
        !! Cubic Lagrange interpolation of the seed map.
        real(dp), intent(in) :: r, z
        real(dp) :: wr(4), wz(4)
        integer :: ir, iz
        ir = min(max(count(seed_r < r) - 1, 1), size(seed_r) - 3)
        iz = min(max(count(seed_z < z) - 1, 1), size(seed_z) - 3)
        call lagrange_weights(4, r, seed_r(ir:ir + 3), wr)
        call lagrange_weights(4, z, seed_z(iz:iz + 3), wz)
        value = dot_product(wr, matmul(seed_psi(ir:ir + 3, iz:iz + 3), wz))
    end function seed_value
end program gs_spline_free_main
