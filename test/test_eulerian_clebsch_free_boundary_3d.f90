program test_eulerian_clebsch_free_boundary_3d
    !! Gate E6b of docs/eulerian-label-plan.md: a three-dimensional free
    !! boundary. The 3D label lane (E4: chi, g per toroidal plane, the axis
    !! curve X_k, the multipliers mu_k) inside a moving plasma boundary Gamma_p,
    !! a vacuum shell between Gamma_p and a fixed conducting wall Gamma_c, and
    !! the pressure balance on Gamma_p.
    !!
    !! Shell: B = B_s + grad Phi with B_s a given divergence-free source field
    !! and Phi single valued in V0 of the shell section (FortFEM O-grid with an
    !! annulus between Gamma_p and Gamma_c) times the trigonometric nodal basis
    !! in phi. Phi minimises int |B_s + grad Phi|^2/2 dV over the shell, so
    !! div B = 0 there and B.n = 0 on Gamma_p and Gamma_c are its natural
    !! boundary conditions (Phi is fixed up to a constant).
    !!
    !! Gamma_p: on plane k, R = c_k + r_k(theta) (cos theta, sin theta) about
    !! the exact axis c_k of that plane, r_k a Fourier series of 2M + 1
    !! coefficients. Rows: per plane and mode,
    !! oint (p + |B_in|^2/2) w_m dl = oint |B_s + grad Phi|^2/2 w_m dl, the
    !! plasma side on the plasma mesh of the plane, the shell side on the
    !! shell mesh, B in 3D (with the phi derivatives of the nodal basis).
    !!
    !! Exact manufactured case: Landreman's sheared family as in E4 (eps 0.6,
    !! S 2.4, delta 0.2), B_s the oracle field, Gamma_c the oracle surface
    !! chi = 1.5, p = (delta - psi)/lambda^2 (zero on Gamma_p). The exact
    !! solution has Gamma_p = {chi = 1}, the oracle labels and Phi = 0; B is
    !! continuous on Gamma_p, so the balance holds.
    !!
    !! Discretisation: a fixed reference geometry (the plane-1 meshes of the
    !! exact state) carries the quadrature, the polar rule at the axis vertex
    !! and the basis; every plane is mapped isoparametrically through its node
    !! coordinates (as clebsch_3d_system's plane_nodes), which come from the
    !! layouts of that plane: the plasma O-grid with wall Gamma_p,k and its core
    !! centred on X_k, the shell O-grid with interface Gamma_p,k and wall
    !! Gamma_c,k. The X and Gamma_p columns of the Newton matrix are central
    !! differences, each rebuilding one plane; the trace rows are
    !! differentiated by the chain rule.
    !!
    !! Arguments: cells (default 4), planes (default 9), Fourier order M
    !! (default 6; at most about a third of the 4n boundary cells, or the
    !! unresolved modes make the Newton matrix nearly singular), perturbation
    !! (default 0.01), 'jacobian'.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate, &
        spline_side_quadrature_t, spline_side_quadrature, spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_t_basis
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve, eu_constrained_t, &
        eu_constrained_factor, eu_constrained_apply, eu_constrained_free
    use milca_eulerian_clebsch, only: clebsch_3d_system, clebsch_3d_shape_system
    use milca_eulerian_clebsch_leaves, only: shell_energy_ref
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_gs_spline, only: gs_spline_project
    use milca_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_b_cart
    use fortfem_spline_eval, only: spline_boundary_nodes
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric
    use milca_nested_oracles, only: NESTED_SHEARED
    use solovev_label_fields, only: triplet_list
    use landreman_sections, only: section_oracle, section_label, section_exact_b, &
        section_sqrt_hessian, section_axis, section_boundary_radius, section_fit_profiles, &
        section_profiles
    implicit none
    integer, parameter :: degree = 2, nsample = 96
    real(dp), parameter :: pi = acos(-1.0_dp), wall_level = 1.5_dp
    real(dp), parameter :: ridge = 1e-12_dp
    integer :: n = 4, nt = 9, nmode = 6, max_steps = 14, failures = 0, nrho
    real(dp) :: perturbation = 0.01_dp
    logical :: jacobian_check = .false.
    character(16) :: arg
    real(dp) :: t_nodes = 0, t_clebsch = 0, t_shell = 0, t_balance = 0, t_unused = 0
    integer :: n_residual = 0, n_nodes = 0

    ! reference spaces (fixed geometry) and work spaces (per-plane geometry)
    type(spline_space_t), target :: ref_p, ref_f, work_p, work_f
    type(spline_derham_t) :: derham_ref
    type(spline_quadrature_t) :: aquad_ref, quad_fref
    ! exact axes c_k, frames, the Gamma_c samples per plane
    real(dp), allocatable :: centre(:, :), frames(:, :, :), wall_s(:, :, :)
    ! half width of the core block per plane, fixed from the exact boundary:
    ! sizing it from the extremes of the current boundary samples makes the
    ! nodes a kinked (max/min) function of Gamma_p and stalls Newton
    real(dp), allocatable :: core_half(:)
    ! per-plane node coordinates of the plasma and the shell meshes
    real(dp), allocatable :: nodes_p(:, :, :), nodes_f(:, :, :)
    ! the shell unknowns: full-space dofs touching the annulus
    integer, allocatable :: shell_dof(:), shell_slot(:)
    integer :: n0p, nphi, nfull
    ! precomputed basis at the reference points: plasma wall side, shell
    ! interface side (block 6), shell volume points (block 6), axis vertex
    type :: point_rows_t
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v0(:), v1(:), v2(:)
    end type point_rows_t
    type(spline_side_quadrature_t) :: side_p, side_f
    ! the volume point of every shell quadrature point (0: not in the annulus)
    integer, allocatable :: vol_index(:)
    ! analytic boundary and axis columns (shape matrices) or differences
    logical :: analytic_columns = .true.
    type(point_rows_t) :: rows_wall, rows_iface, rows_vol, rows_axis
    integer, allocatable :: vol_points(:), iface_points(:)
    real(dp), allocatable :: u_exact(:)
    integer :: ichi, ig, ix, imu, iphi, icirc, irho, nu
    ! wall model: 'ideal' (conducting Gamma_c, the circulation held at 0) or
    ! 'exterior' (axisymmetric Gamma_c, the exterior Steklov-Poincare map per
    ! toroidal mode, the plasma current as the circulation c of a filament)
    character(8) :: wall_model = 'ideal'
    real(dp) :: c_exact = 0, filament_circulation = 1, filament_scale = 1
    type(coil_set_t) :: filament
    ! exterior: wall trace nodes, their rows, the blocks Q(d) of the exterior
    ! energy between planes k and k + d, the wall side quadrature rows
    type(point_rows_t) :: rows_wnodes, rows_wside
    type(spline_side_quadrature_t) :: side_w
    real(dp), allocatable :: qblock(:, :, :)
    integer, allocatable :: wall_p(:)
    ! the exact null directions of the Newton matrix, removed by pins
    ! instead of ridges: Phi -> Phi + a by one Phi value
    integer :: gauge_pins(1)
    ! g -> g + k(chi) for any k is the gauge of the continuum; discretely
    ! k(chi_h) is nearly in the space of g, so every such k is a weak
    ! direction (singular values 1e-7 at 3 cells, more under refinement).
    ! Fixed by orthogonality to k = T_j(2 chi_0 - 1), j < ngauge, with the
    ! reference label chi_0 and the reference measure (linear, constant,
    ! no shape derivative), with multipliers nu_j
    integer :: ignu, ngauge
    real(dp) :: check_src(3, 2)
    real(dp), allocatable :: gauge_rows(:, :)
    ! the coordinate (1: R, 2: Z) and the plane of the projection callbacks
    integer :: current_component = 1, label_plane = 1

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) n
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) nmode
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) perturbation
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        jacobian_check = arg == 'jacobian'
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        wall_model = arg
    end if
    if (wall_model == 'shape') then
        jacobian_check = .true.
    else if (wall_model /= 'ideal' .and. wall_model /= 'exterior') then
        error stop 'wall: ideal or exterior (or shape: check the shape derivative)'
    end if
    nrho = 2*nmode + 1
    section_oracle%family = NESTED_SHEARED
    section_oracle%eps = 0.6_dp
    section_oracle%saxis = 2.4_dp
    section_oracle%delta = 0.2_dp
    if (.not. section_oracle%valid()) error stop 'oracle parameters'
    call section_fit_profiles()
    print '(a,i3,a,i3,a,i3,a,f6.3)', 'cells', n, '  planes', nt, '  Fourier order', nmode, &
        '  perturbation', perturbation
    call setup()
    if (jacobian_check .and. wall_model == 'shape') then
        call shape_check()
        if (failures > 0) error stop 'shape derivative check failed'
        stop
    end if
    call run()
    if (failures > 0) error stop 'free-boundary 3D gate E6b failed'

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

    ! ------------------------------------------------------------------ curves

    pure real(dp) function mode(k, theta)
        !! k = 1: 1; k = 2m: cos(m theta); k = 2m + 1: sin(m theta).
        integer, intent(in) :: k
        real(dp), intent(in) :: theta

        if (k == 1) then
            mode = 1
        else if (mod(k, 2) == 0) then
            mode = cos((k/2)*theta)
        else
            mode = sin((k/2)*theta)
        end if
    end function mode

    pure real(dp) function plane_phi(k)
        integer, intent(in) :: k

        plane_phi = 2*pi*(k - 1)/nt
    end function plane_phi

    function exact_rho(k) result(rho)
        !! Fourier coefficients of the exact Gamma_p of plane k about c_k.
        integer, intent(in) :: k
        real(dp) :: rho(nrho), theta, r
        integer, parameter :: nq = 256
        integer :: j, m

        rho = 0
        do j = 1, nq
            theta = 2*pi*(j - 1)/nq
            r = section_boundary_radius(plane_phi(k), centre(:, k), theta)
            do m = 1, nrho
                rho(m) = rho(m) + r*mode(m, theta)*merge(1.0_dp, 2.0_dp, m == 1)/nq
            end do
        end do
    end function exact_rho

    real(dp) function gamma_error(k, rho) result(e)
        !! The largest radial error of Gamma_p on plane k.
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:)
        real(dp) :: theta, r
        integer :: j, m

        e = 0
        do j = 1, 180
            theta = 2*pi*(j - 1)/180
            r = 0
            do m = 1, nrho
                r = r + rho(m)*mode(m, theta)
            end do
            e = max(e, abs(r - section_boundary_radius(plane_phi(k), centre(:, k), theta)))
        end do
    end function gamma_error

    function curve(k, rho) result(s)
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:)
        real(dp) :: s(2, nsample), theta, r
        integer :: i, m

        do i = 1, nsample
            theta = 2*pi*(i - 1)/nsample
            r = 0
            do m = 1, nrho
                r = r + rho(m)*mode(m, theta)
            end do
            s(:, i) = centre(:, k) + r*[cos(theta), sin(theta)]
        end do
    end function curve

    ! ------------------------------------------------------------------ meshes

    subroutine plasma_layout(k, rho, x, layout)
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:), x(2)
        type(block_layout_t), intent(out) :: layout
        real(dp) :: s(2, nsample)
        integer :: st

        s = curve(k, rho)
        call layout_ogrid(centre(:, k), core_half(k), 0.5_dp, layout, st)
        if (st /= 0) error stop 'plasma layout'
        call layout_set_curve(layout, LAYOUT_WALL, s, st, centre_core=.true.)
        if (st /= 0) error stop 'plasma boundary curve'
        layout%vertices(:, 1:4) = layout%vertices(:, 1:4) + spread(x - centre(:, k), 2, 4)
    end subroutine plasma_layout

    subroutine shell_layout(k, rho, layout)
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:)
        type(block_layout_t), intent(out) :: layout
        integer :: st

        call layout_ogrid(centre(:, k), 0.5_dp*(maxval(wall_s(1, :, k)) - minval(wall_s(1, :, k))), &
                          0.3_dp, layout, st, interface_fraction=0.6_dp)
        if (st /= 0) error stop 'shell layout'
        call layout_set_curve(layout, LAYOUT_WALL, wall_s(:, :, k), st)
        if (st /= 0) error stop 'shell wall curve'
        call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, curve(k, rho), st)
        if (st /= 0) error stop 'shell interface curve'
    end subroutine shell_layout

    subroutine plane_nodes(k, rho, x)
        !! The node coordinates of plane k for its boundary rho and axis x: the
        !! control points of the geometry of the plane's spaces (exactly the
        !! isoparametric V0 values; no projection).
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:), x(2)
        type(block_layout_t) :: layout
        real(dp) :: real_start
        integer :: st

        real_start = tick0()
        n_nodes = n_nodes + 1
        call plasma_layout(k, rho, x, layout)
        call spline_space_update_geometry(work_p, layout, st)
        if (st /= 0) error stop 'plasma plane geometry'
        call control_points(work_p, nodes_p(:, :, k))
        call shell_layout(k, rho, layout)
        call spline_space_update_geometry(work_f, layout, st)
        if (st /= 0) error stop 'shell plane geometry'
        call control_points(work_f, nodes_f(:, :, k))
        t_nodes = t_nodes + tick0()
    end subroutine plane_nodes

    subroutine control_points(space, nodes)
        !! Geometry control points per global function (C0: every broken copy
        !! of a function carries the same point).
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: nodes(:, :)
        integer :: j

        do j = 1, space%broken_count
            nodes(:, space%broken_dof(space%broken_ptr(j))) = space%geometry(:, j)
        end do
    end subroutine control_points

    real(dp) function tick0()
        !! Seconds since the previous call.
        integer(8), save :: last = 0
        integer(8) :: now, rate

        call system_clock(now, rate)
        tick0 = 0
        if (last > 0) tick0 = real(now - last, dp)/rate
        last = now
    end function tick0

    subroutine coordinate_p(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(work_p, block, xi, jet, st)
        values(1) = jet(current_component, 1)
        associate (unused => t)
        end associate
    end subroutine coordinate_p

    subroutine label_ref(block, xi, t, values)
        !! The exact label on plane label_plane at the point of the plane mesh.
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(work_p, block, xi, jet, st)
        values(1) = section_label(jet(1, 1), plane_phi(label_plane), jet(2, 1))
        associate (unused => t)
        end associate
    end subroutine label_ref

    subroutine setup()
        !! Exact axes, frames and walls; reference and work spaces; the
        !! precomputed point rows; the unknown layout and the exact state.
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad
        real(dp) :: theta, xi(2), sk(2, nsample)
        real(dp), allocatable :: rho1(:), v(:)
        logical, allocatable :: in6(:)
        integer :: k, i, j, st, block

        allocate (centre(2, nt), frames(2, 2, nt), wall_s(2, nsample, nt), core_half(nt))
        do k = 1, nt
            centre(:, k) = section_axis(plane_phi(k))
            frames(:, :, k) = section_sqrt_hessian(plane_phi(k), centre(:, k))
            do i = 1, nsample
                theta = 2*pi*(i - 1)/nsample
                wall_s(:, i, k) = centre(:, k) + section_boundary_radius(plane_phi(k), &
                    centre(:, k), theta, wall_level)*[cos(theta), sin(theta)]
            end do
        end do
        do k = 1, nt
            rho1 = exact_rho(k)
            sk = curve(k, rho1)
            core_half(k) = 0.5_dp*(maxval(sk(1, :)) - minval(sk(1, :)))
        end do
        if (wall_model == 'exterior') call axisymmetric_wall()
        print '(a,2f8.4,a,2f8.4)', 'exact axis: R from', minval(centre(1, :)), &
            maxval(centre(1, :)), '  Z from', minval(centre(2, :)), maxval(centre(2, :))
        ! reference plasma space: plane 1 exact
        rho1 = exact_rho(1)
        call plasma_layout(1, rho1, centre(:, 1), layout)
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, ref_p, st)
        if (st /= 0) error stop 'reference plasma space'
        work_p = ref_p
        call spline_derham_build(ref_p, derham_ref, st, 0)
        if (st /= 0) error stop 'reference derham'
        n0p = derham_ref%dimension(0)
        call spline_quadrature_build(ref_p, degree + 8, quad, st)
        if (st /= 0) error stop 'reference quadrature'
        call spline_quadrature_point_singular(ref_p, quad, centre(:, 1), aquad_ref, st)
        if (st /= 0) error stop 'reference polar rule'
        call spline_derham_boundary_dofs(derham_ref, LAYOUT_WALL, 0, wall_p)
        ! reference shell space: plane 1 exact
        call shell_layout(1, rho1, layout)
        deallocate (options%knots)
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        options%knots(1, 6)%element_count = 4*n
        call spline_space_build(layout, options, ref_f, st)
        if (st /= 0) error stop 'reference shell space'
        work_f = ref_f
        nfull = ref_f%dof_count
        call spline_quadrature_build(ref_f, degree + 4, quad_fref, st)
        if (st /= 0) error stop 'reference shell quadrature'
        allocate (in6(nfull))
        in6 = .false.
        associate (blk => ref_f%blocks(6))
            do j = blk%offset + 1, blk%offset + blk%n(1)*blk%n(2)
                do i = ref_f%broken_ptr(j), ref_f%broken_ptr(j + 1) - 1
                    if (ref_f%broken_weight(i) /= 0) in6(ref_f%broken_dof(i)) = .true.
                end do
            end do
        end associate
        shell_dof = pack([(i, i=1, nfull)], in6)
        nphi = size(shell_dof)
        allocate (shell_slot(nfull))
        shell_slot = 0
        shell_slot(shell_dof) = [(i, i=1, nphi)]
        ! precomputed rows
        call spline_side_quadrature(ref_p, LAYOUT_WALL, degree + 3, side_p, st)
        if (st /= 0) error stop 'wall side quadrature'
        call point_rows(ref_p, side_p%block, side_p%xi, rows_wall)
        call spline_side_quadrature(ref_f, LAYOUT_ANNULUS_INTERFACE, degree + 3, side_f, st)
        if (st /= 0) error stop 'interface side quadrature'
        iface_points = pack([(i, i=1, side_f%point_count)], side_f%block == 6)
        call point_rows(ref_f, side_f%block(iface_points), side_f%xi(:, iface_points), &
                        rows_iface)
        vol_points = pack([(i, i=1, quad_fref%point_count)], quad_fref%block == 6)
        allocate (vol_index(quad_fref%point_count))
        vol_index = 0
        vol_index(vol_points) = [(i, i=1, size(vol_points))]
        call point_rows(ref_f, quad_fref%block(vol_points), quad_fref%xi(:, vol_points), &
                        rows_vol)
        if (wall_model == 'exterior') call exterior_setup()
        call filament_setup()
        ! the axis vertex: values and the gradients of one adjacent cell
        call spline_locate(ref_p, 1, centre(:, 1) + [1e-9_dp, 1e-9_dp], block, xi, st)
        if (st /= 0) error stop 'axis vertex'
        call point_rows(ref_p, [block], reshape(xi, [2, 1]), rows_axis)
        ! unknown layout
        ichi = 0
        ig = ichi + n0p*nt
        ix = ig + n0p*nt
        imu = ix + 2*nt
        ! gauge multipliers after the axis multipliers (counted in that block)
        ignu = imu + nt
        ngauge = degree*n + 1
        iphi = ignu + ngauge
        icirc = iphi + nphi*nt + 1
        irho = icirc
        gauge_pins(1) = iphi + 1
        nu = irho + nrho*nt
        allocate (nodes_p(2, n0p, nt), nodes_f(2, nfull, nt))
        print '(a,i6,a,i6,a,i6,a,i6)', 'unknowns', nu, '  (labels per plane', n0p, &
            ', shell potential per plane', nphi, ', boundary', nrho*nt, ')'
        ! the exact state (labels projected on the exact plane meshes)
        allocate (u_exact(nu))
        u_exact = 0
        u_exact(icirc) = c_exact
        allocate (v(n0p))
        do k = 1, nt
            rho1 = exact_rho(k)
            u_exact(irho + (k - 1)*nrho + 1:irho + k*nrho) = rho1
            u_exact(ix + 2*k - 1:ix + 2*k) = centre(:, k)
            call plane_nodes(k, rho1, centre(:, k))
            label_plane = k
            call plasma_layout(k, rho1, centre(:, k), layout)
            call spline_space_update_geometry(work_p, layout, st)
            block
                type(spline_derham_t) :: dw
                call spline_derham_build(work_p, dw, st, 0)
                call spline_derham_project(dw, 0, label_ref, v, st, 12)
            end block
            v(wall_p) = 1
            u_exact(ichi + (k - 1)*n0p + 1:ichi + k*n0p) = v
        end do
        call gauge_setup(quad)
    end subroutine setup

    subroutine gauge_setup(quad)
        !! gauge_rows(j, (k-1) n0p + i) = sum_q w_q N_i(q) T_j(2 chi_0,k(q) - 1)/nt
        !! on the reference plane, rows scaled to unit largest entry.
        type(spline_quadrature_t), intent(in) :: quad
        type(point_rows_t) :: rq
        real(dp) :: chi0, x
        integer :: q, k, j, a, b

        call point_rows(ref_p, quad%block, quad%xi, rq)
        allocate (gauge_rows(ngauge, n0p*nt))
        gauge_rows = 0
        do k = 1, nt
            do q = 1, quad%point_count
                a = rq%rp(q)
                b = rq%rp(q + 1) - 1
                chi0 = dot_product(rq%v0(a:b), u_exact(ichi + (k - 1)*n0p + rq%ci(a:b)))
                x = max(-1.0_dp, min(1.0_dp, 2*chi0 - 1))
                do j = 1, ngauge
                    gauge_rows(j, (k - 1)*n0p + rq%ci(a:b)) = &
                        gauge_rows(j, (k - 1)*n0p + rq%ci(a:b)) + &
                        quad%weight(q)*rq%v0(a:b)*cos((j - 1)*acos(x))/nt
                end do
            end do
        end do
        do j = 1, ngauge
            gauge_rows(j, :) = gauge_rows(j, :)/maxval(abs(gauge_rows(j, :)))
        end do
    end subroutine gauge_setup

    subroutine gauge_kkt(u, r, jt)
        !! The gauge rows and their multiplier terms (symmetric, constant).
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        integer :: j, i

        do j = 1, ngauge
            r(ignu + j) = dot_product(gauge_rows(j, :), u(ig + 1:ix))
            r(ig + 1:ix) = r(ig + 1:ix) + gauge_rows(j, :)*u(ignu + j)
            if (present(jt)) then
                do i = 1, n0p*nt
                    if (gauge_rows(j, i) == 0) cycle
                    call jt%add(ig + i, ignu + j, gauge_rows(j, i))
                    call jt%add(ignu + j, ig + i, gauge_rows(j, i))
                end do
                call jt%add(ignu + j, ignu + j, 0.0_dp)
            end if
        end do
    end subroutine gauge_kkt

    subroutine axisymmetric_wall()
        !! Gamma_c for the exterior map: the ellipse about (R_m, 0) through
        !! 1.15 times the extent of the chi = wall_level sections of every
        !! plane, the same curve on every plane.
        real(dp) :: rmin, rmax, zmax, rm, a, b, theta
        integer :: i, k

        rmin = minval(wall_s(1, :, :))
        rmax = maxval(wall_s(1, :, :))
        zmax = maxval(abs(wall_s(2, :, :)))
        rm = (rmin + rmax)/2
        a = 1.15_dp*(rmax - rmin)/2
        b = 1.15_dp*zmax
        do k = 1, nt
            do i = 1, nsample
                theta = 2*pi*(i - 1)/nsample
                wall_s(:, i, k) = [rm + a*cos(theta), b*sin(theta)]
            end do
        end do
        print '(a,3f8.4)', 'axisymmetric Gamma_c: centre R, semi axes', rm, a, b
    end subroutine axisymmetric_wall

    subroutine exterior_check(wn)
        !! The exterior energy of the wall traces against the exact one for
        !! the decaying harmonic Phi of a ring source with density cos(n phi)
        !! inside Gamma_c, one toroidal mode n <= (nt - 1)/2 at a time:
        !! 1/2 int_ext |grad Phi|^2 = -1/2 oint Phi dPhi/dn dS (n out of the
        !! torus), by a fine tensor rule on the axisymmetric wall.
        real(dp), intent(in) :: wn(:, :)
        real(dp) :: tr(size(wn, 2), nt), eh, ex, x(3), th, ph, rz(2), drz(2), &
                    nrm(3), dth, dph, phi, gphi(3), a, b, rm
        integer :: nw, iq, k, j, d, i, l, set
        integer, parameter :: nth = 512, nph = 512

        nw = size(wn, 2)
        rm = (minval(wall_s(1, :, 1)) + maxval(wall_s(1, :, 1)))/2
        a = (maxval(wall_s(1, :, 1)) - minval(wall_s(1, :, 1)))/2
        b = maxval(wall_s(2, :, 1))
        check_src(1, 1) = rm
        do set = 1, (nt + 1)/2
            do k = 1, nt
                ph = 2*pi*(k - 1)/nt
                do iq = 1, nw
                    x = [wn(1, iq)*cos(ph), wn(1, iq)*sin(ph), wn(2, iq)]
                    call check_sources(x, set, phi, gphi)
                    tr(iq, k) = phi
                end do
            end do
            eh = 0
            do k = 1, nt
                do j = 1, nt
                    d = modulo(k - j, nt)
                    eh = eh + 0.5_dp*dot_product(tr(:, k), matmul(qblock(:, :, d), tr(:, j)))
                end do
            end do
            ex = 0
            dth = 2*pi/nth
            dph = 2*pi/nph
            do i = 1, nth
                th = (i - 0.5_dp)*dth
                rz = [rm + a*cos(th), b*sin(th)]
                drz = [-a*sin(th), b*cos(th)]
                do l = 1, nph
                    ph = (l - 0.5_dp)*dph
                    x = [rz(1)*cos(ph), rz(1)*sin(ph), rz(2)]
                    ! outward normal times |d(R,Z)/dtheta|: (dZ, -dR) rotated to phi
                    nrm = [drz(2)*cos(ph), drz(2)*sin(ph), -drz(1)]
                    call check_sources(x, set, phi, gphi)
                    ex = ex - 0.5_dp*phi*dot_product(gphi, nrm)*rz(1)*dth*dph
                end do
            end do
            print '(a,i2,a,es14.6,a,es14.6,a,es10.2)', 'exterior map check, toroidal mode', set - 1, &
                ': discrete', eh, '  exact', ex, '  relative error', abs(eh - ex)/abs(ex)
            call check(abs(eh - ex) < 1e-2_dp*abs(ex), &
                       'exterior energy of a single-mode ring source within 1 %')
        end do


    end subroutine exterior_check

    subroutine wall_spectrum()
        !! Toroidal spectrum of the source data (B_oracle).n on Gamma_c:
        !! the L2 norm over theta of each mode n, against the planes' n_max.
        integer, parameter :: nth = 128, nph = 128
        real(dp) :: bn(nth, nph), rz(2), drz(2), b(3), a, bb, rm, th, ph, cn(0:nph/2), &
                    re, im
        integer :: i, l, n

        rm = (minval(wall_s(1, :, 1)) + maxval(wall_s(1, :, 1)))/2
        a = (maxval(wall_s(1, :, 1)) - minval(wall_s(1, :, 1)))/2
        bb = maxval(wall_s(2, :, 1))
        do i = 1, nth
            th = 2*pi*(i - 1)/nth
            rz = [rm + a*cos(th), bb*sin(th)]
            drz = [-a*sin(th), bb*cos(th)]
            do l = 1, nph
                ph = 2*pi*(l - 1)/nph
                call section_exact_b([rz(1), ph, rz(2)], b)
                bn(i, l) = (b(1)*drz(2) - b(3)*drz(1))/norm2(drz)
            end do
        end do
        cn = 0
        do n = 0, nph/2
            do i = 1, nth
                re = sum(bn(i, :)*cos(n*2*pi*[(l - 1, l=1, nph)]/nph))/nph
                im = sum(bn(i, :)*sin(n*2*pi*[(l - 1, l=1, nph)]/nph))/nph
                cn(n) = cn(n) + (re**2 + im**2)/nth
            end do
            cn(n) = sqrt(cn(n))
        end do
        print '(a,i3)', 'B.n on Gamma_c, toroidal spectrum (planes resolve n <=', (nt - 1)/2
        print '(a,12es9.1)', '   n = 0..11:', cn(0:11)
    end subroutine wall_spectrum

    subroutine check_sources(x, set, phi, gphi)
        !! Phi of a ring of charge at (R_m, 0) with density cos(n phi'),
        !! n = set - 1 (one toroidal mode, 1024 point charges), and its
        !! gradient.
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: set
        real(dp), intent(out) :: phi, gphi(3)
        real(dp) :: dx(3), r, p
        integer :: l
        integer, parameter :: nl = 1024

        phi = 0
        gphi = 0
        do l = 1, nl
            p = 2*pi*(l - 1)/nl
            dx = x - [check_src(1, 1)*cos(p), check_src(1, 1)*sin(p), 0.0_dp]
            r = norm2(dx)
            phi = phi + cos((set - 1)*p)/r*2*pi/nl
            gphi = gphi - cos((set - 1)*p)*dx/r**3*2*pi/nl
        end do
    end subroutine check_sources

    subroutine exterior_setup()
        !! The exterior energy 1/2 int |grad Phi_ext|^2 of the decaying
        !! continuation of the single-valued Phi beyond the axisymmetric
        !! Gamma_c: with the traces u_k = T Phi_k on the wall nodes of plane k,
        !! 1/2 sum_kj u_k^T Q(k - j) u_j, Q(d) = 2 pi/nt^2 (S_0 + 2 sum_n S_n
        !! cos(2 pi n d/nt)), S_n the Steklov-Poincare map of toroidal mode n
        !! (FortFEM, <S_n u, u> = int R |grad u|^2 + n^2 u^2/R over the
        !! exterior of Gamma_c). The wall side quadrature carries the source
        !! term -oint Phi B_src.n of the exterior field of the sources.
        type(gs_bem_curve_t) :: bcurve
        type(gs_bem_operators_t) :: ops
        integer, allocatable :: blk(:)
        real(dp), allocatable :: xi(:, :), wn(:, :), smat(:, :)
        integer :: st, mn, d, nw

        call spline_boundary_nodes(ref_f, LAYOUT_WALL, 2, blk, xi, wn, st)
        if (st /= 0) error stop 'wall nodes'
        call point_rows(ref_f, blk, xi, rows_wnodes)
        call gs_bem_curve_from_nodes(wn, bcurve, st)
        if (st /= 0) error stop 'wall BEM curve'
        nw = size(wn, 2)
        allocate (qblock(nw, nw, 0:nt - 1))
        qblock = 0
        do mn = 0, (nt - 1)/2
            call tm_bem_assemble(bcurve, mn, 1, ops, st)
            if (st /= 0) error stop 'toroidal-mode BEM'
            call tm_bem_steklov_symmetric(ops, smat, st)
            if (st /= 0) error stop 'Steklov-Poincare map'
            do d = 0, nt - 1
                qblock(:, :, d) = qblock(:, :, d) + merge(1.0_dp, 2.0_dp, mn == 0)* &
                                  cos(2*pi*mn*d/real(nt, dp))*smat
            end do
        end do
        qblock = qblock*2*pi/real(nt, dp)**2
        call exterior_check(wn)
        call wall_spectrum()
        call spline_side_quadrature(ref_f, LAYOUT_WALL, degree + 3, side_w, st)
        if (st /= 0) error stop 'wall side quadrature'
        call point_rows(ref_f, side_w%block, side_w%xi, rows_wside)
        print '(a,i4,a,i3)', 'exterior map: wall nodes', nw, '  toroidal modes', (nt - 1)/2
    end subroutine exterior_setup

    subroutine filament_setup()
        !! A unit-current filament along the axis curve (trigonometric
        !! interpolation of the plane axes); its circulation about the plasma
        !! of plane 1, and in the exterior model the oracle's, so that
        !! B_s = B_oracle - c_exact B_fil and the exact c is c_exact.
        integer, parameter :: nf = 96
        real(dp) :: pts(3, nf), lt(nt, 2), ht(nt, 2), x(2), t, sb(2, 512), bo(3), bf(3), &
                    tng(2), circ_o
        integer :: j, st

        do j = 1, nf
            t = (j - 1)/real(nf, dp)
            call spline_derham_t_basis(nt, t, lt, ht)
            x = matmul(centre, lt(:, 1))
            pts(:, j) = [x(1)*cos(2*pi*t), x(1)*sin(2*pi*t), x(2)]
        end do
        call coil_set_from_polylines(filament, [nf], pts, status=st)
        if (st /= 0) error stop 'filament'
        ! circulations about the exact plane-1 boundary (counter-clockwise)
        filament_circulation = 0
        circ_o = 0
        do j = 1, 512
            sb(:, j) = centre(:, 1) + section_boundary_radius(plane_phi(1), centre(:, 1), &
                2*pi*(j - 1)/512)*[cos(2*pi*(j - 1)/512), sin(2*pi*(j - 1)/512)]
        end do
        do j = 1, 512
            tng = (sb(:, modulo(j, 512) + 1) - sb(:, modulo(j - 2, 512) + 1))/2
            call filament_b([sb(1, j), plane_phi(1), sb(2, j)], bf)
            call section_exact_b([sb(1, j), plane_phi(1), sb(2, j)], bo)
            filament_circulation = filament_circulation + bf(1)*tng(1) + bf(3)*tng(2)
            circ_o = circ_o + bo(1)*tng(1) + bo(3)*tng(2)
        end do
        c_exact = 0
        if (wall_model == 'exterior') then
            ! the filament current in units of the plasma's, so that the
            ! unknown c is O(1) like the others (the Newton matrix had a
            ! singular value 1e-9 of the largest from units alone)
            filament_scale = circ_o/filament_circulation
            filament_circulation = circ_o
            c_exact = 1
        end if
        print '(a,es12.4,a,es12.4)', 'filament circulation', filament_circulation, &
            '  exact c', c_exact
    end subroutine filament_setup

    subroutine filament_b(x, b)
        !! Cylindrical components of the unit-filament field at x = (R, phi, Z).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: bc(3)

        call coil_set_b_cart(filament, [filament_scale], [x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)], bc)
        b = [bc(1)*cos(x(2)) + bc(2)*sin(x(2)), -bc(1)*sin(x(2)) + bc(2)*cos(x(2)), bc(3)]
    end subroutine filament_b

    subroutine source_b(x, b)
        !! The source field B_s = B_oracle - c_exact B_fil.
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: bf(3)

        call section_exact_b(x, b)
        if (c_exact /= 0) then
            call filament_b(x, bf)
            b = b - c_exact*bf
        end if
    end subroutine source_b

    subroutine point_rows(space, block, xi, rows)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block(:)
        real(dp), intent(in) :: xi(:, :)
        type(point_rows_t), intent(out) :: rows
        integer, allocatable :: rp(:), ci(:)
        integer :: st

        call spline_point_matrix(space, block, xi, 0, rows%rp, rows%ci, rows%v0, st)
        if (st /= 0) error stop 'point rows'
        call spline_point_matrix(space, block, xi, 1, rp, ci, rows%v1, st)
        if (any(ci /= rows%ci)) error stop 'point rows: pattern'
        call spline_point_matrix(space, block, xi, 2, rp, ci, rows%v2, st)
        if (any(ci /= rows%ci)) error stop 'point rows: pattern'
    end subroutine point_rows

    ! ------------------------------------------------------- local 3D jets

    subroutine plane_jet(rows, q, nodes, k, coef, ndof, slot, value, grad3, dval, dgrad, s, f, &
                         vel, lt)
        !! At reference point q (rows) on plane k: the physical point s, the
        !! map F = ds/ds0, the mesh velocity vel = sum_j L_j'(t_k) s_j, the
        !! field value and its gradient (R, phi, Z) of the nodal field coef
        !! (ndof per plane, slot maps the row columns to dofs, 0: not a dof),
        !! and the derivatives of value and gradient with respect to the dofs
        !! (dval(local), dgrad(3, local), local = (plane j, column)).
        type(point_rows_t), intent(in) :: rows
        integer, intent(in) :: q, k, ndof
        real(dp), intent(in) :: nodes(:, :, :), coef(:)
        integer, intent(in) :: slot(:)
        real(dp), intent(out) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2)
        real(dp), allocatable, intent(out) :: dval(:, :), dgrad(:, :, :)
        real(dp) :: ht(nt, 2), ft(2, 2), detf, g0(2), gx(2), dt_value, vals(nt)
        integer :: a, b, i, j, col, nloc

        call spline_derham_t_basis(nt, (k - 1)/real(nt, dp), lt, ht)
        a = rows%rp(q)
        b = rows%rp(q + 1) - 1
        nloc = b - a + 1
        s = matmul(nodes(:, rows%ci(a:b), k), rows%v0(a:b))
        f(:, 1) = matmul(nodes(:, rows%ci(a:b), k), rows%v1(a:b))
        f(:, 2) = matmul(nodes(:, rows%ci(a:b), k), rows%v2(a:b))
        vel = 0
        do j = 1, nt
            vel = vel + lt(j, 2)*matmul(nodes(:, rows%ci(a:b), j), rows%v0(a:b))
        end do
        detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
        ! inverse transpose of F
        ft = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/detf
        value = 0
        g0 = 0
        vals = 0
        do i = a, b
            col = slot(rows%ci(i))
            if (col == 0) cycle
            value = value + rows%v0(i)*coef((k - 1)*ndof + col)
            g0 = g0 + [rows%v1(i), rows%v2(i)]*coef((k - 1)*ndof + col)
            do j = 1, nt
                vals(j) = vals(j) + rows%v0(i)*coef((j - 1)*ndof + col)
            end do
        end do
        gx = matmul(ft, g0)
        dt_value = dot_product(lt(:, 2), vals)
        grad3 = [gx(1), (dt_value - dot_product(vel, gx))/(2*pi*s(1)), gx(2)]
        allocate (dval(nt, nloc), dgrad(3, nt, nloc))
        dval = 0
        dgrad = 0
        do i = a, b
            gx = matmul(ft, [rows%v1(i), rows%v2(i)])
            dval(k, i - a + 1) = rows%v0(i)
            dgrad(1, k, i - a + 1) = gx(1)
            dgrad(3, k, i - a + 1) = gx(2)
            do j = 1, nt
                dgrad(2, j, i - a + 1) = rows%v0(i)*lt(j, 2)/(2*pi*s(1))
            end do
            dgrad(2, k, i - a + 1) = dgrad(2, k, i - a + 1) - dot_product(vel, gx)/(2*pi*s(1))
        end do
    end subroutine plane_jet

    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

    ! --------------------------------------------------------------- residual

    subroutine residual(u, r, jt, ferr)
        !! Residual at u; with jt the transposed Newton matrix without the X
        !! and Gamma_p columns; with ferr the field error of the plasma.
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        real(dp), intent(out), optional :: ferr(2)
        real(dp), allocatable :: grad(:)
        real(dp) :: energy, xk(2, nt)
        type(csr_t) :: hessian
        integer :: k, i, p

        allocate (r(nu))
        r = 0
        n_residual = n_residual + 1
        t_unused = tick0()
        xk = reshape(u(ix + 1:ix + 2*nt), [2, nt])
        if (present(ferr)) then
            call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                                   u(ichi + 1:ig), u(ig + 1:ix), energy, grad, &
                                   plane_axis=xk, plane_nodes=nodes_p, plane_frames=frames, &
                                   exact_field=section_exact_b, field_error=ferr)
        else if (present(jt)) then
            call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                                   u(ichi + 1:ig), u(ig + 1:ix), energy, grad, hessian, xk, &
                                   nodes_p, frames)
        else
            call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                                   u(ichi + 1:ig), u(ig + 1:ix), energy, grad, &
                                   plane_axis=xk, plane_nodes=nodes_p, plane_frames=frames)
        end if
        r(1:ix) = grad
        t_clebsch = t_clebsch + tick0()
        call axis_rows(u, r, jt)
        call gauge_kkt(u, r, jt)
        call shell_rows(u, r, jt)
        t_shell = t_shell + tick0()
        call balance_rows(u, r, jt)
        t_balance = t_balance + tick0()
        do k = 1, nt
            r(ichi + (k - 1)*n0p + wall_p) = 0
        end do
        if (.not. present(jt)) return
        do i = 1, hessian%nrow
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call jt%add(hessian%col(p), i, hessian%val(p))
            end do
        end do
    end subroutine residual

    subroutine axis_rows(u, r, jt)
        !! grad chi_k(X_k) = 0, chi_k(X_k) = 0 and the multiplier terms; X_k
        !! is the central vertex of plane k, so N(X_k) is fixed and grad N(X_k)
        !! = F_k^-T grad_0 N.
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        real(dp) :: f(2, 2), ft(2, 2), detf, gx(2)
        integer :: k, a, b, i, col

        a = rows_axis%rp(1)
        b = rows_axis%rp(2) - 1
        do k = 1, nt
            f(:, 1) = matmul(nodes_p(:, rows_axis%ci(a:b), k), rows_axis%v1(a:b))
            f(:, 2) = matmul(nodes_p(:, rows_axis%ci(a:b), k), rows_axis%v2(a:b))
            detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
            ft = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/detf
            do i = a, b
                col = ichi + (k - 1)*n0p + rows_axis%ci(i)
                gx = matmul(ft, [rows_axis%v1(i), rows_axis%v2(i)])
                r(ix + 2*k - 1:ix + 2*k) = r(ix + 2*k - 1:ix + 2*k) + gx*u(col)
                r(imu + k) = r(imu + k) + rows_axis%v0(i)*u(col)
                r(col) = r(col) + u(imu + k)*rows_axis%v0(i)
                if (present(jt)) then
                    call jt%add(col, ix + 2*k - 1, gx(1))
                    call jt%add(col, ix + 2*k, gx(2))
                    call jt%add(col, imu + k, rows_axis%v0(i))
                    call jt%add(imu + k, col, rows_axis%v0(i))
                end if
            end do
            if (present(jt)) call jt%add(imu + k, imu + k, 0.0_dp)
        end do
    end subroutine axis_rows

    subroutine shell_rows(u, r, jt)
        !! d/dPhi of int |B_s + grad Phi|^2/2 dV over the shell, with the
        !! Hessian when jt is present; the constant is held by a ridge.
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        real(dp), allocatable :: dval(:, :), dgrad(:, :, :)
        real(dp) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2), bs(3), b(3), w, &
                    detf, tq, kmax, bf(3), c
        real(dp), allocatable :: e(:, :)
        integer :: q, m, k, j, i, a, nloc, rowi, rowj, jj
        integer, parameter :: per_plane_points = 2
        integer, allocatable :: cols(:)

        kmax = 0
        c = u(icirc)
        if (.not. present(jt) .or. analytic_columns) then
            ! the residual threaded over points; with analytic columns the
            ! Hessian and the c column come from shell_shape (one pass)
            call shell_residual(u, r)
            if (wall_model == 'exterior') call exterior_rows(u, r, jt)
            return
        end if
        do q = 1, size(vol_points)
            a = rows_vol%rp(q)
            nloc = rows_vol%rp(q + 1) - a
            do k = 1, nt
                ! the toroidal quadrature: the planes and the mid-points
                do m = 1, per_plane_points
                    tq = (k - 1 + 0.5_dp*(m - 1))/nt
                    call shell_jet(q, tq, u(iphi + 1:icirc - 1), value, grad3, dval, dgrad, s, f, vel, &
                                   lt)
                    detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
                    call source_b([s(1), 2*pi*tq, s(2)], bs)
                    call filament_b([s(1), 2*pi*tq, s(2)], bf)
                    b = bs + c*bf + grad3
                    w = quad_fref%weight(vol_points(q))*abs(detf)*s(1)*2*pi/(per_plane_points*nt)
                    do j = 1, nt
                        do i = 1, nloc
                            if (shell_slot(rows_vol%ci(a + i - 1)) == 0) cycle
                            rowi = iphi + (j - 1)*nphi + shell_slot(rows_vol%ci(a + i - 1))
                            r(rowi) = r(rowi) + w*dot_product(b, dgrad(:, j, i))
                            if (present(jt)) call jt%add(icirc, rowi, &
                                                         w*dot_product(bf, dgrad(:, j, i)))
                        end do
                    end do
                    if (.not. present(jt)) cycle
                    allocate (e(3, nt*nloc), cols(nt*nloc))
                    jj = 0
                    do j = 1, nt
                        do i = 1, nloc
                            jj = jj + 1
                            e(:, jj) = dgrad(:, j, i)
                            cols(jj) = 0
                            if (shell_slot(rows_vol%ci(a + i - 1)) /= 0) cols(jj) = iphi + &
                                (j - 1)*nphi + shell_slot(rows_vol%ci(a + i - 1))
                        end do
                    end do
                    do rowj = 1, nt*nloc
                        if (cols(rowj) == 0) cycle
                        if (all(e(:, rowj) == 0)) cycle
                        do rowi = 1, nt*nloc
                            if (cols(rowi) == 0) cycle
                            if (all(e(:, rowi) == 0)) cycle
                            call jt%add(cols(rowj), cols(rowi), w*dot_product(e(:, rowi), &
                                                                              e(:, rowj)))
                            kmax = max(kmax, abs(w*dot_product(e(:, rowi), e(:, rowj))))
                        end do
                    end do
                    deallocate (e, cols)
                end do
            end do
        end do
        if (wall_model == 'exterior') call exterior_rows(u, r, jt)
        associate (unused => value)
        end associate
    end subroutine shell_rows

    subroutine shell_residual(u, r)
        !! The shell rows of shell_rows without the Hessian, threaded over the
        !! volume points (the source fields at every point dominate).
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        real(dp), allocatable :: dval(:, :), dgrad(:, :, :)
        real(dp) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2), bs(3), b(3), w, &
                    detf, tq, bf(3), c
        integer :: q, m, k, j, i, a, nloc, rowi
        integer, parameter :: per_plane_points = 2

        c = u(icirc)
        ! (default(private) would privatise the host-associated tables)
        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(q, a, nloc, k, m, tq, value, grad3, dval, dgrad, s, f, vel, lt, detf, bs, &
        !$omp         bf, b, w, j, i, rowi) reduction(+:r)
        do q = 1, size(vol_points)
            a = rows_vol%rp(q)
            nloc = rows_vol%rp(q + 1) - a
            do k = 1, nt
                do m = 1, per_plane_points
                    tq = (k - 1 + 0.5_dp*(m - 1))/nt
                    call shell_jet(q, tq, u(iphi + 1:icirc - 1), value, grad3, dval, dgrad, s, f, &
                                   vel, lt)
                    detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
                    call source_b([s(1), 2*pi*tq, s(2)], bs)
                    call filament_b([s(1), 2*pi*tq, s(2)], bf)
                    b = bs + c*bf + grad3
                    w = quad_fref%weight(vol_points(q))*abs(detf)*s(1)*2*pi/(per_plane_points*nt)
                    do j = 1, nt
                        do i = 1, nloc
                            if (shell_slot(rows_vol%ci(a + i - 1)) == 0) cycle
                            rowi = iphi + (j - 1)*nphi + shell_slot(rows_vol%ci(a + i - 1))
                            r(rowi) = r(rowi) + w*dot_product(b, dgrad(:, j, i))
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine shell_residual

    subroutine exterior_rows(u, r, jt)
        !! The exterior energy 1/2 sum_kj u_k^T Q(k - j) u_j of the wall traces
        !! and the source term -oint Phi (B_s + c B_fil).n dS on Gamma_c (the
        !! exterior field of the sources: its cross term with grad Phi_ext).
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        real(dp), allocatable :: tr(:, :), dval(:, :), dgrad(:, :, :), qu(:)
        real(dp) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2), bs(3), bf(3), nrm(2), &
                    tq, w, tvec(2), c
        integer :: nw, k, j, q, i, a, b, ia, ib, d, m, row, col
        integer, parameter :: per_plane_points = 2

        nw = size(qblock, 1)
        c = u(icirc)
        ! wall traces u_k (nw, nt) of Phi per plane
        allocate (tr(nw, nt))
        tr = 0
        do q = 1, nw
            do i = rows_wnodes%rp(q), rows_wnodes%rp(q + 1) - 1
                if (shell_slot(rows_wnodes%ci(i)) == 0) cycle
                do k = 1, nt
                    tr(q, k) = tr(q, k) + rows_wnodes%v0(i)* &
                               u(iphi + (k - 1)*nphi + shell_slot(rows_wnodes%ci(i)))
                end do
            end do
        end do
        do k = 1, nt
            qu = 0*tr(:, 1)
            do j = 1, nt
                d = modulo(k - j, nt)
                qu = qu + matmul(qblock(:, :, d), tr(:, j))
            end do
            do q = 1, nw
                do i = rows_wnodes%rp(q), rows_wnodes%rp(q + 1) - 1
                    if (shell_slot(rows_wnodes%ci(i)) == 0) cycle
                    row = iphi + (k - 1)*nphi + shell_slot(rows_wnodes%ci(i))
                    r(row) = r(row) + rows_wnodes%v0(i)*qu(q)
                end do
            end do
        end do
        if (present(jt)) then
            do k = 1, nt
                do j = 1, nt
                    d = modulo(k - j, nt)
                    do ia = 1, nw
                        do a = rows_wnodes%rp(ia), rows_wnodes%rp(ia + 1) - 1
                            if (shell_slot(rows_wnodes%ci(a)) == 0) cycle
                            row = iphi + (k - 1)*nphi + shell_slot(rows_wnodes%ci(a))
                            do ib = 1, nw
                                if (qblock(ia, ib, d) == 0) cycle
                                do b = rows_wnodes%rp(ib), rows_wnodes%rp(ib + 1) - 1
                                    if (shell_slot(rows_wnodes%ci(b)) == 0) cycle
                                    col = iphi + (j - 1)*nphi + shell_slot(rows_wnodes%ci(b))
                                    call jt%add(col, row, rows_wnodes%v0(a)*qblock(ia, ib, d)* &
                                                rows_wnodes%v0(b))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if
        ! the source term on Gamma_c (the same wall on every plane)
        do q = 1, side_w%point_count
            do k = 1, nt
                do m = 1, per_plane_points
                    tq = (k - 1 + 0.5_dp*(m - 1))/nt
                    call side_jet(rows_wside, q, tq, value, grad3, dval, dgrad, s, f, vel, lt)
                    tvec = matmul(f, side_w%tangent(:, q))
                    w = side_w%weight(q)*norm2(tvec)*s(1)*2*pi/(per_plane_points*nt)
                    tvec = tvec/norm2(tvec)
                    nrm = [tvec(2), -tvec(1)]
                    call source_b([s(1), 2*pi*tq, s(2)], bs)
                    call filament_b([s(1), 2*pi*tq, s(2)], bf)
                    a = rows_wside%rp(q)
                    do i = a, rows_wside%rp(q + 1) - 1
                        if (shell_slot(rows_wside%ci(i)) == 0) cycle
                        do j = 1, nt
                            row = iphi + (j - 1)*nphi + shell_slot(rows_wside%ci(i))
                            r(row) = r(row) - w*dval(j, i - a + 1)* &
                                     ((bs(1) + c*bf(1))*nrm(1) + (bs(3) + c*bf(3))*nrm(2))
                            if (present(jt)) call jt%add(icirc, row, -w*dval(j, i - a + 1)* &
                                                         (bf(1)*nrm(1) + bf(3)*nrm(2)))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine exterior_rows

    subroutine side_jet(rows, q, tq, value, grad3, dval, dgrad, s, f, vel, lt)
        !! As shell_jet at the point q of other precomputed shell rows.
        type(point_rows_t), intent(in) :: rows
        integer, intent(in) :: q
        real(dp), intent(in) :: tq
        real(dp), intent(out) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2)
        real(dp), allocatable, intent(out) :: dval(:, :), dgrad(:, :, :)
        real(dp) :: ht(nt, 2), sj(2)
        integer :: a, b, i, j

        call spline_derham_t_basis(nt, tq, lt, ht)
        a = rows%rp(q)
        b = rows%rp(q + 1) - 1
        s = 0
        f = 0
        vel = 0
        do j = 1, nt
            sj = matmul(nodes_f(:, rows%ci(a:b), j), rows%v0(a:b))
            s = s + lt(j, 1)*sj
            vel = vel + lt(j, 2)*sj
            f(:, 1) = f(:, 1) + lt(j, 1)*matmul(nodes_f(:, rows%ci(a:b), j), rows%v1(a:b))
            f(:, 2) = f(:, 2) + lt(j, 1)*matmul(nodes_f(:, rows%ci(a:b), j), rows%v2(a:b))
        end do
        allocate (dval(nt, b - a + 1), dgrad(3, nt, b - a + 1))
        dgrad = 0
        do i = a, b
            do j = 1, nt
                dval(j, i - a + 1) = rows%v0(i)*lt(j, 1)
            end do
        end do
        value = 0
        grad3 = 0
    end subroutine side_jet

    subroutine shell_jet(q, tq, phi, value, grad3, dval, dgrad, s, f, vel, lt)
        !! The shell potential at volume point q and toroidal position tq (any
        !! t, not only planes): as plane_jet with L_j(tq).
        integer, intent(in) :: q
        real(dp), intent(in) :: tq, phi(:)
        real(dp), intent(out) :: value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2)
        real(dp), allocatable, intent(out) :: dval(:, :), dgrad(:, :, :)
        real(dp) :: ht(nt, 2), ft(2, 2), detf, gx(2), sj(2)
        integer :: a, b, i, j, col, nloc

        call spline_derham_t_basis(nt, tq, lt, ht)
        a = rows_vol%rp(q)
        b = rows_vol%rp(q + 1) - 1
        nloc = b - a + 1
        s = 0
        f = 0
        vel = 0
        do j = 1, nt
            sj = matmul(nodes_f(:, rows_vol%ci(a:b), j), rows_vol%v0(a:b))
            s = s + lt(j, 1)*sj
            vel = vel + lt(j, 2)*sj
            f(:, 1) = f(:, 1) + lt(j, 1)*matmul(nodes_f(:, rows_vol%ci(a:b), j), rows_vol%v1(a:b))
            f(:, 2) = f(:, 2) + lt(j, 1)*matmul(nodes_f(:, rows_vol%ci(a:b), j), rows_vol%v2(a:b))
        end do
        detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
        ft = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/detf
        value = 0
        allocate (dval(nt, nloc), dgrad(3, nt, nloc))
        dval = 0
        dgrad = 0
        do i = a, b
            col = shell_slot(rows_vol%ci(i))
            gx = matmul(ft, [rows_vol%v1(i), rows_vol%v2(i)])
            do j = 1, nt
                dval(j, i - a + 1) = rows_vol%v0(i)*lt(j, 1)
                dgrad(1, j, i - a + 1) = gx(1)*lt(j, 1)
                dgrad(3, j, i - a + 1) = gx(2)*lt(j, 1)
                dgrad(2, j, i - a + 1) = (rows_vol%v0(i)*lt(j, 2) - dot_product(vel, gx)*lt(j, 1)) &
                                         /(2*pi*s(1))
                if (col /= 0) then
                    value = value + dval(j, i - a + 1)*phi((j - 1)*nphi + col)
                end if
            end do
        end do
        grad3 = 0
        do i = a, b
            col = shell_slot(rows_vol%ci(i))
            if (col == 0) cycle
            do j = 1, nt
                grad3 = grad3 + dgrad(:, j, i - a + 1)*phi((j - 1)*nphi + col)
            end do
        end do
    end subroutine shell_jet

    subroutine balance_rows(u, r, jt)
        !! Per plane k and mode m: oint (p + |B_in|^2/2) w_m dl - oint
        !! |B_s + grad Phi|^2/2 w_m dl, with the chain rule in chi, g, Phi.
        real(dp), intent(in) :: u(:)
        real(dp), intent(inout) :: r(:)
        type(triplet_list), intent(inout), optional :: jt
        real(dp), allocatable :: dvc(:, :), dgc(:, :, :), dvg(:, :), dgg(:, :, :)
        real(dp) :: c, gc(3), gg(3), s(2), f(2, 2), vel(2), lt(nt, 2), hv(3), fp(3), gp(3), &
                    pr(3), bv(3), q, dq_c, dq_gc(3), dq_gg(3), wl, theta, tvec(2), bs(3), &
                    gphi(3), xk(2, nt), mf(2, 2), mdot(2, 2), uu(2), udot(2), rt(2), gx(2), &
                    xdot(2), e(3), db(3), rho2, dqp(3), value, t3(3), dc_c, dc_gc(3), &
                    dc_gg(3), bf(3)
        integer :: k, qq, m, a, i, j, nloc, row, col, b
        integer, allocatable :: pslot(:)

        allocate (pslot(n0p))
        pslot = [(i, i=1, n0p)]
        r(icirc) = u(icirc)*filament_circulation
        if (present(jt)) call jt%add(icirc, icirc, filament_circulation)
        xk = reshape(u(ix + 1:ix + 2*nt), [2, nt])
        do k = 1, nt
            ! plasma side
            do qq = 1, side_p%point_count
                call plane_jet(rows_wall, qq, nodes_p, k, u(ichi + 1:ig), n0p, pslot, c, gc, &
                               dvc, dgc, s, f, vel, lt)
                call plane_jet(rows_wall, qq, nodes_p, k, u(ig + 1:ix), n0p, pslot, value, gg, &
                               dvg, dgg, s, f, vel, lt)
                ! winding form about X(t), frame M(t), at t_k
                mf = frames(:, :, k)
                mdot = 0
                xdot = 0
                do j = 1, nt
                    mdot = mdot + lt(j, 2)*frames(:, :, j)
                    xdot = xdot + lt(j, 2)*xk(:, j)
                end do
                rt = s - xk(:, k)
                uu = matmul(mf, rt)
                udot = matmul(mdot, rt) - matmul(mf, xdot)
                rho2 = sum(uu**2)
                gx = matmul(transpose(mf), [-uu(2), uu(1)])/rho2
                hv = [gx(1), (uu(1)*udot(2) - uu(2)*udot(1))/(rho2*2*pi*s(1)), gx(2)]
                gphi = [0.0_dp, 1/s(1), 0.0_dp]
                call section_profiles(c, fp, gp, pr)
                bv = cross3(gg, gc) + fp(1)*cross3(gc, hv) + gp(1)*cross3(gc, gphi)
                q = pr(1) + sum(bv**2)/2
                dq_c = pr(2) + dot_product(bv, fp(2)*cross3(gc, hv) + gp(2)*cross3(gc, gphi))
                do b = 1, 3
                    e = 0
                    e(b) = 1
                    db = cross3(gg, e) + fp(1)*cross3(e, hv) + gp(1)*cross3(e, gphi)
                    dq_gc(b) = dot_product(bv, db)
                    dq_gg(b) = dot_product(bv, cross3(e, gc))
                end do
                tvec = matmul(f, side_p%tangent(:, qq))
                wl = side_p%weight(qq)*norm2(tvec)
                theta = atan2(s(2) - centre(2, k), s(1) - centre(1, k))
                a = rows_wall%rp(qq)
                nloc = rows_wall%rp(qq + 1) - a
                ! the current row: c oint B_fil.dl = mean_k oint B_in.dl (the
                ! tangent with its length: dl = F t ds0)
                t3 = [tvec(1), 0.0_dp, tvec(2)]*side_p%weight(qq)
                r(icirc) = r(icirc) - dot_product(bv, t3)/nt
                if (present(jt)) then
                    dc_c = dot_product(fp(2)*cross3(gc, hv) + gp(2)*cross3(gc, gphi), t3)
                    do b = 1, 3
                        e = 0
                        e(b) = 1
                        dc_gc(b) = dot_product(cross3(gg, e) + fp(1)*cross3(e, hv) + &
                                              gp(1)*cross3(e, gphi), t3)
                        dc_gg(b) = dot_product(cross3(e, gc), t3)
                    end do
                    do j = 1, nt
                        do i = 1, nloc
                            col = ichi + (j - 1)*n0p + rows_wall%ci(a + i - 1)
                            call jt%add(col, icirc, -(dc_c*dvc(j, i) + &
                                        dot_product(dc_gc, dgc(:, j, i)))/nt)
                            col = ig + (j - 1)*n0p + rows_wall%ci(a + i - 1)
                            call jt%add(col, icirc, -dot_product(dc_gg, dgg(:, j, i))/nt)
                        end do
                    end do
                end if
                do m = 1, nrho
                    row = irho + (k - 1)*nrho + m
                    r(row) = r(row) + wl*q*mode(m, theta)
                    if (.not. present(jt)) cycle
                    do j = 1, nt
                        do i = 1, nloc
                            col = ichi + (j - 1)*n0p + rows_wall%ci(a + i - 1)
                            call jt%add(col, row, wl*mode(m, theta)*(dq_c*dvc(j, i) + &
                                        dot_product(dq_gc, dgc(:, j, i))))
                            col = ig + (j - 1)*n0p + rows_wall%ci(a + i - 1)
                            call jt%add(col, row, wl*mode(m, theta)*dot_product(dq_gg, dgg(:, j, i)))
                        end do
                    end do
                end do
            end do
            ! shell side: B_s + grad Phi on the interface (block 6 points)
            do qq = 1, size(iface_points)
                call plane_jet(rows_iface, qq, nodes_f, k, u(iphi + 1:icirc - 1), nphi, shell_slot, &
                               value, gc, dvc, dgc, s, f, vel, lt)
                call source_b([s(1), plane_phi(k), s(2)], bs)
                call filament_b([s(1), plane_phi(k), s(2)], bf)
                bv = bs + u(icirc)*bf + gc
                q = sum(bv**2)/2
                tvec = matmul(f, side_f%tangent(:, iface_points(qq)))
                wl = side_f%weight(iface_points(qq))*norm2(tvec)
                theta = atan2(s(2) - centre(2, k), s(1) - centre(1, k))
                a = rows_iface%rp(qq)
                nloc = rows_iface%rp(qq + 1) - a
                do m = 1, nrho
                    row = irho + (k - 1)*nrho + m
                    r(row) = r(row) - wl*q*mode(m, theta)
                    if (.not. present(jt)) cycle
                    call jt%add(icirc, row, -wl*mode(m, theta)*dot_product(bv, bf))
                    do j = 1, nt
                        do i = 1, nloc
                            if (shell_slot(rows_iface%ci(a + i - 1)) == 0) cycle
                            col = iphi + (j - 1)*nphi + shell_slot(rows_iface%ci(a + i - 1))
                            dqp = bv
                            call jt%add(col, row, -wl*mode(m, theta)*dot_product(dqp, dgc(:, j, i)))
                        end do
                    end do
                end do
            end do
        end do
        if (present(jt)) then
            do i = irho + 1, nu
                call jt%add(i, i, 0.0_dp)
            end do
        end if
    end subroutine balance_rows

    ! ------------------------------------------------------------------ Newton

    subroutine rebuild(u)
        !! Node coordinates of every plane for u.
        real(dp), intent(in) :: u(:)
        integer :: k

        do k = 1, nt
            call plane_nodes(k, u(irho + (k - 1)*nrho + 1:irho + k*nrho), u(ix + 2*k - 1:ix + 2*k))
        end do
    end subroutine rebuild

    subroutine newton(u, fixed, steps, check_now)
        real(dp), intent(inout) :: u(:)
        integer, intent(in) :: fixed(:)
        integer, intent(out) :: steps
        logical, intent(in) :: check_now
        real(dp), allocatable :: r(:), delta(:), history(:), rt(:), dbar(:)
        real(dp) :: theta, lin0
        real(dp), parameter :: trust = 0.05_dp
        logical :: svd_now
        real(dp), allocatable :: lres(:)
        type(eu_constrained_t) :: fac
        type(triplet_list) :: jt
        type(csr_t) :: jac_t, sp, sf
        logical, allocatable :: is_fixed(:)
        real(dp) :: alpha
        integer :: it, k, st, j, check_at
        character(8) :: env

        allocate (is_fixed(nu))
        is_fixed = .false.
        is_fixed(fixed) = .true.
        check_at = 1
        ! diagnostics: MILCA_CHECK_AT=N checks the Jacobian at step N,
        ! MILCA_SVD also lists the weak modes of the Newton matrix there
        call get_environment_variable('MILCA_SVD', env, status=st)
        svd_now = st == 0
        call get_environment_variable('MILCA_CHECK_AT', env, status=st)
        if (st == 0) read (env, *) check_at
        allocate (history(0))
        steps = 0
        do it = 1, max_steps
            call rebuild(u)
            if (allocated(jt%rows)) deallocate (jt%rows, jt%cols, jt%vals)
            jt%count = 0
            call jt%init(64*nu)
            call residual(u, r, jt)
            r(fixed) = 0
            history = [history, norm2(r)]
            print '(a,i3,a,es10.3,a,5es10.2)', '  step', it, '  |r|', history(it), &
                '  chi, g, axis, shell, balance', norm2(r(ichi + 1:ig)), norm2(r(ig + 1:ix)), &
                norm2(r(ix + 1:iphi)), norm2(r(iphi + 1:irho)), norm2(r(irho + 1:))
            if (history(it) < 1e-7_dp) exit
            ! X and Gamma_p columns: each rebuilds one plane
            if (analytic_columns) then
                call shape_matrices(u, sp, sf, jt)
                do j = ix + 1, ix + 2*nt
                    call analytic_column(u, jt, is_fixed, j, (j - ix + 1)/2, sp, sf)
                end do
                do j = irho + 1, nu
                    call analytic_column(u, jt, is_fixed, j, (j - irho - 1)/nrho + 1, sp, sf)
                end do
            else
                do j = ix + 1, ix + 2*nt
                    call difference_column(u, jt, is_fixed, j, (j - ix + 1)/2)
                end do
                do j = irho + 1, nu
                    call difference_column(u, jt, is_fixed, j, (j - irho - 1)/nrho + 1)
                end do
            end if
            call jt%to_csr(nu, nu, jac_t)
            if (check_now .and. it == check_at) call check_jacobian(u, jac_t)
            if (check_now .and. it == check_at .and. svd_now) then
                if (nu - size(fixed) <= 3000) then
                    call weak_modes(jac_t, is_fixed)
                else
                    call weakest_mode(jac_t, fixed)
                end if
            end if
            allocate (delta(nu), dbar(nu))
            call eu_constrained_factor(jac_t, fixed, fac, st)
            if (st /= 0) error stop 'Newton factorisation'
            call eu_constrained_apply(fac, -r, delta, st)
            if (st /= 0) error stop 'Newton solve'
            ! iterative refinement against the unfactored matrix (the free
            ! block is badly conditioned; the direct factor alone leaves a
            ! linear residual that slows Newton to linear convergence)
            do k = 1, 3
                lres = -r - jac_vec(jac_t, delta)
                lres(fixed) = 0
                if (k == 1) lin0 = norm2(lres)/max(norm2(r), tiny(1.0_dp))
                call eu_constrained_apply(fac, lres, dbar, st)
                delta = delta + dbar
            end do
            lres = -r - jac_vec(jac_t, delta)
            lres(fixed) = 0
            print '(a,2es9.2)', '    linear residual before/after refinement', lin0, &
                norm2(lres)/max(norm2(r), tiny(1.0_dp))
            print '(a,6es10.2)', '    |d| chi, g, axis, mu, shell, boundary', &
                maxval(abs(delta(ichi + 1:ig))), maxval(abs(delta(ig + 1:ix))), &
                maxval(abs(delta(ix + 1:imu))), maxval(abs(delta(imu + 1:iphi))), &
                maxval(abs(delta(iphi + 1:irho))), maxval(abs(delta(irho + 1:)))
            alpha = min(1.0_dp, trust*u(irho + 1)/max(maxval(abs(delta(irho + 1:))), &
                                                      tiny(1.0_dp)), &
                        trust/max(maxval(abs(delta(ix + 1:imu))), tiny(1.0_dp)))
            do k = 1, 12
                call rebuild(u + alpha*delta)
                call residual(u + alpha*delta, rt)
                rt(fixed) = 0
                theta = norm2(rt)/history(it)
                if (theta < 1 - alpha/4) exit
                alpha = alpha/2
            end do
            call eu_constrained_free(fac)
            print '(a,f7.4,a,es9.2)', '    alpha', alpha, '  theta', theta
            u = u + alpha*delta
            deallocate (delta, dbar)
            steps = steps + 1
        end do
        call rebuild(u)
        print '(a,20es9.1)', '  |r| history', history
    end subroutine newton

    subroutine weakest_mode(jac_t, fixed)
        !! The eigenvector of the Newton matrix with the smallest |eigenvalue|
        !! by inverse iteration with its factor, for systems too large for a
        !! dense SVD: block shares, axis shifts and boundary modes per plane.
        type(csr_t), intent(in) :: jac_t
        integer, intent(in) :: fixed(:)
        type(eu_constrained_t) :: c
        real(dp), allocatable :: v(:), w(:)
        real(dp) :: lam
        integer :: it, st, i, k, b, lo(7), hi(7)

        call eu_constrained_factor(jac_t, fixed, c, st)
        allocate (v(nu), w(nu))
        v = [(sin(0.7_dp*i) + 0.3_dp, i=1, nu)]
        v(fixed) = 0
        v = v/norm2(v)
        do it = 1, 60
            call eu_constrained_apply(c, v, w, st)
            lam = 1/norm2(w)
            v = w*lam
        end do
        call eu_constrained_free(c)
        w = jac_vec(jac_t, v)
        w(fixed) = 0
        lo = [ichi + 1, ig + 1, ix + 1, imu + 1, iphi + 1, icirc, irho + 1]
        hi = [ig, ix, imu, iphi, icirc - 1, icirc, nu]
        print '(a,es10.2,a,es10.2)', '  weakest mode: |J v|', norm2(w), '  1/|J^-1 v|', lam
        print '(a,7f8.3)', '    shares chi g axis mu Phi c Gamma_p', &
            (norm2(v(lo(b):hi(b)))**2, b=1, 7)
        do k = 1, nt
            print '(a,i3,a,2es10.2,a,6es10.2)', '    plane', k, '  axis', &
                v(ix + 2*k - 1:ix + 2*k), '  rho', v(irho + (k - 1)*nrho + 1:irho + min(6, nrho) + &
                (k - 1)*nrho)
        end do
    end subroutine weakest_mode

    subroutine weak_modes(jac_t, is_fixed)
        !! The smallest singular values of the Newton matrix on the free
        !! unknowns, each with the share of its right singular vector (the
        !! weak direction of the unknowns) and its left one (the residual
        !! rows that do not see it) per block.
        type(csr_t), intent(in) :: jac_t
        logical, intent(in) :: is_fixed(:)
        interface
            subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, &
                              info)
                import :: dp
                character, intent(in) :: jobu, jobvt
                integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
                integer, intent(out) :: info
            end subroutine dgesvd
        end interface
        real(dp), allocatable :: a(:, :), sv(:), uu(:, :), vt(:, :), work(:), full(:)
        integer, allocatable :: free(:), loc(:)
        integer :: nf, i, p, info, m, b
        integer :: lo(7), hi(7)
        character(8), parameter :: names(7) = ['chi     ', 'g       ', 'axis    ', 'mu      ', &
                                               'Phi     ', 'c       ', 'Gamma_p ']

        free = pack([(i, i=1, nu)], .not. is_fixed)
        nf = size(free)
        allocate (loc(nu))
        loc = 0
        loc(free) = [(i, i=1, nf)]
        allocate (a(nf, nf))
        a = 0
        do i = 1, jac_t%nrow
            if (loc(i) == 0) cycle
            do p = jac_t%row_ptr(i), jac_t%row_ptr(i + 1) - 1
                if (loc(jac_t%col(p)) > 0) a(loc(jac_t%col(p)), loc(i)) = jac_t%val(p)
            end do
        end do
        allocate (sv(nf), uu(nf, nf), vt(nf, nf), work(1))
        call dgesvd('A', 'A', nf, nf, a, nf, sv, uu, nf, vt, nf, work, -1, info)
        deallocate (work)
        allocate (work(int(work_size(a, nf))))
        call dgesvd('A', 'A', nf, nf, a, nf, sv, uu, nf, vt, nf, work, size(work), info)
        if (info /= 0) error stop 'weak_modes: dgesvd'
        lo = [ichi + 1, ig + 1, ix + 1, imu + 1, iphi + 1, icirc, irho + 1]
        hi = [ig, ix, imu, iphi, icirc - 1, icirc, nu]
        print '(a,i6,a,es10.2)', '  weak modes: free unknowns', nf, '  largest singular value', &
            sv(1)
        print '(a,7a9)', '    sigma/sigma1   side ', names
        allocate (full(nu))
        do m = nf, max(1, nf - 11), -1
            full = 0
            full(free) = vt(m, :)
            print '(es14.2,a,7f9.3)', sv(m)/sv(1), '  unknowns', &
                (norm2(full(lo(b):hi(b)))**2, b=1, 7)
            full = 0
            full(free) = uu(:, m)
            print '(14x,a,7f9.3)', '  rows    ', (norm2(full(lo(b):hi(b)))**2, b=1, 7)
        end do
    end subroutine weak_modes

    real(dp) function work_size(a, n)
        !! Workspace for dgesvd with full U and V^T (the query result is not
        !! kept across the reallocation).
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: n

        work_size = max(1.0_dp, real(5*n + n*64, dp))
        if (size(a) < 0) work_size = 0
    end function work_size

    function fixed_wall_rows() result(rows)
        !! The wall rows of chi on every plane (fixed labels).
        integer, allocatable :: rows(:)
        integer :: k, i

        rows = [((ichi + (k - 1)*n0p + wall_p(i), i=1, size(wall_p)), k=1, nt)]
    end function fixed_wall_rows

    function jac_vec(jac_t, v) result(y)
        !! J v with J stored by rows of J^T.
        type(csr_t), intent(in) :: jac_t
        real(dp), intent(in) :: v(:)
        real(dp) :: y(size(v))
        integer :: i, p

        y = 0
        do i = 1, jac_t%nrow
            do p = jac_t%row_ptr(i), jac_t%row_ptr(i + 1) - 1
                y(jac_t%col(p)) = y(jac_t%col(p)) + jac_t%val(p)*v(i)
            end do
        end do
    end function jac_vec

    subroutine shape_matrices(u, sp, sf, jt)
        !! The plasma (clebsch_3d_shape_system) and shell shape matrices at u,
        !! and the shell's Hessian and c column into the Newton triplets jt.
        real(dp), intent(in) :: u(:)
        type(csr_t), intent(out) :: sp, sf
        type(triplet_list), intent(inout) :: jt
        real(dp), allocatable :: grad(:)
        real(dp) :: energy
        type(csr_t) :: hessian

        call clebsch_3d_shape_system(derham_ref, aquad_ref, nt, section_profiles, u(ichi + 1:ig), &
                                     u(ig + 1:ix), reshape(u(ix + 1:ix + 2*nt), [2, nt]), nodes_p, &
                                     frames, energy, grad, hessian, sp)
        call shell_shape(u, sf, jt)
    end subroutine shape_matrices

    subroutine shell_shape(u, sf, jn)
        !! The shape matrix of the shell rows: sf(row, col) = d r_row/d gamma_col
        !! with row a global unknown (the Phi rows) and gamma the shell node
        !! coordinates, column ((k - 1) nfull + i - 1) 2 + c, from the
        !! generated reference leaf shell_energy_ref and the gradient of the
        !! source field at the moving point.
        real(dp), intent(in) :: u(:)
        type(csr_t), intent(out) :: sf
        type(triplet_list), intent(inout) :: jn
        type(triplet_list) :: t, tl, jl
        real(dp), allocatable :: as(:, :, :, :, :), tb(:, :, :), sx(:, :, :), sy(:, :, :), &
                                 ah(:, :, :, :, :), ac(:, :, :)
        integer, allocatable :: dofs(:)
        real(dp) :: lt(nt, 2), ht(nt, 2), tq, lv(3), s(2), fm(2, 2), vel(2), bt(3), bp(3), &
                    bm(3), db(3, 2), w, gl(3), hh(6), hm(33), gmat(3, 11), g8(3, 8), wt, c, &
                    sj(2), bf(3), hmat(3, 3), dcb(3), kmax
        real(dp), parameter :: hb = 1e-6_dp
        integer :: cell, p0, q0, q, p, mt, nq2, nloc, k, i, j, a, ct1, ct2, kr, kc, row, col, cc, &
                   col2
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info

        c = u(icirc)
        nq2 = 2*nt
        kmax = 0
        call t%init(64*nu)
        !$omp parallel default(shared) &
        !$omp private(tl, jl, cell, p0, q0, dofs, nloc, as, ah, ac, tb, sx, sy, mt, tq, lt, ht, p, q, &
        !$omp         a, s, fm, vel, lv, k, sj, i, j, bt, cc, bp, bm, db, w, gl, hh, hm, gmat, g8, &
        !$omp         hmat, dcb, bf, wt, ct1, ct2, kr, kc, row, col, col2) reduction(max:kmax)
        call tl%init(16*nu)
        call jl%init(16*nu)
        !$omp do schedule(dynamic)
        do cell = 1, size(quad_fref%point_offset)
            p0 = quad_fref%point_offset(cell)
            if (vol_index(p0) == 0) cycle
            q0 = vol_index(p0)
            dofs = rows_vol%ci(rows_vol%rp(q0):rows_vol%rp(q0 + 1) - 1)
            nloc = size(dofs)
            allocate (as(nloc, 2*nloc, 2, 2, nq2), tb(nt, 2, nq2), sx(3, nloc, 2), &
                      sy(8, 2*nloc, 2), ah(nloc, nloc, 2, 2, nq2), ac(nloc, 2, nq2))
            as = 0
            ah = 0
            ac = 0
            do mt = 1, nq2
                tq = (mt - 1)/real(nq2, dp)
                call spline_derham_t_basis(nt, tq, lt, ht)
                tb(:, :, mt) = lt
                do p = p0, p0 + quad_fref%points_per_cell - 1
                    q = vol_index(p)
                    a = rows_vol%rp(q)
                    if (any(rows_vol%ci(a:a + nloc - 1) /= dofs)) error stop 'shell shape: cell dofs'
                    s = 0
                    fm = 0
                    vel = 0
                    lv = 0
                    do k = 1, nt
                        sj = matmul(nodes_f(:, dofs, k), rows_vol%v0(a:a + nloc - 1))
                        s = s + lt(k, 1)*sj
                        vel = vel + lt(k, 2)*sj
                        fm(:, 1) = fm(:, 1) + lt(k, 1)*matmul(nodes_f(:, dofs, k), &
                                                             rows_vol%v1(a:a + nloc - 1))
                        fm(:, 2) = fm(:, 2) + lt(k, 1)*matmul(nodes_f(:, dofs, k), &
                                                             rows_vol%v2(a:a + nloc - 1))
                        do i = 1, nloc
                            if (shell_slot(dofs(i)) == 0) cycle
                            j = iphi + (k - 1)*nphi + shell_slot(dofs(i))
                            lv(1) = lv(1) + lt(k, 1)*rows_vol%v1(a + i - 1)*u(j)
                            lv(2) = lv(2) + lt(k, 1)*rows_vol%v2(a + i - 1)*u(j)
                            lv(3) = lv(3) + lt(k, 2)*rows_vol%v0(a + i - 1)*u(j)
                        end do
                    end do
                    call source_total(s, tq, c, bt)
                    do cc = 1, 2
                        call source_total(s + hb*merge([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp], cc == 1), &
                                          tq, c, bp)
                        call source_total(s - hb*merge([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp], cc == 1), &
                                          tq, c, bm)
                        db(:, cc) = (bp - bm)/(2*hb)
                    end do
                    call shell_energy_ref(lv(1), lv(2), lv(3), s(1), s(2), fm(1, 1), fm(2, 1), &
                        fm(1, 2), fm(2, 2), vel(1), vel(2), bt(1), bt(2), bt(3), w, gl(1), &
                        gl(2), gl(3), hh(1), hh(2), hh(3), hh(4), hh(5), hh(6), hm(1), hm(2), &
                        hm(3), hm(4), hm(5), hm(6), hm(7), hm(8), hm(9), hm(10), hm(11), &
                        hm(12), hm(13), hm(14), hm(15), hm(16), hm(17), hm(18), hm(19), hm(20), &
                        hm(21), hm(22), hm(23), hm(24), hm(25), hm(26), hm(27), hm(28), hm(29), &
                        hm(30), hm(31), hm(32), hm(33))
                    gmat = transpose(reshape(hm, [11, 3]))
                    hmat = reshape([hh(1), hh(2), hh(3), hh(2), hh(4), hh(5), hh(3), hh(5), &
                                    hh(6)], [3, 3])
                    call filament_b([s(1), 2*pi*tq, s(2)], bf)
                    dcb = matmul(gmat(:, 9:11), bf)
                    g8 = gmat(:, 1:8)
                    do cc = 1, 2
                        g8(:, cc) = g8(:, cc) + matmul(gmat(:, 9:11), db(:, cc))
                    end do
                    wt = quad_fref%weight(p)/nq2
                    sx = 0
                    sy = 0
                    do i = 1, nloc
                        sx(1, i, 1) = rows_vol%v1(a + i - 1)
                        sx(2, i, 1) = rows_vol%v2(a + i - 1)
                        sx(3, i, 2) = rows_vol%v0(a + i - 1)
                        do cc = 1, 2
                            sy(cc, 2*(i - 1) + cc, 1) = rows_vol%v0(a + i - 1)
                            sy(2 + cc, 2*(i - 1) + cc, 1) = rows_vol%v1(a + i - 1)
                            sy(4 + cc, 2*(i - 1) + cc, 1) = rows_vol%v2(a + i - 1)
                            sy(6 + cc, 2*(i - 1) + cc, 2) = rows_vol%v0(a + i - 1)
                        end do
                    end do
                    do ct1 = 1, 2
                        ac(:, ct1, mt) = ac(:, ct1, mt) + wt*matmul(transpose(sx(:, :, ct1)), dcb)
                        do ct2 = 1, 2
                            as(:, :, ct1, ct2, mt) = as(:, :, ct1, ct2, mt) + &
                                wt*matmul(transpose(sx(:, :, ct1)), matmul(g8, sy(:, :, ct2)))
                            ah(:, :, ct1, ct2, mt) = ah(:, :, ct1, ct2, mt) + &
                                wt*matmul(transpose(sx(:, :, ct1)), matmul(hmat, sx(:, :, ct2)))
                        end do
                    end do
                end do
            end do
            do mt = 1, nq2
                do ct1 = 1, 2
                    do ct2 = 1, 2
                        do kr = 1, nt
                            do kc = 1, nt
                                w = tb(kr, ct1, mt)*tb(kc, ct2, mt)
                                if (w == 0) cycle
                                do i = 1, nloc
                                    if (shell_slot(dofs(i)) == 0) cycle
                                    row = iphi + (kr - 1)*nphi + shell_slot(dofs(i))
                                    do j = 1, 2*nloc
                                        if (as(i, j, ct1, ct2, mt) == 0) cycle
                                        col = ((kc - 1)*nfull + dofs((j - 1)/2 + 1) - 1)*2 + &
                                              mod(j - 1, 2) + 1
                                        call tl%add(row, col, w*as(i, j, ct1, ct2, mt))
                                    end do
                                    do j = 1, nloc
                                        if (shell_slot(dofs(j)) == 0) cycle
                                        if (ah(i, j, ct1, ct2, mt) == 0) cycle
                                        col2 = iphi + (kc - 1)*nphi + shell_slot(dofs(j))
                                        ! the Newton triplets are transposed
                                        call jl%add(col2, row, w*ah(i, j, ct1, ct2, mt))
                                        kmax = max(kmax, abs(w*ah(i, j, ct1, ct2, mt)))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            ! the c column: d r_Phi/d c = int B_fil . d grad Phi (per plane)
            do mt = 1, nq2
                do ct1 = 1, 2
                    do kr = 1, nt
                        if (tb(kr, ct1, mt) == 0) cycle
                        do i = 1, nloc
                            if (shell_slot(dofs(i)) == 0) cycle
                            row = iphi + (kr - 1)*nphi + shell_slot(dofs(i))
                            call jl%add(icirc, row, tb(kr, ct1, mt)*ac(i, ct1, mt))
                        end do
                    end do
                end do
            end do
            deallocate (as, tb, sx, sy, ah, ac)
        end do
        !$omp end do
        !$omp critical (shell_shape_merge)
        do i = 1, tl%count
            call t%add(tl%rows(i), tl%cols(i), tl%vals(i))
        end do
        do i = 1, jl%count
            call jn%add(jl%rows(i), jl%cols(i), jl%vals(i))
        end do
        !$omp end critical (shell_shape_merge)
        !$omp end parallel
        call t%add(nu, 2*nfull*nt, 0.0_dp)
        call csc_from_triplet(2*nfull*nt, nu, t%cols(:t%count), t%rows(:t%count), &
                              t%vals(:t%count), mm, info)
        if (info%code /= 0) error stop 'shell shape triplets'
        sf%nrow = nu
        sf%ncol = 2*nfull*nt
        sf%row_ptr = mm%col_ptr
        sf%col = mm%row_idx
        sf%val = mm%val
        associate (unused => bf)
        end associate
    end subroutine shell_shape

    subroutine source_total(s, tq, c, b)
        !! B_s + c B_fil at the section point s of toroidal position tq.
        real(dp), intent(in) :: s(2), tq, c
        real(dp), intent(out) :: b(3)
        real(dp) :: bf(3)

        call source_b([s(1), 2*pi*tq, s(2)], b)
        call filament_b([s(1), 2*pi*tq, s(2)], bf)
        b = b + c*bf
    end subroutine source_total

    subroutine trace_rows(u, r)
        !! The rows that the volume shape matrices do not carry: the axis rows
        !! and the boundary (balance and current) rows.
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: r(:)

        allocate (r(nu))
        r = 0
        call axis_rows(u, r)
        call balance_rows(u, r)
    end subroutine trace_rows

    subroutine analytic_column(u, jt, is_fixed, col, plane_k, sp, sf)
        !! Column col (an axis coordinate or a boundary coefficient of plane
        !! plane_k): the volume rows through the shape matrices and the node
        !! derivatives (differences of the plane rebuild only, no residual),
        !! the axis point directly; the trace rows by differences.
        real(dp), intent(in) :: u(:)
        type(triplet_list), intent(inout) :: jt
        logical, intent(in) :: is_fixed(:)
        integer, intent(in) :: col, plane_k
        type(csr_t), intent(in) :: sp, sf
        real(dp), allocatable :: up(:), np_p(:, :), np_m(:, :), nf_p(:, :), nf_m(:, :), rp(:), rm(:), &
                                 colv(:)
        real(dp), parameter :: h = 1e-6_dp
        integer :: lo, i, p, j, k, comp, node, axis_col

        if (is_fixed(col)) then
            call jt%add(col, col, 1.0_dp)
            return
        end if
        lo = irho + (plane_k - 1)*nrho
        up = u
        up(col) = u(col) + h
        call plane_nodes(plane_k, up(lo + 1:lo + nrho), up(ix + 2*plane_k - 1:ix + 2*plane_k))
        np_p = nodes_p(:, :, plane_k)
        nf_p = nodes_f(:, :, plane_k)
        call trace_rows(up, rp)
        up(col) = u(col) - h
        call plane_nodes(plane_k, up(lo + 1:lo + nrho), up(ix + 2*plane_k - 1:ix + 2*plane_k))
        np_m = nodes_p(:, :, plane_k)
        nf_m = nodes_f(:, :, plane_k)
        call trace_rows(up, rm)
        call plane_nodes(plane_k, u(lo + 1:lo + nrho), u(ix + 2*plane_k - 1:ix + 2*plane_k))
        allocate (colv(nu))
        colv = (rp - rm)/(2*h)
        np_p = (np_p - np_m)/(2*h)
        nf_p = (nf_p - nf_m)/(2*h)
        ! the axis column of the plasma (the winding form about X)
        axis_col = 0
        if (col > ix .and. col <= ix + 2*nt) axis_col = 2*n0p*nt + (col - ix)
        do i = 1, sp%nrow
            do p = sp%row_ptr(i), sp%row_ptr(i + 1) - 1
                j = sp%col(p)
                if (j == axis_col) then
                    colv(i) = colv(i) + sp%val(p)
                else if (j <= 2*n0p*nt) then
                    k = (j - 1)/(2*n0p) + 1
                    if (k /= plane_k) cycle
                    node = mod((j - 1)/2, n0p) + 1
                    comp = mod(j - 1, 2) + 1
                    colv(i) = colv(i) + sp%val(p)*np_p(comp, node)
                end if
            end do
        end do
        do i = 1, sf%nrow
            do p = sf%row_ptr(i), sf%row_ptr(i + 1) - 1
                j = sf%col(p)
                k = (j - 1)/(2*nfull) + 1
                if (k /= plane_k) cycle
                node = mod((j - 1)/2, nfull) + 1
                comp = mod(j - 1, 2) + 1
                colv(i) = colv(i) + sf%val(p)*nf_p(comp, node)
            end do
        end do
        do k = 1, nt
            colv(ichi + (k - 1)*n0p + wall_p) = 0
        end do
        do i = 1, nu
            if (colv(i) /= 0) call jt%add(col, i, colv(i))
        end do
        call jt%add(col, col, 0.0_dp)
    end subroutine analytic_column

    subroutine difference_column(u, jt, is_fixed, col, plane_k)
        !! Column col of the Newton matrix by central differences; only plane
        !! plane_k is rebuilt.
        real(dp), intent(in) :: u(:)
        type(triplet_list), intent(inout) :: jt
        logical, intent(in) :: is_fixed(:)
        integer, intent(in) :: col, plane_k
        real(dp), allocatable :: up(:), rp(:), rm(:)
        real(dp), parameter :: h = 1e-6_dp
        integer :: i, lo

        if (is_fixed(col)) then
            call jt%add(col, col, 1.0_dp)
            return
        end if
        lo = irho + (plane_k - 1)*nrho
        up = u
        up(col) = u(col) + h
        call plane_nodes(plane_k, up(lo + 1:lo + nrho), up(ix + 2*plane_k - 1:ix + 2*plane_k))
        call residual(up, rp)
        up(col) = u(col) - h
        call plane_nodes(plane_k, up(lo + 1:lo + nrho), up(ix + 2*plane_k - 1:ix + 2*plane_k))
        call residual(up, rm)
        call plane_nodes(plane_k, u(lo + 1:lo + nrho), u(ix + 2*plane_k - 1:ix + 2*plane_k))
        rp = (rp - rm)/(2*h)
        do i = 1, nu
            if (rp(i) /= 0) call jt%add(col, i, rp(i))
        end do
        call jt%add(col, col, 0.0_dp)
    end subroutine difference_column

    subroutine check_jacobian(u, jac_t)
        !! J v against central differences of the residual per block of v.
        real(dp), intent(in) :: u(:)
        type(csr_t), intent(in) :: jac_t
        real(dp), allocatable :: v(:), jv(:), rp(:), rm(:), fd(:)
        integer :: lo(6), hi(6), blk, i, p, k
        character(9), parameter :: names(6) = ['chi      ', 'g        ', 'axis     ', &
                                               'mu       ', 'shell    ', 'boundary ']
        real(dp) :: hv

        lo = [ichi + 1, ig + 1, ix + 1, imu + 1, iphi + 1, irho + 1]
        hi = [ig, ix, imu, iphi, irho, nu]
        allocate (v(nu), jv(nu))
        do blk = 1, 6
            v = 0
            do i = lo(blk), hi(blk)
                v(i) = sin(1.3_dp*i) + 0.5_dp
            end do
            do k = 1, nt
                v(ichi + (k - 1)*n0p + wall_p) = 0
            end do
            hv = 1e-6_dp
            jv = 0
            do i = 1, jac_t%nrow
                do p = jac_t%row_ptr(i), jac_t%row_ptr(i + 1) - 1
                    jv(jac_t%col(p)) = jv(jac_t%col(p)) + jac_t%val(p)*v(i)
                end do
            end do
            call rebuild(u + hv*v)
            call residual(u + hv*v, rp)
            call rebuild(u - hv*v)
            call residual(u - hv*v, rm)
            fd = (rp - rm)/(2*hv)
            do k = 1, nt
                jv(ichi + (k - 1)*n0p + wall_p) = 0
            end do
            print '(a,a9,a,es10.2,a,es10.2)', '  Jacobian check, v in ', names(blk), &
                ': |Jv - FD|/|FD|', norm2(jv - fd)/max(norm2(fd), tiny(1.0_dp)), '  |FD|', &
                norm2(fd)
        end do
        call rebuild(u)
    end subroutine check_jacobian

    subroutine errors(u, eb, ephi, egamma, eaxis)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: eb, ephi, egamma, eaxis
        real(dp), allocatable :: r(:), dval(:, :), dgrad(:, :, :)
        real(dp) :: ferr(2), value, grad3(3), s(2), f(2, 2), vel(2), lt(nt, 2), bs(3), num, &
                    den, w, detf, tq
        integer :: k, q

        call rebuild(u)
        call residual(u, r, ferr=ferr)
        eb = sqrt(ferr(1)/ferr(2))
        num = 0
        den = 0
        do q = 1, size(vol_points)
            do k = 1, nt
                tq = (k - 1)/real(nt, dp)
                call shell_jet(q, tq, u(iphi + 1:icirc - 1), value, grad3, dval, dgrad, s, f, vel, lt)
                call section_exact_b([s(1), 2*pi*tq, s(2)], bs)
                detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
                w = quad_fref%weight(vol_points(q))*abs(detf)*s(1)
                num = num + w*sum(grad3**2)
                den = den + w*sum(bs**2)
            end do
        end do
        ephi = sqrt(num/den)
        egamma = 0
        eaxis = 0
        do k = 1, nt
            egamma = max(egamma, gamma_error(k, u(irho + (k - 1)*nrho + 1:irho + k*nrho)))
            eaxis = max(eaxis, norm2(u(ix + 2*k - 1:ix + 2*k) - centre(:, k)))
        end do
    end subroutine errors

    subroutine run()
        real(dp), allocatable :: u(:)
        real(dp) :: eb, ephi, egamma, eaxis, eb0, egamma0, egamma_start, tk
        integer :: k, steps1, steps2, i
        integer, allocatable :: fixed_wall(:)

        u = u_exact
        ! the exact state on the mesh with g fitted at the exact labels (E4a):
        ! the interpolation level of B
        call rebuild(u)
        call fit_g(u)
        call errors(u, eb0, ephi, egamma0, eaxis)
        print '(a,4es10.2)', 'exact state on the mesh: errors B, grad Phi, Gamma_p, axis', &
            eb0, ephi, egamma0, eaxis
        block
            real(dp), allocatable :: r0(:)
            call rebuild(u)
            call residual(u, r0)
            r0(fixed_wall_rows()) = 0
            print '(a,5es10.2)', 'exact state residual: chi, g, axis, shell, balance', &
                norm2(r0(ichi + 1:ig)), norm2(r0(ig + 1:ix)), norm2(r0(ix + 1:iphi)), &
                norm2(r0(iphi + 1:irho)), norm2(r0(irho + 1:))
        end block
        u = u_exact
        ! perturbed start
        do k = 1, nt
            tk = (k - 1)/real(nt, dp)
            associate (rho => u(irho + (k - 1)*nrho + 1:irho + k*nrho))
                rho(2) = rho(2) + perturbation*rho(1)*(1 + 0.5_dp*cos(2*pi*tk))
                if (nrho >= 3) rho(3) = rho(3) + 0.5_dp*perturbation*rho(1)
                rho(1) = rho(1)*(1 + perturbation)
            end associate
            u(ix + 2*k - 1:ix + 2*k) = centre(:, k) + 0.2_dp*perturbation*[cos(2*pi*tk), &
                                                                             sin(2*pi*tk)]
            associate (c0 => u(ichi + (k - 1)*n0p + 1:ichi + k*n0p))
                c0 = c0*(1 + perturbation*(1 - c0))
            end associate
        end do
        fixed_wall = [((ichi + (k - 1)*n0p + wall_p(i), i=1, size(wall_p)), k=1, nt)]
        ! the circulation: held at 0 with the conducting wall (B_s carries the
        ! plasma current), an unknown with the exterior map
        if (wall_model == 'ideal') then
            fixed_wall = [fixed_wall, icirc]
            u(icirc) = 0
        else
            u(icirc) = c_exact*(1 + perturbation)
        end if
        ! g at fixed labels
        call rebuild(u)
        call fit_g(u)
        call errors(u, eb, ephi, egamma_start, eaxis)
        print '(a,4es10.2)', 'start: errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma_start, &
            eaxis
        ! stage 1: Gamma_p frozen (the fixed-boundary 3D problem of E4, and
        ! the shell), stage 2: all free
        call newton(u, [fixed_wall, gauge_pins, [(irho + i, i=1, nrho*nt)]], steps1, .false.)
        call errors(u, eb, ephi, egamma, eaxis)
        print '(a,i3,a,4es10.2)', 'stage 1 (Gamma_p frozen):', steps1, &
            ' steps; errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma, eaxis
        call newton(u, [fixed_wall, gauge_pins], steps2, jacobian_check)
        call errors(u, eb, ephi, egamma, eaxis)
        print '(a,i3,a,4es10.2)', 'stage 2 (free boundary):', steps2, &
            ' steps; errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma, eaxis
        print '(a,i6,a,4f9.2,a,i6)', 'residuals', n_residual, '  seconds in nodes, clebsch, shell, ' &
            //'balance', t_nodes, t_clebsch, t_shell, t_balance, '  plane rebuilds', n_nodes
        call check(steps2 < max_steps, 'Newton converges')
        ! Gamma_p back to the representation error of the exact boundary by
        ! 2M + 1 modes per plane (its error at the exact state)
        if (wall_model == 'ideal') then
            call check(egamma < 1.2_dp*egamma0 .and. egamma < 0.6_dp*egamma_start, &
                       'Gamma_p returns to the exact boundary within its Fourier truncation')
        else
            ! the exterior state is consistent (shell residual of the exact
            ! state at the discretisation level, like the balance rows), but
            ! on the registered 3-cell mesh Gamma_p is mesh-limited in both
            ! wall models (ideal 1.7e-2, exterior 2.1e-2, truncation 1.6e-2)
            call check(egamma < 1.5_dp*egamma0, &
                       'Gamma_p within 50 % of its Fourier truncation')
        end if
        call check(eb < 1.5_dp*eb0, 'B within 50 % of the exact labels with g fitted')
        if (wall_model == 'exterior') then
            print '(a,es14.6,a,es14.6,a,es10.2)', 'circulation c', u(icirc), '  exact', c_exact, &
                '  relative error', abs(u(icirc)/c_exact - 1)
            call check(abs(u(icirc)/c_exact - 1) < 1e-2_dp, 'the plasma current within 1 %')
        end if
    end subroutine run

    subroutine shape_check()
        !! clebsch_3d_shape_system against clebsch_3d_system (gradient and
        !! Hessian) and its shape matrix against central differences of the
        !! gradient for a perturbation of every node coordinate and axis point.
        real(dp), allocatable :: u(:), g1(:), g2(:), gp(:), gm(:), dn(:, :, :), sd(:), fd(:)
        real(dp) :: e1, e2, xk(2, nt), dx(2, nt), h
        type(csr_t) :: h1, h2, shp
        integer :: k, i, p, j

        u = u_exact
        call rebuild(u)
        call fit_g(u)
        xk = reshape(u(ix + 1:ix + 2*nt), [2, nt])
        call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                               u(ichi + 1:ig), u(ig + 1:ix), e1, g1, h1, xk, nodes_p, frames)
        call clebsch_3d_shape_system(derham_ref, aquad_ref, nt, section_profiles, u(ichi + 1:ig), &
                                     u(ig + 1:ix), xk, nodes_p, frames, e2, g2, h2, shp)
        print '(a,es10.2,a,es10.2,a,es10.2)', 'shape system: |dW|/W', abs(e2 - e1)/abs(e1), &
            '  |d grad|/|grad|', norm2(g2 - g1)/norm2(g1), '  |d H|/|H|', &
            maxval(abs(h2%val - h1%val))/maxval(abs(h1%val))
        call check(abs(e2 - e1) <= 1e-12_dp*abs(e1) .and. norm2(g2 - g1) <= 1e-10_dp*norm2(g1), &
                   'the reference leaf reproduces energy and gradient')
        allocate (dn(2, n0p, nt))
        do k = 1, nt
            do i = 1, n0p
                dn(:, i, k) = [sin(1.7_dp*i + k), cos(2.3_dp*i - k)]
            end do
            dx(:, k) = [cos(1.1_dp*k), sin(0.7_dp*k)]
        end do
        h = 1e-6_dp
        call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                               u(ichi + 1:ig), u(ig + 1:ix), e1, gp, plane_axis=xk + h*dx, &
                               plane_nodes=nodes_p + h*dn, plane_frames=frames)
        call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                               u(ichi + 1:ig), u(ig + 1:ix), e1, gm, plane_axis=xk - h*dx, &
                               plane_nodes=nodes_p - h*dn, plane_frames=frames)
        fd = (gp - gm)/(2*h)
        allocate (sd(shp%nrow))
        sd = 0
        do i = 1, shp%nrow
            do p = shp%row_ptr(i), shp%row_ptr(i + 1) - 1
                j = shp%col(p)
                if (j <= 2*n0p*nt) then
                    sd(i) = sd(i) + shp%val(p)*dn(mod(j - 1, 2) + 1, mod((j - 1)/2, n0p) + 1, &
                                                  (j - 1)/(2*n0p) + 1)
                else
                    sd(i) = sd(i) + shp%val(p)*dx(mod(j - 2*n0p*nt - 1, 2) + 1, &
                                                  (j - 2*n0p*nt - 1)/2 + 1)
                end if
            end do
        end do
        print '(a,es10.2,a,es10.2)', 'shape matrix against differences: |S d - FD|/|FD|', &
            norm2(sd - fd)/norm2(fd), '  |FD|', norm2(fd)
        call check(norm2(sd - fd) <= 1e-6_dp*norm2(fd), 'the shape matrix is the derivative')
    end subroutine shape_check

    subroutine fit_g(u)
        !! g of the current labels, axis and meshes (W is quadratic in g).
        real(dp), intent(inout) :: u(:)
        real(dp), allocatable :: grad(:), work(:)
        real(dp) :: energy
        type(csr_t) :: hessian
        integer :: i, st

        call clebsch_3d_system(derham_ref, aquad_ref, nt, section_profiles, centre(:, 1), &
                               u(ichi + 1:ig), 0*u(ig + 1:ix), energy, grad, hessian, &
                               reshape(u(ix + 1:ix + 2*nt), [2, nt]), nodes_p, frames)
        do i = n0p*nt + 1, 2*n0p*nt
            call add_diagonal(hessian, i, ridge*maxval(abs(hessian%val)))
        end do
        allocate (work(2*n0p*nt))
        call eu_constrained_solve(hessian, -grad, [(i, i=1, n0p*nt)], 0*u(ichi + 1:ig), work, st)
        if (st /= 0) error stop 'g solve'
        u(ig + 1:ix) = work(n0p*nt + 1:)
    end subroutine fit_g

    subroutine add_diagonal(a, i, value)
        type(csr_t), intent(inout) :: a
        integer, intent(in) :: i
        real(dp), intent(in) :: value
        integer :: p

        do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
            if (a%col(p) == i) a%val(p) = a%val(p) + value
        end do
    end subroutine add_diagonal

end program test_eulerian_clebsch_free_boundary_3d
