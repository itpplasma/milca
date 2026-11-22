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
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_3d_system
    use milca_gs_spline, only: gs_spline_project
    use milca_nested_oracles, only: NESTED_SHEARED
    use solovev_label_fields, only: triplet_list
    use landreman_sections, only: section_oracle, section_label, section_exact_b, &
        section_sqrt_hessian, section_axis, section_boundary_radius, section_fit_profiles, &
        section_profiles
    implicit none
    integer, parameter :: degree = 2, nsample = 96
    real(dp), parameter :: pi = acos(-1.0_dp), wall_level = 1.5_dp
    real(dp), parameter :: ridge = 1e-12_dp, newton_ridge = 1e-8_dp
    integer :: n = 4, nt = 9, nmode = 6, max_steps = 14, failures = 0, nrho
    real(dp) :: perturbation = 0.01_dp
    logical :: jacobian_check = .false.
    character(16) :: arg

    ! reference spaces (fixed geometry) and work spaces (per-plane geometry)
    type(spline_space_t), target :: ref_p, ref_f, work_p, work_f
    type(spline_derham_t) :: derham_ref
    type(spline_quadrature_t) :: aquad_ref, quad_fref
    ! exact axes c_k, frames, the Gamma_c samples per plane
    real(dp), allocatable :: centre(:, :), frames(:, :, :), wall_s(:, :, :)
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
    type(point_rows_t) :: rows_wall, rows_iface, rows_vol, rows_axis
    integer, allocatable :: vol_points(:), iface_points(:)
    real(dp), allocatable :: u_exact(:)
    integer :: ichi, ig, ix, imu, iphi, irho, nu
    integer, allocatable :: wall_p(:)
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
        call layout_ogrid(centre(:, k), 0.5_dp*(maxval(s(1, :)) - minval(s(1, :))), 0.5_dp, &
                          layout, st)
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
        !! The node coordinates of plane k for its boundary rho and axis x.
        integer, intent(in) :: k
        real(dp), intent(in) :: rho(:), x(2)
        type(block_layout_t) :: layout
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: v(:), vals(:), pc(:, :)
        integer :: st, q, c

        call plasma_layout(k, rho, x, layout)
        call spline_space_update_geometry(work_p, layout, st)
        if (st /= 0) error stop 'plasma plane geometry'
        call spline_derham_build(work_p, derham, st, 0)
        if (st /= 0) error stop 'plasma plane derham'
        allocate (v(n0p))
        do c = 1, 2
            current_component = c
            call spline_derham_project(derham, 0, coordinate_p, v, st, 12)
            if (st /= 0) error stop 'plasma node projection'
            nodes_p(c, :, k) = v
        end do
        call shell_layout(k, rho, layout)
        call spline_space_update_geometry(work_f, layout, st)
        if (st /= 0) error stop 'shell plane geometry'
        call spline_quadrature_build(work_f, degree + 4, quad, st)
        if (st /= 0) error stop 'shell plane quadrature'
        allocate (vals(quad%point_count), pc(1, nfull))
        do c = 1, 2
            do q = 1, quad%point_count
                vals(q) = quad%s(c, q)
            end do
            call gs_spline_project(work_f, quad, vals, pc, st)
            if (st /= 0) error stop 'shell node projection'
            nodes_f(c, :, k) = pc(1, :)
        end do
    end subroutine plane_nodes

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
        real(dp) :: theta, xi(2)
        real(dp), allocatable :: rho1(:), v(:)
        logical, allocatable :: in6(:)
        integer :: k, i, j, st, block

        allocate (centre(2, nt), frames(2, 2, nt), wall_s(2, nsample, nt))
        do k = 1, nt
            centre(:, k) = section_axis(plane_phi(k))
            frames(:, :, k) = section_sqrt_hessian(plane_phi(k), centre(:, k))
            do i = 1, nsample
                theta = 2*pi*(i - 1)/nsample
                wall_s(:, i, k) = centre(:, k) + section_boundary_radius(plane_phi(k), &
                    centre(:, k), theta, wall_level)*[cos(theta), sin(theta)]
            end do
        end do
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
        call point_rows(ref_f, quad_fref%block(vol_points), quad_fref%xi(:, vol_points), &
                        rows_vol)
        ! the axis vertex: values and the gradients of one adjacent cell
        call spline_locate(ref_p, 1, centre(:, 1) + [1e-9_dp, 1e-9_dp], block, xi, st)
        if (st /= 0) error stop 'axis vertex'
        call point_rows(ref_p, [block], reshape(xi, [2, 1]), rows_axis)
        ! unknown layout
        ichi = 0
        ig = ichi + n0p*nt
        ix = ig + n0p*nt
        imu = ix + 2*nt
        iphi = imu + nt
        irho = iphi + nphi*nt
        nu = irho + nrho*nt
        allocate (nodes_p(2, n0p, nt), nodes_f(2, nfull, nt))
        print '(a,i6,a,i6,a,i6,a,i6)', 'unknowns', nu, '  (labels per plane', n0p, &
            ', shell potential per plane', nphi, ', boundary', nrho*nt, ')'
        ! the exact state (labels projected on the exact plane meshes)
        allocate (u_exact(nu))
        u_exact = 0
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
    end subroutine setup

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
        call axis_rows(u, r, jt)
        call shell_rows(u, r, jt)
        call balance_rows(u, r, jt)
        do k = 1, nt
            r(ichi + (k - 1)*n0p + wall_p) = 0
        end do
        if (.not. present(jt)) return
        do i = 1, hessian%nrow
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call jt%add(hessian%col(p), i, hessian%val(p))
            end do
        end do
        do i = ig + 1, ix
            call jt%add(i, i, newton_ridge*maxval(abs(hessian%val)))
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
                    detf, tq, kmax
        real(dp), allocatable :: e(:, :)
        integer :: q, m, k, j, i, a, nloc, rowi, rowj, jj
        integer, parameter :: per_plane_points = 2
        integer, allocatable :: cols(:)

        kmax = 0
        do q = 1, size(vol_points)
            a = rows_vol%rp(q)
            nloc = rows_vol%rp(q + 1) - a
            do k = 1, nt
                ! the toroidal quadrature: the planes and the mid-points
                do m = 1, per_plane_points
                    tq = (k - 1 + 0.5_dp*(m - 1))/nt
                    call shell_jet(q, tq, u(iphi + 1:irho), value, grad3, dval, dgrad, s, f, vel, &
                                   lt)
                    detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
                    call section_exact_b([s(1), 2*pi*tq, s(2)], bs)
                    b = bs + grad3
                    w = quad_fref%weight(vol_points(q))*abs(detf)*s(1)*2*pi/(per_plane_points*nt)
                    do j = 1, nt
                        do i = 1, nloc
                            if (shell_slot(rows_vol%ci(a + i - 1)) == 0) cycle
                            rowi = iphi + (j - 1)*nphi + shell_slot(rows_vol%ci(a + i - 1))
                            r(rowi) = r(rowi) + w*dot_product(b, dgrad(:, j, i))
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
        if (present(jt)) then
            do i = iphi + 1, irho
                call jt%add(i, i, 1e-10_dp*kmax)
            end do
        end if
        associate (unused => value)
        end associate
    end subroutine shell_rows

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
                    xdot(2), e(3), db(3), rho2, dqp(3), value
        integer :: k, qq, m, a, i, j, nloc, row, col, b
        integer, allocatable :: pslot(:)

        allocate (pslot(n0p))
        pslot = [(i, i=1, n0p)]
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
                call plane_jet(rows_iface, qq, nodes_f, k, u(iphi + 1:irho), nphi, shell_slot, &
                               value, gc, dvc, dgc, s, f, vel, lt)
                call section_exact_b([s(1), plane_phi(k), s(2)], bs)
                bv = bs + gc
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
        real(dp), allocatable :: r(:), delta(:), history(:), rt(:)
        type(triplet_list) :: jt
        type(csr_t) :: jac_t
        logical, allocatable :: is_fixed(:)
        real(dp) :: alpha
        integer :: it, k, st, j

        allocate (is_fixed(nu))
        is_fixed = .false.
        is_fixed(fixed) = .true.
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
            do j = ix + 1, ix + 2*nt
                call difference_column(u, jt, is_fixed, j, (j - ix + 1)/2)
            end do
            do j = irho + 1, nu
                call difference_column(u, jt, is_fixed, j, (j - irho - 1)/nrho + 1)
            end do
            call jt%to_csr(nu, nu, jac_t)
            if (check_now .and. it == 1) call check_jacobian(u, jac_t)
            allocate (delta(nu))
            call eu_constrained_solve(jac_t, -r, fixed, 0*u(fixed), delta, st)
            if (st /= 0) error stop 'Newton solve'
            print '(a,6es10.2)', '    |d| chi, g, axis, mu, shell, boundary', &
                maxval(abs(delta(ichi + 1:ig))), maxval(abs(delta(ig + 1:ix))), &
                maxval(abs(delta(ix + 1:imu))), maxval(abs(delta(imu + 1:iphi))), &
                maxval(abs(delta(iphi + 1:irho))), maxval(abs(delta(irho + 1:)))
            alpha = min(1.0_dp, 0.05_dp*u(irho + 1)/max(maxval(abs(delta(irho + 1:))), &
                                                         tiny(1.0_dp)), &
                        0.05_dp/max(maxval(abs(delta(ix + 1:imu))), tiny(1.0_dp)))
            do k = 1, 12
                call rebuild(u + alpha*delta)
                call residual(u + alpha*delta, rt)
                rt(fixed) = 0
                if (norm2(rt) < (1 - alpha/4)*history(it)) exit
                alpha = alpha/2
            end do
            print '(a,f7.4)', '    alpha', alpha
            u = u + alpha*delta
            deallocate (delta)
            steps = steps + 1
        end do
        call rebuild(u)
        print '(a,20es9.1)', '  |r| history', history
    end subroutine newton

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
                call shell_jet(q, tq, u(iphi + 1:irho), value, grad3, dval, dgrad, s, f, vel, lt)
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
        ! g at fixed labels
        call rebuild(u)
        call fit_g(u)
        call errors(u, eb, ephi, egamma_start, eaxis)
        print '(a,4es10.2)', 'start: errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma_start, &
            eaxis
        ! stage 1: Gamma_p frozen (the fixed-boundary 3D problem of E4, and
        ! the shell), stage 2: all free
        call newton(u, [fixed_wall, [(irho + i, i=1, nrho*nt)]], steps1, .false.)
        call errors(u, eb, ephi, egamma, eaxis)
        print '(a,i3,a,4es10.2)', 'stage 1 (Gamma_p frozen):', steps1, &
            ' steps; errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma, eaxis
        call newton(u, fixed_wall, steps2, jacobian_check)
        call errors(u, eb, ephi, egamma, eaxis)
        print '(a,i3,a,4es10.2)', 'stage 2 (free boundary):', steps2, &
            ' steps; errors B, grad Phi, Gamma_p, axis', eb, ephi, egamma, eaxis
        call check(steps2 < max_steps, 'Newton converges')
        ! Gamma_p back to the representation error of the exact boundary by
        ! 2M + 1 modes per plane (its error at the exact state)
        call check(egamma < 1.2_dp*egamma0 .and. egamma < 0.6_dp*egamma_start, &
                   'Gamma_p returns to the exact boundary within its Fourier truncation')
        call check(eb < 1.5_dp*eb0, 'B within 50 % of the exact labels with g fitted')
    end subroutine run

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
