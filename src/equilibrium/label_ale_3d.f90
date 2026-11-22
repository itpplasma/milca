module milca_label_ale_3d
    !! 3D label lane on an ALE mesh, one magnetic axis (step A1 of
    !! docs/equilibrium-roadmap.md, the per-plane form of milca_label_ale).
    !!
    !! nt planes phi_k = 2 pi (k - 1)/nt carry O-grids of the same topology,
    !! fitted to the wall curve of each plane, with the core block centred on
    !! the axis point X_k (its central vertex, cells per side even). The
    !! section varies with phi by the trigonometric-nodal interpolation of
    !! the node coordinates (clebsch_3d_system, isoparametric), so the polar
    !! rule of the reference plane serves every plane. Newton on
    !!
    !!   grad W + mu_k [N(X_k); 0] = 0,  grad chi_k(X_k) = 0,  chi_k(X_k) = 0,
    !!
    !! with the core of plane k translating with X_k and the rest of the
    !! section following the harmonic blend w (1 on the core, 0 on the wall):
    !! nodes_k = base_k + w (X_k - anchor_k). The X columns are exact: the
    !! shape derivatives of clebsch_3d_shape_system along the node motion w
    !! and the axis curve. The rows N(X_k) and grad N(X_k) do not change with
    !! X (the core moves rigidly). Between passes the meshes are rebuilt with
    !! their cores on X_k, chi and g transferred and the winding frames
    !! (H^(1/2) of chi_h at X_k) recomputed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_quadrature_point_singular, spline_point_matrix, spline_locate
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_boundary_dofs, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_profiles, clebsch_exact_field, &
        clebsch_3d_system, clebsch_3d_shape_system
    use milca_spline_ale, only: spline_geometry_dofs, spline_set_geometry, &
        spline_blend_weights, spline_block_dofs, spline_transfer_v0
    implicit none
    private
    public :: label3d_t, label3d_init, label3d_initial_labels, label3d_solve_g, &
              label3d_newton, label3d_solve, label3d_field_error, label3d_surface, &
              label3d_energy, label3d_frames, label3d_section_at, trig_nodal_weights, label3d_set_walls

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: label3d_t
        integer :: nt = 0, cells = 2, degree = 2, extra_points = 8
        real(dp), allocatable :: samples(:, :, :) ! wall curve per plane (2, m, nt)
        type(spline_space_t) :: space ! reference plane (plane 1)
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad, aquad
        real(dp), allocatable :: base(:, :, :) ! node coordinates (2, n0, nt)
        real(dp), allocatable :: anchor(:, :) ! core centres (2, nt)
        real(dp), allocatable :: w(:) ! core blend
        real(dp), allocatable :: frames(:, :, :) ! winding frames (2, 2, nt)
        real(dp), allocatable :: nrow(:), ngrad(:, :, :)
        integer, allocatable :: wall(:)
        real(dp), allocatable :: chi(:), g(:), x(:, :), mu(:)
        real(dp) :: cell_size = 1, ridge = 1e-12_dp, newton_ridge = 1e-8_dp
        ! Levenberg shift of the chi block (relative to max |H|), adapted:
        ! low magnetic shear leaves near-null label modes (m iota ~ n)
        real(dp) :: levenberg = 0
        logical :: use_levenberg = .false. ! opt-in, capped at 1e-4
        ! gauge g -> g + k(chi) removed by constraints (plan, section 4):
        ! sum_i g_i P_k(chi_i) = 0, k < gauge_k, Legendre P_k on [0, 1];
        ! multipliers in mu(nt + 1:nt + gauge_k). 0: the ridge instead.
        integer :: gauge_k = 0
        integer :: passes = 0, newton_total = 0
        real(dp), allocatable :: history(:)
        logical :: verbose = .false.
    end type

contains

    subroutine label3d_init(s, samples, centre, cells, status)
        !! Meshes of all planes from the wall curves samples(2, m, nt)
        !! (counter-clockwise, first sample on the outboard midplane) with the
        !! cores centred on centre(2, nt); cells per block side (even).
        type(label3d_t), intent(inout) :: s
        real(dp), intent(in) :: samples(:, :, :), centre(:, :)
        integer, intent(in) :: cells
        integer, intent(out) :: status
        integer, allocatable :: core_cells(:), core(:)
        integer :: k
        status = 1
        if (mod(cells, 2) /= 0) return
        s%nt = size(samples, 3)
        s%cells = cells
        s%samples = samples
        call build_meshes(s, centre, status)
        if (status /= 0) return
        call spline_block_dofs(s%space, 1, core_cells, core)
        call spline_blend_weights(s%space, core_cells, wall1_of(s), s%w, status)
        if (status /= 0) return
        if (s%verbose) print '(a,es10.2)', '  max |core vertex - centre|', &
            maxval(norm2(s%anchor - centre, 1))
        s%x = s%anchor
        if (allocated(s%mu)) deallocate (s%mu)
        allocate (s%mu(s%nt + s%gauge_k), source=0.0_dp)
        if (.not. allocated(s%frames)) then
            allocate (s%frames(2, 2, s%nt), source=0.0_dp)
            do k = 1, s%nt
                s%frames(1, 1, k) = 1
                s%frames(2, 2, k) = 1
            end do
        end if
    end subroutine label3d_init

    function wall1_of(s) result(wall1)
        type(label3d_t), intent(in) :: s
        integer, allocatable :: wall1(:)
        call spline_derham_boundary_dofs(s%derham, LAYOUT_WALL, 0, wall1)
    end function wall1_of

    subroutine plane_layout(s, k, centre, layout, status)
        type(label3d_t), intent(in) :: s
        integer, intent(in) :: k
        real(dp), intent(in) :: centre(2)
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        associate (sm => s%samples(:, :, k))
            call layout_ogrid(centre, 0.5_dp*(maxval(sm(1, :)) - minval(sm(1, :))), 0.5_dp, &
                              layout, status)
            if (status /= 0) return
            call layout_set_curve(layout, LAYOUT_WALL, sm, status, centre_core=.true.)
        end associate
    end subroutine plane_layout

    subroutine build_meshes(s, centre, status)
        !! base nodes of every plane with the core on centre(:, k); derham,
        !! quadratures and rows of the reference plane 1.
        type(label3d_t), intent(inout) :: s
        real(dp), intent(in) :: centre(:, :)
        integer, intent(out) :: status
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        real(dp), allocatable :: xg(:, :)
        integer, allocatable :: wall1(:)
        integer :: k, n0
        do k = s%nt, 1, -1
            call plane_layout(s, k, centre(:, k), layout, status)
            if (status /= 0) return
            if (k == s%nt) then
                options%degree = s%degree
                options%regularity = SPLINE_C0
                allocate (options%knots(2, layout%block_count))
                options%knots(:, :)%element_count = s%cells
                call spline_space_build(layout, options, s%space, status)
            else
                call spline_space_update_geometry(s%space, layout, status)
            end if
            if (status /= 0) return
            call spline_geometry_dofs(s%space, xg, status)
            if (status /= 0) return
            n0 = size(xg, 2)
            if (.not. allocated(s%base)) allocate (s%base(2, n0, s%nt))
            if (size(s%base, 2) /= n0) then
                deallocate (s%base)
                allocate (s%base(2, n0, s%nt))
            end if
            s%base(:, :, k) = xg
            call plane_rows(s, k)
        end do
        ! space now holds plane 1; the anchors are the central vertices
        ! (layout_set_curve keeps the core near centre, not exactly on it)
        if (.not. allocated(s%anchor)) allocate (s%anchor(2, s%nt))
        do k = 1, s%nt
            s%anchor(:, k) = matmul(s%base(:, :, k), s%nrow)
        end do
        call spline_derham_build(s%space, s%derham, status, 0)
        if (status /= 0) return
        call spline_quadrature_build(s%space, s%degree + s%extra_points, s%quad, status)
        if (status /= 0) return
        call spline_quadrature_point_singular(s%space, s%quad, s%anchor(:, 1), s%aquad, status)
        if (status /= 0) return
        call spline_derham_boundary_dofs(s%derham, LAYOUT_WALL, 0, wall1)
        s%wall = [((k - 1)*n0 + wall1, k=1, s%nt)]
        s%cell_size = (maxval(s%base(1, :, 1)) - minval(s%base(1, :, 1)))/(3*s%cells)
    end subroutine build_meshes

    subroutine plane_rows(s, k)
        !! N(X_k) (the central vertex function) and grad N(X_k), the mean of
        !! the one-sided gradients of the four cells at the vertex, on the
        !! geometry now in s%space (plane k).
        type(label3d_t), intent(inout) :: s
        integer, intent(in) :: k
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        real(dp), parameter :: e = 1e-9_dp
        real(dp) :: xi(2, 4)
        integer :: st, c, n0, q
        n0 = s%space%dof_count
        if (.not. allocated(s%ngrad)) allocate (s%ngrad(n0, 2, s%nt))
        if (size(s%ngrad, 1) /= n0) then
            deallocate (s%ngrad)
            allocate (s%ngrad(n0, 2, s%nt))
        end if
        xi = reshape([0.5_dp + e, 0.5_dp + e, 0.5_dp - e, 0.5_dp + e, 0.5_dp - e, 0.5_dp - e, &
                      0.5_dp + e, 0.5_dp - e], [2, 4])
        if (allocated(s%nrow)) deallocate (s%nrow)
        allocate (s%nrow(n0), source=0.0_dp)
        call spline_point_matrix(s%space, [1], reshape([0.5_dp, 0.5_dp], [2, 1]), 0, rp, ci, &
                                 v, st)
        if (st /= 0) error stop 'label3d: vertex row'
        s%nrow(ci) = v
        where (abs(s%nrow) < 1e-12_dp) s%nrow = 0
        s%ngrad(:, :, k) = 0
        do c = 1, 2
            call spline_point_matrix(s%space, [1, 1, 1, 1], xi, c, rp, ci, v, st)
            if (st /= 0) error stop 'label3d: vertex gradient rows'
            do q = 1, 4
                s%ngrad(ci(rp(q):rp(q + 1) - 1), c, k) = s%ngrad(ci(rp(q):rp(q + 1) - 1), c, k) + &
                                                          0.25_dp*v(rp(q):rp(q + 1) - 1)
            end do
        end do
    end subroutine plane_rows

    subroutine set_plane(s, k, nodes)
        !! The reference space with the node coordinates of plane k.
        type(label3d_t), intent(inout) :: s
        integer, intent(in) :: k
        real(dp), intent(in) :: nodes(:, :, :)
        call spline_set_geometry(s%space, nodes(:, :, k))
    end subroutine set_plane

    subroutine label3d_initial_labels(s)
        !! chi = (|p - X_k| / rho_wall(angle))^2 at the node points p of each
        !! plane, rho_wall the distance from X_k to the wall curve along the
        !! ray: nested ellipse-like labels, 0 at the axis and 1 on the wall.
        type(label3d_t), intent(inout) :: s
        integer :: k, i, n0
        real(dp) :: d(2), r, rw
        n0 = size(s%base, 2)
        if (allocated(s%chi)) deallocate (s%chi)
        allocate (s%chi(n0*s%nt))
        do k = 1, s%nt
            do i = 1, n0
                d = s%base(:, i, k) - s%x(:, k)
                r = norm2(d)
                if (r < 1e-14_dp) then
                    s%chi((k - 1)*n0 + i) = 0
                    cycle
                end if
                rw = ray_distance(s%samples(:, :, k), s%x(:, k), d/r)
                s%chi((k - 1)*n0 + i) = min(1.0_dp, (r/rw)**2)
            end do
        end do
        s%chi(s%wall) = 1
    end subroutine label3d_initial_labels

    real(dp) function ray_distance(curve, origin, dir) result(t)
        !! Distance from origin along dir to the closed polygon curve.
        real(dp), intent(in) :: curve(:, :), origin(2), dir(2)
        real(dp) :: a(2), b(2), e(2), den, u, tt
        integer :: i, m
        m = size(curve, 2)
        t = huge(1.0_dp)
        do i = 1, m
            a = curve(:, i) - origin
            b = curve(:, mod(i, m) + 1) - origin
            e = b - a
            den = dir(1)*e(2) - dir(2)*e(1)
            if (abs(den) < tiny(1.0_dp)) cycle
            ! origin + tt dir = a + u e
            tt = (a(1)*e(2) - a(2)*e(1))/den
            u = (a(1)*dir(2) - a(2)*dir(1))/den
            if (u >= -1e-12_dp .and. u <= 1 + 1e-12_dp .and. tt > 0) t = min(t, tt)
        end do
    end function ray_distance

    subroutine label3d_frames(s)
        !! Winding frames H^(1/2) of chi_h at X_k (mean of the second
        !! derivative rows of the four cells at the vertex), identity where H
        !! is not positive definite.
        type(label3d_t), intent(inout) :: s
        real(dp), allocatable :: nodes(:, :, :)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        real(dp), parameter :: e = 1e-9_dp
        real(dp) :: xi(2, 4), h(3), hm(2, 2), tr, dt, sq
        integer :: k, c, q, st, n0
        n0 = size(s%base, 2)
        nodes = node_coordinates(s, s%x)
        xi = reshape([0.5_dp + e, 0.5_dp + e, 0.5_dp - e, 0.5_dp + e, 0.5_dp - e, 0.5_dp - e, &
                      0.5_dp + e, 0.5_dp - e], [2, 4])
        do k = 1, s%nt
            call set_plane(s, k, nodes)
            do c = 3, 5
                call spline_point_matrix(s%space, [1, 1, 1, 1], xi, c, rp, ci, v, st)
                h(c - 2) = 0
                do q = 1, 4
                    h(c - 2) = h(c - 2) + 0.25_dp*dot_product(v(rp(q):rp(q + 1) - 1), &
                                                             s%chi((k - 1)*n0 + ci(rp(q):rp(q + 1) - 1)))
                end do
            end do
            hm = reshape([h(1), h(2), h(2), h(3)], [2, 2])
            tr = hm(1, 1) + hm(2, 2)
            dt = hm(1, 1)*hm(2, 2) - hm(1, 2)**2
            if (dt > 0 .and. tr > 0) then
                ! sqrt of a 2x2 SPD matrix: (H + sqrt(det) I)/sqrt(tr + 2 sqrt(det))
                sq = sqrt(dt)
                s%frames(:, :, k) = (hm + sq*reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))/ &
                                    sqrt(tr + 2*sq)
            else
                s%frames(:, :, k) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
            end if
        end do
        call set_plane(s, 1, s%base)
    end subroutine label3d_frames

    function node_coordinates(s, x) result(nodes)
        type(label3d_t), intent(in) :: s
        real(dp), intent(in) :: x(:, :)
        real(dp) :: nodes(2, size(s%w), s%nt)
        integer :: k
        do k = 1, s%nt
            nodes(1, :, k) = s%base(1, :, k) + s%w*(x(1, k) - s%anchor(1, k))
            nodes(2, :, k) = s%base(2, :, k) + s%w*(x(2, k) - s%anchor(2, k))
        end do
    end function node_coordinates

    subroutine label3d_solve_g(s, profiles, status)
        !! g minimising W at fixed labels (W is quadratic and convex in g).
        type(label3d_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        integer, intent(out) :: status
        real(dp), allocatable :: grad(:), delta(:)
        type(csr_t) :: hessian
        real(dp) :: energy
        integer :: ntot, i
        ntot = size(s%chi)
        if (.not. allocated(s%g)) allocate (s%g(ntot), source=0.0_dp)
        call clebsch_3d_system(s%derham, s%aquad, s%nt, profiles, s%x(:, 1), s%chi, s%g, energy, &
                               grad, hessian, s%x, node_coordinates(s, s%x), s%frames)
        if (s%gauge_k == 0) then
            call add_ridge(hessian, ntot, s%ridge)
            allocate (delta(2*ntot))
            call eu_constrained_solve(hessian, -grad, [(i, i=1, ntot)], 0*s%chi, delta, status)
            if (status == 0) s%g = s%g + delta(ntot + 1:)
            return
        end if
        ! [H C; C^T 0] with the gauge rows C(i, k) = P_k(chi_i) on g
        block
            integer, allocatable :: rows(:), cols(:)
            real(dp), allocatable :: vals(:), rhs(:)
            real(dp) :: pk(s%gauge_k), dpk(s%gauge_k)
            type(csc_t) :: m
            type(csr_t) :: a
            type(fortsparse_status_t) :: info
            integer :: nz, p, kk, n
            n = 2*ntot + s%gauge_k
            allocate (rows(size(hessian%val) + 2*ntot*s%gauge_k + s%gauge_k), &
                      cols(size(hessian%val) + 2*ntot*s%gauge_k + s%gauge_k), &
                      vals(size(hessian%val) + 2*ntot*s%gauge_k + s%gauge_k))
            nz = 0
            do i = 1, 2*ntot
                do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                    nz = nz + 1
                    rows(nz) = i; cols(nz) = hessian%col(p); vals(nz) = hessian%val(p)
                end do
            end do
            allocate (rhs(n), source=0.0_dp)
            rhs(1:2*ntot) = -grad
            do i = 1, ntot
                call legendre01(s%chi(i), s%gauge_k, pk, dpk)
                do kk = 1, s%gauge_k
                    nz = nz + 1
                    rows(nz) = ntot + i; cols(nz) = 2*ntot + kk; vals(nz) = pk(kk)
                    nz = nz + 1
                    rows(nz) = 2*ntot + kk; cols(nz) = ntot + i; vals(nz) = pk(kk)
                    rhs(2*ntot + kk) = rhs(2*ntot + kk) - pk(kk)*s%g(i)
                end do
            end do
            do kk = 1, s%gauge_k
                nz = nz + 1
                rows(nz) = 2*ntot + kk; cols(nz) = 2*ntot + kk; vals(nz) = 0
            end do
            call csc_from_triplet(n, n, cols(:nz), rows(:nz), vals(:nz), m, info)
            if (info%code /= 0) error stop 'label3d: g gauge triplets'
            a%nrow = n
            a%ncol = n
            a%row_ptr = m%col_ptr
            a%col = m%row_idx
            a%val = m%val
            allocate (delta(n))
            call eu_constrained_solve(a, rhs, [(i, i=1, ntot)], 0*s%chi, delta, status)
            if (status == 0) s%g = s%g + delta(ntot + 1:2*ntot)
        end block
    end subroutine label3d_solve_g

    real(dp) function label3d_energy(s, profiles) result(energy)
        type(label3d_t), intent(in) :: s
        procedure(clebsch_profiles) :: profiles
        real(dp), allocatable :: grad(:)
        call clebsch_3d_system(s%derham, s%aquad, s%nt, profiles, s%x(:, 1), s%chi, s%g, energy, &
                               grad, plane_axis=s%x, plane_nodes=node_coordinates(s, s%x), &
                               plane_frames=s%frames)
    end function label3d_energy

    subroutine residual(s, profiles, chi, g, x, mu, r, hessian, dgdx)
        !! Newton residual; with hessian and dgdx also the (chi, g) Hessian
        !! and the exact X columns.
        type(label3d_t), intent(in) :: s
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: chi(:), g(:), x(:, :), mu(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), allocatable, intent(out), optional :: dgdx(:, :)
        real(dp), allocatable :: grad(:), dir(:)
        type(csr_t) :: shape
        real(dp) :: energy
        integer :: n0, ntot, nf, k, c, i, j, p
        n0 = size(s%w)
        ntot = n0*s%nt
        nf = 2*ntot
        if (present(dgdx)) then
            call clebsch_3d_shape_system(s%derham, s%aquad, s%nt, profiles, chi, g, x, &
                                         node_coordinates(s, x), s%frames, energy, grad, &
                                         hessian, shape)
            allocate (dgdx(nf + 3*s%nt, 2*s%nt), source=0.0_dp)
            allocate (dir(shape%ncol))
            do k = 1, s%nt
                do c = 1, 2
                    j = 2*(k - 1) + c
                    dir = 0
                    do i = 1, n0
                        dir(((k - 1)*n0 + i - 1)*2 + c) = s%w(i)
                    end do
                    dir(2*n0*s%nt + j) = 1
                    do i = 1, nf
                        do p = shape%row_ptr(i), shape%row_ptr(i + 1) - 1
                            dgdx(i, j) = dgdx(i, j) + shape%val(p)*dir(shape%col(p))
                        end do
                    end do
                end do
            end do
            dgdx(s%wall, :) = 0
        else
            call clebsch_3d_system(s%derham, s%aquad, s%nt, profiles, x(:, 1), chi, g, energy, &
                                   grad, plane_axis=x, plane_nodes=node_coordinates(s, x), &
                                   plane_frames=s%frames)
        end if
        allocate (r(nf + 3*s%nt + s%gauge_k))
        r(1:nf) = grad
        do k = 1, s%nt
            associate (chik => chi((k - 1)*n0 + 1:k*n0))
                r((k - 1)*n0 + 1:k*n0) = r((k - 1)*n0 + 1:k*n0) + mu(k)*s%nrow
                do c = 1, 2
                    r(nf + 2*(k - 1) + c) = dot_product(s%ngrad(:, c, k), chik)
                end do
                r(nf + 2*s%nt + k) = dot_product(s%nrow, chik)
            end associate
        end do
        if (s%gauge_k > 0) then
            block
                real(dp) :: pk(s%gauge_k), dpk(s%gauge_k)
                do i = 1, ntot
                    call legendre01(chi(i), s%gauge_k, pk, dpk)
                    r(ntot + i) = r(ntot + i) + dot_product(pk, mu(s%nt + 1:))
                    r(nf + 3*s%nt + 1:) = r(nf + 3*s%nt + 1:) + pk*g(i)
                end do
            end block
        end if
        r(s%wall) = 0
    end subroutine residual

    pure subroutine legendre01(x, k, p, dp_)
        !! Legendre polynomials P_0..P_(k-1) of 2x - 1 and their x derivatives.
        real(dp), intent(in) :: x
        integer, intent(in) :: k
        real(dp), intent(out) :: p(k), dp_(k)
        real(dp) :: t
        integer :: j
        t = 2*x - 1
        p(1) = 1
        dp_(1) = 0
        if (k >= 2) then
            p(2) = t
            dp_(2) = 2
        end if
        do j = 2, k - 1
            ! (j) P_j = (2j - 1) t P_(j-1) - (j - 1) P_(j-2), index shifted by 1
            p(j + 1) = ((2*j - 1)*t*p(j) - (j - 1)*p(j - 1))/j
            dp_(j + 1) = ((2*j - 1)*(2*p(j) + t*dp_(j)) - (j - 1)*dp_(j - 1))/j
        end do
    end subroutine legendre01

    subroutine label3d_newton(s, profiles, tolerance, max_steps, steps, status)
        !! Newton on (chi, g, X_k, mu_k) with the cores following X_k. Status
        !! 0 converged, 1 linear solve, 2 not converged.
        type(label3d_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: steps, status
        real(dp), allocatable :: r(:), rt(:), delta(:), dgdx(:, :)
        type(csr_t) :: hessian, jac_t
        real(dp) :: alpha, r0, rn, dxmax
        integer :: it, k, st, n0, ntot, nf, nt
        n0 = size(s%w)
        nt = s%nt
        ntot = n0*nt
        nf = 2*ntot
        if (.not. allocated(s%history)) allocate (s%history(0))
        steps = 0
        status = 2
        r0 = -1
        do it = 1, max_steps + 1
            if (it <= max_steps) then
                call residual(s, profiles, s%chi, s%g, s%x, s%mu, r, hessian, dgdx)
            else
                call residual(s, profiles, s%chi, s%g, s%x, s%mu, r)
            end if
            rn = norm2(r)
            s%history = [s%history, rn]
            if (r0 < 0) r0 = rn
            if (s%verbose) print '(a,i3,a,4es10.2)', '    newton', it - 1, '  |r| chi g X mu', &
                norm2(r(1:ntot)), norm2(r(ntot + 1:nf)), norm2(r(nf + 1:nf + 2*nt)), &
                norm2(r(nf + 2*nt + 1:))
            if (it > 1 .and. rn < tolerance*r0) then
                status = 0
                exit
            end if
            if (it > max_steps) exit
            if (s%gauge_k == 0) call add_ridge(hessian, ntot, s%newton_ridge)
            if (s%levenberg > 0) call shift_chi(hessian, ntot, s%levenberg)
            call newton_matrix(s, hessian, dgdx, jac_t)
            if (it == 1 .and. jacobian_check_on()) call jacobian_check(s, profiles, jac_t)
            allocate (delta(nf + 3*nt + s%gauge_k))
            call eu_constrained_solve(jac_t, -r, s%wall, 0*s%chi(s%wall), delta, st)
            if (st /= 0) then
                status = 1
                return
            end if
            ! X moves at most a quarter core cell per step
            dxmax = maxval(abs(delta(nf + 1:nf + 2*nt)))
            alpha = min(1.0_dp, 0.25_dp*s%cell_size/max(dxmax, tiny(1.0_dp)))
            do k = 1, 10
                call residual(s, profiles, s%chi + alpha*delta(1:ntot), &
                              s%g + alpha*delta(ntot + 1:nf), &
                              s%x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt]), &
                              s%mu + alpha*delta(nf + 2*nt + 1:), rt)
                if (norm2(rt) < (1 - alpha/4)*rn) exit
                alpha = alpha/2
            end do
            ! adapt the shift: full steps relax it, cut steps raise it
            if (s%use_levenberg) then
                if (alpha >= 0.99_dp) then
                    s%levenberg = s%levenberg/10
                    if (s%levenberg < 1e-12_dp) s%levenberg = 0
                else if (alpha < 0.1_dp) then
                    s%levenberg = min(max(10*s%levenberg, 1e-8_dp), 1e-4_dp)
                end if
            end if
            if (s%verbose) print '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '      alpha', alpha, &
                '  |dX|', dxmax, '  |dchi|', maxval(abs(delta(1:ntot))), '  trial |r|', norm2(rt)
            s%chi = s%chi + alpha*delta(1:ntot)
            s%g = s%g + alpha*delta(ntot + 1:nf)
            s%x = s%x + alpha*reshape(delta(nf + 1:nf + 2*nt), [2, nt])
            s%mu = s%mu + alpha*delta(nf + 2*nt + 1:)
            deallocate (delta)
            steps = steps + 1
        end do
        s%newton_total = s%newton_total + steps
    end subroutine label3d_newton

    subroutine label3d_solve(s, profiles, tolerance, mesh_tolerance, max_passes, max_steps, &
                             status)
        !! Passes of Newton and remeshing (cores rebuilt on X_k, chi and g
        !! transferred, frames recomputed, g re-solved) until Newton moves no
        !! X_k by more than mesh_tolerance core cells. Status 0 converged,
        !! 1 Newton, 2 mesh, 3 transfer, 4 passes exhausted.
        type(label3d_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: tolerance, mesh_tolerance
        integer, intent(in) :: max_passes, max_steps
        integer, intent(out) :: status
        integer :: pass, steps, st
        status = 4
        do pass = 1, max_passes
            s%passes = pass
            call label3d_newton(s, profiles, tolerance, max_steps, steps, st)
            if (st == 1) then
                status = 1
                return
            end if
            if (s%verbose) print '(a,i2,a,i3,a,es10.2,a,es10.2)', '  pass', pass, '  Newton steps', &
                steps, '  max |X - core|', maxval(norm2(s%x - s%anchor, 1)), '  cell', s%cell_size
            if (st == 0 .and. maxval(norm2(s%x - s%anchor, 1)) < mesh_tolerance*s%cell_size) then
                status = 0
                return
            end if
            if (pass == max_passes) return
            call remesh(s, profiles, st)
            if (st /= 0 .and. s%verbose) print '(a,i0)', '  remesh failed: ', st
            if (st /= 0) then
                status = 2
                return
            end if
        end do
    end subroutine label3d_solve

    subroutine remesh(s, profiles, status)
        !! Cores rebuilt on X_k, chi and g transferred plane by plane.
        type(label3d_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        integer, intent(out) :: status
        real(dp), allocatable :: old_nodes(:, :, :), chi(:), g(:)
        type(spline_derham_t) :: from, to
        integer :: k, n0
        n0 = size(s%w)
        old_nodes = node_coordinates(s, s%x)
        call build_meshes(s, s%x, status)
        if (status /= 0) then
            if (s%verbose) print '(a,i0)', '  build_meshes failed: ', status
            return
        end if
        allocate (chi(size(s%chi)), g(size(s%g)))
        do k = 1, s%nt
            call spline_set_geometry(s%space, old_nodes(:, :, k))
            call spline_derham_build(s%space, from, status, 0)
            if (status /= 0) return
            call spline_set_geometry(s%space, s%base(:, :, k))
            call spline_derham_build(s%space, to, status, 0)
            if (status /= 0) return
            call spline_transfer_v0(from, s%chi((k - 1)*n0 + 1:k*n0), to, &
                                    chi((k - 1)*n0 + 1:k*n0), status)
            if (status == 0) call spline_transfer_v0(from, s%g((k - 1)*n0 + 1:k*n0), to, &
                                                     g((k - 1)*n0 + 1:k*n0), status)
            if (status /= 0) then
                status = 3
                return
            end if
        end do
        call spline_set_geometry(s%space, s%base(:, :, 1))
        chi(s%wall) = 1
        s%chi = chi
        s%g = g
        call label3d_frames(s)
        call label3d_solve_g(s, profiles, status)
    end subroutine remesh

    real(dp) function label3d_field_error(s, profiles, exact) result(rel)
        !! Relative L2 error of B against the exact field.
        type(label3d_t), intent(in) :: s
        procedure(clebsch_profiles) :: profiles
        procedure(clebsch_exact_field) :: exact
        real(dp), allocatable :: grad(:)
        real(dp) :: energy, ferr(2)
        call clebsch_3d_system(s%derham, s%aquad, s%nt, profiles, s%x(:, 1), s%chi, s%g, energy, &
                               grad, plane_axis=s%x, plane_nodes=node_coordinates(s, s%x), &
                               plane_frames=s%frames, exact_field=exact, field_error=ferr)
        rel = sqrt(ferr(1)/ferr(2))
    end function label3d_field_error

    subroutine label3d_surface(s, k, level, npts, rz, status)
        !! Points of the level set chi_h = level on plane k, on npts rays
        !! from X_k (bisection on each ray), counter-clockwise from the
        !! outboard midplane.
        type(label3d_t), intent(inout) :: s
        integer, intent(in) :: k, npts
        real(dp), intent(in) :: level
        real(dp), intent(out) :: rz(2, npts)
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:, :, :)
        type(spline_derham_t) :: d
        real(dp) :: dir(2), lo, hi, mid, rw, ang
        integer :: i, it, n0
        n0 = size(s%w)
        nodes = node_coordinates(s, s%x)
        call spline_set_geometry(s%space, nodes(:, :, k))
        call spline_derham_build(s%space, d, status, 0)
        if (status /= 0) return
        do i = 1, npts
            ang = 2*pi*(i - 1)/npts
            dir = [cos(ang), sin(ang)]
            rw = ray_distance(s%samples(:, :, k), s%x(:, k), dir)
            lo = 0
            hi = rw*(1 - 1e-9_dp)
            if (level >= 1) then
                rz(:, i) = s%x(:, k) + hi*dir
                cycle
            end if
            do it = 1, 60
                mid = 0.5_dp*(lo + hi)
                if (value_at(s%x(:, k) + mid*dir) < level) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            rz(:, i) = s%x(:, k) + 0.5_dp*(lo + hi)*dir
        end do
        call spline_set_geometry(s%space, s%base(:, :, 1))
        status = 0
    contains
        real(dp) function value_at(p) result(v)
            real(dp), intent(in) :: p(2)
            real(dp) :: xi(2), val(1)
            integer :: b, st
            call spline_locate(d%space, 1, p, b, xi, st)
            if (st /= 0) then
                v = 1
                return
            end if
            call spline_derham_evaluate(d, 0, s%chi((k - 1)*n0 + 1:k*n0), b, xi, 0.0_dp, val)
            v = val(1)
        end function value_at
    end subroutine label3d_surface

    subroutine label3d_set_walls(s, samples, status)
        !! New wall curves (continuation in the boundary shape): the meshes
        !! are rebuilt with the cores on X_k, chi and g kept as coefficients
        !! of the same topology (labels in mesh coordinates), the frames
        !! recomputed and g re-solved by the caller.
        type(label3d_t), intent(inout) :: s
        real(dp), intent(in) :: samples(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: x(:, :)
        x = s%x
        s%samples = samples
        call build_meshes(s, x, status)
        if (status /= 0) then
            if (s%verbose) print '(a,i0)', '  set_walls: meshes with the cores on X failed: ', &
                status
            return
        end if
        s%x = s%anchor
    end subroutine label3d_set_walls

    function trig_nodal_weights(nt, t) result(l)
        !! Trigonometric nodal basis L_k(t) of nt (odd) equally spaced nodes
        !! t_k = (k - 1)/nt on the period 1.
        integer, intent(in) :: nt
        real(dp), intent(in) :: t
        real(dp) :: l(nt), d, den
        integer :: k
        do k = 1, nt
            d = t - (k - 1)/real(nt, dp)
            den = nt*sin(pi*d)
            if (abs(den) < 1e-13_dp) then
                l(k) = 1
            else
                l(k) = sin(nt*pi*d)/den
            end if
        end do
        l = l/sum(l)
    end function trig_nodal_weights

    subroutine label3d_section_at(s, t, levels, npts, axis, rz, status)
        !! The axis and the level sets chi_h = levels(j) on the plane
        !! phi = 2 pi t (interpolated nodes, labels and axis), on npts rays
        !! from the axis counter-clockwise from the outboard midplane.
        type(label3d_t), intent(inout) :: s
        real(dp), intent(in) :: t, levels(:)
        integer, intent(in) :: npts
        real(dp), intent(out) :: axis(2), rz(2, npts, size(levels))
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:, :, :), xt(:, :), ct(:)
        real(dp) :: l(s%nt), dir(2), lo, hi, mid, ang, edge
        type(spline_derham_t) :: d
        integer :: i, j, k, it, n0
        n0 = size(s%w)
        l = trig_nodal_weights(s%nt, t)
        nodes = node_coordinates(s, s%x)
        allocate (xt(2, n0), ct(n0), source=0.0_dp)
        axis = 0
        do k = 1, s%nt
            xt = xt + l(k)*nodes(:, :, k)
            ct = ct + l(k)*s%chi((k - 1)*n0 + 1:k*n0)
            axis = axis + l(k)*s%x(:, k)
        end do
        call spline_set_geometry(s%space, xt)
        call spline_derham_build(s%space, d, status, 0)
        if (status /= 0) return
        do i = 1, npts
            ang = 2*pi*(i - 1)/npts
            dir = [cos(ang), sin(ang)]
            ! the wall along the ray: last point inside the mesh
            lo = 0
            hi = 10*maxval(abs(xt(1, :) - axis(1)))
            do it = 1, 60
                mid = 0.5_dp*(lo + hi)
                if (inside(axis + mid*dir)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            edge = lo
            do j = 1, size(levels)
                if (levels(j) >= 1) then
                    rz(:, i, j) = axis + edge*dir
                    cycle
                end if
                lo = 0
                hi = edge
                do it = 1, 60
                    mid = 0.5_dp*(lo + hi)
                    if (value_at(axis + mid*dir) < levels(j)) then
                        lo = mid
                    else
                        hi = mid
                    end if
                end do
                rz(:, i, j) = axis + 0.5_dp*(lo + hi)*dir
            end do
        end do
        call spline_set_geometry(s%space, s%base(:, :, 1))
        status = 0
    contains
        logical function inside(p)
            real(dp), intent(in) :: p(2)
            real(dp) :: xi(2)
            integer :: b, st
            call spline_locate(d%space, 1, p, b, xi, st)
            inside = st == 0
        end function inside
        real(dp) function value_at(p) result(v)
            real(dp), intent(in) :: p(2)
            real(dp) :: xi(2), val(1)
            integer :: b, st
            call spline_locate(d%space, 1, p, b, xi, st)
            if (st /= 0) then
                v = 1
                return
            end if
            call spline_derham_evaluate(d, 0, ct, b, xi, 0.0_dp, val)
            v = val(1)
        end function value_at
    end subroutine label3d_section_at

    logical function jacobian_check_on()
        character(8) :: v
        call get_environment_variable('MILCA_LABEL3D_JACOBIAN', v)
        jacobian_check_on = len_trim(v) > 0
    end function jacobian_check_on

    subroutine jacobian_check(s, profiles, jac_t)
        !! J v against central differences of the residual, per block of v
        !! (chi, g, X, mu); wall entries of v and rows excluded.
        type(label3d_t), intent(in) :: s
        procedure(clebsch_profiles) :: profiles
        type(csr_t), intent(in) :: jac_t
        real(dp), allocatable :: v(:), jv(:), rp(:), rm(:), fd(:)
        real(dp) :: h
        integer :: n0, nt, ntot, nf, blk, lo, hi, j, p
        n0 = size(s%w)
        nt = s%nt
        ntot = n0*nt
        nf = 2*ntot
        allocate (v(nf + 3*nt + s%gauge_k), jv(nf + 3*nt + s%gauge_k))
        do blk = 1, 4
            select case (blk)
            case (1)
                lo = 1; hi = ntot
            case (2)
                lo = ntot + 1; hi = nf
            case (3)
                lo = nf + 1; hi = nf + 2*nt
            case default
                lo = nf + 2*nt + 1; hi = nf + 3*nt + s%gauge_k
            end select
            v = 0
            do j = lo, hi
                v(j) = sin(1.3_dp*j) + 0.5_dp
            end do
            v(s%wall) = 0
            h = merge(1e-6_dp, 1e-5_dp, blk == 3)/maxval(abs(v))
            jv = 0
            do j = 1, jac_t%nrow
                do p = jac_t%row_ptr(j), jac_t%row_ptr(j + 1) - 1
                    jv(jac_t%col(p)) = jv(jac_t%col(p)) + jac_t%val(p)*v(j)
                end do
            end do
            call residual(s, profiles, s%chi + h*v(1:ntot), s%g + h*v(ntot + 1:nf), &
                          s%x + h*reshape(v(nf + 1:nf + 2*nt), [2, nt]), s%mu + h*v(nf + 2*nt + 1:), rp)
            call residual(s, profiles, s%chi - h*v(1:ntot), s%g - h*v(ntot + 1:nf), &
                          s%x - h*reshape(v(nf + 1:nf + 2*nt), [2, nt]), s%mu - h*v(nf + 2*nt + 1:), rm)
            fd = (rp - rm)/(2*h)
            jv(s%wall) = 0
            print '(a,i2,a,es10.2,a,es10.2,a,3es10.2)', '  Jacobian block', blk, '  |Jv-FD|/|FD|', &
                norm2(jv - fd)/norm2(fd), '  |FD|', norm2(fd), '  rows chi g rest', &
                norm2(jv(1:ntot) - fd(1:ntot))/max(norm2(fd(1:ntot)), tiny(1.0_dp)), &
                norm2(jv(ntot + 1:nf) - fd(ntot + 1:nf))/max(norm2(fd(ntot + 1:nf)), tiny(1.0_dp)), &
                norm2(jv(nf + 1:) - fd(nf + 1:))/max(norm2(fd(nf + 1:)), tiny(1.0_dp))
        end do
    end subroutine jacobian_check

    subroutine newton_matrix(s, hessian, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, g, X_k, mu_k), stored by
        !! rows (the layout eu_constrained_solve takes).
        type(label3d_t), intent(in) :: s
        type(csr_t), intent(in) :: hessian
        real(dp), intent(in) :: dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(csc_t) :: m
        type(fortsparse_status_t) :: info
        integer :: n0, nt, nf, nx, i, j, k, c, p, nz, cap
        n0 = size(s%w)
        nt = s%nt
        nf = 2*n0*nt
        nx = 2*nt
        cap = size(hessian%val) + nx*size(dgdx, 1) + 6*n0*nt + 4*nt + &
              s%gauge_k*(3*n0*nt + 1) + n0*nt
        allocate (rows(cap), cols(cap), vals(cap))
        nz = 0
        do i = 1, nf
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call add(i, hessian%col(p), hessian%val(p))
            end do
        end do
        do j = 1, nx
            do i = 1, nf + 3*nt
                if (abs(dgdx(i, j)) > 0) call add(i, nf + j, dgdx(i, j))
            end do
            call add(nf + j, nf + j, 0.0_dp)
        end do
        do k = 1, nt
            do i = 1, n0
                do c = 1, 2
                    if (abs(s%ngrad(i, c, k)) > 0) &
                        call add(nf + 2*(k - 1) + c, (k - 1)*n0 + i, s%ngrad(i, c, k))
                end do
                if (abs(s%nrow(i)) > 0) then
                    call add(nf + nx + k, (k - 1)*n0 + i, s%nrow(i))
                    call add((k - 1)*n0 + i, nf + nx + k, s%nrow(i))
                end if
            end do
            call add(nf + nx + k, nf + nx + k, 0.0_dp)
        end do
        if (s%gauge_k > 0) then
            block
                real(dp) :: pk(s%gauge_k), dpk(s%gauge_k)
                integer :: kk, ntot
                ntot = n0*nt
                do i = 1, ntot
                    call legendre01(s%chi(i), s%gauge_k, pk, dpk)
                    do kk = 1, s%gauge_k
                        call add(ntot + i, nf + 3*nt + kk, pk(kk))
                        call add(nf + 3*nt + kk, ntot + i, pk(kk))
                        call add(nf + 3*nt + kk, i, dpk(kk)*s%g(i))
                    end do
                    call add(ntot + i, i, dot_product(dpk, s%mu(nt + 1:)))
                end do
                do kk = 1, s%gauge_k
                    call add(nf + 3*nt + kk, nf + 3*nt + kk, 0.0_dp)
                end do
            end block
        end if
        call csc_from_triplet(nf + 3*nt + s%gauge_k, nf + 3*nt + s%gauge_k, rows(:nz), cols(:nz), &
                              vals(:nz), m, info)
        if (info%code /= 0) error stop 'label3d: Newton triplets'
        jac_t%nrow = nf + 3*nt + s%gauge_k
        jac_t%ncol = nf + 3*nt + s%gauge_k
        jac_t%row_ptr = m%col_ptr
        jac_t%col = m%row_idx
        jac_t%val = m%val
    contains
        subroutine add(r, c2, v)
            integer, intent(in) :: r, c2
            real(dp), intent(in) :: v
            nz = nz + 1
            rows(nz) = r
            cols(nz) = c2
            vals(nz) = v
        end subroutine add
    end subroutine newton_matrix

    subroutine shift_chi(hessian, ntot, relative)
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: ntot
        real(dp), intent(in) :: relative
        real(dp) :: shift
        integer :: i, p
        shift = relative*maxval(abs(hessian%val))
        do i = 1, ntot
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine shift_chi

    subroutine add_ridge(hessian, ntot, ridge)
        !! g -> g + k(chi) leaves B unchanged: relative ridge on the g block.
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: ntot
        real(dp), intent(in) :: ridge
        real(dp) :: shift
        integer :: i, p
        shift = ridge*maxval(abs(hessian%val))
        do i = ntot + 1, 2*ntot
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine add_ridge

end module milca_label_ale_3d
