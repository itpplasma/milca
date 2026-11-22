module milca_label_ale
    !! Axisymmetric label lane on an ALE mesh, one magnetic axis (step A1 of
    !! docs/equilibrium-roadmap.md). Staggered:
    !!
    !! 1. physics Newton in (chi, g, X, mu):
    !!      grad W + mu [N(X); 0] = 0,  grad chi_h(X) = 0,  chi_h(X) = 0,
    !!    W the direct-Clebsch energy (milca_eulerian_clebsch) with the
    !!    winding form about X. X is a fixed parametric point (the centre of
    !!    one cell of the core block); inside Newton the core block translates
    !!    with X and the rest of the mesh follows by a fixed harmonic blend
    !!    (1 on the core, 0 on the wall), so the polar rule about X and its
    !!    neighbour cells keep their shape and the residual is smooth in X;
    !! 2. mesh equation between solves: the harmonic blend of the core
    !!    translation (default) or, optionally, Winslow between the core and
    !!    the wall. The core is a region, not a point: a point constraint on a
    !!    Winslow (harmonic) map is singular in the continuum, and pinning one
    !!    cell degraded the accuracy under refinement;
    !! 3. transfer of chi and g to the new mesh, g re-solved; back to 1 until
    !!    Newton moves X by less than a fraction of the cell size.
    !!
    !! The mesh coordinates never enter the physics Newton as unknowns: X is
    !! physical, the mesh follows it by a prescribed map.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_boundary_dofs
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_quadrature_point_singular, spline_locate, spline_geometry_at, &
        spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_boundary_dofs, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_profiles, clebsch_axisym_system
    use milca_spline_ale, only: spline_winslow, spline_transfer_v0, spline_geometry_dofs, &
        spline_set_geometry, spline_blend_weights, spline_cell_at, spline_block_dofs, &
        spline_min_jacobian
    implicit none
    private
    public :: label_ale_t, label_ale_init, label_ale_newton, label_ale_solve, &
              label_ale_rebuild, label_ale_initial_g

    type :: label_ale_t
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad, aquad
        real(dp), allocatable :: chi(:), g(:)
        real(dp) :: axis(2) = 0, mu = 0, cell_size = 1
        integer, allocatable :: wall(:), mesh_fixed(:)
        integer :: centre_block = 1, extra_points = 8, pin_cell = 0
        real(dp) :: centre_xi(2) = 0.5_dp
        real(dp) :: ridge = 1e-12_dp
        ! report
        integer :: passes = 0, newton_total = 0, winslow_total = 0
        real(dp), allocatable :: newton_history(:)
        ! the mesh inside Newton: x_ref + blend (X - anchor)
        type(spline_space_t) :: mesh
        real(dp), allocatable :: x_ref(:, :), blend(:)
        real(dp) :: anchor(2) = 0
    end type

contains

    subroutine label_ale_init(s, space, centre_block, centre_xi, extra_points, status)
        !! Mesh from `space` (C0, wall tagged LAYOUT_WALL) with the O-grid
        !! centre vertex at (centre_block, centre_xi).
        type(label_ale_t), intent(out) :: s
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: centre_block, extra_points
        real(dp), intent(in) :: centre_xi(2)
        integer, intent(out) :: status
        real(dp), allocatable :: x(:, :)
        s%extra_points = extra_points
        call spline_geometry_dofs(space, x, status)
        if (status /= 0) then
            status = 2
            return
        end if
        s%centre_block = centre_block
        s%centre_xi = centre_xi
        call spline_boundary_dofs(space, LAYOUT_WALL, s%mesh_fixed)
        s%axis = centre_point(space, s)
        s%cell_size = (maxval(x(1, :)) - minval(x(1, :)))/ &
                      sqrt(real(space%cell_count, dp))
        call label_ale_rebuild(s, space, status)
    end subroutine label_ale_init

    subroutine label_ale_rebuild(s, space, status)
        !! derham, quadratures and wall dofs on the mesh `space`.
        type(label_ale_t), intent(inout) :: s
        type(spline_space_t), intent(in) :: space
        integer, intent(out) :: status
        call spline_derham_build(space, s%derham, status, 0)
        if (status /= 0) then
            if (len_trim(env_debug()) > 0) print *, ' derham build status', status
            return
        end if
        call spline_quadrature_build(s%derham%space, s%derham%degree + s%extra_points, s%quad, &
                                     status)
        if (status /= 0) then
            if (len_trim(env_debug()) > 0) print *, ' quadrature build status', status
            return
        end if
        call spline_derham_boundary_dofs(s%derham, LAYOUT_WALL, 0, s%wall)
        call axis_quadrature(s, s%axis, status)
    end subroutine label_ale_rebuild

    subroutine axis_quadrature(s, x, status)
        type(label_ale_t), intent(inout) :: s
        real(dp), intent(in) :: x(2)
        integer, intent(out) :: status
        call spline_quadrature_point_singular(s%derham%space, s%quad, x, s%aquad, status)
    end subroutine axis_quadrature

    subroutine label_ale_newton(s, profiles, tolerance, max_steps, steps, status)
        !! Newton on (chi, g, X, mu) on the frozen mesh. The X columns are
        !! central differences of the residual at fixed mesh (polar rule about
        !! the shifted X). Status 0 converged, 1 linear solve, 2 not converged.
        type(label_ale_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: steps, status
        real(dp), allocatable :: r(:), rp(:), rm(:), rt(:), delta(:), nrow(:), ngrad(:, :), &
                                 dgdx(:, :)
        type(csr_t) :: hessian, jac_t
        real(dp) :: alpha, r0, rnorm
        real(dp), parameter :: hx = 1e-6_dp
        integer :: n0, it, k, st
        logical :: debug
        character(8) :: env
        call get_environment_variable('MILCA_LABEL_DEBUG', env)
        debug = len_trim(env) > 0
        n0 = s%derham%dimension(0)
        if (.not. allocated(s%newton_history)) allocate (s%newton_history(0))
        steps = 0
        status = 2
        r0 = -1
        do it = 1, max_steps + 1
            call residual(s%axis, s%chi, s%g, s%mu, r, hessian, nrow, ngrad)
            rnorm = norm2(r)
            if (debug) print '(a,4es12.4)', '  r parts chi g X mu', norm2(r(1:n0)), &
                norm2(r(n0 + 1:2*n0)), norm2(r(2*n0 + 1:2*n0 + 2)), r(2*n0 + 3)
            s%newton_history = [s%newton_history, rnorm]
            if (r0 < 0) r0 = rnorm
            if (rnorm < tolerance*max(r0, 1.0_dp) .and. it > 1) then
                status = 0
                exit
            end if
            if (it > max_steps) exit
            allocate (dgdx(2*n0 + 3, 2))
            do k = 1, 2
                call residual(s%axis + hx*unit(k), s%chi, s%g, s%mu, rp)
                call residual(s%axis - hx*unit(k), s%chi, s%g, s%mu, rm)
                dgdx(:, k) = (rp - rm)/(2*hx)
            end do
            call add_ridge(hessian, n0, s%ridge)
            call newton_matrix(hessian, nrow, ngrad, dgdx, jac_t)
            allocate (delta(2*n0 + 3))
            call eu_constrained_solve(jac_t, -r, s%wall, 0*s%chi(s%wall), delta, st)
            if (st /= 0) then
                status = 1
                return
            end if
            if (debug .and. it == 1) then
                block
                    real(dp), allocatable :: jd(:), rpp(:), rmm(:), dd(:)
                    real(dp) :: ee
                    integer :: blk
                    ee = 1e-6_dp
                    dd = delta
                    do blk = 1, 4
                        if (blk == 1) dd = [delta(1:n0), 0*delta(n0 + 1:)]
                        if (blk == 2) dd = [0*delta(1:n0), delta(n0 + 1:2*n0), 0*delta(2*n0 + 1:)]
                        if (blk == 3) dd = [0*delta(1:2*n0), delta(2*n0 + 1:2*n0 + 2), 0.0_dp]
                        if (blk == 4) dd = [0*delta(1:2*n0 + 2), delta(2*n0 + 3)]
                        call residual(s%axis + ee*dd(2*n0 + 1:2*n0 + 2), s%chi + ee*dd(1:n0), &
                                      s%g + ee*dd(n0 + 1:2*n0), s%mu + ee*dd(2*n0 + 3), rpp)
                        call residual(s%axis - ee*dd(2*n0 + 1:2*n0 + 2), s%chi - ee*dd(1:n0), &
                                      s%g - ee*dd(n0 + 1:2*n0), s%mu - ee*dd(2*n0 + 3), rmm)
                        jd = jmul(jac_t, dd)
                        jd(s%wall) = 0
                        print '(a,i2,a,es10.2,a,es10.2)', '  block', blk, ' |J d - fd| ', &
                            norm2(jd - (rpp - rmm)/(2*ee)), ' |fd| ', norm2((rpp - rmm)/(2*ee))
                    end do
                end block
            end if
            ! X moves at most a quarter cell per step (the tied mesh)
            alpha = min(1.0_dp, 0.25_dp*s%cell_size/max(norm2(delta(2*n0 + 1:2*n0 + 2)), &
                                                         tiny(1.0_dp)))
            do k = 1, 12
                call residual(s%axis + alpha*delta(2*n0 + 1:2*n0 + 2), &
                              s%chi + alpha*delta(1:n0), s%g + alpha*delta(n0 + 1:2*n0), &
                              s%mu + alpha*delta(2*n0 + 3), rt)
                if (norm2(rt) < (1 - alpha/4)*rnorm) exit
                alpha = alpha/2
            end do
            if (debug) print '(a,i3,a,es10.2,a,2es11.3,a,es10.2,a,es10.2)', '  newton', it, ' |r|', &
                rnorm, ' dX', delta(2*n0 + 1:2*n0 + 2), ' alpha', alpha, ' |dchi|', &
                maxval(abs(delta(1:n0)))
            s%chi = s%chi + alpha*delta(1:n0)
            s%g = s%g + alpha*delta(n0 + 1:2*n0)
            s%axis = s%axis + alpha*delta(2*n0 + 1:2*n0 + 2)
            s%mu = s%mu + alpha*delta(2*n0 + 3)
            deallocate (delta, dgdx)
            steps = steps + 1
        end do
        call tie_mesh(s, s%axis, st)
        s%newton_total = s%newton_total + steps

    contains

        function jmul(jt, v) result(y)
            !! J v from the CSR of J^T
            type(csr_t), intent(in) :: jt
            real(dp), intent(in) :: v(:)
            real(dp), allocatable :: y(:)
            integer :: i, p
            allocate (y(jt%ncol), source=0.0_dp)
            do i = 1, jt%nrow
                do p = jt%row_ptr(i), jt%row_ptr(i + 1) - 1
                    y(jt%col(p)) = y(jt%col(p)) + jt%val(p)*v(i)
                end do
            end do
        end function jmul

        function unit(k) result(e)
            integer, intent(in) :: k
            real(dp) :: e(2)
            e = 0
            e(k) = 1
        end function unit

        subroutine residual(x, chi, g, mu, rr, hm, nr, ng)
            real(dp), intent(in) :: x(2), chi(:), g(:), mu
            real(dp), allocatable, intent(out) :: rr(:)
            type(csr_t), intent(out), optional :: hm
            real(dp), allocatable, intent(out), optional :: nr(:), ng(:, :)
            type(csr_t) :: htmp
            real(dp), allocatable :: nrw(:), ngr(:, :), gr(:)
            real(dp) :: e
            integer :: stq
            call tie_mesh(s, x, stq)
            if (stq /= 0) then
                ! the tied mesh does not hold X: an infinite residual
                rr = [(huge(1.0_dp)/1e10_dp, stq=1, 2*size(chi) + 3)]
                if (present(hm)) error stop 'label_ale: tied mesh at the current X'
                return
            end if
            call vertex_rows(s, nrw, ngr)
            call clebsch_axisym_system(s%derham, s%aquad, profiles, x, chi, g, e, gr, htmp)
            rr = [gr + mu*[nrw, 0*nrw], matmul(chi, ngr), dot_product(nrw, chi)]
            rr(s%wall) = 0
            if (present(hm)) hm = htmp
            if (present(nr)) nr = nrw
            if (present(ng)) ng = ngr
            associate (unused => e)
            end associate
        end subroutine residual

    end subroutine label_ale_newton

    subroutine label_ale_solve(s, space, profiles, tolerance, mesh_tolerance, max_passes, &
                               max_steps, status, winslow)
        !! The staggered ALE loop: physics Newton on the frozen mesh, then the
        !! Winslow mesh with its centre at X and the transfer of chi and g,
        !! until |centre - X| < mesh_tolerance * cell size. `space` is the
        !! current mesh on entry and the final one on return. winslow
        !! (default false) relaxes the mesh between core and wall by Winslow
        !! after the blend. Status 0 converged, 1 Newton, 2 mesh, 3 transfer,
        !! 4 passes exhausted.
        type(label_ale_t), intent(inout) :: s
        type(spline_space_t), intent(inout) :: space
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: tolerance, mesh_tolerance
        integer, intent(in) :: max_passes, max_steps
        integer, intent(out) :: status
        logical, intent(in), optional :: winslow
        type(spline_derham_t) :: old
        real(dp), allocatable :: x(:, :), chi(:), g(:)
        integer :: pass, steps, st, wsteps
        ! the pin: centre of the cell that holds X, kept for all passes (X
        ! stays at a fixed parametric point inside one cell, away from edges)
        block
            real(dp) :: xi(2)
            integer :: b
            call spline_locate(space, 1, s%axis, b, xi, st)
            if (st /= 0) then
                status = 2
                return
            end if
            s%centre_block = b
            call spline_cell_at(space, b, xi, s%pin_cell, s%centre_xi)
        end block
        status = 4
        do pass = 1, max_passes
            s%passes = pass
            if (len_trim(env_debug()) > 0) print '(a,i3,2es12.4)', '  pass, |centre - X|, tol', pass, &
                norm2(centre_point(space, s) - s%axis), mesh_tolerance*s%cell_size
            if (pass > 1 .or. norm2(centre_point(space, s) - s%axis) > 1e-12_dp*s%cell_size) then
                ! mesh equation: the core block translated so that the pin is
                ! at X, the harmonic blend elsewhere (optionally Winslow
                ! between core and wall), then transfer
                old = s%derham
                block
                    integer, allocatable :: core_cells(:), core(:)
                    real(dp), allocatable :: w(:), x0(:, :)
                    real(dp) :: shift(2), frac, xi(2)
                    type(spline_quadrature_t) :: probe
                    integer :: i, k, b
                    call spline_block_dofs(space, s%centre_block, core_cells, core)
                    shift = s%axis - centre_point(space, s)
                    call spline_blend_weights(space, core_cells, s%mesh_fixed, w, st)
                    call spline_geometry_dofs(space, x0, st)
                    ! the full centring, halved while the moved mesh folds at
                    ! the quadrature points of the physics (a coarse ring of
                    ! cells cannot take a large core translation)
                    frac = 1
                    do k = 1, 6
                        x = x0
                        do i = 1, size(w)
                            x(:, i) = x(:, i) + frac*w(i)*shift
                        end do
                        call spline_set_geometry(space, x)
                        call spline_quadrature_build(space, space%degree + s%extra_points, probe, &
                                                     st)
                        if (st == 0) exit
                        frac = frac/2
                    end do
                    if (st /= 0) then
                        frac = 0
                        call spline_set_geometry(space, x0)
                    end if
                    if (frac < 1) then
                        ! partial centring: X stays at its own point of the core
                        call spline_locate(space, 1, s%axis, b, xi, st)
                        if (st /= 0 .or. b /= s%centre_block) then
                            st = 1
                        else
                            s%centre_xi = xi
                            call spline_cell_at(space, b, xi, s%pin_cell, xi)
                        end if
                    end if
                    wsteps = 0
                    if (st == 0 .and. present(winslow)) then
                        if (winslow) call spline_winslow(space, [s%mesh_fixed, core], 0, &
                                                         s%centre_xi, s%axis, 1e-10_dp, 50, &
                                                         wsteps, st)
                    end if
                end block
                s%winslow_total = s%winslow_total + wsteps
                if (st /= 0) then
                    status = 2
                    return
                end if
                call label_ale_rebuild(s, space, st)
                if (st /= 0) then
                    if (len_trim(env_debug()) > 0) print '(a,i0,a,es10.2)', &
                        '  rebuild of the moved mesh failed: ', st, ', min det ', &
                        spline_min_jacobian(space)
                    status = 2
                    return
                end if
                allocate (chi(size(s%chi)), g(size(s%g)))
                call spline_transfer_v0(old, s%chi, s%derham, chi, st)
                if (st == 0) call spline_transfer_v0(old, s%g, s%derham, g, st)
                if (st /= 0) then
                    status = 3
                    return
                end if
                chi(s%wall) = 1
                call move_alloc(chi, s%chi)
                call move_alloc(g, s%g)
                call axis_quadrature(s, s%axis, st)
                ! g is slaved to chi (W convex in g): re-solve on the new mesh
                call label_ale_initial_g(s, profiles, st)
            end if
            call label_ale_anchor(s, space, st)
            if (st /= 0) then
                status = 2
                return
            end if
            ! physics Newton; the singular patch follows X by the fixed blend
            call label_ale_newton(s, profiles, tolerance, max_steps, steps, st)
            if (st == 1) then
                status = 1
                return
            end if
            space = s%mesh
            if (st == 0 .and. norm2(s%axis - s%anchor) < mesh_tolerance*s%cell_size) then
                status = 0
                return
            end if
            if (st == 0) then
                ! converged, but the blend moved the centre far: remesh once more
                continue
            end if
        end do
    end subroutine label_ale_solve

    subroutine label_ale_initial_g(s, profiles, status)
        !! g minimising W at fixed labels (W is quadratic and convex in g).
        type(label_ale_t), intent(inout) :: s
        procedure(clebsch_profiles) :: profiles
        integer, intent(out) :: status
        real(dp), allocatable :: grad(:), delta(:)
        type(csr_t) :: hessian
        real(dp) :: energy
        integer :: n0, i
        n0 = s%derham%dimension(0)
        call clebsch_axisym_system(s%derham, s%aquad, profiles, s%axis, s%chi, s%g, energy, &
                                   grad, hessian)
        call add_ridge(hessian, n0, s%ridge)
        allocate (delta(2*n0))
        call eu_constrained_solve(hessian, -grad, [(i, i=1, n0)], 0*s%chi, delta, status)
        if (status == 0) s%g = s%g + delta(n0 + 1:)
    end subroutine label_ale_initial_g

    subroutine label_ale_anchor(s, space, status)
        !! The Newton mesh: `space` translated by blend (X - anchor), blend 1
        !! on the centre block and 0 on the wall, anchor = X now.
        type(label_ale_t), intent(inout) :: s
        type(spline_space_t), intent(in) :: space
        integer, intent(out) :: status
        s%mesh = space
        call spline_geometry_dofs(space, s%x_ref, status)
        if (status /= 0) return
        block
            integer, allocatable :: core_cells(:), core(:)
            call spline_block_dofs(space, s%centre_block, core_cells, core)
            call spline_blend_weights(space, core_cells, s%mesh_fixed, s%blend, status)
        end block
        s%anchor = centre_point(space, s)
    end subroutine label_ale_anchor

    subroutine tie_mesh(s, x, status)
        type(label_ale_t), intent(inout) :: s
        real(dp), intent(in) :: x(2)
        integer, intent(out) :: status
        real(dp), allocatable :: xg(:, :)
        integer :: i
        xg = s%x_ref
        do i = 1, size(xg, 2)
            xg(:, i) = xg(:, i) + s%blend(i)*(x - s%anchor)
        end do
        call spline_set_geometry(s%mesh, xg)
        call label_ale_rebuild(s, s%mesh, status)
        if (status /= 0) return
        call axis_quadrature(s, x, status)
    end subroutine tie_mesh

    function env_debug() result(v)
        character(8) :: v
        call get_environment_variable('MILCA_LABEL_DEBUG', v)
    end function env_debug

    function centre_point(space, s) result(p)
        type(spline_space_t), intent(in) :: space
        type(label_ale_t), intent(in) :: s
        real(dp) :: p(2), jet(2, 6)
        integer :: st
        call spline_geometry_at(space, s%centre_block, s%centre_xi, jet, st)
        p = jet(:, 1)
    end function centre_point

    subroutine newton_matrix(hessian, nrow, ngrad, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, g, X_R, X_Z, mu), stored by
        !! rows (the layout eu_constrained_solve takes).
        type(csr_t), intent(in) :: hessian
        real(dp), intent(in) :: nrow(:), ngrad(:, :), dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(csc_t) :: m
        type(fortsparse_status_t) :: info
        integer :: n0, nf, i, k, p, nz, cap
        n0 = size(nrow)
        nf = 2*n0
        cap = size(hessian%val) + 2*size(dgdx, 1) + 4*n0 + 8
        allocate (rows(cap), cols(cap), vals(cap))
        nz = 0
        ! J(i, j) entries; J = H on the (chi, g) block
        do i = 1, nf
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call add(i, hessian%col(p), hessian%val(p))
            end do
        end do
        do k = 1, 2
            do i = 1, nf + 3
                if (dgdx(i, k) /= 0) call add(i, nf + k, dgdx(i, k))
            end do
            call add(nf + k, nf + k, 0.0_dp)
            do i = 1, n0
                if (ngrad(i, k) /= 0) call add(nf + k, i, ngrad(i, k))
            end do
        end do
        do i = 1, n0
            if (nrow(i) /= 0) then
                call add(i, nf + 3, nrow(i))
                call add(nf + 3, i, nrow(i))
            end if
        end do
        call add(nf + 3, nf + 3, 0.0_dp)
        ! CSC of J = CSR of J^T
        call csc_from_triplet(nf + 3, nf + 3, rows(:nz), cols(:nz), vals(:nz), m, info)
        if (info%code /= 0) error stop 'label_ale: Newton triplets'
        jac_t%nrow = nf + 3
        jac_t%ncol = nf + 3
        jac_t%row_ptr = m%col_ptr
        jac_t%col = m%row_idx
        jac_t%val = m%val
    contains
        subroutine add(r, c, v)
            integer, intent(in) :: r, c
            real(dp), intent(in) :: v
            nz = nz + 1
            rows(nz) = r
            cols(nz) = c
            vals(nz) = v
        end subroutine add
    end subroutine newton_matrix

    subroutine vertex_rows(s, nrow, grows)
        !! Rows of chi_h(X) and grad chi_h(X), X at the pin (an interior
        !! point of one cell), analytic (spline_point_matrix).
        type(label_ale_t), intent(in) :: s
        real(dp), allocatable, intent(out) :: nrow(:), grows(:, :)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        integer :: st, k
        allocate (nrow(s%derham%dimension(0)), grows(s%derham%dimension(0), 2), source=0.0_dp)
        do k = 0, 2
            call spline_point_matrix(s%derham%space, [s%centre_block], &
                                     reshape(s%centre_xi, [2, 1]), k, rp, ci, v, st)
            if (st /= 0) error stop 'label_ale: pin rows'
            if (k == 0) then
                nrow(ci) = v
            else
                grows(ci, k) = v
            end if
        end do
    end subroutine vertex_rows

    subroutine add_ridge(hessian, n0, ridge)
        !! g -> g + k(chi) leaves B unchanged: relative ridge on the g block.
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: n0
        real(dp), intent(in) :: ridge
        real(dp) :: shift
        integer :: i, p
        shift = ridge*maxval(abs(hessian%val))
        do i = n0 + 1, 2*n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine add_ridge

end module milca_label_ale
