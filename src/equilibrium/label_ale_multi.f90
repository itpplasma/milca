module milca_label_ale_multi
    !! Axisymmetric label lane on an ALE mesh with several magnetic axes
    !! (step A2 of docs/equilibrium-roadmap.md), in the Grad-Shafranov form
    !! of the label: psi = psi_O + s chi with the flux scale s given, the
    !! poloidal current f(chi) prescribed (R B_phi = f), so that
    !!
    !!   W(chi) = int [ s^2 |grad chi|^2/(2 R^2) - p(chi) - f(chi)^2/(2 R^2) ] dV
    !!
    !! is stationary at Grad-Shafranov equilibria. Unknowns chi (V0, chi = 2
    !! on the wall), the O-points O_j and multipliers mu_j:
    !!
    !!   grad W + sum_j mu_j N(O_j) = 0,  grad chi_h(O_j) = 0,  chi_h(O_j) = 0.
    !!
    !! The toroidal field needs no winding form here, so the separatrix carries
    !! no log-singular profile (the open part of gate E2, the direct-Clebsch
    !! toroidal representation, stays separate). Each O_j is a fixed
    !! parametric point of its own core block; inside Newton core j translates
    !! with O_j and the mesh follows the sum of harmonic blends w_j (1 on core
    !! j; 0 on the wall, on the other cores and on the separatrix sides), as in
    !! milca_label_ale. Between passes the X-point vertex follows the discrete
    !! saddle of chi_h by its own blend (1 at X, 0 on the wall and the cores).
    !! The mesh coordinates are never unknowns.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_locate, spline_geometry_at, spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_boundary_dofs, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_spline_ale, only: spline_geometry_dofs, spline_set_geometry, &
        spline_blend_weights, spline_cell_at, spline_block_dofs, spline_transfer_v0
    implicit none
    private
    public :: label_multi_t, gs_profiles, label_multi_init, label_multi_newton, &
              label_multi_solve, label_multi_saddle, label_multi_rebuild

    abstract interface
        subroutine gs_profiles(chi, p, ff)
            !! p = (p, p', p''), ff = (f^2/2, (f^2/2)', (f^2/2)'') at chi.
            import :: dp
            real(dp), intent(in) :: chi
            real(dp), intent(out) :: p(3), ff(3)
        end subroutine gs_profiles
    end interface

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: label_multi_t
        type(spline_derham_t) :: derham
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: chi(:), axes(:, :), mu(:), centre_xi(:, :), anchor(:, :)
        integer, allocatable :: core_block(:), wall(:), fixed(:)
        real(dp) :: scale = 1, cell_size = 1, wall_value = 2
        integer :: points = 4, xpoint_dof = 0
        real(dp) :: xpoint(2) = 0, xanchor(2) = 0
        ! the mesh inside Newton: x_ref + sum_j blend(:, j) (O_j - anchor_j)
        ! + xblend (X - xanchor)
        type(spline_space_t) :: mesh, orig
        real(dp), allocatable :: x_ref(:, :), blend(:, :), xblend(:)
        integer :: newton_total = 0, passes = 0
        real(dp), allocatable :: newton_history(:)
    end type

contains

    subroutine label_multi_init(s, space, core_block, axes, fixed, xpoint_dof, scale, points, &
                                status, xpoint_start)
        !! Mesh `space` (C0, wall tagged LAYOUT_WALL), core block per axis,
        !! initial O-points `axes` (inside their cores), dofs `fixed` held by
        !! every axis blend (the separatrix sides), the X-point vertex dof
        !! (0: no X-point), flux scale s, Gauss points per direction, and
        !! optionally a start for the X vertex (default: where it is).
        !! Each O_j is pinned at the centre of the core cell holding the initial
        !! O_j (else the central cell), and the cores are moved so that the
        !! pins are at the initial O_j
        !! (status 2: the moved mesh folds, 3: no central cell).
        type(label_multi_t), intent(out) :: s
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: core_block(:), fixed(:), xpoint_dof, points
        real(dp), intent(in) :: axes(:, :), scale
        integer, intent(out) :: status
        real(dp), intent(in), optional :: xpoint_start(2)
        integer, allocatable :: cells(:), dofs(:), zero(:)
        real(dp), allocatable :: wj(:)
        real(dp) :: jet(2, 6), xi(2)
        integer :: m, j, k, cell, st, b
        logical, allocatable :: mark(:)

        m = size(core_block)
        s%core_block = core_block
        s%axes = axes
        s%scale = scale
        s%points = points
        s%fixed = fixed
        s%xpoint_dof = xpoint_dof
        allocate (s%mu(m), s%centre_xi(2, m), s%anchor(2, m))
        s%mu = 0
        s%mesh = space
        s%orig = space
        call spline_geometry_dofs(space, s%x_ref, status)
        if (status /= 0) return
        s%cell_size = (maxval(s%x_ref(1, :)) - minval(s%x_ref(1, :)))/ &
                      sqrt(real(space%cell_count, dp))
        ! pins: the centre of a core cell; the core is then moved so that
        ! the pin is at O_j
        do j = 1, m
            ! the core cell that holds the initial O_j, else the central one
            call spline_locate(space, 1, axes(:, j), b, xi, st)
            if (st /= 0 .or. b /= core_block(j)) xi = 0.5_dp - 1e-9_dp
            call spline_cell_at(space, core_block(j), xi, cell, s%centre_xi(:, j))
            call spline_geometry_at(space, core_block(j), s%centre_xi(:, j), jet, st)
            if (st /= 0 .or. cell == 0) then
                status = 3
                return
            end if
            s%anchor(:, j) = jet(:, 1)
        end do
        ! blends: core j to 1, the wall, the other cores, `fixed` and X to 0
        call spline_derham_build(space, s%derham, status, 0)
        if (status /= 0) return
        call spline_derham_boundary_dofs(s%derham, LAYOUT_WALL, 0, s%wall)
        allocate (s%blend(space%dof_count, m), mark(space%dof_count))
        do j = 1, m
            mark = .false.
            mark(s%wall) = .true.
            mark(fixed) = .true.
            if (xpoint_dof > 0) mark(xpoint_dof) = .true.
            do k = 1, m
                if (k == j) cycle
                call spline_block_dofs(space, core_block(k), cells, dofs)
                mark(dofs) = .true.
            end do
            call spline_block_dofs(space, core_block(j), cells, dofs)
            mark(dofs) = .false.
            zero = pack([(k, k=1, space%dof_count)], mark)
            call spline_blend_weights(space, cells, zero, wj, status)
            if (status /= 0) return
            s%blend(:, j) = wj
        end do
        allocate (s%xblend(space%dof_count), source=0.0_dp)
        if (xpoint_dof > 0) then
            s%xpoint = s%x_ref(:, xpoint_dof)
            s%xanchor = s%xpoint
            ! X to 1, wall and cores to 0: blend about the single vertex X
            mark = .false.
            mark(s%wall) = .true.
            do k = 1, m
                call spline_block_dofs(space, core_block(k), cells, dofs)
                mark(dofs) = .true.
            end do
            zero = pack([(k, k=1, space%dof_count)], mark)
            call vertex_blend(space, xpoint_dof, zero, s%xblend, status)
            if (status /= 0) return
            if (present(xpoint_start)) s%xpoint = xpoint_start
        end if
        ! move the cores to the initial O_j, halved while the mesh folds (a
        ! coarse ring cannot take a large core translation): the start is
        ! then the partially moved pin
        block
            real(dp) :: frac, xtarget(2)
            real(dp), allocatable :: target(:, :)
            target = axes
            xtarget = s%xpoint
            frac = 1
            do k = 1, 8
                s%axes = s%anchor + frac*(target - s%anchor)
                s%xpoint = s%xanchor + frac*(xtarget - s%xanchor)
                call tie_mesh(s, s%axes, s%xpoint, status)
                if (status == 0) exit
                frac = frac/2
            end do
            if (status /= 0) then
                s%axes = s%anchor
                s%xpoint = s%xanchor
                call tie_mesh(s, s%axes, s%xpoint, status)
            end if
        end block
        if (status /= 0) status = 2
    end subroutine label_multi_init

    subroutine vertex_blend(space, dof, zero, w, status)
        !! Reference-harmonic weights, 1 at one dof, 0 on `zero`.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: dof, zero(:)
        real(dp), allocatable, intent(out) :: w(:)
        integer, intent(out) :: status
        integer, allocatable :: none(:)
        allocate (none(0))
        call spline_blend_weights(space, none, zero, w, status, inner_dofs=[dof])
    end subroutine vertex_blend

    subroutine label_multi_rebuild(s, status)
        !! derham and quadrature on the current mesh.
        type(label_multi_t), intent(inout) :: s
        integer, intent(out) :: status
        call spline_derham_build(s%mesh, s%derham, status, 0)
        if (status /= 0) return
        call spline_quadrature_build(s%derham%space, s%points, s%quad, status)
    end subroutine label_multi_rebuild

    subroutine tie_mesh(s, axes, xpoint, status)
        type(label_multi_t), intent(inout) :: s
        real(dp), intent(in) :: axes(:, :), xpoint(2)
        integer, intent(out) :: status
        real(dp), allocatable :: xg(:, :)
        integer :: i, j
        xg = s%x_ref
        do i = 1, size(xg, 2)
            do j = 1, size(axes, 2)
                xg(:, i) = xg(:, i) + s%blend(i, j)*(axes(:, j) - s%anchor(:, j))
            end do
            xg(:, i) = xg(:, i) + s%xblend(i)*(xpoint - s%xanchor)
        end do
        call spline_set_geometry(s%mesh, xg)
        call label_multi_rebuild(s, status)
    end subroutine tie_mesh

    subroutine energy_system(s, profiles, chi, grad, hessian)
        !! Gradient and Hessian of W in the V0 node values of chi.
        type(label_multi_t), intent(in) :: s
        procedure(gs_profiles) :: profiles
        real(dp), intent(in) :: chi(:)
        real(dp), allocatable, intent(out) :: grad(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), vals(:), kl(:, :)
        integer, allocatable :: d0(:), ri(:), ci(:)
        real(dp) :: c, gc(2), w, r, pp(3), ff(3), a, b2
        integer :: cell, q, p, i, k, n, nz, nloc, ncell
        type(csc_t) :: m
        type(fortsparse_status_t) :: info
        n = size(chi)
        ncell = size(s%quad%point_offset)
        allocate (grad(n), source=0.0_dp)
        nz = 0
        do cell = 1, ncell
            call clebsch_node_gradients(s%derham, s%quad, cell, v0, d0, gn)
            nloc = size(d0)
            if (.not. allocated(vals)) allocate (ri(ncell*nloc**2), ci(ncell*nloc**2), &
                                                 vals(ncell*nloc**2), kl(nloc, nloc))
            kl = 0
            do q = 1, s%quad%points_per_cell
                p = s%quad%point_offset(cell) + q - 1
                r = s%quad%s(1, p)
                w = 2*pi*r*s%quad%weight(p)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                call profiles(c, pp, ff)
                ! dW/dchi_i = int s^2 grad N_i . grad chi / R^2 - (p' + ff'/R^2) N_i
                a = s%scale**2/r**2
                b2 = -(pp(3) + ff(3)/r**2)
                grad(d0) = grad(d0) + w*(a*matmul(gc, gn(:, :, q)) - &
                                         (pp(2) + ff(2)/r**2)*v0(1, :, q))
                if (present(hessian)) then
                    do k = 1, nloc
                        kl(:, k) = kl(:, k) + w*(a*matmul(gn(:, k, q), gn(:, :, q)) + &
                                                 b2*v0(1, :, q)*v0(1, k, q))
                    end do
                end if
            end do
            if (present(hessian)) then
                do k = 1, nloc
                    do i = 1, nloc
                        nz = nz + 1
                        ri(nz) = d0(i)
                        ci(nz) = d0(k)
                        vals(nz) = kl(i, k)
                    end do
                end do
            end if
        end do
        if (present(hessian)) then
            call csc_from_triplet(n, n, ri(:nz), ci(:nz), vals(:nz), m, info)
            if (info%code /= 0) error stop 'label_ale_multi: Hessian triplets'
            hessian%nrow = n
            hessian%ncol = n
            hessian%row_ptr = m%col_ptr
            hessian%col = m%row_idx
            hessian%val = m%val
        end if
    end subroutine energy_system

    subroutine pin_rows(s, j, nrow, grows)
        !! Rows of chi_h and grad chi_h at the pin of axis j (analytic).
        type(label_multi_t), intent(in) :: s
        integer, intent(in) :: j
        real(dp), allocatable, intent(out) :: nrow(:), grows(:, :)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        integer :: st, k
        allocate (nrow(s%derham%dimension(0)), grows(s%derham%dimension(0), 2), source=0.0_dp)
        do k = 0, 2
            call spline_point_matrix(s%derham%space, [s%core_block(j)], &
                                     reshape(s%centre_xi(:, j), [2, 1]), k, rp, ci, v, st)
            if (st /= 0) error stop 'label_ale_multi: pin rows'
            if (k == 0) then
                nrow(ci) = v
            else
                grows(ci, k) = v
            end if
        end do
    end subroutine pin_rows

    subroutine label_multi_newton(s, profiles, tolerance, max_steps, steps, status)
        !! Newton on (chi, O_j, mu_j) with the cores tied to O_j. The O columns
        !! are central differences of the residual (the mesh moves with O).
        !! Status 0 converged, 1 linear solve, 2 not converged.
        type(label_multi_t), intent(inout) :: s
        procedure(gs_profiles) :: profiles
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: steps, status
        real(dp), allocatable :: r(:), rp(:), rm(:), rt(:), delta(:), dgdx(:, :), xp(:, :), &
                                 nrows(:, :), grows(:, :, :)
        type(csr_t) :: hessian, jac_t
        real(dp) :: alpha, r0, rnorm, dmax
        real(dp), parameter :: hx = 1e-6_dp
        integer :: n0, m, it, k, st, nt, j
        n0 = size(s%chi)
        m = size(s%core_block)
        nt = n0 + 3*m
        if (.not. allocated(s%newton_history)) allocate (s%newton_history(0))
        steps = 0
        status = 2
        r0 = -1
        do it = 1, max_steps + 1
            call residual(s%axes, s%chi, s%mu, r, hessian, nrows, grows)
            rnorm = norm2(r)
            s%newton_history = [s%newton_history, rnorm]
            if (r0 < 0) r0 = rnorm
            if (rnorm < tolerance*max(r0, 1.0_dp) .and. it > 1) then
                status = 0
                exit
            end if
            if (it > max_steps) exit
            allocate (dgdx(nt, 2*m))
            do k = 1, 2*m
                xp = s%axes
                j = (k + 1)/2
                xp(mod(k - 1, 2) + 1, j) = s%axes(mod(k - 1, 2) + 1, j) + hx
                call residual(xp, s%chi, s%mu, rp)
                xp(mod(k - 1, 2) + 1, j) = s%axes(mod(k - 1, 2) + 1, j) - hx
                call residual(xp, s%chi, s%mu, rm)
                dgdx(:, k) = (rp - rm)/(2*hx)
            end do
            call newton_matrix(hessian, nrows, grows, dgdx, jac_t)
            allocate (delta(nt))
            call eu_constrained_solve(jac_t, -r, s%wall, 0*s%chi(s%wall), delta, st)
            if (st /= 0) then
                status = 1
                return
            end if
            if (len_trim(env_debug()) > 0) then
                block
                    real(dp), allocatable :: rpp(:), rmm(:), jd(:), dd(:)
                    real(dp) :: ee
                    integer :: blk, ii, pp
                    ee = 1e-7_dp
                    do blk = 1, 3
                        dd = 0*delta
                        if (blk == 1) dd(1:n0) = delta(1:n0)
                        if (blk == 2) dd(n0 + 1:n0 + 2*m) = delta(n0 + 1:n0 + 2*m)
                        if (blk == 3) dd(n0 + 2*m + 1:) = delta(n0 + 2*m + 1:)
                        call residual(s%axes + ee*reshape(dd(n0 + 1:n0 + 2*m), [2, m]), &
                                      s%chi + ee*dd(1:n0), s%mu + ee*dd(n0 + 2*m + 1:), rpp)
                        call residual(s%axes - ee*reshape(dd(n0 + 1:n0 + 2*m), [2, m]), &
                                      s%chi - ee*dd(1:n0), s%mu - ee*dd(n0 + 2*m + 1:), rmm)
                        allocate (jd(nt), source=0.0_dp)
                        do ii = 1, nt
                            do pp = jac_t%row_ptr(ii), jac_t%row_ptr(ii + 1) - 1
                                jd(jac_t%col(pp)) = jd(jac_t%col(pp)) + jac_t%val(pp)*dd(ii)
                            end do
                        end do
                        jd(s%wall) = 0
                        print '(a,i2,a,es10.2,a,es10.2)', '  J check block', blk, ' |J d - fd|', &
                            norm2(jd - (rpp - rmm)/(2*ee)), ' |fd|', norm2((rpp - rmm)/(2*ee))
                        deallocate (jd)
                    end do
                end block
            end if
            ! each O moves at most one cell per step (the tied mesh; a fold is
            ! rejected by the line search)
            dmax = 0
            do j = 1, m
                dmax = max(dmax, norm2(delta(n0 + 2*j - 1:n0 + 2*j)))
            end do
            alpha = min(1.0_dp, s%cell_size/max(dmax, tiny(1.0_dp)))
            do k = 1, 12
                call residual(s%axes + alpha*reshape(delta(n0 + 1:n0 + 2*m), [2, m]), &
                              s%chi + alpha*delta(1:n0), s%mu + alpha*delta(n0 + 2*m + 1:), rt)
                if (norm2(rt) < (1 - alpha/4)*rnorm) exit
                alpha = alpha/2
            end do
            if (len_trim(env_debug()) > 0) print '(a,i3,a,es10.2,a,es10.2,a,es10.2,a,es10.2)', &
                '  newton', it, ' |r|', rnorm, ' |dO|', dmax, ' cap', s%cell_size, &
                ' alpha', alpha
            s%chi = s%chi + alpha*delta(1:n0)
            s%axes = s%axes + alpha*reshape(delta(n0 + 1:n0 + 2*m), [2, m])
            s%mu = s%mu + alpha*delta(n0 + 2*m + 1:)
            deallocate (delta, dgdx)
            steps = steps + 1
        end do
        call tie_mesh(s, s%axes, s%xpoint, st)
        s%newton_total = s%newton_total + steps

    contains

        subroutine residual(x, chi, mu, rr, hm, nr, gr)
            real(dp), intent(in) :: x(:, :), chi(:), mu(:)
            real(dp), allocatable, intent(out) :: rr(:)
            type(csr_t), intent(out), optional :: hm
            real(dp), allocatable, intent(out), optional :: nr(:, :), gr(:, :, :)
            real(dp), allocatable :: g(:), nrow(:), grow(:, :)
            integer :: stq, jj
            call tie_mesh(s, x, s%xpoint, stq)
            if (stq /= 0) then
                ! the tied mesh folds: an infinite residual rejects the trial
                rr = [(huge(1.0_dp)/1e10_dp, stq=1, size(chi) + 3*size(mu))]
                if (present(hm)) error stop 'label_ale_multi: tied mesh folds at the current O'
                return
            end if
            if (present(hm)) then
                call energy_system(s, profiles, chi, g, hm)
                allocate (nr(size(chi), size(mu)), gr(size(chi), 2, size(mu)))
            else
                call energy_system(s, profiles, chi, g)
            end if
            allocate (rr(size(chi) + 3*size(mu)))
            rr(1:size(chi)) = g
            do jj = 1, size(mu)
                call pin_rows(s, jj, nrow, grow)
                rr(1:size(chi)) = rr(1:size(chi)) + mu(jj)*nrow
                rr(size(chi) + 2*jj - 1:size(chi) + 2*jj) = matmul(chi, grow)
                rr(size(chi) + 2*size(mu) + jj) = dot_product(nrow, chi)
                if (present(nr)) then
                    nr(:, jj) = nrow
                    gr(:, :, jj) = grow
                end if
            end do
            rr(s%wall) = 0
        end subroutine residual

    end subroutine label_multi_newton

    subroutine newton_matrix(hessian, nrows, grows, dgdx, jac_t)
        !! Transpose of the Newton matrix in (chi, O_1, ..., O_m, mu_1..mu_m),
        !! stored by rows (the layout eu_constrained_solve takes).
        type(csr_t), intent(in) :: hessian
        real(dp), intent(in) :: nrows(:, :), grows(:, :, :), dgdx(:, :)
        type(csr_t), intent(out) :: jac_t
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info
        integer :: n0, m, nt, i, j, k, p, nz, cap
        n0 = size(nrows, 1)
        m = size(nrows, 2)
        nt = n0 + 3*m
        cap = size(hessian%val) + size(dgdx) + 6*n0*m + 4*m
        allocate (rows(cap), cols(cap), vals(cap))
        nz = 0
        do i = 1, n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call add(i, hessian%col(p), hessian%val(p))
            end do
        end do
        do k = 1, 2*m
            do i = 1, nt
                if (abs(dgdx(i, k)) > 0) call add(i, n0 + k, dgdx(i, k))
            end do
            call add(n0 + k, n0 + k, 0.0_dp)
        end do
        do j = 1, m
            do i = 1, n0
                do k = 1, 2
                    if (abs(grows(i, k, j)) > 0) call add(n0 + 2*(j - 1) + k, i, grows(i, k, j))
                end do
                if (abs(nrows(i, j)) > 0) then
                    call add(i, n0 + 2*m + j, nrows(i, j))
                    call add(n0 + 2*m + j, i, nrows(i, j))
                end if
            end do
            call add(n0 + 2*m + j, n0 + 2*m + j, 0.0_dp)
        end do
        call csc_from_triplet(nt, nt, rows(:nz), cols(:nz), vals(:nz), mm, info)
        if (info%code /= 0) error stop 'label_ale_multi: Newton triplets'
        jac_t%nrow = nt
        jac_t%ncol = nt
        jac_t%row_ptr = mm%col_ptr
        jac_t%col = mm%row_idx
        jac_t%val = mm%val
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

    subroutine label_multi_saddle(s, x0, x, status)
        !! Saddle of chi_h near x0, from the quadratic least-squares fit of
        !! chi_h at the quadrature points within 1.5 cells of the estimate
        !! (weights (1 - r^2/rho^2)^2),
        !! recentred until it moves less than 1e-8 cells. chi_h is only C0:
        !! at a vertex its gradient jumps by O(h^p), and Newton on grad chi_h
        !! circles the vertex; the fit is smooth in chi_h and exact for a
        !! quadratic label, so its saddle is within O(h^p) of the exact one.
        !! Status 1: too few points or no saddle (det of the fit >= 0).
        type(label_multi_t), intent(in) :: s
        real(dp), intent(in) :: x0(2)
        real(dp), intent(out) :: x(2)
        integer, intent(out) :: status
        real(dp) :: a(6, 6), rhs(6), phi(6), y(2), v(1), rho, c(2, 2), det, step(2), wk
        integer :: it, p, np, info, ipiv(6)
        status = 1
        x = x0
        rho = 1.5_dp*s%cell_size
        do it = 1, 20
            a = 0
            rhs = 0
            np = 0
            do p = 1, s%quad%point_count
                y = s%quad%s(:, p) - x
                if (norm2(y) > rho) cycle
                call spline_derham_evaluate(s%derham, 0, s%chi, s%quad%block(p), &
                                            s%quad%xi(:, p), 0.0_dp, v)
                y = y/rho
                ! a smooth kernel: the fit is continuous as points enter
                wk = (1 - sum(y**2))**2
                phi = [1.0_dp, y(1), y(2), 0.5_dp*y(1)**2, y(1)*y(2), 0.5_dp*y(2)**2]
                a = a + wk*spread(phi, 2, 6)*spread(phi, 1, 6)
                rhs = rhs + wk*phi*v(1)
                np = np + 1
            end do
            if (np < 12) return
            call dgesv(6, 1, a, 6, ipiv, rhs, 6, info)
            if (info /= 0) return
            c = reshape([rhs(4), rhs(5), rhs(5), rhs(6)], [2, 2])
            det = c(1, 1)*c(2, 2) - c(1, 2)**2
            if (det >= 0) return
            step = -[c(2, 2)*rhs(2) - c(1, 2)*rhs(3), c(1, 1)*rhs(3) - c(2, 1)*rhs(2)]/det*rho
            step = step*min(1.0_dp, s%cell_size/max(norm2(step), tiny(1.0_dp)))
            x = x + step
            if (len_trim(env_debug()) > 0) print '(a,i3,2es14.6,a,es10.2,a,i5)', '  saddle', it, x, &
                ' |step|', norm2(step), ' points', np
            if (norm2(step) < 1e-8_dp*s%cell_size) then
                status = 0
                return
            end if
        end do
    end subroutine label_multi_saddle

    subroutine label_multi_solve(s, profiles, tolerance, xpoint_tolerance, max_passes, &
                                 max_steps, status)
        !! Newton with the cores tied to the O-points; then, with an X-point
        !! vertex, the vertex moved to the fitted saddle of chi_h
        !! (label_multi_saddle) by its blend (the
        !! node values of chi ride with the mesh, the next Newton corrects
        !! them), until the saddle is within xpoint_tolerance cell sizes of
        !! the vertex. Before every later pass the pins are renewed
        !! (repin): each O_j is pinned at the centre of the cell of the
        !! original mesh that holds it, so a core never moves by more than
        !! about half a cell from its place in the original mesh, and chi is
        !! transferred. Status 0 converged, 1 Newton, 2 saddle or mesh,
        !! 4 passes.
        type(label_multi_t), intent(inout) :: s
        procedure(gs_profiles) :: profiles
        real(dp), intent(in) :: tolerance, xpoint_tolerance
        integer, intent(in) :: max_passes, max_steps
        integer, intent(out) :: status
        real(dp) :: saddle(2)
        integer :: pass, steps, st
        status = 4
        do pass = 1, max_passes
            s%passes = pass
            if (pass > 1) then
                call repin(s, st)
                if (st /= 0) then
                    status = 2
                    return
                end if
            end if
            call label_multi_newton(s, profiles, tolerance, max_steps, steps, st)
            if (st == 1) then
                status = 1
                return
            end if
            ! not converged (a step blocked by a fold of the tied mesh):
            ! re-pin and go on
            if (st /= 0) cycle
            if (s%xpoint_dof == 0) then
                status = 0
                return
            end if
            call label_multi_saddle(s, s%xpoint, saddle, st)
            if (st /= 0) then
                status = 2
                return
            end if
            if (norm2(saddle - s%xpoint) < xpoint_tolerance*s%cell_size) then
                status = 0
                return
            end if
            s%xpoint = saddle
        end do
    end subroutine label_multi_solve

    function env_debug() result(v)
        character(8) :: v
        call get_environment_variable('MILCA_LABEL_DEBUG', v)
    end function env_debug

    subroutine repin(s, status)
        !! Pins at the centres of the cells of the original mesh that hold the
        !! O_j (unchanged where O_j left its core), the mesh tied again and
        !! chi transferred to it (Greville interpolation).
        type(label_multi_t), intent(inout) :: s
        integer, intent(out) :: status
        type(spline_derham_t) :: old
        real(dp), allocatable :: chi(:)
        real(dp) :: xi(2), jet(2, 6)
        integer :: j, b, cell, st
        old = s%derham
        do j = 1, size(s%core_block)
            call spline_locate(s%orig, 1, s%axes(:, j), b, xi, st)
            if (st /= 0 .or. b /= s%core_block(j)) cycle
            call spline_cell_at(s%orig, b, xi, cell, s%centre_xi(:, j))
            call spline_geometry_at(s%orig, b, s%centre_xi(:, j), jet, st)
            s%anchor(:, j) = jet(:, 1)
        end do
        call tie_mesh(s, s%axes, s%xpoint, status)
        if (status /= 0) return
        allocate (chi(size(s%chi)))
        call spline_transfer_v0(old, s%chi, s%derham, chi, status)
        if (status /= 0) return
        chi(s%wall) = s%wall_value
        s%chi = chi
    end subroutine repin

end module milca_label_ale_multi
