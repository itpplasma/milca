module milca_eulerian_clebsch
    !! Axisymmetric label lane (docs/eulerian-label-plan.md): energy
    !! W = 2 pi int w dR dZ of the direct Clebsch representation with the
    !! label chi and g in V0 of a spline section, given profiles F(chi),
    !! G(chi), p(chi) and the winding form about a given axis. Returns W, its
    !! gradient in (chi, g) node values and the Hessian, from the generated
    !! leaf clebsch_energy_axisym by the chain rule.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_cell_basis, &
        spline_derham_t_basis
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_clebsch_leaves, only: clebsch_energy_axisym, clebsch_energy_3d, &
        clebsch_energy_3d_ref
    implicit none
    private

    abstract interface
        subroutine clebsch_profiles(chi, f, g, p)
            !! f = (F', F'', F'''), g = (G', G'', G'''), p = (p, p', p'') at chi.
            import :: dp
            real(dp), intent(in) :: chi
            real(dp), intent(out) :: f(3), g(3), p(3)
        end subroutine clebsch_profiles

        subroutine clebsch_exact_field(x, b)
            !! Reference B in cylindrical components at x = (R, phi, Z).
            import :: dp
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: b(3)
        end subroutine clebsch_exact_field
    end interface

    !! Field periods of the 3D systems: the planes span one period, phi =
    !! 2 pi t / clebsch_nfp, so d/dphi = clebsch_nfp d/dt / (2 pi); the
    !! energy is that of the whole torus (one period times clebsch_nfp).
    integer, public :: clebsch_nfp = 1
    public :: clebsch_profiles, clebsch_axisym_system, clebsch_node_gradients, &
              clebsch_3d_system, clebsch_exact_field, clebsch_3d_shape_system

contains

    subroutine clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
        !! Values v0(1, i, q) and physical (R, Z) gradients gn(:, i, q) of the
        !! node functions of a section cell, the gradients from the incidence
        !! grad N_n = sum_e G_en E_e and the inverse transposed Jacobian.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell
        real(dp), allocatable, intent(out) :: v0(:, :, :), gn(:, :, :)
        integer, allocatable, intent(out) :: d0(:)
        real(dp), allocatable :: v1(:, :, :)
        integer, allocatable :: d1(:), s0(:), s1(:)
        real(dp) :: jet(2, 6), jinv(2, 2), det
        integer :: q, p, j, e, kk, i, st

        call spline_derham_cell_basis(derham, quad, cell, 0, v0, d0, s0)
        call spline_derham_cell_basis(derham, quad, cell, 1, v1, d1, s1)
        allocate (gn(2, size(d0), quad%points_per_cell))
        gn = 0
        do q = 1, quad%points_per_cell
            p = quad%point_offset(cell) + q - 1
            do j = 1, size(d1)
                e = d1(j)
                do kk = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                    do i = 1, size(d0)
                        if (d0(i) == derham%grad%colidx(kk)) gn(:, i, q) = gn(:, i, q) + &
                            derham%grad%val(kk)*s1(j)*v1(:, j, q)
                    end do
                end do
            end do
            call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
            det = jet(1, 2)*jet(2, 3) - jet(1, 3)*jet(2, 2)
            jinv = reshape([jet(2, 3), -jet(2, 2), -jet(1, 3), jet(1, 2)], [2, 2])/det
            gn(:, :, q) = matmul(transpose(jinv), gn(:, :, q))
        end do
    end subroutine clebsch_node_gradients

    subroutine clebsch_axisym_system(derham, quad, profiles, axis, chi, g, energy, &
                                     gradient, hessian, angle_shift, frame)
        !! Unknowns ordered (chi(1:n0), g(1:n0)). With angle_shift (V0
        !! coefficients of a single-valued psi) the winding form is
        !! grad theta_X + grad psi_h: a closed form of the same class, chosen
        !! near the straight-field-line angle so that g stays small (the
        !! toroidal field then needs no cancellation between g and the form).
        !! With frame (M, typically H^(1/2) of chi at X) the angle is atan2 of
        !! u = M (x - X): at an elliptic axis it follows the straight-field-line
        !! angle, so g has no O(F') jump at X.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: axis(2), chi(:), g(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: gradient(:)
        type(csr_t), intent(out) :: hessian
        real(dp), intent(in), optional :: angle_shift(:), frame(2, 2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), vals(:)
        integer, allocatable :: d0(:), rows(:), cols(:), idx(:)
        real(dp) :: uf(2)
        real(dp) :: s(2), wq, c, gc(2), gg(2), th(2), rho2, f(3), gp(3), pr(3), w, &
                    gu(5), h(15), hm(5, 5), jet(2, 6)
        real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
        integer :: n0, cell, q, p, i, j, a, b, k, nloc, nz, st
        type(csc_t) :: m
        type(fortsparse_status_t) :: info

        n0 = derham%dimension(0)
        allocate (gradient(2*n0))
        gradient = 0
        energy = 0
        nz = quad%point_count*4*(derham%degree + 1)**4
        allocate (rows(nz), cols(nz), vals(nz))
        nz = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            if (allocated(du)) deallocate (du, idx)
            allocate (du(5, 2*nloc), idx(2*nloc))
            idx = [d0, n0 + d0]
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
                s = jet(:, 1)
                wq = two_pi*quad%weight(p)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                gg = matmul(gn(:, :, q), g(d0))
                if (present(frame)) then
                    uf = matmul(frame, s - axis)
                    th = matmul(transpose(frame), [-uf(2), uf(1)])/(uf(1)**2 + uf(2)**2)
                else
                    rho2 = (s(1) - axis(1))**2 + (s(2) - axis(2))**2
                    th = [-(s(2) - axis(2)), s(1) - axis(1)]/rho2
                end if
                if (present(angle_shift)) th = th + matmul(gn(:, :, q), angle_shift(d0))
                call profiles(c, f, gp, pr)
                call clebsch_energy_axisym(s(1), c, gc(1), gc(2), gg(1), gg(2), th(1), th(2), &
                                           f(1), f(2), f(3), gp(1), gp(2), gp(3), pr(1), &
                                           pr(2), pr(3), w, gu(1), gu(2), gu(3), gu(4), &
                                           gu(5), h(1), h(2), h(3), h(4), h(5), h(6), h(7), &
                                           h(8), h(9), h(10), h(11), h(12), h(13), h(14), &
                                           h(15))
                energy = energy + wq*w
                ! du(a, local): d u_a / d unknown, u = (c, c_R, c_Z, g_R, g_Z)
                du = 0
                du(1, 1:nloc) = v0(1, :, q)
                du(2, 1:nloc) = gn(1, :, q)
                du(3, 1:nloc) = gn(2, :, q)
                du(4, nloc + 1:) = gn(1, :, q)
                du(5, nloc + 1:) = gn(2, :, q)
                gradient(idx) = gradient(idx) + wq*matmul(gu, du)
                k = 0
                do a = 1, 5
                    do b = a, 5
                        k = k + 1
                        hm(a, b) = h(k)
                        hm(b, a) = h(k)
                    end do
                end do
                do j = 1, 2*nloc
                    do i = 1, 2*nloc
                        nz = nz + 1
                        rows(nz) = idx(i)
                        cols(nz) = idx(j)
                        vals(nz) = wq*dot_product(du(:, i), matmul(hm, du(:, j)))
                    end do
                end do
            end do
        end do
        call csc_from_triplet(2*n0, 2*n0, cols(:nz), rows(:nz), vals(:nz), m, info)
        if (info%code /= 0) error stop 'Clebsch Hessian triplets'
        hessian%nrow = 2*n0
        hessian%ncol = 2*n0
        hessian%row_ptr = m%col_ptr
        hessian%col = m%row_idx
        hessian%val = m%val
    end subroutine clebsch_axisym_system


    subroutine clebsch_3d_system(derham, quad, nt, profiles, axis, chi, g, energy, &
                                 gradient, hessian, plane_axis, plane_nodes, plane_frames, &
                                 exact_field, field_error)
        !! 3D label lane on the section x trigonometric-nodal complex (nt =
        !! 2M + 1 planes, cylindrical chart, phi = 2 pi t). Unknowns chi and g
        !! as node values per plane, section index fastest, ordered (chi, g).
        !! Physical gradient components (R, phi, Z): section gradients times
        !! L_k(t) and N_i L_k'(t)/(2 pi R). The winding form is about the
        !! axis X (independent of phi). Per section cell the local Hessian is
        !! accumulated densely over all quadrature points and planes.
        !! With plane_nodes (V0 node values of (R, Z) on each plane, the same
        !! section topology) the section varies with t: s(xi, t) = sum_k
        !! L_k(t) s_k(xi), so a wall that changes with phi and a core that
        !! follows the axis are both isoparametric. With plane_axis (the
        !! axis X_k on each plane) the winding form is about the axis curve
        !! X(t) = sum_k L_k(t) X_k. The quadrature stays that of derham.
        !! With plane_frames (M_k, typically H^(1/2) of chi at X_k) the angle
        !! is theta = atan2 of u = M(t)(x - X(t)), M(t) = sum_k L_k(t) M_k: at
        !! an elliptic axis it follows the straight-field-line angle, so g
        !! stays regular there.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: nt
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: axis(2), chi(:), g(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), intent(in), optional :: plane_axis(:, :), plane_nodes(:, :, :), &
                                           plane_frames(:, :, :)
        procedure(clebsch_exact_field), optional :: exact_field
        !! sum of w |B - B_exact|^2 and of w |B_exact|^2 over the quadrature
        real(dp), intent(out), optional :: field_error(2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), gloc(:), hloc(:, :), &
                                 vals(:)
        real(dp) :: dt(2), xt(2), sm(2), fm(2, 2), ft(2, 2), detf, gx(2), v(2), rt(2), &
                    sk(2, nt), gk(2, 2, nt), mf(2, 2), mdot(2, 2), uu(2), udot(2), bv(3), &
                    bx(3)
        integer :: kk, ka, kb, la, lb, il, jl, blk
        real(dp) :: ferr(2)
        integer, allocatable :: d0(:), idx(:), rows(:), cols(:)
        real(dp), allocatable :: sx(:, :, :), am(:, :, :, :, :), tb(:, :, :)
        real(dp) :: lt(nt, 2), ht(nt, 2), s(2), jet(2, 6), wq, wt, c, gc(3), gg(3), hv(3), &
                    rho2, f(3), gp(3), pr(3), w, gu(7), h(28), hm(7, 7), t
        real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
        integer :: n0, ntot, nqt, cell, q, p, m, i, k, a, b, nloc, nl, st, nz, col
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info

        n0 = derham%dimension(0)
        ntot = n0*nt
        nqt = 2*nt + 2
        allocate (gradient(2*ntot))
        gradient = 0
        energy = 0
        if (present(field_error)) field_error = 0
        ! one block of triplets per cell (zeros where a cell has fewer nodes)
        blk = (2*(derham%degree + 1)**2*nt)**2
        if (present(hessian)) then
            allocate (rows(size(quad%point_offset)*blk), cols(size(quad%point_offset)*blk), &
                      vals(size(quad%point_offset)*blk))
            rows = 1
            cols = 1
            vals = 0
        end if
        ferr = 0
        !$omp parallel do schedule(dynamic) default(private) &
        !$omp shared(derham, quad, nt, chi, g, axis, plane_axis, plane_nodes, plane_frames, &
        !$omp        gradient, rows, cols, vals, n0, ntot, nqt, blk, hessian, field_error, &
        !$omp        clebsch_nfp) &
        !$omp reduction(+:energy, ferr)
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            nl = nloc*nt
            if (allocated(du)) deallocate (du, gloc, hloc, idx, sx, am, tb)
            allocate (du(7, 2*nl), gloc(2*nl), hloc(2*nl, 2*nl), idx(2*nl), &
                      sx(7, 2*nloc, 2), am(2*nloc, 2*nloc, 2, 2, nqt), tb(nt, 2, nqt))
            am = 0
            do k = 1, nt
                idx((k - 1)*nloc + 1:k*nloc) = d0 + n0*(k - 1)
                idx(nl + (k - 1)*nloc + 1:nl + k*nloc) = ntot + d0 + n0*(k - 1)
            end do
            gloc = 0
            hloc = 0
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
                s = jet(:, 1)
                wq = quad%weight(p)
                if (present(plane_nodes)) then
                    do kk = 1, nt
                        sk(:, kk) = matmul(plane_nodes(:, d0, kk), v0(1, :, q))
                        gk(:, :, kk) = matmul(plane_nodes(:, d0, kk), transpose(gn(:, :, q)))
                    end do
                end if
                do m = 1, nqt
                    t = (m - 1)/real(nqt, dp)
                    call spline_derham_t_basis(nt, t, lt, ht)
                    lt(:, 2) = clebsch_nfp*lt(:, 2)
                    ! the point, F = ds/ds0, the mesh velocity ds/dt and the axis
                    sm = s
                    fm = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
                    v = 0
                    if (present(plane_nodes)) then
                        sm = matmul(sk, lt(:, 1))
                        fm = 0
                        do kk = 1, nt
                            fm = fm + lt(kk, 1)*gk(:, :, kk)
                        end do
                        v = matmul(sk, lt(:, 2))
                    end if
                    xt = axis
                    dt = 0
                    if (present(plane_axis)) then
                        xt = matmul(plane_axis, lt(:, 1))
                        dt = matmul(plane_axis, lt(:, 2))
                    end if
                    detf = fm(1, 1)*fm(2, 2) - fm(1, 2)*fm(2, 1)
                    ft = reshape([fm(2, 2), -fm(1, 2), -fm(2, 1), fm(1, 1)], [2, 2])/detf
                    rt = sm - xt
                    mf = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
                    mdot = 0
                    if (present(plane_frames)) then
                        mf = 0
                        do kk = 1, nt
                            mf = mf + lt(kk, 1)*plane_frames(:, :, kk)
                            mdot = mdot + lt(kk, 2)*plane_frames(:, :, kk)
                        end do
                    end if
                    ! u = M (x - X); d/dt at fixed x: u' = M' (x - X) - M X'
                    uu = matmul(mf, rt)
                    udot = matmul(mdot, rt) - matmul(mf, dt)
                    rho2 = uu(1)**2 + uu(2)**2
                    gx = matmul(transpose(mf), [-uu(2), uu(1)])/rho2
                    hv = [gx(1), (uu(1)*udot(2) - uu(2)*udot(1))/(rho2*two_pi*sm(1)), gx(2)]
                    du = 0
                    do k = 1, nt
                        do i = 1, nloc
                            a = (k - 1)*nloc + i
                            ! F^-T grad_0; d/dphi at fixed (R, Z) removes the mesh velocity
                            gx = matmul(ft, gn(:, i, q))
                            du(1, a) = v0(1, i, q)*lt(k, 1)
                            du(2, a) = gx(1)*lt(k, 1)
                            du(3, a) = (v0(1, i, q)*lt(k, 2) - dot_product(v, gx)*lt(k, 1)) &
                                       /(two_pi*sm(1))
                            du(4, a) = gx(2)*lt(k, 1)
                            du(5, nl + a) = du(2, a)
                            du(6, nl + a) = du(3, a)
                            du(7, nl + a) = du(4, a)
                        end do
                    end do
                    c = dot_product(du(1, 1:nl), chi(idx(1:nl)))
                    gc = matmul(du(2:4, 1:nl), chi(idx(1:nl)))
                    gg = matmul(du(5:7, nl + 1:), g(idx(nl + 1:) - ntot))
                    call profiles(c, f, gp, pr)
                    call clebsch_energy_3d(sm(1), c, gc(1), gc(2), gc(3), gg(1), gg(2), gg(3), &
                                           hv(1), hv(2), hv(3), f(1), f(2), f(3), gp(1), &
                                           gp(2), gp(3), pr(1), pr(2), pr(3), w, gu(1), gu(2), &
                                           gu(3), gu(4), gu(5), gu(6), gu(7), h(1), h(2), &
                                           h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10), &
                                           h(11), h(12), h(13), h(14), h(15), h(16), h(17), &
                                           h(18), h(19), h(20), h(21), h(22), h(23), h(24), &
                                           h(25), h(26), h(27), h(28))
                    wt = wq*abs(detf)*two_pi/nqt
                    energy = energy + wt*w
                    if (present(exact_field)) then
                        ! B = grad g x grad chi + F' grad chi x h + G' grad chi x grad phi
                        bv = cross3(gg, gc) + f(1)*cross3(gc, hv) &
                             + gp(1)*cross3(gc, [0.0_dp, 1/sm(1), 0.0_dp])
                        call exact_field([sm(1), two_pi*t/clebsch_nfp, sm(2)], bx)
                        ferr = ferr + wt*[sum((bv - bx)**2), sum(bx**2)]
                    end if
                    gloc = gloc + wt*matmul(gu, du)
                    k = 0
                    do a = 1, 7
                        do b = a, 7
                            k = k + 1
                            hm(a, b) = h(k)
                            hm(b, a) = h(k)
                        end do
                    end do
                    if (present(hessian)) then
                        ! du = sum_c sx(:, :, c) (x) T_c(t) with T_1 = L_k, T_2 = L_k':
                        ! accumulate the section factors per toroidal point, expand
                        ! over the planes once per cell
                        sx = 0
                        do i = 1, nloc
                            gx = matmul(ft, gn(:, i, q))
                            sx(1, i, 1) = v0(1, i, q)
                            sx(2, i, 1) = gx(1)
                            sx(3, i, 1) = -dot_product(v, gx)/(two_pi*sm(1))
                            sx(3, i, 2) = v0(1, i, q)/(two_pi*sm(1))
                            sx(4, i, 1) = gx(2)
                            sx(5:7, nloc + i, :) = sx(2:4, i, :)
                        end do
                        tb(:, :, m) = lt
                        do b = 1, 2
                            do a = 1, 2
                                am(:, :, a, b, m) = am(:, :, a, b, m) + &
                                    wt*matmul(transpose(sx(:, :, a)), matmul(hm, sx(:, :, b)))
                            end do
                        end do
                    end if
                end do
            end do
            if (present(hessian)) then
                do m = 1, nqt
                    do b = 1, 2
                        do a = 1, 2
                            do kb = 1, 2*nloc
                                do lb = 1, nt
                                    jl = merge((lb - 1)*nloc + kb, nl + (lb - 1)*nloc + kb - nloc, kb <= nloc)
                                    do ka = 1, 2*nloc
                                        if (am(ka, kb, a, b, m) == 0) cycle
                                        do la = 1, nt
                                            il = merge((la - 1)*nloc + ka, nl + (la - 1)*nloc + ka - nloc, ka <= nloc)
                                            hloc(il, jl) = hloc(il, jl) + tb(la, a, m)*tb(lb, b, m) &
                                                           *am(ka, kb, a, b, m)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            do i = 1, 2*nl
                !$omp atomic
                gradient(idx(i)) = gradient(idx(i)) + gloc(i)
            end do
            if (.not. present(hessian)) cycle
            nz = (cell - 1)*blk
            do col = 1, 2*nl
                do i = 1, 2*nl
                    nz = nz + 1
                    rows(nz) = idx(i)
                    cols(nz) = idx(col)
                    vals(nz) = hloc(i, col)
                end do
            end do
        end do
        !$omp end parallel do
        if (present(field_error)) field_error = ferr
        if (.not. present(hessian)) return
        call csc_from_triplet(2*ntot, 2*ntot, cols, rows, vals, mm, info)
        if (info%code /= 0) error stop 'Clebsch 3D Hessian triplets'
        hessian%nrow = 2*ntot
        hessian%ncol = 2*ntot
        hessian%row_ptr = mm%col_ptr
        hessian%col = mm%row_idx
        hessian%val = mm%val
    end subroutine clebsch_3d_system

    subroutine clebsch_3d_shape_system(derham, quad, nt, profiles, chi, g, plane_axis, &
                                       plane_nodes, plane_frames, energy, gradient, hessian, &
                                       shape)
        !! The 3D label lane of clebsch_3d_system (with plane_axis, plane_nodes
        !! and plane_frames) from the reference-coordinate leaf
        !! clebsch_energy_3d_ref, and its exact derivatives with respect to the
        !! geometry: shape(i, j) = d gradient_i / d gamma_j with gamma the node
        !! coordinates of every plane, column ((k - 1) n0 + i - 1) 2 + c for
        !! component c of node i on plane k, then the axis points, column
        !! 2 n0 nt + 2 (k - 1) + c. The field locals (value, reference gradient,
        !! nodal t derivative) do not depend on the geometry, so the shape
        !! derivative is the mixed second derivative of the leaf.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: nt
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: chi(:), g(:), plane_axis(:, :), plane_nodes(:, :, :), &
                                plane_frames(:, :, :)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: gradient(:)
        type(csr_t), intent(out) :: hessian, shape
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), sx(:, :, :), sy(:, :, :), &
                                 am(:, :, :, :, :), as(:, :, :, :, :), ax(:, :, :, :, :), &
                                 tb(:, :, :), gloc(:, :), hv(:), sv(:), hvl(:, :), svl(:, :)
        integer, allocatable :: d0(:), rh(:), ch(:), rs(:), cs(:)
        real(dp) :: lt(nt, 2), ht(nt, 2), lv(7), qv(12), mf(2, 2), md(2, 2), f(3), gp(3), &
                    pr(3), w, gl(7), hh(28), hm(84), hmat(7, 7), gmat(7, 12), t, wt, s(2), &
                    fm(2, 2), vel(2), xa(2), xd(2), ck(nt), gk(nt), cak(2, nt), gak(2, nt)
        real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
        integer :: n0, ntot, nqt, cell, q, p, m, k, i, a, b, nloc, nl2, nc2, ct1, ct2, &
                   col, nz_h, nz_s, blk_h, blk_s, c, ncol, ki, kj, idx_col
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info

        n0 = derham%dimension(0)
        ntot = n0*nt
        nqt = 2*nt + 2
        ncol = 2*n0*nt + 2*nt
        allocate (gradient(2*ntot))
        gradient = 0
        energy = 0
        blk_h = (2*(derham%degree + 1)**2*nt)**2
        blk_s = (2*(derham%degree + 1)**2*nt)*(2*(derham%degree + 1)**2*nt + 2*nt)
        allocate (rh(size(quad%point_offset)*blk_h), ch(size(quad%point_offset)*blk_h), &
                  hv(size(quad%point_offset)*blk_h), rs(size(quad%point_offset)*blk_s), &
                  cs(size(quad%point_offset)*blk_s), sv(size(quad%point_offset)*blk_s))
        rh = 1
        ch = 1
        hv = 0
        rs = 1
        cs = 1
        sv = 0
        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(cell, v0, gn, d0, nloc, nl2, nc2, sx, sy, am, as, ax, tb, gloc, hvl, svl, &
        !$omp         m, t, lt, ht, q, p, k, ck, gk, cak, gak, lv, s, fm, vel, xa, xd, mf, md, qv, &
        !$omp         f, gp, pr, w, gl, hh, hm, hmat, gmat, wt, i, a, b, c, ct1, ct2, col, &
        !$omp         ki, kj, nz_h, nz_s, idx_col) reduction(+:energy)
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            nl2 = 2*nloc
            nc2 = 2*nloc
            allocate (sx(7, nl2, 2), sy(12, nc2, 2), am(nl2, nl2, 2, 2, nqt), &
                      as(nl2, nc2, 2, 2, nqt), ax(nl2, 2, 2, 2, nqt), tb(nt, 2, nqt), &
                      gloc(nt, nl2))
            am = 0
            as = 0
            ax = 0
            gloc = 0
            do m = 1, nqt
                t = (m - 1)/real(nqt, dp)
                call spline_derham_t_basis(nt, t, lt, ht)
                lt(:, 2) = clebsch_nfp*lt(:, 2)
                tb(:, :, m) = lt
                xa = matmul(plane_axis, lt(:, 1))
                xd = matmul(plane_axis, lt(:, 2))
                mf = 0
                md = 0
                do k = 1, nt
                    mf = mf + lt(k, 1)*plane_frames(:, :, k)
                    md = md + lt(k, 2)*plane_frames(:, :, k)
                end do
                do q = 1, quad%points_per_cell
                    p = quad%point_offset(cell) + q - 1
                    ! field and geometry locals per plane, then in t
                    do k = 1, nt
                        ck(k) = dot_product(v0(1, :, q), chi(d0 + n0*(k - 1)))
                        gk(k) = dot_product(v0(1, :, q), g(d0 + n0*(k - 1)))
                        cak(:, k) = matmul(gn(:, :, q), chi(d0 + n0*(k - 1)))
                        gak(:, k) = matmul(gn(:, :, q), g(d0 + n0*(k - 1)))
                    end do
                    lv = [dot_product(lt(:, 1), ck), dot_product(lt(:, 1), cak(1, :)), &
                          dot_product(lt(:, 1), cak(2, :)), dot_product(lt(:, 2), ck), &
                          dot_product(lt(:, 1), gak(1, :)), dot_product(lt(:, 1), gak(2, :)), &
                          dot_product(lt(:, 2), gk)]
                    s = 0
                    fm = 0
                    vel = 0
                    do k = 1, nt
                        s = s + lt(k, 1)*matmul(plane_nodes(:, d0, k), v0(1, :, q))
                        fm = fm + lt(k, 1)*matmul(plane_nodes(:, d0, k), transpose(gn(:, :, q)))
                        vel = vel + lt(k, 2)*matmul(plane_nodes(:, d0, k), v0(1, :, q))
                    end do
                    qv = [s(1), s(2), fm(1, 1), fm(2, 1), fm(1, 2), fm(2, 2), vel(1), vel(2), &
                          xa(1), xa(2), xd(1), xd(2)]
                    call profiles(lv(1), f, gp, pr)
                    call clebsch_energy_3d_ref(lv(1), lv(2), lv(3), lv(4), lv(5), lv(6), lv(7), &
                        qv(1), qv(2), qv(3), qv(4), qv(5), qv(6), qv(7), qv(8), qv(9), qv(10), &
                        qv(11), qv(12), mf(1, 1), mf(2, 1), mf(1, 2), mf(2, 2), md(1, 1), &
                        md(2, 1), md(1, 2), md(2, 2), f(1), f(2), f(3), gp(1), gp(2), gp(3), &
                        pr(1), pr(2), pr(3), w, gl(1), gl(2), gl(3), gl(4), gl(5), gl(6), &
                        gl(7), hh(1), hh(2), hh(3), hh(4), hh(5), hh(6), hh(7), hh(8), hh(9), &
                        hh(10), hh(11), hh(12), hh(13), hh(14), hh(15), hh(16), hh(17), hh(18), &
                        hh(19), hh(20), hh(21), hh(22), hh(23), hh(24), hh(25), hh(26), hh(27), &
                        hh(28), hm(1), hm(2), hm(3), hm(4), hm(5), hm(6), hm(7), hm(8), hm(9), &
                        hm(10), hm(11), hm(12), hm(13), hm(14), hm(15), hm(16), hm(17), hm(18), &
                        hm(19), hm(20), hm(21), hm(22), hm(23), hm(24), hm(25), hm(26), hm(27), &
                        hm(28), hm(29), hm(30), hm(31), hm(32), hm(33), hm(34), hm(35), hm(36), &
                        hm(37), hm(38), hm(39), hm(40), hm(41), hm(42), hm(43), hm(44), hm(45), &
                        hm(46), hm(47), hm(48), hm(49), hm(50), hm(51), hm(52), hm(53), hm(54), &
                        hm(55), hm(56), hm(57), hm(58), hm(59), hm(60), hm(61), hm(62), hm(63), &
                        hm(64), hm(65), hm(66), hm(67), hm(68), hm(69), hm(70), hm(71), hm(72), &
                        hm(73), hm(74), hm(75), hm(76), hm(77), hm(78), hm(79), hm(80), hm(81), &
                        hm(82), hm(83), hm(84))
                    wt = quad%weight(p)/nqt
                    energy = energy + wt*w
                    k = 0
                    do a = 1, 7
                        do b = a, 7
                            k = k + 1
                            hmat(a, b) = hh(k)
                            hmat(b, a) = hh(k)
                        end do
                    end do
                    gmat = transpose(reshape(hm, [12, 7]))
                    ! section factors of the field locals: column 1 with L_k,
                    ! column 2 with L_k'
                    sx = 0
                    sy = 0
                    do i = 1, nloc
                        sx(1, i, 1) = v0(1, i, q)
                        sx(2, i, 1) = gn(1, i, q)
                        sx(3, i, 1) = gn(2, i, q)
                        sx(4, i, 2) = v0(1, i, q)
                        sx(5, nloc + i, 1) = gn(1, i, q)
                        sx(6, nloc + i, 1) = gn(2, i, q)
                        sx(7, nloc + i, 2) = v0(1, i, q)
                        do c = 1, 2
                            ! node coordinate (i, c): s_c, f_c1, f_c2 with L_k, v_c with L_k'
                            sy(c, 2*(i - 1) + c, 1) = v0(1, i, q)
                            sy(2 + c, 2*(i - 1) + c, 1) = gn(1, i, q)
                            sy(4 + c, 2*(i - 1) + c, 1) = gn(2, i, q)
                            sy(6 + c, 2*(i - 1) + c, 2) = v0(1, i, q)
                        end do
                    end do
                    do ct1 = 1, 2
                        do i = 1, nl2
                            do k = 1, nt
                                gloc(k, i) = gloc(k, i) + wt*dot_product(gl, sx(:, i, ct1))*lt(k, ct1)
                            end do
                        end do
                        do ct2 = 1, 2
                            am(:, :, ct1, ct2, m) = am(:, :, ct1, ct2, m) + &
                                wt*matmul(transpose(sx(:, :, ct1)), matmul(hmat, sx(:, :, ct2)))
                            as(:, :, ct1, ct2, m) = as(:, :, ct1, ct2, m) + &
                                wt*matmul(transpose(sx(:, :, ct1)), matmul(gmat, sy(:, :, ct2)))
                        end do
                        ! axis: x_c with L_k, xd_c with L_k'
                        do c = 1, 2
                            ax(:, c, ct1, 1, m) = ax(:, c, ct1, 1, m) + &
                                wt*matmul(transpose(sx(:, :, ct1)), gmat(:, 8 + c))
                            ax(:, c, ct1, 2, m) = ax(:, c, ct1, 2, m) + &
                                wt*matmul(transpose(sx(:, :, ct1)), gmat(:, 10 + c))
                        end do
                    end do
                end do
            end do
            ! expand over the planes: rows (plane ki, local l), columns
            allocate (hvl(nl2*nt, nl2*nt), svl(nl2*nt, nc2*nt + 2*nt))
            hvl = 0
            svl = 0
            do m = 1, nqt
                do ct1 = 1, 2
                    do ct2 = 1, 2
                        do kj = 1, nt
                            do ki = 1, nt
                                hvl((ki - 1)*nl2 + 1:ki*nl2, (kj - 1)*nl2 + 1:kj*nl2) = &
                                    hvl((ki - 1)*nl2 + 1:ki*nl2, (kj - 1)*nl2 + 1:kj*nl2) + &
                                    tb(ki, ct1, m)*tb(kj, ct2, m)*am(:, :, ct1, ct2, m)
                                svl((ki - 1)*nl2 + 1:ki*nl2, (kj - 1)*nc2 + 1:kj*nc2) = &
                                    svl((ki - 1)*nl2 + 1:ki*nl2, (kj - 1)*nc2 + 1:kj*nc2) + &
                                    tb(ki, ct1, m)*tb(kj, ct2, m)*as(:, :, ct1, ct2, m)
                                do c = 1, 2
                                    svl((ki - 1)*nl2 + 1:ki*nl2, nc2*nt + 2*(kj - 1) + c) = &
                                        svl((ki - 1)*nl2 + 1:ki*nl2, nc2*nt + 2*(kj - 1) + c) + &
                                        tb(ki, ct1, m)*tb(kj, ct2, m)*ax(:, c, ct1, ct2, m)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            do ki = 1, nt
                do i = 1, nl2
                    !$omp atomic
                    gradient(global_row(i, ki, nloc, d0, n0, ntot)) = gradient(global_row(i, ki, nloc, d0, n0, ntot)) + gloc(ki, i)
                end do
            end do
            nz_h = (cell - 1)*blk_h
            do kj = 1, nt
                do b = 1, nl2
                    do ki = 1, nt
                        do a = 1, nl2
                            nz_h = nz_h + 1
                            rh(nz_h) = global_row(a, ki, nloc, d0, n0, ntot)
                            ch(nz_h) = global_row(b, kj, nloc, d0, n0, ntot)
                            hv(nz_h) = hvl((ki - 1)*nl2 + a, (kj - 1)*nl2 + b)
                        end do
                    end do
                end do
            end do
            nz_s = (cell - 1)*blk_s
            do col = 1, nc2*nt + 2*nt
                if (col <= nc2*nt) then
                    kj = (col - 1)/nc2 + 1
                    i = mod(col - 1, nc2) + 1
                    idx_col = ((kj - 1)*n0 + d0((i - 1)/2 + 1) - 1)*2 + mod(i - 1, 2) + 1
                else
                    idx_col = 2*n0*nt + col - nc2*nt
                end if
                do ki = 1, nt
                    do a = 1, nl2
                        nz_s = nz_s + 1
                        rs(nz_s) = global_row(a, ki, nloc, d0, n0, ntot)
                        cs(nz_s) = idx_col
                        sv(nz_s) = svl((ki - 1)*nl2 + a, col)
                    end do
                end do
            end do
            deallocate (sx, sy, am, as, ax, tb, gloc, hvl, svl)
        end do
        !$omp end parallel do
        call csc_from_triplet(2*ntot, 2*ntot, ch, rh, hv, mm, info)
        if (info%code /= 0) error stop 'Clebsch 3D shape: Hessian triplets'
        hessian%nrow = 2*ntot
        hessian%ncol = 2*ntot
        hessian%row_ptr = mm%col_ptr
        hessian%col = mm%row_idx
        hessian%val = mm%val
        ! shape(i, j): stored by rows (CSR) from the transposed CSC
        call csc_from_triplet(ncol, 2*ntot, cs, rs, sv, mm, info)
        if (info%code /= 0) error stop 'Clebsch 3D shape: shape triplets'
        shape%nrow = 2*ntot
        shape%ncol = ncol
        shape%row_ptr = mm%col_ptr
        shape%col = mm%row_idx
        shape%val = mm%val
    end subroutine clebsch_3d_shape_system

    pure integer function global_row(l, k, nloc, d0, n0, ntot)
        !! Global unknown of local field dof l (chi 1..nloc, g nloc + 1..2 nloc)
        !! on plane k.
        integer, intent(in) :: l, k, nloc, d0(:), n0, ntot

        if (l <= nloc) then
            global_row = d0(l) + n0*(k - 1)
        else
            global_row = ntot + d0(l - nloc) + n0*(k - 1)
        end if
    end function global_row

    pure function cross3(a, b) result(c)
        !! Cross product in right-handed (R, phi, Z) components.
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

end module milca_eulerian_clebsch
