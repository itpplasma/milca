module milca_flux_map_3d
    !! Cartesian flux map in 3D (docs/flux-map.md): u = (x, y) maps the torus
    !! onto the unit disc on every toroidal plane, s = x^2 + y^2 is the flux
    !! label and
    !!
    !!   B = -[tau(s) grad x X grad y + psi'(s) grad phi X grad s],
    !!   W = int [|B|^2/2 - p(s)] dV   (whole torus).
    !!
    !! Discretisation as in clebsch_3d_system: a spline section times the
    !! trigonometric nodal basis of nt = 2M + 1 planes over one field period
    !! (clebsch_nfp), phi = 2 pi t/clebsch_nfp, the section of plane k given by
    !! its node coordinates plane_nodes(:, :, k) (isoparametric in t). The mesh
    !! knows nothing about the magnetic axis; there are no axis unknowns,
    !! winding forms or singular quadratures. On the wall u = (cos Theta, sin
    !! Theta) holds at every point and every t: with sigma the sum of the wall
    !! functions, the wall part of u is sigma (cos, sin)(Theta_h), Theta_h =
    !! sum_k L_k(t) sum_b N_b Theta_bk/sigma.
    !!
    !! Unknowns z(1:2 n0 nt), section index fastest: z((k-1) n0 + i) = x_ik
    !! (interior function) or Theta_ik (wall function), z(n0 nt + (k-1) n0 + i)
    !! = y_ik (interior; unused and fixed for wall functions). Kernel:
    !! generated leaf flux_map_energy_3d. B enters pointwise: every per-cell
    !! projection tried in 3D was unstable (docs/flux-map.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_spline_eval, only: spline_quadrature_t
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_t_basis
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_clebsch, only: clebsch_node_gradients, clebsch_nfp, clebsch_exact_field
    use milca_eulerian_vacuum, only: eu_constrained_t, eu_constrained_refactor, &
        eu_constrained_apply, eu_constrained_free, eu_constrained_inertia
    use milca_flux_map, only: flux_map_profiles
    use milca_flux_map_leaves, only: flux_map_energy_3d
    implicit none
    private

    type, public :: flux_map3d_t
        integer :: n0 = 0, nt = 0
        logical, allocatable :: is_wall(:) ! section functions with a wall trace
        integer, allocatable :: fixed(:) ! y slots of the wall functions
        real(dp), allocatable :: z(:)
        !! toroidal quadrature points (trapezoid) per plane
        integer :: toroidal_points = 2
        !! Gauge theta* -> theta* + k(s) by a penalty on linear rows, exact for
        !! stellarator-symmetric states with rows y_h = 0 on Z = 0 of the
        !! plane phi = 0. Without rows the solver's shift holds the gauge; the
        !! general case needs constraints orthogonal to the gauge.
        integer, allocatable :: gauge_ptr(:), gauge_col(:)
        real(dp), allocatable :: gauge_val(:), gauge_target(:)
        real(dp) :: gauge_weight = 0
        !! Solver: .false. minimises W (steps accepted on energy decrease);
        !! .true. seeks the stationary point next to the start by Newton
        !! with a line search on the scaled gradient norm, whatever the
        !! signature of the Hessian (an equilibrium need not be a minimum).
        logical :: stationary = .false.
        integer :: steps = 0
        real(dp) :: energy = 0, gradient_norm = 0
        logical :: converged = .false., verbose = .false.
        !! Newton steps by the symmetric LDL^T factorisation of the upper
        !! triangle of the Hessian (FortSparse); .false.: LU of the mirrored
        !! full matrix (same triangle assembly), which gives no inertia.
        logical :: symmetric_solver = .true.
        !! Inertia of the free block of the Hessian (gauge penalty included)
        !! at the converged state: numbers of negative and zero eigenvalues
        !! (-1: not computed, i.e. not converged or LU), and whether every
        !! pivot was clear of its rounding threshold.
        integer :: negative_eigenvalues = -1, zero_eigenvalues = -1
        logical :: inertia_reliable = .false.
    end type flux_map3d_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    public :: flux_map3d_init, flux_map3d_gauge_row, flux_map3d_system, flux_map3d_solve, &
              flux_map3d_point, flux_map3d_expand_symmetric

contains

    subroutine flux_map3d_init(fm, n0, nt, wall)
        !! wall: the section V0 functions with a trace on the wall.
        type(flux_map3d_t), intent(out) :: fm
        integer, intent(in) :: n0, nt, wall(:)
        integer :: k

        fm%n0 = n0
        fm%nt = nt
        allocate (fm%is_wall(n0), source=.false.)
        fm%is_wall(wall) = .true.
        fm%fixed = [(n0*nt + (k - 1)*n0 + wall, k=1, nt)]
        allocate (fm%z(2*n0*nt), source=0.0_dp)
        allocate (fm%gauge_ptr(1), source=1)
        allocate (fm%gauge_col(0), fm%gauge_val(0), fm%gauge_target(0))
    end subroutine flux_map3d_init

    subroutine flux_map3d_gauge_row(fm, cols, vals, target)
        type(flux_map3d_t), intent(inout) :: fm
        integer, intent(in) :: cols(:)
        real(dp), intent(in) :: vals(:), target

        fm%gauge_col = [fm%gauge_col, cols]
        fm%gauge_val = [fm%gauge_val, vals]
        fm%gauge_target = [fm%gauge_target, target]
        fm%gauge_ptr = [fm%gauge_ptr, size(fm%gauge_col) + 1]
    end subroutine flux_map3d_gauge_row

    pure subroutine locals(is_w, ix, iy, nv, gx, vel, rr, ax, ay, bx, by, u, sx, sig, th, &
                           gsig, gth, lam, mw)
        !! u = (x, y, grad x, grad y) with gradient components (R, phi, Z) at
        !! a point, from the amplitudes ax, ay (values in t) and bx, by (t
        !! derivatives times clebsch_nfp) of the section functions with values
        !! nv and (R, Z) gradients gx; vel the mesh velocity, rr = R. sx(:, a,
        !! 1:2): derivatives of u with respect to the amplitude and the t
        !! derivative of slot a. For the second derivatives of the wall part:
        !! sig, th, gsig, gth (3 components) and per function lam, mw.
        logical, intent(in) :: is_w(:)
        integer, intent(in) :: ix(:), iy(:)
        real(dp), intent(in) :: nv(:), gx(:, :), vel(2), rr, ax(:), ay(:), bx(:), by(:)
        real(dp), intent(out) :: u(8), sx(:, :, :), sig, th, gsig(3), gth(3), lam(:), mw(:, :)
        real(dp) :: c, g3(3), acc, gacc(3), tacc, cs, sn
        integer :: i, a

        c = 1/(two_pi*rr)
        u = 0
        sx = 0
        sig = 0
        gsig = 0
        acc = 0
        gacc = 0
        tacc = 0
        th = 0
        gth = 0
        lam = 0
        mw = 0
        do i = 1, size(nv)
            g3 = [gx(1, i), -c*dot_product(vel, gx(:, i)), gx(2, i)]
            if (is_w(i)) then
                sig = sig + nv(i)
                gsig = gsig + g3
                acc = acc + nv(i)*ax(i)
                gacc = gacc + g3*ax(i)
                tacc = tacc + nv(i)*bx(i)
            else
                u(1) = u(1) + nv(i)*ax(i)
                u(2) = u(2) + nv(i)*ay(i)
                u(3:5) = u(3:5) + g3*ax(i) + [0.0_dp, c*nv(i)*bx(i), 0.0_dp]
                u(6:8) = u(6:8) + g3*ay(i) + [0.0_dp, c*nv(i)*by(i), 0.0_dp]
                sx(1, ix(i), 1) = nv(i)
                sx(3:5, ix(i), 1) = g3
                sx(4, ix(i), 2) = c*nv(i)
                sx(2, iy(i), 1) = nv(i)
                sx(6:8, iy(i), 1) = g3
                sx(7, iy(i), 2) = c*nv(i)
            end if
        end do
        if (sig <= 0) return
        th = acc/sig
        gth = (gacc - th*gsig)/sig
        gth(2) = gth(2) + c*tacc/sig
        cs = cos(th)
        sn = sin(th)
        u(1) = u(1) + sig*cs
        u(2) = u(2) + sig*sn
        u(3:5) = u(3:5) + gsig*cs - sig*sn*gth
        u(6:8) = u(6:8) + gsig*sn + sig*cs*gth
        do i = 1, size(nv)
            if (.not. is_w(i)) cycle
            g3 = [gx(1, i), -c*dot_product(vel, gx(:, i)), gx(2, i)]
            lam(i) = nv(i)/sig
            mw(:, i) = (g3 - lam(i)*gsig)/sig
            a = ix(i)
            sx(1, a, 1) = -sig*sn*lam(i)
            sx(2, a, 1) = sig*cs*lam(i)
            sx(3:5, a, 1) = -gsig*sn*lam(i) - sig*cs*lam(i)*gth - sig*sn*mw(:, i)
            sx(6:8, a, 1) = gsig*cs*lam(i) - sig*sn*lam(i)*gth + sig*cs*mw(:, i)
            sx(4, a, 2) = -sig*sn*lam(i)*c
            sx(7, a, 2) = sig*cs*lam(i)*c
        end do
    end subroutine locals

    subroutine flux_map3d_system(fm, derham, quad, profiles, plane_nodes, z, energy, gradient, &
                                 hessian, exact_field, field_error)
        !! Energy at z; with gradient and hessian also the derivatives. The
        !! Hessian is symmetric and only its UPPER triangle (columns >= row,
        !! diagonal included) is assembled, by rows (CSR): symmetric by
        !! construction, a third of the work of the full matrix.
        !! flux_map3d_expand_symmetric mirrors it where a full matrix is
        !! needed. derham and quad belong to the section
        !! geometry the gradients of clebsch_node_gradients refer to (plane 1
        !! of the callers here); plane_nodes(2, n0, nt) are the (R, Z) node
        !! coordinates of every plane. field_error: sum of w |B - B_exact|^2
        !! and of w |B_exact|^2 over the quadrature.
        type(flux_map3d_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        procedure(clebsch_exact_field), optional :: exact_field
        real(dp), intent(out), optional :: field_error(2)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), gx(:, :), sx(:, :, :), am(:, :, :, :), &
                                 tb(:, :, :), gloc(:), hloc(:, :), vals(:), zx(:, :), zy(:, :), &
                                 ax(:), ay(:), bx(:), by(:), lam(:), mw(:, :), gsec(:, :), &
                                 hs(:, :), lm(:), ltab(:, :, :)
        integer, allocatable :: d0(:), ix(:), iy(:), idx(:), rows(:), cols(:)
        logical, allocatable :: is_w(:)
        real(dp) :: lt(fm%nt, 2), ht(fm%nt, 2), sk(2, fm%nt), gk(2, 2, fm%nt), sm(2), f(2, 2), &
                    ft(2, 2), detf, vel(2), t, wq, wt, u(8), sig, th, gsig(3), gth(3), s, &
                    tau(3), gp(3), pr(3), w, wg(8), wh(36), h8(8, 8), ref, shift, k11, k12, &
                    vm(3), et(3), c3(3), gs(3), bv(3), bx3(3), ferr(2), resid, sn, cs
        integer :: n0, nt, ntot, nqt, cell, q, p, m, i, j, k, a, b, c1, c2, nloc, nsec, nl, &
                   nz, blk, il, jl, la, lb, row, gi, gj, nlmax, ia, ib
        logical :: derivs
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info

        n0 = fm%n0
        nt = fm%nt
        ntot = n0*nt
        nqt = fm%toroidal_points*nt + 2
        derivs = present(gradient)
        energy = 0
        ferr = 0
        nlmax = 2*(derham%degree + 1)**2*nt
        ! triplets of the upper triangle of one cell's local Hessian
        blk = nlmax*(nlmax + 1)/2
        ! the toroidal basis at the nqt points t_m, the same for every section
        ! point (one evaluation costs O(nt^2) trigonometric functions)
        allocate (ltab(nt, 2, nqt))
        do m = 1, nqt
            call spline_derham_t_basis(nt, (m - 1)/real(nqt, dp), lt, ht)
            ltab(:, 1, m) = lt(:, 1)
            ltab(:, 2, m) = clebsch_nfp*lt(:, 2)
        end do
        if (derivs) then
            allocate (gradient(2*ntot), source=0.0_dp)
            if (present(hessian)) then
                nz = size(quad%point_offset)*blk
                allocate (rows(nz + size(fm%gauge_col)**2 + size(fm%fixed)), &
                          cols(nz + size(fm%gauge_col)**2 + size(fm%fixed)), &
                          vals(nz + size(fm%gauge_col)**2 + size(fm%fixed)))
                rows = 1
                cols = 1
                vals = 0
            end if
        end if
        !$omp parallel do schedule(dynamic) default(private) &
        !$omp shared(fm, derham, quad, plane_nodes, z, gradient, rows, cols, vals, n0, nt, &
        !$omp        ntot, nqt, blk, derivs, clebsch_nfp, hessian, field_error, ltab) &
        !$omp reduction(+:energy, ferr)
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            if (allocated(ix)) deallocate (ix, iy, is_w, zx, zy, ax, ay, bx, by, lam, mw, gx, lm)
            allocate (ix(nloc), iy(nloc), is_w(nloc), zx(nloc, nt), zy(nloc, nt), ax(nloc), &
                      ay(nloc), bx(nloc), by(nloc), lam(nloc), mw(3, nloc), gx(2, nloc), lm(nloc))
            is_w = fm%is_wall(d0)
            nsec = 0
            do i = 1, nloc
                nsec = nsec + 1
                ix(i) = nsec
                iy(i) = 0
                if (.not. is_w(i)) then
                    nsec = nsec + 1
                    iy(i) = nsec
                end if
            end do
            nl = nsec*nt
            if (allocated(sx)) deallocate (sx, gsec, hs, idx, gloc)
            allocate (sx(8, nsec, 2), gsec(nsec, 2), hs(nsec, nsec), idx(nl), gloc(nl))
            if (present(hessian)) then
                if (allocated(am)) deallocate (am, tb, hloc)
                allocate (am(nsec, nsec, 3, nqt), tb(nt, 2, nqt), hloc(nl, nl))
                am = 0
                hloc = 0
            end if
            gloc = 0
            ! values per plane, the wall angles of the cell on one branch
            ref = 0
            do i = 1, nloc
                if (is_w(i)) then
                    ref = z(d0(i))
                    exit
                end if
            end do
            do i = 1, nloc
                shift = 0
                if (is_w(i)) shift = two_pi*nint((ref - z(d0(i)))/two_pi)
                do k = 1, nt
                    zx(i, k) = z((k - 1)*n0 + d0(i)) + shift
                    zy(i, k) = z(ntot + (k - 1)*n0 + d0(i))
                    idx((k - 1)*nsec + ix(i)) = (k - 1)*n0 + d0(i)
                    if (iy(i) > 0) idx((k - 1)*nsec + iy(i)) = ntot + (k - 1)*n0 + d0(i)
                end do
            end do
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                wq = quad%weight(p)
                do k = 1, nt
                    sk(:, k) = matmul(plane_nodes(:, d0, k), v0(1, :, q))
                    gk(:, :, k) = matmul(plane_nodes(:, d0, k), transpose(gn(:, :, q)))
                end do
                do m = 1, nqt
                    t = (m - 1)/real(nqt, dp)
                    lt = ltab(:, :, m)
                    sm = matmul(sk, lt(:, 1))
                    vel = matmul(sk, lt(:, 2))
                    f = 0
                    do k = 1, nt
                        f = f + lt(k, 1)*gk(:, :, k)
                    end do
                    detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
                    ft = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/detf
                    gx = matmul(ft, gn(:, :, q))
                    ax = matmul(zx, lt(:, 1))
                    ay = matmul(zy, lt(:, 1))
                    bx = matmul(zx, lt(:, 2))
                    by = matmul(zy, lt(:, 2))
                    call locals(is_w, ix, iy, v0(1, :, q), gx, vel, sm(1), ax, ay, bx, by, u, &
                                sx, sig, th, gsig, gth, lam, mw)
                    s = u(1)**2 + u(2)**2
                    call profiles(s, tau, gp, pr)
                    call flux_map_energy_3d(sm(1), u(1), u(2), u(3), u(4), u(5), u(6), u(7), &
                        u(8), tau(1), tau(2), tau(3), gp(1), gp(2), gp(3), pr(1), pr(2), pr(3), &
                        w, wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8), wh(1), &
                        wh(2), wh(3), wh(4), wh(5), wh(6), wh(7), wh(8), wh(9), wh(10), &
                        wh(11), wh(12), wh(13), wh(14), wh(15), wh(16), wh(17), wh(18), &
                        wh(19), wh(20), wh(21), wh(22), wh(23), wh(24), wh(25), wh(26), &
                        wh(27), wh(28), wh(29), wh(30), wh(31), wh(32), wh(33), wh(34), &
                        wh(35), wh(36))
                    wt = wq*abs(detf)*two_pi/nqt
                    energy = energy + wt*w
                    if (present(exact_field)) then
                        c3 = [u(4)*u(8) - u(5)*u(7), u(5)*u(6) - u(3)*u(8), u(3)*u(7) - u(4)*u(6)]
                        gs = 2*(u(1)*u(3:5) + u(2)*u(6:8))
                        bv = -[tau(1)*c3(1) + gp(1)*gs(3)/sm(1), tau(1)*c3(2), &
                               tau(1)*c3(3) - gp(1)*gs(1)/sm(1)]
                        call exact_field([sm(1), two_pi*t/clebsch_nfp, sm(2)], bx3)
                        ferr = ferr + wt*sm(1)*[sum((bv - bx3)**2), sum(bx3**2)]
                    end if
                    if (.not. derivs) cycle
                    do c1 = 1, 2
                        gsec(:, c1) = matmul(wg, sx(:, :, c1))
                    end do
                    do k = 1, nt
                        gloc((k - 1)*nsec + 1:k*nsec) = gloc((k - 1)*nsec + 1:k*nsec) + &
                            wt*(gsec(:, 1)*lt(k, 1) + gsec(:, 2)*lt(k, 2))
                    end do
                    if (.not. present(hessian)) cycle
                    a = 0
                    do i = 1, 8
                        do j = i, 8
                            a = a + 1
                            h8(i, j) = wh(a)
                            h8(j, i) = wh(a)
                        end do
                    end do
                    tb(:, :, m) = lt
                    ! blocks (1,1), (1,2), (2,2); (2,1) is the transpose of (1,2)
                    am(:, :, 1, m) = am(:, :, 1, m) + &
                        wt*matmul(transpose(sx(:, :, 1)), matmul(h8, sx(:, :, 1)))
                    am(:, :, 2, m) = am(:, :, 2, m) + &
                        wt*matmul(transpose(sx(:, :, 1)), matmul(h8, sx(:, :, 2)))
                    am(:, :, 3, m) = am(:, :, 3, m) + &
                        wt*matmul(transpose(sx(:, :, 2)), matmul(h8, sx(:, :, 2)))
                    if (sig > 0) then
                        ! second derivatives of the wall part, contracted with wg
                        sn = sin(th)
                        cs = cos(th)
                        et = [0.0_dp, 1/(two_pi*sm(1)), 0.0_dp]
                        k11 = -sig*(wg(1)*cs + wg(2)*sn) &
                              + dot_product(wg(3:5), -gsig*cs + sig*sn*gth) &
                              + dot_product(wg(6:8), -gsig*sn - sig*cs*gth)
                        vm = -sig*(cs*wg(3:5) + sn*wg(6:8))
                        k12 = dot_product(vm, et)
                        lm = matmul(vm, mw)
                        do j = 1, nloc
                            if (.not. is_w(j)) cycle
                            b = ix(j)
                            do i = 1, nloc
                                if (.not. is_w(i)) cycle
                                a = ix(i)
                                am(a, b, 1, m) = am(a, b, 1, m) + wt*(k11*lam(i)*lam(j) &
                                                 + lam(i)*lm(j) + lam(j)*lm(i))
                                am(a, b, 2, m) = am(a, b, 2, m) + wt*k12*lam(i)*lam(j)
                            end do
                        end do
                    end if
                end do
            end do
            if (.not. derivs) cycle
            do i = 1, nl
                !$omp atomic
                gradient(idx(i)) = gradient(idx(i)) + gloc(i)
            end do
            if (.not. present(hessian)) cycle
            ! upper block triangle of the local Hessian: block (la, lb), la <= lb,
            ! is sum_m [t1a t1b A11 + t1a t2b A12 + t2a t1b A12^T + t2a t2b A22]
            do m = 1, nqt
                do lb = 1, nt
                    do la = 1, lb
                        il = (la - 1)*nsec
                        jl = (lb - 1)*nsec
                        w = tb(la, 1, m)*tb(lb, 1, m)
                        if (w /= 0) hloc(il + 1:il + nsec, jl + 1:jl + nsec) = &
                            hloc(il + 1:il + nsec, jl + 1:jl + nsec) + w*am(:, :, 1, m)
                        w = tb(la, 1, m)*tb(lb, 2, m)
                        if (w /= 0) hloc(il + 1:il + nsec, jl + 1:jl + nsec) = &
                            hloc(il + 1:il + nsec, jl + 1:jl + nsec) + w*am(:, :, 2, m)
                        w = tb(la, 2, m)*tb(lb, 1, m)
                        if (w /= 0) hloc(il + 1:il + nsec, jl + 1:jl + nsec) = &
                            hloc(il + 1:il + nsec, jl + 1:jl + nsec) + &
                            w*transpose(am(:, :, 2, m))
                        w = tb(la, 2, m)*tb(lb, 2, m)
                        if (w /= 0) hloc(il + 1:il + nsec, jl + 1:jl + nsec) = &
                            hloc(il + 1:il + nsec, jl + 1:jl + nsec) + w*am(:, :, 3, m)
                    end do
                end do
            end do
            ! local entries (ia, ib), ia <= ib (a pair once; the diagonal blocks
            ! hold both orders, only ia <= ib is taken), to global (min, max)
            nz = (cell - 1)*blk
            do ib = 1, nl
                do ia = 1, ib
                    nz = nz + 1
                    rows(nz) = min(idx(ia), idx(ib))
                    cols(nz) = max(idx(ia), idx(ib))
                    vals(nz) = hloc(ia, ib)
                end do
            end do
        end do
        !$omp end parallel do
        if (present(field_error)) field_error = ferr
        nz = size(quad%point_offset)*blk
        ! gauge penalty
        do row = 1, size(fm%gauge_target)
            associate (cg => fm%gauge_col(fm%gauge_ptr(row):fm%gauge_ptr(row + 1) - 1), &
                       vg => fm%gauge_val(fm%gauge_ptr(row):fm%gauge_ptr(row + 1) - 1))
                resid = dot_product(vg, z(cg)) - fm%gauge_target(row)
                energy = energy + 0.5_dp*fm%gauge_weight*resid**2
                if (derivs) gradient(cg) = gradient(cg) + fm%gauge_weight*resid*vg
                if (present(hessian)) then
                    do b = 1, size(cg)
                        do a = 1, b
                            nz = nz + 1
                            rows(nz) = min(cg(a), cg(b))
                            cols(nz) = max(cg(a), cg(b))
                            vals(nz) = fm%gauge_weight*vg(a)*vg(b)
                        end do
                    end do
                end if
            end associate
        end do
        if (.not. present(hessian)) return
        shift = maxval(abs(vals(:nz)))
        do k = 1, size(fm%fixed)
            nz = nz + 1
            rows(nz) = fm%fixed(k)
            cols(nz) = fm%fixed(k)
            vals(nz) = shift
        end do
        call csc_from_triplet(2*ntot, 2*ntot, cols(:nz), rows(:nz), vals(:nz), mm, info)
        if (info%code /= 0) error stop 'flux map 3D Hessian triplets'
        hessian%nrow = 2*ntot
        hessian%ncol = 2*ntot
        hessian%row_ptr = mm%col_ptr
        hessian%col = mm%row_idx
        hessian%val = mm%val
    end subroutine flux_map3d_system

    subroutine flux_map3d_point(fm, plane_nodes, d0, nv, gn, t, profiles, rz, u, b)
        !! State at one point of the section and the toroidal parameter t:
        !! d0, nv, gn are the section functions that do not vanish there,
        !! their values and their gradients in the geometry of the reference
        !! plane (as from spline_point_matrix). Returns the position (R, Z),
        !! u = (x, y, grad x, grad y) and B in components (R, phi, Z).
        type(flux_map3d_t), intent(in) :: fm
        real(dp), intent(in) :: plane_nodes(:, :, :), nv(:), gn(:, :), t
        integer, intent(in) :: d0(:)
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: rz(2), u(8), b(3)
        real(dp) :: lt(fm%nt, 2), ht(fm%nt, 2), f(2, 2), ft(2, 2), detf, vel(2), sig, th, &
                    gsig(3), gth(3), tau(3), gp(3), pr(3), c3(3), gs(3), ref, shift
        real(dp) :: gx(2, size(d0)), ax(size(d0)), ay(size(d0)), bx(size(d0)), by(size(d0)), &
                    lam(size(d0)), mw(3, size(d0)), sx(8, 2*size(d0), 2), zx(size(d0), fm%nt), &
                    zy(size(d0), fm%nt)
        integer :: ix(size(d0)), iy(size(d0)), i, k, nsec, ntot
        logical :: is_w(size(d0))

        ntot = fm%n0*fm%nt
        call spline_derham_t_basis(fm%nt, t, lt, ht)
        lt(:, 2) = clebsch_nfp*lt(:, 2)
        is_w = fm%is_wall(d0)
        nsec = 0
        ref = 0
        do i = size(d0), 1, -1
            if (is_w(i)) ref = fm%z(d0(i))
        end do
        do i = 1, size(d0)
            nsec = nsec + 1
            ix(i) = nsec
            iy(i) = 0
            if (.not. is_w(i)) then
                nsec = nsec + 1
                iy(i) = nsec
            end if
            shift = 0
            if (is_w(i)) shift = two_pi*nint((ref - fm%z(d0(i)))/two_pi)
            do k = 1, fm%nt
                zx(i, k) = fm%z((k - 1)*fm%n0 + d0(i)) + shift
                zy(i, k) = fm%z(ntot + (k - 1)*fm%n0 + d0(i))
            end do
        end do
        rz = 0
        vel = 0
        f = 0
        do k = 1, fm%nt
            rz = rz + lt(k, 1)*matmul(plane_nodes(:, d0, k), nv)
            vel = vel + lt(k, 2)*matmul(plane_nodes(:, d0, k), nv)
            f = f + lt(k, 1)*matmul(plane_nodes(:, d0, k), transpose(gn))
        end do
        detf = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
        ft = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/detf
        gx = matmul(ft, gn)
        ax = matmul(zx, lt(:, 1))
        ay = matmul(zy, lt(:, 1))
        bx = matmul(zx, lt(:, 2))
        by = matmul(zy, lt(:, 2))
        call locals(is_w, ix, iy, nv, gx, vel, rz(1), ax, ay, bx, by, u, sx(:, :nsec, :), sig, &
                    th, gsig, gth, lam, mw)
        call profiles(u(1)**2 + u(2)**2, tau, gp, pr)
        c3 = [u(4)*u(8) - u(5)*u(7), u(5)*u(6) - u(3)*u(8), u(3)*u(7) - u(4)*u(6)]
        gs = 2*(u(1)*u(3:5) + u(2)*u(6:8))
        b = -[tau(1)*c3(1) + gp(1)*gs(3)/rz(1), tau(1)*c3(2), tau(1)*c3(3) - gp(1)*gs(1)/rz(1)]
    end subroutine flux_map3d_point

    subroutine flux_map3d_solve(fm, derham, quad, profiles, plane_nodes, tolerance, max_steps, &
                                status)
        !! Minimises W from fm%z by Newton steps with an identity shift, as
        !! flux_map_solve: a step is accepted when the energy decreases as the
        !! quadratic model predicts; a rejected step raises the shift.
        !! Converged when the gradient scaled by the Hessian diagonal is below
        !! tolerance. Status 1: no acceptable step, 2: step limit.
        type(flux_map3d_t), intent(inout) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: status
        type(csr_t) :: hessian, shifted, full
        type(eu_constrained_t) :: factor
        real(dp), allocatable :: grad(:), diag(:), delta(:), hd(:), zt(:), zero(:), gt(:)
        real(dp) :: w, wt, lam, dmin, pred, alpha
        integer :: it, i, p, n, st, nneg, nzero, npos
        logical :: accepted, sure
        logical, allocatable :: free(:)
        integer(int64) :: c0, c1, rate
        real(dp) :: t_hess, t_energy, t_solve

        n = 2*fm%n0*fm%nt
        t_hess = 0
        t_energy = 0
        t_solve = 0
        call system_clock(c0, rate)
        allocate (diag(n), delta(n), hd(n), free(n), zero(size(fm%fixed)))
        zero = 0
        free = .true.
        free(fm%fixed) = .false.
        lam = 1e-9_dp
        fm%steps = 0
        fm%converged = .false.
        fm%negative_eigenvalues = -1
        fm%zero_eigenvalues = -1
        fm%inertia_reliable = .false.
        status = 2
        if (fm%gauge_weight == 0 .and. size(fm%gauge_target) > 0) then
            call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, fm%z, w, grad, hessian)
            call diagonal(hessian, diag)
            fm%gauge_weight = sum(abs(diag), free)/count(free)
        end if
        do it = 0, max_steps
            call tick(c0)
            call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, fm%z, w, grad, hessian)
            call tock(c0, t_hess)
            call diagonal(hessian, diag)
            grad(fm%fixed) = 0
            fm%energy = w
            fm%gradient_norm = norm2(pack(grad/sqrt(abs(diag) + tiny(1.0_dp)), free))
            if (fm%verbose) print '(a,i4,a,es20.12,a,es10.2,a,es9.1)', '    step', it, '  W', w, &
                '  |g|', fm%gradient_norm, '  shift', lam
            if (fm%gradient_norm < tolerance) then
                fm%converged = .true.
                status = 0
                call inertia_at_convergence()
                call report()
                return
            end if
            if (it == max_steps) then
                call report()
                return
            end if
            dmin = minval(abs(diag), free)
            accepted = .false.
            if (fm%stationary) then
                call tick(c0)
                call constrained_solve(hessian, -grad, delta, st)
                call tock(c0, t_solve)
                if (st /= 0) then
                    status = 1
                    call report()
                    return
                end if
                alpha = 1
                do while (alpha > 1e-4_dp)
                    zt = fm%z + alpha*delta
                    call tick(c0)
                    call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, zt, wt, gt)
                    call tock(c0, t_energy)
                    gt(fm%fixed) = 0
                    if (norm2(pack(gt/sqrt(abs(diag) + tiny(1.0_dp)), free)) < &
                        (1 - alpha/4)*fm%gradient_norm) then
                        accepted = .true.
                        exit
                    end if
                    alpha = alpha/2
                end do
                if (.not. accepted) then
                    status = 1
                    call report()
                    return
                end if
                fm%z = zt
                fm%steps = fm%steps + 1
                cycle
            end if
            do while (lam < 1e12_dp)
                shifted = hessian
                do i = 1, n
                    do p = shifted%row_ptr(i), shifted%row_ptr(i + 1) - 1
                        if (shifted%col(p) == i) shifted%val(p) = shifted%val(p) + lam*dmin
                    end do
                end do
                call tick(c0)
                call constrained_solve(shifted, -grad, delta, st)
                call tock(c0, t_solve)
                if (st == 0) then
                    call upper_matvec(hessian, delta, hd)
                    pred = -(dot_product(grad, delta) + 0.5_dp*dot_product(delta, hd))
                    zt = fm%z + delta
                    call tick(c0)
                    call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, zt, wt)
                    call tock(c0, t_energy)
                    if (pred > 0 .and. wt == wt .and. &
                        w - wt > 1e-4_dp*pred - 1e-14_dp*abs(w)) then
                        fm%z = zt
                        fm%steps = fm%steps + 1
                        lam = max(lam/100, 1e-12_dp)
                        accepted = .true.
                        exit
                    end if
                end if
                lam = 10*lam
            end do
            if (.not. accepted) then
                status = 1
                call report()
                return
            end if
        end do

    contains

        subroutine tick(c)
            integer(int64), intent(out) :: c
            call system_clock(c)
        end subroutine tick

        subroutine tock(c, acc)
            integer(int64), intent(in) :: c
            real(dp), intent(inout) :: acc
            integer(int64) :: now
            call system_clock(now)
            acc = acc + real(now - c, dp)/rate
        end subroutine tock

        subroutine constrained_solve(k, rhs, x, status)
            !! K x = rhs on the free slots, x = 0 on the fixed ones; one factor
            !! object for the whole solve, so the ordering and the symbolic
            !! analysis of the fixed pattern are computed once.
            type(csr_t), intent(in) :: k
            real(dp), intent(in) :: rhs(:)
            real(dp), intent(out) :: x(:)
            integer, intent(out) :: status

            call factorise(k, status)
            if (status == 0) call eu_constrained_apply(factor, rhs, x, status)
        end subroutine constrained_solve

        subroutine factorise(k, status)
            !! k: upper triangle of a symmetric matrix of the pattern of the
            !! Hessian: LDL^T of it, or LU of its mirror.
            type(csr_t), intent(in) :: k
            integer, intent(out) :: status

            if (fm%symmetric_solver) then
                call eu_constrained_refactor(k, fm%fixed, factor, status, symmetric=.true.)
            else
                call flux_map3d_expand_symmetric(k, full)
                call eu_constrained_refactor(full, fm%fixed, factor, status)
            end if
        end subroutine factorise

        subroutine inertia_at_convergence()
            !! Inertia of the Hessian itself (no shift) at the final state.
            if (.not. fm%symmetric_solver) return
            call tick(c0)
            call factorise(hessian, st)
            call tock(c0, t_solve)
            ! an exactly singular block still has its counts
            call eu_constrained_inertia(factor, nneg, nzero, npos, sure, st)
            if (st /= 0) return
            fm%negative_eigenvalues = nneg
            fm%zero_eigenvalues = nzero
            fm%inertia_reliable = sure
            if (fm%verbose) print '(a,i0,a,i0,a,i0,a,l1)', &
                '    Hessian inertia at the solution (negative, zero, positive): ', nneg, ', ', &
                nzero, ', ', npos, '  reliable ', sure
        end subroutine inertia_at_convergence

        subroutine report()
            call eu_constrained_free(factor)
            if (fm%verbose) print '(a,3f9.2)', '    time: energy+gradient+Hessian, energy(+gradient), linear solves (s)', &
                t_hess, t_energy, t_solve
        end subroutine report
    end subroutine flux_map3d_solve

    subroutine flux_map3d_expand_symmetric(upper, full)
        !! The full symmetric matrix of a Hessian given as its upper
        !! triangle by rows: every off-diagonal entry mirrored, exactly
        !! symmetric, rows ascending.
        type(csr_t), intent(in) :: upper
        type(csr_t), intent(out) :: full
        integer, allocatable :: cnt(:), fill(:)
        integer :: n, i, p, j, q

        n = upper%nrow
        allocate (cnt(n), fill(n))
        cnt = 0
        do i = 1, n
            do p = upper%row_ptr(i), upper%row_ptr(i + 1) - 1
                j = upper%col(p)
                cnt(i) = cnt(i) + 1
                if (j /= i) cnt(j) = cnt(j) + 1
            end do
        end do
        full%nrow = n
        full%ncol = n
        allocate (full%row_ptr(n + 1))
        full%row_ptr(1) = 1
        do i = 1, n
            full%row_ptr(i + 1) = full%row_ptr(i) + cnt(i)
        end do
        allocate (full%col(full%row_ptr(n + 1) - 1), full%val(full%row_ptr(n + 1) - 1))
        fill = full%row_ptr(1:n)
        ! rows ascending: row j receives the mirrored entries of the rows i < j
        ! (in the order of i) before its own, which have columns >= j
        do i = 1, n
            do p = upper%row_ptr(i), upper%row_ptr(i + 1) - 1
                j = upper%col(p)
                q = fill(i)
                full%col(q) = j
                full%val(q) = upper%val(p)
                fill(i) = q + 1
                if (j /= i) then
                    q = fill(j)
                    full%col(q) = i
                    full%val(q) = upper%val(p)
                    fill(j) = q + 1
                end if
            end do
        end do
    end subroutine flux_map3d_expand_symmetric

    pure subroutine upper_matvec(a, x, y)
        !! y = A x for the symmetric A given by its upper triangle by rows.
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, p, j

        y = 0
        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col(p)
                y(i) = y(i) + a%val(p)*x(j)
                if (j /= i) y(j) = y(j) + a%val(p)*x(i)
            end do
        end do
    end subroutine upper_matvec

    pure subroutine diagonal(a, d)
        type(csr_t), intent(in) :: a
        real(dp), intent(out) :: d(:)
        integer :: i, p

        d = 0
        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(p) == i) d(i) = d(i) + a%val(p)
            end do
        end do
    end subroutine diagonal

end module milca_flux_map_3d
