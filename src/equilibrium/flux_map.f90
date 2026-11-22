module milca_flux_map
    !! Axisymmetric Cartesian flux map (docs/flux-map.md): the unknown is the
    !! map u = (x, y) from the section to the unit disc, s = x^2 + y^2 the
    !! flux label, with
    !!
    !!   B_phi = tau(s) det(grad u),  B_pol = psi'(s) grad s x grad phi,
    !!   W = 2 pi int [B_phi^2/2 + |B_pol|^2/2 - p(s)] R dR dZ.
    !!
    !! x and y are in V0 of a spline section whose mesh knows nothing about
    !! the magnetic axis: the axis is the zero of u. On the wall u = (cos
    !! Theta, sin Theta) holds along the whole wall, not only at nodes: with
    !! sigma the sum of the wall functions, the wall part of u is sigma (cos,
    !! sin)(Theta_h), Theta_h = sum_b N_b Theta_b/sigma. The toroidal flux is
    !! then the area of the disc weighted by tau and cannot leak.
    !!
    !! Unknowns z(1:2 n0): z(i) = x_i (interior function) or Theta_i (wall
    !! function), z(n0 + i) = y_i (interior; unused and fixed for wall
    !! functions). B_phi enters the energy pointwise (projection < 0) or as
    !! its L2 projection in the area measure onto the polynomials of degree
    !! projection in each parameter of a cell, which removes the locking at a
    !! strong guide field. Kernels: generated leaves flux_map_toroidal and
    !! flux_map_poloidal.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t
    use fortsparse_csr, only: csr_t, csr_matvec
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_flux_map_leaves, only: flux_map_toroidal, flux_map_poloidal
    implicit none
    private

    abstract interface
        subroutine flux_map_profiles(s, tau, gp, p)
            !! tau = (tau, tau', tau''), gp = (psi', psi'', psi'''),
            !! p = (p, p', p'') at the label s; tau = Phi'(s)/pi with Phi the
            !! toroidal flux inside s.
            import :: dp
            real(dp), intent(in) :: s
            real(dp), intent(out) :: tau(3), gp(3), p(3)
        end subroutine flux_map_profiles
    end interface

    type, public :: flux_map_t
        integer :: n0 = 0
        logical, allocatable :: is_wall(:)
        integer, allocatable :: fixed(:) ! y slots of the wall functions
        real(dp), allocatable :: z(:)
        !! degree of the per-cell projection of B_phi; negative: pointwise
        integer :: projection = -1
        !! Gauge theta* -> theta* + k(s) by a penalty on linear rows
        !! sum_j gauge_val z(gauge_col) = gauge_target. Exact for up-down
        !! symmetric states with rows y_h = 0 on the symmetry line; the
        !! general case needs constraints orthogonal to the gauge instead.
        integer, allocatable :: gauge_ptr(:), gauge_col(:)
        real(dp), allocatable :: gauge_val(:), gauge_target(:)
        real(dp) :: gauge_weight = 0
        integer :: steps = 0
        real(dp) :: energy = 0, gradient_norm = 0
        logical :: converged = .false., verbose = .false.
    end type flux_map_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    public :: flux_map_profiles, flux_map_init, flux_map_gauge_row, flux_map_system, &
              flux_map_solve

contains

    subroutine flux_map_init(fm, n0, wall)
        !! wall: the V0 functions with a trace on the wall.
        type(flux_map_t), intent(out) :: fm
        integer, intent(in) :: n0, wall(:)

        fm%n0 = n0
        allocate (fm%is_wall(n0), source=.false.)
        fm%is_wall(wall) = .true.
        fm%fixed = n0 + wall
        allocate (fm%z(2*n0), source=0.0_dp)
        allocate (fm%gauge_ptr(1), source=1)
        allocate (fm%gauge_col(0), fm%gauge_val(0), fm%gauge_target(0))
    end subroutine flux_map_init

    subroutine flux_map_gauge_row(fm, cols, vals, target)
        !! Adds the gauge row sum vals z(cols) = target.
        type(flux_map_t), intent(inout) :: fm
        integer, intent(in) :: cols(:)
        real(dp), intent(in) :: vals(:), target

        fm%gauge_col = [fm%gauge_col, cols]
        fm%gauge_val = [fm%gauge_val, vals]
        fm%gauge_target = [fm%gauge_target, target]
        fm%gauge_ptr = [fm%gauge_ptr, size(fm%gauge_col) + 1]
    end subroutine flux_map_gauge_row

    subroutine flux_map_system(fm, derham, quad, profiles, z, energy, gradient, hessian, &
                               label, grad_label, bphi, bphi_projected)
        !! Energy at z; with gradient and hessian also the derivatives (the
        !! Hessian by rows, symmetric). The optional arrays return, per
        !! quadrature point, s, grad s in (R, Z), and B_phi pointwise and
        !! projected.
        type(flux_map_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), intent(out), optional :: label(:), grad_label(:, :), bphi(:), &
                                           bphi_projected(:)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), d2(:, :, :), &
                                 hloc(:, :), gloc(:), cq(:), dc(:, :), d2c(:, :, :), &
                                 proj(:, :), cb(:), lam(:), wr(:), pdc(:, :), vals(:), &
                                 lw(:), mw(:, :)
        integer, allocatable :: d0(:), idx(:), ix(:), iy(:), rows(:), cols(:)
        real(dp) :: jet(2, 6), wq, u(6), sig, gsig(2), th, gth(2), acc, gacc(2), ref, tb, &
                    tau(3), gp(3), pr(3), c, cg(6), ch(21), w, wg(6), wh(21), h6(6, 6), s, &
                    resid, shift
        integer :: n0, cell, q, p, i, j, a, b, k, nloc, nu, nq, nz, st, row
        logical :: derivs
        type(csc_t) :: m
        type(fortsparse_status_t) :: info

        n0 = fm%n0
        derivs = present(gradient)
        nq = quad%points_per_cell
        energy = 0
        if (derivs) then
            allocate (gradient(2*n0), source=0.0_dp)
            nz = size(quad%point_offset)*(2*(derham%degree + 1)**2)**2 + size(fm%gauge_col)**2 &
                 + size(fm%fixed)
            allocate (rows(nz), cols(nz), vals(nz))
            nz = 0
        end if
        allocate (cq(nq), cb(nq), lam(nq), wr(nq))
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            if (allocated(idx)) deallocate (idx, ix, iy, du, d2, hloc, gloc, dc, d2c, pdc, lw, mw)
            allocate (ix(nloc), iy(nloc), idx(2*nloc), lw(nloc), mw(2, nloc))
            nu = 0
            do i = 1, nloc
                nu = nu + 1
                ix(i) = nu
                idx(nu) = d0(i)
                iy(i) = 0
                if (.not. fm%is_wall(d0(i))) then
                    nu = nu + 1
                    iy(i) = nu
                    idx(nu) = n0 + d0(i)
                end if
            end do
            allocate (du(6, nu), d2(6, nu, nu), hloc(nu, nu), gloc(nu), dc(nu, nq), &
                      d2c(nu, nu, nq), pdc(nu, nq))
            hloc = 0
            gloc = 0
            ! wall angles of this cell on one branch
            ref = 0
            do i = 1, nloc
                if (fm%is_wall(d0(i))) then
                    ref = z(d0(i))
                    exit
                end if
            end do
            do q = 1, nq
                p = quad%point_offset(cell) + q - 1
                call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
                wq = two_pi*quad%weight(p)
                wr(q) = wq*jet(1, 1)
                u = 0
                du = 0
                d2 = 0
                sig = 0
                gsig = 0
                acc = 0
                gacc = 0
                do i = 1, nloc
                    if (fm%is_wall(d0(i))) then
                        tb = z(d0(i)) + two_pi*nint((ref - z(d0(i)))/two_pi)
                        sig = sig + v0(1, i, q)
                        gsig = gsig + gn(:, i, q)
                        acc = acc + v0(1, i, q)*tb
                        gacc = gacc + gn(:, i, q)*tb
                    else
                        u(1) = u(1) + v0(1, i, q)*z(d0(i))
                        u(2) = u(2) + v0(1, i, q)*z(n0 + d0(i))
                        u(3:4) = u(3:4) + gn(:, i, q)*z(d0(i))
                        u(5:6) = u(5:6) + gn(:, i, q)*z(n0 + d0(i))
                        du(1, ix(i)) = v0(1, i, q)
                        du(3:4, ix(i)) = gn(:, i, q)
                        du(2, iy(i)) = v0(1, i, q)
                        du(5:6, iy(i)) = gn(:, i, q)
                    end if
                end do
                if (sig > 0) then
                    th = acc/sig
                    gth = (gacc - th*gsig)/sig
                    u(1) = u(1) + sig*cos(th)
                    u(2) = u(2) + sig*sin(th)
                    u(3:4) = u(3:4) + gsig*cos(th) - sig*sin(th)*gth
                    u(5:6) = u(5:6) + gsig*sin(th) + sig*cos(th)*gth
                    do i = 1, nloc
                        if (.not. fm%is_wall(d0(i))) cycle
                        lw(i) = v0(1, i, q)/sig
                        mw(:, i) = (gn(:, i, q) - lw(i)*gsig)/sig
                        a = ix(i)
                        du(1, a) = -sig*sin(th)*lw(i)
                        du(2, a) = sig*cos(th)*lw(i)
                        du(3:4, a) = -gsig*sin(th)*lw(i) - sig*cos(th)*lw(i)*gth - sig*sin(th)*mw(:, i)
                        du(5:6, a) = gsig*cos(th)*lw(i) - sig*sin(th)*lw(i)*gth + sig*cos(th)*mw(:, i)
                    end do
                    if (derivs) then
                        do i = 1, nloc
                            if (.not. fm%is_wall(d0(i))) cycle
                            do j = 1, nloc
                                if (.not. fm%is_wall(d0(j))) cycle
                                a = ix(i)
                                b = ix(j)
                                d2(1, a, b) = -sig*cos(th)*lw(i)*lw(j)
                                d2(2, a, b) = -sig*sin(th)*lw(i)*lw(j)
                                d2(3:4, a, b) = (-gsig*cos(th) + sig*sin(th)*gth)*lw(i)*lw(j) &
                                                - sig*cos(th)*(lw(i)*mw(:, j) + lw(j)*mw(:, i))
                                d2(5:6, a, b) = (-gsig*sin(th) - sig*cos(th)*gth)*lw(i)*lw(j) &
                                                - sig*sin(th)*(lw(i)*mw(:, j) + lw(j)*mw(:, i))
                            end do
                        end do
                    end if
                end if
                s = u(1)**2 + u(2)**2
                call profiles(s, tau, gp, pr)
                call flux_map_toroidal(u(1), u(2), u(3), u(4), u(5), u(6), tau(1), tau(2), &
                                       tau(3), c, cg(1), cg(2), cg(3), cg(4), cg(5), cg(6), &
                                       ch(1), ch(2), ch(3), ch(4), ch(5), ch(6), ch(7), ch(8), &
                                       ch(9), ch(10), ch(11), ch(12), ch(13), ch(14), ch(15), &
                                       ch(16), ch(17), ch(18), ch(19), ch(20), ch(21))
                call flux_map_poloidal(jet(1, 1), u(1), u(2), u(3), u(4), u(5), u(6), gp(1), &
                                       gp(2), gp(3), pr(1), pr(2), pr(3), w, wg(1), wg(2), &
                                       wg(3), wg(4), wg(5), wg(6), wh(1), wh(2), wh(3), wh(4), &
                                       wh(5), wh(6), wh(7), wh(8), wh(9), wh(10), wh(11), &
                                       wh(12), wh(13), wh(14), wh(15), wh(16), wh(17), wh(18), &
                                       wh(19), wh(20), wh(21))
                energy = energy + wq*w
                cq(q) = c
                if (present(label)) then
                    label(p) = s
                    grad_label(:, p) = 2*(u(1)*u(3:4) + u(2)*u(5:6))
                    bphi(p) = c
                end if
                if (.not. derivs) cycle
                gloc = gloc + wq*matmul(wg, du)
                call unpack6(wh, h6)
                hloc = hloc + wq*matmul(transpose(du), matmul(h6, du))
                dc(:, q) = matmul(cg, du)
                call unpack6(ch, h6)
                d2c(:, :, q) = matmul(transpose(du), matmul(h6, du))
                do a = 1, 6
                    hloc = hloc + wq*wg(a)*d2(a, :, :)
                    d2c(:, :, q) = d2c(:, :, q) + cg(a)*d2(a, :, :)
                end do
            end do
            ! toroidal energy sum_q wr cb^2/2, cb = proj c
            if (fm%projection >= 0) then
                call cell_projection(quad, cell, fm%projection, proj)
                cb = matmul(proj, cq)
            else
                cb = cq
            end if
            energy = energy + 0.5_dp*sum(wr*cb**2)
            if (present(bphi_projected)) bphi_projected(quad%point_offset(cell): &
                                                        quad%point_offset(cell) + nq - 1) = cb
            if (derivs) then
                if (fm%projection >= 0) then
                    lam = matmul(wr*cb, proj)
                    pdc = matmul(dc, transpose(proj))
                else
                    lam = wr*cb
                    pdc = dc
                end if
                gloc = gloc + matmul(dc, lam)
                do q = 1, nq
                    hloc = hloc + lam(q)*d2c(:, :, q)
                    do b = 1, nu
                        hloc(:, b) = hloc(:, b) + wr(q)*pdc(:, q)*pdc(b, q)
                    end do
                end do
                gradient(idx(:nu)) = gradient(idx(:nu)) + gloc
                do b = 1, nu
                    do a = 1, nu
                        nz = nz + 1
                        rows(nz) = idx(a)
                        cols(nz) = idx(b)
                        vals(nz) = hloc(a, b)
                    end do
                end do
            end if
        end do
        ! gauge penalty
        do row = 1, size(fm%gauge_target)
            associate (cg_ => fm%gauge_col(fm%gauge_ptr(row):fm%gauge_ptr(row + 1) - 1), &
                       vg => fm%gauge_val(fm%gauge_ptr(row):fm%gauge_ptr(row + 1) - 1))
                resid = dot_product(vg, z(cg_)) - fm%gauge_target(row)
                energy = energy + 0.5_dp*fm%gauge_weight*resid**2
                if (derivs) then
                    gradient(cg_) = gradient(cg_) + fm%gauge_weight*resid*vg
                    do b = 1, size(cg_)
                        do a = 1, size(cg_)
                            nz = nz + 1
                            rows(nz) = cg_(a)
                            cols(nz) = cg_(b)
                            vals(nz) = fm%gauge_weight*vg(a)*vg(b)
                        end do
                    end do
                end if
            end associate
        end do
        if (.not. derivs) return
        ! unit diagonal for the unused y slots of the wall functions
        shift = maxval(abs(vals(:nz)))
        do k = 1, size(fm%fixed)
            nz = nz + 1
            rows(nz) = fm%fixed(k)
            cols(nz) = fm%fixed(k)
            vals(nz) = shift
        end do
        if (.not. present(hessian)) return
        call csc_from_triplet(2*n0, 2*n0, cols(:nz), rows(:nz), vals(:nz), m, info)
        if (info%code /= 0) error stop 'flux map Hessian triplets'
        hessian%nrow = 2*n0
        hessian%ncol = 2*n0
        hessian%row_ptr = m%col_ptr
        hessian%col = m%row_idx
        hessian%val = m%val
    end subroutine flux_map_system

    pure subroutine unpack6(h, hm)
        real(dp), intent(in) :: h(21)
        real(dp), intent(out) :: hm(6, 6)
        integer :: a, b, k

        k = 0
        do a = 1, 6
            do b = a, 6
                k = k + 1
                hm(a, b) = h(k)
                hm(b, a) = h(k)
            end do
        end do
    end subroutine unpack6

    subroutine cell_projection(quad, cell, degree, proj)
        !! proj(q, q'): L2 projection in the area measure of point values in
        !! a cell onto the polynomials of the given degree in each parameter.
        !! It reproduces constants, so the flux sum_q weight c is unchanged.
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell, degree
        real(dp), allocatable, intent(out) :: proj(:, :)
        real(dp), allocatable :: pm(:, :), mass(:, :), rhs(:, :)
        real(dp) :: lo(2), hi(2), t(2), f
        integer :: nq, nt, q, p, i, j, k, a, b

        nq = quad%points_per_cell
        nt = (degree + 1)**2
        allocate (pm(nq, nt), mass(nt, nt), rhs(nt, nq), proj(nq, nq))
        p = quad%point_offset(cell)
        lo = minval(quad%xi(:, p:p + nq - 1), 2)
        hi = maxval(quad%xi(:, p:p + nq - 1), 2)
        do q = 1, nq
            t = (2*quad%xi(:, p + q - 1) - lo - hi)/(hi - lo)
            k = 0
            do i = 0, degree
                do j = 0, degree
                    k = k + 1
                    pm(q, k) = t(1)**i*t(2)**j
                end do
            end do
            rhs(:, q) = pm(q, :)*quad%weight(p + q - 1)
        end do
        mass = matmul(rhs, pm)
        ! Gauss elimination without pivoting (mass is symmetric positive definite)
        do a = 1, nt
            do b = a + 1, nt
                f = mass(b, a)/mass(a, a)
                mass(b, a:) = mass(b, a:) - f*mass(a, a:)
                rhs(b, :) = rhs(b, :) - f*rhs(a, :)
            end do
        end do
        do a = nt, 1, -1
            do b = a + 1, nt
                rhs(a, :) = rhs(a, :) - mass(a, b)*rhs(b, :)
            end do
            rhs(a, :) = rhs(a, :)/mass(a, a)
        end do
        proj = matmul(pm, rhs)
    end subroutine cell_projection

    subroutine flux_map_solve(fm, derham, quad, profiles, tolerance, max_steps, status)
        !! Minimises W from fm%z by Newton steps with an identity shift: a
        !! step is accepted when the energy decreases as the quadratic model
        !! predicts; a rejected step raises the shift. Converged when the
        !! gradient scaled by the Hessian diagonal is below tolerance. Status
        !! 1: no acceptable step, 2: step limit.
        type(flux_map_t), intent(inout) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: status
        type(csr_t) :: hessian, shifted
        real(dp), allocatable :: grad(:), diag(:), delta(:), hd(:), zt(:), zero(:)
        real(dp) :: w, wt, lam, dmin, pred
        integer :: it, i, p, n, st
        logical :: accepted
        logical, allocatable :: free(:)

        n = 2*fm%n0
        allocate (diag(n), delta(n), hd(n), free(n), zero(size(fm%fixed)))
        zero = 0
        free = .true.
        free(fm%fixed) = .false.
        lam = 1e-9_dp
        fm%steps = 0
        fm%converged = .false.
        status = 2
        if (fm%gauge_weight == 0 .and. size(fm%gauge_target) > 0) then
            ! weight of the gauge penalty at the scale of the energy Hessian
            call flux_map_system(fm, derham, quad, profiles, fm%z, w, grad, hessian)
            call diagonal(hessian, diag)
            fm%gauge_weight = sum(abs(diag), free)/count(free)
        end if
        do it = 0, max_steps
            call flux_map_system(fm, derham, quad, profiles, fm%z, w, grad, hessian)
            call diagonal(hessian, diag)
            grad(fm%fixed) = 0
            fm%energy = w
            fm%gradient_norm = norm2(pack(grad/sqrt(abs(diag) + tiny(1.0_dp)), free))
            if (fm%verbose) print '(a,i4,a,es20.12,a,es10.2,a,es9.1)', '    step', it, '  W', w, &
                '  |g|', fm%gradient_norm, '  shift', lam
            if (fm%gradient_norm < tolerance) then
                fm%converged = .true.
                status = 0
                return
            end if
            if (it == max_steps) return
            dmin = minval(abs(diag), free)
            accepted = .false.
            do while (lam < 1e12_dp)
                shifted = hessian
                do i = 1, n
                    do p = shifted%row_ptr(i), shifted%row_ptr(i + 1) - 1
                        if (shifted%col(p) == i) shifted%val(p) = shifted%val(p) + lam*dmin
                    end do
                end do
                call eu_constrained_solve(shifted, -grad, fm%fixed, zero, delta, st)
                if (st == 0) then
                    call csr_matvec(hessian, delta, hd)
                    pred = -(dot_product(grad, delta) + 0.5_dp*dot_product(delta, hd))
                    zt = fm%z + delta
                    call flux_map_system(fm, derham, quad, profiles, zt, wt)
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
                return
            end if
        end do
    end subroutine flux_map_solve

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

end module milca_flux_map
