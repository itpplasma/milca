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
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use fortfem_constrained_newton, only: newton_constraints_t, newton_options_t, &
        newton_report_t, newton_system, newton_constraint_update, &
        newton_constraints_init, constrained_newton, &
        NEWTON_ACCEPT_ENERGY, NEWTON_ACCEPT_ENERGY_RELAXED, NEWTON_ACCEPT_GRADIENT
    use milca_eulerian_clebsch, only: clebsch_node_gradients
    use milca_gauge_moments, only: gauge_moment_constraints
    use milca_flux_map_leaves, only: flux_map_toroidal, flux_map_poloidal, &
        flux_map_gauge_moment
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
        !! Gauge. Of theta* -> theta* + k(s) the discrete energy keeps one
        !! exact symmetry, the rigid rotation k = const (interior values
        !! rotate, wall angles shift). It is removed by one moment constraint
        !! orthogonal to it (milca_gauge_moments), which follows the state:
        !! every Newton step keeps the mean wall angle, weighted by the
        !! moment of the rotation generator (-y, x) with each wall function
        !! (flux_map_gauge). .false.: no constraint (the Hessian is singular
        !! at a stationary point).
        logical :: gauge_constraint = .true.
        type(newton_constraints_t) :: gauge
        !! Solver: .false. minimises W (steps accepted on energy decrease);
        !! .true. seeks the stationary point next to the start by Newton
        !! with a line search on the scaled gradient norm, whatever the
        !! signature of the Hessian (an equilibrium need not be a minimum).
        logical :: stationary = .false.
        !! Minimisation with relaxed steps (fortfem_constrained_newton,
        !! NEWTON_ACCEPT_ENERGY_RELAXED): fewer steps along a curved soft
        !! valley (measured: the cold start on a mesh without symmetry), more
        !! on the standard decks (docs/flux-map.md).
        logical :: relaxed_steps = .false.
        integer :: steps = 0
        real(dp) :: energy = 0, gradient_norm = 0
        logical :: converged = .false., verbose = .false.
        !! Inertia of the Hessian on the tangent space of the gauge
        !! constraints at the converged state: numbers of negative and zero
        !! eigenvalues (-1: not converged), and whether every pivot was clear
        !! of its rounding threshold.
        integer :: negative_eigenvalues = -1, zero_eigenvalues = -1
        logical :: inertia_reliable = .false.
    end type flux_map_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    public :: flux_map_profiles, flux_map_init, flux_map_gauge, flux_map_system, &
              flux_map_solve, flux_map_gauge_row, flux_map_newton

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
        call newton_constraints_init(fm%gauge)
    end subroutine flux_map_init

    subroutine flux_map_gauge_row(moment, wall_slots, z, gauge)
        !! The gauge constraint at the state z from the moment of the
        !! rotation generator with the wall functions: the row dz -> sum_b
        !! moment_b dTheta_b through z. Its Gram number with the generator
        !! (dTheta_b = 1) is the sum of the moments. An earlier gauge keeps
        !! its pattern.
        real(dp), intent(in) :: moment(:), z(:)
        integer, intent(in) :: wall_slots(:)
        type(newton_constraints_t), intent(inout) :: gauge
        real(dp) :: gram(1, 1)
        logical :: kept(1)

        gram(1, 1) = sum(moment(wall_slots))
        call gauge_moment_constraints(reshape(moment, [size(moment), 1]), gram, z, gauge, kept)
        if (.not. kept(1)) error stop 'flux map: the gauge moment vanishes'
    end subroutine flux_map_gauge_row

    subroutine flux_map_gauge(fm, derham, quad, profiles, z, gauge)
        !! The gauge constraint at the state z (the reference of
        !! milca_gauge_moments).
        type(flux_map_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: z(:)
        type(newton_constraints_t), intent(inout) :: gauge
        real(dp) :: moment(2*fm%n0), energy

        if (.not. fm%gauge_constraint) then
            call newton_constraints_init(gauge)
            return
        end if
        call flux_map_system(fm, derham, quad, profiles, z, energy, gauge_moment=moment)
        call flux_map_gauge_row(moment, pack(fm%fixed - fm%n0, .true.), z, gauge)
    end subroutine flux_map_gauge

    subroutine flux_map_system(fm, derham, quad, profiles, z, energy, gradient, hessian, &
                               label, grad_label, bphi, bphi_projected, gauge_moment)
        !! Energy at z; with gradient and hessian also the derivatives (the
        !! Hessian as its UPPER triangle by rows). The optional arrays
        !! return, per quadrature point, s, grad s in (R, Z), and B_phi
        !! pointwise and projected. gauge_moment: per wall function b the
        !! moment int (-y, x) . du/dTheta_b dV of the rotation generator at
        !! z, zero in the other slots (flux_map_gauge).
        type(flux_map_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), intent(out), optional :: label(:), grad_label(:, :), bphi(:), &
                                           bphi_projected(:), gauge_moment(:)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), d2(:, :, :), &
                                 hloc(:, :), gloc(:), cq(:), dc(:, :), d2c(:, :, :), &
                                 proj(:, :), cb(:), lam(:), wr(:), pdc(:, :), vals(:), &
                                 lw(:), mw(:, :), turn(:)
        integer, allocatable :: d0(:), idx(:), ix(:), iy(:), rows(:), cols(:)
        real(dp) :: jet(2, 6), wq, u(6), sig, gsig(2), th, gth(2), acc, gacc(2), ref, tb, &
                    tau(3), gp(3), pr(3), c, cg(6), ch(21), w, wg(6), wh(21), h6(6, 6), s, &
                    shift
        integer :: n0, cell, q, p, i, j, a, b, k, nloc, nu, nq, nz, st
        logical :: derivs
        type(csc_t) :: m
        type(fortsparse_status_t) :: info

        n0 = fm%n0
        derivs = present(gradient)
        nq = quad%points_per_cell
        energy = 0
        if (present(gauge_moment)) gauge_moment = 0
        if (derivs) then
            allocate (gradient(2*n0), source=0.0_dp)
            nz = size(quad%point_offset)*(2*(derham%degree + 1)**2)**2 + size(fm%fixed)
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
                if (present(gauge_moment)) then
                    turn = [(0.0_dp, i=1, nu)]
                    do i = 1, nloc
                        if (fm%is_wall(d0(i))) call flux_map_gauge_moment(u(1), u(2), &
                            du(1, ix(i)), du(2, ix(i)), turn(ix(i)))
                    end do
                    gauge_moment(idx(:nu)) = gauge_moment(idx(:nu)) + wr(q)*turn
                end if
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
                ! a pair once, to the upper triangle
                do b = 1, nu
                    do a = 1, b
                        nz = nz + 1
                        rows(nz) = min(idx(a), idx(b))
                        cols(nz) = max(idx(a), idx(b))
                        vals(nz) = hloc(a, b)
                    end do
                end do
            end if
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

    subroutine flux_map_newton(system, z, fixed, gauge, update, stationary, relaxed, &
                               verbose, tolerance, max_steps, report, status, model_agreement_damping)
        !! The solver loop of the flux map, axisymmetric and 3D: Newton on
        !! the KKT system of W with the gauge constraints, by LDL^T with its
        !! inertia (fortfem_constrained_newton). update rebuilds the gauge
        !! row at every iterate. Not stationary: W is minimised (a step is
        !! accepted when W decreases as the quadratic model predicts, with
        !! relaxed steps if relaxed; a rejected step raises an identity
        !! shift).
        !! stationary: the stationary point next to z by a line search on
        !! the scaled gradient norm. Converged when the gradient of the Lagrangian
        !! scaled by the Hessian diagonal is below tolerance. Status 1: no
        !! acceptable step, 2: step limit.
        procedure(newton_system) :: system
        real(dp), intent(inout) :: z(:)
        integer, intent(in) :: fixed(:), max_steps
        type(newton_constraints_t), intent(inout) :: gauge
        procedure(newton_constraint_update) :: update
        logical, intent(in) :: stationary, relaxed, verbose
        real(dp), intent(in) :: tolerance
        type(newton_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: model_agreement_damping
        type(newton_options_t) :: options

        options%tolerance = tolerance
        options%max_steps = max_steps
        options%acceptance = NEWTON_ACCEPT_ENERGY
        if (relaxed) options%acceptance = NEWTON_ACCEPT_ENERGY_RELAXED
        if (stationary) options%acceptance = NEWTON_ACCEPT_GRADIENT
        options%verbose = verbose
        if (present(model_agreement_damping)) options%model_agreement_damping = model_agreement_damping
        call constrained_newton(system, z, fixed, gauge, options, report, status, update)
        if (verbose) print '(a,3f9.2)', &
            '    time: energy+gradient+Hessian, energy(+gradient), linear solves (s)', &
            report%time_hessian, report%time_energy, report%time_solve
    end subroutine flux_map_newton

    subroutine flux_map_solve(fm, derham, quad, profiles, tolerance, max_steps, status)
        !! Solves from fm%z by flux_map_newton, with the gauge constraints
        !! built at fm%z.
        type(flux_map_t), intent(inout) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: status
        type(newton_report_t) :: report
        real(dp), allocatable :: z(:)

        call flux_map_gauge(fm, derham, quad, profiles, fm%z, fm%gauge)
        z = fm%z
        call flux_map_newton(energy_system, z, fm%fixed, fm%gauge, gauge_at, fm%stationary, &
                             fm%relaxed_steps, fm%verbose, tolerance, max_steps, report, status)
        fm%z = z
        fm%steps = report%steps
        fm%energy = report%energy
        fm%gradient_norm = report%gradient_norm
        fm%converged = report%converged
        fm%negative_eigenvalues = report%negative_eigenvalues
        fm%zero_eigenvalues = report%zero_eigenvalues
        fm%inertia_reliable = report%inertia_reliable

    contains

        subroutine gauge_at(zz, gauge)
            real(dp), intent(in) :: zz(:)
            type(newton_constraints_t), intent(inout) :: gauge

            call flux_map_gauge(fm, derham, quad, profiles, zz, gauge)
        end subroutine gauge_at

        subroutine energy_system(zz, energy, gradient, hessian)
            real(dp), intent(in) :: zz(:)
            real(dp), intent(out) :: energy
            real(dp), intent(out), optional :: gradient(:)
            type(csr_t), intent(out), optional :: hessian
            real(dp), allocatable :: g(:)

            if (present(hessian)) then
                call flux_map_system(fm, derham, quad, profiles, zz, energy, g, hessian)
            else if (present(gradient)) then
                call flux_map_system(fm, derham, quad, profiles, zz, energy, g)
            else
                call flux_map_system(fm, derham, quad, profiles, zz, energy)
            end if
            if (present(gradient)) gradient = g
        end subroutine energy_system
    end subroutine flux_map_solve

end module milca_flux_map
