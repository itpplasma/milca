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
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_cell_xi_basis
    use fortfem_spline_space, only: spline_cell_broken
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_t_basis
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet, csc_triplet_plan_t, &
        csc_triplet_plan_build, csc_triplet_plan_apply
    use fortsparse_status, only: fortsparse_status_t, status_ok
    use milca_eulerian_clebsch, only: clebsch_nfp, clebsch_exact_field
    use fortfem_constrained_newton, only: newton_constraints_t, newton_report_t, &
        newton_constraints_init
    use milca_flux_map, only: flux_map_profiles, flux_map_newton, flux_map_gauge_row
    use milca_flux_map_leaves, only: flux_map_energy_3d, flux_map_energy_3d_value, &
        flux_map_energy_3d_gradient, flux_map_gauge_moment
    implicit none
    private

    type :: axis_section_cell_t
        real(dp), allocatable :: basis(:, :, :), gradient(:, :, :), radius(:), weight(:)
        integer, allocatable :: dofs(:)
    end type axis_section_cell_t

    type :: flux_map3d_hessian_plan_t
        type(csc_triplet_plan_t) :: triplets
        ! The section geometry stays fixed across Newton iterations.
        type(axis_section_cell_t), allocatable :: section(:)
        logical :: ready = .false.
    end type flux_map3d_hessian_plan_t

    type, public :: flux_map3d_t
        integer :: n0 = 0, nt = 0
        logical, allocatable :: is_wall(:) ! section functions with a wall trace
        integer, allocatable :: fixed(:) ! y slots of the wall functions
        real(dp), allocatable :: z(:)
        !! toroidal quadrature points (trapezoid) per plane
        integer :: toroidal_points = 2
        !! Gauge, as in flux_map_t: the rigid rotation of the disc is the one
        !! exact symmetry of the discrete energy and is removed by one moment
        !! constraint that follows the state. Every Newton step keeps the
        !! mean wall angle of plane 1, weighted by the moment of the rotation
        !! generator with each wall function of that plane (the row touches
        !! one plane, so the factorisation stays sparse). .false.: no
        !! constraint.
        logical :: gauge_constraint = .true.
        type(newton_constraints_t) :: gauge
        !! Solver: .false. minimises W (steps accepted on energy decrease);
        !! .true. seeks the stationary point next to the start by Newton
        !! with a line search on the scaled gradient norm, whatever the
        !! signature of the Hessian (an equilibrium need not be a minimum).
        logical :: stationary = .false.
        !! minimisation with relaxed steps, as in flux_map_t
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
    contains
        procedure :: toroidal_quadrature_count
    end type flux_map3d_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    public :: flux_map3d_init, flux_map3d_gauge, flux_map3d_system, flux_map3d_solve, &
              flux_map3d_point, flux_map3d_expand_symmetric, flux_map3d_section_basis

contains

    pure integer function toroidal_quadrature_count(fm, reference_field) result(count)
        class(flux_map3d_t), intent(in) :: fm
        logical, intent(in), optional :: reference_field

        count = fm%toroidal_points*fm%nt + 2
        if (present(reference_field)) then
            if (reference_field) return
        end if
        ! The model integrands are constant when only the zero mode is present.
        if (fm%nt == 1) count = 1
    end function toroidal_quadrature_count

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
        call newton_constraints_init(fm%gauge)
    end subroutine flux_map3d_init

    subroutine flux_map3d_gauge(fm, derham, quad, profiles, plane_nodes, z, gauge, hessian_plan)
        !! The gauge constraint at the state z, as flux_map_gauge, from the
        !! wall functions of plane 1.
        type(flux_map3d_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:)
        type(newton_constraints_t), intent(inout) :: gauge
        type(flux_map3d_hessian_plan_t), intent(inout), optional :: hessian_plan
        real(dp) :: moment(2*fm%n0*fm%nt), energy
        integer :: i

        if (.not. fm%gauge_constraint) then
            call newton_constraints_init(gauge)
            return
        end if
        call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, z, energy, &
                               gauge_moment=moment, hessian_plan=hessian_plan)
        call flux_map_gauge_row(moment, pack([(i, i=1, fm%n0)], fm%is_wall), z, gauge)
    end subroutine flux_map3d_gauge

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
        integer, intent(in), optional :: ix(:), iy(:)
        real(dp), intent(in) :: nv(:), gx(:, :), vel(2), rr, ax(:), ay(:), bx(:), by(:)
        real(dp), intent(out) :: u(8), sig, th, gsig(3), gth(3)
        real(dp), intent(out), optional :: sx(:, :, :), lam(:), mw(:, :)
        real(dp) :: c, g3(3), acc, gacc(3), tacc, cs, sn
        integer :: i, a

        c = 1/(two_pi*rr)
        u = 0
        if (present(sx)) sx = 0
        sig = 0
        gsig = 0
        acc = 0
        gacc = 0
        tacc = 0
        th = 0
        gth = 0
        if (present(lam)) lam = 0
        if (present(mw)) mw = 0
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
                if (present(sx)) then
                    sx(1, ix(i), 1) = nv(i)
                    sx(3:5, ix(i), 1) = g3
                    sx(4, ix(i), 2) = c*nv(i)
                    sx(2, iy(i), 1) = nv(i)
                    sx(6:8, iy(i), 1) = g3
                    sx(7, iy(i), 2) = c*nv(i)
                end if
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
        if (.not. present(sx)) return
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

    subroutine flux_map3d_section_basis(derham, quad, cell, v0, d0, gn)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell
        real(dp), allocatable, intent(out) :: v0(:, :, :), gn(:, :, :)
        integer, allocatable, intent(out) :: d0(:)

        real(dp), allocatable :: jets(:, :, :)
        integer :: nloc, nq, q, p

        nloc = (derham%degree + 1)**2
        nq = quad%points_per_cell
        allocate (d0(nloc), v0(1, nloc, nq), gn(2, nloc, nq), jets(nloc, 6, nq))
        call spline_cell_broken(derham%space, cell, d0)
        d0 = derham%space%broken_dof(derham%space%broken_ptr(d0))
        call spline_cell_xi_basis(derham%space, quad, cell, jets)
        v0(1, :, :) = jets(:, 1, :)
        do q = 1, nq
            p = quad%point_offset(cell) + q - 1
            gn(:, :, q) = matmul(quad%chain(2:3, 2:3, p), transpose(jets(:, 2:3, q)))
        end do
    end subroutine flux_map3d_section_basis

    subroutine cache_axisymmetric_section(derham, quad, nodes, plan)
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: nodes(:, :)
        type(flux_map3d_hessian_plan_t), intent(inout) :: plan
        real(dp) :: f(2, 2), inverse_transpose(2, 2), determinant, position(2)
        integer :: cell, q, nq, point

        nq = quad%points_per_cell
        allocate (plan%section(size(quad%point_offset)))
        !$omp parallel do default(none) shared(plan, derham, quad, nodes, nq) &
        !$omp private(q, point, f, inverse_transpose, determinant, position)
        do cell = 1, size(plan%section)
            associate (table => plan%section(cell))
                call flux_map3d_section_basis(derham, quad, cell, table%basis, table%dofs, table%gradient)
                allocate (table%radius(nq), table%weight(nq))
                do q = 1, nq
                    point = quad%point_offset(cell) + q - 1
                    position = matmul(nodes(:, table%dofs), table%basis(1, :, q))
                    table%radius(q) = position(1)
                    f = matmul(nodes(:, table%dofs), transpose(table%gradient(:, :, q)))
                    call section_inverse_transpose(f, inverse_transpose, determinant)
                    table%gradient(:, :, q) = matmul(inverse_transpose, table%gradient(:, :, q))
                    table%weight(q) = quad%weight(point)*abs(determinant)*two_pi
                end do
            end associate
        end do
        !$omp end parallel do
    end subroutine cache_axisymmetric_section

    pure subroutine section_inverse_transpose(f, inverse_transpose, determinant)
        real(dp), intent(in) :: f(2, 2)
        real(dp), intent(out) :: inverse_transpose(2, 2), determinant

        determinant = f(1, 1)*f(2, 2) - f(1, 2)*f(2, 1)
        inverse_transpose = reshape([f(2, 2), -f(1, 2), -f(2, 1), f(1, 1)], [2, 2])/determinant
    end subroutine section_inverse_transpose

    subroutine flux_map3d_system(fm, derham, quad, profiles, plane_nodes, z, energy, gradient, &
                                 hessian, exact_field, field_error, gauge_moment, hessian_plan)
        !! Energy at z; with gradient and hessian also the derivatives. The
        !! Hessian is symmetric and only its UPPER triangle (columns >= row,
        !! diagonal included) is assembled, by rows (CSR): symmetric by
        !! construction, a third of the work of the full matrix.
        !! flux_map3d_expand_symmetric mirrors it where a full matrix is
        !! needed. derham and quad belong to the section
        !! geometry (plane 1 of the callers here).
        !! plane_nodes(2, n0, nt) are the (R, Z) node
        !! coordinates of every plane. field_error: sum of w |B - B_exact|^2
        !! and of w |B_exact|^2 over the quadrature. gauge_moment: per wall
        !! function b of plane 1 the moment int (-y, x) . du/dTheta_b1 dV of
        !! the rotation generator at z, zero in the other slots
        !! (flux_map3d_gauge).
        type(flux_map3d_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        procedure(clebsch_exact_field), optional :: exact_field
        real(dp), intent(out), optional :: field_error(2), gauge_moment(:)
        type(flux_map3d_hessian_plan_t), target, intent(inout), optional :: hessian_plan
        real(dp), allocatable :: gn(:, :, :), gx(:, :), sx(:, :, :), am(:, :, :, :), &
                                 products(:, :, :), gloc(:), blocks(:, :, :), vals(:), zx(:, :), zy(:, :), &
                                 ax(:), ay(:), bx(:), by(:), lam(:), mw(:, :), gsec(:, :), &
                                 hs(:, :), lm(:), ltab(:, :, :), mloc(:), turn(:), &
                                 trial_jets(:, :), weighted_jets(:, :)
        integer, allocatable :: ix(:), iy(:), idx(:), rows(:), cols(:)
        logical, allocatable :: is_w(:)
        real(dp) :: lt(fm%nt, 2), ht(fm%nt, 2), sk(2, fm%nt), gk(2, 2, fm%nt), sm(2), f(2, 2), &
                    ft(2, 2), detf, vel(2), t, wq, wt, u(8), sig, th, gsig(3), gth(3), s, &
                    tau(3), gp(3), pr(3), w, wg(8), wh(36), h8(8, 8), ref, shift, k11, k12, &
                    vm(3), et(3), c3(3), gs(3), bv(3), bx3(3), ferr(2), sn, cs
        integer :: n0, nt, ntot, nqt, cell, q, p, m, i, j, k, a, b, c1, nloc, nsec, nl, &
                   nz, blk, la, lb, nlmax, ia, ib, ng, nb, pair
        real(dp), pointer, contiguous :: v0(:, :, :)
        integer, pointer, contiguous :: d0(:)
        real(dp), allocatable, target :: v_work(:, :, :)
        integer, allocatable, target :: d_work(:)
        logical :: derivs, cached, cached_section
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info
        interface
            subroutine dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
                import dp
                character, intent(in) :: ta, tb
                integer, intent(in) :: m, n, k, lda, ldb, ldc
                real(dp), intent(in) :: alpha, a(lda, *), b(ldb, *), beta
                real(dp), intent(inout) :: c(ldc, *)
            end subroutine dgemm
        end interface

        n0 = fm%n0
        nt = fm%nt
        ntot = n0*nt
        ! The external field can vary in angle even for an axisymmetric map.
        nqt = fm%toroidal_quadrature_count(reference_field=present(exact_field))
        nb = nt*(nt + 1)/2
        derivs = present(gradient)
        cached = .false.
        cached_section = .false.
        if (present(hessian_plan)) then
            cached = hessian_plan%ready
            cached_section = nt == 1 .and. nqt == 1 .and. allocated(hessian_plan%section)
        end if
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
        ng = 0
        if (present(gauge_moment)) then
            ng = 1
            gauge_moment = 0
        end if
        if (derivs) then
            allocate (gradient(2*ntot), source=0.0_dp)
            if (present(hessian)) then
                nz = size(quad%point_offset)*blk
                allocate (vals(nz + size(fm%fixed)))
                if (.not. cached) then
                    allocate (rows(size(vals)), cols(size(vals)))
                    rows = 1
                    cols = 1
                end if
                vals = 0
                allocate (products(nqt, 4, nb))
                pair = 0
                do lb = 1, nt
                    do la = 1, lb
                        pair = pair + 1
                        products(:, 1, pair) = ltab(la, 1, :)*ltab(lb, 1, :)
                        products(:, 2, pair) = ltab(la, 1, :)*ltab(lb, 2, :)
                        products(:, 3, pair) = ltab(la, 2, :)*ltab(lb, 1, :)
                        products(:, 4, pair) = ltab(la, 2, :)*ltab(lb, 2, :)
                    end do
                end do
            end if
        end if
        !$omp parallel do schedule(dynamic) default(private) &
        !$omp shared(fm, derham, quad, plane_nodes, z, gradient, rows, cols, vals, n0, nt, &
        !$omp        ntot, nqt, blk, derivs, clebsch_nfp, hessian, field_error, ltab, ng, &
        !$omp        gauge_moment, nb, products, cached, cached_section, hessian_plan) &
        !$omp reduction(+:energy, ferr)
        do cell = 1, size(quad%point_offset)
            if (cached_section) then
                v0 => hessian_plan%section(cell)%basis
                d0 => hessian_plan%section(cell)%dofs
            else
                call flux_map3d_section_basis(derham, quad, cell, v_work, d_work, gn)
                v0 => v_work
                d0 => d_work
            end if
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
            allocate (sx(8, nsec, 2), gsec(nsec, 2), hs(8, nsec), idx(nl), gloc(nl))
            if (present(hessian)) then
                if (allocated(am)) deallocate (am, blocks)
                allocate (am(nsec, nsec, nqt, 4), blocks(nsec, nsec, nb))
                am = 0
                if (cached_section .and. .not. allocated(trial_jets)) then
                    allocate (trial_jets(8*quad%points_per_cell, 2*nloc), &
                              weighted_jets(8*quad%points_per_cell, 2*nloc))
                end if
            end if
            gloc = 0
            if (ng > 0) then
                if (allocated(mloc)) deallocate (mloc, turn)
                allocate (mloc(nsec), turn(nsec))
                mloc = 0
            end if
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
                if (.not. cached_section) then
                    p = quad%point_offset(cell) + q - 1
                    wq = quad%weight(p)
                    do k = 1, nt
                        sk(:, k) = matmul(plane_nodes(:, d0, k), v0(1, :, q))
                        gk(:, :, k) = matmul(plane_nodes(:, d0, k), transpose(gn(:, :, q)))
                    end do
                end if
                do m = 1, nqt
                    t = (m - 1)/real(nqt, dp)
                    lt = ltab(:, :, m)
                    if (cached_section) then
                        sm(1) = hessian_plan%section(cell)%radius(q)
                        gx = hessian_plan%section(cell)%gradient(:, :, q)
                        wt = hessian_plan%section(cell)%weight(q)
                        vel = 0
                    else
                        sm = matmul(sk, lt(:, 1))
                        vel = matmul(sk, lt(:, 2))
                        f = 0
                        do k = 1, nt
                            f = f + lt(k, 1)*gk(:, :, k)
                        end do
                        call section_inverse_transpose(f, ft, detf)
                        gx = matmul(ft, gn(:, :, q))
                        wt = wq*abs(detf)*two_pi/nqt
                    end if
                    ax = matmul(zx, lt(:, 1))
                    ay = matmul(zy, lt(:, 1))
                    bx = matmul(zx, lt(:, 2))
                    by = matmul(zy, lt(:, 2))
                    if (derivs .or. ng > 0) then
                        call locals(is_w, ix, iy, v0(1, :, q), gx, vel, sm(1), ax, ay, bx, by, u, &
                                    sx, sig, th, gsig, gth, lam, mw)
                    else
                        call locals(is_w, nv=v0(1, :, q), gx=gx, vel=vel, rr=sm(1), &
                                    ax=ax, ay=ay, bx=bx, by=by, u=u, sig=sig, th=th, gsig=gsig, gth=gth)
                    end if
                    s = u(1)**2 + u(2)**2
                    call profiles(s, tau, gp, pr)
                    if (present(hessian)) then
                        call flux_map_energy_3d(sm(1), u(1), u(2), u(3), u(4), u(5), u(6), u(7), &
                            u(8), tau(1), tau(2), tau(3), gp(1), gp(2), gp(3), pr(1), pr(2), pr(3), &
                            w, wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8), wh(1), &
                            wh(2), wh(3), wh(4), wh(5), wh(6), wh(7), wh(8), wh(9), wh(10), &
                            wh(11), wh(12), wh(13), wh(14), wh(15), wh(16), wh(17), wh(18), &
                            wh(19), wh(20), wh(21), wh(22), wh(23), wh(24), wh(25), wh(26), &
                            wh(27), wh(28), wh(29), wh(30), wh(31), wh(32), wh(33), wh(34), &
                            wh(35), wh(36))
                    else if (derivs) then
                        call flux_map_energy_3d_gradient(sm(1), u(1), u(2), u(3), u(4), &
                            u(5), u(6), u(7), u(8), tau(1), tau(2), tau(3), gp(1), gp(2), gp(3), &
                            pr(1), pr(2), pr(3), w, wg(1), wg(2), wg(3), wg(4), wg(5), wg(6), wg(7), wg(8))
                    else
                        call flux_map_energy_3d_value(sm(1), u(1), u(2), u(3), u(4), &
                            u(5), u(6), u(7), u(8), tau(1), gp(1), pr(1), w)
                    end if
                    energy = energy + wt*w
                    if (ng > 0) then
                        ! the wall functions of plane 1
                        turn = 0
                        do i = 1, nloc
                            if (is_w(i)) call flux_map_gauge_moment(u(1), u(2), &
                                sx(1, ix(i), 1), sx(2, ix(i), 1), turn(ix(i)))
                        end do
                        mloc = mloc + wt*sm(1)*lt(1, 1)*turn
                    end if
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
                    ! blocks (1,1), (1,2), (2,2); (2,1) is the transpose of (1,2)
                    hs = matmul(h8, sx(:, :, 1))
                    if (cached_section) then
                        trial_jets(8*q - 7:8*q, :nsec) = sx(:, :, 1)
                        weighted_jets(8*q - 7:8*q, :nsec) = wt*hs
                    else
                        call dgemm('T', 'N', nsec, nsec, 8, wt, sx(1, 1, 1), 8, &
                                   hs, 8, 1.0_dp, am(1, 1, m, 1), nsec)
                    end if
                    ! The zero Fourier mode has no toroidal derivatives.
                    if (nt > 1) then
                        hs = matmul(h8, sx(:, :, 2))
                        call dgemm('T', 'N', nsec, nsec, 8, wt, sx(1, 1, 1), 8, &
                                   hs, 8, 1.0_dp, am(1, 1, m, 2), nsec)
                        call dgemm('T', 'N', nsec, nsec, 8, wt, sx(1, 1, 2), 8, &
                                   hs, 8, 1.0_dp, am(1, 1, m, 4), nsec)
                    end if
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
                                am(a, b, m, 1) = am(a, b, m, 1) + wt*(k11*lam(i)*lam(j) &
                                                 + lam(i)*lm(j) + lam(j)*lm(i))
                                if (nt > 1) am(a, b, m, 2) = am(a, b, m, 2) + wt*k12*lam(i)*lam(j)
                            end do
                        end do
                    end if
                end do
            end do
            if (ng > 0) then
                !$omp critical (flux_map3d_moments)
                gauge_moment(idx(:nsec)) = gauge_moment(idx(:nsec)) + mloc
                !$omp end critical (flux_map3d_moments)
            end if
            if (.not. derivs) cycle
            do i = 1, nl
                !$omp atomic
                gradient(idx(i)) = gradient(idx(i)) + gloc(i)
            end do
            if (.not. present(hessian)) cycle
            if (cached_section) then
                ! Stack the point jets to assemble the cell with one matrix product.
                blocks(:, :, 1) = am(:, :, 1, 1)
                call dgemm('T', 'N', nsec, nsec, size(trial_jets, 1), 1.0_dp, &
                           trial_jets, size(trial_jets, 1), weighted_jets, size(weighted_jets, 1), &
                           1.0_dp, blocks(1, 1, 1), nsec)
            else if (nt == 1) then
                blocks(:, :, 1) = sum(am(:, :, :, 1), dim=3)
            else
                ! Batch the contractions over section and upper-triangular plane pairs.
                do m = 1, nqt
                    am(:, :, m, 3) = transpose(am(:, :, m, 2))
                end do
                call dgemm('N', 'N', nsec**2, nb, 4*nqt, 1.0_dp, am(1, 1, 1, 1), nsec**2, &
                           products(1, 1, 1), 4*nqt, 0.0_dp, blocks(1, 1, 1), nsec**2)
            end if
            ! local entries (ia, ib), ia <= ib (a pair once; the diagonal blocks
            ! hold both orders, only ia <= ib is taken), to global (min, max)
            nz = (cell - 1)*blk
            do lb = 1, nt
                do j = 1, nsec
                    ib = (lb - 1)*nsec + j
                    do la = 1, lb
                        pair = lb*(lb - 1)/2 + la
                        do i = 1, merge(j, nsec, la == lb)
                            ia = (la - 1)*nsec + i
                            nz = nz + 1
                            if (.not. cached) then
                                rows(nz) = min(idx(ia), idx(ib))
                                cols(nz) = max(idx(ia), idx(ib))
                            end if
                            vals(nz) = blocks(i, j, pair)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if (present(field_error)) field_error = ferr
        nz = size(quad%point_offset)*blk
        if (.not. present(hessian)) return
        shift = maxval(abs(vals(:nz)))
        do k = 1, size(fm%fixed)
            nz = nz + 1
            if (.not. cached) then
                rows(nz) = fm%fixed(k)
                cols(nz) = fm%fixed(k)
            end if
            vals(nz) = shift
        end do
        if (present(hessian_plan)) then
            if (.not. cached) then
                call csc_triplet_plan_build(2*ntot, 2*ntot, cols, rows, hessian_plan%triplets, info)
                if (.not. status_ok(info)) error stop 'flux map 3D Hessian pattern'
                deallocate (rows, cols)
                hessian_plan%ready = .true.
            end if
            call csc_triplet_plan_apply(hessian_plan%triplets, vals, mm, info)
        else
            call csc_from_triplet(2*ntot, 2*ntot, cols, rows, vals, mm, info)
        end if
        if (.not. status_ok(info)) error stop 'flux map 3D Hessian triplets'
        hessian%nrow = 2*ntot
        hessian%ncol = 2*ntot
        call move_alloc(mm%col_ptr, hessian%row_ptr)
        call move_alloc(mm%row_idx, hessian%col)
        call move_alloc(mm%val, hessian%val)
    end subroutine flux_map3d_system

    subroutine flux_map3d_point(fm, plane_nodes, d0, nv, gn, t, profiles, rz, u, b, label_jacobian)
        !! State at one point of the section and the toroidal parameter t:
        !! d0, nv, gn are the section functions that do not vanish there,
        !! their values and their gradients in the geometry of the reference
        !! plane (as from spline_point_matrix). Returns the position (R, Z),
        !! u = (x, y, grad x, grad y) and B in components (R, phi, Z).
        !! label_jacobian gives the derivatives of (x, y) with respect to the
        !! coordinates (R, Z) of the reference plane at fixed t.
        type(flux_map3d_t), intent(in) :: fm
        real(dp), intent(in) :: plane_nodes(:, :, :), nv(:), gn(:, :), t
        integer, intent(in) :: d0(:)
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out) :: rz(2), u(8), b(3)
        real(dp), intent(out), optional :: label_jacobian(2, 2)
        real(dp) :: lt(fm%nt, 2), ht(fm%nt, 2), f(2, 2), ft(2, 2), detf, vel(2), sig, th, &
                    gsig(3), gth(3), tau(3), gp(3), pr(3), c3(3), gs(3), ref, shift
        real(dp) :: gx(2, size(d0)), ax(size(d0)), ay(size(d0)), bx(size(d0)), by(size(d0)), &
                    zx(size(d0), fm%nt), &
                    zy(size(d0), fm%nt)
        integer :: i, k, ntot
        logical :: is_w(size(d0))

        ntot = fm%n0*fm%nt
        call spline_derham_t_basis(fm%nt, t, lt, ht)
        lt(:, 2) = clebsch_nfp*lt(:, 2)
        is_w = fm%is_wall(d0)
        ref = 0
        do i = size(d0), 1, -1
            if (is_w(i)) ref = fm%z(d0(i))
        end do
        do i = 1, size(d0)
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
        call locals(is_w, nv=nv, gx=gx, vel=vel, rr=rz(1), &
                    ax=ax, ay=ay, bx=bx, by=by, u=u, sig=sig, th=th, gsig=gsig, gth=gth)
        if (present(label_jacobian)) then
            label_jacobian(1, :) = matmul(u([3, 5]), f)
            label_jacobian(2, :) = matmul(u([6, 8]), f)
        end if
        call profiles(u(1)**2 + u(2)**2, tau, gp, pr)
        c3 = [u(4)*u(8) - u(5)*u(7), u(5)*u(6) - u(3)*u(8), u(3)*u(7) - u(4)*u(6)]
        gs = 2*(u(1)*u(3:5) + u(2)*u(6:8))
        b = -[tau(1)*c3(1) + gp(1)*gs(3)/rz(1), tau(1)*c3(2), tau(1)*c3(3) - gp(1)*gs(1)/rz(1)]
    end subroutine flux_map3d_point

    subroutine flux_map3d_solve(fm, derham, quad, profiles, plane_nodes, tolerance, max_steps, &
                                status)
        !! Solves from fm%z by flux_map_newton (the one solver loop of the
        !! flux map), with the gauge constraints built at fm%z. Status 1: no
        !! acceptable step, 2: step limit.
        type(flux_map3d_t), intent(inout) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), tolerance
        integer, intent(in) :: max_steps
        integer, intent(out) :: status
        type(newton_report_t) :: report
        type(flux_map3d_hessian_plan_t) :: hessian_plan
        real(dp), allocatable :: z(:)

        if (fm%nt == 1) call cache_axisymmetric_section(derham, quad, plane_nodes(:, :, 1), hessian_plan)
        call flux_map3d_gauge(fm, derham, quad, profiles, plane_nodes, fm%z, fm%gauge, hessian_plan)
        z = fm%z
        call flux_map_newton(energy_system, z, fm%fixed, fm%gauge, gauge_at, fm%stationary, &
                             fm%relaxed_steps, fm%verbose, tolerance, max_steps, report, status, &
                             model_agreement_damping=.true., initial_damping=0.1_dp)
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

            call flux_map3d_gauge(fm, derham, quad, profiles, plane_nodes, zz, gauge, hessian_plan)
        end subroutine gauge_at

        subroutine energy_system(zz, energy, gradient, hessian)
            real(dp), intent(in) :: zz(:)
            real(dp), intent(out) :: energy
            real(dp), intent(out), optional :: gradient(:)
            type(csr_t), intent(out), optional :: hessian
            real(dp), allocatable :: g(:)

            if (present(hessian)) then
                call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, zz, energy, g, &
                                       hessian, hessian_plan=hessian_plan)
            else if (present(gradient)) then
                call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, zz, energy, g, hessian_plan=hessian_plan)
            else
                call flux_map3d_system(fm, derham, quad, profiles, plane_nodes, zz, energy, hessian_plan=hessian_plan)
            end if
            if (present(gradient)) gradient = g
        end subroutine energy_system
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

end module milca_flux_map_3d
