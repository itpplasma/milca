!> Monoenergetic drift-kinetic coefficients D11, D31, D13, D33 on an
!! axisymmetric surface, the equation of milca_dk_mono, by continuous
!! spectral elements on a mesh of the (theta, xi) cylinder that is fitted to
!! the orbit classes.
!!
!! lambda = (1 - xi^2)/b is constant along orbits, so the layers of the
!! solution lie on level curves of lambda, above all the trapped-passing
!! boundary lambda_c = 1/b_max, whose width shrinks like sqrt(nu). Element
!! rows are bounded by such curves: in the passing region by levels
!! xi = sigma sqrt(1 - Lambda b(theta)), in the trapped region by scaled
!! separatrices xi = zeta xi_c(theta), xi_c = sqrt(1 - b/b_max), so that
!! zeta = +-1 is the separatrix itself. Columns are theta intervals of one
!! period that starts and ends at the maximum theta_x of b, where the
!! separatrix closes in the X point. Trapped elements next to theta_x
!! collapse to that point; its nodes are one unknown, which also makes
!! the discrete functions regular there. Rows are graded geometrically
!! towards the separatrix from both sides and columns towards theta_x.
!! Reflection points need no treatment in (theta, xi): the solution is
!! smooth there, unlike in (theta, lambda).
!!
!! Each element carries tensor Lagrange polynomials of degree p on
!! Gauss-Lobatto nodes in the reference square, mapped to the element by
!! theta linear in the first and xi linear between the bounding curves in
!! the second coordinate. The weak form uses the measure (1/b^2) dtheta dxi/2
!! of the flux-surface and pitch averages, in which streaming is
!! divergence free; there is no inflow boundary, so no upwinding. The
!! density null space is bordered as in milca_dk_mono, and the system is
!! factorised once by a sparse direct solver; adjoint solves reuse the
!! factors through the transpose.
module milca_dk_hp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre, gauss_lobatto_legendre
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use fortnum_special_legendre, only: legendre_p
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, sparse_solver_t, sparse_factor, &
        sparse_solve_op, sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_LU_SUPERNODAL, FORTSPARSE_ORDER_AMD
    use milca_dk_hp_leaves, only: dk_hp_point, dk_hp_stream, dk_hp_ntv
    use milca_dk_mono, only: mono_coefficients_t
    implicit none
    private
    public :: hp_geometry_t, hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, hp_system_t, &
        hp_assemble, hp_factor, hp_solve, hp_coefficients, hp_free, hp_estimate, &
        hp_refine, hp_operator_t, hp_ntv_kernel, hp_ntv_estimate, hp_ntv_factor_t, &
        hp_ntv_factor, hp_ntv_solve, hp_ntv_free, hp_terms_t, hp_ntv_terms, &
        hp_node_matrix, hp_mass, hp_mesh_join, hp_locate, hp_lines_t, hp_lines, &
        hp_line_ops
    integer, parameter, public :: HP_OK = 0, HP_BAD_INPUT = 1, HP_SINGULAR = 2, &
        HP_NO_CONVERGENCE = 3

    real(dp), parameter :: pi = 4*atan(1.0_dp)


    !> Operator per unit speed: b.grad theta = k b**k_power, nu = nu_D/v; for
    !! one toroidal harmonic n of the NTV problem (milca_dk_ntv) also q, the
    !! speed v, omega_E, the precession prefactor pre and offset, the drift
    !! weight sign wsign and the perturbation c_m exp(i m theta).
    type :: hp_operator_t
        real(dp) :: k = 0, nu = 0, q = 0, v = 1, omega_e = 0, pre = 0, offset = 0, &
            wsign = -1
        integer :: k_power = 0, n = 0
        integer, allocatable :: m(:)
        complex(dp), allocatable :: c(:)
        ! The two further NEO-2-QL sources of milca_dk_ntv (thin fields only):
        ! the parallel-derivative part of the perturbed radial drift and the
        ! perturbed operator acting on the axisymmetric solution f_10.
        logical :: source_dpar = .false., source_dlf10 = .false.
    end type

    !> The orbit geometry of the mesh: on the (theta, xi) cylinder of one
    !> surface, the stream function of the orbits (its bracket is the
    !> streaming term in the measure of the weak form), the row curves
    !> (passing levels of sign kind = +-1 at level val, scaled separatrices
    !> kind = 0 at val = zeta), and the pointwise coefficients of the weak form
    !> per unit dtheta dxi. The separatrix closes at theta_x at the passing
    !> level lambda_c.
    type, abstract :: hp_geometry_t
        real(dp) :: kappa = 0, theta_x = 0, lambda_c = 0
    contains
        procedure(stream_i), deferred :: stream
        procedure(level_i), deferred :: level
        procedure(point_i), deferred :: point
        procedure :: perturbation => op_perturbation
        ! perturbation and point at one (theta, xi), for geometries that
        ! share work between them.
        procedure :: point_source => field_point_source
        ! The flux label s of the orbit point at (theta, xi), or -1 where the
        ! point lies on the surface itself (thin orbits).
        procedure :: position => field_position
    end type

    abstract interface
        real(dp) function stream_i(g, th, xi)
            import :: dp, hp_geometry_t
            class(hp_geometry_t), intent(in) :: g
            real(dp), intent(in) :: th, xi
        end function stream_i

        subroutine level_i(g, kind, val, th, xi, dxi)
            import :: dp, hp_geometry_t
            class(hp_geometry_t), intent(in) :: g
            integer, intent(in) :: kind
            real(dp), intent(in) :: val, th
            real(dp), intent(out) :: xi, dxi
        end subroutine level_i

        !> collide: coefficient of d_xi w d_xi g; drive_r, drive_e: the radial
        !> and parallel drives of the axisymmetric problem; density: the
        !> measure; phase, precess, drive, torque: the NTV terms of
        !> milca_dk_hp_leaves (dk_hp_ntv).
        subroutine point_i(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
                precess, drive, torque)
            import :: dp, hp_geometry_t, hp_operator_t
            class(hp_geometry_t), intent(in) :: g
            real(dp), intent(in) :: th, xi
            type(hp_operator_t), intent(in) :: op
            real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, &
                drive, torque
        end subroutine point_i
    end interface

    !> b(theta) = a_0 + sum_m a_m cos(m theta) + c_m sin(m theta), the
    !! trigonometric interpolant of the samples, with its single maximum
    !! b_max at theta_x and minimum b_min: the thin-orbit geometry, whose
    !! stream function is -lambda/2, lambda = (1 - xi^2)/b.
    type, extends(hp_geometry_t) :: hp_field_t
        integer :: m = 0
        real(dp) :: b_max = 0, b_min = 0
        real(dp), allocatable :: a(:), c(:)
    contains
        procedure :: stream => field_stream
        procedure :: level => curve
        procedure :: point => field_point
    end type

    !> Column breakpoints theta(1) = theta_x - 2 pi < ... < theta(end) =
    !! theta_x; passing levels lam(1) = 0 < ... < lam(end) = lambda_c;
    !! trapped rows zeta(1) = -1 < ... < zeta(end) = 1; degree p.
    type :: hp_mesh_t
        integer :: p = 0
        real(dp), allocatable :: theta(:), lam(:), zeta(:)
    end type


    !> Assembled system of n pitch unknowns with its element matrices
    !! (columns of mat), loads and density weights for residual estimates;
    !! loc maps local element nodes to unknowns.
    !> Equilibrated sparse factor of an NTV system: A, its row and column
    !> scales and the factor of R A C.
    type :: hp_ntv_factor_t
        type(csc_z_t) :: a
        real(dp), allocatable :: rs(:), cs(:)
        type(sparse_solver_t) :: solver
    end type

    interface hp_ntv_factor
        module procedure factor_elements, factor_matrix
    end interface

    type :: hp_system_t
        integer :: n = 0
        real(dp) :: norm = 0 ! int dtheta dxi/(b b.grad theta)
        real(dp), allocatable :: density(:), mat(:, :), dens(:, :)
        complex(dp), allocatable :: drive(:, :), load(:, :, :)
        integer, allocatable :: loc(:, :)
        type(csc_t) :: a
        type(sparse_solver_t) :: solver
        logical :: factored = .false.
    end type

    !> The NTV system of one mesh with its terms apart, so that systems at
    !> other nu, v, omega_E and pre need no new quadrature. Element matrices
    !> (column major, columns of the arrays): the streaming bracket s, pitch
    !> scattering per unit nu c, the factor ph of -i n q, the mass m in the
    !> measure and the factor pr of i n pre/v. On the n unknowns: the drive
    !> load3 without its 1/v, the torque weight load4, the radial drive load1
    !> of the axisymmetric problem and the density weight; theta and xi of
    !> the nodes. With source_dlf10 the element matrices mir and ls give the
    !> load (mir + nu ls) f_10 of f_10 on the same mesh. The quadrature of the
    !> mass: angles q_theta(t), t = (column - 1) nq + i, and at each the points
    !> q_xi(j, t), j = (row - 1) nq + i, with weights q_w (measure included),
    !> nq = 2 p + 3; m is the sum of q_w phi_a phi_b over them.
    type :: hp_terms_t
        integer :: n = 0, p = 0
        real(dp) :: norm = 0
        integer, allocatable :: loc(:, :)
        real(dp), allocatable :: s(:, :), c(:, :), ph(:, :), m(:, :), pr(:, :), load1(:), &
            density(:), theta(:), xi(:), q_theta(:), q_xi(:, :), q_w(:, :)
        complex(dp), allocatable :: load3(:), load4(:), mir(:, :), ls(:, :)
    end type

    !> Angles theta with weights w, sum w = 2 pi, on which meshes of one
    !> surface exchange pitch moments.
    type :: hp_lines_t
        real(dp), allocatable :: theta(:), w(:)
    end type

contains

    !> Field from b on n (odd) uniform points theta_j = 2 pi (j - 1)/n and
    !! kappa = b.grad theta / b; refuses fields with more than one maximum.
    subroutine hp_field(b, kappa, f, status)
        real(dp), intent(in) :: b(:), kappa
        type(hp_field_t), intent(out) :: f
        integer, intent(out) :: status
        integer :: n, k, j, n_fine, n_max, n_min
        real(dp) :: th, d0, d1
        real(dp), allocatable :: theta(:)

        status = HP_BAD_INPUT
        n = size(b)
        if (mod(n, 2) == 0 .or. n < 3 .or. any(b <= 0) .or. &
            abs(kappa) <= tiny(kappa)) return
        f%m = (n - 1)/2
        f%kappa = kappa
        allocate (f%a(0:f%m), f%c(f%m), theta(n))
        theta = [(2*pi*(j - 1)/n, j=1, n)]
        f%a(0) = sum(b)/n
        do k = 1, f%m
            f%a(k) = 2*sum(b*cos(k*theta))/n
            f%c(k) = 2*sum(b*sin(k*theta))/n
        end do
        ! Extrema from sign changes of b' on a fine grid, polished by Newton.
        n_fine = 64*n
        n_max = 0; n_min = 0
        d0 = field_db(f, 0.0_dp)
        do j = 1, n_fine
            th = 2*pi*j/n_fine
            d1 = field_db(f, th)
            if (d0 > 0 .and. d1 <= 0) then
                n_max = n_max + 1
                f%theta_x = newton_extremum(f, th - pi/n_fine)
            else if (d0 < 0 .and. d1 >= 0) then
                n_min = n_min + 1
                f%b_min = field_b(f, newton_extremum(f, th - pi/n_fine))
            end if
            d0 = d1
        end do
        if (n_max /= 1 .or. n_min /= 1) return
        f%b_max = field_b(f, f%theta_x)
        f%lambda_c = 1/f%b_max
        status = HP_OK
    end subroutine hp_field

    !> Geometric mesh: n_bulk passing levels and trapped rows in the bulk,
    !! layers towards the separatrix down to the relative width delta with
    !! ratio q, and n_col columns graded towards theta_x down to sqrt(delta).
    subroutine hp_default_mesh(f, p, n_col, n_bulk, delta, q, mesh)
        class(hp_geometry_t), intent(in) :: f
        integer, intent(in) :: p, n_col, n_bulk
        real(dp), intent(in) :: delta, q
        type(hp_mesh_t), intent(out) :: mesh
        real(dp), allocatable :: layer(:), col(:)
        real(dp) :: lc, h
        integer :: j, n_l

        mesh%p = p
        lc = f%lambda_c
        ! Relative distances to the separatrix: bulk spacing, then q^k.
        n_l = max(0, ceiling(log(delta/0.25_dp)/log(q)))
        layer = [(0.25_dp*q**j, j=0, n_l)]
        mesh%lam = [(lc*(1 - 0.25_dp)*j/n_bulk, j=0, n_bulk), &
            (lc*(1 - layer(j)), j=2, size(layer)), lc]
        mesh%zeta = [-1.0_dp, (-(1 - layer(j)), j=size(layer), 2, -1), &
            (-(1 - 0.25_dp)*(n_bulk - j)/real(n_bulk, dp), j=0, n_bulk), &
            ((1 - 0.25_dp)*j/real(n_bulk, dp), j=1, n_bulk), &
            ((1 - layer(j)), j=2, size(layer)), 1.0_dp]
        ! Columns: uniform interior, geometric towards both ends.
        h = 2*pi/n_col
        n_l = max(0, ceiling(log(0.1_dp*sqrt(delta)/0.5_dp)/log(q)))
        col = [(h*0.5_dp*q**j, j=n_l, 1, -1)]
        mesh%theta = f%theta_x - 2*pi + [0.0_dp, col, &
            (h*(0.5_dp + j), j=0, n_col - 2), 2*pi - h*0.5_dp, &
            (2*pi - col(j), j=size(col), 1, -1), 2*pi]
    end subroutine hp_default_mesh

    !> Element matrices, loads and the map of local nodes to unknowns:
    !! periodic in theta, the trapped nodes on theta_x (the X point) one
    !! unknown. mat holds each element matrix column major.
    subroutine elements(f, mesh, op, sys, mat, status, f10, terms)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_system_t), intent(inout) :: sys
        complex(dp), allocatable, intent(out) :: mat(:, :)
        integer, intent(out) :: status
        ! f_10 per unit speed on the unknowns of the same mesh (source_dlf10).
        real(dp), intent(in), optional :: f10(:)
        ! The split of mat and the f_10 operators (never together with f10).
        type(hp_terms_t), intent(inout), optional :: terms
        integer :: p, n_col, n_row, n_th, ec, er, e, a1, a2, n_loc, n_el, c, k, i
        integer, allocatable :: dof(:, :), curve_kind(:)
        real(dp), allocatable :: curve_val(:), node(:), wn(:), xq(:), wq(:), xi_n(:, :), &
            xi_q(:, :), parts(:, :), qx(:, :), qw(:, :)
        complex(dp), allocatable :: f10_ops(:, :)
        real(dp) :: th_a, th, d

        status = HP_BAD_INPUT
        p = mesh%p
        if (p < 1 .or. op%nu <= 0 .or. op%k == 0) return
        n_col = size(mesh%theta) - 1
        call row_curves(mesh, curve_kind, curve_val)
        n_row = size(curve_val) - 1
        n_th = n_col*p
        call numbering(mesh, dof, sys%n)

        n_loc = (p + 1)**2
        n_el = n_col*n_row
        allocate (sys%loc(n_loc, n_el), sys%load(n_loc, 4, n_el), &
            sys%dens(n_loc, n_el), mat(n_loc**2, n_el))
        if (present(terms)) then
            allocate (terms%s(n_loc**2, n_el), terms%c(n_loc**2, n_el), &
                terms%ph(n_loc**2, n_el), terms%m(n_loc**2, n_el), &
                terms%pr(n_loc**2, n_el), terms%theta(sys%n), terms%xi(sys%n))
            if (op%source_dlf10) allocate (terms%mir(n_loc**2, n_el), &
                terms%ls(n_loc**2, n_el))
            allocate (terms%q_theta(n_col*(2*p + 3)), terms%q_xi(n_row*(2*p + 3), &
                n_col*(2*p + 3)), terms%q_w(n_row*(2*p + 3), n_col*(2*p + 3)))
        end if
        allocate (node(p + 1), wn(p + 1), xq(2*p + 3), wq(2*p + 3))
        call gauss_lobatto_legendre(p + 1, node, wn)
        call gauss_legendre(2*p + 3, xq, wq)
        ! Column by column: every row curve once at the column's nodes and
        ! quadrature angles, shared by the elements above and below it.
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(e, ec, er, a1, a2, c, k, th_a, th, d, xi_n, xi_q, parts, f10_ops, qx, qw)
        do ec = 1, n_col
            allocate (xi_n(n_row + 1, p + 1), xi_q(n_row + 1, 2*p + 3))
            if (present(terms)) allocate (parts(n_loc**2, 5), f10_ops(n_loc**2, 2), &
                qx(2*p + 3, 2*p + 3), qw(2*p + 3, 2*p + 3))
            th_a = (mesh%theta(ec + 1) - mesh%theta(ec))/2
            if (present(terms)) terms%q_theta((ec - 1)*(2*p + 3) + 1:ec*(2*p + 3)) = &
                mesh%theta(ec) + th_a*(1 + xq)
            do c = 1, n_row + 1
                do a1 = 1, p + 1
                    th = mesh%theta(ec) + th_a*(1 + node(a1))
                    call f%level(curve_kind(c), curve_val(c), th, xi_n(c, a1), d)
                end do
                do a1 = 1, 2*p + 3
                    th = mesh%theta(ec) + th_a*(1 + xq(a1))
                    call f%level(curve_kind(c), curve_val(c), th, xi_q(c, a1), d)
                end do
            end do
            do er = 1, n_row
                e = ec + (er - 1)*n_col
                do a2 = 0, p
                    do a1 = 0, p
                        sys%loc(1 + a1 + (p + 1)*a2, e) = dof(mod((ec - 1)*p + a1, n_th), &
                            (er - 1)*p + a2)
                    end do
                end do
                if (present(terms)) then
                    call element(f, mesh, op, xi_n(er:er + 1, :), xi_q(er:er + 1, :), &
                        ec, sys%load(:, :, e), sys%dens(:, e), mat(:, e), parts=parts, &
                        f10_ops=f10_ops, qxi=qx, qw=qw)
                    terms%q_xi((er - 1)*(2*p + 3) + 1:er*(2*p + 3), (ec - 1)*(2*p + 3) &
                        + 1:ec*(2*p + 3)) = transpose(qx)
                    terms%q_w((er - 1)*(2*p + 3) + 1:er*(2*p + 3), (ec - 1)*(2*p + 3) &
                        + 1:ec*(2*p + 3)) = transpose(qw)
                    terms%s(:, e) = parts(:, 1)
                    terms%c(:, e) = parts(:, 2)
                    terms%ph(:, e) = parts(:, 3)
                    terms%m(:, e) = parts(:, 4)
                    terms%pr(:, e) = parts(:, 5)
                    if (op%source_dlf10) then
                        terms%mir(:, e) = f10_ops(:, 1)
                        terms%ls(:, e) = f10_ops(:, 2)
                    end if
                    ! Each column owns its first p node angles, so every
                    ! unknown is written by one thread.
                    do a2 = 0, p
                        do a1 = 0, p - 1
                            k = sys%loc(1 + a1 + (p + 1)*a2, e)
                            terms%theta(k) = mesh%theta(ec) + th_a*(1 + node(a1 + 1))
                            terms%xi(k) = xi_n(er, a1 + 1) + (xi_n(er + 1, a1 + 1) &
                                - xi_n(er, a1 + 1))*(1 + node(a2 + 1))/2
                        end do
                    end do
                else if (present(f10)) then
                    call element(f, mesh, op, xi_n(er:er + 1, :), xi_q(er:er + 1, :), &
                        ec, sys%load(:, :, e), sys%dens(:, e), mat(:, e), f10(sys%loc(:, e)))
                else
                    call element(f, mesh, op, xi_n(er:er + 1, :), xi_q(er:er + 1, :), &
                        ec, sys%load(:, :, e), sys%dens(:, e), mat(:, e))
                end if
            end do
            deallocate (xi_n, xi_q)
            if (present(terms)) deallocate (parts, f10_ops, qx, qw)
        end do
        !$omp end parallel do
        allocate (sys%drive(sys%n, 4), sys%density(sys%n))
        sys%drive = 0; sys%density = 0
        ! Node by node: at the X point several nodes of one element share an
        ! unknown, and a vector-subscript update would keep only one of them.
        do e = 1, n_el
            do i = 1, n_loc
                sys%drive(sys%loc(i, e), :) = sys%drive(sys%loc(i, e), :) + sys%load(i, :, e)
                sys%density(sys%loc(i, e)) = sys%density(sys%loc(i, e)) + sys%dens(i, e)
            end do
        end do
        sys%norm = sum(sys%density)
        status = HP_OK
    end subroutine elements

    !> The row curves from xi = -1 to xi = 1: passing levels (sigma = -1),
    !> the trapped rows, passing levels (sigma = +1).
    pure subroutine row_curves(mesh, kind, val)
        type(hp_mesh_t), intent(in) :: mesh
        integer, allocatable, intent(out) :: kind(:)
        real(dp), allocatable, intent(out) :: val(:)
        integer :: n_lev, i
        n_lev = size(mesh%lam)
        val = [mesh%lam(1:n_lev - 1), mesh%zeta, mesh%lam(n_lev - 1:1:-1)]
        kind = [(-1, i=1, n_lev - 1), (0, i=1, size(mesh%zeta)), (1, i=1, n_lev - 1)]
    end subroutine row_curves

    !> The unknown of node (i, j), i the node angle of the period and j the
    !> row node: the trapped nodes on theta_x (i = 0) are one unknown.
    pure subroutine numbering(mesh, dof, n)
        type(hp_mesh_t), intent(in) :: mesh
        integer, allocatable, intent(out) :: dof(:, :)
        integer, intent(out) :: n
        integer :: n_th, n_jr, j_lo, j_hi, i, j
        n_th = (size(mesh%theta) - 1)*mesh%p
        j_lo = (size(mesh%lam) - 1)*mesh%p
        j_hi = j_lo + (size(mesh%zeta) - 1)*mesh%p
        n_jr = j_lo + j_hi
        allocate (dof(0:n_th - 1, 0:n_jr))
        n = 0
        do j = 0, n_jr
            do i = 0, n_th - 1
                if (i == 0 .and. j > j_lo .and. j <= j_hi) then
                    dof(i, j) = dof(0, j_lo)
                else
                    n = n + 1
                    dof(i, j) = n
                end if
            end do
        end do
    end subroutine numbering

    !> Real system of the axisymmetric monoenergetic problem, b.grad theta =
    !! kappa b; drives 1 (radial) and 2 (parallel).
    subroutine hp_assemble(f, mesh, nu, sys, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        type(hp_system_t), intent(out) :: sys
        integer, intent(out) :: status
        type(hp_operator_t) :: op
        complex(dp), allocatable :: mat(:, :)
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: e, i, j, n_loc, nt
        type(fortsparse_status_t) :: st

        op%k = f%kappa
        op%k_power = 1
        op%nu = nu
        call elements(f, mesh, op, sys, mat, status)
        if (status /= HP_OK) return
        status = HP_BAD_INPUT
        sys%mat = real(mat, dp)
        n_loc = (mesh%p + 1)**2
        nt = size(mat)
        allocate (rows(nt + 1), cols(nt + 1))
        do e = 1, size(mat, 2)
            do j = 1, n_loc
                do i = 1, n_loc
                    rows((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(i, e)
                    cols((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(j, e)
                end do
            end do
        end do
        vals = [reshape(sys%mat, [nt]), 0.0_dp]
        ! Constants span the null spaces of A and A^T; the equation of the
        ! first unknown, a combination of the others, pins it instead. A
        ! dense border row would defeat the fill-reducing ordering.
        where (rows(:nt) == 1) vals(:nt) = 0
        rows(nt + 1) = 1; cols(nt + 1) = 1; vals(nt + 1) = maxval(abs(vals(:nt)))
        call csc_from_triplet(sys%n, sys%n, rows, cols, vals, sys%a, st)
        if (.not. status_ok(st)) return
        status = HP_OK
    end subroutine hp_assemble

    !> Galerkin matrix of one element by tensor Gauss quadrature, column
    !> major in s, with its loads (radial and parallel drive, NTV drive and
    !> torque functional) and density weights.
    !> xi_n, xi_q: xi of the element's lower and upper row curve at its nodal
    !> and quadrature angles.
    subroutine element(f, mesh, op, xi_n, xi_q, ec, load, dens, s, f10, parts, f10_ops, &
            qxi, qw)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(in) :: xi_n(:, :), xi_q(:, :)
        integer, intent(in) :: ec
        complex(dp), intent(out) :: load(:, :), s(:)
        real(dp), intent(out) :: dens(:)
        real(dp), intent(in), optional :: f10(:)
        ! The terms of s apart (bracket, collide, phase, mass, precess) and
        ! the operators mir and ls that f10 enters the drive with.
        real(dp), intent(out), optional :: parts(:, :)
        complex(dp), intent(out), optional :: f10_ops(:, :)
        ! xi and weight (measure included) of the quadrature point (q1, q2).
        real(dp), intent(out), optional :: qxi(:, :), qw(:, :)
        integer :: p, nq, q1, q2, a1, a2, ia, ic, n_loc, j, lo, hi
        real(dp), allocatable :: node(:), wn(:), xq(:), wq(:), phi(:, :), dphi(:, :), &
            b_ab(:), d_ab(:), da(:), db_(:), psi(:, :), psi1(:)
        real(dp) :: th, th_a, xi_lo, xi_hi, xi, xi_b, &
            col, dr, de, den, ww, psi_a, psi_b, phase, precess, drive, torque
        complex(dp) :: loc_term, src, s_dpar, s_dlf, s_db
        real(dp) :: bb, bb_t, f10_xi

        p = mesh%p
        n_loc = (p + 1)**2
        nq = 2*p + 3
        allocate (node(p + 1), wn(p + 1), xq(nq), wq(nq), phi(p + 1, nq), &
            dphi(p + 1, nq), b_ab(n_loc), d_ab(n_loc), da(n_loc), db_(n_loc), &
            psi(p + 1, p + 1), psi1(n_loc))
        call gauss_lobatto_legendre(p + 1, node, wn)
        call gauss_legendre(nq, xq, wq)
        do q1 = 1, nq
            call lagrange_weights(p + 1, xq(q1), node, phi(:, q1))
            call lagrange_deriv_weights(p + 1, xq(q1), node, dphi(:, q1))
        end do
        s = 0; load = 0; dens = 0
        if (present(parts)) parts = 0
        if (present(f10_ops)) f10_ops = 0
        th_a = (mesh%theta(ec + 1) - mesh%theta(ec))/2
        ! Streaming as the bracket of the nodal interpolant of the stream
        ! function: the bracket is invariant under the element map and
        ! polynomial, so the quadrature is exact, constants are exact left
        ! and right null vectors and the form is skew to rounding.
        do a1 = 0, p
            th = mesh%theta(ec) + th_a*(1 + node(a1 + 1))
            xi_lo = xi_n(1, a1 + 1)
            xi_hi = xi_n(2, a1 + 1)
            do a2 = 0, p
                psi(a1 + 1, a2 + 1) = f%stream(th, xi_lo + (xi_hi - xi_lo)*(1 + node(a2 + 1))/2)
            end do
        end do
        psi1 = reshape(psi, [n_loc])
        do q1 = 1, nq
            th = mesh%theta(ec) + th_a*(1 + xq(q1))
            xi_lo = xi_q(1, q1)
            xi_hi = xi_q(2, q1)
            xi_b = (xi_hi - xi_lo)/2
            do q2 = 1, nq
                xi = xi_lo + xi_b*(1 + xq(q2))
                ww = wq(q1)*wq(q2)*th_a*xi_b
                call f%point_source(th, xi, op, src, col, dr, de, den, phase, precess, drive, &
                    torque)
                if (present(qxi)) qxi(q1, q2) = xi
                if (present(qw)) qw(q1, q2) = ww*den
                ! -i n q xi k + i (n/v)(omega_E + pre (1 + wsign xi^2)(cos - offset)/b)
                loc_term = (0, 1)*op%n*(-op%q*phase + (op%omega_e*den + op%pre*precess) &
                    /op%v)
                do a2 = 0, p
                    do a1 = 0, p
                        ia = 1 + a1 + (p + 1)*a2
                        b_ab(ia) = phi(a1 + 1, q1)*phi(a2 + 1, q2)
                        da(ia) = dphi(a1 + 1, q1)*phi(a2 + 1, q2)
                        db_(ia) = phi(a1 + 1, q1)*dphi(a2 + 1, q2)
                        d_ab(ia) = db_(ia)/xi_b
                    end do
                end do
                psi_a = sum(psi1*da)
                psi_b = sum(psi1*db_)
                do ic = 1, n_loc
                    s((ic - 1)*n_loc + 1:ic*n_loc) = s((ic - 1)*n_loc + 1:ic*n_loc) &
                        + wq(q1)*wq(q2)*b_ab*(psi_b*da(ic) - psi_a*db_(ic)) &
                        + ww*(col*d_ab*d_ab(ic) + loc_term*b_ab*b_ab(ic))
                end do
                if (present(parts)) then
                    do ic = 1, n_loc
                        lo = (ic - 1)*n_loc + 1
                        hi = ic*n_loc
                        parts(lo:hi, 1) = parts(lo:hi, 1) &
                            + wq(q1)*wq(q2)*b_ab*(psi_b*da(ic) - psi_a*db_(ic))
                        parts(lo:hi, 2) = parts(lo:hi, 2) + ww*col*d_ab*d_ab(ic)
                        parts(lo:hi, 3) = parts(lo:hi, 3) + ww*phase*b_ab*b_ab(ic)
                        parts(lo:hi, 4) = parts(lo:hi, 4) + ww*den*b_ab*b_ab(ic)
                        parts(lo:hi, 5) = parts(lo:hi, 5) + ww*precess*b_ab*b_ab(ic)
                    end do
                end if
                load(:, 1) = load(:, 1) + ww*dr*b_ab
                load(:, 2) = load(:, 2) + ww*de*b_ab
                load(:, 3) = load(:, 3) + ww*drive*src/op%v*b_ab
                if (op%source_dpar .or. present(f10) .or. &
                    (op%source_dlf10 .and. present(f10_ops))) then
                    call field_value(f, th, bb, bb_t)
                    call perturbed_sources(op, th, bb_t/bb, s_dpar, s_dlf, s_db)
                    if (op%source_dpar) load(:, 3) = load(:, 3) + ww*drive*s_dpar/op%v*b_ab
                    if (present(f10)) then
                        ! Mirror force of the perturbed parallel gradient on
                        ! f_10, (k/b) s_dlf (1 - xi^2)/2 d_xi f_10, and minus
                        ! the perturbed collision term, (dB/B)(nu/2) L f_10
                        ! with L = -d_xi (1 - xi^2) d_xi, by parts in xi.
                        f10_xi = sum(f10*d_ab)
                        load(:, 3) = load(:, 3) + ww*den*(op%k*bb**op%k_power/bb)*s_dlf &
                            *(1 - xi*xi)/2*f10_xi*b_ab &
                            + ww*den*(s_db/bb)*op%nu/2*(1 - xi*xi)*f10_xi*d_ab
                    end if
                    if (op%source_dlf10 .and. present(f10_ops)) then
                        ! The same two terms as operators on f_10, the second
                        ! per unit nu.
                        do ic = 1, n_loc
                            lo = (ic - 1)*n_loc + 1
                            hi = ic*n_loc
                            f10_ops(lo:hi, 1) = f10_ops(lo:hi, 1) + ww*den &
                                *(op%k*bb**op%k_power/bb)*s_dlf*(1 - xi*xi)/2*b_ab*d_ab(ic)
                            f10_ops(lo:hi, 2) = f10_ops(lo:hi, 2) + ww*den*(s_db/bb) &
                                *(1 - xi*xi)/2*d_ab*d_ab(ic)
                        end do
                    end if
                end if
                load(:, 4) = load(:, 4) + ww*torque*conjg(src)*b_ab
                dens = dens + ww*den*b_ab
            end do
        end do
    end subroutine element

    subroutine hp_factor(sys, status)
        type(hp_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        status = HP_SINGULAR
        sys%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        sys%solver%options%ordering = FORTSPARSE_ORDER_AMD
        sys%solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(sys%solver, sys%a, st)
        if (.not. status_ok(st)) return
        sys%factored = .true.
        status = HP_OK
    end subroutine hp_factor

    !> Solve A x + mu d = rhs, <d, x> = 0, or the same with the transpose of
    !! A (the adjoint); d is the density weight, mu the density source that
    !! makes the right-hand side compatible.
    subroutine hp_solve(sys, rhs, x, mu, status, adjoint)
        type(hp_system_t), intent(inout) :: sys
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:), mu
        integer, intent(out) :: status
        logical, intent(in), optional :: adjoint
        type(fortsparse_status_t) :: st
        integer :: op, it
        real(dp), allocatable :: r(:), dr(:), dx(:)
        op = FORTSPARSE_OP_NORMAL
        if (present(adjoint)) then
            if (adjoint) op = FORTSPARSE_OP_TRANSPOSE
        end if
        status = HP_SINGULAR
        mu = sum(rhs)/sum(sys%density)
        r = rhs - mu*sys%density
        r(1) = 0
        call sparse_solve_op(sys%solver, r, x, op, st)
        if (.not. status_ok(st)) return
        ! The supernodal factors use static pivots; refine on the true residual.
        allocate (dx(size(x)))
        do it = 1, 3
            dr = r - apply(sys%a, x, op == FORTSPARSE_OP_TRANSPOSE)
            call sparse_solve_op(sys%solver, dr, dx, op, st)
            if (.not. status_ok(st)) return
            x = x + dx
        end do
        x = x - dot_product(sys%density, x)/sys%norm
        status = HP_OK
    end subroutine hp_solve

    !> J = <s_j g_i> (D_ji) on the mesh and the dual-weighted residual
    !! estimate eta of J - J_h: the residual of g_h in the space of degree
    !! p + 1 weighted by the adjoint z+ - I z_h of that space, where
    !! A^T z = s_j. eta_el holds element indicators for refinement.
    !! This is an indicator, not a bound.
    subroutine hp_estimate(f, mesh, nu, i, j, value, eta, eta_el, status, n_dof)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        integer, intent(in) :: i, j
        real(dp), intent(out) :: value, eta
        real(dp), allocatable, intent(out) :: eta_el(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: n_dof
        type(hp_system_t) :: sys, fine
        type(hp_mesh_t) :: mesh_fine
        real(dp), allocatable :: g(:), z(:), zf(:), gf(:), izf(:), node(:), &
            node_f(:), w(:), lag(:, :), r(:)
        real(dp) :: mu, mu_z, mu_f
        integer :: p, e, a, n_loc

        call hp_assemble(f, mesh, nu, sys, status)
        if (status /= HP_OK) return
        call hp_factor(sys, status)
        if (status /= HP_OK) return
        if (present(n_dof)) n_dof = sys%n
        allocate (g(sys%n), z(sys%n))
        call hp_solve(sys, real(sys%drive(:, i), dp), g, mu, status)
        if (status == HP_OK) call hp_solve(sys, real(sys%drive(:, j), dp), z, mu_z, status, &
            adjoint=.true.)
        call hp_free(sys)
        if (status /= HP_OK) return
        value = dot_product(real(sys%drive(:, j), dp), g)/sys%norm

        mesh_fine = mesh
        mesh_fine%p = mesh%p + 1
        call hp_assemble(f, mesh_fine, nu, fine, status)
        if (status /= HP_OK) return
        call hp_factor(fine, status)
        if (status /= HP_OK) return
        allocate (zf(fine%n))
        call hp_solve(fine, real(fine%drive(:, j), dp), zf, mu_f, status, adjoint=.true.)
        call hp_free(fine)
        if (status /= HP_OK) return

        ! Interpolation of degree p into degree p + 1 per element.
        p = mesh%p
        allocate (node(p + 1), node_f(p + 2), w(p + 2), lag(p + 2, p + 1))
        call gauss_lobatto_legendre(p + 1, node, w(:p + 1))
        call gauss_lobatto_legendre(p + 2, node_f, w)
        do a = 1, p + 2
            call lagrange_weights(p + 1, node_f(a), node, lag(a, :))
        end do
        n_loc = (p + 2)**2
        allocate (gf(fine%n), izf(fine%n), eta_el(size(fine%loc, 2)), r(n_loc))
        eta = 0
        do e = 1, size(fine%loc, 2)
            gf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(g(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
            izf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(z(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
        end do
        do e = 1, size(fine%loc, 2)
            r = real(fine%load(:, i, e), dp) - mu*fine%dens(:, e) &
                - matmul(reshape(fine%mat(:, e), [n_loc, n_loc]), gf(fine%loc(:, e)))
            eta = eta + dot_product(zf(fine%loc(:, e)) - izf(fine%loc(:, e)), r)/fine%norm
            ! Marking uses node-wise magnitudes; signed element sums cancel.
            eta_el(e) = sum(abs((zf(fine%loc(:, e)) - izf(fine%loc(:, e)))*r))/fine%norm
        end do
    end subroutine hp_estimate

    !> Bisect the rows and columns that carry the largest indicator sums
    !! until they hold the fraction of the total (Doerfler marking). Rows
    !! are bisected in their curve parameter, mirrored in sigma.
    subroutine hp_refine(mesh, eta_el, fraction, out)
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: eta_el(:), fraction
        type(hp_mesh_t), intent(out) :: out
        integer :: n_col, n_lev, n_z, n_row, e, er, ec, m
        real(dp), allocatable :: s_col(:), s_lam(:), s_z(:), all(:)
        real(dp) :: total, acc, cut

        n_col = size(mesh%theta) - 1
        n_lev = size(mesh%lam); n_z = size(mesh%zeta)
        n_row = 2*(n_lev - 1) + (n_z - 1)
        allocate (s_col(n_col), s_lam(n_lev - 1), s_z(n_z - 1))
        s_col = 0; s_lam = 0; s_z = 0
        do e = 1, size(eta_el)
            ec = 1 + mod(e - 1, n_col)
            er = 1 + (e - 1)/n_col
            s_col(ec) = s_col(ec) + abs(eta_el(e))
            if (er <= n_lev - 1) then
                s_lam(er) = s_lam(er) + abs(eta_el(e))
            else if (er <= n_lev - 1 + n_z - 1) then
                s_z(er - (n_lev - 1)) = s_z(er - (n_lev - 1)) + abs(eta_el(e))
            else
                s_lam(n_row + 1 - er) = s_lam(n_row + 1 - er) + abs(eta_el(e))
            end if
        end do
        ! Trapped rows mirror in xi.
        s_z = s_z + s_z(n_z - 1:1:-1)
        ! Smallest threshold whose larger entries hold the fraction.
        all = [s_col, s_lam, s_z]
        total = sum(all)
        call sort_desc(all)
        acc = 0; cut = all(1)
        do m = 1, size(all)
            acc = acc + all(m)
            cut = all(m)
            if (acc >= fraction*total) exit
        end do
        out%p = mesh%p
        out%theta = bisect(mesh%theta, s_col >= cut)
        out%lam = bisect(mesh%lam, s_lam >= cut)
        out%zeta = bisect(mesh%zeta, s_z >= cut)
    end subroutine hp_refine

    pure function bisect(x, mark) result(y)
        real(dp), intent(in) :: x(:)
        logical, intent(in) :: mark(:)
        real(dp), allocatable :: y(:)
        integer :: k
        y = [x(1)]
        do k = 1, size(mark)
            if (mark(k)) y = [y, (x(k) + x(k + 1))/2]
            y = [y, x(k + 1)]
        end do
    end function bisect

    pure subroutine sort_desc(x)
        real(dp), intent(inout) :: x(:)
        integer :: k, m
        real(dp) :: t
        do k = 2, size(x)
            t = x(k)
            m = k - 1
            do while (m >= 1)
                if (x(m) >= t) exit
                x(m + 1) = x(m)
                m = m - 1
            end do
            x(m + 1) = t
        end do
    end subroutine sort_desc

    !> Flux kernel R = int dtheta dxi (2/b)(1 + wsign xi^2)/b Re(S^* F) of one
    !! toroidal harmonic (the kernel of milca_dk_ntv), op%k the constant
    !! b.grad theta (k_power = 0) or k b (k_power = 1). Without a null space
    !! the complex system is solved directly. With edges (increasing flux
    !! labels), dep(b) is the part of R whose quadrature points lie at f%position
    !! in bin b (below edges(1), between edges, above the last); sum(dep) = R.
    subroutine hp_ntv_kernel(f, mesh, op, r, status, n_dof, edges, dep)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: r
        integer, intent(out) :: status
        integer, intent(out), optional :: n_dof
        real(dp), intent(in), optional :: edges(:)
        real(dp), intent(out), optional :: dep(:)
        type(hp_system_t) :: sys
        type(hp_ntv_factor_t) :: fac
        complex(dp), allocatable :: mat(:, :), x(:)

        call ntv_system(f, mesh, op, sys, mat, status)
        if (status /= HP_OK) return
        if (present(n_dof)) n_dof = sys%n
        call hp_ntv_factor(sys, mesh%p, mat, fac, status)
        if (status /= HP_OK) return
        allocate (x(sys%n))
        call hp_ntv_solve(fac, sys%drive(:, 3), x, .false., status)
        call hp_ntv_free(fac)
        if (status /= HP_OK) return
        r = real(sum(sys%drive(:, 4)*x), dp)
        if (present(edges) .and. present(dep)) call ntv_deposit(f, mesh, op, sys, x, edges, dep)
    end subroutine hp_ntv_kernel

    !> The torque functional of hp_ntv_kernel split by where its quadrature
    !! points lie: the same tensor Gauss rule as the load of drive 4 in element,
    !! so that the bins sum to Re sum(drive(:, 4) x) up to rounding. A finite
    !! orbit spends time in proportion to the measure at each of its points, so
    !! this is the torque deposited where the orbit is, not on its psi*.
    subroutine ntv_deposit(f, mesh, op, sys, x, edges, dep)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x(:)
        real(dp), intent(in) :: edges(:)
        real(dp), intent(out) :: dep(:)
        integer :: p, nq, n_col, n_lev, n_row, ec, er, e, c, q1, q2, a1, a2
        integer, allocatable :: curve_kind(:)
        real(dp), allocatable :: curve_val(:), node(:), wn(:), xq(:), wq(:), phi(:, :), &
            xi_q(:, :)
        real(dp) :: th_a, th, d, xi_lo, xi_hi, xi_b, xi, ww, col, dr, de, den, phase, &
            precess, drive, torque, s
        complex(dp) :: src, xh
        complex(dp), allocatable :: xl(:, :)

        dep = 0
        p = mesh%p
        nq = 2*p + 3
        n_col = size(mesh%theta) - 1
        n_lev = size(mesh%lam)
        curve_val = [mesh%lam(1:n_lev - 1), mesh%zeta, mesh%lam(n_lev - 1:1:-1)]
        curve_kind = [(-1, c=1, n_lev - 1), (0, c=1, size(mesh%zeta)), (1, c=1, n_lev - 1)]
        n_row = size(curve_val) - 1
        allocate (node(p + 1), wn(p + 1), xq(nq), wq(nq), phi(p + 1, nq), &
            xi_q(n_row + 1, nq), xl(p + 1, p + 1))
        call gauss_lobatto_legendre(p + 1, node, wn)
        call gauss_legendre(nq, xq, wq)
        do q1 = 1, nq
            call lagrange_weights(p + 1, xq(q1), node, phi(:, q1))
        end do
        do ec = 1, n_col
            th_a = (mesh%theta(ec + 1) - mesh%theta(ec))/2
            do c = 1, n_row + 1
                do q1 = 1, nq
                    th = mesh%theta(ec) + th_a*(1 + xq(q1))
                    call f%level(curve_kind(c), curve_val(c), th, xi_q(c, q1), d)
                end do
            end do
            do er = 1, n_row
                e = ec + (er - 1)*n_col
                xl = reshape(x(sys%loc(:, e)), [p + 1, p + 1])
                do q1 = 1, nq
                    th = mesh%theta(ec) + th_a*(1 + xq(q1))
                    xi_lo = xi_q(er, q1)
                    xi_hi = xi_q(er + 1, q1)
                    xi_b = (xi_hi - xi_lo)/2
                    do q2 = 1, nq
                        xi = xi_lo + xi_b*(1 + xq(q2))
                        ww = wq(q1)*wq(q2)*th_a*xi_b
                        call f%point_source(th, xi, op, src, col, dr, de, den, phase, &
                            precess, drive, torque)
                        xh = 0
                        do a2 = 1, p + 1
                            do a1 = 1, p + 1
                                xh = xh + phi(a1, q1)*phi(a2, q2)*xl(a1, a2)
                            end do
                        end do
                        s = f%position(th, xi)
                        c = 1 + count(edges <= s)
                        dep(c) = dep(c) + real(ww*torque*conjg(src)*xh, dp)
                    end do
                end do
            end do
        end do
    end subroutine ntv_deposit

    !> hp_ntv_kernel with the dual-weighted residual estimate eta of the
    !! error of R: the primal and adjoint solutions of the mesh are
    !! interpolated into degree p + 1, where the adjoint is solved again, and
    !! eta = Re sum_e (z_f - I z)^T (b_f - A_f I x). R + eta is the
    !! corrected value; eta_el(e) = sum |(z_f - I z) (b_f - A_f I x)| over the
    !! nodes of element e marks it for refinement (hp_refine).
    subroutine hp_ntv_estimate(f, mesh, op, r, eta, eta_el, status, n_dof)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: r, eta
        real(dp), allocatable, intent(out) :: eta_el(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: n_dof
        type(hp_system_t) :: sys, fine
        type(hp_mesh_t) :: mesh_fine
        type(hp_ntv_factor_t) :: fac
        complex(dp), allocatable :: mat(:, :), matf(:, :), x(:), z(:), zf(:), xf(:), izf(:), &
            res(:)
        real(dp), allocatable :: node(:), node_f(:), w(:), lag(:, :)
        complex(dp) :: acc
        integer :: p, e, a, n_loc

        call ntv_system(f, mesh, op, sys, mat, status)
        if (status /= HP_OK) return
        if (present(n_dof)) n_dof = sys%n
        call hp_ntv_factor(sys, mesh%p, mat, fac, status)
        if (status /= HP_OK) return
        allocate (x(sys%n), z(sys%n))
        call hp_ntv_solve(fac, sys%drive(:, 3), x, .false., status)
        if (status == HP_OK) call hp_ntv_solve(fac, sys%drive(:, 4), z, .true., status)
        call hp_ntv_free(fac)
        if (status /= HP_OK) return
        r = real(sum(sys%drive(:, 4)*x), dp)

        mesh_fine = mesh
        mesh_fine%p = mesh%p + 1
        call ntv_system(f, mesh_fine, op, fine, matf, status)
        if (status /= HP_OK) return
        call hp_ntv_factor(fine, mesh_fine%p, matf, fac, status)
        if (status /= HP_OK) return
        allocate (zf(fine%n))
        call hp_ntv_solve(fac, fine%drive(:, 4), zf, .true., status)
        call hp_ntv_free(fac)
        if (status /= HP_OK) return

        p = mesh%p
        allocate (node(p + 1), node_f(p + 2), w(p + 2), lag(p + 2, p + 1))
        call gauss_lobatto_legendre(p + 1, node, w(:p + 1))
        call gauss_lobatto_legendre(p + 2, node_f, w)
        do a = 1, p + 2
            call lagrange_weights(p + 1, node_f(a), node, lag(a, :))
        end do
        n_loc = (p + 2)**2
        allocate (xf(fine%n), izf(fine%n), eta_el(size(fine%loc, 2)), res(n_loc))
        do e = 1, size(fine%loc, 2)
            xf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(x(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
            izf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(z(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
        end do
        acc = 0
        do e = 1, size(fine%loc, 2)
            res = fine%load(:, 3, e) - matmul(reshape(matf(:, e), [n_loc, n_loc]), &
                xf(fine%loc(:, e)))
            acc = acc + sum((zf(fine%loc(:, e)) - izf(fine%loc(:, e)))*res)
            eta_el(e) = sum(abs((zf(fine%loc(:, e)) - izf(fine%loc(:, e)))*res))
        end do
        eta = real(acc, dp)
    end subroutine hp_ntv_estimate

    !> The NTV system of hp_ntv_kernel on mesh: element matrices mat, loads
    !! (drive 3) and output weights (drive 4); with source_dlf10, f_10 of the
    !! axisymmetric problem on the same mesh enters the loads.
    subroutine ntv_system(f, mesh, op, sys, mat, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_system_t), intent(out) :: sys
        complex(dp), allocatable, intent(out) :: mat(:, :)
        integer, intent(out) :: status

        status = HP_BAD_INPUT
        if (op%n == 0 .or. .not. allocated(op%m)) return
        if (op%source_dlf10) then
            ! f_10 of the axisymmetric problem on this mesh, in the kernel's
            ! normalisation per unit speed (milca_dk_ntv ntv_f10).
            block
                type(hp_system_t) :: ax
                real(dp), allocatable :: f10(:)
                real(dp) :: mu
                select type (f)
                class is (hp_field_t)
                class default
                    return
                end select
                if (op%q == 0 .or. f%kappa == 0) return
                call hp_assemble(f, mesh, op%nu, ax, status)
                if (status /= HP_OK) return
                call hp_factor(ax, status)
                if (status /= HP_OK) return
                allocate (f10(ax%n))
                call hp_solve(ax, real(ax%drive(:, 1), dp), f10, mu, status)
                call hp_free(ax)
                if (status /= HP_OK) return
                f10 = -2/(op%v*op%q*f%kappa)*f10
                call elements(f, mesh, op, sys, mat, status, f10)
            end block
        else
            call elements(f, mesh, op, sys, mat, status)
        end if
    end subroutine ntv_system

    !> The terms of the NTV system of hp_ntv_kernel on mesh, apart (hp_terms_t);
    !! nu and v of op are not used.
    subroutine hp_ntv_terms(f, mesh, op, terms, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_terms_t), intent(out) :: terms
        integer, intent(out) :: status
        type(hp_operator_t) :: unit
        type(hp_system_t) :: sys
        complex(dp), allocatable :: mat(:, :)

        status = HP_BAD_INPUT
        if (op%n == 0 .or. .not. allocated(op%m)) return
        if (op%source_dlf10) then
            ! The f_10 operators need b of a thin field, as in ntv_system.
            select type (f)
            class is (hp_field_t)
            class default
                return
            end select
        end if
        ! collide and the drive are linear in nu and 1/v.
        unit = op
        unit%nu = 1
        unit%v = 1
        call elements(f, mesh, unit, sys, mat, status, terms=terms)
        if (status /= HP_OK) return
        terms%n = sys%n
        terms%p = mesh%p
        terms%norm = sys%norm
        call move_alloc(sys%loc, terms%loc)
        call move_alloc(sys%density, terms%density)
        terms%load1 = real(sys%drive(:, 1), dp)
        terms%load3 = sys%drive(:, 3)
        terms%load4 = sys%drive(:, 4)
    end subroutine hp_ntv_terms

    !> The matrix S + nu C + i n (-q Ph + (omega_e M + pre Pr)/v) of the NTV
    !! system at collisionality nu = nu_D/v and speed v from the terms of its
    !! mesh; a is left empty when terms are not those of hp_ntv_terms.
    subroutine hp_node_matrix(terms, nu, v, omega_e, pre, n, q, a)
        type(hp_terms_t), intent(in) :: terms
        real(dp), intent(in) :: nu, v, omega_e, pre, q
        integer, intent(in) :: n
        type(csc_z_t), intent(out) :: a
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        integer :: nt

        if (.not. allocated(terms%loc)) return
        nt = size(terms%s)
        call pattern(terms%loc, rows, cols)
        vals = reshape(terms%s + nu*terms%c + (0, 1)*n*(-q*terms%ph &
            + (omega_e*terms%m + pre*terms%pr)/v), [nt])
        call csc_from_triplet(terms%n, terms%n, rows, cols, vals, a, st)
        if (.not. status_ok(st)) a = csc_z_t()
    end subroutine hp_node_matrix

    !> The mass matrix int measure phi_i phi_j of the mesh of terms; m is
    !! left empty when terms are not those of hp_ntv_terms.
    subroutine hp_mass(terms, m)
        type(hp_terms_t), intent(in) :: terms
        type(csc_t), intent(out) :: m
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:)

        if (.not. allocated(terms%loc)) return
        call pattern(terms%loc, rows, cols)
        call csc_from_triplet(terms%n, terms%n, rows, cols, &
            reshape(terms%m, [size(terms%m)]), m, st)
        if (.not. status_ok(st)) m = csc_t()
    end subroutine hp_mass

    !> Adds the breakpoints of m to mesh.
    pure subroutine hp_mesh_join(mesh, m)
        type(hp_mesh_t), intent(inout) :: mesh
        type(hp_mesh_t), intent(in) :: m
        mesh%theta = merged(mesh%theta, m%theta)
        mesh%lam = merged(mesh%lam, m%lam)
        mesh%zeta = merged(mesh%zeta, m%zeta)
    end subroutine hp_mesh_join

    !> Sorted union of two increasing breakpoint lists; points closer than
    !> 1e-12 count as one.
    pure function merged(u, v) result(w)
        real(dp), intent(in) :: u(:), v(:)
        real(dp), allocatable :: w(:)
        real(dp) :: tmp(size(u) + size(v)), x
        integer :: i, j, n
        i = 1
        j = 1
        n = 0
        do while (i <= size(u) .or. j <= size(v))
            if (j > size(v)) then
                x = u(i)
                i = i + 1
            else if (i > size(u)) then
                x = v(j)
                j = j + 1
            else if (u(i) <= v(j)) then
                x = u(i)
                i = i + 1
            else
                x = v(j)
                j = j + 1
            end if
            if (n > 0) then
                if (abs(x - tmp(n)) <= 1e-12_dp) cycle
            end if
            n = n + 1
            tmp(n) = x
        end do
        w = tmp(1:n)
    end function merged

    !> Element elem and reference coordinates (r, s) in [-1, 1] of the points
    !! (th(t), xi(j, t)) on mesh of geometry f, numbered j + (t - 1) size(xi, 1):
    !! theta taken into the period of the mesh, xi placed linearly between
    !! the row curves of f at that angle, as the element map does. The mesh
    !! interpolant at a point is then sum phi_a(r) phi_b(s) over the nodes of
    !! its element, exact on the row curves whatever their width, and points
    !! of one angle share the row curves. The points are physical: for the
    !! mesh of another energy at finite orbit width they lie elsewhere
    !! relative to its separatrix.
    subroutine hp_locate(f, mesh, th, xi, elem, r, s, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: th(:), xi(:, :)
        integer, allocatable, intent(out) :: elem(:)
        real(dp), allocatable, intent(out) :: r(:), s(:)
        integer, intent(out) :: status
        integer, allocatable :: kind(:)
        real(dp), allocatable :: val(:), xc(:)
        real(dp) :: t, rr, d
        integer :: n_col, n_row, nx, it, j, ec, er, c, k

        status = HP_BAD_INPUT
        if (mesh%p < 1 .or. size(xi, 2) /= size(th)) return
        call row_curves(mesh, kind, val)
        n_col = size(mesh%theta) - 1
        n_row = size(val) - 1
        nx = size(xi, 1)
        allocate (elem(nx*size(th)), r(nx*size(th)), s(nx*size(th)))
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(it, t, ec, rr, xc, c, d, j, er, k)
        do it = 1, size(th)
            allocate (xc(n_row + 1))
            t = in_period(mesh, th(it))
            call interval(mesh%theta, t, ec, rr)
            do c = 1, n_row + 1
                call f%level(kind(c), val(c), t, xc(c), d)
            end do
            do j = 1, nx
                k = j + (it - 1)*nx
                call interval(xc, xi(j, it), er, s(k))
                elem(k) = ec + (er - 1)*n_col
                r(k) = rr
            end do
            deallocate (xc)
        end do
        !$omp end parallel do
        status = HP_OK
    end subroutine hp_locate

    !> theta in the period [theta(1), theta(end)] of the mesh columns.
    pure real(dp) function in_period(mesh, th) result(t)
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: th
        t = mesh%theta(1) + modulo(th - mesh%theta(1), 2*pi)
        t = min(t, mesh%theta(size(mesh%theta)))
    end function in_period

    !> The interval k of the increasing points x that holds x0 and the
    !> reference coordinate r in [-1, 1] of x0 in it; -1 in an interval of no
    !> width (the rows collapsed on the X point).
    pure subroutine interval(x, x0, k, r)
        real(dp), intent(in) :: x(:), x0
        integer, intent(out) :: k
        real(dp), intent(out) :: r
        real(dp) :: d
        k = max(1, min(count(x(2:) <= x0) + 1, size(x) - 1))
        d = x(k + 1) - x(k)
        r = -1
        if (d > 0) r = max(-1.0_dp, min(1.0_dp, 2*(x0 - x(k))/d - 1))
    end subroutine interval

    !> 2 p + 3 Gauss angles, p the highest degree, on every column of the
    !! union of the theta breakpoints of the meshes of one surface, taken into
    !! the period of the first mesh: at finite orbit width the meshes of other
    !! energies start at their own X point.
    subroutine hp_lines(meshes, lines)
        type(hp_mesh_t), intent(in) :: meshes(:)
        type(hp_lines_t), intent(out) :: lines
        real(dp), allocatable :: theta(:), xq(:), wq(:), v(:)
        integer :: k, nq, ec

        lines%theta = [real(dp) ::]
        lines%w = [real(dp) ::]
        if (size(meshes) < 1) return
        theta = meshes(1)%theta
        do k = 2, size(meshes)
            v = [(in_period(meshes(1), meshes(k)%theta(ec)), ec=1, size(meshes(k)%theta))]
            v = -v
            call sort_desc(v)
            theta = merged(theta, -v)
        end do
        nq = 2*maxval(meshes%p) + 3
        allocate (xq(nq), wq(nq))
        call gauss_legendre(nq, xq, wq)
        lines%theta = [((theta(ec) + (theta(ec + 1) - theta(ec))*(1 + xq(k))/2, k=1, nq), &
            ec=1, size(theta) - 1)]
        lines%w = [(((theta(ec + 1) - theta(ec))/2*wq(k), k=1, nq), ec=1, size(theta) - 1)]
    end subroutine hp_lines

    !> Pitch moments on lines and their loads on mesh, l = 0..l_max:
    !!   g(l) (lines x unknowns)  G[t, i] = (2 l + 1)/2 int P_l(xi) phi_i dxi,
    !!   d(l) (unknowns x lines)  D[i, t] = w_t int measure P_l(xi) phi_i dxi
    !! at theta_t, with xi the pitch of this geometry. The rule is exact for
    !! the moments, whose integrand is a polynomial of degree p + l between
    !! the row curves. With measure and measure_mesh (the speed node of
    !! another energy at finite orbit width), d integrates that geometry's
    !! measure over the rows of measure_mesh at the same (theta, xi), 2 p + 3
    !! Gauss points per row, phi_i evaluated there through the element map of
    !! mesh, and P_l(xi) of l >= 1 less its mean in that measure along the
    !! line: 1^T D does not depend on mesh and vanishes for l >= 1, so the
    !! field part conserves the density of the orbit measure on any set of
    !! meshes and at any width; for thin orbits the mean is zero.
    subroutine hp_line_ops(f, mesh, op, lines, l_max, g, d, status, measure, measure_mesh)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_lines_t), intent(in) :: lines
        integer, intent(in) :: l_max
        type(csc_t), allocatable, intent(out) :: g(:), d(:)
        integer, intent(out) :: status
        class(hp_geometry_t), intent(in), optional :: measure
        type(hp_mesh_t), intent(in), optional :: measure_mesh
        type(fortsparse_status_t) :: st
        type(hp_operator_t) :: o
        integer, allocatable :: dof(:, :), kind(:), col(:, :), row(:, :), kind_m(:)
        real(dp), allocatable :: val(:), node(:), wn(:), xq(:), wq(:), phi(:, :), pa(:), &
            xc(:), gj(:, :), dj(:, :), gv(:, :, :), dv(:, :, :), pl(:), val_m(:), xm(:), &
            xqm(:), wqm(:), ps(:), pts(:, :), cl(:)
        real(dp) :: th, r, xi, xi_b, dxi, collide, dr, de, den, phase, precess, drive, &
            torque, sr
        complex(dp) :: src
        integer :: p, n, n_col, n_row, n_th, n_jr, n_t, nq, nqm, n_line, t, ec, er, c, q, a1, &
            a2, j, k, l, ek, np_m
        logical :: other

        status = HP_BAD_INPUT
        p = mesh%p
        n_t = size(lines%theta)
        if (p < 1 .or. l_max < 0 .or. n_t < 1 .or. size(lines%w) /= n_t) return
        other = present(measure)
        if (other .neqv. present(measure_mesh)) return
        n_col = size(mesh%theta) - 1
        ! The measure does not depend on the perturbation.
        o = op
        if (.not. allocated(o%m)) then
            o%m = [0]
            o%c = [(0.0_dp, 0.0_dp)]
        end if
        call row_curves(mesh, kind, val)
        call numbering(mesh, dof, n)
        n_row = size(val) - 1
        n_th = n_col*p
        n_jr = n_row*p
        nq = (p + l_max + 2)/2
        n_line = (p + 1)*(n_jr + 1)
        allocate (node(p + 1), wn(p + 1), xq(nq), wq(nq), phi(p + 1, nq), &
            col(n_line, n_t), row(n_line, n_t), gv(n_line, 0:l_max, n_t), &
            dv(n_line, 0:l_max, n_t))
        call gauss_lobatto_legendre(p + 1, node, wn)
        call gauss_legendre(nq, xq, wq)
        do q = 1, nq
            call lagrange_weights(p + 1, xq(q), node, phi(:, q))
        end do
        nqm = 0
        if (other) then
            call row_curves(measure_mesh, kind_m, val_m)
            nqm = 2*measure_mesh%p + 3
            allocate (xqm(nqm), wqm(nqm))
            call gauss_legendre(nqm, xqm, wqm)
        end if
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(t, th, ec, r, pa, xc, xm, gj, dj, pl, ps, c, er, q, a1, a2, j, k, l, xi, &
        !$omp xi_b, dxi, src, collide, dr, de, den, phase, precess, drive, torque, ek, sr, &
        !$omp pts, cl, np_m)
        do t = 1, n_t
            allocate (pa(p + 1), xc(n_row + 1), gj(0:l_max, 0:n_jr), dj(0:l_max, 0:n_jr), &
                pl(0:l_max), ps(p + 1))
            th = in_period(mesh, lines%theta(t))
            ec = min(count(mesh%theta(2:) <= th) + 1, n_col)
            r = 2*(th - mesh%theta(ec))/(mesh%theta(ec + 1) - mesh%theta(ec)) - 1
            call lagrange_weights(p + 1, r, node, pa)
            do c = 1, n_row + 1
                call f%level(kind(c), val(c), th, xc(c), dxi)
            end do
            gj = 0
            dj = 0
            do er = 1, n_row
                xi_b = (xc(er + 1) - xc(er))/2
                do q = 1, nq
                    xi = xc(er) + xi_b*(1 + xq(q))
                    pl = [(legendre_p(l, 0, xi), l=0, l_max)]
                    if (.not. other) call f%point_source(th, xi, o, src, collide, dr, de, &
                        den, phase, precess, drive, torque)
                    do a2 = 0, p
                        j = (er - 1)*p + a2
                        gj(:, j) = gj(:, j) + wq(q)*xi_b*phi(a2 + 1, q)*pl
                        if (.not. other) dj(:, j) = dj(:, j) + wq(q)*xi_b*phi(a2 + 1, q)*den*pl
                    end do
                end do
            end do
            if (other) then
                allocate (xm(size(val_m)), pts(3, (size(val_m) - 1)*nqm), cl(0:l_max))
                do c = 1, size(val_m)
                    call measure%level(kind_m(c), val_m(c), th, xm(c), dxi)
                end do
                ! Points (xi, weight times measure) on the rows of measure_mesh.
                np_m = 0
                do er = 1, size(val_m) - 1
                    xi_b = (xm(er + 1) - xm(er))/2
                    if (xi_b <= 0) cycle
                    do q = 1, nqm
                        xi = xm(er) + xi_b*(1 + xqm(q))
                        call measure%point_source(th, xi, o, src, collide, dr, de, den, &
                            phase, precess, drive, torque)
                        np_m = np_m + 1
                        pts(:, np_m) = [xi, wqm(q)*xi_b*den, 0.0_dp]
                    end do
                end do
                ! For l >= 1, P_l less its mean in the measure along the line:
                ! the measure of the orbits depends on xi, and int m P_l dxi
                ! would make the field part change the density at finite
                ! width. Thin orbits have a constant measure and mean zero.
                cl = 0
                do q = 1, np_m
                    cl = cl + pts(2, q)*[(legendre_p(l, 0, pts(1, q)), l=0, l_max)]
                end do
                cl(1:) = cl(1:)/cl(0)
                cl(0) = 0
                do q = 1, np_m
                    pl = [(legendre_p(l, 0, pts(1, q)), l=0, l_max)] - cl
                    call interval(xc, pts(1, q), ek, sr)
                    call lagrange_weights(p + 1, sr, node, ps)
                    do a2 = 0, p
                        j = (ek - 1)*p + a2
                        dj(:, j) = dj(:, j) + pts(2, q)*ps(a2 + 1)*pl
                    end do
                end do
                deallocate (xm, pts, cl)
            end if
            do j = 0, n_jr
                do a1 = 0, p
                    k = 1 + a1 + (p + 1)*j
                    row(k, t) = t
                    col(k, t) = dof(mod((ec - 1)*p + a1, n_th), j)
                    gv(k, :, t) = [((2*l + 1)/2.0_dp, l=0, l_max)]*pa(a1 + 1)*gj(:, j)
                    dv(k, :, t) = lines%w(t)*pa(a1 + 1)*dj(:, j)
                end do
            end do
            deallocate (pa, xc, gj, dj, pl, ps)
        end do
        !$omp end parallel do
        allocate (g(0:l_max), d(0:l_max))
        do l = 0, l_max
            call csc_from_triplet(n_t, n, reshape(row, [n_line*n_t]), &
                reshape(col, [n_line*n_t]), reshape(gv(:, l, :), [n_line*n_t]), g(l), st)
            if (.not. status_ok(st)) return
            call csc_from_triplet(n, n_t, reshape(col, [n_line*n_t]), &
                reshape(row, [n_line*n_t]), reshape(dv(:, l, :), [n_line*n_t]), d(l), st)
            if (.not. status_ok(st)) return
        end do
        status = HP_OK
    end subroutine hp_line_ops

    !> Equilibrated rows and columns keep the static-pivot growth of the
    !> supernodal LU in bounds on graded meshes.
    subroutine factor_elements(sys, p, mat, fac, status)
        type(hp_system_t), intent(in) :: sys
        integer, intent(in) :: p
        complex(dp), intent(in) :: mat(:, :)
        type(hp_ntv_factor_t), intent(out) :: fac
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        type(csc_z_t) :: a
        integer, allocatable :: rows(:), cols(:)

        status = HP_BAD_INPUT
        if (size(sys%loc, 1) /= (p + 1)**2) return
        call pattern(sys%loc, rows, cols)
        call csc_from_triplet(sys%n, sys%n, rows, cols, reshape(mat, [size(mat)]), a, st)
        if (.not. status_ok(st)) return
        call factor_matrix(a, fac, status)
    end subroutine factor_elements

    !> The same factor of an assembled square matrix (hp_node_matrix).
    subroutine factor_matrix(a, fac, status)
        type(csc_z_t), intent(in) :: a
        type(hp_ntv_factor_t), intent(out) :: fac
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        type(csc_z_t) :: scaled

        status = HP_BAD_INPUT
        if (a%nrow < 1 .or. a%nrow /= a%ncol) return
        status = HP_SINGULAR
        fac%a = a
        scaled = a
        call equilibrate(scaled, fac%rs, fac%cs)
        fac%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        fac%solver%options%ordering = FORTSPARSE_ORDER_AMD
        fac%solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(fac%solver, scaled, st)
        if (.not. status_ok(st)) return
        status = HP_OK
    end subroutine factor_matrix

    !> Rows and columns of the entries of element matrices stored column
    !> major, element after element.
    pure subroutine pattern(loc, rows, cols)
        integer, intent(in) :: loc(:, :)
        integer, allocatable, intent(out) :: rows(:), cols(:)
        integer :: e, i, j, n_loc
        n_loc = size(loc, 1)
        allocate (rows(n_loc**2*size(loc, 2)), cols(n_loc**2*size(loc, 2)))
        do e = 1, size(loc, 2)
            do j = 1, n_loc
                do i = 1, n_loc
                    rows((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = loc(i, e)
                    cols((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = loc(j, e)
                end do
            end do
        end do
    end subroutine pattern

    !> A x = b, or A^T x = b with transpose, by the equilibrated factor
    !> S = R A C and refine (default three) refinement steps against A.
    subroutine hp_ntv_solve(fac, b, x, transpose, status, refine)
        type(hp_ntv_factor_t), intent(inout) :: fac
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        logical, intent(in) :: transpose
        integer, intent(out) :: status
        integer, intent(in), optional :: refine
        type(fortsparse_status_t) :: st
        complex(dp), allocatable :: res(:), dx(:)
        integer :: it, op, n_refine
        real(dp), allocatable :: left(:), right(:)

        status = HP_SINGULAR
        n_refine = 3
        if (present(refine)) n_refine = refine
        if (transpose) then
            ! A^T = C^-1 S^T R^-1: x = R S^-T C b.
            left = fac%cs
            right = fac%rs
            op = FORTSPARSE_OP_TRANSPOSE
        else
            left = fac%rs
            right = fac%cs
            op = FORTSPARSE_OP_NORMAL
        end if
        allocate (dx(size(b)))
        res = left*b
        call sparse_solve_op(fac%solver, res, x, op, st)
        x = right*x
        do it = 1, n_refine
            if (.not. status_ok(st)) exit
            if (transpose) then
                res = left*(b - apply_zt(fac%a, x))
            else
                res = left*(b - apply_z(fac%a, x))
            end if
            call sparse_solve_op(fac%solver, res, dx, op, st)
            x = x + right*dx
        end do
        if (status_ok(st)) status = HP_OK
    end subroutine hp_ntv_solve

    subroutine hp_ntv_free(fac)
        type(hp_ntv_factor_t), intent(inout) :: fac
        call sparse_destroy(fac%solver)
    end subroutine hp_ntv_free

    !> Row and column scalings rs, cs with rs A cs of unit row and column
    !> maxima, by three sweeps of Ruiz's iteration; a is replaced.
    subroutine equilibrate(a, rs, cs)
        type(csc_z_t), intent(inout) :: a
        real(dp), allocatable, intent(out) :: rs(:), cs(:)
        real(dp), allocatable :: r(:), c(:)
        integer :: sweep, j, k
        allocate (rs(a%nrow), cs(a%ncol), r(a%nrow), c(a%ncol))
        rs = 1; cs = 1
        do sweep = 1, 3
            r = 0; c = 0
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                    r(a%row_idx(k)) = max(r(a%row_idx(k)), abs(a%val(k)))
                    c(j) = max(c(j), abs(a%val(k)))
                end do
            end do
            where (r > 0) r = 1/sqrt(r)
            where (c > 0) c = 1/sqrt(c)
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                    a%val(k) = a%val(k)*r(a%row_idx(k))*c(j)
                end do
            end do
            rs = rs*r
            cs = cs*c
        end do
    end subroutine equilibrate

    pure function apply_zt(a, x) result(y)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(a%ncol)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(j) = y(j) + a%val(k)*x(a%row_idx(k))
            end do
        end do
    end function apply_zt

    pure function apply_z(a, x) result(y)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
            end do
        end do
    end function apply_z

    !> A x or A^T x of a CSC matrix.
    pure function apply(a, x, transpose) result(y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        logical, intent(in) :: transpose
        real(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (transpose) then
                    y(j) = y(j) + a%val(k)*x(a%row_idx(k))
                else
                    y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
                end if
            end do
        end do
    end function apply

    subroutine hp_free(sys)
        type(hp_system_t), intent(inout) :: sys
        if (sys%factored) call sparse_destroy(sys%solver)
        sys%factored = .false.
    end subroutine hp_free

    !> D_ij = <s_i g_j> on the given mesh; border holds the density-source
    !! multipliers of the two drives.
    subroutine hp_coefficients(f, mesh, nu, c, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        type(mono_coefficients_t), intent(out) :: c
        integer, intent(out) :: status
        type(hp_system_t) :: sys
        real(dp), allocatable :: g1(:), g3(:)
        integer :: n

        call hp_assemble(f, mesh, nu, sys, status)
        if (status /= HP_OK) return
        call hp_factor(sys, status)
        if (status /= HP_OK) return
        n = sys%n
        allocate (g1(n), g3(n))
        call hp_solve(sys, real(sys%drive(:, 1), dp), g1, c%border(1), status)
        if (status == HP_OK) call hp_solve(sys, real(sys%drive(:, 2), dp), g3, c%border(2), status)
        call hp_free(sys)
        if (status /= HP_OK) return
        c%d11 = dot_product(real(sys%drive(:, 1), dp), g1)/sys%norm
        c%d31 = dot_product(real(sys%drive(:, 2), dp), g1)/sys%norm
        c%d13 = dot_product(real(sys%drive(:, 1), dp), g3)/sys%norm
        c%d33 = dot_product(real(sys%drive(:, 2), dp), g3)/sys%norm
    end subroutine hp_coefficients

    !> xi and d xi/d theta on a row curve: a passing level of sign kind
    !! (kind = +-1, val = Lambda) or a scaled separatrix (kind = 0, val =
    !! zeta). Differences of b are formed without cancellation near the X
    !! point.
    subroutine curve(g, kind, val, th, xi, dxi)
        class(hp_field_t), intent(in) :: g
        integer, intent(in) :: kind
        real(dp), intent(in) :: val, th
        real(dp), intent(out) :: xi, dxi
        real(dp) :: v, q, db
        v = separatrix2(g, th)
        db = field_db(g, th)
        if (kind /= 0) then
            q = g%b_max*((1/g%b_max - val) + val*v)
            xi = kind*sqrt(q)
            dxi = -kind*val*db/(2*sqrt(q))
        else if (v > 0) then
            xi = val*sqrt(v)
            dxi = -val*db/(2*g%b_max*sqrt(v))
        else
            xi = 0; dxi = 0
        end if
    end subroutine curve

    !> The drive of one toroidal harmonic at (theta, xi), i n sum_m c_m
    !> exp(i m theta) from the operator's amplitudes; a geometry whose orbits
    !> leave the surface overrides it with the amplitudes at the orbit point.
    !> b and db/dtheta of a thin field (zero slope otherwise).
    subroutine field_value(f, th, b, b_t)
        class(hp_geometry_t), intent(in) :: f
        real(dp), intent(in) :: th
        real(dp), intent(out) :: b, b_t
        select type (f)
        class is (hp_field_t)
            b = field_b(f, th)
            b_t = field_db(f, th)
        class default
            b = 1
            b_t = 0
        end select
    end subroutine field_value

    !> At theta with d_theta ln b = slope: the parallel-derivative source
    !> sum_m (i (m - n q) - 2 slope)/q c_m e^{i m theta}, the f_10 mirror
    !> factor sum_m (i (m - n q) - slope) c_m e^{i m theta} and dB/B.
    pure subroutine perturbed_sources(op, th, slope, s_dpar, s_dlf, s_db)
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(in) :: th, slope
        complex(dp), intent(out) :: s_dpar, s_dlf, s_db
        complex(dp) :: e, z(size(op%m))
        integer :: j
        s_dpar = 0; s_dlf = 0; s_db = 0
        z = phases(op%m, th)
        do j = 1, size(op%m)
            e = op%c(j)*z(j)
            s_dpar = s_dpar + ((0, 1)*(op%m(j) - op%n*op%q) - 2*slope)/op%q*e
            s_dlf = s_dlf + ((0, 1)*(op%m(j) - op%n*op%q) - slope)*e
            s_db = s_db + e
        end do
    end subroutine perturbed_sources

    subroutine field_point_source(g, th, xi, op, src, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque)
        class(hp_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(out) :: src
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        src = g%perturbation(th, xi, op)
        call g%point(th, xi, op, collide, drive_r, drive_e, density, phase, precess, drive, &
            torque)
    end subroutine field_point_source

    real(dp) function field_position(g, th, xi) result(s)
        class(hp_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        s = -1
    end function field_position

    complex(dp) function op_perturbation(g, th, xi, op) result(src)
        class(hp_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp) :: z(size(op%m))
        integer :: j
        src = 0
        if (.not. allocated(op%m)) return
        z = phases(op%m, th)
        do j = 1, size(op%m)
            src = src + (0, 1)*op%n*op%c(j)*z(j)
        end do
    end function op_perturbation

    !> exp(i m th) by angle addition from one sincos (rounding grows like |m|
    !> epsilon).
    pure function phases(m, th) result(z)
        integer, intent(in) :: m(:)
        real(dp), intent(in) :: th
        complex(dp) :: z(size(m)), w(0:max(0, maxval(abs(m)))), e
        integer :: k
        e = cmplx(cos(th), sin(th), dp)
        w(0) = 1
        do k = 1, ubound(w, 1)
            w(k) = w(k - 1)*e
        end do
        do k = 1, size(m)
            z(k) = w(abs(m(k)))
            if (m(k) < 0) z(k) = conjg(z(k))
        end do
    end function phases

    real(dp) function field_stream(g, th, xi) result(psi)
        class(hp_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        call dk_hp_stream(field_b(g, th), xi, psi)
    end function field_stream

    subroutine field_point(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
            precess, drive, torque)
        class(hp_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp) :: bb, k
        bb = field_b(g, th)
        k = op%k*bb**op%k_power
        call dk_hp_point(bb, field_db(g, th), xi, k, op%nu, collide, drive_r, drive_e, density)
        call dk_hp_ntv(bb, xi, k, cos(th), op%offset, op%wsign, phase, precess, drive, torque)
        ! The flux is weighted with the surface measure, which streaming
        ! b.grad theta = k b**k_power divides by b**k_power (Boozer: 2/b^3).
        torque = torque/bb**op%k_power
    end subroutine field_point

    !> xi_c^2 = 1 - b/b_max.
    pure real(dp) function separatrix2(f, th) result(s)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        real(dp) :: d
        d = modulo(th - f%theta_x + pi, 2*pi) - pi
        s = max(0.0_dp, -field_diff(f, f%theta_x, d)/f%b_max)
    end function separatrix2

    pure real(dp) function field_b(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        complex(dp) :: e, z
        integer :: k
        b = f%a(0)
        e = cmplx(cos(th), sin(th), dp)
        z = 1
        do k = 1, f%m
            z = z*e
            b = b + f%a(k)*real(z, dp) + f%c(k)*aimag(z)
        end do
    end function field_b

    pure real(dp) function field_db(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        complex(dp) :: e, z
        integer :: k
        b = 0
        e = cmplx(cos(th), sin(th), dp)
        z = 1
        do k = 1, f%m
            z = z*e
            b = b + k*(f%c(k)*real(z, dp) - f%a(k)*aimag(z))
        end do
    end function field_db

    pure real(dp) function field_d2b(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        complex(dp) :: e, z
        integer :: k
        b = 0
        e = cmplx(cos(th), sin(th), dp)
        z = 1
        do k = 1, f%m
            z = z*e
            b = b - k*k*(f%a(k)*real(z, dp) + f%c(k)*aimag(z))
        end do
    end function field_d2b

    !> b(y + d) - b(y) from products of sines, exact to rounding for small d.
    pure real(dp) function field_diff(f, y, d) result(s)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: y, d
        integer :: k
        s = 0
        do k = 1, f%m
            s = s + 2*sin(k*d/2)*(f%c(k)*cos(k*(y + d/2)) - f%a(k)*sin(k*(y + d/2)))
        end do
    end function field_diff

    real(dp) function newton_extremum(f, th0) result(th)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th0
        integer :: it
        real(dp) :: step
        th = th0
        do it = 1, 50
            step = field_db(f, th)/field_d2b(f, th)
            th = th - step
            if (abs(step) < 1e-15_dp) exit
        end do
        th = modulo(th, 2*pi)
    end function newton_extremum

end module milca_dk_hp
