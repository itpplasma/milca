module milca_ideal_fb_gamma_p
    !! Vacuum of the 3D free boundary as a boundary integral on the moving
    !! plasma boundary Gamma_p, with no vacuum mesh.
    !! Gamma_p is the chain of the plasma map: nedge edges of degree order
    !! (Lagrange in t, shared end nodes, nchain = nedge*order nodes) whose
    !! chart positions (q1, q2) carry real Fourier modes in phi-hat (n = k nfp,
    !! mode 1 constant, 2k cos, 2k + 1 sin), x = H(q1, q2, phi-hat).
    !! n points out of the plasma. With the coil field B_e at fixed currents
    !! and the circulation c (mu0 I_p) carried by a fixed closed filament
    !! inside the plasma (field B_f, potential A_f per unit circulation), the
    !! vacuum field is B = B_x + grad Phi, B_x = B_e + c B_f, with Phi single
    !! valued, decaying and B.n = 0. Up to a constant, per field period,
    !!   mu0 E_v = -1/2 int_p |B_e|^2 dV + c X + c^2/2 G - 1/2 h^T S^-1 h,
    !!   X = -int n.(A_f x B_e) dA,  G = -int n.(A_f x B_f) dA,
    !!   h_i = int B_x.n N_i dA,
    !! S the Galerkin exterior Dirichlet-to-Neumann map of Gamma_p (FortFEM
    !! P2 BEM). X and G follow from int_ext B.curl A_f = surface terms plus
    !! int_ext A_f.J_coils, which is constant while the coils stay outside.
    !! The volume term is the surface integral of F = (1/R) int_Rref^R
    !! |B_e|^2 R' dR' e_R, so E_v depends on Gamma_p only. Its shape gradient
    !! is the Hadamard form dE_v = -1/(2 mu0) int |B|^2 V_n dA with B the
    !! discrete field on Gamma_p (gamma_p_traction): consistent with the
    !! discrete E_v to the order of the BEM, not to roundoff.
    !! Two exteriors: GAMMA_P_AXISYMMETRIC (cylindrical chart, nmode = 1,
    !! the n = 0 toroidal-mode BEM of FortFEM on P2 panels through the chain)
    !! and GAMMA_P_CIRCULANT (any nfp-symmetric chart, parametric P2 BEM
    !! assembled block-circulant over the nfp field periods; only the zero
    !! sector harmonic acts on the periodic Phi).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t, bem_operators_t, &
        bem_exterior_dtn_t, build_exterior_dtn
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t, assemble_parametric_bem_p2_circulant
    use fortfem_bem_circulant, only: bem_circulant_t
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_exterior_harmonic_fields, only: evaluate_current_loops
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_panel_point
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric
    use fortnum_quadrature, only: gauss_legendre
    use milca_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t
    use milca_fb_gamma_p_leaves, only: gamma_p_field, gamma_p_traction, &
        gamma_p_traction_jvp, gamma_p_flux, gamma_p_edge_p2, gamma_p_edge_p3
    implicit none
    private
    public :: gamma_p_t, gamma_p_field_t, gamma_p_build, gamma_p_vacuum
    public :: gamma_p_surface_t, gamma_p_edge_point, gamma_p_hessian
    public :: gamma_p_edge_basis
    integer, parameter, public :: GAMMA_P_AXISYMMETRIC = 1, GAMMA_P_CIRCULANT = 2
    integer, parameter, public :: GAMMA_P_OK = 0, GAMMA_P_BAD_INPUT = 1, &
        GAMMA_P_BEM_FAILED = 2, GAMMA_P_FIELD_FAILED = 3, GAMMA_P_SOLVE_FAILED = 4
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! External (coil) field at Cartesian points, Cartesian components.
    ! axisymmetric: (B_R, B_phi, B_Z) of the n = 0 harmonic at meridional
    ! points (R, Z); the default evaluates points at phi = 0, exact for
    ! axisymmetric fields only.
    type, abstract :: gamma_p_field_t
    contains
        procedure(field_points), deferred :: points
        procedure :: axisymmetric => field_axisymmetric
    end type
    abstract interface
        subroutine field_points(self, x, b, status)
            import :: dp, gamma_p_field_t
            class(gamma_p_field_t), intent(in) :: self
            real(dp), intent(in) :: x(:, :)
            real(dp), intent(out) :: b(:, :)
            integer, intent(out) :: status
        end subroutine
    end interface

    ! Parametric Gamma_p for the 3D BEM: panel p lies on chain edge
    ! panel_edge(p) with (t, phi-hat) corners parameters(:, :, p).
    type, extends(bem_surface_t) :: gamma_p_surface_t
        class(ideal_chart_t), allocatable :: chart
        integer :: order = 0, nfp = 1, nmode = 1
        real(dp), allocatable :: coef(:, :, :, :), parameters(:, :, :)
        integer, allocatable :: panel_edge(:)
    contains
        procedure :: evaluate => surface_evaluate
    end type

    type :: gamma_p_t
        integer :: kind = 0, nfp = 1, nmode = 1, order = 0, nedge = 0, nchain = 0
        ! Quadrature over one field period: nq points with Cartesian x, the
        ! tangents a, b (a x b = n dA per du dv, out of the plasma), weight w
        ! (du dv and the period factor), chart angle phi-hat, nb trace basis
        ! values and (u, v) derivatives on trace dofs, the chain edge, the
        ! weights of its order + 1 nodes in the displacement with their
        ! (u, v) derivatives, the chart columns dx/dq1, dx/dq2, the chart
        ! point q with dq/du, dq/dv (qd) and (dphi/du, dphi/dv) of the modes.
        ! The axisymmetric exterior takes u = phi, v = the panel parameter.
        integer :: nq = 0, nb = 0, ntrace = 0, nline = 16
        real(dp), allocatable :: x(:, :), a(:, :), b(:, :), w(:), angle(:)
        real(dp), allocatable :: basis(:, :), basis_u(:, :), basis_v(:, :)
        integer, allocatable :: dof(:, :), edge(:)
        real(dp), allocatable :: shape(:, :), shape_u(:, :), shape_v(:, :)
        real(dp), allocatable :: d1(:, :, :), q(:, :), qd(:, :, :), angle_uv(:, :)
        class(ideal_chart_t), allocatable :: chart
        ! Filament polygon and its potential and field per unit circulation.
        real(dp), allocatable :: filament(:, :), af(:, :), bf(:, :)
        ! Exterior DtN of one period on the trace dofs and its Cholesky
        ! factor; Cholesky factor of the trace mass matrix (Hessian).
        real(dp), allocatable :: s(:, :), factor(:, :), mass_factor(:, :)
        real(dp) :: r_ref = 0
        type(bem_pair_report_t) :: report
        ! Last gamma_p_vacuum: trace potential Phi, field B at the points,
        ! load h and the energy parts (bem, volume, filament; T^2 m^3).
        ! B_x = B_e + c B_f at the points and its Cartesian Jacobian
        ! (gamma_p_hessian, central differences of the field).
        real(dp), allocatable :: phi(:), field(:, :), load(:), bx(:, :), dbx(:, :, :)
        real(dp) :: energy_parts(3) = 0, circulation = 0, mu0 = 1
    end type
contains
    subroutine field_axisymmetric(self, rz, b, status)
        class(gamma_p_field_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        real(dp), intent(out) :: b(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: x(:, :)
        allocate (x(3, size(rz, 2)))
        x(1, :) = rz(1, :); x(2, :) = 0; x(3, :) = rz(2, :)
        call self%points(x, b, status)
    end subroutine

    pure subroutine gamma_p_edge_basis(order, t, w, dw)
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        call edge_basis(order, t, w, dw)
    end subroutine

    pure subroutine edge_basis(order, t, w, dw)
        ! Lagrange basis of the chain edge on its order + 1 equispaced nodes,
        ! the edge restriction of the plasma map's triangle basis (checked
        ! against triangle_basis_jet by test_ideal_fb_gamma_p edge-basis).
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        if (order == 3) then
            call gamma_p_edge_p3(t, w(0), w(1), w(2), w(3), dw(0), dw(1), dw(2), dw(3))
        else
            call gamma_p_edge_p2(t, w(0), w(1), w(2), dw(0), dw(1), dw(2))
        end if
    end subroutine

    pure subroutine trig_modes(nmode, nfp, phi, trig, dtrig)
        ! One sincos per point; harmonics by the angle-addition recurrence.
        integer, intent(in) :: nmode, nfp
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: trig(nmode), dtrig(nmode)
        real(dp) :: c1, s1, c, s, cn
        integer :: m, k
        trig(1) = 1; dtrig(1) = 0
        c1 = cos(nfp*phi); s1 = sin(nfp*phi); c = 1; s = 0
        do m = 2, nmode, 2
            cn = c*c1 - s*s1; s = s*c1 + c*s1; c = cn
            k = (m/2)*nfp
            trig(m) = c; dtrig(m) = -k*s
            if (m + 1 <= nmode) then
                trig(m + 1) = s; dtrig(m + 1) = k*c
            end if
        end do
    end subroutine

    pure subroutine gamma_p_edge_point(coef, order, nfp, t, phi, q, qt, qp)
        !! Chart point (q1, q2) of edge coefficients coef(2, nmode, 0:order) at
        !! (t, phi-hat) and its t and phi-hat derivatives.
        real(dp), intent(in) :: coef(:, :, 0:)
        integer, intent(in) :: order, nfp
        real(dp), intent(in) :: t, phi
        real(dp), intent(out) :: q(2), qt(2), qp(2)
        real(dp) :: w(0:3), dw(0:3), trig(size(coef, 2)), dtrig(size(coef, 2))
        integer :: l
        call edge_basis(order, t, w(:order), dw(:order))
        call trig_modes(size(coef, 2), nfp, phi, trig, dtrig)
        q = 0; qt = 0; qp = 0
        do l = 0, order
            q = q + w(l)*matmul(coef(:, :, l), trig)
            qt = qt + dw(l)*matmul(coef(:, :, l), trig)
            qp = qp + w(l)*matmul(coef(:, :, l), dtrig)
        end do
    end subroutine

    pure subroutine surface_evaluate(self, panel, uv, x, a, b, status)
        class(gamma_p_surface_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        real(dp) :: p(2), deriv(2, 2), q(3), qt(2), qp(2), qd(3, 2)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), jac, normal(3)
        integer :: k
        status = GAMMA_P_BAD_INPUT; x = 0; a = 0; b = 0
        if (panel < 1 .or. panel > size(self%panel_edge)) return
        if (.not. all(ieee_is_finite(uv))) return
        deriv(:, 1) = self%parameters(:, 2, panel) - self%parameters(:, 1, panel)
        deriv(:, 2) = self%parameters(:, 3, panel) - self%parameters(:, 1, panel)
        p = self%parameters(:, 1, panel) + uv(1)*deriv(:, 1) + uv(2)*deriv(:, 2)
        call gamma_p_edge_point(self%coef(:, :, :, self%panel_edge(panel)), &
            self%order, self%nfp, p(1), p(2), q(1:2), qt, qp)
        q(3) = p(2)
        do k = 1, 2
            qd(1:2, k) = qt*deriv(1, k) + qp*deriv(2, k)
            qd(3, k) = deriv(2, k)
        end do
        call self%chart%jet(q, x, d1, d2, d3)
        a = matmul(d1, qd(:, 1)); b = matmul(d1, qd(:, 2))
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(a)) .or. &
            .not. all(ieee_is_finite(b))) return
        call bem_frame(a(1), a(2), a(3), b(1), b(2), b(3), jac, &
                       normal(1), normal(2), normal(3))
        if (.not. ieee_is_finite(jac) .or. jac <= 0) return
        status = GAMMA_P_OK
    end subroutine

    subroutine gamma_p_build(gp, chart, nfp, order, chain, kind, status, &
                             panels_per_edge, nphi, degree, options, filament, r_ref)
        !! chain(2, nmode, nchain): chart positions of the chain nodes per
        !! real Fourier mode, nchain = nedge*order, node (e - 1) order + l + 1
        !! is local node l of edge e (the last one wraps to node 1).
        !! AXISYMMETRIC: panels_per_edge P2 panels per edge through the chain,
        !! degree Gauss points per panel (default 8). CIRCULANT: nphi strips
        !! per full turn (multiple of nfp), panels_per_edge subdivisions in t,
        !! triangle rule of polynomial degree `degree` (default 8) and the
        !! FortFEM pair options. filament(3, nv): closed polygon inside the
        !! plasma carrying the circulation (CIRCULANT: nfp symmetric).
        !! r_ref: start of the R lines of the volume term (default: mean R).
        type(gamma_p_t), intent(out) :: gp
        class(ideal_chart_t), intent(in) :: chart
        integer, intent(in) :: nfp, order, kind
        real(dp), intent(in) :: chain(:, :, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: panels_per_edge, nphi, degree
        type(bem_pair_options_t), intent(in), optional :: options
        real(dp), intent(in), optional :: filament(:, :)
        real(dp), intent(in), optional :: r_ref
        integer :: ppe, nd
        status = GAMMA_P_BAD_INPUT
        ppe = 2
        if (present(panels_per_edge)) ppe = panels_per_edge
        nd = 8
        if (present(degree)) nd = degree
        if (ppe < 1 .or. nd < 2 .or. nfp < 1) return
        if (order /= 2 .and. order /= 3) return
        if (size(chain, 1) /= 2 .or. mod(size(chain, 3), order) /= 0) return
        if (mod(size(chain, 2), 2) /= 1) return
        if (.not. all(ieee_is_finite(chain))) return
        gp%kind = kind; gp%nfp = nfp; gp%order = order
        gp%nmode = size(chain, 2); gp%nchain = size(chain, 3)
        gp%nedge = gp%nchain/order
        if (gp%nedge < 2) return
        select case (kind)
        case (GAMMA_P_AXISYMMETRIC)
            if (gp%nmode /= 1) return
            select type (chart)
            type is (cylindrical_chart_t)
            class default
                return
            end select
            call build_axisymmetric(gp, chain, ppe, nd, status)
            allocate (gp%chart, source=chart)
        case (GAMMA_P_CIRCULANT)
            if (.not. present(nphi)) return
            if (nphi < 4 .or. mod(nphi, nfp) /= 0) return
            if (nphi <= 2*(gp%nmode/2)*nfp) return
            call build_circulant(gp, chart, chain, ppe, nphi, nd, options, status)
            allocate (gp%chart, source=chart)
        case default
            return
        end select
        if (status /= GAMMA_P_OK) return
        status = GAMMA_P_BAD_INPUT
        gp%r_ref = sum(sqrt(gp%x(1, :)**2 + gp%x(2, :)**2)*gp%w)/sum(gp%w)
        if (present(r_ref)) gp%r_ref = r_ref
        allocate (gp%af(3, gp%nq), gp%bf(3, gp%nq))
        gp%af = 0; gp%bf = 0
        if (present(filament)) then
            if (size(filament, 1) /= 3 .or. size(filament, 2) < 3) return
            gp%filament = filament
            call filament_fields(gp, filament, status)
            if (status /= GAMMA_P_OK) return
        end if
        call factor_dtn(gp, status)
        if (status /= GAMMA_P_OK) return
        call factor_mass(gp, status)
    end subroutine

    subroutine filament_fields(gp, filament, status)
        type(gamma_p_t), intent(inout) :: gp
        real(dp), intent(in) :: filament(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: potential(:, :, :), field(:, :, :)
        integer :: st
        status = GAMMA_P_FIELD_FAILED
        allocate (potential(3, gp%nq, 1), field(3, gp%nq, 1))
        call evaluate_current_loops(filament, [1, size(filament, 2) + 1], gp%x, &
                                    potential, field, st)
        if (st /= 0) return
        gp%af = potential(:, :, 1); gp%bf = field(:, :, 1)
        status = GAMMA_P_OK
    end subroutine

    subroutine factor_dtn(gp, status)
        type(gamma_p_t), intent(inout) :: gp
        integer, intent(out) :: status
        integer :: info
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        status = GAMMA_P_SOLVE_FAILED
        gp%factor = gp%s
        call dpotrf('L', gp%ntrace, gp%factor, gp%ntrace, info)
        if (info /= 0) return
        status = GAMMA_P_OK
    end subroutine

    subroutine build_axisymmetric(gp, chain, ppe, nd, status)
        ! P2 panels through the chain at t = (j - 1)/ppe, (j - 1/2)/ppe of each
        ! edge, counter-clockwise in (R, Z). Points: a = dx/dphi, b = dx/dxi at
        ! phi = 0, so a x b points out; w carries 2 pi/nfp.
        type(gamma_p_t), intent(inout) :: gp
        real(dp), intent(in) :: chain(:, :, :)
        integer, intent(in) :: ppe, nd
        integer, intent(out) :: status
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: nodes(:, :), tnode(:, :), xg(:), wg(:), s(:, :)
        integer, allocatable :: pedge(:)
        real(dp) :: area, q(2), qt(2), qp(2), w(0:3), dw(0:3), psi(3), dpsi(3)
        real(dp) :: xr(2), tau(2), speed, t
        integer :: npan, p, j, e, k, g, iq, st, l
        status = GAMMA_P_BAD_INPUT
        npan = gp%nedge*ppe
        allocate (nodes(2, 2*npan), tnode(2, npan), pedge(npan))
        do e = 1, gp%nedge
            do j = 1, ppe
                p = (e - 1)*ppe + j
                pedge(p) = e
                tnode(:, p) = [real(j - 1, dp), real(j, dp)]/ppe
                do k = 1, 2
                    t = (real(j - 1, dp) + 0.5_dp*(k - 1))/ppe
                    call gamma_p_edge_point(edge_coef(chain, gp%order, e), gp%order, &
                        gp%nfp, t, 0.0_dp, q, qt, qp)
                    nodes(:, 2*p - 2 + k) = q
                end do
            end do
        end do
        area = 0
        do k = 1, 2*npan
            j = mod(k, 2*npan) + 1
            area = area + nodes(1, k)*nodes(2, j) - nodes(1, j)*nodes(2, k)
        end do
        if (area < 0) then
            ! Reverse: panel p runs backwards along edge pedge.
            nodes = nodes(:, [(mod(2*npan - k + 1, 2*npan) + 1, k=1, 2*npan)])
            pedge = pedge(npan:1:-1)
            tnode = tnode([2, 1], npan:1:-1)
        end if
        call gs_bem_curve_from_nodes(nodes, curve, st)
        if (st /= 0) return
        status = GAMMA_P_BEM_FAILED
        call tm_bem_assemble(curve, 0, 2, ops, st)
        if (st /= 0) return
        call tm_bem_steklov_symmetric(ops, s, st)
        if (st /= 0) return
        gp%ntrace = 2*npan
        gp%s = (2*pi/gp%nfp)*0.5_dp*(s + transpose(s))
        allocate (xg(nd), wg(nd))
        call gauss_legendre(nd, xg, wg)
        xg = 0.5_dp*(xg + 1); wg = 0.5_dp*wg
        gp%nq = npan*nd; gp%nb = 3
        call allocate_points(gp)
        do p = 1, npan
            do g = 1, nd
                iq = (p - 1)*nd + g
                call gs_bem_panel_point(curve, p, xg(g), xr, tau, speed)
                psi = [(1 - xg(g))*(1 - 2*xg(g)), 4*xg(g)*(1 - xg(g)), &
                       xg(g)*(2*xg(g) - 1)]
                dpsi = [4*xg(g) - 3, 4 - 8*xg(g), 4*xg(g) - 1]
                gp%x(:, iq) = [xr(1), 0.0_dp, xr(2)]
                gp%a(:, iq) = [0.0_dp, xr(1), 0.0_dp]
                gp%b(:, iq) = speed*[tau(1), 0.0_dp, tau(2)]
                gp%w(iq) = wg(g)*2*pi/gp%nfp
                gp%angle(iq) = 0
                gp%basis(:, iq) = psi; gp%basis_u(:, iq) = 0
                gp%basis_v(:, iq) = dpsi
                gp%dof(:, iq) = [(mod(2*p - 3 + k, 2*npan) + 1, k=1, 3)]
                gp%edge(iq) = pedge(p)
                gp%shape(:, iq) = 0; gp%shape_u(:, iq) = 0; gp%shape_v(:, iq) = 0
                do k = 1, 3
                    t = tnode(1, p) + 0.5_dp*(k - 1)*(tnode(2, p) - tnode(1, p))
                    call edge_basis(gp%order, t, w(:gp%order), dw(:gp%order))
                    do l = 0, gp%order
                        gp%shape(l, iq) = gp%shape(l, iq) + psi(k)*w(l)
                        gp%shape_v(l, iq) = gp%shape_v(l, iq) + dpsi(k)*w(l)
                    end do
                end do
                gp%d1(:, 1, iq) = [1.0_dp, 0.0_dp, 0.0_dp]
                gp%d1(:, 2, iq) = [0.0_dp, 0.0_dp, 1.0_dp]
                gp%q(:, iq) = [xr(1), xr(2), 0.0_dp]
                gp%qd(:, 1, iq) = [0.0_dp, 0.0_dp, 1.0_dp]
                gp%qd(:, 2, iq) = speed*[tau(1), tau(2), 0.0_dp]
                gp%angle_uv(:, iq) = [1.0_dp, 0.0_dp]
            end do
        end do
        status = GAMMA_P_OK
    end subroutine

    pure function edge_coef(chain, order, e) result(coef)
        real(dp), intent(in) :: chain(:, :, :)
        integer, intent(in) :: order, e
        real(dp) :: coef(2, size(chain, 2), 0:order)
        integer :: l
        do l = 0, order
            coef(:, :, l) = chain(:, :, chain_node(size(chain, 3), order, e, l))
        end do
    end function

    pure integer function chain_node(nchain, order, e, l)
        integer, intent(in) :: nchain, order, e, l
        chain_node = mod((e - 1)*order + l, nchain) + 1
    end function

    subroutine allocate_points(gp)
        type(gamma_p_t), intent(inout) :: gp
        allocate (gp%x(3, gp%nq), gp%a(3, gp%nq), gp%b(3, gp%nq), gp%w(gp%nq), &
                  gp%angle(gp%nq), gp%basis(gp%nb, gp%nq), gp%basis_u(gp%nb, gp%nq), &
                  gp%basis_v(gp%nb, gp%nq), gp%dof(gp%nb, gp%nq), gp%edge(gp%nq), &
                  gp%shape(0:gp%order, gp%nq), gp%shape_u(0:gp%order, gp%nq), &
                  gp%shape_v(0:gp%order, gp%nq), gp%d1(3, 2, gp%nq), gp%q(3, gp%nq), &
                  gp%qd(3, 2, gp%nq), gp%angle_uv(2, gp%nq))
    end subroutine

    subroutine build_circulant(gp, chart, chain, refine, nphi, nd, options, status)
        type(gamma_p_t), intent(inout) :: gp
        class(ideal_chart_t), intent(in) :: chart
        real(dp), intent(in) :: chain(:, :, :)
        integer, intent(in) :: refine, nphi, nd
        type(bem_pair_options_t), intent(in), optional :: options
        integer, intent(out) :: status
        type(gamma_p_surface_t) :: surface
        type(bem_circulant_t) :: cop
        type(bem_operators_t) :: op0
        type(bem_exterior_dtn_t) :: exterior
        type(bem_pair_options_t) :: pair
        integer, allocatable :: faces(:, :), td(:, :), fd(:, :), rep(:)
        real(dp), allocatable :: xi(:), eta(:), qw(:)
        real(dp) :: volume, x(3), a(3), b(3), q(3), qt(2), qp(2), uv(2), p(2)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), w(0:3), dw(0:3)
        real(dp) :: hess(2, 2), g(2), deriv(2, 2)
        integer :: np, nf, i, k, t, iq, st, r, orientation
        status = GAMMA_P_BAD_INPUT
        pair = bem_pair_options_t()
        if (present(options)) pair = options
        allocate (surface%chart, source=chart)
        surface%order = gp%order; surface%nfp = gp%nfp; surface%nmode = gp%nmode
        allocate (surface%coef(2, gp%nmode, 0:gp%order, gp%nedge))
        do i = 1, gp%nedge
            surface%coef(:, :, :, i) = edge_coef(chain, gp%order, i)
        end do
        call triangle_duffy_quadrature(nd, xi, eta, qw, st)
        if (st /= 0) return
        orientation = 1
        call panels(surface, gp%nedge, refine, nphi, orientation, faces)
        ! (1/3) int x.n dA is the enclosed volume; negative means inward.
        volume = 0
        do i = 1, size(faces, 2)
            do t = 1, size(qw)
                call surface%evaluate(i, [xi(t), eta(t)], x, a, b, st)
                if (st /= 0) return
                volume = volume + qw(t)*dot_product(x, cross(a, b))/3
            end do
        end do
        if (volume < 0) then
            orientation = -1
            call panels(surface, gp%nedge, refine, nphi, orientation, faces)
        end if
        status = GAMMA_P_BEM_FAILED
        call assemble_parametric_bem_p2_circulant(surface, faces, gp%nfp, pair, &
            cop, td, fd, gp%report, st)
        if (st /= 0) return
        op0%v = sum(cop%v%blocks, dim=3)
        op0%k = sum(cop%k%blocks, dim=3)
        op0%w = sum(cop%w%blocks, dim=3)
        op0%mass = sum(cop%mass%blocks, dim=3)
        op0%genus = cop%genus
        call build_exterior_dtn(op0, exterior, st)
        if (st /= 0) return
        call exterior%matrix(gp%s, st)
        if (st /= 0) return
        gp%ntrace = size(gp%s, 1)
        allocate (rep(maxval(td)))
        rep = 0
        do r = 1, size(cop%w%row_dofs, 2)
            rep(cop%w%row_dofs(:, r)) = [(i, i=1, size(cop%w%row_dofs, 1))]
        end do
        if (any(rep < 1)) return
        status = GAMMA_P_BAD_INPUT
        nf = size(faces, 2); np = nf/gp%nfp
        gp%nq = np*size(qw); gp%nb = 6
        call allocate_points(gp)
        do i = 1, np
            do t = 1, size(qw)
                iq = (i - 1)*size(qw) + t
                uv = [xi(t), eta(t)]
                call surface%evaluate(i, uv, x, a, b, st)
                if (st /= 0) return
                gp%x(:, iq) = x; gp%a(:, iq) = a; gp%b(:, iq) = b; gp%w(iq) = qw(t)
                do k = 1, 6
                    call generated_p2_jet(k, uv(1), uv(2), gp%basis(k, iq), g, hess)
                    gp%basis_u(k, iq) = g(1); gp%basis_v(k, iq) = g(2)
                    gp%dof(k, iq) = rep(td(k, i))
                end do
                p = surface%parameters(:, 1, i) + &
                    uv(1)*(surface%parameters(:, 2, i) - surface%parameters(:, 1, i)) + &
                    uv(2)*(surface%parameters(:, 3, i) - surface%parameters(:, 1, i))
                gp%angle(iq) = p(2); gp%edge(iq) = surface%panel_edge(i)
                deriv(:, 1) = surface%parameters(:, 2, i) - surface%parameters(:, 1, i)
                deriv(:, 2) = surface%parameters(:, 3, i) - surface%parameters(:, 1, i)
                call edge_basis(gp%order, p(1), w(:gp%order), dw(:gp%order))
                gp%shape(:, iq) = w(:gp%order)
                gp%shape_u(:, iq) = dw(:gp%order)*deriv(1, 1)
                gp%shape_v(:, iq) = dw(:gp%order)*deriv(1, 2)
                call gamma_p_edge_point(surface%coef(:, :, :, gp%edge(iq)), gp%order, &
                    gp%nfp, p(1), p(2), q(1:2), qt, qp)
                q(3) = p(2)
                call chart%jet(q, x, d1, d2, d3)
                gp%d1(:, :, iq) = d1(:, 1:2)
                gp%q(:, iq) = q
                do k = 1, 2
                    gp%qd(1:2, k, iq) = qt*deriv(1, k) + qp*deriv(2, k)
                    gp%qd(3, k, iq) = deriv(2, k)
                end do
                gp%angle_uv(:, iq) = deriv(2, :)
            end do
        end do
        status = GAMMA_P_OK
    end subroutine

    subroutine panels(surface, nedge, refine, nphi, orientation, faces)
        ! Strips in phi-hat outermost, so the first nf/nfp faces are the first
        ! field period (the circulant sector order).
        type(gamma_p_surface_t), intent(inout) :: surface
        integer, intent(in) :: nedge, refine, nphi, orientation
        integer, allocatable, intent(out) :: faces(:, :)
        integer :: ntheta, nf, p, j, after, e, sub, segment, next, ids(4), perm(3)
        real(dp) :: p0, p1, t0, t1
        ntheta = nedge*refine; nf = 2*ntheta*nphi
        perm = [1, 2, 3]
        if (orientation < 0) perm = [1, 3, 2]
        if (allocated(surface%parameters)) deallocate (surface%parameters, &
                                                       surface%panel_edge)
        allocate (surface%parameters(2, 3, nf), surface%panel_edge(nf), faces(3, nf))
        p = 0
        do j = 1, nphi
            after = mod(j, nphi) + 1
            p0 = 2*pi*real(j - 1, dp)/nphi; p1 = 2*pi*real(j, dp)/nphi
            do e = 1, nedge
                do sub = 1, refine
                    segment = (e - 1)*refine + sub
                    next = mod(segment, ntheta) + 1
                    t0 = real(sub - 1, dp)/refine; t1 = real(sub, dp)/refine
                    ids = [segment + (j - 1)*ntheta, segment + (after - 1)*ntheta, &
                           next + (j - 1)*ntheta, next + (after - 1)*ntheta]
                    p = p + 1; surface%panel_edge(p) = e
                    faces(:, p) = ids(perm)
                    surface%parameters(:, :, p) = reshape([t0, p0, t0, p1, t1, p0], [2, 3])
                    surface%parameters(:, :, p) = surface%parameters(:, perm, p)
                    p = p + 1; surface%panel_edge(p) = e
                    faces(:, p) = ids([2, 4, 3])
                    faces(:, p) = faces(perm, p)
                    surface%parameters(:, :, p) = reshape([t0, p1, t1, p1, t1, p0], [2, 3])
                    surface%parameters(:, :, p) = surface%parameters(:, perm, p)
                end do
            end do
        end do
    end subroutine

    pure function cross(u, v) result(c)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: c(3)
        c = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function

    subroutine external_field(gp, field, x, b, status)
        ! B_e at points of the period: Cartesian, or the n = 0 harmonic in
        ! (R, phi, Z) = Cartesian at phi = 0 for the axisymmetric exterior.
        type(gamma_p_t), intent(in) :: gp
        class(gamma_p_field_t), intent(in) :: field
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :)
        integer, intent(out) :: status
        if (gp%kind == GAMMA_P_AXISYMMETRIC) then
            call field%axisymmetric(x([1, 3], :), b, status)
        else
            call field%points(x, b, status)
        end if
        if (status == 0 .and. .not. all(ieee_is_finite(b))) status = 1
    end subroutine

    subroutine gamma_p_vacuum(gp, field, circulation, mu0, energy, gradient, status)
        !! E_v per field period (J, up to a constant) and its Hadamard shape
        !! gradient gradient(2, nmode, nchain) with respect to the chain
        !! coefficients. Stores Phi, B on the points and the energy parts.
        type(gamma_p_t), intent(inout) :: gp
        class(gamma_p_field_t), intent(in) :: field
        real(dp), intent(in) :: circulation, mu0
        real(dp), intent(out) :: energy, gradient(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: be(:, :), bx(:, :), line(:, :), bl(:, :), xg(:), wg(:)
        real(dp) :: n(3), pu, pv, s(3), trig(gp%nmode), dtrig(gp%nmode), r, er(3)
        real(dp) :: flux, volume, filament, f, rr, fsum
        integer :: iq, k, l, c, m, j, st, info
        interface
            subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(in) :: a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        status = GAMMA_P_BAD_INPUT; energy = 0
        if (.not. allocated(gp%factor)) return
        if (any(shape(gradient) /= [2, gp%nmode, gp%nchain])) return
        if (.not. ieee_is_finite(circulation) .or. .not. ieee_is_finite(mu0)) return
        if (mu0 <= 0) return
        gradient = 0
        status = GAMMA_P_FIELD_FAILED
        allocate (be(3, gp%nq))
        call external_field(gp, field, gp%x, be, st)
        if (st /= 0) return
        bx = be + circulation*gp%bf
        gp%bx = bx; gp%circulation = circulation; gp%mu0 = mu0
        if (allocated(gp%dbx)) deallocate (gp%dbx)
        ! Load h and the filament terms.
        if (allocated(gp%load)) deallocate (gp%load)
        allocate (gp%load(gp%ntrace))
        gp%load = 0; filament = 0
        do iq = 1, gp%nq
            call gamma_p_flux(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                gp%b(2, iq), gp%b(3, iq), bx(1, iq), bx(2, iq), bx(3, iq), flux)
            gp%load(gp%dof(:, iq)) = gp%load(gp%dof(:, iq)) + gp%w(iq)*flux*gp%basis(:, iq)
            n = cross(gp%a(:, iq), gp%b(:, iq))
            filament = filament - gp%w(iq)*dot_product(n, &
                circulation*cross(gp%af(:, iq), be(:, iq)) &
                + 0.5_dp*circulation**2*cross(gp%af(:, iq), gp%bf(:, iq)))
        end do
        status = GAMMA_P_SOLVE_FAILED
        gp%phi = gp%load
        call dpotrs('L', gp%ntrace, 1, gp%factor, gp%ntrace, gp%phi, gp%ntrace, info)
        if (info /= 0) return
        ! Volume term along R lines from r_ref at fixed (Z, phi).
        status = GAMMA_P_FIELD_FAILED
        allocate (xg(gp%nline), wg(gp%nline), line(3, gp%nline*gp%nq), &
                  bl(3, gp%nline*gp%nq))
        call gauss_legendre(gp%nline, xg, wg)
        xg = 0.5_dp*(xg + 1); wg = 0.5_dp*wg
        do iq = 1, gp%nq
            r = sqrt(gp%x(1, iq)**2 + gp%x(2, iq)**2)
            do l = 1, gp%nline
                rr = gp%r_ref + xg(l)*(r - gp%r_ref)
                line(:, (iq - 1)*gp%nline + l) = [rr*gp%x(1, iq)/r, rr*gp%x(2, iq)/r, &
                                                  gp%x(3, iq)]
            end do
        end do
        call external_field(gp, field, line, bl, st)
        if (st /= 0) return
        volume = 0
        if (allocated(gp%field)) deallocate (gp%field)
        allocate (gp%field(3, gp%nq))
        do iq = 1, gp%nq
            r = sqrt(gp%x(1, iq)**2 + gp%x(2, iq)**2)
            er = [gp%x(1, iq)/r, gp%x(2, iq)/r, 0.0_dp]
            fsum = 0
            do l = 1, gp%nline
                k = (iq - 1)*gp%nline + l
                rr = gp%r_ref + xg(l)*(r - gp%r_ref)
                fsum = fsum + wg(l)*rr*sum(bl(:, k)**2)
            end do
            f = fsum*(r - gp%r_ref)/r
            n = cross(gp%a(:, iq), gp%b(:, iq))
            volume = volume + gp%w(iq)*f*dot_product(er, n)
            pu = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_u(:, iq))
            pv = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_v(:, iq))
            call gamma_p_field(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                gp%b(2, iq), gp%b(3, iq), bx(1, iq), bx(2, iq), bx(3, iq), pu, pv, &
                gp%field(1, iq), gp%field(2, iq), gp%field(3, iq))
            call gamma_p_traction(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                gp%b(2, iq), gp%b(3, iq), bx(1, iq), bx(2, iq), bx(3, iq), pu, pv, &
                s(1), s(2), s(3))
            call trig_modes(gp%nmode, gp%nfp, gp%angle(iq), trig, dtrig)
            do l = 0, gp%order
                j = chain_node(gp%nchain, gp%order, gp%edge(iq), l)
                do m = 1, gp%nmode
                    do c = 1, 2
                        gradient(c, m, j) = gradient(c, m, j) + gp%w(iq)*gp%shape(l, iq) &
                            *trig(m)*dot_product(s, gp%d1(:, c, iq))
                    end do
                end do
            end do
        end do
        gp%energy_parts = [-0.5_dp*dot_product(gp%load, gp%phi), -0.5_dp*volume, filament]
        energy = sum(gp%energy_parts)/mu0
        gradient = gradient/mu0
        status = GAMMA_P_OK
        if (.not. ieee_is_finite(energy) .or. .not. all(ieee_is_finite(gradient))) &
            status = GAMMA_P_SOLVE_FAILED
    end subroutine

    subroutine factor_mass(gp, status)
        type(gamma_p_t), intent(inout) :: gp
        integer, intent(out) :: status
        integer :: iq, i, j, info
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        status = GAMMA_P_SOLVE_FAILED
        allocate (gp%mass_factor(gp%ntrace, gp%ntrace))
        gp%mass_factor = 0
        do iq = 1, gp%nq
            do j = 1, gp%nb
                do i = 1, gp%nb
                    gp%mass_factor(gp%dof(i, iq), gp%dof(j, iq)) = &
                        gp%mass_factor(gp%dof(i, iq), gp%dof(j, iq)) + gp%w(iq) &
                        *norm2(cross(gp%a(:, iq), gp%b(:, iq)))*gp%basis(i, iq)*gp%basis(j, iq)
                end do
            end do
        end do
        call dpotrf('L', gp%ntrace, gp%mass_factor, gp%ntrace, info)
        if (info /= 0) return
        status = GAMMA_P_OK
    end subroutine

    subroutine field_jacobian(gp, field, status)
        ! dbx(i, k, q) = dB_x,i/dx_k by central differences (step 1e-6 of
        ! the largest coordinate). The axisymmetric exterior has an n = 0
        ! field, so its y column at phi = 0 is the rotation (-B_phi, B_R, 0)/R.
        type(gamma_p_t), intent(inout) :: gp
        class(gamma_p_field_t), intent(in) :: field
        integer, intent(out) :: status
        real(dp), allocatable :: xs(:, :), bs(:, :), potential(:, :, :), bfs(:, :, :)
        real(dp) :: step
        integer :: iq, k, sgn, col, st
        status = GAMMA_P_FIELD_FAILED
        allocate (xs(3, 6*gp%nq), bs(3, 6*gp%nq), gp%dbx(3, 3, gp%nq))
        step = 1.0e-6_dp*maxval(abs(gp%x))
        do iq = 1, gp%nq
            do k = 1, 3
                do sgn = 1, 2
                    col = 6*(iq - 1) + 2*(k - 1) + sgn
                    xs(:, col) = gp%x(:, iq)
                    xs(k, col) = xs(k, col) + (3 - 2*sgn)*step
                end do
            end do
        end do
        if (gp%kind == GAMMA_P_AXISYMMETRIC) xs(2, :) = 0
        call external_field(gp, field, xs, bs, st)
        if (st /= 0) return
        if (allocated(gp%filament) .and. abs(gp%circulation) > 0) then
            allocate (potential(3, 6*gp%nq, 1), bfs(3, 6*gp%nq, 1))
            call evaluate_current_loops(gp%filament, [1, size(gp%filament, 2) + 1], &
                                        xs, potential, bfs, st)
            if (st /= 0) return
            bs = bs + gp%circulation*bfs(:, :, 1)
        end if
        do iq = 1, gp%nq
            do k = 1, 3
                col = 6*(iq - 1) + 2*(k - 1)
                gp%dbx(:, k, iq) = (bs(:, col + 1) - bs(:, col + 2))/(2*step)
            end do
            if (gp%kind == GAMMA_P_AXISYMMETRIC) gp%dbx(:, 2, iq) = &
                [-gp%bx(2, iq), gp%bx(1, iq), 0.0_dp]/gp%x(1, iq)
        end do
        status = GAMMA_P_OK
    end subroutine

    subroutine point_displacement(gp, iq, direction, dq, dq_u, dq_v)
        ! Chart increment (dq1, dq2) of a chain direction at point iq and its
        ! u, v derivatives.
        type(gamma_p_t), intent(in) :: gp
        integer, intent(in) :: iq
        real(dp), intent(in) :: direction(:, :, :)
        real(dp), intent(out) :: dq(2), dq_u(2), dq_v(2)
        real(dp) :: trig(gp%nmode), dtrig(gp%nmode), c(2), ct(2)
        integer :: l, j
        call trig_modes(gp%nmode, gp%nfp, gp%angle(iq), trig, dtrig)
        dq = 0; dq_u = 0; dq_v = 0
        do l = 0, gp%order
            j = chain_node(gp%nchain, gp%order, gp%edge(iq), l)
            c = matmul(direction(:, :, j), trig)
            ct = matmul(direction(:, :, j), dtrig)
            dq = dq + gp%shape(l, iq)*c
            dq_u = dq_u + gp%shape_u(l, iq)*c + gp%shape(l, iq)*ct*gp%angle_uv(1, iq)
            dq_v = dq_v + gp%shape_v(l, iq)*c + gp%shape(l, iq)*ct*gp%angle_uv(2, iq)
        end do
    end subroutine

    subroutine gamma_p_hessian(gp, field, direction, product, status)
        !! Hessian product of E_v (derivative of the gradient of
        !! gamma_p_vacuum, called before at the same state) along a chain
        !! direction(2, nmode, nchain). The trace potential follows its
        !! material derivative Phi' + dx.grad Phi: the shape derivative Phi'
        !! solves the exterior Neumann problem dPhi'/dn = div_Gamma(V_n B),
        !! V_n = dx.n, weakly S Phi' = int V_n B.grad_Gamma N_i dA with the
        !! same S, and the transport dx.grad Phi = dx.grad_Gamma Phi -
        !! (dx.n)(B_x.n) is projected on the trace space. Consistent with the
        !! discrete gradient to the order of the BEM, like that gradient.
        type(gamma_p_t), intent(inout) :: gp
        class(gamma_p_field_t), intent(in) :: field
        real(dp), intent(in) :: direction(:, :, :)
        real(dp), intent(out) :: product(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: dx(:, :), da(:, :), db(:, :), dqs(:, :)
        real(dp), allocatable :: shift(:), transport(:)
        real(dp) :: dq(2), dq_u(2), dq_v(2), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: n(3), nh(3), gn(3), gphi(3), pu, pv, dpu, dpv, s(3), ds(3), de(3)
        real(dp) :: trig(gp%nmode), dtrig(gp%nmode)
        integer :: iq, i, k, l, j, m, c, st, info
        interface
            subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(in) :: a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        status = GAMMA_P_BAD_INPUT
        if (.not. allocated(gp%phi) .or. .not. allocated(gp%bx)) return
        if (any(shape(direction) /= [2, gp%nmode, gp%nchain])) return
        if (any(shape(product) /= [2, gp%nmode, gp%nchain])) return
        if (.not. allocated(gp%dbx)) then
            call field_jacobian(gp, field, st)
            status = GAMMA_P_FIELD_FAILED
            if (st /= GAMMA_P_OK) return
        end if
        allocate (dx(3, gp%nq), da(3, gp%nq), db(3, gp%nq), dqs(2, gp%nq), &
                  shift(gp%ntrace), transport(gp%ntrace))
        shift = 0; transport = 0
        do iq = 1, gp%nq
            call point_displacement(gp, iq, direction, dq, dq_u, dq_v)
            call gp%chart%jet(gp%q(:, iq), x, d1, d2, d3)
            dqs(:, iq) = dq
            dx(:, iq) = matmul(d1(:, 1:2), dq)
            da(:, iq) = matmul(d1(:, 1:2), dq_u)
            db(:, iq) = matmul(d1(:, 1:2), dq_v)
            do k = 1, 3
                do j = 1, 2
                    da(:, iq) = da(:, iq) + d2(:, j, k)*dq(j)*gp%qd(k, 1, iq)
                    db(:, iq) = db(:, iq) + d2(:, j, k)*dq(j)*gp%qd(k, 2, iq)
                end do
            end do
            n = cross(gp%a(:, iq), gp%b(:, iq)); nh = n/norm2(n)
            pu = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_u(:, iq))
            pv = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_v(:, iq))
            call gamma_p_field(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                gp%b(2, iq), gp%b(3, iq), 0.0_dp, 0.0_dp, 0.0_dp, pu, pv, &
                gphi(1), gphi(2), gphi(3))
            do i = 1, gp%nb
                call gamma_p_field(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                    gp%b(2, iq), gp%b(3, iq), 0.0_dp, 0.0_dp, 0.0_dp, gp%basis_u(i, iq), &
                    gp%basis_v(i, iq), gn(1), gn(2), gn(3))
                shift(gp%dof(i, iq)) = shift(gp%dof(i, iq)) + gp%w(iq) &
                    *dot_product(dx(:, iq), n)*dot_product(gp%field(:, iq), gn)
            end do
            transport(gp%dof(:, iq)) = transport(gp%dof(:, iq)) + gp%w(iq)*norm2(n) &
                *gp%basis(:, iq)*(dot_product(dx(:, iq), gphi) &
                - dot_product(dx(:, iq), nh)*dot_product(gp%bx(:, iq), nh))
        end do
        status = GAMMA_P_SOLVE_FAILED
        call dpotrs('L', gp%ntrace, 1, gp%factor, gp%ntrace, shift, gp%ntrace, info)
        if (info /= 0) return
        call dpotrs('L', gp%ntrace, 1, gp%mass_factor, gp%ntrace, transport, &
                    gp%ntrace, info)
        if (info /= 0) return
        shift = shift + transport
        product = 0
        do iq = 1, gp%nq
            call gp%chart%jet(gp%q(:, iq), x, d1, d2, d3)
            pu = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_u(:, iq))
            pv = dot_product(gp%phi(gp%dof(:, iq)), gp%basis_v(:, iq))
            dpu = dot_product(shift(gp%dof(:, iq)), gp%basis_u(:, iq))
            dpv = dot_product(shift(gp%dof(:, iq)), gp%basis_v(:, iq))
            de = matmul(gp%dbx(:, :, iq), dx(:, iq))
            call gamma_p_traction(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), gp%b(1, iq), &
                gp%b(2, iq), gp%b(3, iq), gp%bx(1, iq), gp%bx(2, iq), gp%bx(3, iq), pu, &
                pv, s(1), s(2), s(3))
            call gamma_p_traction_jvp(gp%a(1, iq), gp%a(2, iq), gp%a(3, iq), &
                gp%b(1, iq), gp%b(2, iq), gp%b(3, iq), gp%bx(1, iq), gp%bx(2, iq), &
                gp%bx(3, iq), pu, pv, da(1, iq), da(2, iq), da(3, iq), db(1, iq), &
                db(2, iq), db(3, iq), de(1), de(2), de(3), dpu, dpv, ds(1), ds(2), ds(3))
            call trig_modes(gp%nmode, gp%nfp, gp%angle(iq), trig, dtrig)
            do l = 0, gp%order
                j = chain_node(gp%nchain, gp%order, gp%edge(iq), l)
                do m = 1, gp%nmode
                    do c = 1, 2
                        ! d2 term: the chart column dx/dq_c moves with q.
                        product(c, m, j) = product(c, m, j) + gp%w(iq)*gp%shape(l, iq) &
                            *trig(m)*(dot_product(ds, d1(:, c)) + dot_product(s, &
                            d2(:, c, 1)*dqs(1, iq) + d2(:, c, 2)*dqs(2, iq)))
                    end do
                end do
            end do
        end do
        product = product/gp%mu0
        status = GAMMA_P_OK
        if (.not. all(ieee_is_finite(product))) status = GAMMA_P_SOLVE_FAILED
    end subroutine
end module milca_ideal_fb_gamma_p
