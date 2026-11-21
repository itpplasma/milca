module milca_ideal_fb_bem
    ! Axisymmetric free-boundary vacuum on the plasma boundary Gamma_p
    ! (theory/derivations/free-boundary-plasma-boundary-bem.md).  The vacuum
    ! is the field of the coils plus a sheet current tau/mu0 on Gamma_p; its
    ! contribution to the Lagrangian (per full turn) is
    !
    !   L_v = -(2 pi/mu0) (tau' V tau/2 + tau' g) + (psi_hat/mu0) (s' tau - c),
    !
    ! with V_ij = int int G phi_i phi_j, g_i = int psi_c phi_i, s_i = int phi_i
    ! on the curved panels through the plasma boundary nodes, c = mu0 I_phi
    ! (toroidal current along +phi; MILCA's circulation is -c).
    ! The coil flux psi_c is a fit of ring sources to the mgrid field.  This
    ! module sees Gamma_p only as the coordinates of its chain nodes; the
    ! plasma map and the unknown layout belong to milca_ideal_free_boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_gs_single_layer_shape_2d, only: gs_sl_curve_t, gs_sl_plan_t, &
        gs_sl_curve_create, gs_sl_plan_create, gs_sl_assemble, gs_sl_load
    use fortfem_axisymmetric_gs_kernel, only: gs_green_hessian
    implicit none
    real(dp), public :: fb_bem_seconds(0:2) = 0
    integer, public :: fb_bem_calls(0:2) = 0
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    ! psi_c(x) = sum_j coef(j) G(x, source(:, j)) (Wb/rad).
    type, public :: fb_coil_flux_t
        real(dp), allocatable :: source(:, :), coef(:)
        ! relative l2 misfit of R B on the samples
        real(dp) :: misfit = 0
    end type

    type, public :: fb_bem_t
        integer :: degree = 0, flux_degree = 0, nb = 0, nt = 0
        ! plasma node of each chain node, counter-clockwise in (R, Z)
        integer, allocatable :: node(:)
        type(gs_sl_curve_t) :: curve
        type(gs_sl_plan_t) :: plan
        type(fb_coil_flux_t) :: coil
        real(dp), allocatable :: tau(:)
        real(dp) :: psi_hat = 0
    end type

    public :: fb_coil_flux_fit, fb_coil_flux_value
    public :: fb_bem_setup, fb_bem_replan, fb_bem_terms, fb_bem_exact

contains

    subroutine fb_coil_flux_fit(points, field, sources, coil, status)
        ! Least-squares fit of psi_c to the poloidal field samples:
        ! d psi/dR = R B_Z, d psi/dZ = -R B_R at points(:, i), field(:, i) =
        ! (B_R, B_Z).  The sources must lie outside the region of the points.
        real(dp), intent(in) :: points(:, :), field(:, :), sources(:, :)
        type(fb_coil_flux_t), intent(out) :: coil
        integer, intent(out) :: status
        real(dp), allocatable :: a(:, :), b(:), work(:), sv(:)
        real(dp) :: g, grad(4), hess(4, 4), rcond
        integer :: i, j, m, n, rank, lwork, info
        interface
            subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, &
                    lwork, info)
                import :: dp
                integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                real(dp), intent(out) :: s(*), work(*)
                real(dp), intent(in) :: rcond
                integer, intent(out) :: rank, info
            end subroutine dgelss
        end interface
        status = 1
        n = size(sources, 2)
        m = 2*size(points, 2)
        if (m < n .or. n < 1) return
        allocate (a(m, n), b(max(m, n)), sv(n))
        do i = 1, size(points, 2)
            do j = 1, n
                call gs_green_hessian(points(1, i), points(2, i), sources(1, j), &
                    sources(2, j), g, grad, hess)
                a(2*i - 1, j) = grad(1)
                a(2*i, j) = grad(2)
            end do
            b(2*i - 1) = points(1, i)*field(2, i)
            b(2*i) = -points(1, i)*field(1, i)
        end do
        rcond = 1.0e-13_dp
        allocate (work(1))
        call dgelss(m, n, 1, a, m, b, max(m, n), sv, rcond, rank, work, -1, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))
        a = 0
        do i = 1, size(points, 2)
            do j = 1, n
                call gs_green_hessian(points(1, i), points(2, i), sources(1, j), &
                    sources(2, j), g, grad, hess)
                a(2*i - 1, j) = grad(1)
                a(2*i, j) = grad(2)
            end do
        end do
        call dgelss(m, n, 1, a, m, b, max(m, n), sv, rcond, rank, work, lwork, info)
        if (info /= 0) return
        coil%source = sources
        coil%coef = b(1:n)
        coil%misfit = fit_misfit(coil, points, field)
        status = 0
    end subroutine fb_coil_flux_fit

    function fit_misfit(coil, points, field) result(misfit)
        type(fb_coil_flux_t), intent(in) :: coil
        real(dp), intent(in) :: points(:, :), field(:, :)
        real(dp) :: misfit, v, gr(2), h(3), num, den
        integer :: i

        num = 0; den = 0
        do i = 1, size(points, 2)
            call fb_coil_flux_value(coil, points(:, i), v, gr, h)
            num = num + (gr(1) - points(1, i)*field(2, i))**2 &
                + (gr(2) + points(1, i)*field(1, i))**2
            den = den + points(1, i)**2*sum(field(:, i)**2)
        end do
        misfit = sqrt(num/max(den, tiny(1.0_dp)))
    end function fit_misfit

    pure subroutine fb_coil_flux_value(coil, x, value, gradient, hessian)
        ! psi_c with its gradient and Hessian (RR, RZ, ZZ) at x = (R, Z).
        type(fb_coil_flux_t), intent(in) :: coil
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: value, gradient(2), hessian(3)
        real(dp) :: g, grad(4), hess(4, 4)
        integer :: j

        value = 0; gradient = 0; hessian = 0
        do j = 1, size(coil%coef)
            call gs_green_hessian(x(1), x(2), coil%source(1, j), coil%source(2, j), &
                g, grad, hess)
            value = value + coil%coef(j)*g
            gradient = gradient + coil%coef(j)*grad(1:2)
            hessian = hessian + coil%coef(j)*[hess(1, 1), hess(1, 2), hess(2, 2)]
        end do
    end subroutine fb_coil_flux_value

    subroutine fb_bem_setup(bem, node, coords, degree, flux_degree, status)
        ! Chain of the plasma boundary (node, in panel order) with its
        ! coordinates (R, Z); a clockwise chain is reversed.
        type(fb_bem_t), intent(inout) :: bem
        integer, intent(in) :: node(:), degree, flux_degree
        real(dp), intent(in) :: coords(:, :)
        integer, intent(out) :: status
        real(dp) :: area
        integer :: b, nb

        status = 1
        nb = size(node)
        if (mod(nb, degree) /= 0) return
        area = 0
        do b = 1, nb
            area = area + coords(1, b)*coords(2, modulo(b, nb) + 1) &
                - coords(1, modulo(b, nb) + 1)*coords(2, b)
        end do
        bem%degree = degree
        bem%flux_degree = flux_degree
        bem%nb = nb
        if (area > 0) then
            bem%node = node
        else
            ! keep node 1, reverse the rest: panel order stays valid
            bem%node = [node(1), node(nb:2:-1)]
        end if
        bem%nt = (flux_degree + 1)*(nb/degree)
        if (.not. allocated(bem%tau)) then
            allocate (bem%tau(bem%nt), source=0.0_dp)
        end if
        call fb_bem_replan(bem, chain_coords(bem, coords, node), status)
    end subroutine fb_bem_setup

    function chain_coords(bem, coords, node) result(x)
        ! coords given in the order of node, returned in bem%node order.
        type(fb_bem_t), intent(in) :: bem
        real(dp), intent(in) :: coords(:, :)
        integer, intent(in) :: node(:)
        real(dp) :: x(2, size(node))
        integer :: b, k

        do b = 1, size(node)
            k = findloc(node, bem%node(b), dim=1)
            x(:, b) = coords(:, k)
        end do
    end function chain_coords

    subroutine fb_bem_replan(bem, x, status)
        ! Curve through the chain coordinates x (bem%node order) and a
        ! quadrature plan frozen at this geometry.
        type(fb_bem_t), intent(inout) :: bem
        real(dp), intent(in) :: x(:, :)
        integer, intent(out) :: status

        call gs_sl_curve_create(bem%degree, x, bem%curve, status)
        if (status /= 0) return
        call gs_sl_plan_create(bem%curve, bem%flux_degree, bem%plan, status)
    end subroutine fb_bem_replan

    subroutine fb_bem_terms(bem, x, mu0, current, energy, gx, gt, gp, hxx, hxt, &
            hxp, htt, htp)
        ! L_v and its derivatives at the chain coordinates x (bem%node order)
        ! and the present tau, psi_hat; current = mu0 I_p.  gx, hxx, hxt, hxp
        ! are in the flattened coordinates X(2(b-1)+c).
        type(fb_bem_t), intent(in) :: bem
        real(dp), intent(in) :: x(:, :), mu0, current
        real(dp), intent(out) :: energy, gx(:), gt(:), gp
        real(dp), intent(out), optional :: hxx(:, :), hxt(:, :), hxp(:), htt(:, :), &
            htp(:)
        type(gs_sl_curve_t) :: curve
        real(dp), allocatable :: v(:, :), dvt(:, :), gv(:), hv(:, :)
        real(dp), allocatable :: g(:), dg(:, :), gg(:), hg(:, :)
        real(dp), allocatable :: s(:), ds(:, :), gs(:), hs(:, :)
        real(dp) :: k
        integer :: nt, nx

        nt = bem%nt
        nx = 2*bem%nb
        curve = bem%curve
        curve%nodes = x
        allocate (v(nt, nt), dvt(nt, nx), gv(nx), hv(nx, nx), g(nt), dg(nt, nx), &
            gg(nx), hg(nx, nx), s(nt), ds(nt, nx), gs(nx), hs(nx, nx))
        block
            integer(8) :: c0, c1, cr
            call system_clock(c0, cr)
            if (present(hxx)) then
                call gs_sl_assemble(curve, bem%plan, v, bem%tau, dvt, gv, hv)
            else
                call gs_sl_assemble(curve, bem%plan, v, bem%tau, dvt, gv)
                hv = 0
            end if
            call system_clock(c1)
            fb_bem_seconds(merge(2, 1, present(hxx))) = fb_bem_seconds(merge(2, 1, present(hxx))) + real(c1 - c0, dp)/cr
            fb_bem_calls(merge(2, 1, present(hxx))) = fb_bem_calls(merge(2, 1, present(hxx))) + 1
        end block
        call gs_sl_load(curve, bem%flux_degree, coil_field, g, bem%tau, dg, gg, hg)
        call gs_sl_load(curve, bem%flux_degree, unit_field, s, bem%tau, ds, gs, hs)
        k = 2*pi/mu0
        energy = -k*(0.5_dp*dot_product(bem%tau, matmul(v, bem%tau)) &
            + dot_product(bem%tau, g)) &
            + bem%psi_hat/mu0*(dot_product(s, bem%tau) - current)
        gx = -k*(gv + gg) + bem%psi_hat/mu0*gs
        gt = -k*(matmul(v, bem%tau) + g) + bem%psi_hat/mu0*s
        gp = (dot_product(s, bem%tau) - current)/mu0
        if (present(hxx)) hxx = -k*(hv + hg) + bem%psi_hat/mu0*hs
        if (present(hxt)) hxt = transpose(-k*(dvt + dg) + bem%psi_hat/mu0*ds)
        if (present(hxp)) hxp = gs/mu0
        if (present(htt)) htt = -k*v
        if (present(htp)) htp = s/mu0
    contains
        subroutine coil_field(p, value, gradient, hessian)
            real(dp), intent(in) :: p(2)
            real(dp), intent(out) :: value, gradient(2), hessian(3)
            call fb_coil_flux_value(bem%coil, p, value, gradient, hessian)
        end subroutine coil_field
    end subroutine fb_bem_terms

    subroutine unit_field(p, value, gradient, hessian)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: value, gradient(2), hessian(3)

        value = 1 + 0*p(1)
        gradient = 0
        hessian = 0
    end subroutine unit_field

    subroutine fb_bem_exact(bem, x, current, status)
        ! tau and psi_hat stationary at fixed geometry:
        !   2 pi V tau - psi_hat s = -2 pi g,  s' tau = current.
        type(fb_bem_t), intent(inout) :: bem
        real(dp), intent(in) :: x(:, :), current
        integer, intent(out) :: status
        type(gs_sl_curve_t) :: curve
        real(dp), allocatable :: v(:, :), g(:), s(:), a(:, :), rhs(:)
        integer, allocatable :: ipiv(:)
        integer :: nt, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface
        nt = bem%nt
        curve = bem%curve
        curve%nodes = x
        allocate (v(nt, nt), g(nt), s(nt), a(nt + 1, nt + 1), rhs(nt + 1), &
            ipiv(nt + 1))
        block
            integer(8) :: c0, c1, cr
            call system_clock(c0, cr)
            call gs_sl_assemble(curve, bem%plan, v)
            call system_clock(c1)
            fb_bem_seconds(0) = fb_bem_seconds(0) + real(c1 - c0, dp)/cr
            fb_bem_calls(0) = fb_bem_calls(0) + 1
        end block
        call gs_sl_load(curve, bem%flux_degree, coil_field, g)
        call gs_sl_load(curve, bem%flux_degree, unit_field, s)
        a = 0
        a(1:nt, 1:nt) = 2*pi*0.5_dp*(v + transpose(v))
        a(1:nt, nt + 1) = -s
        a(nt + 1, 1:nt) = s
        rhs(1:nt) = -2*pi*g
        rhs(nt + 1) = current
        call dgesv(nt + 1, 1, a, nt + 1, ipiv, rhs, nt + 1, info)
        status = info
        if (info /= 0) return
        bem%tau = rhs(1:nt)
        bem%psi_hat = rhs(nt + 1)
    contains
        subroutine coil_field(p, value, gradient, hessian)
            real(dp), intent(in) :: p(2)
            real(dp), intent(out) :: value, gradient(2), hessian(3)
            call fb_coil_flux_value(bem%coil, p, value, gradient, hessian)
        end subroutine coil_field
    end subroutine fb_bem_exact

end module milca_ideal_fb_bem
