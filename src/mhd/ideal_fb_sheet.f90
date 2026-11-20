module milca_ideal_fb_sheet
    ! Vacuum of a 3D free boundary on the plasma boundary Gamma_p as a surface
    ! current (theory/derivations/free-boundary-plasma-boundary-sheet-3d.md).
    ! Its contribution to the Lagrangian (per full turn) is L_v = -e/mu0,
    !
    !   e = E_self[Phi] - E_self[Phi_F] - int phi B_c . N du dv + c <Psi_t(u)>_u
    !       - (j/2 pi) int Psi~_p(v) dv,
    !
    ! Phi = phi + c u + j v/(2 pi) (u in [0, 1) along the boundary chain, v
    ! in [0, 2 pi)), E_self the energy of K = n x grad Phi, c = mu0 I_phi,
    ! j = -2 pi F o (o = 1 for a counter-clockwise chain in (R, Z)), Psi_t(u)
    ! the coil flux through the ribbon from a fixed curve L inside the plasma
    ! to the toroidal loop at u and Psi_p(v) through the fan from L(v) to the
    ! section at v. Every term is periodic: the net currents enter through
    ! the fluxes averaged over the loops, not through a cut. The poloidal
    ! sheet current Phi_F = j v/(2 pi) that cancels the coil toroidal field
    ! F/R in the plasma contributes E_self[Phi_F] - (j/2 pi) int Psi_p^F dv =
    ! -F^2/2 int_D |grad phi|^2, two large terms that cancel to this volume
    ! integral; it is left to the plasma side (exact in the plasma
    ! quadrature), E_self[Phi_F] is subtracted on the same grid and Psi~_p
    ! holds only the part of B_phi beyond F/R.
    !
    ! Discretisation: the boundary nodes lie on the Fourier curve x_b = F y
    ! (kb harmonics in the chain parameter u_b = (b - 1)/nb, the plasma
    ! holds them there); the surface is that curve, sampled on a uniform grid
    ! of nu = 4 kb + 1 points and nv per field period, X = U F^+ x_b with exact
    ! u derivatives, times the map's Fourier modes in v. E_self is the corrected
    ! trapezoidal rule (fortfem_toroidal_sheet_trapezoid), the load and the
    ! fluxes the trapezoidal rule on the same grid (Gauss along s). The grid
    ! cells are at most twice as long toroidally as poloidally (the rule
    ! loses accuracy on longer cells). phi(k, n) has Fourier modes in u up to
    ! about nu/4, so the grid resolves it (the rule is accurate only for
    ! resolved densities; with nodal phi its block is indefinite), and the
    ! map's toroidal modes; phi(1, 1) is the gauge. e is quadratic in phi:
    ! phi is eliminated exactly with the Cholesky factor of its block.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_toroidal_sheet_trapezoid, only: tst_grid_t, tst_energy, tst_hvp, &
        tst_current_modes
    use fortfem_generated_toroidal_sheet_load_jet, only: generated_toroidal_sheet_load_jet
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gauss_unit
    use milca_sheet_flux_leaves, only: fan_flux_jet, ribbon_flux_jet
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    abstract interface
        subroutine sheet_field(r, z, v, b, db, d2b)
            !! (B_R, B_phi, B_Z) at (R, Z, phi = v), derivatives in (R, Z) and
            !! second derivatives (RR, RZ, ZZ).
            import :: dp
            real(dp), intent(in) :: r, z, v
            real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        end subroutine sheet_field
    end interface

    type, public :: fb_sheet_t
        integer :: nb = 0, nmode = 1, nfp = 1, nu = 0, nv = 0
        ! plasma node of each chain node
        integer, allocatable :: node(:)
        type(tst_grid_t) :: grid
        ! the surface at the u samples from the node values, and its u
        ! derivative: smat = U F^+, dsmat = U' F^+ (nu x nb)
        real(dp), allocatable :: smat(:, :), dsmat(:, :)
        ! spectral derivative matrix on the u samples (period 1)
        real(dp), allocatable :: du(:, :)
        ! Fourier modes and their v derivative at v_j over the torus
        real(dp), allocatable :: t(:, :), dt(:, :)
        ! current potential phi(k, n): u mode k (values ub(i, k) at u_i),
        ! toroidal mode n
        integer :: nph = 0
        real(dp), allocatable :: ub(:, :), dub(:, :), phi(:, :)
        ! L: (R, Z) per geometry mode; o: orientation of the chain
        real(dp), allocatable :: axis(:, :)
        real(dp) :: o = 1
        integer :: ng_s = 8
        ! Cholesky factor of the phi block of e (gauge removed)
        real(dp), allocatable :: chol(:, :)
    end type

    public :: sheet_field
    public :: fb_sheet_setup, fb_sheet_terms, fb_sheet_hvp, fb_sheet_exact
    public :: fb_sheet_reduced_hvp, fb_sheet_fourier

contains

    subroutine fb_sheet_setup(sheet, node, x, axis, nfp, ntor, kb, status)
        ! node: plasma node of each chain node (nb = size(node)); x(2, nb,
        ! nmode) their (R, Z) coefficients, on the Fourier curve of kb
        ! harmonics (2 kb + 1 <= nb); axis(2, nmode) the fixed curve L. nu =
        ! 4 kb + 1 resolves the curve and phi; nv makes the grid cells at most
        ! twice as long toroidally (the rule loses accuracy on longer cells) and
        ! resolves the toroidal modes.
        type(fb_sheet_t), intent(out) :: sheet
        integer, intent(in) :: node(:), nfp, ntor, kb
        real(dp), intent(in) :: x(:, :, :), axis(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: f(:, :), fp(:, :), u(:, :), du(:, :)
        real(dp) :: area, perimeter, v
        integer :: nb, b, i, j, k
        nb = size(node)
        status = 1
        if (kb < 1 .or. 2*kb + 1 > nb) return
        if (any(shape(x) /= [2, nb, 2*ntor + 1])) return
        if (any(shape(axis) /= [2, 2*ntor + 1])) return
        sheet%nb = nb; sheet%nmode = 2*ntor + 1; sheet%nfp = nfp
        sheet%node = node
        sheet%axis = axis
        area = 0; perimeter = 0
        do b = 1, nb
            area = area + x(1, b, 1)*x(2, 1 + mod(b, nb), 1) &
                - x(1, 1 + mod(b, nb), 1)*x(2, b, 1)
            perimeter = perimeter + norm2(x(:, 1 + mod(b, nb), 1) - x(:, b, 1))
        end do
        sheet%o = sign(1.0_dp, area)
        sheet%nu = 4*kb + 1
        sheet%nv = max(4*ntor + 2, ceiling(2*pi*maxval(x(1, :, 1))*sheet%nu &
            /(2*nfp*perimeter)))
        sheet%grid = tst_grid_t(nu=sheet%nu, nv=sheet%nv, nfp=nfp)
        ! F at the nodes, F^+ = diag(1/F'F) F' (orthogonal columns), U and U'
        call fb_sheet_fourier(nb, kb, f)
        fp = transpose(f)
        do k = 1, 2*kb + 1
            fp(k, :) = fp(k, :)/sum(f(:, k)**2)
        end do
        allocate (u(sheet%nu, 2*kb + 1), du(sheet%nu, 2*kb + 1))
        do i = 1, sheet%nu
            u(i, 1) = 1; du(i, 1) = 0
            do k = 1, kb
                u(i, 2*k) = cos(2*pi*k*(i - 1)/sheet%nu)
                u(i, 2*k + 1) = sin(2*pi*k*(i - 1)/sheet%nu)
                du(i, 2*k) = -2*pi*k*u(i, 2*k + 1)
                du(i, 2*k + 1) = 2*pi*k*u(i, 2*k)
            end do
        end do
        sheet%smat = matmul(u, fp)
        sheet%dsmat = matmul(du, fp)
        sheet%du = spectral_derivative(sheet%nu, 2*pi)
        ! phi: the same kb harmonics in u, resolved by the grid
        sheet%nph = 2*kb + 1
        sheet%ub = u
        sheet%dub = du
        allocate (sheet%t(sheet%nmode, sheet%nv*nfp), sheet%dt(sheet%nmode, sheet%nv*nfp))
        do j = 1, sheet%nv*nfp
            v = 2*pi*(j - 1)/(nfp*sheet%nv)
            call trig(sheet%nmode, nfp, v, sheet%t(:, j), sheet%dt(:, j))
        end do
        allocate (sheet%phi(sheet%nph, sheet%nmode), source=0.0_dp)
        status = 0
    end subroutine fb_sheet_setup

    pure subroutine fb_sheet_fourier(nb, kb, f)
        ! The real Fourier basis of kb harmonics at u_b = (b - 1)/nb.
        integer, intent(in) :: nb, kb
        real(dp), allocatable, intent(out) :: f(:, :)
        integer :: b, k
        allocate (f(nb, 2*kb + 1))
        do b = 1, nb
            f(b, 1) = 1
            do k = 1, kb
                f(b, 2*k) = cos(2*pi*k*(b - 1)/nb)
                f(b, 2*k + 1) = sin(2*pi*k*(b - 1)/nb)
            end do
        end do
    end subroutine fb_sheet_fourier

    pure function spectral_derivative(n, scale) result(d)
        ! d/dt of the trigonometric interpolant through n equidistant samples
        ! of a period 2 pi/scale (the Nyquist mode of an even n dropped).
        integer, intent(in) :: n
        real(dp), intent(in) :: scale
        real(dp) :: d(n, n), x
        integer :: i, k
        d = 0
        do k = 1, n
            do i = 1, n
                if (i == k) cycle
                x = (i - k)*pi/n
                if (mod(n, 2) == 0) then
                    d(i, k) = scale*0.5_dp*(-1)**mod(abs(i - k), 2)/tan(x)
                else
                    d(i, k) = scale*0.5_dp*(-1)**mod(abs(i - k), 2)/sin(x)
                end if
            end do
        end do
    end function spectral_derivative

    pure subroutine trig(nmode, nfp, v, t, dt)
        ! T_1 = 1, T_2k = cos(k nfp v), T_2k+1 = sin(k nfp v) and dT/dv.
        integer, intent(in) :: nmode, nfp
        real(dp), intent(in) :: v
        real(dp), intent(out) :: t(nmode), dt(nmode)
        integer :: k
        t(1) = 1; dt(1) = 0
        do k = 1, (nmode - 1)/2
            t(2*k) = cos(k*nfp*v); t(2*k + 1) = sin(k*nfp*v)
            dt(2*k) = -k*nfp*t(2*k + 1); dt(2*k + 1) = k*nfp*t(2*k)
        end do
    end subroutine trig

    subroutine grid_inputs(sheet, x, phi, current, jnet, u)
        ! u(8, nu, nv) = (R, Z, R_u, Z_u, R_v, Z_v, Phi_u, Phi_v) on the grid
        ! of one period from x(2, nb, nmode) and phi(nph, nmode) (with the
        ! ramps).
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), phi(:, :), current, jnet
        real(dp), intent(out) :: u(:, :, :)
        real(dp) :: gm(sheet%nu, sheet%nmode)
        integer :: c
        do c = 1, 2
            gm = matmul(sheet%smat, x(c, :, :))
            u(c, :, :) = matmul(gm, sheet%t(:, :sheet%nv))
            u(c + 4, :, :) = matmul(gm, sheet%dt(:, :sheet%nv))
            u(c + 2, :, :) = matmul(matmul(sheet%dsmat, x(c, :, :)), sheet%t(:, :sheet%nv))
        end do
        u(7, :, :) = matmul(matmul(sheet%dub, phi), sheet%t(:, :sheet%nv)) + current
        u(8, :, :) = matmul(matmul(sheet%ub, phi), sheet%dt(:, :sheet%nv)) + jnet/(2*pi)
    end subroutine grid_inputs

    pure function toroidal_current(u, jnet) result(u0)
        ! The inputs with Phi = j v/(2 pi) only: the poloidal sheet current that
        ! cancels the coil toroidal field.
        real(dp), intent(in) :: u(:, :, :), jnet
        real(dp) :: u0(size(u, 1), size(u, 2), size(u, 3))
        u0 = u
        u0(7, :, :) = 0
        u0(8, :, :) = jnet/(2*pi)
    end function toroidal_current

    subroutine grid_adjoint(sheet, g, gx, gphi)
        ! Transpose of grid_inputs (without the ramps): g(8, nu, nv) ->
        ! gx(2, nb, nmode), gphi(nph, nmode).
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: g(:, :, :)
        real(dp), intent(out) :: gx(:, :, :), gphi(:, :)
        real(dp) :: tt(sheet%nv, sheet%nmode), dtt(sheet%nv, sheet%nmode)
        integer :: c
        tt = transpose(sheet%t(:, :sheet%nv))
        dtt = transpose(sheet%dt(:, :sheet%nv))
        do c = 1, 2
            gx(c, :, :) = matmul(transpose(sheet%smat), matmul(g(c, :, :), tt) &
                + matmul(g(c + 4, :, :), dtt)) &
                + matmul(transpose(sheet%dsmat), matmul(g(c + 2, :, :), tt))
        end do
        gphi = matmul(transpose(sheet%dub), matmul(g(7, :, :), tt)) &
            + matmul(transpose(sheet%ub), matmul(g(8, :, :), dtt))
    end subroutine grid_adjoint

    subroutine coupling(sheet, u, phi, field, current, jnet, order, e, g, gl, du, dphi, &
            hg, hl)
        ! e += -load + c Psi_t - j Psi_p with its gradient in u (g) and phi
        ! (gl); order 2: the Hessian applied to (du, dphi) into hg, hl.
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: u(:, :, :), phi(:, :), current, jnet
        procedure(sheet_field) :: field
        integer, intent(in) :: order
        real(dp), intent(inout) :: e, g(:, :, :), gl(:, :)
        real(dp), intent(in), optional :: du(:, :, :), dphi(:, :)
        real(dp), intent(inout), optional :: hg(:, :, :), hl(:, :)
        real(dp) :: hu, hv, w, v, ph, out(28), b(3), db(3, 2), d2b(3, 3), gr(6), he(6, 6)
        real(dp) :: gs(sheet%ng_s), ws(sheet%ng_s), ta(sheet%nmode), dta(sheet%nmode)
        real(dp) :: a(4), q(4), s, f(15), g4(4), h4(4, 4), c0(2)
        real(dp) :: pg(sheet%nu, sheet%nv), dpg(sheet%nu, sheet%nv), gl2(sheet%nu, sheet%nv)
        real(dp) :: hl2(sheet%nu, sheet%nv), fv, rp
        integer :: i, j, p, jt, k, l, m, idx(4)
        hu = 1.0_dp/sheet%nu
        hv = 2*pi/(sheet%nfp*sheet%nv)
        w = hu*hv
        fv = -jnet/(2*pi*sheet%o)
        ! phi on the grid (its modes are periodic over a period)
        pg = matmul(phi, sheet%t(:, :sheet%nv))
        if (order == 2) dpg = matmul(dphi, sheet%t(:, :sheet%nv))
        gl2 = 0; hl2 = 0
        ! load over the torus: -w Phi B.N
        do p = 0, sheet%nfp - 1
            do j = 1, sheet%nv
                jt = j + p*sheet%nv
                v = (jt - 1)*hv
                do i = 1, sheet%nu
                    call field(u(1, i, j), u(2, i, j), v, b, db, d2b)
                    call generated_toroidal_sheet_load_jet(u(1, i, j), u(3, i, j), &
                        u(4, i, j), u(5, i, j), u(6, i, j), b(1), db(1, 1), db(1, 2), &
                        d2b(1, 1), d2b(1, 2), d2b(1, 3), b(2), db(2, 1), db(2, 2), &
                        d2b(2, 1), d2b(2, 2), d2b(2, 3), b(3), db(3, 1), db(3, 2), &
                        d2b(3, 1), d2b(3, 2), d2b(3, 3), out(1), out(2), out(3), out(4), &
                        out(5), out(6), out(7), out(8), out(9), out(10), out(11), &
                        out(12), out(13), out(14), out(15), out(16), out(17), out(18), &
                        out(19), out(20), out(21), out(22), out(23), out(24), out(25), &
                        out(26), out(27), out(28))
                    ph = pg(i, j)
                    e = e - w*ph*out(1)
                    gr = out(2:7)
                    g(1:6, i, j) = g(1:6, i, j) - w*ph*gr
                    gl2(i, j) = gl2(i, j) - w*out(1)
                    if (order < 2) cycle
                    k = 7
                    do l = 1, 6
                        do m = l, 6
                            k = k + 1
                            he(l, m) = out(k); he(m, l) = out(k)
                        end do
                    end do
                    hg(1:6, i, j) = hg(1:6, i, j) - w*(ph*matmul(he, du(1:6, i, j)) &
                        + gr*dpg(i, j))
                    hl2(i, j) = hl2(i, j) - w*dot_product(gr, du(1:6, i, j))
                end do
            end do
        end do
        gl = gl + matmul(gl2, transpose(sheet%t(:, :sheet%nv)))
        if (order == 2) hl = hl + matmul(hl2, transpose(sheet%t(:, :sheet%nv)))
        call gauss_unit(sheet%ng_s, gs, ws)
        ! ribbons from L to the toroidal loops at every u_i, averaged over u
        idx = [1, 2, 5, 6]
        do p = 0, sheet%nfp - 1
            do j = 1, sheet%nv
                jt = j + p*sheet%nv
                v = (jt - 1)*hv
                ta = sheet%t(:, jt); dta = sheet%dt(:, jt)
                a = [dot_product(sheet%axis(1, :), ta), dot_product(sheet%axis(2, :), ta), &
                    dot_product(sheet%axis(1, :), dta), dot_product(sheet%axis(2, :), dta)]
                do i = 1, sheet%nu
                    q = u(idx, i, j)
                    do l = 1, sheet%ng_s
                        s = gs(l)
                        call field(a(1) + s*(q(1) - a(1)), a(2) + s*(q(2) - a(2)), v, b, &
                            db, d2b)
                        call ribbon_flux_jet(q(1), q(2), q(3), q(4), s, a(1), a(2), a(3), &
                            a(4), b(1), db(1, 1), db(1, 2), d2b(1, 1), d2b(1, 2), d2b(1, 3), &
                            b(2), db(2, 1), db(2, 2), d2b(2, 1), d2b(2, 2), d2b(2, 3), b(3), &
                            db(3, 1), db(3, 2), d2b(3, 1), d2b(3, 2), d2b(3, 3), f(1), f(2), &
                            f(3), f(4), f(5), f(6), f(7), f(8), f(9), f(10), f(11), f(12), &
                            f(13), f(14), f(15))
                        call add_leaf(current*ws(l)*hv*hu, f, [1, 2, 5, 6], i, j)
                    end do
                end do
            end do
        end do
        ! fans from L(v) to the sections at every v, integrated over v
        do p = 0, sheet%nfp - 1
            do j = 1, sheet%nv
                jt = j + p*sheet%nv
                v = (jt - 1)*hv
                c0 = [dot_product(sheet%axis(1, :), sheet%t(:, jt)), &
                    dot_product(sheet%axis(2, :), sheet%t(:, jt))]
                do i = 1, sheet%nu
                    q = u([1, 2, 3, 4], i, j)
                    do l = 1, sheet%ng_s
                        s = gs(l)
                        call field(c0(1) + s*(q(1) - c0(1)), c0(2) + s*(q(2) - c0(2)), v, &
                            b, db, d2b)
                        ! B_phi beyond F/R (the F/R part is the volume term)
                        rp = c0(1) + s*(q(1) - c0(1))
                        b(2) = b(2) - fv/rp
                        db(2, 1) = db(2, 1) + fv/rp**2
                        d2b(2, 1) = d2b(2, 1) - 2*fv/rp**3
                        call fan_flux_jet(q(1), q(2), q(3), q(4), s, c0(1), c0(2), b(2), &
                            db(2, 1), db(2, 2), d2b(2, 1), d2b(2, 2), d2b(2, 3), f(1), &
                            f(2), f(3), f(4), f(5), f(6), f(7), f(8), f(9), f(10), f(11), &
                            f(12), f(13), f(14), f(15))
                        call add_leaf(-jnet/(2*pi)*ws(l)*hu*hv, f, [1, 2, 3, 4], i, j)
                    end do
                end do
            end do
        end do
    contains
        subroutine add_leaf(wt, f, idx, i, j)
            real(dp), intent(in) :: wt, f(15)
            integer, intent(in) :: idx(4), i, j
            integer :: r, c, k
            e = e + wt*f(1)
            g4 = f(2:5)
            g(idx, i, j) = g(idx, i, j) + wt*g4
            if (order < 2) return
            k = 5
            do r = 1, 4
                do c = r, 4
                    k = k + 1
                    h4(r, c) = f(k); h4(c, r) = f(k)
                end do
            end do
            hg(idx, i, j) = hg(idx, i, j) + wt*matmul(h4, du(idx, i, j))
        end subroutine add_leaf
    end subroutine coupling

    subroutine fb_sheet_terms(sheet, x, field, current, f_vac, mu0, energy, gx, gphi)
        ! L_v at x(2, nb, nmode) and the present phi, with its gradient in x
        ! (same shape) and phi (nph, nmode; zero at the gauge).
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0
        procedure(sheet_field) :: field
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gx(:, :, :), gphi(:, :)
        real(dp) :: u(8, sheet%nu, sheet%nv), g(8, sheet%nu, sheet%nv), gl(sheet%nu, sheet%nmode)
        real(dp) :: g0(8, sheet%nu, sheet%nv), e0
        real(dp) :: e, jnet, ax(2, sheet%nb, sheet%nmode), ap(sheet%nph, sheet%nmode)
        real(dp) :: pn(sheet%nu, sheet%nmode)
        jnet = -2*pi*f_vac*sheet%o
        pn = matmul(sheet%ub, sheet%phi)
        call grid_inputs(sheet, x, sheet%phi, current, jnet, u)
        if (present(gx)) then
            call tst_energy(sheet%grid, u, e, g)
            call tst_energy(sheet%grid, toroidal_current(u, jnet), e0, g0)
            g(1:6, :, :) = g(1:6, :, :) - g0(1:6, :, :)
        else
            call tst_energy(sheet%grid, u, e)
            call tst_energy(sheet%grid, toroidal_current(u, jnet), e0)
            g = 0
        end if
        e = e - e0
        gl = 0
        call coupling(sheet, u, pn, field, current, jnet, 1, e, g, gl)
        energy = -e/mu0
        if (.not. present(gx)) return
        call grid_adjoint(sheet, g, ax, ap)
        gx = -ax/mu0
        if (present(gphi)) then
            gphi = -(ap + matmul(transpose(sheet%ub), gl))/mu0
            gphi(1, 1) = 0
        end if
    end subroutine fb_sheet_terms

    subroutine fb_sheet_hvp(sheet, x, field, current, f_vac, mu0, dx, dphi, hx, hphi)
        ! Hessian of L_v in (x, phi) applied to (dx, dphi).
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0, dx(:, :, :), dphi(:, :)
        procedure(sheet_field) :: field
        real(dp), intent(out) :: hx(:, :, :), hphi(:, :)
        call hvp_e(sheet, x, field, current, f_vac, dx, dphi, hx, hphi)
        hx = -hx/mu0
        hphi = -hphi/mu0
    end subroutine fb_sheet_hvp

    subroutine hvp_e(sheet, x, field, current, f_vac, dx, dphi, hx, hphi)
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, dx(:, :, :), dphi(:, :)
        procedure(sheet_field) :: field
        real(dp), intent(out) :: hx(:, :, :), hphi(:, :)
        real(dp) :: u(8, sheet%nu, sheet%nv), du(8, sheet%nu, sheet%nv)
        real(dp) :: hg(8, sheet%nu, sheet%nv), g(8, sheet%nu, sheet%nv)
        real(dp) :: gl(sheet%nu, sheet%nmode), hl(sheet%nu, sheet%nmode)
        real(dp) :: dp0(sheet%nph, sheet%nmode), pn(sheet%nu, sheet%nmode)
        real(dp) :: dpn(sheet%nu, sheet%nmode)
        real(dp) :: e, jnet
        jnet = -2*pi*f_vac*sheet%o
        pn = matmul(sheet%ub, sheet%phi)
        call grid_inputs(sheet, x, sheet%phi, current, jnet, u)
        dp0 = dphi
        dp0(1, 1) = 0
        dpn = matmul(sheet%ub, dp0)
        call grid_inputs(sheet, dx, dp0, 0.0_dp, 0.0_dp, du)
        call tst_hvp(sheet%grid, u, du, hg)
        block
            real(dp) :: d0(8, sheet%nu, sheet%nv), h0(8, sheet%nu, sheet%nv)
            d0 = du
            d0(7:8, :, :) = 0
            call tst_hvp(sheet%grid, toroidal_current(u, jnet), d0, h0)
            hg(1:6, :, :) = hg(1:6, :, :) - h0(1:6, :, :)
        end block
        g = 0; gl = 0; hl = 0; e = 0
        call coupling(sheet, u, pn, field, current, jnet, 2, e, g, gl, du, dpn, hg, hl)
        call grid_adjoint(sheet, hg, hx, hphi)
        hphi = hphi + matmul(transpose(sheet%ub), hl)
        hphi(1, 1) = 0
    end subroutine hvp_e

    subroutine fb_sheet_exact(sheet, x, field, current, f_vac, mu0, status)
        ! phi stationary at x: e is quadratic in phi with the block M = L^T F L
        ! (F the current block of E_self, L the spectral derivatives); its
        ! Cholesky factor is kept for fb_sheet_reduced_hvp.
        type(fb_sheet_t), intent(inout) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0
        procedure(sheet_field) :: field
        integer, intent(out) :: status
        real(dp), allocatable :: m(:, :), mn(:, :)
        real(dp), allocatable :: rhs(:), gphi(:, :), gx(:, :, :), ph(:)
        real(dp) :: u(8, sheet%nu, sheet%nv), jnet, en
        integer :: n, info, b1, b2
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine dpotrf
        end interface
        n = sheet%nph*sheet%nmode
        jnet = -2*pi*f_vac*sheet%o
        call grid_inputs(sheet, x, sheet%phi, current, jnet, u)
        allocate (mn(sheet%nu*sheet%nmode, sheet%nu*sheet%nmode), m(n, n))
        call tst_current_modes(sheet%grid, u, sheet%du, sheet%t(:, :sheet%nv), &
            sheet%dt(:, :sheet%nv), mn)
        ! to the u modes: blocks ub^T M_nn' ub
        do b2 = 1, sheet%nmode
            do b1 = 1, sheet%nmode
                m((b1 - 1)*sheet%nph + 1:b1*sheet%nph, (b2 - 1)*sheet%nph + 1:b2*sheet%nph) = &
                    matmul(transpose(sheet%ub), matmul(mn((b1 - 1)*sheet%nu + 1:b1*sheet%nu, &
                    (b2 - 1)*sheet%nu + 1:b2*sheet%nu), sheet%ub))
            end do
        end do
        deallocate (mn)
        m = 0.5_dp*(m + transpose(m))
        ! drop the gauge (point 1) and factor
        sheet%chol = m(2:, 2:)
        call dpotrf('L', n - 1, sheet%chol, n - 1, info)
        status = info
        if (info /= 0) return
        allocate (gphi(sheet%nph, sheet%nmode), gx(2, sheet%nb, sheet%nmode))
        call fb_sheet_terms(sheet, x, field, current, f_vac, mu0, en, gx, gphi)
        ! gphi is the gradient of L_v = -e/mu0: d e/d phi = -mu0 gphi
        rhs = reshape(mu0*gphi, [n])
        ph = rhs(2:)
        call chol_solve(sheet, ph)
        sheet%phi = sheet%phi + reshape([0.0_dp, ph], [sheet%nph, sheet%nmode])
    end subroutine fb_sheet_exact

    subroutine chol_solve(sheet, b)
        ! b <- M^-1 b (gauge removed) with the stored factor.
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(inout) :: b(:)
        integer :: info
        interface
            subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(in) :: a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine dpotrs
        end interface
        call dpotrs('L', size(b), 1, sheet%chol, size(b), b, size(b), info)
    end subroutine chol_solve

    subroutine fb_sheet_reduced_hvp(sheet, x, field, current, f_vac, mu0, dx, hx)
        ! Hessian of the reduced L_v (phi eliminated at x by fb_sheet_exact)
        ! applied to dx: H_xx dx - H_xphi M^-1 H_phix dx.
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0, dx(:, :, :)
        procedure(sheet_field) :: field
        real(dp), intent(out) :: hx(:, :, :)
        real(dp) :: hp(sheet%nph, sheet%nmode), zero(sheet%nph, sheet%nmode)
        real(dp) :: flat(sheet%nph*sheet%nmode)
        real(dp) :: hb(2, sheet%nb, sheet%nmode), zx(2, sheet%nb, sheet%nmode)
        real(dp), allocatable :: y(:)
        zero = 0
        call hvp_e(sheet, x, field, current, f_vac, dx, zero, hx, hp)
        flat = reshape(hp, [size(flat)])
        y = -flat(2:)
        call chol_solve(sheet, y)
        zx = 0
        call hvp_e(sheet, x, field, current, f_vac, zx, reshape([0.0_dp, y], &
            [sheet%nph, sheet%nmode]), hb, hp)
        hx = -(hx + hb)/mu0
    end subroutine fb_sheet_reduced_hvp


end module milca_ideal_fb_sheet
