module milca_ideal_fb_sheet
    ! Vacuum of a 3D free boundary on the plasma boundary Gamma_p as a surface
    ! current (theory/derivations/free-boundary-plasma-boundary-sheet-3d.md).
    ! Its contribution to the Lagrangian (per full turn) is L_v = -e/mu0,
    !
    !   e = E_self[Phi] - int Phi B_c . N du dv + c Psi_t - j Psi_p,
    !
    ! Phi = phi + c u + j v/(2 pi) on the boundary chain (u in [0, 1) from
    ! chain node 1, v in [0, 2 pi)), E_self the energy of K = n x grad Phi,
    ! c = mu0 I_phi, j = -2 pi F o (o = 1 for a counter-clockwise chain in
    ! (R, Z)), Psi_t the coil flux through the ribbon from a fixed curve L
    ! inside the plasma to the loop u = 0, and Psi_p the flux through the fan
    ! from L(0) to the section v = 0. The flux through L itself is constant and
    ! dropped. The surface is the chain of P_k panels through the plasma
    ! boundary nodes times the map's real Fourier modes; phi has the same
    ! panels and modes, its value at chain node 1 in mode 1 is the gauge.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_toroidal_sheet_current_3d, only: ts_surface_t, ts_plan_t, &
        ts_field, ts_surface_create, ts_size, ts_energy, ts_load
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gauss_unit
    use milca_sheet_flux_leaves, only: fan_flux_jet, ribbon_flux_jet
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    type, public :: fb_sheet_t
        integer :: nb = 0, nmode = 1, nfp = 1
        ! plasma node of each chain node
        integer, allocatable :: node(:)
        type(ts_surface_t) :: surface
        type(ts_plan_t) :: plan
        ! L: (R, Z) per geometry mode; o: orientation of the chain
        real(dp), allocatable :: axis(:, :)
        real(dp) :: o = 1
        ! Gauss points along s and along u or v for the fluxes
        integer :: ng_s = 8, ng_line = 8
    end type

    public :: fb_sheet_setup, fb_sheet_size, fb_sheet_terms, fb_sheet_exact
    public :: fb_sheet_geometry_index, fb_sheet_phi_index

contains

    subroutine fb_sheet_setup(sheet, node, x, axis, degree, nfp, ntor, status)
        ! node: plasma node of each chain node (nb = size(node), a multiple of
        ! degree); x(2, nb, nmode) their (R, Z) coefficients; axis(2, nmode)
        ! the fixed curve L.
        type(fb_sheet_t), intent(out) :: sheet
        integer, intent(in) :: node(:), degree, nfp, ntor
        real(dp), intent(in) :: x(:, :, :), axis(:, :)
        integer, intent(out) :: status
        integer :: nb, b
        real(dp) :: area, perimeter

        nb = size(node)
        status = 1
        if (mod(nb, degree) /= 0 .or. nb < 3*degree) return
        call ts_surface_create(degree, nb/degree, nfp, ntor, ntor, sheet%surface, status)
        if (status /= 0) return
        status = 1
        if (any(shape(x) /= [2, nb, sheet%surface%nmode])) return
        if (any(shape(axis) /= [2, sheet%surface%nmode])) return
        sheet%nb = nb; sheet%nmode = sheet%surface%nmode; sheet%nfp = nfp
        sheet%node = node
        sheet%axis = axis
        sheet%surface%rz = x
        sheet%surface%phi = 0
        ! toroidal intervals at most 3 poloidal panel lengths long (relative
        ! energy error near 1e-6 on an ellipse, exponential in nv)
        area = 0; perimeter = 0
        do b = 1, nb
            area = area + x(1, b, 1)*x(2, 1 + mod(b, nb), 1) &
                - x(1, 1 + mod(b, nb), 1)*x(2, b, 1)
            perimeter = perimeter + norm2(x(:, 1 + mod(b, nb), 1) - x(:, b, 1))
        end do
        sheet%plan%nv = max(4, 2*ntor, ceiling(2*pi*maxval(x(1, :, 1)) &
            /(3*nfp*perimeter/(nb/degree))))
        sheet%o = sign(1.0_dp, area)
        status = 0
    end subroutine fb_sheet_setup

    pure integer function fb_sheet_size(sheet) result(n)
        ! Geometry coefficients and phi without its gauge.
        type(fb_sheet_t), intent(in) :: sheet
        n = 3*sheet%nb*sheet%nmode - 1
    end function fb_sheet_size

    pure integer function fb_sheet_geometry_index(sheet, c, m, b) result(i)
        ! Component c (R, Z) of mode m at chain node b.
        type(fb_sheet_t), intent(in) :: sheet
        integer, intent(in) :: c, m, b
        i = b + sheet%nb*(m - 1) + (c - 1)*sheet%nb*sheet%nmode
    end function fb_sheet_geometry_index

    pure integer function fb_sheet_phi_index(sheet, m, b) result(i)
        ! phi of mode m at chain node b; 0 for the gauge.
        type(fb_sheet_t), intent(in) :: sheet
        integer, intent(in) :: m, b
        i = 2*sheet%nb*sheet%nmode + b + sheet%nb*(m - 1) - 1
        if (m == 1 .and. b == 1) i = 0
    end function fb_sheet_phi_index

    subroutine fb_sheet_terms(sheet, x, field, current, f_vac, mu0, energy, &
            gradient, hessian)
        ! L_v at the chain coefficients x(2, nb, nmode) and the present phi,
        ! with its gradient and Hessian in (geometry, phi without the gauge);
        ! field gives B_c and its jets, current = mu0 I_phi, f_vac = R B_phi
        ! of the coils.
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0
        procedure(ts_field) :: field
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gradient(:), hessian(:, :)
        type(ts_surface_t) :: s
        real(dp), allocatable :: g(:), h(:, :), gl(:), hl(:, :)
        real(dp) :: e, el, j
        integer, allocatable :: keep(:)
        integer :: n, ng, k
        logical :: second

        s = sheet%surface
        s%rz = x
        s%i_t = current
        j = -2*pi*f_vac*sheet%o
        s%i_p = j
        n = ts_size(s)
        ng = 2*sheet%nb*sheet%nmode
        second = present(hessian)
        if (present(gradient)) then
            allocate (g(n), gl(n))
            if (second) then
                allocate (h(n, n), hl(n, n))
                call ts_energy(s, sheet%plan, e, g, h)
                call ts_load(s, sheet%plan, field, el, gl, hl)
                h = h - hl
            else
                call ts_energy(s, sheet%plan, e, g)
                call ts_load(s, sheet%plan, field, el, gl)
            end if
            g = g - gl
        else
            call ts_energy(s, sheet%plan, e)
            call ts_load(s, sheet%plan, field, el)
        end if
        e = e - el
        call ribbon_flux(sheet, x, field, current, e, g, h)
        call fan_flux(sheet, x, field, -j, e, g, h)
        energy = -e/mu0
        ! drop the gauge, FortFEM index ng + 1
        keep = [(k, k=1, ng), (k, k=ng + 2, n)]
        if (present(gradient)) gradient = -g(keep)/mu0
        if (second) hessian = -h(keep, keep)/mu0
    end subroutine fb_sheet_terms

    subroutine ribbon_flux(sheet, x, field, weight, e, g, h)
        ! e += weight Psi_t with its derivatives in the geometry of chain
        ! node 1 (FortFEM coefficient order in g, h when allocated).
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), weight
        procedure(ts_field) :: field
        real(dp), intent(inout) :: e
        real(dp), allocatable, intent(inout) :: g(:), h(:, :)
        real(dp) :: gs(sheet%ng_s), ws(sheet%ng_s), gv(sheet%ng_line), wv(sheet%ng_line)
        real(dp) :: t(sheet%nmode), dt(sheet%nmode), v, w, q(4), a(4), s, out(15)
        real(dp) :: b(3), db(3, 2), d2b(3, 3), jac(4, 2*sheet%nmode)
        integer :: k, i, l, m, idx(2*sheet%nmode), nv

        call gauss_unit(sheet%ng_s, gs, ws)
        call gauss_unit(sheet%ng_line, gv, wv)
        do m = 1, sheet%nmode
            idx(m) = fb_sheet_geometry_index(sheet, 1, m, 1)
            idx(sheet%nmode + m) = fb_sheet_geometry_index(sheet, 2, m, 1)
        end do
        nv = sheet%nfp*sheet%plan%nv
        do k = 1, nv
            do i = 1, sheet%ng_line
                v = 2*pi*(k - 1 + gv(i))/nv
                call trig(sheet%nmode, sheet%nfp, v, t, dt)
                ! boundary loop (r, z, r_v, z_v) and L (r, z, r_v, z_v)
                q = [dot_product(x(1, 1, :), t), dot_product(x(2, 1, :), t), &
                    dot_product(x(1, 1, :), dt), dot_product(x(2, 1, :), dt)]
                a = [dot_product(sheet%axis(1, :), t), dot_product(sheet%axis(2, :), t), &
                    dot_product(sheet%axis(1, :), dt), dot_product(sheet%axis(2, :), dt)]
                jac = 0
                jac(1, :sheet%nmode) = t; jac(3, :sheet%nmode) = dt
                jac(2, sheet%nmode + 1:) = t; jac(4, sheet%nmode + 1:) = dt
                do l = 1, sheet%ng_s
                    s = gs(l)
                    w = weight*ws(l)*wv(i)*2*pi/nv
                    call field(a(1) + s*(q(1) - a(1)), a(2) + s*(q(2) - a(2)), v, b, db, d2b)
                    call ribbon_flux_jet(q(1), q(2), q(3), q(4), s, a(1), a(2), a(3), &
                        a(4), b(1), db(1, 1), db(1, 2), d2b(1, 1), d2b(1, 2), d2b(1, 3), &
                        b(2), db(2, 1), db(2, 2), d2b(2, 1), d2b(2, 2), d2b(2, 3), b(3), &
                        db(3, 1), db(3, 2), d2b(3, 1), d2b(3, 2), d2b(3, 3), out(1), &
                        out(2), out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
                        out(10), out(11), out(12), out(13), out(14), out(15))
                    call accumulate(w, out, jac, idx, e, g, h)
                end do
            end do
        end do
    end subroutine ribbon_flux

    subroutine fan_flux(sheet, x, field, weight, e, g, h)
        ! e += weight Psi_p with its derivatives in the geometry at v = 0.
        type(fb_sheet_t), intent(in) :: sheet
        real(dp), intent(in) :: x(:, :, :), weight
        procedure(ts_field) :: field
        real(dp), intent(inout) :: e
        real(dp), allocatable, intent(inout) :: g(:), h(:, :)
        type(ts_surface_t) :: sf
        real(dp) :: gs(sheet%ng_s), ws(sheet%ng_s), gu(sheet%ng_line), wu(sheet%ng_line)
        real(dp) :: t(sheet%nmode), dt(sheet%nmode), c(2), q(4), s, w, out(15)
        real(dp) :: b(3), db(3, 2), d2b(3, 3), xi
        real(dp) :: nb_(sheet%surface%degree + 1), dnb(sheet%surface%degree + 1)
        real(dp), allocatable :: jac(:, :)
        integer, allocatable :: idx(:)
        integer :: p, i, l, a, m, node, k, deg, ncol

        sf = sheet%surface
        deg = sf%degree
        call gauss_unit(sheet%ng_s, gs, ws)
        call gauss_unit(sheet%ng_line, gu, wu)
        call trig(sheet%nmode, sheet%nfp, 0.0_dp, t, dt)
        c = matmul(sheet%axis, t)
        ncol = 2*(deg + 1)*sheet%nmode
        allocate (jac(4, ncol), idx(ncol))
        do p = 1, sf%n_panels
            do i = 1, sheet%ng_line
                xi = gu(i)
                call lagrange(sf, xi, nb_, dnb)
                q = 0; jac = 0; k = 0
                do a = 1, deg + 1
                    node = (p - 1)*deg + a
                    if (node > sheet%nb) node = node - sheet%nb
                    do m = 1, sheet%nmode
                        q(1) = q(1) + nb_(a)*t(m)*x(1, node, m)
                        q(2) = q(2) + nb_(a)*t(m)*x(2, node, m)
                        q(3) = q(3) + dnb(a)*t(m)*x(1, node, m)
                        q(4) = q(4) + dnb(a)*t(m)*x(2, node, m)
                        k = k + 1
                        idx(k) = fb_sheet_geometry_index(sheet, 1, m, node)
                        jac(1, k) = nb_(a)*t(m); jac(3, k) = dnb(a)*t(m)
                        k = k + 1
                        idx(k) = fb_sheet_geometry_index(sheet, 2, m, node)
                        jac(2, k) = nb_(a)*t(m); jac(4, k) = dnb(a)*t(m)
                    end do
                end do
                do l = 1, sheet%ng_s
                    s = gs(l)
                    w = weight*ws(l)*wu(i)
                    call field(c(1) + s*(q(1) - c(1)), c(2) + s*(q(2) - c(2)), 0.0_dp, &
                        b, db, d2b)
                    call fan_flux_jet(q(1), q(2), q(3), q(4), s, c(1), c(2), b(2), &
                        db(2, 1), db(2, 2), d2b(2, 1), d2b(2, 2), d2b(2, 3), out(1), &
                        out(2), out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
                        out(10), out(11), out(12), out(13), out(14), out(15))
                    call accumulate(w, out, jac, idx, e, g, h)
                end do
            end do
        end do
    end subroutine fan_flux

    subroutine accumulate(w, out, jac, idx, e, g, h)
        ! e += w f, g += w J' grad f, h += w J' hess f J for a leaf of the four
        ! inputs (r, z, a, b) = jac c(idx).
        real(dp), intent(in) :: w, out(15), jac(:, :)
        integer, intent(in) :: idx(:)
        real(dp), intent(inout) :: e
        real(dp), allocatable, intent(inout) :: g(:), h(:, :)
        real(dp) :: he(4, 4), hj(4, size(idx))
        integer :: i, j, k

        e = e + w*out(1)
        if (.not. allocated(g)) return
        g(idx) = g(idx) + w*matmul(out(2:5), jac)
        if (.not. allocated(h)) return
        k = 5
        do i = 1, 4
            do j = i, 4
                k = k + 1
                he(i, j) = out(k); he(j, i) = out(k)
            end do
        end do
        hj = matmul(he, jac)
        h(idx, idx) = h(idx, idx) + w*matmul(transpose(jac), hj)
    end subroutine accumulate

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

    pure subroutine lagrange(s, xi, n, dn)
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: n(:), dn(:)
        integer :: i, j
        n = 0; dn = 0
        do i = 1, s%degree + 1
            do j = s%degree + 1, 1, -1
                n(i) = n(i)*xi + s%lagrange(i, j)
                if (j > 1) dn(i) = dn(i)*xi + (j - 1)*s%lagrange(i, j)
            end do
        end do
    end subroutine lagrange

    subroutine fb_sheet_exact(sheet, x, field, current, f_vac, mu0, status)
        ! phi stationary at the chain coefficients x: L_v is quadratic and
        ! concave in phi, one Newton step on the phi block.
        type(fb_sheet_t), intent(inout) :: sheet
        real(dp), intent(in) :: x(:, :, :), current, f_vac, mu0
        procedure(ts_field) :: field
        integer, intent(out) :: status
        real(dp), allocatable :: g(:), h(:, :), a(:, :), rhs(:), phi(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: e
        integer :: n, ng, np, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        n = fb_sheet_size(sheet)
        ng = 2*sheet%nb*sheet%nmode
        np = n - ng
        allocate (g(n), h(n, n), ipiv(np))
        call fb_sheet_terms(sheet, x, field, current, f_vac, mu0, e, g, h)
        a = h(ng + 1:, ng + 1:)
        rhs = -g(ng + 1:)
        call dgesv(np, 1, a, np, ipiv, rhs, np, info)
        status = info
        if (info /= 0) return
        phi = reshape(sheet%surface%phi, [sheet%nb*sheet%nmode])
        phi(2:) = phi(2:) + rhs
        sheet%surface%phi = reshape(phi, [sheet%nb, sheet%nmode])
    end subroutine fb_sheet_exact

end module milca_ideal_fb_sheet
