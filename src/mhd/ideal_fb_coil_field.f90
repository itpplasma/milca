module milca_ideal_fb_coil_field
    ! Smooth harmonic representation of a coil field inside a coil-free
    ! region, for the 3D vacuum on the plasma boundary
    ! (theory/derivations/free-boundary-plasma-boundary-sheet-3d.md):
    !
    !   B = F grad(phi) + grad(chi),
    !   chi = sum_j sum_k [ac(j,k) C_jk(x) + as(j,k) S_jk(x)],
    !   C_jk = sum_d cos(n phi_d)/|x - y_jd| dphi/(4 pi),  n = k nfp,
    !
    ! with S_jk the same with sin(n phi_d), y_jd = (R_j cos phi_d, R_j sin
    ! phi_d, Z_j) the nd points of ring j at phi_d = 2 pi (d - 1)/nd. Every
    ! term is harmonic, so chi is exactly harmonic where no ring lies, and its
    ! jets to third order come from the generated leaf point_source_jet. The
    ! ring sum approximates the toroidal-mode kernel with an error of order
    ! exp(-nd delta/R), delta the distance of the rings from the region; nd is
    ! chosen to make it negligible. The coefficients are a least-squares fit,
    ! harmonic by harmonic, to the field at points of the region (mgrid
    ! nodes); the misfit is reported.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_coil_field_leaves, only: point_source_jet
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    type, public :: fb_coil_field_t
        integer :: nfp = 1, kmax = 0, nd = 0
        real(dp) :: f_vac = 0, misfit = 0
        real(dp), allocatable :: ring(:, :), ac(:, :), as(:, :)
    end type

    public :: fb_coil_field_fit, fb_coil_field_eval, fb_coil_field_nd

contains

    pure function fb_coil_field_nd(rings, points, nfp) result(nd)
        ! Points per ring: nd delta/R_max >= 40 over the region, a multiple of
        ! nfp, at least 64.
        real(dp), intent(in) :: rings(:, :), points(:, :)
        integer, intent(in) :: nfp
        integer :: nd, i, j
        real(dp) :: delta, rmax

        delta = huge(1.0_dp)
        rmax = max(maxval(rings(1, :)), maxval(points(1, :)))
        do j = 1, size(rings, 2)
            do i = 1, size(points, 2)
                delta = min(delta, norm2(rings(:, j) - points(:, i)))
            end do
        end do
        nd = max(64, ceiling(40*rmax/max(delta, 1.0e-6_dp)))
        nd = min(nd, 8192)
        nd = nfp*((nd + nfp - 1)/nfp)
    end function fb_coil_field_nd

    subroutine fb_coil_field_fit(points, bcos, bsin, nfp, kmax, rings, nd, field, &
            status)
        ! points(2, np) = (R, Z); bcos(3, 0:kmax, np), bsin: the cylindrical
        ! components (B_R, B_phi, B_Z) of the harmonics cos/sin(n phi), n =
        ! k nfp, at the points. rings(2, nr): ring positions outside the
        ! region.
        real(dp), intent(in) :: points(:, :), bcos(:, 0:, :), bsin(:, 0:, :)
        integer, intent(in) :: nfp, kmax, nd
        real(dp), intent(in) :: rings(:, :)
        type(fb_coil_field_t), intent(out) :: field
        integer, intent(out) :: status
        real(dp), allocatable :: amp(:, :, :, :), a(:, :), b(:), work(:), sv(:)
        real(dp) :: num, den, rcond
        integer :: np, nr, k, j, q, m, rank, info, lwork, kind
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
        np = size(points, 2)
        nr = size(rings, 2)
        if (np < 1 .or. nr < 1 .or. kmax < 0 .or. nd < 1) return
        field%nfp = nfp; field%kmax = kmax; field%nd = nd
        field%ring = rings
        allocate (field%ac(nr, 0:kmax), field%as(nr, 0:kmax))
        field%ac = 0; field%as = 0
        ! F: the n = 0 part of R B_phi (grad chi has no n = 0 phi component)
        field%f_vac = sum(points(1, :)*bcos(2, 0, :))/np
        ! amp(c, kind, j, q): cylindrical harmonic amplitudes at point q of the
        ! cos (kind 1) and sin (kind 2) basis function of ring j, mode k
        allocate (amp(3, 2, nr, np))
        rcond = 1.0e-13_dp
        num = 0; den = 0
        do k = 0, kmax
            call basis_amplitudes(field, k, points, amp)
            do kind = 1, merge(1, 2, k == 0)
                ! cos basis: (B_R, B_Z) cos, B_phi sin; sin basis: the reverse
                m = 3*np
                if (k == 0) m = 2*np
                allocate (a(m, nr), b(max(m, nr)))
                do q = 1, np
                    do j = 1, nr
                        a(q, j) = amp(1, kind, j, q)
                        a(np + q, j) = amp(3, kind, j, q)
                        if (k > 0) a(2*np + q, j) = amp(2, kind, j, q)
                    end do
                    if (kind == 1) then
                        b(q) = bcos(1, k, q); b(np + q) = bcos(3, k, q)
                        if (k > 0) b(2*np + q) = bsin(2, k, q)
                    else
                        b(q) = bsin(1, k, q); b(np + q) = bsin(3, k, q)
                        b(2*np + q) = bcos(2, k, q)
                    end if
                end do
                den = den + sum(b(1:m)**2)
                allocate (sv(nr), work(1))
                call dgelss(m, nr, 1, a, m, b, max(m, nr), sv, rcond, rank, work, -1, &
                    info)
                lwork = int(work(1))
                deallocate (work)
                allocate (work(lwork))
                call dgelss(m, nr, 1, a, m, b, max(m, nr), sv, rcond, rank, work, &
                    lwork, info)
                if (info /= 0) return
                if (kind == 1) then
                    field%ac(:, k) = b(1:nr)
                else
                    field%as(:, k) = b(1:nr)
                end if
                if (m > nr) num = num + sum(b(nr + 1:m)**2)
                deallocate (a, b, sv, work)
            end do
        end do
        field%misfit = sqrt(num/max(den, tiny(1.0_dp)))
        status = 0
    end subroutine fb_coil_field_fit

    subroutine basis_amplitudes(field, k, points, amp)
        ! Harmonic amplitudes of the ring basis functions of mode k: the cos
        ! function gives (B_R, B_Z) along cos(n phi) and B_phi along
        ! sin(n phi), read at phi = 0 and phi = pi/(2n); the sin function the
        ! reverse. amp(c, kind, j, q), c = (R, phi, Z).
        type(fb_coil_field_t), intent(in) :: field
        integer, intent(in) :: k
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: amp(:, :, :, :)
        real(dp) :: phi0, phi1, gc0(3), gc1(3), gs0(3), gs1(3)
        integer :: q, j, n

        n = k*field%nfp
        phi0 = 0
        phi1 = 0
        if (n > 0) phi1 = pi/(2*n)
        amp = 0
        !$omp parallel do default(shared) private(q, j, gc0, gc1, gs0, gs1) &
        !$omp& schedule(dynamic, 4)
        do q = 1, size(points, 2)
            do j = 1, size(field%ring, 2)
                call ring_gradient(field%ring(:, j), field%nd, n, points(:, q), &
                    phi0, gc0, gs0)
                amp(1, 1, j, q) = gc0(1)
                amp(3, 1, j, q) = gc0(3)
                amp(2, 2, j, q) = gs0(2)
                if (n > 0) then
                    call ring_gradient(field%ring(:, j), field%nd, n, points(:, q), &
                        phi1, gc1, gs1)
                    amp(2, 1, j, q) = gc1(2)
                    amp(1, 2, j, q) = gs1(1)
                    amp(3, 2, j, q) = gs1(3)
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine basis_amplitudes

    subroutine ring_gradient(ring, nd, n, rz, phi, gc, gs)
        ! Cylindrical gradient at (R, Z, phi) of the cos and sin ring basis
        ! functions of harmonic n.
        real(dp), intent(in) :: ring(2), rz(2), phi
        integer, intent(in) :: nd, n
        real(dp), intent(out) :: gc(3), gs(3)
        real(dp) :: x(3), y(3), o(20), ang, er(3), ep(3), gcart(3)
        integer :: d

        x = [rz(1)*cos(phi), rz(1)*sin(phi), rz(2)]
        gc = 0; gs = 0
        er = [cos(phi), sin(phi), 0.0_dp]
        ep = [-sin(phi), cos(phi), 0.0_dp]
        do d = 1, nd
            ang = 2*pi*(d - 1)/nd
            y = [ring(1)*cos(ang), ring(1)*sin(ang), ring(2)]
            call point_source_jet(x(1) - y(1), x(2) - y(2), x(3) - y(3), o(1), o(2), &
                o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), o(11), o(12), &
                o(13), o(14), o(15), o(16), o(17), o(18), o(19), o(20))
            gcart = o(2:4)/(2*nd)
            gc = gc + cos(n*ang)*[dot_product(er, gcart), dot_product(ep, gcart), &
                gcart(3)]
            gs = gs + sin(n*ang)*[dot_product(er, gcart), dot_product(ep, gcart), &
                gcart(3)]
        end do
    end subroutine ring_gradient

    subroutine fb_coil_field_eval(field, r, z, v, b, db, d2b)
        ! (B_R, B_phi, B_Z) at (R, Z, phi = v), their derivatives in (R, Z)
        ! (db(:, 1:2)) and second derivatives (RR, RZ, ZZ) at fixed v.
        type(fb_coil_field_t), intent(in) :: field
        real(dp), intent(in) :: r, z, v
        real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        real(dp) :: g(3), h(3, 3), t(3, 3, 3)
        real(dp) :: e(3, 3)

        call chi_jet(field, r, z, v, g, h, t)
        ! e(:, 1) = e_R, e(:, 2) = e_phi, e(:, 3) = e_Z
        e(:, 1) = [cos(v), sin(v), 0.0_dp]
        e(:, 2) = [-sin(v), cos(v), 0.0_dp]
        e(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
        call project(g, h, t, e, b, db, d2b)
        b(2) = b(2) + field%f_vac/r
        db(2, 1) = db(2, 1) - field%f_vac/r**2
        d2b(2, 1) = d2b(2, 1) + 2*field%f_vac/r**3
    end subroutine fb_coil_field_eval

    pure subroutine project(g, h, t, e, b, db, d2b)
        ! Components c of grad chi along e(:, c), differentiated along e_R
        ! (e(:, 1)) and e_Z (e(:, 3)).
        real(dp), intent(in) :: g(3), h(3, 3), t(3, 3, 3), e(3, 3)
        real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        integer, parameter :: dir(2) = [1, 3]
        integer :: c, i, j, p, q, s

        do c = 1, 3
            b(c) = dot_product(e(:, c), g)
            do i = 1, 2
                db(c, i) = dot_product(e(:, c), matmul(h, e(:, dir(i))))
            end do
            s = 0
            do i = 1, 2
                do j = i, 2
                    s = s + 1
                    d2b(c, s) = 0
                    do p = 1, 3
                        do q = 1, 3
                            d2b(c, s) = d2b(c, s) + e(p, c)*sum(t(p, q, :) &
                                *e(:, dir(j)))*e(q, dir(i))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine project

    subroutine chi_jet(field, r, z, v, g, h, t)
        ! Gradient, Hessian and third derivatives (Cartesian) of chi at
        ! (R cos v, R sin v, Z).
        type(fb_coil_field_t), intent(in) :: field
        real(dp), intent(in) :: r, z, v
        real(dp), intent(out) :: g(3), h(3, 3), t(3, 3, 3)
        real(dp) :: x(3), y(3), wc, ws, w, phi, o(20)
        integer :: j, d, k, n

        x = [r*cos(v), r*sin(v), z]
        g = 0; h = 0; t = 0
        do j = 1, size(field%ring, 2)
            do d = 1, field%nd
                phi = 2*pi*(d - 1)/field%nd
                y = [field%ring(1, j)*cos(phi), field%ring(1, j)*sin(phi), &
                    field%ring(2, j)]
                w = 0
                do k = 0, field%kmax
                    n = k*field%nfp
                    wc = cos(n*phi); ws = sin(n*phi)
                    w = w + field%ac(j, k)*wc + field%as(j, k)*ws
                end do
                if (w == 0) cycle
                w = w/(2*field%nd)
                call point_source_jet(x(1) - y(1), x(2) - y(2), x(3) - y(3), o(1), &
                    o(2), o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), o(11), &
                    o(12), o(13), o(14), o(15), o(16), o(17), o(18), o(19), o(20))
                call accumulate(w, o, g, h, t)
            end do
        end do
    end subroutine chi_jet

    pure subroutine accumulate(w, o, g, h, t)
        ! o = (v, g(3), Hessian 6, third 10) in the leaf's order; the
        ! symmetric Hessian and third-derivative tensors are filled.
        real(dp), intent(in) :: w, o(20)
        real(dp), intent(inout) :: g(3), h(3, 3), t(3, 3, 3)
        integer, parameter :: hi(2, 6) = reshape([1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], &
            [2, 6])
        integer, parameter :: ti(3, 10) = reshape([1, 1, 1, 1, 1, 2, 1, 1, 3, 1, 2, 2, &
            1, 2, 3, 1, 3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3], [3, 10])
        integer :: m, a, b, c, i, j, l

        g = g + w*o(2:4)
        do m = 1, 6
            a = hi(1, m); b = hi(2, m)
            h(a, b) = h(a, b) + w*o(4 + m)
            if (a /= b) h(b, a) = h(b, a) + w*o(4 + m)
        end do
        do i = 1, 3
            do j = 1, 3
                do l = 1, 3
                    ! sorted index triple of (i, j, l)
                    a = min(i, j, l); c = max(i, j, l); b = i + j + l - a - c
                    do m = 1, 10
                        if (ti(1, m) == a .and. ti(2, m) == b .and. ti(3, m) == c) &
                            t(i, j, l) = t(i, j, l) + w*o(10 + m)
                    end do
                end do
            end do
        end do
    end subroutine accumulate

end module milca_ideal_fb_coil_field
