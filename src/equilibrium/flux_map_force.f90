module milca_flux_map_force
    !! Force balance of a flux-map state, independent of the energy that is
    !! minimised (docs/flux-map.md, "Force-balance residual"): with
    !!
    !!   B = -[tau(s) grad x X grad y + psi'(s) grad phi X grad s],
    !!
    !! the current J = curl B is taken from the second derivatives of the
    !! map u = (x, y) itself (generated leaf flux_map_force_3d), not from the
    !! gradient of W, and F = J X B - grad p is integrated over the volume
    !! and per band of the label s. F . B = 0 holds identically, so F has a
    !! part along grad s (the radial balance) and a part along B X grad s
    !! (which is |B| |J . grad s|/|grad s|: current crossing the surfaces).
    !!
    !! The measure takes the 2-jet of u in the cylindrical coordinates
    !! (R, phi, Z), whatever produced it: flux_map3d_force evaluates the
    !! discrete state of milca_flux_map_3d cell by cell (second derivatives
    !! exist inside the cells; across the C0 block interfaces B_tangential
    !! may jump, which this volume measure does not see), and a state given
    !! in other variables (a VMEC state) is added through
    !! flux_map_jet_to_cylindrical and flux_map_force_add.
    !!
    !! 2-jets are stored as (f, f_1, f_2, f_3, f_11, f_12, f_13, f_22, f_23,
    !! f_33) in three parameters.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map, only: flux_map_profiles
    use milca_flux_map_3d, only: flux_map3d_t
    use milca_flux_map_leaves, only: flux_map_force_3d
    implicit none
    private

    integer, parameter, public :: FORCE_VOLUME = 1, FORCE_F = 2, FORCE_GRAD_P = 3, &
        FORCE_JXB = 4, FORCE_GRAD_B2 = 5, FORCE_J = 6, FORCE_B = 7, FORCE_F_RADIAL = 8, &
        FORCE_F_BINORMAL = 9

    type, public :: flux_map_force_t
        !! sums(q, b): volume integrals over the band b of s (b = 0: the
        !! whole volume; bands of equal width in s, s >= 1 in the last) of 1
        !! and of the squares of |F|, |grad p|, |J X B|, |grad(B^2/2)|, |J|,
        !! |B|, F . grad s/|grad s| and F . (B X grad s)/|B X grad s|.
        real(dp), allocatable :: sums(:, :)
        !! range of grad phi . (grad x X grad y) over the points (one sign:
        !! the map does not fold and B_phi keeps its sign) and their number
        real(dp) :: jacobian_min = huge(1.0_dp), jacobian_max = -huge(1.0_dp)
        integer :: points = 0
    end type flux_map_force_t

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    ! position of the second derivative (a, b) in a 2-jet
    integer, parameter :: second(3, 3) = reshape([5, 6, 7, 6, 8, 9, 7, 9, 10], [3, 3])

    public :: flux_map_force_reset, flux_map_force_add, flux_map_force_merge, &
              flux_map_jet_to_cylindrical, flux_map_jet_product, flux_map_jet_cos_sin, &
              flux_map_jet_quotient, flux_map3d_force, flux_map3d_force_point, flux_map_force_norm, &
              flux_map_force_print

contains

    pure subroutine flux_map_force_reset(acc, bands)
        type(flux_map_force_t), intent(out) :: acc
        integer, intent(in) :: bands

        allocate (acc%sums(9, 0:bands), source=0.0_dp)
    end subroutine flux_map_force_reset

    pure subroutine flux_map_force_merge(acc, part)
        type(flux_map_force_t), intent(inout) :: acc
        type(flux_map_force_t), intent(in) :: part

        acc%sums = acc%sums + part%sums
        acc%jacobian_min = min(acc%jacobian_min, part%jacobian_min)
        acc%jacobian_max = max(acc%jacobian_max, part%jacobian_max)
        acc%points = acc%points + part%points
    end subroutine flux_map_force_merge

    subroutine flux_map_force_add(acc, weight, r, xj, yj, profiles, field)
        !! One point: weight = its volume, xj, yj the 2-jets of x and y in
        !! (R, phi, Z) (coordinate derivatives).
        type(flux_map_force_t), intent(inout) :: acc
        real(dp), intent(in) :: weight, r, xj(10), yj(10)
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(out), optional :: field(3)
        real(dp) :: s, tau(3), gp(3), pr(3), b(3), j(3), f(3), m(3), jac, gs(3), jxb(3), &
                    bn(3), q(9)
        integer :: band

        s = xj(1)**2 + yj(1)**2
        call profiles(s, tau, gp, pr)
        call flux_map_force_3d(r, xj(1), yj(1), xj(2), xj(3), xj(4), xj(5), xj(6), xj(7), xj(8), &
            xj(9), xj(10), yj(2), yj(3), yj(4), yj(5), yj(6), yj(7), yj(8), yj(9), yj(10), &
            tau(1), tau(2), gp(1), gp(2), pr(2), b(1), b(2), b(3), j(1), j(2), j(3), f(1), f(2), &
            f(3), m(1), m(2), m(3), jac)
        if (present(field)) field = b
        gs = 2*[xj(1)*xj(2) + yj(1)*yj(2), (xj(1)*xj(3) + yj(1)*yj(3))/r, &
                xj(1)*xj(4) + yj(1)*yj(4)]
        jxb = f + pr(2)*gs
        bn = [b(2)*gs(3) - b(3)*gs(2), b(3)*gs(1) - b(1)*gs(3), b(1)*gs(2) - b(2)*gs(1)]
        q = [1.0_dp, sum(f**2), pr(2)**2*sum(gs**2), sum(jxb**2), sum(m**2), sum(j**2), &
             sum(b**2), dot_product(f, gs)**2/max(sum(gs**2), tiny(1.0_dp)), &
             dot_product(f, bn)**2/max(sum(bn**2), tiny(1.0_dp))]
        band = min(int(s*ubound(acc%sums, 2)) + 1, ubound(acc%sums, 2))
        acc%sums(:, 0) = acc%sums(:, 0) + weight*q
        acc%sums(:, band) = acc%sums(:, band) + weight*q
        acc%jacobian_min = min(acc%jacobian_min, jac)
        acc%jacobian_max = max(acc%jacobian_max, jac)
        acc%points = acc%points + 1
    end subroutine flux_map_force_add

    pure real(dp) function flux_map_force_norm(acc, quantity, band) result(v)
        !! Root mean square of a quantity over a band (0: the volume).
        type(flux_map_force_t), intent(in) :: acc
        integer, intent(in) :: quantity, band

        v = sqrt(acc%sums(quantity, band)/max(acc%sums(FORCE_VOLUME, band), tiny(1.0_dp)))
    end function flux_map_force_norm

    subroutine flux_map_force_print(acc, label)
        !! The residual relative to |grad(B^2/2)| and to |grad p| (ratios of
        !! L2 norms over the volume, and over the bands outside s = 0.1),
        !! then per band the root mean squares.
        type(flux_map_force_t), intent(in) :: acc
        character(*), intent(in) :: label
        real(dp) :: f, scale_b, scale_p, outer(9)
        integer :: b, n, first

        n = ubound(acc%sums, 2)
        f = flux_map_force_norm(acc, FORCE_F, 0)
        scale_b = flux_map_force_norm(acc, FORCE_GRAD_B2, 0)
        scale_p = flux_map_force_norm(acc, FORCE_GRAD_P, 0)
        first = ceiling(0.1_dp*n - 1e-9_dp) + 1
        outer = sum(acc%sums(:, first:n), 2)
        print '(3a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)', '  force balance (', label, &
            '): |F|/|grad B^2/2|', f/scale_b, '  |F|/|grad p|', ratio(f, scale_p), &
            '  radial part', flux_map_force_norm(acc, FORCE_F_RADIAL, 0)/scale_b, &
            '  binormal part', flux_map_force_norm(acc, FORCE_F_BINORMAL, 0)/scale_b, &
            '  |J|/|B| (1/m)', &
            flux_map_force_norm(acc, FORCE_J, 0)/flux_map_force_norm(acc, FORCE_B, 0), &
            '  volume', acc%sums(FORCE_VOLUME, 0)
        print '(a,es10.3,a,es10.3,a,i0,a,2es11.3)', '    outside s = 0.1: |F|/|grad B^2/2|', &
            ratio(sqrt(outer(FORCE_F)), sqrt(outer(FORCE_GRAD_B2))), '  |F|/|grad p|', &
            ratio(sqrt(outer(FORCE_F)), sqrt(outer(FORCE_GRAD_P))), '  points ', acc%points, &
            '  range of grad phi . (grad x X grad y)', acc%jacobian_min, acc%jacobian_max
        print '(a)', '    band of s      rms |F|   |grad p|     |J X B|  |grad B^2/2|  ' &
            //'|F|/|grad B^2/2| (volume scale)'
        do b = 1, n
            if (acc%sums(FORCE_VOLUME, b) <= 0) cycle
            print '(4x,f5.2,a,f5.2,5es12.3)', (b - 1)/real(n, dp), ' -', b/real(n, dp), &
                flux_map_force_norm(acc, FORCE_F, b), flux_map_force_norm(acc, FORCE_GRAD_P, b), &
                flux_map_force_norm(acc, FORCE_JXB, b), flux_map_force_norm(acc, FORCE_GRAD_B2, b), &
                flux_map_force_norm(acc, FORCE_F, b)/scale_b
        end do

    contains

        pure real(dp) function ratio(a, b)
            real(dp), intent(in) :: a, b
            ratio = 0
            if (b > 0) ratio = a/b
        end function ratio
    end subroutine flux_map_force_print

    pure function flux_map_jet_product(a, b) result(c)
        real(dp), intent(in) :: a(10), b(10)
        real(dp) :: c(10)
        integer :: i, j

        c(1) = a(1)*b(1)
        c(2:4) = a(2:4)*b(1) + a(1)*b(2:4)
        do j = 1, 3
            do i = 1, j
                c(second(i, j)) = a(second(i, j))*b(1) + a(1 + i)*b(1 + j) + a(1 + j)*b(1 + i) &
                                  + a(1)*b(second(i, j))
            end do
        end do
    end function flux_map_jet_product

    pure function flux_map_jet_quotient(a, b) result(c)
        real(dp), intent(in) :: a(10), b(10)
        real(dp) :: c(10), inv(10), v
        integer :: i, j

        v = 1/b(1)
        inv(1) = v
        inv(2:4) = -v**2*b(2:4)
        do j = 1, 3
            do i = 1, j
                inv(second(i, j)) = 2*v**3*b(1 + i)*b(1 + j) - v**2*b(second(i, j))
            end do
        end do
        c = flux_map_jet_product(a, inv)
    end function flux_map_jet_quotient

    pure subroutine flux_map_jet_cos_sin(a, c, s)
        !! 2-jets of cos and sin of the function with 2-jet a.
        real(dp), intent(in) :: a(10)
        real(dp), intent(out) :: c(10), s(10)
        real(dp) :: c0, s0
        integer :: i, j

        c0 = cos(a(1))
        s0 = sin(a(1))
        c(1) = c0
        s(1) = s0
        c(2:4) = -s0*a(2:4)
        s(2:4) = c0*a(2:4)
        do j = 1, 3
            do i = 1, j
                c(second(i, j)) = -c0*a(1 + i)*a(1 + j) - s0*a(second(i, j))
                s(second(i, j)) = -s0*a(1 + i)*a(1 + j) + c0*a(second(i, j))
            end do
        end do
    end subroutine flux_map_jet_cos_sin

    pure subroutine flux_map_jet_to_cylindrical(rj, zj, phi_slope, fq, fc)
        !! 2-jet fc of a function in (R, phi, Z) from its 2-jet fq in the
        !! parameters q, with R(q), Z(q) given by their 2-jets and
        !! phi = phi_slope q_3: grad_q f = A^T grad f and
        !! H_q f = A^T (H f) A + f_R H_q R + f_Z H_q Z, A = d(R, phi, Z)/dq.
        real(dp), intent(in) :: rj(10), zj(10), phi_slope, fq(10)
        real(dp), intent(out) :: fc(10)
        real(dp) :: a(3, 3), ai(3, 3), det, g(3), h(3, 3)
        integer :: i, j

        a(1, :) = rj(2:4)
        a(2, :) = [0.0_dp, 0.0_dp, phi_slope]
        a(3, :) = zj(2:4)
        det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) - a(1, 2)*(a(2, 1)*a(3, 3) - &
              a(2, 3)*a(3, 1)) + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
        ai(1, :) = [a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2), a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3), &
                    a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)]
        ai(2, :) = [a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3), a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1), &
                    a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)]
        ai(3, :) = [a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1), a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2), &
                    a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)]
        ai = ai/det
        g = matmul(transpose(ai), fq(2:4))
        do j = 1, 3
            do i = 1, 3
                h(i, j) = fq(second(i, j)) - g(1)*rj(second(i, j)) - g(3)*zj(second(i, j))
            end do
        end do
        h = matmul(transpose(ai), matmul(h, ai))
        fc(1) = fq(1)
        fc(2:4) = g
        do j = 1, 3
            do i = 1, j
                fc(second(i, j)) = h(i, j)
            end do
        end do
    end subroutine flux_map_jet_to_cylindrical

    pure subroutine toroidal_jet(nt, t, l)
        !! Trigonometric nodal basis of nt = 2M + 1 planes at t with its
        !! first and second t derivatives.
        integer, intent(in) :: nt
        real(dp), intent(in) :: t
        real(dp), intent(out) :: l(nt, 3)
        real(dp) :: x, w
        integer :: k, m

        do k = 1, nt
            x = two_pi*(t - (k - 1)/real(nt, dp))
            l(k, :) = [1.0_dp, 0.0_dp, 0.0_dp]
            do m = 1, (nt - 1)/2
                w = two_pi*m
                l(k, :) = l(k, :) + 2*[cos(m*x), -w*sin(m*x), -w**2*cos(m*x)]
            end do
        end do
        l = l/nt
    end subroutine toroidal_jet

    pure function parameter_jet(nj, l, c) result(j)
        !! 2-jet in (X_1, X_2, t) of sum_ik c(i, k) N_i(X) L_k(t) from the
        !! section jets nj(:, i) = (N, N_1, N_2, N_11, N_12, N_22).
        real(dp), intent(in) :: nj(:, :), l(:, :), c(:, :)
        real(dp) :: j(10), w(size(c, 1), 3)

        w = matmul(c, l)
        j(1) = dot_product(nj(1, :), w(:, 1))
        j(2) = dot_product(nj(2, :), w(:, 1))
        j(3) = dot_product(nj(3, :), w(:, 1))
        j(4) = dot_product(nj(1, :), w(:, 2))
        j(5) = dot_product(nj(4, :), w(:, 1))
        j(6) = dot_product(nj(5, :), w(:, 1))
        j(7) = dot_product(nj(2, :), w(:, 2))
        j(8) = dot_product(nj(6, :), w(:, 1))
        j(9) = dot_product(nj(3, :), w(:, 2))
        j(10) = dot_product(nj(1, :), w(:, 3))
    end function parameter_jet

    subroutine state_jets(fm, plane_nodes, z, d0, nj, l, rj, zj, xc, yc)
        !! At one point of the section (functions d0 with section jets nj)
        !! and one t (toroidal basis jets l): the 2-jets of R and Z in the
        !! parameters (X_1, X_2, t) and of x, y of the state z in
        !! (R, phi, Z), the wall part of u included.
        type(flux_map3d_t), intent(in) :: fm
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:), nj(:, :), l(:, :)
        integer, intent(in) :: d0(:)
        real(dp), intent(out) :: rj(10), zj(10), xc(10), yc(10)
        real(dp), dimension(size(d0), fm%nt) :: cr, cz, cx, cy, cw
        real(dp) :: cs(size(d0), 1), sj(10), wj(10), tj(10), cj(10), nj2(10), xq(10), yq(10), &
                    ref, shift
        integer :: n0, ntot, i, k
        logical :: wall

        n0 = fm%n0
        ntot = n0*fm%nt
        cx = 0
        cy = 0
        cw = 0
        cs = 0
        wall = .false.
        ref = 0
        do i = 1, size(d0)
            if (fm%is_wall(d0(i)) .and. .not. wall) then
                ref = z(d0(i))
                wall = .true.
            end if
        end do
        do i = 1, size(d0)
            cr(i, :) = plane_nodes(1, d0(i), :)
            cz(i, :) = plane_nodes(2, d0(i), :)
            if (fm%is_wall(d0(i))) then
                ! the wall angles of the point on one branch
                shift = two_pi*nint((ref - z(d0(i)))/two_pi)
                cs(i, 1) = 1
                do k = 1, fm%nt
                    cw(i, k) = z((k - 1)*n0 + d0(i)) + shift
                end do
            else
                do k = 1, fm%nt
                    cx(i, k) = z((k - 1)*n0 + d0(i))
                    cy(i, k) = z(ntot + (k - 1)*n0 + d0(i))
                end do
            end if
        end do
        rj = parameter_jet(nj, l, cr)
        zj = parameter_jet(nj, l, cz)
        xq = parameter_jet(nj, l, cx)
        yq = parameter_jet(nj, l, cy)
        if (wall) then
            ! sigma does not depend on t
            sj = parameter_jet(nj, reshape([1.0_dp, 0.0_dp, 0.0_dp], [1, 3]), cs)
            wj = parameter_jet(nj, l, cw)
            tj = flux_map_jet_quotient(wj, sj)
            call flux_map_jet_cos_sin(tj, cj, nj2)
            xq = xq + flux_map_jet_product(sj, cj)
            yq = yq + flux_map_jet_product(sj, nj2)
        end if
        call flux_map_jet_to_cylindrical(rj, zj, two_pi/clebsch_nfp, xq, xc)
        call flux_map_jet_to_cylindrical(rj, zj, two_pi/clebsch_nfp, yq, yc)
    end subroutine state_jets

    subroutine flux_map3d_force_point(fm, derham, profiles, plane_nodes, z, block, xi, t, rz, b, &
                                      j, f, status, label, grad_b2, volume)
        !! Position (R, Z), B, J = curl B and F = J X B - grad p, components
        !! (R, phi, Z), of the state z at the block parameters xi and the
        !! toroidal parameter t.
        type(flux_map3d_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:), xi(2), t
        integer, intent(in) :: block
        real(dp), intent(out) :: rz(2), b(3), j(3), f(3)
        integer, intent(out) :: status
        !! optional: the label s, grad(B^2/2) and the volume element per unit
        !! of the block parameters and of t
        real(dp), intent(out), optional :: label, grad_b2(3), volume
        integer, allocatable :: rowptr(:), colidx(:), d0(:)
        real(dp), allocatable :: v(:), nj(:, :)
        real(dp) :: l(fm%nt, 3), rj(10), zj(10), xc(10), yc(10), tau(3), gp(3), pr(3), m(3), jac
        integer :: c

        do c = 0, 5
            call spline_point_matrix(derham%space, [block], reshape(xi, [2, 1]), c, rowptr, &
                                     colidx, v, status)
            if (status /= 0) return
            if (c == 0) then
                d0 = colidx
                allocate (nj(6, size(d0)))
            else
                status = 1
                if (size(colidx) /= size(d0)) return
            end if
            nj(c + 1, :) = v
        end do
        status = 0
        call toroidal_jet(fm%nt, t, l)
        call state_jets(fm, plane_nodes, z, d0, nj, l, rj, zj, xc, yc)
        rz = [rj(1), zj(1)]
        call profiles(xc(1)**2 + yc(1)**2, tau, gp, pr)
        call flux_map_force_3d(rj(1), xc(1), yc(1), xc(2), xc(3), xc(4), xc(5), xc(6), xc(7), &
            xc(8), xc(9), xc(10), yc(2), yc(3), yc(4), yc(5), yc(6), yc(7), yc(8), yc(9), yc(10), &
            tau(1), tau(2), gp(1), gp(2), pr(2), b(1), b(2), b(3), j(1), j(2), j(3), f(1), f(2), &
            f(3), m(1), m(2), m(3), jac)
        if (present(label)) label = xc(1)**2 + yc(1)**2
        if (present(grad_b2)) grad_b2 = m
        if (present(volume)) volume = abs(rj(2)*zj(3) - rj(3)*zj(2))*rj(1)*two_pi/clebsch_nfp
    end subroutine flux_map3d_force_point

    subroutine flux_map3d_force(fm, derham, quad, profiles, plane_nodes, z, bands, acc, status)
        !! The force balance of the state z over the quadrature points of
        !! quad and the toroidal quadrature of flux_map3d_system.
        !! Status 1: the section jets could not be evaluated.
        type(flux_map3d_t), intent(in) :: fm
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: plane_nodes(:, :, :), z(:)
        integer, intent(in) :: bands
        type(flux_map_force_t), intent(out) :: acc
        integer, intent(out) :: status
        type(flux_map_force_t) :: part
        integer, allocatable :: rowptr(:), colidx(:), rp(:), ci(:)
        real(dp), allocatable :: val(:, :), v(:), ltab(:, :, :)
        real(dp) :: rj(10), zj(10), xc(10), yc(10)
        integer :: nt, nqt, p, m, c

        nt = fm%nt
        nqt = fm%toroidal_quadrature_count()
        call flux_map_force_reset(acc, bands)
        do c = 0, 5
            call spline_point_matrix(derham%space, quad%block, quad%xi, c, rp, ci, v, status)
            if (status /= 0) then
                status = 1
                return
            end if
            if (c == 0) then
                rowptr = rp
                colidx = ci
                allocate (val(6, size(v)))
            else
                status = 1
                if (size(ci) /= size(colidx)) return
                if (any(ci /= colidx)) return
            end if
            val(c + 1, :) = v
        end do
        status = 0
        allocate (ltab(nt, 3, nqt))
        do m = 1, nqt
            call toroidal_jet(nt, (m - 1)/real(nqt, dp), ltab(:, :, m))
        end do
        !$omp parallel default(private) shared(fm, quad, plane_nodes, z, acc, rowptr, colidx, &
        !$omp     val, ltab, nqt, bands)
        call flux_map_force_reset(part, bands)
        !$omp do schedule(dynamic, 16)
        do p = 1, quad%point_count
            do m = 1, nqt
                call state_jets(fm, plane_nodes, z, colidx(rowptr(p):rowptr(p + 1) - 1), &
                                val(:, rowptr(p):rowptr(p + 1) - 1), ltab(:, :, m), rj, zj, xc, yc)
                call flux_map_force_add(part, quad%weight(p)*abs(rj(2)*zj(3) - rj(3)*zj(2)) &
                                        *two_pi/nqt*rj(1), rj(1), xc, yc, profiles)
            end do
        end do
        !$omp end do
        !$omp critical (flux_map_force_merge_parts)
        call flux_map_force_merge(acc, part)
        !$omp end critical (flux_map_force_merge_parts)
        !$omp end parallel
    end subroutine flux_map3d_force

end module milca_flux_map_force
