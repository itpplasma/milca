module milca_ideal_fb_mgrid
    !! Coil-group fields tabulated on an mgrid (MAKEGRID's cylindrical grid)
    !! as a field basis of the nested-map free boundary
    !! (docs/free-boundary-nested-map-3d.md, section 2). Field k of the basis
    !! is coil group k at unit EXTCUR, so the vacuum field of a deck is
    !! sum_k EXTCUR_k b_k (ideal_fb_set_currents).
    !!
    !! Grid: R_i = rmin + (i - 1) dR, i = 1..ir; Z_j = zmin + (j - 1) dZ,
    !! j = 1..jz; phi_l = 2 pi (l - 1)/(nfp kp), l = 1..kp, over one field
    !! period of the grid. In (R, Z) every plane is interpolated by the
    !! tensor Lagrange polynomial of degree p on the (p + 1)^2 nearest nodes
    !! (error O(h^(p+1)) in B, p = 3 or 5); across the planes the toroidal
    !! harmonics are the exact discrete Fourier coefficients, which is the
    !! trigonometric interpolant through the kp planes. Harmonics at or
    !! above the grid's Nyquist index kp/2 are not resolved and are zero.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_polynomial, only: lagrange_weights
    use milca_ideal_free_boundary, only: ideal_fb_field_basis_t, IDEAL_FB_OK, &
        IDEAL_FB_BAD_INPUT
    implicit none
    private
    public :: mgrid_basis_t, mgrid_basis_init, mgrid_current_nodes
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, extends(ideal_fb_field_basis_t) :: mgrid_basis_t
        integer :: ir = 0, jz = 0, kp = 0, nfp = 1, ngroup = 0, degree = 5
        real(dp) :: rmin = 0, rmax = 0, zmin = 0, zmax = 0
        !> (B_R, B_phi, B_Z) per node, plane and group at unit EXTCUR.
        real(dp), allocatable :: b(:, :, :, :, :)
    contains
        procedure :: nfield => mgrid_nfield
        procedure :: sample => mgrid_sample
    end type

contains

    subroutine mgrid_basis_init(basis, b, nfp, rmin, rmax, zmin, zmax, scale, degree, &
                                status)
        !! b(3, ir, jz, kp, ngroup) as stored in the file; scale(k) converts
        !! group k to unit EXTCUR (1 for scaled grids, 1/raw_coil_cur for raw).
        type(mgrid_basis_t), intent(out) :: basis
        real(dp), intent(in) :: b(:, :, :, :, :), rmin, rmax, zmin, zmax, scale(:)
        integer, intent(in) :: nfp, degree
        integer, intent(out) :: status
        integer :: k
        status = IDEAL_FB_BAD_INPUT
        if (size(b, 1) /= 3 .or. nfp < 1) return
        if (degree /= 1 .and. degree /= 3 .and. degree /= 5) return
        if (size(b, 2) < degree + 1 .or. size(b, 3) < degree + 1) return
        if (size(b, 4) < 1 .or. size(b, 5) < 1 .or. size(scale) /= size(b, 5)) return
        if (.not. (rmax > rmin .and. zmax > zmin .and. rmin > 0)) return
        if (.not. all(ieee_is_finite(scale))) return
        basis%ir = size(b, 2); basis%jz = size(b, 3); basis%kp = size(b, 4)
        basis%ngroup = size(b, 5); basis%nfp = nfp; basis%degree = degree
        basis%rmin = rmin; basis%rmax = rmax; basis%zmin = zmin; basis%zmax = zmax
        allocate (basis%b, mold=b)
        do k = 1, basis%ngroup
            basis%b(:, :, :, :, k) = scale(k)*b(:, :, :, :, k)
        end do
        if (.not. all(ieee_is_finite(basis%b))) return
        status = IDEAL_FB_OK
    end subroutine

    integer function mgrid_nfield(self)
        class(mgrid_basis_t), intent(in) :: self
        mgrid_nfield = self%ngroup
    end function

    subroutine mgrid_sample(self, rz, nfp, ntor, bcos, bsin, status)
        !! Harmonics n = k nfp, k = 0..ntor, of every group field at each
        !! meridional point (the basis contract). The plasma's nfp must be a
        !! multiple of the grid's, so n is a grid harmonic l = n/nfp_grid.
        !! Points are independent (OpenMP over points).
        class(mgrid_basis_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: ck(:, :), sk(:, :)
        real(dp) :: plane(3, self%kp, self%ngroup), wr(self%degree + 1)
        real(dp) :: wz(self%degree + 1), x, y, dr, dz
        integer :: q, k, l, g, i0, j0, a, c, np, harmonic, failed
        status = IDEAL_FB_BAD_INPUT
        if (size(rz, 1) /= 2 .or. nfp < 1 .or. ntor < 0) return
        if (mod(nfp, self%nfp) /= 0 .and. self%kp > 1) return
        if (size(bcos, 1) /= 3 .or. ubound(bcos, 2) /= ntor) return
        if (size(bcos, 3) /= self%ngroup .or. size(bcos, 4) /= size(rz, 2)) return
        if (any(shape(bsin) /= shape(bcos))) return
        np = self%degree + 1
        dr = (self%rmax - self%rmin)/(self%ir - 1)
        dz = (self%zmax - self%zmin)/(self%jz - 1)
        ! Discrete Fourier weights of grid harmonic l = k nfp/nfp_grid over
        ! the kp planes; unresolved harmonics (2 l >= kp) stay zero.
        allocate (ck(self%kp, 0:ntor), sk(self%kp, 0:ntor))
        ck = 0; sk = 0
        do k = 0, ntor
            harmonic = k*nfp/self%nfp
            if (2*harmonic >= self%kp .and. k > 0) cycle
            do l = 1, self%kp
                ck(l, k) = merge(1.0_dp, 2.0_dp, k == 0) &
                    *cos(2*pi*mod(harmonic*(l - 1), self%kp)/self%kp)/self%kp
                sk(l, k) = merge(0.0_dp, 2.0_dp, k == 0) &
                    *sin(2*pi*mod(harmonic*(l - 1), self%kp)/self%kp)/self%kp
            end do
        end do
        failed = 0
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(self, rz, ntor, np, dr, dz, ck, sk, bcos, bsin) &
        !$omp private(q, k, l, g, i0, j0, a, c, x, y, wr, wz, plane) &
        !$omp reduction(+:failed)
        do q = 1, size(rz, 2)
            bcos(:, :, :, q) = 0; bsin(:, :, :, q) = 0
            ! Grid coordinates from 0; the stencil is centred on the point and
            ! shifted inward at the grid edges.
            x = (rz(1, q) - self%rmin)/dr; y = (rz(2, q) - self%zmin)/dz
            if (.not. (x >= 0 .and. x <= self%ir - 1 .and. y >= 0 &
                .and. y <= self%jz - 1)) then
                failed = failed + 1
                cycle
            end if
            i0 = min(max(int(x) - (np/2 - 1), 0), self%ir - np)
            j0 = min(max(int(y) - (np/2 - 1), 0), self%jz - np)
            call lagrange_weights(np, x - i0, [(real(a, dp), a=0, np - 1)], wr)
            call lagrange_weights(np, y - j0, [(real(a, dp), a=0, np - 1)], wz)
            plane = 0
            do g = 1, self%ngroup
                do l = 1, self%kp
                    do c = 1, np
                        do a = 1, np
                            plane(:, l, g) = plane(:, l, g) + wr(a)*wz(c) &
                                *self%b(:, i0 + a, j0 + c, l, g)
                        end do
                    end do
                end do
            end do
            do g = 1, self%ngroup
                do k = 0, ntor
                    bcos(:, k, g, q) = matmul(plane(:, :, g), ck(:, k))
                    bsin(:, k, g, q) = matmul(plane(:, :, g), sk(:, k))
                end do
            end do
        end do
        !$omp end parallel do
        if (failed == 0) status = IDEAL_FB_OK
    end subroutine
    subroutine mgrid_current_nodes(basis, currents, points)
        !! (R, Z) of the grid nodes, in any plane, where the field of the
        !! groups at these currents carries current: a coil inside the grid.
        !! A node is flagged when curl B, by second-order central differences
        !! in (R, Z) and the derivative of the trigonometric interpolant in
        !! phi, exceeds 5 % of |B|/h, h the smaller grid step (|B_pol| for the
        !! toroidal current): it is O(h^2) in a vacuum field and O(|B|/h)
        !! within about three cells of a filament. The flags are
        !! widened by the half-width of the interpolation stencil, so that no
        !! interpolated value outside them uses a node next to a coil.
        class(mgrid_basis_t), intent(in) :: basis
        real(dp), intent(in) :: currents(:)
        real(dp), allocatable, intent(out) :: points(:, :)
        real(dp), allocatable :: b(:, :, :, :)
        logical, allocatable :: flagged(:, :), core(:, :)
        real(dp), allocatable :: dphi(:, :)
        real(dp) :: dr, dz, h, r, curl(3), pol, full, d(3, 3), delta
        integer :: i, j, l, g, m, k, n
        dr = (basis%rmax - basis%rmin)/(basis%ir - 1)
        dz = (basis%zmax - basis%zmin)/(basis%jz - 1)
        h = min(dr, dz)
        ! dphi(l, m): derivative at plane l of the trigonometric interpolant
        ! through the planes m (the Nyquist harmonic of an even kp dropped)
        allocate (dphi(basis%kp, basis%kp))
        dphi = 0
        do m = 1, basis%kp
            do l = 1, basis%kp
                delta = 2*pi*(l - m)/basis%kp
                do k = 1, (basis%kp - 1)/2
                    dphi(l, m) = dphi(l, m) - 2.0_dp/basis%kp*k*basis%nfp*sin(k*delta)
                end do
            end do
        end do
        allocate (b(3, basis%ir, basis%jz, basis%kp), flagged(basis%ir, basis%jz))
        b = 0
        do g = 1, min(basis%ngroup, size(currents))
            b = b + currents(g)*basis%b(:, :, :, :, g)
        end do
        flagged = .false.
        do l = 1, basis%kp
            do j = 2, basis%jz - 1
                do i = 2, basis%ir - 1
                    r = basis%rmin + (i - 1)*dr
                    ! d(:, c): derivatives of component c in R, phi, Z
                    d(1, :) = (b(:, i + 1, j, l) - b(:, i - 1, j, l))/(2*dr)
                    d(3, :) = (b(:, i, j + 1, l) - b(:, i, j - 1, l))/(2*dz)
                    d(2, :) = matmul(b(:, i, j, :), dphi(l, :))
                    curl = [d(2, 3)/r - d(3, 2), d(3, 1) - d(1, 3), &
                        d(1, 2) + b(2, i, j, l)/r - d(2, 1)/r]
                    pol = norm2(b([1, 3], i, j, l)); full = norm2(b(:, i, j, l))
                    if (abs(curl(2)) > 0.05_dp*pol/h .or. &
                        norm2(curl([1, 3])) > 0.05_dp*full/h) flagged(i, j) = .true.
                end do
            end do
        end do
        core = flagged
        n = (basis%degree + 1)/2
        do j = 1, basis%jz
            do i = 1, basis%ir
                if (core(i, j)) flagged(max(1, i - n):min(basis%ir, i + n), &
                    max(1, j - n):min(basis%jz, j + n)) = .true.
            end do
        end do
        n = count(flagged)
        allocate (points(2, n))
        n = 0
        do j = 1, basis%jz
            do i = 1, basis%ir
                if (.not. flagged(i, j)) cycle
                n = n + 1
                points(:, n) = [basis%rmin + (i - 1)*dr, basis%zmin + (j - 1)*dz]
            end do
        end do
    end subroutine
end module milca_ideal_fb_mgrid
