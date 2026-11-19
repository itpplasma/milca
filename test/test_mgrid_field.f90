program test_mgrid_field
    ! mgrid field basis against an analytic vacuum field written to synthetic
    ! mgrid files. Groups (unit EXTCUR), each B = grad Phi with Phi harmonic:
    !   1: Phi = phi                          (toroidal field 1/R)
    !   2: Phi = Z log R                      (B_R = Z/R, B_Z = log R)
    !   3: Phi = Z (R0/R)^N cos(N phi)        (N = nfp, grid harmonic 1)
    !      + Z (R0/R)^(2N) sin(2 N phi)/2     (grid harmonic 2)
    ! Checks: file round trip, harmonics against the exact ones at points off
    ! the nodes with the error falling as h^(p+1) for p = 3 and 5, the raw
    ! mode scaling by raw_coil_cur, unresolved harmonics zero, and rejection
    ! of points outside the grid; no coil found in the vacuum field and a
    ! current sheet found where it is.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_mgrid_file, only: mgrid_t, mgrid_read, mgrid_write
    use milca_ideal_fb_mgrid, only: mgrid_basis_t, mgrid_basis_init, mgrid_current_nodes
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 1.5_dp
    integer, parameter :: nfp = 3, kp = 8, npoint = 40, ntor = 3
    real(dp) :: points(2, npoint), error(2, 2), ratio
    integer :: q, degree, level
    character(256) :: mode, path
    if (command_argument_count() == 2) then
        call get_command_argument(1, mode)
        call get_command_argument(2, path)
        if (trim(mode) /= 'tokamak') error stop 'usage: test_mgrid_field [tokamak FILE]'
        call tokamak(trim(path))
        stop
    end if
    do q = 1, npoint
        ! deterministic points off the nodes, away from the grid edge
        points(:, q) = [1.05_dp + 0.9_dp*frac(0.6180339887_dp*q), &
            -0.45_dp + 0.9_dp*frac(0.7548776662_dp*q)]
    end do
    do degree = 3, 5, 2
        do level = 1, 2
            error(level, (degree - 1)/2) = sample_error(16*2**level, degree, level == 1)
        end do
        ratio = error(1, (degree - 1)/2)/error(2, (degree - 1)/2)
        print '(a,i0,a,2es11.3,a,f6.2)', 'degree ', degree, ': max error ', &
            error(:, (degree - 1)/2), ' T, ratio ', ratio
        ! h^(p+1): 16 and 64 for p = 3 and 5, allowed down to 70 %
        if (ratio < 0.7_dp*2.0_dp**(degree + 1)) error stop 'interpolation order'
    end do
    ! |B| reaches 43 T on this grid: about 2e-6 and 2e-8 relative
    if (error(2, 1) > 8e-5_dp .or. error(2, 2) > 8e-7_dp) &
        error stop 'interpolation error'
    call raw_mode_and_bounds()
    call coil_detection()
    print '(a)', 'mgrid field basis: order, harmonics, scaling and bounds verified'
contains
    pure real(dp) function frac(x)
        real(dp), intent(in) :: x
        frac = x - floor(x)
    end function

    pure subroutine exact(r, z, phi, b)
        ! b(:, group) = (B_R, B_phi, B_Z)
        real(dp), intent(in) :: r, z, phi
        real(dp), intent(out) :: b(3, 3)
        real(dp) :: u1, u2
        integer, parameter :: n = nfp
        b(:, 1) = [0.0_dp, 1/r, 0.0_dp]
        b(:, 2) = [z/r, 0.0_dp, log(r)]
        u1 = (r0/r)**n; u2 = (r0/r)**(2*n)/2
        b(:, 3) = [-n*z*u1/r*cos(n*phi) - 2*n*z*u2/r*sin(2*n*phi), &
            -n*z*u1/r*sin(n*phi) + 2*n*z*u2/r*cos(2*n*phi), &
            u1*cos(n*phi) + u2*sin(2*n*phi)]
    end subroutine

    subroutine synthetic(grid, nr, mode, raw)
        type(mgrid_t), intent(out) :: grid
        integer, intent(in) :: nr
        character, intent(in) :: mode
        real(dp), intent(in) :: raw
        real(dp) :: r, z, phi
        integer :: i, j, l
        grid%ir = nr; grid%jz = nr; grid%kp = kp; grid%nfp = nfp; grid%nextcur = 3
        grid%rmin = 0.9_dp; grid%rmax = 2.1_dp; grid%zmin = -0.6_dp; grid%zmax = 0.6_dp
        grid%mode = mode
        allocate (grid%raw_coil_cur(3), grid%b(3, nr, nr, kp, 3))
        grid%raw_coil_cur = raw
        do l = 1, kp
            phi = 2*pi*(l - 1)/(nfp*kp)
            do j = 1, nr
                z = grid%zmin + (j - 1)*(grid%zmax - grid%zmin)/(nr - 1)
                do i = 1, nr
                    r = grid%rmin + (i - 1)*(grid%rmax - grid%rmin)/(nr - 1)
                    call exact(r, z, phi, grid%b(:, i, j, l, :))
                end do
            end do
        end do
        grid%b = raw*grid%b
    end subroutine

    real(dp) function sample_error(nr, degree, round_trip) result(worst)
        integer, intent(in) :: nr, degree
        logical, intent(in) :: round_trip
        type(mgrid_t) :: grid, back
        type(mgrid_basis_t) :: basis
        real(dp) :: bcos(3, 0:ntor, 3, npoint), bsin(3, 0:ntor, 3, npoint)
        real(dp) :: b(3, 3), sum(3, 3), phi
        integer :: st, q, j, k
        character(256) :: message
        call synthetic(grid, nr, 'S', 1.0_dp)
        call mgrid_write('mgrid_synthetic.nc', grid, st)
        if (st /= 0) error stop 'mgrid write'
        call mgrid_read('mgrid_synthetic.nc', back, st, message)
        if (st /= 0) error stop 'mgrid read'
        if (round_trip) then
            if (any([back%ir, back%jz, back%kp, back%nfp, back%nextcur] &
                /= [grid%ir, grid%jz, grid%kp, grid%nfp, grid%nextcur])) &
                error stop 'round trip dimensions'
            if (any(back%b /= grid%b) .or. back%rmax /= grid%rmax) &
                error stop 'round trip data'
        end if
        call mgrid_basis_init(basis, back%b, back%nfp, back%rmin, back%rmax, &
            back%zmin, back%zmax, back%unit_scale(), degree, st)
        if (st /= 0) error stop 'basis init'
        call basis%sample(points, nfp, ntor, bcos, bsin, st)
        if (st /= 0) error stop 'basis sample'
        ! grid harmonic 3 of the field period is unresolved by kp = 8 planes
        ! (Nyquist 4) only from k = 4 on; the field has none at k = 3
        worst = 0
        do q = 1, npoint
            do j = 0, 6
                phi = 2*pi*j/(7*nfp)
                sum = 0
                do k = 0, ntor
                    sum = sum + bcos(:, k, :, q)*cos(k*nfp*phi) + bsin(:, k, :, q) &
                        *sin(k*nfp*phi)
                end do
                call exact(points(1, q), points(2, q), phi, b)
                worst = max(worst, maxval(abs(sum - b)))
            end do
        end do
    end function

    subroutine tokamak(path)
        ! Synthetic axisymmetric coil field of the free-boundary app test,
        ! per unit EXTCUR: 1 the toroidal field 1/R, 2 a uniform vertical
        ! field, 3 grad(Z log R), whose B_Z = log R gives the vertical field
        ! a decay index. 101 x 101 nodes on [0.55, 1.45] x [-0.45, 0.45].
        character(*), intent(in) :: path
        type(mgrid_t) :: grid
        real(dp) :: r, z
        integer :: i, j, st, nr
        nr = 101
        grid%ir = nr; grid%jz = nr; grid%kp = 1; grid%nfp = 1; grid%nextcur = 3
        grid%rmin = 0.55_dp; grid%rmax = 1.45_dp
        grid%zmin = -0.45_dp; grid%zmax = 0.45_dp
        allocate (grid%raw_coil_cur(3), grid%b(3, nr, nr, 1, 3))
        grid%raw_coil_cur = 1
        do j = 1, nr
            z = grid%zmin + (j - 1)*(grid%zmax - grid%zmin)/(nr - 1)
            do i = 1, nr
                r = grid%rmin + (i - 1)*(grid%rmax - grid%rmin)/(nr - 1)
                grid%b(:, i, j, 1, 1) = [0.0_dp, 1/r, 0.0_dp]
                grid%b(:, i, j, 1, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
                grid%b(:, i, j, 1, 3) = [z/r, 0.0_dp, log(r)]
            end do
        end do
        call mgrid_write(path, grid, st)
        if (st /= 0) error stop 'mgrid write'
    end subroutine

    subroutine coil_detection()
        ! The analytic groups are curl free on the 3D grid (kp = 8): no node
        ! may be flagged. Adding the sheet B_Z = tanh((R - 1.5)/0.02), whose
        ! toroidal current density -B_Z'(R) sits at R = 1.5, flags nodes there
        ! only.
        type(mgrid_t) :: grid
        type(mgrid_basis_t) :: basis
        real(dp), allocatable :: points(:, :), b(:, :, :, :, :)
        real(dp) :: r
        integer :: st, i
        call synthetic(grid, 64, 'S', 1.0_dp)
        call mgrid_basis_init(basis, grid%b, nfp, grid%rmin, grid%rmax, grid%zmin, &
            grid%zmax, grid%unit_scale(), 5, st)
        call mgrid_current_nodes(basis, [1.0_dp, 1.0_dp, 1.0_dp], points)
        if (size(points, 2) /= 0) error stop 'coil found in a vacuum field'
        allocate (b(3, 64, 64, kp, 4))
        b(:, :, :, :, :3) = grid%b
        b(:, :, :, :, 4) = 0
        do i = 1, 64
            r = grid%rmin + (i - 1)*(grid%rmax - grid%rmin)/63
            b(3, i, :, :, 4) = tanh((r - 1.5_dp)/0.02_dp)
        end do
        call mgrid_basis_init(basis, b, nfp, grid%rmin, grid%rmax, grid%zmin, &
            grid%zmax, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 5, st)
        call mgrid_current_nodes(basis, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], points)
        if (size(points, 2) == 0) error stop 'current sheet not found'
        if (maxval(abs(points(1, :) - 1.5_dp)) > 0.2_dp) &
            error stop 'coil nodes misplaced'
        print '(a,i0,a,2f7.3)', 'current sheet: ', size(points, 2), ' nodes in R ', &
            minval(points(1, :)), maxval(points(1, :))
    end subroutine

    subroutine raw_mode_and_bounds()
        ! raw grids store the field at raw_coil_cur; the basis is per unit EXTCUR
        type(mgrid_t) :: grid, back
        type(mgrid_basis_t) :: scaled, raw
        real(dp) :: c1(3, 0:ntor, 3, npoint), s1(3, 0:ntor, 3, npoint)
        real(dp) :: c2(3, 0:ntor, 3, npoint), s2(3, 0:ntor, 3, npoint)
        real(dp) :: outside(2, 1), co(3, 0:ntor, 3, 1), so(3, 0:ntor, 3, 1)
        integer :: st
        character(256) :: message
        call synthetic(grid, 32, 'S', 1.0_dp)
        call mgrid_basis_init(scaled, grid%b, nfp, grid%rmin, grid%rmax, grid%zmin, &
            grid%zmax, grid%unit_scale(), 5, st)
        call synthetic(grid, 32, 'R', 2.5_dp)
        call mgrid_write('mgrid_synthetic_raw.nc', grid, st)
        call mgrid_read('mgrid_synthetic_raw.nc', back, st, message)
        if (st /= 0 .or. back%mode /= 'R') error stop 'raw mode read'
        call mgrid_basis_init(raw, back%b, nfp, back%rmin, back%rmax, back%zmin, &
            back%zmax, back%unit_scale(), 5, st)
        call scaled%sample(points, nfp, ntor, c1, s1, st)
        call raw%sample(points, nfp, ntor, c2, s2, st)
        if (maxval(abs(c1 - c2)) + maxval(abs(s1 - s2)) > 1e-14_dp) &
            error stop 'raw scaling'
        ! grid harmonics 1 and 2 only: k = 3 is absent from the field, zero
        if (maxval(abs(c1(:, 3, :, :))) + maxval(abs(s1(:, 3, :, :))) > 1e-13_dp) &
            error stop 'absent harmonic'
        outside(:, 1) = [0.5_dp, 0.0_dp]
        call scaled%sample(outside, nfp, ntor, co, so, st)
        if (st == 0) error stop 'point outside the grid accepted'
        ! a plasma nfp that is not a multiple of the grid's is rejected
        call scaled%sample(points, 2, ntor, c1, s1, st)
        if (st == 0) error stop 'incompatible nfp accepted'
    end subroutine
end program test_mgrid_field
