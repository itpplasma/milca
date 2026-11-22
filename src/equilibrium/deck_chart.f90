module milca_deck_chart
    !! A curve chart for a VMEC-style deck: the guide curve is the m = 0 part
    !! of the boundary (its centroid curve), and the boundary is expressed in
    !! the normal-plane coordinates (q1, q2) of that chart on every plane
    !! phi-hat, as a Fourier series in the geometric angle about the curve.
    !! On a strongly helical deck the sections then stay close to a fixed
    !! disc in the chart, so a coupling surface Gamma_c that is a fixed
    !! contour in the normal plane follows the plasma at a uniform distance.
    !!
    !! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp phi),
    !! the layout of nested_equilibrium and milca_import.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_chart, only: curve_chart_t
    implicit none
    private
    public :: deck_axis_chart, deck_boundary_to_chart

    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    subroutine deck_axis_chart(boundary, ntor, nfp, chart, status)
        !! The chart along the m = 0 curve of the boundary, R(phi), Z(phi),
        !! as full-turn Cartesian Fourier modes k <= nfp ntor + 1 (x and y
        !! carry the k = j nfp +- 1 harmonics, z the k = j nfp ones).
        integer, intent(in) :: ntor, nfp
        real(dp), intent(in) :: boundary(:, -ntor:, 0:)
        type(curve_chart_t), intent(out) :: chart
        integer, intent(out) :: status
        integer :: kmax, nphi, j, k, n
        real(dp), allocatable :: ac(:, :), as(:, :), xyz(:, :)
        real(dp) :: phi, r, z, arg

        kmax = nfp*ntor + 1
        nphi = 8*(kmax + 1)
        allocate (xyz(3, nphi), ac(3, 0:kmax), as(3, 0:kmax))
        do j = 1, nphi
            phi = 2*pi*(j - 1)/nphi
            r = 0; z = 0
            do n = -ntor, ntor
                arg = -n*nfp*phi
                r = r + boundary(1, n, 0)*cos(arg) + boundary(2, n, 0)*sin(arg)
                z = z + boundary(3, n, 0)*cos(arg) + boundary(4, n, 0)*sin(arg)
            end do
            xyz(:, j) = [r*cos(phi), r*sin(phi), z]
        end do
        ac = 0; as = 0
        do k = 0, kmax
            do j = 1, nphi
                phi = 2*pi*(j - 1)/nphi
                ac(:, k) = ac(:, k) + xyz(:, j)*cos(k*phi)
                as(:, k) = as(:, k) + xyz(:, j)*sin(k*phi)
            end do
        end do
        ac = ac*2/nphi; as = as*2/nphi
        ac(:, 0) = ac(:, 0)/2; as(:, 0) = 0
        ! keep the harmonics a period-nfp curve can carry
        do k = 0, kmax
            if (mod(k - 1, nfp) /= 0 .and. mod(k + 1, nfp) /= 0) then
                ac(1:2, k) = 0; as(1:2, k) = 0
            end if
            if (mod(k, nfp) /= 0) then
                ac(3, k) = 0; as(3, k) = 0
            end if
        end do
        where (abs(ac) < 1e-13_dp*maxval(abs(ac))) ac = 0
        where (abs(as) < 1e-13_dp*maxval(abs(ac))) as = 0
        call chart%init(ac, as, nfp, status)
    end subroutine deck_axis_chart

    subroutine deck_boundary_to_chart(boundary, ntor, nfp, chart, mpol_out, ntor_out, qb, &
                                      status)
        !! The boundary in chart coordinates: on nphi planes phi-hat of a
        !! period the points of the surface in the normal plane, (q1, q2),
        !! resampled at uniform geometric angle about the curve and transformed
        !! in (alpha, phi-hat) to qb(1:4, n, m) (the layout of boundary with
        !! (q1, q2) in place of (R, Z)), m <= mpol_out, |n| <= ntor_out (the
        !! frame turns with the curve, so the chart sections need more toroidal
        !! modes than the input). status 1 when a point
        !! of the surface does not cross a normal plane near its own plane.
        integer, intent(in) :: ntor, nfp, mpol_out, ntor_out
        real(dp), intent(in) :: boundary(:, -ntor:, 0:)
        type(curve_chart_t), intent(in) :: chart
        real(dp), allocatable, intent(out) :: qb(:, :, :)
        integer, intent(out) :: status
        integer, parameter :: nth = 1024
        integer :: mpol, nphi, j, i, it, m, n, k
        real(dp) :: ph, th, f, df, x0(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), &
                    nrm(3), p(3), dp3(3), q(2), w, arg, target, phase
        real(dp), allocatable :: qs(:, :, :), al(:), qq(:, :), resampled(:, :, :), ph0(:)

        status = 1
        mpol = ubound(boundary, 3)
        nphi = max(1, 4*ntor_out + 1)
        allocate (qs(2, nth + 1, nphi), al(nth + 1), qq(2, nth + 1), &
                  resampled(2, nth, nphi), ph0(nphi))
        do j = 1, nphi
            ph = 2*pi*(j - 1)/(nphi*nfp)
            call chart%jet([0.0_dp, 0.0_dp, ph], x0, d1, d2, d3)
            nrm = cross(d1(:, 1), d1(:, 2))
            nrm = nrm/norm2(nrm)
            do i = 1, nth
                th = 2*pi*(i - 1)/nth
                ! the cylindrical angle where the surface line theta meets the
                ! normal plane: Newton on (X(theta, phi) - x0) . n = 0
                f = ph
                do it = 1, 40
                    call surface(th, f, p, dp3)
                    df = dot_product(dp3, nrm)
                    if (abs(df) < tiny(1.0_dp)) return
                    w = dot_product(p - x0, nrm)/df
                    f = f - w
                    if (abs(w) < 1e-14_dp) exit
                end do
                if (abs(f - ph) > pi/nfp) return
                call surface(th, f, p, dp3)
                q = [dot_product(p - x0, d1(:, 1))/dot_product(d1(:, 1), d1(:, 1)), &
                     dot_product(p - x0, d1(:, 2))/dot_product(d1(:, 2), d1(:, 2))]
                qs(:, i, j) = q
                al(i) = atan2(q(2), q(1))
            end do
            ! unwrap; resample at uniform angle from al(1)
            do i = 2, nth
                al(i) = al(i) + 2*pi*nint((al(i - 1) - al(i))/(2*pi))
            end do
            if (al(nth) < al(1)) then
                ! clockwise in (q1, q2): reverse
                al(1:nth) = al(nth:1:-1)
                qs(:, 1:nth, j) = qs(:, nth:1:-1, j)
            end if
            al(nth + 1) = al(1) + 2*pi
            qs(:, nth + 1, j) = qs(:, 1, j)
            k = 1
            do i = 1, nth
                target = al(1) + 2*pi*(i - 1)/nth
                do while (al(k + 1) < target .and. k < nth)
                    k = k + 1
                end do
                w = (target - al(k))/max(al(k + 1) - al(k), tiny(1.0_dp))
                resampled(:, i, j) = (1 - w)*qs(:, k, j) + w*qs(:, k + 1, j)
            end do
            ph0(j) = al(1)
        end do
        allocate (qb(4, -ntor_out:ntor_out, 0:mpol_out), source=0.0_dp)
        do m = 0, mpol_out
            do n = -ntor_out, ntor_out
                if (m == 0 .and. n < 0) cycle
                w = merge(1.0_dp, 2.0_dp, m == 0 .and. n == 0)/(nth*nphi)
                do j = 1, nphi
                    ph = 2*pi*(j - 1)/(nphi*nfp)
                    do i = 1, nth
                        phase = ph0(j) + 2*pi*(i - 1)/nth
                        arg = m*phase - n*nfp*ph
                        qb(1, n, m) = qb(1, n, m) + w*resampled(1, i, j)*cos(arg)
                        qb(2, n, m) = qb(2, n, m) + w*resampled(1, i, j)*sin(arg)
                        qb(3, n, m) = qb(3, n, m) + w*resampled(2, i, j)*cos(arg)
                        qb(4, n, m) = qb(4, n, m) + w*resampled(2, i, j)*sin(arg)
                    end do
                end do
            end do
        end do
        status = 0

    contains

        subroutine surface(theta, phi, x, dx)
            !! Cartesian point of the boundary at (theta, cylindrical phi) and
            !! its phi derivative.
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: x(3), dx(3)
            real(dp) :: r, z, dr, dz, a
            integer :: mm, nn
            r = 0; z = 0; dr = 0; dz = 0
            do mm = 0, mpol
                do nn = -ntor, ntor
                    a = mm*theta - nn*nfp*phi
                    r = r + boundary(1, nn, mm)*cos(a) + boundary(2, nn, mm)*sin(a)
                    z = z + boundary(3, nn, mm)*cos(a) + boundary(4, nn, mm)*sin(a)
                    dr = dr + nn*nfp*(boundary(1, nn, mm)*sin(a) - boundary(2, nn, mm)*cos(a))
                    dz = dz + nn*nfp*(boundary(3, nn, mm)*sin(a) - boundary(4, nn, mm)*cos(a))
                end do
            end do
            x = [r*cos(phi), r*sin(phi), z]
            dx = [dr*cos(phi) - r*sin(phi), dr*sin(phi) + r*cos(phi), dz]
        end subroutine surface

    end subroutine deck_boundary_to_chart

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

end module milca_deck_chart
