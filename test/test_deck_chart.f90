program test_deck_chart
    !! milca_deck_chart on a helical boundary (nfp 5, 5 cm helical shift of a
    !! 0.13 x 0.16 m ellipse, like VMEC++'s CTH-like deck): the chart of its
    !! m = 0 curve initialises, the boundary in chart coordinates maps back to
    !! points of the original surface (round trip), and in the chart the
    !! sections stay centred: the helical swing seen in the cylindrical
    !! planes (5 cm) is gone.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_chart, only: curve_chart_t
    use milca_deck_chart, only: deck_axis_chart, deck_boundary_to_chart
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 0.78_dp, d = 0.05_dp, &
                           a = 0.13_dp, b = 0.16_dp
    integer, parameter :: nfp = 5, ntor = 2, mpol = 2, nout = 8
    real(dp) :: boundary(4, -ntor:ntor, 0:mpol), x(3), d1(3, 3), d2(3, 3, 3), &
                d3(3, 3, 3, 3), ph, al, q(2), p(3), err, swing, rr, rz(2), dist, best
    real(dp), allocatable :: qb(:, :, :)
    type(curve_chart_t) :: chart
    integer :: st, j, i, m, n, failures

    failures = 0
    boundary = 0
    ! R = r0 + d cos(nfp phi) + a cos theta, Z = d sin(nfp phi) + b sin theta
    boundary(1, 0, 0) = r0
    boundary(1, -1, 0) = d
    boundary(4, -1, 0) = d
    boundary(1, 0, 1) = a
    boundary(4, 0, 1) = b
    call deck_axis_chart(boundary, ntor, nfp, chart, st)
    call check(st == 0, 'the chart of the m = 0 curve initialises')
    call deck_boundary_to_chart(boundary, ntor, nfp, chart, 16, nout, qb, st)
    call check(st == 0, 'the boundary converts to chart coordinates')
    ! round trip: chart points lie on the surface
    err = 0
    swing = 0
    do j = 1, 7
        ph = 2*pi*(j - 1)/(7*nfp)
        call chart%jet([0.0_dp, 0.0_dp, ph], x, d1, d2, d3)
        rz = 0
        do i = 1, 64
            al = 2*pi*(i - 1)/64
            q = 0
            do m = 0, 16
                do n = -nout, nout
                    q(1) = q(1) + qb(1, n, m)*cos(m*al - n*nfp*ph) + qb(2, n, m)*sin(m*al - n*nfp*ph)
                    q(2) = q(2) + qb(3, n, m)*cos(m*al - n*nfp*ph) + qb(4, n, m)*sin(m*al - n*nfp*ph)
                end do
            end do
            rz = rz + q/64
            call chart%jet([q(1), q(2), ph], p, d1, d2, d3)
            ! distance to the surface: the section at the point's cylindrical
            ! angle is an ellipse; its level-set defect times the local scale
            best = section_defect(p)
            err = max(err, best)
        end do
        swing = max(swing, norm2(rz))
    end do
    print '(a,es10.2,a,es10.2)', 'round trip distance (m)', err, &
        '  section centroid off the curve (m)', swing
    call check(err < 1e-5_dp, 'chart points lie on the surface')
    call check(swing < 0.2_dp*d, 'the sections stay centred in the chart')
    if (failures > 0) error stop 'deck chart'
    print '(a)', 'PASS deck chart'

contains

    real(dp) function section_defect(p) result(dist)
        !! approximate distance of p from the surface: on the plane phi of p
        !! the section is the ellipse ((R - Rc)/a)^2 + ((Z - Zc)/b)^2 = 1
        real(dp), intent(in) :: p(3)
        real(dp) :: phi, r, u, v, g
        phi = atan2(p(2), p(1))
        r = hypot(p(1), p(2))
        u = (r - r0 - d*cos(nfp*phi))/a
        v = (p(3) - d*sin(nfp*phi))/b
        g = sqrt(u**2 + v**2)
        ! radial distance along the ray from the centre of the ellipse
        dist = abs(g - 1)*hypot(u*a, v*b)/max(g, tiny(1.0_dp))
    end function section_defect

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (ok) return
        print '(2a)', 'FAIL ', what
        failures = failures + 1
    end subroutine check

end program test_deck_chart
