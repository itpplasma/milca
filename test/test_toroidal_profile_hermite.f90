program test_toroidal_profile_hermite
    ! Piecewise cubic Hermite profiles on nonuniform knots reproduce cubic
    ! polynomial profiles exactly: values, the iota integral from the axis, and
    ! the shifted iota. The polynomial profile is the independent oracle.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_seed_3d, only: toroidal_profile_t
    implicit none
    real(dp), parameter :: knots(5) = [0.0_dp, 0.07_dp, 0.3_dp, 0.81_dp, 1.0_dp]
    real(dp), parameter :: ci(4) = [0.4_dp, 0.2_dp, -0.35_dp, 0.12_dp]
    real(dp), parameter :: cp(4) = [2.0_dp, -1.5_dp, 0.25_dp, -0.6_dp]
    type(toroidal_profile_t) :: exact, hermite
    real(dp) :: s, iota(2), pressure(2), integral(2), worst
    integer :: j
    exact%iota = ci; exact%pressure = cp
    hermite%iota = [0.0_dp]; hermite%pressure = [0.0_dp]
    call hermite%set_hermite_iota(knots, poly(ci, knots), slope(ci, knots))
    call hermite%set_hermite_pressure(knots, poly(cp, knots), slope(cp, knots))
    if (.not. hermite%valid()) error stop 'Hermite profile rejected'
    worst = compare()
    print '(a,es10.3)', 'largest difference to the cubic profile: ', worst
    if (worst > 1e-14_dp) error stop 'Hermite profile differs from the cubic'
    call exact%shift_iota(0.05_dp); call hermite%shift_iota(0.05_dp)
    worst = compare()
    print '(a,es10.3)', 'after an iota shift: ', worst
    if (worst > 1e-14_dp) error stop 'shifted Hermite iota differs'
    hermite%knots(3) = hermite%knots(2)
    if (hermite%valid()) error stop 'repeated knot accepted'
    print '(a)', 'PASS Hermite toroidal profiles'
contains
    real(dp) function compare() result(worst)
        worst = 0
        do j = 0, 200
            s = j/200.0_dp
            call exact%evaluate(s, iota(1), pressure(1), integral(1))
            call hermite%evaluate(s, iota(2), pressure(2), integral(2))
            worst = max(worst, abs(iota(2) - iota(1)), abs(pressure(2) - pressure(1)), &
                abs(integral(2) - integral(1)))
        end do
    end function

    function poly(c, x) result(values)
        real(dp), intent(in) :: c(:), x(:)
        real(dp) :: values(size(x))
        values = c(1) + x*(c(2) + x*(c(3) + x*c(4)))
    end function

    function slope(c, x) result(values)
        real(dp), intent(in) :: c(:), x(:)
        real(dp) :: values(size(x))
        values = c(2) + x*(2*c(3) + x*3*c(4))
    end function
end program
