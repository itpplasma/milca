program test_strong_field
    ! Quantum synchrotron emission against derived limits: the Bessel
    ! identity int_y K_{5/3} = 2 K_{2/3}(y) - int_y K_{1/3}, which makes the
    ! spectrum classical at delta << 1; the Gaunt factor g(chi) -> 1, the
    ! Landau-Lifshitz power, as 1 - (55 sqrt(3)/16) chi + 48 chi^2 (Erber,
    ! Rev. Mod. Phys. 38, 626, 1966); the emitted power falls below the
    ! classical one.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_strong_field
    implicit none
    real(dp), parameter :: c1 = 55*sqrt(3.0_dp)/16
    real(dp) :: y, lhs, rhs, s1, s2, g
    integer :: i

    do i = 1, 3
        y = 10.0_dp**(i - 3)
        lhs = bessel_k_tail(5.0_dp/3, y)
        rhs = 2*bessel_k_frac(2.0_dp/3, y) - bessel_k_tail(1.0_dp/3, y)
        call require(abs(lhs/rhs - 1) < 1e-10_dp, 'Bessel identity of the classical spectrum')
    end do

    s1 = (1 - quantum_gaunt(1e-4_dp, 400))/1e-4_dp
    s2 = (1 - quantum_gaunt(2e-4_dp, 400))/2e-4_dp
    print '(a, 2f12.7)', 'first quantum correction and 55 sqrt(3)/16:', 2*s1 - s2, c1
    call require(abs((2*s1 - s2)/c1 - 1) < 1e-5_dp, 'first quantum correction')
    print '(a, f10.4)', 'second coefficient:', (s1 - s2)/1e-4_dp
    call require(abs((s1 - s2)/1e-4_dp/48 - 1) < 1e-2_dp, 'second quantum correction')
    g = quantum_gaunt(1.0_dp, 400)
    print '(a, f10.6)', 'g(chi = 1):', g
    call require(g > 0 .and. g < 1, 'quantum emission below the classical power')
    print '(a)', 'quantum synchrotron emission matches its derived limits'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_strong_field
