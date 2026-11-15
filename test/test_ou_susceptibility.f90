!> Susceptibility functions of the Ornstein-Uhlenbeck response
!! (milca_ou_susceptibility) against an independent solution of the same
!! kinetic equation: with f = phi(u) g, phi the unit normal density, the
!! collision operator acts on g as -N with N He_n = n He_n, u acts as
!! u He_n = He_(n+1) + n He_(n-1), and the integral parts of the conserving
!! operators add the projections onto He_2 = u^2 - 1 and He_1 = u, so the
!! equation (-i x2 + i x1 u + N - P) g = u^l is tridiagonal in the Hermite
!! coefficients and I^{kl} = E[u^k g]. Further checks: the collisionless limit
!! I^{00} -> -i Z(x2/(sqrt(2) x1))/(sqrt(2) x1), the pure collisional limit
!! x1 = 0, the symmetry I^{kl} = I^{lk}, both evaluation paths on either side
!! of their switch, and values of the independent closed forms of KiLCA
!! (calc_Imn_array, KAMEL 266210c5).
program test_ou_susceptibility
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ou_susceptibility, only: ou_susceptibility, ou_conserving, OU_OK, &
        OU_SINGULAR, OU_PARTICLES, OU_ENERGY, OU_MOMENTUM, OU_ENERGY_MOMENTUM
    use fortnum_special_faddeeva, only: plasma_dispersion_z
    implicit none
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    integer :: failures = 0

    call hermite_oracle()
    call collisional_limit()
    call collisionless_limit()
    call kilca_values()
    call singular_point()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    subroutine hermite_oracle()
        real(dp), parameter :: x1s(7) = [0.05_dp, 0.3_dp, 1.0_dp, 2.9_dp, 3.1_dp, 7.0_dp, 12.0_dp]
        real(dp), parameter :: x2s(5) = [-5.0_dp, 0.0_dp, 0.2_dp, 40.0_dp, 400.0_dp]
        complex(dp) :: got(0:3, 0:3), want(0:3, 0:3), cons(0:3, 0:3)
        integer :: i, j, model, status
        character(64) :: tag
        do i = 1, size(x1s)
            do j = 1, size(x2s)
                call ou_susceptibility(x1s(i), x2s(j), got, status)
                write (tag, '(a,f6.2,a,f7.1)') 'x1 =', x1s(i), ' x2 =', x2s(j)
                call check('status '//trim(tag), status == OU_OK)
                do model = OU_PARTICLES, OU_ENERGY_MOMENTUM
                    call hermite(x1s(i), x2s(j), model, want)
                    call ou_conserving(got, model, cons)
                    ! The Kummer sum keeps about 1e-15 x1^6 of the scale near
                    ! the switch to the time integral; a conserving correction
                    ! I u v^T I/d amplifies the error of I by max|I|/|d|, which
                    ! is large where the conserved moment is barely damped
                    ! (x1, x2 -> 0).
                    call close_all('Hermite '//trim(tag), cons, want, &
                        2e-12_dp*max(1.0_dp, x1s(i))**3*condition(got, model))
                end do
                call check('symmetric '//trim(tag), &
                    maxval(abs(got - transpose(got))) <= 1e-12_dp*maxval(abs(got)))
            end do
        end do
    end subroutine hermite_oracle

    !> max|I|/|d| over the denominators d of the conserving correction.
    real(dp) function condition(imn, model) result(c)
        complex(dp), intent(in) :: imn(0:, 0:)
        integer, intent(in) :: model
        complex(dp) :: de, dm
        de = 1 - (imn(2, 2) - imn(2, 0) - imn(0, 2) + imn(0, 0))
        dm = 1 - imn(1, 1)
        c = 1
        if (model == OU_ENERGY .or. model == OU_ENERGY_MOMENTUM) &
            c = max(c, maxval(abs(imn))/abs(de))
        if (model == OU_MOMENTUM .or. model == OU_ENERGY_MOMENTUM) &
            c = max(c, maxval(abs(imn))/abs(dm))
    end function condition

    !> Hermite coefficients of g for u^l, l = 0..3, and I^{kl} = E[u^k g].
    subroutine hermite(x1, x2, model, imn)
        real(dp), intent(in) :: x1, x2
        integer, intent(in) :: model
        complex(dp), intent(out) :: imn(0:3, 0:3)
        ! u^0..u^3 in He_0..He_3, and n! E[He_n He_n].
        real(dp), parameter :: power(0:3, 0:3) = reshape([1, 0, 0, 0, 0, 1, 0, 0, &
            1, 0, 1, 0, 0, 3, 0, 1], [4, 4])
        complex(dp), allocatable :: lower(:), diag(:), upper(:), rhs(:, :)
        real(dp), allocatable :: norm(:)
        integer :: n, m, k, l, info
        n = max(600, ceiling(60*x1*x1))
        allocate (lower(n), diag(0:n), upper(n), rhs(0:n, 0:3), norm(0:n))
        ! Row m: (-i x2 + m) g_m + i x1 (g_(m-1) + (m + 1) g_(m+1)).
        do m = 0, n
            diag(m) = cmplx(m, -x2, dp)
        end do
        do m = 1, n
            lower(m) = iu*x1
            upper(m) = iu*x1*m
        end do
        ! Energy: -He_2 <He_2, g>/... with <u^2 - 1, g> = 2 g_2; momentum: g_1.
        if (model == OU_ENERGY .or. model == OU_ENERGY_MOMENTUM) diag(2) = diag(2) - 2
        if (model == OU_MOMENTUM .or. model == OU_ENERGY_MOMENTUM) diag(1) = diag(1) - 1
        rhs = 0
        rhs(0:3, :) = power
        call zgtsv(n + 1, 4, lower, diag, upper, rhs, n + 1, info)
        call check('Hermite solve', info == 0)
        norm(0) = 1
        do m = 1, n
            norm(m) = norm(m - 1)*m
        end do
        do l = 0, 3
            do k = 0, 3
                imn(k, l) = sum(power(0:3, k)*norm(0:3)*rhs(0:3, l))
            end do
        end do
    end subroutine hermite

    !> x1 = 0: g = u^l is a sum of eigenfunctions, I^{00} = i/x2 (x2 /= 0),
    !> I^{11} = 1/(1 - i x2), I^{22} = 2/(2 - i x2) + 1/(-i x2).
    subroutine collisional_limit()
        complex(dp) :: got(0:3, 0:3)
        integer :: status
        call ou_susceptibility(0.0_dp, 2.5_dp, got, status)
        call close_to('x1 = 0: I00', got(0, 0), iu/2.5_dp, 1e-14_dp)
        call close_to('x1 = 0: I11', got(1, 1), 1/(1 - iu*2.5_dp), 1e-14_dp)
        call close_to('x1 = 0: I22', got(2, 2), 2/(2 - iu*2.5_dp) + 1/(-iu*2.5_dp), 1e-14_dp)
        call close_to('x1 = 0: I01', got(0, 1), (0.0_dp, 0.0_dp), 1e-14_dp)
    end subroutine collisional_limit

    !> nu -> 0 at fixed zeta = x2/(sqrt(2) x1): the collisional correction
    !> is of order x1^(-2/3) near the resonance.
    subroutine collisionless_limit()
        real(dp), parameter :: zetas(3) = [0.3_dp, 1.1_dp, 2.5_dp]
        complex(dp) :: got(0:3, 0:3), want
        real(dp) :: x1
        integer :: i, status
        x1 = 3000
        do i = 1, size(zetas)
            call ou_susceptibility(x1, sqrt(2.0_dp)*x1*zetas(i), got, status)
            want = -iu*plasma_dispersion_z(cmplx(zetas(i), 0.0_dp, dp))/(sqrt(2.0_dp)*x1)
            call close_to('collisionless I00', got(0, 0), want, 3*x1**(-2.0_dp/3))
        end do
    end subroutine collisionless_limit

    !> KiLCA calc_Imn_array (closed forms in 1F1(1; x1^2 - i x2 + 3; x1^2)).
    subroutine kilca_values()
        complex(dp) :: got(0:3, 0:3)
        integer :: status
        call ou_susceptibility(0.3_dp, -5.0_dp, got, status)
        call close_to('KiLCA I00 (0.3, -5)', got(0, 0), &
            (1.4276208656880432e-04_dp, -2.0069793989049009e-01_dp), 1e-12_dp)
        call close_to('KiLCA I01 (0.3, -5)', got(0, 1), &
            (-2.3793681094797835e-03_dp, 1.1632331508168140e-02_dp), 1e-12_dp)
        call close_to('KiLCA I02 (0.3, -5)', got(0, 2), &
            (8.8169679743591251e-04_dp, -2.0180341883440164e-01_dp), 1e-12_dp)
    end subroutine kilca_values

    subroutine singular_point()
        complex(dp) :: got(0:3, 0:3)
        integer :: status
        call ou_susceptibility(0.0_dp, 0.0_dp, got, status)
        call check('x1 = x2 = 0 is singular', status == OU_SINGULAR)
    end subroutine singular_point

    subroutine close_all(name, got, want, tol)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: got(0:, 0:), want(0:, 0:)
        real(dp), intent(in) :: tol
        if (maxval(abs(got - want)) <= tol*maxval(abs(want))) return
        failures = failures + 1
        print '(a,a,es10.2,a,es10.2)', 'FAIL ', name, maxval(abs(got - want))/maxval(abs(want)), &
            ' > ', tol
    end subroutine close_all

    subroutine close_to(name, got, want, tol)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: tol
        if (abs(got - want) <= tol*max(abs(want), tiny(1.0_dp))) return
        if (abs(want) == 0 .and. abs(got) <= tol) return
        failures = failures + 1
        print '(a,a,4es24.15)', 'FAIL ', name, got, want
    end subroutine close_to

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_ou_susceptibility
