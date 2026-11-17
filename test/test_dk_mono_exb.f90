! Independent continuum oracles for the incompressible electric drift.
program test_dk_mono_exb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_mono, only: mono_factor_t, mono_coefficients_t, mono_surface_point, &
        mono_factor_rates, mono_apply, mono_solve, mono_coefficients_3d, MONO_OK, MONO_SINGULAR
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    character(32) :: check
    logical :: ok
    call get_command_argument(1, check)
    select case (trim(check))
    case ('manufactured'); ok = manufactured()
    case ('structure'); ok = structure()
    case ('tangent'); ok = tangent()
    case ('reciprocity'); ok = reciprocity()
    case ('rational'); ok = rational()
    case ('sfincs'); ok = sfincs()
    case default
        error stop 'unknown electric drift check'
    end select
    if (.not. ok) error stop 1
contains
    logical function sfincs() result(ok)
        real(dp), parameter :: kap = 1.25_dp/6.375_dp, q = 0.8_dp
        type(mono_coefficients_t) :: c
        real(dp) :: estar, expected(4), actual(4), nu, c_deflection
        character(1024) :: filename
        integer :: unit, st, j
        call get_command_argument(2, filename)
        open (newunit=unit, file=trim(filename), status='old')
        read (unit, *)
        c_deflection = 3*sqrt(pi)/4*(erf(1.0_dp) &
            - (erf(1.0_dp)-2*exp(-1.0_dp)/sqrt(pi))/2)
        nu = 0.1_dp*c_deflection/6.375_dp
        ok = .true.
        do j = 1, 2
            read (unit, *) estar, expected
            call mono_coefficients_3d(field(15, 9, 3), 3, q, 0.05_dp, kap, nu, 32, &
                c, st, kap*estar)
            if (st /= MONO_OK) error stop 'SFINCS comparison solve'
            actual = [-16*6.375_dp/(sqrt(pi)*1.25_dp**2)*c%d11, &
                -8/(sqrt(pi)*1.25_dp)*c%d31, 8/(sqrt(pi)*1.25_dp)*c%d13, &
                4/(sqrt(pi)*6.375_dp)*c%d33]
            ok = small('SFINCS matrix', maxval(abs((actual-expected)/expected)), &
                1.0_dp, 1e-2_dp) .and. ok
        end do
        close (unit)
    end function sfincs

    logical function rational() result(ok)
        integer, parameter :: nt = 25, nz = 11, n = nt*nz, nl = 8
        type(mono_factor_t) :: f
        real(dp) :: x(n, 0:nl-1), rate(0:nl-1), th, phase
        real(dp), allocatable :: action(:, :)
        integer :: i, j, st
        call mono_surface_point(field(nt, nz, 3), 3, 0.8_dp, 0.05_dp, 0.3_dp, nl, f, st)
        if (st /= MONO_OK) error stop 'rational surface setup'
        x = 0
        do j = 1, nz
            phase = 2*pi*(j - 1)/nz
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                x(i + nt*(j - 1), 0) = cos(12*th - 5*phase)
            end do
        end do
        action = mono_apply(f, 0.1_dp, x)
        ok = small('extra density null mode', maxval(abs(action)), 1.0_dp, 1e-12_dp)
        rate = [(0.1_dp*i*(i + 1)/2, i=0, nl-1)]
        call mono_factor_rates(f, rate, st)
        ok = (st == MONO_SINGULAR) .and. ok
        print '(a,i0)', 'rational factor status ', st
    end function rational

    function field(nt, nz, nper) result(b)
        integer, intent(in) :: nt, nz, nper
        real(dp) :: b(nt, nz), th, ze
        integer :: i, j
        do j = 1, nz
            ze = 2*pi*(j - 1)/(nper*nz)
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                b(i, j) = 1 - 0.12_dp*cos(th) - 0.04_dp*cos(th - nper*ze) &
                    + 0.02_dp*sin(2*th - nper*ze)
            end do
        end do
    end function field

    logical function small(name, error, scale, tolerance) result(ok)
        character(*), intent(in) :: name
        real(dp), intent(in) :: error, scale, tolerance
        ok = abs(error) <= tolerance*scale
        print '(a,es12.4,a,es12.4,a,l1)', name//' error ', abs(error), &
            ' tolerance ', tolerance*scale, ' pass ', ok
    end function small

    logical function manufactured() result(ok)
        integer, parameter :: nt = 9, nz = 7, nper = 3, n = nt*nz, nl = 6
        real(dp), parameter :: q = 0.8_dp, ig = 0.07_dp, kap = 0.3_dp, &
            nu = 0.2_dp, omega = 0.13_dp
        type(mono_factor_t) :: f
        real(dp) :: b(nt, nz), x(n, 0:nl-1), rhs(n, 0:nl-1, 1), rate(0:nl-1)
        real(dp) :: th, chi, db, k, mirror, ve, mean_b2, constraint(1), mu(1)
        real(dp), allocatable :: actual(:, :), solution(:, :, :)
        real(dp), allocatable :: before(:, :)
        character(1024) :: trace
        integer :: i, j, a, l, st, unit
        do j = 1, nz
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                chi = th - 2*pi*(j - 1)/nz
                b(i, j) = 1 - 0.12_dp*cos(th) - 0.04_dp*cos(chi)
            end do
        end do
        mean_b2 = n/sum(1/b**2)
        x = 0
        rhs = 0
        do j = 1, nz
            do i = 1, nt
                a = i + nt*(j - 1)
                th = 2*pi*(i - 1)/nt
                chi = th - 2*pi*(j - 1)/nz
                x(a, 0:2) = [cos(th), sin(th), cos(chi)]
                db = 0.12_dp*sin(th) + 0.04_dp*(1 - q*nper)*sin(chi)
                k = kap*b(i, j)
                mirror = kap*db/2
                ve = omega*b(i, j)**2/mean_b2
                ! Projection of xi d_parallel g - mirror (1-xi^2) d_xi g
                ! + v_E d_perp g - nu L_xi g, with analytic angle derivatives.
                rhs(a, 0, 1) = k*cos(th)/3 - 2*mirror*sin(th)/3 - ve*sin(th)
                rhs(a, 1, 1) = -k*sin(th) - 2*k*(1 - q*nper)*sin(chi)/5 &
                    - 6*mirror*cos(chi)/5 + ve*cos(th) + nu*sin(th)
                rhs(a, 2, 1) = 2*k*cos(th)/3 + 2*mirror*sin(th)/3 &
                    - ve*(1 + ig*nper)*sin(chi) + 3*nu*cos(chi)
                rhs(a, 3, 1) = -3*k*(1 - q*nper)*sin(chi)/5 + 6*mirror*cos(chi)/5
            end do
        end do
        call mono_surface_point(b, nper, q, ig, kap, nl, f, st, omega)
        if (st /= MONO_OK) error stop 'surface setup'
        allocate (actual(n, 0:nl-1), before(n, 0:nl-1))
        actual = mono_apply(f, nu, x)
        ok = small('continuum action', maxval(abs(actual - rhs(:, :, 1))), &
            maxval(abs(rhs)), 1e-12_dp)
        call get_command_argument(2, trace)
        if (len_trim(trace) > 0) then
            f%omega_e = 0
            before = mono_apply(f, nu, x)
            f%omega_e = omega
            open (newunit=unit, file=trim(trace), status='replace')
            write (unit, '(a)') 'theta_rad,analytic_per_m,computed_per_m,zero_electric_per_m'
            do i = 1, nt
                write (unit, '(*(g0,:,","))') 2*pi*(i-1)/nt, rhs(i, 0, 1), &
                    actual(i, 0), before(i, 0)
            end do
            close (unit)
        end if
        rate = [(nu*l*(l + 1)/2, l=0, nl-1)]
        call mono_factor_rates(f, rate, st)
        if (st /= MONO_OK) error stop 'factorization'
        constraint = sum(x(:, 0)*f%p%weight)
        call mono_solve(f, rhs, constraint, solution, mu, st)
        if (st /= MONO_OK) error stop 'manufactured solve'
        ok = small('continuum inverse', maxval(abs(solution(:, :, 1) - x)), &
            1.0_dp, 2e-11_dp) .and. ok
        ok = small('density source', mu(1), 1.0_dp, 1e-12_dp) .and. ok
    end function manufactured

    logical function structure() result(ok)
        type(mono_factor_t) :: f
        real(dp), allocatable :: w(:, :), a(:, :)
        integer :: st, n
        call mono_surface_point(field(9, 7, 3), 3, 0.8_dp, 0.07_dp, 0.3_dp, 6, &
            f, st, 0.13_dp)
        if (st /= MONO_OK) error stop 'structure setup'
        n = f%n
        w = spread(f%p%weight, 2, n)
        a = w*f%exb
        ok = small('weighted skew', maxval(abs(a + transpose(a))), maxval(abs(a)), 1e-13_dp)
        ok = small('constant nullspace', maxval(abs(sum(f%exb, dim=2))), &
            maxval(abs(f%exb)), 1e-13_dp) .and. ok
    end function structure

    function entries(c) result(d)
        type(mono_coefficients_t), intent(in) :: c
        real(dp) :: d(4)
        d = [c%d11, c%d31, c%d13, c%d33]
    end function entries

    logical function tangent() result(ok)
        real(dp), parameter :: omega = 0.013_dp, h = 1e-6_dp
        type(mono_coefficients_t) :: c, dc, cp, cm
        real(dp) :: b(9, 7), fd(4), derivative(4)
        integer :: st
        b = field(9, 7, 3)
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 16, c, st, &
            omega, dc)
        if (st /= MONO_OK) error stop 'tangent solve'
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 16, cp, st, omega+h)
        if (st /= MONO_OK) error stop 'positive difference solve'
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 16, cm, st, omega-h)
        if (st /= MONO_OK) error stop 'negative difference solve'
        fd = (entries(cp) - entries(cm))/(2*h)
        derivative = entries(dc)
        ok = all(abs(fd - derivative) < 2e-6_dp*max(abs(derivative), 1e-5_dp))
        print '(a,4es14.5)', 'implicit derivative ', derivative
        print '(a,4es14.5)', 'central difference  ', fd
    end function tangent

    logical function reciprocity() result(ok)
        type(mono_coefficients_t) :: c, cr, c0, cold
        real(dp) :: b(15, 9), e
        integer :: st
        b = field(15, 9, 3)
        e = 0.03_dp
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 28, c, st, e)
        if (st /= MONO_OK) error stop 'positive field solve'
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 28, cr, st, -e)
        if (st /= MONO_OK) error stop 'negative field solve'
        ok = small('D11 field reversal', c%d11-cr%d11, abs(c%d11), 2e-4_dp)
        ok = small('D33 field reversal', c%d33-cr%d33, abs(c%d33), 2e-4_dp) .and. ok
        ok = small('D13 Casimir', c%d13+cr%d31, max(abs(c%d13), abs(cr%d31)), 2e-4_dp) .and. ok
        ok = small('radial density source', maxval(abs(c%border)), 1.0_dp, 1e-10_dp) .and. ok
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 28, c0, st, 0.0_dp)
        if (st /= MONO_OK) error stop 'zero field solve'
        call mono_coefficients_3d(b, 3, 0.8_dp, 0.07_dp, 0.3_dp, 0.1_dp, 28, cold, st)
        if (st /= MONO_OK) error stop 'default field solve'
        ok = small('zero field limit', maxval(abs(entries(c0)-entries(cold))), &
            maxval(abs(entries(cold))), 1e-13_dp) .and. ok
    end function reciprocity
end program test_dk_mono_exb
