!> Window background of the cylinder kinetic layer response
!> (milca_kinetic_layer_window): without pressure and with constant q the
!> force balance du/dr = -2 r u/(q^2 R0^2 + r^2) has the solution
!> u (q^2 R0^2 + r^2) = const, so B0 and the wavenumbers are known in closed
!> form; on the core the profiles are those given; the continuation is
!> continuous across the window edge; and with a pressure gradient B0 drops
!> by the diamagnetic amount 8 pi Delta p/(2 B0) to first order.
program test_kinetic_layer_window
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_window, KLW_OK
    use milca_kinetic_layer_cylinder, only: kl_background_t
    implicit none
    real(dp), parameter :: r0 = 165, btor = -18000, qc = -3, rm = 55, a = 0.5_dp, w = 3.0_dp
    real(dp), parameter :: pi = 4*atan(1.0_dp), ev = 1.602176634e-12_dp
    integer :: failures = 0

    call vacuum_constant_q()
    call core_and_edge()
    call diamagnetic_drop()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    subroutine profiles(p, pressure)
        type(kl_profiles_t), intent(out) :: p
        logical, intent(in) :: pressure
        integer, parameter :: n = 401
        integer :: i
        allocate (p%r(n), p%ne(n), p%te(n), p%ti(n), p%q(n), p%er(n))
        do i = 1, n
            p%r(i) = 45 + 20*(i - 1)/real(n - 1, dp)
        end do
        p%ne = 3e13_dp*exp(-(p%r - 45)/30)
        p%te = 800*exp(-(p%r - 45)/25)
        p%ti = 900*exp(-(p%r - 45)/20)
        if (.not. pressure) then
            p%ne = 3e13_dp; p%te = 800; p%ti = 900
        end if
        p%q = qc
        p%er = -0.1_dp + 0.01_dp*(p%r - 55)
        p%btor = btor; p%r0 = r0; p%m = 8; p%n = 3; p%zi = 1; p%ai = 2
    end subroutine profiles

    subroutine vacuum_constant_q()
        type(kl_profiles_t) :: p
        type(kl_background_t) :: bg
        real(dp) :: u0, u, b0z, bth, b0, hz, hth, r
        integer :: st, j
        call profiles(p, .false.)
        call kl_window(p, rm, a, w, 64, bg, st)
        call check('window', st == KLW_OK)
        if (failures > 0) return
        u0 = btor**2*(1 + p%r(1)**2/(r0**2*qc**2))*(qc**2*r0**2 + p%r(1)**2)
        do j = 1, size(bg%r)
            r = bg%r(j)
            if (abs(r - rm) > a) cycle
            u = u0/(qc**2*r0**2 + r**2)
            b0z = sign(sqrt(u/(1 + r**2/(r0**2*qc**2))), btor)
            bth = b0z*r/(qc*r0)
            b0 = sqrt(u)
            hz = b0z/b0; hth = bth/b0
            call close_to('B0', bg%b0(j), b0, 1e-12_dp)
            call close_to('k_par', bg%kpar(j), hth*8/r + hz*3/r0, 1e-9_dp)
            call close_to('k_s', bg%ks(j), hz*8/r - hth*3/r0, 1e-12_dp)
            call close_to('omega_E', bg%omega_e(j), &
                -2.99792458e10_dp*bg%er(j)*(hz*8/r - hth*3/r0)/b0, 1e-12_dp)
        end do
    end subroutine vacuum_constant_q

    subroutine core_and_edge()
        type(kl_profiles_t) :: p
        type(kl_background_t) :: bg
        integer :: st, j, n
        call profiles(p, .true.)
        call kl_window(p, rm, a, w, 256, bg, st)
        do j = 1, size(bg%r)
            if (abs(bg%r(j) - rm) > a) cycle
            call close_to('core n_e', bg%species(1)%n(j), 3e13_dp*exp(-(bg%r(j) - 45)/30), 1e-9_dp)
            call close_to('core T_i', bg%species(2)%t(j), 900*ev*exp(-(bg%r(j) - 45)/20), 1e-9_dp)
            call close_to('core dn/dr', bg%species(1)%dn(j), -3e13_dp/30*exp(-(bg%r(j) - 45)/30), &
                1e-7_dp)
        end do
        ! Last point and the first point of the next period differ by one step.
        n = size(bg%r)
        call check('continuous n_e at the edge', &
            abs(bg%species(1)%n(n) - bg%species(1)%n(1)) < 1e-2_dp*bg%species(1)%n(1))
        call check('continuous B0 at the edge', abs(bg%b0(n) - bg%b0(1)) < 1e-3_dp*bg%b0(1))
    end subroutine core_and_edge

    !> u(r) - u_vacuum(r) ~ -8 pi (p(r) - p(r_1)) to first order in beta.
    subroutine diamagnetic_drop()
        type(kl_profiles_t) :: p, v
        type(kl_background_t) :: bp, bv
        real(dp) :: dp_, du, pr
        integer :: st, j
        call profiles(p, .true.)
        v = p
        v%ne = p%ne(1); v%te = p%te(1); v%ti = p%ti(1)
        call kl_window(p, rm, a, w, 64, bp, st)
        call kl_window(v, rm, a, w, 64, bv, st)
        j = 33
        pr = 3e13_dp*exp(-(bp%r(j) - 45)/30)*(800*exp(-(bp%r(j) - 45)/25) &
            + 900*exp(-(bp%r(j) - 45)/20))*ev
        dp_ = pr - 3e13_dp*1700*ev
        du = bp%b0(j)**2 - bv%b0(j)**2
        call close_to('diamagnetic drop', du, -8*pi*dp_, 5e-2_dp)
    end subroutine diamagnetic_drop

    subroutine close_to(name, got, want, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        if (abs(got - want) <= tol*abs(want)) return
        failures = failures + 1
        print '(a,a,2es24.15)', 'FAIL ', name, got, want
    end subroutine close_to

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_kinetic_layer_window
