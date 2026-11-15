program test_nested_flux
    ! The sheared exact equilibrium (Landreman, arXiv:2609.26742, the
    ! eps = 1.2 example of docs/nested-oracles.md) as a flux-coordinate
    ! equilibrium. Checks: the field rebuilt from psi', iota and the stream
    ! function equals the exact field; iota equals the oracle's flux-ratio
    ! transform; the finite-radius magnetic well approaches the near-axis
    ! value from the exact volume and flux densities as rho^2; both Mercier
    ! forms agree, which needs force balance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_flux, only: flux_point_t, flux_space_t, flux_space_init, &
        flux_assemble, flux_energy_terms, flux_normal_spectrum
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK
    use milca_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, s0_solve
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    use milca_ideal_flux_leaves, only: flux_jets_array
    use milca_ideal_mercier, only: mercier_t, flux_mercier
    use milca_nested_flux, only: nested_flux_t
    use milca_nested_oracles, only: nested_oracle_t, NESTED_SHEARED
    use milca_nested_oracle_leaves, only: nested_sheared_axis_well
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(nested_flux_t) :: eq
    character(32) :: mode
    eq%oracle = nested_oracle_t(family=NESTED_SHEARED, eps=1.2_dp, tau=0.0_dp, &
        delta=0.14_dp, saxis=1.4_dp, lam=1.6_dp)
    if (.not. eq%valid()) error stop 'invalid equilibrium'
    mode = 'all'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (mode == 'all' .or. mode == 'field') call field_and_iota()
    if (mode == 'all' .or. mode == 'well') call near_axis_well()
    if (mode == 'study') call study()
    print '(a)', 'PASS nested flux ('//trim(mode)//')'
contains
    subroutine field_and_iota()
        type(flux_point_t) :: pt
        real(dp) :: jac, b(3), g(3, 3), q(3), t(12, 12), v(20), x(3), xj, errb, &
            errf, iota, iref, cur(3), force(3)
        integer :: i, j, ok
        errb = 0; errf = 0
        do j = 1, 5
            do i = 1, 6
                call eq%evaluate(0.15_dp*i, 1.1_dp*j, 0.7_dp*j - 0.3_dp*i, pt)
                call flux_jets_array(pt%xu, pt%xuu, pt%psi1, pt%psi2, pt%iota, &
                    pt%iota1, pt%pp, pt%lam, jac, b, g, q, t)
                call eq%oracle%surface(0.15_dp*i*sqrt(2*eq%oracle%delta), &
                    1.1_dp*j, 0.7_dp*j - 0.3_dp*i, x, xj)
                call eq%oracle%point(x, v, ok)
                if (ok /= 0) error stop 'oracle point'
                errb = max(errb, norm2(b - v(3:5))/norm2(v(3:5)))
                errb = max(errb, maxval(abs(g - reshape(v(12:20), [3, 3]))) &
                    /maxval(abs(v(12:20))))
                cur = [g(3, 2) - g(2, 3), g(1, 3) - g(3, 1), g(2, 1) - g(1, 2)]
                force = [cur(2)*b(3) - cur(3)*b(2), cur(3)*b(1) - cur(1)*b(3), &
                    cur(1)*b(2) - cur(2)*b(1)] - q
                errf = max(errf, norm2(force)/norm2(q))
            end do
        end do
        print '(a,2es10.2)', 'field: rebuilt B and grad B, force residual', errb, errf
        if (errb > 1e-11_dp) error stop 'rebuilt field differs from the exact field'
        if (errf > 1e-10_dp) error stop 'force balance'
        do i = 1, 4
            call eq%evaluate(0.25_dp*i, 0.3_dp, 0.2_dp, pt)
            call eq%oracle%rotation(eq%oracle%delta*(0.25_dp*i)**2, 1, 1e-12_dp, &
                iref, ok)
            if (ok /= 0) error stop 'oracle rotation'
            iota = abs(pt%iota)
            print '(a,f5.2,a,2f12.8)', 'iota at rho', 0.25_dp*i, ': flux chart, oracle', &
                iota, iref
            if (abs(iota - iref) > 1e-10_dp) error stop 'iota differs from the oracle'
        end do
    end subroutine field_and_iota

    subroutine near_axis_well()
        ! With psi_l = rho_o^2/2 the oracle label, dV/dpsi_l and dPhi/dpsi_l are
        ! even in rho_o; their O(rho_o^2) coefficients are the surface means of
        ! half the second rho-derivatives of the densities at the axis. With
        ! Psi = Phi/(2 pi): d2V/dPsi2 = (V_ll Psi_l - V_l Psi_ll)/Psi_l^3.
        integer, parameter :: nq = 64
        type(mercier_t) :: mc
        real(dp) :: v1, v2, t1, t2, vol, volrr, tor, torrr, axis, dev(3), angle, &
            rho
        integer :: i, j, ok
        associate(o => eq%oracle)
            v1 = 0; v2 = 0; t1 = 0; t2 = 0
            do j = 1, nq
                do i = 1, nq
                    call nested_sheared_axis_well(0.0_dp, 2*pi*(i - 1)/nq, &
                        2*pi*(j - 1)/nq, o%eps, o%tau, o%saxis, o%lam, vol, volrr, &
                        tor, torrr)
                    v1 = v1 + vol; v2 = v2 + volrr
                end do
                angle = 2*pi*(j - 1)/nq
                call nested_sheared_axis_well(0.0_dp, angle, 0.0_dp, o%eps, o%tau, &
                    o%saxis, o%lam, vol, volrr, tor, torrr)
                t1 = t1 + tor; t2 = t2 + torrr
            end do
        end associate
        v1 = v1*(2*pi/nq)**2; v2 = v2*(2*pi/nq)**2
        t1 = t1*(2*pi/nq)/(2*pi); t2 = t2*(2*pi/nq)/(2*pi)
        axis = (v2*t1 - v1*t2)/t1**3
        do i = 1, 3
            rho = 0.02_dp*2**(i - 1)
            call flux_mercier(eq, rho, 2, 48, 24, 1.0_dp, mc, ok)
            if (ok /= 0) error stop 'mercier'
            dev(i) = abs(mc%vpp/axis - 1)
            print '(a,f5.2,a,2es15.7,a,es10.2)', 'well rho', rho, &
                ': d2V/dPsi2, near-axis value', mc%vpp, axis, ', (4.1) - D', &
                (mc%total_alt - mc%total)/abs(mc%total)
            if (abs(mc%total_alt - mc%total) > 1e-8_dp*abs(mc%total)) &
                error stop 'Mercier forms disagree'
        end do
        print '(a,3es10.2)', 'well: relative deviation from the axis value', dev
        if (dev(1) > 1e-2_dp) error stop 'well off the near-axis value'
        if (dev(2) > 0.3_dp*dev(3) .or. dev(1) > 0.3_dp*dev(2)) &
            error stop 'well does not approach the axis as rho^2'
    end subroutine near_axis_well

    subroutine study()
        ! The stability study of the sheared case: beta, the transform at the
        ! axis and edge, Mercier and the magnetic well across the radius, and
        ! for each mode family the lowest eigenvalue of the fixed-boundary
        ! problem (gamma = 5/3, kinetic norm rho |xi|^2), its intuitive-form
        ! split and the delta psi spectrum in the straight-field-line angle.
        ! The axisymmetric control eps = 0 keeps all other parameters.
        type(mercier_t) :: mc
        type(nested_flux_t) :: control
        type(flux_point_t) :: pt
        real(dp) :: avg(4), lowest, helical
        integer :: i, ok
        call eq%oracle%averages(48, avg)
        call eq%evaluate(1e-3_dp, 0.0_dp, 0.0_dp, pt)
        print '(a,f8.4,a,f10.6)', 'study: <beta> [%]', 100*avg(4), ', iota axis', &
            abs(pt%iota)
        call eq%evaluate(1.0_dp, 0.0_dp, 0.0_dp, pt)
        print '(a,f10.6)', 'study: iota edge', abs(pt%iota)
        print '(a)', 'study: rho, D_shear, D_curr, D_well, D_geod, D_Merc, d2V/dPsi2'
        do i = 1, 18
            call flux_mercier(eq, 0.05_dp*(i + 1), 2, 64, 32, 1.0_dp, mc, ok)
            if (ok /= 0) error stop 'mercier'
            print '(f6.3,6es12.4)', 0.05_dp*(i + 1), mc%d_shear, mc%d_curr, &
                mc%d_well, mc%d_geod, mc%total, mc%vpp
            if (mc%total >= 0) error stop 'Mercier-stable surface'
            if (mc%vpp <= 0) error stop 'magnetic well'
        end do
        call family(eq, 1, 2, 16, 3.5_dp, 'antiperiodic (odd n)', lowest, helical)
        if (lowest >= 0 .or. helical < 0.9_dp) error stop 'antiperiodic (1, -3) mode'
        call family(eq, 1, 2, 24, 3.5_dp, 'antiperiodic (odd n)', lowest, helical)
        call family(eq, 1, 2, 16, 5.5_dp, 'antiperiodic (odd n)', lowest, helical)
        call family(eq, 0, 2, 16, 3.5_dp, 'periodic (even n)', lowest, helical)
        call family(eq, 0, 2, 16, 5.5_dp, 'periodic (even n)', lowest, helical)
        control%oracle = eq%oracle
        control%oracle%eps = 0
        call control%oracle%averages(48, avg)
        print '(a,f8.4)', 'study: control eps = 0, <beta> [%]', 100*avg(4)
        call family(control, -3, 0, 16, 3.5_dp, 'control eps = 0, n = -3', lowest, &
            helical)
        if (lowest >= 0 .or. helical < 0.9_dp) error stop 'axisymmetric control'
    end subroutine study

    subroutine family(equilibrium, n0, nfp, ne, margin, label, lowest, helical)
        ! Harmonics m = -1..8 and n = n0 mod nfp (only n = n0 when nfp = 0, the
        ! axisymmetric case) with |m iota + n| < margin, cubic hybrid
        ! elements on ne cells. The lowest eigenvalue comes from bisection of
        ! the LDL^H inertia; up to three eigenpairs nearest to it follow by
        ! shift-invert Lanczos.
        type(nested_flux_t), intent(in) :: equilibrium
        integer, intent(in) :: n0, nfp, ne
        real(dp), intent(in) :: margin
        character(*), intent(in) :: label
        real(dp), intent(out) :: lowest
        !> Largest delta psi power fraction in (m, n) = (1, -3) among the pairs.
        real(dp), intent(out) :: helical
        type(flux_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: st
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        type(flux_point_t) :: pt
        integer, allocatable :: ms(:), ns(:)
        complex(dp), allocatable :: start(:), vectors(:, :), omega2(:)
        real(dp), allocatable :: residuals(:)
        real(dp) :: iota, lo, hi, mid, terms(5), total, power(-4:4, -12:12)
        integer :: m, n, j, k, ok, below, nev, top(2)
        call equilibrium%evaluate(0.5_dp, 0.0_dp, 0.0_dp, pt)
        iota = pt%iota
        allocate(ms(0), ns(0))
        do m = -1, 8
            do n = -40, 40
                if (nfp == 0) then
                    if (n /= n0) cycle
                else if (modulo(n - n0, nfp) /= 0) then
                    cycle
                end if
                if (abs(m*iota + n) < margin) then
                    ms = [ms, m]; ns = [ns, n]
                end if
            end do
        end do
        call flux_space_init(space, [(real(j, dp)/ne, j = 0, ne)], 3, ms, ns, &
            max(nfp, 1), 32, merge(24, 1, nfp > 0), ok)
        if (ok /= 0) error stop 'family space'
        call flux_assemble(space, equilibrium, 5.0_dp/3, 1.0_dp, stiffness, mass, st)
        if (st%code /= FORTSPARSE_OK) error stop 'family assembly'
        lo = -1; hi = -1e-9_dp
        call s0_prepare(pencil, stiffness, mass, cmplx(lo, 0, dp), report)
        if (report%below /= 0) error stop 'inertia at -1'
        call s0_prepare(pencil, stiffness, mass, cmplx(hi, 0, dp), report)
        below = int(report%below)
        print '(a,a,a,i3,a,f4.1,a,i4,a,i7,a,i5)', 'study: ', label, ', ne', ne, &
            ', |m iota + n| <', margin, ',', size(ms), ' harmonics,', space%ndof, &
            ' dofs; eigenvalues below -1e-9:', below
        lowest = 0; helical = 0
        if (below <= 0) return
        nev = min(3, below)
        do j = 1, 60
            mid = -sqrt(lo*hi)
            call s0_prepare(pencil, stiffness, mass, cmplx(mid, 0, dp), report)
            if (report%below > 0) then
                hi = mid
            else
                lo = mid
            end if
            if (abs(hi - lo) < 1e-6_dp*abs(lo)) exit
        end do
        call s0_prepare(pencil, stiffness, mass, cmplx(1.01_dp*lo, 0, dp), report)
        allocate(start(space%ndof), vectors(space%ndof, nev), omega2(nev), &
            residuals(nev))
        do j = 1, space%ndof
            start(j) = cmplx(sin(1.3_dp*j), cos(0.7_dp*j), dp)
        end do
        options%max_basis = 32; options%max_restarts = 80
        options%which = EIGEN_NEAREST
        call s0_solve(pencil, start, 1e-10_dp, options, omega2, vectors, &
            residuals, report)
        if (report%status /= EIGEN_OK) error stop 'family eigensolve'
        lowest = minval(real(omega2))
        do k = 1, nev
            call flux_energy_terms(space, equilibrium, 5.0_dp/3, 1.0_dp, &
                vectors(:, k), terms, total, ok)
            if (ok /= 0) error stop 'family energy terms'
            call flux_normal_spectrum(space, equilibrium, vectors(:, k), 0.5_dp, &
                48, 48, 4, 12, power, ok)
            if (ok /= 0) error stop 'family spectrum'
            top = maxloc(power) - [5, 13]
            helical = max(helical, power(1, -3))
            print '(a,es13.5,a,es9.2)', 'study:   omega^2', real(omega2(k)), &
                ', residual', residuals(k)
            print '(a,5f9.4,a,es9.2)', 'study:   2 dW/|2 dW| bend, comp, sound, '// &
                'current, curv:', terms/abs(total), '; sum - total', &
                (sum(terms) - total)/abs(total)
            print '(a,2i4,a,f7.4,a,f7.4,a,f7.4)', 'study:   delta psi at rho = '// &
                '0.5: top (m, n)', top, ' fraction', maxval(power), &
                ', (1, -3)', power(1, -3), ', window', sum(power)
        end do
    end subroutine family
end program test_nested_flux
