!> Kinetic layers matched to an outer solution (milca_kinetic_layer_match),
!! the three limits of the KIM-1L design: a vacuum-cylinder outer region
!! returns the single-surface torque of milca_kinetic_layer_torque with the
!! layer response of kl_delta, ideal shielding (|Delta| -> infinity) removes
!! flux and torque, and without layer response the outer field
!! psi = -E^{-1} C passes, with no torque from a purely reactive layer.
program test_kinetic_layer_match
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_match, only: klm_solve, KLM_OK
    use milca_kinetic_layer_torque, only: klt_torque
    use milca_kinetic_layer_delta, only: kl_delta
    use milca_kinetic_layer_window, only: kl_profiles_t
    use milca_kinetic_layer_cylinder, only: KL_OK
    implicit none
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    integer :: failures = 0

    call vacuum_cylinder()
    call ideal_shielding()
    call no_layer()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    !> One layer of kl_delta on the profiles of test_kinetic_layer_torque in a
    !> vacuum cylinder, E = -2 |m| and C = 2 |m| r_s B_vac, against klt_torque
    !> at B_vac (#488), which fixes the sign and normalisation of Delta in
    !> Eq. (1). klt_torque takes the binormal wavenumber as k_s = m s/r_s,
    !> kl_delta the exact k_s = (m h_z - k_z r h_theta)/r_s, which differs by
    !> about 3/2 h_theta^2, h_theta ~ r_s/(q R0): the torques agree as R0 grows,
    !> with the difference falling as R0^-2 (1.2e-2 at R0 = 165 cm). A second,
    !> uncoupled layer (diagonal E) without drive carries no flux.
    subroutine vacuum_cylinder()
        real(dp) :: d10, d100
        d10 = matched_difference(1650.0_dp)
        d100 = matched_difference(16500.0_dp)
        print '(a,2es10.2)', 'matched against single-surface torque at R0 x 10, x 100 ', &
            d10, d100
        call check('vacuum cylinder: single-surface torque', d100 < 1e-5_dp)
        call check('vacuum cylinder: difference falls as R0^-2', &
            d10/d100 > 50 .and. d10/d100 < 200)
    end subroutine vacuum_cylinder

    real(dp) function matched_difference(r0) result(rel)
        real(dp), intent(in) :: r0
        type(kl_profiles_t) :: p
        complex(dp) :: d, ipar, psi(2)
        real(dp) :: t(2), t1, rm, bvac
        integer :: st, i, np
        np = 401
        allocate (p%r(np), p%ne(np), p%te(np), p%ti(np), p%q(np), p%er(np))
        do i = 1, np
            p%r(i) = 40 + 30*(i - 1)/real(np - 1, dp)
        end do
        rm = 55
        p%ne = 3e13_dp*exp(-(p%r - 45)/30)
        p%te = 800*exp(-(p%r - 45)/25)
        p%ti = 900*exp(-(p%r - 45)/20)
        p%q = -3*(1 + 0.01_dp*(p%r - rm))
        p%er = -0.3_dp + 0.01_dp*(p%r - rm)
        p%btor = -18000; p%r0 = r0; p%m = 6; p%n = 2; p%zi = 1; p%ai = 2
        call kl_delta(p, rm, 1.5_dp, 3.0_dp, 256, 96, 1, d, ipar, st)
        call check('layer response', st == KL_OK)
        bvac = 2.0_dp
        call klm_solve(reshape([(-12.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
            (-8.0_dp, 0.0_dp)], [2, 2]), [12*rm*bvac*(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)], &
            [d, (0.01_dp, 0.0_dp)], [ipar, (0.0_dp, 0.0_dp)], [rm, 40.0_dp], [6, 4], 2, &
            r0, psi, t, st)
        call check('vacuum cylinder solve', st == KLM_OK)
        t1 = klt_torque(ipar, rm, r0, 6, 2, -1.0_dp, bvac)
        call check('vacuum cylinder: torque nonzero', abs(t1) > 0)
        call check('vacuum cylinder: uncoupled layer without drive', abs(psi(2)) == 0)
        rel = abs(t(1) - t1)/abs(t1)
    end function matched_difference

    !> |Delta| -> infinity at fixed phase, sigma scaling with Delta: psi and
    !> the torque fall like 1/|Delta|.
    subroutine ideal_shielding()
        complex(dp) :: e(2, 2), c(2), d(2), s(2), psi(2), psi_big(2)
        real(dp) :: t(2), t_big(2), r(2)
        integer :: st
        call outer(e, c, r)
        d = [(0.3_dp, 0.7_dp), (-0.2_dp, 0.4_dp)]
        s = 1e10_dp*d
        call klm_solve(e, c, d, s, r, [2, 3], 1, 165.0_dp, psi, t, st)
        call check('shielding: solve', st == KLM_OK)
        call klm_solve(e, c, 1e8_dp*d, 1e8_dp*s, r, [2, 3], 1, 165.0_dp, psi_big, t_big, st)
        call check('shielding: strong layer solve', st == KLM_OK)
        print '(a,2es12.4,a,2es12.4)', 'shielding: |psi| ', abs(psi_big), &
            ', torque ', t_big
        call check('shielding: flux removed', all(abs(psi_big) < 1e-6_dp*abs(psi)))
        call check('shielding: torque removed', all(abs(t_big) < 1e-6_dp*abs(t)))
    end subroutine ideal_shielding

    !> Delta = 0: psi = -E^{-1} C of the coupled outer region; a purely
    !> reactive layer (real Delta, sigma = i c r Delta/(2 k) imaginary) exerts
    !> no torque.
    subroutine no_layer()
        complex(dp) :: e(2, 2), c(2), psi(2), want(2)
        real(dp) :: t(2), r(2), det
        integer :: st
        call outer(e, c, r)
        call klm_solve(e, c, [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)], &
            [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)], r, [2, 3], 1, 165.0_dp, psi, t, st)
        call check('no layer: solve', st == KLM_OK)
        det = real(e(1, 1)*e(2, 2) - e(1, 2)*e(2, 1))
        want = -[e(2, 2)*c(1) - e(1, 2)*c(2), -e(2, 1)*c(1) + e(1, 1)*c(2)]/det
        call check('no layer: psi = -E^-1 C', all(abs(psi - want) <= 1e-13_dp*abs(want)))
        call check('no layer: no torque', all(t == 0))
        call klm_solve(e, c, [(0.05_dp, 0.0_dp), (-0.02_dp, 0.0_dp)], &
            iu*[3e9_dp, -1e9_dp], r, [2, 3], 1, 165.0_dp, psi, t, st)
        call check('reactive layer: solve', st == KLM_OK)
        call check('reactive layer: no torque', all(t == 0) .and. all(abs(psi) > 0))
    end subroutine no_layer

    !> A coupled, symmetric outer region of two surfaces (q = 2, 3) with the
    !> cylinder diagonal -2 m and coupling 0.8, driven at both.
    subroutine outer(e, c, r)
        complex(dp), intent(out) :: e(2, 2), c(2)
        real(dp), intent(out) :: r(2)
        e = reshape([(-4.0_dp, 0.0_dp), (0.8_dp, 0.0_dp), (0.8_dp, 0.0_dp), &
            (-6.0_dp, 0.0_dp)], [2, 2])
        c = [(30.0_dp, 5.0_dp), (-12.0_dp, 20.0_dp)]
        r = [25.0_dp, 40.0_dp]
    end subroutine outer

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_kinetic_layer_match
