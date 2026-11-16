!> Torque balance at one resonant surface (milca_kinetic_layer_torque): the
!! fold of a manufactured layer torque against its closed form, the torque
!! of a purely reactive layer, and the response scan against a direct
!! kl_delta call at the unshifted E x B velocity.
program test_kinetic_layer_torque
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_torque, only: klt_sigma, klt_torque, klt_branch, klt_fold_t
    use milca_kinetic_layer_delta, only: kl_delta
    use milca_kinetic_layer_window, only: kl_profiles_t
    use milca_kinetic_layer_cylinder, only: kl_light, KL_OK
    implicit none
    integer :: failures = 0

    call fold_closed_form()
    call reactive_layer()
    call scan_matches_delta()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    !> tau_1(V) = tau0/(1 + ((V - Vc)/w)^2) with V_u - Vc = 3 w: along the path
    !> y = K (V_u - V)(1 + u^2)/(h tau0), u = (V - Vc)/w, has its first
    !> maximum at 3 u^2 - 6 u + 1 = 0, u = 1 + sqrt(2/3).
    subroutine fold_closed_form()
        real(dp), parameter :: vc = 0.0_dp, w = 1.0_dp, vu = 3.0_dp, tau0 = 2.0_dp, &
            kv = 0.5_dp, h = 0.25_dp
        real(dp), allocatable :: ve(:), t(:), v(:), y(:)
        type(klt_fold_t), allocatable :: f(:)
        logical :: locked
        real(dp) :: u, yf
        integer :: i
        ve = [(-6.0_dp + 1.0e-3_dp*i, i = 0, 12000)]
        t = tau0/(1 + ((ve - vc)/w)**2)
        call klt_branch(ve, t, vu, kv, h, 30001, v, y, f, locked)
        u = 1 + sqrt(2.0_dp/3)
        yf = kv*(vu - vc - u*w)*1.0e3_dp*(1 + u*u)/(h*tau0)
        print '(a,f9.5,a,f9.5,a,es12.5,a,es12.5)', 'fold at V ', f(1)%v_before, ' (closed form ', &
            vc + u*w, '), y ', f(1)%y, ' (', yf, ')'
        call check('one fold found', size(f) >= 1)
        call check('fold position', abs(f(1)%v_before - (vc + u*w)) < 1.0e-3_dp)
        call check('fold value', abs(f(1)%y - yf) < 1.0e-6_dp*yf)
        call check('path not locked (tau1 keeps its sign)', .not. locked)
    end subroutine fold_closed_form

    !> A purely reactive layer (Re sigma = 0) exerts no torque; for small
    !> |sigma| the field is the vacuum field and tau = -pi R0 r n Re(sigma) B^2/(m c).
    subroutine reactive_layer()
        real(dp) :: t, t0
        t = klt_torque(cmplx(0.0_dp, 1.0e12_dp, dp), 50.0_dp, 165.0_dp, 6, 2, -1.0_dp, 1.0_dp)
        call check('reactive layer has no torque', abs(t) == 0)
        t = klt_torque(cmplx(1.0e3_dp, 0.0_dp, dp), 50.0_dp, 165.0_dp, 6, 2, -1.0_dp, 2.0_dp)
        t0 = -4*atan(1.0_dp)*165.0_dp*50.0_dp*2*1.0e3_dp*4.0_dp/(6*kl_light)/1.0e7_dp
        call check('weak layer torque', abs(t - t0) < 1.0e-9_dp*abs(t0))
    end subroutine reactive_layer

    !> sigma at the E x B velocity of the unshifted profiles equals I_par of
    !> kl_delta on those profiles.
    subroutine scan_matches_delta()
        type(kl_profiles_t) :: p
        complex(dp) :: d, ipar, s(1)
        real(dp) :: ve(1), rm, er_rm
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
        p%btor = -18000; p%r0 = 165; p%m = 6; p%n = 2; p%zi = 1; p%ai = 2
        call kl_delta(p, rm, 1.5_dp, 3.0_dp, 256, 96, 1, d, ipar, st)
        call check('direct solve', st == KL_OK)
        er_rm = -0.3_dp
        ve(1) = -er_rm*kl_light/abs(p%btor)*1.0e-5_dp
        call klt_sigma(p, rm, 1.5_dp, 3.0_dp, 256, 96, 1, ve, s, st)
        call check('scan solve', st == KL_OK)
        call check('scan at the profile V_E equals kl_delta', abs(s(1) - ipar) < 1.0e-9_dp*abs(ipar))
    end subroutine scan_matches_delta

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_kinetic_layer_torque
