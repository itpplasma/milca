program test_radiation_reaction
    ! The gyro-averaged synchrotron friction against the Landau-Lifshitz
    ! force it reduces: an electron pushed through 40 gyrations in a uniform
    ! magnetic field drifts in momentum and pitch at the rates of
    ! synchrotron_friction, d p/dt = -nu gamma p (1 - xi^2) and
    ! d xi/dt = nu xi (1 - xi^2)/gamma, nu = kappa B^2; without the
    ! reaction term the integrator's energy drift is negligible.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_radiation_reaction_leaves, only: landau_lifshitz_force, synchrotron_friction
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), kappa = 1e-7_dp, bz = 1
    real(dp) :: u(3), u0(3), p0, xi0, p1, xi1, t, dt, gam, pdot, xidot, pm, xim, drift
    integer :: step, nsteps

    p0 = 3
    xi0 = 0.4_dp
    u0 = p0*[sqrt(1 - xi0**2), 0.0_dp, xi0]
    gam = sqrt(1 + p0**2)
    dt = 2*pi*gam/bz/400
    nsteps = 40*400
    ! the integrator's own energy drift must be negligible against the
    ! radiative loss it is to resolve
    u = push(u0, 0.0_dp)
    drift = abs(norm2(u) - p0)
    u = push(u0, kappa)
    t = nsteps*dt
    p1 = norm2(u)
    print '(a, 2es12.4)', 'integrator drift and radiative loss in p:', drift, p0 - p1
    call require(drift < 1e-4_dp*(p0 - p1), 'integrator drift below the radiative loss')
    xi1 = u(3)/p1
    pm = (p0 + p1)/2
    xim = (xi0 + xi1)/2
    call synchrotron_friction(pm, xim, kappa*bz**2, pdot, xidot)
    print '(a, 2es12.4)', 'measured and averaged dp/dt:', (p1 - p0)/t, pdot
    print '(a, 2es12.4)', 'measured and averaged dxi/dt:', (xi1 - xi0)/t, xidot
    call require(abs((p1 - p0)/t/pdot - 1) < 1e-3_dp, 'momentum loss of the gyro-average')
    call require(abs((xi1 - xi0)/t/xidot - 1) < 1e-3_dp, 'pitch drift of the gyro-average')
    print '(a)', 'synchrotron friction is the gyro-average of the Landau-Lifshitz force'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    ! electron (charge -e): the leaf is for charge +e, so the fields enter
    ! with reversed sign; the reaction terms are even in the fields
    function rhs(v, k) result(f)
        real(dp), intent(in) :: v(3), k
        real(dp) :: f(3)
        call landau_lifshitz_force(v(1), v(2), v(3), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            -bz, k, f(1), f(2), f(3))
    end function rhs

    function push(start, k) result(v)
        real(dp), intent(in) :: start(3), k
        real(dp) :: v(3), k1(3), k2(3), k3(3), k4(3)
        v = start
        do step = 1, nsteps
            k1 = rhs(v, k)
            k2 = rhs(v + dt/2*k1, k)
            k3 = rhs(v + dt/2*k2, k)
            k4 = rhs(v + dt*k3, k)
            v = v + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        end do
    end function push

end program test_radiation_reaction
