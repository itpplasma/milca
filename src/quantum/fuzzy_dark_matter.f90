!> Fuzzy dark matter: the ground state of the Schrodinger-Poisson system, the
!> solitonic core of a self-gravitating ultralight boson field. Units
!> hbar = m = G = 1 (lengths hbar^2/(G m^3 M_unit) follow by scaling):
!>   i d psi/dt = -(1/2) lap psi + Phi psi,   lap Phi = 4 pi |psi|^2.
!> A stationary state psi(r) exp(-i eps t) with u = Phi - eps solves
!>   psi'' + 2 psi'/r = 2 u psi,   u'' + 2 u'/r = 4 pi psi^2,
!> shot from psi(0) = psi0, psi'(0) = u'(0) = 0 on u(0) by bisection, which
!> separates nodes (psi crossing zero) from growing tails (psi' turning
!> positive). The system is invariant under psi -> lambda^2 psi(lambda r),
!> M -> lambda M. Derived checks on the result: the virial theorem
!> 2 K + W = 0 and eps M = K + 2 W, with K = (1/2) int |grad psi|^2 and
!> W = (1/2) int rho Phi = -int rho M(r)/r.
module milca_fuzzy_dark_matter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, ode_integrate
    implicit none
    private
    public :: soliton_t, soliton_ground_state

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    type :: soliton_t
        !> Mass, eigenvalue, kinetic and potential energy, half-mass radius.
        real(dp) :: mass = 0, eps = 0, kinetic = 0, potential = 0, r_half = 0
        logical :: ok = .false.
    end type

contains

    !> Ground state with central amplitude psi0.
    function soliton_ground_state(psi0) result(s)
        real(dp), intent(in) :: psi0
        type(soliton_t) :: s
        type(ode_solution_t) :: sol
        real(dp) :: lo, hi, mid
        integer :: it, k, kbest
        logical :: node
        ! u(0) = -|u0|: too deep gives a node, too shallow a growing tail
        lo = -10*psi0
        hi = 0
        do it = 1, 200
            mid = (lo + hi)/2
            call shoot(psi0, mid, sol, node)
            if (node) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 4*epsilon(1.0_dp)*abs(lo)) exit
        end do
        call shoot(psi0, (lo + hi)/2, sol, node)
        ! the tail is trustworthy up to its smallest amplitude
        kbest = minloc(abs(sol%y(1, :)), dim=1)
        k = kbest
        s%mass = sol%y(5, k)
        s%kinetic = sol%y(6, k)
        s%potential = sol%y(7, k)
        ! outside the mass Phi = -M/r, so eps = -(u + M/r)
        s%eps = -(sol%y(3, k) + s%mass/sol%t(k))
        do k = 2, kbest
            if (sol%y(5, k) >= s%mass/2) then
                s%r_half = sol%t(k - 1) + (sol%t(k) - sol%t(k - 1))*(s%mass/2 - sol%y(5, k - 1))/ &
                    (sol%y(5, k) - sol%y(5, k - 1))
                exit
            end if
        end do
        ! near the eigenvalue the integration may stop on the decayed tail;
        ! the state is resolved once the tail has fallen far below psi0
        s%ok = abs(sol%y(1, kbest)) < 1e-7_dp*psi0
    end function soliton_ground_state

    subroutine shoot(psi0, u0, sol, node)
        real(dp), intent(in) :: psi0, u0
        type(ode_solution_t), intent(out) :: sol
        logical, intent(out) :: node
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: work
        type(fortnum_status_t) :: status
        real(dp) :: r0
        ! regular series at the origin
        r0 = 1e-6_dp/sqrt(psi0)
        problem%rhs => rhs
        problem%event => turn
        problem%t0 = r0
        problem%t1 = 200/sqrt(psi0)
        problem%y0 = [psi0*(1 + u0*r0**2/3), psi0*2*u0*r0/3, u0 + 4*pi*psi0**2*r0**2/6, &
            4*pi*psi0**2*r0/3, 4*pi*psi0**2*r0**3/3, 0.0_dp, 0.0_dp]
        problem%rtol = 1e-12_dp
        problem%atol = 1e-14_dp*psi0
        call ode_integrate(problem, work, sol, status)
        sol%status = status
        node = sol%y(1, size(sol%t)) <= 0
    end subroutine shoot

    !> y = [psi, psi', u, u', M(r), K(r), W(r)]
    subroutine rhs(r, y, dy, ctx)
        real(dp), intent(in) :: r, y(:)
        real(dp), intent(out) :: dy(:)
        class(*), intent(in), optional :: ctx
        dy(1) = y(2)
        dy(2) = 2*y(3)*y(1) - 2*y(2)/r
        dy(3) = y(4)
        dy(4) = 4*pi*y(1)**2 - 2*y(4)/r
        dy(5) = 4*pi*r**2*y(1)**2
        dy(6) = 2*pi*r**2*y(2)**2
        dy(7) = -4*pi*r*y(1)**2*y(5)
    end subroutine rhs

    !> Zero when psi crosses zero or psi' turns positive.
    function turn(r, y, ctx) result(g)
        real(dp), intent(in) :: r, y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        g = min(y(1), -y(2))
    end function turn

end module milca_fuzzy_dark_matter
