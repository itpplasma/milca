!> Torque balance of a rotating plasma at one resonant surface with the
!! cylinder kinetic layer response, Gaussian units inside, torques in N m.
!!
!! The layer response sigma(V_E) = I_par/B^r [statA/G] is computed with
!! kl_delta at E x B velocities V_E = -c E_r(r_s)/|B| set by a uniform shift
!! of E_r; the response depends on (Omega, k) only through V_E(r_s) for the
!! broad features (uncertainty-map pre-study). The resonant field follows
!! from the matching to the vacuum cylinder,
!!     |B_m| = |B_vac| / |1 - i sigma/(c r_s s)|,  s = sign(m) sign(B_z),
!! and the toroidal torque of the layer current on the plasma is
!!     tau = -pi R0 r_s n Re(sigma) |B_m|^2/(m c).
!! A viscous stiffness K balances it, K (Omega - Omega_u) = y tau_1(V_E),
!! with V_E = V_u - h Omega and y = (I/I_ref)^2, so that along the path from
!! the unperturbed state the steady states satisfy
!!     y(V) = K (V_u - V)/(h tau_1(V)),
!! and a fold (loss of the steady state on a ramp of y) is a local maximum of
!! y along that path that a later point exceeds.
module milca_kinetic_layer_torque
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_window, only: kl_profiles_t
    use milca_kinetic_layer_delta, only: kl_delta
    use milca_kinetic_layer_cylinder, only: kl_light, KL_OK
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, cubic_spline_eval
    implicit none
    private
    public :: klt_sigma, klt_torque, klt_branch, klt_fold_t

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    real(dp), parameter :: dyncm_per_nm = 1.0e7_dp

    !> One fold: y at the fold, V_E before it and where the ramp lands
    !> (NaN if no later point of the path exceeds it).
    type :: klt_fold_t
        real(dp) :: y, v_before, v_after
    end type klt_fold_t

contains

    !> sigma(V_E) [statA/G] at the E x B velocities ve [km/s] of r_s.
    subroutine klt_sigma(p, r_s, a, w, nwin, mmax, model, ve, sigma, status, &
            boltzmann_gyroaverage)
        type(kl_profiles_t), intent(in) :: p
        real(dp), intent(in) :: r_s, a, w, ve(:)
        integer, intent(in) :: nwin, mmax, model
        complex(dp), intent(out) :: sigma(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: boltzmann_gyroaverage
        type(kl_profiles_t) :: q
        type(cubic_spline_t) :: sp
        real(dp) :: e0(1), de(1), er_s
        complex(dp) :: delta
        integer :: i

        sigma = 0
        call cubic_spline_init(p%r, reshape(p%er, [size(p%r), 1]), sp, status)
        if (status /= 0) return
        call cubic_spline_eval(sp, r_s, e0, de)
        do i = 1, size(ve)
            q = p
            er_s = -ve(i)*1.0e5_dp*abs(p%btor)/kl_light
            q%er = p%er + (er_s - e0(1))
            call kl_delta(q, r_s, a, w, nwin, mmax, model, delta, sigma(i), status, &
                boltzmann_gyroaverage)
            if (status /= KL_OK) return
        end do
    end subroutine klt_sigma

    !> Toroidal torque [N m] of the layer current on the plasma, co-current
    !> positive, at vacuum resonant field bvac [G].
    elemental real(dp) function klt_torque(sigma, r_s, r0, m, n, s, bvac) result(tau)
        complex(dp), intent(in) :: sigma
        real(dp), intent(in) :: r_s, r0, bvac, s
        integer, intent(in) :: m, n
        real(dp) :: bm
        bm = bvac/abs(1 - iu*sigma/(kl_light*r_s*s))
        tau = -pi*r0*r_s*n*real(sigma)*bm**2/(m*kl_light)/dyncm_per_nm
    end function klt_torque


    !> Steady states and folds along the path from the unperturbed V_E = v_u
    !> [km/s]: tau1 [N m] tabulated on the ascending grid ve at the reference
    !> current, viscous stiffness k_visc [N m s/rad], lever h [m] (V_E [km/s]
    !> changes by -1e-3 h per rad/s of Omega). npts points on the path; the
    !> path ends where tau1 changes sign or the grid ends.
    subroutine klt_branch(ve, tau1, v_u, k_visc, h, npts, v, y, folds, locked)
        real(dp), intent(in) :: ve(:), tau1(:), v_u, k_visc, h
        integer, intent(in) :: npts
        real(dp), allocatable, intent(out) :: v(:), y(:)
        type(klt_fold_t), allocatable, intent(out) :: folds(:)
        logical, intent(out) :: locked
        real(dp), allocatable :: t(:)
        real(dp) :: sgn, v_end, runmax
        integer :: i, j, iend
        type(klt_fold_t) :: f

        sgn = sign(1.0_dp, interp(v_u))
        v_end = merge(ve(1), ve(size(ve)), sgn > 0)
        allocate (v(npts), t(npts))
        do i = 1, npts
            v(i) = v_u + (v_end - v_u)*(i - 1)/real(npts - 1, dp)
            t(i) = interp(v(i))
        end do
        iend = npts
        do i = 1, npts
            if (sign(1.0_dp, t(i)) /= sgn .or. t(i) == 0) then
                iend = i - 1
                exit
            end if
        end do
        locked = iend < npts
        v = v(:iend)
        t = t(:iend)
        allocate (y(iend))
        y(1) = 0
        do i = 2, iend
            y(i) = k_visc*(v_u - v(i))*1.0e3_dp/h/t(i)
        end do
        allocate (folds(0))
        runmax = 0
        i = 2
        do while (i < iend)
            if (y(i) >= runmax .and. y(i + 1) < y(i)) then
                f = klt_fold_t(y(i), v(i), ieee_nan())
                j = i + 1
                do while (j <= iend)
                    if (y(j) > y(i)) exit
                    j = j + 1
                end do
                if (j <= iend) f%v_after = v(j)
                folds = [folds, f]
                if (j > iend) exit
                i = j
            end if
            runmax = max(runmax, y(i))
            i = i + 1
        end do

    contains

        real(dp) function interp(x)
            real(dp), intent(in) :: x
            integer :: k
            if (x <= ve(1)) then
                interp = tau1(1)
                return
            end if
            if (x >= ve(size(ve))) then
                interp = tau1(size(ve))
                return
            end if
            k = 1
            do while (ve(k + 1) < x)
                k = k + 1
            end do
            interp = tau1(k) + (tau1(k + 1) - tau1(k))*(x - ve(k))/(ve(k + 1) - ve(k))
        end function interp

    end subroutine klt_branch

    real(dp) function ieee_nan()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        ieee_nan = ieee_value(1.0_dp, ieee_quiet_nan)
    end function ieee_nan

end module milca_kinetic_layer_torque
