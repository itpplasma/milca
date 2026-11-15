!> Layer response Delta of the cylinder kinetic layer for the matching of an
!! outer solution (KL closure of docs/kinetic-layer-response.md), Gaussian
!! units.
!!
!! With psi = r B^r and B^r = i k_s A_par, Ampere's law across the layer gives
!! the jump of d psi/dr in terms of the net parallel current I_par through the
!! layer, so for a radially constant B^r (constant psi)
!!     Delta = [d psi/dr]/psi(r_s) = -2 i k_s(r_s) I_par/(c r_s B^r)  [1/cm].
!! I_par = 2 pi int r j_par dr is integrated analytically from the Fourier
!! series of j_par over the core |r - r_s| <= a of the window, where the
!! profiles are those given.
module milca_kinetic_layer_delta
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_cylinder, only: kl_background_t, kl_kernels, kl_solve, &
        kl_current, kl_light, KL_OK
    use milca_kinetic_layer_window, only: kl_profiles_t, kl_window, KLW_OK
    implicit none
    private
    public :: kl_delta
    integer, parameter, public :: KL_BAD_WINDOW = 100

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

contains

    !> Delta [1/cm] and I_par per gauss of B^r [statA/G] of the resonant
    !> surface r_s, from a window of nwin points (even) with core half-width a
    !> and transition width w, Fourier modes |l| <= mmax and conservation
    !> model `model` of the collision operator. status is KL_OK on success,
    !> otherwise the status of the failing stage (window, kernels, solve).
    subroutine kl_delta(p, r_s, a, w, nwin, mmax, model, delta, ipar, status, &
            boltzmann_gyroaverage)
        type(kl_profiles_t), intent(in) :: p
        real(dp), intent(in) :: r_s, a, w
        integer, intent(in) :: nwin, mmax, model
        complex(dp), intent(out) :: delta, ipar
        integer, intent(out) :: status
        logical, intent(in), optional :: boltzmann_gyroaverage
        type(kl_background_t) :: bg
        complex(dp), allocatable :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :), phi(:), j(:)
        complex(dp), parameter :: br = (1.0_dp, 0.0_dp)
        real(dp) :: k, lo, hi
        integer :: l

        delta = 0
        ipar = 0
        status = KL_BAD_WINDOW
        if (mod(nwin, 2) /= 0) return
        call kl_window(p, r_s, a, w, nwin, bg, status)
        if (status /= KLW_OK) return
        call kl_kernels(bg, mmax, model, kphi, kb, kjphi, kjb, status, &
            boltzmann_gyroaverage=boltzmann_gyroaverage)
        if (status /= KL_OK) return
        call kl_solve(kphi, kb, bg%period, mmax, br, phi, status)
        if (status /= KL_OK) return
        j = kl_current(kjphi, kjb, phi, br, mmax)
        lo = r_s - a
        hi = r_s + a
        do l = -mmax, mmax
            k = 2*pi*l/bg%period
            ipar = ipar + j(l + mmax + 1)*radial_moment(k)
        end do
        ipar = 2*pi*ipar
        ! the window is centred on r_s, which is its point nwin/2 + 1
        delta = -2*iu*bg%ks(nwin/2 + 1)*ipar/(kl_light*r_s*br)

    contains

        !> int_lo^hi r exp(i k r) dr
        complex(dp) function radial_moment(k)
            real(dp), intent(in) :: k
            if (k == 0) then
                radial_moment = (hi*hi - lo*lo)/2
            else
                radial_moment = prim(hi, k) - prim(lo, k)
            end if
        end function radial_moment

        complex(dp) function prim(r, k)
            real(dp), intent(in) :: r, k
            prim = exp(iu*k*r)*(r/(iu*k) + 1/(k*k))
        end function prim

    end subroutine kl_delta

end module milca_kinetic_layer_delta
