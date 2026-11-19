!> Kinetic layers at N resonant surfaces matched to a toroidal outer solution
!! (KIM-1L), Gaussian units inside, torques in N m.
!!
!! With psi = r B^r at the surfaces r_j, the outer (ideal, toroidal) region
!! gives the jumps of r dpsi/dr across the layers from the reconnected
!! fluxes psi_k and the coil drive C_j,
!!     [r dpsi/dr]_j = sum_k E_jk psi_k + C_j,
!! E the outer tearing matrix (dimensionless) and C [G cm], as in
!! R. Fitzpatrick, Nucl. Fusion 33, 1049 (1993). Each layer closes the jump
!! by its own response, [r dpsi/dr]_j = r_j Delta_j psi_j, with Delta_j
!! [1/cm] and the current per gauss sigma_j = I_par/B^r [statA/G] of
!! kl_delta at the profiles of r_j, so that
!!     sum_k (r_j Delta_j delta_jk - E_jk) psi_k = C_j,                (1)
!! and the toroidal torque of the layer current on the plasma is the
!! single-surface torque of milca_kinetic_layer_torque with B^r = psi/r,
!!     T_j = -pi R0 r_j n Re(sigma_j) |psi_j/r_j|^2/(m_j c).          (2)
!! Without layers (Delta = 0) psi = -E^{-1} C is the field the outer region
!! lets through, and for a vacuum cylinder E_jj = -2 |m_j| and C_j = 2 |m_j|
!! r_j B_vac,j, with which (2) is klt_torque at the vacuum field B_vac,j up to
!! the binormal wavenumber: kl_delta takes the exact k_s of the helical sheet,
!! klt_torque k_s = m s/r_s, a difference of about 3/2 (r_s/(q R0))^2.
module milca_kinetic_layer_match
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_kinetic_layer_cylinder, only: kl_light
    implicit none
    private
    public :: klm_solve
    integer, parameter, public :: KLM_OK = 0, KLM_BAD_INPUT = 1, KLM_SINGULAR = 2

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    real(dp), parameter :: dyncm_per_nm = 1.0e7_dp

contains

    !> psi [G cm] and torque [N m] of N layers from Eq. (1) and (2): e (N x N),
    !> c [G cm], delta [1/cm] and sigma [statA/G] of each layer at r_s [cm]
    !> with mode numbers m and n, major radius r0 [cm].
    subroutine klm_solve(e, c, delta, sigma, r_s, m, n, r0, psi, torque, status)
        complex(dp), intent(in) :: e(:, :), c(:), delta(:), sigma(:)
        real(dp), intent(in) :: r_s(:), r0
        integer, intent(in) :: m(:), n
        complex(dp), intent(out) :: psi(:)
        real(dp), intent(out) :: torque(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: a(:, :), b(:, :)
        integer, allocatable :: piv(:)
        integer :: nr, j, info

        psi = 0
        torque = 0
        status = KLM_BAD_INPUT
        nr = size(c)
        if (any(shape(e) /= [nr, nr]) .or. size(delta) /= nr .or. size(sigma) /= nr .or. &
            size(r_s) /= nr .or. size(m) /= nr .or. size(psi) /= nr .or. &
            size(torque) /= nr) return
        if (any(r_s <= 0) .or. any(m == 0)) return
        a = -e
        do j = 1, nr
            a(j, j) = a(j, j) + r_s(j)*delta(j)
        end do
        b = reshape(c, [nr, 1])
        allocate (piv(nr))
        call zgesv(nr, 1, a, nr, piv, b, nr, info)
        status = KLM_SINGULAR
        if (info /= 0) return
        psi = b(:, 1)
        torque = -pi*r0*r_s*n*real(sigma)*abs(psi/r_s)**2/(m*kl_light)/dyncm_per_nm
        status = KLM_OK
    end subroutine klm_solve

end module milca_kinetic_layer_match
