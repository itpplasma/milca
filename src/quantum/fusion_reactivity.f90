! Maxwellian fusion reactivities from S factors (N3 of
! docs/nuclear-reactions.md):
!   <sigma v> = c sqrt(8/(pi m_r c^2)) T^(-3/2) int sigma(E) E exp(-E/T) dE,
!   sigma(E) = S(E)/E exp(-B_G/sqrt(E)),
! E and T in keV, S in keV b, <sigma v> in cm^3/s. The integral runs in ln E
! on Gauss-Legendre panels no wider than a set width, with extra panel edges
! around narrow resonances, so the 148 keV resonance of p-11B is resolved.
module milca_fusion_reactivity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_fusion_reactivity_leaves, only: pb11_s_low, pb11_s_mid, pb11_s_high
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: maxwellian_reactivity, s_factor_t, pb11_s_factor, pb11_breaks
    public :: PB11_GAMOW_SQRT_KEV, PB11_MRC2_KEV

    real(dp), parameter :: light_cm = 2.99792458e10_dp, barn_cm2 = 1e-24_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    !> Gamow energy 22.589 MeV (Wang et al. 2026) and reduced mass of p-11B
    !> from CODATA 2018 masses (11B nuclear mass: atomic mass less five
    !> electrons; electron binding below 1 keV neglected).
    real(dp), parameter :: PB11_GAMOW_SQRT_KEV = sqrt(22589.0_dp)
    real(dp), parameter :: PB11_MRC2_KEV = 859604.53_dp
    !> Panel edges around the narrow resonances of the p-11B fit, keV.
    real(dp), parameter :: pb11_breaks(6) = [140.0_dp, 146.0_dp, 150.0_dp, 156.0_dp, &
        400.0_dp, 700.0_dp]

    abstract interface
        !> S factor in keV b at centre-of-mass energy e in keV.
        real(dp) function s_factor_t(e)
            import :: dp
            real(dp), intent(in) :: e
        end function s_factor_t
    end interface

contains

    !> The p-11B S factor of Wang et al. (2026) in keV b, for e <= 10 MeV.
    real(dp) function pb11_s_factor(e) result(s)
        real(dp), intent(in) :: e
        if (e <= 400) then
            call pb11_s_low(e, s)
        else if (e <= 700) then
            call pb11_s_mid(e, s)
        else
            call pb11_s_high(e, s)
        end if
        s = 1000*s
    end function pb11_s_factor

    !> <sigma v> in cm^3/s at ion temperature t (keV) for the S factor sf,
    !> Gamow constant bg = sqrt(E_G/keV), reduced mass mrc2 (keV), over
    !> energies up to emax (keV); breaks are extra panel edges; panels are at
    !> most dz wide in ln E.
    real(dp) function maxwellian_reactivity(t, bg, mrc2, sf, emax, breaks, dz) result(rate)
        real(dp), intent(in) :: t, bg, mrc2, emax, breaks(:), dz
        procedure(s_factor_t) :: sf
        integer, parameter :: nq = 16
        real(dp) :: x(nq), w(nq), edges(size(breaks) + 2), zlo, zhi, z, e, a, b
        real(dp), allocatable :: cuts(:)
        integer :: i, j, k, n, npan
        call gauss_legendre(nq, x, w)
        ! below 1e-3 of the Gamow peak energy the integrand is negligible
        zlo = log(1e-3_dp*(bg**2*t**2/4)**(1.0_dp/3))
        zhi = log(min(emax, 80*t + (bg**2*t**2/4)**(1.0_dp/3)*20))
        n = 0
        edges(1) = zlo
        n = 1
        do i = 1, size(breaks)
            if (log(breaks(i)) > zlo .and. log(breaks(i)) < zhi) then
                n = n + 1
                edges(n) = log(breaks(i))
            end if
        end do
        n = n + 1
        edges(n) = zhi
        call sort(edges(1:n))
        rate = 0
        do i = 1, n - 1
            npan = max(1, ceiling((edges(i + 1) - edges(i))/dz))
            allocate (cuts(npan + 1))
            cuts = [(edges(i) + (edges(i + 1) - edges(i))*k/npan, k=0, npan)]
            do k = 1, npan
                a = cuts(k)
                b = cuts(k + 1)
                do j = 1, nq
                    z = (a + b)/2 + (b - a)/2*x(j)
                    e = exp(z)
                    ! sigma E exp(-E/T) dE with dE = E dz
                    rate = rate + w(j)*(b - a)/2*sf(e)*exp(-bg/sqrt(e) - e/t)*e
                end do
            end do
            deallocate (cuts)
        end do
        rate = light_cm*sqrt(8/(pi*mrc2))*t**(-1.5_dp)*rate*barn_cm2
    end function maxwellian_reactivity

    pure subroutine sort(v)
        real(dp), intent(inout) :: v(:)
        real(dp) :: tmp
        integer :: i, j
        do i = 2, size(v)
            tmp = v(i)
            j = i - 1
            do while (j >= 1)
                if (v(j) <= tmp) exit
                v(j + 1) = v(j)
                j = j - 1
            end do
            v(j + 1) = tmp
        end do
    end subroutine sort

end module milca_fusion_reactivity
