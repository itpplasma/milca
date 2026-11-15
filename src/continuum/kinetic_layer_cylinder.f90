!> Linear kinetic response of a straight periodic cylinder to a static helical
!! perturbation near its resonant surface, with the Larmor radius treated by
!! integral kernels (KL1 of docs/kinetic-layer-response.md; the model of Markl,
!! Kasilov and Albert, J. Plasma Phys., accepted 2026, local branch KIM-1L).
!!
!! Gaussian units throughout (cm, g, s, erg, statV, G), as in that paper.
!!
!! The background is given on one period [r_lo, r_lo + L) of a radially
!! periodised window around the resonant surface (N equidistant points,
!! endpoint excluded). For one poloidal and toroidal harmonic m, the charge
!! and parallel current densities respond to the potential Phi and to a
!! prescribed radial field B_r through kernels in the Fourier basis
!! exp(i k_l r), k_l = 2 pi l/L, |l| <= M:
!!     K(l, l') = (2 pi/N) sum_j exp(-i (k_l - k_l') r_j)/(8 pi^2)
!!                 sum_species core(k_l, k_l', r_j),
!! with the gyroaverage harmonic only (higher cyclotron harmonics are of
!! relative order omega_E/omega_c). Per species, with b+ = rho_T^2 (k_perp^2
!! + k_perp'^2)/2, bx = rho_T^2 k_perp k_perp', k_perp^2 = k_s^2 + k_r^2,
!! S0 = exp(-b+) I_0(bx), S1 = exp(-b+) I_1(bx) and the susceptibility
!! functions I^{kl}(x1, x2) of milca_ou_susceptibility (x1 = k_par v_T/nu,
!! x2 = -omega_E/nu):
!!     F(I^a, I^b) = I^a (S0 (A1 + A2 (1 - b+)) + A2 bx S1) + A2 S0 I^b/2,
!!     rho-Phi:  (1/lambda_D^2) (-G_B + i v_T^2 k_s/(omega_c nu) F(I00, I02)),
!!     rho-B:   -(1/lambda_D^2) v_T^3/(omega_c nu c) F(I01, I03),
!!     j-Phi:    (1/lambda_D^2) i v_T^3 k_s/(omega_c nu) F(I10, I12),
!!     j-B:     -(1/lambda_D^2) v_T^4/(omega_c nu c) F(I11, I13).
!! The u_par moments follow from counting: the output moment supplies k
!! powers, the source (u_par for B_r, and u_par^2 of the energy factor of the
!! A2 drive) supplies l. G_B = exp(-rho_T^2 (k - k')^2/2) is the gyroaverage of
!! the Boltzmann part at the particle position (paper Eq. 2.41); KIM sets it
!! to 1, available here as boltzmann_gyroaverage = .false.
!! Poisson's equation reads (4 pi K^rhoPhi - diag(k_l^2)) Phi = -4 pi B_r
!! K^rhoB(:, 0) for a radially constant B_r (constant-psi), and the current
!! harmonics are j = K^jPhi Phi + B_r K^jB(:, 0).
module milca_kinetic_layer_cylinder
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_bessel, only: bessel_in_scaled
    use milca_ou_susceptibility, only: ou_susceptibility, ou_conserving, OU_OK, OU_ENERGY
    implicit none
    private
    public :: kl_species_t, kl_background_t, kl_species_point_t, kl_kernels, kl_solve, &
        kl_series, kl_current, kl_point, kl_core, scaled_bessel
    integer, parameter, public :: KL_OK = 0, KL_SUSCEPTIBILITY = 1, KL_SINGULAR = 2, &
        KL_BAD_INPUT = 3

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    real(dp), parameter, public :: kl_charge = 4.80320471257e-10_dp, kl_light = 2.99792458e10_dp
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    !> One species on the window grid: charge number, mass [g], density
    !> [cm^-3], temperature [erg], their radial derivatives and the collision
    !> frequency [1/s] of the Ornstein-Uhlenbeck model.
    type :: kl_species_t
        real(dp) :: z = 0, mass = 0
        real(dp), allocatable :: n(:), t(:), dn(:), dt(:), nu(:)
    end type

    !> The window: equidistant radii r over one period L [cm], the parallel and
    !> perpendicular wavenumbers k_par, k_s [1/cm] of the harmonic, the E x B
    !> frequency omega_E = k_s V_E [1/s], the field B0 [G] and the radial
    !> electric field E_r [statV/cm] of the equilibrium.
    type :: kl_background_t
        real(dp) :: period = 0
        real(dp), allocatable :: r(:), kpar(:), ks(:), omega_e(:), b0(:), er(:)
        type(kl_species_t), allocatable :: species(:)
    end type

    !> Derived quantities of one species at one point.
    type :: kl_species_point_t
        real(dp) :: vt = 0, omega_c = 0, rho = 0, debye2 = 0, a1 = 0, a2 = 0, nu = 0, ks = 0
        complex(dp) :: imn(0:3, 0:3) = 0
    end type

contains

    !> Species a at point j: v_T = sqrt(T/m), omega_c = Z e B0/(m c),
    !> rho_T = v_T/|omega_c|, lambda_D^2 = T/(4 pi Z^2 e^2 n),
    !> A1 = n'/n - Z e E_r/T - 3 T'/(2 T), A2 = T'/T, and the conserving
    !> susceptibility functions of the given model.
    subroutine kl_point(bg, a, j, model, p, status)
        type(kl_background_t), intent(in) :: bg
        integer, intent(in) :: a, j, model
        type(kl_species_point_t), intent(out) :: p
        integer, intent(out) :: status
        complex(dp) :: raw(0:3, 0:3)
        associate (s => bg%species(a))
            p%vt = sqrt(s%t(j)/s%mass)
            p%omega_c = s%z*kl_charge*bg%b0(j)/(s%mass*kl_light)
            p%rho = p%vt/abs(p%omega_c)
            p%debye2 = s%t(j)/(4*pi*(s%z*kl_charge)**2*s%n(j))
            p%a2 = s%dt(j)/s%t(j)
            p%a1 = s%dn(j)/s%n(j) - s%z*kl_charge*bg%er(j)/s%t(j) - 1.5_dp*p%a2
            p%nu = s%nu(j)
            p%ks = bg%ks(j)
            call ou_susceptibility(bg%kpar(j)*p%vt/s%nu(j), -bg%omega_e(j)/s%nu(j), raw, status)
        end associate
        if (status /= OU_OK) then
            status = KL_SUSCEPTIBILITY
            return
        end if
        call ou_conserving(raw, model, p%imn)
        status = KL_OK
    end subroutine kl_point

    !> The four per-species cores at radial wavenumbers k (row) and k' (column):
    !> c(1) rho-Phi, c(2) rho-B, c(3) j-Phi, c(4) j-B.
    pure function kl_core(p, k, kp, boltzmann_gyroaverage) result(c)
        type(kl_species_point_t), intent(in) :: p
        real(dp), intent(in) :: k, kp
        logical, intent(in) :: boltzmann_gyroaverage
        complex(dp) :: c(4)
        real(dp) :: r2, bplus, bcross, s0, s1, gb
        r2 = p%rho**2
        bplus = r2*(2*p%ks**2 + k**2 + kp**2)/2
        bcross = r2*sqrt(p%ks**2 + k**2)*sqrt(p%ks**2 + kp**2)
        call scaled_bessel(bplus, bcross, s0, s1)
        gb = 1
        if (boltzmann_gyroaverage) gb = exp(-r2*(k - kp)**2/2)
        c(1) = (-gb + iu*p%vt**2*p%ks/(p%omega_c*p%nu)*f(p%imn(0, 0), p%imn(0, 2)))/p%debye2
        c(2) = -p%vt**3/(p%omega_c*p%nu*kl_light)*f(p%imn(0, 1), p%imn(0, 3))/p%debye2
        c(3) = iu*p%vt**3*p%ks/(p%omega_c*p%nu)*f(p%imn(1, 0), p%imn(1, 2))/p%debye2
        c(4) = -p%vt**4/(p%omega_c*p%nu*kl_light)*f(p%imn(1, 1), p%imn(1, 3))/p%debye2
    contains
        pure complex(dp) function f(ia, ib)
            complex(dp), intent(in) :: ia, ib
            f = ia*(s0*(p%a1 + p%a2*(1 - bplus)) + p%a2*bcross*s1) + p%a2*s0*ib/2
        end function f
    end function kl_core

    !> s0 = exp(-bplus) I_0(bcross), s1 = exp(-bplus) I_1(bcross), with
    !> 0 <= bcross <= bplus, overflow-free through exp(-x) I_n(x).
    elemental subroutine scaled_bessel(bplus, bcross, s0, s1)
        real(dp), intent(in) :: bplus, bcross
        real(dp), intent(out) :: s0, s1
        s0 = exp(-(bplus - bcross))*bessel_in_scaled(0, bcross)
        s1 = exp(-(bplus - bcross))*bessel_in_scaled(1, bcross)
    end subroutine scaled_bessel

    !> Kernel matrices (rows l, columns l', index l + mmax + 1) on the window.
    subroutine kl_kernels(bg, mmax, model, kphi, kb, kjphi, kjb, status, &
            boltzmann_gyroaverage, only_species)
        type(kl_background_t), intent(in) :: bg
        integer, intent(in) :: mmax, model
        complex(dp), allocatable, intent(out) :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :)
        integer, intent(out) :: status
        logical, intent(in), optional :: boltzmann_gyroaverage
        ! Restrict the sum to one species (diagnostics).
        integer, intent(in), optional :: only_species
        type(kl_species_point_t), allocatable :: pts(:, :)
        complex(dp), allocatable :: acc(:, :, :)
        complex(dp) :: c(4), ph
        real(dp) :: k(-mmax:mmax)
        logical :: gyro
        integer :: n, ns, j, a, l, lp, st
        n = size(bg%r)
        ns = size(bg%species)
        status = KL_BAD_INPUT
        if (n < 1 .or. ns < 1 .or. bg%period <= 0 .or. mmax < 0) return
        gyro = .true.
        if (present(boltzmann_gyroaverage)) gyro = boltzmann_gyroaverage
        allocate (pts(n, ns))
        status = KL_OK
        do a = 1, ns
            do j = 1, n
                call kl_point(bg, a, j, model, pts(j, a), st)
                if (st /= KL_OK) status = st
            end do
        end do
        if (status /= KL_OK) return
        do l = -mmax, mmax
            k(l) = 2*pi*l/bg%period
        end do
        allocate (acc(4, -mmax:mmax, -mmax:mmax))
        acc = 0
        !$omp parallel do default(shared) private(lp, l, j, a, c, ph) schedule(static)
        do lp = -mmax, mmax
            do l = -mmax, mmax
                do j = 1, n
                    c = 0
                    do a = 1, ns
                        if (present(only_species)) then
                            if (a /= only_species) cycle
                        end if
                        c = c + kl_core(pts(j, a), k(l), k(lp), gyro)
                    end do
                    ph = exp(-iu*(k(l) - k(lp))*bg%r(j))
                    acc(:, l, lp) = acc(:, l, lp) + ph*c
                end do
            end do
        end do
        !$omp end parallel do
        ! (2 pi/N) of the periodic trapezoidal rule and 1/(8 pi^2).
        acc = acc/(4*pi*n)
        kphi = acc(1, :, :)
        kb = acc(2, :, :)
        kjphi = acc(3, :, :)
        kjb = acc(4, :, :)
    end subroutine kl_kernels

    !> Poisson's equation for the Fourier coefficients of Phi driven by a
    !> radially constant B_r [G].
    subroutine kl_solve(kphi, kb, period, mmax, br, phi, status)
        complex(dp), intent(in) :: kphi(:, :), kb(:, :)
        real(dp), intent(in) :: period
        integer, intent(in) :: mmax
        complex(dp), intent(in) :: br
        complex(dp), allocatable, intent(out) :: phi(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: a(:, :)
        integer, allocatable :: piv(:)
        integer :: l, info, d
        d = 2*mmax + 1
        a = 4*pi*kphi
        do l = -mmax, mmax
            a(l + mmax + 1, l + mmax + 1) = a(l + mmax + 1, l + mmax + 1) &
                - (2*pi*l/period)**2
        end do
        phi = -4*pi*br*kb(:, mmax + 1)
        allocate (piv(d))
        call zgesv(d, 1, a, d, piv, phi, d, info)
        status = merge(KL_OK, KL_SINGULAR, info == 0)
    end subroutine kl_solve

    !> Parallel current harmonics [statA/cm^2] from Phi and B_r.
    pure function kl_current(kjphi, kjb, phi, br, mmax) result(j)
        complex(dp), intent(in) :: kjphi(:, :), kjb(:, :), phi(:), br
        integer, intent(in) :: mmax
        complex(dp) :: j(size(phi))
        j = matmul(kjphi, phi) + br*kjb(:, mmax + 1)
    end function kl_current

    !> sum_l c_l exp(i k_l r) at the radii r.
    pure function kl_series(c, period, mmax, r) result(v)
        complex(dp), intent(in) :: c(:)
        real(dp), intent(in) :: period, r(:)
        integer, intent(in) :: mmax
        complex(dp) :: v(size(r))
        integer :: l
        v = 0
        do l = -mmax, mmax
            v = v + c(l + mmax + 1)*exp(iu*(2*pi*l/period)*r)
        end do
    end function kl_series

end module milca_kinetic_layer_cylinder
