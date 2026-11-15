!> Background of the cylinder kinetic layer response on a radially periodised
!! window around a resonant surface (KL1 of docs/kinetic-layer-response.md),
!! built from radial profiles of a large-aspect-ratio tokamak mapped to an
!! effective radius. Gaussian units, temperatures in eV on input.
!!
!! The equilibrium field follows from radial force balance of the straight
!! cylinder with B_theta = r B_z/(q R0): with u = B0^2 = B_z^2 (1 + r^2/(q^2 R0^2)),
!!     du/dr = -2 r u/(q^2 R0^2 + r^2) - 8 pi dp/dr,
!! integrated over the profile radii from the vacuum value at the first one.
!! The harmonic (m, n) has k_par = h_theta m/r + h_z n/R0,
!! k_s = h_z m/r - h_theta n/R0 and omega_E = -c E_r k_s/B0.
!!
!! On the window [r_m - L/2, r_m + L/2), L = 2 (a + w), the profiles n_e,
!! T_e, T_i, E_r and the geometry B0, k_par, k_s, omega_E of the true
!! equilibrium are continued periodically as functions (fortnum_periodize):
!! unchanged on |r - r_m| <= a, blended across bands of width w. Integrating
!! the force balance through the blended profiles instead would add the
!! artificial pressure gradient of the bands to B0 in the core.
!!
!! Collision frequencies [1/s] of the Ornstein-Uhlenbeck model (T in eV, n in
!! cm^-3) with the Coulomb logarithms of the NRL formulary, as KIM uses them
!! (KAMEL after #399), for one ion species with n_i = n_e/Z:
!!     nu_e = (5.8e-6 n_e L_ee + 7.7e-6 n_i Z^2 L_ei)/T_e^(3/2),
!!     nu_i = 3.2e-9 n_e Z^2 L_ei/(A T_e^(1/2) T_i) + 1.4e-7 n_i Z^4 L_ii/(A^(1/2) T_i^(3/2)),
!!     L_ee = 23.5 - ln(n_e^(1/2) T_e^(-5/4)) - (1e-5 + (ln T_e - 2)^2/16)^(1/2),
!!     L_ei = 24 - ln(n_e^(1/2)/T_e),  L_ii = 23 - ln(Z^3 (2 n_i)^(1/2)/T_i^(3/2)).
module milca_kinetic_layer_window
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, cubic_spline_eval
    use fortnum_periodize, only: periodize_points
    use milca_kinetic_layer_cylinder, only: kl_background_t, kl_species_t, kl_light
    implicit none
    private
    public :: kl_profiles_t, kl_window
    integer, parameter, public :: KLW_OK = 0, KLW_BAD_INPUT = 1, KLW_SPLINE = 2

    real(dp), parameter :: pi = 4*atan(1.0_dp), ev = 1.602176634e-12_dp
    real(dp), parameter :: me = 9.1093837015e-28_dp, mp = 1.67262192369e-24_dp

    !> Profiles against the effective radius r [cm]: n_e [cm^-3], T_e, T_i [eV],
    !> q and E_r [statV/cm]; toroidal field btor [G] and major radius r0 [cm];
    !> poloidal and toroidal mode numbers; ion charge and mass numbers.
    type :: kl_profiles_t
        real(dp), allocatable :: r(:), ne(:), te(:), ti(:), q(:), er(:)
        real(dp) :: btor = 0, r0 = 0, zi = 1, ai = 2
        integer :: m = 0, n = 0
    end type

contains

    !> Window of nwin points centred on r_m with core half-width a and
    !> transition width w.
    subroutine kl_window(p, r_m, a, w, nwin, bg, status)
        type(kl_profiles_t), intent(in) :: p
        real(dp), intent(in) :: r_m, a, w
        integer, intent(in) :: nwin
        type(kl_background_t), intent(out) :: bg
        integer, intent(out) :: status
        type(cubic_spline_t) :: sp, spu
        real(dp), allocatable :: y(:, :), us(:, :)
        real(dp) :: period, u, v(5), dv(5), f(5), df(5), g(4), xs(3), ws(3), dws(3), ni
        integer :: j, k, st

        status = KLW_BAD_INPUT
        if (size(p%r) < 4 .or. nwin < 2 .or. a <= 0 .or. w <= 0 .or. p%r0 <= 0) return
        y = reshape([p%ne, p%te, p%ti, p%q, p%er], [size(p%r), 5])
        call cubic_spline_init(p%r, y, sp, st)
        status = KLW_SPLINE
        if (st /= 0) return
        ! B0^2 on the profile radii, then its spline.
        allocate (us(size(p%r), 1))
        u = p%btor**2*(1 + p%r(1)**2/(p%r0**2*p%q(1)**2))
        us(1, 1) = u
        do j = 2, size(p%r)
            call integrate(p%r(j - 1), p%r(j))
            us(j, 1) = u
        end do
        call cubic_spline_init(p%r, us, spu, st)
        if (st /= 0) return

        period = 2*(a + w)
        bg%period = period
        allocate (bg%r(nwin), bg%kpar(nwin), bg%ks(nwin), bg%omega_e(nwin), bg%b0(nwin), &
            bg%er(nwin), bg%species(2))
        call alloc_species(bg%species(1), -1.0_dp, me, nwin)
        call alloc_species(bg%species(2), p%zi, p%ai*mp, nwin)
        do j = 1, nwin
            bg%r(j) = r_m - period/2 + period*(j - 1)/nwin
            call periodize_points(bg%r(j), r_m, a, w, xs, ws, dws)
            v = 0
            dv = 0
            g = 0
            do k = 1, 3
                call cubic_spline_eval(sp, xs(k), f, df)
                v = v + ws(k)*f
                dv = dv + ws(k)*df + dws(k)*f
                g = g + ws(k)*geometry(xs(k))
            end do
            bg%b0(j) = g(1)
            bg%kpar(j) = g(2)
            bg%ks(j) = g(3)
            bg%omega_e(j) = g(4)
            bg%er(j) = v(5)
            ni = v(1)/p%zi
            bg%species(1)%n(j) = v(1)
            bg%species(1)%t(j) = v(2)*ev
            bg%species(1)%dn(j) = dv(1)
            bg%species(1)%dt(j) = dv(2)*ev
            bg%species(2)%n(j) = ni
            bg%species(2)%t(j) = v(3)*ev
            bg%species(2)%dn(j) = dv(1)/p%zi
            bg%species(2)%dt(j) = dv(3)*ev
            call rates(v(1), v(2), ni, v(3), p%zi, p%ai, bg%species(1)%nu(j), &
                bg%species(2)%nu(j))
        end do
        status = KLW_OK

    contains

        !> B0, k_par, k_s and omega_E of the true equilibrium at r.
        function geometry(r) result(gm)
            real(dp), intent(in) :: r
            real(dp) :: gm(4), vv(5), dvv(5), uu(1), du(1), b0z, bth, b0, hz, hth
            call cubic_spline_eval(sp, r, vv, dvv)
            call cubic_spline_eval(spu, r, uu, du)
            b0z = sign(1.0_dp, p%btor)*sqrt(uu(1)/(1 + r**2/(p%r0**2*vv(4)**2)))
            bth = b0z*r/(vv(4)*p%r0)
            b0 = sqrt(b0z**2 + bth**2)
            hz = b0z/b0
            hth = bth/b0
            gm(1) = b0
            gm(2) = hth*p%m/r + hz*p%n/p%r0
            gm(3) = hz*p%m/r - hth*p%n/p%r0
            gm(4) = -kl_light*vv(5)*gm(3)/b0
        end function geometry

        !> du/dr = -2 r u/(q^2 R0^2 + r^2) - 8 pi dp/dr by RK4 from r1 to r2.
        subroutine integrate(r1, r2)
            real(dp), intent(in) :: r1, r2
            integer, parameter :: steps = 16
            real(dp) :: h, r, k1, k2, k3, k4
            integer :: i
            h = (r2 - r1)/steps
            r = r1
            do i = 1, steps
                k1 = rhs(r, u)
                k2 = rhs(r + h/2, u + h*k1/2)
                k3 = rhs(r + h/2, u + h*k2/2)
                k4 = rhs(r + h, u + h*k3)
                u = u + h*(k1 + 2*k2 + 2*k3 + k4)/6
                r = r + h
            end do
        end subroutine integrate

        real(dp) function rhs(r, uu)
            real(dp), intent(in) :: r, uu
            real(dp) :: vv(5), dvv(5), dpdr
            call cubic_spline_eval(sp, r, vv, dvv)
            ! p = n_e T_e + n_i T_i [erg/cm^3]
            dpdr = ev*(dvv(1)*vv(2) + vv(1)*dvv(2) + (dvv(1)*vv(3) + vv(1)*dvv(3))/p%zi)
            rhs = -2*r*uu/(vv(4)**2*p%r0**2 + r*r) - 8*pi*dpdr
        end function rhs

    end subroutine kl_window

    subroutine alloc_species(s, z, mass, n)
        type(kl_species_t), intent(out) :: s
        real(dp), intent(in) :: z, mass
        integer, intent(in) :: n
        s%z = z
        s%mass = mass
        allocate (s%n(n), s%t(n), s%dn(n), s%dt(n), s%nu(n))
    end subroutine alloc_species

    !> Collision frequencies of electrons and ions (header formulas).
    pure subroutine rates(ne, te, ni, ti, z, a, nue, nui)
        real(dp), intent(in) :: ne, te, ni, ti, z, a
        real(dp), intent(out) :: nue, nui
        real(dp) :: lee, lei, lii
        lee = 23.5_dp - log(sqrt(ne)/te**1.25_dp) - sqrt(1e-5_dp + (log(te) - 2)**2/16)
        lei = 24 - log(sqrt(ne)/te)
        lii = 23 - log(z**3/ti**1.5_dp*sqrt(2*ni))
        nue = (5.8e-6_dp*ne*lee + 7.7e-6_dp*ni*z*z*lei)/te**1.5_dp
        nui = 3.2e-9_dp*ne*z*z*lei/(a*sqrt(te)*ti) + 1.4e-7_dp*ni*z**4*lii/(sqrt(a)*ti**1.5_dp)
    end subroutine rates

end module milca_kinetic_layer_window
