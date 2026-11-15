!> Cylinder kinetic layer response (milca_kinetic_layer_cylinder): the
!> overflow-free Bessel products match the direct ones where both apply; the
!> local rho-Phi kernel reduces to the second-order Larmor-radius form of
!> Markl, Kasilov and Albert (Eq. 2.66) with an error of order b^2; on a
!> uniform background the kernels are diagonal and equal to the local core;
!> and far from the resonance (k_par v_T >> nu, omega_E) the solution tends to
!> the ideal one, constant potential on perturbed flux surfaces,
!> Phi = -i E_r B_r/(k_par B0) (Eq. 2.68).
program test_kinetic_layer_cylinder
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_bessel, only: bessel_in
    use milca_kinetic_layer_cylinder, only: kl_species_t, kl_background_t, kl_species_point_t, &
        kl_kernels, kl_solve, kl_point, kl_core, scaled_bessel, KL_OK, kl_charge
    use milca_ou_susceptibility, only: OU_ENERGY
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-9_dp
    real(dp), parameter :: me = 9.1093837e-28_dp, md = 3.3435837768e-24_dp
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    integer :: failures = 0

    call bessel_branches()
    call flr2_limit()
    call uniform_diagonal()
    call ideal_limit()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    subroutine bessel_branches()
        real(dp), parameter :: xs(3) = [610.0_dp, 650.0_dp, 690.0_dp]
        real(dp) :: s0, s1, x
        integer :: i
        do i = 1, 3
            x = xs(i)
            call scaled_bessel(x, x, s0, s1)
            call close_to('large-argument I0', s0, exp(-x)*bessel_in(0, x), 1e-13_dp)
            call close_to('large-argument I1', s1, exp(-x)*bessel_in(1, x), 1e-13_dp)
        end do
        call scaled_bessel(5.0_dp, 3.0_dp, s0, s1)
        call close_to('product I0', s0, exp(-5.0_dp)*bessel_in(0, 3.0_dp), 1e-14_dp)
    end subroutine bessel_branches

    !> Uniform deuterium-electron plasma near an AUG pedestal top.
    subroutine uniform(bg, n, kpar, ks)
        type(kl_background_t), intent(out) :: bg
        integer, intent(in) :: n
        real(dp), intent(in) :: kpar
        real(dp), intent(in), optional :: ks
        integer :: j
        bg%period = 4
        bg%r = [(55 - 2 + 4*(j - 1)/real(n, dp), j=1, n)]
        bg%kpar = [(kpar, j=1, n)]
        bg%ks = [(0.12_dp, j=1, n)]
        if (present(ks)) bg%ks = ks
        bg%b0 = [(2e4_dp, j=1, n)]
        ! E_r = -0.3 statV/cm (about -9 kV/m); V_E = -c E_r/B0.
        bg%er = [(-0.3_dp, j=1, n)]
        bg%omega_e = bg%ks*(-2.99792458e10_dp*bg%er/bg%b0)
        allocate (bg%species(2))
        call species(bg%species(1), -1.0_dp, me, 3e13_dp, 1.0_dp*kev, 2e5_dp, n)
        call species(bg%species(2), 1.0_dp, md, 3e13_dp, 1.0_dp*kev, 2e3_dp, n)
    end subroutine uniform

    subroutine species(s, z, m, dens, t, nu, n)
        type(kl_species_t), intent(out) :: s
        real(dp), intent(in) :: z, m, dens, t, nu
        integer, intent(in) :: n
        integer :: j
        s%z = z; s%mass = m
        s%n = [(dens, j=1, n)]
        s%t = [(t, j=1, n)]
        s%dn = s%n*(-0.5_dp)
        s%dt = s%t*(-0.4_dp)
        s%nu = [(nu, j=1, n)]
    end subroutine species

    !> Taylor expansion of exp(-b) I_0(b), exp(-b) b I_1(b) and (1 - b) to first
    !> order: the rho-Phi core of Eq. (2.66) is
    !>   (1/lambda_D^2) (-1 + i v^2 k_s/(w nu) ((A1 + A2) I00 + A2 I02/2
    !>                   - b ((A1 + 2 A2) I00 + A2 I02/2))).
    subroutine flr2_limit()
        type(kl_background_t) :: bg
        type(kl_species_point_t) :: p
        complex(dp) :: c(4), want
        real(dp), parameter :: bs(2) = [4e-3_dp, 2e-3_dp]
        real(dp) :: err(2), k, b
        integer :: st, i
        call uniform(bg, 1, 1e-3_dp)
        call kl_point(bg, 2, 1, OU_ENERGY, p, st)
        call check('ion point', st == KL_OK)
        do i = 1, 2
            ! k_perp^2 rho^2 = b with k_s included.
            b = bs(i)
            k = sqrt(b/p%rho**2 - p%ks**2)
            c = kl_core(p, k, k, .true.)
            want = (-1 + iu*p%vt**2*p%ks/(p%omega_c*p%nu)*((p%a1 + p%a2)*p%imn(0, 0) &
                + p%a2*p%imn(0, 2)/2 - b*((p%a1 + 2*p%a2)*p%imn(0, 0) + p%a2*p%imn(0, 2)/2))) &
                /p%debye2
            err(i) = abs(c(1) - want)/abs(want)
        end do
        ! Second-order remainder: halving b quarters the error.
        call check('FLR2 limit is second order', err(1) < 1e-4_dp .and. &
            abs(err(1)/err(2) - 4) < 0.1_dp)
    end subroutine flr2_limit

    !> Uniform background: no coupling between radial harmonics, and the
    !> diagonal is the local core over 4 pi.
    subroutine uniform_diagonal()
        type(kl_background_t) :: bg
        type(kl_species_point_t) :: pe, pi_
        complex(dp), allocatable :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :)
        complex(dp) :: want(4)
        real(dp) :: k
        integer :: st
        call uniform(bg, 16, 1e-3_dp)
        call kl_kernels(bg, 3, OU_ENERGY, kphi, kb, kjphi, kjb, st)
        call check('uniform kernels', st == KL_OK)
        if (failures > 0) return
        call kl_point(bg, 1, 1, OU_ENERGY, pe, st)
        call kl_point(bg, 2, 1, OU_ENERGY, pi_, st)
        k = 2*pi*2/bg%period
        want = (kl_core(pe, k, k, .true.) + kl_core(pi_, k, k, .true.))/(4*pi)
        call close_c('diagonal rho-Phi', kphi(6, 6), want(1), 1e-13_dp)
        call close_c('diagonal j-B', kjb(6, 6), want(4), 1e-13_dp)
        call check('no coupling', abs(kphi(6, 5)) <= 1e-13_dp*abs(kphi(6, 6)) .and. &
            abs(kb(2, 7)) <= 1e-13_dp*maxval(abs(kb)))
    end subroutine uniform_diagonal

    !> Far from the resonance the electrons keep the perturbed field lines
    !> equipotential. The remaining deviation is the ion Larmor-radius
    !> correction, of order (rho_i k_s)^2: it falls by 100 when k_s falls by 10.
    subroutine ideal_limit()
        type(kl_background_t) :: bg
        complex(dp), allocatable :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :), phi(:)
        complex(dp) :: br, ideal
        real(dp), parameter :: kss(2) = [0.12_dp, 0.012_dp]
        real(dp) :: dev(2)
        integer :: st, i
        br = (1.0_dp, 0.0_dp)
        do i = 1, 2
            call uniform(bg, 4, 0.3_dp, kss(i))
            call kl_kernels(bg, 0, OU_ENERGY, kphi, kb, kjphi, kjb, st)
            call kl_solve(kphi, kb, bg%period, 0, br, phi, st)
            call check('ideal solve', st == KL_OK)
            ideal = -iu*bg%er(1)*br/(bg%kpar(1)*bg%b0(1))
            dev(i) = abs(phi(1) - ideal)/abs(ideal)
        end do
        print '(a,2es10.2)', 'ideal-limit deviation ', dev
        call check('ideal limit approached', dev(2) < 1e-4_dp .and. &
            abs(dev(1)/dev(2)/100 - 1) < 0.05_dp)
    end subroutine ideal_limit

    subroutine close_to(name, got, want, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        if (abs(got - want) <= tol*abs(want)) return
        failures = failures + 1
        print '(a,a,2es24.15)', 'FAIL ', name, got, want
    end subroutine close_to

    subroutine close_c(name, got, want, tol)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: tol
        if (abs(got - want) <= tol*abs(want)) return
        failures = failures + 1
        print '(a,a,4es24.15)', 'FAIL ', name, got, want
    end subroutine close_c

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a,a)', 'FAIL ', name
    end subroutine check

end program test_kinetic_layer_cylinder
