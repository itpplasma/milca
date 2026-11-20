!> Local 3D transport on a synthetic helical-core `.bc` file (one field
!! period, up-down asymmetric helical harmonic, iota = 1.2):
!!   geometry  the sampled |B|/B0, q, I/G and kappa equal the closed form of
!!             the file's harmonics;
!!   helical   as the helical amplitude delta goes to zero at fixed
!!             collisionality the transport matrix and the bootstrap current
!!             approach the axisymmetric ones like delta^2: delta -> -delta is
!!             the shift phi_B -> phi_B + pi, so the coefficients are even in
!!             delta, and the differences at delta and delta/2 have ratio
!!             4 (1 + O(delta^2)): 3.982 and 3.960 at delta = 0.01, 3.999 and
!!             3.998 at 0.0025 (matrix and bootstrap current);
!!   sfincs    particle flux, kinetic-energy flux and <B u_par> for a density
!!             gradient on the delta = 0.03 file against SFINCS
!!             (test/data/dk_surface_helical_sfincs.md);
!!   sfincs_phi1  the part of <B u_par> and of the local u_par that is linear
!!             in a prescribed Phi_1, for the inductive drive, relative to
!!             <B u_par> of that drive, against SFINCS
!!             (test/data/dk_surface_helical_ephi1_sfincs.md).
program test_dk_surface_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, boozer_surface_t, load_boozer_file, &
        boozer_surface, BOOZER_OK
    use milca_dk_surface_device, only: dk_surface_input_t, dk_surface_result_t, &
        dk_surface_harmonic_t, dk_surface_geometry, dk_surface_profile, DKS_OK
    use milca_dk_surface, only: surface_transport, surface_fluxes, SURFACE_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, COLLOP_OK
    implicit none
    real(dp), parameter :: twopi = 8*atan(1.0_dp), mu0 = 2*twopi*1e-7_dp
    real(dp), parameter :: b0 = 2, r0 = 1.7_dp, a = 0.5_dp, eps = 0.15_dp, iota = 1.2_dp
    real(dp), parameter :: g_cov = 3.4_dp, s_eval = 0.5_dp
    real(dp), parameter :: elementary = 1.602176634e-19_dp
    character(64) :: check
    character(1024) :: path
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('geometry'); ok = geometry()
    case ('helical'); ok = helical()
    case ('sfincs')
        call get_command_argument(2, path)
        ok = sfincs(trim(path))
    case ('sfincs_phi1')
        call get_command_argument(2, path)
        ok = sfincs_phi1(trim(path))
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    !> |B| = b0 (1 - eps r cos theta + delta r (cos(theta - phi) + 0.5 sin(theta - phi))),
    !> r = sqrt(s), on circular surfaces of minor radius a r.
    subroutine write_file(name, delta)
        character(*), intent(in) :: name
        real(dp), intent(in) :: delta
        integer :: unit, k
        real(dp), parameter :: ss(9) = [(0.1_dp*k, k=1, 9)]
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic helical-core test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 1 9 1 0.785 0.5 1.7'
        do k = 1, 9
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') ss(k), iota, twopi*g_cov/mu0, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write (unit, '(2i5,8es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, b0, 0.0_dp
            write (unit, '(2i5,8es24.16)') 1, 0, a*sqrt(ss(k)), 0.0_dp, 0.0_dp, &
                a*sqrt(ss(k)), 0.0_dp, 0.0_dp, -eps*b0*sqrt(ss(k)), 0.0_dp
            write (unit, '(2i5,8es24.16)') 1, 1, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, delta*b0*sqrt(ss(k)), 0.5_dp*delta*b0*sqrt(ss(k))
        end do
        close (unit)
    end subroutine write_file

    logical function geometry() result(pass)
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        real(dp), allocatable :: b(:, :)
        real(dp) :: bb0, q, ig, kappa, th, ph, exact, err
        integer :: st, i, j
        real(dp), parameter :: delta = 0.03_dp
        call write_file('dk_surface_geometry.bc', delta)
        call load_boozer_file('dk_surface_geometry.bc', bf, st)
        pass = st == BOOZER_OK
        if (.not. pass) return
        call boozer_surface(bf, s_eval, surf, st)
        call dk_surface_geometry(surf, 9, 7, b, bb0, q, ig, kappa, st)
        pass = st == DKS_OK
        if (.not. pass) return
        err = 0
        do j = 1, 7
            do i = 1, 9
                th = twopi*(i - 1)/9
                ph = twopi*(j - 1)/7
                exact = 1 - eps*sqrt(s_eval)*cos(th) &
                    + delta*sqrt(s_eval)*(cos(th + ph) + 0.5_dp*sin(th + ph))
                err = max(err, abs(b(i, j) - exact))
            end do
        end do
        err = max(err, abs(bb0/b0 - 1), abs(q*iota + 1), abs(ig), &
            abs(kappa*g_cov/(iota*b0) + 1))
        print '(a,es10.2)', 'helical-core geometry error ', err
        pass = err < 1e-12_dp
    end function geometry

    logical function helical() result(pass)
        type(dk_surface_result_t) :: r(3)
        real(dp) :: dd(3), db(3), ratio_d, ratio_b
        integer :: k
        real(dp), parameter :: delta(3) = [0.01_dp, 0.005_dp, 0.0_dp]
        pass = .true.
        do k = 1, 3
            call solve(delta(k), r(k), pass)
            if (.not. pass) return
        end do
        dd = [(maxval(abs(r(k)%d - r(3)%d)), k=1, 3)]
        db = [(abs(r(k)%bootstrap - r(3)%bootstrap), k=1, 3)]
        ratio_d = dd(1)/dd(2)
        ratio_b = db(1)/db(2)
        print '(a,2es11.3,a,f8.4)', 'matrix change at delta, delta/2 ', dd(1:2), &
            ', ratio ', ratio_d
        print '(a,2es11.3,a,f8.4)', 'bootstrap change [T A/m^2] ', db(1:2), ', ratio ', ratio_b
        pass = abs(ratio_d/4 - 1) < 0.02_dp .and. abs(ratio_b/4 - 1) < 0.02_dp
    end function helical

    !> SFINCS reads the file as its mirror image (psi, iota and n flipped); the
    !> fluxes per psi and <B u_par> agree in sign with this reading.
    logical function sfincs(path) result(pass)
        character(*), intent(in) :: path
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp), allocatable :: b(:, :)
        real(dp) :: bb0, q, ig, kappa, d(1, 1, 3, 3), lnl(1, 1), particle(1), energy(1), &
            flow(1), expected(3), computed(3), err, collisionality
        integer :: st, its, unit, cases, c, model
        call write_file('dk_surface_sfincs.bc', 0.03_dp)
        call load_boozer_file('dk_surface_sfincs.bc', bf, st)
        pass = st == BOOZER_OK
        if (.not. pass) return
        call boozer_surface(bf, s_eval, surf, st)
        call dk_surface_geometry(surf, 11, 7, b, bb0, q, ig, kappa, st)
        pass = st == DKS_OK
        if (.not. pass) return
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2000*elementary)
        call collop_make_grid(6.0_dp, 12, 6, grid, st)
        open (newunit=unit, file=path, status='old')
        read (unit, *) cases
        do c = 1, cases
            read (unit, *) model, collisionality, expected
            lnl = collisionality
            call surface_transport(b, surf%nper, q, ig, kappa, sp, lnl, grid, model, 24, d, &
                its, st)
            pass = pass .and. st == SURFACE_OK
            call surface_fluxes(d, q, ig, sp, [-1/surf%psi_prime], [0.0_dp], 0.0_dp, 0.0_dp, &
                particle, energy, flow, st)
            computed = [particle(1), energy(1), bb0*flow(1)]
            err = maxval(abs(computed/expected - 1))
            print '(a,i0,a,f8.1,a,3es12.4,a,es10.2)', 'SFINCS collision ', model, &
                ', lnlambda ', collisionality, ': Gamma, Q, <B u_par> ', computed, &
                ', relative error ', err
            pass = pass .and. st == SURFACE_OK .and. err < 0.01_dp
        end do
        close (unit)
    end function sfincs

    !> eps = Z e Phi_1/T = 1e-4 [0.4 + cos th + 0.6 sin th + 0.5 cos(2 th + ze
    !> + 0.3)] on the delta = 0.03 file, deuterium at 2 keV; the force
    !> normalisation of the inductive drive cancels in both ratios. The mean
    !> of eps carries the Boltzmann factor of the drive: the varying part of
    !> eps B is a parallel gradient, whose response is a density without flow.
    logical function sfincs_phi1(path) result(pass)
        character(*), intent(in) :: path
        integer, parameter :: nt = 11, nz = 7
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp), allocatable :: b(:, :)
        real(dp) :: bb0, q, ig, kappa, d(1, 1, 3, 3), d_phi(1, 1, 3, 3), lnl(1, 1), &
            u(nt*nz, 1, 1, 3), u_phi(nt*nz, 1, 1, 3), phi1(nt, nz), th, ze, &
            expected_ratio, expected(nt*nz), ratio, local(nt*nz), err_mean, err_local
        integer :: st, its, unit, cases, c, model, i, j
        call write_file('dk_surface_ephi1.bc', 0.03_dp)
        call load_boozer_file('dk_surface_ephi1.bc', bf, st)
        pass = st == BOOZER_OK
        if (.not. pass) return
        call boozer_surface(bf, s_eval, surf, st)
        call dk_surface_geometry(surf, nt, nz, b, bb0, q, ig, kappa, st)
        pass = st == DKS_OK
        if (.not. pass) return
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2000*elementary)
        call collop_make_grid(6.0_dp, 12, 6, grid, st)
        do j = 1, nz
            do i = 1, nt
                th = twopi*(i - 1)/nt
                ze = twopi*(j - 1)/nz
                phi1(i, j) = 2000*1e-4_dp*(0.4_dp + cos(th) + 0.6_dp*sin(th) &
                    + 0.5_dp*cos(2*th + ze + 0.3_dp))
            end do
        end do
        open (newunit=unit, file=path, status='old')
        read (unit, *) cases
        do c = 1, cases
            read (unit, *) model
            read (unit, *) expected_ratio
            read (unit, *) expected
            lnl = 136.0478_dp
            call surface_transport(b, surf%nper, q, ig, kappa, sp, lnl, grid, model, 24, d, &
                its, st, phi1=phi1, d_phi=d_phi, flow_local=u, flow_local_phi=u_phi)
            pass = pass .and. st == SURFACE_OK
            ratio = d_phi(1, 1, 3, 3)/d(1, 1, 3, 3)
            ! <B u_par> = B0 <b U> 2 s f and u_par = 2 s U f, so u_par/<B u_par>
            ! = U/(B0 d(3, 3)).
            local = u_phi(:, 1, 1, 3)/(bb0*d(1, 1, 3, 3))
            err_mean = abs(ratio - expected_ratio)/abs(expected_ratio)
            err_local = maxval(abs(local - expected))/maxval(abs(expected))
            print '(a,i0,a,es12.4,a,es12.4,a,2es10.2)', 'SFINCS collision ', model, &
                ': Phi_1 part of <B u_par> over <B u_par> ', ratio, ' (SFINCS ', &
                expected_ratio, '), errors mean and local ', err_mean, err_local
            pass = pass .and. st == SURFACE_OK .and. max(err_mean, err_local) < 2e-3_dp
        end do
        close (unit)
    end function sfincs_phi1

    subroutine solve(delta, r, pass)
        real(dp), intent(in) :: delta
        type(dk_surface_result_t), intent(out) :: r
        logical, intent(inout) :: pass
        type(boozer_file_t) :: bf
        type(dk_surface_input_t) :: inp
        type(dk_surface_result_t), allocatable :: res(:)
        integer :: st
        call write_file('dk_surface_helical.bc', delta)
        call load_boozer_file('dk_surface_helical.bc', bf, st)
        if (st /= BOOZER_OK) then
            pass = .false.
            return
        end if
        inp%n_theta = 11
        inp%n_zeta = 7
        inp%n_leg = 20
        inp%collisions = 1
        inp%x_max = 5
        inp%nelem = 3
        inp%order = 4
        inp%z = [1.0_dp, -1.0_dp]
        inp%mass = [3.3435837768e-27_dp, 9.1093837015e-31_dp]
        inp%s = [s_eval]
        inp%omega_te = [0.0_dp]
        inp%density = reshape([3e19_dp, 3e19_dp], [1, 2])
        inp%temperature = reshape([2000*elementary, 2000*elementary], [1, 2])
        inp%density_s = reshape([-3e19_dp, -3e19_dp], [1, 2])
        inp%temperature_s = reshape([-2000*elementary, -2000*elementary], [1, 2])
        inp%kappa = reshape([1e-2_dp, 1e-2_dp], [1, 2])
        allocate (inp%phi(0))
        call dk_surface_profile(bf, inp, res, st)
        pass = st == DKS_OK
        if (pass) r = res(1)
    end subroutine solve
end program test_dk_surface_device
