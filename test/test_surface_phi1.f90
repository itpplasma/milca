!> Prescribed in-surface potential Phi_1 in the 3D surface solver
!! (surface_transport with phi1) and the local parallel flow:
!!   closure  <b U> of the local flow equals the parallel-flow moment of every
!!            drive, including the phi1 drive, and a constant phi1 acts as
!!            the density rescaling exp(-Z e Phi_1/T) of F (pitch-angle
!!            scattering);
!!   continuity  particle conservation on the surface fixes the local parallel
!!            flow of every drive up to a constant: with psi the toroidal flux
!!            over 2 pi, Boozer B.grad = B^2 (iota d_theta + d_zeta)/(G + iota I)
!!            and d_perp = d_theta - (I/G) d_zeta,
!!              B^2 (iota d_theta + d_zeta)(u/B) + G S = const,
!!            S = T d_perp B (2 f1 + 5 f2)/(Z e B) for the radial magnetic drift
!!            and S = d_perp Phi_1 (f1 + 3 f2/2) - Phi_1 (d_perp B/B)(2 f1 + 3 f2)
!!            for the Phi_1 drive: the radial E x B drift on f_M and the magnetic
!!            drift on the Phi_1 part of dF/dpsi at fixed total energy,
!!            Maxwellian moments of x^2 (1 + xi^2)(f1 + x^2 f2), (f1 + x^2 f2) and
!!            x^2 (1 + xi^2)(-f1 + (1 - x^2) f2); the parallel electric
!!            acceleration conserves density. Independent of the collision
!!            operator.
program test_surface_phi1
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_surface, only: surface_transport, SURFACE_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        COLLOP_OK, COLLOP_LORENTZ, COLLOP_FULL
    use fortnum_fourier, only: fourier_derivative
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), elementary = 1.602176634e-19_dp
    ! SFINCS reference case: Boozer field of three helicities, deuterium.
    real(dp), parameter :: iota = 1.25_dp, g_cov = 6, i_cov = 0.3_dp
    real(dp), parameter :: q = 1/iota, ig = i_cov/g_cov, kappa = iota/(g_cov + iota*i_cov)
    real(dp), parameter :: lnlambda = 136.04778278483653_dp, phi_ref = 2000
    integer, parameter :: nt = 11, nz = 7, nper = 3, n_leg = 24
    character(1024) :: check, path
    logical :: ok

    call get_command_argument(1, check)
    call get_command_argument(2, path)
    select case (trim(check))
    case ('closure'); ok = closure()
    case ('sfincs'); ok = sfincs(trim(path))
    case ('continuity')
        ok = continuity(COLLOP_LORENTZ)
        ok = continuity(COLLOP_FULL) .and. ok
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    subroutine plasma(sp, lnl, grid)
        type(collop_species_t), intent(out) :: sp(1)
        real(dp), intent(out) :: lnl(1, 1)
        type(collop_grid_t), intent(out) :: grid
        integer :: st
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2000*elementary)
        lnl = lnlambda
        call collop_make_grid(6.0_dp, 12, 6, grid, st)
        if (st /= COLLOP_OK) error stop 'grid'
    end subroutine plasma

    function field() result(b)
        real(dp) :: b(nt, nz), th, ze
        integer :: i, j
        do j = 1, nz
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                ze = 2*pi*(j - 1)/(nz*nper)
                b(i, j) = 1 - 0.12_dp*cos(th) - 0.04_dp*cos(th - 3*ze) &
                    + 0.02_dp*sin(2*th - 3*ze)
            end do
        end do
    end function field

    !> Phi_1 in volts: the fixture's Phi1Hat times the reference potential.
    function potential() result(phi)
        real(dp) :: phi(nt, nz), th, ze
        integer :: i, j
        do j = 1, nz
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                ze = 2*pi*(j - 1)/(nz*nper)
                phi(i, j) = phi_ref*1e-4_dp*(cos(th) + 0.6_dp*sin(th) &
                    + 0.5_dp*cos(2*th - 3*ze + 0.3_dp))
            end do
        end do
    end function potential

    logical function closure() result(pass)
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), b(nt, nz), d(1, 1, 3, 3), d_phi(1, 1, 3, 2), &
            u(nt*nz, 1, 1, 3), u_phi(nt*nz, 1, 1, 2), w(nt*nz), err, null, scale
        integer :: its, st, jf
        call plasma(sp, lnl, grid)
        b = field()
        call surface_transport(b, nper, q, ig, kappa, sp, lnl, grid, COLLOP_LORENTZ, n_leg, &
            d, its, st, phi1=potential(), d_phi=d_phi, flow_local=u, flow_local_phi=u_phi)
        pass = st == SURFACE_OK
        if (.not. pass) return
        w = reshape(1/b, [nt*nz])
        err = 0
        do jf = 1, 3
            err = max(err, abs(sum(w*u(:, 1, 1, jf))/sum(w**2) - d(1, 1, 3, jf)) &
                /maxval(abs(d(1, 1, 3, :))))
        end do
        do jf = 1, 2
            err = max(err, abs(sum(w*u_phi(:, 1, 1, jf))/sum(w**2) - d_phi(1, 1, 3, jf)) &
                /maxval(abs(d_phi(1, 1, 3, :))))
        end do
        ! A constant c rescales F by exp(-eps), eps = Z e c/T, and shifts ln n'
        ! by eps T'/T: the response is -eps g1 for force 1 and eps (g1 - g2) for
        ! force 2 (the constant has zero spectral derivative up to rounding).
        call surface_transport(b, nper, q, ig, kappa, sp, lnl, grid, COLLOP_LORENTZ, n_leg, &
            d, its, st, phi1=spread(spread(3.0_dp, 1, nt), 2, nz), d_phi=d_phi)
        scale = 3.0_dp/(sp(1)%t/(sp(1)%z*elementary))
        null = max(maxval(abs(d_phi(1, 1, :, 1) + scale*d(1, 1, :, 1))), &
            maxval(abs(d_phi(1, 1, :, 2) - scale*(d(1, 1, :, 1) - d(1, 1, :, 2))))) &
            /(scale*maxval(abs(d(1, 1, :, 1:2))))
        print '(a,es10.2,a,es10.2)', 'local-to-average closure ', err, &
            ', constant potential against density rescaling ', null
        pass = st == SURFACE_OK .and. err < 1e-12_dp .and. null < 1e-10_dp
    end function closure

    !> u_par(i) = 2 s U(i) f per unit force f, s = q/(1 + (I/G)/q), as the
    !> parallel flow of surface_fluxes; <B u_par> with B = b in tesla. The
    !> SFINCS response is the part linear in Phi_1 of runs at 0, Phi_1, 2 Phi_1.
    logical function sfincs(path) result(pass)
        character(*), intent(in) :: path
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), b(nt, nz), d(1, 1, 3, 3), d_phi(1, 1, 3, 2), &
            u_phi(nt*nz, 1, 1, 2), expected(nt*nz), expected_mean, computed(nt*nz), &
            computed_mean, w(nt*nz), s, err_local, err_mean
        integer :: unit, cases, c, model, jf, its, st
        call plasma(sp, lnl, grid)
        b = field()
        w = reshape(1/b**2, [nt*nz])
        s = q/(1 + ig/q)
        pass = .true.
        open (newunit=unit, file=path, status='old')
        read (unit, *) cases
        do c = 1, cases
            read (unit, *) model, jf
            read (unit, *) expected_mean
            read (unit, *) expected
            call surface_transport(b, nper, q, ig, kappa, sp, lnl, grid, model, n_leg, d, &
                its, st, phi1=potential(), d_phi=d_phi, flow_local_phi=u_phi)
            if (st /= SURFACE_OK) then
                pass = .false.
                cycle
            end if
            computed = 2*s*u_phi(:, 1, 1, jf)
            computed_mean = sum(w*reshape(b, [nt*nz])*computed)/sum(w)
            err_local = maxval(abs(computed - expected))/maxval(abs(expected))
            err_mean = abs(computed_mean - expected_mean)/maxval(abs(expected))
            print '(a,i0,a,i0,a,es12.4,a,es12.4,a)', 'SFINCS collision ', model, ', force ', &
                jf, ': <B u_par> ', computed_mean, ' T m/s (SFINCS ', expected_mean, ')'
            print '(a,es10.2,a,es10.2)', '  local u_par error ', err_local, &
                ', <B u_par> error / max |u_par| ', err_mean
            ! The full operator carries its speed-grid density defect (1e-3).
            pass = pass .and. max(err_local, err_mean) &
                < merge(1e-3_dp, 1e-2_dp, model == COLLOP_LORENTZ)
        end do
        close (unit)
    end function sfincs

    !> Residual of the continuity identity, relative to max |G S|, for the
    !> magnetic drives (force 1, 2) and the phi1 drives (force 1, 2).
    logical function continuity(model) result(pass)
        integer, intent(in) :: model
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), b(nt, nz), d(1, 1, 3, 3), d_phi(1, 1, 3, 2), &
            u(nt*nz, 1, 1, 3), u_phi(nt*nz, 1, 1, 2), dt(nt, nt), dz(nz, nz), &
            bv(nt*nz), dpar_b(nt*nz), dperp_b(nt*nz), dperp_phi(nt*nz), s, tz, err(4)
        integer :: its, st, jf
        call plasma(sp, lnl, grid)
        b = field()
        bv = reshape(b, [nt*nz])
        call surface_transport(b, nper, q, ig, kappa, sp, lnl, grid, model, n_leg, d, its, &
            st, phi1=potential(), d_phi=d_phi, flow_local=u, flow_local_phi=u_phi)
        pass = st == SURFACE_OK
        if (.not. pass) return
        call fourier_derivative(nt, dt)
        call fourier_derivative(nz, dz)
        dz = nper*dz
        s = q/(1 + ig/q)
        tz = sp(1)%t/(sp(1)%z*elementary)
        dperp_b = perp(b, dt, dz)
        dperp_phi = perp(potential(), dt, dz)
        do jf = 1, 2
            err(jf) = residual(2*s*u(:, 1, 1, jf), bv, dt, dz, &
                g_cov*tz*dperp_b/bv*merge(2.0_dp, 5.0_dp, jf == 1))
            err(2 + jf) = residual(2*s*u_phi(:, 1, 1, jf), bv, dt, dz, &
                g_cov*(dperp_phi*merge(1.0_dp, 1.5_dp, jf == 1) &
                - reshape(potential(), [nt*nz])*dperp_b/bv*merge(2.0_dp, 3.0_dp, jf == 1)))
        end do
        print '(a,i0,a,2es10.2,a,2es10.2)', 'continuity, collision ', model, &
            ': magnetic drive ', err(1:2), ', phi1 drive ', err(3:4)
        ! The full operator conserves density only to its speed-element width
        ! (3.7e-3 at 12 elements for the phi1 drive). The parallel electric
        ! acceleration conserves density in the Galerkin form, but the dropped
        ! x = 0 row leaves 1.9e-6, 1.6e-8, 6.0e-11 at 12, 24, 48 elements.
        if (model == COLLOP_LORENTZ) then
            pass = all(err(1:2) < 1e-10_dp) .and. all(err(3:4) < 1e-5_dp)
        else
            pass = all(err < 1e-2_dp)
        end if
    end function continuity

    function perp(f, dt, dz) result(v)
        real(dp), intent(in) :: f(nt, nz), dt(nt, nt), dz(nz, nz)
        real(dp) :: v(nt*nz), w(nt, nz)
        w = matmul(dt, f) - ig*transpose(matmul(dz, transpose(f)))
        v = reshape(w, [nt*nz])
    end function perp

    !> max |r - mean r|/max |G S| with r = B^2 (iota d_theta + d_zeta)(u/B) + G S.
    real(dp) function residual(flow, bv, dt, dz, source)
        real(dp), intent(in) :: flow(nt*nz), bv(nt*nz), dt(nt, nt), dz(nz, nz), &
            source(nt*nz)
        real(dp) :: r(nt*nz)
        ! B^2 D(u/B) = B D u - u D B, in the product form the collocation keeps.
        r = bv*par(flow, dt, dz) - flow*par(bv, dt, dz) - source
        residual = maxval(abs(r - sum(r)/size(r)))/maxval(abs(source))
    end function residual

    !> (iota d_theta + d_zeta) f by collocation.
    function par(f, dt, dz) result(v)
        real(dp), intent(in) :: f(nt*nz), dt(nt, nt), dz(nz, nz)
        real(dp) :: v(nt*nz), w(nt, nz)
        w = reshape(f, [nt, nz])
        w = iota*matmul(dt, w) + transpose(matmul(dz, transpose(w)))
        v = reshape(w, [nt*nz])
    end function par

end program test_surface_phi1
