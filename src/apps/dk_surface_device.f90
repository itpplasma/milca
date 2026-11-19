!> Local 3D transport on the surfaces of a Boozer `.bc` file: the transport
!! matrix of surface_transport, bootstrap current and, for a prescribed
!! in-surface potential, its parallel current, averaged and local. The file
!! may be up-down asymmetric and have one field period with a helical axis
!! (helical-core tokamak, #222).
!!
!! Input (dk_surface_read_input), free format:
!!   path of the .bc file
!!   n_theta n_zeta n_leg collisions x_max nelem order
!!   n_species, then per species: Z  mass [kg]
!!   n_surfaces, then per surface: s  Omega_tE [rad/s], and per species
!!     n [m^-3]  T [J]  dn/ds  dT/ds  kappa [1/m] (NEO-2 collisionality)
!!   optional: n_phi, then per harmonic: m  n  Phi_c  Phi_s [V], with
!!     Phi_1 = sum Phi_c cos(m theta - n N_p phi_B) + Phi_s sin(...),
!!     n without the field period N_p as in the .bc file.
!! n_theta and n_zeta are odd; zeta spans one field period. Omega_tE is
!! dPhi_0/dpsi_pol (NEO-2), so dPhi_0/dpsi = iota Omega_tE with psi the
!! toroidal flux over 2 pi. collisions: 1 pitch-angle scattering, 3 full
!! linearised operator (milca_dk_collop).
module milca_dk_surface_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_boozer_file, only: boozer_file_t, boozer_surface_t, boozer_surface, &
        boozer_point, boozer_point_t, BOOZER_OK
    use milca_dk_surface, only: surface_transport, surface_fluxes, SURFACE_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_vth, collop_nuhat, COLLOP_OK
    use milca_dk_mono_leaves, only: dk_surface_force_point
    implicit none
    private
    public :: dk_surface_input_t, dk_surface_result_t, dk_surface_harmonic_t
    public :: dk_surface_read_input, dk_surface_geometry, dk_surface_profile
    integer, parameter, public :: DKS_OK = 0, DKS_BAD_INPUT = 1, DKS_OUTSIDE = 2, &
        DKS_SOLVE = 3
    real(dp), parameter :: twopi = 8*atan(1.0_dp), pi = twopi/2
    real(dp), parameter :: elementary = 1.602176634e-19_dp

    type :: dk_surface_harmonic_t
        integer :: m = 0, n = 0
        real(dp) :: c = 0, s = 0
    end type

    type :: dk_surface_input_t
        character(1024) :: eq_path = ''
        integer :: n_theta = 0, n_zeta = 0, n_leg = 0, collisions = 0, nelem = 0, &
            order = 0
        real(dp) :: x_max = 0
        real(dp), allocatable :: z(:), mass(:)
        real(dp), allocatable :: s(:), omega_te(:)
        ! (surface, species)
        real(dp), allocatable :: density(:, :), temperature(:, :), density_s(:, :), &
            temperature_s(:, :), kappa(:, :)
        type(dk_surface_harmonic_t), allocatable :: phi(:)
    end type

    !> Per surface: geometry (q = 1/iota, I/G, kappa = iota B0/(G + iota I),
    !> B0 the mean of |B| on the grid [T], dpsi/ds [T m^2]); d(a, b, i, j) of
    !> surface_transport; <B u_par> per species [T m/s] and the bootstrap
    !> current <j_par B> [T A/m^2]; the local bootstrap current j_bs [A/m^2]
    !> on the grid (point i_theta + n_theta (i_zeta - 1)); with Phi_1 also
    !> d_phi, its <B u_par> and <j_par B> and its local current j_phi.
    type :: dk_surface_result_t
        integer :: status = DKS_BAD_INPUT
        real(dp) :: q = 0, i_over_g = 0, kappa = 0, b0 = 0, psi_prime = 0
        real(dp), allocatable :: d(:, :, :, :), flow(:), j_bs(:)
        real(dp) :: bootstrap = 0
        real(dp), allocatable :: d_phi(:, :, :, :), flow_phi(:), j_phi(:)
        real(dp) :: current_phi = 0
    end type

contains

    subroutine dk_surface_read_input(path, inp, status)
        character(*), intent(in) :: path
        type(dk_surface_input_t), intent(out) :: inp
        integer, intent(out) :: status
        integer :: u, nsp, ns, i, a, nphi, io
        status = DKS_BAD_INPUT
        open (newunit=u, file=path, status='old', action='read', iostat=io)
        if (io /= 0) return
        read (u, '(a)') inp%eq_path
        read (u, *) inp%n_theta, inp%n_zeta, inp%n_leg, inp%collisions, inp%x_max, &
            inp%nelem, inp%order
        read (u, *) nsp
        allocate (inp%z(nsp), inp%mass(nsp))
        do a = 1, nsp
            read (u, *) inp%z(a), inp%mass(a)
        end do
        read (u, *) ns
        allocate (inp%s(ns), inp%omega_te(ns), inp%density(ns, nsp), &
            inp%temperature(ns, nsp), inp%density_s(ns, nsp), inp%temperature_s(ns, nsp), &
            inp%kappa(ns, nsp))
        do i = 1, ns
            read (u, *) inp%s(i), inp%omega_te(i)
            do a = 1, nsp
                read (u, *) inp%density(i, a), inp%temperature(i, a), inp%density_s(i, a), &
                    inp%temperature_s(i, a), inp%kappa(i, a)
            end do
        end do
        read (u, *, iostat=io) nphi
        if (io /= 0) nphi = 0
        allocate (inp%phi(nphi))
        do i = 1, nphi
            read (u, *) inp%phi(i)%m, inp%phi(i)%n, inp%phi(i)%c, inp%phi(i)%s
        end do
        close (u)
        status = DKS_OK
    end subroutine dk_surface_read_input

    !> b = |B|/B0 on n_theta x n_zeta points of theta in [0, 2 pi) and phi_B in
    !> one field period, with the coefficients surface_transport takes.
    subroutine dk_surface_geometry(surf, n_theta, n_zeta, b, b0, q, i_over_g, kappa, status)
        type(boozer_surface_t), intent(in) :: surf
        integer, intent(in) :: n_theta, n_zeta
        real(dp), allocatable, intent(out) :: b(:, :)
        real(dp), intent(out) :: b0, q, i_over_g, kappa
        integer, intent(out) :: status
        type(boozer_point_t) :: pt
        integer :: i, j
        status = DKS_BAD_INPUT
        if (mod(n_theta, 2) == 0 .or. mod(n_zeta, 2) == 0 .or. n_theta < 3) return
        if (surf%iota == 0 .or. surf%g_cov == 0) return
        allocate (b(n_theta, n_zeta))
        do j = 1, n_zeta
            do i = 1, n_theta
                pt = boozer_point(surf, twopi*(i - 1)/n_theta, twopi*(j - 1)/(n_zeta*surf%nper))
                b(i, j) = pt%bmod
            end do
        end do
        b0 = sum(b)/size(b)
        if (.not. (b0 > 0) .or. any(b <= 0)) return
        b = b/b0
        q = 1/surf%iota
        i_over_g = surf%i_cov/surf%g_cov
        kappa = surf%iota*b0/(surf%g_cov + surf%iota*surf%i_cov)
        status = DKS_OK
    end subroutine dk_surface_geometry

    subroutine dk_surface_profile(eq, inp, result, status)
        type(boozer_file_t), intent(in) :: eq
        type(dk_surface_input_t), intent(in) :: inp
        type(dk_surface_result_t), allocatable, intent(out) :: result(:)
        integer, intent(out) :: status
        type(collop_grid_t) :: grid
        integer :: i, st
        status = DKS_BAD_INPUT
        allocate (result(size(inp%s)))
        if (size(inp%s) < 1 .or. size(inp%z) < 1) return
        call collop_make_grid(inp%x_max, inp%nelem, inp%order, grid, st)
        if (st /= COLLOP_OK) return
        do i = 1, size(inp%s)
            call one_surface(eq, inp, i, grid, result(i))
        end do
        status = merge(DKS_OK, DKS_SOLVE, all(result%status == DKS_OK))
    end subroutine dk_surface_profile

    subroutine one_surface(eq, inp, i, grid, r)
        type(boozer_file_t), intent(in) :: eq
        type(dk_surface_input_t), intent(in) :: inp
        integer, intent(in) :: i
        type(collop_grid_t), intent(in) :: grid
        type(dk_surface_result_t), intent(out) :: r
        type(boozer_surface_t) :: surf
        type(collop_species_t), allocatable :: sp(:)
        real(dp), allocatable :: b(:, :), lnl(:, :), phi1(:, :), u(:, :, :, :), &
            u_phi(:, :, :, :), ln_n(:), ln_t(:), force(:, :), particle(:), energy(:), &
            upar(:)
        real(dp) :: phi_prime, rate, s
        integer :: nsp, n, a, bb, st, its
        logical :: with_phi
        nsp = size(inp%z)
        call boozer_surface(eq, inp%s(i), surf, st)
        if (st /= BOOZER_OK) then
            r%status = DKS_OUTSIDE
            return
        end if
        call dk_surface_geometry(surf, inp%n_theta, inp%n_zeta, b, r%b0, r%q, r%i_over_g, &
            r%kappa, st)
        if (st /= DKS_OK) return
        r%psi_prime = surf%psi_prime
        n = size(b)
        allocate (sp(nsp), lnl(nsp, nsp), ln_n(nsp), ln_t(nsp), force(3, nsp), &
            particle(nsp), energy(nsp), upar(n))
        do a = 1, nsp
            sp(a) = collop_species_t(inp%mass(a), inp%z(a), inp%density(i, a), &
                inp%temperature(i, a))
        end do
        do a = 1, nsp
            lnl(a, :) = 3*sqrt(pi)/8*collop_vth(sp(a))*inp%kappa(i, a) &
                /collop_nuhat(sp(a), sp(a), 1.0_dp)
            ln_n(a) = inp%density_s(i, a)/(inp%density(i, a)*r%psi_prime)
            ln_t(a) = inp%temperature_s(i, a)/(inp%temperature(i, a)*r%psi_prime)
        end do
        phi_prime = surf%iota*inp%omega_te(i)
        ! Incompressible E x B rate Phi' G/(G + iota I) of surface_transport.
        rate = phi_prime/(1 + surf%iota*r%i_over_g)
        with_phi = size(inp%phi) > 0
        allocate (r%d(nsp, nsp, 3, 3), u(n, nsp, nsp, 3), r%flow(nsp), r%j_bs(n))
        if (with_phi) then
            phi1 = potential(inp%phi, inp%n_theta, inp%n_zeta, surf%nper)
            allocate (r%d_phi(nsp, nsp, 3, 2), u_phi(n, nsp, nsp, 2), r%flow_phi(nsp), &
                r%j_phi(n))
            call surface_transport(b, surf%nper, r%q, r%i_over_g, r%kappa, sp, lnl, grid, &
                inp%collisions, inp%n_leg, r%d, its, st, electric_rate=rate, phi1=phi1, &
                d_phi=r%d_phi, flow_local=u, flow_local_phi=u_phi)
        else
            call surface_transport(b, surf%nper, r%q, r%i_over_g, r%kappa, sp, lnl, grid, &
                inp%collisions, inp%n_leg, r%d, its, st, electric_rate=rate, flow_local=u)
        end if
        if (st /= SURFACE_OK) then
            r%status = DKS_SOLVE
            return
        end if
        call surface_fluxes(r%d, r%q, r%i_over_g, sp, ln_n, ln_t, phi_prime, 0.0_dp, &
            particle, energy, r%flow, st)
        if (st /= SURFACE_OK) then
            r%status = DKS_SOLVE
            return
        end if
        ! surface_fluxes returns <b u_par>, b = B/B0.
        r%flow = r%b0*r%flow
        r%bootstrap = sum(sp%z*elementary*sp%n*r%flow)
        do a = 1, nsp
            call dk_surface_force_point(ln_n(a), ln_t(a), phi_prime, 0.0_dp, sp(a)%z, &
                elementary, sp(a)%t, force(1, a), force(2, a), force(3, a))
        end do
        s = r%q/(1 + r%i_over_g/r%q)
        r%j_bs = 0
        do a = 1, nsp
            upar = 0
            do bb = 1, nsp
                upar = upar + 2*s*(u(:, a, bb, 1)*force(1, bb) + u(:, a, bb, 2)*force(2, bb))
            end do
            r%j_bs = r%j_bs + sp(a)%z*elementary*sp(a)%n*upar
        end do
        if (with_phi) then
            r%j_phi = 0
            do a = 1, nsp
                upar = 0
                r%flow_phi(a) = 0
                do bb = 1, nsp
                    upar = upar + 2*s*(u_phi(:, a, bb, 1)*force(1, bb) &
                        + u_phi(:, a, bb, 2)*force(2, bb))
                    r%flow_phi(a) = r%flow_phi(a) + r%b0*2*s*(r%d_phi(a, bb, 3, 1)*force(1, bb) &
                        + r%d_phi(a, bb, 3, 2)*force(2, bb))
                end do
                r%j_phi = r%j_phi + sp(a)%z*elementary*sp(a)%n*upar
            end do
            r%current_phi = sum(sp%z*elementary*sp%n*r%flow_phi)
        end if
        if (.not. all(ieee_is_finite(r%flow)) .or. .not. all(ieee_is_finite(r%j_bs))) then
            r%status = DKS_SOLVE
            return
        end if
        r%status = DKS_OK
    end subroutine one_surface

    function potential(h, n_theta, n_zeta, nper) result(phi)
        type(dk_surface_harmonic_t), intent(in) :: h(:)
        integer, intent(in) :: n_theta, n_zeta, nper
        real(dp) :: phi(n_theta, n_zeta), a, th, ph
        integer :: i, j, k
        phi = 0
        do j = 1, n_zeta
            ph = twopi*(j - 1)/(n_zeta*nper)
            do i = 1, n_theta
                th = twopi*(i - 1)/n_theta
                do k = 1, size(h)
                    a = h(k)%m*th - h(k)%n*nper*ph
                    phi(i, j) = phi(i, j) + h(k)%c*cos(a) + h(k)%s*sin(a)
                end do
            end do
        end do
    end function potential
end module milca_dk_surface_device
