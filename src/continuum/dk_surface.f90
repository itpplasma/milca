!> Transport matrix of a 3D surface with the energy-coupled linearised
!! collision operator in real arithmetic (the n = 0 problem has no complex
!! part), the faster solver of the operator of ntv_surface_coefficients: same
!! drives, responses and units (see there), same collision rows
!! (ntv_collision_rows), same Legendre x Fourier surface discretisation, with
!! the density null space of each speed node bordered as in dk_mono (a
!! density source and <w h_0> = 0 per node instead of the unweighted mean).
!!
!! Unknowns h(i, l, m): surface point i, Legendre degree l, speed-species
!! index m = (a - 1) n_x + k. After the row scaling the operator is
!! S x R + E x W - I x C_l (S streaming and mirror force per unit speed, R the
!! streaming coefficient of each node, C_l the collision rows of degree l);
!! S acts on all m at once as two dense surface matrices per degree.
!! FGMRES is right-preconditioned by one dk_mono factorisation per eigenmode
!! of R^-1 C_4 (complex pairs by their real and imaginary parts), with the
!! diagonal of V^-1 R^-1 C_l V as the rate of degree l: exact for pitch-angle
!! scattering. Memory is that of n_x n_s monoenergetic
!! surface problems; the modes are factorised and solved in parallel.
module milca_dk_surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_legendre, only: legendre_p, legendre_p_derivative
    use milca_dk_mono, only: mono_factor_t, mono_surface_point, mono_factor_rates, &
        mono_solve, mono_coupling, mono_surface_operators, MONO_OK
    use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_moments
    use milca_dk_mono_leaves, only: dk_surface_electric_scale, dk_surface_force_point, &
        dk_surface_flux_point, dk_surface_torque_point, dk_surface_moment_weights, &
        dk_surface_density_point, dk_surface_quasineutral_point, dk_surface_potential_point
    implicit none
    private
    public :: surface_transport, surface_fluxes, surface_density, surface_quasineutral_phi1
    public :: surface_full_exb_t, surface_full_exb_setup, surface_full_exb_apply
    integer, parameter, public :: SURFACE_OK = 0, SURFACE_BAD_INPUT = 1, &
        SURFACE_SINGULAR = 2, SURFACE_NO_CONVERGENCE = 3

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    ! Relative residual of the FGMRES solves.
    real(dp), parameter :: tolerance = 1e-10_dp
    integer, parameter :: restart = 100, max_iterations = 5000

    type :: surface_full_exb_t
        private
        real(dp), allocatable :: d_perp(:, :), measure(:), metric(:), &
            mass(:, :), speed_derivative(:, :), speed_pitch(:, :), pitch_derivative(:, :)
    end type

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, &
                lwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine dgeev
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    !> Weak full electrostatic drift on a fixed surface with f = exp(-x^2) h.
    !! The velocity coordinate x is in thermal-speed units and the pitch
    !! coordinate is v_parallel/v. Integrating the conservative phase-volume
    !! form by parts makes the density moment vanish in the discrete action.
    subroutine surface_full_exb_setup(b, nper, q, i_over_g, grid, n_leg, operator, status)
        real(dp), intent(in) :: b(:, :), q, i_over_g
        integer, intent(in) :: nper, n_leg
        type(collop_grid_t), intent(in) :: grid
        type(surface_full_exb_t), intent(out) :: operator
        integer, intent(out) :: status
        real(dp), allocatable :: parallel(:, :), b_parallel(:), b_perp(:), flat_b(:), &
            unused_mass(:, :), unused_moment(:), pitch(:), weight(:), p(:, :), dp_dxi(:, :)
        integer :: l, count

        status = SURFACE_BAD_INPUT
        if (grid%n < 3 .or. n_leg < 3) return
        if (.not. all(ieee_is_finite(b)) .or. .not. all(ieee_is_finite([q, i_over_g]))) return
        call mono_surface_operators(b, nper, q, i_over_g, parallel, b_parallel, b_perp, &
            status, operator%d_perp)
        if (status /= MONO_OK) then
            status = SURFACE_BAD_INPUT
            return
        end if
        flat_b = reshape(b, [size(b)])
        operator%measure = 1/flat_b**2
        operator%metric = operator%measure*b_perp/(2*flat_b)
        allocate (operator%mass(grid%n, grid%n), operator%speed_derivative(grid%n, grid%n), &
            unused_mass(grid%n, grid%n), unused_moment(grid%n))
        call collop_moments(grid, 2, 2, operator%mass, unused_moment)
        call collop_moments(grid, 3, 2, unused_mass, unused_moment, operator%speed_derivative)
        count = n_leg + 1
        allocate (pitch(count), weight(count), p(count, 0:n_leg-1), dp_dxi(count, 0:n_leg-1))
        call gauss_legendre(count, pitch, weight)
        do l = 0, n_leg-1
            p(:, l) = legendre_p(l, 0, pitch)
            dp_dxi(:, l) = legendre_p_derivative(l, 0, pitch)
        end do
        operator%speed_pitch = matmul(transpose(p), spread(weight*(1+pitch**2), 2, n_leg)*p)
        operator%pitch_derivative = matmul(transpose(dp_dxi), &
            spread(weight*pitch*(1-pitch**2), 2, n_leg)*p)
        status = SURFACE_OK
    end subroutine surface_full_exb_setup

    !> Full electric contribution at unit angular rate to the Galerkin
    !! residual. h(point, degree, speed) and the output include x=0. The
    !! output carries the unnormalised phase-volume weight on the test basis.
    subroutine surface_full_exb_apply(operator, h, action, status)
        type(surface_full_exb_t), intent(in) :: operator
        real(dp), intent(in) :: h(:, 0:, :)
        real(dp), intent(out) :: action(:, 0:, :)
        integer, intent(out) :: status
        integer :: lmax, n, nx

        status = SURFACE_BAD_INPUT
        if (.not. allocated(operator%mass)) return
        n = size(operator%measure)
        nx = size(operator%mass, 1)
        lmax = size(operator%speed_pitch, 1)-1
        if (any(shape(h) /= [n, lmax+1, nx]) .or. any(shape(h) /= shape(action))) return
        call full_exb_action(operator, h, action)
        status = SURFACE_OK
    end subroutine surface_full_exb_apply

    subroutine full_exb_action(operator, h, action)
        type(surface_full_exb_t), intent(in) :: operator
        real(dp), intent(in) :: h(:, 0:, :)
        real(dp), intent(out) :: action(:, 0:, :)
        real(dp), allocatable :: mass(:, :, :), derivative(:, :, :)
        integer :: l, j, lmax, n, nx
        n = size(h, 1)
        lmax = ubound(h, 2)
        nx = size(h, 3)
        allocate (mass(n, 0:lmax, nx), derivative(n, 0:lmax, nx))
        do l = 0, lmax
            mass(:, l, :) = matmul(h(:, l, :), transpose(operator%mass))
            derivative(:, l, :) = matmul(h(:, l, :), transpose(operator%speed_derivative))
        end do
        do l = 0, lmax
            action(:, l, :) = -2.0_dp/(2*l+1)*matmul(transpose(operator%d_perp), &
                spread(operator%measure, 2, nx)*mass(:, l, :))
            do j = l-2, l+2, 2
                if (j < 0 .or. j > lmax) cycle
                action(:, l, :) = action(:, l, :) - spread(operator%metric, 2, nx) &
                    *(operator%speed_pitch(l+1, j+1)*derivative(:, j, :) &
                    +operator%pitch_derivative(l+1, j+1)*mass(:, j, :))
            end do
        end do
    end subroutine full_exb_action

    !> SI moments of a surface_transport matrix. Logarithmic gradients and
    !> phi_prime are derivatives in toroidal flux psi_t (Wb/radian), and
    !> E_parallel = parallel_field*b. Outputs are <Gamma.grad psi_t>, the
    !> kinetic-energy flux, and <b*u_parallel>. The optional torque is the
    !> toroidal Lorentz torque on each species for
    !> B = grad psi_t x grad theta + iota grad zeta x grad psi_t.
    subroutine surface_fluxes(d, q, i_over_g, species, ln_n_prime, ln_t_prime, &
            phi_prime, parallel_field, particle_flux, energy_flux, parallel_flow, &
            status, torque, radial_current)
        real(dp), intent(in) :: d(:, :, :, :), q, i_over_g, ln_n_prime(:), &
            ln_t_prime(:), phi_prime, parallel_field
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(out) :: particle_flux(:), energy_flux(:), parallel_flow(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: torque(:), radial_current
        real(dp), parameter :: elementary = 1.602176634e-19_dp
        real(dp) :: forces(3, size(species)), contribution(3), current, moment
        integer :: ns, a, sb
        status = SURFACE_BAD_INPUT
        ns = size(species)
        particle_flux = 0
        energy_flux = 0
        parallel_flow = 0
        if (present(torque)) torque = 0
        if (present(radial_current)) radial_current = 0
        if (ns < 1 .or. any(shape(d) /= [ns, ns, 3, 3])) return
        if (size(ln_n_prime) /= ns .or. size(ln_t_prime) /= ns) return
        if (size(particle_flux) /= ns .or. size(energy_flux) /= ns .or. &
            size(parallel_flow) /= ns) return
        if (present(torque)) then
            if (size(torque) /= ns) return
        end if
        if (.not. all(ieee_is_finite([q, i_over_g, phi_prime, parallel_field]))) return
        if (q == 0 .or. q+i_over_g == 0) return
        if (.not. all(ieee_is_finite(d))) return
        if (.not. all(ieee_is_finite(ln_n_prime)) .or. &
            .not. all(ieee_is_finite(ln_t_prime))) return
        do a = 1, ns
            if (.not. all(ieee_is_finite([species(a)%n, species(a)%t, species(a)%z]))) return
            if (species(a)%n <= 0 .or. species(a)%t <= 0 .or. species(a)%z == 0) return
            call dk_surface_force_point(ln_n_prime(a), ln_t_prime(a), phi_prime, &
                parallel_field, species(a)%z, elementary, species(a)%t, &
                forces(1, a), forces(2, a), forces(3, a))
        end do
        do a = 1, ns
            do sb = 1, ns
                call dk_surface_flux_point(q, i_over_g, species(a)%n, species(a)%t, &
                    d(a, sb, 1, 1), d(a, sb, 2, 1), d(a, sb, 3, 1), &
                    d(a, sb, 1, 2), d(a, sb, 2, 2), d(a, sb, 3, 2), &
                    d(a, sb, 1, 3), d(a, sb, 2, 3), d(a, sb, 3, 3), &
                    forces(1, sb), forces(2, sb), forces(3, sb), &
                    contribution(1), contribution(2), contribution(3))
                particle_flux(a) = particle_flux(a)+contribution(1)
                energy_flux(a) = energy_flux(a)+contribution(2)
                parallel_flow(a) = parallel_flow(a)+contribution(3)
            end do
            call dk_surface_torque_point(particle_flux(a), q, species(a)%z, elementary, &
                current, moment)
            if (.not. all(ieee_is_finite([current, moment]))) return
            if (present(torque)) torque(a) = moment
            if (present(radial_current)) radial_current = radial_current+current
        end do
        if (.not. all(ieee_is_finite(particle_flux)) .or. &
            .not. all(ieee_is_finite(energy_flux)) .or. &
            .not. all(ieee_is_finite(parallel_flow))) return
        if (present(radial_current)) then
            if (.not. ieee_is_finite(radial_current)) return
        end if
        status = SURFACE_OK
    end subroutine surface_fluxes

    !> Local density perturbation delta_n(i, a) [m^-3] of each species at the
    !> points i of the grid, from density_local of surface_transport and the
    !> gradients and fields of surface_fluxes.
    subroutine surface_density(density_local, q, i_over_g, species, ln_n_prime, &
            ln_t_prime, phi_prime, parallel_field, delta_n, status)
        real(dp), intent(in) :: density_local(:, :, :, :), q, i_over_g, ln_n_prime(:), &
            ln_t_prime(:), phi_prime, parallel_field
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(out) :: delta_n(:, :)
        integer, intent(out) :: status
        real(dp), parameter :: elementary = 1.602176634e-19_dp
        real(dp) :: forces(3, size(species))
        real(dp), allocatable :: contribution(:)
        integer :: ns, n, a, sb
        status = SURFACE_BAD_INPUT
        ns = size(species)
        n = size(density_local, 1)
        delta_n = 0
        if (ns < 1 .or. any(shape(density_local) /= [n, ns, ns, 3])) return
        if (any(shape(delta_n) /= [n, ns])) return
        if (size(ln_n_prime) /= ns .or. size(ln_t_prime) /= ns) return
        if (.not. all(ieee_is_finite([q, i_over_g, phi_prime, parallel_field]))) return
        if (q == 0 .or. q+i_over_g == 0) return
        if (.not. all(ieee_is_finite(density_local))) return
        do a = 1, ns
            if (species(a)%n <= 0 .or. species(a)%t <= 0 .or. species(a)%z == 0) return
            call dk_surface_force_point(ln_n_prime(a), ln_t_prime(a), phi_prime, &
                parallel_field, species(a)%z, elementary, species(a)%t, &
                forces(1, a), forces(2, a), forces(3, a))
        end do
        allocate (contribution(n))
        do a = 1, ns
            do sb = 1, ns
                call dk_surface_density_point(q, i_over_g, species(a)%n, &
                    density_local(:, a, sb, 1), density_local(:, a, sb, 2), &
                    density_local(:, a, sb, 3), forces(1, sb), forces(2, sb), forces(3, sb), &
                    contribution)
                delta_n(:, a) = delta_n(:, a)+contribution
            end do
        end do
        if (.not. all(ieee_is_finite(delta_n))) return
        status = SURFACE_OK
    end subroutine surface_density

    !> Linear quasineutral potential phi1 [V] on the grid of b: the charge
    !> sum_a Z_a e delta_n_a of the kinetic species equals Phi_1 sum Z^2 e^2 n/T
    !> over the kinetic species (Boltzmann factor of F) and the adiabatic
    !> ones, in the gauge <Phi_1> = 0 with the surface measure 1/b^2.
    !> delta_n(i, a) [m^-3] as from surface_density, point i = i_theta +
    !> n_theta (i_zeta - 1).
    subroutine surface_quasineutral_phi1(b, species, delta_n, phi1, status, adiabatic)
        real(dp), intent(in) :: b(:, :), delta_n(:, :)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(out) :: phi1(:, :)
        integer, intent(out) :: status
        type(collop_species_t), intent(in), optional :: adiabatic(:)
        real(dp), parameter :: elementary = 1.602176634e-19_dp
        real(dp), allocatable :: charge(:), contribution(:), parts(:), weight(:), potential(:)
        real(dp) :: shielding, part, unused
        integer :: n, a
        status = SURFACE_BAD_INPUT
        n = size(b)
        phi1 = 0
        if (size(species) < 1 .or. any(shape(delta_n) /= [n, size(species)])) return
        if (any(shape(phi1) /= shape(b)) .or. any(b <= 0)) return
        if (.not. all(ieee_is_finite(delta_n))) return
        allocate (charge(n), contribution(n), parts(n), potential(n))
        charge = 0
        shielding = 0
        do a = 1, size(species)
            if (species(a)%n <= 0 .or. species(a)%t <= 0 .or. species(a)%z == 0) return
            call dk_surface_quasineutral_point(species(a)%z, species(a)%n, species(a)%t, &
                elementary, delta_n(:, a), contribution, parts)
            charge = charge+contribution
            shielding = shielding+parts(1)
        end do
        if (present(adiabatic)) then
            do a = 1, size(adiabatic)
                if (adiabatic(a)%n <= 0 .or. adiabatic(a)%t <= 0) return
                if (adiabatic(a)%z == 0) return
                call dk_surface_quasineutral_point(adiabatic(a)%z, adiabatic(a)%n, &
                    adiabatic(a)%t, elementary, 0.0_dp, unused, part)
                shielding = shielding+part
            end do
        end if
        weight = reshape(1/b**2, [n])
        charge = charge-sum(weight*charge)/sum(weight)
        call dk_surface_potential_point(charge, shielding, potential)
        phi1 = reshape(potential, shape(b))
        if (.not. all(ieee_is_finite(phi1))) return
        status = SURFACE_OK
    end subroutine surface_quasineutral_phi1

    !> Arguments and d as in ntv_surface_coefficients; iterations the FGMRES
    !> iterations summed over the 3 n_s drives and requested tangents.
    !> electric_rate = Phi'_t G/(G+iota I) in 1/s selects incompressible
    !> electric drift. full_local_electric includes electric speed and pitch
    !> advection with the same per-speed distribution constraint. d_drate
    !> holds field, species and collision inputs fixed.
    !>
    !> phi1(n_theta, n_zeta) is a prescribed in-surface potential in volts on
    !> the grid of b. With F = f_M(psi, H), H = m v^2/2 + Z e Phi_1, the
    !> unknown is f - F, and the response linear in Phi_1 has three sources,
    !> eps = Z e Phi_1/T:
    !>   (a) the radial E x B drift on f_M, whose ratio to the magnetic drift is
    !>       Z e B d_perp Phi_1 / (T x^2 (1 + xi^2) d_perp B), d_perp = d_theta -
    !>       (I/G) d_zeta: the radial-drive row with (1 + xi^2) x^2 tz b_drift/b
    !>       replaced by d_perp Phi_1 (P_0 only);
    !>   (b) the magnetic drift on the Phi_1 part of d F/d psi at fixed H,
    !>       eps f_M (-f1 + (1 - x^2) f2);
    !>   (c) the parallel electric field accelerating the first-order response g,
    !>       x' = -(v_t/2) xi E, xi' = -(v_t/(2 x)) (1 - xi^2) E,
    !>       E = b.grad eps, in the conservative weak form.
    !> The inductive drive is proportional to F, so its response has (c) and
    !>   (d) the Boltzmann factor of F, -eps times the inductive source.
    !> d_phi(a, sb, i, jf) holds the moments of d for the response of species
    !> a to these sources with force jf = 1..3 of species sb (per unit force).
    !> flow_local(i, a, sb, jf) is the local parallel flow U at point i of the
    !> response to drive jf = 1..3, normalised so that <b U> = d(a, sb, 3, jf)
    !> with <.> of weight 1/b^2; flow_local_phi(i, a, sb, jf) the same for the
    !> phi1 drive, <b U> = d_phi(a, sb, 3, jf). Point i = i_theta +
    !> n_theta (i_zeta - 1). condense_isotropic (default true) condenses the
    !> isotropic x = 0 node of the collision rows, which makes the collisions
    !> conserve density exactly; false sets h_0(0) = 0 as the NTV energy
    !> solvers do. density_local(i, a, sb, jf) is the local density moment
    !> N of the response to drive jf, delta n/n = 2 (s (N_1 f_1 + N_2 f_2) +
    !> N_3 f_3) (surface_density); its average with weight 1/b^2 is zero.
    subroutine surface_transport(b, nper, q, i_over_g, kappa, species, lnlambda, grid, &
            model, n_leg, d, iterations, status, electric_rate, d_drate, full_local_electric, &
            phi1, d_phi, flow_local, flow_local_phi, condense_isotropic, density_local)
        real(dp), intent(in) :: b(:, :), q, i_over_g, kappa
        integer, intent(in) :: nper, n_leg
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model
        real(dp), intent(out) :: d(:, :, :, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(in), optional :: electric_rate
        real(dp), intent(out), optional :: d_drate(:, :, :, :)
        logical, intent(in), optional :: full_local_electric
        real(dp), intent(in), optional :: phi1(:, :)
        real(dp), intent(out), optional :: d_phi(:, :, :, :), flow_local(:, :, :, :), &
            flow_local_phi(:, :, :, :)
        logical, intent(in), optional :: condense_isotropic
        real(dp), intent(out), optional :: density_local(:, :, :, :)
        type(mono_factor_t) :: geo
        type(mono_factor_t), allocatable :: modes(:)
        type(surface_full_exb_t) :: full_exb
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), d3(:), d4(:), d6(:), &
            vt(:), tz(:), coll(:, :, :), lo(:, :, :), up(:, :, :), v(:, :), vu(:, :), &
            rinv(:), rhs(:, :, :), x(:, :, :), radial(:), measure(:), bv(:), &
            electric_unit(:), electric_rows(:), dx(:, :, :), d_par(:, :), b_par(:), &
            b_drift(:), d_perp(:, :), dphi(:), xg(:, :, :, :, :), eps(:), epar(:), &
            pitch_a(:, :), pitch_b(:, :), d2(:, :), m1(:, :), unused(:, :), unused_v(:)
        real(dp) :: j0, rate
        integer :: n, lmax, ns, nx, nm, a, sb, jf, k, l, st, its
        logical :: full, condense

        status = SURFACE_BAD_INPUT
        condense = .true.
        if (present(condense_isotropic)) condense = condense_isotropic
        iterations = 0
        ns = size(species)
        if (any(shape(d) /= [ns, ns, 3, 3]) .or. mod(n_leg, 2) /= 0) return
        if (present(d_phi) .or. present(flow_local_phi)) then
            if (.not. present(phi1)) return
        end if
        if (present(phi1)) then
            if (any(shape(phi1) /= shape(b))) return
            if (.not. all(ieee_is_finite(phi1))) return
        end if
        if (present(d_phi)) then
            if (any(shape(d_phi) /= [ns, ns, 3, 3])) return
            d_phi = 0
        end if
        if (present(flow_local)) then
            if (any(shape(flow_local) /= [size(b), ns, ns, 3])) return
            flow_local = 0
        end if
        if (present(flow_local_phi)) then
            if (any(shape(flow_local_phi) /= [size(b), ns, ns, 3])) return
            flow_local_phi = 0
        end if
        if (present(density_local)) then
            if (any(shape(density_local) /= [size(b), ns, ns, 3])) return
            density_local = 0
        end if
        rate = 0
        full = .false.
        if (present(full_local_electric)) full = full_local_electric
        if (present(electric_rate)) rate = electric_rate
        if (.not. ieee_is_finite(rate)) return
        if (present(d_drate)) then
            if (any(shape(d_drate) /= shape(d))) return
            d_drate = 0
        end if
        if (present(electric_rate) .or. present(d_drate)) then
            call mono_surface_point(b, nper, q, i_over_g, kappa, n_leg, geo, st, 0.0_dp)
        else
            call mono_surface_point(b, nper, q, i_over_g, kappa, n_leg, geo, st)
        end if
        if (st /= MONO_OK) return
        if (full .and. (rate /= 0 .or. present(d_drate))) then
            call surface_full_exb_setup(b, nper, q, i_over_g, grid, n_leg, full_exb, st)
            if (st /= SURFACE_OK) return
        end if
        call ntv_collision_rows(species, lnlambda, grid, model, n_leg, rows, ratio, d3, &
            d4, d6, vt, tz, st, condense_isotropic=condense)
        if (st /= NTV_OK) return
        n = geo%n
        lmax = geo%lmax
        nx = size(ratio, 1)
        nm = nx*ns
        allocate (electric_unit(nm), electric_rows(nm))
        do a = 1, ns
            call dk_surface_electric_scale(1.0_dp, vt(a), grid%x(2:), ratio(:, a), &
                                           electric_unit((a - 1)*nx + 1:a*nx))
        end do
        electric_rows = rate*electric_unit
        ! Collision rows over (node, species), node fastest.
        allocate (coll(nm, nm, 0:lmax), lo(n, n, 0:lmax), up(n, n, 0:lmax))
        do l = 0, lmax
            do sb = 1, ns
                do a = 1, ns
                    do k = 1, nx
                        coll((a - 1)*nx + k, (sb - 1)*nx + 1:sb*nx, l) = rows(k, :, a, sb, l)
                    end do
                end do
            end do
            call mono_coupling(geo, l, lo(:, :, l), up(:, :, l))
        end do
        rinv = 1/reshape(ratio, [nm])
        status = SURFACE_SINGULAR
        call modal_factors(st)
        if (st /= SURFACE_OK) return

        bv = reshape(b, [n])
        ! Radial drive (1 + xi^2) S/b, S = b_drift/q, in terms of the dk_mono
        ! drift: P_0 and P_2 components -2 drift0/(kappa q), -2 drift2/(kappa q).
        radial = -2/(kappa*q)*geo%p%drift0
        measure = 2/bv**3
        j0 = sum(1/bv**2)
        allocate (rhs(n, 0:lmax, nm))
        if (present(phi1) .and. (present(d_phi) .or. present(flow_local_phi))) &
            allocate (xg(n, 0:lmax, nm, ns, 3))
        do sb = 1, ns
            do jf = 1, 3
                rhs = 0
                do k = 1, nx
                    associate (m => (sb - 1)*nx + k, row => ratio(k, sb)/(vt(sb)*d3(k)))
                        if (jf < 3) then
                            rhs(:, 0, m) = row*tz(sb)*merge(d4(k), d6(k), jf == 1)*radial
                            rhs(:, 2, m) = row*tz(sb)*merge(d4(k), d6(k), jf == 1) &
                                *(-2/(kappa*q))*geo%p%drift2
                        else
                            rhs(:, 1, m) = row*vt(sb)*d3(k)*geo%p%drive_e
                        end if
                    end associate
                end do
                call solve(rhs, x, its, st)
                iterations = iterations + its
                if (st /= SURFACE_OK) then
                    status = st
                    return
                end if
                do a = 1, ns
                    d(a, sb, :, jf) = moments(x, a)
                    if (present(flow_local)) flow_local(:, a, sb, jf) = local_flow(x, a)
                    if (present(density_local)) density_local(:, a, sb, jf) = &
                        local_density(x, a)
                end do
                if (allocated(xg)) xg(:, :, :, sb, jf) = x
                if (.not. present(d_drate)) cycle
                if (full) then
                    call full_electric_action(x, rhs)
                    rhs = -rhs
                else
                    do k = 1, nm
                        rhs(:, :, k) = -electric_unit(k)*matmul(geo%exb, x(:, :, k))
                    end do
                end if
                call solve(rhs, dx, its, st)
                iterations = iterations + its
                if (st /= SURFACE_OK) then
                    status = st
                    return
                end if
                do a = 1, ns
                    d_drate(a, sb, :, jf) = moments(dx, a)
                end do
            end do
        end do
        if (present(phi1) .and. (present(d_phi) .or. present(flow_local_phi))) then
            call mono_surface_operators(b, nper, q, i_over_g, d_par, b_par, b_drift, st, d_perp)
            if (st /= MONO_OK) then
                status = SURFACE_BAD_INPUT
                return
            end if
            dphi = matmul(d_perp, reshape(phi1, [n]))
            ! b.grad Phi_1 = kappa b (d_theta + q d_zeta) Phi_1.
            epar = kappa*bv*matmul(d_par, reshape(phi1, [n]))
            call parallel_field_operators(st)
            if (st /= SURFACE_OK) then
                status = st
                return
            end if
            do sb = 1, ns
                eps = reshape(phi1, [n])/tz(sb)
                do jf = 1, 3
                    call parallel_field_action(xg(:, :, :, sb, jf), rhs)
                    rhs = -rhs
                    do k = 1, nx
                        associate (m => (sb - 1)*nx + k, row => ratio(k, sb)/(vt(sb)*d3(k)), &
                                w => merge(-d4(k), d4(k) - d6(k), jf == 1))
                            if (jf < 3) then
                                rhs(:, 0, m) = rhs(:, 0, m) + row*merge(d4(k), d6(k), jf == 1) &
                                    /grid%x(k + 1)**2*dphi/q + row*tz(sb)*w*radial*eps
                                rhs(:, 2, m) = rhs(:, 2, m) + row*tz(sb)*w &
                                    *(-2/(kappa*q))*geo%p%drift2*eps
                            else
                                rhs(:, 1, m) = rhs(:, 1, m) - row*vt(sb)*d3(k)*geo%p%drive_e*eps
                            end if
                        end associate
                    end do
                    call solve(rhs, x, its, st)
                    iterations = iterations + its
                    if (st /= SURFACE_OK) then
                        status = st
                        return
                    end if
                    do a = 1, ns
                        if (present(d_phi)) d_phi(a, sb, :, jf) = moments(x, a)
                        if (present(flow_local_phi)) flow_local_phi(:, a, sb, jf) = &
                            local_flow(x, a)
                    end do
                end do
            end do
        end if
        status = SURFACE_OK
    contains
        !> Pitch and speed matrices of the parallel electric acceleration:
        !> pitch_a(l + 1, j + 1) = int xi P_l P_j, pitch_b(l + 1, j + 1) =
        !> int (1 - xi^2) P_l' P_j,
        !> d2(i, j) = int x^2 exp(-x^2) phi_i' phi_j, m1(i, j) = int x exp(-x^2)
        !> phi_i phi_j on the speed grid with the x = 0 node.
        subroutine parallel_field_operators(status)
            integer, intent(out) :: status
            real(dp), allocatable :: xi(:), wt(:), p(:, :), dp_dxi(:, :)
            integer :: count, l
            status = SURFACE_BAD_INPUT
            count = n_leg + 1
            allocate (xi(count), wt(count), p(count, 0:lmax), dp_dxi(count, 0:lmax), &
                d2(grid%n, grid%n), m1(grid%n, grid%n), unused(grid%n, grid%n), &
                unused_v(grid%n))
            call gauss_legendre(count, xi, wt)
            do l = 0, lmax
                p(:, l) = legendre_p(l, 0, xi)
                dp_dxi(:, l) = legendre_p_derivative(l, 0, xi)
            end do
            pitch_a = matmul(transpose(p), spread(wt*xi, 2, lmax + 1)*p)
            pitch_b = matmul(transpose(dp_dxi), spread(wt*(1 - xi**2), 2, lmax + 1)*p)
            call collop_moments(grid, 2, 2, unused, unused_v, d2)
            call collop_moments(grid, 1, 1, m1, unused_v)
            if (.not. all(ieee_is_finite(epar))) return
            status = SURFACE_OK
        end subroutine parallel_field_operators

        !> Rows of (c) on h: the weak residual (v_t/2) E [d2 x pitch_a + m1 x
        !> pitch_b] h, E = b.grad Phi_1/tz, mapped to the solver rows like the
        !> full electric action; v_t and the surface measure cancel.
        subroutine parallel_field_action(h, action)
            real(dp), intent(in) :: h(:, 0:, :)
            real(dp), intent(out) :: action(:, 0:, :)
            real(dp), allocatable :: samples(:, :, :), weak(:, :, :)
            integer :: a, k, l, j, m
            allocate (samples(n, 0:lmax, grid%n), weak(n, 0:lmax, grid%n))
            do a = 1, ns
                samples(:, :, 1) = 0
                samples(:, :, 2:) = h(:, :, (a - 1)*nx + 1:a*nx)
                weak = 0
                do l = 0, lmax
                    do j = max(0, l - 1), min(lmax, l + 1)
                        weak(:, l, :) = weak(:, l, :) &
                            + pitch_a(l + 1, j + 1)*matmul(samples(:, j, :), transpose(d2)) &
                            + pitch_b(l + 1, j + 1)*matmul(samples(:, j, :), transpose(m1))
                    end do
                end do
                ! The x = 0 test function joins node 2 for l = 0 (density).
                weak(:, 0, 2) = weak(:, 0, 2) + weak(:, 0, 1)
                do k = 1, nx
                    m = (a - 1)*nx + k
                    do l = 0, lmax
                        action(:, l, m) = epar/(2*tz(a))*weak(:, l, k + 1)*ratio(k, a) &
                            *(2*l + 1)/(2*d3(k))
                    end do
                end do
            end do
        end subroutine parallel_field_action

        !> U(i) with <b U> = moments(h, a)(3).
        function local_flow(h, species_index) result(value)
            real(dp), intent(in) :: h(:, 0:, :)
            integer, intent(in) :: species_index
            real(dp) :: value(n), density_weight, flow_weight
            integer :: k, m
            value = 0
            do k = 1, nx
                m = (species_index - 1)*nx + k
                call dk_surface_moment_weights(d3(k), grid%x(k + 1), vt(species_index), &
                    density_weight, flow_weight)
                value = value + flow_weight*h(:, 1, m)
            end do
        end function local_flow

        !> Local density of h, normalised like local_flow: delta n/n =
        !> 2 (s (N_1 f_1 + N_2 f_2) + N_3 f_3), see dk_surface_density_point.
        function local_density(h, species_index) result(value)
            real(dp), intent(in) :: h(:, 0:, :)
            integer, intent(in) :: species_index
            real(dp) :: value(n), density_weight, flow_weight
            integer :: k, m
            value = 0
            do k = 1, nx
                m = (species_index - 1)*nx + k
                call dk_surface_moment_weights(d3(k), grid%x(k + 1), vt(species_index), &
                    density_weight, flow_weight)
                value = value + density_weight*h(:, 0, m)
            end do
        end function local_density

        function moments(h, species_index) result(value)
            real(dp), intent(in) :: h(:, 0:, :)
            integer, intent(in) :: species_index
            real(dp) :: value(3)
            integer :: k, m
            value = 0
            do k = 1, nx
                m = (species_index - 1)*nx + k
                value(1) = value(1) + d4(k)*flux(h(:, :, m))
                value(2) = value(2) + d6(k)*flux(h(:, :, m))
                value(3) = value(3) + d3(k)*sum(measure*bv**2*2/3.0_dp*h(:, 1, m))
            end do
            value(1:2) = tz(species_index)/(sqrt(pi)*2*j0)*value(1:2)
            value(3) = vt(species_index)/(sqrt(pi)*2*j0)*value(3)
        end function moments

        !> sum measure S (8/3 h_0 + 4/15 h_2), S = b_drift/q = -3 b drift0/(2 kappa q).
        real(dp) function flux(h)
            real(dp), intent(in) :: h(:, 0:)
            flux = sum(measure*(-1.5_dp*bv*geo%p%drift0/(kappa*q)) &
                *(8/3.0_dp*h(:, 0) + 4/15.0_dp*h(:, 2)))
        end function flux

        !> V, V^-1 and one dk_mono factorisation per eigenmode of R^-1 C_4.
        subroutine modal_factors(status)
            integer, intent(out) :: status
            real(dp), allocatable :: w(:, :), wr(:), wi(:), work(:), vl(:, :), collision_rate(:, :), &
                mode_electric(:)
            integer, allocatable :: piv(:), mode_status(:)
            integer :: m, info, a, k, first, last
            status = SURFACE_SINGULAR
            allocate (w(nm, nm), wr(nm), wi(nm), work(8*nm), vl(1, 1), v(nm, nm), &
                vu(nm, nm), piv(nm), collision_rate(0:lmax, nm), mode_status(nm))
            ! Eigenvectors of degree 4: fewest FGMRES iterations among degrees 1
            ! to 9 on the DIII-D helical core (2402 for degree 1, 857 for 4).
            do m = 1, nm
                w(m, :) = rinv(m)*coll(m, :, min(4, lmax))
            end do
            call dgeev('N', 'V', nm, w, nm, wr, wi, vl, 1, v, nm, work, size(work), info)
            if (info /= 0) return
            ! A complex pair enters V as its real and imaginary parts, a real
            ! basis in which the mode block is 2 x 2 and its diagonal the real part.
            w = v
            vu = 0
            do m = 1, nm
                vu(m, m) = 1
            end do
            call dgesv(nm, nm, w, nm, piv, vu, nm, info)
            if (info /= 0) return
            do l = 0, lmax
                do m = 1, nm
                    w(m, :) = rinv(m)*coll(m, :, l)
                end do
                w = matmul(vu, matmul(w, v))
                do m = 1, nm
                    collision_rate(l, m) = -w(m, m)
                end do
            end do
            if (full .and. rate /= 0) then
                w = 0
                do a = 1, ns
                    first = (a-1)*nx+1
                    last = a*nx
                    do k = 1, nx
                        w(first+k-1, first:last) = rate*full_exb%mass(k+1, 2:)/(vt(a)*d3(k))
                    end do
                end do
                w = matmul(vu, matmul(w, v))
            else
                w = matmul(vu, spread(rinv*electric_rows, 2, nm)*v)
            end if
            mode_electric = [(w(m, m), m=1, nm)]
            allocate (modes(nm))
            !$omp parallel do schedule(dynamic, 1) default(shared) private(m)
            do m = 1, nm
                modes(m) = geo
                if (full .and. rate /= 0) modes(m)%exb = full_exb%d_perp
                modes(m)%omega_e = mode_electric(m)
                call mono_factor_rates(modes(m), collision_rate(:, m), mode_status(m))
            end do
            !$omp end parallel do
            if (any(mode_status /= MONO_OK)) return
            status = SURFACE_OK
        end subroutine modal_factors

        !> FGMRES on the bordered operator: unknowns h(n, 0:lmax, nm) and the
        !> density sources mu(nm).
        subroutine solve(rhs, x, its, status)
            real(dp), intent(in) :: rhs(:, 0:, :)
            real(dp), allocatable, intent(out) :: x(:, :, :)
            integer, intent(out) :: its, status
            real(dp), allocatable :: bvec(:), u(:)
            real(dp) :: scale, residual
            integer :: ntot, info
            ntot = n*(lmax + 1)*nm + nm
            allocate (bvec(ntot), u(ntot))
            bvec = 0
            bvec(1:ntot - nm) = reshape(rhs, [ntot - nm])
            scale = norm2(bvec)
            its = 0
            allocate (x(n, 0:lmax, nm))
            x = 0
            status = SURFACE_OK
            if (scale == 0) return
            bvec = bvec/scale
            u = 0
            call real_fgmres_operator(apply, precondition, bvec, u, tolerance, &
                max_iterations, restart, info, its, residual)
            status = SURFACE_NO_CONVERGENCE
            if (info /= KRYLOV_OK) return
            x = scale*reshape(u(1:ntot - nm), [n, lmax + 1, nm])
            status = SURFACE_OK
        end subroutine solve

        !> y = A x with the border: h-rows plus mu at l = 0, then <w h_0> rows.
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: h(:, :, :), y(:, :, :), full_electric_rows(:, :, :)
            integer :: l, m, nh
            nh = n*(lmax + 1)*nm
            h = reshape(input(1:nh), [n, lmax + 1, nm])
            allocate (y(n, 0:lmax, nm))
            !$omp parallel do schedule(static) default(shared) private(l, m)
            do l = 0, lmax
                block
                    real(dp) :: electric(n, nm)
                    y(:, l, :) = 0
                    if (l > 0) call dgemm('N', 'N', n, nm, n, 1.0_dp, lo(:, :, l), n, &
                        h(:, l, :), n, 1.0_dp, y(:, l, :), n)
                    if (l < lmax) call dgemm('N', 'N', n, nm, n, 1.0_dp, up(:, :, l), n, &
                        h(:, l + 2, :), n, 1.0_dp, y(:, l, :), n)
                    do m = 1, nm
                        y(:, l, m) = ratio_flat(m)*y(:, l, m)
                    end do
                    if (allocated(geo%exb) .and. .not. full) then
                        call dgemm('N', 'N', n, nm, n, 1.0_dp, geo%exb, n, &
                                   h(:, l + 1, :), n, 0.0_dp, electric, n)
                        do m = 1, nm
                            y(:, l, m) = y(:, l, m) + electric_rows(m)*electric(:, m)
                        end do
                    end if
                    call dgemm('N', 'T', n, nm, nm, -1.0_dp, h(:, l + 1, :), n, coll(:, :, l), &
                        nm, 1.0_dp, y(:, l, :), n)
                end block
            end do
            !$omp end parallel do
            if (full .and. rate /= 0) then
                allocate (full_electric_rows(n, 0:lmax, nm))
                call full_electric_action(h, full_electric_rows)
                y = y + rate*full_electric_rows
            end if
            do m = 1, nm
                y(:, 0, m) = y(:, 0, m) + input(nh + m)
                output(nh + m) = sum(geo%p%weight*h(:, 1, m))
            end do
            output(1:nh) = reshape(y, [nh])
        end subroutine apply

        !> Map the electric weak form to the transport solver's scaled rows.
        !! Its basis fixes the coefficient at x=0 to zero.
        subroutine full_electric_action(h, action)
            real(dp), intent(in) :: h(:, 0:, :)
            real(dp), intent(out) :: action(:, 0:, :)
            real(dp), allocatable :: samples(:, :, :), weak(:, :, :)
            integer :: a, k, l, m
            allocate (samples(n, 0:lmax, nx+1), weak(n, 0:lmax, nx+1))
            samples(:, :, 1) = 0
            do a = 1, ns
                samples(:, :, 2:) = h(:, :, (a-1)*nx+1:a*nx)
                call full_exb_action(full_exb, samples, weak)
                ! The x = 0 test function joins node 2 for l = 0, so that the
                ! kept test functions sum to one and the density is conserved.
                weak(:, 0, 2) = weak(:, 0, 2) + weak(:, 0, 1)
                do k = 1, nx
                    m = (a-1)*nx+k
                    do l = 0, lmax
                        action(:, l, m) = weak(:, l, k+1)*ratio(k, a)*(2*l+1) &
                            /(2*full_exb%measure*vt(a)*d3(k))
                    end do
                end do
            end do
        end subroutine full_electric_action

        pure real(dp) function ratio_flat(m)
            integer, intent(in) :: m
            ratio_flat = 1/rinv(m)
        end function ratio_flat

        !> x = P^-1 input: modal coordinates V^-1 R^-1 for the rows and V^-1
        !> for the constraints, one bordered dk_mono solve per mode, back by V
        !> (h) and R V (mu).
        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: r(:, :, :), z(:, :, :), c(:), cm(:), mu(:)
            integer :: l, m, nh
            nh = n*(lmax + 1)*nm
            r = reshape(input(1:nh), [n, lmax + 1, nm])
            allocate (z(n, lmax + 1, nm), mu(nm))
            do l = 1, lmax + 1
                do m = 1, nm
                    r(:, l, m) = rinv(m)*r(:, l, m)
                end do
                z(:, l, :) = matmul(r(:, l, :), transpose(vu))
            end do
            c = input(nh + 1:nh + nm)
            cm = matmul(vu, c)
            !$omp parallel do schedule(dynamic, 1) default(shared) private(m)
            do m = 1, nm
                block
                    real(dp), allocatable :: xm(:, :, :)
                    real(dp) :: mum(1)
                    integer :: st
                    call mono_solve(modes(m), reshape(z(:, :, m), [n, lmax + 1, 1]), &
                        [cm(m)], xm, mum, st)
                    z(:, :, m) = xm(:, :, 1)
                    mu(m) = mum(1)
                end block
            end do
            !$omp end parallel do
            do l = 1, lmax + 1
                r(:, l, :) = matmul(z(:, l, :), transpose(v))
            end do
            output(1:nh) = reshape(r, [nh])
            output(nh + 1:nh + nm) = matmul(v, mu)/rinv
        end subroutine precondition
    end subroutine surface_transport
end module milca_dk_surface
