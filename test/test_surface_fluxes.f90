program test_surface_fluxes
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_surface, only: surface_fluxes, SURFACE_OK
    use milca_dk_collop, only: collop_species_t
    use milca_dk_mono_leaves, only: dk_surface_force_point, dk_surface_force_point_jvp, &
        dk_surface_force_point_vjp, dk_surface_flux_point, dk_surface_flux_point_jvp, &
        dk_surface_flux_point_vjp, dk_surface_torque_point, dk_surface_torque_point_jvp, &
        dk_surface_torque_point_vjp
    implicit none
    real(dp), parameter :: elementary = 1.602176634e-19_dp
    character(1024) :: check, path
    call get_command_argument(1, check)
    call get_command_argument(2, path)
    select case (trim(check))
    case ('sfincs'); call sfincs_reference(trim(path))
    case ('equilibrium'); call equilibrium()
    case ('ad'); call derivatives()
    case default; error stop 'unknown surface flux check'
    end select
contains
    subroutine sfincs_reference(path)
        character(*), intent(in) :: path
        type(collop_species_t) :: species(1)
        real(dp), parameter :: amplitude = 1e-6_dp
        real(dp) :: d(1, 1, 3, 3), expected(3, 3), mean_b2, ln_n(1), ln_t(1), field, &
            particle(1), energy(1), flow(1), torque(1), current, error
        integer :: unit, cases, k, model, j, st
        species(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2000*elementary)
        open (newunit=unit, file=path, status='old')
        read (unit, *) cases
        do k = 1, cases
            read (unit, *) model, mean_b2
            read (unit, *) d
            read (unit, *) expected
            do j = 1, 3
                ln_n = 0; ln_t = 0; field = 0
                if (j == 1) ln_n = amplitude
                if (j == 2) then
                    ln_n = 1.5_dp*amplitude
                    ln_t = amplitude
                end if
                if (j == 3) field = 2000*amplitude/mean_b2
                call surface_fluxes(d, 0.8_dp, 0.05_dp, species, ln_n, ln_t, 0.0_dp, &
                    field, particle, energy, flow, st, torque, current)
                if (st /= SURFACE_OK) error stop 'SI response status'
                error = maxval(abs([particle(1), energy(1), flow(1)]/expected(:, j)-1))
                print '(a,i0,a,i0,a,es12.4)', 'SFINCS collision ', model, ', drive ', j, &
                    ': SI moment relative error ', error
                if (error > 0.01_dp) error stop 'SFINCS SI moments'
                if (abs(current-elementary*particle(1)) > 1e-13_dp*abs(current)) &
                    error stop 'radial current units'
                if (abs(torque(1)*0.8_dp-current) > 1e-13_dp*abs(current)) &
                    error stop 'toroidal torque units'
            end do
        end do
        close (unit)
    end subroutine sfincs_reference

    subroutine equilibrium()
        type(collop_species_t) :: species(2)
        real(dp) :: d(2, 2, 3, 3), particle(2), energy(2), flow(2), torque(2), current
        real(dp) :: reference, scale
        integer :: i, st
        species(1) = collop_species_t(6.64e-27_dp, 2.0_dp, 2e19_dp, 2000*elementary)
        species(2) = collop_species_t(9.109e-31_dp, -1.0_dp, 4e19_dp, 2000*elementary)
        d = reshape([(real(i, dp)/10, i=1, size(d))], shape(d))
        call surface_fluxes(d, 0.8_dp, 0.05_dp, species, [-2.0_dp, 1.0_dp], &
            [0.0_dp, 0.0_dp], 2000.0_dp, 0.0_dp, particle, energy, flow, st)
        if (st /= SURFACE_OK) error stop 'Boltzmann response status'
        if (maxval(abs(particle)/species%n) > 1e-12_dp .or. &
            maxval(abs(energy)/(species%n*species%t)) > 1e-12_dp .or. &
            maxval(abs(flow)) > 1e-12_dp) error stop 'Boltzmann equilibrium'
        d = 0
        d(1, 2, 1, 1) = 3
        call surface_fluxes(d, 0.8_dp, 0.05_dp, species, [0.0_dp, 1.0_dp], &
            [0.0_dp, 0.0_dp], 0.0_dp, 0.0_dp, particle, energy, flow, st, torque, current)
        if (st /= SURFACE_OK) error stop 'cross-species response status'
        scale = 0.8_dp/(1+0.05_dp/0.8_dp)
        reference = -2*species(1)%n*scale**2*3
        if (abs(particle(1)/reference-1) > 1e-13_dp .or. particle(2) /= 0) &
            error stop 'receiving-species density'
        if (abs(sum(torque)*0.8_dp/current-1) > 1e-13_dp) &
            error stop 'species-current torque ledger'
    end subroutine equilibrium

    subroutine derivatives()
        real(dp) :: p(16), direction(16), seed(3), tangent(3), gradient(16), fd(3), &
            error, dual, lhs, rhs
        real(dp), parameter :: step = 1e-5_dp
        integer :: i, kind
        p = [(real(i+1, dp)/10, i=1, 16)]
        direction = [(real((-1)**i, dp)*i/30, i=1, 16)]
        seed = [0.7_dp, -0.2_dp, 0.3_dp]
        do kind = 1, 3
            call tangent_point(kind, p, direction, tangent)
            call reverse_point(kind, p, seed, gradient)
            fd = (point(kind, p+step*direction)-point(kind, p-step*direction))/(2*step)
            error = maxval(abs(fd-tangent))/max(1e-12_dp, maxval(abs(fd)))
            lhs = dot_product(seed, tangent)
            rhs = dot_product(gradient, direction)
            dual = abs(lhs-rhs)/max(1e-12_dp, abs(lhs)+abs(rhs))
            print '(a,i0,a,2es12.4)', 'surface kernel ', kind, ': tangent and dual errors ', error, dual
            if (error > 1e-8_dp .or. dual > 1e-12_dp) error stop 'generated surface derivatives'
        end do
    end subroutine derivatives

    function point(kind, p) result(f)
        integer, intent(in) :: kind
        real(dp), intent(in) :: p(16)
        real(dp) :: f(3)
        f = 0
        select case (kind)
        case (1)
            call dk_surface_force_point(p(1), p(2), p(3), p(4), p(5), p(6), p(7), f(1), f(2), f(3))
        case (2)
            call dk_surface_flux_point(p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), &
                p(9), p(10), p(11), p(12), p(13), p(14), p(15), p(16), f(1), f(2), f(3))
        case (3)
            call dk_surface_torque_point(p(1), p(2), p(3), p(4), f(1), f(2))
        end select
    end function point

    subroutine tangent_point(kind, p, v, f)
        integer, intent(in) :: kind
        real(dp), intent(in) :: p(16), v(16)
        real(dp), intent(out) :: f(3)
        f = 0
        select case (kind)
        case (1)
            call dk_surface_force_point_jvp(p(1), p(2), p(3), p(4), p(5), p(6), p(7), &
                v(1), v(2), v(3), v(4), v(5), v(6), v(7), f(1), f(2), f(3))
        case (2)
            call dk_surface_flux_point_jvp(p(1), p(2), p(3), p(4), p(5), p(6), p(7), &
                p(8), p(9), p(10), p(11), p(12), p(13), p(14), p(15), p(16), &
                v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11), &
                v(12), v(13), v(14), v(15), v(16), f(1), f(2), f(3))
        case (3)
            call dk_surface_torque_point_jvp(p(1), p(2), p(3), p(4), v(1), v(2), v(3), &
                v(4), f(1), f(2))
        end select
    end subroutine tangent_point

    subroutine reverse_point(kind, p, seed, gradient)
        integer, intent(in) :: kind
        real(dp), intent(in) :: p(16), seed(3)
        real(dp), intent(out) :: gradient(16)
        gradient = 0
        select case (kind)
        case (1)
            call dk_surface_force_point_vjp(p(1), p(2), p(3), p(4), p(5), p(6), p(7), &
                seed(1), seed(2), seed(3), gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5), gradient(6), gradient(7))
        case (2)
            call dk_surface_flux_point_vjp(p(1), p(2), p(3), p(4), p(5), p(6), p(7), &
                p(8), p(9), p(10), p(11), p(12), p(13), p(14), p(15), p(16), &
                seed(1), seed(2), seed(3), gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5), gradient(6), gradient(7), gradient(8), gradient(9), gradient(10), &
                gradient(11), gradient(12), gradient(13), gradient(14), gradient(15), gradient(16))
        case (3)
            call dk_surface_torque_point_vjp(p(1), p(2), p(3), p(4), seed(1), seed(2), &
                gradient(1), gradient(2), gradient(3), gradient(4))
        end select
    end subroutine reverse_point
end program test_surface_fluxes
