!> G2: Newtonian self-gravitating point masses as the shared N-body
!> Hamiltonian of milca_nbody with the Coulomb coupling replaced by -G and
!> the charges by the masses (one pair-force definition for both). No
!> softening, as for point charges. The Plummer sphere (Plummer 1911) is the
!> equilibrium: density (3 M/(4 pi a^3)) (1 + r^2/a^2)^(-5/2), sampled from its
!> distribution function (Aarseth, Henon and Wielen 1974) in units G = M = a = 1,
!> where the kinetic energy is 3 pi/64 and the potential energy -3 pi/32.
module milca_self_gravity
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, status_set
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use milca_nbody, only: nbody_system_t, nbody_isolated
    use milca_species, only: units_t, species_t
    implicit none
    private
    public :: gravity_system, plummer_sample

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    !> An isolated self-gravitating system of the given masses.
    subroutine gravity_system(masses, g, system, ok)
        real(dp), intent(in) :: masses(:), g
        type(nbody_system_t), intent(out) :: system
        logical, intent(out) :: ok
        type(species_t) :: bodies(size(masses))
        type(units_t) :: units
        units%coulomb = -g
        bodies%name = 'mass'
        bodies%charge = masses
        bodies%mass = masses
        call nbody_isolated(system, units, bodies, ok)
    end subroutine gravity_system

    !> n equal masses 1/n of a Plummer sphere (G = M = a = 1) in the centre-of-
    !> mass frame; state(1:3n) positions, state(3n+1:6n) velocities. Radii
    !> beyond 100 a (a mass fraction 1.5e-4) are redrawn.
    subroutine plummer_sample(n, seed, state, status)
        integer, intent(in) :: n
        integer(int64), intent(in) :: seed
        real(dp), intent(out) :: state(:)
        type(fortnum_status_t), intent(out) :: status
        type(rng_t) :: g
        real(dp) :: u, r, q, y, v, d(3), mean(3)
        integer :: i
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'state size must be 6 n')
        if (size(state) /= 6*n .or. n < 1) return
        call rng_seed(g, seed, status)
        if (status%code /= FORTNUM_OK) return
        do i = 1, n
            ! enclosed mass u = r^3/(1 + r^2)^(3/2) inverted
            do
                call rng_uniform(g, u)
                if (u <= 0) cycle
                r = 1/sqrt(u**(-2.0_dp/3) - 1)
                if (r <= 100) exit
            end do
            state(3*i - 2:3*i) = r*direction()
            ! speed fraction q of the escape speed, density q^2 (1 - q^2)^(7/2)
            do
                call rng_uniform(g, q)
                call rng_uniform(g, y)
                if (0.1_dp*y < q**2*(1 - q**2)**3.5_dp) exit
            end do
            v = q*sqrt(2.0_dp)*(1 + r**2)**(-0.25_dp)
            state(3*n + 3*i - 2:3*n + 3*i) = v*direction()
        end do
        do i = 0, 1
            mean = sum(reshape(state(3*n*i + 1:3*n*(i + 1)), [3, n]), dim=2)/n
            state(3*n*i + 1:3*n*(i + 1)) = state(3*n*i + 1:3*n*(i + 1)) - &
                reshape(spread(mean, 2, n), [3*n])
        end do
    contains
        function direction() result(e)
            real(dp) :: e(3), z, phi
            call rng_uniform(g, d(1))
            call rng_uniform(g, d(2))
            z = 2*d(1) - 1
            phi = 2*pi*d(2)
            e = [sqrt(1 - z**2)*cos(phi), sqrt(1 - z**2)*sin(phi), z]
        end function direction
    end subroutine plummer_sample

end module milca_self_gravity
