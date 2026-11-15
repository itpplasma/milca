program test_plummer
    ! A Plummer sphere of point masses on the shared N-body Hamiltonian with
    ! coupling -G (G2), units G = M = a = 1: the sampled kinetic and potential
    ! energies of 1000 bodies against 3 pi/64 and -3 pi/32 within the
    ! 1/sqrt(N) sampling noise; 300 bodies over five time units keep their
    ! energy, a virial ratio 2K/|W| near 1 and the half-mass radius near
    ! (2^(2/3) - 1)^(-1/2), as an equilibrium does.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use milca_nbody, only: nbody_system_t, nbody_step, nbody_energy
    use milca_self_gravity, only: gravity_system, plummer_sample
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), dt = 1e-3_dp
    type(nbody_system_t) :: system
    type(fortnum_status_t) :: status
    real(dp), allocatable :: state(:), m(:)
    real(dp) :: k, e0, e, w, rh, virial
    logical :: ok
    integer :: step, n, samples

    call setup(1000)
    call energies(k, w)
    print '(a, 2f9.5, a, 2f9.5)', 'K, W:', k, w, '   Plummer:', 3*pi/64, -3*pi/32
    call require(abs(k/(3*pi/64) - 1) < 0.08_dp .and. abs(w/(-3*pi/32) - 1) < 0.08_dp, &
        'Plummer energies')

    call setup(300)
    call energies(k, w)
    e0 = k + w
    virial = 0
    samples = 0
    do step = 1, 5000
        call nbody_step(system, state, dt, status)
        call require(status%code == FORTNUM_OK, 'gravity step')
        if (step > 2500 .and. mod(step, 50) == 0) then
            call energies(k, w)
            virial = virial + 2*k/abs(w)
            samples = samples + 1
        end if
    end do
    virial = virial/samples
    call energies(k, w)
    e = k + w
    rh = half_mass_radius()
    print '(a, es10.2, a, f7.4, a, f7.4)', 'energy drift:', e/e0 - 1, '  mean virial ratio:', &
        virial, '  half-mass radius:', rh
    call require(abs(e/e0 - 1) < 2e-3_dp, 'energy conservation')
    call require(abs(virial - 1) < 0.15_dp, 'virial equilibrium')
    call require(abs(rh/1.3048_dp - 1) < 0.15_dp, 'half-mass radius of the equilibrium')
    print '(a)', 'Plummer sphere stays in equilibrium'

contains

    subroutine setup(bodies)
        integer, intent(in) :: bodies
        n = bodies
        if (allocated(state)) deallocate (state, m)
        allocate (state(6*n), m(n))
        m = 1.0_dp/n
        call plummer_sample(n, 20261115_int64, state, status)
        call require(status%code == FORTNUM_OK, 'Plummer sample')
        call gravity_system(m, 1.0_dp, system, ok)
        call require(ok, 'gravity system')
    end subroutine setup

    subroutine energies(kin, pot)
        real(dp), intent(out) :: kin, pot
        real(dp) :: total
        kin = 0.5_dp*sum(spread(m, 1, 3)*reshape(state(3*n + 1:6*n), [3, n])**2)
        call nbody_energy(system, state, total, status)
        pot = total - kin
    end subroutine energies

    real(dp) function half_mass_radius() result(r_half)
        real(dp) :: r(n), t
        integer :: i, j
        r = norm2(reshape(state(1:3*n), [3, n]), dim=1)
        do i = 2, n
            t = r(i)
            j = i - 1
            do while (j >= 1)
                if (r(j) <= t) exit
                r(j + 1) = r(j)
                j = j - 1
            end do
            r(j + 1) = t
        end do
        r_half = (r(n/2) + r(n/2 + 1))/2
    end function half_mass_radius

    subroutine require(cond, what)
        logical, intent(in) :: cond
        character(*), intent(in) :: what
        if (.not. cond) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_plummer
