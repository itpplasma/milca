program test_fuzzy_soliton
    ! The fuzzy-dark-matter soliton (Schrodinger-Poisson ground state,
    ! hbar = m = G = 1) against derived identities: the virial theorem
    ! 2 K + W = 0, the eigenvalue relation eps M = K + 2 W, and the scaling
    ! symmetry M -> lambda M, r -> r/lambda for psi0 -> lambda^2 psi0; the
    ! invariant r_half M against the 3.925 of the literature.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_fuzzy_dark_matter
    implicit none
    type(soliton_t) :: s, t

    s = soliton_ground_state(1.0_dp)
    call require(s%ok, 'shooting')
    print '(a, 4f12.7)', 'M, eps, K, W:', s%mass, s%eps, s%kinetic, s%potential
    print '(a, 2es11.3)', 'virial and eigenvalue defects:', (2*s%kinetic + s%potential)/s%potential, &
        (s%eps*s%mass - s%kinetic - 2*s%potential)/s%potential
    call require(abs((2*s%kinetic + s%potential)/s%potential) < 1e-6_dp, 'virial theorem')
    call require(abs((s%eps*s%mass - s%kinetic - 2*s%potential)/s%potential) < 1e-6_dp, &
        'eigenvalue relation')
    t = soliton_ground_state(4.0_dp)
    call require(abs(t%mass/(2*s%mass) - 1) < 1e-6_dp .and. abs(t%r_half/(s%r_half/2) - 1) < 1e-6_dp, &
        'scaling symmetry')
    print '(a, f10.5)', 'r_half M:', s%r_half*s%mass
    call require(abs(s%r_half*s%mass/3.925_dp - 1) < 2e-3_dp, 'half-mass radius times mass')
    print '(a)', 'fuzzy-dark-matter soliton satisfies its identities'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_fuzzy_soliton
