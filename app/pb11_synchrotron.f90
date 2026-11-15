! Synchrotron power escaping a uniform slab of thickness L between walls of
! reflectivity R, against its optically thin value (the Juttner-averaged
! Larmor power), for electrons at temperature T in a field B; CSV on standard
! output.
!   pb11_synchrotron T_keV n_e[m^-3] B[T] L[m] R
program pb11_synchrotron
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_synchrotron_transport, only: slab_synchrotron_power
    use milca_qed_thermal, only: scaled_bessel_k
    implicit none
    ! CODATA 2018
    real(dp), parameter :: e = 1.602176634e-19_dp, me = 9.1093837015e-31_dp
    real(dp), parameter :: eps0 = 8.8541878128e-12_dp, c = 299792458.0_dp
    real(dp), parameter :: mc2_kev = 510.99895_dp, pi = acos(-1.0_dp)
    real(dp) :: t, ne, b, l, r, theta, wc, wpe2, lambda, unit, esc, larmor
    character(64) :: arg(5)
    integer :: i

    do i = 1, 5
        call get_command_argument(i, arg(i))
    end do
    read (arg(1), *) t
    read (arg(2), *) ne
    read (arg(3), *) b
    read (arg(4), *) l
    read (arg(5), *) r
    theta = t/mc2_kev
    wc = e*b/me
    wpe2 = ne*e**2/(eps0*me)
    lambda = wpe2*l/(wc*c)
    ! per electron, K omega_c^2 with K = e^2/(8 pi^2 epsilon_0 c)
    unit = e**2/(8*pi**2*eps0*c)*wc**2
    esc = slab_synchrotron_power(theta, lambda, r, 24, 8, 1e-9_dp)
    larmor = 8*pi/3*theta*scaled_bessel_k(3, 1/theta)/scaled_bessel_k(2, 1/theta)
    print '(a)', 'T_keV,n_e,B,L,R,Lambda,P_thin_Wm3,P_escape_Wm3,escape_fraction'
    print '(f0.1, 4(",", es11.4), 4(",", es12.5))', t, ne, b, l, r, lambda, ne*unit*larmor, &
        ne*unit*esc, esc/larmor
end program pb11_synchrotron
