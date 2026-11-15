! Electron-nucleus bremsstrahlung rates from the Born amplitudes of
! milca_qed_tree: the photon spectrum d sigma/d omega by angular quadrature,
! the energy an electron radiates per unit path and density of nuclei, and
! the radiated power of a relativistic Maxwellian (Juttner) plasma. Units
! m = c = 1; cross sections in Z^2 alpha r_e^2. The photon angle is
! integrated in ln(1 - beta cos theta_k) and the electron angle in
! ln |p - p' - k|^2, which follow the forward and Coulomb peaks.
module milca_qed_bremsstrahlung
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_tree, only: four_vector_t, on_shell, nucleus_bremsstrahlung_spin_sum
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: nucleus_photon_spectrum, nucleus_energy_loss, nucleus_thermal_power

    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    !> d sigma/d omega of an electron of total energy e for photon energy
    !> omega < e - 1, with n quadrature points per angle.
    real(dp) function nucleus_photon_spectrum(e, omega, n) result(total)
        real(dp), intent(in) :: e, omega
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), p, pp, beta, lo, hi, q2, jac, slo, shi, sk, ck, jk, phi, a(3)
        type(four_vector_t) :: pin, pout, k
        integer :: i, j, l
        call gauss_legendre(n, x, w)
        p = sqrt(e**2 - 1)
        pp = sqrt((e - omega)**2 - 1)
        beta = p/e
        pin = on_shell(1.0_dp, [0.0_dp, 0.0_dp, p])
        slo = log(1 - beta)
        shi = log(1 + beta)
        total = 0
        do j = 1, n
            sk = slo + (shi - slo)*(x(j) + 1)/2
            ck = (1 - exp(sk))/beta
            jk = (shi - slo)/2*exp(sk)/beta
            k%c = [omega, omega*sqrt(max(0.0_dp, 1 - ck**2)), 0.0_dp, omega*ck]
            a = pin%c(1:3) - k%c(1:3)
            lo = log(max((norm2(a) - pp)**2, tiny(1.0_dp)))
            hi = log((norm2(a) + pp)**2)
            do i = 1, n
                q2 = exp(lo + (hi - lo)*(x(i) + 1)/2)
                jac = (hi - lo)/2*q2/(2*norm2(a)*pp)
                do l = 1, n
                    phi = pi*(x(l) + 1)
                    pout = on_sphere(a, pp, q2, phi)
                    total = total + w(i)*w(j)*w(l)*jac*jk*pi*2*pi* &
                        nucleus_bremsstrahlung_spin_sum(1.0_dp, pin, pout, k)
                end do
            end do
        end do
        ! d sigma/(d omega dOmega' dOmega_k) = X p' omega/(4 pi^2 p)
        total = total*pp*omega/(4*pi**2*p)
    end function nucleus_photon_spectrum

    ! p' of modulus pp at azimuth phi around a, with |a - p'|^2 = q2
    function on_sphere(a, pp, q2, phi) result(v)
        real(dp), intent(in) :: a(3), pp, q2, phi
        type(four_vector_t) :: v
        real(dp) :: n(3), e1(3), e2(3), c, s, t(3)
        n = a/norm2(a)
        t = [1.0_dp, 0.0_dp, 0.0_dp]
        if (abs(n(1)) > 0.9_dp) t = [0.0_dp, 1.0_dp, 0.0_dp]
        e1 = t - dot_product(t, n)*n
        e1 = e1/norm2(e1)
        e2 = [n(2)*e1(3) - n(3)*e1(2), n(3)*e1(1) - n(1)*e1(3), n(1)*e1(2) - n(2)*e1(1)]
        c = max(-1.0_dp, min(1.0_dp, (dot_product(a, a) + pp**2 - q2)/(2*norm2(a)*pp)))
        s = sqrt(1 - c**2)
        v = on_shell(1.0_dp, pp*(c*n + s*(cos(phi)*e1 + sin(phi)*e2)))
    end function on_sphere

    !> Radiated energy per unit path and unit density of nuclei,
    !> integral of omega d sigma/d omega over 0 < omega < e - 1. With
    !> omega = (e - 1) sin^2(pi y/2) the logarithm at omega -> 0 and the
    !> square root of p' at the upper end become smooth in y.
    real(dp) function nucleus_energy_loss(e, nang, nomega) result(total)
        real(dp), intent(in) :: e
        integer, intent(in) :: nang, nomega
        real(dp) :: x(nomega), w(nomega), omega, y, jac
        integer :: i
        call gauss_legendre(nomega, x, w)
        total = 0
        do i = 1, nomega
            y = (x(i) + 1)/2
            omega = (e - 1)*sin(pi*y/2)**2
            jac = (e - 1)*pi*sin(pi*y/2)*cos(pi*y/2)/2
            total = total + w(i)*jac*omega*nucleus_photon_spectrum(e, omega, nang)
        end do
    end function nucleus_energy_loss

    !> Radiated power per n_e n_i of a Juttner plasma of temperature
    !> theta = kT/(m c^2), in Z^2 alpha r_e^2 m c^3; the kinetic energy is
    !> integrated up to umax theta with nmom points in its square root.
    real(dp) function nucleus_thermal_power(theta, nang, nomega, nmom, umax) result(power)
        real(dp), intent(in) :: theta, umax
        integer, intent(in) :: nang, nomega, nmom
        real(dp) :: x(nmom), w(nmom), u, t, gam, p, weight, norm
        integer :: i
        call gauss_legendre(nmom, x, w)
        power = 0
        norm = 0
        do i = 1, nmom
            ! u = t^2 keeps p ~ t smooth at the origin
            t = sqrt(umax)*(x(i) + 1)/2
            u = t**2
            gam = 1 + theta*u
            p = sqrt(gam**2 - 1)
            ! p^2 dp exp(-(gam-1)/theta) = p gam theta exp(-u) du, du = 2 t dt
            weight = w(i)*sqrt(umax)/2*2*t*p*gam*theta*exp(-u)
            norm = norm + weight
            power = power + weight*(p/gam)*nucleus_energy_loss(gam, nang, nomega)
        end do
        power = power/norm
    end function nucleus_thermal_power

end module milca_qed_bremsstrahlung
