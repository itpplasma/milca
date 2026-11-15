! Electron-nucleus bremsstrahlung rates from the Born amplitudes of
! milca_qed_tree: the photon spectrum d sigma/d omega by angular quadrature,
! the energy an electron radiates per unit path and density of nuclei, and
! the radiated power of a relativistic Maxwellian (Juttner) plasma. Units
! m = c = 1; cross sections in Z^2 alpha r_e^2. The photon angle is
! integrated in ln(1 - beta cos theta_k) and the electron angle in
! ln |p - p' - k|^2, which follow the forward and Coulomb peaks.
! An optional nuclear charge z applies the Coulomb correction of Elwert (Ann.
! Phys. 426, 178, 1939), the ratio of the Coulomb-wave densities at the
! nucleus, (beta/beta') (1 - exp(-2 pi nu))/(1 - exp(-2 pi nu')), nu = z
! alpha/beta. It keeps the spectrum finite at the tip, where the Born result
! vanishes, and matches Sommerfeld's exact nonrelativistic tip up to the
! factor (pi nu/sinh(pi nu))^2.
module milca_qed_bremsstrahlung
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_tree, only: four_vector_t, on_shell, minkowski, &
        nucleus_bremsstrahlung_spin_sum, electron_bremsstrahlung_spin_sum
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    public :: nucleus_photon_spectrum, nucleus_energy_loss, nucleus_thermal_power
    public :: electron_energy_radiated, electron_photon_spectrum

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: alpha = 7.2973525693e-3_dp

contains

    !> d sigma/d omega of an electron of total energy e for photon energy
    !> omega < e - 1, with n quadrature points per angle; Born unless the
    !> nuclear charge z is present.
    real(dp) function nucleus_photon_spectrum(e, omega, n, z) result(total)
        real(dp), intent(in) :: e, omega
        integer, intent(in) :: n
        real(dp), intent(in), optional :: z
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
        if (present(z)) total = total*elwert(z, beta, pp/(e - omega))
    end function nucleus_photon_spectrum

    pure real(dp) function elwert(z, beta, beta_out)
        real(dp), intent(in) :: z, beta, beta_out
        elwert = beta/beta_out*(-expm1(-2*pi*z*alpha/beta))/(-expm1(-2*pi*z*alpha/beta_out))
    end function elwert

    pure real(dp) function expm1(x)
        real(dp), intent(in) :: x
        if (abs(x) < 1e-5_dp) then
            expm1 = x*(1 + x/2*(1 + x/3))
        else
            expm1 = exp(x) - 1
        end if
    end function expm1

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
    real(dp) function nucleus_energy_loss(e, nang, nomega, z) result(total)
        real(dp), intent(in) :: e
        integer, intent(in) :: nang, nomega
        real(dp), intent(in), optional :: z
        real(dp) :: x(nomega), w(nomega), omega, y, jac
        integer :: i
        call gauss_legendre(nomega, x, w)
        total = 0
        do i = 1, nomega
            y = (x(i) + 1)/2
            omega = (e - 1)*sin(pi*y/2)**2
            jac = (e - 1)*pi*sin(pi*y/2)*cos(pi*y/2)/2
            total = total + w(i)*jac*omega*nucleus_photon_spectrum(e, omega, nang, z)
        end do
    end function nucleus_energy_loss

    !> Radiated power per n_e n_i of a Juttner plasma of temperature
    !> theta = kT/(m c^2), in Z^2 alpha r_e^2 m c^3; the kinetic energy is
    !> integrated up to umax theta with nmom points in its square root.
    real(dp) function nucleus_thermal_power(theta, nang, nomega, nmom, umax, z) result(power)
        real(dp), intent(in) :: theta, umax
        integer, intent(in) :: nang, nomega, nmom
        real(dp), intent(in), optional :: z
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
            power = power + weight*(p/gam)*nucleus_energy_loss(gam, nang, nomega, z)
        end do
        power = power/norm
    end function nucleus_thermal_power

    !> Photon energy radiated per electron-electron collision in the
    !> centre-of-mass frame, integral of omega d sigma/d omega over the
    !> photon energy, in alpha r_e^2 m c^2, for invariant s > 4.
    real(dp) function electron_energy_radiated(s, n) result(total)
        real(dp), intent(in) :: s
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), omax, y, omega, jw
        integer :: i
        call gauss_legendre(n, x, w)
        omax = (s - 4)/(2*sqrt(s))
        total = 0
        do i = 1, n
            y = (x(i) + 1)/2
            omega = omax*sin(pi*y/2)**2
            jw = omax*pi*sin(pi*y/2)*cos(pi*y/2)/2
            if (omega <= 0) cycle
            total = total + w(i)*jw*omega*electron_photon_spectrum(s, omega, n)
        end do
    end function electron_energy_radiated

    !> d sigma/d omega of electron-electron bremsstrahlung in the
    !> centre-of-mass frame, in alpha r_e^2, integrated over the photon
    !> direction and the electron directions. The electron direction in the
    !> rest frame of the outgoing pair is integrated over one hemisphere
    !> (exchanging the two electrons covers the other, which cancels the
    !> identical-particle 1/2) in ln |t| about the incoming electron's
    !> direction, which follows the Coulomb peak.
    real(dp) function electron_photon_spectrum(s, omega, n) result(total)
        real(dp), intent(in) :: s, omega
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), root_s, e, p, ck, mass, q, beta(3), gam
        real(dp) :: n1(3), e1(3), e2(3), tr(3), p1s(0:3), tlo, thi, tt, ca, ja, phi, flux
        real(dp) :: b2, bp, e3(0:3)
        type(four_vector_t) :: p1, p2, p3, p4, k
        integer :: j, a, b
        call gauss_legendre(n, x, w)
        root_s = sqrt(s)
        e = root_s/2
        p = sqrt(e**2 - 1)
        p1 = on_shell(1.0_dp, [0.0_dp, 0.0_dp, p])
        p2 = on_shell(1.0_dp, [0.0_dp, 0.0_dp, -p])
        ! F = sqrt((p1.p2)^2 - 1) = sqrt(s (s - 4))/2
        flux = sqrt(s*(s - 4))/2
        total = 0
        do j = 1, n
            ck = x(j)
            k%c = [omega, omega*sqrt(max(0.0_dp, 1 - ck**2)), 0.0_dp, omega*ck]
            ! the outgoing pair: total four-momentum P - k, mass, boost
            mass = sqrt((root_s - omega)**2 - omega**2)
            q = sqrt(max(mass**2/4 - 1, 0.0_dp))
            if (q <= 0) cycle
            beta = -k%c(1:3)/(root_s - omega)
            b2 = dot_product(beta, beta)
            gam = 1/sqrt(1 - b2)
            ! p1 in the pair rest frame
            bp = dot_product(beta, p1%c(1:3))
            p1s(0) = gam*(p1%c(0) - bp)
            p1s(1:3) = p1%c(1:3) + ((gam - 1)*bp/b2 - gam*p1%c(0))*beta
            n1 = p1s(1:3)/norm2(p1s(1:3))
            tr = [1.0_dp, 0.0_dp, 0.0_dp]
            if (abs(n1(1)) > 0.9_dp) tr = [0.0_dp, 1.0_dp, 0.0_dp]
            e1 = tr - dot_product(tr, n1)*n1
            e1 = e1/norm2(e1)
            e2 = [n1(2)*e1(3) - n1(3)*e1(2), n1(3)*e1(1) - n1(1)*e1(3), &
                n1(1)*e1(2) - n1(2)*e1(1)]
            ! |t| = 2 (E1* E3* - p1* q cos a) - 2, linear in cos a on [0, 1]
            tlo = log(max(2*(p1s(0)*mass/2 - norm2(p1s(1:3))*q) - 2, tiny(1.0_dp)))
            thi = log(2*(p1s(0)*mass/2) - 2)
            do a = 1, n
                tt = exp(tlo + (thi - tlo)*(x(a) + 1)/2)
                ca = (p1s(0)*mass/2 - (tt + 2)/2)/(norm2(p1s(1:3))*q)
                ja = (thi - tlo)/2*tt/(2*norm2(p1s(1:3))*q)
                do b = 1, n
                    phi = pi*(x(b) + 1)
                    e3(1:3) = q*(ca*n1 + sqrt(max(0.0_dp, 1 - ca**2))* &
                        (cos(phi)*e1 + sin(phi)*e2))
                    e3(0) = mass/2
                    p3 = boost(e3, beta, gam)
                    p4%c = p1%c + p2%c - k%c - p3%c
                    ! dPhi3 = omega d omega dOmega_k/(16 pi^3) q dOmega*/(16 pi^2 M);
                    ! d sigma/(alpha r_e^2) = 64 pi^3 X/(4F) dPhi3; photon
                    ! azimuth 2 pi, electron azimuth pi dx
                    total = total + w(j)*w(a)*ja*w(b)*pi*2*pi* &
                        64*pi**3/(4*flux)*omega/(16*pi**3)*q/(16*pi**2*mass)* &
                        electron_bremsstrahlung_spin_sum(1.0_dp, p1, p2, p3, p4, k)
                end do
            end do
        end do
    end function electron_photon_spectrum

    ! the four-vector v of the frame moving with -beta seen from the lab:
    ! boost from the pair rest frame back to the centre of mass
    function boost(v, beta, gam) result(u)
        real(dp), intent(in) :: v(0:3), beta(3), gam
        type(four_vector_t) :: u
        real(dp) :: b2, bp
        b2 = dot_product(beta, beta)
        bp = dot_product(beta, v(1:3))
        u%c(0) = gam*(v(0) + bp)
        if (b2 == 0) then
            u%c(1:3) = v(1:3)
        else
            u%c(1:3) = v(1:3) + ((gam - 1)*bp/b2 + gam*v(0))*beta
        end if
    end function boost

end module milca_qed_bremsstrahlung
