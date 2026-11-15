program test_qed_bremsstrahlung
    ! Bremsstrahlung from the Feynman rules (milca_qed_tree): gauge
    ! invariance of the electron-nucleus and electron-electron amplitudes;
    ! the soft-photon theorem, which factorizes both into elastic Mott and
    ! Moller scattering times the eikonal emission factor; the
    ! nonrelativistic Born photon spectrum of Bethe and Heitler,
    ! d sigma/d omega = (16/3) Z^2 alpha r_e^2 ln((p+p')/(p-p'))/(p^2 omega),
    ! by angular quadrature, with a defect shrinking like p^2. Units m = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_tree
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: m = 1, pi = acos(-1.0_dp)
    real(dp) :: defect(2)

    call gauge()
    call soft_nucleus()
    call soft_electron()
    defect(1) = born_spectrum(0.05_dp)
    defect(2) = born_spectrum(0.02_dp)
    print '(a, 2es10.2)', 'nonrelativistic Born defect at p = 0.05, 0.02:', defect
    call require(abs(defect(1)) < 1e-2_dp, 'Bethe-Heitler nonrelativistic Born spectrum')
    call require(abs(defect(2)) < abs(defect(1))*(0.02_dp/0.05_dp)**2*1.5_dp + 1e-7_dp, &
        'Born defect shrinks like p^2')
    print '(a)', 'bremsstrahlung: gauge invariant, soft limits and Born spectrum hold'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    pure function direction(theta, phi) result(n)
        real(dp), intent(in) :: theta, phi
        real(dp) :: n(3)
        n = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    end function direction

    pure function photon(e, n) result(k)
        real(dp), intent(in) :: e, n(3)
        type(four_vector_t) :: k
        k%c = [e, e*n]
    end function photon

    ! e(p1) e(p2) -> e(p3) e(p4) gamma(k) in the centre-of-mass frame: the
    ! pair recoiling against k decays isotropically in its rest frame along n
    subroutine two_electrons_photon(p, omega, nk, n, p1, p2, p3, p4, k)
        real(dp), intent(in) :: p, omega, nk(3), n(3)
        type(four_vector_t), intent(out) :: p1, p2, p3, p4, k
        real(dp) :: w, rest(0:3), mass, q, beta(3), gam, bp, e3(0:3), e4(0:3)
        p1 = on_shell(m, [0.0_dp, 0.0_dp, p])
        p2 = on_shell(m, [0.0_dp, 0.0_dp, -p])
        w = 2*p1%c(0)
        k = photon(omega, nk)
        rest = [w - omega, -omega*nk]
        mass = sqrt(rest(0)**2 - sum(rest(1:3)**2))
        q = sqrt(mass**2/4 - m**2)
        e3 = [mass/2, q*n]
        e4 = [mass/2, -q*n]
        beta = rest(1:3)/rest(0)
        gam = rest(0)/mass
        bp = dot_product(beta, e3(1:3))
        p3%c(0) = gam*(e3(0) + bp)
        p3%c(1:3) = e3(1:3) + ((gam - 1)*bp/dot_product(beta, beta) + gam*e3(0))*beta
        bp = dot_product(beta, e4(1:3))
        p4%c(0) = gam*(e4(0) + bp)
        p4%c(1:3) = e4(1:3) + ((gam - 1)*bp/dot_product(beta, beta) + gam*e4(0))*beta
    end subroutine two_electrons_photon

    subroutine gauge()
        type(four_vector_t) :: p, pp, k, p1, p2, p3, p4
        real(dp) :: e, omega, x
        e = on_shell_energy(2.0_dp)
        omega = 0.4_dp*(e - m)
        p = on_shell(m, [0.0_dp, 0.0_dp, 2.0_dp])
        k = photon(omega, direction(0.7_dp, 0.4_dp))
        pp = on_shell(m, direction(1.9_dp, 2.3_dp)*sqrt((e - omega)**2 - m**2))
        x = nucleus_bremsstrahlung_spin_sum(m, p, pp, k)
        call require(ward_nucleus_bremsstrahlung(m, p, pp, k)/(sqrt(x)*sum((pp%c(1:3) + &
            k%c(1:3) - p%c(1:3))**2)) < 1e-12_dp, 'nucleus bremsstrahlung gauge invariance')
        call two_electrons_photon(1.5_dp, 0.6_dp, direction(0.9_dp, 1.2_dp), &
            direction(2.2_dp, 0.3_dp), p1, p2, p3, p4, k)
        x = electron_bremsstrahlung_spin_sum(m, p1, p2, p3, p4, k)
        call require(maxval(abs(p1%c + p2%c - p3%c - p4%c - k%c)) < 1e-12_dp, &
            'electron-electron kinematics conserve four-momentum')
        call require(ward_electron_bremsstrahlung(m, p1, p2, p3, p4, k) < 1e-11_dp*sqrt(x), &
            'electron-electron bremsstrahlung gauge invariance')
    end subroutine gauge

    pure real(dp) function on_shell_energy(p)
        real(dp), intent(in) :: p
        on_shell_energy = sqrt(m**2 + p**2)
    end function on_shell_energy

    ! sum over polarizations of |eps . J|^2 for a real current J
    pure real(dp) function eikonal(k, j)
        type(four_vector_t), intent(in) :: k
        real(dp), intent(in) :: j(0:3)
        complex(dp) :: eps(0:3, 2)
        integer :: l
        eps = transverse_polarizations(k)
        eikonal = 0
        do l = 1, 2
            eikonal = eikonal + abs(sum(real(eps(1:3, l), dp)*j(1:3)))**2
        end do
    end function eikonal

    subroutine soft_nucleus()
        type(four_vector_t) :: p, pp, pe, k
        real(dp) :: e, omega, pmag, ppmag, x, j(0:3)
        pmag = 0.8_dp
        e = on_shell_energy(pmag)
        omega = 1e-7_dp
        p = on_shell(m, [0.0_dp, 0.0_dp, pmag])
        k = photon(omega, direction(1.1_dp, 0.5_dp))
        ppmag = sqrt((e - omega)**2 - m**2)
        pp = on_shell(m, direction(0.9_dp, 2.0_dp)*ppmag)
        pe = on_shell(m, direction(0.9_dp, 2.0_dp)*pmag)
        x = nucleus_bremsstrahlung_spin_sum(m, p, pp, k)*ppmag/pmag
        j = pp%c/minkowski(pp, k) - p%c/minkowski(p, k)
        call require(abs(x/(mott_spin_sum(m, p, pe)*eikonal(k, j)) - 1) < 1e-5_dp, &
            'soft-photon theorem, electron on a nucleus')
    end subroutine soft_nucleus

    subroutine soft_electron()
        type(four_vector_t) :: p1, p2, p3, p4, k
        real(dp) :: x, j(0:3)
        call two_electrons_photon(0.9_dp, 1e-7_dp, direction(1.3_dp, 0.2_dp), &
            direction(0.8_dp, 1.7_dp), p1, p2, p3, p4, k)
        x = electron_bremsstrahlung_spin_sum(m, p1, p2, p3, p4, k)
        j = p3%c/minkowski(p3, k) + p4%c/minkowski(p4, k) - p1%c/minkowski(p1, k) - &
            p2%c/minkowski(p2, k)
        call require(abs(x/(moller_spin_sum(m, p1, p2, p3, p4)*eikonal(k, j)) - 1) < 1e-5_dp, &
            'soft-photon theorem, electron-electron')
    end subroutine soft_electron

    ! relative defect of the quadrature spectrum against the nonrelativistic
    ! Born formula at momentum p and omega = 0.3 of the kinetic energy; the
    ! electron angle is integrated in ln q^2 to follow the Coulomb peak
    real(dp) function born_spectrum(p) result(defect)
        real(dp), intent(in) :: p
        integer, parameter :: n = 48
        real(dp) :: x(n), w(n), e, omega, pp, total, lo, hi, q2, cosine, ck, phi, jac
        type(four_vector_t) :: pin, pout, k
        integer :: a, b, c
        call gauss_legendre(n, x, w)
        e = on_shell_energy(p)
        omega = 0.3_dp*(e - m)
        pp = sqrt((e - omega)**2 - m**2)
        pin = on_shell(m, [0.0_dp, 0.0_dp, p])
        lo = log((p - pp)**2)
        hi = log((p + pp)**2)
        total = 0
        do a = 1, n
            q2 = exp(lo + (hi - lo)*(x(a) + 1)/2)
            cosine = (p**2 + pp**2 - q2)/(2*p*pp)
            jac = (hi - lo)/2*q2/(2*p*pp)
            do b = 1, n
                ck = x(b)
                k = photon(omega, [sqrt(1 - ck**2), 0.0_dp, ck])
                do c = 1, n
                    phi = pi*(x(c) + 1)
                    pout = on_shell(m, pp*[sqrt(1 - cosine**2)*cos(phi), &
                        sqrt(1 - cosine**2)*sin(phi), cosine])
                    total = total + w(a)*jac*w(b)*w(c)*pi*2*pi* &
                        nucleus_bremsstrahlung_spin_sum(m, pin, pout, k)
                end do
            end do
        end do
        ! d sigma/(d omega dOmega dOmega_k) = Z^2 alpha r_e^2 X p' omega/(4 pi^2 p)
        total = total*pp*omega/(4*pi**2*p)
        defect = total/(16.0_dp/3/(p**2*omega)*log((p + pp)/(p - pp))) - 1
    end function born_spectrum

end program test_qed_bremsstrahlung
