program test_qed_tree
    ! Tree-level QED two ways: the Feynman-rule amplitudes of milca_qed_tree
    ! against the generated closed forms of milca_qed_leaves at scattered
    ! kinematics; spinor completeness, gauge invariance and crossing; total
    ! cross sections by angular quadrature of the amplitudes against the
    ! Klein-Nishina, Breit-Wheeler and Dirac formulas; Thomson and Rutherford
    ! limits. Units m = 1, cross sections in r_e^2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_qed_tree
    use milca_qed_leaves, only: qed_moller_msq, qed_bhabha_msq, qed_compton_msq, &
        qed_annihilation_msq, qed_mott_dcs, qed_klein_nishina_dcs, &
        qed_klein_nishina_total, qed_klein_nishina_total_low, qed_breit_wheeler_total, &
        qed_dirac_annihilation_total
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: m = 1, pi = acos(-1.0_dp), tol = 1e-11_dp
    integer, parameter :: nq = 256
    real(dp) :: xq(nq), wq(nq)
    real(dp) :: momenta(4) = [0.05_dp, 0.7_dp, 3.0_dp, 40.0_dp]
    real(dp) :: angles(2, 3) = reshape([0.4_dp, 0.3_dp, 1.3_dp, 2.1_dp, 2.6_dp, 5.0_dp], [2, 3])
    integer :: i, j

    call gauss_legendre(nq, xq, wq)
    do i = 1, size(momenta)
        call require(spinor_completeness(m, on_shell(m, [0.3_dp, -0.2_dp, 1.0_dp]*momenta(i))) < &
            1e-12_dp*(1 + momenta(i)), 'spinor completeness')
        do j = 1, size(angles, 2)
            call two_to_two(momenta(i), angles(1, j), angles(2, j))
        end do
    end do
    call totals()
    call limits()
    print '(a)', 'tree-level QED: amplitudes match closed forms and total cross sections'

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

    pure real(dp) function sq(a, b)
        type(four_vector_t), intent(in) :: a, b
        type(four_vector_t) :: d
        d%c = a%c - b%c
        sq = minkowski(d, d)
    end function sq

    pure logical function close(a, b)
        real(dp), intent(in) :: a, b
        close = abs(a - b) <= tol*max(abs(a), abs(b))
    end function close

    ! centre-of-mass kinematics of momentum p, scattering direction (theta, phi)
    subroutine two_to_two(p, theta, phi)
        real(dp), intent(in) :: p, theta, phi
        type(four_vector_t) :: p1, p2, p3, p4, k1, k2, kk
        real(dp) :: n(3), z(3), e, s, t, u, msq, ref, ew
        z = [0.0_dp, 0.0_dp, 1.0_dp]
        n = direction(theta, phi)
        p1 = on_shell(m, p*z)
        p2 = on_shell(m, -p*z)
        p3 = on_shell(m, p*n)
        p4 = on_shell(m, -p*n)
        e = p1%c(0)
        s = 4*e**2
        t = sq(p1, p3)
        u = sq(p1, p4)
        call qed_moller_msq(s, t, u, m, ref)
        call require(close(moller_spin_sum(m, p1, p2, p3, p4), ref), 'Moller closed form')
        call qed_bhabha_msq(s, t, u, m, ref)
        call require(close(bhabha_spin_sum(m, p1, p2, p3, p4), ref), 'Bhabha closed form')
        ! Compton: photon p along z, electron against it
        k1 = photon(p, z)
        p1 = on_shell(m, -p*z)
        k2 = photon(p, n)
        p3 = on_shell(m, -p*n)
        s = minkowski(sum4(p1, k1), sum4(p1, k1))
        u = sq(p1, k2)
        call qed_compton_msq(s, u, m, ref)
        msq = compton_spin_sum(m, p1, k1, p3, k2)
        call require(close(msq, ref), 'Compton closed form')
        call require(ward_compton(m, p1, k1, p3, k2) <= 1e-12_dp*sqrt(msq)*(1 + p), &
            'Compton gauge invariance')
        ! annihilation e- e+ -> gamma gamma and its crossing gamma gamma -> e- e+
        p1 = on_shell(m, p*z)
        p2 = on_shell(m, -p*z)
        ew = p1%c(0)
        k1 = photon(ew, n)
        k2 = photon(ew, -n)
        t = sq(p1, k1)
        u = sq(p1, k2)
        call qed_annihilation_msq(t, u, m, ref)
        msq = annihilation_spin_sum(m, p1, p2, k1, k2)
        call require(close(msq, ref), 'annihilation closed form')
        kk = photon(ew, z)
        k1 = kk
        k2 = photon(ew, -z)
        p1 = on_shell(m, p*n)
        p2 = on_shell(m, -p*n)
        call qed_annihilation_msq(sq(p1, k1), sq(p1, k2), m, ref)
        call require(close(breit_wheeler_spin_sum(m, k1, k2, p1, p2), ref), &
            'Breit-Wheeler by crossing')
    end subroutine two_to_two

    pure function sum4(a, b) result(c)
        type(four_vector_t), intent(in) :: a, b
        type(four_vector_t) :: c
        c%c = a%c + b%c
    end function sum4

    subroutine totals()
        real(dp), parameter :: xs(5) = [1e-3_dp, 0.01_dp, 0.3_dp, 2.0_dp, 50.0_dp]
        real(dp), parameter :: betas(3) = [0.1_dp, 0.6_dp, 0.95_dp]
        real(dp) :: x, cosine, omega_p, dcs, ref, total, beta, e, p, s, gamma_lab
        type(four_vector_t) :: pe, k, kp, pp
        integer :: i, q
        ! Klein-Nishina in the electron rest frame
        do i = 1, size(xs)
            x = xs(i)
            total = 0
            do q = 1, nq
                cosine = xq(q)
                omega_p = x/(1 + x*(1 - cosine))
                pe = on_shell(m, [0.0_dp, 0.0_dp, 0.0_dp])
                k = photon(x, [0.0_dp, 0.0_dp, 1.0_dp])
                kp = photon(omega_p, [sqrt(1 - cosine**2), 0.0_dp, cosine])
                pp%c = pe%c + k%c - kp%c
                dcs = (omega_p/x)**2*compton_spin_sum(m, pe, k, pp, kp)/4
                call qed_klein_nishina_dcs(x, cosine, ref)
                call require(close(dcs, ref), 'Klein-Nishina differential')
                total = total + 2*pi*wq(q)*dcs
            end do
            if (x < 0.2_dp) then
                call qed_klein_nishina_total_low(x, ref)
            else
                call qed_klein_nishina_total(x, ref)
            end if
            call require(abs(total - ref) <= 1e-13_dp*ref, 'Klein-Nishina total by quadrature')
        end do
        ! Breit-Wheeler and Dirac annihilation in the centre-of-mass frame
        do i = 1, 3
            beta = betas(i)
            e = 1/sqrt(1 - beta**2)
            p = beta*e
            s = 4*e**2
            total = 0
            do q = 1, nq
                cosine = xq(q)
                total = total + 2*pi*wq(q)*breit_wheeler_spin_sum(m, photon(e, [0.0_dp, 0.0_dp, 1.0_dp]), &
                    photon(e, [0.0_dp, 0.0_dp, -1.0_dp]), &
                    on_shell(m, p*[sqrt(1 - cosine**2), 0.0_dp, cosine]), &
                    on_shell(m, -p*[sqrt(1 - cosine**2), 0.0_dp, cosine]))/(4*s)*beta
            end do
            call qed_breit_wheeler_total(beta, ref)
            call require(abs(total - ref) <= 1e-10_dp*ref, 'Breit-Wheeler total by quadrature')
            total = 0
            do q = 1, nq
                cosine = xq(q)
                ! identical photons: half the solid angle
                total = total + pi*wq(q)*annihilation_spin_sum(m, on_shell(m, [0.0_dp, 0.0_dp, p]), &
                    on_shell(m, [0.0_dp, 0.0_dp, -p]), photon(e, [sqrt(1 - cosine**2), 0.0_dp, cosine]), &
                    photon(e, [-sqrt(1 - cosine**2), 0.0_dp, -cosine]))/(4*s)/beta
            end do
            gamma_lab = s/2 - 1
            call qed_dirac_annihilation_total(gamma_lab, ref)
            call require(abs(total - ref) <= 1e-10_dp*ref, 'Dirac annihilation total by quadrature')
        end do
    end subroutine totals

    subroutine limits()
        real(dp) :: x, ref, beta, theta, half, rutherford, mott
        type(four_vector_t) :: p, pp
        ! Thomson: sigma_KN = (8 pi/3)(1 - 2x + O(x^2)); the two branches
        ! agree where they meet
        x = 1e-4_dp
        call qed_klein_nishina_total_low(x, ref)
        call require(abs(ref/(8*pi/3) - (1 - 2*x)) < 10*x**2, 'Thomson limit')
        call qed_klein_nishina_total_low(0.2_dp, ref)
        call qed_klein_nishina_total(0.2_dp, x)
        call require(abs(ref - x) <= 1e-13_dp*ref, 'Klein-Nishina branches agree at x = 0.2')
        ! Mott: amplitude against closed form, and Rutherford as beta -> 0
        beta = 0.5_dp
        theta = 1.1_dp
        p = on_shell(m, [0.0_dp, 0.0_dp, beta/sqrt(1 - beta**2)])
        pp = on_shell(m, direction(theta, 0.7_dp)*beta/sqrt(1 - beta**2))
        call qed_mott_dcs(beta, theta, 1.0_dp, ref)
        call require(close(mott_spin_sum(m, p, pp), ref), 'Mott closed form')
        beta = 1e-3_dp
        half = sin(theta/2)**2
        rutherford = 1/(4*beta**4*half**2)
        call qed_mott_dcs(beta, theta, 1.0_dp, mott)
        call require(abs(mott/rutherford - 1) < 2*beta**2, 'Rutherford limit')
    end subroutine limits

end program test_qed_tree
