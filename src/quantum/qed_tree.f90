! Tree-level QED amplitudes from the Feynman rules, natural units with the
! electron mass m, metric (+,-,-,-), Dirac representation, e = 1. Each
! *_spin_sum returns (1/4) times the sum over all spins and photon
! polarizations of |M|^2, the quantity milca_qed_leaves gives in closed form.
! Photon polarizations are the two real transverse vectors; gauge invariance
! is a checked property (ward_compton), not an assumption.
module milca_qed_tree
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: four_vector_t, on_shell, minkowski, dirac_u, dirac_v, slash
    public :: transverse_polarizations, moller_spin_sum, bhabha_spin_sum
    public :: compton_spin_sum, annihilation_spin_sum, breit_wheeler_spin_sum
    public :: mott_spin_sum, ward_compton, spinor_completeness

    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp), one = (1.0_dp, 0.0_dp), &
        iu = (0.0_dp, 1.0_dp)
    real(dp), parameter :: metric(0:3) = [1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]

    !> Contravariant components p^mu.
    type :: four_vector_t
        real(dp) :: c(0:3) = 0
    end type four_vector_t

contains

    !> On-shell four-momentum of mass m and three-momentum p.
    pure function on_shell(m, p) result(v)
        real(dp), intent(in) :: m, p(3)
        type(four_vector_t) :: v
        v%c(0) = sqrt(m**2 + sum(p**2))
        v%c(1:3) = p
    end function on_shell

    pure real(dp) function minkowski(a, b)
        type(four_vector_t), intent(in) :: a, b
        minkowski = sum(metric*a%c*b%c)
    end function minkowski

    pure function gamma_matrix(mu) result(g)
        integer, intent(in) :: mu
        complex(dp) :: g(4, 4), sigma(2, 2)
        g = zero
        if (mu == 0) then
            g(1, 1) = one
            g(2, 2) = one
            g(3, 3) = -one
            g(4, 4) = -one
            return
        end if
        sigma = pauli(mu)
        g(1:2, 3:4) = sigma
        g(3:4, 1:2) = -sigma
    end function gamma_matrix

    pure function pauli(i) result(sigma)
        integer, intent(in) :: i
        complex(dp) :: sigma(2, 2)
        sigma = zero
        select case (i)
        case (1)
            sigma(1, 2) = one
            sigma(2, 1) = one
        case (2)
            sigma(1, 2) = -iu
            sigma(2, 1) = iu
        case (3)
            sigma(1, 1) = one
            sigma(2, 2) = -one
        end select
    end function pauli

    !> gamma^mu a_mu for a complex four-vector of contravariant components.
    pure function slash(a) result(g)
        complex(dp), intent(in) :: a(0:3)
        complex(dp) :: g(4, 4)
        integer :: mu
        g = zero
        do mu = 0, 3
            g = g + metric(mu)*a(mu)*gamma_matrix(mu)
        end do
    end function slash

    pure function sigma_dot(p) result(s)
        real(dp), intent(in) :: p(3)
        complex(dp) :: s(2, 2)
        s = p(1)*pauli(1) + p(2)*pauli(2) + p(3)*pauli(3)
    end function sigma_dot

    !> Positive-energy spinor, spin up (s = 1) or down (s = 2) along z in
    !> the rest frame; u-bar u = 2m.
    pure function dirac_u(m, p, s) result(u)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p
        integer, intent(in) :: s
        complex(dp) :: u(4), chi(2)
        chi = zero
        chi(s) = one
        u(1:2) = chi
        u(3:4) = matmul(sigma_dot(p%c(1:3)), chi)/(p%c(0) + m)
        u = sqrt(p%c(0) + m)*u
    end function dirac_u

    !> Negative-energy spinor; v-bar v = -2m.
    pure function dirac_v(m, p, s) result(v)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p
        integer, intent(in) :: s
        complex(dp) :: v(4), eta(2)
        eta = zero
        eta(3 - s) = one
        v(1:2) = matmul(sigma_dot(p%c(1:3)), eta)/(p%c(0) + m)
        v(3:4) = eta
        v = sqrt(p%c(0) + m)*v
    end function dirac_v

    pure function bar(psi) result(b)
        complex(dp), intent(in) :: psi(4)
        complex(dp) :: b(4)
        b = conjg(psi)
        b(3:4) = -b(3:4)
    end function bar

    !> Current a-bar gamma^mu b, contravariant.
    pure function current(a, b) result(j)
        complex(dp), intent(in) :: a(4), b(4)
        complex(dp) :: j(0:3)
        integer :: mu
        do mu = 0, 3
            j(mu) = dot_product(conjg(bar(a)), matmul(gamma_matrix(mu), b))
        end do
    end function current

    pure complex(dp) function contract(a, b)
        complex(dp), intent(in) :: a(0:3), b(0:3)
        contract = sum(metric*a*b)
    end function contract

    !> The two real transverse polarization four-vectors of a photon of
    !> momentum k.
    pure function transverse_polarizations(k) result(eps)
        type(four_vector_t), intent(in) :: k
        complex(dp) :: eps(0:3, 2)
        real(dp) :: n(3), e1(3), e2(3), trial(3)
        n = k%c(1:3)/norm2(k%c(1:3))
        trial = [1.0_dp, 0.0_dp, 0.0_dp]
        if (abs(n(1)) > 0.9_dp) trial = [0.0_dp, 1.0_dp, 0.0_dp]
        e1 = trial - dot_product(trial, n)*n
        e1 = e1/norm2(e1)
        e2 = [n(2)*e1(3) - n(3)*e1(2), n(3)*e1(1) - n(1)*e1(3), n(1)*e1(2) - n(2)*e1(1)]
        eps = zero
        eps(1:3, 1) = e1
        eps(1:3, 2) = e2
    end function transverse_polarizations

    ! (q-slash + m)/(q^2 - m^2)
    pure function propagator(m, q) result(sprop)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: q
        complex(dp) :: sprop(4, 4)
        integer :: i
        sprop = slash(cmplx(q%c, 0.0_dp, dp))
        do i = 1, 4
            sprop(i, i) = sprop(i, i) + m
        end do
        sprop = sprop/(minkowski(q, q) - m**2)
    end function propagator

    pure function sub(a, b) result(c)
        type(four_vector_t), intent(in) :: a, b
        type(four_vector_t) :: c
        c%c = a%c - b%c
    end function sub

    pure function add(a, b) result(c)
        type(four_vector_t), intent(in) :: a, b
        type(four_vector_t) :: c
        c%c = a%c + b%c
    end function add

    !> e-(p1) e-(p2) -> e-(p3) e-(p4): t and u channels with the relative
    !> sign of exchanged identical fermions.
    pure real(dp) function moller_spin_sum(m, p1, p2, p3, p4) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p1, p2, p3, p4
        complex(dp) :: amp, u1(4), u2(4), u3(4), u4(4)
        real(dp) :: t, u
        integer :: s1, s2, s3, s4
        t = minkowski(sub(p1, p3), sub(p1, p3))
        u = minkowski(sub(p1, p4), sub(p1, p4))
        total = 0
        do s1 = 1, 2
            u1 = dirac_u(m, p1, s1)
            do s2 = 1, 2
                u2 = dirac_u(m, p2, s2)
                do s3 = 1, 2
                    u3 = dirac_u(m, p3, s3)
                    do s4 = 1, 2
                        u4 = dirac_u(m, p4, s4)
                        amp = contract(current(u3, u1), current(u4, u2))/t - &
                            contract(current(u4, u1), current(u3, u2))/u
                        total = total + abs(amp)**2
                    end do
                end do
            end do
        end do
        total = total/4
    end function moller_spin_sum

    !> e-(p1) e+(p2) -> e-(p3) e+(p4): scattering (t) and annihilation (s)
    !> channels.
    pure real(dp) function bhabha_spin_sum(m, p1, p2, p3, p4) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p1, p2, p3, p4
        complex(dp) :: amp, u1(4), v2(4), u3(4), v4(4)
        real(dp) :: s, t
        integer :: s1, s2, s3, s4
        s = minkowski(add(p1, p2), add(p1, p2))
        t = minkowski(sub(p1, p3), sub(p1, p3))
        total = 0
        do s1 = 1, 2
            u1 = dirac_u(m, p1, s1)
            do s2 = 1, 2
                v2 = dirac_v(m, p2, s2)
                do s3 = 1, 2
                    u3 = dirac_u(m, p3, s3)
                    do s4 = 1, 2
                        v4 = dirac_v(m, p4, s4)
                        amp = contract(current(u3, u1), current(v2, v4))/t - &
                            contract(current(v2, u1), current(u3, v4))/s
                        total = total + abs(amp)**2
                    end do
                end do
            end do
        end do
        total = total/4
    end function bhabha_spin_sum

    ! u-bar(a) [eps2-slash S(q1) eps1-slash + eps1-slash S(q2) eps2-slash] b
    pure complex(dp) function two_photon(m, a, b, eps1, eps2, q1, q2)
        real(dp), intent(in) :: m
        complex(dp), intent(in) :: a(4), b(4), eps1(0:3), eps2(0:3)
        type(four_vector_t), intent(in) :: q1, q2
        complex(dp) :: op(4, 4)
        op = matmul(slash(eps2), matmul(propagator(m, q1), slash(eps1))) + &
            matmul(slash(eps1), matmul(propagator(m, q2), slash(eps2)))
        two_photon = dot_product(conjg(bar(a)), matmul(op, b))
    end function two_photon

    !> e(p) gamma(k) -> e(pp) gamma(kp).
    pure real(dp) function compton_spin_sum(m, p, k, pp, kp) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p, k, pp, kp
        complex(dp) :: eps(0:3, 2), epsp(0:3, 2), u(4), up(4)
        integer :: s, sp, l, lp
        eps = transverse_polarizations(k)
        epsp = transverse_polarizations(kp)
        total = 0
        do s = 1, 2
            u = dirac_u(m, p, s)
            do sp = 1, 2
                up = dirac_u(m, pp, sp)
                do l = 1, 2
                    do lp = 1, 2
                        total = total + abs(two_photon(m, up, u, eps(:, l), conjg(epsp(:, lp)), &
                            add(p, k), sub(p, kp)))**2
                    end do
                end do
            end do
        end do
        total = total/4
    end function compton_spin_sum

    !> The Compton amplitude with the incoming polarization replaced by k:
    !> zero by gauge invariance. Returns the largest modulus over spins.
    pure real(dp) function ward_compton(m, p, k, pp, kp) result(worst)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p, k, pp, kp
        complex(dp) :: epsp(0:3, 2), u(4), up(4)
        integer :: s, sp, lp
        epsp = transverse_polarizations(kp)
        worst = 0
        do s = 1, 2
            u = dirac_u(m, p, s)
            do sp = 1, 2
                up = dirac_u(m, pp, sp)
                do lp = 1, 2
                    worst = max(worst, abs(two_photon(m, up, u, cmplx(k%c, 0.0_dp, dp), &
                        conjg(epsp(:, lp)), add(p, k), sub(p, kp))))
                end do
            end do
        end do
    end function ward_compton

    !> e-(p1) e+(p2) -> gamma(k1) gamma(k2).
    pure real(dp) function annihilation_spin_sum(m, p1, p2, k1, k2) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p1, p2, k1, k2
        complex(dp) :: e1(0:3, 2), e2(0:3, 2), u1(4), v2(4)
        integer :: s1, s2, l1, l2
        e1 = transverse_polarizations(k1)
        e2 = transverse_polarizations(k2)
        total = 0
        do s1 = 1, 2
            u1 = dirac_u(m, p1, s1)
            do s2 = 1, 2
                v2 = dirac_v(m, p2, s2)
                do l1 = 1, 2
                    do l2 = 1, 2
                        total = total + abs(two_photon(m, v2, u1, conjg(e1(:, l1)), &
                            conjg(e2(:, l2)), sub(p1, k1), sub(p1, k2)))**2
                    end do
                end do
            end do
        end do
        total = total/4
    end function annihilation_spin_sum

    !> gamma(k1) gamma(k2) -> e-(p1) e+(p2).
    pure real(dp) function breit_wheeler_spin_sum(m, k1, k2, p1, p2) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: k1, k2, p1, p2
        complex(dp) :: e1(0:3, 2), e2(0:3, 2), u1(4), v2(4)
        integer :: s1, s2, l1, l2
        e1 = transverse_polarizations(k1)
        e2 = transverse_polarizations(k2)
        total = 0
        do s1 = 1, 2
            u1 = dirac_u(m, p1, s1)
            do s2 = 1, 2
                v2 = dirac_v(m, p2, s2)
                do l1 = 1, 2
                    do l2 = 1, 2
                        ! photon k1 at the electron end: u-bar eps1 S(p1-k1) eps2 v
                        total = total + abs(two_photon(m, u1, v2, e2(:, l2), e1(:, l1), &
                            sub(p1, k1), sub(p1, k2)))**2
                    end do
                end do
            end do
        end do
        total = total/4
    end function breit_wheeler_spin_sum

    !> Electron p -> pp on a static point charge Z = 1: (1/2) sum over spins
    !> of |u-bar(pp) gamma^0 u(p)|^2 / |q|^4.
    pure real(dp) function mott_spin_sum(m, p, pp) result(total)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p, pp
        complex(dp) :: j(0:3)
        real(dp) :: q2
        integer :: s, sp
        q2 = sum((pp%c(1:3) - p%c(1:3))**2)
        total = 0
        do s = 1, 2
            do sp = 1, 2
                j = current(dirac_u(m, pp, sp), dirac_u(m, p, s))
                total = total + abs(j(0))**2
            end do
        end do
        total = total/(2*q2**2)
    end function mott_spin_sum

    !> Largest deviation of sum_s u u-bar from p-slash + m and of
    !> sum_s v v-bar from p-slash - m.
    pure real(dp) function spinor_completeness(m, p) result(worst)
        real(dp), intent(in) :: m
        type(four_vector_t), intent(in) :: p
        complex(dp) :: su(4, 4), sv(4, 4), target_u(4, 4), target_v(4, 4), a(4), b(4)
        integer :: s, i, j
        su = zero
        sv = zero
        do s = 1, 2
            a = dirac_u(m, p, s)
            b = bar(a)
            do j = 1, 4
                do i = 1, 4
                    su(i, j) = su(i, j) + a(i)*b(j)
                end do
            end do
            a = dirac_v(m, p, s)
            b = bar(a)
            do j = 1, 4
                do i = 1, 4
                    sv(i, j) = sv(i, j) + a(i)*b(j)
                end do
            end do
        end do
        target_u = slash(cmplx(p%c, 0.0_dp, dp))
        target_v = target_u
        do i = 1, 4
            target_u(i, i) = target_u(i, i) + m
            target_v(i, i) = target_v(i, i) - m
        end do
        worst = max(maxval(abs(su - target_u)), maxval(abs(sv - target_v)))
    end function spinor_completeness

end module milca_qed_tree
