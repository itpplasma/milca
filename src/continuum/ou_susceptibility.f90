!> Susceptibility functions of the linear drift-kinetic response with the
!! Ornstein-Uhlenbeck collision model (Heyn et al. 2014; Markl et al. 2023,
!! 2026), the generalised plasma dispersion functions of the integral
!! cylinder response.
!!
!! In units of the thermal speed v_T and the collision frequency nu, the
!! velocity-space Green's function G of
!!     (-i x2 + i x1 u - L) G(u, u') = delta(u - u'),
!!     L g = d/du (dg/du + u g),
!! with x1 = k_par v_T/nu and x2 = (omega - omega_E)/nu, gives
!!     I^{kl}(x1, x2) = int du du' u^k G(u, u') u'^l exp(-u'^2/2)/sqrt(2 pi).
!! L is the generator of the stationary Ornstein-Uhlenbeck process with unit
!! variance and correlation exp(-|t|), so by Feynman-Kac
!!     I^{kl} = int_0^inf dt exp(i x2 t) E[u(t)^k u(0)^l exp(-i x1 X)],
!!     X = int_0^t u ds.
!! The process is Gaussian: Var X = 2 (t - 1 + s) with s = exp(-t), and
!! tilting by exp(-i x1 X) shifts u(0) and u(t) by mu = -i x1 (1 - s) at
!! unchanged covariance (correlation s). Hence
!!     I^{kl} = int_0^inf dt exp(i x2 t - x1^2 (t - 1 + s)) P_kl(s),
!!     P_kl(s) = E[(b + mu)^k (a + mu)^l],
!! (a, b) standard normal with correlation s. With p = x1^2 - i x2 and
!! z = x1^2 the substitution s = exp(-t) and Kummer's integral give for each
!! power s^j of P_kl
!!     int_0^1 s^(p + j - 1) exp(z (1 - s)) ds = M(1, p + j + 1, z)/(p + j).
!! That sum cancels like x1^(k + l) for large x1, where the time integral is
!! short (the Gaussian decay exp(-x1^2 t^2/2)) and is evaluated directly; where
!! it would need too many panels (|x2| >> x1^2) and x1 is moderate, the Kummer
!! sum is kept, losing about x1^(k + l) of the double precision.
module milca_ou_susceptibility
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1_a1
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    public :: ou_susceptibility, ou_conserving, ou_polynomial
    integer, parameter, public :: OU_OK = 0, OU_SINGULAR = 1, OU_NOT_CONVERGED = 2
    ! Conservation added by the integral part of the collision operator.
    integer, parameter, public :: OU_PARTICLES = 0, OU_ENERGY = 1, OU_MOMENTUM = 2, &
        OU_ENERGY_MOMENTUM = 3

    ! Kummer sums below, the time integral from this |x1| on, and the Kummer
    ! sum as fallback below x1_kummer.
    real(dp), parameter :: x1_direct = 3, x1_kummer = 30
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

contains

    !> I^{kl}(x1, x2) for 0 <= k, l <= ubound(imn); OU_SINGULAR at x1 = x2 = 0,
    !> where the response has no damping.
    subroutine ou_susceptibility(x1, x2, imn, status)
        real(dp), intent(in) :: x1, x2
        complex(dp), intent(out) :: imn(0:, 0:)
        integer, intent(out) :: status
        complex(dp), allocatable :: c(:, :, :)
        integer :: n
        imn = 0
        status = OU_SINGULAR
        if (x1 == 0 .and. x2 == 0) return
        n = ubound(imn, 1)
        if (abs(x1) < x1_direct) then
            call ou_polynomial(x1, n, c)
            call kummer(x1, x2, c, imn, status)
        else
            call direct(x1, x2, imn, status)
            if (status == OU_NOT_CONVERGED .and. abs(x1) < x1_kummer) then
                call ou_polynomial(x1, n, c)
                call kummer(x1, x2, c, imn, status)
            end if
        end if
    end subroutine ou_susceptibility

    !> Coefficients c(j, k, l) of P_kl(s) = sum_j c(j, k, l) s^j, 0 <= k, l <= n,
    !> from the binomial expansion in mu = -i x1 (1 - s) and the moments of the
    !> correlated normal pair,
    !>   E[a^i b^j] = sum_r C(i, r) C(j, r) r! s^r E[Z^(i-r)] E[Z^(j-r)].
    pure subroutine ou_polynomial(x1, n, c)
        real(dp), intent(in) :: x1
        integer, intent(in) :: n
        complex(dp), allocatable, intent(out) :: c(:, :, :)
        complex(dp) :: mu(0:2*n, 0:2*n), m(0:2*n), t(0:2*n)
        integer :: k, l, i, j, r, p
        allocate (c(0:2*n, 0:n, 0:n))
        c = 0
        ! mu^p = (-i x1)^p (1 - s)^p as polynomials in s.
        mu = 0
        do p = 0, 2*n
            do j = 0, p
                mu(j, p) = (-iu*x1)**p*binomial(p, j)*(-1)**j
            end do
        end do
        do k = 0, n
            do l = 0, n
                do i = 0, l
                    do j = 0, k
                        ! E[a^i b^j] in s.
                        m = 0
                        do r = 0, min(i, j)
                            m(r) = m(r) + binomial(i, r)*binomial(j, r)*factorial(r) &
                                *normal_moment(i - r)*normal_moment(j - r)
                        end do
                        t = product_poly(mu(:, l - i + k - j), m)
                        c(:, k, l) = c(:, k, l) + binomial(l, i)*binomial(k, j)*t
                    end do
                end do
            end do
        end do
    end subroutine ou_polynomial

    !> I^{kl} = sum_j c_j M(1, p + j + 1, z)/(p + j), p = x1^2 - i x2, z = x1^2.
    subroutine kummer(x1, x2, c, imn, status)
        real(dp), intent(in) :: x1, x2
        complex(dp), intent(in) :: c(0:, 0:, 0:)
        complex(dp), intent(out) :: imn(0:, 0:)
        integer, intent(out) :: status
        complex(dp) :: f(0:ubound(c, 1)), p, z, m
        type(fortnum_status_t) :: st
        integer :: j, k, l
        status = OU_NOT_CONVERGED
        z = cmplx(x1*x1, 0.0_dp, dp)
        p = cmplx(x1*x1, -x2, dp)
        do j = 0, ubound(c, 1)
            if (j == 0 .and. abs(p) == 0) then
                status = OU_SINGULAR
                return
            end if
            call hyperg_1f1_a1(p + j + 1, z, m, st)
            if (.not. status_ok(st)) return
            f(j) = m/(p + j)
        end do
        do l = 0, ubound(imn, 2)
            do k = 0, ubound(imn, 1)
                imn(k, l) = sum(c(:, k, l)*f)
            end do
        end do
        status = OU_OK
    end subroutine kummer

    !> The time integral on Gauss-Legendre panels, each bisected until its
    !> 16- and 32-point values agree to 1e-14 of the running total; it ends
    !> where exp(-x1^2 (t - 1 + exp(-t))) < 1e-18.
    subroutine direct(x1, x2, imn, status)
        real(dp), intent(in) :: x1, x2
        complex(dp), intent(out) :: imn(0:, 0:)
        integer, intent(out) :: status
        integer, parameter :: max_panels = 50000, max_evaluations = 400000
        real(dp) :: x16(16), w16(16), x32(32), w32(32), t_end, lo(max_panels), hi(max_panels)
        complex(dp) :: q16(0:ubound(imn, 1), 0:ubound(imn, 2)), &
            q32(0:ubound(imn, 1), 0:ubound(imn, 2))
        real(dp) :: scale, h
        integer :: n, top, np, evaluations
        status = OU_NOT_CONVERGED
        evaluations = 0
        call gauss_legendre(16, x16, w16)
        call gauss_legendre(32, x32, w32)
        ! Gaussian decay up to t ~ 1, then exponential at rate x1^2.
        t_end = sqrt(2*log(1e18_dp))/abs(x1)
        do while (x1*x1*(t_end - 1 + exp(-t_end)) < log(1e18_dp))
            t_end = 2*t_end
        end do
        imn = 0
        ! Panels of at most a quarter oscillation and a tenth of the decay length.
        h = min(0.1_dp/abs(x1), 1.0_dp)
        if (x2 /= 0) h = min(h, 1.5_dp/abs(x2))
        np = ceiling(t_end/h)
        if (np > max_panels) return
        do n = 1, np
            lo(np - n + 1) = (n - 1)*t_end/np
            hi(np - n + 1) = n*t_end/np
        end do
        top = np
        scale = 0
        do while (top > 0)
            evaluations = evaluations + 1
            if (evaluations > max_evaluations) return
            call panel(lo(top), hi(top), x16, w16, q16)
            call panel(lo(top), hi(top), x32, w32, q32)
            scale = max(scale, maxval(abs(q32)))
            if (maxval(abs(q32 - q16)) <= 1e-14_dp*max(scale, maxval(abs(imn))) &
                    .or. hi(top) - lo(top) < 1e-12_dp) then
                imn = imn + q32
                top = top - 1
            else
                if (top + 1 > max_panels) return
                lo(top + 1) = (lo(top) + hi(top))/2
                hi(top + 1) = hi(top)
                hi(top) = lo(top + 1)
                ! Work on the upper half first, keep the lower half below it.
                call swap(top, top + 1)
                top = top + 1
            end if
        end do
        status = OU_OK
    contains
        subroutine panel(a, b, x, w, q)
            real(dp), intent(in) :: a, b, x(:), w(:)
            complex(dp), intent(out) :: q(0:, 0:)
            real(dp) :: t, s
            complex(dp) :: e
            integer :: i
            q = 0
            do i = 1, size(x)
                t = (a + b)/2 + (b - a)/2*x(i)
                s = exp(-t)
                ! 1 - s and t - 1 + s without cancellation at small t.
                e = exp(cmplx(-x1*x1*tail(t), x2*t, dp))*w(i)*(b - a)/2
                q = q + e*moments(-iu*x1*(t - tail(t)), s, ubound(q, 1))
            end do
        end subroutine panel
        subroutine swap(i, j)
            integer, intent(in) :: i, j
            real(dp) :: a, b
            a = lo(i); b = hi(i)
            lo(i) = lo(j); hi(i) = hi(j)
            lo(j) = a; hi(j) = b
        end subroutine swap
    end subroutine direct

    !> Adds the integral part of the collision operator that restores the
    !> conservation laws the Ornstein-Uhlenbeck part breaks. Particles are
    !> conserved by L itself. The energy part nu phi(u) (u^2 - 1) <u^2 - 1, f>
    !> and the momentum part nu phi(u) u <u, f> are rank-one sources in the
    !> kinetic equation, so the conserving functions follow from the Woodbury
    !> identity on I: with e = (u^2 - 1) and d = u,
    !>   energy:    I + I e e^T I/(1 - e^T I e),
    !>   momentum:  I + I d d^T I/(1 - d^T I d),
    !> and both together with the 2 x 2 capacitance matrix of (e, d). Needs
    !> ubound(imn) >= 2.
    pure subroutine ou_conserving(imn, model, out)
        complex(dp), intent(in) :: imn(0:, 0:)
        integer, intent(in) :: model
        complex(dp), intent(out) :: out(0:, 0:)
        complex(dp) :: ie(0:ubound(imn, 1)), ei(0:ubound(imn, 2)), id(0:ubound(imn, 1)), &
            di(0:ubound(imn, 2)), eie, did, eid, die, det
        integer :: k, l
        ! (I e)_k = I^{k2} - I^{k0}, (e^T I)_l = I^{2l} - I^{0l}.
        ie = imn(:, 2) - imn(:, 0)
        ei = imn(2, :) - imn(0, :)
        id = imn(:, 1)
        di = imn(1, :)
        eie = imn(2, 2) - imn(2, 0) - imn(0, 2) + imn(0, 0)
        did = imn(1, 1)
        eid = imn(2, 1) - imn(0, 1)
        die = imn(1, 2) - imn(1, 0)
        out = imn
        select case (model)
        case (OU_ENERGY)
            do l = 0, ubound(out, 2)
                do k = 0, ubound(out, 1)
                    out(k, l) = imn(k, l) + ie(k)*ei(l)/(1 - eie)
                end do
            end do
        case (OU_MOMENTUM)
            do l = 0, ubound(out, 2)
                do k = 0, ubound(out, 1)
                    out(k, l) = imn(k, l) + id(k)*di(l)/(1 - did)
                end do
            end do
        case (OU_ENERGY_MOMENTUM)
            ! (1 - V^T I U)^-1 with U = V = [e, d].
            det = (1 - eie)*(1 - did) - eid*die
            do l = 0, ubound(out, 2)
                do k = 0, ubound(out, 1)
                    out(k, l) = imn(k, l) + ((1 - did)*ie(k)*ei(l) + eid*ie(k)*di(l) &
                        + die*id(k)*ei(l) + (1 - eie)*id(k)*di(l))/det
                end do
            end do
        end select
    end subroutine ou_conserving

    !> P_kl = E[(b + mu)^k (a + mu)^l] for 0 <= k, l <= n at correlation s,
    !> evaluated directly so that no powers of x1 cancel.
    pure function moments(mu, s, n) result(p)
        complex(dp), intent(in) :: mu
        real(dp), intent(in) :: s
        integer, intent(in) :: n
        complex(dp) :: p(0:n, 0:n), mp(0:2*n)
        real(dp) :: m(0:n, 0:n)
        integer :: i, j, k, l, r
        mp(0) = 1
        do i = 1, 2*n
            mp(i) = mp(i - 1)*mu
        end do
        do i = 0, n
            do j = 0, n
                m(i, j) = 0
                do r = 0, min(i, j)
                    m(i, j) = m(i, j) + binomial(i, r)*binomial(j, r)*factorial(r) &
                        *normal_moment(i - r)*normal_moment(j - r)*s**r
                end do
            end do
        end do
        p = 0
        do l = 0, n
            do k = 0, n
                do i = 0, l
                    do j = 0, k
                        p(k, l) = p(k, l) + binomial(l, i)*binomial(k, j)*mp(l - i + k - j) &
                            *m(i, j)
                    end do
                end do
            end do
        end do
    end function moments

    !> t - 1 + exp(-t), by its series below t = 0.1.
    pure real(dp) function tail(t)
        real(dp), intent(in) :: t
        real(dp) :: term
        integer :: k
        if (t > 0.1_dp) then
            tail = t - 1 + exp(-t)
            return
        end if
        ! t^2/2 - t^3/6 + t^4/24 - ...
        tail = 0
        term = t
        do k = 2, 20
            term = -term*t/k
            tail = tail - term
        end do
    end function tail

    pure function product_poly(a, b) result(r)
        complex(dp), intent(in) :: a(0:), b(0:)
        complex(dp) :: r(0:ubound(a, 1))
        integer :: i, j
        r = 0
        do i = 0, ubound(a, 1)
            do j = 0, ubound(b, 1) - i
                r(i + j) = r(i + j) + a(i)*b(j)
            end do
        end do
    end function product_poly

    pure real(dp) function binomial(n, k)
        integer, intent(in) :: n, k
        binomial = factorial(n)/(factorial(k)*factorial(n - k))
    end function binomial

    pure real(dp) function factorial(n)
        integer, intent(in) :: n
        integer :: i
        factorial = 1
        do i = 2, n
            factorial = factorial*i
        end do
    end function factorial

    !> E[Z^n] of a standard normal: (n - 1)!! for even n, zero for odd n.
    pure real(dp) function normal_moment(n)
        integer, intent(in) :: n
        integer :: i
        normal_moment = 0
        if (mod(n, 2) /= 0) return
        normal_moment = 1
        do i = n - 1, 1, -2
            normal_moment = normal_moment*i
        end do
    end function normal_moment

end module milca_ou_susceptibility
