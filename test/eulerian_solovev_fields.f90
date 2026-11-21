module solovev_label_fields
    !! Exact Solov'ev equilibrium for gate G2a, mu0 = 1:
    !! psi = alpha (R^2 - r0^2)^2 + beta R^2 Z^2, Delta* psi = (8 alpha + 2 beta) R^2,
    !! B = grad psi x grad phi + f0 grad phi, A = psi grad phi - f0 ln R e_z,
    !! J_phi = -(8 alpha + 2 beta) R, dp/dpsi = -(8 alpha + 2 beta).
    !! Labels chi = psi/psi_b, axis (r0, 0), wall chi = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_map
    use fortfem_spline_eval, only: spline_quadrature_t
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_evaluate
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    implicit none
    private

    real(dp), parameter, public :: r0 = 3.0_dp, alpha = 0.25_dp, beta = 1.0_dp, &
                                   f0 = 3.0_dp, psi_b = 6.25_dp
    real(dp), parameter, public :: pi = acos(-1.0_dp), source = 8*alpha + 2*beta
    integer, parameter :: rays = 128
    logical, public :: nonlinear_label = .false.
    type(spline_space_t), public :: space
    type(curve_chart_t), public :: chart

    type, public :: triplet_list
        !! Coordinate entries, summed when converted to CSR.
        integer :: count = 0
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
    contains
        procedure :: init => triplet_init
        procedure :: add => triplet_add
        procedure :: add_csr => triplet_add_csr
        procedure :: to_csr => triplet_to_csr
    end type triplet_list

    public :: psi, grad_psi, full, section_jet, potential, exact_b, metric2, det3, &
              ray_radius, enclosed_current, winding, label, label_inverse, &
              boundary_samples, l2_error, apply, curl

contains

    pure real(dp) function psi(r, z)
        real(dp), intent(in) :: r, z

        psi = alpha*(r**2 - r0**2)**2 + beta*r**2*z**2
    end function psi

    pure subroutine label_inverse(chi, psi_n, dpsi_n)
        !! psi_N(chi) and d psi_N/d chi for chi = psi_N or the nonlinear label
        !! chi = (psi_N + psi_N^2/2)/1.5.
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: psi_n, dpsi_n

        if (nonlinear_label) then
            psi_n = -1 + sqrt(1 + 3*chi)
            dpsi_n = 1.5_dp/sqrt(1 + 3*chi)
        else
            psi_n = chi
            dpsi_n = 1
        end if
    end subroutine label_inverse

    subroutine label(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: s(2), js(2, 2), pn

        call section_jet(block, xi, s, js)
        pn = psi(s(1), s(2))/psi_b
        values(1) = pn
        if (nonlinear_label) values(1) = (pn + pn**2/2)/1.5_dp
        associate (unused => t)
        end associate
    end subroutine label

    pure function grad_psi(r, z) result(g)
        real(dp), intent(in) :: r, z
        real(dp) :: g(2)

        g = [4*alpha*r*(r**2 - r0**2) + 2*beta*r*z**2, 2*beta*r**2*z]
    end function grad_psi

    pure real(dp) function ray_radius(level, theta) result(rho)
        !! Root of psi(axis + rho e_theta) = level by bisection (psi grows
        !! along every ray for these parameters).
        real(dp), intent(in) :: level, theta
        real(dp) :: lo, hi
        integer :: k

        lo = 0
        hi = 2
        do k = 1, 60
            rho = 0.5_dp*(lo + hi)
            if (psi(r0 + rho*cos(theta), rho*sin(theta)) < level) then
                lo = rho
            else
                hi = rho
            end if
        end do
    end function ray_radius

    pure real(dp) function enclosed_current(chi) result(current)
        !! I(chi) = int_{S(chi)} J_phi dA with J_phi = -source R, by rays from
        !! the axis (trapezoidal in theta, spectral).
        real(dp), intent(in) :: chi
        real(dp) :: rho, theta, area_r
        integer :: k

        area_r = 0
        do k = 1, rays
            theta = 2*pi*(k - 1)/rays
            rho = ray_radius(chi*psi_b, theta)
            area_r = area_r + r0*rho**2/2 + rho**3*cos(theta)/3
        end do
        current = -source*area_r*2*pi/rays
    end function enclosed_current

    pure function winding(r, z) result(h)
        !! grad theta about the axis in (R, Z) components.
        real(dp), intent(in) :: r, z
        real(dp) :: h(2)

        h = [-z, r - r0]/((r - r0)**2 + z**2)
    end function winding

    subroutine section_jet(block, xi, s, js)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: s(2), js(2, 2)
        real(dp) :: jet(2, 6)
        integer :: status

        call spline_geometry_at(space, block, xi, jet, status)
        s = jet(:, 1)
        js = jet(:, 2:3)
    end subroutine section_jet

    subroutine full(block, xi, t, x, j3)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: x(3), j3(3, 3)
        real(dp) :: s(2), js(2, 2), jc(3, 3)

        call section_jet(block, xi, s, js)
        call curve_chart_map(chart, s, t, x, jc)
        j3(:, 1:2) = matmul(jc(:, 1:2), js)
        j3(:, 3) = jc(:, 3)
    end subroutine full

    subroutine potential(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3), r, a(3)

        call full(block, xi, t, x, j3)
        r = hypot(x(1), x(2))
        a = psi(r, x(3))*[-x(2), x(1), 0.0_dp]/r**2 - f0*log(r)*[0.0_dp, 0.0_dp, 1.0_dp]
        values(1:3) = matmul(transpose(j3), a)
    end subroutine potential

    pure function exact_b(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), r, g(2), er(3), ephi(3)

        r = hypot(x(1), x(2))
        er = [x(1), x(2), 0.0_dp]/r
        ephi = [-x(2), x(1), 0.0_dp]/r
        g = grad_psi(r, x(3))
        b = (g(1)*[0.0_dp, 0.0_dp, 1.0_dp] - g(2)*er)/r + f0*ephi/r
    end function exact_b

    subroutine metric2(block, xi, t, weight)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: weight(:, :)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        weight(1:3, 1:3) = matmul(transpose(j3), j3)/abs(det3(j3))
    end subroutine metric2

    pure real(dp) function det3(j)
        real(dp), intent(in) :: j(3, 3)

        det3 = j(1, 1)*(j(2, 2)*j(3, 3) - j(3, 2)*j(2, 3)) - &
               j(1, 2)*(j(2, 1)*j(3, 3) - j(3, 1)*j(2, 3)) + &
               j(1, 3)*(j(2, 1)*j(3, 2) - j(3, 1)*j(2, 2))
    end function det3

    subroutine boundary_samples(samples)
        real(dp), allocatable, intent(out) :: samples(:, :)
        integer, parameter :: n = 256
        real(dp) :: theta, rho
        integer :: k

        allocate (samples(2, n))
        do k = 1, n
            theta = 2*pi*(k - 1)/n
            rho = ray_radius(psi_b, theta)
            samples(:, k) = [r0 + rho*cos(theta), rho*sin(theta)]
        end do
    end subroutine boundary_samples


    real(dp) function l2_error(derham, quad, b2) result(rel)
        !! || B_h - B || / || B ||; dV = w 2 pi R (axisymmetric, quad%weight
        !! holds the section Jacobian).
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: b2(:)
        real(dp) :: x(3), j3(3, 3), b(3), bh(3), proxy(3), dv, num, den
        integer :: q

        num = 0
        den = 0
        do q = 1, quad%point_count
            call full(quad%block(q), quad%xi(:, q), 0.0_dp, x, j3)
            b = exact_b(x)
            call spline_derham_evaluate(derham, 2, b2, quad%block(q), quad%xi(:, q), &
                                        0.0_dp, proxy)
            bh = matmul(j3, proxy)/det3(j3)
            dv = quad%weight(q)*2*pi*hypot(x(1), x(2))
            num = num + sum((bh - b)**2)*dv
            den = den + sum(b**2)*dv
        end do
        rel = sqrt(num/den)
    end function l2_error


    function apply(n, rowptr, colidx, val, x) result(y)
        integer, intent(in) :: n, rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:), x(:)
        real(dp), allocatable :: y(:)
        integer :: i, p

        allocate (y(n))
        do i = 1, n
            y(i) = 0
            do p = rowptr(i), rowptr(i + 1) - 1
                y(i) = y(i) + val(p)*x(colidx(p))
            end do
        end do
    end function apply


    function curl(derham, a) result(b)
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: b(:)
        integer :: i, p

        allocate (b(derham%dimension(2)))
        do i = 1, derham%curl%n
            b(i) = 0
            do p = derham%curl%rowptr(i), derham%curl%rowptr(i + 1) - 1
                b(i) = b(i) + derham%curl%val(p)*a(derham%curl%colidx(p))
            end do
        end do
    end function curl


    subroutine triplet_init(self, capacity)
        class(triplet_list), intent(inout) :: self
        integer, intent(in) :: capacity

        self%count = 0
        allocate (self%rows(max(capacity, 1)), self%cols(max(capacity, 1)), &
                  self%vals(max(capacity, 1)))
    end subroutine triplet_init

    subroutine triplet_add(self, row, col, val)
        class(triplet_list), intent(inout) :: self
        integer, intent(in) :: row, col
        real(dp), intent(in) :: val

        if (self%count == size(self%rows)) then
            self%rows = [self%rows, self%rows]
            self%cols = [self%cols, self%cols]
            self%vals = [self%vals, self%vals]
        end if
        self%count = self%count + 1
        self%rows(self%count) = row
        self%cols(self%count) = col
        self%vals(self%count) = val
    end subroutine triplet_add

    subroutine triplet_add_csr(self, a, row0, col0)
        class(triplet_list), intent(inout) :: self
        type(csr_t), intent(in) :: a
        integer, intent(in) :: row0, col0
        integer :: i, p

        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                call self%add(row0 + i, col0 + a%col(p), a%val(p))
            end do
        end do
    end subroutine triplet_add_csr

    subroutine triplet_to_csr(self, nrow, ncol, a)
        !! The CSC form of the transpose is the CSR form of the matrix.
        class(triplet_list), intent(in) :: self
        integer, intent(in) :: nrow, ncol
        type(csr_t), intent(out) :: a
        type(csc_t) :: m
        type(fortsparse_status_t) :: info

        call csc_from_triplet(ncol, nrow, self%cols(:self%count), self%rows(:self%count), &
                              self%vals(:self%count), m, info)
        if (info%code /= 0) error stop 'triplets'
        a%nrow = nrow
        a%ncol = ncol
        a%row_ptr = m%col_ptr
        a%col = m%row_idx
        a%val = m%val
    end subroutine triplet_to_csr

end module solovev_label_fields
