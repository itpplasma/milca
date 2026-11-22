module solovev_doublet
    !! Exact up-down symmetric Solov'ev doublet of the label-lane gates E2a and
    !! E2b (docs/eulerian-label-plan.md): psi = R^4/8 + sum_k a_k psi_k with the
    !! Cerfon-Freidberg homogeneous solutions (Phys. Plasmas 17, 032502
    !! (2010)), R0 = 1, f = f0. Magnetic axes O1, O2 at (1, +-y0), the X-point
    !! at (1, 0), label chi = (psi - psi_O)/(psi_X - psi_O): 0 at the axes, 1 on
    !! the separatrix, 2 on the wall. The toroidal flux per unit chi of each
    !! region, dPhi/dchi = f0 oint dl/(R |grad chi|), is tabulated in
    !! log10 |chi - 1| by level-curve tracing; the separatrix lobes are traced
    !! from X for separatrix-aligned layouts.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter, public :: doublet_f0 = 1.0_dp, doublet_y0 = 0.17_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! the X-point at (1, 0) and the axes at (1, +-y0) fix a(1:3) (depth
    ! psi_X - psi_O = 9.3e-4)
    real(dp), parameter :: a(6) = [1.2373456637_dp, -3.3088858713_dp, -2.0930442997_dp, &
                                   0.4209_dp, 0.1550_dp, 0.00694_dp]
    integer, parameter :: table_n = 400
    ! step bound from the label's second derivative (ray and level tracing)
    real(dp), parameter :: curvature = 3000

    type, public :: solovev_doublet_t
        real(dp) :: psi_o = 0, psi_x = 0, xo(2, 2) = 0, xx(2) = 0
        real(dp) :: tab_s(table_n) = 0, tab_lobe(table_n) = 0, tab_out(table_n) = 0
        real(dp) :: slope_o(2) = 0
        !> largest arc-length step of the level-curve tracing
        real(dp) :: max_step = 2e-3_dp
    contains
        procedure :: init => doublet_init
        procedure :: label => doublet_label
        procedure :: label_gradient => doublet_label_gradient
        procedure :: lobe_slope => doublet_lobe_slope
        procedure :: outer_slope => doublet_outer_slope
        procedure :: flux_slope => doublet_flux_slope
        procedure :: ray_radius => doublet_ray_radius
        procedure :: wall_samples => doublet_wall_samples
        procedure :: lobe_samples => doublet_lobe_samples
    end type solovev_doublet_t

    public :: doublet_psi

contains

    subroutine doublet_init(self)
        class(solovev_doublet_t), intent(inout) :: self
        integer :: k

        self%xx = critical([1.0_dp, 0.0_dp])
        self%xo(:, 1) = critical([1.0_dp, doublet_y0])
        self%xo(:, 2) = critical([1.0_dp, -doublet_y0])
        self%psi_x = doublet_psi(self%xx(1), self%xx(2))
        self%psi_o = doublet_psi(self%xo(1, 1), self%xo(2, 1))
        do k = 1, table_n
            self%tab_s(k) = -10 + 10.1_dp*(k - 1)/(table_n - 1)
            self%tab_lobe(k) = self%flux_slope(max(1 - 10**self%tab_s(k), 1e-3_dp), .true.)
            self%tab_out(k) = self%flux_slope(1 + 10**self%tab_s(k), .false.)
        end do
        self%slope_o = [self%flux_slope(1e-3_dp, .true.), self%flux_slope(2e-3_dp, .true.)]
    end subroutine doublet_init

    pure real(dp) function doublet_psi(x, y) result(psi)
        real(dp), intent(in) :: x, y
        real(dp) :: l

        l = log(x)
        psi = x**4/8 + a(1)*x**2 + a(2)*(y**2 - x**2*l) + a(3)*(x**4 - 4*x**2*y**2) &
              + a(4)*(2*y**4 - 9*y**2*x**2 + 3*x**4*l - 12*x**2*y**2*l) &
              + a(5)*(x**6 - 12*x**4*y**2 + 8*x**2*y**4) &
              + a(6)*(8*y**6 - 140*y**4*x**2 + 75*y**2*x**4 - 15*x**6*l &
                      + 180*x**4*y**2*l - 120*x**2*y**4*l)
    end function doublet_psi

    pure function grad_psi(x, y) result(gp)
        real(dp), intent(in) :: x, y
        real(dp) :: gp(2)
        real(dp), parameter :: h = 1e-6_dp

        gp = [doublet_psi(x + h, y) - doublet_psi(x - h, y), &
              doublet_psi(x, y + h) - doublet_psi(x, y - h)]/(2*h)
    end function grad_psi

    function critical(x0) result(x)
        !! Newton for grad psi = 0 from x0 (finite-difference Hessian).
        real(dp), intent(in) :: x0(2)
        real(dp) :: x(2), g(2), hm(2, 2), det
        real(dp), parameter :: h = 1e-5_dp
        integer :: it

        x = x0
        do it = 1, 30
            g = grad_psi(x(1), x(2))
            hm(:, 1) = (grad_psi(x(1) + h, x(2)) - grad_psi(x(1) - h, x(2)))/(2*h)
            hm(:, 2) = (grad_psi(x(1), x(2) + h) - grad_psi(x(1), x(2) - h))/(2*h)
            det = hm(1, 1)*hm(2, 2) - hm(1, 2)*hm(2, 1)
            x = x - [hm(2, 2)*g(1) - hm(1, 2)*g(2), -hm(2, 1)*g(1) + hm(1, 1)*g(2)]/det
        end do
    end function critical

    real(dp) function doublet_label(self, x, y) result(chi)
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: x, y

        chi = (doublet_psi(x, y) - self%psi_o)/(self%psi_x - self%psi_o)
    end function doublet_label

    function doublet_label_gradient(self, x, y) result(g)
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: x, y
        real(dp) :: g(2)

        g = grad_psi(x, y)/(self%psi_x - self%psi_o)
    end function doublet_label_gradient

    real(dp) function doublet_ray_radius(self, c, theta, chi) result(r)
        !! First radius along the ray from c where the label reaches chi. Each
        !! step solves v + G s + K s^2/2 = chi with G = |grad chi| and K a bound
        !! of the label's second derivative, so no crossing is passed; near the
        !! X-point the steps shrink to the width of the band.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: c(2), theta, chi
        real(dp) :: v, gl, e(2), x(2), step
        integer :: j

        e = [cos(theta), sin(theta)]
        r = 0
        do j = 1, 2000
            x = c + r*e
            v = self%label(x(1), x(2))
            if (chi - v < 1e-13_dp) exit
            gl = norm2(self%label_gradient(x(1), x(2)))
            step = (-gl + sqrt(gl**2 + 2*curvature*(chi - v)))/curvature
            r = r + step
        end do
    end function doublet_ray_radius

    real(dp) function doublet_flux_slope(self, chi, lobe) result(slope)
        !! dPhi/dchi = f0 oint dl/(R |grad chi|) along the level chi, traced
        !! with RK4 in arc length from the point on the ray straight away
        !! from the X-point; the step is limited by |grad chi|/K, so it
        !! shrinks near the X-point where the integrand peaks.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: chi
        logical, intent(in) :: lobe
        real(dp) :: c(2), y(3), k1(3), k2(3), k3(3), k4(3), h, turned, a0, a1, x0(2), g(2)
        integer :: it

        if (lobe) then
            c = self%xo(:, 1)
            x0 = c + self%ray_radius(c, pi/2, chi)*[0.0_dp, 1.0_dp]
        else
            c = self%xx
            x0 = c + self%ray_radius(c, 0.0_dp, chi)*[1.0_dp, 0.0_dp]
        end if
        y = [x0(1), x0(2), 0.0_dp]
        turned = 0
        a0 = atan2(y(2) - c(2), y(1) - c(1))
        do it = 1, 200000
            h = min(self%max_step, 0.05_dp*norm2(self%label_gradient(y(1), y(2)))/curvature)
            h = max(h, 1e-9_dp)
            k1 = level_rhs(y)
            k2 = level_rhs(y + h/2*k1)
            k3 = level_rhs(y + h/2*k2)
            k4 = level_rhs(y + h*k3)
            y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6
            a1 = atan2(y(2) - c(2), y(1) - c(1))
            turned = turned + modulo(a1 - a0 + pi, 2*pi) - pi
            a0 = a1
            if (abs(turned) >= 2*pi) exit
        end do
        ! the last step overshoots the start: remove the excess arc (straight)
        g = self%label_gradient(x0(1), x0(2))
        slope = doublet_f0*y(3) - doublet_f0*norm2(y(1:2) - x0)/(x0(1)*norm2(g))

    contains

        function level_rhs(y) result(f)
            !! d/dl of (x, integral of 1/(R |grad chi|)) along the level curve.
            real(dp), intent(in) :: y(3)
            real(dp) :: f(3), g(2), gl

            g = self%label_gradient(y(1), y(2))
            gl = norm2(g)
            f = [-g(2)/gl, g(1)/gl, 1/(y(1)*gl)]
        end function level_rhs

    end function doublet_flux_slope

    real(dp) function table(self, values, chi)
        !! Linear interpolation in log10 |chi - 1|, clamped at the far end.
        !! Below the first entry (|chi - 1| < 1e-10) the slope is extrapolated
        !! linearly in log10 |chi - 1| from the entries at 1e-10 and 1e-8: at a
        !! separatrix dPhi/dchi = A log|chi - 1| + B + o(1) (the passage time
        !! past the X-point), and the table follows that law to 1e-3 relative
        !! over 1e-8 ... 1e-2. Graded cells at X reach |chi - 1| < 1e-10.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: values(:), chi
        real(dp) :: s, u
        integer :: k, k8

        s = log10(max(abs(chi - 1), tiny(1.0_dp)))
        if (s < self%tab_s(1)) then
            k8 = minloc(abs(self%tab_s + 8), 1)
            table = values(1) + (s - self%tab_s(1))*(values(k8) - values(1)) &
                    /(self%tab_s(k8) - self%tab_s(1))
            return
        end if
        u = (s - self%tab_s(1))/(self%tab_s(2) - self%tab_s(1))
        k = min(max(int(u) + 1, 1), table_n - 1)
        u = min(max(u - (k - 1), 0.0_dp), 1.0_dp)
        table = (1 - u)*values(k) + u*values(k + 1)
    end function table

    real(dp) function doublet_lobe_slope(self, chi) result(slope)
        !! dPhi/dchi of a lobe; below the tabulated labels near O linear in chi.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: chi

        if (chi > 1e-3_dp) then
            slope = table(self, self%tab_lobe, chi)
        else
            slope = self%slope_o(1) + (chi - 1e-3_dp)*(self%slope_o(2) - self%slope_o(1))/1e-3_dp
        end if
    end function doublet_lobe_slope

    real(dp) function doublet_outer_slope(self, chi) result(slope)
        !! dPhi/dchi of the region between the separatrix and the wall.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), intent(in) :: chi

        slope = table(self, self%tab_out, chi)
    end function doublet_outer_slope

    subroutine doublet_wall_samples(self, samples)
        !! The level chi = 2, counter-clockwise from the outboard midplane, as
        !! rays from X.
        class(solovev_doublet_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: samples(:, :)
        integer, parameter :: m = 256
        real(dp) :: theta
        integer :: i

        allocate (samples(2, m))
        do i = 1, m
            theta = 2*pi*(i - 1)/m
            samples(:, i) = self%xx + self%ray_radius(self%xx, theta, 2.0_dp)* &
                            [cos(theta), sin(theta)]
        end do
    end subroutine doublet_wall_samples

    subroutine doublet_lobe_samples(self, lobe, step, samples)
        !! Separatrix loop of lobe 1 (above X) or 2 (below), from X
        !! counter-clockwise about its O-point back to X, first and last
        !! sample X. Lobe 1 leaves X along the upper null direction of the
        !! label Hessian with dR > 0 and follows the level chi = 1 by RK4 in
        !! arc length, each step projected back onto chi = 1 by Newton; it
        !! ends where it returns to within one step of X. Lobe 2 is its mirror
        !! image in Z = Z_X, reversed so that it runs counter-clockwise.
        class(solovev_doublet_t), intent(in) :: self
        integer, intent(in) :: lobe
        real(dp), intent(in) :: step
        real(dp), allocatable, intent(out) :: samples(:, :)
        real(dp), allocatable :: path(:, :)
        real(dp) :: hm(2, 2), d(2), x(2), k1(2), k2(2), k3(2), k4(2), beta, travelled
        real(dp), parameter :: h = 1e-5_dp
        integer :: n, i, it

        hm(:, 1) = (self%label_gradient(self%xx(1) + h, self%xx(2)) - &
                    self%label_gradient(self%xx(1) - h, self%xx(2)))/(2*h)
        hm(:, 2) = (self%label_gradient(self%xx(1), self%xx(2) + h) - &
                    self%label_gradient(self%xx(1), self%xx(2) - h))/(2*h)
        if (hm(1, 1)*hm(2, 2) - hm(1, 2)**2 >= 0) error stop 'doublet: X is not a saddle'
        ! null directions (cos b, sin b): H11 + 2 H12 tan b + H22 tan^2 b = 0
        beta = atan((-hm(1, 2) + sqrt(hm(1, 2)**2 - hm(1, 1)*hm(2, 2)))/hm(2, 2))
        d = [cos(beta), sin(beta)]
        if (d(2) < 0) d = -d
        if (d(1) < 0) then
            beta = atan((-hm(1, 2) - sqrt(hm(1, 2)**2 - hm(1, 1)*hm(2, 2)))/hm(2, 2))
            d = [cos(beta), sin(beta)]
            if (d(2) < 0) d = -d
        end if
        if (d(1) <= 0 .or. d(2) <= 0) error stop 'doublet: no upper-right separatrix branch'
        allocate (path(2, 100000))
        n = 1
        path(:, 1) = self%xx
        x = self%xx + step*d
        if (dot_product(tangent(x), d) <= 0) error stop 'doublet: lobe 1 is not counter-clockwise'
        travelled = step
        do it = 1, size(path, 2) - 2
            n = n + 1
            path(:, n) = x
            k1 = tangent(x)
            k2 = tangent(x + step/2*k1)
            k3 = tangent(x + step/2*k2)
            k4 = tangent(x + step*k3)
            x = x + step*(k1 + 2*k2 + 2*k3 + k4)/6
            x = on_separatrix(x)
            travelled = travelled + step
            if (travelled > 10*step .and. norm2(x - self%xx) < 1.5_dp*step) exit
        end do
        n = n + 1
        path(:, n) = self%xx
        if (lobe == 1) then
            samples = path(:, 1:n)
        else
            allocate (samples(2, n))
            do i = 1, n
                samples(:, i) = [path(1, n + 1 - i), 2*self%xx(2) - path(2, n + 1 - i)]
            end do
        end if

    contains

        function tangent(y) result(t)
            !! Unit tangent of the level curve, counter-clockwise about a
            !! label minimum.
            real(dp), intent(in) :: y(2)
            real(dp) :: t(2), g(2)

            g = self%label_gradient(y(1), y(2))
            t = [-g(2), g(1)]/norm2(g)
        end function tangent

        function on_separatrix(y0) result(y)
            real(dp), intent(in) :: y0(2)
            real(dp) :: y(2), g(2)
            integer :: j

            y = y0
            do j = 1, 4
                g = self%label_gradient(y(1), y(2))
                y = y - (self%label(y(1), y(2)) - 1)*g/sum(g**2)
            end do
        end function on_separatrix

    end subroutine doublet_lobe_samples

end module solovev_doublet
