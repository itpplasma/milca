module ideal_fb_gamma_p_test_support
    ! Fields of test_ideal_fb_gamma_p. source_field_t: grad u of a ring of
    ! nring point sources on R = rs, Z = 0 inside the plasma, u decaying and
    ! harmonic outside, so the exact vacuum response is Phi = -u and B = 0.
    ! coil_field_t: F/R e_phi + bz e_z + two external polygon loops.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_exterior_harmonic_fields, only: evaluate_current_loops
    use milca_ideal_fb_gamma_p, only: gamma_p_field_t
    implicit none
    private
    public :: ring, source_potential, torus_chain
    real(dp), parameter :: pi = acos(-1.0_dp)
    type, extends(gamma_p_field_t), public :: source_field_t
        real(dp) :: rs = 1, zs = 0
        integer :: nring = 360
    contains
        procedure :: points => source_points
    end type
    type, extends(gamma_p_field_t), public :: coil_field_t
        real(dp) :: f = 0.8_dp, bz = 0.05_dp, current(2) = [0.4_dp, -0.25_dp]
        real(dp), allocatable :: loops(:, :)
    contains
        procedure :: points => coil_points
    end type
contains
    function ring(r, z, n) result(v)
        real(dp), intent(in) :: r, z
        integer, intent(in) :: n
        real(dp) :: v(3, n)
        integer :: k
        do k = 1, n
            v(:, k) = [r*cos(2*pi*(k - 0.5_dp)/n), r*sin(2*pi*(k - 0.5_dp)/n), z]
        end do
    end function

    real(dp) function source_potential(field, x) result(u)
        type(source_field_t), intent(in) :: field
        real(dp), intent(in) :: x(3)
        real(dp) :: y(3, field%nring)
        integer :: k
        y = ring(field%rs, field%zs, field%nring)
        u = 0
        do k = 1, field%nring
            u = u + 1/(4*pi*norm2(x - y(:, k)))
        end do
        u = u/field%nring
    end function

    subroutine source_points(self, x, b, status)
        class(source_field_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :)
        integer, intent(out) :: status
        real(dp) :: y(3, self%nring), d(3)
        integer :: q, k
        y = ring(self%rs, self%zs, self%nring)
        b = 0
        do q = 1, size(x, 2)
            do k = 1, self%nring
                d = x(:, q) - y(:, k)
                b(:, q) = b(:, q) - d/(4*pi*norm2(d)**3)
            end do
        end do
        b = b/self%nring
        status = 0
    end subroutine

    subroutine coil_points(self, x, b, status)
        class(coil_field_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: potential(:, :, :), field(:, :, :)
        integer :: q, n
        n = size(self%loops, 2)/2
        allocate (potential(3, size(x, 2), 2), field(3, size(x, 2), 2))
        call evaluate_current_loops(self%loops, [1, n + 1, 2*n + 1], x, potential, &
                                    field, status)
        if (status /= 0) return
        do q = 1, size(x, 2)
            b(:, q) = self%current(1)*field(:, q, 1) + self%current(2)*field(:, q, 2) &
                + self%f*[-x(2, q), x(1, q), 0.0_dp]/(x(1, q)**2 + x(2, q)**2) &
                + [0.0_dp, 0.0_dp, self%bz]
        end do
    end subroutine

    subroutine torus_chain(nfp, nmode, nchain, chain, d)
        ! R = 1 + 0.3 cos t + d cos(t + nfp phi), Z = 0.39 sin t + d sin(t + nfp phi).
        integer, intent(in) :: nfp, nmode, nchain
        real(dp), allocatable, intent(out) :: chain(:, :, :)
        real(dp), intent(in) :: d
        real(dp) :: t
        integer :: j
        if (nfp < 1) error stop 'nfp'
        allocate (chain(2, nmode, nchain))
        chain = 0
        do j = 1, nchain
            t = 2*pi*(j - 1)/nchain
            chain(:, 1, j) = [1 + 0.3_dp*cos(t), 0.39_dp*sin(t)]
            if (nmode >= 3) then
                chain(:, 2, j) = d*[cos(t), sin(t)]
                chain(:, 3, j) = d*[-sin(t), cos(t)]
            end if
        end do
    end subroutine
end module ideal_fb_gamma_p_test_support

program test_ideal_fb_gamma_p
    ! Vacuum of the free boundary on Gamma_p (milca_ideal_fb_gamma_p).
    ! oracle-*: grad u of an interior ring source as B_e. The exact response
    ! cancels it, so E_bem = -1/2 int_ext |grad u|^2 = 1/2 int u du/dn dA,
    ! Phi = -u and B = 0 on Gamma_p; errors and orders under refinement.
    ! gradient-*: Hadamard gradient against central differences of the
    ! discrete E_v (coils, toroidal field and a filament circulation).
    ! hessian-*: Hessian product against differences of that gradient.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_triangle_basis, only: triangle_basis_jet
    use milca_ideal_fb_gamma_p
    use ideal_fb_gamma_p_test_support
    implicit none
    character(64) :: check
    call get_command_argument(1, check)
    select case (trim(check))
    case ('oracle-axisymmetric')
        call oracle(GAMMA_P_AXISYMMETRIC)
    case ('oracle-circulant')
        call oracle(GAMMA_P_CIRCULANT)
    case ('gradient-axisymmetric')
        call gradient(GAMMA_P_AXISYMMETRIC)
    case ('gradient-circulant')
        call gradient(GAMMA_P_CIRCULANT)
    case ('edge-basis')
        call edge_check()
    case ('hessian-axisymmetric')
        call hessian(GAMMA_P_AXISYMMETRIC)
    case ('hessian-circulant')
        call hessian(GAMMA_P_CIRCULANT)
    case default
        error stop 'unknown check'
    end select
contains
    subroutine build(kind, level, chain, nfp, gp, filament, r_ref)
        integer, intent(in) :: kind, level, nfp
        real(dp), intent(in) :: chain(:, :, :)
        type(gamma_p_t), intent(out) :: gp
        real(dp), intent(in), optional :: filament(:, :), r_ref
        type(cylindrical_chart_t) :: chart
        type(bem_pair_options_t) :: options
        integer :: st
        ! Two fixed pair rules, the second always accepted: the discrete
        ! operator is then a smooth function of the geometry, as the
        ! difference checks need, and the cost is bounded.
        options%initial_degree = 6
        options%maximum_degree = 10
        options%degree_increment = 4
        options%relative_tolerance = huge(1.0_dp)
        if (kind == GAMMA_P_AXISYMMETRIC) then
            call gamma_p_build(gp, chart, nfp, 3, chain, kind, st, &
                panels_per_edge=level, filament=filament, r_ref=r_ref)
        else
            call gamma_p_build(gp, chart, nfp, 3, chain, kind, st, &
                panels_per_edge=level, nphi=4*nfp*level, options=options, &
                filament=filament, r_ref=r_ref)
        end if
        if (st /= GAMMA_P_OK) then
            write (*, '(a,i0)') 'gamma_p_build status ', st
            error stop 1
        end if
    end subroutine

    subroutine oracle(kind)
        integer, intent(in) :: kind
        type(source_field_t) :: field
        type(gamma_p_t) :: gp
        real(dp), allocatable :: chain(:, :, :), g(:, :, :)
        real(dp) :: energy, exact, e_rel(3), b_rel(3), p_rel(3), u, du(3, 1), n(3)
        real(dp) :: bb, uu, pp, ph, order(3)
        integer :: level, nl, iq, st, nfp, nmode, lv(3), first
        nfp = 1; nmode = 1
        if (kind == GAMMA_P_CIRCULANT) then
            nfp = 3; nmode = 3; field%nring = 3*120
        end if
        lv = [1, 2, 4]
        first = 1
        if (kind == GAMMA_P_CIRCULANT) first = 2
        if (kind == GAMMA_P_CIRCULANT) lv = [1, 1, 2]
        do nl = first, 3
            level = lv(nl)
            call torus_chain(nfp, nmode, 12, chain, 0.04_dp)
            call build(kind, level, chain, nfp, gp)
            allocate (g(2, nmode, 12))
            call gamma_p_vacuum(gp, field, 0.0_dp, 1.0_dp, energy, g, st)
            if (st /= GAMMA_P_OK) error stop 'gamma_p_vacuum'
            exact = 0; bb = 0; uu = 0; pp = 0
            do iq = 1, gp%nq
                n = cross(gp%a(:, iq), gp%b(:, iq))
                u = source_potential(field, gp%x(:, iq))
                call field%points(gp%x(:, iq:iq), du, st)
                exact = exact + 0.5_dp*gp%w(iq)*u*dot_product(du(:, 1), n)
                bb = bb + gp%w(iq)*norm2(n)*sum(gp%field(:, iq)**2)
                uu = uu + gp%w(iq)*norm2(n)*(sum(du**2) + u**2)
                ph = dot_product(gp%phi(gp%dof(:, iq)), gp%basis(:, iq))
                pp = pp + gp%w(iq)*norm2(n)*(ph + u)**2
            end do
            e_rel(nl) = abs(gp%energy_parts(1) - exact)/abs(exact)
            b_rel(nl) = sqrt(bb/uu); p_rel(nl) = sqrt(pp/uu)
            write (*, '(a,i2,a,i6,a,es12.4,a,es10.3,a,es10.3,a,es10.3)') 'level ', level, &
                ' dofs ', gp%ntrace, ' E_exact ', exact, ' E rel err ', e_rel(nl), &
                ' |B| rel ', b_rel(nl), ' Phi+u rel ', p_rel(nl)
            deallocate (g)
        end do
        order = log([e_rel(2)/e_rel(3), b_rel(2)/b_rel(3), p_rel(2)/p_rel(3)]) &
            /log(real(lv(3), dp)/lv(2))
        write (*, '(a,3f7.2)') 'observed orders (energy, |B|, Phi) ', order
        if (e_rel(3) > 1.0e-3_dp .or. order(1) < 2) error stop 'energy not converging'
        if (b_rel(3) > 2.0e-2_dp .or. order(2) < 1) error stop 'field not converging'
    end subroutine

    subroutine gradient(kind)
        integer, intent(in) :: kind
        type(coil_field_t) :: field
        type(gamma_p_t) :: gp, plus, minus
        real(dp), allocatable :: chain(:, :, :), g(:, :, :), gs(:, :, :), filament(:, :)
        real(dp) :: energy, ep, em, h, fd, worst(3), scale
        integer :: level, nl, st, nfp, nmode, k, c, m, j, picks(3, 4), lv(3), first
        real(dp), parameter :: circulation = 0.25_dp
        nfp = 1; nmode = 1
        if (kind == GAMMA_P_CIRCULANT) then
            nfp = 3; nmode = 3
        end if
        field%loops = reshape([ring(1.9_dp, 0.5_dp, 600), ring(0.45_dp, -0.6_dp, 600)], &
                              [3, 1200])
        filament = ring(1.0_dp, 0.0_dp, 120)
        picks = reshape([1, 1, 1, 2, 1, 4, 1, 1, 7, 2, 1, 10], [3, 4])
        if (kind == GAMMA_P_CIRCULANT) picks(:, 3:4) = reshape([1, 2, 7, 2, 3, 10], [3, 2])
        lv = [1, 2, 4]
        first = 1
        if (kind == GAMMA_P_CIRCULANT) first = 2
        if (kind == GAMMA_P_CIRCULANT) lv = [1, 1, 2]
        h = 1.0e-4_dp
        do nl = first, 3
            level = lv(nl)
            call torus_chain(nfp, nmode, 12, chain, 0.04_dp)
            call build(kind, level, chain, nfp, gp, filament)
            allocate (g(2, nmode, 12), gs(2, nmode, 12))
            call gamma_p_vacuum(gp, field, circulation, 1.0_dp, energy, g, st)
            if (st /= GAMMA_P_OK) error stop 'gamma_p_vacuum'
            scale = maxval(abs(g))
            worst(nl) = 0
            do k = 1, size(picks, 2)
                c = picks(1, k); m = picks(2, k); j = picks(3, k)
                chain(c, m, j) = chain(c, m, j) + h
                call build(kind, level, chain, nfp, plus, filament, gp%r_ref)
                call gamma_p_vacuum(plus, field, circulation, 1.0_dp, ep, gs, st)
                chain(c, m, j) = chain(c, m, j) - 2*h
                call build(kind, level, chain, nfp, minus, filament, gp%r_ref)
                call gamma_p_vacuum(minus, field, circulation, 1.0_dp, em, gs, st)
                chain(c, m, j) = chain(c, m, j) + h
                fd = (ep - em)/(2*h)
                write (*, '(a,i2,a,3i3,a,es14.6,a,es14.6,a,es10.3)') 'level ', level, &
                    ' (c,m,j)', c, m, j, ' Hadamard ', g(c, m, j), ' FD ', fd, &
                    ' rel ', abs(fd - g(c, m, j))/scale
                worst(nl) = max(worst(nl), abs(fd - g(c, m, j))/scale)
            end do
            write (*, '(a,i2,a,i6,a,es10.3)') 'level ', level, ' dofs ', gp%ntrace, &
                ' worst |FD - Hadamard|/max|g| ', worst(nl)
            deallocate (g, gs)
        end do
        write (*, '(a,f7.2)') 'observed order ', log(worst(2)/worst(3)) &
            /log(real(lv(3), dp)/lv(2))
        if (worst(3) > 1.0e-2_dp .or. worst(3) > worst(first)) &
            error stop 'Hadamard gradient inconsistent with the discrete energy'
    end subroutine

    subroutine hessian(kind)
        ! Hessian product against central differences of the Hadamard gradient
        ! along a smooth chain direction.
        integer, intent(in) :: kind
        type(coil_field_t) :: field
        type(gamma_p_t) :: gp, plus, minus
        real(dp), allocatable :: chain(:, :, :), g(:, :, :), gp_(:, :, :), gm(:, :, :)
        real(dp), allocatable :: filament(:, :), direction(:, :, :), product(:, :, :)
        real(dp) :: energy, h, worst(3), t
        integer :: level, nl, st, nfp, nmode, j, m, lv(3), first
        real(dp), parameter :: circulation = 0.25_dp
        nfp = 1; nmode = 1
        if (kind == GAMMA_P_CIRCULANT) then
            nfp = 3; nmode = 3
        end if
        field%loops = reshape([ring(1.9_dp, 0.5_dp, 600), ring(0.45_dp, -0.6_dp, 600)], &
                              [3, 1200])
        filament = ring(1.0_dp, 0.0_dp, 120)
        allocate (direction(2, nmode, 12))
        do j = 1, 12
            t = 2*acos(-1.0_dp)*(j - 1)/12
            do m = 1, nmode
                direction(:, m, j) = [cos(2*t + m), 0.7_dp*sin(t - m)]/m
            end do
        end do
        lv = [1, 2, 4]
        first = 1
        if (kind == GAMMA_P_CIRCULANT) first = 2
        if (kind == GAMMA_P_CIRCULANT) lv = [1, 1, 2]
        h = 1.0e-4_dp
        do nl = first, 3
            level = lv(nl)
            call torus_chain(nfp, nmode, 12, chain, 0.04_dp)
            call build(kind, level, chain, nfp, gp, filament)
            allocate (g(2, nmode, 12), gp_(2, nmode, 12), gm(2, nmode, 12), &
                      product(2, nmode, 12))
            call gamma_p_vacuum(gp, field, circulation, 1.0_dp, energy, g, st)
            if (st /= GAMMA_P_OK) error stop 'gamma_p_vacuum'
            call gamma_p_hessian(gp, field, direction, product, st)
            if (st /= GAMMA_P_OK) error stop 'gamma_p_hessian'
            call build(kind, level, chain + h*direction, nfp, plus, filament, gp%r_ref)
            call gamma_p_vacuum(plus, field, circulation, 1.0_dp, energy, gp_, st)
            call build(kind, level, chain - h*direction, nfp, minus, filament, gp%r_ref)
            call gamma_p_vacuum(minus, field, circulation, 1.0_dp, energy, gm, st)
            worst(nl) = maxval(abs((gp_ - gm)/(2*h) - product))/maxval(abs(product))
            write (*, '(a,i2,a,i6,a,es12.4,a,es10.3)') 'level ', level, ' dofs ', &
                gp%ntrace, ' max|Hv| ', maxval(abs(product)), &
                ' max|FD - Hv|/max|Hv| ', worst(nl)
            deallocate (g, gp_, gm, product)
        end do
        write (*, '(a,f7.2)') 'observed order ', log(worst(2)/worst(3)) &
            /log(real(lv(3), dp)/lv(2))
        if (worst(3) > 2.0e-2_dp .or. worst(3) > worst(first)) &
            error stop 'Hessian product inconsistent with the gradient'
    end subroutine

    subroutine edge_check()
        ! The chain edge basis is the edge restriction of the plasma map's
        ! triangle basis (local nodes 1, 4, (5,) 2 along s2 = 0).
        real(dp) :: w(0:3), dw(0:3), values(10), gradient(2, 10), t, worst
        integer :: order, k, index(4), nloc
        worst = 0
        do order = 2, 3
            nloc = 6; index = [1, 4, 2, 0]
            if (order == 3) then
                nloc = 10; index = [1, 4, 5, 2]
            end if
            do k = 0, 10
                t = 0.1_dp*k
                call gamma_p_edge_basis(order, t, w(:order), dw(:order))
                call triangle_basis_jet([t, 0.0_dp], values(:nloc), gradient(:, :nloc))
                worst = max(worst, maxval(abs(w(:order) - values(index(:order + 1)))), &
                            maxval(abs(dw(:order) - gradient(1, index(:order + 1)))))
            end do
        end do
        write (*, '(a,es10.3)') 'edge basis against the triangle basis ', worst
        if (worst > 1.0e-12_dp) error stop 'edge basis'
    end subroutine

    pure function cross(u, v) result(c)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: c(3)
        c = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function
end program test_ideal_fb_gamma_p
