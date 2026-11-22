program test_eulerian_clebsch_free_boundary
    !! Gate E6a of docs/eulerian-label-plan.md: the axisymmetric label lane
    !! with a free plasma boundary Gamma_p, coupled to a vacuum shell up to a
    !! fixed surface Gamma_c and to the exterior through the Steklov-Poincare
    !! map on Gamma_c (FortFEM BEM, milca_gs_spline vacuum terms).
    !!
    !! Manufactured exact solution: the Solov'ev field of E1 (test
    !! eulerian_solovev_fields, mu0 = 1) everywhere inside Gamma_c,
    !! psi = alpha (R^2 - r0^2)^2 + beta R^2 Z^2, B = grad psi x grad phi +
    !! f0 grad phi. Inside Gamma_p = {psi = psi_b} it is the label-lane
    !! equilibrium of E1 (p = source psi_b (1 - chi) vanishes on Gamma_p).
    !! Outside, psi is not harmonic: the shell carries the known toroidal
    !! current density of a distributed coil, the source f = -source R of the
    !! shell energy int |grad psi|^2/(2R) - f psi, and the exterior sees the
    !! coil data psi_c = psi, d_n psi_c = d_n psi on Gamma_c. B is continuous
    !! across Gamma_p, so the boundary balance p + B^2/2 holds there. Every
    !! unknown below has its exact value.
    !!
    !! Unknowns: the plasma labels chi (1 on Gamma_p) and g in V0 of an
    !! O-grid whose wall is Gamma_p and whose core follows the axis X (E1),
    !! the multiplier mu of chi(X) = 0; the shell flux psi in V0 of an O-grid
    !! with an annulus between Gamma_p (interface) and Gamma_c (wall), whose
    !! functions inside Gamma_p are tied to one value psi_b (psi is constant
    !! there, so they add no energy); Gamma_p as R = r0 + r(theta) (cos, sin)
    !! with r a Fourier series of 2M + 1 coefficients rho. Equations:
    !!
    !!   label lane (E1): grad W/(2 pi) + mu N(X) = 0, grad chi(X) = 0,
    !!     chi(X) = 0;
    !!   shell: the stationarity of its energy plus the exterior terms
    !!     T^T S T psi - T^T S psi_c - oint d_n psi_c/R v (free functions);
    !!   psi_b: the current through Gamma_p from both sides,
    !!     oint G'(chi) d_n chi/R dl = oint d_n psi/R dl;
    !!   rho: the pressure balance tested with the Fourier modes w_m(theta),
    !!     oint (p + |B_in|^2/2) w_m dl = oint |B_out|^2/2 w_m dl,
    !!     |B_out|^2 = (|grad psi|^2 + f0^2)/R^2, each side on its own mesh.
    !!
    !! The X and rho columns of the Newton matrix are central differences with
    !! both meshes rebuilt (ALE: fields move with the mesh); the trace rows are
    !! differentiated by the chain rule. B, psi in the shell, the boundary and
    !! the axis are compared with the exact solution from a perturbed start.
    !!
    !! Arguments: levels (default 2: 6, 12 cells), Fourier order M (default 10),
    !! perturbation eps of the start (default 0.01: Gamma_p radius and psi_b by
    !! eps, the m = 1 boundary modes by eps and eps/2 of the radius, the axis
    !! by eps (0.6, 0.3), the labels by 1.5 eps chi (1 - chi)), 'jacobian', the
    !! cells of the first level (default 6, even: the axis is a vertex).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_geometry_at, spline_quadrature_point_singular, spline_locate, &
        spline_side_quadrature_t, spline_side_quadrature, spline_interface_quadrature, &
        spline_point_matrix
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_evaluate
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortsparse_csr, only: csr_t
    use milca_eulerian_vacuum, only: eu_constrained_solve
    use milca_eulerian_clebsch, only: clebsch_axisym_system, clebsch_node_gradients
    use milca_gs_spline, only: gs_spline_kernel_t, gs_spline_vacuum_t, gs_spline_vacuum_setup, &
        gs_spline_project
    use solovev_label_fields, only: r0, alpha, beta, f0, psi_b, source, pi, psi, grad_psi, &
        ray_radius, solovev_profiles, triplet_list
    implicit none
    integer, parameter :: degree = 2, nsample = 128
    real(dp), parameter :: ridge = 1e-12_dp, newton_ridge = 1e-8_dp
    integer :: levels = 2, nmode = 10, max_steps = 16, level, failures, nrho
    integer :: cells(3) = [6, 12, 24]
    real(dp) :: perturbation = 0.01_dp
    logical :: jacobian_check = .false.
    real(dp), allocatable :: err_b(:), err_psi(:), err_gamma(:), err_axis(:)
    integer, allocatable :: steps(:)
    character(16) :: arg

    ! the current meshes (rebuilt for every axis and boundary)
    type(spline_space_t), target :: space_p, space_f
    type(spline_derham_t) :: derham_p
    type(spline_quadrature_t) :: aquad_p, quad_f
    type(gs_spline_vacuum_t) :: vac
    integer :: cells_built = 0
    ! dof sets: plasma wall (chi = 1), shell functions tied to psi_b, free ones
    integer, allocatable :: wall_p(:), tied(:), free(:)
    real(dp), allocatable :: psi_c(:), dn_psi_c(:)
    integer :: n0p, nfull, nfree
    ! the previous (coarse) level: its spaces and solution, for the start of
    ! the next level (mesh continuation)
    type(spline_space_t), target :: coarse_p, coarse_f
    real(dp), allocatable :: coarse_u(:), coarse_psi(:)
    integer :: coarse_n0p = 0, coarse_nfree = 0
    real(dp) :: coarse_x(2)
    integer :: transfer_field = 1

    failures = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nmode
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) perturbation
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        jacobian_check = arg == 'jacobian'
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) cells(1)
        cells(2:) = [2*cells(1), 4*cells(1)]
    end if
    nrho = 2*nmode + 1
    print '(a,i3,a,f6.3)', 'Fourier order of Gamma_p', nmode, '  boundary perturbation', &
        perturbation
    print '(a,es10.2)', 'Fourier truncation of the exact Gamma_p (max radius error)', &
        truncation_error()
    allocate (err_b(levels), err_psi(levels), err_gamma(levels), err_axis(levels), &
              steps(levels))
    do level = 1, levels
        call solve(cells(level), err_b(level), err_psi(level), err_gamma(level), &
                   err_axis(level), steps(level))
        print '(a,i3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,i3)', 'cells', cells(level), &
            '  B (plasma)', err_b(level), '  psi (shell)', err_psi(level), '  Gamma_p', &
            err_gamma(level), '  axis', err_axis(level), '  Newton steps', steps(level)
        call check(steps(level) <= 10, 'Newton converges in at most 10 steps')
    end do
    do level = 2, levels
        print '(a,4f7.2)', 'rates of B, psi, Gamma_p, axis', &
            log(err_b(level - 1)/err_b(level))/log(2.0_dp), &
            log(err_psi(level - 1)/err_psi(level))/log(2.0_dp), &
            log(err_gamma(level - 1)/err_gamma(level))/log(2.0_dp), &
            log(err_axis(level - 1)/err_axis(level))/log(2.0_dp)
        call check(log(err_b(level - 1)/err_b(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'B in the plasma converges at order p')
        ! psi is measured on the discrete shell, whose inner side Gamma_p
        ! converges at order p, and Gamma_p in the maximum norm
        call check(log(err_psi(level - 1)/err_psi(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'psi in the shell converges at order p')
        call check(log(err_gamma(level - 1)/err_gamma(level))/log(2.0_dp) > degree - 0.5_dp, &
                   'the plasma boundary converges at order p')
        call check(log(err_axis(level - 1)/err_axis(level))/log(2.0_dp) > degree - 0.3_dp, &
                   'the axis converges at order p')
    end do
    if (failures > 0) error stop 'free-boundary label-lane gate E6a failed'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            print '(2a)', 'PASS ', what
        else
            print '(2a)', 'FAIL ', what
            failures = failures + 1
        end if
    end subroutine check

    ! ---------------------------------------------------------------- Gamma_p

    pure real(dp) function mode(k, theta)
        !! k = 1: 1; k = 2m: cos(m theta); k = 2m + 1: sin(m theta).
        integer, intent(in) :: k
        real(dp), intent(in) :: theta

        if (k == 1) then
            mode = 1
        else if (mod(k, 2) == 0) then
            mode = cos((k/2)*theta)
        else
            mode = sin((k/2)*theta)
        end if
    end function mode

    real(dp) function radius(rho, theta)
        real(dp), intent(in) :: rho(:), theta
        integer :: k

        radius = 0
        do k = 1, size(rho)
            radius = radius + rho(k)*mode(k, theta)
        end do
    end function radius

    function exact_rho() result(rho)
        !! Fourier coefficients of the exact Gamma_p (psi = psi_b) about (r0, 0).
        real(dp) :: rho(nrho), theta, r
        integer, parameter :: nq = 512
        integer :: j, k

        rho = 0
        do j = 1, nq
            theta = 2*pi*(j - 1)/nq
            r = ray_radius(psi_b, theta)
            do k = 1, nrho
                rho(k) = rho(k) + r*mode(k, theta)*merge(1.0_dp, 2.0_dp, k == 1)/nq
            end do
        end do
    end function exact_rho

    real(dp) function truncation_error() result(e)
        real(dp) :: rho(nrho), theta
        integer :: j

        rho = exact_rho()
        e = 0
        do j = 1, 720
            theta = 2*pi*(j - 1)/720
            e = max(e, abs(radius(rho, theta) - ray_radius(psi_b, theta)))
        end do
    end function truncation_error

    function curve_samples(rho) result(s)
        real(dp), intent(in) :: rho(:)
        real(dp) :: s(2, nsample), theta
        integer :: i

        do i = 1, nsample
            theta = 2*pi*(i - 1)/nsample
            s(:, i) = [r0, 0.0_dp] + radius(rho, theta)*[cos(theta), sin(theta)]
        end do
    end function curve_samples

    function wall_samples() result(s)
        !! Gamma_c: the level psi = 2 psi_b.
        real(dp) :: s(2, nsample), theta
        integer :: i

        do i = 1, nsample
            theta = 2*pi*(i - 1)/nsample
            s(:, i) = [r0, 0.0_dp] + ray_radius(2*psi_b, theta)*[cos(theta), sin(theta)]
        end do
    end function wall_samples

    ! ------------------------------------------------------------------ meshes

    subroutine build_meshes(n, x, rho)
        !! Plasma O-grid with wall Gamma_p and its core moved by x - (r0, 0)
        !! (E1); full O-grid with the annulus between Gamma_p and Gamma_c. The
        !! spaces are built at the first call for n, then only moved.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2), rho(:)
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad_p
        real(dp) :: sp(2, nsample), sw(2, nsample), rp, rw
        integer :: st

        sp = curve_samples(rho)
        sw = wall_samples()
        rp = 0.5_dp*(maxval(sp(1, :)) - minval(sp(1, :)))
        rw = 0.5_dp*(maxval(sw(1, :)) - minval(sw(1, :)))
        ! plasma
        call layout_ogrid([r0, 0.0_dp], rp, 0.5_dp, layout, st)
        if (st /= 0) error stop 'plasma layout'
        ! the core centred on (r0, 0) whatever Gamma_p is, then moved to X: X
        ! is always the central vertex of the core block (n even), so the
        ! polar cells never change and the residual is smooth in Gamma_p and X
        ! (with the default core, X on a cell edge switches cells under a
        ! boundary deformation)
        call layout_set_curve(layout, LAYOUT_WALL, sp, st, centre_core=.true.)
        if (st /= 0) error stop 'plasma boundary curve'
        layout%vertices(:, 1:4) = layout%vertices(:, 1:4) + spread(x - [r0, 0.0_dp], 2, 4)
        if (cells_built /= n) then
            options%degree = degree
            options%regularity = SPLINE_C0
            allocate (options%knots(2, layout%block_count))
            options%knots(:, :)%element_count = n
            call spline_space_build(layout, options, space_p, st)
            if (st /= 0) error stop 'plasma space'
            deallocate (options%knots)
        else
            call spline_space_update_geometry(space_p, layout, st)
            if (st /= 0) error stop 'plasma geometry'
        end if
        call spline_derham_build(space_p, derham_p, st, 0)
        if (st /= 0) error stop 'plasma derham'
        call spline_quadrature_build(space_p, degree + 8, quad_p, st)
        if (st /= 0) error stop 'plasma quadrature'
        call spline_quadrature_point_singular(space_p, quad_p, x, aquad_p, st)
        if (st /= 0) error stop 'axis quadrature'
        ! full: core and ruled blocks inside Gamma_p, annulus (block 6) to Gamma_c
        call layout_ogrid([r0, 0.0_dp], rw, 0.3_dp, layout, st, &
                          interface_fraction=0.6_dp)
        if (st /= 0) error stop 'full layout'
        call layout_set_curve(layout, LAYOUT_WALL, sw, st)
        if (st /= 0) error stop 'wall curve'
        call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, sp, st)
        if (st /= 0) error stop 'interface curve'
        if (cells_built /= n) then
            options%degree = degree
            options%regularity = SPLINE_C0
            allocate (options%knots(2, layout%block_count))
            options%knots(:, :)%element_count = n
            options%knots(1, 6)%element_count = 4*n
            call spline_space_build(layout, options, space_f, st)
            if (st /= 0) error stop 'full space'
        else
            call spline_space_update_geometry(space_f, layout, st)
            if (st /= 0) error stop 'full geometry'
        end if
        call spline_quadrature_build(space_f, degree + 4, quad_f, st)
        if (st /= 0) error stop 'full quadrature'
        if (cells_built /= n) call setup_shell()
        cells_built = n
    end subroutine build_meshes

    subroutine setup_shell()
        !! Dof sets, the exterior terms on the fixed Gamma_c and the coil data
        !! (the exact psi there); none of them moves with Gamma_p.
        logical, allocatable :: inside(:)
        integer :: k, j, b, st
        real(dp) :: gr(2)

        nfull = space_f%dof_count
        allocate (inside(nfull))
        inside = .false.
        do b = 1, 5
            associate (blk => space_f%blocks(b))
                do j = blk%offset + 1, blk%offset + blk%n(1)*blk%n(2)
                    do k = space_f%broken_ptr(j), space_f%broken_ptr(j + 1) - 1
                        if (space_f%broken_weight(k) /= 0) inside(space_f%broken_dof(k)) = .true.
                    end do
                end do
            end associate
        end do
        tied = pack([(k, k=1, nfull)], inside)
        free = pack([(k, k=1, nfull)], .not. inside)
        nfree = size(free)
        call gs_spline_vacuum_setup(space_f, LAYOUT_WALL, 2, 1.0_dp, vac, st)
        if (st /= 0) error stop 'exterior setup'
        if (any(inside(vac%dofs))) error stop 'Gamma_c touches the tied functions'
        if (allocated(psi_c)) deallocate (psi_c, dn_psi_c)
        allocate (psi_c(size(vac%nodes, 2)), dn_psi_c(vac%side%point_count))
        do k = 1, size(psi_c)
            psi_c(k) = psi(vac%nodes(1, k), vac%nodes(2, k))
        end do
        do k = 1, size(dn_psi_c)
            gr = grad_psi(vac%side%s(1, k), vac%side%s(2, k))
            dn_psi_c(k) = dot_product(gr, vac%normal(:, k))
        end do
    end subroutine setup_shell

    function full_psi(psi_free, psi_bv) result(v)
        real(dp), intent(in) :: psi_free(:), psi_bv
        real(dp) :: v(nfull)

        v(free) = psi_free
        v(tied) = psi_bv
    end function full_psi

    ! ---------------------------------------------------------------- residual

    subroutine residual(n, u, r, jac_t)
        !! The residual at u = (chi, g, X, mu, psi_free, psi_b, rho); with jac_t
        !! the transposed Newton matrix without the X and rho columns (they are
        !! added by the caller).
        integer, intent(in) :: n
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: r(:)
        type(triplet_list), intent(inout), optional :: jac_t
        real(dp), allocatable :: chi(:), g(:), grad(:), nrow(:), ngrad(:, :), psi_full(:), &
                                 gpsi(:, :), coef(:, :)
        real(dp) :: x(2), mu, energy, psi_bv
        type(csr_t) :: hessian
        type(csr_matrix_t) :: hv
        type(gs_spline_kernel_t) :: kernel
        integer :: ix, imu, ipsi, ipb, irho, k, i, p, st, j

        call unpack(u, chi, g, x, mu, psi_bv)
        call build_meshes(n, x, u(irho_of() + 1:))
        ix = 2*n0p
        imu = ix + 2
        ipsi = imu + 1
        ipb = ipsi + nfree
        irho = ipb + 1
        allocate (r(size(u)))
        r = 0
        ! label lane
        call clebsch_axisym_system(derham_p, aquad_p, solovev_profiles, x, chi, g, energy, &
                                   grad, hessian)
        call axis_jets(x, nrow, ngrad)
        r(1:2*n0p) = grad/(2*pi)
        r(1:n0p) = r(1:n0p) + mu*nrow
        r(ix + 1:ix + 2) = matmul(chi, ngrad)
        r(imu + 1) = dot_product(nrow, chi)
        ! shell
        psi_full = full_psi(u(ipsi + 1:ipb), psi_bv)
        allocate (kernel%r(quad_f%point_count), kernel%f(quad_f%point_count))
        kernel%r = quad_f%s(1, :)
        kernel%f = merge(-source*quad_f%s(1, :), 0.0_dp, quad_f%block == 6)
        allocate (coef(1, nfull), gpsi(1, nfull))
        coef(1, :) = psi_full
        call spline_energy(space_f, quad_f, kernel, coef, energy, st, gpsi)
        if (st /= 0) error stop 'shell energy'
        gpsi(1, vac%dofs) = gpsi(1, vac%dofs) + matmul(vac%block, psi_full(vac%dofs)) &
                            - matmul(transpose(vac%trace), matmul(vac%steklov, psi_c))
        do k = 1, vac%side%point_count
            do j = vac%side_rowptr(k), vac%side_rowptr(k + 1) - 1
                gpsi(1, vac%side_colidx(j)) = gpsi(1, vac%side_colidx(j)) - &
                    vac%side%weight(k)*dn_psi_c(k)/vac%side%s(1, k)*vac%side_val(j)
            end do
        end do
        r(ipsi + 1:ipb) = gpsi(1, free)
        ! Gamma_p rows: current and pressure balance
        call boundary_rows(chi, g, x, psi_full, r(ipb + 1), r(irho + 1:), jac_t, ipsi, ipb, irho)
        r(wall_p) = 0
        if (.not. present(jac_t)) return
        ! the label-lane block, the axis rows, the multiplier
        do i = 1, hessian%nrow
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                call jac_t%add(hessian%col(p), i, hessian%val(p)/(2*pi))
            end do
        end do
        do i = n0p + 1, 2*n0p
            call jac_t%add(i, i, newton_ridge*maxval(abs(hessian%val))/(2*pi))
        end do
        do i = 1, n0p
            call jac_t%add(i, ix + 1, ngrad(i, 1))
            call jac_t%add(i, ix + 2, ngrad(i, 2))
            call jac_t%add(i, imu + 1, nrow(i))
            call jac_t%add(imu + 1, i, nrow(i))
        end do
        ! the shell block: H(free, free), H(free, tied) 1, the exterior block
        call spline_energy_hessian_csr(space_f, quad_f, kernel, coef, hv, st)
        if (st /= 0) error stop 'shell Hessian'
        block
            integer, allocatable :: slot(:)
            logical, allocatable :: is_tied(:)
            allocate (slot(nfull), is_tied(nfull))
            slot = 0
            slot(free) = [(k, k=1, nfree)]
            is_tied = .false.
            is_tied(tied) = .true.
            do i = 1, hv%n
                if (slot(i) == 0) cycle
                do p = hv%rowptr(i), hv%rowptr(i + 1) - 1
                    j = hv%colidx(p)
                    if (slot(j) /= 0) then
                        call jac_t%add(ipsi + slot(j), ipsi + slot(i), hv%val(p))
                    else if (is_tied(j)) then
                        call jac_t%add(ipb + 1, ipsi + slot(i), hv%val(p))
                    end if
                end do
            end do
            do i = 1, size(vac%dofs)
                do k = 1, size(vac%dofs)
                    call jac_t%add(ipsi + slot(vac%dofs(k)), ipsi + slot(vac%dofs(i)), &
                                   vac%block(i, k))
                end do
            end do
        end block
        do i = ipsi + 1, ipb + 1
            call jac_t%add(i, i, 0.0_dp)
        end do
    end subroutine residual

    integer function irho_of() result(i)
        i = 2*n0p + 3 + nfree + 1
    end function irho_of

    subroutine unpack(u, chi, g, x, mu, psi_bv)
        real(dp), intent(in) :: u(:)
        real(dp), allocatable, intent(out) :: chi(:), g(:)
        real(dp), intent(out) :: x(2), mu, psi_bv

        chi = u(1:n0p)
        g = u(n0p + 1:2*n0p)
        x = u(2*n0p + 1:2*n0p + 2)
        mu = u(2*n0p + 3)
        psi_bv = u(2*n0p + 3 + nfree + 1)
    end subroutine unpack

    subroutine axis_jets(x, nrow, ngrad)
        !! Node function values at X and their X derivatives.
        real(dp), intent(in) :: x(2)
        real(dp), allocatable, intent(out) :: nrow(:), ngrad(:, :)
        real(dp), parameter :: h = 1e-5_dp

        nrow = values_at(x)
        allocate (ngrad(size(nrow), 2))
        ngrad(:, 1) = (values_at(x + [h, 0.0_dp]) - values_at(x - [h, 0.0_dp]))/(2*h)
        ngrad(:, 2) = (values_at(x + [0.0_dp, h]) - values_at(x - [0.0_dp, h]))/(2*h)
    end subroutine axis_jets

    function values_at(x) result(row)
        real(dp), intent(in) :: x(2)
        real(dp), allocatable :: row(:)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: val(:)
        real(dp) :: xi(2)
        integer :: block, st

        call spline_locate(space_p, 1, x, block, xi, st)
        if (st /= 0) error stop 'axis location'
        call spline_point_matrix(space_p, [block], reshape(xi, [2, 1]), 0, rp, ci, val, st)
        allocate (row(n0p))
        row = 0
        row(ci) = val
    end function values_at

    subroutine boundary_rows(chi, g, x, psi_full, r_current, r_balance, jac_t, ipsi, ipb, irho)
        !! The current row and the pressure-balance rows, each side on its own
        !! representation of Gamma_p; with jac_t their derivatives in chi, g,
        !! psi_free and psi_b by the chain rule.
        real(dp), intent(in) :: chi(:), g(:), x(2), psi_full(:)
        real(dp), intent(out) :: r_current, r_balance(:)
        type(triplet_list), intent(inout), optional :: jac_t
        integer, intent(in) :: ipsi, ipb, irho
        type(spline_side_quadrature_t) :: side
        integer, allocatable :: rp0(:), ci0(:), rp1(:), ci1(:), rp2(:), ci2(:), slot(:)
        real(dp), allocatable :: v0(:), v1(:), v2(:)
        real(dp) :: s(2), nrm(2), c, gc(2), gg(2), th(2), f(3), gp(3), pr(3), bphi, bpol(2), &
                    q, dq(5), wl, theta, rr, gpsi(2), ip
        integer :: k, j, m, st, col

        r_current = 0
        r_balance = 0
        ! plasma side: the wall of the plasma mesh
        call spline_side_quadrature(space_p, LAYOUT_WALL, degree + 3, side, st)
        if (st /= 0) error stop 'plasma side quadrature'
        call spline_point_matrix(space_p, side%block, side%xi, 0, rp0, ci0, v0, st)
        call spline_point_matrix(space_p, side%block, side%xi, 1, rp1, ci1, v1, st)
        call spline_point_matrix(space_p, side%block, side%xi, 2, rp2, ci2, v2, st)
        do k = 1, side%point_count
            s = side%s(:, k)
            nrm = [side%tangent(2, k), -side%tangent(1, k)]
            c = dot_product(v0(rp0(k):rp0(k + 1) - 1), chi(ci0(rp0(k):rp0(k + 1) - 1)))
            gc = [dot_product(v1(rp1(k):rp1(k + 1) - 1), chi(ci1(rp1(k):rp1(k + 1) - 1))), &
                  dot_product(v2(rp2(k):rp2(k + 1) - 1), chi(ci2(rp2(k):rp2(k + 1) - 1)))]
            gg = [dot_product(v1(rp1(k):rp1(k + 1) - 1), g(ci1(rp1(k):rp1(k + 1) - 1))), &
                  dot_product(v2(rp2(k):rp2(k + 1) - 1), g(ci2(rp2(k):rp2(k + 1) - 1)))]
            rr = (s(1) - x(1))**2 + (s(2) - x(2))**2
            th = [-(s(2) - x(2)), s(1) - x(1)]/rr
            call solovev_profiles(c, f, gp, pr)
            bphi = gg(2)*gc(1) - gg(1)*gc(2) - f(1)*(gc(1)*th(2) - gc(2)*th(1))
            bpol = gp(1)*[-gc(2), gc(1)]/s(1)
            q = pr(1) + (bphi**2 + sum(bpol**2))/2
            ! dq/d(c, c_R, c_Z, g_R, g_Z)
            dq(1) = pr(2) - bphi*f(2)*(gc(1)*th(2) - gc(2)*th(1)) + &
                    gp(2)*(bpol(1)*(-gc(2)) + bpol(2)*gc(1))/s(1)
            dq(2) = bphi*(gg(2) - f(1)*th(2)) + bpol(2)*gp(1)/s(1)
            dq(3) = bphi*(-gg(1) + f(1)*th(1)) - bpol(1)*gp(1)/s(1)
            dq(4) = -bphi*gc(2)
            dq(5) = bphi*gc(1)
            wl = side%weight(k)
            theta = atan2(s(2), s(1) - r0)
            ip = gp(1)*dot_product(gc, nrm)/s(1)
            r_current = r_current + wl*ip
            do m = 1, size(r_balance)
                r_balance(m) = r_balance(m) + wl*q*mode(m, theta)
            end do
            if (.not. present(jac_t)) cycle
            ! d(current)/d chi_j and d(balance_m)/d(chi_j, g_j)
            do j = rp0(k), rp0(k + 1) - 1
                col = ci0(j)
                call jac_t%add(col, ipb + 1, wl*gp(2)*dot_product(gc, nrm)/s(1)*v0(j))
                do m = 1, size(r_balance)
                    call jac_t%add(col, irho + m, wl*mode(m, theta)*dq(1)*v0(j))
                end do
            end do
            do j = rp1(k), rp1(k + 1) - 1
                col = ci1(j)
                call jac_t%add(col, ipb + 1, wl*gp(1)*nrm(1)/s(1)*v1(j))
                do m = 1, size(r_balance)
                    call jac_t%add(col, irho + m, wl*mode(m, theta)*dq(2)*v1(j))
                    call jac_t%add(n0p + col, irho + m, wl*mode(m, theta)*dq(4)*v1(j))
                end do
            end do
            do j = rp2(k), rp2(k + 1) - 1
                col = ci2(j)
                call jac_t%add(col, ipb + 1, wl*gp(1)*nrm(2)/s(1)*v2(j))
                do m = 1, size(r_balance)
                    call jac_t%add(col, irho + m, wl*mode(m, theta)*dq(3)*v2(j))
                    call jac_t%add(n0p + col, irho + m, wl*mode(m, theta)*dq(5)*v2(j))
                end do
            end do
        end do
        ! shell side: the annulus (block 6) at the interface
        ! (spline_interface_quadrature finds no interface next to the periodic
        ! annulus; the side quadrature of the tag holds both sides, of which
        ! the points of block 6 are the shell side)
        call spline_side_quadrature(space_f, LAYOUT_ANNULUS_INTERFACE, degree + 3, side, st)
        if (st /= 0) error stop 'interface quadrature'
        if (count(side%block == 6) == 0) error stop 'interface quadrature: no shell points'
        call spline_point_matrix(space_f, side%block, side%xi, 1, rp1, ci1, v1, st)
        call spline_point_matrix(space_f, side%block, side%xi, 2, rp2, ci2, v2, st)
        allocate (slot(nfull))
        slot = 0
        slot(free) = [(k, k=1, nfree)]
        do k = 1, side%point_count
            if (side%block(k) /= 6) cycle
            s = side%s(:, k)
            nrm = [side%tangent(2, k), -side%tangent(1, k)]
            if (dot_product(nrm, s - x) < 0) nrm = -nrm
            gpsi = [dot_product(v1(rp1(k):rp1(k + 1) - 1), psi_full(ci1(rp1(k):rp1(k + 1) - 1))), &
                    dot_product(v2(rp2(k):rp2(k + 1) - 1), psi_full(ci2(rp2(k):rp2(k + 1) - 1)))]
            wl = side%weight(k)
            theta = atan2(s(2), s(1) - r0)
            q = (sum(gpsi**2) + f0**2)/(2*s(1)**2)
            r_current = r_current - wl*dot_product(gpsi, nrm)/s(1)
            do m = 1, size(r_balance)
                r_balance(m) = r_balance(m) - wl*q*mode(m, theta)
            end do
            if (.not. present(jac_t)) cycle
            do j = rp1(k), rp1(k + 1) - 1
                call shell_entry(jac_t, ci1(j), slot, ipsi, ipb, irho, size(r_balance), theta, &
                                 -wl*nrm(1)/s(1)*v1(j), -wl*gpsi(1)/s(1)**2*v1(j))
            end do
            do j = rp2(k), rp2(k + 1) - 1
                call shell_entry(jac_t, ci2(j), slot, ipsi, ipb, irho, size(r_balance), theta, &
                                 -wl*nrm(2)/s(1)*v2(j), -wl*gpsi(2)/s(1)**2*v2(j))
            end do
        end do
    end subroutine boundary_rows

    subroutine shell_entry(jac_t, col, slot, ipsi, ipb, irho, nbalance, theta, dcurrent, &
                           dbalance)
        !! d/d psi_col of the current row and of the balance rows (dbalance
        !! times mode m), routed to psi_free or, for a tied function, psi_b.
        type(triplet_list), intent(inout) :: jac_t
        integer, intent(in) :: col, slot(:), ipsi, ipb, irho, nbalance
        real(dp), intent(in) :: theta, dcurrent, dbalance
        integer :: target, m

        if (slot(col) /= 0) then
            target = ipsi + slot(col)
        else
            target = ipb + 1
        end if
        call jac_t%add(target, ipb + 1, dcurrent)
        do m = 1, nbalance
            call jac_t%add(target, irho + m, dbalance*mode(m, theta))
        end do
    end subroutine shell_entry

    ! ------------------------------------------------------------------- solve

    subroutine solve(n, err_b, err_psi, err_gamma, err_axis, steps)
        !! Perturbed start, then Newton in two stages: Gamma_p frozen (the
        !! fixed-boundary problem of E1, with the axis free, which makes labels,
        !! g, the axis and the shell consistent with the start), then
        !! everything free. (With the axis frozen too, chi(X) = 0 holds away
        !! from the label minimum and stage 1 stalls on fine meshes.)
        integer, intent(in) :: n
        real(dp), intent(out) :: err_b, err_psi, err_gamma, err_axis
        integer, intent(out) :: steps
        real(dp), allocatable :: u(:), r(:), chi(:), g(:), grad(:), work(:)
        real(dp) :: x(2), rho(nrho), energy, eb, ep, eg, ea
        type(csr_t) :: hessian
        integer :: st, k, nu, irho, stage1
        integer, allocatable :: wall2(:)

        rho = exact_rho()
        rho(1) = rho(1)*(1 + perturbation)
        rho(2) = rho(2) + perturbation*rho(1)
        rho(3) = rho(3) + 0.5_dp*perturbation*rho(1)
        x = [r0, 0.0_dp] + perturbation*[0.6_dp, 0.3_dp]
        cells_built = 0
        n0p = 0
        call build_meshes(n, x, rho)
        n0p = derham_p%dimension(0)
        call spline_derham_boundary_dofs(derham_p, LAYOUT_WALL, 0, wall2)
        wall_p = wall2
        nu = 2*n0p + 3 + nfree + 1 + nrho
        irho = irho_of()
        allocate (u(nu), chi(n0p), g(n0p))
        u = 0
        call spline_derham_project(derham_p, 0, label_cb, chi, st, 12)
        if (st /= 0) error stop 'label projection'
        chi = chi*(1 + 1.5_dp*perturbation*(1 - chi))
        chi(wall_p) = 1
        ! g at fixed labels
        g = 0
        call clebsch_axisym_system(derham_p, aquad_p, solovev_profiles, x, chi, g, energy, &
                                   grad, hessian)
        call add_ridge(hessian, n0p)
        allocate (work(2*n0p))
        call eu_constrained_solve(hessian, -grad, [(k, k=1, n0p)], 0*chi, work, st)
        if (st /= 0) error stop 'g solve'
        g = work(n0p + 1:)
        u(1:n0p) = chi
        u(n0p + 1:2*n0p) = g
        u(2*n0p + 1:2*n0p + 2) = x
        u(irho) = psi_b*(1 + perturbation)
        u(irho + 1:) = rho
        if (allocated(coarse_u)) then
            ! mesh continuation: the previous level's solution on this mesh
            call continuation_start(n, u)
            call errors(u, eb, ep, eg, ea)
            print '(a,4es10.2)', '  start from the coarse level: errors B, psi, Gamma_p, axis', &
                eb, ep, eg, ea
        else
            call errors(u, eb, ep, eg, ea)
            print '(a,4es10.2)', '  start: errors B, psi (shell), Gamma_p, axis', eb, ep, eg, ea
            call newton(n, u, [wall_p, [(irho + k, k=1, nrho)]], stage1, .false.)
            call errors(u, eb, ep, eg, ea)
            print '(a,i3,a,4es10.2)', '  stage 1 (Gamma_p frozen):', stage1, &
                ' steps; errors B, psi, Gamma_p, axis', eb, ep, eg, ea
        end if
        call newton(n, u, wall_p, steps, jacobian_check)
        call residual(n, u, r)
        call errors(u, err_b, err_psi, err_gamma, err_axis)
        ! keep this level for the start of the next
        coarse_p = space_p
        coarse_f = space_f
        coarse_u = u
        coarse_psi = full_psi(u(2*n0p + 4:2*n0p + 3 + nfree), u(irho))
        coarse_n0p = n0p
        coarse_nfree = nfree
        coarse_x = u(2*n0p + 1:2*n0p + 2)
    end subroutine solve

    subroutine continuation_start(n, u)
        !! Axis, boundary, psi_b and mu of the coarse level; chi and g projected
        !! from the coarse plasma mesh, psi from the coarse shell mesh (at the
        !! fine quadrature points, the tied functions at psi_b).
        integer, intent(in) :: n
        real(dp), intent(inout) :: u(:)
        real(dp), allocatable :: chi(:), g(:), vals(:), pc(:, :)
        integer :: st, k, irho, cirho

        irho = irho_of()
        cirho = 2*coarse_n0p + 4 + coarse_nfree
        u(2*n0p + 1:2*n0p + 3) = coarse_u(2*coarse_n0p + 1:2*coarse_n0p + 3)
        u(irho:) = coarse_u(cirho:)
        call build_meshes(n, u(2*n0p + 1:2*n0p + 2), u(irho + 1:))
        allocate (chi(n0p), g(n0p))
        transfer_field = 1
        call spline_derham_project(derham_p, 0, transfer_cb, chi, st, 12)
        if (st /= 0) error stop 'label transfer'
        transfer_field = 2
        call spline_derham_project(derham_p, 0, transfer_cb, g, st, 12)
        if (st /= 0) error stop 'g transfer'
        chi(wall_p) = 1
        u(1:n0p) = chi
        u(n0p + 1:2*n0p) = g
        allocate (vals(quad_f%point_count), pc(1, nfull))
        do k = 1, quad_f%point_count
            if (quad_f%block(k) == 6) then
                vals(k) = coarse_value(coarse_f, coarse_psi, quad_f%s(:, k), coarse_x)
            else
                vals(k) = u(irho)
            end if
        end do
        call gs_spline_project(space_f, quad_f, vals, pc, st)
        if (st /= 0) error stop 'psi transfer'
        u(2*n0p + 4:2*n0p + 3 + nfree) = pc(1, free)
    end subroutine continuation_start

    subroutine transfer_cb(block, xi, t, values)
        !! The coarse chi (transfer_field 1) or g (2) at the fine point.
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(space_p, block, xi, jet, st)
        if (transfer_field == 1) then
            values(1) = coarse_value(coarse_p, coarse_u(1:coarse_n0p), jet(:, 1), coarse_x)
        else
            values(1) = coarse_value(coarse_p, coarse_u(coarse_n0p + 1:2*coarse_n0p), &
                                     jet(:, 1), coarse_x)
        end if
        associate (unused => t)
        end associate
    end subroutine transfer_cb

    real(dp) function coarse_value(space, coef, s, centre) result(v)
        !! A field of the coarse space at the physical point s; a point just
        !! outside the coarse boundary (the two meshes represent Gamma_p to
        !! O(h^3)) is moved towards centre until it is located.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:), s(2), centre(2)
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: val(:)
        real(dp) :: xi(2), y(2)
        integer :: block, st, k

        y = s
        do k = 1, 30
            call spline_locate(space, 1, y, block, xi, st)
            if (st == 0) exit
            y = centre + (1 - 1e-6_dp*2**k)*(s - centre)
        end do
        if (st /= 0) error stop 'coarse location'
        call spline_point_matrix(space, [block], reshape(xi, [2, 1]), 0, rp, ci, val, st)
        v = dot_product(val, coef(ci))
    end function coarse_value

    subroutine newton(n, u, fixed, steps, check)
        !! Newton with the unknowns `fixed` held. The X and rho columns that are
        !! free are central differences with the meshes rebuilt. Monotone
        !! backtracking within a trust region (Gamma_p by at most a tenth of
        !! its mean radius, the axis by at most 0.1).
        integer, intent(in) :: n, fixed(:)
        real(dp), intent(inout) :: u(:)
        integer, intent(out) :: steps
        logical, intent(in) :: check
        real(dp), allocatable :: r(:), rp(:), rm(:), delta(:), history(:), up(:), rt(:)
        real(dp) :: alpha_step, h
        type(triplet_list) :: jt
        type(csr_t) :: jac_t
        logical, allocatable :: is_fixed(:)
        integer :: st, k, it, j, nu, irho, i

        nu = size(u)
        irho = irho_of()
        allocate (is_fixed(nu))
        is_fixed = .false.
        is_fixed(fixed) = .true.
        allocate (history(0))
        steps = 0
        do it = 1, max_steps
            call reset(jt)
            call jt%init(64*nu)
            call residual(n, u, r, jt)
            r(fixed) = 0
            history = [history, norm2(r)]
            ! just above the floor of the difference columns
            if (history(it) < 5e-8_dp) exit
            if (it > 1) then
                if (history(it) < 1e-11_dp*history(1)) exit
            end if
            do j = 1, 2 + nrho
                if (j <= 2) then
                    k = 2*n0p + j
                else
                    k = irho + j - 2
                end if
                if (is_fixed(k)) then
                    call jt%add(k, k, 1.0_dp)
                    cycle
                end if
                h = 1e-6_dp
                up = u
                up(k) = u(k) + h
                call residual(n, up, rp)
                up(k) = u(k) - h
                call residual(n, up, rm)
                rp = (rp - rm)/(2*h)
                do i = 1, nu
                    if (rp(i) /= 0) call jt%add(k, i, rp(i))
                end do
                call jt%add(k, k, 0.0_dp)
            end do
            call jt%to_csr(nu, nu, jac_t)
            if (check .and. it == 1) call check_jacobian(n, u, r, jac_t)
            allocate (delta(nu))
            call eu_constrained_solve(jac_t, -r, fixed, 0*u(fixed), delta, st)
            if (st /= 0) error stop 'Newton solve'
            print '(a,7es10.2)', '    step |d| chi, g, X, mu, psi, psi_b, rho', &
                maxval(abs(delta(1:n0p))), maxval(abs(delta(n0p + 1:2*n0p))), &
                maxval(abs(delta(2*n0p + 1:2*n0p + 2))), abs(delta(2*n0p + 3)), &
                maxval(abs(delta(2*n0p + 4:irho - 1))), abs(delta(irho)), &
                maxval(abs(delta(irho + 1:)))
            alpha_step = min(1.0_dp, &
                             0.1_dp*u(irho + 1)/max(maxval(abs(delta(irho + 1:))), tiny(1.0_dp)), &
                             0.1_dp/max(maxval(abs(delta(2*n0p + 1:2*n0p + 2))), tiny(1.0_dp)))
            do k = 1, 12
                call residual(n, u + alpha_step*delta, rt)
                rt(fixed) = 0
                if (norm2(rt) < (1 - alpha_step/4)*history(it)) exit
                alpha_step = alpha_step/2
            end do
            print '(a,i3,a,es10.3,a,f6.3,a,4es10.2)', '  step', it, '  |r|', history(it), &
                '  alpha', alpha_step, '  |r| plasma, shell, current, balance', &
                norm2(r(1:2*n0p + 3)), norm2(r(2*n0p + 4:irho - 1)), abs(r(irho)), &
                norm2(r(irho + 1:))
            u = u + alpha_step*delta
            deallocate (delta)
            steps = steps + 1
        end do
        print '(a,20es9.1)', '  |r| history', history
    end subroutine newton

    subroutine check_jacobian(n, u, r, jac_t)
        !! J v against central differences of the residual per block of v,
        !! the plasma wall rows excluded.
        integer, intent(in) :: n
        real(dp), intent(in) :: u(:), r(:)
        type(csr_t), intent(in) :: jac_t
        real(dp), allocatable :: v(:), jv(:), rp(:), rm(:), fd(:)
        integer :: blk, lo(6), hi(6), i, p
        character(8), parameter :: names(6) = ['chi     ', 'g       ', 'X, mu   ', &
                                               'psi     ', 'psi_b   ', 'rho     ']
        real(dp) :: hv

        lo = [1, n0p + 1, 2*n0p + 1, 2*n0p + 4, 2*n0p + 4 + nfree, 2*n0p + 5 + nfree]
        hi = [n0p, 2*n0p, 2*n0p + 3, 2*n0p + 3 + nfree, 2*n0p + 4 + nfree, size(u)]
        allocate (v(size(u)), jv(size(u)))
        do blk = 1, 6
            v = 0
            do i = lo(blk), hi(blk)
                v(i) = sin(1.3_dp*i) + 0.5_dp
            end do
            v(wall_p) = 0
            hv = 1e-6_dp/maxval(abs(v))
            jv = 0
            do i = 1, jac_t%nrow
                do p = jac_t%row_ptr(i), jac_t%row_ptr(i + 1) - 1
                    jv(jac_t%col(p)) = jv(jac_t%col(p)) + jac_t%val(p)*v(i)
                end do
            end do
            call residual(n, u + hv*v, rp)
            call residual(n, u - hv*v, rm)
            fd = (rp - rm)/(2*hv)
            jv(wall_p) = 0
            print '(a,a8,a,es10.2,a,es10.2)', '  Jacobian check, v in ', names(blk), &
                ': |Jv - FD|/|FD|', norm2(jv - fd)/max(norm2(fd), tiny(1.0_dp)), '  |FD|', &
                norm2(fd)
        end do
        ! single rho modes at three steps: smoothness of the residual in rho
        block
            real(dp) :: hs(4) = [1e-4_dp, 1e-5_dp, 1e-6_dp, 1e-7_dp]
            real(dp), allocatable :: d(:, :)
            integer :: mk, ih, kk(3)
            kk = [1, 2, 2*nmode + 1]
            do mk = 1, 3
                allocate (d(size(u), 4))
                do ih = 1, 4
                    v = 0
                    v(lo(6) + kk(mk) - 1) = 1
                    call residual(n, u + hs(ih)*v, rp)
                    call residual(n, u - hs(ih)*v, rm)
                    d(:, ih) = (rp - rm)/(2*hs(ih))
                    d(wall_p, ih) = 0
                end do
                print '(a,i3,a,es10.2,a,3es10.2)', '  rho mode', kk(mk), ': |FD|', norm2(d(:, 3)), &
                    '  |FD(h) - FD(1e-6)|/|FD| for h = 1e-4, 1e-5, 1e-7', &
                    norm2(d(:, 1) - d(:, 3))/norm2(d(:, 3)), norm2(d(:, 2) - d(:, 3))/norm2(d(:, 3)), &
                    norm2(d(:, 4) - d(:, 3))/norm2(d(:, 3))
                ! the largest entries and where
                print '(a,3i6,a,3es10.2)', '    largest rows', maxloc(abs(d(:, 3))), 0, 0, &
                    '  value', maxval(abs(d(:, 3))), 0.0_dp, 0.0_dp
                deallocate (d)
            end do
        end block
        associate (unused => r)
        end associate
    end subroutine check_jacobian

    subroutine errors(u, err_b, err_psi, err_gamma, err_axis)
        !! Relative L2 errors of B in the plasma and of psi in the shell, the
        !! largest radial error of Gamma_p and the axis error.
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: err_b, err_psi, err_gamma, err_axis
        real(dp), allocatable :: chi(:), g(:), v0(:, :, :), gn(:, :, :), psi_full(:)
        integer, allocatable :: d0(:), rp(:), ci(:)
        real(dp), allocatable :: val(:), pv(:)
        real(dp) :: x(2), mu, psi_bv, rho(nrho), s(2), w, c, gc(2), gg(2), th(2), f(3), gp(3), &
                    pr(3), bphi, bpol(2), gr(2), num, den, theta, jet(2, 6)
        integer :: cell, q, p, j, st

        call unpack(u, chi, g, x, mu, psi_bv)
        rho = u(irho_of() + 1:)
        call build_meshes(cells_built, x, rho)
        num = 0
        den = 0
        do cell = 1, size(aquad_p%point_offset)
            call clebsch_node_gradients(derham_p, aquad_p, cell, v0, d0, gn)
            do q = 1, aquad_p%points_per_cell
                p = aquad_p%point_offset(cell) + q - 1
                call spline_geometry_at(space_p, aquad_p%block(p), aquad_p%xi(:, p), jet, st)
                s = jet(:, 1)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                gg = matmul(gn(:, :, q), g(d0))
                th = [-(s(2) - x(2)), s(1) - x(1)]/sum((s - x)**2)
                call solovev_profiles(c, f, gp, pr)
                bphi = gg(2)*gc(1) - gg(1)*gc(2) - f(1)*(gc(1)*th(2) - gc(2)*th(1))
                bpol = gp(1)*[-gc(2), gc(1)]/s(1)
                gr = grad_psi(s(1), s(2))
                w = 2*pi*s(1)*aquad_p%weight(p)
                num = num + w*((bphi - f0/s(1))**2 + (bpol(1) + gr(2)/s(1))**2 + &
                               (bpol(2) - gr(1)/s(1))**2)
                den = den + w*((f0/s(1))**2 + sum(gr**2)/s(1)**2)
            end do
        end do
        err_b = sqrt(num/den)
        psi_full = full_psi(u(2*n0p + 4:2*n0p + 3 + nfree), psi_bv)
        call spline_point_matrix(space_f, quad_f%block, quad_f%xi, 0, rp, ci, val, st)
        num = 0
        den = 0
        do q = 1, quad_f%point_count
            if (quad_f%block(q) /= 6) cycle
            pv = val(rp(q):rp(q + 1) - 1)
            c = dot_product(pv, psi_full(ci(rp(q):rp(q + 1) - 1)))
            s = quad_f%s(:, q)
            num = num + quad_f%weight(q)*(c - psi(s(1), s(2)))**2
            den = den + quad_f%weight(q)*psi(s(1), s(2))**2
        end do
        err_psi = sqrt(num/den)
        err_gamma = 0
        do j = 1, 720
            theta = 2*pi*(j - 1)/720
            err_gamma = max(err_gamma, abs(radius(rho, theta) - ray_radius(psi_b, theta)))
        end do
        err_axis = norm2(x - [r0, 0.0_dp])
        print '(a,es12.4,a,es10.2)', '  psi_b', psi_bv, '  (exact 6.25), mu', mu
    end subroutine errors

    subroutine label_cb(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6)
        integer :: st

        call spline_geometry_at(space_p, block, xi, jet, st)
        values(1) = psi(jet(1, 1), jet(2, 1))/psi_b
        associate (unused => t)
        end associate
    end subroutine label_cb

    subroutine reset(t)
        type(triplet_list), intent(inout) :: t

        if (allocated(t%rows)) deallocate (t%rows, t%cols, t%vals)
        t%count = 0
    end subroutine reset

    subroutine add_ridge(hessian, n0)
        type(csr_t), intent(inout) :: hessian
        integer, intent(in) :: n0
        real(dp) :: shift
        integer :: i, p

        shift = ridge*maxval(abs(hessian%val))
        do i = n0 + 1, 2*n0
            do p = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                if (hessian%col(p) == i) hessian%val(p) = hessian%val(p) + shift
            end do
        end do
    end subroutine add_ridge

end program test_eulerian_clebsch_free_boundary
