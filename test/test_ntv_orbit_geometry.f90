!> Finite-orbit geometry of the orbit-class spectral elements
!! (milca_ntv_orbit_geometry) in its thin-orbit limit. On a synthetic Boozer
!! file with B = 2.3 (1 - 0.1 sqrt(s) cos(theta)), constant G, I = 0 and no
!! magnetic shear (with shear the orbit also precesses through -q' <theta
!! psi_dot>, which the local thin field lacks) the
!! orbit has theta_dot = v_par iota B/G + H_psi and zeta_dot = v_par B/G, so in
!! the frame exp(-i n zeta) it precesses at n q H_psi with
!!   H_psi = (m v^2/(2 Z e B)) (1 + xi^2) dB/dpsi.
!! The thin field (hp_field_t) with Boozer streaming, the full drift weight
!! and the precession prefactor
!!   pre = q m v^2 dB/ds(cos)/(2 Z e B_ref psi'),
!! (offset 0) solves the same problem to first order in the orbit width, so
!! the NTV kernels of the two geometries agree when the orbit is thin, for
!! either sign of the charge.
!!
!! The energy_* checks run the energy-coupled kernel (milca_dk_hp_energy)
!! with every speed node on its own orbit geometry and mesh against the same
!! kernel on the thin field, the per-node pitch-angle kernel, its own
!! conservation laws and a manufactured solution:
!!   test_ntv_orbit_geometry [factor]       thin limit of one kernel
!!   test_ntv_orbit_geometry thin_negative  the same for charge -e
!!   test_ntv_orbit_geometry energy_<name>  width, lorentz, conservation,
!!                                          manufactured, electrons

program test_ntv_orbit_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_state
    use milca_dk_hp, only: hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, hp_operator_t, &
        hp_ntv_kernel, HP_OK
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t, hp_orbit_geometry, &
        hp_orbit_perturbation
    use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, hp_energy_coefficients, &
        hp_energy_collide, hp_energy_apply, hp_energy_solve, hp_energy_free, hp_node_geometry_t
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_nuhat, collop_vth, COLLOP_OK, COLLOP_LORENTZ, COLLOP_TEST_PARTICLE, COLLOP_FULL
    implicit none
    !> The energy-coupled kernel of one surface with its nodes: per species a
    !> and speed node k the geometry (orbit geometry, or the thin field f),
    !> mesh meshes(mid(k, a)) and operator.
    type :: coupled_t
        type(hp_orbit_geometry_t), allocatable :: go(:, :)
        type(hp_field_t) :: f
        type(hp_mesh_t), allocatable :: meshes(:)
        type(hp_operator_t), allocatable :: ops(:, :)
        integer, allocatable :: mid(:, :)
        type(collop_species_t), allocatable :: sp(:)
        type(collop_grid_t) :: grid
        real(dp), allocatable :: lnl(:, :)
        logical :: thin = .true.
        type(hp_energy_system_t) :: sys
    end type
    real(dp), parameter :: pi = 4*atan(1.0_dp), e = 1.602176634e-19_dp
    real(dp), parameter :: md = 3.3435837768e-27_dp, me = 9.1093837015e-31_dp, &
        s_star = 0.45_dp, b_ref = 2.3_dp, kev = 1.602176634e-16_dp
    character(*), parameter :: file = 'ntv_orbit_geometry.bc'
    type(boozer_file_t), target :: bf
    character(32) :: check
    real(dp) :: factor
    integer :: status
    logical :: ok

    call write_file(file)
    call load_boozer_file(file, bf, status)
    if (status /= BOOZER_OK) error stop 'cannot read the test file'
    check = 'thin'
    factor = 1e-3_dp
    if (command_argument_count() > 0) then
        call get_command_argument(1, check)
        read (check, *, iostat=status) factor
        if (status == 0) check = 'thin'
        if (status /= 0) factor = 1e-3_dp
    end if
    select case (trim(check))
    case ('thin'); ok = thin(factor, 1.0_dp)
    case ('thin_negative'); ok = thin(1e-3_dp, -1.0_dp)
    case ('energy_width'); ok = energy_width()
    case ('energy_mixed'); ok = energy_mixed()
    case ('energy_lorentz'); ok = energy_lorentz()
    case ('energy_conservation'); ok = energy_conservation()
    case ('energy_manufactured'); ok = energy_manufactured()
    case ('energy_electrons'); ok = energy_electrons()
    case default; error stop 'unknown check'
    end select
    if (.not. ok) error stop 1

contains

    !> The monoenergetic kernel on the orbit geometry against the thin field
    !> at mass factor factor (header).
    logical function thin(factor, z) result(ok)
        real(dp), intent(in) :: factor, z
        type(boozer_file_t), target, save :: pert
        type(orbit_particle_t) :: p
        type(hp_orbit_geometry_t) :: g
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh_g, mesh_f
        type(hp_operator_t) :: op_g, op_f
        real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, v, r_g, r_f, samples(33), &
            omega_d
        integer :: status, j
        integer(8) :: t0, t1, rate

        p%mass = md*factor
        p%charge = z*e
        v = sqrt(2*10e3_dp*e/md)
        p%energy = p%mass*v**2/2
        p%s_ref = s_star
        call hp_orbit_geometry(bf, p, s_star, b_ref, g, status)
        if (status /= 0) error stop 'orbit geometry failed'
        call orbit_state(bf, s_star, 0.0_dp, b, b_s, gg, i, g_s, i_s, iota, psi_prime, status)
        do j = 1, 33
            call orbit_state(bf, s_star, 2*pi*(j - 1)/33, b, b_s, gg, i, g_s, i_s, iota, &
                psi_prime, status)
            samples(j) = b/b_ref
        end do
        call hp_field(samples, iota*b_ref/gg, f, status)
        if (status /= HP_OK) error stop 'hp_field failed'

        ! dB/ds = -0.23/(2 sqrt(s)) cos(theta) for this file.
        op_f%pre = (1/iota)*p%mass*v**2*(-0.23_dp/(2*sqrt(s_star)))/(2*z*e*b_ref*psi_prime)
        omega_d = abs(op_f%pre)/v
        op_f%k = iota*b_ref/gg
        op_f%k_power = 1
        op_f%q = 1/iota
        op_f%n = 3
        op_f%m = [3]
        op_f%c = [(1e-3_dp, 0.0_dp)]
        op_f%v = v
        op_f%wsign = 1
        op_f%offset = 0
        op_f%nu = 0.1_dp*op_f%n*omega_d
        op_g = op_f
        op_g%pre = 0

        call hp_default_mesh(f, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_f)
        call hp_ntv_kernel(f, mesh_f, op_f, r_f, status)
        if (status /= HP_OK) error stop 'thin kernel failed'
        call system_clock(t0, rate)
        call hp_default_mesh(g, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_g)
        call hp_ntv_kernel(g, mesh_g, op_g, r_g, status)
        call system_clock(t1)
        if (status /= HP_OK) error stop 'orbit kernel failed'
        print '(a,2es20.11,a,es9.2,a,f8.2)', 'R thin field, orbit geometry ', r_f, r_g, ' rel ', &
            abs(r_g/r_f - 1), ' orbit s ', real(t1 - t0, dp)/rate
        ok = abs(r_g/r_f - 1) <= 1e-3_dp

        ! The drive taken along the orbit from a perturbation file whose m = 3
        ! amplitude at s* equals the operator's: the same kernel for thin orbits.
        call write_pert('ntv_orbit_pert.bc')
        call load_boozer_file('ntv_orbit_pert.bc', pert, status)
        if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
        call hp_orbit_perturbation(g, pert, 3, b_ref)
        op_g%c = [cmplx(1e-3_dp*b_ref*s_star**1.5_dp, 0.0_dp, dp)/b_ref]
        op_g%m = [-3]
        call hp_ntv_kernel(g, mesh_g, op_g, r_f, status)
        if (status /= HP_OK) error stop 'orbit kernel with perturbation file failed'
        g%pert => null()
        call hp_ntv_kernel(g, mesh_g, op_g, r_g, status)
        print '(a,2es20.11,a,es9.2)', 'R drive along orbit, drive at s* ', r_f, r_g, ' rel ', &
            abs(r_f/r_g - 1)
        ok = ok .and. abs(r_f/r_g - 1) <= 1e-4_dp
    end function thin


    !> The coupled kernel of species with base masses mass and charge
    !> numbers z, all at 10 keV and 3e19 m^-3 at the base mass, the mass and
    !> temperature scaled by factor (fixed thermal speed, orbit width linear
    !> in factor). Speed grid of 3 elements of order 4 on [0, 4]. Orbit
    !> geometry and mesh per node, or with thin the thin field with the
    !> matching precession and one mesh. nu_D/v at x = 1 is a tenth of n
    !> times the drift frequency of species 1 there, at every factor.
    subroutine build(mass, z, factor, thin, model, c, status, l_field, apart, slow)
        real(dp), intent(in) :: mass(:), z(:), factor
        logical, intent(in) :: thin
        integer, intent(in) :: model
        type(coupled_t), intent(out) :: c
        integer, intent(out) :: status
        integer, intent(in), optional :: l_field
        ! Thin nodes that name distinct geometries: the shared points of the
        ! finite-width coupling instead of the common mass.
        logical, intent(in), optional :: apart
        ! With finite width: the slowest nodes of each species on the thin
        ! field, each on its own mesh.
        integer, intent(in), optional :: slow
        type(orbit_particle_t) :: p
        type(hp_node_geometry_t), allocatable :: nodes(:, :)
        integer :: n_slow
        real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, v, samples(33), kappa, pre1
        integer, allocatable :: gid(:, :)
        type(hp_field_t), allocatable :: geom(:, :)
        integer :: ns, nx, a, k, j

        ns = size(mass)
        c%thin = thin
        n_slow = 0
        if (present(slow)) n_slow = slow
        call collop_make_grid(4.0_dp, 3, 4, c%grid, status)
        if (status /= COLLOP_OK) error stop 'collop_make_grid failed'
        nx = c%grid%n - 1
        allocate (c%sp(ns), c%lnl(ns, ns), c%ops(nx, ns), c%mid(nx, ns), gid(nx, ns))
        do a = 1, ns
            c%sp(a) = collop_species_t(mass(a)*factor, z(a), 3e19_dp, 10*kev*factor)
        end do
        do j = 1, 33
            call orbit_state(bf, s_star, 2*pi*(j - 1)/33, b, b_s, gg, i, g_s, i_s, iota, &
                psi_prime, status)
            samples(j) = b/b_ref
        end do
        call hp_field(samples, iota*b_ref/gg, c%f, status)
        if (status /= HP_OK) error stop 'hp_field failed'
        ! dB/ds = -0.23/(2 sqrt(s)) cos(theta) for this file.
        pre1 = (1/iota)*c%sp(1)%m*collop_vth(c%sp(1))**2*(-0.23_dp/(2*sqrt(s_star))) &
            /(2*z(1)*e*b_ref*psi_prime)
        kappa = 0.1_dp*3*abs(pre1)/collop_vth(c%sp(1))
        c%lnl = 3*sqrt(pi)/8*collop_vth(c%sp(1))*kappa/collop_nuhat(c%sp(1), c%sp(1), 1.0_dp)
        if (thin) then
            allocate (c%meshes(1))
            call hp_default_mesh(c%f, 4, 16, 6, 1e-3_dp, 0.3_dp, c%meshes(1))
        else
            allocate (c%go(nx, ns), c%meshes(nx*ns))
        end if
        do a = 1, ns
            do k = 1, nx
                v = collop_vth(c%sp(a))*c%grid%x(k + 1)
                associate (op => c%ops(k, a))
                    op%k = iota*b_ref/gg
                    op%k_power = 1
                    op%q = 1/iota
                    op%n = 3
                    op%m = [3]
                    op%c = [(1e-3_dp, 0.0_dp)]
                    op%v = v
                    op%wsign = 1
                    op%offset = 0
                    op%pre = 0
                    if (thin .or. k <= n_slow) op%pre = (1/iota)*c%sp(a)%m*v**2*(-0.23_dp/(2*sqrt(s_star))) &
                        /(2*z(a)*e*b_ref*psi_prime)
                end associate
                if (thin) then
                    c%mid(k, a) = 1
                    gid(k, a) = 1
                    if (present(apart)) then
                        if (apart) gid(k, a) = (a - 1)*nx + k
                    end if
                    cycle
                end if
                c%mid(k, a) = (a - 1)*nx + k
                gid(k, a) = c%mid(k, a)
                if (k <= n_slow) then
                    call hp_default_mesh(c%f, 4, 16, 6, 1e-3_dp, 0.3_dp, c%meshes(c%mid(k, a)))
                    cycle
                end if
                p%mass = c%sp(a)%m
                p%charge = z(a)*e
                p%energy = p%mass*v**2/2
                p%s_ref = s_star
                call hp_orbit_geometry(bf, p, s_star, b_ref, c%go(k, a), status)
                if (status /= 0) error stop 'orbit geometry failed'
                call hp_default_mesh(c%go(k, a), 4, 16, 6, 1e-3_dp, 0.3_dp, c%meshes(c%mid(k, a)))
            end do
        end do
        if (thin) then
            allocate (c%go(0, 0), geom(nx, ns))
            geom = c%f
            call hp_energy_setup(geom, c%ops, c%meshes, c%mid, c%grid, c%sp, c%lnl, model, &
                huge(1_8), c%sys, status, l_field=l_field, gid=gid)
        else if (n_slow > 0) then
            allocate (nodes(nx, ns))
            do a = 1, ns
                do k = 1, nx
                    if (k <= n_slow) then
                        allocate (nodes(k, a)%g, source=c%f)
                    else
                        allocate (nodes(k, a)%g, source=c%go(k, a))
                    end if
                end do
            end do
            call hp_energy_setup(nodes, c%ops, c%meshes, c%mid, c%grid, c%sp, c%lnl, model, &
                huge(1_8), c%sys, status, l_field=l_field, gid=gid)
        else
            call hp_energy_setup(c%go, c%ops, c%meshes, c%mid, c%grid, c%sp, c%lnl, model, &
                huge(1_8), c%sys, status, l_field=l_field, gid=gid)
        end if
    end subroutine build

    !> K and q of the coupled kernel of build.
    subroutine coefficients(mass, z, factor, thin, model, kc, qc, node_r, its, l_field, apart, &
            slow)
        real(dp), intent(in) :: mass(:), z(:), factor
        logical, intent(in) :: thin
        integer, intent(in) :: model
        real(dp), intent(out) :: kc(:, :, :), qc(:, :, :)
        real(dp), allocatable, intent(out) :: node_r(:, :, :, :)
        integer, intent(out) :: its
        integer, intent(in), optional :: l_field
        logical, intent(in), optional :: apart
        integer, intent(in), optional :: slow
        type(coupled_t) :: c
        integer :: status
        integer(8) :: t0, t1, rate
        call system_clock(t0, rate)
        call build(mass, z, factor, thin, model, c, status, l_field, apart, slow)
        if (status /= HP_OK) error stop 'hp_energy_setup failed'
        call hp_energy_coefficients(c%sys, 2*pi, kc, qc, node_r, its, status)
        if (status /= HP_OK) error stop 'hp_energy_coefficients failed'
        call system_clock(t1)
        print '(a,l2,a,es9.2,a,i0,a,i0,a,f8.2)', 'thin', thin, ' factor', factor, &
            ' unknowns ', c%sys%n_dof, ' GMRES ', its, ' s', real(t1 - t0, dp)/rate
        call hp_energy_free(c%sys)
    end subroutine coefficients

    !> Largest relative difference of K and q of a against b.
    pure real(dp) function rel(ka, qa, kb, qb)
        real(dp), intent(in) :: ka(:, :, :), qa(:, :, :), kb(:, :, :), qb(:, :, :)
        rel = max(maxval(abs(ka - kb)/abs(kb)), maxval(abs(qa - qb)/abs(qb)))
    end function rel

    logical function below(name, got, bound) result(ok)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, bound
        ok = got <= bound
        print '(a,es10.2,a,es10.2,a,l2)', name//' ', got, ' bound ', bound, ' pass', ok
    end function below

    !> Zero orbit width: deuterium, full operator, with the mass scaled by
    !> factor at fixed thermal speed, against the thin kernel at the same
    !> factor whose nodes are coupled through the same shared points: the
    !> difference falls at least linearly with the orbit width (measured: 1.5e-6
    !> and 5.9e-7 at 4e-3 and 2e-3). On one geometry the shared points give
    !> the kernel of the common mass to rounding.
    logical function energy_width() result(ok)
        real(dp), parameter :: factors(2) = [4e-3_dp, 2e-3_dp]
        real(dp) :: kt(1, 1, 2), qt(1, 1, 2), kf(1, 1, 2), qf(1, 1, 2), ka(1, 1, 2), &
            qa(1, 1, 2), d(2), shared
        real(dp), allocatable :: node_r(:, :, :, :)
        integer :: its, i
        do i = 1, 2
            if (i == 1) call coefficients([md], [1.0_dp], factors(i), .true., COLLOP_FULL, &
                kt, qt, node_r, its)
            call coefficients([md], [1.0_dp], factors(i), .true., COLLOP_FULL, ka, qa, &
                node_r, its, apart=.true.)
            call coefficients([md], [1.0_dp], factors(i), .false., COLLOP_FULL, kf, qf, &
                node_r, its)
            d(i) = rel(kf, qf, ka, qa)
            if (i == 1) shared = rel(ka, qa, kt, qt)
            print '(a,es9.2,a,4es15.7)', 'factor ', factors(i), ' thin, shared points', ka, qa
            print '(a,es9.2,a,4es15.7)', 'factor ', factors(i), ' finite width       ', kf, qf
            print '(a,es9.2,a,es10.3)', 'factor ', factors(i), &
                ' finite width against thin', d(i)
        end do
        ok = below('shared points against the common mass, thin', shared, 1e-9_dp)
        ok = below('relative difference at the smaller width', d(2), 1e-5_dp) .and. ok
        ok = below('smaller over larger width, halved width', d(2)/d(1), 0.55_dp) .and. ok
    end function energy_width

    !> The slowest nodes on the thin field inside the finite-width kernel, as
    !> the device solves nodes whose orbits have no pitch point: at small
    !> width the kernel equals the finite-width one to the bound of
    !> energy_width, so both geometries carry the same measure and
    !> normalisation into the energy coupling.
    logical function energy_mixed() result(ok)
        real(dp), parameter :: factor = 2e-3_dp
        real(dp) :: kf(1, 1, 2), qf(1, 1, 2), km(1, 1, 2), qm(1, 1, 2)
        real(dp), allocatable :: node_r(:, :, :, :)
        integer :: its
        call coefficients([md], [1.0_dp], factor, .false., COLLOP_FULL, kf, qf, node_r, its)
        call coefficients([md], [1.0_dp], factor, .false., COLLOP_FULL, km, qm, node_r, its, &
            slow=3)
        print '(a,4es15.7)', 'finite width       ', kf, qf
        print '(a,4es15.7)', 'three thin nodes   ', km, qm
        ok = below('three thin slow nodes against finite width', rel(km, qm, kf, qf), 1e-5_dp)
    end function energy_mixed

    !> Pitch-angle scattering alone on orbit geometries: every node is the
    !> monoenergetic finite-width kernel at its speed and deflection rate.
    logical function energy_lorentz() result(ok)
        real(dp), parameter :: factor = 1e-2_dp
        type(coupled_t) :: c
        real(dp), allocatable :: node_r(:, :, :, :), r(:), x(:)
        real(dp) :: kc(1, 1, 2), qc(1, 1, 2), y, g, vt, tz
        type(hp_operator_t) :: op
        integer :: its, status, k, nx
        call build([md], [1.0_dp], factor, .false., COLLOP_LORENTZ, c, status)
        if (status /= HP_OK) error stop 'hp_energy_setup failed'
        call hp_energy_coefficients(c%sys, 2*pi, kc, qc, node_r, its, status)
        if (status /= HP_OK) error stop 'hp_energy_coefficients failed'
        nx = c%grid%n - 1
        x = c%grid%x(2:)
        vt = collop_vth(c%sp(1))
        tz = c%sp(1)%t/(c%sp(1)%z*e)
        allocate (r(nx))
        do k = 1, nx
            y = x(k)
            g = (erf(y) - 2*y*exp(-y*y)/sqrt(pi))/(2*y*y)
            op = c%ops(k, 1)
            op%nu = collop_nuhat(c%sp(1), c%sp(1), c%lnl(1, 1))*(erf(y) - g)/x(k)**3/(vt*x(k))
            op%v = vt*x(k)
            call hp_ntv_kernel(c%go(k, 1), c%meshes(c%mid(k, 1)), op, r(k), status)
            if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
        end do
        ok = below('node responses A_1, relative', maxval(abs(node_r(:, 1, 1, 1) - tz*x**2*r)) &
            /maxval(abs(tz*x**2*r)), 1e-7_dp)
        ok = below('node responses A_2, relative', maxval(abs(node_r(:, 1, 1, 2) - tz*x**4*r)) &
            /maxval(abs(tz*x**4*r)), 1e-7_dp) .and. ok
        call hp_energy_free(c%sys)
    end function energy_lorentz

    !> Manufactured solution: b = N h on the orbit-geometry nodes of both
    !> species (full operator), h a smooth vector plus random nodal values so
    !> that it lies outside the coarse space of the preconditioner; the
    !> solve returns h, and N u - b is at the GMRES tolerance, both in the
    !> norm GMRES works in, weighted by the speed mass v_a d3_k of each node.
    !> Weakly damped orbit-constant modes amplify the residual into the
    !> error (2.9e-9 into 2.6e-4 here), so the error bound is wider.
    logical function energy_manufactured() result(ok)
        type(coupled_t) :: c
        complex(dp), allocatable :: h(:), b(:), u(:)
        real(dp), allocatable :: re(:), im(:)
        integer :: its, status
        call build([md, me], [1.0_dp, -1.0_dp], 1e-2_dp, .false., COLLOP_FULL, c, status)
        if (status /= HP_OK) error stop 'hp_energy_setup failed'
        call smooth_vector(c, h)
        allocate (b(size(h)), u(size(h)), re(size(h)), im(size(h)))
        call random_number(re)
        call random_number(im)
        h = h + cmplx(re - 0.5_dp, im - 0.5_dp, dp)
        call hp_energy_apply(c%sys, h, b)
        call hp_energy_solve(c%sys, b, u, its, status)
        if (status /= HP_OK) error stop 'hp_energy_solve failed'
        print '(a,i0)', 'GMRES iterations ', its
        call speed_mass(c, re)
        print '(a,es10.2)', 'largest nodal error, relative ', maxval(abs(u - h))/maxval(abs(h))
        ok = below('solution error in the speed-mass norm, relative', norm2(re*abs(u - h)) &
            /norm2(re*abs(h)), 1e-3_dp)
        call hp_energy_apply(c%sys, u, h)
        ok = below('residual N u - b in the speed-mass norm, relative', norm2(re*abs(h - b)) &
            /norm2(re*abs(b)), 1e-7_dp) .and. ok
        call hp_energy_free(c%sys)
    end function energy_manufactured

    !> Electrons at their own mass and 10 keV, full operator: the orbit
    !> width is a sixtieth of that of deuterium, and the finite-width kernel
    !> agrees with the thin one.
    logical function energy_electrons() result(ok)
        real(dp) :: kt(1, 1, 2), qt(1, 1, 2), kf(1, 1, 2), qf(1, 1, 2)
        real(dp), allocatable :: node_r(:, :, :, :)
        integer :: its
        call coefficients([me], [-1.0_dp], 1.0_dp, .true., COLLOP_FULL, kt, qt, node_r, its)
        call coefficients([me], [-1.0_dp], 1.0_dp, .false., COLLOP_FULL, kf, qf, node_r, its)
        print '(a,4es14.6)', 'thin K q', kt, qt
        print '(a,4es14.6)', 'fow  K q', kf, qf
        ok = below('largest relative difference', rel(kf, qf, kt, qt), 1e-4_dp)
    end function energy_electrons

    !> v_a d3_k of the node of every unknown of c.
    subroutine speed_mass(c, w)
        type(coupled_t), intent(in) :: c
        real(dp), intent(out) :: w(:)
        real(dp) :: x
        integer :: a, k
        do a = 1, size(c%sp)
            do k = 1, c%grid%n - 1
                x = c%grid%x(k + 1)
                associate (nd => c%sys%node(k, a))
                    w(nd%off + 1:nd%off + nd%n) = collop_vth(c%sp(a))*c%grid%w(k + 1)*x**3 &
                        *exp(-x*x)
                end associate
            end do
        end do
    end subroutine speed_mass

    !> Random multiples of (1 + c x^2) xi^l exp(i m theta), l, m <= 2, on the
    !> nodes of c, fixed seed.
    subroutine smooth_vector(c, h)
        type(coupled_t), intent(in) :: c
        complex(dp), allocatable, intent(out) :: h(:)
        real(dp) :: cf(9, 4), x
        integer :: a, k, q, n
        call random_seed(size=n)
        call random_seed(put=[(17*k + 3, k=1, n)])
        allocate (h(c%sys%n_dof))
        do a = 1, size(c%sp)
            call random_number(cf)
            cf = cf - 0.5_dp
            do k = 1, c%grid%n - 1
                x = c%grid%x(k + 1)
                associate (nd => c%sys%node(k, a))
                    h(nd%off + 1:nd%off + nd%n) = 0
                    do q = 0, 8
                        h(nd%off + 1:nd%off + nd%n) = h(nd%off + 1:nd%off + nd%n) &
                            + cmplx(cf(q + 1, 1) + cf(q + 1, 2)*x**2, cf(q + 1, 3) &
                            + cf(q + 1, 4)*x**2, dp)*nd%terms%xi**(q/3) &
                            *exp(cmplx(0.0_dp, mod(q, 3)*nd%terms%theta, dp))
                    end do
                end associate
            end do
        end do
    end subroutine smooth_vector

    !> Moments of hp_energy_collide on a smooth vector, deuterium and
    !> electrons, thin and at two orbit widths. Density of each species,
    !> sum_k v_a d3_k 1^T y_ak, against the defect of the dropped x = 0 row
    !> (its energy and field rows on the points and lines of the measure
    !> nodes): exact at every width. Total energy sum_a n_a T_a sum_k x_k^2
    !> v_a d3_k 1^T y_ak and momentum sum_a n_a m_a v_a sum_k x_k v_a d3_k
    !> xi^T y_ak relative to the same sums of |y|: exact (energy) in the thin
    !> limit, and at finite width a defect of the orbit measure's speed
    !> dependence that falls with the width.
    logical function energy_conservation() result(ok)
        real(dp), parameter :: factors(3) = [0.0_dp, 4e-3_dp, 2e-3_dp]
        real(dp) :: dens(2), en(3), mom(3)
        integer :: i
        ok = .true.
        do i = 1, 3
            call moments(factors(i), dens, en(i), mom(i))
            print '(a,es9.2,a,2es10.2,a,es10.2,a,es10.2)', 'factor ', factors(i), &
                ' density ', dens, ' energy ', en(i), ' momentum ', mom(i)
            ok = below('density minus defect, relative', maxval(dens), 1e-9_dp) .and. ok
        end do
        ok = below('energy, thin, relative', en(1), 1e-10_dp) .and. ok
        ok = below('energy at half the width over the width', en(3)/en(2), 0.7_dp) .and. ok
    end function energy_conservation

    !> The moments of energy_conservation at mass factor (0: thin field).
    subroutine moments(factor, dens, en, mom)
        use milca_dk_hp, only: hp_terms_t, hp_lines_t, hp_lines, hp_line_ops, hp_geometry_t
        use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
        use fortsparse, only: csc_t
        real(dp), intent(in) :: factor
        real(dp), intent(out) :: dens(2), en, mom
        integer, parameter :: l_field = 3
        type(coupled_t) :: c
        type(hp_lines_t) :: lines
        type(csc_t), allocatable :: g(:), d(:)
        complex(dp), allocatable :: h(:), y(:)
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), d3(:), d4(:), d6(:), &
            vt(:), tz(:), e0w(:, :, :), f0w(:, :, :, :, :), x(:), hr(:), wp(:), col(:, :)
        real(dp) :: want, got, size_, scale, en_sum, mom_sum, momv
        integer :: status, a, b, k, j, q, l, nx, pk
        logical :: thin
        thin = factor == 0
        call build([md, me], [1.0_dp, -1.0_dp], merge(1e-3_dp, factor, thin), thin, &
            COLLOP_FULL, c, status, l_field)
        if (status /= HP_OK) error stop 'hp_energy_setup failed'
        nx = c%grid%n - 1
        x = c%grid%x(2:)
        call smooth_vector(c, h)
        allocate (y(size(h)))
        call hp_energy_collide(c%sys, h, y)
        call ntv_collision_rows(c%sp, c%lnl, c%grid, COLLOP_FULL, l_field + 1, rows, ratio, &
            d3, d4, d6, vt, tz, status, energy0_w=e0w, field0_w=f0w)
        if (status /= NTV_OK) error stop 'ntv_collision_rows failed'
        d3 = c%grid%w(2:)*x**3*exp(-x**2)
        call hp_lines(c%meshes, lines)
        en = 0
        en_sum = 0
        momv = 0
        mom_sum = 0
        do a = 1, 2
            got = 0
            size_ = 0
            do k = 1, nx
                associate (nd => c%sys%node(k, a))
                    associate (yk => y(nd%off + 1:nd%off + nd%n))
                        got = got + vt(a)*d3(k)*real(sum(yk), dp)
                        size_ = size_ + vt(a)*d3(k)*sum(abs(yk))
                        scale = c%sp(a)%n*c%sp(a)%t*x(k)**2*vt(a)*d3(k)
                        en = en + scale*real(sum(yk), dp)
                        en_sum = en_sum + scale*sum(abs(yk))
                        scale = c%sp(a)%n*c%sp(a)%m*vt(a)*x(k)*vt(a)*d3(k)
                        momv = momv + scale*real(sum(nd%terms%xi*yk), dp)
                        mom_sum = mom_sum + scale*sum(abs(nd%terms%xi*yk))
                    end associate
                end associate
            end do
            ! The x = 0 row: energy scattering on the points of measure node q,
            ! field part on the lines with its measure.
            want = 0
            do q = 1, c%grid%n
                pk = max(q - 1, 1)
                wp = reshape(c%sys%node(pk, a)%terms%q_w, [size(c%sys%node(pk, a)%terms%q_w)])
                do j = 1, nx
                    if (e0w(q, j, a) == 0) cycle
                    hr = real(h(c%sys%node(j, a)%off + 1:c%sys%node(j, a)%off &
                        + c%sys%node(j, a)%n), dp)
                    if (thin) then
                        hr = at_points(c%f, c%meshes(c%mid(j, a)), c%sys%node(j, a)%terms, &
                            c%sys%node(pk, a)%terms, hr)
                    else
                        hr = at_points(c%go(j, a), c%meshes(c%mid(j, a)), &
                            c%sys%node(j, a)%terms, c%sys%node(pk, a)%terms, hr)
                    end if
                    want = want - e0w(q, j, a)*sum(wp*hr)
                end do
                if (thin) then
                    call hp_line_ops(c%f, c%meshes(1), c%ops(1, 1), lines, l_field, g, d, &
                        status)
                else
                    call hp_line_ops(c%go(pk, a), c%meshes(c%mid(pk, a)), c%ops(pk, a), lines, &
                        l_field, g, d, status, measure=c%go(pk, a), &
                        measure_mesh=c%meshes(c%mid(pk, a)))
                end if
                if (status /= HP_OK) error stop 'hp_line_ops failed'
                allocate (col(size(lines%theta), 0:l_field))
                do l = 0, l_field
                    col(:, l) = column_sums(d(l))
                end do
                do b = 1, 2
                    do j = 1, nx
                        associate (nj => c%sys%node(j, b))
                            hr = real(h(nj%off + 1:nj%off + nj%n), dp)
                            do l = 0, l_field
                                if (f0w(q, j, a, b, l) == 0) cycle
                                want = want - f0w(q, j, a, b, l)*sum(col(:, l) &
                                    *csc_times(nj%g(l), hr))
                            end do
                        end associate
                    end do
                end do
                deallocate (col)
            end do
            dens(a) = abs(got - want)/size_
            print '(a,i0,a,es14.6,a,es14.6)', 'species ', a, ': density moment ', got, &
                ' x = 0 defect ', want
        end do
        en = abs(en)/en_sum
        mom = abs(momv)/mom_sum
        call hp_energy_free(c%sys)
    end subroutine moments

    function column_sums(m) result(v)
        use fortsparse, only: csc_t
        type(csc_t), intent(in) :: m
        real(dp) :: v(m%ncol)
        integer :: j
        do j = 1, m%ncol
            v(j) = sum(m%val(m%col_ptr(j):m%col_ptr(j + 1) - 1))
        end do
    end function column_sums

    function csc_times(m, x) result(v)
        use fortsparse, only: csc_t
        type(csc_t), intent(in) :: m
        real(dp), intent(in) :: x(:)
        real(dp) :: v(m%nrow)
        integer :: j, i
        v = 0
        do j = 1, m%ncol
            do i = m%col_ptr(j), m%col_ptr(j + 1) - 1
                v(m%row_idx(i)) = v(m%row_idx(i)) + m%val(i)*x(j)
            end do
        end do
    end function csc_times

    !> The interpolant of x, on mesh of geometry f with the node map of terms,
    !> at the mass quadrature points of terms_at (hp_locate).
    function at_points(f, mesh, terms, terms_at, x) result(v)
        use milca_dk_hp, only: hp_terms_t, hp_locate, hp_geometry_t
        use fortnum_quadrature, only: gauss_lobatto_legendre
        use fortnum_polynomial, only: lagrange_weights
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_terms_t), intent(in) :: terms, terms_at
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: v(:), r(:), s(:)
        integer, allocatable :: elem(:)
        real(dp) :: node(mesh%p + 1), wn(mesh%p + 1), w1(mesh%p + 1), w2(mesh%p + 1)
        integer :: status, i, b1, b2, np
        call hp_locate(f, mesh, terms_at%q_theta, terms_at%q_xi, elem, r, s, status)
        if (status /= HP_OK) error stop 'hp_locate failed'
        np = mesh%p + 1
        call gauss_lobatto_legendre(np, node, wn)
        allocate (v(size(elem)))
        v = 0
        do i = 1, size(elem)
            call lagrange_weights(np, r(i), node, w1)
            call lagrange_weights(np, s(i), node, w2)
            do b2 = 1, np
                do b1 = 1, np
                    v(i) = v(i) + w1(b1)*w2(b2)*x(terms%loc(b1 + np*(b2 - 1), elem(i)))
                end do
            end do
        end do
    end function at_points

    !> One harmonic m = 3, n = 3 with amplitude 1e-3 b0 s^1.5.
    subroutine write_pert(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 3 3 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.765_dp, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 3, 3, 0.0_dp, 0.0_dp, 0.0_dp, 1e-3_dp*b_ref*s**1.5_dp
        end do
        close (unit)
    end subroutine write_pert

    subroutine write_file(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            ! iota = +0.765 in the right-handed chart of this left-handed file
            write (unit, '(6es24.16)') s, -0.765_dp, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, 1.7_dp, 0.0_dp, 0.0_dp, 2.3_dp
            write (unit, '(2i5,4es24.16)') 1, 0, 0.5_dp*sqrt(s), 0.5_dp*sqrt(s), 0.0_dp, &
                -0.23_dp*sqrt(s)
        end do
        close (unit)
    end subroutine write_file
end program test_ntv_orbit_geometry
