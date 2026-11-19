program test_equilibrium_import
    ! Importers of foreign equilibrium inputs (docs/equilibrium-import.md):
    ! every fixture is written here, its expected boundary and profiles follow
    ! from the source convention by hand or from published inputs.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_import, only: imported_equilibrium_t, import_write, profile_value, &
        boundary_point, IMPORT_OK, IMPORT_UNSUPPORTED
    use milca_import_any, only: import_equilibrium
    use fortio_hdf5_writer, only: hdf5_writer_t
    use fortio_status, only: fortio_status_t
    use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
        nf90_close, NF90_CLOBBER, NF90_DOUBLE, NF90_INT, NF90_NOERR
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi
    character(4096) :: work
    integer :: failures = 0
    if (command_argument_count() /= 1) error stop 'usage: test_equilibrium_import WORKDIR'
    call get_command_argument(1, work)
    call execute_command_line('mkdir -p '//trim(work))
    call gvec_cylindrical()
    call gvec_bspline()
    call gvec_frenet_circle()
    call gvec_gframe()
    call wout()
    call desc_dshape()
    call spectre_currents()
    call vmec_variants()
    if (failures > 0) then
        print '(i0,a)', failures, ' import checks failed'
        error stop 1
    end if
    print '(a)', 'equilibrium import checks passed'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(a)', 'FAIL: '//what
        end if
    end subroutine

    logical function near(a, b, tol)
        real(dp), intent(in) :: a, b, tol
        near = abs(a - b) <= tol
    end function

    function path(name) result(out)
        character(*), intent(in) :: name
        character(len=:), allocatable :: out
        out = trim(work)//'/'//name
    end function

    subroutine write_text(name, lines)
        character(*), intent(in) :: name, lines(:)
        integer :: unit, i
        open (newunit=unit, file=path(name), status='replace', action='write')
        do i = 1, size(lines)
            write (unit, '(a)') trim(lines(i))
        end do
        close (unit)
    end subroutine

    subroutine gvec_cylindrical()
        ! GVEC hmap 1: R = X1, Z = X2, phi = -zeta, so X1_b_cos(m;n) is
        ! rbc(-n, m); iota and the flux change sign (GVEC util flip_parameters_zeta).
        ! R = 4 + cos t + 0.4 cos(t - 2 zeta), Z = sin t - 0.4 sin(t - 2 zeta)
        ! runs counterclockwise, so no poloidal flip follows.
        character(80) :: ini(12), toml(17), yaml(18)
        type(imported_equilibrium_t), allocatable :: eq, other
        character(1024) :: label
        integer :: k
        logical :: written
        ini = [character(80) :: 'which_hmap = 1  ! cylindrical', 'nfp = 2', &
            'X1_mn_max = (/ 2, 2/)', 'X1_sin_cos = _cos_', 'X2_sin_cos = _sin_', &
            'X1_b_cos( 0; 0) = 4.0', 'X1_b_cos( 1; 0) = 1.0 # comment', 'X2_b_sin(1;0) = 1.0', &
            'X1_b_cos( 1; 1) = 0.4', 'X2_b_sin( 1; 1) = -0.4', &
            'iota_type = polynomial', 'iota_coefs = (/ 0.9, &']
        call write_text('gvec.ini', [ini, [character(80) :: '  0.1 /)', 'iota_scale = -1.', &
            'pres_type = polynomial', 'pres_coefs = (/1.0, -1.0/)', 'pres_scale = 1500.', &
            'PHIEDGE = -1.5', 'X1_b_sin(1;1) = 7.0  ! outside X1_sin_cos']])
        toml = [character(80) :: 'which_hmap = 1', 'nfp = 2', 'PHIEDGE = -1.5', &
            'X1_sin_cos = "_cos_"', 'X2_sin_cos = "_sin_"', '[iota]', 'type = "polynomial"', &
            'coefs = [0.9,', '  0.1]', 'scale = -1.0', '[pres]', 'coefs = [1500.0, -1500.0]', &
            '[X1_b_cos]', '"(0, 0)" = 4.0', '"(1, 0)" = 1.0', '"(1, 1)" = 0.4', '[X2_b_sin]']
        call write_text('gvec.toml', [toml, [character(80) :: '"(1, 0)" = 1.0', '"(1, 1)" = -0.4']])
        yaml = [character(80) :: '# GVEC YAML', 'which_hmap: 1', 'nfp: 2', 'PHIEDGE: -1.5', &
            'stages:', '- minimize_tol: 0.001', '  sgrid:', '    nElems: 3', 'iota:', &
            '  type: polynomial', '  coefs: [0.9, 0.1]', '  scale: -1.0', 'pres:', &
            '  coefs: [1500.0, -1500.0]', 'X1_b_cos:', '  (0, 0): 4.0', '  (1, 0): 1.0', '  (1, 1): 0.4']
        call write_text('gvec.yaml', [yaml, [character(80) :: 'X2_b_sin:', '  (1, 0): 1.0', &
            '  (1, 1): -0.4']])
        allocate (eq, other)
        call import_equilibrium(path('gvec.ini'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'GVEC INI imported: '//trim(eq%reason))
        call check(eq%nfp == 2 .and. near(eq%phiedge, 1.5_dp, 0.0_dp), 'GVEC nfp and |phiedge|')
        call check(near(eq%rbc(0, 0), 4.0_dp, 1e-14_dp) .and. near(eq%rbc(0, 1), 1.0_dp, 1e-14_dp) &
            .and. near(eq%zbs(0, 1), 1.0_dp, 1e-14_dp), 'GVEC axisymmetric modes')
        call check(near(eq%rbc(-1, 1), 0.4_dp, 1e-14_dp) .and. near(eq%zbs(-1, 1), -0.4_dp, 1e-14_dp) &
            .and. eq%rbc(1, 1) == 0 .and. eq%rbs(-1, 1) == 0, 'GVEC (m, n) -> VMEC (-n, m), sin of a cos basis dropped')
        call check(near(profile_value(eq%iota, 0.5_dp), 0.95_dp, 1e-14_dp), 'GVEC iota sign')
        call check(near(profile_value(eq%pressure, 0.25_dp), 1125.0_dp, 1e-10_dp), 'GVEC pressure')
        call import_write(eq, path('gvec.nml'), written)
        call check(written .and. readable(path('gvec.nml')), 'namelist readable as nested_input')
        do k = 1, 2
            if (k == 1) call import_equilibrium(path('gvec.toml'), 'auto', other)
            if (k == 2) call import_equilibrium(path('gvec.yaml'), 'auto', other)
            label = merge('TOML', 'YAML', k == 1)
            call check(other%status == IMPORT_OK .and. all(other%rbc == eq%rbc) .and. &
                all(other%zbs == eq%zbs) .and. near(profile_value(other%iota, 0.3_dp), &
                profile_value(eq%iota, 0.3_dp), 1e-15_dp) .and. near(profile_value(other%pressure, &
                0.3_dp), profile_value(eq%pressure, 0.3_dp), 1e-10_dp), 'GVEC '//trim(label)//' = INI')
        end do
    end subroutine

    logical function readable(file)
        ! the declarations of app/nested_equilibrium.f90
        character(*), intent(in) :: file
        integer, parameter :: mmax = 64, nmax = 32, kmax = 32, knot_max = 1024
        integer :: nfp, mpol, ntor, angle_mpol, max_iterations, n_knots, unit, ios
        real(dp) :: phiedge, rbc(-nmax:nmax, 0:mmax), rbs(-nmax:nmax, 0:mmax)
        real(dp) :: zbc(-nmax:nmax, 0:mmax), zbs(-nmax:nmax, 0:mmax)
        real(dp) :: iota(0:kmax), pressure(0:kmax), current(0:kmax)
        real(dp) :: profile_knots(knot_max), iota_values(knot_max), iota_slopes(knot_max)
        real(dp) :: pressure_values(knot_max), pressure_slopes(knot_max)
        real(dp) :: current_values(knot_max), current_slopes(knot_max)
        character(16) :: iota_type, pressure_type, current_type
        logical :: prescribe_current
        namelist /nested_input/ nfp, mpol, ntor, rbc, rbs, zbc, zbs, phiedge, iota, pressure, &
            angle_mpol, max_iterations, prescribe_current, current, iota_type, pressure_type, &
            current_type, n_knots, profile_knots, iota_values, iota_slopes, pressure_values, &
            pressure_slopes, current_values, current_slopes
        rbc = 0; zbs = 0; iota = 0
        open (newunit=unit, file=file, status='old', action='read')
        read (unit, nml=nested_input, iostat=ios)
        close (unit)
        readable = ios == 0 .and. nfp == 2 .and. abs(rbc(-1, 1) - 0.4_dp) < 1e-15_dp .and. &
            abs(iota(1) - 0.1_dp) < 1e-15_dp
    end function

    subroutine gvec_bspline()
        ! a clamped cubic B-spline on one span is the Bernstein polynomial
        ! sum c_i C(3,i) s^i (1-s)^(3-i); a C1 cubic is held exactly
        type(imported_equilibrium_t), allocatable :: eq
        real(dp) :: s, expected
        allocate (eq)
        call write_text('bspline.ini', [character(80) :: 'X1_b_cos(0;0) = 3.0', 'X1_b_cos(1;0) = 1.0', &
            'X2_b_sin(1;0) = 1.0', 'iota_type = bspline', 'iota_knots = (/0,0,0,0,0.5,1,1,1,1/)', &
            'iota_coefs = (/0.3, 0.5, 0.2, 0.9, 0.4/)', 'pres_type = interpolation', &
            'pres_rho2 = (/0.0, 0.25, 0.5, 0.75, 1.0/)', 'pres_vals = (/4.0, 3.0, 2.0, 1.0, 0.0/)', &
            'pres_scale = 10.0'])
        call import_equilibrium(path('bspline.ini'), 'gvec', eq)
        call check(eq%status == IMPORT_OK .and. eq%iota%kind == 'hermite', 'GVEC B-spline imported')
        ! s in [0, 0.5]: knots 0,0,0,0,0.5,1 -> de Boor by hand at s = 0.2
        s = 0.2_dp
        expected = bspline_reference(s)
        call check(near(profile_value(eq%iota, s), -expected, 1e-14_dp), 'GVEC cubic B-spline exact')
        call check(near(profile_value(eq%pressure, 0.6_dp), 10*(4 - 4*0.6_dp), 1e-12_dp), &
            'GVEC not-a-knot spline reproduces a line')
    end subroutine

    real(dp) function bspline_reference(s)
        ! knots (0,0,0,0,0.5,1,1,1,1): on [0, 0.5] with u = 2 s the basis is
        ! B1 = (1-u)^3, B2 = u (1-u)^2 (3 - 3u/2)... evaluated by Cox-de Boor
        real(dp), intent(in) :: s
        real(dp) :: t(9), c(5), b(8), left, right
        integer :: i, k
        t = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        c = [0.3_dp, 0.5_dp, 0.2_dp, 0.9_dp, 0.4_dp]
        b = 0
        b(4) = 1
        do k = 1, 3
            do i = 1, 8 - k
                left = 0; right = 0
                if (t(i + k) > t(i)) left = (s - t(i))/(t(i + k) - t(i))*b(i)
                if (t(i + k + 1) > t(i + 1)) right = (t(i + k + 1) - s)/(t(i + k + 1) - t(i + 1))*b(i + 1)
                b(i) = left + right
            end do
        end do
        bspline_reference = sum(c*b(:5))
    end function

    subroutine gvec_frenet_circle()
        ! hmap 20 on the circle R0 = 3: T = e_phi, B = e_Z, N = -e_R, phi = zeta,
        ! so R = 3 - X1, Z = X2. X1 = 0.5 cos t + 0.05 cos(t - 2 zeta),
        ! X2 = 0.6 sin t is clockwise; after t -> -t, rbc(0,1) = -0.5,
        ! rbc(-1,1) = -0.05, zbs(0,1) = -0.6, and iota changes sign twice.
        type(imported_equilibrium_t), allocatable :: eq
        allocate (eq)
        call write_text('frenet.ini', [character(80) :: 'which_hmap = 20', 'nfp = 2', 'hmap_nfp = 2', &
            'hmap_n_max = 1', 'hmap_rc = (/3.0, 0.0/)', 'hmap_zs = (/0.0, 0.0/)', &
            'X1_b_cos(1;0) = 0.5', 'X1_b_cos(1;1) = 0.05', 'X2_b_sin(1;0) = 0.6', &
            'iota_type = polynomial', 'iota_coefs = (/0.4/)', 'PHIEDGE = 0.7'])
        call import_equilibrium(path('frenet.ini'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'Frenet imported: '//trim(eq%reason))
        call check(near(eq%rbc(0, 0), 3.0_dp, 1e-12_dp) .and. near(eq%rbc(0, 1), -0.5_dp, 1e-12_dp) &
            .and. near(eq%rbc(-1, 1), -0.05_dp, 1e-12_dp) .and. near(eq%zbs(0, 1), -0.6_dp, 1e-12_dp), &
            'Frenet circle refitted on cylindrical planes')
        call check(eq%geometry_error < 1e-10_dp, 'Frenet refit error')
        call check(near(profile_value(eq%iota, 0.5_dp), -0.4_dp, 1e-15_dp) .and. &
            near(eq%phiedge, 0.7_dp, 1e-15_dp), 'Frenet: phi = zeta keeps iota, the theta flip reverses it')
    end subroutine

    subroutine gvec_gframe()
        ! G-frame (hmap 21) about the circle R0 = 2 with N tilted towards
        ! e_phi by 0.3 rad, so a cross-section is not a phi plane. The raw
        ! grid points must lie on the imported surface.
        integer, parameter :: nfp = 3, nzeta = 15, ntheta = 9
        type(hdf5_writer_t) :: w
        type(fortio_status_t) :: status
        type(imported_equilibrium_t), allocatable :: eq
        real(dp) :: xyz(nfp*nzeta, 3), nv(nfp*nzeta, 3), bv(nfp*nzeta, 3)
        real(dp) :: x(nzeta, ntheta), y(nzeta, ntheta), zeta(nzeta), theta(ntheta), z, p(3)
        real(dp) :: worst, alpha
        integer :: i, j
        alpha = 0.3_dp
        do j = 1, nfp*nzeta
            z = 2*pi*(j - 1)/(nfp*nzeta)
            xyz(j, :) = [2*cos(z), 2*sin(z), 0.0_dp]
            nv(j, :) = cos(alpha)*[cos(z), sin(z), 0.0_dp] + sin(alpha)*[-sin(z), cos(z), 0.0_dp]
            bv(j, :) = [0.0_dp, 0.0_dp, 1.0_dp]
        end do
        do j = 1, nzeta
            zeta(j) = 2*pi*(j - 1)/(nfp*nzeta)
        end do
        do i = 1, ntheta
            theta(i) = 2*pi*(i - 1)/ntheta
            x(:, i) = 0.4_dp*cos(theta(i)) + 0.03_dp*cos(theta(i) - nfp*zeta)
            y(:, i) = 0.5_dp*sin(theta(i)) + 0.02_dp*sin(2*theta(i) - nfp*zeta)
        end do
        call w%create(path('gframe.nc'), status)
        call w%add_i32_scalar('NFP', nfp, status)
        call w%define_group('axis', status)
        call w%add_i32_scalar('axis/nzeta', nzeta, status)
        call w%add_i32_scalar('axis/n_max', 4, status)
        call w%add_r64_1('axis/zeta(:)', zeta, status)
        call w%add_r64_2('axis/xyz(::)', xyz, status)
        call w%add_r64_2('axis/Nxyz(::)', nv, status)
        call w%add_r64_2('axis/Bxyz(::)', bv, status)
        call w%define_group('boundary', status)
        call w%add_i32_scalar('boundary/ntheta', ntheta, status)
        call w%add_i32_scalar('boundary/nzeta', nzeta, status)
        call w%add_i32_scalar('boundary/m_max', 4, status)
        call w%add_i32_scalar('boundary/n_max', 7, status)
        call w%add_i32_scalar('boundary/lasym', 0, status)
        call w%add_r64_1('boundary/theta(:)', theta, status)
        call w%add_r64_1('boundary/zeta(:)', zeta, status)
        call w%add_r64_2('boundary/X(::)', x, status)
        call w%add_r64_2('boundary/Y(::)', y, status)
        call w%close(status)
        call check(status%ok(), 'G-frame fixture written')
        call write_text('gframe.ini', [character(80) :: 'which_hmap = 21', 'nfp = 3', &
            'hmap_ncfile = gframe.nc', 'getBoundaryFromFile = 1', 'boundary_filename = gframe.nc', &
            'iota_type = polynomial', 'iota_coefs = (/0.5/)', 'PHIEDGE = 1.0'])
        allocate (eq)
        call import_equilibrium(path('gframe.ini'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'G-frame imported: '//trim(eq%reason))
        worst = 0
        do j = 1, nzeta
            do i = 1, ntheta
                p = xyz(j, :) + x(j, i)*nv(j, :) + y(j, i)*bv(j, :)
                worst = max(worst, distance(eq, hypot(p(1), p(2)), atan2(p(2), p(1)), p(3)))
            end do
        end do
        call check(worst < 1e-8_dp, 'G-frame grid points on the imported surface')
    end subroutine

    real(dp) function distance(eq, r, phi, z)
        ! distance of (r, z) from the section phi of the surface, by a
        ! dense scan in theta refined by golden sections
        type(imported_equilibrium_t), intent(in) :: eq
        real(dp), intent(in) :: r, phi, z
        real(dp) :: a, b, c, d, t0
        integer :: k, best
        real(dp) :: dk(720)
        do k = 1, 720
            dk(k) = gap(eq, 2*pi*(k - 1)/720, r, phi, z)
        end do
        best = minloc(dk, 1)
        t0 = 2*pi*(best - 1)/720
        a = t0 - 2*pi/720; b = t0 + 2*pi/720
        do k = 1, 80
            c = b - (b - a)*0.618033988749895_dp
            d = a + (b - a)*0.618033988749895_dp
            if (gap(eq, c, r, phi, z) < gap(eq, d, r, phi, z)) then
                b = d
            else
                a = c
            end if
        end do
        distance = gap(eq, (a + b)/2, r, phi, z)
    end function

    real(dp) function gap(eq, t, r, phi, z)
        type(imported_equilibrium_t), intent(in) :: eq
        real(dp), intent(in) :: t, r, phi, z
        real(dp) :: rr, zz
        call boundary_point(eq, t, phi, rr, zz)
        gap = hypot(rr - r, zz - z)
    end function

    subroutine wout()
        ! iotaf, presf cubic in s are held exactly by the not-a-knot spline
        integer, parameter :: ns = 7, mnmax = 3
        type(imported_equilibrium_t), allocatable :: eq
        real(dp) :: s(ns), rmnc(mnmax, ns), zmns(mnmax, ns)
        integer :: ncid, dr, dm, v(10), j, code
        do j = 1, ns
            s(j) = (j - 1.0_dp)/(ns - 1)
        end do
        ! clockwise: Z = -0.5 sin t; a helical (m, n) = (1, 1) term
        rmnc = 0; zmns = 0
        rmnc(1, :) = 5; rmnc(2, :) = sqrt(s); rmnc(3, :) = 0.1_dp*s
        zmns(2, :) = -0.5_dp*sqrt(s); zmns(3, :) = 0.1_dp*s
        code = nf90_create(path('wout_test.nc'), NF90_CLOBBER, ncid)
        code = nf90_def_dim(ncid, 'radius', ns, dr)
        code = nf90_def_dim(ncid, 'mn_mode', mnmax, dm)
        code = nf90_def_var(ncid, 'ns', NF90_INT, v(1))
        code = nf90_def_var(ncid, 'nfp', NF90_INT, v(2))
        code = nf90_def_var(ncid, 'mnmax', NF90_INT, v(3))
        code = nf90_def_var(ncid, 'lasym__logical__', NF90_INT, v(4))
        code = nf90_def_var(ncid, 'xm', NF90_DOUBLE, [dm], v(5))
        code = nf90_def_var(ncid, 'xn', NF90_DOUBLE, [dm], v(6))
        code = nf90_def_var(ncid, 'iotaf', NF90_DOUBLE, [dr], v(7))
        code = nf90_def_var(ncid, 'presf', NF90_DOUBLE, [dr], v(8))
        code = nf90_def_var(ncid, 'phi', NF90_DOUBLE, [dr], v(9))
        code = nf90_def_var(ncid, 'rmnc', NF90_DOUBLE, [dm, dr], v(10))
        code = nf90_def_var(ncid, 'zmns', NF90_DOUBLE, [dm, dr], j)
        code = nf90_enddef(ncid)
        code = nf90_put_var(ncid, v(1), ns)
        code = nf90_put_var(ncid, v(2), 4)
        code = nf90_put_var(ncid, v(3), mnmax)
        code = nf90_put_var(ncid, v(4), 0)
        code = nf90_put_var(ncid, v(5), [0.0_dp, 1.0_dp, 1.0_dp])
        code = nf90_put_var(ncid, v(6), [0.0_dp, 0.0_dp, 4.0_dp])
        code = nf90_put_var(ncid, v(7), 0.4_dp + 0.2_dp*s - 0.3_dp*s**3)
        code = nf90_put_var(ncid, v(8), 1e3_dp*(1 - s)**2)
        code = nf90_put_var(ncid, v(9), -0.8_dp*s)
        code = nf90_put_var(ncid, v(10), rmnc)
        code = nf90_put_var(ncid, j, zmns)
        code = nf90_close(ncid)
        call check(code == NF90_NOERR, 'wout fixture written')
        allocate (eq)
        call import_equilibrium(path('wout_test.nc'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'wout imported: '//trim(eq%reason))
        ! the theta flip maps (n, m) = (1, 1) to (-1, 1) and negates sines and iota
        call check(eq%nfp == 4 .and. near(eq%rbc(0, 1), 1.0_dp, 1e-15_dp) .and. &
            near(eq%zbs(0, 1), 0.5_dp, 1e-15_dp) .and. near(eq%rbc(-1, 1), 0.1_dp, 1e-15_dp) .and. &
            near(eq%zbs(-1, 1), -0.1_dp, 1e-15_dp), 'wout boundary, orientation flipped')
        call check(near(profile_value(eq%iota, 0.37_dp), -(0.4_dp + 0.2_dp*0.37_dp - 0.3_dp*0.37_dp**3), &
            1e-14_dp) .and. near(profile_value(eq%pressure, 0.61_dp), 1e3_dp*0.39_dp**2, 1e-10_dp), &
            'wout profiles exact for cubics')
        call check(near(eq%phiedge, 0.8_dp, 1e-15_dp), 'wout |phiedge|')
    end subroutine

    subroutine desc_dshape()
        ! DESC's DSHAPE input (DESC tests/inputs/input.DSHAPE_desc) against the
        ! VMEC deck it was generated from (Hirshman and Whitson 1983):
        ! RBC(0,0..2) = 3.51, 1.0, 0.106, ZBS(0,1..2) = 1.47, -0.16,
        ! AM = 1600 (1, -2, 1), AI = 1, -0.67
        type(imported_equilibrium_t), allocatable :: eq
        allocate (eq)
        call write_text('input.DSHAPE_desc', [character(80) :: '# DESC input', 'sym = 1', 'NFP =   1', &
            'Psi = 1.00000000e+00', 'M_pol = 11', 'bdry_mode = lcfs', &
            'l:   0  p =   1.60000000E+03  i =   1.00000000E+00', &
            'l:   2  p =  -3.20000000E+03  i =  -6.70000000E-01', &
            'l:   4  p =   1.60000000E+03  i =   0.00000000E+00', &
            'l:   0  m:  -2  n:   0  R1 =   0.00000000E+00  Z1 =  -1.60000000E-01', &
            'l:   0  m:  -1  n:   0  R1 =   0.00000000E+00  Z1 =   1.47000000E+00', &
            'l:   0  m:   0  n:   0  R1 =   3.51000000E+00  Z1 =   0.00000000E+00', &
            'l:   0  m:   1  n:   0  R1 =   1.00000000E+00  Z1 =   0.00000000E+00', &
            'l:   0  m:   2  n:   0  R1 =   1.06000000E-01  Z1 =   0.00000000E+00', &
            'n:   0  R0 =   3.51000000E+00  Z0 =   0.00000000E+00'])
        call import_equilibrium(path('input.DSHAPE_desc'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'DESC imported: '//trim(eq%reason))
        call check(near(eq%rbc(0, 0), 3.51_dp, 1e-14_dp) .and. near(eq%rbc(0, 1), 1.0_dp, 1e-14_dp) &
            .and. near(eq%rbc(0, 2), 0.106_dp, 1e-14_dp) .and. near(eq%zbs(0, 1), 1.47_dp, 1e-14_dp) &
            .and. near(eq%zbs(0, 2), -0.16_dp, 1e-14_dp), 'DESC boundary = VMEC DSHAPE')
        call check(near(profile_value(eq%pressure, 0.3_dp), 1600*0.49_dp, 1e-10_dp) .and. &
            near(profile_value(eq%iota, 0.3_dp), 1 - 0.67_dp*0.3_dp, 1e-14_dp), 'DESC rho^l profiles in s')
        ! a 3-D DESC boundary with a current: G_-1(t) G_1(nfp phi) = sin t cos(nfp phi)
        ! = (sin(t + nfp phi) + sin(t - nfp phi))/2; counterclockwise, so
        ! the current is c_DESC (signgs = -1, nested I = -sign(Psi) I_V)
        call write_text('input.helix_desc', [character(80) :: 'NFP = 3', 'Psi = 0.5', &
            'l: 0  p = 0.0  c = 0.0', 'l: 2  p = 0.0  c = 2.0E+03', 'l: 1  p = 1.0', &
            'm: 0  n: 0  R1 = 4.0  Z1 = 0.0', 'm: 1  n: 0  R1 = 1.0  Z1 = 0.0', &
            'm: -1  n: 0  R1 = 0.0  Z1 = 1.0', 'm: -1  n: 1  R1 = 0.0  Z1 = 0.2'])
        call import_equilibrium(path('input.helix_desc'), 'auto', eq)
        call check(eq%status == IMPORT_OK .and. eq%prescribe_current, 'DESC current input')
        call check(near(eq%zbs(1, 1), 0.1_dp, 1e-14_dp) .and. near(eq%zbs(-1, 1), 0.1_dp, 1e-14_dp), &
            'DESC product basis to VMEC basis')
        call check(near(profile_value(eq%current, 0.5_dp), 1e3_dp, 1e-9_dp), 'DESC current sign')
        call check(eq%pressure%kind == 'hermite' .and. near(profile_value(eq%pressure, 0.49_dp), 0.7_dp, &
            1e-6_dp), 'DESC odd rho power as Hermite')
    end subroutine

    subroutine spectre_currents()
        ! SPECTRE "(m, n)" keys, theta reversed against VMEC: the clockwise
        ! boundary is flipped, (m, n) = (1, 1) lands on VMEC (n, m) = (-1, 1);
        ! the enclosed current is (ivolume + sum isurf)/mu0 at each interface
        type(imported_equilibrium_t), allocatable :: eq
        allocate (eq)
        call write_text('spectre.toml', [character(80) :: '[physics]', 'igeometry = 3', 'nfp = 5', &
            'nvol = 3', 'lconstraint = 3', 'lfreebound = 0', 'phiedge = 1.0', 'pscale = 1.0', &
            'tflux = [0.5, 1.0, 2.0]', 'pressure = [0.0, 0.0, 0.0]', 'iota = [0.27, 0.27, 0.272, 0.275]', &
            'ivolume = [-1e-5, -2e-5,', ' -4e-5]', 'isurf = [-1e-6, -2e-6, -3e-6]', '', '[physics.rbc]', &
            '"(0, 0)" = 10.0', '"(1, 0)" = 1.0', '"(1, 1)" = 0.25', '[physics.zbs]', '"(1, 0)" = -1.0', &
            '"(1, 1)" = 0.25', '[physics.rbs]', '[numeric]', 'linitialize = 1'])
        call import_equilibrium(path('spectre.toml'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'SPECTRE imported: '//trim(eq%reason))
        call check(near(eq%rbc(0, 1), 1.0_dp, 1e-15_dp) .and. near(eq%zbs(0, 1), 1.0_dp, 1e-15_dp) .and. &
            near(eq%rbc(-1, 1), 0.25_dp, 1e-15_dp) .and. near(eq%zbs(-1, 1), -0.25_dp, 1e-15_dp), &
            'SPECTRE boundary flipped')
        call check(near(profile_value(eq%current, 1.0_dp), (-4e-5_dp - 6e-6_dp)/mu0, 1e-9_dp) .and. &
            near(profile_value(eq%current, 0.25_dp), (-1e-5_dp - 1e-6_dp)/mu0, 1e-9_dp), &
            'SPECTRE interface currents')
        call write_text('spectre_free.toml', [character(80) :: '[physics]', 'igeometry = 3', &
            'lfreebound = 1'])
        call import_equilibrium(path('spectre_free.toml'), 'spectre', eq)
        call check(eq%status == IMPORT_UNSUPPORTED, 'SPECTRE free boundary refused')
    end subroutine

    subroutine vmec_variants()
        type(imported_equilibrium_t), allocatable :: eq
        allocate (eq)
        call write_text('input.member', [character(80) :: '&INDATA', ' NFP = 2, NCURR = 1,', &
            ' PHIEDGE = -2.0  CURTOR = 1.0E5', ' AC = 1.0 1.0', ' AM = 2*1.0E3', &
            ' AI = 0.1', ' AI =', ' RBC(0,0) = 5.0  RBC(0,1) = 1.0', ' ZBS(0,1) = 1.0', &
            ' RBC(1,1) = 0.1 ZBS(1,1) = 0.1', '/'])
        call write_text('input.sequence', [character(80) :: '&VSEQ NSEQ = 1', &
            ' EXTENSION = "other", "member"', '&END'])
        call import_equilibrium(path('input.sequence'), 'auto', eq)
        call check(eq%status == IMPORT_OK, 'VSEQ follows to its last member: '//trim(eq%reason))
        ! I(s) = CURTOR (s + s^2/2)/(3/2); nested I = -sign(PHIEDGE) I
        call check(near(profile_value(eq%current, 0.5_dp), 1e5_dp*(0.5_dp + 0.125_dp)/1.5_dp, 1e-8_dp) &
            .and. eq%prescribe_current, 'INDATA current')
        call check(near(profile_value(eq%pressure, 0.5_dp), 1.5e3_dp, 1e-10_dp), 'INDATA repeat count')
        call check(near(eq%rbc(1, 1), 0.1_dp, 1e-15_dp), 'INDATA boundary')
        call write_text('input.free', [character(80) :: '&INDATA', " MGRID_FILE = 'mgrid.x'", &
            ' RBC(0,0) = 5.0', '/'])
        call import_equilibrium(path('input.free'), 'vmec', eq)
        call check(eq%status == IMPORT_UNSUPPORTED, 'VMEC free boundary refused')
        call write_text('input.mirror', [character(80) :: '&MIRROR', ' NS = 7', '/'])
        call import_equilibrium(path('input.mirror'), 'auto', eq)
        call check(eq%status == IMPORT_UNSUPPORTED, 'VMEX mirror refused')
        call write_text('cyl.ini', [character(80) :: 'which_hmap = 3', 'X1_b_cos(1;0) = 1.0'])
        call import_equilibrium(path('cyl.ini'), 'auto', eq)
        call check(eq%status == IMPORT_UNSUPPORTED, 'GVEC straight cylinder refused')
    end subroutine
end program test_equilibrium_import
