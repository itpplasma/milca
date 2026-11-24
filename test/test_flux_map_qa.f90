program test_flux_map_qa
    !! Opt-in driver (not a registered test): the 3D flux map
    !! (milca_flux_map_3d) on a fixed-boundary VMEC deck with a prescribed
    !! rotational transform and the deck's power-series pressure (zero for the
    !! default, the Landreman-Paul QA deck of benchmarks/label-qa). s is the
    !! normalised toroidal flux, B = PHIEDGE/(2 pi) [2 grad x X grad y +
    !! iota(s) grad phi X grad s], W = int (|B|^2/2 - mu0 p) dV (B in T).
    !! The planes span one field period; the wall curves are the only
    !! geometric input (meshes on their area centroids, polar start after
    !! the second-moment normalisation). Stages are solved in turn, each from
    !! the previous one (fm3_prolong).
    !!
    !! Writes axis.txt, surface_s0.25.txt, surface_s0.50.txt and summary.json
    !! in the layout benchmarks/label-qa/compare.py reads (s = 1 is the given
    !! wall and is not written).
    !!
    !! Arguments: deck, output directory, stages "degree:cells:planes,..."
    !! (default 2:2:9), iota scale (1), verbose.
    !!
    !! With MILCA_FM3_VMEC=<wout file> every stage also evaluates the VMEC
    !! state of that file on its mesh (flux_map_vmec, both orientations of
    !! theta*): W and the scaled gradient of the flux map there, the gradient
    !! by bands of s, and the state Newton reaches from it; before the stages
    !! the force balance J X B - grad p of the VMEC state itself
    !! (milca_flux_map_force; with the stages argument "vmec" only this).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_import, only: imported_equilibrium_t, IMPORT_OK, profile_value
    use milca_import_vmec, only: import_vmec
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_solve
    use fortsparse_csr, only: csr_t
    use fortfem_spline_eval, only: spline_point_matrix
    use flux_map_sections
    use flux_map_vmec, only: fm3_from_vmec, fm3_from_vmec_failures, fm3_vmec_force
    use milca_flux_map_force, only: flux_map_force_t, flux_map3d_force, flux_map_force_print, &
        flux_map3d_force_point
    use milca_vmec_state, only: vmec_state_t, vmec_state_read
    use milca_nested_oracles, only: NESTED_SHEARED
    use landreman_sections, only: section_oracle, section_exact_b
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nout = 64, nray = 128, max_stages = 8
    integer :: nsamples = 256
    type(imported_equilibrium_t) :: eq
    type(fm3_case_t), target :: c, previous, axis_seed
    character(512) :: deck = '../benchmarks/label-qa/input.landreman-paul-qa-iota', &
                      outdir = 'flux-map-qa-out', stages = '2:2:9', arg
    integer :: eval_planes, stage_spec(3, max_stages), nstages, istage, k, st, u, newton_total, j, solve_status
    real(dp) :: phiedge, kappa = 1, t0, t1, ts, energy, phi, axis(2), x0(2)
    real(dp) :: deck_angle_sense, deck_angle_origin
    real(dp), allocatable :: iota_c(:), pres_c(:), samples(:, :, :), quality(:), flux(:)
    real(dp), allocatable :: axes(:, :), rings(:, :, :, :), bring(:, :, :, :)
    logical, allocatable :: fring(:, :, :)
    integer, allocatable :: axis_status(:)
    logical :: pres_hermite = .false.
    real(dp) :: bmin, bmax
    logical :: verbose
    integer :: skipped = 0
    type(vmec_state_t) :: vstate
    character(512) :: vpath, vreason
    integer :: vlen
    logical :: manufactured, sheared
    real(dp) :: magnetic_error(2)
    real(dp) :: pscale = 1 ! diagnostic: MILCA_FM3_PSCALE multiplies the pressure
    real(dp) :: grading = 0 ! MILCA_FM3_GRADING: curvature grading of the section meshes
    ! diagnostic: MILCA_FM3_ORIGIN moves the first wall sample (disc angle 0,
    ! and with it the block corners on the wall) by this fraction of the
    ! deck angle
    real(dp) :: origin_shift = 0
    ! convergence basin: MILCA_FM3_START_AXIS = d moves the start u = G to
    ! u + d (1 - |u|^2) e_x (the axis of the start by about d minor radii);
    ! MILCA_FM3_START_ANGLE = a turns it by the angle a sin(theta), the wall
    ! angle included
    real(dp) :: start_axis = 0, start_angle = 0
    type(flux_map_force_t) :: force
    real(dp), allocatable :: winding(:), least_step(:)
    real(dp) :: jump

    if (command_argument_count() >= 1) call get_command_argument(1, deck)
    if (command_argument_count() >= 2) call get_command_argument(2, outdir)
    if (command_argument_count() >= 3) call get_command_argument(3, stages)
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) kappa
    end if
    verbose = command_argument_count() >= 5
    if (trim(stages) /= 'vmec') call parse_stages()
    call get_environment_variable('MILCA_FM3_PSCALE', arg, vlen)
    if (vlen > 0) read (arg, *) pscale
    call get_environment_variable('MILCA_FM3_GRADING', arg, vlen)
    if (vlen > 0) read (arg, *) grading
    call get_environment_variable('MILCA_FM3_START_AXIS', arg, vlen)
    if (vlen > 0) read (arg, *) start_axis
    call get_environment_variable('MILCA_FM3_START_ANGLE', arg, vlen)
    if (vlen > 0) read (arg, *) start_angle
    call get_environment_variable('MILCA_FM3_ORIGIN', arg, vlen)
    if (vlen > 0) read (arg, *) origin_shift
    call get_environment_variable('MILCA_FM3_EXACT', arg, vlen)
    sheared = arg == 'landreman-sheared'
    manufactured = arg == 'landreman-integer' .or. sheared
    if (vlen > 0 .and. .not. manufactured) error stop 'unknown exact field'
    if (sheared) then
        section_oracle%family = NESTED_SHEARED
        section_oracle%eps = 0.6_dp
        section_oracle%saxis = 2.4_dp
        section_oracle%delta = 0.2_dp
    else
        section_oracle%eps = 0.5_dp
        section_oracle%delta = 1.0_dp/64
    end if
    call get_environment_variable('MILCA_FM3_VMEC', vpath, vlen)
    if (vlen > 0) then
        call vmec_state_read(trim(vpath), vstate, st, vreason)
        if (st /= 0) then
            print '(2a)', 'VMEC state: ', trim(vreason)
            error stop 'VMEC state'
        end if
    end if
    call wall_time(t0)
    call import_vmec(trim(deck), eq)
    if (eq%status /= IMPORT_OK) then
        print '(2a)', 'import: ', trim(eq%reason)
        error stop 'deck'
    end if
    phiedge = eq%phiedge
    clebsch_nfp = eq%nfp
    iota_c = eq%iota%coefficients
    if (trim(eq%iota%kind) /= 'power_series') error stop 'iota: power series required'
    select case (trim(eq%pressure%kind))
    case ('none')
    case ('power_series')
        if (any(eq%pressure%coefficients /= 0)) pres_c = eq%pressure%coefficients
    case ('hermite')
        pres_hermite = .true.
    case default
        error stop 'pressure: power series or Hermite spline required'
    end select
    print '(a,i0,a,f10.6,a,*(1x,es16.8))', 'NFP ', eq%nfp, '  PHIEDGE', phiedge, '  iota', iota_c
    if (vlen > 0) call vmec_force()
    if (trim(stages) == 'vmec') then
        if (vlen == 0) error stop 'VMEC force diagnostic needs MILCA_FM3_VMEC'
        stop
    end if
    call exact_origin(0.0_dp, deck_angle_origin, deck_angle_sense)
    print '(a,f4.0)', 'poloidal deck winding ', deck_angle_sense
    newton_total = 0
    do istage = 1, nstages
        call wall_time(ts)
        associate (deg => stage_spec(1, istage), cells => stage_spec(2, istage), &
                   nt => stage_spec(3, istage))
            if (allocated(quality)) deallocate (quality, flux)
            allocate (quality(nt), flux(nt))
            ! wall samples doubled until the harmonic section maps resolve the curves
            nsamples = 256
            do
                if (allocated(samples)) deallocate (samples)
                allocate (samples(2, nsamples, nt))
                do k = 1, nt
                    call plane_curve(2*pi*(k - 1)/(nt*eq%nfp), samples(:, :, k))
                end do
                call fm3_build(c, samples, cells, deg, 4, st, preserve_area=keep_area(), &
                               grading=grading)
                if (st == 0 .and. c%map_residual < 1e-9_dp) exit
                if (st /= 0 .and. st /= 13) exit
                if (nsamples >= 2048) exit
                nsamples = 2*nsamples
            end do
            if (st /= 0) then
                print '(a,i0)', 'fm3_build status ', st
                error stop 'meshes'
            end if
            print '(a,i0,a,es8.1,a,es8.1,a,es9.1,a,es8.1)', '  meshes: wall samples ', nsamples, &
                '  harmonic map residual', c%map_residual, '  wall distance', &
                maxval(c%wall_distance), ' m  area error', maxval(abs(c%area_error)), &
                '  certified Jacobian ratio', minval(c%jacobian_ratio)
            call fm3_mesh_quality(c, quality)
            if (istage == 1) then
                if (sheared) then
                    call fm3_start_polar(c)
                else
                    call fm3_start(c)
                end if
                call perturb_start()
            else
                call fm3_prolong(previous, c, st)
                if (st /= 0) error stop 'prolongation'
                c%fm%stationary = sheared
            end if
            c%fm%verbose = verbose
            call flux_map3d_solve(c%fm, c%derham, c%quad, profiles, c%base, 1e-9_dp, 80, solve_status)
            if (manufactured) then
                call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, energy, &
                    exact_field=section_exact_b, field_error=magnetic_error)
                print '(a,es12.4)', '  volume relative L2 error of B ', sqrt(magnetic_error(1)/magnetic_error(2))
            else
                call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, energy)
            end if
            call fm3_dump_hessian(c, profiles)
            call fm3_flux(c, profiles, flux, bmin, bmax)
            newton_total = newton_total + c%fm%steps
            call wall_time(t1)
            print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es9.2,a,f6.3,a,es17.10,a,es9.2,a,f8.2,a,i0)', &
                'degree ', deg, '  cells ', cells, '  planes ', nt, '  unknowns ', 2*c%n0*nt, &
                '  status ', solve_status, '  Newton ', c%fm%steps, '  |g|', c%fm%gradient_norm, &
                '  worst Jacobian ratio', minval(quality), '  W', energy, '  flux defect', &
                maxval(abs(flux - phiedge))/abs(phiedge), '  ', t1 - ts, ' s  negative eigenvalues ', &
                c%fm%negative_eigenvalues
            call flux_map3d_force(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, 10, force, st)
            if (st /= 0) error stop 'force balance'
            call flux_map_force_print(force, 'CEMV')
            call force_dump()
            call wall_misfit()
            allocate (winding(nt), least_step(nt))
            call fm3_wall_winding(c, profiles, winding, least_step, st)
            if (st /= 0) error stop 'wall winding'
            call fm3_interface_jump(c, profiles, jump, st)
            if (st /= 0) error stop 'interface jump'
            print '(a,2f10.6,a,es10.2,a,es10.2)', '  winding of the wall angle per plane, min max', &
                minval(winding), maxval(winding), '  least step of the wall angle / mean step', &
                minval(least_step), '  jump of B^2/2 at the block interfaces (relative rms)', jump
            deallocate (winding, least_step)
            if (solve_status /= 0) error stop 'CEMV Newton did not converge'
            if (vlen > 0) call vmec_check(energy)
        end associate
        if (sheared .and. istage == 1 .and. nstages > 1) axis_seed = c
        previous = c
    end do
    call wall_time(t1)

    call execute_command_line('mkdir -p '//trim(outdir))
    ! The one-plane map has the same section at every output angle.
    eval_planes = merge(1, nout, c%nt == 1)
    allocate (axes(2, 0:eval_planes - 1), rings(2, nray, 2, 0:eval_planes - 1), &
              bring(3, nray, 2, 0:eval_planes - 1), fring(nray, 2, 0:eval_planes - 1), &
              axis_status(0:eval_planes - 1))
    !$omp parallel do schedule(dynamic) private(x0, st)
    do j = 0, eval_planes - 1
        call fm3_axis(c, j/real(nout, dp), profiles, axes(:, j), axis_status(j), x0)
        if (axis_status(j) /= 0) cycle
        call ring(0.25_dp, j/real(nout, dp), x0, rings(:, :, 1, j), bring(:, :, 1, j), fring(:, 1, j))
        call ring(0.50_dp, j/real(nout, dp), x0, rings(:, :, 2, j), bring(:, :, 2, j), fring(:, 2, j))
    end do
    !$omp end parallel do
    if (any(axis_status /= 0)) error stop 'axis'
    skipped = count(.not. fring)*(nout/eval_planes)
    open (newunit=u, file=trim(outdir)//'/axis.txt', status='replace')
    write (u, '(a)') '# phi R Z'
    do j = 0, nout - 1
        write (u, '(3es24.16)') 2*pi*j/(nout*eq%nfp), axes(:, min(j, eval_planes - 1))
    end do
    close (u)
    call write_surface('0.25', 1)
    call write_surface('0.50', 2)
    call write_currents()
    open (newunit=u, file=trim(outdir)//'/summary.json', status='replace')
    write (u, '(a,es24.16,a,a,a,i0,a,f10.3,a)') '{"energy_T2m3": ', energy, ', "stages": "', &
        trim(stages), '", "newton_steps": ', newton_total, ', "time_s": ', t1 - t0, '}'
    close (u)
    print '(a,i0,a,f8.2,3a,i0)', 'Newton steps ', newton_total, '  time ', t1 - t0, ' s  wrote ', &
        trim(outdir), '  surface points skipped ', skipped

contains

    subroutine parse_stages()
        integer :: p, q, r
        character(512) :: rest

        rest = stages
        nstages = 0
        do while (len_trim(rest) > 0)
            p = index(rest, ',')
            if (p == 0) p = len_trim(rest) + 1
            nstages = nstages + 1
            q = index(rest(:p - 1), ':')
            r = q + index(rest(q + 1:p - 1), ':')
            read (rest(:q - 1), *) stage_spec(1, nstages)
            read (rest(q + 1:r - 1), *) stage_spec(2, nstages)
            read (rest(r + 1:p - 1), *) stage_spec(3, nstages)
            rest = rest(p + 1:)
        end do
    end subroutine parse_stages

    subroutine wall_time(t)
        real(dp), intent(out) :: t
        integer(8) :: cnt, rate
        call system_clock(cnt, rate)
        t = real(cnt, dp)/real(rate, dp)
    end subroutine wall_time

    subroutine ring(level, t, x0, rz, b, found)
        !! The level set s_h = level on the plane t at nray equally spaced
        !! angles theta* of the map (the preimage of the circle), with B.
        real(dp), intent(in) :: level, t, x0(2)
        real(dp), intent(out) :: rz(2, nray), b(3, nray)
        logical, intent(out) :: found(nray)
        real(dp) :: e(2), x(2), xprev(2)
        integer :: q, it, stl, nc
        logical :: have

        have = .false.
        do q = 1, nray
            e = [cos(2*pi*(q - 1)/nray), sin(2*pi*(q - 1)/nray)]
            ! from the previous ray's point (one Newton solve); on failure
            ! continuation in the radius of the disc from the axis, then
            ! finer, then the ray is skipped
            stl = 1
            if (have) then
                x = xprev
                call fm3_preimage(c, t, profiles, sqrt(level)*e, x, rz(:, q), b(:, q), stl)
            end if
            if (stl /= 0) then
                do nc = 16, 256, 240
                    x = x0
                    do it = 1, nc
                        call fm3_preimage(c, t, profiles, sqrt(level)*it/nc*e, x, rz(:, q), &
                                          b(:, q), stl)
                        if (stl /= 0) exit
                    end do
                    if (stl == 0) exit
                end do
            end if
            found(q) = stl == 0
            have = found(q)
            if (have) xprev = x
        end do
    end subroutine ring

    subroutine write_surface(tag, l)
        !! surface_s<tag>.txt from the rings of level index l, with |B|.
        character(*), intent(in) :: tag
        integer, intent(in) :: l
        integer :: v, q, jj, plane

        open (newunit=v, file=trim(outdir)//'/surface_s'//tag//'.txt', status='replace')
        write (v, '(a)') '# phi angle R Z |B|'
        do jj = 0, nout - 1
            plane = min(jj, eval_planes - 1)
            do q = 1, nray
                if (fring(q, l, plane)) write (v, '(5es24.16)') 2*pi*jj/(nout*eq%nfp), &
                    2*pi*(q - 1)/nray, rings(:, q, l, plane), norm2(bring(:, q, l, plane))
            end do
        end do
        close (v)
    end subroutine write_surface

    subroutine enclosed_current(level, t, x0, circulation, reference, ok)
        !! circulation = the loop integral of B along the surface s_h = level
        !! in the plane t (mu0 times the toroidal current through the section
        !! inside it; zero for a vacuum field), reference = the integral of
        !! |B.dl| along the same loop. The tangent is the derivative of the
        !! trigonometric interpolant in theta*.
        real(dp), intent(in) :: level, t, x0(2)
        real(dp), intent(out) :: circulation, reference
        logical, intent(out) :: ok
        real(dp) :: rz(2, nray), b(3, nray), d(2), bt
        complex(dp) :: cf(2, -(nray/2 - 1):nray/2 - 1)
        logical :: found(nray)
        integer :: q, m, j

        circulation = 0
        reference = 0
        call ring(level, t, x0, rz, b, found)
        ok = all(found)
        if (.not. ok) return
        cf = 0
        do m = -(nray/2 - 1), nray/2 - 1
            do j = 1, nray
                cf(:, m) = cf(:, m) + rz(:, j)*exp(cmplx(0.0_dp, -2*pi*m*(j - 1)/real(nray, dp), dp))
            end do
        end do
        do q = 1, nray
            d = 0
            do m = -(nray/2 - 1), nray/2 - 1
                d = d + real(cmplx(0.0_dp, real(m, dp), dp)*cf(:, m)/nray* &
                             exp(cmplx(0.0_dp, 2*pi*m*(q - 1)/real(nray, dp), dp)), dp)
            end do
            bt = b(1, q)*d(1) + b(3, q)*d(2)
            circulation = circulation + bt*2*pi/nray
            reference = reference + abs(bt)*2*pi/nray
        end do
    end subroutine enclosed_current

    subroutine write_currents()
        !! current.txt: per plane (4 of the output planes) and surface the
        !! circulation of B and its ratio to the integral of |B.dl|; the
        !! plane means are printed.
        real(dp), parameter :: levels(6) = [0.05_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.9_dp]
        real(dp) :: circ(6), ref(6), mean(6), ax(2), xx(2)
        integer :: v, jj, l, stc, cnt(6)
        logical :: ok(6)

        open (newunit=v, file=trim(outdir)//'/current.txt', status='replace')
        write (v, '(a)') '# phi s circulation_of_B[Tm] ratio_to_integral_of_|B.dl|'
        mean = 0
        cnt = 0
        do jj = 0, nout - 1, nout/4
            if (jj == 0 .or. c%nt > 1) then
                call fm3_axis(c, jj/real(nout, dp), profiles, ax, stc, xx)
                ok = .false.
                if (stc == 0) then
                    do l = 1, 6
                        call enclosed_current(levels(l), jj/real(nout, dp), xx, circ(l), ref(l), ok(l))
                    end do
                end if
            end if
            do l = 1, 6
                if (.not. ok(l)) cycle
                write (v, '(4es24.16)') 2*pi*jj/(nout*eq%nfp), levels(l), circ(l), circ(l)/ref(l)
                mean(l) = mean(l) + circ(l)/ref(l)
                cnt(l) = cnt(l) + 1
            end do
        end do
        close (v)
        print '(a,6f6.2)', 'enclosed toroidal current / poloidal circulation at s =', levels
        print '(a,6es10.2)', '   plane mean', mean/max(cnt, 1)
    end subroutine write_currents

    subroutine force_dump()
        !! diagnostic: MILCA_FM3_FORCE_DUMP=prefix writes, per quadrature
        !! point of the section and 2 nt values of t, block, xi, t, R, Z, s,
        !! the volume weight, B, J, F and grad(B^2/2) to
        !! prefix_<degree>_<cells>_<planes>.txt.
        character(512) :: prefix, path
        integer :: length, v, p, m, stp
        real(dp) :: rz(2), b(3), jc(3), f(3), sl, gb(3), vol, t

        call get_environment_variable('MILCA_FM3_FORCE_DUMP', prefix, length)
        if (length == 0) return
        write (path, '(a,a,i0,a,i0,a,i0,a)') trim(prefix), '_', c%degree, '_', c%cells, '_', &
            c%nt, '.txt'
        open (newunit=v, file=trim(path), status='replace')
        write (v, '(a)') '# block xi1 xi2 t R Z s weight B(3) J(3) F(3) gradB2(3)'
        do p = 1, c%quad%point_count
            do m = 1, 2*c%nt
                t = (m - 1)/real(2*c%nt, dp)
                call flux_map3d_force_point(c%fm, c%derham, profiles, c%base, c%fm%z, &
                    c%quad%block(p), c%quad%xi(:, p), t, rz, b, jc, f, stp, sl, gb, vol)
                if (stp /= 0) error stop 'force dump'
                write (v, '(i3,19es16.8)') c%quad%block(p), c%quad%xi(:, p), t, rz, sl, &
                    c%quad%weight(p)*vol/(2*c%nt), b, jc, f, gb
            end do
        end do
        close (v)
    end subroutine force_dump

    subroutine perturb_start()
        !! The start of the first stage moved off u = G (see start_axis,
        !! start_angle).
        integer :: i, kk, ix, iy
        real(dp) :: x, y, th, da

        if (start_axis == 0 .and. start_angle == 0) return
        do kk = 1, c%nt
            do i = 1, c%n0
                ix = (kk - 1)*c%n0 + i
                iy = c%n0*c%nt + ix
                if (c%fm%is_wall(i)) then
                    c%fm%z(ix) = c%fm%z(ix) + start_angle*sin(c%fm%z(ix))
                else
                    x = c%fm%z(ix)
                    y = c%fm%z(iy)
                    th = atan2(y, x)
                    da = start_angle*sin(th)
                    c%fm%z(ix) = cos(da)*x - sin(da)*y + start_axis*(1 - x**2 - y**2)
                    c%fm%z(iy) = sin(da)*x + cos(da)*y
                end if
            end do
        end do
    end subroutine perturb_start

    subroutine wall_misfit()
        !! diagnostic: MILCA_FM3_MISFIT=prefix writes per plane, outer block
        !! and wall cell the largest and the signed mean distance of the mesh
        !! wall from the wall curve (8 points per cell) to
        !! prefix_<degree>_<cells>_<planes>.txt.
        integer, parameter :: nd = 4096
        character(512) :: prefix, path
        integer :: length, v, k, b, e, q, i, ib
        real(dp) :: dense(2, nd), pos(2), d, best, dmax, a(2), ab(2), tt, th0, sense, phi, pr(2)

        call get_environment_variable('MILCA_FM3_MISFIT', prefix, length)
        if (length == 0) return
        write (path, '(a,a,i0,a,i0,a,i0,a)') trim(prefix), '_', c%degree, '_', c%cells, '_', &
            c%nt, '.txt'
        open (newunit=v, file=trim(path), status='replace')
        write (v, '(a)') '# plane block cell max_distance R Z'
        do k = 1, c%nt
            phi = 2*pi*(k - 1)/(c%nt*eq%nfp)
            call exact_origin(phi, th0, sense)
            do i = 1, nd
                call wall_point(2*pi*(i - 1)/nd, phi, dense(:, i))
            end do
            do b = 2, 5
                do e = 1, c%cells
                    dmax = 0
                    do q = 1, 8
                        call wpoint(b, (e - 1 + (q - 0.5_dp)/8)/c%cells, k, pos)
                        best = huge(1.0_dp)
                        do i = 1, nd
                            ib = mod(i, nd) + 1
                            a = dense(:, i)
                            ab = dense(:, ib) - a
                            tt = max(0.0_dp, min(1.0_dp, dot_product(pos - a, ab)/dot_product(ab, ab)))
                            d = norm2(pos - a - tt*ab)
                            best = min(best, d)
                        end do
                        dmax = max(dmax, best)
                        if (q == 4) pr = pos
                    end do
                    write (v, '(3i5,3es16.8)') k, b, e, dmax, pr
                end do
            end do
        end do
        close (v)
    end subroutine wall_misfit

    subroutine wpoint(b, x, kk, pos)
        !! Point of the mesh wall of plane kk: outer block b, wall parameter x.
        integer, intent(in) :: b, kk
        real(dp), intent(in) :: x
        real(dp), intent(out) :: pos(2)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:)
        integer :: st2
        call spline_point_matrix(c%derham%space, [b], reshape([x, 0.0_dp], [2, 1]), 0, rowptr, &
                                 colidx, val, st2)
        if (st2 /= 0) error stop 'wall point'
        pos = matmul(c%base(:, colidx, kk), val)
    end subroutine wpoint

    subroutine vmec_force()
        !! VMEC's straight angle obeys d(theta + lambda)/d(zeta) = iota.
        !! Reflecting it would reverse the transform for the same profile.
        type(flux_map_force_t) :: fv
        integer :: bands, length, radial
        character(16) :: v

        ! diagnostic: MILCA_FM3_BANDS = number of bands of s (10)
        bands = 10
        call get_environment_variable('MILCA_FM3_BANDS', v, length)
        if (length > 0) read (v, *) bands
        if (bands < 10 .or. mod(bands, 10) /= 0) error stop 'force bands must be a multiple of 10'
        radial = max(200, 4*bands)
        call get_environment_variable('MILCA_FM3_RADIAL', v, length)
        if (length > 0) read (v, *) radial
        if (radial < bands) error stop 'radial quadrature needs at least one point per band'
        if (manufactured) then
            call fm3_vmec_force(vstate, deck_profiles, 1.0_dp, bands, radial, fv, &
                                section_exact_b, magnetic_error)
            print '(a,es12.4)', '  volume relative L2 error of B ', sqrt(magnetic_error(1)/magnetic_error(2))
        else
            call fm3_vmec_force(vstate, deck_profiles, 1.0_dp, bands, radial, fv)
        end if
        call flux_map_force_print(fv, 'VMEC state')
    end subroutine vmec_force

    subroutine vmec_check(w_ours)
        !! The VMEC state on the mesh of c: W and gradient of the flux map.
        real(dp), intent(in) :: w_ours
        type(fm3_case_t) :: cv
        type(csr_t) :: h
        real(dp), allocatable :: g(:), d(:), sc(:)
        real(dp) :: wv, orient, edges(6), band(5), sv
        integer :: o, i, kk, n0, nt, b, stv


        edges = [0.0_dp, 0.05_dp, 0.25_dp, 0.5_dp, 0.8_dp, huge(1.0_dp)]
        do o = 1, 2
            orient = merge(1.0_dp, -1.0_dp, o == 1)
            cv = c
            call fm3_from_vmec(cv, vstate, orient, stv)
            if (stv /= 0) error stop 'fm3_from_vmec'
            call flux_map3d_system(cv%fm, cv%derham, cv%quad, profiles, cv%base, cv%fm%z, wv, g, h)
            d = diag_of(h)
            sc = g/sqrt(abs(d) + tiny(1.0_dp))
            sc(cv%fm%fixed) = 0
            n0 = cv%n0
            nt = cv%nt
            band = 0
            do kk = 1, nt
                do i = 1, n0
                    if (cv%fm%is_wall(i)) then
                        band(5) = band(5) + sc((kk - 1)*n0 + i)**2
                        cycle
                    end if
                    sv = cv%fm%z((kk - 1)*n0 + i)**2 + cv%fm%z(n0*nt + (kk - 1)*n0 + i)**2
                    b = count(sv >= edges(2:5)) + 1
                    band(b) = band(b) + sc((kk - 1)*n0 + i)**2 + sc(n0*nt + (kk - 1)*n0 + i)**2
                end do
            end do
            print '(a,f5.1,a,es17.10,a,es10.3,a,es9.2,a,5es9.2,a,i0)', '  VMEC state, orientation', &
                orient, ':  W', wv, '  W - W_flux map', wv - w_ours, '  |g| scaled', &
                norm2(sc), '  by s bands [0,.05,.25,.5,.8,wall]', sqrt(band), &
                '  inversion failures ', fm3_from_vmec_failures()
            if (o == 1 .or. wv - w_ours < 1e-2_dp*abs(w_ours)) then
                cv%fm%verbose = .false.
                call flux_map3d_solve(cv%fm, cv%derham, cv%quad, profiles, cv%base, 1e-9_dp, 80, stv)
                call flux_map3d_system(cv%fm, cv%derham, cv%quad, profiles, cv%base, cv%fm%z, wv)
                print '(a,i0,a,i0,a,es17.10,a,es10.3,a,es10.3)', '    Newton from it: status ', stv, &
                    '  steps ', cv%fm%steps, '  W', wv, '  W - W_flux map', wv - w_ours, &
                    '  max |z - z_flux map| (up to a rotation of the disc)', &
                    maxval(abs(cv%fm%z - c%fm%z))
            end if
        end do
    end subroutine vmec_check

    function diag_of(a) result(dg)
        type(csr_t), intent(in) :: a
        real(dp), allocatable :: dg(:)
        integer :: i, q
        allocate (dg(a%nrow), source=0.0_dp)
        do i = 1, a%nrow
            do q = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(q) == i) dg(i) = a%val(q)
            end do
        end do
    end function diag_of

    subroutine plane_curve(phi, curve)
        !! Wall curve of the deck in the plane phi, counter-clockwise, equally
        !! spaced in the deck angle from the outboard point at the height of
        !! the area centroid (a smooth function of phi).
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: curve(:, :)
        real(dp) :: th0, sense, coarse_axis(2)
        integer :: ii, ns, local

        ns = size(curve, 2)
        if (axis_seed%n0 == 0) then
            call exact_origin(phi, th0, sense)
        else
            ! Fix the section anchor across the stationary continuation.
            call fm3_axis(axis_seed, phi*eq%nfp/(2*pi), profiles, coarse_axis, local)
            if (local /= 0) error stop 'cannot locate the coarse magnetic axis'
            call exact_origin(phi, th0, sense, coarse_axis(2))
        end if
        th0 = th0 + sense*2*pi*origin_shift
        do ii = 1, ns
            call wall_point(th0 + sense*2*pi*(ii - 1)/ns, phi, curve(:, ii))
        end do
    end subroutine plane_curve

    subroutine wall_point(th, phi, p)
        real(dp), intent(in) :: th, phi
        real(dp), intent(out) :: p(2)
        real(dp) :: a
        integer :: m, n
        p = 0
        do m = lbound(eq%rbc, 2), ubound(eq%rbc, 2)
            do n = lbound(eq%rbc, 1), ubound(eq%rbc, 1)
                a = m*th - n*eq%nfp*phi
                p(1) = p(1) + eq%rbc(n, m)*cos(a) + eq%rbs(n, m)*sin(a)
                p(2) = p(2) + eq%zbc(n, m)*cos(a) + eq%zbs(n, m)*sin(a)
            end do
        end do
    end subroutine wall_point

    subroutine exact_origin(phi, th0, sense, height)
        !! Outboard wall angle at the requested height or the area centroid.
        !! The winding sense is -1 for a clockwise deck angle.
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: th0, sense
        real(dp), intent(in), optional :: height
        integer, parameter :: nq = 1024
        real(dp) :: p(2), q(2), area, cen(2), a, lo, hi, best, w
        integer :: i, it, ib
        area = 0
        cen = 0
        call wall_point(0.0_dp, phi, p)
        do i = 1, nq
            call wall_point(2*pi*i/nq, phi, q)
            a = p(1)*q(2) - q(1)*p(2)
            area = area + a/2
            cen = cen + (p + q)*a/6
            p = q
        end do
        cen = cen/area
        if (present(height)) cen(2) = height
        sense = sign(1.0_dp, area)
        best = huge(1.0_dp)
        ib = 0
        do i = 0, nq - 1
            call wall_point(2*pi*i/nq, phi, p)
            if (p(1) > cen(1) .and. abs(p(2) - cen(2)) < best) then
                best = abs(p(2) - cen(2))
                ib = i
            end if
        end do
        lo = 2*pi*(ib - 1)/nq
        hi = 2*pi*(ib + 1)/nq
        do it = 1, 60
            th0 = (lo + hi)/2
            call wall_point(th0, phi, p)
            call wall_point(lo, phi, q)
            w = (p(2) - cen(2))*(q(2) - cen(2))
            if (w <= 0) then
                hi = th0
            else
                lo = th0
            end if
        end do
    end subroutine exact_origin

    logical function keep_area()
        !! diagnostic: MILCA_FM3_FREE_AREA set = the mesh wall interpolates the
        !! wall curve without the area correction
        integer :: length
        call get_environment_variable('MILCA_FM3_FREE_AREA', length=length)
        keep_area = length == 0
    end function keep_area

    subroutine profiles(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)

        call deck_profiles(s, tau, gp, p)
        ! CEMV uses counterclockwise labels in (R, Z) and signed toroidal flux.
        tau = -tau
        gp = -deck_angle_sense*gp
    end subroutine profiles

    subroutine deck_profiles(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: io(3)
        integer :: m

        io = 0
        do m = 1, size(iota_c)
            io(1) = io(1) + iota_c(m)*s**(m - 1)
            if (m >= 2) io(2) = io(2) + (m - 1)*iota_c(m)*s**(m - 2)
            if (m >= 3) io(3) = io(3) + (m - 1)*(m - 2)*iota_c(m)*s**(m - 3)
        end do
        tau = [-phiedge/pi, 0.0_dp, 0.0_dp]
        gp = -kappa*phiedge/(2*pi)*io
        ! mu0 p(s), p' and p'' of the deck's power series (Pa -> T^2)
        p = 0
        if (allocated(pres_c)) then
            do m = 1, size(pres_c)
                p(1) = p(1) + pres_c(m)*s**(m - 1)
                if (m >= 2) p(2) = p(2) + (m - 1)*pres_c(m)*s**(m - 2)
                if (m >= 3) p(3) = p(3) + (m - 1)*(m - 2)*pres_c(m)*s**(m - 3)
            end do
            p = 4e-7_dp*pi*p
        else if (pres_hermite) then
            ! sampled profiles (e.g. two_power): p', p'' by central differences
            block
                real(dp), parameter :: h = 1e-4_dp
                real(dp) :: pm, p0, pp
                p0 = profile_value(eq%pressure, s)
                pm = profile_value(eq%pressure, s - h)
                pp = profile_value(eq%pressure, s + h)
                p = 4e-7_dp*pi*[p0, (pp - pm)/(2*h), (pp - 2*p0 + pm)/h**2]
            end block
        end if
        p = pscale*p
    end subroutine deck_profiles

end program test_flux_map_qa
