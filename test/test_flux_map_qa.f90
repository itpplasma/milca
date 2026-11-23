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
    !! by bands of s, and the state Newton reaches from it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_import, only: imported_equilibrium_t, IMPORT_OK, profile_value
    use milca_import_vmec, only: import_vmec
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_solve
    use fortsparse_csr, only: csr_t
    use flux_map_sections
    use flux_map_vmec, only: fm3_from_vmec, fm3_from_vmec_failures
    use milca_vmec_state, only: vmec_state_t, vmec_state_read
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nout = 64, nray = 128, max_stages = 8
    integer :: nsamples = 256
    type(imported_equilibrium_t) :: eq
    type(fm3_case_t), target :: c, previous
    character(512) :: deck = '../benchmarks/label-qa/input.landreman-paul-qa-iota', &
                      outdir = 'flux-map-qa-out', stages = '2:2:9', arg
    integer :: stage_spec(3, max_stages), nstages, istage, k, st, u, newton_total, j
    real(dp) :: phiedge, kappa = 1, t0, t1, ts, energy, phi, axis(2), x0(2), sym
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
    real(dp) :: pscale = 1 ! diagnostic: MILCA_FM3_PSCALE multiplies the pressure

    if (command_argument_count() >= 1) call get_command_argument(1, deck)
    if (command_argument_count() >= 2) call get_command_argument(2, outdir)
    if (command_argument_count() >= 3) call get_command_argument(3, stages)
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) kappa
    end if
    verbose = command_argument_count() >= 5
    call parse_stages()
    call get_environment_variable('MILCA_FM3_PSCALE', arg, vlen)
    if (vlen > 0) read (arg, *) pscale
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
    print '(a,i0,a,f10.6,a,7f10.6)', 'NFP ', eq%nfp, '  PHIEDGE', phiedge, '  iota', iota_c
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
                call fm3_build(c, samples, cells, deg, 4, st, preserve_area=keep_area())
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
                call fm3_start(c)
            else
                call fm3_prolong(previous, c, st)
                if (st /= 0) error stop 'prolongation'
            end if
            sym = maxval(abs(samples(2, 2:nsamples/2, 1) + samples(2, nsamples:nsamples/2 + 2:-1, 1)))
            if (sym < 1e-8_dp) then
                call fm3_gauge_symmetric(c, st)
                if (st /= 0) error stop 'gauge rows'
            end if
            c%fm%verbose = verbose
            call flux_map3d_solve(c%fm, c%derham, c%quad, profiles, c%base, 1e-9_dp, 80, st)
            call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, energy)
            call fm3_dump_hessian(c, profiles)
            call fm3_flux(c, profiles, flux, bmin, bmax)
            newton_total = newton_total + c%fm%steps
            call wall_time(t1)
            print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es9.2,a,f6.3,a,es17.10,a,es9.2,a,f8.2,a,i0)', &
                'degree ', deg, '  cells ', cells, '  planes ', nt, '  unknowns ', 2*c%n0*nt, &
                '  status ', st, '  Newton ', c%fm%steps, '  |g|', c%fm%gradient_norm, &
                '  worst Jacobian ratio', minval(quality), '  W', energy, '  flux defect', &
                maxval(abs(abs(flux) - abs(phiedge)))/abs(phiedge), '  ', t1 - ts, ' s  negative eigenvalues ', &
                c%fm%negative_eigenvalues
            if (vlen > 0) call vmec_check(energy)
        end associate
        previous = c
    end do
    call wall_time(t1)

    call execute_command_line('mkdir -p '//trim(outdir))
    ! axis and the surfaces s = 0.25, 0.5 of every output plane (independent
    ! planes in parallel), then written in plane order
    allocate (axes(2, 0:nout - 1), rings(2, nray, 2, 0:nout - 1), bring(3, nray, 2, 0:nout - 1), &
              fring(nray, 2, 0:nout - 1), axis_status(0:nout - 1))
    !$omp parallel do schedule(dynamic) private(x0, st)
    do j = 0, nout - 1
        call fm3_axis(c, j/real(nout, dp), profiles, axes(:, j), axis_status(j), x0)
        if (axis_status(j) /= 0) cycle
        call ring(0.25_dp, j/real(nout, dp), x0, rings(:, :, 1, j), bring(:, :, 1, j), fring(:, 1, j))
        call ring(0.50_dp, j/real(nout, dp), x0, rings(:, :, 2, j), bring(:, :, 2, j), fring(:, 2, j))
    end do
    !$omp end parallel do
    if (any(axis_status /= 0)) error stop 'axis'
    skipped = count(.not. fring)
    open (newunit=u, file=trim(outdir)//'/axis.txt', status='replace')
    write (u, '(a)') '# phi R Z'
    do j = 0, nout - 1
        write (u, '(3es24.16)') 2*pi*j/(nout*eq%nfp), axes(:, j)
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
        integer :: v, q, jj

        open (newunit=v, file=trim(outdir)//'/surface_s'//tag//'.txt', status='replace')
        write (v, '(a)') '# phi angle R Z |B|'
        do jj = 0, nout - 1
            do q = 1, nray
                if (fring(q, l, jj)) write (v, '(5es24.16)') 2*pi*jj/(nout*eq%nfp), &
                    2*pi*(q - 1)/nray, rings(:, q, l, jj), norm2(bring(:, q, l, jj))
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
        real(dp) :: circ, ref, mean(6), ax(2), xx(2)
        integer :: v, jj, l, stc, cnt(6)
        logical :: ok

        open (newunit=v, file=trim(outdir)//'/current.txt', status='replace')
        write (v, '(a)') '# phi s circulation_of_B[Tm] ratio_to_integral_of_|B.dl|'
        mean = 0
        cnt = 0
        do jj = 0, nout - 1, nout/4
            call fm3_axis(c, jj/real(nout, dp), profiles, ax, stc, xx)
            if (stc /= 0) cycle
            do l = 1, 6
                call enclosed_current(levels(l), jj/real(nout, dp), xx, circ, ref, ok)
                if (.not. ok) cycle
                write (v, '(4es24.16)') 2*pi*jj/(nout*eq%nfp), levels(l), circ, circ/ref
                mean(l) = mean(l) + circ/ref
                cnt(l) = cnt(l) + 1
            end do
        end do
        close (v)
        print '(a,6f6.2)', 'enclosed toroidal current / poloidal circulation at s =', levels
        print '(a,6es10.2)', '   plane mean', mean/max(cnt, 1)
    end subroutine write_currents

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
                    '  max |z - z_flux map|', maxval(abs(cv%fm%z - c%fm%z))
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
        real(dp) :: th0, sense
        integer :: ii, ns

        ns = size(curve, 2)
        call exact_origin(phi, th0, sense)
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

    subroutine exact_origin(phi, th0, sense)
        !! Deck angle of the outboard wall point at the height of the area
        !! centroid (bisection), a smooth function of phi; sense = -1 for a
        !! clockwise deck angle.
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: th0, sense
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
    end subroutine profiles

end program test_flux_map_qa
