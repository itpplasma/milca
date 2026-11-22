program test_flux_map_qa
    !! Opt-in driver (not a registered test): the 3D flux map
    !! (milca_flux_map_3d) on a fixed-boundary VMEC deck with p = 0 and a
    !! prescribed rotational transform, by default the Landreman-Paul QA deck
    !! of benchmarks/label-qa. s is the normalised toroidal flux,
    !! B = PHIEDGE/(2 pi) [2 grad x X grad y + iota(s) grad phi X grad s].
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
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_import, only: imported_equilibrium_t, IMPORT_OK
    use milca_import_vmec, only: import_vmec
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_system, flux_map3d_solve
    use flux_map_sections
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nsamples = 256, nout = 64, nray = 128, max_stages = 8
    type(imported_equilibrium_t) :: eq
    type(fm3_case_t), target :: c, previous
    character(512) :: deck = '../benchmarks/label-qa/input.landreman-paul-qa-iota', &
                      outdir = 'flux-map-qa-out', stages = '2:2:9', arg
    integer :: stage_spec(3, max_stages), nstages, istage, k, st, u, newton_total, j
    real(dp) :: phiedge, kappa = 1, t0, t1, ts, energy, phi, axis(2), x0(2), sym
    real(dp), allocatable :: iota_c(:), samples(:, :, :), quality(:), flux(:)
    real(dp) :: bmin, bmax
    logical :: verbose
    integer :: skipped = 0

    if (command_argument_count() >= 1) call get_command_argument(1, deck)
    if (command_argument_count() >= 2) call get_command_argument(2, outdir)
    if (command_argument_count() >= 3) call get_command_argument(3, stages)
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) kappa
    end if
    verbose = command_argument_count() >= 5
    call parse_stages()
    call wall_time(t0)
    call import_vmec(trim(deck), eq)
    if (eq%status /= IMPORT_OK) then
        print '(2a)', 'import: ', trim(eq%reason)
        error stop 'deck'
    end if
    phiedge = eq%phiedge
    clebsch_nfp = eq%nfp
    iota_c = eq%iota%coefficients
    print '(a,i0,a,f10.6,a,7f10.6)', 'NFP ', eq%nfp, '  PHIEDGE', phiedge, '  iota', iota_c
    newton_total = 0
    do istage = 1, nstages
        call wall_time(ts)
        associate (deg => stage_spec(1, istage), cells => stage_spec(2, istage), &
                   nt => stage_spec(3, istage))
            if (allocated(samples)) deallocate (samples, quality, flux)
            allocate (samples(2, nsamples, nt), quality(nt), flux(nt))
            do k = 1, nt
                call plane_curve(2*pi*(k - 1)/(nt*eq%nfp), samples(:, :, k))
            end do
            call fm3_build(c, samples, cells, deg, 4, st)
            if (st /= 0) error stop 'meshes'
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
        end associate
        previous = c
    end do
    call wall_time(t1)

    call execute_command_line('mkdir -p '//trim(outdir))
    open (newunit=u, file=trim(outdir)//'/axis.txt', status='replace')
    write (u, '(a)') '# phi R Z'
    do j = 0, nout - 1
        phi = 2*pi*j/(nout*eq%nfp)
        call fm3_axis(c, j/real(nout, dp), profiles, axis, st, x0)
        if (st /= 0) error stop 'axis'
        write (u, '(3es24.16)') phi, axis
        call surface(0.25_dp, '0.25', phi, j/real(nout, dp), x0, j == 0)
        call surface(0.50_dp, '0.50', phi, j/real(nout, dp), x0, j == 0)
    end do
    close (u)
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

    subroutine surface(level, tag, phi, t, x0, first)
        !! The level set s_h = level on the plane t, parametrised by the
        !! angle theta* of the map (the preimage of the circle), with |B|.
        real(dp), intent(in) :: level, phi, t, x0(2)
        character(*), intent(in) :: tag
        logical, intent(in) :: first
        real(dp) :: e(2), x(2), rz(2), b(3)
        integer :: v, q, it, stl, nc

        if (first) then
            open (newunit=v, file=trim(outdir)//'/surface_s'//tag//'.txt', status='replace')
            write (v, '(a)') '# phi angle R Z |B|'
        else
            open (newunit=v, file=trim(outdir)//'/surface_s'//tag//'.txt', position='append', &
                  status='old')
        end if
        do q = 1, nray
            e = [cos(2*pi*(q - 1)/nray), sin(2*pi*(q - 1)/nray)]
            ! continuation in the radius of the disc from the axis; a ray
            ! that fails is retried with finer steps, then skipped
            do nc = 16, 256, 240
                x = x0
                do it = 1, nc
                    call fm3_preimage(c, t, profiles, sqrt(level)*it/nc*e, x, rz, b, stl)
                    if (stl /= 0) exit
                end do
                if (stl == 0) exit
            end do
            if (stl /= 0) then
                skipped = skipped + 1
                cycle
            end if
            write (v, '(5es24.16)') phi, 2*pi*(q - 1)/nray, rz, norm2(b)
        end do
        close (v)
    end subroutine surface

    subroutine plane_curve(phi, curve)
        !! Wall curve of the deck in the plane phi, counter-clockwise, first
        !! sample on the outboard side at the centroid height.
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: curve(:, :)
        real(dp), allocatable :: raw(:, :)
        real(dp) :: th, area, a, best, cen(2)
        integer :: ii, m, n, first, ns

        ns = size(curve, 2)
        allocate (raw(2, ns))
        do ii = 1, ns
            th = 2*pi*(ii - 1)/ns
            raw(:, ii) = 0
            do m = lbound(eq%rbc, 2), ubound(eq%rbc, 2)
                do n = lbound(eq%rbc, 1), ubound(eq%rbc, 1)
                    a = m*th - n*eq%nfp*phi
                    raw(1, ii) = raw(1, ii) + eq%rbc(n, m)*cos(a) + eq%rbs(n, m)*sin(a)
                    raw(2, ii) = raw(2, ii) + eq%zbc(n, m)*cos(a) + eq%zbs(n, m)*sin(a)
                end do
            end do
        end do
        area = 0
        cen = 0
        do ii = 1, ns
            associate (p => raw(:, ii), q => raw(:, mod(ii, ns) + 1))
                a = p(1)*q(2) - q(1)*p(2)
                area = area + a/2
                cen = cen + (p + q)*a/6
            end associate
        end do
        cen = cen/area
        if (area < 0) raw = raw(:, ns:1:-1)
        best = huge(1.0_dp)
        first = 1
        do ii = 1, ns
            if (raw(1, ii) > cen(1) .and. abs(raw(2, ii) - cen(2)) < best) then
                best = abs(raw(2, ii) - cen(2))
                first = ii
            end if
        end do
        curve = cshift(raw, first - 1, dim=2)
    end subroutine plane_curve

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
        p = 0
    end subroutine profiles

end program test_flux_map_qa
