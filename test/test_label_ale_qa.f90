program test_label_ale_qa
    !! The 3D label lane on a real stellarator (docs/equilibrium-roadmap.md,
    !! step 4): Landreman–Paul QA, fixed boundary, p = 0, the rotational
    !! transform of the current-free state prescribed
    !! (benchmarks/label-qa/input.landreman-paul-qa-iota). chi = s
    !! (normalized toroidal flux), F(chi) = fs (-PHIEDGE chi/(2 pi)),
    !! G'(chi) = gs (-iota(chi) PHIEDGE/(2 pi)). Planes over one field period
    !! carry O-grids fitted to the boundary with their cores on the boundary
    !! centroid; labels start from the geometric guess. Writes the layout of
    !! benchmarks/label-qa/reference (axis.txt, surface_s*.txt, summary.json)
    !! to the output directory for benchmarks/label-qa/compare.py.
    !!
    !! Arguments: deck, output directory, cells per block side (2), planes
    !! per period (21), fs (1), gs (1), continuation steps in the shaping (1: none),
    !! added shear (0), degree (2), iota stages (1: none; k: iota 2^(k-1) down to 1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_import, only: imported_equilibrium_t, IMPORT_OK
    use milca_import_vmec, only: import_vmec
    use milca_label_ale_3d
    use milca_eulerian_clebsch, only: clebsch_nfp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nsamples = 256, nout = 64, nray = 128
    type(imported_equilibrium_t) :: eq
    type(label3d_t) :: s
    character(512) :: deck = '../benchmarks/label-qa/input.landreman-paul-qa-iota', &
                      outdir = 'label-qa-out', arg
    integer :: cells = 2, nt = 21, k, i, st, j, u
    real(dp) :: fs = 1, gs = 1, phiedge, t0, t1, w0, energy, shaping = 1
    real(dp), allocatable :: steps_lambda(:)
    integer :: istep, nsteps = 1
    real(dp) :: extra_shear = 0, iota_scale = 1
    integer :: degree = 2, iota_stages = 1, ist
    real(dp) :: newton_ridge = 1e-8_dp
    integer :: gauge_k = 0
    real(dp), allocatable :: iota_c(:), samples(:, :, :), centre(:, :)
    real(dp) :: axis(2), rz(2, nray, 3), levels(3) = [0.25_dp, 0.5_dp, 1.0_dp], phi
    character(4), parameter :: tags(3) = ['0.25', '0.50', '1.00']

    if (command_argument_count() >= 1) call get_command_argument(1, deck)
    if (command_argument_count() >= 2) call get_command_argument(2, outdir)
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) cells
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) nt
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(5, arg)
        read (arg, *) fs
        call get_command_argument(6, arg)
        read (arg, *) gs
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) nsteps
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, arg)
        read (arg, *) extra_shear
    end if
    if (command_argument_count() >= 9) then
        call get_command_argument(9, arg)
        read (arg, *) degree
    end if
    if (command_argument_count() >= 10) then
        call get_command_argument(10, arg)
        read (arg, *) iota_stages
    end if
    if (command_argument_count() >= 11) then
        call get_command_argument(11, arg)
        read (arg, *) newton_ridge
    end if
    if (command_argument_count() >= 12) then
        call get_command_argument(12, arg)
        read (arg, *) gauge_k
    end if
    allocate (steps_lambda(nsteps))
    steps_lambda = [(real(istep, dp)/nsteps, istep=1, nsteps)]
    if (nsteps > 1) steps_lambda = [(real(istep - 1, dp)/(nsteps - 1), istep=1, nsteps)]
    call wall_time(t0)
    call import_vmec(trim(deck), eq)
    if (eq%status /= IMPORT_OK) then
        print '(2a)', 'import: ', trim(eq%reason)
        error stop 'deck'
    end if
    phiedge = eq%phiedge
    ! the planes span one field period
    clebsch_nfp = eq%nfp
    iota_c = eq%iota%coefficients
    print '(a,i0,a,f10.6,a,7f10.6)', 'NFP ', eq%nfp, '  PHIEDGE', phiedge, '  iota', iota_c

    ! wall curves of the planes phi_k = 2 pi (k - 1)/nt, counter-clockwise,
    ! first sample on the outboard side at the centroid height
    allocate (samples(2, nsamples, nt), centre(2, nt))
    ! continuation in the shaping: boundary = axisymmetric ellipse (m <= 1,
    ! n = 0) + lambda (rest), lambda over steps_lambda
    ! continuation: shaping at the largest rotational transform (iota scaled
    ! by 2^(iota_stages - 1), B_phi/B_pol smaller), then iota down to the
    ! deck's; g restarts from 0 at every stage (no gauge leftovers)
    do ist = 1, nsteps + iota_stages - 1
        istep = min(ist, nsteps)
        shaping = steps_lambda(istep)
        iota_scale = 2.0_dp**max(0, iota_stages - 1 - max(0, ist - nsteps))
        do k = 1, nt
            phi = 2*pi*(k - 1)/(nt*eq%nfp)
            call plane_curve(phi, samples(:, :, k), centre(:, k))
        end do
        if (ist == 1) then
            s%degree = degree
            s%newton_ridge = newton_ridge
            s%gauge_k = gauge_k
            call label3d_init(s, samples, centre, cells, st)
            if (st /= 0) error stop 'meshes'
            call label3d_initial_labels(s)
        else if (ist <= nsteps) then
            call label3d_set_walls(s, samples, st)
            if (st /= 0) error stop 'continuation meshes'
        end if
        call label3d_frames(s)
        if (allocated(s%g)) s%g = 0
        call label3d_solve_g(s, profiles, st)
        if (st /= 0) error stop 'initial g'
        if (ist == 1) w0 = label3d_energy(s, profiles)
        s%verbose = .true.
        call label3d_solve(s, profiles, 1e-9_dp, 0.05_dp, 6, 16, st)
        print '(a,f6.3,a,f6.2,a,i0,a,i0,a,es14.7)', 'shaping', shaping, '  iota scale', &
            iota_scale, '  status ', st, '  Newton ', s%newton_total, '  W', &
            label3d_energy(s, profiles)
        if (st /= 0) exit
    end do
    energy = label3d_energy(s, profiles)
    call wall_time(t1)
    print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,es14.7,a,es14.7,a,f8.2,a)', 'cells ', cells, '  planes ', nt, &
        '  status ', st, '  passes ', s%passes, '  Newton ', s%newton_total, '  W', energy, &
        '  (start', w0, ')  ', t1 - t0, ' s'

    ! output on the reference planes of one period
    call execute_command_line('mkdir -p '//trim(outdir))
    open (newunit=u, file=trim(outdir)//'/axis.txt', status='replace')
    write (u, '(a)') '# phi R Z'
    do j = 0, nout - 1
        phi = 2*pi*j/(nout*eq%nfp)
        call label3d_section_at(s, phi*eq%nfp/(2*pi), levels, nray, axis, rz, st)
        if (st /= 0) error stop 'section output'
        write (u, '(3es24.16)') phi, axis
        do i = 1, 3
            call append_surface(i, phi, rz(:, :, i), j == 0)
        end do
    end do
    close (u)
    open (newunit=u, file=trim(outdir)//'/summary.json', status='replace')
    write (u, '(a,es24.16,a,i0,a,i0,a,i0,a,f10.3,a)') '{"energy_T2m3": ', energy, &
        ', "cells": ', cells, ', "planes": ', nt, ', "status": ', st, ', "time_s": ', &
        t1 - t0, '}'
    close (u)
    print '(2a)', 'wrote ', trim(outdir)

contains

    subroutine wall_time(t)
        real(dp), intent(out) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end subroutine wall_time

    subroutine append_surface(i, phi, pts, first)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi, pts(:, :)
        logical, intent(in) :: first
        integer :: v, q
        if (first) then
            open (newunit=v, file=trim(outdir)//'/surface_s'//tags(i)//'.txt', status='replace')
            write (v, '(a)') '# phi theta R Z |B| (|B| not computed: nan)'
        else
            open (newunit=v, file=trim(outdir)//'/surface_s'//tags(i)//'.txt', &
                  position='append', status='old')
        end if
        do q = 1, size(pts, 2)
            write (v, '(4es24.16,a)') phi, 2*pi*(q - 1)/size(pts, 2), pts(:, q), '  nan'
        end do
        close (v)
    end subroutine append_surface

    subroutine plane_curve(phi, curve, c)
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: curve(:, :), c(2)
        real(dp), allocatable :: raw(:, :)
        real(dp) :: th, area, a, best, w
        integer :: i, m, n, first, ns
        ns = size(curve, 2)
        allocate (raw(2, ns))
        do i = 1, ns
            th = 2*pi*(i - 1)/ns
            raw(:, i) = 0
            do m = lbound(eq%rbc, 2), ubound(eq%rbc, 2)
                do n = lbound(eq%rbc, 1), ubound(eq%rbc, 1)
                    a = m*th - n*eq%nfp*phi
                    w = 1
                    if (n /= 0 .or. m > 1) w = shaping
                    raw(1, i) = raw(1, i) + w*(eq%rbc(n, m)*cos(a) + eq%rbs(n, m)*sin(a))
                    raw(2, i) = raw(2, i) + w*(eq%zbc(n, m)*cos(a) + eq%zbs(n, m)*sin(a))
                end do
            end do
        end do
        ! orientation and area centroid
        area = 0
        c = 0
        do i = 1, ns
            associate (p => raw(:, i), q => raw(:, mod(i, ns) + 1))
                a = p(1)*q(2) - q(1)*p(2)
                area = area + a/2
                c = c + (p + q)*a/6
            end associate
        end do
        c = c/area
        if (area < 0) raw = raw(:, ns:1:-1)
        ! start on the outboard side nearest the centroid height
        best = huge(1.0_dp)
        first = 1
        do i = 1, ns
            if (raw(1, i) > c(1) .and. abs(raw(2, i) - c(2)) < best) then
                best = abs(raw(2, i) - c(2))
                first = i
            end if
        end do
        curve = cshift(raw, first - 1, dim=2)
    end subroutine plane_curve

    subroutine profiles(chi, f, g, p)
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: f(3), g(3), p(3)
        real(dp) :: io(3)
        integer :: k
        io = 0
        do k = 1, size(iota_c)
            io(1) = io(1) + iota_c(k)*chi**(k - 1)
            if (k >= 2) io(2) = io(2) + (k - 1)*iota_c(k)*chi**(k - 2)
            if (k >= 3) io(3) = io(3) + (k - 1)*(k - 2)*iota_c(k)*chi**(k - 3)
        end do
        ! optional added shear (diagnosis): iota + extra_shear chi
        io(1) = io(1) + extra_shear*chi
        io(2) = io(2) + extra_shear
        f = [fs*(-phiedge/(2*pi)), 0.0_dp, 0.0_dp]
        g = gs*(-phiedge/(2*pi))*io*iota_scale
        p = 0
    end subroutine profiles

end program test_label_ale_qa
