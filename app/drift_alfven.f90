program drift_alfven_main
    !! Electromagnetic drift-Alfven turbulence with one parallel wave number on
    !! a doubly periodic box (docs/drift-turbulence.md, stage DA, homogeneous
    !! step): time series of the ledger, time averages with standard errors and
    !! the fluctuation spectrum of the saturated state. Units: lengths rho_s,
    !! time L_n/c_s, fluctuations times L_n/rho_s.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    use milca_drift_spectral, only: drift_to_spectral
    use milca_drift_alfven, only: da_parameters_t, da_state_t, da_ledger_t, da_initialize, &
        da_step, da_ledger
    use milca_turbulence_statistics, only: turbulence_accumulator_t, turbulence_spectrum_t, &
        turbulence_scalar_t, turbulence_accumulator_init, turbulence_accumulate, &
        turbulence_spectrum, turbulence_series_scalar
    use milca_result_io, only: run_t, result_file_t
    use milca_turbulence_output, only: turbulence_units_t, write_turbulence_spectrum, &
        write_turbulence_scalar
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(run_t) :: run
    type(result_file_t) :: result
    type(da_parameters_t) :: par
    type(da_state_t) :: state
    type(da_ledger_t) :: ledger
    type(turbulence_accumulator_t) :: acc
    type(turbulence_spectrum_t) :: spectrum
    type(turbulence_units_t) :: units
    type(turbulence_scalar_t) :: flux, sink, energy, enstrophy, zonal, diffusivity, &
        hyper_e, hyper_w
    !> beta = (beta_e/2)(L_par/L_n)^2, mu = (m_e/m_i)(L_par/L_n)^2,
    !> c = 0.51 (nu_e L_n/c_s) mu, kpar in 1/L_par, kappa the drive.
    real(dp) :: beta = 1.0_dp, mu = 1.0_dp, c = 1.0_dp, kpar = 1.0_dp, kappa = 1.0_dp, &
        d = 1.0e-4_dp
    integer :: order = 2
    integer :: nx = 128, ny = 128
    !> Box lx = ly = 2 pi/k0.
    real(dp) :: k0 = 0.15_dp
    real(dp) :: dt = 0.05_dp, time_end = 400.0_dp, average_from = 200.0_dp
    !> Steps between samples of the ledger series and of the spectrum input.
    integer :: sample_steps = 5
    !> Samples per Welch segment and largest retained mode number per direction.
    integer :: segment = 256, modes = 16
    real(dp) :: amplitude = 1.0e-2_dp
    integer :: seed = 1
    !> Size of rho_s, L_n/c_s and T_e/e in cm, s and statvolt.
    real(dp) :: rho_s_cm = 1.0_dp, ln_over_cs_s = 1.0_dp, te_over_e_statvolt = 1.0_dp
    namelist /drift_alfven/ beta, mu, c, kpar, kappa, d, order, nx, ny, k0, dt, time_end, &
        average_from, sample_steps, segment, modes, amplitude, seed, rho_s_cm, &
        ln_over_cs_s, te_over_e_statvolt
    real(dp), allocatable :: series(:, :), kx(:), ky(:)
    complex(dp), allocatable :: phi(:, :), dens(:, :)
    integer, allocatable :: ix(:)
    real(dp) :: box, sample_dt, balance(2)
    integer :: unit, ios, steps, step, samples, first, k, i, j
    character(len=256) :: message
    character(len=*), parameter :: names(8) = [character(len=13) :: 'energy', 'enstrophy', &
        'flux', 'd_c', 'd_e', 'd_w', 'kinetic', 'zonal_kinetic']

    call run%start('drift_alfven')
    ios = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=drift_alfven, iostat=ios, iomsg=message)
    call run%close_config(unit, 'drift_alfven', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: drift_alfven [--config FILE]'
    if (.not. (k0 > 0 .and. dt > 0 .and. time_end > average_from .and. average_from >= 0)) &
        error stop 'drift_alfven: need k0, dt > 0 and time_end > average_from >= 0'
    if (sample_steps < 1 .or. modes < 1 .or. 3*modes > min(nx, ny) - 1) &
        error stop 'drift_alfven: need sample_steps >= 1 and 1 <= modes <= (n - 1)/3'
    call run%open_resolved(unit)
    write (unit, nml=drift_alfven)
    call run%close_resolved(unit)

    par = da_parameters_t(beta, mu, c, kpar, kappa, d, order)
    box = 2*pi/k0
    call da_initialize(state, par, nx, ny, box, box, dt)
    call seed_noise()
    steps = nint(time_end/dt)
    samples = steps/sample_steps
    sample_dt = dt*sample_steps
    first = int(average_from/sample_dt) + 1
    if (samples - first + 1 < segment) &
        error stop 'drift_alfven: averaging window shorter than one segment'
    allocate (series(samples, 9))
    kx = [(k0*i, i=-modes, modes)]
    ky = [(k0*j, j=0, modes)]
    ix = [(modulo(i, nx) + 1, i=-modes, modes)]
    allocate (phi(2*modes + 1, modes + 1), dens(2*modes + 1, modes + 1))
    call turbulence_accumulator_init(acc, kx, ky, sample_dt, segment)

    do step = 1, steps
        call da_step(state)
        if (mod(step, sample_steps) /= 0) cycle
        k = step/sample_steps
        call da_ledger(state, ledger)
        if (.not. ledger%energy < huge(1.0_dp)) then
            write (*, '(a,es12.5)') 'non-finite energy at time ', state%time
            error stop 'drift_alfven: run diverged'
        end if
        series(k, :) = [state%time, ledger%energy, ledger%enstrophy, ledger%flux, &
                        ledger%d_c, ledger%d_e, ledger%d_w, ledger%kinetic, &
                        ledger%zonal_kinetic]
        if (k < first) cycle
        do j = 1, modes + 1
            do i = 1, 2*modes + 1
                phi(i, j) = 0
                if (state%grid%k2(ix(i), j) > 0) &
                    phi(i, j) = -state%zeta(ix(i), j)/state%grid%k2(ix(i), j)
                dens(i, j) = state%n(ix(i), j)
            end do
        end do
        call turbulence_accumulate(acc, phi, dens)
    end do
    call turbulence_spectrum(acc, spectrum)

    flux = turbulence_series_scalar(series(first:, 4), sample_dt)
    sink = turbulence_series_scalar(series(first:, 5), sample_dt)
    hyper_e = turbulence_series_scalar(series(first:, 6), sample_dt)
    hyper_w = turbulence_series_scalar(series(first:, 7), sample_dt)
    energy = turbulence_series_scalar(series(first:, 2), sample_dt)
    enstrophy = turbulence_series_scalar(series(first:, 3), sample_dt)
    zonal = turbulence_series_scalar(series(first:, 9)/series(first:, 8), sample_dt)
    ! Gamma_n = kappa <n v_x> and the gradient of n0 is -kappa.
    diffusivity = flux
    diffusivity%value = flux%value/kappa**2
    diffusivity%error = flux%error/kappa**2
    ! Mean rate of change over the averaging window against the mean ledger.
    balance(1) = (series(samples, 2) - series(first, 2))/(series(samples, 1) - series(first, 1)) &
                 - (flux%value - sink%value - hyper_e%value)
    balance(2) = (series(samples, 3) - series(first, 3))/(series(samples, 1) - series(first, 1)) &
                 - (flux%value - hyper_w%value)

    units%length_cm = rho_s_cm
    units%time_s = ln_over_cs_s
    units%potential_statvolt = te_over_e_statvolt
    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('time', series(:, 1), 'L_n/c_s', 'time', dims='time')
    do k = 1, 8
        call result%write(trim(names(k)), series(:, k + 1), 'T_e^2/e^2 box average', &
            'ledger term '//trim(names(k)), dims='time', status='candidate')
    end do
    call write_turbulence_scalar(result, 'Gamma_n', flux, '1', &
        'particle flux kappa <n v_x>, time average')
    call write_turbulence_scalar(result, 'D_C', sink, '1', 'resistive sink, time average')
    call write_turbulence_scalar(result, 'D_E', hyper_e, '1', &
        'hyperviscous energy sink, time average')
    call write_turbulence_scalar(result, 'D_W', hyper_w, '1', &
        'hyperviscous enstrophy sink, time average')
    call write_turbulence_scalar(result, 'E', energy, '1', 'energy, time average')
    call write_turbulence_scalar(result, 'W', enstrophy, '1', &
        'potential enstrophy, time average')
    call write_turbulence_scalar(result, 'zonal_fraction', zonal, '1', &
        'zonal share of the kinetic energy, time average')
    call write_turbulence_scalar(result, 'D_n', diffusivity, 'rho_s^2 c_s/L_n', &
        'particle diffusivity <n v_x>/kappa')
    call result%write('energy_balance_residual', balance(1), '1', &
        'mean dE/dt minus mean (Gamma_n - D_C - D_E)', status='indicator')
    call result%write('enstrophy_balance_residual', balance(2), '1', &
        'mean dW/dt minus mean (Gamma_n - D_W)', status='indicator')
    call result%attribute('derivative_status', &
        'unavailable: time averages of a chaotic state')
    call write_turbulence_spectrum(result, spectrum, units, &
        'drift-Alfven, one parallel wave number', 'nonlinear')
    call result%close()

    write (*, '(a,es12.5,a,es9.2)') 'Gamma_n = ', flux%value, ' +- ', flux%error
    write (*, '(a,es12.5,a,es9.2)') 'D_C     = ', sink%value, ' +- ', sink%error
    write (*, '(a,es12.5,a,es9.2)') 'E       = ', energy%value, ' +- ', energy%error
    write (*, '(a,es12.5,a,es9.2)') 'W       = ', enstrophy%value, ' +- ', enstrophy%error
    write (*, '(a,es12.5,a,es9.2)') 'zonal   = ', zonal%value, ' +- ', zonal%error
    write (*, '(a,2es11.3)') 'balance residuals (E, W): ', balance
    write (*, '(a,es12.5,a,es9.2)') 'D_E     = ', hyper_e%value, ' +- ', hyper_e%error
    write (*, '(a,es12.5,a,es9.2)') 'tau_eul = ', spectrum%tau_euler%value, ' +- ', &
        spectrum%tau_euler%error
    write (*, '(a,es12.5,a,es9.2)') 'tau_mod = ', spectrum%tau_mode%value, ' +- ', &
        spectrum%tau_mode%error
    write (*, '(a,es12.5,a,es9.2)') 'ell_x   = ', spectrum%ell_x%value, ' +- ', &
        spectrum%ell_x%error
    write (*, '(a,es12.5,a,es9.2)') 'ell_y   = ', spectrum%ell_y%value, ' +- ', &
        spectrum%ell_y%error
    write (*, '(a,es12.5,a,es9.2)') 'kubo    = ', spectrum%kubo%value, ' +- ', &
        spectrum%kubo%error
contains
    !> Band-limited Gaussian noise of standard deviation `amplitude` per grid
    !> point in zeta and n.
    subroutine seed_noise()
        type(rng_t) :: rng
        type(fortnum_status_t) :: status
        real(dp), allocatable :: field(:, :)
        integer :: i, j, pass

        call rng_seed(rng, int(seed, int64), status)
        allocate (field(nx, ny))
        do pass = 1, 2
            do j = 1, ny
                do i = 1, nx
                    call rng_normal(rng, field(i, j))
                end do
            end do
            field = amplitude*(field - sum(field)/size(field))
            if (pass == 1) call drift_to_spectral(state%grid, state%work, field, state%zeta)
            if (pass == 2) call drift_to_spectral(state%grid, state%work, field, state%n)
        end do
    end subroutine
end program drift_alfven_main
