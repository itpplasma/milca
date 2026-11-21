program test_frozen_map_spline
    ! Frozen-in energy map on the multipatch spline space (milca_frozen_map_spline):
    ! energy gradient and Hessian against central differences, including the
    ! relabelling potential f, and the fixed-boundary Newton on a circular
    ! tokamak to the requested residual.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_energy, only: spline_energy_hessian_csr
    use fortfem_element_to_csr, only: csr_matrix_t
    use milca_frozen_map_spline, only: frozen_spline_map_t, frozen_spline_init, &
        frozen_spline_energy, frozen_spline_solve, FROZEN_SPLINE_OK
    use milca_ideal_seed_3d, only: toroidal_seed_t
    implicit none
    real(dp), parameter :: mu0 = 4.0e-7_dp*acos(-1.0_dp)
    type(frozen_spline_map_t) :: m
    type(toroidal_seed_t) :: seed
    type(csr_matrix_t) :: h
    real(dp) :: boundary(4, 1, 0:1), e0, ep, em, step, err, residual
    real(dp), allocatable :: g(:), gp(:), gm(:), v(:), hv(:), c0(:, :)
    integer :: st, i, p, n, iterations

    seed%flux = 5.0_dp
    seed%centre = [0.0_dp, 0.0_dp]; seed%axes = [1.0_dp, 1.0_dp]
    seed%profile%iota = [0.9_dp, -0.4_dp]
    seed%profile%pressure = mu0*[2.0e4_dp, -2.0e4_dp]
    boundary = 0
    boundary(1, 1, 0) = 3; boundary(1, 1, 1) = 1; boundary(4, 1, 1) = 1.3_dp
    call frozen_spline_init(m, 3, 2, 0.5_dp, seed, boundary, mu0, st)
    if (st /= FROZEN_SPLINE_OK) error stop 'init'
    m%kernel%f_gauge = 0; m%kernel%slide_gauge = 0
    n = size(m%coef)
    allocate (v(n), hv(n))
    do i = 1, n
        v(i) = cos(0.71_dp*i)
        m%coef(mod(i - 1, 3) + 1, (i - 1)/3 + 1) = m%coef(mod(i - 1, 3) + 1, (i - 1)/3 + 1) &
            + 1.0e-3_dp*sin(1.3_dp*i)
    end do
    call frozen_spline_energy(m, e0, st, g)
    if (st /= 0) error stop 'base energy'
    c0 = m%coef
    step = 1.0e-6_dp
    m%coef = c0 + step*reshape(v, shape(c0))
    call frozen_spline_energy(m, ep, st, gp)
    m%coef = c0 - step*reshape(v, shape(c0))
    call frozen_spline_energy(m, em, st, gm)
    m%coef = c0
    err = abs((ep - em)/(2*step) - dot_product(g, v))/sum(abs(g*v))
    print '(a,es9.2)', 'gradient against differences', err
    if (err > 1.0e-7_dp) error stop 'gradient mismatch'
    call spline_energy_hessian_csr(m%space, m%quad, m%kernel, m%coef, h, st)
    if (st /= 0) error stop 'hessian'
    hv = 0
    do i = 1, n
        do p = h%rowptr(i), h%rowptr(i + 1) - 1
            hv(i) = hv(i) + h%val(p)*v(h%colidx(p))
        end do
    end do
    err = maxval(abs((gp - gm)/(2*step) - hv))/maxval(abs(hv))
    print '(a,es9.2)', 'Hessian against differences', err
    if (err > 1.0e-6_dp) error stop 'Hessian mismatch'
    err = maxval(abs(hv(3::3)))/maxval(abs(hv))
    print '(a,es9.2)', 'f block share of H v', err
    if (.not. err > 1.0e-6_dp) error stop 'the relabelling potential does not enter'

    call frozen_spline_init(m, 3, 3, 0.5_dp, seed, boundary, mu0, st)
    if (st /= FROZEN_SPLINE_OK) error stop 'init 2'
    m%kernel%f_gauge = 1.0e-2_dp; m%kernel%slide_gauge = 1.0e-2_dp
    call frozen_spline_solve(m, 1.0e-10_dp, 60, iterations, residual, st)
    print '(a,i0,a,i0,a,es9.2)', 'Newton status ', st, ', iterations ', iterations, &
        ', residual ', residual
    if (st /= FROZEN_SPLINE_OK) error stop 'fixed-boundary Newton'
    print '(a)', 'frozen map on splines: derivatives and Newton pass'
end program test_frozen_map_spline
