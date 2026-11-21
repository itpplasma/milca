module vacuum_torus_fields
    !! Exact vacuum field of gate G1 in Cartesian x: A = a_z(x) e_z with
    !! a_z = -b0 r0 ln R + 2 eps x y, so B = b0 r0 e_phi/R + eps grad(x^2 - y^2),
    !! curl-free and divergence-free, with toroidal and n = 2 parts.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_map
    implicit none
    private

    real(dp), parameter, public :: b0 = 1.0_dp, r0 = 3.0_dp, eps = 0.1_dp
    type(spline_space_t), public :: space
    type(curve_chart_t), public :: chart

    public :: full, potential, field_b, metric2, det3, jacobian3

contains

    subroutine full(block, xi, t, x, j3)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: x(3), j3(3, 3)
        real(dp) :: jet(2, 6), jc(3, 3)
        integer :: status

        call spline_geometry_at(space, block, xi, jet, status)
        call curve_chart_map(chart, jet(:, 1), t, x, jc)
        j3(:, 1:2) = matmul(jc(:, 1:2), jet(:, 2:3))
        j3(:, 3) = jc(:, 3)
    end subroutine full

    pure function exact_b(x) result(b)
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), rr

        rr = x(1)**2 + x(2)**2
        b = b0*r0*[-x(2), x(1), 0.0_dp]/rr + eps*[2*x(1), -2*x(2), 0.0_dp]
    end function exact_b

    subroutine potential(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: x(3), j3(3, 3), az

        call full(block, xi, t, x, j3)
        az = -b0*r0*0.5_dp*log(x(1)**2 + x(2)**2) + 2*eps*x(1)*x(2)
        values(1:3) = j3(3, :)*az
    end subroutine potential

    subroutine field_b(block, xi, t, x, b)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: x(3), b(3)
        real(dp) :: j3(3, 3)

        call full(block, xi, t, x, j3)
        b = exact_b(x)
    end subroutine field_b

    subroutine jacobian3(block, xi, t, jacobian)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: jacobian(:, :)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        jacobian(1:3, 1:3) = j3
    end subroutine jacobian3

    subroutine metric2(block, xi, t, weight)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: weight(:, :)
        real(dp) :: x(3), j3(3, 3)

        call full(block, xi, t, x, j3)
        weight(1:3, 1:3) = matmul(transpose(j3), j3)/abs(det3(j3))
    end subroutine metric2

    pure real(dp) function det3(j)
        real(dp), intent(in) :: j(3, 3)

        det3 = j(1, 1)*(j(2, 2)*j(3, 3) - j(3, 2)*j(2, 3)) - &
               j(1, 2)*(j(2, 1)*j(3, 3) - j(3, 1)*j(2, 3)) + &
               j(1, 3)*(j(2, 1)*j(3, 2) - j(3, 1)*j(2, 2))
    end function det3

end module vacuum_torus_fields

program test_vacuum_torus_3d
    !! Gate G1 of docs/eulerian-3d.md: the vacuum field in a torus from the
    !! ungauged vector potential with n x A fixed on the wall. The solution
    !! reproduces the exact field to the approximation error, converging at
    !! order p in L2, and has a lower energy than the interpolant with the
    !! same boundary data.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_options_t, spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build
    use fortfem_curve_chart, only: curve_chart_cylindrical
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_project, spline_derham_boundary_dofs, spline_derham_tree, spline_derham_evaluate, &
        spline_derham_interpolation
    use fortsparse_csr, only: csr_t, csr_matvec
    use fortfem_element_to_csr, only: csr_matrix_t
    use milca_eulerian_vacuum, only: eu_curl_curl, eu_constrained_solve, eu_ams_solve
    use vacuum_torus_fields, only: space, chart, r0, potential, field_b, metric2, full, det3, &
        jacobian3
    implicit none
    integer, parameter :: modes = 2, degree = 2, levels = 3
    integer :: level, cells(levels) = [2, 4, 8], failures
    real(dp) :: error(levels), error_interp, energy, energy_interp, rate, error_ams
    integer :: iterations(levels)

    failures = 0
    call curve_chart_cylindrical(chart, 1)
    do level = 1, levels
        call solve(cells(level), error(level), error_interp, energy, energy_interp, &
                   error_ams, iterations(level))
        print '(a,i3,a,2es10.3,a,2es22.14)', 'cells', cells(level), &
            ' L2 error of B (Galerkin, interpolant)', error(level), error_interp, &
            ' energy', energy, energy_interp
        print '(a,i0,a,es10.3)', 'ungauged PCG-AMS iterations ', iterations(level), &
            ' L2 error of B ', error_ams
        call check(abs(error_ams - error(level)) < 1e-3_dp*error(level), &
                   'ungauged AMS solve gives the field of the gauged direct solve')
        call check(energy <= energy_interp*(1 + 1e-12_dp), &
                   'Galerkin energy below the interpolant energy')
    end do
    do level = 2, levels
        rate = log(error(level - 1)/error(level))/log(2.0_dp)
        print '(a,f6.2)', 'rate', rate
        call check(rate > degree - 0.3_dp, 'L2 error of B converges at order p')
    end do
    call check(error(levels) < 1e-3_dp, 'relative L2 error of B below 1e-3')
    if (failures > 0) error stop 'vacuum torus gate G1 failed'

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

    subroutine solve(n, err, e_p, e_h, e_i, err_ams, its)
        integer, intent(in) :: n
        real(dp), intent(out) :: err, e_p, e_h, e_i, err_ams
        integer, intent(out) :: its
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        type(spline_quadrature_t) :: quad
        type(spline_derham_t) :: derham
        type(csr_t) :: k
        type(csr_matrix_t) :: interpolation
        integer, allocatable :: nodes(:), edges(:), tree(:), fixed(:)
        real(dp), allocatable :: a_i(:), a_h(:), kv(:), a_ams(:)
        logical, allocatable :: tree_mask(:)
        integer :: st, i

        call layout_ogrid([r0, 0.0_dp], 1.0_dp, 0.5_dp, layout, st)
        if (st /= 0) error stop 'layout'
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = n
        call spline_space_build(layout, options, space, st)
        if (st /= 0) error stop 'space'
        call spline_quadrature_build(space, degree + 2, quad, st)
        if (st /= 0) error stop 'quadrature'
        call spline_derham_build(space, derham, st, modes)
        if (st /= 0) error stop 'derham'
        allocate (a_i(derham%dimension(1)), a_h(derham%dimension(1)))
        call spline_derham_project(derham, 1, potential, a_i, st, 12)
        if (st /= 0) error stop 'projection'
        call eu_curl_curl(derham, quad, metric2, k, st, 2*modes + 4)
        if (st /= 0) error stop 'curl-curl'
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 0, nodes)
        call spline_derham_boundary_dofs(derham, LAYOUT_WALL, 1, edges)
        call spline_derham_tree(derham, tree_mask, roots=nodes)
        tree = pack([(i, i=1, size(tree_mask))], tree_mask)
        fixed = [edges, tree]
        call eu_constrained_solve(k, 0*a_i, fixed, [a_i(edges), 0*a_i(tree)], a_h, st)
        if (st /= 0) error stop 'solve'
        call check(size(tree) == derham%dimension(0) - size(nodes), &
                   'one gauge edge per interior node')
        allocate (kv(size(a_h)))
        call csr_matvec(k, a_h, kv)
        e_h = dot_product(a_h, kv)
        call csr_matvec(k, a_i, kv)
        e_i = dot_product(a_i, kv)
        err = l2_error(derham, quad, curl(derham, a_h))
        e_p = l2_error(derham, quad, curl(derham, a_i))
        call spline_derham_interpolation(derham, jacobian3, interpolation, st, 1e-4_dp)
        if (st /= 0) error stop 'interpolation'
        allocate (a_ams(size(a_h)))
        call eu_ams_solve(k, derham%grad, interpolation, derham%dimension(0), edges, &
                          a_i(edges), nodes, a_ams, 1e-10_dp, its, st)
        if (st /= 0) print *, 'AMS status', st, ' iterations', its
        err_ams = l2_error(derham, quad, curl(derham, a_ams))
    end subroutine solve

    function curl(derham, a) result(b)
        type(spline_derham_t), intent(in) :: derham
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: b(:)
        integer :: i, p

        allocate (b(derham%dimension(2)))
        do i = 1, derham%curl%n
            b(i) = 0
            do p = derham%curl%rowptr(i), derham%curl%rowptr(i + 1) - 1
                b(i) = b(i) + derham%curl%val(p)*a(derham%curl%colidx(p))
            end do
        end do
    end function curl

    real(dp) function l2_error(derham, quad, b2) result(rel)
        !! || B_h - B || / || B || by the section quadrature times a 16-point
        !! uniform rule in t (exact for the trigonometric degree present).
        !! B_h = J b/det J from the 2-form proxy b, and dV = w 2 pi R dt
        !! because quad%weight holds the section Jacobian.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: b2(:)
        integer, parameter :: nt = 16
        real(dp) :: x(3), j3(3, 3), b(3), bh(3), proxy(3), t, dv, num, den
        integer :: q, m

        num = 0
        den = 0
        do q = 1, quad%point_count
            do m = 1, nt
                t = (m - 1)/real(nt, dp)
                call field_b(quad%block(q), quad%xi(:, q), t, x, b)
                call full(quad%block(q), quad%xi(:, q), t, x, j3)
                call spline_derham_evaluate(derham, 2, b2, quad%block(q), quad%xi(:, q), &
                                            t, proxy)
                bh = matmul(j3, proxy)/det3(j3)
                dv = quad%weight(q)*2*acos(-1.0_dp)*hypot(x(1), x(2))/nt
                num = num + sum((bh - b)**2)*dv
                den = den + sum(b**2)*dv
            end do
        end do
        rel = sqrt(num/den)
    end function l2_error

end program test_vacuum_torus_3d
