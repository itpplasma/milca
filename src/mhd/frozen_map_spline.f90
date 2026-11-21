module milca_frozen_map_spline
    ! Frozen-in energy map (docs/ideal-map-3d.md) on FortFEM's multipatch spline
    ! spaces instead of the C0 P3 triangles: the unknowns are the chart
    ! coordinates q1, q2 of the map and the relabelling potential f of the
    ! frozen form, b -> b + grad f x grad chi (docs/free-boundary-ale.md (6)),
    ! all in one spline space over the reference plane s. Axisymmetric (q3 =
    ! phi, Q(:, 3) = e3), fixed boundary. The pointwise mathematics are the
    ! generated leaves of milca_ideal_map_leaves, as for the triangles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_BOUNDARY
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_boundary_dofs, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build
    use fortfem_spline_energy, only: pointwise_energy_t, spline_energy, &
        spline_energy_hessian_csr, SPLINE_HESSIAN_DENSE
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_solve_factored, sparse_direct_free
    use milca_ideal_map, only: ideal_frozen_t, hessian_packed_pair
    use milca_ideal_map_leaves, only: ideal_point_gradient_jets, ideal_point_hessian_jets, &
        ideal_point_field_jets, ideal_point_field_gradient_jets
    use milca_ideal_seed_3d, only: ideal_seed_boundary_coefficients
    use milca_equilibrium_chart, only: cylindrical_chart_t
    implicit none
    private
    public :: frozen_spline_kernel_t, frozen_spline_map_t
    public :: frozen_spline_init, frozen_spline_energy, frozen_spline_solve
    public :: FROZEN_SPLINE_OK, FROZEN_SPLINE_FOLDED, FROZEN_SPLINE_FAILED

    integer, parameter :: FROZEN_SPLINE_OK = 0, FROZEN_SPLINE_FOLDED = 1, &
        FROZEN_SPLINE_FAILED = 2
    integer, parameter :: NF = 3, NC = 9
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Leaf variables v = (q(1:3), Q(k,i) at 3 + k + 3(i - 1)) that are jet
    ! components c = field + NF*d (d = 0 value, 1 d/ds1, 2 d/ds2) of q1, q2;
    ! for these v and c coincide (q3 = phi and the toroidal column are fixed).
    integer, parameter :: free_v(6) = [1, 2, 4, 5, 7, 8]

    type, extends(pointwise_energy_t) :: frozen_spline_kernel_t
        ! frozen data per quadrature point: b0 = (-d chi/ds2, d chi/ds1, b3),
        ! grad chi and p in solver units
        real(dp), allocatable :: b0(:, :), dchi(:, :), p(:)
        type(cylindrical_chart_t) :: chart
        real(dp) :: mu0 = 1, sigma = 1, measure = 2*pi
        ! Hessian-only gauge of the exact nulls: f(chi) (weight on the
        ! derivative of f across the surfaces) and tangential slides of q
        real(dp) :: f_gauge = 0, slide_gauge = 0
    contains
        procedure :: evaluate => kernel_evaluate
    end type

    type, extends(pointwise_energy_t) :: projection_kernel_t
        ! w = |u - target|^2/2 on the first two fields (L2 projection)
        real(dp), allocatable :: target(:, :)
    contains
        procedure :: evaluate => projection_evaluate
    end type

    type :: frozen_spline_map_t
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(frozen_spline_kernel_t) :: kernel
        real(dp), allocatable :: coef(:, :) ! (NF, dof)
        logical, allocatable :: fixed(:) ! (NF*dof)
    end type

contains

    subroutine frozen_spline_init(m, degree, cells, core_fraction, field, boundary, &
            mu0, status)
        ! O-grid over the unit reference disc, frozen data at the quadrature
        ! points, the map projected from the harmonic extension of the boundary
        ! (ideal_seed_boundary_coefficients), boundary positions and one f fixed.
        type(frozen_spline_map_t), intent(out) :: m
        integer, intent(in) :: degree, cells
        real(dp), intent(in) :: core_fraction, mu0
        class(ideal_frozen_t), intent(in) :: field
        real(dp), intent(in) :: boundary(:, :, 0:)
        integer, intent(out) :: status
        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        integer, allocatable :: dofs(:)
        real(dp) :: chi, grad(2), b3, p, values(2, 1), e
        real(dp), allocatable :: target(:, :)
        integer :: k, c, n
        status = FROZEN_SPLINE_FAILED
        call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, core_fraction, layout, k)
        if (k /= 0) return
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate (options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, m%space, k)
        if (k /= 0) return
        call spline_quadrature_build(m%space, degree + 1, m%quad, k)
        if (k /= 0) return
        n = m%quad%point_count
        allocate (m%kernel%b0(3, n), m%kernel%dchi(2, n), m%kernel%p(n), target(2, n))
        do k = 1, n
            call field%evaluate(m%quad%s(:, k), chi, grad, b3, p)
            m%kernel%b0(:, k) = [-grad(2), grad(1), b3]
            m%kernel%dchi(:, k) = grad
            m%kernel%p(k) = p
            call ideal_seed_boundary_coefficients(m%quad%s(:, k), boundary, values, &
                centre=[0.0_dp, 0.0_dp], axes=[1.0_dp, 1.0_dp])
            target(:, k) = values(:, 1)
        end do
        m%kernel%field_count = NF
        m%kernel%order = 1
        m%kernel%hessian_storage = SPLINE_HESSIAN_DENSE
        m%kernel%mu0 = mu0
        call project_map(m, target, status)
        if (status /= FROZEN_SPLINE_OK) return
        ! orientation fixed by the initial map, as in ideal_seed_3d_init (the
        ! boundary may run clockwise in theta); mixed signs still fail
        m%kernel%sigma = 1
        call frozen_spline_energy(m, e, k)
        if (k == FROZEN_SPLINE_FOLDED) then
            m%kernel%sigma = -1
            call frozen_spline_energy(m, e, k)
        end if
        status = FROZEN_SPLINE_FOLDED
        if (k /= 0) return
        status = FROZEN_SPLINE_FAILED
        allocate (m%fixed(NF*m%space%dof_count))
        m%fixed = .false.
        do c = 1, 2
            call spline_boundary_dofs(m%space, LAYOUT_BOUNDARY, dofs, c, NF)
            m%fixed(dofs) = .true.
        end do
        ! f + const is an exact null
        m%fixed(3) = .true.
        status = FROZEN_SPLINE_OK
    end subroutine

    subroutine project_map(m, target, status)
        ! L2 projection of target (R, Z) at the quadrature points onto q1, q2;
        ! f = 0. Mass matrix and load through the energy hooks.
        type(frozen_spline_map_t), intent(inout) :: m
        real(dp), intent(in) :: target(:, :)
        integer, intent(out) :: status
        type(projection_kernel_t) :: pk
        type(csr_matrix_t) :: mass
        real(dp), allocatable :: zero(:, :), load(:, :), x(:)
        real(dp) :: e
        integer :: st
        status = FROZEN_SPLINE_FAILED
        pk%field_count = NF; pk%order = 0; pk%target = target
        allocate (zero(NF, m%space%dof_count), load(NF, m%space%dof_count))
        zero = 0
        call spline_energy(m%space, m%quad, pk, zero, e, st, load)
        if (st /= 0) return
        call spline_energy_hessian_csr(m%space, m%quad, pk, zero, mass, st)
        if (st /= 0) return
        ! f has no projection term: its rows get a unit diagonal
        call solve_csr(mass, -reshape(load, [size(load)]), x, st, field3=.true.)
        if (st /= 0) return
        m%coef = reshape(x, [NF, m%space%dof_count])
        m%coef(3, :) = 0
        status = FROZEN_SPLINE_OK
    end subroutine

    subroutine kernel_evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(frozen_spline_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: b(3), db(NC), jac, gv(12), hp(78), hv(12, 12), gb(3), hbb(6), hvb(12, 3)
        real(dp) :: h(NC, NC), n(2), t(2), tq(2), dn
        integer :: k, pt, i, j, ci, cj
        status = 0
        do k = 1, npoint
            pt = first_point + k - 1
            q = [jet(1, k), jet(2, k), 0.0_dp]
            qd(:, 1) = [jet(4, k), jet(5, k), 0.0_dp]
            qd(:, 2) = [jet(7, k), jet(8, k), 0.0_dp]
            qd(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
            ! d b3 / d jet: through d f/ds1 (comp 6) and d f/ds2 (comp 9)
            db = 0
            db(6) = self%dchi(2, pt); db(9) = -self%dchi(1, pt)
            b = self%b0(:, pt)
            b(3) = b(3) + jet(6, k)*db(6) + jet(9, k)*db(9)
            call self%chart%jet(q, pos, d1, d2, d3)
            if (need < 2) then
                call ideal_point_gradient_jets(d1, d2, d3, qd, b, self%p(pt), self%mu0, &
                    self%sigma, w(k), jac, gv)
            else
                call ideal_point_hessian_jets(d1, d2, d3, qd, b, self%p(pt), self%mu0, &
                    self%sigma, w(k), jac, gv, hp)
            end if
            if (.not. jac > 0) status = FROZEN_SPLINE_FOLDED
            w(k) = self%measure*w(k)
            if (need < 1) cycle
            call ideal_point_field_gradient_jets(d1, d2, d3, qd, b, self%p(pt), self%mu0, &
                self%sigma, gb)
            gradient(:, k) = 0
            gradient(free_v, k) = gv(free_v)
            gradient(:, k) = self%measure*(gradient(:, k) + gb(3)*db)
            if (need < 2) cycle
            call ideal_point_field_jets(d1, d2, d3, qd, b, self%p(pt), self%mu0, &
                self%sigma, gb, hbb, hvb)
            do i = 1, 78
                call hessian_packed_pair(i, ci, cj)
                hv(ci, cj) = hp(i); hv(cj, ci) = hp(i)
            end do
            h = 0
            h(free_v, free_v) = hv(free_v, free_v)
            ! d2w/dv db3 times d b3/d jet, both orders
            do i = 1, size(free_v)
                j = free_v(i)
                h(j, :) = h(j, :) + hvb(j, 3)*db
                h(:, j) = h(:, j) + hvb(j, 3)*db
            end do
            h = h + hbb(6)*spread(db, 2, NC)*spread(db, 1, NC)
            ! Hessian-only gauges: f across the surfaces, q along them
            ! weights relative to the block each gauge complements
            dn = norm2(self%dchi(:, pt))
            if (dn > 0) then
                n = self%dchi(:, pt)/dn
                h([6, 9], [6, 9]) = h([6, 9], [6, 9]) + self%f_gauge*abs(hbb(6))*dn**2 &
                    *spread(n, 2, 2)*spread(n, 1, 2)
                t = [-n(2), n(1)]
                tq = qd(1:2, 1)*t(1) + qd(1:2, 2)*t(2)
                if (norm2(tq) > 0) tq = tq/norm2(tq)
                h(1:2, 1:2) = h(1:2, 1:2) + self%slide_gauge &
                    *(abs(h(4, 4)) + abs(h(5, 5)) + abs(h(7, 7)) + abs(h(8, 8))) &
                    *spread(tq, 2, 2)*spread(tq, 1, 2)
            end if
            hessian(:, k) = self%measure*reshape(h, [NC*NC])
        end do
    end subroutine

    subroutine projection_evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(projection_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: r(2)
        integer :: k, pt
        status = 0
        do k = 1, npoint
            pt = first_point + k - 1
            r = jet(1:2, k) - self%target(:, pt)
            w(k) = dot_product(r, r)/2
            if (need >= 1) then
                gradient(:, k) = 0
                gradient(1:2, k) = r
            end if
            if (need >= 2) then
                hessian(:, k) = 0
                hessian(1, k) = 1
                hessian(1 + NF + 1, k) = 1
            end if
        end do
    end subroutine

    subroutine frozen_spline_energy(m, energy, status, gradient)
        ! energy per field period; gradient flattened as coef(field, dof)
        type(frozen_spline_map_t), intent(in) :: m
        real(dp), intent(out) :: energy
        integer, intent(out) :: status
        real(dp), allocatable, intent(out), optional :: gradient(:)
        real(dp), allocatable :: g(:, :)
        if (present(gradient)) then
            allocate (g(NF, m%space%dof_count))
            call spline_energy(m%space, m%quad, m%kernel, m%coef, energy, status, g)
            gradient = reshape(g, [size(g)])
        else
            call spline_energy(m%space, m%quad, m%kernel, m%coef, energy, status)
        end if
        ! FortFEM: 1 invalid call, 2 a kernel status (here only the fold)
        select case (status)
        case (0)
            status = FROZEN_SPLINE_OK
        case (2)
            status = FROZEN_SPLINE_FOLDED
        case default
            status = FROZEN_SPLINE_FAILED
        end select
    end subroutine

    subroutine frozen_spline_solve(m, tolerance, max_iterations, iterations, residual, status)
        ! Newton on the free coefficients with a Levenberg shift on the
        ! diagonal and backtracking on the energy; stops when the free
        ! gradient falls below tolerance times its initial value.
        type(frozen_spline_map_t), intent(inout) :: m
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: iterations, status
        real(dp), intent(out) :: residual
        type(csr_matrix_t) :: h
        real(dp), allocatable :: g(:), du(:), c0(:, :)
        real(dp) :: e0, e1, g0, mu, alpha
        integer :: st, tries
        status = FROZEN_SPLINE_FAILED
        mu = 0
        call frozen_spline_energy(m, e0, st, g)
        if (st /= 0) then
            status = FROZEN_SPLINE_FOLDED; return
        end if
        where (m%fixed) g = 0
        g0 = norm2(g)
        residual = 1
        do iterations = 0, max_iterations
            residual = norm2(g)/max(g0, tiny(1.0_dp))
            if (residual <= tolerance) then
                status = FROZEN_SPLINE_OK; return
            end if
            if (iterations == max_iterations) exit
            call spline_energy_hessian_csr(m%space, m%quad, m%kernel, m%coef, h, st)
            if (st /= 0) exit
            c0 = m%coef
            do tries = 1, 20
                call solve_csr(h, -g, du, st, fixed=m%fixed, shift=mu)
                if (st /= 0) then
                    mu = max(4*mu, 1.0e-8_dp); cycle
                end if
                alpha = 1
                m%coef = c0 + alpha*reshape(du, shape(c0))
                call frozen_spline_energy(m, e1, st)
                if (st == 0 .and. e1 <= e0 + 1.0e-4_dp*alpha*dot_product(g, du)) exit
                m%coef = c0
                mu = max(4*mu, 1.0e-8_dp)
            end do
            if (tries > 20) exit
            mu = mu/4
            e0 = e1
            call frozen_spline_energy(m, e0, st, g)
            where (m%fixed) g = 0
        end do
    end subroutine

    subroutine solve_csr(a, b, x, status, fixed, shift, field3)
        ! x = A^-1 b on the rows not fixed (fixed rows: x = 0), with an
        ! optional relative shift of the diagonal; field3 replaces the rows of
        ! the third field by unit rows (projection without an f term).
        type(csr_matrix_t), intent(in) :: a
        real(dp), intent(in) :: b(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: fixed(:)
        real(dp), intent(in), optional :: shift
        logical, intent(in), optional :: field3
        type(sparse_direct_factor_t) :: factor
        logical, allocatable :: drop(:)
        integer, allocatable :: ptr(:), idx(:)
        real(dp), allocatable :: val(:), rhs(:)
        real(dp) :: dmax
        integer :: n, i, p, cnt
        n = size(b)
        allocate (drop(n))
        drop = .false.
        if (present(fixed)) drop = fixed
        if (present(field3)) then
            if (field3) drop(3::NF) = .true.
        end if
        dmax = 0
        do i = 1, n
            do p = a%rowptr(i), a%rowptr(i + 1) - 1
                if (a%colidx(p) == i) dmax = max(dmax, abs(a%val(p)))
            end do
        end do
        allocate (ptr(n + 1), idx(size(a%val) + n), val(size(a%val) + n))
        cnt = 0; ptr(1) = 1
        do i = 1, n
            if (drop(i)) then
                cnt = cnt + 1; idx(cnt) = i; val(cnt) = 1
            else
                do p = a%rowptr(i), a%rowptr(i + 1) - 1
                    if (drop(a%colidx(p))) cycle
                    cnt = cnt + 1; idx(cnt) = a%colidx(p); val(cnt) = a%val(p)
                    if (present(shift) .and. a%colidx(p) == i) &
                        val(cnt) = val(cnt) + shift*dmax
                end do
            end if
            ptr(i + 1) = cnt + 1
        end do
        rhs = b
        where (drop) rhs = 0
        ! symmetric pattern: the CSR rows are the CSC columns
        call sparse_direct_factor_csc(factor, n, ptr, idx(:cnt), val(:cnt), status)
        if (status /= 0) return
        allocate (x(n))
        call sparse_direct_solve_factored(factor, rhs, x, status)
        call sparse_direct_free(factor)
    end subroutine

end module milca_frozen_map_spline
