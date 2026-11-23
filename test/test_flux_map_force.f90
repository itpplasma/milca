program test_flux_map_force
    !! Force-balance measure of the flux map (milca_flux_map_force) against
    !! oracles that do not share its derivation:
    !!
    !! 1. the generated leaf on analytic fields: the vacuum field
    !!    B = e_phi/R (u = (ln R, Z), tau = -1, psi' = 0) has J = 0, and with
    !!    psi' = g the current is J_phi = g Delta* s/R, Delta* s = s_RR -
    !!    s_R/R + s_ZZ (curl(grad s X grad phi) = -Delta* s grad phi);
    !! 2. the 2-jet change of variables on an analytic function under a
    !!    curved parametrisation;
    !! 3. J and F of a discrete 3D state (Landreman sections, 13 planes, the
    !!    polar start: any state serves) against the curl of B_h of
    !!    flux_map3d_point by central differences in the block parameters and
    !!    t, at points of every block, next to the wall included. The same
    !!    state checks the assembled energy derivatives by differences
    !!    in the interior coefficients and wall angles;
    !! 4. the same measure through the VMEC route (vmec_state_jet,
    !!    fm3_vmec_force) on an exact state: the vacuum field B0 R0/R e_phi
    !!    in a torus of circular section, whose flux coordinates are closed
    !!    forms (concentric circles of radius r(s) with 2 pi B0 R0 (R0 -
    !!    sqrt(R0^2 - r^2)) = s Phi_edge, and theta* = 2 atan(sqrt((R0 - r)/
    !!    (R0 + r)) tan(theta/2)), i.e. lambda = 2 sum_m (-q)^m sin(m theta)/m
    !!    with q = (R0 - sqrt(R0^2 - r^2))/r). J = 0, and the residual must
    !!    fall with the radial grid at the order of the cubic splines.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_point_matrix
    use fortsparse_csr, only: csr_t
    use milca_nested_oracles, only: NESTED_SHEARED
    use milca_eulerian_clebsch, only: clebsch_nfp
    use milca_flux_map_3d, only: flux_map3d_point, flux_map3d_system
    use milca_flux_map_leaves, only: flux_map_force_3d
    use milca_flux_map_force, only: flux_map_jet_to_cylindrical, flux_map3d_force_point
    use flux_map_sections
    use flux_map_vmec, only: fm3_vmec_force
    use milca_vmec_state, only: vmec_state_t, vmec_state_set
    use milca_flux_map_force, only: flux_map_force_t, flux_map_force_norm, FORCE_J, FORCE_B, &
        FORCE_F, FORCE_GRAD_B2
    use landreman_sections, only: section_oracle, section_axis, section_boundary_radius, &
        section_fit_profiles, section_profiles
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: q0(3) = [0.3_dp, -0.2_dp, 0.5_dp], slope = 1.3_dp, hq = 1e-4_dp
    type(fm3_case_t) :: c
    integer :: failures = 0

    call leaf_on_analytic_fields()
    call change_of_variables()
    call discrete_state()
    call vmec_route()
    if (failures > 0) error stop 'flux map force balance test failed'
    print '(a)', 'PASS flux map force balance'

contains

    subroutine check(ok, what, value)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        real(dp), intent(in) :: value

        print '(a,1x,a,es10.2)', merge('PASS', 'FAIL', ok), what, value
        if (.not. ok) failures = failures + 1
    end subroutine check

    subroutine leaf_on_analytic_fields()
        real(dp), parameter :: r = 1.7_dp, z = 0.4_dp, g = 0.6_dp
        real(dp) :: b(3), j(3), f(3), m(3), jac, x, dstar, s_r

        x = log(r)
        ! vacuum: B = e_phi/R
        call flux_map_force_3d(r, x, z, 1/r, 0.0_dp, 0.0_dp, -1/r**2, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, b(1), b(2), b(3), j(1), j(2), j(3), &
            f(1), f(2), f(3), m(1), m(2), m(3), jac)
        call check(norm2(b - [0.0_dp, -1/r, 0.0_dp]) < 1e-14_dp, 'vacuum field B = -e_phi/R', &
                   norm2(b - [0.0_dp, -1/r, 0.0_dp]))
        call check(norm2(j) < 1e-14_dp, 'vacuum field: J = 0', norm2(j))
        call check(norm2(m - [-1/r**3, 0.0_dp, 0.0_dp]) < 1e-14_dp, 'vacuum field: grad B^2/2', &
                   norm2(m - [-1/r**3, 0.0_dp, 0.0_dp]))
        call check(abs(jac + 1/r**2) < 1e-14_dp, 'grad phi . (grad x X grad y)', jac + 1/r**2)
        ! poloidal field of s = ln^2 R + Z^2 with psi' = g: B_pol = g grad s X
        ! grad phi, J_phi = -g Delta* s/R
        s_r = 2*x/r
        dstar = (2 - 2*x)/r**2 - s_r/r + 2
        call flux_map_force_3d(r, x, z, 1/r, 0.0_dp, 0.0_dp, -1/r**2, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, -1.0_dp, 0.0_dp, g, 0.0_dp, 0.0_dp, b(1), b(2), b(3), j(1), j(2), j(3), &
            f(1), f(2), f(3), m(1), m(2), m(3), jac)
        call check(norm2(b - [-g*2*z/r, -1/r, g*s_r/r]) < 1e-14_dp, 'poloidal field', &
                   norm2(b - [-g*2*z/r, -1/r, g*s_r/r]))
        call check(norm2(j - [0.0_dp, -g*dstar/r, 0.0_dp]) < 1e-13_dp, &
                   'toroidal current -psi'' Delta* s/R', norm2(j - [0.0_dp, -g*dstar/r, 0.0_dp]))
    end subroutine leaf_on_analytic_fields

    subroutine change_of_variables()
        real(dp) :: rj(10), zj(10), fq(10), fc(10), exact(10), r, p, z

        call numerical_jet(1, rj)
        call numerical_jet(2, zj)
        call numerical_jet(3, fq)
        call flux_map_jet_to_cylindrical(rj, zj, slope, fq, fc)
        r = rj(1)
        z = zj(1)
        p = slope*q0(3)
        ! f = R^2 Z + R sin phi + Z^3
        exact = [r**2*z + r*sin(p) + z**3, 2*r*z + sin(p), r*cos(p), r**2 + 3*z**2, 2*z, cos(p), &
                 2*r, -r*sin(p), 0.0_dp, 6*z]
        call check(maxval(abs(fc - exact)) < 1e-6_dp, '2-jet in (R, phi, Z) from the parameters', &
                   maxval(abs(fc - exact)))
    end subroutine change_of_variables

    real(dp) function fun(which, q) result(v)
        integer, intent(in) :: which
        real(dp), intent(in) :: q(3)
        real(dp) :: rr, zz

        rr = 1.5_dp + q(1) + 0.3_dp*q(2)**2 + 0.1_dp*q(3) + 0.2_dp*q(1)*q(3)**2
        zz = q(2) + 0.2_dp*q(1)*q(3) + 0.1_dp*q(1)**2
        select case (which)
        case (1)
            v = rr
        case (2)
            v = zz
        case default
            v = rr**2*zz + rr*sin(slope*q(3)) + zz**3
        end select
    end function fun

    subroutine numerical_jet(which, jet)
        !! central differences of the fourth order
        integer, intent(in) :: which
        real(dp), intent(out) :: jet(10)
        integer, parameter :: pos(3, 3) = reshape([5, 6, 7, 6, 8, 9, 7, 9, 10], [3, 3])
        real(dp) :: e(3, 3)
        integer :: a, b

        e = 0
        do a = 1, 3
            e(a, a) = hq
        end do
        jet(1) = fun(which, q0)
        do a = 1, 3
            jet(1 + a) = (8*(fun(which, q0 + e(:, a)) - fun(which, q0 - e(:, a))) &
                          - (fun(which, q0 + 2*e(:, a)) - fun(which, q0 - 2*e(:, a))))/(12*hq)
            jet(pos(a, a)) = (-fun(which, q0 + 2*e(:, a)) + 16*fun(which, q0 + e(:, a)) &
                - 30*jet(1) + 16*fun(which, q0 - e(:, a)) - fun(which, q0 - 2*e(:, a)))/(12*hq**2)
            do b = a + 1, 3
                jet(pos(a, b)) = (fun(which, q0 + e(:, a) + e(:, b)) &
                    - fun(which, q0 + e(:, a) - e(:, b)) - fun(which, q0 - e(:, a) + e(:, b)) &
                    + fun(which, q0 - e(:, a) - e(:, b)))/(4*hq**2)
            end do
        end do
    end subroutine numerical_jet

    subroutine vacuum_profiles(s, tau, gp, p)
        !! Phi_edge = 1, no transform, no pressure
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)

        tau = [-1/pi, 0.0_dp, 0.0_dp]
        gp = 0
        p = 0
        associate (unused => s)
        end associate
    end subroutine vacuum_profiles

    subroutine vmec_route()
        integer, parameter :: mmax = 14, grids(2) = [26, 51]
        real(dp), parameter :: r0 = 3.0_dp, a = 0.9_dp
        type(vmec_state_t) :: st
        type(flux_map_force_t) :: acc(2)
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), lmns(:, :)
        real(dp) :: xm(0:mmax), xn(0:mmax), b0, s, r, q, res(2), best
        integer :: level, ns, j, m, o

        ! Phi_edge = 1 = 2 pi B0 R0 (R0 - sqrt(R0^2 - a^2))
        b0 = 1/(2*pi*r0*(r0 - sqrt(r0**2 - a**2)))
        xm = [(real(m, dp), m=0, mmax)]
        xn = 0
        do level = 1, 2
            ns = grids(level)
            if (allocated(rmnc)) deallocate (rmnc, zmns, lmns)
            allocate (rmnc(0:mmax, ns), zmns(0:mmax, ns), lmns(0:mmax, ns), source=0.0_dp)
            do j = 1, ns
                s = (j - 1.0_dp)/(ns - 1)
                r = sqrt(r0**2 - (r0 - s/(2*pi*b0*r0))**2)
                rmnc(0, j) = r0
                rmnc(1, j) = r
                zmns(1, j) = r
                if (j == 1) cycle
                ! lambda on the half grid
                s = (j - 1.5_dp)/(ns - 1)
                r = sqrt(r0**2 - (r0 - s/(2*pi*b0*r0))**2)
                q = (r0 - sqrt(r0**2 - r**2))/r
                do m = 1, mmax
                    lmns(m, j) = 2*(-q)**m/m
                end do
            end do
            call vmec_state_set(st, 1, xm, xn, rmnc, zmns, lmns, 1.0_dp)
            best = huge(1.0_dp)
            do o = 1, 2
                call fm3_vmec_force(st, vacuum_profiles, merge(1.0_dp, -1.0_dp, o == 1), 10, 60, &
                                    acc(o))
                best = min(best, flux_map_force_norm(acc(o), FORCE_F, 0))
            end do
            o = merge(1, 2, flux_map_force_norm(acc(1), FORCE_F, 0) <= best)
            ! outside s = 0.1, where the radial splines do not feel the axis
            res(level) = sqrt(sum(acc(o)%sums(FORCE_F, 2:))/sum(acc(o)%sums(FORCE_GRAD_B2, 2:)))
            print '(a,i3,a,es10.3,a,es10.3,a,es10.3)', '  exact vacuum state through the VMEC route, ns', &
                ns, ': |F|/|grad B^2/2| outside s = 0.1', res(level), '  whole volume', &
                flux_map_force_norm(acc(o), FORCE_F, 0)/flux_map_force_norm(acc(o), FORCE_GRAD_B2, 0), &
                '  |B| error', abs(flux_map_force_norm(acc(o), FORCE_B, 0) &
                                   - b0*r0*sqrt(sum_inverse_r2(a, r0)))
        end do
        call check(res(2) < 1e-4_dp, 'VMEC route: J X B = 0 of the exact vacuum state', res(2))
        call check(log(res(1)/res(2))/log(2.0_dp) > 3.0_dp, &
                   'VMEC route: order above 3 in the radial grid', log(res(1)/res(2))/log(2.0_dp))
    end subroutine vmec_route

    real(dp) function sum_inverse_r2(a, r0) result(v)
        !! Volume mean of 1/R^2 over the torus of circular section:
        !! int dA/R over int R dA = 2 pi (R0 - sqrt(R0^2 - a^2))/(pi a^2 R0).
        real(dp), intent(in) :: a, r0

        v = 2*(r0 - sqrt(r0**2 - a**2))/(a**2*r0)
    end function sum_inverse_r2

    subroutine profiles(s, tau, gp, p)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: tau(3), gp(3), p(3)
        real(dp) :: f(3)

        call section_profiles(s, f, gp, p)
        tau = -2*f
    end subroutine profiles

    subroutine discrete_state()
        integer, parameter :: nsamples = 128, nt = 13
        real(dp), parameter :: h = 1e-5_dp, t0 = 0.283_dp
        real(dp), parameter :: points(2, 3) = reshape([0.37_dp, 0.41_dp, 0.37_dp, 0.93_dp, &
                                                       0.93_dp, 0.37_dp], [2, 3])
        real(dp) :: samples(2, nsamples, nt), a(2), tk, phi, rz(2), b(3), j(3), f(3), &
                    bp(3), bm(3), rp(2), rm(2), db(3, 3), amat(3, 3), ainv(3, 3), gb(3, 3), &
                    jfd(3), ffd(3), u(8), gs(3), tau(3), gp(3), pr(3), q(3), worst_j, worst_f, &
                    scale_j
        integer :: st, k, i, block, ip, d, wall_points

        clebsch_nfp = 1
        section_oracle%family = NESTED_SHEARED
        section_oracle%eps = 0.6_dp
        section_oracle%saxis = 2.4_dp
        section_oracle%delta = 0.2_dp
        call section_fit_profiles()
        do k = 1, nt
            phi = 2*pi*(k - 1)/nt
            a = section_axis(phi)
            do i = 1, nsamples
                tk = 2*pi*(i - 1)/nsamples
                samples(:, i, k) = a + section_boundary_radius(phi, a, tk)*[cos(tk), sin(tk)]
            end do
        end do
        call fm3_build(c, samples, 2, 3, 4, st)
        if (st /= 0) error stop 'meshes'
        call fm3_start_polar(c)
        worst_j = 0
        worst_f = 0
        scale_j = 0
        wall_points = 0
        do block = 1, c%space%layout_block_count
            do ip = 1, 3
                call flux_map3d_force_point(c%fm, c%derham, profiles, c%base, c%fm%z, block, &
                                            points(:, ip), t0, rz, b, j, f, st)
                if (st /= 0) error stop 'force point'
                q = [points(:, ip), t0]
                do d = 1, 3
                    q(d) = q(d) + h
                    call field(block, q, rp, bp, u)
                    q(d) = q(d) - 2*h
                    call field(block, q, rm, bm, u)
                    q(d) = q(d) + h
                    db(:, d) = (bp - bm)/(2*h)
                    amat(1, d) = (rp(1) - rm(1))/(2*h)
                    amat(3, d) = (rp(2) - rm(2))/(2*h)
                end do
                amat(2, :) = [0.0_dp, 0.0_dp, 2*pi/clebsch_nfp]
                call invert(amat, ainv)
                ! gb(c, m) = d B_c/d(R, phi, Z)_m
                gb = matmul(db, ainv)
                call field(block, q, rp, bp, u)
                if (any(c%fm%is_wall(wall_functions(block, points(:, ip))))) &
                    wall_points = wall_points + 1
                jfd = [gb(3, 2)/rz(1) - gb(2, 3), gb(1, 3) - gb(3, 1), &
                       bp(2)/rz(1) + gb(2, 1) - gb(1, 2)/rz(1)]
                call profiles(u(1)**2 + u(2)**2, tau, gp, pr)
                gs = 2*(u(1)*u(3:5) + u(2)*u(6:8))
                ffd = [jfd(2)*bp(3) - jfd(3)*bp(2), jfd(3)*bp(1) - jfd(1)*bp(3), &
                       jfd(1)*bp(2) - jfd(2)*bp(1)] - pr(2)*gs
                call check(norm2(b - bp) < 1e-11_dp*norm2(bp), &
                           'B of the 2-jet path equals flux_map3d_point', norm2(b - bp)/norm2(bp))
                worst_j = max(worst_j, norm2(j - jfd))
                worst_f = max(worst_f, norm2(f - ffd))
                scale_j = max(scale_j, norm2(jfd))
            end do
        end do
        call check(wall_points >= 4, 'points next to the wall are included', real(wall_points, dp))
        call check(worst_j < 1e-6_dp*scale_j, 'J = curl B_h against central differences of B_h', &
                   worst_j/scale_j)
        call check(worst_f < 1e-6_dp*scale_j, 'F = J X B - grad p against the same', &
                   worst_f/scale_j)
        call discrete_energy_derivatives()
    end subroutine discrete_state

    subroutine discrete_energy_derivatives()
        real(dp), parameter :: h = 1e-5_dp
        type(csr_t) :: hessian
        real(dp), allocatable :: gradient(:), plus(:), minus(:), direction(:), product(:)
        real(dp) :: energy, energy_plus, energy_minus, derivative, relative_error
        integer :: i, j, entry, n

        n = size(c%fm%z)
        direction = [(sin(real(i, dp)), i=1, n)]
        direction(c%fm%fixed) = 0
        direction = direction/norm2(direction)
        call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z, energy, gradient, hessian)
        call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z + h*direction, energy_plus, plus)
        call flux_map3d_system(c%fm, c%derham, c%quad, profiles, c%base, c%fm%z - h*direction, energy_minus, minus)
        derivative = (energy_plus - energy_minus)/(2*h)
        relative_error = abs(dot_product(gradient, direction) - derivative)/max(abs(derivative), 1e-5_dp)
        call check(relative_error < 1e-6_dp, 'energy gradient against central differences', relative_error)

        allocate (product(n), source=0.0_dp)
        do i = 1, n
            do entry = hessian%row_ptr(i), hessian%row_ptr(i + 1) - 1
                j = hessian%col(entry)
                product(i) = product(i) + hessian%val(entry)*direction(j)
                if (i /= j) product(j) = product(j) + hessian%val(entry)*direction(i)
            end do
        end do
        plus = (plus - minus)/(2*h)
        relative_error = norm2(product - plus)/norm2(plus)
        call check(relative_error < 1e-7_dp, 'assembled Hessian action against central differences', relative_error)
    end subroutine discrete_energy_derivatives

    function wall_functions(blk, xi) result(d0)
        integer, intent(in) :: blk
        real(dp), intent(in) :: xi(2)
        integer, allocatable :: d0(:), rowptr(:)
        real(dp), allocatable :: val(:)
        integer :: s

        call spline_point_matrix(c%space, [blk], reshape(xi, [2, 1]), 0, rowptr, d0, val, s)
        if (s /= 0) error stop 'point'
    end function wall_functions

    subroutine field(blk, qq, pos, bb, uu)
        integer, intent(in) :: blk
        real(dp), intent(in) :: qq(3)
        real(dp), intent(out) :: pos(2), bb(3), uu(8)
        integer, allocatable :: rowptr(:), colidx(:), d0(:)
        real(dp), allocatable :: val(:), nv(:), gn(:, :)
        integer :: comp, s

        do comp = 0, 2
            call spline_point_matrix(c%space, [blk], reshape(qq(1:2), [2, 1]), comp, rowptr, &
                                     colidx, val, s)
            if (s /= 0) error stop 'point'
            if (comp == 0) then
                d0 = colidx
                nv = val
                allocate (gn(2, size(d0)))
            else
                gn(comp, :) = val
            end if
        end do
        call flux_map3d_point(c%fm, c%base, d0, nv, gn, qq(3), profiles, pos, uu, bb)
    end subroutine field

    subroutine invert(m, mi)
        real(dp), intent(in) :: m(3, 3)
        real(dp), intent(out) :: mi(3, 3)
        real(dp) :: det

        det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - &
              m(2, 3)*m(3, 1)) + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))
        mi(1, :) = [m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2), m(1, 3)*m(3, 2) - m(1, 2)*m(3, 3), &
                    m(1, 2)*m(2, 3) - m(1, 3)*m(2, 2)]
        mi(2, :) = [m(2, 3)*m(3, 1) - m(2, 1)*m(3, 3), m(1, 1)*m(3, 3) - m(1, 3)*m(3, 1), &
                    m(1, 3)*m(2, 1) - m(1, 1)*m(2, 3)]
        mi(3, :) = [m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1), m(1, 2)*m(3, 1) - m(1, 1)*m(3, 2), &
                    m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)]
        mi = mi/det
    end subroutine invert

end program test_flux_map_force
