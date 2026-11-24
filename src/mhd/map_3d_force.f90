module milca_map_3d_force
    !! Force balance J X B - grad p of a Lagrangian-map state (ideal_map_t with
    !! the toroidal seed), on the measure of the flux map
    !! (milca_flux_map_force, docs/flux-map.md "Force-balance residual"), so
    !! that both directions of the ideal map are judged by one kernel.
    !!
    !! The map X: reference disc -> space carries the labels with it: with the
    !! reference coordinates (u, v) of the seed's disc, s = u^2 + v^2 and the
    !! field-line angle theta* = atan2(v, u) + eta (eta: the boundary
    !! field-line angle of the seed, zero without it), the inverse flux map
    !! is x + i y = (u + i v) exp(i eta) as a function of position. Its 2-jet
    !! in (R, phi, Z) follows inside every cell from the 2-jets of R, Z (the
    !! element functions times the toroidal Fourier modes) by the chain rule
    !! of flux_map_jet_to_cylindrical; J = curl B comes from these second
    !! derivatives, not from the gradient of the energy. Nothing is fitted
    !! or smoothed: the derivatives are those of the discrete map inside its
    !! cells (jumps of the tangential field across cell edges are not seen by
    !! this volume measure, as for the flux map).
    !!
    !! Cylindrical chart with the toroidal label physical (no phi
    !! displacement) only.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_ideal_map, only: ideal_map_t
    use milca_ideal_seed_3d, only: toroidal_seed_t
    use milca_triangle_basis, only: triangle_basis_jet, TRIANGLE_MAX_NODES
    use milca_curve_chart_leaves, only: curve_fourier_mode
    use milca_flux_map, only: flux_map_profiles
    use milca_flux_map_leaves, only: flux_map_force_3d
    use milca_flux_map_force, only: flux_map_force_t, flux_map_force_reset, flux_map_force_add, &
        flux_map_force_merge, flux_map_jet_to_cylindrical, flux_map_jet_product, &
        flux_map_jet_cos_sin
    implicit none
    private
    public :: map_3d_force, map_3d_force_jets, map_3d_force_field

    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)

contains

    subroutine map_3d_force_jets(map, seed, cell, xi, phi, orientation, rj, zj, xq, yq, status)
        !! At the point xi of a cell and the toroidal angle phi: the 2-jets
        !! (f, f_1, f_2, f_3, f_11, f_12, f_13, f_22, f_23, f_33) in the
        !! parameters (xi_1, xi_2, phi) of R, Z and of the flux-map
        !! components x, y = orientation * Im. Status 1: not a cylindrical
        !! chart or a map that displaces phi.
        type(ideal_map_t), intent(in) :: map
        type(toroidal_seed_t), intent(in) :: seed
        integer, intent(in) :: cell
        real(dp), intent(in) :: xi(2), phi, orientation
        real(dp), intent(out) :: rj(10), zj(10), xq(10), yq(10)
        integer, intent(out) :: status
        real(dp) :: shape(TRIANGLE_MAX_NODES), gradient(2, TRIANGLE_MAX_NODES), &
                    hessian(2, 2, TRIANGLE_MAX_NODES), section(6, 3), f(5), c, s, phase, &
                    uj(10), vj(10), ej(10), re(10), im(10), next(10), cj(10), sj(10), &
                    cosj(10), sinj(10), unit(10)
        real(dp), allocatable :: coefficient(:, :)
        integer :: k, n, a, i, m, level

        status = 1
        select type (chart => map%chart)
        type is (cylindrical_chart_t)
        class default
            return
        end select
        if (any(map%x(3, :, :) /= 0)) return
        n = map%nloc
        call triangle_basis_jet(xi, shape(:n), gradient(:, :n), hessian(:, :, :n))
        allocate (coefficient(3, n))
        rj = 0
        zj = 0
        do k = 1, map%nmode
            f = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
            if (k > 1) then
                phase = real((k/2)*map%nfp, dp)
                c = merge(1.0_dp, 0.0_dp, mod(k, 2) == 0)
                s = 1 - c
                call curve_fourier_mode(c, s, phase, phi, f(1), f(2), f(3), f(4), f(5))
            end if
            coefficient = map%x(:, k, map%conn(:, cell))
            do a = 1, 2
                call section_jet(coefficient(a, :), section(:, a))
            end do
            rj = rj + toroidal_jet(section(:, 1), f)
            zj = zj + toroidal_jet(section(:, 2), f)
        end do
        ! reference coordinates of the seed's disc: functions of xi only
        call section_jet((map%reference(1, map%conn(:, cell)) - seed%centre(1))/seed%axes(1), &
                         section(:, 1))
        call section_jet((map%reference(2, map%conn(:, cell)) - seed%centre(2))/seed%axes(2), &
                         section(:, 2))
        f = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        uj = toroidal_jet(section(:, 1), f)
        vj = toroidal_jet(section(:, 2), f)
        ! eta = sum_i ca_i (c Re + s Im) + cb_i (c Im - s Re), Re + i Im =
        ! (u + i v)^m, c + i s = exp(i n nfp phi) (leaf field_angle_mode)
        ej = 0
        unit = 0
        unit(1) = 1
        if (map%nangle > 0) then
            level = 0
            re = unit
            im = 0
            do m = 0, map%angle_mpol
                do while (level < m)
                    next = flux_map_jet_product(re, uj) - flux_map_jet_product(im, vj)
                    im = flux_map_jet_product(re, vj) + flux_map_jet_product(im, uj)
                    re = next
                    level = level + 1
                end do
                do i = 1, size(map%angle_m)
                    if (map%angle_m(i) /= m) cycle
                    phase = real(map%angle_n(i)*map%nfp, dp)
                    call curve_fourier_mode(1.0_dp, 0.0_dp, phase, phi, f(1), f(2), f(3), f(4), &
                                            f(5))
                    cosj = 0
                    cosj([1, 4, 10]) = f(1:3)
                    call curve_fourier_mode(0.0_dp, 1.0_dp, phase, phi, f(1), f(2), f(3), f(4), &
                                            f(5))
                    sinj = 0
                    sinj([1, 4, 10]) = f(1:3)
                    ej = ej + map%angle(2*i - 1)*(flux_map_jet_product(cosj, re) &
                                                  + flux_map_jet_product(sinj, im)) &
                         + map%angle(2*i)*(flux_map_jet_product(cosj, im) &
                                           - flux_map_jet_product(sinj, re))
                end do
            end do
        end if
        call flux_map_jet_cos_sin(ej, cj, sj)
        xq = flux_map_jet_product(uj, cj) - flux_map_jet_product(vj, sj)
        yq = orientation*(flux_map_jet_product(uj, sj) + flux_map_jet_product(vj, cj))
        status = 0

    contains

        subroutine section_jet(nodal, jet)
            !! (f, f_1, f_2, f_11, f_12, f_22) in xi of an element function.
            real(dp), intent(in) :: nodal(:)
            real(dp), intent(out) :: jet(6)

            jet(1) = dot_product(nodal, shape(:n))
            jet(2) = dot_product(nodal, gradient(1, :n))
            jet(3) = dot_product(nodal, gradient(2, :n))
            jet(4) = dot_product(nodal, hessian(1, 1, :n))
            jet(5) = dot_product(nodal, hessian(1, 2, :n))
            jet(6) = dot_product(nodal, hessian(2, 2, :n))
        end subroutine section_jet

        pure function toroidal_jet(sec, mode) result(jet)
            !! 2-jet in (xi_1, xi_2, phi) of sec(xi) mode(phi), mode = (f, f',
            !! f'', ...).
            real(dp), intent(in) :: sec(6), mode(5)
            real(dp) :: jet(10)

            jet = [sec(1)*mode(1), sec(2)*mode(1), sec(3)*mode(1), sec(1)*mode(2), &
                   sec(4)*mode(1), sec(5)*mode(1), sec(2)*mode(2), sec(6)*mode(1), &
                   sec(3)*mode(2), sec(1)*mode(3)]
        end function toroidal_jet
    end subroutine map_3d_force_jets

    subroutine map_3d_force_field(map, seed, profiles, orientation, cell, xi, phi, b, j, f, status)
        !! B, J = curl B and F = J X B - grad p, components (R, phi, Z), of
        !! the flux-map form at one point: B is to agree with the frozen-in
        !! field of the map (map_3d_sample), which is computed another way.
        type(ideal_map_t), intent(in) :: map
        type(toroidal_seed_t), intent(in) :: seed
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: orientation, xi(2), phi
        integer, intent(in) :: cell
        real(dp), intent(out) :: b(3), j(3), f(3)
        integer, intent(out) :: status
        real(dp) :: rj(10), zj(10), xq(10), yq(10), xc(10), yc(10), tau(3), gp(3), pr(3), m(3), &
                    jac

        call map_3d_force_jets(map, seed, cell, xi, phi, orientation, rj, zj, xq, yq, status)
        if (status /= 0) return
        call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, xq, xc)
        call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, yq, yc)
        call profiles(xc(1)**2 + yc(1)**2, tau, gp, pr)
        call flux_map_force_3d(rj(1), xc(1), yc(1), xc(2), xc(3), xc(4), xc(5), xc(6), xc(7), &
            xc(8), xc(9), xc(10), yc(2), yc(3), yc(4), yc(5), yc(6), yc(7), yc(8), yc(9), yc(10), &
            tau(1), tau(2), gp(1), gp(2), pr(2), b(1), b(2), b(3), j(1), j(2), j(3), f(1), f(2), &
            f(3), m(1), m(2), m(3), jac)
    end subroutine map_3d_force_field

    subroutine map_3d_force(map, seed, profiles, orientation, bands, acc, status)
        !! The force balance of the map over Gauss points of every cell
        !! (collapsed square, 8 per direction) and max(16, 4 ntor + 4) planes
        !! of one field period. orientation +1 or -1: the sense of theta*
        !! relative to the reference angle that matches the sign convention
        !! of profiles (checked against the frozen-in field by the caller). profiles:
        !! those of the flux map, tau = Phi'/pi, psi' and p of the label s.
        !! Status 1 includes a failed geometry sample or a singular or
        !! inconsistently oriented flux label at the quadrature points.
        type(ideal_map_t), intent(in) :: map
        type(toroidal_seed_t), intent(in) :: seed
        procedure(flux_map_profiles) :: profiles
        real(dp), intent(in) :: orientation
        integer, intent(in) :: bands
        type(flux_map_force_t), intent(out) :: acc
        integer, intent(out) :: status
        integer, parameter :: nq = 8
        type(flux_map_force_t) :: part
        real(dp) :: x(nq), w(nq), xi(2), weight, phi, rj(10), zj(10), xq(10), yq(10), xc(10), &
                    yc(10), volume_jacobian
        integer :: e, i, j, k, nphi, local, failed

        call gauss_legendre(nq, x, w)
        x = (x + 1)/2
        w = w/2
        nphi = max(16, 4*map%ntor + 4)
        call flux_map_force_reset(acc, bands)
        failed = 0
        !$omp parallel default(private) shared(map, seed, acc, x, w, nphi, bands, orientation) &
        !$omp     reduction(+:failed)
        call flux_map_force_reset(part, bands)
        !$omp do schedule(dynamic, 4)
        do e = 1, map%cells
            do k = 0, nphi - 1
                phi = two_pi*k/(map%nfp*nphi)
                do j = 1, nq
                    do i = 1, nq
                        xi = [x(i)*(1 - x(j)), x(j)]
                        weight = w(i)*w(j)*(1 - x(j))
                        call map_3d_force_jets(map, seed, e, xi, phi, orientation, rj, zj, xq, &
                                               yq, local)
                        if (local /= 0) then
                            failed = failed + 1
                            cycle
                        end if
                        volume_jacobian = -map%sigma*rj(1)*(rj(2)*zj(3) - rj(3)*zj(2))
                        if (.not. (rj(1) > 0 .and. volume_jacobian > 0)) then
                            failed = failed + 1
                            cycle
                        end if
                        call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, xq, xc)
                        call flux_map_jet_to_cylindrical(rj, zj, 1.0_dp, yq, yc)
                        call flux_map_force_add(part, weight*volume_jacobian &
                                                *two_pi/nphi, rj(1), xc, yc, profiles)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp critical (map_3d_force_merge_parts)
        call flux_map_force_merge(acc, part)
        !$omp end critical (map_3d_force_merge_parts)
        !$omp end parallel
        status = 1
        if (failed /= 0 .or. acc%points == 0) return
        if (acc%jacobian_min > 0 .or. acc%jacobian_max < 0) status = 0
    end subroutine map_3d_force

end module milca_map_3d_force
