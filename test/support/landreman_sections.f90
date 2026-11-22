module landreman_sections
    !! Toroidal sections of the exact nested 3D equilibria of Landreman
    !! (docs/nested-oracles.md) for the label-lane gates E4a and E4: the
    !! magnetic axis in a plane phi, the label chi = psi/delta, the boundary
    !! psi = delta along rays from the axis, H^(1/2) of chi at a point and
    !! the O-grid of a plane, centred on its axis and fitted to its boundary.
    !! Both families (integer and sheared) through nested_oracle_t; the
    !! oracle is module state because the exact-field callback of
    !! milca_eulerian_clebsch takes no context.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    use milca_nested_oracles, only: nested_oracle_t
    use milca_nested_oracle_leaves, only: nested_ray
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: samples_n = 128
    type(nested_oracle_t), public :: section_oracle

    public :: section_axis, section_label, section_exact_b, section_boundary_radius, &
              section_sqrt_hessian, section_layout

contains

    function section_axis(phi) result(x)
        !! (R, Z) of the magnetic axis in the plane phi: the surface map at
        !! rho = 0, its parameter t found where the ray angle matches phi.
        real(dp), intent(in) :: phi
        real(dp) :: x(2), lo, hi, t, c(3), jac, ray, radius
        integer :: k

        lo = phi - pi/2
        hi = phi + pi/2
        do k = 1, 60
            t = (lo + hi)/2
            call section_oracle%surface(0.0_dp, 0.0_dp, t, c, jac)
            call nested_ray(c(1), c(2), phi, ray, radius)
            if (ray < 0) then
                lo = t
            else
                hi = t
            end if
        end do
        x = [sqrt(c(1)**2 + c(2)**2), c(3)]
    end function section_axis

    real(dp) function section_label(r, phi, z) result(chi)
        !! chi = psi/delta at the cylindrical point (R, phi, Z).
        real(dp), intent(in) :: r, phi, z
        real(dp) :: v(20)
        integer :: st

        call section_oracle%point([r*cos(phi), r*sin(phi), z], v, st)
        if (st /= 0) error stop 'landreman_sections: point outside the domain'
        chi = v(1)/section_oracle%delta
    end function section_label

    subroutine section_exact_b(x, b)
        !! Cylindrical (R, phi, Z) components of the exact field at x = (R, phi, Z).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: v(20)
        integer :: st

        call section_oracle%point([x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)], v, st)
        if (st /= 0) error stop 'landreman_sections: exact field outside the domain'
        b = [v(3)*cos(x(2)) + v(4)*sin(x(2)), -v(3)*sin(x(2)) + v(4)*cos(x(2)), v(5)]
    end subroutine section_exact_b

    real(dp) function section_boundary_radius(phi, centre, angle) result(r)
        !! Distance from centre to chi = 1 along a ray in the plane phi.
        real(dp), intent(in) :: phi, centre(2), angle
        real(dp) :: lo, hi, e(2)
        integer :: k

        e = [cos(angle), sin(angle)]
        lo = 0
        hi = 0.05_dp
        do while (section_label(centre(1) + hi*e(1), phi, centre(2) + hi*e(2)) < 1)
            lo = hi
            hi = 1.5_dp*hi
        end do
        do k = 1, 60
            r = (lo + hi)/2
            if (section_label(centre(1) + r*e(1), phi, centre(2) + r*e(2)) < 1) then
                lo = r
            else
                hi = r
            end if
        end do
    end function section_boundary_radius

    function section_sqrt_hessian(phi, x) result(m)
        !! H^(1/2) of the exact chi in the plane phi at x:
        !! (H + sqrt(det H) I)/sqrt(tr H + 2 sqrt(det H)).
        real(dp), intent(in) :: phi, x(2)
        real(dp) :: m(2, 2), hm(2, 2), sd
        real(dp), parameter :: h = 1e-3_dp

        hm(1, 1) = (section_label(x(1) + h, phi, x(2)) - 2*section_label(x(1), phi, x(2)) &
                    + section_label(x(1) - h, phi, x(2)))/h**2
        hm(2, 2) = (section_label(x(1), phi, x(2) + h) - 2*section_label(x(1), phi, x(2)) &
                    + section_label(x(1), phi, x(2) - h))/h**2
        hm(1, 2) = (section_label(x(1) + h, phi, x(2) + h) - section_label(x(1) + h, phi, x(2) - h) &
                    - section_label(x(1) - h, phi, x(2) + h) + &
                    section_label(x(1) - h, phi, x(2) - h))/(4*h**2)
        hm(2, 1) = hm(1, 2)
        sd = sqrt(hm(1, 1)*hm(2, 2) - hm(1, 2)**2)
        m = hm
        m(1, 1) = m(1, 1) + sd
        m(2, 2) = m(2, 2) + sd
        m = m/sqrt(hm(1, 1) + hm(2, 2) + 2*sd)
    end function section_sqrt_hessian

    subroutine section_layout(phi, centre, layout, core_centre)
        !! O-grid of the plane phi centred on its axis (returned in centre)
        !! and fitted to its boundary chi = 1, sampled by rays from the axis
        !! counter-clockwise from the outboard midplane. With core_centre the
        !! core is centred there instead (the same boundary): the mesh of an
        !! axis guess, whose core follows the axis (ALE).
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: centre(2)
        type(block_layout_t), intent(out) :: layout
        real(dp), intent(in), optional :: core_centre(2)
        real(dp) :: samples(2, samples_n), angle, c(2)
        integer :: i, st

        centre = section_axis(phi)
        do i = 1, samples_n
            angle = 2*pi*(i - 1)/samples_n
            samples(:, i) = centre + section_boundary_radius(phi, centre, angle)* &
                            [cos(angle), sin(angle)]
        end do
        c = centre
        if (present(core_centre)) c = core_centre
        call layout_ogrid(c, 0.5_dp*(maxval(samples(1, :)) - minval(samples(1, :))), &
                          0.5_dp, layout, st)
        if (st /= 0) error stop 'landreman_sections: layout'
        ! the core centred on c, so the axis is the central vertex of the
        ! core block on every plane
        call layout_set_curve(layout, LAYOUT_WALL, samples, st, centre_core=.true.)
        if (st /= 0) error stop 'landreman_sections: boundary curve'
    end subroutine section_layout

end module landreman_sections
