module milca_eulerian_clebsch
    !! Axisymmetric label lane (docs/eulerian-label-plan.md): energy
    !! W = 2 pi int w dR dZ of the direct Clebsch representation with the
    !! label chi and g in V0 of a spline section, given profiles F(chi),
    !! G(chi), p(chi) and the winding form about a given axis. Returns W, its
    !! gradient in (chi, g) node values and the Hessian, from the generated
    !! leaf clebsch_energy_axisym by the chain rule.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_cell_basis
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_clebsch_leaves, only: clebsch_energy_axisym
    implicit none
    private

    abstract interface
        subroutine clebsch_profiles(chi, f, g, p)
            !! f = (F', F'', F'''), g = (G', G'', G'''), p = (p, p', p'') at chi.
            import :: dp
            real(dp), intent(in) :: chi
            real(dp), intent(out) :: f(3), g(3), p(3)
        end subroutine clebsch_profiles
    end interface

    public :: clebsch_profiles, clebsch_axisym_system, clebsch_node_gradients

contains

    subroutine clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
        !! Values v0(1, i, q) and physical (R, Z) gradients gn(:, i, q) of the
        !! node functions of a section cell, the gradients from the incidence
        !! grad N_n = sum_e G_en E_e and the inverse transposed Jacobian.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell
        real(dp), allocatable, intent(out) :: v0(:, :, :), gn(:, :, :)
        integer, allocatable, intent(out) :: d0(:)
        real(dp), allocatable :: v1(:, :, :)
        integer, allocatable :: d1(:), s0(:), s1(:)
        real(dp) :: jet(2, 6), jinv(2, 2), det
        integer :: q, p, j, e, kk, i, st

        call spline_derham_cell_basis(derham, quad, cell, 0, v0, d0, s0)
        call spline_derham_cell_basis(derham, quad, cell, 1, v1, d1, s1)
        allocate (gn(2, size(d0), quad%points_per_cell))
        gn = 0
        do q = 1, quad%points_per_cell
            p = quad%point_offset(cell) + q - 1
            do j = 1, size(d1)
                e = d1(j)
                do kk = derham%grad%rowptr(e), derham%grad%rowptr(e + 1) - 1
                    do i = 1, size(d0)
                        if (d0(i) == derham%grad%colidx(kk)) gn(:, i, q) = gn(:, i, q) + &
                            derham%grad%val(kk)*s1(j)*v1(:, j, q)
                    end do
                end do
            end do
            call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
            det = jet(1, 2)*jet(2, 3) - jet(1, 3)*jet(2, 2)
            jinv = reshape([jet(2, 3), -jet(2, 2), -jet(1, 3), jet(1, 2)], [2, 2])/det
            gn(:, :, q) = matmul(transpose(jinv), gn(:, :, q))
        end do
    end subroutine clebsch_node_gradients

    subroutine clebsch_axisym_system(derham, quad, profiles, axis, chi, g, energy, &
                                     gradient, hessian)
        !! Unknowns ordered (chi(1:n0), g(1:n0)).
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: axis(2), chi(:), g(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: gradient(:)
        type(csr_t), intent(out) :: hessian
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), vals(:)
        integer, allocatable :: d0(:), rows(:), cols(:), idx(:)
        real(dp) :: s(2), wq, c, gc(2), gg(2), th(2), rho2, f(3), gp(3), pr(3), w, &
                    gu(5), h(15), hm(5, 5), jet(2, 6)
        real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
        integer :: n0, cell, q, p, i, j, a, b, k, nloc, nz, st
        type(csc_t) :: m
        type(fortsparse_status_t) :: info

        n0 = derham%dimension(0)
        allocate (gradient(2*n0))
        gradient = 0
        energy = 0
        nz = quad%point_count*4*(derham%degree + 1)**4
        allocate (rows(nz), cols(nz), vals(nz))
        nz = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            if (allocated(du)) deallocate (du, idx)
            allocate (du(5, 2*nloc), idx(2*nloc))
            idx = [d0, n0 + d0]
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
                s = jet(:, 1)
                wq = two_pi*quad%weight(p)
                c = dot_product(v0(1, :, q), chi(d0))
                gc = matmul(gn(:, :, q), chi(d0))
                gg = matmul(gn(:, :, q), g(d0))
                rho2 = (s(1) - axis(1))**2 + (s(2) - axis(2))**2
                th = [-(s(2) - axis(2)), s(1) - axis(1)]/rho2
                call profiles(c, f, gp, pr)
                call clebsch_energy_axisym(s(1), c, gc(1), gc(2), gg(1), gg(2), th(1), th(2), &
                                           f(1), f(2), f(3), gp(1), gp(2), gp(3), pr(1), &
                                           pr(2), pr(3), w, gu(1), gu(2), gu(3), gu(4), &
                                           gu(5), h(1), h(2), h(3), h(4), h(5), h(6), h(7), &
                                           h(8), h(9), h(10), h(11), h(12), h(13), h(14), &
                                           h(15))
                energy = energy + wq*w
                ! du(a, local): d u_a / d unknown, u = (c, c_R, c_Z, g_R, g_Z)
                du = 0
                du(1, 1:nloc) = v0(1, :, q)
                du(2, 1:nloc) = gn(1, :, q)
                du(3, 1:nloc) = gn(2, :, q)
                du(4, nloc + 1:) = gn(1, :, q)
                du(5, nloc + 1:) = gn(2, :, q)
                gradient(idx) = gradient(idx) + wq*matmul(gu, du)
                k = 0
                do a = 1, 5
                    do b = a, 5
                        k = k + 1
                        hm(a, b) = h(k)
                        hm(b, a) = h(k)
                    end do
                end do
                do j = 1, 2*nloc
                    do i = 1, 2*nloc
                        nz = nz + 1
                        rows(nz) = idx(i)
                        cols(nz) = idx(j)
                        vals(nz) = wq*dot_product(du(:, i), matmul(hm, du(:, j)))
                    end do
                end do
            end do
        end do
        call csc_from_triplet(2*n0, 2*n0, cols(:nz), rows(:nz), vals(:nz), m, info)
        if (info%code /= 0) error stop 'Clebsch Hessian triplets'
        hessian%nrow = 2*n0
        hessian%ncol = 2*n0
        hessian%row_ptr = m%col_ptr
        hessian%col = m%row_idx
        hessian%val = m%val
    end subroutine clebsch_axisym_system


end module milca_eulerian_clebsch
