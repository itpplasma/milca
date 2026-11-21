module milca_eulerian_clebsch
    !! Axisymmetric label lane (docs/eulerian-label-plan.md): energy
    !! W = 2 pi int w dR dZ of the direct Clebsch representation with the
    !! label chi and g in V0 of a spline section, given profiles F(chi),
    !! G(chi), p(chi) and the winding form about a given axis. Returns W, its
    !! gradient in (chi, g) node values and the Hessian, from the generated
    !! leaf clebsch_energy_axisym by the chain rule.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_eval, only: spline_quadrature_t, spline_geometry_at
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_cell_basis, &
        spline_derham_t_basis
    use fortsparse_csr, only: csr_t
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use milca_eulerian_clebsch_leaves, only: clebsch_energy_axisym, clebsch_energy_3d
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

    public :: clebsch_profiles, clebsch_axisym_system, clebsch_node_gradients, &
              clebsch_3d_system

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


    subroutine clebsch_3d_system(derham, quad, nt, profiles, axis, chi, g, energy, &
                                 gradient, hessian, plane_axis, plane_nodes)
        !! 3D label lane on the section x trigonometric-nodal complex (nt =
        !! 2M + 1 planes, cylindrical chart, phi = 2 pi t). Unknowns chi and g
        !! as node values per plane, section index fastest, ordered (chi, g).
        !! Physical gradient components (R, phi, Z): section gradients times
        !! L_k(t) and N_i L_k'(t)/(2 pi R). The winding form is about the
        !! axis X (independent of phi). Per section cell the local Hessian is
        !! accumulated densely over all quadrature points and planes.
        !! With plane_nodes (V0 node values of (R, Z) on each plane, the same
        !! section topology) the section varies with t: s(xi, t) = sum_k
        !! L_k(t) s_k(xi), so a wall that changes with phi and a core that
        !! follows the axis are both isoparametric. With plane_axis (the
        !! axis X_k on each plane) the winding form is about the axis curve
        !! X(t) = sum_k L_k(t) X_k. The quadrature stays that of derham.
        type(spline_derham_t), intent(in) :: derham
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: nt
        procedure(clebsch_profiles) :: profiles
        real(dp), intent(in) :: axis(2), chi(:), g(:)
        real(dp), intent(out) :: energy
        real(dp), allocatable, intent(out) :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp), intent(in), optional :: plane_axis(:, :), plane_nodes(:, :, :)
        real(dp), allocatable :: v0(:, :, :), gn(:, :, :), du(:, :), gloc(:), hloc(:, :), &
                                 vals(:)
        real(dp) :: dt(2), xt(2), sm(2), fm(2, 2), ft(2, 2), detf, gx(2), v(2), rt(2), &
                    sk(2, nt), gk(2, 2, nt)
        integer :: kk
        integer, allocatable :: d0(:), idx(:), rows(:), cols(:)
        real(dp) :: lt(nt, 2), ht(nt, 2), s(2), jet(2, 6), wq, wt, c, gc(3), gg(3), hv(3), &
                    rho2, f(3), gp(3), pr(3), w, gu(7), h(28), hm(7, 7), t
        real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
        integer :: n0, ntot, nqt, cell, q, p, m, i, k, a, b, nloc, nl, st, nz, col
        type(csc_t) :: mm
        type(fortsparse_status_t) :: info

        n0 = derham%dimension(0)
        ntot = n0*nt
        nqt = 2*nt + 2
        allocate (gradient(2*ntot))
        gradient = 0
        energy = 0
        nz = size(quad%point_offset)*(2*(derham%degree + 1)**2*nt)**2
        allocate (rows(nz), cols(nz), vals(nz))
        nz = 0
        do cell = 1, size(quad%point_offset)
            call clebsch_node_gradients(derham, quad, cell, v0, d0, gn)
            nloc = size(d0)
            nl = nloc*nt
            if (allocated(du)) deallocate (du, gloc, hloc, idx)
            allocate (du(7, 2*nl), gloc(2*nl), hloc(2*nl, 2*nl), idx(2*nl))
            do k = 1, nt
                idx((k - 1)*nloc + 1:k*nloc) = d0 + n0*(k - 1)
                idx(nl + (k - 1)*nloc + 1:nl + k*nloc) = ntot + d0 + n0*(k - 1)
            end do
            gloc = 0
            hloc = 0
            do q = 1, quad%points_per_cell
                p = quad%point_offset(cell) + q - 1
                call spline_geometry_at(derham%space, quad%block(p), quad%xi(:, p), jet, st)
                s = jet(:, 1)
                wq = quad%weight(p)
                if (present(plane_nodes)) then
                    do kk = 1, nt
                        sk(:, kk) = matmul(plane_nodes(:, d0, kk), v0(1, :, q))
                        gk(:, :, kk) = matmul(plane_nodes(:, d0, kk), transpose(gn(:, :, q)))
                    end do
                end if
                do m = 1, nqt
                    t = (m - 1)/real(nqt, dp)
                    call spline_derham_t_basis(nt, t, lt, ht)
                    ! the point, F = ds/ds0, the mesh velocity ds/dt and the axis
                    sm = s
                    fm = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
                    v = 0
                    if (present(plane_nodes)) then
                        sm = matmul(sk, lt(:, 1))
                        fm = 0
                        do kk = 1, nt
                            fm = fm + lt(kk, 1)*gk(:, :, kk)
                        end do
                        v = matmul(sk, lt(:, 2))
                    end if
                    xt = axis
                    dt = 0
                    if (present(plane_axis)) then
                        xt = matmul(plane_axis, lt(:, 1))
                        dt = matmul(plane_axis, lt(:, 2))
                    end if
                    detf = fm(1, 1)*fm(2, 2) - fm(1, 2)*fm(2, 1)
                    ft = reshape([fm(2, 2), -fm(1, 2), -fm(2, 1), fm(1, 1)], [2, 2])/detf
                    rt = sm - xt
                    rho2 = rt(1)**2 + rt(2)**2
                    hv = [-rt(2), (rt(2)*dt(1) - rt(1)*dt(2))/(two_pi*sm(1)), rt(1)]/rho2
                    du = 0
                    do k = 1, nt
                        do i = 1, nloc
                            a = (k - 1)*nloc + i
                            ! F^-T grad_0; d/dphi at fixed (R, Z) removes the mesh velocity
                            gx = matmul(ft, gn(:, i, q))
                            du(1, a) = v0(1, i, q)*lt(k, 1)
                            du(2, a) = gx(1)*lt(k, 1)
                            du(3, a) = (v0(1, i, q)*lt(k, 2) - dot_product(v, gx)*lt(k, 1)) &
                                       /(two_pi*sm(1))
                            du(4, a) = gx(2)*lt(k, 1)
                            du(5, nl + a) = du(2, a)
                            du(6, nl + a) = du(3, a)
                            du(7, nl + a) = du(4, a)
                        end do
                    end do
                    c = dot_product(du(1, 1:nl), chi(idx(1:nl)))
                    gc = matmul(du(2:4, 1:nl), chi(idx(1:nl)))
                    gg = matmul(du(5:7, nl + 1:), g(idx(nl + 1:) - ntot))
                    call profiles(c, f, gp, pr)
                    call clebsch_energy_3d(sm(1), c, gc(1), gc(2), gc(3), gg(1), gg(2), gg(3), &
                                           hv(1), hv(2), hv(3), f(1), f(2), f(3), gp(1), &
                                           gp(2), gp(3), pr(1), pr(2), pr(3), w, gu(1), gu(2), &
                                           gu(3), gu(4), gu(5), gu(6), gu(7), h(1), h(2), &
                                           h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10), &
                                           h(11), h(12), h(13), h(14), h(15), h(16), h(17), &
                                           h(18), h(19), h(20), h(21), h(22), h(23), h(24), &
                                           h(25), h(26), h(27), h(28))
                    wt = wq*abs(detf)*two_pi/nqt
                    energy = energy + wt*w
                    gloc = gloc + wt*matmul(gu, du)
                    k = 0
                    do a = 1, 7
                        do b = a, 7
                            k = k + 1
                            hm(a, b) = h(k)
                            hm(b, a) = h(k)
                        end do
                    end do
                    if (present(hessian)) hloc = hloc + wt*matmul(transpose(du), matmul(hm, du))
                end do
            end do
            gradient(idx) = gradient(idx) + gloc
            if (.not. present(hessian)) cycle
            do col = 1, 2*nl
                do i = 1, 2*nl
                    if (hloc(i, col) == 0) cycle
                    nz = nz + 1
                    rows(nz) = idx(i)
                    cols(nz) = idx(col)
                    vals(nz) = hloc(i, col)
                end do
            end do
        end do
        if (.not. present(hessian)) return
        call csc_from_triplet(2*ntot, 2*ntot, cols(:nz), rows(:nz), vals(:nz), mm, info)
        if (info%code /= 0) error stop 'Clebsch 3D Hessian triplets'
        hessian%nrow = 2*ntot
        hessian%ncol = 2*ntot
        hessian%row_ptr = mm%col_ptr
        hessian%col = mm%row_idx
        hessian%val = mm%val
    end subroutine clebsch_3d_system

end module milca_eulerian_clebsch
