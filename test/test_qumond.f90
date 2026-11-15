program test_qumond
    ! QUMOND (G = 1) against derived limits, on a fine polar mesh and on one
    ! with twice its spacing:
    ! 1. spherical symmetry (Plummer sphere, M = a = a0 = 1): the curl term
    !    vanishes and g = nu(g_N/a0) g_N holds exactly, so Phi(r) - Phi(0) is
    !    the integral of nu g_N; for nu_1 and nu_2, second order in the mesh;
    ! 2. deep MOND far from the mass: v^4/(G M a0) = (y nu_1(y))^2/y
    !    = 1 + sqrt(y) + O(y), y = G M/(r^2 a0), so v^4 -> G M a0;
    ! 3. Newtonian limit: with nu = 1 the projection returns Phi_N; for
    !    g_N >> a0 the correction Phi - Phi_N is the integral of (nu - 1) g_N,
    !    about a0 r;
    ! 4. Kuzmin disk: outside the midplane g_N is a point-mass field about the
    !    image point, so nu g_N is curl-free and the algebraic relation is
    !    exact (Brada and Milgrom 1995) although the system is not spherical;
    ! 5. two separated Plummer spheres: nu g_N has curl, and the solution
    !    departs from it by a mesh-independent amount, while the flux identity
    !    (the curl part carries no flux through a closed surface) and the
    !    projection identity ||g||^2 + ||g - gt||^2 = ||gt||^2 hold.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_qumond
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    type(qumond_t) :: fine, coarse
    type(qumond_field_t) :: f
    real(dp) :: pair_d = 1, pair_a = 0.25_dp, a0 = 1
    real(dp) :: err(2), dep_sphere(2), dep_disk(2), dep_pair(2), circ(2)
    integer :: n_nu, k

    call qumond_setup(fine, 0.02_dp, 1e4_dp, 270, 64)
    call require(fine%ok, 'fine mesh')
    call qumond_setup(coarse, 0.02_dp, 1e4_dp, 135, 32)
    call require(coarse%ok, 'coarse mesh')

    print '(a)', '1. Plummer sphere: max |Phi - int nu g_N dr| / max |Phi|, r < 50, coarse and fine'
    do n_nu = 1, 2
        call qumond_solve(coarse, 1.0_dp, a0, n_nu, plummer, f)
        call require(f%ok, 'Plummer solve')
        err(1) = sphere_error(coarse, f, 50.0_dp, n_nu)
        call qumond_solve(fine, 1.0_dp, a0, n_nu, plummer, f)
        call require(f%ok, 'Plummer solve')
        err(2) = sphere_error(fine, f, 50.0_dp, n_nu)
        print '(a, i1, a, 2es10.3)', '   nu_', n_nu, ':', err
        call require(err(2) < 1.5e-3_dp .and. err(1)/err(2) > 3, 'spherical algebraic relation')
    end do

    print '(a)', '2. deep MOND, Plummer sphere:'
    call qumond_solve(fine, 1.0_dp, a0, 1, plummer, f)
    print '(a)', '      r      G M/(r^2 a0)   v^4/(G M a0)   exact      1 + sqrt(y)'
    do k = 1, 3
        call deep(fine, f, minloc(abs(fine%r - 10.0_dp**k), dim=1), k == 3)
    end do

    print '(a)', '3. Newtonian limit, Plummer sphere, r < 10:'
    a0 = 1e-30_dp
    call qumond_solve(fine, 1.0_dp, a0, 1, plummer, f)
    call newton_limit(fine, f)
    a0 = 1e-4_dp
    call qumond_solve(fine, 1.0_dp, a0, 1, plummer, f)
    call newton_limit(fine, f)

    a0 = 1
    print '(a)', '4. Kuzmin disk: max |Phi - Phi_exact| / max |Phi_exact|, r < 30, coarse and fine'
    call qumond_solve(coarse, 1.0_dp, a0, 1, zero, f, surface=kuzmin)
    err(1) = kuzmin_error(coarse, f, 30.0_dp)
    call qumond_solve(fine, 1.0_dp, a0, 1, zero, f, surface=kuzmin)
    call require(f%ok, 'Kuzmin solve')
    err(2) = kuzmin_error(fine, f, 30.0_dp)
    print '(a, 2es10.3)', '  ', err
    call require(err(2) < 1.5e-3_dp .and. err(1)/err(2) > 3, 'Kuzmin disk algebraic relation')

    print '(a)', '5. circulation of gt around R, z in [0.5, 1.5] over int |gt| dl, Kuzmin and two masses'
    circ = [circulation(kuzmin_gn), circulation(pair_gn)]
    print '(a, 2es11.3)', '  ', circ
    call require(abs(circ(1)) < 1e-12_dp .and. abs(circ(2)) > 1e-2_dp, 'curl of the algebraic field')
    print '(a)', '   departure ||g - gt|| / ||gt||, r < 10, coarse and fine:'
    call departures(coarse, dep_sphere(1), dep_disk(1), dep_pair(1))
    call departures(fine, dep_sphere(2), dep_disk(2), dep_pair(2))
    print '(a, 2es10.3)', '   Plummer sphere:', dep_sphere
    print '(a, 2es10.3)', '   Kuzmin disk:   ', dep_disk
    print '(a, 2es10.3)', '   two masses:    ', dep_pair
    call require(dep_sphere(2) < 2e-4_dp .and. dep_disk(2) < 5e-3_dp .and. &
        dep_sphere(1)/dep_sphere(2) > 1.8_dp .and. dep_disk(1)/dep_disk(2) > 1.8_dp, &
        'curl-free cases converge to the algebraic field')
    call require(abs(dep_pair(1)/dep_pair(2) - 1) < 0.1_dp .and. dep_pair(2) > 5*dep_disk(2), &
        'two masses depart from the algebraic field')
    call qumond_solve(fine, 1.0_dp, a0, 1, pair, f)
    call require(f%ok, 'two-mass solve')
    call pair_identities(fine, f)

    print '(a)', 'QUMOND satisfies its limits'

contains

    function zero(r_cyl, z) result(rho)
        real(dp), intent(in) :: r_cyl, z
        real(dp) :: rho
        rho = 0*(r_cyl + z)
    end function zero

    function plummer(r_cyl, z) result(rho)
        real(dp), intent(in) :: r_cyl, z
        real(dp) :: rho
        rho = 3/(4*pi)*(1 + r_cyl**2 + z**2)**(-2.5_dp)
    end function plummer

    function plummer_gn(r) result(g)
        real(dp), intent(in) :: r
        real(dp) :: g
        g = r/(r**2 + 1)**1.5_dp
    end function plummer_gn

    !> Two Plummer spheres of mass 1/2 and radius pair_a at z = +-pair_d.
    function pair(r_cyl, z) result(rho)
        real(dp), intent(in) :: r_cyl, z
        real(dp) :: rho
        rho = 0.5_dp*3/(4*pi*pair_a**3)*((1 + (r_cyl**2 + (z - pair_d)**2)/pair_a**2)**(-2.5_dp) &
            + (1 + (r_cyl**2 + (z + pair_d)**2)/pair_a**2)**(-2.5_dp))
    end function pair

    function pair_gn(x) result(g)
        real(dp), intent(in) :: x(2)
        real(dp) :: g(2), s(2)
        integer :: k
        g = 0
        do k = -1, 1, 2
            s = x - [0.0_dp, k*pair_d]
            g = g - 0.5_dp*s/(sum(s**2) + pair_a**2)**1.5_dp
        end do
    end function pair_gn

    !> Kuzmin disk, M = a = 1.
    function kuzmin(r_cyl, z) result(sigma)
        real(dp), intent(in) :: r_cyl, z
        real(dp) :: sigma
        sigma = 1/(2*pi)*(r_cyl**2 + 1 + 0*z)**(-1.5_dp)
    end function kuzmin

    !> Point mass at (0, -sign(z)) seen from either side of the disk.
    function kuzmin_gn(x) result(g)
        real(dp), intent(in) :: x(2)
        real(dp) :: g(2), s(2)
        s = [x(1), abs(x(2)) + 1]
        g = -s/norm2(s)**3
        g(2) = sign(g(2), -x(2))
    end function kuzmin_gn

    !> y nu_n(y), written independently of the module.
    function nu_y(n, y) result(f)
        integer, intent(in) :: n
        real(dp), intent(in) :: y
        real(dp) :: f
        f = (0.5_dp*y**n + sqrt(0.25_dp*y**(2*n) + y**n))**(1.0_dp/n)
    end function nu_y

    function algebraic(gn_of, x) result(g)
        procedure(pair_gn) :: gn_of
        real(dp), intent(in) :: x(2)
        real(dp) :: g(2), gn(2)
        gn = gn_of(x)
        g = a0*nu_y(1, norm2(gn)/a0)*gn/norm2(gn)
    end function algebraic

    function sphere_error(q, f, r_out, n_nu) result(err)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        real(dp), intent(in) :: r_out
        integer, intent(in) :: n_nu
        real(dp) :: err, exact(0:q%nr), u(16), w(16), r0, r1, scale
        integer :: i, j, k
        call gauss_legendre(16, u, w)
        exact(0) = 0
        do i = 1, q%nr
            r0 = merge(0.0_dp, q%r(max(i - 1, 1)), i == 1)
            r1 = q%r(i)
            exact(i) = exact(i - 1)
            do k = 1, 16
                exact(i) = exact(i) + w(k)*(r1 - r0)/2* &
                    a0*nu_y(n_nu, plummer_gn(r0 + (u(k) + 1)*(r1 - r0)/2)/a0)
            end do
        end do
        err = 0
        scale = 0
        do i = 1, q%nr
            if (q%r(i) > r_out) exit
            scale = max(scale, abs(exact(i)))
            do j = 0, q%nth
                err = max(err, abs(f%phi(qumond_node(q, i, j)) - exact(i)))
            end do
        end do
        err = err/scale
    end function sphere_error

    !> Rotation speed from the secant slope of Phi along the equator.
    function speed4(q, f, i) result(v4)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        integer, intent(in) :: i
        real(dp) :: v4
        v4 = ((q%r(i) + q%r(i + 1))/2*(f%phi(qumond_node(q, i + 1, q%nth/2)) - &
            f%phi(qumond_node(q, i, q%nth/2)))/(q%r(i + 1) - q%r(i)))**2/a0
    end function speed4

    subroutine deep(q, f, i, far)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        integer, intent(in) :: i
        logical, intent(in) :: far
        real(dp) :: r, v4, exact
        r = (q%r(i) + q%r(i + 1))/2
        v4 = speed4(q, f, i)
        exact = r**2*a0*nu_y(1, plummer_gn(r)/a0)**2
        print '(f10.1, es12.3, 3f13.7)', r, 1/(r**2*a0), v4, exact, 1 + sqrt(1/(r**2*a0))
        call require(abs(v4 - exact) < 1.5e-3_dp, 'deep-MOND rotation speed')
        if (far) call require(abs(v4 - 1) < 1.5e-3_dp, 'flat rotation curve v^4 = G M a0')
    end subroutine deep

    subroutine newton_limit(q, f)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        real(dp) :: err, scale, newton, u(16), w(16), r0, r1, exact(0:q%nr)
        integer :: i, j, k
        ! exact correction int_0^r (nu - 1) g_N ds
        call gauss_legendre(16, u, w)
        exact(0) = 0
        do i = 1, q%nr
            r0 = merge(0.0_dp, q%r(max(i - 1, 1)), i == 1)
            r1 = q%r(i)
            exact(i) = exact(i - 1)
            do k = 1, 16
                associate (gn => plummer_gn(r0 + (u(k) + 1)*(r1 - r0)/2))
                    exact(i) = exact(i) + w(k)*(r1 - r0)/2*(a0*nu_y(1, gn/a0) - gn)
                end associate
            end do
        end do
        err = 0
        scale = 0
        newton = 0
        do i = 1, q%nr
            if (q%r(i) > 10) exit
            scale = max(scale, exact(i))
            do j = 0, q%nth
                k = qumond_node(q, i, j)
                err = max(err, abs(f%phi(k) - (f%phi_n(k) - f%phi_n(1)) - exact(i)))
                newton = max(newton, abs(f%phi_n(k) + 1/sqrt(q%r(i)**2 + 1)))
            end do
        end do
        print '(a, es8.1, a, es10.3, a, es10.3, a, es10.3)', '   a0 =', a0, &
            ': max |Phi_N + G M/sqrt(r^2 + a^2)| =', newton, &
            ', max |(Phi - Phi_N) - int (nu - 1) g_N dr| =', err, ', max correction =', scale
        call require(newton < 1.5e-3_dp, 'Newtonian Plummer potential')
        ! the floor is the round-off of the two sparse solves
        call require(err < max(1e-7_dp, 1e-2_dp*scale), 'Newtonian limit')
    end subroutine newton_limit

    subroutine departures(q, dsphere, ddisk, dpair)
        type(qumond_t), intent(inout) :: q
        real(dp), intent(out) :: dsphere, ddisk, dpair
        type(qumond_field_t) :: f
        call qumond_solve(q, 1.0_dp, a0, 1, plummer, f)
        dsphere = departure(q, f, 10.0_dp)
        call qumond_solve(q, 1.0_dp, a0, 1, zero, f, surface=kuzmin)
        ddisk = departure(q, f, 10.0_dp)
        call qumond_solve(q, 1.0_dp, a0, 1, pair, f)
        dpair = departure(q, f, 10.0_dp)
    end subroutine departures

    function departure(q, f, r_out) result(d)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        real(dp), intent(in) :: r_out
        real(dp) :: d, num, den
        integer :: t
        num = 0
        den = 0
        do t = 1, q%ntri
            if (norm2(q%centroid(:, t)) > r_out) cycle
            num = num + q%weight(t)*sum((f%g(:, t) - f%g_alg(:, t))**2)
            den = den + q%weight(t)*sum(f%g_alg(:, t)**2)
        end do
        d = sqrt(num/den)
    end function departure

    !> Phi(x) - Phi(0) = int_a^s nu(G M/(s^2 a0)) G M/s^2 ds, s the distance
    !> from the image point.
    function kuzmin_error(q, f, r_out) result(err)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        real(dp), intent(in) :: r_out
        real(dp) :: err, scale, u(32), w(32), s, exact
        integer :: i, j, k, node
        call gauss_legendre(32, u, w)
        err = 0
        scale = 0
        do i = 1, q%nr
            if (q%r(i) > r_out) exit
            do j = 0, q%nth
                node = qumond_node(q, i, j)
                s = norm2([q%x(1, node), abs(q%x(2, node)) + 1])
                exact = 0
                do k = 1, 32
                    exact = exact + w(k)*(s - 1)/2*a0*nu_y(1, 1/((1 + (u(k) + 1)*(s - 1)/2)**2*a0))
                end do
                scale = max(scale, abs(exact))
                err = max(err, abs(f%phi(node) - exact))
            end do
        end do
        err = err/scale
    end function kuzmin_error

    !> Circulation of the algebraic field around the square R, z in
    !> [0.5, 1.5], relative to the line integral of its magnitude.
    function circulation(gn_of) result(c)
        procedure(pair_gn) :: gn_of
        real(dp) :: c, total, corner(2, 5), u(32), w(32), x(2), dl(2), g(2)
        integer :: side, k
        call gauss_legendre(32, u, w)
        corner = reshape([0.5_dp, 0.5_dp, 1.5_dp, 0.5_dp, 1.5_dp, 1.5_dp, 0.5_dp, 1.5_dp, &
            0.5_dp, 0.5_dp], [2, 5])
        c = 0
        total = 0
        do side = 1, 4
            dl = corner(:, side + 1) - corner(:, side)
            do k = 1, 32
                x = corner(:, side) + (u(k) + 1)/2*dl
                g = algebraic(gn_of, x)
                c = c + w(k)/2*dot_product(g, dl)
                total = total + w(k)/2*norm2(g)*norm2(dl)
            end do
        end do
        c = c/total
    end function circulation

    subroutine pair_identities(q, f)
        type(qumond_t), intent(in) :: q
        type(qumond_field_t), intent(in) :: f
        real(dp) :: g2, h2, t2, flux_h, flux_a, p0(2), p1(2), n(2), u(8), w(8), x(2), len, v4
        integer :: t, i, j, k
        g2 = sum(q%weight*sum(f%g**2, dim=1))
        h2 = sum(q%weight*sum((f%g - f%g_alg)**2, dim=1))
        t2 = sum(q%weight*sum(f%g_alg**2, dim=1))
        print '(a, es10.3)', '   (||g||^2 + ||g - gt||^2) / ||gt||^2 - 1 =', (g2 + h2)/t2 - 1
        call require(abs((g2 + h2)/t2 - 1) < 1e-10_dp, 'projection identity')
        ! flux of the QUMOND field (mean of the elements on both sides of
        ! the polygonal sphere) against that of the analytic algebraic field
        call gauss_legendre(8, u, w)
        do t = 1, 2
            i = minloc(abs(q%r - merge(2, 10, t == 1)), dim=1)
            flux_h = 0
            flux_a = 0
            do j = 1, q%nth
                p0 = q%x(:, qumond_node(q, i, j - 1))
                p1 = q%x(:, qumond_node(q, i, j))
                len = norm2(p1 - p0)
                n = [p1(2) - p0(2), p0(1) - p1(1)]/len
                if (dot_product(n, p0 + p1) < 0) n = -n
                flux_h = flux_h + pi*(p0(1) + p1(1))*len*0.5_dp* &
                    dot_product(f%g(:, qumond_cell(q, i, j, .false.)) + &
                    f%g(:, qumond_cell(q, i - 1, j, .true.)), n)
                do k = 1, 8
                    x = p0 + (u(k) + 1)/2*(p1 - p0)
                    flux_a = flux_a + w(k)/2*len*2*pi*x(1)*dot_product(algebraic(pair_gn, x), n)
                end do
            end do
            print '(a, f5.2, a, 2f13.7, es10.2)', '   flux through r =', q%r(i), &
                ' of g and of gt:', flux_h, flux_a, flux_h/flux_a - 1
            call require(abs(flux_h/flux_a - 1) < 5e-4_dp, 'flux identity')
        end do
        i = minloc(abs(q%r - 1e3_dp), dim=1)
        v4 = speed4(q, f, i)
        print '(a, f7.1, a, f11.7)', '   v^4/(G M a0) at r =', q%r(i), ':', v4
        call require(abs(v4 - 1) < 1.5e-3_dp, 'two-mass deep-MOND rotation speed')
    end subroutine pair_identities

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_qumond
