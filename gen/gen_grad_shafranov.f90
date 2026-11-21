program gen_grad_shafranov
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_expr, only: expr_t, sym, rat, log, sqrt, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_products, only: jvp
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_precision, only: PRECISION_REAL32, PRECISION_REAL64
    use milca_gen_support, only: leaf_writer_t, fortsym_revision
    implicit none
    type(leaf_writer_t) :: w, wb, ws
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    type(expr_t) :: r, r0, mu0, pp, ff, pr, pz, f, b0, k
    type(expr_t) :: psi, br, bz, bt
    type(expr_t) :: width, derivative, slope
    type(expr_t) :: hrr, hrz, hzz, nr, nz, jump_r, jump_z
    type(expr_t) :: x, a0, a1, a2, a3, axis, cubic, primitive, coordinate
    type(expr_t) :: rank_a, rank_b, rank_t
    type(expr_t) :: tt, level, seg_v(0:3), seg_node(0:3), seg_basis(0:3), seg_value
    type(expr_t) :: seg_slope
    integer :: i, ik
    type(expr_t) :: dpp, dff, app, aff, source, df, af, tangent(1)
    type(expr_t) :: pres, energy, g_psi, g_r, g_z, h_psi
    type(expr_t) :: rho, slope_r, a1c, a2c, sig, phie, dpds, dids, dsdu, ppn, ffn
    type(expr_t) :: hh, y0, y1, m0, m1, c0, c1, c2, c3, herm, tt2, h0, h1, ym, yp, mslope
    type(expr_t) :: xe, x1e, x2e, y1e, y2e, lin, mend, moment, uproj, vproj, wproj
    engine = make_native_engine()
    r = sym('r'); r0 = sym('r0')
    mu0 = sym('mu0')
    pp = sym('pprime'); ff = sym('ffprime')
    pr = sym('psi_r'); pz = sym('psi_z'); f = sym('f')
    b0 = sym('b0'); k = sym('k')
    width = sym('width'); derivative = sym('derivative')
    slope = sym('slope')
    hrr = sym('h_rr'); hrz = sym('h_rz'); hzz = sym('h_zz')
    nr = sym('n_r'); nz = sym('n_z')
    jump_r = sym('jump_r'); jump_z = sym('jump_z')
    x = sym('x'); axis = sym('axis')
    a0 = sym('a0'); a1 = sym('a1')
    a2 = sym('a2'); a3 = sym('a3')
    identity = engine%zero_test(diff(width*derivative + width**2*slope/2, width) &
        - derivative - width*slope)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS profile primitive unproved'
    cubic = a0 + x*(a1 + x*(a2 + x*a3))
    primitive = x*(a0 + x*(a1/2 + x*(a2/3 + x*a3/4)))
    identity = engine%zero_test(diff(primitive, x) - cubic)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS cubic primitive unproved'
    identity = engine%zero_test(diff(cubic, x) - (a1 + 2*a2*x + 3*a3*x*x))
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS cubic source slope unproved'
    ! Physical pprime and FFprime remain derivatives with respect to psi,
    ! even when their polynomial argument is normalized poloidal flux.
    psi = sym('psi')
    coordinate = 1 - psi/axis
    identity = engine%zero_test(diff(-axis*(a0*coordinate &
        + a1*coordinate**2/2 + a2*coordinate**3/3 + a3*coordinate**4/4), &
        psi) - (a0 + a1*coordinate + a2*coordinate**2 &
        + a3*coordinate**3))
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS normalized primitive unproved'
    ! J=B-u*v^T, B*z=r, B*w=u: a=v^T*z, b=v^T*w, with 1-b nonzero.
    ! Then J*(z+w*t)-r=u*(t-a-b*t); swapping u/v proves the transpose.
    rank_a = sym('rank_a'); rank_b = sym('rank_b')
    rank_t = rank_a/(1-rank_b)
    identity = engine%zero_test(rank_t-rank_a-rank_b*rank_t)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS rank-one response unproved'
    ! Smooth-cell source batch: the source and its two directional derivatives
    ! (Newton response d/dpsi and axis response d/dpsi_axis of the profile law).
    ! The tangent expressions are checked against the symbolic JVP.
    dpp = sym('dpprime'); dff = sym('dffprime')
    app = sym('apprime'); aff = sym('affprime')
    source = mu0*r*pp + ff/r
    df = mu0*r*dpp + dff/r
    af = mu0*r*app + aff/r
    tangent = jvp([source], [pp, ff], [dpp, dff])
    identity = engine%zero_test(tangent(1) - df)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS source tangent unproved'
    tangent = jvp([source], [pp, ff], [app, aff])
    identity = engine%zero_test(tangent(1) - af)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS source axis tangent unproved'
    ! Energy density of the Grad-Shafranov energy on spline spaces,
    ! w = |grad psi|^2/(2 mu0 R) - R p(psi) - F(psi)^2/(2 mu0 R): its
    ! derivative along psi, by the chain rule with dp/dpsi = pprime and
    ! d(F^2/2)/dpsi = ffprime, is -(R pprime + ffprime/(mu0 R)) = -source/mu0;
    ! the second derivative follows with the source tangent.
    pres = sym('pressure')
    energy = (pr**2 + pz**2)/(2*mu0*r) - r*pres - f**2/(2*mu0*r)
    g_psi = -(r*pp + ff/(mu0*r))
    tangent = jvp([energy], [pres, f], [pp, ff/f])
    identity = engine%zero_test(tangent(1) - g_psi)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS energy chain rule unproved'
    identity = engine%zero_test(g_psi + source/mu0)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS energy source unproved'
    g_r = diff(energy, pr)
    g_z = diff(energy, pz)
    h_psi = -(r*dpp + dff/(mu0*r))
    tangent = jvp([g_psi], [pp, ff], [dpp, dff])
    identity = engine%zero_test(tangent(1) - h_psi)
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS energy Hessian unproved'
    ! Current-profile mode (gs_free_current): coarea weights along a ray,
    ! profile update from dI/du = a1 p' + a2 FF'/mu0, cubic Hermite pieces.
    rho = sym('rho'); slope_r = sym('slope')
    a1c = sym('a1'); a2c = sym('a2'); sig = sym('sigma'); phie = sym('phi_edge')
    dpds = sym('dpds'); dids = sym('dids')
    dsdu = f*a2c/phie
    ppn = -sig*dpds*dsdu
    ffn = mu0*(sig*dids*dsdu - a1c*ppn)/a2c
    identity = engine%zero_test(a1c*ppn + a2c*ffn/mu0 - sig*dids*dsdu)
    if (identity%verdict /= VERDICT_TRUE) error stop 'current balance of the profile update unproved'
    hh = sym('h'); y0 = sym('y0'); y1 = sym('y1'); m0 = sym('m0'); m1 = sym('m1')
    tt2 = sym('t')
    c0 = y0; c1 = m0
    c2 = (3*(y1 - y0)/hh - 2*m0 - m1)/hh
    c3 = (m0 + m1 - 2*(y1 - y0)/hh)/hh**2
    herm = c0 + tt2*(c1 + tt2*(c2 + tt2*c3))
    identity = engine%zero_test(subs(herm, tt2, hh) - y1)
    if (identity%verdict /= VERDICT_TRUE) error stop 'Hermite end value unproved'
    identity = engine%zero_test(subs(diff(herm, tt2), tt2, hh) - m1)
    if (identity%verdict /= VERDICT_TRUE) error stop 'Hermite end slope unproved'
    h0 = sym('h0'); h1 = sym('h1'); ym = sym('ym'); yp = sym('yp')
    mslope = (h0**2*(yp - y0) + h1**2*(y0 - ym))/(h0*h1*(h0 + h1))
    ! exact for every quadratic y(t) = y0 + m0 t + c t^2 about the middle node
    identity = engine%zero_test(subs(subs(mslope, yp, y0 + m0*h1 + k*h1**2), ym, &
        y0 - m0*h0 + k*h0**2) - m0)
    if (identity%verdict /= VERDICT_TRUE) error stop 'three-point slope unproved'
    ! End slope of the Hermite law: the slope at the end node of the quadratic
    ! through the end values and the slope m1 at the neighbour.
    mend = 2*(y1 - y0)/hh - m1
    identity = engine%zero_test(subs(subs(mend, y1, y0 + m0*hh + (m1 - m0)*hh/2), m1, m1) &
        - m0)
    if (identity%verdict /= VERDICT_TRUE) error stop 'Hermite end slope unproved'
    xe = sym('x'); x1e = sym('x1'); x2e = sym('x2'); y1e = sym('y1'); y2e = sym('y2')
    lin = y1e + (y2e - y1e)*(xe - x1e)/(x2e - x1e)
    identity = engine%zero_test(subs(lin, xe, x2e) - y2e)
    if (identity%verdict /= VERDICT_TRUE) error stop 'linear extrapolation unproved'
    moment = sig*(a1c*pp + a2c*ff/mu0)
    identity = engine%zero_test(subs(subs(moment, pp, ppn), ff, ffn) - sig**2*dids*dsdu)
    if (identity%verdict /= VERDICT_TRUE) error stop 'current moment inverse unproved'
    call w%open('gen_grad_shafranov', &
                'src/generated/milca_grad_shafranov_leaves.f90')
    call w%line('! Generated by gen/gen_grad_shafranov.f90; do not edit.')
    call w%line('! FortSym revision: '//fortsym_revision)
    call w%line('! SI, CCW phi, psi flux per radian, B=grad psi x grad phi+F grad phi.')
    call w%line('module milca_grad_shafranov_leaves')
    call w%line('implicit none')
    call w%line('contains')
    call w%emit('gs_source', ['r      ', 'mu0    ', 'pprime ', 'ffprime'], &
                ['source'], [mu0*r*pp + ff/r], ad=.true., interval=.true.)
    br = -pz/r; bz = pr/r; bt = f/r
    call w%emit('gs_field', ['r    ', 'psi_r', 'psi_z', 'f    '], &
                ['b_r  ', 'b_z  ', 'b_phi'], [br, bz, bt], ad=.true.)
    ! Classical current inside a smooth element and distributional face current.
    call w%emit('gs_current_density', &
        ['r     ', 'mu0   ', 'psi_r ', 'h_rr  ', 'h_zz  '], &
        ['j_phi'], [-(hrr+hzz-pr/r)/(mu0*r)], ad=.true.)
    call w%emit('gs_current_sheet', &
        ['r     ', 'mu0   ', 'jump_r', 'jump_z', 'n_r   ', 'n_z   '], &
        ['k_phi'], [-(jump_r*nr+jump_z*nz)/(mu0*r)], ad=.true.)
    ! Positive-outward psi, elliptic stationary minimum and physical CCW phi.
    call w%emit('gs_smooth_axis_q', &
        ['r   ', 'f   ', 'h_rr', 'h_rz', 'h_zz'], &
        ['q'], [f/(r*sqrt(hrr*hzz-hrz**2))], ad=.true.)
    ! A straight cylindrical model; r is distance from its axis, not tokamak R.
    bz = b0/(1 + k**2*r**2); bt = k*r*bz
    psi = r0*b0*log(1 + k**2*r**2)/(2*k)
    call w%emit('gs_gold_hoyle', ['r ', 'r0', 'b0', 'k '], &
        ['b_z  ', 'b_t  ', 'psi  '], [bz, bt, psi], ad=.true.)
    call w%emit('gs_segment_primitive', ['width     ', 'derivative', 'slope     '], &
        ['integral'], [width*derivative + width**2*slope/2], ad=.true.)
    call w%emit('gs_cubic', ['x ', 'a0', 'a1', 'a2', 'a3'], &
        ['value'], [cubic], ad=.true.)
    call w%emit('gs_cubic_primitive', ['x ', 'a0', 'a1', 'a2', 'a3'], &
        ['value'], [primitive], ad=.true.)
    call w%emit('gs_normalized_flux', ['x   ', 'axis'], &
        ['s'], [1 - x/axis], ad=.true.)
    ! Source-cut topology: discriminants of the flux polynomial in the inner
    ! coordinate. Their outer roots include closed contours inside a cell.
    call w%emit('gs_cut_discriminant', ['a0', 'a1', 'a2', 'a3'], &
        ['value'], [18*a3*a2*a1*a0 - 4*a2**3*a0 + a2**2*a1**2 &
        - 4*a3*a1**3 - 27*a3**2*a0**2])
    call w%emit('gs_cut_quadratic_discriminant', ['a0', 'a1', 'a2'], &
        ['value'], [a1**2 - 4*a2*a0])
    call w%emit('gs_spline_energy', &
        ['r       ', 'mu0     ', 'psi_r   ', 'psi_z   ', 'pressure', 'f       ', &
        'pprime  ', 'ffprime ', 'dpprime ', 'dffprime'], &
        ['w      ', 'g_psi  ', 'g_r    ', 'g_z    ', 'h_psi  ', 'h_grad '], &
        [energy, g_psi, g_r, g_z, h_psi, diff(g_r, pr)])
    call w%emit('gs_ray_coarea_weights', ['r    ', 'rho  ', 'slope'], ['w1', 'w2'], &
        [r*rho/slope_r, rho/(r*slope_r)])
    call w%emit('gs_current_profile_update', &
        ['sigma   ', 'mu0     ', 'f       ', 'a1      ', 'a2      ', 'phi_edge', 'dpds    ', &
        'dids    '], ['pprime ', 'ffprime'], [ppn, ffn])
    call w%emit('gs_hermite_cubic', ['h ', 'y0', 'y1', 'm0', 'm1'], &
        ['c0', 'c1', 'c2', 'c3'], [c0, c1, c2, c3])
    call w%emit('gs_three_point_slope', ['h0', 'h1', 'ym', 'y0', 'yp'], ['slope'], [mslope])
    ! L2 projection of a sampled field onto a spline space: w = (u - v)^2/2
    uproj = sym('u'); vproj = sym('v')
    wproj = (uproj - vproj)**2/2
    call w%emit('gs_l2_projection', ['u', 'v'], ['w     ', 'g_u   ', 'h_uu  '], &
        [wproj, diff(wproj, uproj), diff(diff(wproj, uproj), uproj)])
    call w%emit('gs_hermite_end_slope', ['h ', 'y0', 'y1', 'm1'], ['slope'], [mend])
    call w%emit('gs_linear_extrapolation', ['x ', 'x1', 'x2', 'y1', 'y2'], ['y'], [lin])
    call w%emit('gs_current_moment', ['sigma  ', 'mu0    ', 'a1     ', 'a2     ', 'pprime ', &
        'ffprime'], ['didu'], [moment])
    call w%line('end module milca_grad_shafranov_leaves')
    call w%close()
    ! Cubic Lagrange segment on the nodes 0, 1/3, 2/3, 1 (flux samples along a
    ! P3 edge): value at t of (sum_i L_i(t) (v_i - level)) and the slopes of
    ! the interpolant at t = 0, 1/2, 1. Consumed by milca_triangle_basis.
    tt = sym('t'); level = sym('level')
    do i = 0, 3
        seg_v(i) = sym('v'//achar(48 + i))
        seg_node(i) = rat(int(i, int64), 3_int64)
    end do
    do i = 0, 3
        seg_basis(i) = rat(1_int64, 1_int64)
        do ik = 0, 3
            if (ik == i) cycle
            seg_basis(i) = seg_basis(i)*(tt - seg_node(ik)) &
                /(seg_node(i) - seg_node(ik))
        end do
    end do
    seg_value = seg_basis(0)*(seg_v(0) - level)
    do i = 1, 3
        seg_value = seg_value + seg_basis(i)*(seg_v(i) - level)
    end do
    identity = engine%zero_test(subs(seg_value, tt, seg_node(2)) - (seg_v(2) - level))
    if (identity%verdict /= VERDICT_TRUE) error stop 'GS segment cubic node unproved'
    seg_slope = seg_basis(1)*(seg_v(1) - seg_v(0)) + seg_basis(2)*(seg_v(2) - seg_v(0)) &
        + seg_basis(3)*(seg_v(3) - seg_v(0))
    seg_slope = diff(seg_slope, tt)
    call ws%open('gen_grad_shafranov', &
                 'src/generated/milca_segment_leaves.f90')
    ws%simplify_roots = .true.
    call ws%line('! Generated by gen/gen_grad_shafranov.f90; do not edit.')
    call ws%line('! FortSym revision: '//fortsym_revision)
    call ws%line('! Cubic Lagrange segment leaves (nodes 0, 1/3, 2/3, 1).')
    call ws%line('module milca_segment_leaves')
    call ws%line('implicit none')
    call ws%line('contains')
    call ws%emit('gs_segment_cubic_value', &
        ['t    ', 'level', 'v0   ', 'v1   ', 'v2   ', 'v3   '], ['value'], [seg_value])
    call ws%emit('gs_segment_cubic_slopes', ['v0', 'v1', 'v2', 'v3'], &
        ['g0    ', 'g_half', 'g1    '], &
        [subs(seg_slope, tt, seg_node(0)), &
        subs(seg_slope, tt, rat(1_int64, 2_int64)), &
        subs(seg_slope, tt, seg_node(3))])
    call ws%line('end module milca_segment_leaves')
    call ws%close()
    call wb%open('gen_grad_shafranov', &
                 'src/generated/milca_grad_shafranov_batch.f90')
    call wb%line('! Generated by gen/gen_grad_shafranov.f90; do not edit.')
    call wb%line('! FortSym revision: '//fortsym_revision)
    call wb%line('! SIMD batch drivers (FortSym FY1) for the smooth-cell load.')
    call wb%emit_batch('gs_source', ['r      ', 'mu0    ', 'pprime ', 'ffprime'], &
                       ['source'], [source], ['mu0'])
    call wb%emit_batch('gs_source_tangents', &
        ['r       ', 'mu0     ', 'pprime  ', 'ffprime ', 'dpprime ', 'dffprime', &
        'apprime ', 'affprime'], ['source', 'dsrc  ', 'asrc  '], &
        [source, df, af], ['mu0'])
    ! Profile law: cubic piece with its slope, and the normalized flux with its
    ! two partial derivatives (psi and axis flux).
    call wb%emit_batch('gs_cubic_slope', ['x ', 'a0', 'a1', 'a2', 'a3'], &
                       ['value', 'slope'], [cubic, diff(cubic, x)], &
                       ['a0', 'a1', 'a2', 'a3'])
    call wb%emit_batch('gs_normalized_flux_jets', ['x   ', 'axis'], &
                       ['s       ', 'ds_dx   ', 'ds_daxis'], &
                       [1 - x/axis, diff(1 - x/axis, x), diff(1 - x/axis, axis)], &
                       ['axis'])
    call wb%close()
    ! Experimental CPU benchmark seam: use the upstream emitter's precision.
    ! No candidate or estimator policy is changed by these opt-in leaves.
    call wb%open('gen_grad_shafranov', &
                 'src/generated/milca_gs_cpu_precision.f90')
    call wb%emit_batch('gs_cpu_fp64', &
        ['r       ', 'mu0     ', 'pprime  ', 'ffprime ', 'dpprime ', 'dffprime', &
        'apprime ', 'affprime'], ['source', 'dsrc  ', 'asrc  '], &
        [source, df, af], ['mu0'], precision=PRECISION_REAL64)
    call wb%emit_batch('gs_cpu_fp32', &
        ['r       ', 'mu0     ', 'pprime  ', 'ffprime ', 'dpprime ', 'dffprime', &
        'apprime ', 'affprime'], ['source', 'dsrc  ', 'asrc  '], &
        [source, df, af], ['mu0'], precision=PRECISION_REAL32)
    call wb%close()
end program gen_grad_shafranov
