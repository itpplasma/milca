# Full local electric drift

The full local trajectory includes the speed and pitch changes associated with
radial electrostatic work. This note gives the discretization used by the
optional full-electric surface transport model.

Let `theta` and `zeta` be Boozer angles, `b` the field magnitude divided by
its surface Fourier constant, and `G` and `I` the covariant toroidal and
poloidal magnetic-field components. Let `x = v/v_t` be speed divided by the
thermal speed, `xi = v_parallel/v` the pitch, and `E` the electric angular
rate in inverse seconds. Define

\[
D_\perp=\partial_\theta-(I/G)\partial_\zeta,\quad
g=\frac{D_\perp b}{2b},\quad
F=e^{-x^2}.
\tag{1}
\]

For the perturbation \(f=Fh\), the electric terms in the local characteristic
are [Landreman et al., 2014](https://arxiv.org/abs/1312.6058)

\[
L_E
=E\left[D_\perp+gx(1+\xi^2)\partial_x+
              g\xi(1-\xi^2)\partial_\xi\right].
\tag{2}
\]

The phase-volume density is \(J=x^2/b^2\). Together with the streaming and
mirror terms, this trajectory conserves the magnetic moment \(x^2(1-\xi^2)/b\).
The electric and streaming-plus-mirror terms each satisfy phase-volume
conservation. The electric identity is

\[
\begin{aligned}
D_\perp(JE)&=-4JgE,\\
\partial_x\big[J Egx(1+\xi^2)\big]&=3JgE(1+\xi^2),\\
\partial_\xi\big[J Eg\xi(1-\xi^2)\big]&=JgE(1-3\xi^2),\\
\nabla\cdot(J\dot{\boldsymbol z}_E)&=0.
\end{aligned}
\tag{3}
\]

Writing \(u\) for a test function, integration by parts gives

\[
\int J u L_E(Fh)\,d\theta\,d\zeta\,dx\,d\xi
=-\int J Fh L_Eu\,d\theta\,d\zeta\,dx\,d\xi .
\tag{4}
\]

The angle boundary terms vanish by periodicity. The speed flux vanishes at
\(x=0\) and at infinity because of \(x^3F\). The pitch flux vanishes at
\(\xi=\pm1\). This form includes the Maxwellian derivative and gives zero
density moment for the constant test function.

The angle derivative is Fourier collocation. The speed basis is the collision
operator's continuous Gauss–Lobatto–Legendre basis with its polynomial tail.
For test index \(i\) and trial index \(j\), the speed matrices are
\(\int x^2F\phi_i\phi_j\,dx\) and
\(\int x^3F\phi_i'\phi_j\,dx\). Pitch integrals use Gauss–Legendre quadrature
with Legendre test and trial polynomials.

The transport solver uses this action with its existing speed-row scaling and
fixes the coefficient at zero speed to zero. The per-speed surface-average
constraint is retained. The preconditioner uses the modal diagonal of the
full-electric angular advection and the collision blocks. The rate tangent
uses the same electric action as the primal solve. The separate particle and
pressure constraints and tangential magnetic drifts are open.
