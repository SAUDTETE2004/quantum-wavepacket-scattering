# From scattering states to a time-dependent experiment

## 1. Coordinates, units and initial state

The barrier is centered at the origin, with $a=L/2$:

$$V(x)=\begin{cases}U_0,&|x|<a,\\0,&|x|>a.\end{cases}$$

The code uses $\hbar=m=1$. Choose a physical length $\ell$ to restore units:
$x_{phys}=\ell x$, $E_{phys}=\hbar^2 E/(m\ell^2)$ and
$t_{phys}=m\ell^2t/\hbar$. Here $m$ in these conversions is the physical mass.

The normalized incident packet is

$$\psi(x,0)=(2\pi\sigma_x^2)^{-1/4}
\exp[-(x-x_0)^2/(4\sigma_x^2)+ik_0(x-x_0)].$$

Thus $\Delta x=\sigma_x$, $\Delta k=1/(2\sigma_x)$,
$p_0=\hbar k_0$ and $E_0=\hbar^2k_0^2/(2m)$. $E_0$ is the carrier energy;
the mean kinetic energy is $E_0+\hbar^2/(8m\sigma_x^2)$.
Under $\psi=(2\pi)^{-1/2}\int f(k)e^{ikx}dk$,

$$f(k)=(2\sigma_x^2/\pi)^{1/4}e^{-\sigma_x^2(k-k_0)^2}e^{-ikx_0}.$$

Using $e^{-i(k-k_0)x_0}$ instead changes only a constant phase. **Dropping
the linear phase would incorrectly drop the initial position.** The Gaussian
has a negative-momentum weight $\Phi(-2k_0\sigma_x)$; left-incidence theory
conditions on $k>0$. This distinction is negligible for the supplied cases.
Also require $-a-x_0\gg\sigma_x$ so the initial overlap with the barrier is tiny.

## 2. Stationary barrier amplitudes

For left incidence and $k>0$, use global coordinates:

$$\phi_k(x)=\begin{cases}
e^{ikx}+r(k)e^{-ikx},&x<-a,\\
A e^{iqx}+B e^{-iqx},&|x|<a,\\
t(k)e^{ikx},&x>a,
\end{cases}\qquad q=\sqrt{k^2-2mU_0/\hbar^2}.$$

Take $q=i\kappa$ below the barrier. At threshold the interior solution is
$A+Bx$, the continuous $q\to0$ limit. Match both the function and its derivative
at each edge. A convenient way to eliminate the interior coefficients is

$$\binom{\phi(a)}{\phi'(a)}=
\begin{pmatrix}\cos(qL)&\sin(qL)/q\\-q\sin(qL)&\cos(qL)\end{pmatrix}
\binom{\phi(-a)}{\phi'(-a)}.$$

Insert the exterior waves and solve the resulting two equations for $r,t$:

$$D=\cos(qL)-i\frac{k^2+q^2}{2k}\frac{\sin(qL)}q,$$
$$t=\frac{e^{-ikL}}D,\qquad
r=\frac{-i(k^2-q^2)\sin(qL)/(2kq)\ e^{-ikL}}D.$$

The reflected phase here is for a *centered* barrier. Translating a barrier by
$b$ multiplies $r$ by $e^{2ikb}$ while leaving $t$ unchanged. Using the phase of
a barrier on $[0,L]$ in the centered formulas would shift reflected centers.

Equal exterior wave numbers imply $T_k=|t|^2$, $R_k=|r|^2$, $T_k+R_k=1$.
For $E>U_0$,

$$T_k=\left[1+\frac{U_0^2\sin^2(qL)}{4E(E-U_0)}\right]^{-1};$$

below threshold replace the ratio by
$U_0^2\sinh^2(\kappa L)/[4E(U_0-E)]$. At threshold,
$T_k=[1+(kL/2)^2]^{-1}$. Above-barrier resonances occur at $qL=n\pi$.
The implementation uses scaled exponentials to avoid overflow for large
$\kappa L$, and evaluates $\sin(qL)/q\to L$ at threshold.

### Connection to the parity derivation

Right exterior parity states are proportional to $\cos(kx+\delta_e)$ and
$\sin(kx+\delta_o)$. Continue these to the left with even/odd parity; do not
use the same right-exterior formula through the barrier. In this convention,
$t=(e^{2i\delta_e}+e^{2i\delta_o})/2$ and
$r=(e^{2i\delta_e}-e^{2i\delta_o})/2$.
Parity does not require the sine coefficient on the right exterior to vanish,
and matching does not quantize continuum $k$.

## 3. Probability for a packet

The outgoing transmitted amplitude is $f(k)t(k)$; the reflected one is
$f(k)r(k)$ with spatial factor $e^{-ikx}$. Parseval's identity gives

$$T=\int_0^\infty |f_+(k)|^2|t(k)|^2dk,\qquad
R=\int_0^\infty |f_+(k)|^2|r(k)|^2dk=1-T,$$

where $f_+$ is normalized on the positive half-line. Only when $|t|^2$ varies
little across $\Delta k$ may we replace this by $|t(k_0)|^2$.

In the simulation define $P_{right}(t)=\int_a^\infty|\psi|^2dx$,
$P_{left}(t)=\int_{-\infty}^{-a}|\psi|^2dx$ and $P_{barrier}(t)$.
These sum to one. $P_{left}$ initially includes the incident packet and is
not a reflection probability during collision. After separation,
$P_{right}\to T$ and $P_{left}\to R$. A finite periodic box has no literal
$t\to\infty$ scattering limit: measure its plateau before boundary return.
The continuity equation gives $dP_{right}/dt=j(a)$ on the infinite line.

## 4. Gaussian asymptotics, including position and phase curvature

Write $k=k_0+u$ and $t=|t|e^{i\theta}$. With constant magnitude and
$\theta=\theta_0+\theta'_0u+\theta''_0u^2/2$, the integral is Gaussian:

$$\psi_T\propto A^{-1/2}\exp[-(x-x_0-v_0t+\theta'_0)^2/(4A)],$$
$$A=\sigma_x^2+i\left(\frac{\hbar t}{2m}-\frac{\theta''_0}{2}\right),
\qquad v_0=\frac{\hbar k_0}m.$$

The omitted carrier phase has unit modulus. For $A=\alpha+i\beta$,
$\operatorname{Re}(1/A)=\alpha/(\alpha^2+\beta^2)$, hence

$$x_T=x_0+v_0t-\theta'_0,\quad
\Delta x_T^2=\sigma_x^2+\frac{(\hbar t/(2m)-\theta''_0/2)^2}{\sigma_x^2}.$$

For reflection replace $\theta$ by $\arg r$ and use
$x_R=-x_0-v_0t+(\arg r)'_0$; the width formula is unchanged.
The relative transmission group delay is
$\tau=\hbar\,d\theta/dE=\theta'_0/v_0$. This is a shift relative to free
propagation, not a measurement of time spent inside the barrier.

The free-energy expansion is exactly quadratic and its quadratic term grows
with time. Dropping phase curvature requires the separate condition
$|\theta''_0|(\Delta k)^2\ll1$; keeping free dispersion does not establish it.
For second-order phase theory, neglected cubic terms must be small over the
occupied bandwidth. Curvature can cause initial focusing: the width need not
increase monotonically for every phase-modified packet, although eventually
it grows linearly. These are exterior, asymptotic formulas, not global
solutions inside the barrier. Extending a fitted Gaussian there is not a
calculation of probability leaking through the barrier.

## 5. The barrier also filters momentum

Expand $\log c(k)=\log c_0+(g_1+i\theta_1)u+
(g_2+i\theta_2)u^2/2$, for $c=t$ or $r$. In code units define

$$\alpha=\sigma_x^2-g_2/2,\quad\beta=(t-\theta_2)/2,\quad
B=sx-x_0-k_0t+\theta_1,\quad s=+1\ (T),\ -1\ (R).$$

The integral is proportional to
$A^{-1/2}\exp[(g_1+iB)^2/(4A)]$, with $A=\alpha+i\beta$.
Its maximum occurs at $B=\beta g_1/\alpha$, giving

$$x_c=s[x_0+k_0t-\theta_1+\beta g_1/\alpha],\qquad
\Delta x^2=(\alpha^2+\beta^2)/\alpha,$$
$$v_c=s[k_0+g_1/(2\alpha)].$$

This includes a small shift in outgoing mean momentum. Require $\alpha>0$
and smooth, nonzero amplitude across the packet. Near a resonance zero of
$r$, or strong energy filtering, a Gaussian need not be a good approximation.
We fit the density as $M\exp[-(x-\mu)^2/(2s_x^2)]/(\sqrt{2\pi}s_x)$,
report $M,\mu,s_x$, peak density and relative L2 residual, and compare a
history of centers and widths with the formula above. Analytic density peaks
use the spectral channel mass. A density fit does not verify the complex
wavefunction's phase; this project's Gaussian claim concerns the envelope.

## References

- [Richard Fitzpatrick, Square Potential Barrier](https://phys.libretexts.org/Bookshelves/Quantum_Mechanics/Introductory_Quantum_Mechanics_%28Fitzpatrick%29/04%3A_One-Dimensional_Potentials/4.02%3A_Square_Potential_Barrier): stationary matching and probabilities.
- [Algorithm Archive, Split-Operator Method](https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html): position/momentum splitting implementation.
- [ETH computational quantum physics lecture notes](https://edu.itp.phys.ethz.ch/fs11/cqp/cqp.pdf): numerical time propagation and split operators.

The derivation above fixes its own centered-barrier phase convention; compare
conventions before copying phase-delay expressions from another source.
