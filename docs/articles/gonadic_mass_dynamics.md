# Gonadic mass dynamics

## Deriving the PDE for gonadic mass

Standard mizer describes the consumer number density \\N(w,t)\\ by the
McKendrick-von Foerster equation with an added diffusion term: \\
\frac{\partial N}{\partial t} + \frac{\partial}{\partial w}\left( g\\N -
\frac{1}{2}\frac{\partial(d\\N)}{\partial w} \right) = -\mu\\N,
\tag{1}\\ where \\g(w,t)\\ is the somatic growth rate of an individual
of size \\w\\, \\d(w,t)\\ is its diffusion coefficient, and \\\mu(w,t)\\
is its total mortality rate. The bracketed term inside the spatial
derivative is the total consumer flux \\J_N(w,t) = J_N^{\text{adv}} +
J_N^{\text{diff}}\\, consisting of the advective growth flux \\g\\N\\
and the diffusive flux \\-\frac{1}{2}\frac{\partial(d\\N)}{\partial
w}\\.

Instead of the number density \\N\\, we could describe the population by
the somatic biomass density \\B = w\\N\\. In the advection-only case
(\\d=0\\), differentiating \\B = w\\N\\ gives: \\\begin{split}
\frac{\partial B}{\partial t} &= w\frac{\partial N}{\partial t} =
-w\frac{\partial}{\partial w}\left(g\\\frac{B}{w}\right)-\mu\\B,\\
&=-\frac{\partial}{\partial w}\left(g\\B\right)+g\\\frac{B}{w}-\mu\\B.
\end{split}\\ This is intuitive: in addition to the advective transport
term \\-\frac{\partial}{\partial w}(g B)\\ and the mortality loss \\-\mu
B\\, there is a distributed source term \\g B/w = g N\\ representing
somatic growth occurring within the size class.

### Total gonadic biomass density \\Q(w,t)\\

Now consider the gonadic mass carried by the population. Let \\Q(w,t)\\
denote the total gonadic biomass density across size \\w\\ at time
\\t\\, and let \\G(w,t)\\ be the gonadic investment rate of an
individual of size \\w\\ (i.e. the energy allocated to reproduction).

In mizer, the energy available for growth and reproduction \\E_r\\ is
split into somatic growth \\g = E_r(1-\psi)\\ and reproduction
investment \\G= E_r\psi = g\\\psi / (1 - \psi)\\, where \\\psi(w)\\ is
the proportion of energy allocated to reproduction. In mizer code,
\\G(w,t)\\ corresponds to `e_repro(w,t)`.

Because gonadic mass is an extensive quantity carried by the individuals
\\N\\, \\Q(w,t)\\ undergoes the exact same advective transport and
stochastic size diffusion as \\N(w,t)\\. In addition, \\Q(w,t)\\
receives a distributed source from reproductive investment \\G\\N\\,
loses mass to natural/fishing mortality at rate \\\mu\\Q\\, and releases
mass during spawning at a mass-specific release rate \\r(w,t)\\: \\
\frac{\partial Q}{\partial t} + \frac{\partial}{\partial w}\left( g\\Q-
\frac{1}{2}\frac{\partial(d\\Q)}{\partial w} \right) = G\\N - \mu\\Q-
r\\Q. \tag{2}\\

The release rate \\r(w,t)\\ encodes the seasonality of reproduction. All
terms in [Equation 2](#eq-PDEQ) are species-specific, so there is one
equation for each species \\i\\.

### Per-capita gonadic mass \\q(w,t) = Q/ N\\

Instead of the total gonadic biomass density \\Q(w,t)\\, one might
consider the average gonadic mass \\q(w,t)\\ of an individual fish of
size \\w\\: \\ q(w,t) = \frac{Q(w,t)}{N(w,t)}. \\ Let us derive the
evolution equation for \\q(w,t)\\ by applying the quotient rule
\\\frac{\partial q}{\partial t} = \frac{1}{N}\frac{\partial Q}{\partial
t} - \frac{Q}{N^2}\frac{\partial N}{\partial t}\\.

Substituting [Equation 2](#eq-PDEQ) and [Equation 1](#eq-PDEN) into the
quotient rule: \\\begin{split} \frac{\partial q}{\partial t} &=
\frac{1}{N}\left( -\frac{\partial J_Q}{\partial w} + G\\N - \mu\\Q-
r\\Q\right) - \frac{q}{N}\left( -\frac{\partial J_N}{\partial w} -
\mu\\N \right)\\ &= G- r\\q- \frac{1}{N}\left( \frac{\partial
J_Q}{\partial w} - q\frac{\partial J_N}{\partial w} \right), \end{split}
\tag{3}\\ where the mortality terms \\-\mu q\\ and \\+\mu q\\ cancel out
exactly.

Let us evaluate the spatial flux divergence term \\\frac{1}{N}\left(
\frac{\partial J_Q}{\partial w} - q\frac{\partial J_N}{\partial w}
\right)\\:

1.  **Advective contribution:** \\\frac{\partial(g\\qN)}{\partial w} -
    q\frac{\partial(g\\N)}{\partial w} = g\\N\frac{\partial q}{\partial
    w}.\\ Dividing by \\N\\ gives \\g\\\frac{\partial q}{\partial w}\\.

2.  **Diffusive contribution:** Using the product rule on
    \\\frac{\partial^2(d\\qN)}{\partial w^2} -
    q\frac{\partial^2(d\\N)}{\partial w^2}\\:
    \\\frac{\partial^2(d\\qN)}{\partial w^2} -
    q\frac{\partial^2(d\\N)}{\partial w^2} = d\\N\frac{\partial^2
    q}{\partial w^2} + 2\frac{\partial(d\\N)}{\partial w}\frac{\partial
    q}{\partial w}.\\ Dividing by \\N\\ and multiplying by \\-1/2\\
    gives \\-\frac{1}{2}d\frac{\partial^2 q}{\partial w^2} -
    \frac{1}{N}\frac{\partial(d\\N)}{\partial w}\frac{\partial
    q}{\partial w}\\.

Combining these terms gives the full PDE for individual gonadic mass
\\q(w,t)\\: \\ \frac{\partial q}{\partial t} + g\\\frac{\partial
q}{\partial w} - \frac{1}{2}d\\\frac{\partial^2 q}{\partial w^2} -
\frac{1}{N}\frac{\partial(d\\N)}{\partial w}\frac{\partial q}{\partial
w} = G- r\\q. \tag{4}\\

### The deterministic limit (\\d = 0\\) vs. diffusion (\\d \> 0\\)

In the absence of diffusion (\\d = 0\\), [Equation 4](#eq-wg-diff)
simplifies to: \\ \frac{\partial q}{\partial t} + g\\\frac{\partial
q}{\partial w} = G- r\\q. \tag{5}\\ Because \\\frac{dq(w(t),t)}{dt} =
\frac{\partial q}{\partial t} + \frac{\partial q}{\partial
w}\frac{dw}{dt}\\ and \\\frac{dw}{dt} = g\\, [Equation 5](#eq-wg)
reduces along the deterministic growth trajectory of an individual to
the ordinary differential equation: \\ \frac{dq}{dt} = G- r\\q.
\tag{6}\\ Without diffusion, the evolution of per-capita gonadic mass
\\q\\ is completely decoupled from the abundance \\N(w,t)\\ and
mortality \\\mu(w,t)\\.

However, when diffusion is present (\\d \> 0\\),
[Equation 4](#eq-wg-diff) reveals several key difficulties:

1.  **Coupling to population gradients:** The effective drift velocity
    of \\q\\ acquires an additional term
    \\-\frac{1}{N}\frac{\partial(d\\N)}{\partial w}\\. This arises
    because individual fish undergo random walks in size; the average
    gonadic mass at size \\w\\ depends on the relative abundance of fish
    diffusing upwards from smaller sizes versus downwards from larger
    sizes.
2.  **Singularities in unpopulated size classes:** Wherever \\N(w,t) \to
    0\\ (such as sizes below maturity, sizes above \\w\_{\max}\\, or
    gaps between discrete cohorts), the term
    \\\frac{1}{N}\frac{\partial(d\\N)}{\partial w}\\ becomes
    indeterminate or singular (\\0/0\\).
3.  **Non-conservative form:** Equation [Equation 4](#eq-wg-diff) is not
    in divergence (flux) form \\\frac{\partial q}{\partial t} +
    \frac{\partial J}{\partial w} = \dots\\, making it impossible to
    apply standard finite-volume conservation schemes and flux limiters
    (such as van Leer TVD) without violating global conservation of
    total gonadic biomass.

**Conclusion:** For size-spectrum models with diffusion and higher-order
spatial schemes, the total gonadic biomass density \\Q(w,t)\\ is the
proper fundamental state variable. Once \\Q(w,t)\\ and \\N(w,t)\\ are
known, the per-capita gonadic mass is simply recovered as \\q(w,t) =
Q(w,t) / N(w,t)\\ for all sizes where \\N(w,t) \> 0\\.

## Reproduction rate

The mass-specific release rate \\r(w,t)\\ determines the total rate of
egg production \\R_p(t)\\ (density-independent reproduction rate
\\R\_{di}\\): \\ R_p(t) = \frac{\epsilon}{2\\w_0} \int
Q(w,t)\\r(w,t)\\dw = \frac{\epsilon}{2\\w_0} \int
N(w,t)\\q(w,t)\\r(w,t)\\dw, \tag{7}\\ where \\w_0 = w\_{\min}\\ is the
egg size, \\\epsilon\\ is the reproductive efficiency (`erepro`), and
the factor \\1/2\\ accounts for the female proportion of the population.

The release rate \\r(w,t)\\ is typically periodic with a period of one
year. It reaches a large peak once a year during the spawning season so
that the accumulated gonadic biomass is released to produce the seasonal
cohort.

## Numerical implementation

To discretise the PDE [Equation 2](#eq-PDEQ) for total gonadic biomass
density \\Q(w,t)\\, we use the finite-volume framework and higher-order
schemes described in mizer’s numerical details vignette.

### The weight grid and finite-volume averages

We discretise the logarithmic weight grid \\w_1, w_2, \ldots, w_K\\ with
grid ratio \\\beta = w\_{j+1}/w_j = 10^{\Delta x}\\. The \\j\\-th bin is
\\\[w_j, w\_{j+1})\\ with width \\\Delta w_j = w\_{j+1} - w_j\\ and
geometric centre \\w_j^c = \sqrt{w_j w\_{j+1}}\\.

In the finite-volume method, discrete state values represent **bin
averages**: \\ Q_j = \frac{1}{\Delta w_j}\int\_{w_j}^{w\_{j+1}}
Q(w)\\dw, \qquad N_j = \frac{1}{\Delta w_j}\int\_{w_j}^{w\_{j+1}}
N(w)\\dw. \\ Integrating [Equation 2](#eq-PDEQ) over bin \\j\\ and
dividing by \\\Delta w_j\\ gives: \\ \frac{\partial Q_j}{\partial t} +
\frac{J\_{Q, j+1} - J\_{Q, j}}{\Delta w_j} = S\_{Q, j} - (\mu_j +
r_j)\\Q_j, \tag{8}\\ where \\J\_{Q, j}\\ is the total flux of gonadic
biomass across the lower bin boundary \\w_j\\, and \\S\_{Q, j}\\ is the
bin-averaged source from reproductive investment.

### Discretisation of spatial fluxes

The flux \\J\_{Q, j} = J\_{Q, j}^{\text{adv}} + J\_{Q,
j}^{\text{diff}}\\ at boundary \\w_j\\ consists of:

1.  **Advective flux:** \\ J\_{Q, j}^{\text{adv}} = g_j\\Q(w_j), \\
    where \\g_j = g(w_j)\\ is the point growth velocity at the boundary
    \\w_j\\. Under the **first-order upwind scheme**, the boundary
    density is taken from the bin below: \\ J\_{Q, j}^{\text{adv}} =
    g_j\\Q\_{j-1}. \\ Under the **second-order scheme**
    (`second_order_w$flux`), the boundary density is reconstructed using
    the slope-limited formula: \\ J\_{Q, j}^{\text{adv}} = g_j\left\[
    Q\_{j-1} + \frac{1}{2}\chi_j(Q_j - Q\_{j-1}) \right\], \\ where
    \\\chi_j = 1\\ gives the centred second-order reconstruction and
    \\\chi_j = \chi(r_j)\\ with the **van Leer** limiter \\ \chi(r) =
    \frac{r + \|r\|}{1 + \|r\|}, \qquad r_j = \frac{Q\_{j-1} -
    Q\_{j-2}}{Q_j - Q\_{j-1}} \\ gives the total-variation diminishing
    (TVD) scheme that avoids spurious oscillations while preserving
    non-negativity (\\Q\ge 0\\).

2.  **Diffusive flux:** The diffusive flux is evaluated by central
    differencing the bin-averaged products \\d_j Q_j\\ between adjacent
    bin centres: \\ J\_{Q, j}^{\text{diff}} =
    -\frac{1}{2}\\\frac{d_j\\Q_j - d\_{j-1}\\Q\_{j-1}}{w_j^c -
    w\_{j-1}^c}. \\

### Bin averaging of source and sink terms (`second_order_w$bin_average`)

- **Reproduction source term:** The total investment in bin \\j\\ is
  \\\int\_{w_j}^{w\_{j+1}} G(w)\\N(w)\\dw\\. When `bin_average = TRUE`,
  \\G(w) = e\_{\text{repro}}(w)\\ is trapezoidal bin-averaged, giving
  the discrete source: \\ S\_{Q, j} = \overline{G}\_j\\N_j =
  \overline{e\_{\text{repro}}}\_j\\N_j. \\
- **Loss terms:** The total sink rate is \\(\mu_j + r_j)\\, where
  \\\mu_j\\ is the bin-averaged total mortality and \\r_j\\ is the
  bin-averaged release rate (equal to \\r(t)\\ when the release rate is
  size-independent).

### Boundary conditions

- **Lower boundary (\\j = j\_{\min}\\):** Newborn larvae enter the
  smallest size class with zero gonadic mass. Thus, there is no influx
  of gonadic mass across the egg boundary: \\ J\_{Q, j\_{\min}} = 0. \\
  The sub-diagonal coefficient \\A\_{j\_{\min}}\\ is set to 0, and the
  lower-face advection/diffusion terms drop out of the diagonal
  coefficient \\B\_{j\_{\min}}\\.
- **Upper boundary (\\j = j\_{\max}\\):** The grid is truncated at
  \\w\_{\max}\\. The super-diagonal coefficient \\C\_{j\_{\max}}\\ is
  set to 0, decoupling the retained spectrum from sizes above
  \\w\_{\max}\\, and \\Q_j\\ is held at zero for \\j \> j\_{\max}\\.

### Semi-implicit linear system

Using semi-implicit time discretisation with time step \\\Delta t\\,
evaluating \\Q\\ at time \\t+1\\ and rates at time \\t\\ yields the
tridiagonal linear system: \\ A_j\\Q\_{j-1}^{t+1} + B_j\\Q_j^{t+1} +
C_j\\Q\_{j+1}^{t+1} = S_j, \tag{9}\\ where the right-hand side source
vector is: \\ S_j = Q_j^t + \Delta t\\S\_{Q, j}^t = Q_j^t + \Delta
t\\\overline{e\_{\text{repro}}}\_{j}^t\\N_j^t. \\ The tridiagonal
coefficients \\A_j, B_j, C_j\\ are identical to those built by mizer’s
transport operator `get_transport_coefs()`, with total mortality
replaced by \\\mu_j + r_j^t\\ and the recruitment influx set to 0.

### Higher-order time-stepping schemes

1.  **Predictor-Corrector method:**
    - **Predict:** Solve [Equation 9](#eq-tridiagonal-Q) using
      start-of-step rates \\r^t\\ to obtain predicted end-of-step
      densities \\\hat{N}^{t+1}\\ and \\\hat{Q}^{t+1}\\.
    - **Midpoint rates:** Evaluate provisional rates \\\hat{r}^{t+1}\\
      and form second-order midpoint rates \\r^{t+1/2} =
      \frac{1}{2}(r^t + \hat{r}^{t+1})\\.
    - **Correct:** Advance \\Q\\ using a Crank-Nicolson corrector with
      midpoint rates: \\ \frac{Q_j^{t+1} - Q_j^t}{\Delta t} +
      \frac{1}{2\Delta w_j}\left\[ (J\_{Q, j+1}^{t+1} - J\_{Q,
      j}^{t+1}) + (J\_{Q, j+1}^t - J\_{Q, j}^t) \right\] = S\_{Q,
      j}^{t+1/2} - \frac{1}{2}(\mu_j^{t+1/2} + r_j^{t+1/2})(Q_j^{t+1} +
      Q_j^t). \\
2.  **TR-BDF2 method:** An L-stable, two-stage second-order method with
    stage fraction \\\gamma = 2 - \sqrt{2}\\:
    - **Trapezoidal (TR) stage** over \\\[t, t+\gamma\Delta t\]\\: \\
      \frac{Q^{t+\gamma} - Q^t}{\gamma\Delta t} = \frac{1}{2}\left(
      LQ^t + LQ^{t+\gamma} \right) + S_Q^{t+1/2}. \\
    - **BDF2 stage** over the whole step: \\ Q^{t+1} =
      \frac{1}{\gamma(2-\gamma)}Q^{t+\gamma} -
      \frac{(1-\gamma)^2}{\gamma(2-\gamma)}Q^t +
      \frac{1-\gamma}{2-\gamma}\Delta t\left( LQ^{t+1} + S_Q^{t+1/2}
      \right). \\ Both stages share the exact same tridiagonal matrix
      \\(I - \alpha\Delta t\\L)\\ with \\\alpha = 1 - 1/\sqrt{2}\\,
      providing second-order accuracy while damping stiff high-frequency
      modes without the numerical ringing of Crank-Nicolson.

### Relationship with the legacy per-capita scheme

In the first-order upwind limit without diffusion (\\d = 0\\, \\\chi =
0\\), [Equation 9](#eq-tridiagonal-Q) for \\Q\\ reduces to: \\
Q_j^{t+1}\left(1 + g_j^t\frac{\Delta t}{\Delta w_j} + (\mu_j^t +
r_j^t)\Delta t\right) = Q\_{j-1}^{t+1}\left(g\_{j-1}^t\frac{\Delta
t}{\Delta w_j}\right) + Q_j^t + e\_{\text{repro}, j}^t\\N_j^t\\\Delta t.
\\ Dividing by \\N_j^{t+1}\\ and using the stationary approximation
\\N_j^{t+1} \approx N_j^t\\ and \\g\_{j-1} N\_{j-1} \approx g_j N_j -
\mu_j N_j \Delta w_j\\ reproduces the legacy per-capita recurrence
relation: \\ q_j^{t+1}\left(1 + g_j^t\frac{\Delta t}{\Delta w_j} +
r_j^t\Delta t\right) = q\_{j-1}^{t+1}\left(g_j^t\frac{\Delta t}{\Delta
w_j}\right) + q_j^t + G_j^t\\\Delta t. \\ Solving for \\Q(w,t)\\
directly removes all stationary approximations, handles size diffusion
\\d(w) \> 0\\ rigorously, and integrates seamlessly with mizer’s
higher-order spatial limiters and time-stepping solvers.
