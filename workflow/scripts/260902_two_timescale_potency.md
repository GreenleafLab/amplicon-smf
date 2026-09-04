---
title: "Two-Timescale Potency: Splitting Potency into Instantaneous and Integrated Components"
author: "Design discussion notes"
date: "2026-09-02"
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
  - \usepackage{booktabs}
  - \usepackage{microtype}
---

# Scope

This note develops a model in which **potency** (Nature 2024, Fig. 3d) is decomposed into an
*instantaneous* component --- how much a TF bound **right now, on this molecule** promotes promoter
activity --- and an *integrated* component --- how much the molecule's **time-averaged** binding
history promotes activity. A single mixing parameter $w \in [0,1]$ sets the split.

The motivation is a specific objection to an earlier proposal (in
`260717_ising_intuition.md` §7) that potency could be read off as an equilibrium coupling
$J_\text{prom}$ between promoter state and same-molecule TF occupancy. **That proposal is wrong for
this system, for two independent reasons**, both visible in Fig. 3:

1. The functional form is wrong (Fig. 3d): an equilibrium energy coupling predicts a *convex*
   occupancy--activity curve; the data are straight-to-concave (§2.1).
2. It has no representation of history (Fig. 3h): promoter activity depends on the number of
   *available* sites and not only on the number *currently bound* (§2.2).

The model below fixes both, reproduces the published potency exactly as a limiting case, and makes
$w$ measurable from data already in hand.

---

# 1. Notation

| symbol | meaning |
|---|---|
| $K$ | TetO copy number of the construct (0--8, or 0--9 for opJS5) |
| $n$ | TF occupancy **on this molecule, in this snapshot** (v5's `n_tf_teto`) |
| $\langle n\rangle$ | mean occupancy for this molecule's *condition* $(K, \text{dox}, \text{AD})$ = Fig. 3c/3i |
| $P$ | probability this molecule's promoter is active (v5 region-anchored promoter call) |
| $\bar P$ | condition-level mean of $P$ = "fraction of promoters active", Fig. 3b/3d |
| $p_0$ | basal activity at zero occupancy (intercept) |
| $A$ | odds of being active, $A=P/(1-P)$; used only for the hyperbolic-link alternative |
| $k_\text{pot}$ | potency (units: drive per bound TF) |
| $w$ | integrated fraction of potency ($0$ = purely instantaneous, $1$ = purely historical) |
| $\tau_\text{int}$ | promoter integration timescale |
| $\tau_c$ | correlation (residence) timescale of TF occupancy on one molecule |

---

# 2. What Fig. 3 constrains

## 2.1 Panel 3d fixes the *link function*: line-like, definitely not energy-additive

Panel 3d plots $\bar P$ against $\langle n\rangle$: $\bar P$ runs from $\approx 0.03$ at
$\langle n\rangle = 0$ to $\approx 0.47$ at $\langle n\rangle \approx 5$.

Consider three candidate couplings.

**Energy-additive (Ising / Boltzmann).** Each bound TF adds a fixed free energy stabilizing the
active state, so the *log-odds* are linear:
$$
\operatorname{logit}\bar P = h_\text{prom} + J\,\langle n\rangle .
$$
Fit to the endpoints: $\operatorname{logit}(0.03) = -3.48$, $\operatorname{logit}(0.47) = -0.12$,
so $J \approx 3.36/5 \approx 0.67$. Evaluate at the midpoint $\langle n\rangle = 2.5$:
$\operatorname{logit}\bar P = -1.80 \Rightarrow \bar P \approx 0.14$. A straight line between the
endpoints would give $0.25$. **So this model predicts pronounced upward (convex) curvature ---
the midpoint should sit far *below* the chord.**

**Rate-additive (the paper's Fig. 3f additive activation model).** Each bound TF independently adds
to the promoter's opening *rate*, $k_\text{on} \propto \text{TF}_\text{occ}$, giving
$$
\bar P = \frac{1}{1 + \dfrac{k_\text{off}}{k_\text{on}\!\cdot\!\text{TF}_\text{occ}}}
\;=\; \frac{A}{1+A}, \qquad A = k_0 + k_\text{pot}\langle n\rangle ,
$$
i.e. the *odds* (not the log-odds) are linear in occupancy. This is hyperbolic in
$\langle n\rangle$: mildly **concave**, essentially straight while $A \lesssim 1$, and it is what
`compute_potency.py`'s `saturation_model_mechanistic` already implements.

**Identity (a straight line).** $\bar P = p_0 + k_\text{pot}\langle n\rangle$, with no saturation
over the observed range.

Ranking these against the data:

| link | $\bar P$ vs $\langle n\rangle$ | verdict |
|---|---|---|
| energy-additive (logistic) | strongly **convex** | **rejected** --- predicts $\approx0.14$ at mid-occupancy where the data sit near $0.25$--$0.3$ |
| rate-additive (hyperbolic) | mildly concave | fits ($r^2=0.96$), but *worse than a straight line* |
| identity (linear) | straight | **best empirical fit**; adopted here |

Two conclusions, at different confidence levels.

**Robust:** the curve is *not convex*, which rejects the energy-additive/Boltzmann coupling
regardless of how the concave-vs-straight question resolves. This is the conclusion that kills the
earlier $J_\text{prom}$ proposal, and it does not depend on distinguishing the other two links.

**Empirical:** a straight line fits better than the hyperbola, so **this note uses the identity
link** and treats eventual saturation as outside the measured range. Nothing below depends on
reproducing the hyperbola exactly; the requirement is only that the relation be *line-like*, which
it is. Conveniently, the identity link also makes the central result of §4 *exact* rather than
approximate.

> **The linearity itself remains the open mechanistic puzzle**, exactly as in the original analysis.
> A probability cannot be linear in a driving variable forever, so linearity means the system sits
> far from saturation across the entire accessible range --- which is a statement worth explaining,
> not an explanation. Diagnostic for future ADs: a *stronger* AD pushes $\bar P$ higher and should
> therefore be the first place linearity breaks. **If a new AD reaches high $\bar P$ and is still
> straight, that is a real result**; if it bends, the bend's shape identifies the link function that
> the current ADs are too weak to reveal. Either outcome is informative, so this is worth checking
> per-AD rather than pooling.

## 2.2 Panel 3h requires *two* occupancy terms

Panel 3h stratifies "fraction of promoters active" by the number of TFs bound (0--5) as a function
of TetO number. Two effects are simultaneously present:

- At fixed $K$, activity rises with $n$ (the curves are ordered by shade).
- **At fixed $n$, activity rises with $K$** (each curve slopes upward). A molecule with 4 of 8 sites
  bound is more likely active than one with 4 of 4 bound.

The second effect is the one that no single-snapshot equilibrium model can produce. In equilibrium,
conditioning on the configuration is sufficient: $P(\text{active}\mid\sigma_\text{TF})$ is a
function of $\sigma_\text{TF}$ alone, so it cannot depend on how many *other* sites merely exist.

Two further facts sharpen this into strong evidence, both established from the construct design
(all opJS4 variants are exactly 630 bp; operator coordinates are shared and absolute --- 4x is the
*first four* positions of 8x; the promoter sits at low coordinates so the proximal operator is at
282 in every construct; spacer blocks are CG-free scramble):

- **A cell-to-cell "permissiveness" confound predicts the opposite sign.** If some cells are
  globally more permissive (more binding *and* more activity), then in a 4x construct "4 bound" is
  the top of the occupancy distribution and so is enriched for permissive cells, whereas in 8x
  "4 bound" is typical. That predicts $P(\text{active}\mid 4/4) > P(\text{active}\mid 4/8)$ ---
  the reverse of the observation.
- **A distance-decay coupling also predicts the opposite sign.** In 4x, four bound sites *must* be
  the four promoter-proximal ones; in 8x they are on average more distal. If proximal binding is
  more effective, again $4/4 > 4/8$.

Both of the obvious alternatives work *against* the observed effect, so controlling for them should
make it larger, not smaller. The one equilibrium mechanism with the right sign --- a distal gap
acting as a sink for nucleosomes that would otherwise cover the promoter --- is excluded by the
equal-length design: conditioned on *which* sites are bound, a 4x and an 8x molecule present the
same physical DNA (the sole loophole being any intrinsic difference in nucleosome affinity between
TetO and scramble sequence).

## 2.3 Panels 3i--3k constrain how potency may depend on condition

Changing dox changes occupancy (3i), but all dox levels **collapse onto one occupancy--activity
curve** (3j), and potency is constant across dox while effective TF concentration varies about
two-fold (3k).

So the map (occupancy $\to$ activity) is a function of occupancy alone: it must not acquire a
separate explicit dependence on dox or on $K$. §4 shows the model satisfies this automatically ---
which is why panel 3j, despite being the panel that *defines* potency, carries **no information
about $w$ at all**.

---

# 3. The model

## 3.1 The drive, and why the mixture form is derived rather than assumed

Let $x(t)$ be the (stochastic) TF occupancy on one molecule over time: stationary, mean
$\langle n\rangle$, autocorrelation $e^{-s/\tau_c}$. Suppose the promoter responds not to the
instantaneous occupancy but to an exponentially filtered version of it with integration time
$\tau_\text{int}$:
$$
D(t) \;=\; \frac{1}{\tau_\text{int}}\int_0^\infty e^{-s/\tau_\text{int}}\,x(t-s)\,ds .
$$
SMF gives us a *snapshot*: we observe $x(t) = n$, not $D$. The relevant quantity is therefore
$D$'s conditional expectation given what we measured:
$$
\mathbb{E}[D \mid x(t) = n]
= \langle n\rangle + (n - \langle n\rangle)\,
  \frac{1}{\tau_\text{int}}\int_0^\infty e^{-s/\tau_\text{int}}e^{-s/\tau_c}\,ds
= \frac{\tau_c}{\tau_\text{int}+\tau_c}\,n \;+\; \frac{\tau_\text{int}}{\tau_\text{int}+\tau_c}\,\langle n\rangle .
$$
Defining
$$
\boxed{\;w \;=\; \frac{\tau_\text{int}}{\tau_\text{int}+\tau_c}\;,\qquad\text{equivalently}\qquad
\frac{w}{1-w} = \frac{\tau_\text{int}}{\tau_c}\;}
$$
gives exactly the mixture
$$
\boxed{\;I \;=\; (1-w)\,n \;+\; w\,\langle n\rangle\;}
$$
So the two-timescale form is **not an ad hoc interpolation**: it is the snapshot-conditional
expectation of a time-integrated drive, and $w$ is a direct readout of the ratio of the promoter's
integration time to the TF's residence time. $w \to 0$ means the promoter reads instantaneous
occupancy; $w \to 1$ means it reads only the long-run average, which is a property of the
*condition*, not of the molecule.

(Approximation used: $\mathbb{E}[f(D)\mid n] \approx f(\mathbb{E}[D\mid n])$, i.e. linear response
of the promoter to the drive. This is exactly the rate-additive assumption of §2.1, so the two
approximations are the same one.)

## 3.2 The full model

Combining §2.1's link with §3.1's drive:
$$
\boxed{\;P(\text{active}\mid n, \langle n\rangle) = p_0 + k_\text{pot}\,I
= p_0 + k_\text{pot}\big[(1-w)\,n + w\,\langle n\rangle\big]\;}
$$
(If a future AD turns out to saturate (§2.1), swap the link for $P = A/(1+A)$ with
$A = p_0' + k_\text{pot}I$. Everything structural below is unchanged; only §4's cancellation
weakens from exact to first-order and §7's estimator changes from least squares to a bounded
nonlinear fit.)

Four parameters, all $\ge 0$:

| parameter | meaning | identified by |
|---|---|---|
| $p_0$ | basal activity at zero occupancy | intercept, Fig. 3d at $\langle n\rangle=0$ |
| $k_\text{pot}$ | **potency** --- total drive added per unit occupancy | the bulk slope (Fig. 3d) |
| $w$ | integrated fraction of potency | *within*- vs *between*-condition contrast (Fig. 3h) |
| $\langle n\rangle$ | not fitted --- supplied per condition | measured (Fig. 3c/3i) or predicted by the array model |

Equivalently, in the linearly-parametrized coordinates
$k_\text{fast} = k_\text{pot}(1-w)$ and $k_\text{slow} = k_\text{pot}w$:
$$
P = p_0 + k_\text{fast}\,n + k_\text{slow}\,\langle n\rangle ,
\qquad k_\text{pot} = k_\text{fast}+k_\text{slow},
\qquad w = \frac{k_\text{slow}}{k_\text{fast}+k_\text{slow}} .
$$
Fit in $(p_0, k_\text{fast}, k_\text{slow})$ --- linear in the parameters, so the fit is
well-behaved and $k_\text{pot}$, $w$ are read off afterwards. (Fitting $w$ and $k_\text{pot}$
directly introduces a product of unknowns for no benefit.)

**Special cases.**

- $w = 0$: the pure equilibrium/instantaneous model. Promoter state depends only on this molecule's
  own current binding. This is what the earlier $J_\text{prom}$ proposal assumed, and Fig. 3h
  rejects it.
- $w = 1$: pure history. $A$ depends only on the condition, so **all** molecules in a condition have
  the same activation probability regardless of their own occupancy --- potency is fully real and
  fully invisible to any same-molecule analysis. This is the limit the objection correctly
  identified.
- $k_\text{slow} = 0$ *and* logistic link: recovers the original $J_\text{prom}$ proposal, now
  doubly falsified.

---

# 4. Why $w$ is invisible in bulk: the exact cancellation

This section formalizes the central intuition: *if the promoter integrates over a long time,
potency is real but does not appear in any same-molecule coupling.*

Average over all molecules in one condition. Since $\mathbb{E}[n] = \langle n\rangle$ by definition,
and $P$ is **linear** in $n$:
$$
\bar P = \mathbb{E}[P]
= p_0 + k_\text{pot}\big[(1-w)\,\mathbb{E}[n] + w\,\langle n\rangle\big]
= p_0 + k_\text{pot}\big[(1-w)\langle n\rangle + w\langle n\rangle\big]
= p_0 + k_\text{pot}\,\langle n\rangle .
$$
**The $w$ terms cancel identically, and with the identity link this is exact --- no approximation.**
The bulk occupancy--activity curve depends only on $p_0$ and $k_\text{pot}$, for *any* value of $w$.

The intuition in one sentence: a molecule's own $n$ is a noisy sample of its condition's
$\langle n\rangle$, and averaging that sample over molecules returns $\langle n\rangle$ again --- so
bulk data cannot tell whether the promoter read the sample or read the mean. Panel h can, because
conditioning on $n$ deliberately examines molecules that *differ* from their condition mean, which
is exactly the information averaging destroys.

An equivalent way to see it, which also gives the shape of the answer:
$$
\underbrace{k_\text{pot}}_{\text{Fig. 3d slope}}
= \underbrace{k_\text{pot}\,w}_{\substack{\text{panel h: slope of each}\\\text{curve in }\langle n\rangle\text{ at fixed }n}}
+ \underbrace{k_\text{pot}(1-w)}_{\substack{\text{panel h: spacing between}\\\text{consecutive }n\text{ curves}}}
$$
**Potency is the conserved quantity; $w$ only decides how it is split between the two axes of panel
h.** Fig. 3d measures the sum and is therefore blind to $w$; panel h resolves the decomposition.

Three consequences:

1. **The published potency is $k_\text{pot}$**, cleanly, with no $w$ correction. The Fig. 3d slope
   and `compute_potency.py`'s $k_\text{tf}$ estimate $k_\text{pot}$, and they remain valid
   regardless of how the fast/slow split turns out. Nothing in the published analysis needs
   revision.
2. **Panel 3j cannot measure $w$.** Its collapse across dox is a genuine and important check that
   potency is a property of the AD/promoter rather than of how occupancy was achieved (§2.3) --- and
   the model satisfies it automatically --- but it is silent on $w$.
3. **Any purely same-molecule coupling estimates $k_\text{fast} = k_\text{pot}(1-w)$, not
   potency.** It *underestimates* potency by exactly the factor $(1-w)$, and reports zero when
   $w = 1$. Hence: a same-molecule coupling is not potency, and potency is not a same-molecule
   coupling. They are two different, both-meaningful projections of the same $k_\text{pot}$.

**If the link is ever found to be nonlinear** (§2.1, e.g. a strong new AD that saturates), the
cancellation degrades from exact to first-order: $\mathbb{E}[f(I)] \ne f(\mathbb{E}[I])$ by Jensen,
and larger $w$ means less molecule-to-molecule variance in $I$ and hence a smaller Jensen
correction. So $w$ would leave a faint bulk signature --- in the *curvature* of Fig. 3d, never its
slope. Far too weak to estimate $w$ from, but it means a joint bulk + single-molecule fit should
then average over the molecule-level distribution of $n$ rather than substituting
$\langle n\rangle$.

## 4.1 Numerical demonstration

`260902_w_regimes_cartoon.py` simulates molecules under $w = 0$, $0.25$, $1$ (20,000 molecules per
copy number, occupancy drawn Binomial$(K, \langle n\rangle_K/K)$ so the condition mean matches by
construction) and confirms both halves of the claim:

| $w$ true | bulk intercept | bulk slope | $\hat w$ recovered | $\hat k_\text{pot}$ |
|---|---|---|---|---|
| 0.00 | 0.031 | 0.090 | 0.004 | 0.090 |
| 0.25 | 0.030 | 0.090 | 0.254 | 0.091 |
| 1.00 | 0.030 | 0.091 | 1.006 | 0.091 |

The bulk line is indistinguishable across regimes (true $k_\text{pot} = 0.09$), while $w$ is
recovered accurately from the molecule-level regression. The script also emits a figure contrasting
the three regimes in bulk (identical) and in panel-h coordinates (qualitatively different).

---

# 5. Where $w$ *is* visible: three tests

All three condition on $n$ and vary $\langle n\rangle$, which is the only way to break the §4
cancellation.

## 5.1 Copy number at fixed $n$ (Fig. 3h --- already measured)

Rearranged, the model predicts that at fixed $n$,
$$
P(n, K) = \underbrace{\big[p_0 + k_\text{pot}(1-w)\,n\big]}_{\text{intercept in }\langle n\rangle}
\;+\; k_\text{pot}w\cdot\langle n\rangle_K ,
$$
i.e. **plotting $P$ against $\langle n\rangle_K$ at fixed $n$ should give a family of parallel
straight lines**, one per $n$, with common slope $k_\text{pot}w$ and spacing $k_\text{pot}(1-w)$
between consecutive $n$. That single plot yields both parameters.

> **Correction (2026-09-04): the "sum rule" is NOT a test.** An earlier version of this note and
> of §6 treated $k_\text{pot}w + k_\text{pot}(1-w) = k_\text{pot}^\text{bulk}$ as a falsifiable
> check. When both sides are obtained from **one joint fit to the same data it is an algebraic
> identity**, verified to hold to machine precision on random data. The reason is the §4
> cancellation reappearing in the estimator: within an amplicon the $N$-weighted mean of $n$ *is*
> $\langle n\rangle$, so writing $n = \langle n\rangle + d$ makes $d$ within-amplicon centred and
> therefore $N$-orthogonal to every amplicon-level function. The design splits into orthogonal
> between/within blocks, $k_\text{pot}$ is estimated purely from between-amplicon variation (i.e.
> it *is* the bulk fit) and $k_\text{fast}$ purely from within. They cannot disagree.
> (§6's version compares two *independent hand readings* off two published panels, which is
> informative about the figures' mutual consistency — but it is not a test of the model, and it is
> not what any joint fit computes.)
>
> **The real structural test is PARALLELISM**: add an interaction $n\times\langle n\rangle$ and
> check its coefficient is zero. Non-zero means the lines fan out and the two-timescale form is
> wrong. `fit_two_timescale_potency.py` reports `k_interact` and its $z$-score for exactly this. Panel 3h contains this information already --- it plots $P$ against $K$, and
the model wants $\langle n\rangle_K$ on $x$, which is just re-plotting through Fig. 3c. Panel 3h
lacks a statistical test; this gives it one.

The three regimes are qualitatively distinguishable by eye in this plot (simulated in
`260902_w_regimes_cartoon.py`):

| regime | panel h appearance | interpretation |
|---|---|---|
| $w = 0$ | curves **flat** in $\langle n\rangle$, maximally separated by $n$ | promoter reads instantaneous occupancy only |
| $0 < w < 1$ | **parallel sloping** lines, partially separated | mixture; slope/spacing ratio gives $w/(1-w) = \tau_\text{int}/\tau_c$ |
| $w = 1$ | curves **collapse** onto a single sloping line | promoter reads history only; potency fully real, single-molecule signal zero |

Note that $w=1$ is the regime in which potency exists and is measurable in bulk while *no*
same-molecule analysis can detect any coupling whatsoever --- the case that motivated this note.

## 5.2 Dox at fixed $K$ and fixed $n$ --- the clean test, not yet run

This is the sharpest available experiment and the data already exist. Hold copy number fixed, so
the construct is *physically identical*; vary dox, which changes $\langle n\rangle$ (Fig. 3i) with
no change whatsoever in geometry, site positions, or sequence. Then compare molecules with the
*same* $n$ across dox levels:

$$
\text{if } w > 0:\quad P(\text{active}\mid n,\ \text{high dox}) \;>\; P(\text{active}\mid n,\ \text{low dox})
$$
$$
\text{if } w = 0:\quad \text{these are equal.}
$$

Why this beats §5.1: the copy-number comparison must contend with *which* sites are bound differing
between constructs (§2.2), even though that confound has the wrong sign. The dox comparison has no
structural confound at all --- same molecule design, same sites, only the occupancy statistics
differ. A monotone increase in $P(\text{active}\mid n)$ with dox at fixed $K$ and fixed $n$ is
essentially unimpeachable evidence for $w > 0$.

## 5.3 Fixed exact configuration (v5 per-site calls)

Strengthen §5.1 by conditioning on *which* sites are bound rather than only how many: compare
"operators 282/322/362/402 bound" in an 8x construct against the same four sites bound in a 4x
construct. Given the equal-length shared-coordinate design (§2.2), those molecules are physically
near-identical, so a residual difference cannot be equilibrium physics of a single snapshot. v5
emits per-site occupancy (`site_{name}` columns), so this is a re-tabulation, not new machinery.
The old count-based analysis could not run this test.

---

# 6. Estimating $w$ from the published figure (illustrative only)

Approximate values read by eye off Fig. 3h/3c, purely to show the arithmetic and check that the
model is identifiable and gives a plausible answer. **These are not measurements** --- they should
be recomputed from the molecule tables.

Take $\langle n\rangle \approx 0.15$ at $K=2$ and $\approx 4.1$ at $K=9$, so
$\Delta\langle n\rangle \approx 4.0$.

*Between-condition contrast (the slope), at fixed $n = 0$:* $P \approx 0.01 \to 0.10$ from $K=2$ to
$K=9$, so $\Delta P \approx 0.09$ and
$$k_\text{pot}\,w \;\approx\; 0.09/4.0 \;\approx\; 0.023 .$$

*Within-condition contrast (the spacing), at fixed $K = 9$:* $P \approx 0.10 \to 0.47$ from $n=0$ to
$n=5$, so $\Delta P \approx 0.37$ and
$$k_\text{pot}(1-w) \;\approx\; 0.37/5 \;\approx\; 0.074 .$$

$$
k_\text{pot} \approx 0.097, \qquad w \approx \frac{0.023}{0.097} \approx 0.23,
\qquad \frac{\tau_\text{int}}{\tau_c} = \frac{w}{1-w} \approx 0.31 .
$$

So roughly a quarter of potency would be integrated and three quarters instantaneous, with the
promoter integrating over about a third of a TF residence time.

**A genuine consistency check falls out of this.** $k_\text{pot}$ estimated from panel h
($\approx 0.097$) is an *independent* route to the same quantity as the Fig. 3d slope
($\approx (0.47-0.03)/4.9 \approx 0.090$). Those agree to about 7%, which is non-trivial support for
the model: nothing forced two different panels, analysed two different ways, to give the same
potency. Confirming this properly on the real molecule tables is worth doing early --- if the two
disagree badly, the decomposition is wrong and no amount of fitting will fix it.

**$w$ is link-dependent, so §2.1 comes first.** The identical arithmetic gives $w \approx 0.14$
under the hyperbolic link and $w \approx 0.59$ under the logistic --- a qualitatively different
conclusion (majority-integrated, $\tau_\text{int} \sim \tau_c$) from the same readings. The linear
link is the empirically favoured one, so $w \approx 0.23$ is the current best guess, but the
interpretation of $w$ is not robust to the link, and the link is settled on bulk data.

---

# 7. Fitting

**Data.** One row per molecule: $n$ (v5 `n_tf_teto`), the promoter activity call (v5
region-anchored `{region}_open` / `{region}_nuc`), and the molecule's condition key
$(K, \text{dox}, \text{AD}, \text{background})$. Join $\langle n\rangle$ on by condition.

**Objective.** With the identity link the model is a **linear probability model**,
$$
P = p_0 + k_\text{fast}\,n + k_\text{slow}\,\langle n\rangle ,
$$
so it fits by ordinary (or molecule-count-weighted) least squares --- a closed form, no optimizer,
no convergence question, and the same WLS family `compute_potency.py` already uses. This is a
side benefit of §2.1's empirical finding: the link that fits best is also the one that makes the
estimator trivial.

Two standard caveats of linear probability models apply and are both benign here: predictions are
not guaranteed to lie in $[0,1]$ (harmless while $\bar P \lesssim 0.5$, but check, and revisit if a
strong new AD pushes activity high), and the residuals are heteroscedastic, so use robust or
$\sqrt{P(1-P)}$-weighted standard errors rather than the naive OLS ones. If either becomes a
problem, switch to the hyperbolic link and a bounded nonlinear fit --- note that link's likelihood
is *not* concave ($\partial^2/\partial A^2$ of $y\log A - \log(1+A)$ is $-y/A^2 + 1/(1+A)^2$,
negative for $y=1$ but positive for $y=0$), so it needs multi-start.

**Identifiability.** $k_\text{fast}$ comes from *within*-condition variation in $n$;
$k_\text{slow}$ from *between*-condition variation in $\langle n\rangle$. Both are needed. A single
condition identifies only $k_\text{fast}$; bulk data alone identifies only $k_\text{pot}$. The
copy-number series and the dox series each supply the between-condition axis, and having both
allows a consistency check: **$w$ estimated from the copy-number series and from the dox series
should agree.** If they do not, the "history" variable is not simply time-averaged occupancy.

**Hierarchy.** Same discipline as elsewhere: $p_0$ per promoter (sequence differences are real);
$k_\text{pot}$ per AD (that is the scientific quantity); $w$ shared across ADs to start --- it is a
timescale ratio, so per-AD variation is a claim needing evidence. Freeing $w$ per AD is the
interesting extension (an AD acting through a stably recruited cofactor should have larger $w$ than
one acting by transient PIC recruitment), but it should be earned by model comparison, not assumed.

**Do not use `--background None`.** The runs of record pooled opJS4 (b0) and opJS5 (b1/b2) into one
regression per sample. Fig. 3b shows the backgrounds differ substantially in activity at matched
copy number, so $p_0$ at least must be indexed by background here.

---

# 8. Mechanism: what is doing the integrating?

If $w > 0$, some slow variable carries memory of past binding. The leading candidate is
**chromatin**, and it is testable with existing data.

The old `3param_nuc` model already encodes the relevant coupling: it shifts the nucleosome energy
whenever any TF is bound,
```
num_tfs * tf_e + num_nuc * (nuc_e - delta_e * (num_tfs > 0))
```
i.e. a recruited-remodeler term. If a remodeler evicts a nucleosome and the depleted state persists
longer than TF residence, then the nucleosome *is* the integrator, and $\tau_\text{int}$ is a
nucleosome re-formation time.

**Test:** run the §5.1/§5.2 analysis on **nucleosome occupancy** rather than promoter activity.

- If $P(\text{nuc}\mid n)$ shows the same $\langle n\rangle$-dependence, the memory is in the
  chromatin layer --- and $\delta_\text{remodel}$, fitted in the equilibrium array model, becomes
  the mechanistic parameter behind $w$.
- If nucleosome occupancy shows *no* $\langle n\rangle$-dependence while promoter activity does, the
  integrator is downstream of chromatin (cofactor residence, PIC assembly, a modification).

Either outcome identifies the physical carrier of the memory, which is the mechanistic basis for
linear potency that the original analysis could not supply.

---

# 9. What this model costs

Stated plainly, because it is a real conceptual price:

- **It is not a thermodynamic model.** $\langle n\rangle$ is a property of the *ensemble*, not of
  the configuration, so no Boltzmann distribution over single-molecule states can produce this
  conditional. $k_\text{pot}$ and $w$ are not free energies and should not be reported in $k_BT$.
  The reciprocity argument (TF stabilizes the promoter and the promoter stabilizes TF binding via
  one symmetric coupling) applies only to the $w=0$ limit, which the data reject.
- **The equilibrium array model is unaffected.** TetO occupancy, cooperativity $J$, nucleosome
  fugacity $\mu$, and $\delta_\text{remodel}$ all remain a legitimate equilibrium lattice-gas
  problem (`260717_ising_intuition.md` §§1--6, §8). Only the *promoter* node leaves equilibrium.
  The clean architecture is therefore two stages: an equilibrium array model that predicts
  $\langle n\rangle$ and the per-molecule distribution of $n$, and a kinetic promoter model on top
  of it. §7 of the intuition doc --- the promoter as one more Ising node --- should be treated as
  superseded by this note.
- **Snapshots cannot see timescales directly.** $\tau_\text{int}/\tau_c$ is inferred from the
  *ratio* $w/(1-w)$ under the linear-filter assumption of §3.1, not measured. Converting it to
  seconds requires an independent estimate of rTetR residence time.

---

# 10. Checklist

1. **Confirm the link per AD, not pooled** (§2.1). $\bar P$ vs $\langle n\rangle$ should be
   straight; the interesting case is a strong new AD that reaches high $\bar P$ and either stays
   straight (a real result) or bends (the bend identifies the link). Pooling ADs averages this away.
2. **Re-plot 3h as the model wants it** (§5.1): $P$ vs $\langle n\rangle_K$ at fixed $n$. Check the
   lines are parallel. Read off $k_\text{pot}w$ (common slope) and $k_\text{pot}(1-w)$ (spacing),
   and check that $k_\text{pot}$ agrees with the Fig. 3d slope (§6).
3. **Run the dox test** (§5.2): $P(\text{active}\mid n)$ vs dox at fixed $K$. This is the
   confound-free test for $w > 0$.
4. **Cross-check** that $w$ from the copy-number and dox series agree (§7).
5. **Condition on exact configuration** (§5.3) using v5 per-site calls.
6. **Ask whether the memory is chromatin** (§8) by repeating 2--3 on nucleosome occupancy.
7. Only then fit the full model, hierarchically, with backgrounds separated (§7).

Steps 1--3 need no new machinery --- they are re-tabulations of existing molecule tables --- and
their outcomes determine whether steps 4--7 are worth building.
