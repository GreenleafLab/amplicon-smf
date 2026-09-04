---
title: "The Ising / MaxEnt Reframe of the Partition-Function Model --- Intuition"
author: "Design discussion notes"
date: "2026-07-17"
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
  - \usepackage{booktabs}
  - \usepackage{microtype}
---

# Scope

Intuition behind recasting the thermodynamic (partition-function) model as an Ising / maximum-entropy
model to consume the `v5` HSMM binding calls. Companion to `260715_v5_downstream_design.md` (which
has the formal derivation); this file is the conceptual walkthrough: what an Ising model is, why it
is the natural model for this system, its exact relationship to the old
`fit_partition_function_model_v3.py`, whether the old approach was wrong, and how to think about the
one genuine modeling decision (discretizing nucleosomes).

---

# 1. What an Ising model is

The Ising model is the canonical physics model of a set of interacting **binary** variables. Each
site $i$ has a variable $\sigma_i \in \{0,1\}$ with

- a **field** $h_i$: the energy of that site being "on" by itself, and
- a **coupling** $J_{ij}$: the extra energy of a *pair* of sites being on together.

The energy of a whole configuration $\sigma$ is **linear** in these pieces,

$$
E(\sigma) = \sum_i h_i\,\sigma_i \;+\; \sum_{ij} J_{ij}\,\sigma_i\sigma_j ,
$$

and configurations are Boltzmann-distributed,

$$
P(\sigma) = \frac{1}{Z}\,e^{-E(\sigma)}, \qquad Z = \sum_\sigma e^{-E(\sigma)} .
$$

The statistics view is the more useful one here: this is exactly a **maximum-entropy (MaxEnt)**
model. If the only things you are willing to commit to about your system are a set of *average*
quantities --- mean per-site occupancy, mean neighbor co-occupancy, mean number of nucleosomes ---
then the least-committal distribution that reproduces those averages (maximum entropy subject to
matching them) is provably

$$
P(\sigma) \propto \exp\!\Big(-\sum_k \theta_k\,\phi_k(\sigma)\Big),
$$

i.e. Boltzmann with energy linear in those feature functions $\phi_k$. So **"Ising model,"
"pairwise MaxEnt model," and "exponential-family Gibbs distribution" are three names for the same
object.** Fields $h_i$ are single-site energies; couplings $J_{ij}$ are cooperativities.

---

# 2. Why it is the natural model here

The system is almost a textbook Ising setup:

- **Binary variables on a 1D template.** Each TetO is bound or not; each region is
  nucleosome-occupied or not.
- **Fields.** The intrinsic affinity of rTetR for an operator (set by dox, the AD, etc.) is a
  per-site field $h_i$.
- **Couplings.** rTetR binding may be **cooperative** (one bound operator helps its neighbor)
  $\Rightarrow$ a $J$ between adjacent sites. Nucleosomes **compete** with TF binding for the same
  DNA $\Rightarrow$ a TF--nucleosome antagonism coupling.
- **The data is exactly right for it.** SMF gives *single-molecule joint configurations*, not just
  bulk marginals. Bulk methylation would give only the marginal occupancies $\langle\sigma_i\rangle$,
  which cannot separate "high affinity" from "cooperativity." Single molecules give the
  co-occupancies $\langle\sigma_i\sigma_j\rangle$ directly, and those are precisely what identify the
  coupling $J$. An Ising model is the tool that turns *which sites are co-bound on the same molecule*
  into an energy decomposition.

So Ising is not a fancy alternative --- the biology **is** a 1D interacting-occupancy system, and
SMF measures the joint states such a model needs.

---

# 3. Relationship to the old model, and: was the old approach wrong?

**No. And the thing that surprised you --- that the fit was "just counting total TFs and nucs" ---
is the tell that you had it right.**

## 3.1 The old model already *is* an Ising / MaxEnt model

`fit_partition_function_model_v3.py` computes, per configuration,

$$
E(\sigma) = n_{\text{tf}}\,e_{\text{tf}} + n_{\text{nuc}}\,e_{\text{nuc}} \;(+\ \text{coupling terms}),
$$

energy linear in feature counts (`assign_energy_*`), then makes it Boltzmann
(`convert_to_probabilities`). That is a MaxEnt/Ising model already. The *only* thing that was
"enumeration-based" is **how it computed $Z$**: it summed $e^{-E}$ over the rows of
`valid_states_table.txt` (brute force over microstates) instead of using a transfer matrix. That is
a computational choice about $Z$, **not** a different model.

## 3.2 Why the likelihood collapsed to counts (the surprising part)

The negative log-likelihood in that script is

$$
\text{NLL}(\theta) = \sum_m E(\sigma_m) + N\log Z(\theta)
= N\sum_k \theta_k \langle\phi_k\rangle_{\text{data}} + N\log Z(\theta),
$$

where

$$
\langle\phi_k\rangle_{\text{data}} = \frac{1}{N}\sum_m \phi_k(\sigma_m)
$$

is just the **average of feature $k$ over your molecules** --- mean number of TFs bound, mean number
of nucleosomes, mean neighbor co-occupancy. The individual molecules enter the likelihood *only*
through those averages: two datasets with the same mean $n_{\text{tf}}$ and mean $n_{\text{nuc}}$
give an identical likelihood and an identical fit.

That is **not a bug** --- it is the defining property of an exponential family: the feature averages
are **sufficient statistics** (Pitman--Koopman--Darmois). It is the same reason the MLE of a Gaussian
needs only the sample mean and variance, not the individual points. Noticing "this is just counting
totals" means you correctly spotted that the likelihood factors through sufficient statistics ---
which is exactly why the MLE is **moment matching**,

$$
\boxed{\ \langle\phi_k\rangle_{\text{model}}(\theta) = \langle\phi_k\rangle_{\text{data}}\ }
$$

and why the reparametrization can throw away the per-molecule `idx` join and just read the moments
off v5's calls.

## 3.3 What the reframe actually changes (almost nothing conceptual)

The model and the objective are the *same*. The reframe buys purely mechanical wins:

1. **Compute $Z$ by transfer matrix / lattice-gas DP** instead of materializing every microstate ---
   scales to any copy number, no enumeration explosion, and it is what lets you consume v5 (which has
   no enumerated state table).
2. **Convexity.** $\log Z$ is log-sum-exp $\Rightarrow$ convex $\Rightarrow$ NLL is convex
   $\Rightarrow$ unique optimum; gradient ascent instead of Nelder--Mead groping.
3. **Free error bars.** The Hessian is $N\,\mathrm{Cov}_{\text{model}}(\phi)$, giving parameter
   confidence intervals directly --- you currently have none.

## 3.4 Should you worry about past results?

Two honest caveats, both about *modeling choices baked into the old enumeration*, not the fitting
logic:

- **The TF energetics are trustworthy.** TF occupancy is clean binary-per-TetO; the sufficient
  statistic (mean bound count / co-occupancy) is unambiguous, and $Z$ over TF configs is exactly
  right whether enumerated or transfer-matrix'd.
- **The nucleosome side is the soft spot.** The old code counts a nucleosome only if its midpoint
  lands in a hardcoded `fixed_e_bounds = (200,600)` window --- i.e. it *already* silently made the
  "coarse region indicator" discretization choice, with an arbitrary boundary. And brute enumeration
  over a discrete nuc-state list means the *positional degeneracy* of nucleosomes (how many ways a
  nuc can sit) is however the state table happened to encode it. If that encoding under/over-counted
  placements, $Z$ and therefore the fitted **nucleosome** energy (and the TF--nuc antagonism $\delta$)
  are biased. The TF field and cooperativity are robust to this; the nucleosome energy is the number
  to treat with suspicion.

## 3.5 A clarification on "mean field"

This is **not** a mean-field model, and that is good news. Mean-field is an *approximation* that
replaces couplings with an averaged field and pretends sites are independent. The reason the 2-param
fit reduced to total counts is that with no coupling ($J=0$) the sites genuinely *are* independent
--- total count is exactly sufficient, no approximation involved. And once you *do* add coupling, the
1D transfer matrix computes $Z$ **exactly** --- strictly better than mean-field, which you would only
reach for in 2D/3D where exact $Z$ is intractable. You never need the mean-field approximation here.

## 3.6 What the fit is *for*

If you only want "how much binding," the empirical moments (mean occupancy per site) already answer
that --- no partition function needed. The Ising fit earns its keep when you want **mechanism**:
separating intrinsic affinity ($h$) from cooperativity ($J$) from nucleosome antagonism, reported in
$k_BT$ units comparable across ADs and doses, plus a generative model you can extrapolate (e.g.
predict occupancy at an unmeasured dox dose). That is the scientific payload.

---

# 4. The nucleosome-binning decision

This is the one real modeling choice, and it is real precisely because --- unlike TFs ---
nucleosomes have no natural discrete coordinate. v5 emits them as continuous-position segments; the
partition function needs a discrete state variable *and* a rule for how $Z$ counts nucleosome
configurations.

## 4.1 Coarse: region indicators

A handful of binary "slots" --- promoter-nuc yes/no, one or two array-region slots --- each a feature
coupled to TF / promoter state. This is essentially what the old code did (nuc-midpoint-in-window).

- **Pros:** few parameters, trivial $Z$, directly comparable to the old fits, and the boundary can be
  made principled instead of `(200,600)`.
- **Cons:** throws away position; a 200--400 bp "slot" can hold one or two nucs and the indicator
  cannot tell; slot boundaries are somewhat arbitrary.

## 4.2 Fine: 1D hard-rod lattice gas

Bin the template; a nucleosome occupies a ${\sim}147$ bp run of bins with hard-core exclusion; one
fugacity / energy per nucleosome. $Z$ is a forward DP over positions.

- **Pros:** the physically correct thermodynamic nucleosome model; the DP **integrates over
  positional degeneracy automatically** --- the "many places a nuc could sit" multiplicity that the
  old enumeration handled ad hoc is now done correctly and cheaply, $O(\text{length})$.
- **Cons:** more machinery, a bin-size choice, and it must be reconciled with v5's notion of a
  nucleosome.

## 4.3 Recommendation: the hard-rod lattice gas is the model; region indicators are a sanity check

The primary formulation is the **1D hard-rod lattice gas** (§4.2), for the reasons developed in
detail in Section 5: it is the physically correct thermodynamic nucleosome model, it integrates over
positional degeneracy automatically, and --- decisively --- it is the *only* one of the two that can
represent the excluded-volume packing that generates the site-dependent occupancy gradient across
identical operators. Region indicators cannot represent a reservoir or phasing at all (§5.5), so they
are not a candidate for the real model; they are a scaffold.

1. **Region indicators are a one-off sanity check, not a stepping stone.** Use them only to
   **reproduce the old fit** and confirm the TF field + cooperativity come out where the old
   enumeration put them, now with real convexity and error bars. Those TF parameters are the paper's
   headline numbers (AD potency lives in the TF energetics) and are clean regardless of how
   nucleosomes are binned --- so this check is cheap and worth doing once.
2. **Build the hard-rod DP as the actual nucleosome model.** Everything downstream that touches
   nucleosome energetics --- the TF--nucleosome antagonism $\delta$, any nucleosome-eviction
   $\Delta\Delta G$ per AD, and the site-dependent occupancy prediction of Section 5 --- should come
   from the lattice gas, not the indicators. Do not treat the coarse model as the default and the rod
   as an optional upgrade; the rod is the default and the coarse model is a disposable cross-check.

## 4.4 Two consistency points, either way

- **Match v5's notion of a nucleosome.** The moment you feed the fit
  ($\langle n_{\text{nuc}}\rangle$, per-region occupancy) is already processed with v5's duration
  prior (mode 147, min 100 / max 210) baked in. If the lattice gas assumes a different footprint
  length, that is a silent mismatch. Cleanest framing: the partition function models the distribution
  of **v5's discrete calls**, treating those calls as the observed configuration. Marginalizing over
  v5's posterior instead of hard calls is the rigorous upgrade if observation noise ever bites, but
  hard calls are the right start.
- **Keep the operator field shared; do NOT fit free per-site $h_i$.** The TetOs are
  sequence-identical, so their intrinsic affinity should be one shared $h$. The observed
  position-dependent occupancy (edges vs. middle) is an emergent excluded-volume / packing effect,
  not intrinsic heterogeneity --- fitting it with free $h_i$ would launder a collective effect into
  fake affinity differences. See Section 5. (This refines an earlier draft suggestion to use per-site
  fields --- that was wrong for identical operators.)

---

# 5. Site-dependent occupancy across identical operators

## 5.1 The puzzle

The TetO sites are sequence-identical, yet occupancy varies across the array (edge sites lower,
interior sites higher). This should **not** be read as differing intrinsic affinity. Encoding it as
site-specific fields $h_i$ would be a mechanistic error: it fits a collective effect as if it were a
property of individual sites.

## 5.2 It is statistical positioning (Kornberg--Stryer)

The right picture is a **1D lattice gas of mixed-size hard rods** competing for the template. Large
objects (nucleosomes, ${\sim}147$ bp) and small objects (TFs, one operator) compete for the same DNA
with excluded volume. Given some interior operators are TF-bound, a nucleosome can only fit in the
flanking gaps, and it is geometrically hard to occlude an interior site without conflicting with its
bound neighbors. The result --- nucleosomes toward the flanks, TFs in the interior --- is the
equilibrium of such a system. It is emergent, collective, and entropic: **uniform affinities suffice.**

## 5.3 Uniform affinity recovers the gradient

With a shared TF field $h$, a single nucleosome fugacity, and hard-rod excluded volume, the per-site
TF occupancy gradient (middle $>$ edge) **emerges** from the partition function. Mechanism:
conditional on interior TFs bound, a nucleosome covering an interior site must conflict with
neighboring bound TFs and with rod-packing, whereas a nucleosome covering an edge site can hang into
the free exterior with no conflict --- so edge sites face more nucleosome competition and lose TF
occupancy. Equivalently: for a given nucleosome load there are more (lower-free-energy) ways to seat
nucleosomes toward the flanks/exterior than threaded through the TF-dense interior.

The statistical payoff: a $K$-number occupancy profile is **predicted** by ${\sim}2$ parameters
($h$, fugacity), rather than fit by $K$ free $h_i$ (saturated, non-mechanistic). If the prediction
matches, that is a real result --- "site-dependent occupancy is fully accounted for by nucleosome
packing against identical sites, zero intrinsic heterogeneity."

## 5.4 The infinity / boundary subtlety --- this is where the *imposition* matters

The nucleosome array effectively extends beyond the amplicon on both sides, but only a small window
is measured, so the boundary conditions on $Z$ must be chosen with care.

**Key fact:** an infinite, uniform, barrier-free nucleosome array has *flat* occupancy
(translational symmetry). Positional structure requires a **barrier**. In this system the barriers
are the **bound TFs** (and the promoter / any real boundary element) --- **not** the amplicon edge,
which is a measurement window, not a physical wall.

Consequence: **do not** compute $Z$ on a finite lattice with hard walls at the amplicon boundary ---
that injects a *spurious* barrier and manufactures fake edge phasing/depletion. Correct treatment:

- **Grand-canonical**, not fixed count: a nucleosome fugacity $\mu$ (chemical potential); the count
  in the window fluctuates, so nucleosomes enter/leave from outside.
- **Pad the lattice by $\ge$ one footprint each side** and allow dyads *outside* the window
  (nucleosomes overhanging inward), so edge accessibility is not artificially starved. This is
  exactly the "ghost-segment" room v5's decoder already reserves
  (`amp_width = last GpC + max(nuc_max/2, ...) + 10`).
- **Exterior = bulk nucleosomal array at ambient density**, entering the transfer matrix only as a
  boundary vector.

**Reasoning about what you cannot see:** you do *not* model individual unseen nucleosomes --- you
integrate them out. The entire infinite exterior is summarized by just two things: the fugacity
$\mu$ (set by ambient nucleosome density) and the boundary vector ("what is just outside is bulk
chromatin"). The measured window's occupancy profile is then the exact marginal of that infinite
system, computed by the transfer matrix / DP over the padded window. That is the whole advantage of
the grand-canonical lattice gas: the bath is one number.

## 5.5 Consequences for the model design

- Strong argument for the **fine hard-rod lattice gas** over coarse region indicators (Section 4):
  indicators cannot represent a reservoir or phasing at all.
- The nucleosome fugacity $\mu$ is a fit parameter, constrained by the window's overall nucleosome
  occupancy.
- Keep the operator field $h$ **shared**; the per-site gradient is a prediction to test, not a
  parameter. A residual *beyond* the packing prediction is where a site-specific term becomes
  evidence of something real (a genuine end effect, a promoter-proximal asymmetry) --- earned, not
  assumed. This is the same discipline UNID footprints get: explanation of last resort.

## 5.6 The lattice-gas DP, concretely

Discretize the padded template into bins of size $b$ (start $b=1$ bp for correctness; coarsen later
for speed --- bin size affects resolution/cost, not the model). A nucleosome is a hard rod of length
$L\approx147$ bp that, if it starts at bin $i$, occupies $i\dots i+L-1$ and **excludes** any other
rod or bound TF from those bins. There is **one** nucleosome parameter, a fugacity
$z=e^{-\varepsilon_{\text{nuc}}}$ (equivalently a chemical potential $\mu$) --- the ambient-density
knob. Position enters *only* through geometry (where a rod fits without overlap), so positional
degeneracy is summed automatically.

Let $Z_i$ be the summed Boltzmann weight of all valid configurations up to bin $i$. The forward
recursion has two moves at each bin:

$$
Z_i = \underbrace{Z_{i-1}\,w_{\text{empty}}(i)}_{\text{bin }i\text{ not a rod start}}
\;+\;
\underbrace{Z_{i-L}\,z\,w_{\text{TFcompat}}(i)}_{\text{a nucleosome starts at }i}.
$$

- **First term:** leave bin $i$ un-nucleosomed. It may instead be a TF site --- that is the coupled
  TetO transfer matrix carrying the shared field $h$ and cooperativity $J$, run on the same lattice.
- **Second term:** place a rod over $i\dots i+L-1$, pay $z$, allowed only if the footprint does not
  collide with a bound TF in those bins. That collision term is the **TF--nucleosome antagonism**
  $\delta$ (hard exclusion, or a soft penalty if TFs and nucleosomes are only partly mutually
  exclusive).

This is $O(\text{length})$, $\log Z$ is convex, and $Z$ is exact --- no enumeration, no state table.
Faithful footprint length: let the "start a rod" term sum over $L\in[100,210]$ to match v5's duration
prior instead of a rigid 147 (a few extra terms; do this only if a nucleosome-energetics claim needs
it).

**Boundary conditions (§5.4, restated operationally):** pad the lattice by $\ge L$ each side beyond
the first/last scored GpC and allow rod dyads in the pad (nucleosomes overhanging inward); seed the
recursion's initial/terminal vector with the **bulk stationary occupancy at fugacity $z$**, not an
empty/occupied wall. Correctness check: a barrier-free uniform stretch must return **flat**
occupancy. Any phasing then present is caused only by an explicitly modeled barrier (a bound TF, the
promoter) --- never by the amplicon frame.

**The whole model is four numbers**, all in $k_BT$ with Hessian error bars: $h$ (shared TetO field),
$J$ (TetO--TetO cooperativity), $\varepsilon_{\text{nuc}}$/$z$ (nucleosome fugacity), $\delta$
(TF--nucleosome antagonism). The middle$>$edge TF-occupancy gradient must **emerge** from these four
via excluded volume, not be fit with per-site $h_i$.

---

# 6. One-paragraph summary

Your old partition-function fit was already an Ising / MaxEnt model; its likelihood depended on the
data only through average feature counts because those averages are the model's sufficient statistics
--- a correct hallmark of the exponential family, not an error. The reparametrization keeps the model
and objective identical and only (a) computes $Z$ with a transfer matrix / lattice-gas DP so it needs
no enumerated state table (and thus consumes v5), (b) exploits convexity for a unique optimum, and
(c) delivers Hessian-based error bars. The TF energetics from the old fit are trustworthy; the
nucleosome energy is the one number to re-examine, because the nucleosome discretization (the sole
genuine modeling choice) was made implicitly by a hardcoded midpoint window. The nucleosome model
should be a proper **1D hard-rod lattice gas**; coarse region indicators are worth running once only
as a disposable sanity check that reproduces the old fit and locks the TF field + cooperativity, not
as the working model, because indicators cannot represent nucleosome packing or a reservoir at all.
Because the operators are sequence-identical, keep their field shared: the site-dependent occupancy
across the array is statistical positioning (an excluded-volume packing effect against the bound-TF
barriers), so it should *emerge* from a uniform-affinity rod lattice gas as a prediction, not be fit
with per-site energies --- provided $Z$ treats the amplicon as a grand-canonical window on an
effectively infinite array (fugacity $\mu$ + a bulk boundary), never a hard-walled finite box.
The promoter can then be added as one extra node coupled to TF occupancy (Section 7) --- the point
where "potency" becomes a fitted coupling in $k_BT$.

---

# 7. Incorporating the promoter (optional extension)

The array model above (TetOs + nucleosomes) reproduces the current repo's partition-function scope.
The promoter is the natural next node: it is where binding is *read out* as activation, and encoding
it in the same $Z$ turns the paper's central quantity --- potency --- from a downstream curve-fit
into a single energetic parameter.

## 7.1 Staging

1. **Reproduce first.** Array only (TetOs + nucleosomes); promoter stays the separate
   `assign_promoter_state_from_model_new.py` pipeline. This is the conditional factorization the repo
   already uses and the sanity target for the new code.
2. **One promoter node.** Add a single binary $\sigma_{\text{prom}}\in\{\text{closed, open}\}$ with
   its own field and a coupling to TF occupancy (this section).
3. **Multi-state promoter.** Promote it to a small categorical (Potts) node
   $\{\text{nuc, open, TBP/PIC, paused}\}$, each state with its own energy, fed by v5 promoter
   footprints or the promoter-state caller. Only when a claim needs the substates.

## 7.2 The promoter as one more node

Give the promoter a field $h_{\text{prom}}$ (intrinsic openness) and, crucially, a coupling
$J_{\text{prom}}$ to the bound TFs. Two things make it well-behaved:

- **It is a *real* barrier**, unlike the amplicon edge or the reservoir nucleosomes: it sits at a
  fixed, sequence-defined position, so a position-specific promoter term is legitimate (§5.4 lists
  the promoter alongside bound TFs as an allowed barrier).
- **Per-promoter heterogeneity is *earned*.** Different minimal promoters differ in sequence, so
  fitting a separate $h_{\text{prom}}$ (and $J_{\text{prom}}$, and any substate energies) per
  promoter is real, not a laundered collective effect --- the exact opposite of the per-site $h_i$
  we forbid for the *sequence-identical* operators. The clean structure is **hierarchical**: share
  the array parameters ($h,J,\mu$ --- same TetO array everywhere) and vary the promoter parameters
  per promoter.

## 7.3 The coupling term --- the one modeling choice

Write the coupling as

$$
E_{\text{coupling}}(\sigma) = -\,J_{\text{prom}}\; f(\sigma_{\text{TF}})\;\sigma_{\text{prom}},
$$

so an open promoter is stabilized by $J_{\text{prom}}\,f$ where $f$ is *some function of the TF
configuration*. **Choosing $f$ is a mechanistic hypothesis about how binding talks to the
promoter**, and --- the payoff of single-molecule data --- SMF can *discriminate* between the
choices, because each molecule reports the TF configuration *and* the promoter state jointly:

| $f(\sigma_{\text{TF}})$ | mechanism | identifying joint moment |
|---|---|---|
| $n_{\text{tf}}$ (count / dosage) | each bound TF adds equally to opening | $\langle n_{\text{tf}}\,\sigma_{\text{prom}}\rangle$ |
| $\mathbf{1}[n_{\text{tf}}\ge 1]$ (switch) | any binding flips the promoter | $\langle \mathbf{1}[n_{\text{tf}}\ge1]\,\sigma_{\text{prom}}\rangle$ |
| $\mathbf{1}[n_{\text{tf}}\ge k]$ (threshold) | needs $k$ bound to fire | $P(\text{open}\mid n_{\text{tf}})$ curve shape |
| $\sigma_{K}$ (nearest operator) | only the promoter-proximal site acts | $\langle \sigma_K\,\sigma_{\text{prom}}\rangle$ |
| $\sum_i e^{-d_i/\lambda}\sigma_i$ (distance decay) | communication has a range $\lambda$ | per-site $\langle\sigma_i\,\sigma_{\text{prom}}\rangle$ profile |

The most agnostic version fits the per-operator profile $\langle\sigma_i\,\sigma_{\text{prom}}\rangle$
directly ($K$ features) and *reads off* whether it is flat (dosage), edge-weighted (proximity), or
threshold-like. The distance-decay form is the principled 2-parameter reduction ($J_{\text{prom}}$,
range $\lambda$) if the full profile is too noisy.

## 7.4 Why this is exact and cheap: the field-shift identity

Because $\sigma_{\text{prom}}$ is a single binary variable, sum over its two states explicitly:

$$
Z \;=\; \underbrace{Z_{\text{array}}(h,J,\mu)}_{\sigma_{\text{prom}}=0}
\;+\;
e^{\,h_{\text{prom}}}\;\underbrace{Z_{\text{array}}\!\big(h_i \to h_i + J_{\text{prom}} f_i,\; J,\mu\big)}_{\sigma_{\text{prom}}=1}.
$$

**Conditioning the promoter open is exactly equivalent to boosting each operator's field by its
coupling weight $J_{\text{prom}} f_i$.** So no new DP is needed --- it is two evaluations of the
*same* array transfer matrix, the second with a per-operator field bump (the machinery a per-site
occupancy readout already uses). The count model boosts every operator uniformly by $J_{\text{prom}}$;
the nearest/decay models boost a profile. This also makes the reciprocity explicit: TF binding
stabilizes the open promoter and the open promoter stabilizes TF binding, by the same
$J_{\text{prom}}$ --- one number, symmetric, as thermodynamics requires.

## 7.5 Potency as a coupling

$J_{\text{prom}}$ **is** potency, in $k_BT$: the free-energy coupling between occupancy and
activation, comparable across ADs and doses, with Hessian error bars --- replacing the current
`compute_potency.py` curve-fit. Note a linear-in-energy coupling already yields a *saturating*
$P(\text{open}\mid n_{\text{tf}})$ (a logistic in $h_{\text{prom}} + J_{\text{prom}} n_{\text{tf}}$),
so it recovers the mechanistic saturation curve's shape without extra machinery --- the curve was
always a shadow of this coupling.

## 7.6 Deferred refinements

- **Promoter--nucleosome exclusion.** An open promoter physically forbids a nucleosome over the
  promoter window --- a *local geometric* coupling to the rod lattice. Dropped from the abstract node
  above (which couples only to TFs); add it when the promoter is placed on the same lattice as the
  rods.
- **Marginalizing v5 posterior.** As with the array, start from v5's hard promoter-state calls;
  integrate over the decoder's posterior only if observation noise bites.

---

# 8. Joint fitting across amplicons and samples

The old partition-function fit was, by design, **joint across the 0--8$\times$ TetO copy-number
amplicons**: one $h_{\text{tf}}$ and one $h_{\text{nuc}}$ shared across all of them, learned together
so a single parameter set had to explain every copy number at once. The MaxEnt reframe makes this the
base case and extends cleanly to many samples with sample-specific parameters. This is where the
reparametrization pays off most.

## 8.1 Joint across the copy-number series (the base case)

Each copy-number amplicon is just a different `Geometry` (different number/positions of operators),
same parameters. The objective is a sum over amplicons:

$$
\text{NLL}(\theta) = \sum_a N_a\Big[\log Z_a(\theta) - \theta\cdot\langle\phi\rangle_{\text{data},a}\Big],
\qquad
\frac{\partial\text{NLL}}{\partial\theta_k} = \sum_a N_a\Big(\langle\phi_k\rangle_{\text{model},a}
- \langle\phi_k\rangle_{\text{data},a}\Big).
$$

Each amplicon has its own $\log Z_a$ (its own transfer-matrix DP, because geometry differs) but all
share one $h$, one $\mu$, one $J$. This is exactly the old joint fit --- now **convex** (a sum of
convex functions), so a unique optimum and Hessian error bars instead of Nelder--Mead. The
copy-number series is what makes $h$, $J$, $\mu$ **identifiable**: different copy numbers give
different $\langle n_{\text{tf}}\rangle$/$\langle n_{\text{pairs}}\rangle$ that one shared
$(h,J,\mu)$ must reproduce simultaneously.

## 8.2 Sample-specific parameters (the same objective, indexed)

Split parameters into shared vs per-sample and sum over samples too:

$$
\text{NLL} = \sum_s\sum_a N_{s,a}\Big[\log Z_{s,a}(\theta_{\text{shared}},\,\theta_s)
- (\theta_{\text{shared}},\theta_s)\cdot\langle\phi\rangle_{s,a}\Big].
$$

Still convex in **all** parameters jointly. The gradient does the bookkeeping automatically: a
*shared* parameter's gradient pools over all data; a *sample-specific* parameter's gradient sees only
that sample's moments. One optimizer call fits everything. The natural split (same "earned
heterogeneity" discipline as Section 5 and 7):

- **Shared** (physical constants of the system): nucleosome fugacity $\mu$; array cooperativity $J$;
  intrinsic TetO affinity at fixed dox.
- **Sample-specific** (what the biology varies): the effective TF field $h_s$ (dox- and
  AD-dependent), and especially the promoter coupling $J_{\text{prom},s}$ --- *the potency*, which is
  the whole reason to index by AD.

This is Section 7.2's promoter hierarchy generalized: any parameter can be shared or indexed, in a
single convex objective.

## 8.3 Covariate-dependent couplings (a design matrix / GLM)

Because the energy is **linear** in $\theta$, any parameter can be made a linear function of sample
covariates and the fit stays convex (composition of convex $\log Z$ with an affine map):

- $h_s = h_0 + \beta\cdot[\text{dox}_s]$ --- fit a **titration/dose-response** with a couple of
  coefficients instead of one $h$ per dox.
- $J_{\text{prom},s}$ = per-AD fixed effect, or a function $g(\text{AD features})$.
- interaction terms (AD $\times$ copy number, dox $\times$ AD, ...).

It becomes a **GLM on top of the lattice gas**: define the feature vector, declare which coefficients
are shared vs indexed, sum the per-$(s,a)$ log-likelihoods, optimize. The Hessian gives confidence
intervals on the covariate coefficients too. To shrink noisy low-$N$ samples toward a group mean, add
a Gaussian prior on $\{h_s\}$ (a ridge penalty, still convex) --- a standard hierarchical /
random-effects model. Optional.

## 8.4 The one caveat: identifiability budget

Floating both $h_s$ and $J_{\text{prom},s}$ per sample needs enough constraints per sample. The
**copy-number series within a sample is what separates $h$ from $J$** --- a sample measured at only
one copy number cannot identify both. Practical recipe: share as much as is physically defensible
($\mu$, and usually $J$), float only what genuinely varies ($h$ with dox, $J_{\text{prom}}$ with AD),
and lean on the 0--8$\times$ series for identification. This is the same discipline as refusing
per-site $h_i$: never spend a parameter the data cannot constrain.

## 8.5 Why this is easier than under the old enumeration

The old code *could* fit jointly across amplicons (and did), but sample-specific and
covariate-dependent parameters would have been painful: everything was tied to enumerating a
per-configuration state table and a bespoke non-convex optimizer. The MaxEnt form turns the whole
thing into a standard convex regression --- define which coefficients are shared vs indexed, sum the
per-$(s,a)$ log-likelihoods, hand to L-BFGS --- with error bars on every parameter, including the
sample-specific ones and the covariate slopes.
