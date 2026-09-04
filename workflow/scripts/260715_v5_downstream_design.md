---
title: "v5 Binding Model — Downstream Output Format & Partition-Function Design"
author: "Design discussion notes"
date: "2026-07-15"
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
  - \usepackage{booktabs}
  - \usepackage{microtype}
---

# Scope

Working notes from a design discussion about how the `v5` HSMM binding classifier
(`classify_single_molecule_binding_v5_hsmm.py`) output should feed downstream analyses, and how
to reparametrize the thermodynamic (partition-function) model to consume it. Covers: the v5
output schema and its compatibility with old code, how to make single recurrent footprints
callable, a plan for a v5-native aggregator, the nucleosome meta-plot denominator, the role of
UNID footprints, and — the centerpiece — an Ising/MaxEnt reformulation of the partition function.

---

# 1. v5 output schema

Two files per (sample, amplicon).

**Main per-molecule file** (`{sample}.{amplicon}.single_molecule_classification.txt`, indexed by
`read_id`):

- `tfbs_1 … tfbs_K` — boolean, one per motif/TetO site
- `n_tf`, `n_nuc`, `n_unid` — counts
- `nucs` — semicolon-joined `start:end:dyad` intervals
- `unids` — semicolon-joined `start:end` intervals
- `log_likelihood` — total Viterbi path score

**Segments file** (`.segments.txt`, tidy long — one row per segment):
`read_id, seg_index, type, start, end, dyad_or_motif, seg_loglik`.

This is deliberately **not** coerced back to the old bin-based schema. There is no `idx` column
and no companion `valid_states_table.txt`, because the HSMM does not enumerate a global
microstate space.

**Old schema** (`v2` / `v2_carosversion`), for contrast: `tfbs_{k}`; ragged wide
`nuc{N}_present/start/end`; and `idx` into `valid_states_table.txt`.

---

# 2. Compatibility with existing downstream code

| Analysis | Works as-is? | Why |
|---|---|---|
| Binding histogram, `n_bound`, `tf_bound` | **Yes** | Both formats share `tfbs_{k}`; `aggregate_binding_model.py` uses only `df.filter(like='tfbs_').sum(axis=1)`. None of v5's extra columns contain the substring `tfbs_`. |
| Peripheral accessibility / remodeling | **Yes** | Computed from the matrix, not the classification — classifier-independent. |
| Nucleosome occupancy | **Adapt** | Old `nuc{N}_*` columns gone; use the `nucs` string or `.segments.txt`. |
| Partition function fit | **No** | `fit_partition_function_model_v3.py` needs `valid_states_table.txt` + the `idx` join; v5 produces neither by design. |
| UNID footprints | **New** | No equivalent in the old format. |

Decision: **do not** coerce v5 back to the old wide schema. The tidy `.segments.txt` is a strictly
better substrate; write new, non-legacy analyses against it.

---

# 3. Making a single recurrent footprint callable

Why two adjacent ATF1 sites get called UNID but a single site does not. To open a footprint
segment, Viterbi must accumulate enough protection evidence to pay the state's **entry cost**.
With the v5 emission parameters, one protected GpC is worth roughly

$$
\log\frac{P(\text{protected}\mid\text{UNID/TF})}{P(\text{protected}\mid\text{OPEN})}
\;\approx\;
\log\frac{0.95\cdot 0.9 + 0.15\cdot 0.1}{0.95\cdot 0.05 + 0.15\cdot 0.95}
\;=\;\log\frac{0.87}{0.19}\;\approx\;+1.5 .
$$

The entry costs are `start_cost_unid = 6.0` and `start_cost_tf = 1.0`. So:

- **Two adjacent sites** (about 3–4 protected GpCs, $\approx +4.5$ to $+6.0$) clear the UNID hurdle of $6.0$.
- **A single site** (about 1–2 GpCs, $\approx +1.5$ to $+3.0$) does not — it stays OPEN.

Two ways to fix, both on the same "position-specific prior" lever:

1. **Discovery pass** (`--discovered_footprints_file`): lowers the UNID entry cost *at that locus*
   from $6.0$ to $\sim2$, so a single site's evidence now clears it. Automated, no naming.
2. **Name it as a motif $\Rightarrow$ it calls as TF, not UNID** (preferred for a site you care
   about): entry cost drops to $1.0$, so **one protected GpC ($+1.5$) already clears it**; and its
   duration is the exact motif window rather than being bounded by `unid_min = 15 bp`
   (a single footprint narrower than 15 bp cannot form a legal UNID segment at all, but can be a TF
   of any width).

The "name everything eventually" endgame and the "make single sites callable" need are therefore
the *same* lever.

---

# 4. A v5-native aggregator

New script (operating on the main file + `.segments.txt`), tidy at all three levels, losing no
columns. Proposed per-molecule superset:

- **Keys:** `read_id, sample, amplicon, background, n_tfbs`
- **Binding (verbatim):** `tfbs_1…K`, `n_tf`, `tf_bound`
- **Nucleosome:** `n_nuc`, `nuc_bases_covered`, `frac_nuc`, raw `nucs`
- **UNID:** `n_unid`, `unid_bases_covered`, `frac_unid`, raw `unids`
- **Open:** `open_bases`, `frac_open` (sanity: fractions sum to 1)
- **Any footprint (the "is there a footprint here" query):**
  `footprint_bases_covered` $= \text{NUC}\cup\text{TF}\cup\text{UNID}$, `frac_footprint`
- **Kept from old:** `peripheral_accessibility` (from matrix), `log_likelihood`

Per-amplicon: means + count histograms of `n_tf`/`n_nuc`/`n_unid` (reuse the `n:freq`
histogram-string trick). Per-sample: as today.

Decision to make explicit: `frac_*` needs a per-molecule denominator. Use decoded span
(first-to-last GpC), not raw amplicon length, so molecules with different coverage compare fairly.

---

# 5. Nucleosome meta-plot denominator

Because Viterbi assigns a definite state to **every base of every molecule** (it fills through GpC
gaps), nucleosome occupancy is *defined everywhere* — no `-1`/no-info denominator correction like
on the raw accessibility matrix. So, per base $b$:

$$
\text{occupancy}(b) = \frac{\#\{\text{molecules with a NUC segment spanning } b\}}{\#\{\text{total molecules}\}} .
$$

The denominator comes straight from the segments file: every molecule has a full path, so every
`read_id` appears (with its OPEN/TF segments even if it has zero NUCs). One footnote: a molecule
with no GpC coverage in a region still decodes (toward OPEN) and counts as "not nucleosome,"
slightly deflating occupancy there — mask per-molecule uninformative regions only if it matters.

---

# 6. Role of UNID footprints

Endgame: name recurrent footprints (like TetOs — "this NRF, this ATF") so each becomes a named
entry in the TF file. But UNID is **not** made obsolete:

- It is the **residual/discovery channel** — "recurrent protection I have not accounted for yet."
  Deleting it forces the model to explain unexpected protection as OPEN (missed) or NUC
  (miscalled). It is how you find the next thing to name.
- The **"is there a footprint here" query** is not a UNID query — it is the per-bp occupancy of
  $\text{TF}\cup\text{UNID}$ ($\cup\text{NUC}$ if nucleosomes count). Named sites tell you *which*
  footprint; the union track tells you *whether* one is there. You want both, permanently.

So: build bulk/union tracks now, name sites incrementally into TF, and let UNID asymptote toward
"residual" rather than "primary output." It never has to reach zero to be done.

---

# 7. Partition function, Ising/MaxEnt style

## 7.1 The reframe

The literal microstate-enumeration approach *is* an Ising/MaxEnt model — just parametrized by
enumeration. A microstate $\sigma$ (which TetOs are TF-bound, where nucleosomes sit) is assigned an
energy that is **linear in a set of interaction features** $\phi_k(\sigma)$:

$$
E(\sigma;\theta) = \sum_k \theta_k\,\phi_k(\sigma),
$$

with, for example,

$$
\begin{aligned}
\phi_{\text{tf},i}(\sigma) &= \mathbb{1}[\text{site } i \text{ TF-bound}] && \text{(per-TF energy / field } h_i) \\
\phi_{\text{coop},i}(\sigma) &= \mathbb{1}[\text{sites } i, i{+}1 \text{ both bound}] && \text{(cooperativity } J) \\
\phi_{\text{nuc}}(\sigma) &= \#\{\text{nucleosomes}\}\ \text{or a region indicator} && \text{(per-nuc energy)} \\
\phi_{\text{tf}\cdot\text{nuc}}(\sigma) &= \mathbb{1}[\text{TF and nuc overlap}] && \text{(TF–nucleosome antagonism)}.
\end{aligned}
$$

The probability is exactly Boltzmann:

$$
P(\sigma\mid\theta) = \frac{1}{Z(\theta)}\,e^{-E(\sigma;\theta)},
\qquad Z(\theta) = \sum_{\sigma} e^{-E(\sigma;\theta)} .
$$

This is identical to "energy per TF, energy per nuc, optional cooperative term" — a Gibbs
distribution, i.e. the same object as a MaxEnt model and a generalized Ising model. Fields $h_i$
are single-site energies; couplings $J_{ij}$ are cooperativities.

## 7.2 The theorem that removes enumeration from *fitting*

The log-likelihood of observed molecules $\sigma_1,\dots,\sigma_N$ is

$$
\ell(\theta) = -N\sum_k \theta_k\,\langle\phi_k\rangle_{\text{data}} \;-\; N\log Z(\theta),
\qquad \langle\phi_k\rangle_{\text{data}} = \frac{1}{N}\sum_m \phi_k(\sigma_m).
$$

Using $\partial \log Z / \partial\theta_k = -\langle\phi_k\rangle_{\text{model}}$, the gradient is

$$
\frac{\partial\ell}{\partial\theta_k} = N\big(\langle\phi_k\rangle_{\text{model}} - \langle\phi_k\rangle_{\text{data}}\big),
$$

so the MLE is **moment matching**:

$$
\boxed{\;\langle\phi_k\rangle_{\text{model}}(\theta) = \langle\phi_k\rangle_{\text{data}}\;}
$$

Three consequences:

1. **The data enters only through feature averages** $\langle\phi_k\rangle_{\text{data}}$ —
   per-site occupancy, adjacent co-occupancy, mean nucleosome count — read directly off the v5
   per-molecule calls. No microstate assignment, no `idx`.
2. **The fit is convex.** $\log Z$ is convex (log-sum-exp), so $\ell$ is concave in $\theta$:
   unique optimum, gradient ascent suffices. The Hessian is $-N\,\mathrm{Cov}_{\text{model}}(\phi)$,
   giving parameter error bars for free.
3. **No multiplicity factors.** $P(\sigma)$ is per-configuration; no $\binom{N}{n}$ degeneracy
   bookkeeping. The "SMF measures microstates directly" intuition is preserved exactly.

## 7.3 Where configurations enter — $Z$ — and why 1D structure is cheap

Enumeration would only have been needed to compute $Z$ and $\langle\phi_k\rangle_{\text{model}}$.
The system is one-dimensional, and 1D lattice models have exact, cheap partition functions:

- **TetO array $\to$ transfer matrix.** With $E = \sum_i h_i n_i + \sum_i J\,n_i n_{i+1}$, $Z$ is a
  product of $2\times2$ matrices — $O(L)$, exact, for any number of sites. Marginals $\langle n_i\rangle$
  and co-occupancies $\langle n_i n_{i+1}\rangle$ (i.e. all $\langle\phi\rangle_{\text{model}}$) come
  from forward–backward on the chain. The combinatorial sum is done by the matrix product, correctly
  weighted, with no multiplicity bookkeeping.
- **Nucleosomes $\to$ 1D lattice gas of hard rods.** Extended excluded-volume objects on a 1D
  template — the same math as thermodynamic nucleosome-positioning models. $Z$ is a forward DP over
  positions with a fugacity per nucleosome and hard-core exclusion. Exact, $O(\text{length})$.
- Brute enumeration is still permitted when a model is small; the point is $Z$ is now a decoupled
  computational detail, not a prerequisite for touching the data.

You keep the literal partition function — still $\sum_\sigma e^{-E}$ — computed with a transfer
matrix / DP instead of materializing every microstate.

## 7.4 The one real modeling decision: discretizing nucleosomes

TF occupancy is clean (binary per TetO, from `tfbs_{k}`). Nucleosomes are emitted by v5 as
continuous-position segments and must be mapped to a discrete state variable:

- **Coarse (recommended to start): region indicators.** A few nucleosome "slots" (promoter nuc
  yes/no, one or two array-region slots) as binary features, coupled to TF/promoter state. Few
  parameters, trivial $Z$.
- **Fine: a proper lattice gas.** Bin the amplicon; a nucleosome occupies a run of about 147 bp of
  bins with exclusion; one energy per nucleosome. The positional degeneracy (a nuc can sit at many
  nearby positions) is integrated over automatically by the DP.

## 7.5 Minimal first model

1. **Features:** one TF field $h$ (shared, or per-site), one nearest-neighbor cooperativity $J$,
   one TF–promoter-nucleosome antagonism term. About 3 parameters — close to the old
   `3param_nuc` / `3param_tfcoop` flavors, but cleanly derived.
2. **Empirical moments** from v5 calls: per-site occupancy, adjacent co-occupancy, TF/nuc
   co-occurrence.
3. **$Z$ and model moments** via a transfer matrix over the array (+ a coarse promoter-nuc slot).
4. **Fit** by concave gradient ascent to moment-match; report energies with Hessian-derived CIs.

Two subtleties to revisit later: v5 calls are themselves inferred (could marginalize over v5's
posterior instead of hard calls if observation noise matters); and molecules with partial coverage
should contribute to a feature's moment only where they are informative.
