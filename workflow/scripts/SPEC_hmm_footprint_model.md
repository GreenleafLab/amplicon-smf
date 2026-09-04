# SPEC: HSMM footprint segmentation for single-molecule SMF (v5 / clean-room)

**Status:** design spec, not implemented. Hand this to Claude Code to build a new classifier as an
alternative to `classify_single_molecule_binding_v4.py`. Read `NOTES_v4_status.md`,
`NOTES_v4_nuc_model_changes.md`, and this repo's `CLAUDE.md` for background first.

---

## 0. Why (motivation)

The current v4 classifier **enumerates every global microstate** (all nucleosome-dyad × TF-binding
combinations), prunes overlaps, collapses observationally-equivalent states, and scores each molecule
against all of them. To fix mis-calls we've been bolting on penalty terms (TF-boundary penalty,
run-length Occam cap). This accretes hacks onto a model that is **site-independent Bernoulli with
channels combined by `max`** — it has no native representation of *runs* or *boundaries*, so
"a TF is a punctate protected blob bounded by accessibility" and "usually-but-not-always accessible
between adjacent TFs" cannot be expressed except as ad-hoc penalties.

A **hidden semi-Markov model (HSMM)** over each molecule expresses all of this natively:
segments (footprints) with explicit durations, position-dependent emissions, and transition
priors, decoded by Viterbi in **linear time per molecule with no global enumeration**. It replaces
enumeration + pruning + collapse + boundary-penalty + run-cap with one generative model.

**The new capability that makes this worth doing now (user's driver):** an **UNIDENTIFIED** footprint
state — a protected segment with *no* spatial anchoring (not nucleosome-length, not at a known TetO
motif). New data has real footprints in uncertain places that the enumerate-microstates approach
**cannot represent at all**. The HSMM's UNID state calls "there is a footprint here of unknown
identity," which is qualitatively new, not just a cleaner reimplementation.

---

## 1. Observations

Per molecule, from the `*.dedup.full_unclustered.matrix` (see `CLAUDE.md` matrix conventions):
a sequence of GpC calls `o_j` at genomic positions `p_j` (bp), `o_j ∈ {1 protected, 0 accessible,
NaN no-info}`. Emissions exist **only at GpC positions**; the latent segmentation is over **bp**
(footprint geometry is in bp). TFBS motif positions come from the `positions.txt` file
(`load_tfbs_positions`, as in v4) — TF segments are only allowed there.

---

## 2. Latent segment types (states)

| State | Where allowed | Duration | Emission (p_U = P(truly protected)) |
|---|---|---|---|
| `OPEN` | anywhere | geometric (variable linker) | low: `prob_unmeth_given_open` (~0.05); promoter window uses `prob_unmeth_given_open_promoter` |
| `NUC` | dyad anywhere (bp or 10-bp bin) | ~fixed 147 bp (see Q1) | position-dependent: `expit((d_edge − |pos − dyad|)/softness)` — reuse v4 `build_nuc_protection_matrix` math |
| `TF` | **only spanning an annotated TFBS motif** | fixed ≈ motif width + margin (~21–30 bp) | flat high: `prob_unmeth_given_tf` (~0.9) |
| `UNID` (new) | anywhere | flexible (broad; see Q2) | flat high: `prob_unmeth_given_unid` (~0.9) |

All emissions pass through the conversion model to get P(observed protected):
`p(obs=1 | state,j) = p_U · p_t_given_unmeth + (1 − p_U) · p_t_given_meth` (reuse v4 constants /
`invert_t_fraction`). NaN observations contribute nothing (skipped) to a segment's emission.

`UNID` is the explanation of last resort: it must carry a **higher start cost** than `NUC`/`TF`
(transition prior) so specific explanations win when they fit, and `UNID` only wins when neither a
nucleosome (wrong length/position) nor a TF (not at a motif) explains a protected stretch.

---

## 3. Transitions / structure (the priors)

- Allowed: `OPEN ↔ {NUC, TF, UNID}`. A protected segment is normally entered from and exited to
  `OPEN`.
- **Segment-start costs** encode parsimony (a per-segment penalty for `TF`, `NUC`, `UNID`), so the
  model does not over-segment. `UNID` start cost > `NUC`/`TF` start cost.
- **TF ↔ TF adjacency** (directly abutting, no intervening `OPEN`): allowed but penalized relative
  to `TF → OPEN → TF`. This is how "adjacent TFs are real, and the base between them is *usually*
  (not always) methylated/accessible" is captured — the linker `OPEN` path is usually higher
  likelihood, but direct `TF→TF` is available when the between-GpC reads protected. **This replaces
  the v4 boundary-penalty hack**; do NOT hard-code a bilateral-accessibility rule.
- **NUC ↔ TF adjacency**: allowed, governed by a tunable transition cost. Whether the model places a
  TF against a nucleosome becomes a cost-vs-likelihood tradeoff, not a hard rule. **Calibrate this
  cost using the `protection_streak_histogram.py` result** (see `NOTES_v4_nuc_model_changes.md`): if
  the ±TF streak analysis shows nuc-TF adjacency is real, keep this cost low/permissive; if it's an
  artifact, make it high.
- Hard constraint: a `TF` segment may only occupy an annotated TFBS motif position.
- Off-read / ghost segments: allow `NUC` (and `UNID`) segments to begin before the first observed
  GpC or end after the last (a nucleosome extending off the amplicon), as in the v4 ghost-position
  idea. Emission only counts observed GpCs, so a partially-observed segment is scored on what's seen.

---

## 4. Algorithm (HSMM Viterbi)

Standard segmental Viterbi over bp positions:
```
V[i] = max over (type, dur) of  V[i − dur] + trans_cost(prev_type → type)
                                 + start_cost(type) + dur_logprob(type, dur)
                                 + segment_emission(i − dur, i, type)
```
- Complexity `O(L · S · D_max)` per molecule (`L`≈600 bp, `S`=4, `D_max`≈200) → trivial; molecules
  are independent → embarrassingly parallel. **No global state enumeration, no collapse, no pruning.**
- **GpC-sparse emission precompute:** for flat-emission states (`OPEN`/`TF`/`UNID`), a segment's
  emission = `n_protected_in_[a,b] · log p(1|type) + n_accessible_in_[a,b] · log p(0|type)`, O(1)
  from cumulative counts of protected/accessible GpCs. For `NUC`, emission depends on the dyad
  center; precompute per-candidate-center or evaluate the sigmoid over the segment's GpCs.
- Backtrace → per-molecule segmentation (ordered list of segments with type, start, end, and for
  `NUC` the dyad, for `TF` the motif index).

---

## 5. Parameters & fitting

Emission params (`p_t_given_unmeth`, `p_t_given_meth`, `nuc_d_edge`, `nuc_softness`,
`prob_unmeth_given_tf`, `prob_unmeth_given_open`, `prob_unmeth_given_unid`) and structural params
(start costs, transition costs, duration distributions) are `ModelParams`-style config with CLI
overrides + optional YAML (mirror v4's two-pass argparse). Ship with sensible fixed defaults first;
add **EM/Baum-Welch** later (M-step estimates emissions by inverting observed T-fractions per assigned
segment type — reuse v4 `invert_t_fraction` / `run_em_mstep` logic; the E-step is
forward-backward or Viterbi-EM). Gate EM behind `--do_em` like v4. Do NOT block v1 on EM.

---

## 6. I/O & integration

- **Input:** identical to v4 — a `.dedup.full_unclustered.matrix`, a `positions.txt` (TFBS motifs),
  `--amplicon_name`. Reuse v4's `load_single_molecule_matrix`, `get_methyl_positions`,
  `load_tfbs_positions`, `filter_all_converted_reads`, `fix_missing_data`, `adjust_gcgs`.
- **Output:** one row per molecule. To stay drop-in for downstream (`aggregate_binding_model.py`,
  `fit_partition_function_model_v3.py`), emit the v4-compatible columns where they map
  (`nuc{k}_present/start/end`, `tfbs_{k}`, `log_likelihood`) **plus new** `unid_{n}_start/end`
  columns for unidentified footprints. Also emit a tidy long-format segmentation
  (`read_id, seg_index, type, start, end, dyad_or_motif, seg_loglik`) since UNID footprints don't fit
  the wide schema cleanly. Keep the per-read plotting (`plot_single_read` / `decorate_single_read_plot2`)
  and add UNID rendering (a distinct color/shape from nuc gray ovals and TF red ovals).
- Standalone CLI first (like v4); a Snakemake rule can come later.

---

## 7. Testing

- **Unit:** synthetic molecules with known segmentation (a lone TF with accessible flanks; two
  adjacent TFs with/without an accessible base between; a 147-bp nucleosome; a nucleosome + adjacent
  TF; a long protected stretch that is neither → should call UNID). Assert the Viterbi path.
- **Regression vs v4** on `data/binding_model_code_updating/260708_v4updates_fromclaude/` (L650,
  opJS4_6x): reproduce GOOD `...1102:21757:18698` (4 TF + 3 nuc); check the two hard molecules
  `...1119:19194:19277` and `...1109:4975:13836` behave sensibly given the histogram result.
- **UNID validation** on the new data with off-motif footprints (user to provide the dataset).

---

## 8. Open design questions (decide during implementation)

1. **Nucleosome duration:** fixed 147 bp, a small discrete set {147, 180}, or a distribution? (This
   is the Option-A/Option-B tension from `NOTES_v4_nuc_model_changes.md`, now expressible as an HSMM
   duration model.) Start fixed 147; make it a config knob.
2. **UNID duration & identifiability:** how to keep UNID from swallowing real nucleosomes/TFs. Lever:
   UNID start cost + a duration prior that disfavors nucleosome-length UNID (a nuc-length protected
   blob should be called NUC, not UNID) and disfavors TF-length UNID *at* a motif (should be TF).
   UNID should preferentially explain protected stretches of *atypical* length or *off-motif* position.
3. **NUC-TF adjacency transition cost** — set from the `protection_streak_histogram.py` result.
4. **Promoter region:** carry over v4's promoter-specific open probability (`--promoter_positions`,
   `prob_unmeth_given_open_promoter`).
5. **Missing data (NaN):** confirmed skipped in emissions; make sure long NaN runs don't let a
   segment span implausibly (cap via duration prior).
6. **Coordinate off-by-one:** the matrix drops the first reference base (known repo bug, see README
   To-Do) — keep the same coordinate convention as v4 so TFBS/positions line up.

---

## 9. Suggested layout

- `classify_single_molecule_binding_v5_hsmm.py` (new; do not modify v4).
- Reuse v4 helpers by import or copy: matrix loading, positions parsing, conversion/emission math,
  plotting. Factor the emission functions so v4 and v5 share them if practical.
- `ModelParams` dataclass mirroring v4 + new fields (`prob_unmeth_given_unid`, start/transition
  costs, duration params).
- Keep it standalone-runnable; add a smoke test like `test_v4_smoke.py`.

**Guiding principle:** the model should be *generative and local* — every rule we hand-coded as a
penalty in v4 (bilateral TF boundaries, run-length cap, parsimony) should instead fall out of
emissions + durations + transition costs. If you find yourself adding a special-case penalty, ask
whether it belongs in the transition/duration structure instead.
