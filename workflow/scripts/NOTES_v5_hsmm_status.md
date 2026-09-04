# Notes: `classify_single_molecule_binding_v5_hsmm.py` — status as of 2026-07-10

Conceptual walkthrough of the model (intuition + math): [`EXPLAINER_hsmm_v5.md`](./EXPLAINER_hsmm_v5.md).
Edge-case tuning log + how to diagnose a bad molecule: [`NOTES_v5_calibration_log.md`](./NOTES_v5_calibration_log.md)
(uses [`inspect_molecule_v5.py`](./inspect_molecule_v5.py)).

Implementation of [`SPEC_hmm_footprint_model.md`](./SPEC_hmm_footprint_model.md). Read that spec,
[`NOTES_v4_status.md`](./NOTES_v4_status.md), and [`NOTES_v4_nuc_model_changes.md`](./NOTES_v4_nuc_model_changes.md)
for background. This is the HSMM alternative to v4's enumerate-microstates approach.

## What was built (2026-07-10)

- **`classify_single_molecule_binding_v5_hsmm.py`** — segmental (semi-Markov) Viterbi over bp
  space, decoded per molecule, vectorized across all molecules at once. No global state
  enumeration / pruning / collapse. Four states: `OPEN`, `NUC`, `TF`, `UNID`.
- **`test_v5_hsmm.py`** — 8 synthetic known-answer tests. All pass (lone TF, two adjacent TFs w/
  OPEN linker, single nuc, nuc+TF, off-motif short UNID, di-nucleosome→2 NUCs, all-open, missing-
  data-skipped).

## Design decisions locked with user (2026-07-10)

- **NUC duration** = flexible Gaussian prior (`nuc_mode=147`, `nuc_sigma=25`, `nuc_min=100`,
  `nuc_max=210`). Adjacent nucs = two `NUC` segments (`NUC→NUC` transition), NOT one long NUC;
  linker-vs-no-linker decided by the inter-dyad GpC emissions. This is the dead-zone cure — a
  ~280bp protected stretch decodes as 2 NUCs (test confirms).
- **UNID = off-motif AND non-nucleosome-length only** (hard constraints): a UNID segment may not
  overlap any motif interval, and `unid_max` is clamped ≤ `nuc_min`. Higher `start_cost_unid`
  (6.0) makes it the explanation of last resort.
- **Output = new clean schema only** (NOT coerced to v4's bin-based `nuc{bin}_present` columns —
  that mapping is lossy: continuous dyads collide on bins, and the HSMM per-segment likelihood is
  a different quantity from v4's global-microstate likelihood). Three outputs:
  - `<output>`: one row/molecule — `tfbs_{k}` binary calls, `n_tf/n_nuc/n_unid`, `nucs`
    (`start:end:dyad;...`), `unids` (`start:end;...`), `log_likelihood`.
  - `<output>.segments.txt`: tidy long path — `read_id, seg_index, type, start, end, dyad_or_motif`.
  - `<plot>`: per-read traces, NUC gray / TF red / UNID purple.
- **EM deferred** (fixed params first). `ModelParams` is structured so a Viterbi-EM/Baum-Welch
  M-step (reuse v4 `invert_t_fraction`) drops in behind a future `--do_em`.

## Emission model (reused from v4, verbatim math)

`p(obs=protected | state,j) = p_U·p_t_given_unmeth + (1−p_U)·p_t_given_meth`. NaN/no-info GpCs
(matrix `-1`) contribute nothing (NOT the v4 `fix_missing_data`→0.5 soft-count; v5 skips them per
spec). Per-state p_U: OPEN=`prob_unmeth_given_open` (promoter override), TF/UNID flat high, NUC =
`expit((d_edge−|pos−dyad|)/softness)` (v4 `build_nuc_protection_matrix` math).

## Algorithm / perf notes

- Boundary grid = regular `grid_step` (default 5bp) lattice ∪ every motif edge (so TF segments
  align to motifs exactly). Flat-state emissions via cumulative-over-GpC arrays (O(1) window).
  NUC emission computed on the fly per (dyad, window).
- DP is `V[boundary, state, molecule]` with per-molecule backpointers; transitions depend only on
  previous segment type. Backtrace per molecule in a python loop.
- Smoke run (300 reads, L650 6x): seconds. Full-run timing being recorded (see below).

## Validation status

- Synthetic unit tests: 8/8 PASS.
- Real data smoke (L650 opJS4_6x, 300 reads): sensible distributions — most molecules 3–4 nucs,
  per-motif TF occupancy 5–26%, UNID rarely fires (last resort). Ghost nucs at read edges work.
- **TODO next session:**
  1. Inspect the 3 spec diagnostic molecules on full L650 6x run — GOOD `...1102:21757:18698`
     (expect 4 TF + 3 nuc), missed-TF `...1119:19194:19277`, extra-TF `...1109:4975:13836`.
     Compare to v4 outputs in `data/binding_model_code_updating/260708_v4updates_fromclaude/`.
  2. Calibrate `trans_nuc_tf` from `protection_streak_histogram.py` (SPEC §8.3).
  3. Tune start/transition costs against real molecules (defaults are first-guess).
  4. UNID validation on the new off-motif-footprint dataset (user to provide).
  5. Add EM (`--do_em`) once fixed-param behavior is trusted.
  6. Snakemake rule (later; standalone CLI for now, like v4).

## Parameter reference (all CLI flags / `ModelParams` fields)

Defaults shown. YAML `--config` supplies defaults; CLI flags override (two-pass argparse, like v4).
Run driver lives at `data/binding_model_code_updating/260711_v5_hsmm/log.sh`.

### I/O (required unless noted)
| Flag | Meaning |
|---|---|
| `--input` | `.dedup.full_unclustered.matrix` |
| `--output` | per-molecule wide table |
| `--segments_output` | tidy per-segment path (default `<output>.segments.txt`) |
| `--plot` | per-read trace PDF |
| `--positions` | `positions.long.txt` (TFBS motifs per amplicon) |
| `--amplicon_name` | which amplicon block in the positions file |
| `--config` | optional YAML of defaults |

### Preprocessing
| Flag | Default | Meaning |
|---|---|---|
| `--reads_to_use` | 0 (all) | subsample N reads (seeded) |
| `--filter_threshold` | 1.0 | drop reads with mean protection above this (v4 `filter_all_converted_reads`) |
| `--convert_ambiguous_gcgs` | '' | comma-separated `gcg_pos,fix_pos,...` pairs for GCG imputation (v4 `adjust_gcgs`; opt-in, amplicon-specific — see the GCG section of the pipeline CLAUDE.md) |

### Grid / geometry
| Flag | Default | Meaning |
|---|---|---|
| `--amp_width` | auto | bp length; auto = last GpC + max(nuc_max/2, d_edge) + 10 (ghost-segment room) |
| `--grid_step` | 5 | bp spacing of segment boundaries. Motif edges always added exactly. Smaller = finer dyad/boundary placement (dyad res ≈ step/2), quadratically slower. |

### Conversion model (shared by all emissions; same meaning as v4)
| Flag | Default | Meaning |
|---|---|---|
| `--p_t_given_unmeth` | 0.95 | P(observe protected \| truly protected). **Pipeline runs pass 0.99.** |
| `--p_t_given_meth` | 0.15 | P(observe protected \| truly accessible). **Pipeline runs pass 0.05.** |

### Per-state "truly protected" probabilities (`p_U`)
| Flag | Default | Meaning |
|---|---|---|
| `--prob_unmeth_given_open` | 0.05 | OPEN (accessible linker) |
| `--prob_unmeth_given_tf` | 0.9 | TF footprint (flat) |
| `--prob_unmeth_given_unid` | 0.9 | UNID footprint (flat) |
| `--promoter_positions` | 0,0 | `lo,hi` promoter window (0,0 = off). Pipeline passes `75,175`. |
| `--prob_unmeth_given_open_promoter` | 0.5 | OPEN prob inside the promoter window (chromatin there is more open) |

### Nucleosome footprint + duration
| Flag | Default | Meaning |
|---|---|---|
| `--nuc_d_edge` | 65 | protection half-max distance from dyad (~half the footprint); v4 sigmoid |
| `--nuc_softness` | 5 | logistic edge softness (larger = softer/broader edge) |
| `--nuc_mode` | 147 | peak of the Gaussian duration prior |
| `--nuc_sigma` | 25 | width of the duration prior |
| `--nuc_min` | 100 | min NUC segment length (also the floor on abutting di-nuc dyad spacing) |
| `--nuc_max` | 210 | max NUC segment length |
| `--nuc_min_prot_gpcs` | 3 | min observed-protected GpCs per NUC segment (evidence floor; blocks spurious 1–2-GpC edge/+1 nucs enabled by the low-edge ghost pad). 0 disables. Added 2026-07-19. |

Note: `nuc_d_edge`/`nuc_softness` (the *emission* footprint) are decoupled from `nuc_min/max`
(the *segment duration*). A long protected stretch is tiled by multiple NUC segments (di-nuc) rather
than one over-wide footprint — this is why the flexible duration cures the v4 dead zone.

### TF / UNID geometry
| Flag | Default | Meaning |
|---|---|---|
| `--tf_margin` | 2 | bp added each side of a motif for the TF segment/emission window (matches v4 −2/+2) |
| `--unid_min` | 15 | min UNID length |
| `--unid_max` | 90 | max UNID length; **clamped in `__post_init__` to ≤ `nuc_min`** so UNID can never be nucleosome-length |

### Structural costs (all POSITIVE numbers, subtracted as log-penalties)
Emission scale for calibration intuition: a *correctly* explained GpC contributes ≈ `log(0.99)` ≈
−0.01; a *mis*-explained one ≈ `log(0.01)` ≈ −4.6 (with 0.99/0.05). So a cost of ~2–6 trades off
against roughly 1–2 mis-explained GpCs.

| Flag | Default | Meaning |
|---|---|---|
| `--start_cost_open` | 0.5 | per-OPEN-segment cost |
| `--start_cost_nuc` | 2.0 | per-NUC parsimony |
| `--start_cost_tf` | 1.0 | per-TF parsimony (was 2.0 — lowered so an anchored nucleosome stops swallowing bound TFs; see `NOTES_v5_calibration_log.md` 2026-07-10a) |
| `--start_cost_unid` | 6.0 | per-UNID parsimony — **highest, makes UNID the explanation of last resort** |
| `--trans_nuc_nuc` | 1.0 | adjacent nucleosomes with no linker (di-nucleosome) |
| `--trans_tf_tf` | 1.0 | directly adjacent TFs (no linker) |
| `--trans_nuc_tf` | 1.0 | nuc directly against a TF — **the knob to CALIBRATE from `protection_streak_histogram.py`** (SPEC §8.3) |
| `--trans_unid_adj` | 1.0 | UNID directly against a nuc or TF |

Transition rules baked in (not flags): `OPEN→OPEN` and `UNID→UNID` are disallowed (a segment is
already maximal); all other same-type/cross-type transitions are allowed at the costs above; TF
segments may only span an annotated motif; UNID may not overlap any motif.

### Plotting
| Flag | Default | Meaning |
|---|---|---|
| `--reads_to_plot` | 10 | random reads to draw (seeded) |
| `--individual_reads_to_plot` | '' | comma-separated read_ids to draw explicitly |
| `--bulk_reads` | 1000 | reads in the bulk sanity panel (first PDF page): observed data (top) vs predicted state (bottom, 4 colors), same reads, sorted by prediction, aligned on the bp axis. 0 disables. |

### EM (Viterbi / hard EM) — implemented 2026-07-10
| Flag | Default | Meaning |
|---|---|---|
| `--do_em` | off | fit params from data by hard EM (decode ⇄ re-estimate) |
| `--em_max_iters` | 8 | max EM iterations |
| `--em_min_obs` | 100 | min observations to update a param |
| `--em_log` | `<output>.em_log.tsv` | per-iter convergence log (applied values + `est_*` raw estimates) |
| `--em_fit` | safe set | comma list of params EM may update. Default = `prob_unmeth_given_tf,prob_unmeth_given_unid,nuc_mode` |

**Only the safe set is fit by default.** `nuc_d_edge`, `nuc_softness`, `prob_unmeth_given_open`,
`nuc_sigma` are EXCLUDED — hard (Viterbi) EM degenerates on them (footprint → wide top-hat, sigma →
collapse); see `NOTES_v5_calibration_log.md` 2026-07-10b. They're still estimated and logged as
`est_*` for diagnostics but not applied. Proper footprint-shape fitting needs soft/forward-backward
EM (future) or external 0x-streak calibration. Nucleosome footprint from 0x data ≈ **140 bp**.

## Updates 2026-07-13

**0. Footprint discovery / position-specific UNID prior (empirical-Bayes across molecules).** A
per-molecule Viterbi can't see that a 2-GpC (or even 1-GpC) footprint is REAL because it recurs
across molecules. New optional pass borrows that cross-molecule strength cheaply (no joint fit):
- `--discover_footprints off|2pass|3pass` (default `off` = unchanged 1-pass scalar UNID cost).
  `2pass` = discover from raw bulk. **`3pass` (recommended)** = decode once → nucleosome(+TF)-
  DECONTAMINATED bulk (`per_gpc_protection(paths, exclude_states=(NUC,TF))`) → discover → re-decode.
  Decontamination matters: raw bulk is nucleosome-contaminated (on BAX_7x, raw flagged 5 sites, 4
  nucleosomal; decontam flagged only the real ones).
- Discovery (`discover_footprints_from_bulk`) is **local-background-relative** (an absolute 0.5 is
  too strict): flag a GpC if decontam frac ≥ `--footprint_abs_floor` (0.25) AND exceeds the local
  bg (`bg_pct`-ile within ±`--footprint_bg_window` 80bp) by ≥ `--footprint_rel_delta` (0.20).
  `min_gpcs=1` (a lone strong-over-background spike counts), off-motif only, span ≤ unid_max, padded
  to ≥ unid_min. Each interval carries an evidence **score = Σ per-GpC excess-over-bg** (folds in
  strength AND #GpCs).
- The discovered UNID start cost is **graded** smoothly by score: `np.interp(score, [rel_delta,
  --footprint_score_hi(1.0)], [start_cost_unid(6), --unid_discovered_start_cost(2)])`. So a strong
  site gets the full discount, a marginal one barely any. Applied per-segment in `viterbi_decode`
  via `discovered_footprints=[(lo,hi,cost),...]` (a UNID whose CENTER is in an interval uses `cost`).
- **Always-on diagnostic** (`plot_footprint_diagnostic`, 1st page of `--plot`): per-bp NUC/TF/UNID
  occupancy (top) + raw & decontaminated bulk with 0.3/0.4/0.5 guides + discovered shading (bottom).
- Smoke (BAX_7x, 3pass): discovers `(171,188,cost3.5)` + `(251,267,cost2.0)` — the ~0.40 footprint
  the absolute 0.5 missed is now caught at a moderate discount; the strong 0.8 one at full discount.
- **POOLED per-promoter footprint vocabulary (IMPLEMENTED 2026-07-13).** Per-(sample x amplicon)
  discovery gave noisy small-N bulks. New `discover_footprint_vocabulary.py` POOLS the
  nucleosome(+TF)-decontaminated bulk across ALL samples AND copy-number variants sharing a promoter
  (e.g. all BAX_opJS4_{0..8}xTetO), discovers the off-motif footprint vocabulary ONCE per promoter,
  writes `lo,hi,cost` blocks; the classifier applies it via `--discovered_footprints_file` (no
  per-amplicon discovery). Mirrors the TF-motif design (fixed vocabulary + per-sample occupancy ->
  comparable across conditions). COORDINATE-SAFE without alignment: pools ONLY positions `< the
  promoter's first-TetO` (upstream promoter is byte-identical at the same coords across variants --
  RC encoding + array grows downstream; verified BAX first-TetO=481). Wired as STAGE 1 of
  `260521_P110_moreprom_vp48_drugs/260713_modelrunning_v5_attempt1/log.sh` (STAGE 2 = per-amplicon
  decode with `--discovered_footprints_file`). Test (BAX, 9 variants x 2 samples, 9162 reads pooled):
  footprints (36,51,c5.9)(171,188,c3.6)(216,231,c5.9)(251,267,c2.0) -- strong sites full discount,
  marginal ones near-full cost (won't over-call). Discipline: verify RC-conserved coords per
  promoter by sequence; condition-blurring acceptable (prior = WHERE, decode = HOW MUCH).
  Details: [[project_pooled_footprint_vocabulary_todo]].
- **OPEN (reminder for tomorrow):** densely-bound TetO ARRAYS are undervalued — bound TetOs with
  protected linkers get swallowed into nucleosomes (contiguous-protection identifiability limit).
  UNID can't fix it (motif exclusion). Fix = symmetric prior that favors TF-over-NUC at
  recurrently-bound TetOs (discover from decontam-over-motif bulk); carries over-calling risk,
  needs its own diagnostic. Not yet built.

**1. Nucleosome footprint tied to segment length (removed `nuc_d_edge`).** The emission sigmoid's
half-max is now tied to the NUC segment half-length (`p_U = sigmoid((half_len - |pos-dyad|)/
softness)`, `half_len=(end-start)/2`), so footprint FWHM == segment duration. Extent/size
variability is carried ENTIRELY by the duration prior; `nuc_softness` is a fixed breathing
constant (not EM-fit). This removed the old `d_edge`<->`nuc_sigma` degeneracy (two knobs both
encoding single-nucleosome extent) and the "footprint narrower than segment -> spurious tail
segment" bug (2026-07-10b). `nuc_d_edge` deleted (field, CLI, EM). EM safe set unchanged
(`p_tf,p_unid,nuc_mode`); `fit_nuc_sigmoid` removed. 8/8 synthetic tests still pass. Validation on
the 3 diagnostic molecules: 4x lone-nuc (07-10b) now a clean single NUC (fixed); 6x GOOD unchanged
(4 TFs, TetO1 has an accessible gap so unambiguous); **7x (07-10a) TetO1 now folded into the
upstream nucleosome instead of called as a TF -- user CONFIRMED this is the biologically correct
call** (TetO1 shouldn't be bound there; the contiguous-protection no-gap case). Net: the change is
kept. Image: `data/binding_model_code_updating/260711_v5_hsmm/260713_revertant_failure_7x_TetO1_swallowed.png`.

**2. Positions-file right-flank off-by-one (verified, harmless, documented).** In
`positions.long.txt`, each TetO window `[start,end)` drops the operator's RIGHT-flanking GpC
readout C (it lands on the half-open exclusive `end`); the left flank is kept. Rescued by
`tf_margin=2` at runtime (shortfall 1bp < margin 2bp), so results are unaffected -- keep
`tf_margin>=2`. Cause is `convert_fa_to_positions_for_script.py`'s RC + right-exclusive end. For
new positions files use `--l_offset 0 --r_offset 0` and re-verify flank GpCs land in
`[start-tf_margin, end+tf_margin)`. Full detail in the memory
`project_positions_file_rightflank_offbyone`.

## TetO footprint + hemi/flank-spike investigation (2026-07-13)

Pseudobulk protection around TetOs (L648+L650, opJS4 arrays), aligned to operator center. Scripts
in the session scratchpad; plots in `data/binding_model_code_updating/260711_v5_hsmm/260713_*.png`.

- **TF footprint is TIGHT — no real extension beyond the immediate flanks.** Selecting molecules by
  "both flanks protected" showed protection extending to +/-20bp, BUT that is nucleosome
  contamination (a nuc protects both flanks symmetrically, mimicking binding). Conditioning instead
  on the MODEL's TF call (excludes NUC-explained molecules), and especially on TF-called +
  neighbors-not-bound, the protection OUTSIDE the operator flanks collapses to the unbound
  baseline. So the genuine rTetR footprint = operator + the two immediate flank GpCs (~+/-10.5bp),
  which `tf_margin=2` (+/-11.5bp window) already captures. **Do NOT widen the TF window / add a
  decaying flank zone** — an earlier suggestion to do so was based on the nuc-contaminated view.
- **Hemi is common + flank spikes are enriched (nuc-immune).** 76.7% both-flanks / 6.3% neither /
  **17.0% hemi** (one flank only). Isolated single-GpC protection spikes (a lone 1 with both
  GpC-neighbors accessible -> nucleosome-excluded by construction) occur at **1.85x** the background
  rate exactly at operator flanks (4.07% vs 2.20%; NEAR==FAR), highly significant. The conversion
  asymmetry (false-protection 0.05 >> false-accessible 0.01) already makes the model correctly lean
  hemi -> not-bound.
- **OPEN QUESTION gated on a control:** is the flank-spike excess real partial/asymmetric binding, or
  flank enzymatic/sequence bias? Decisive test = rerun the FLANK/NEAR/FAR spike test + the
  model-conditioned pseudobulk on a **-dox / no-rTetR sample** (user locating). If the enrichment
  vanishes -> bias -> keep discarding hemis (current behavior). If it persists -> real -> consider a
  soft per-flank TF emission (probabilistic flanks) that can call a real hemi without over-calling.
  No emission change until then.

## Updates 2026-07-14

**Motif overlay on the plots (`--motif_track_file`).** New OPTIONAL cosmetic overlay of known
TF/promoter motifs onto the always-on diagnostic page (top strip on the occupancy axis) and the
per-read trace pages. Fully backward-compatible: default `None` => no strip, decode byte-identical;
an amplicon absent from the track file draws no strip. **Never affects decoding** — plotting only.
- File format (block layout like `positions.long.txt`): `>amplicon` then `lo,hi,label[,color]` in
  the SAME matrix-column/bp frame the model uses (color optional hex, default gray).
- Loader `load_motif_track()`; drawer `draw_motif_strip()` (blended-transform bar strip above the
  axis + faint in-plot guide band + vertical family-colored labels).
- PRODUCER (for the opoBD9 panel): `export_motif_track.py` in the P110 project dir — calls that
  project's `motif_annotations.py` (curated view = manual/literature/encode) and writes a
  per-amplicon track. Coordinate frame VERIFIED (motif side does `col = prom_end - pos`, the RC
  reversal; landmark: BAX prom_end=400 abuts the array at col 481; both sides land in the same
  post-RC bottom-strand column frame, so RC is handled consistently — one reversal per side).
  opoBD9-specific; a different amplicon panel (opJS45/ad_smf) would need its own motif source.
- Wired into the P110 v5 driver (`260521_P110.../260713_modelrunning_v5_attempt1/log.sh`) and a new
  P116 driver (`260629_P116_8ADsmf/260714_modelrunning_v5hsmm/log.sh`, same opoBD9 library so it
  shares P110's positions + track). Bonus finding: discovered UNID footprints (171-188, 251-267 on
  BAX) coincide with ATF1/NRF1 motifs — direct validation of the empirical-Bayes discovery.
  Details: memory `project_motif_overlay_v5_plots`.

## Updates 2026-07-19 — symmetric LOW-edge ghost padding + NUC evidence floor (+1 nuc fix)

**Problem (confirmed by inspection).** The auto `amp_width` padded only the HIGH side (`max_gpc +
nuc_max//2+10`); the tiling floor was bp 0 (~the first GpC). So a **+1 nucleosome** sitting just
downstream of the TSS — i.e. at the LOW-numbered edge after RC, contiguous with bp 0 — could not be
placed as a NUC: once the edge truncates its footprint below `nuc_min=100`, no NUC segment fits and
the protection fell to **UNID**. That UNID then survived NUC/TF decontamination and got promoted
into the pooled footprint vocabulary — hence a spurious low-edge footprint (`~26–51`) in ~9/16
opoBD9 promoters (`260713_modelrunning_v5_attempt1/footprint_vocabulary.txt`; H4C5's ONLY entry was
`32,47`). Measured on H4C5 0x (12,335 edge-protected molecules): **30.5% had the edge called UNID,
and 100% of those were <100 bp truncated** (median observed extent 76 bp) — the call flips UNID↔NUC
purely on whether the truncated stretch clears `nuc_min`.

**Fix 1 — symmetric ghost padding, fully encapsulated in `viterbi_decode`.** Decoding now happens
in an internal frame shifted RIGHT by `left_pad = nuc_max//2+10` (mirrors the existing high-side
ghost room), so a NUC can anchor its dyad in the `[0,left_pad)` ghost zone (no GpCs there → no
emission cost) and represent a nucleosome whose center is off the low edge. GpC / motif / promoter
(`promoter_lo/hi` via `dc_replace`) / discovered-footprint coords are all shifted in; on backtrace,
segment coords are shifted back and **CLAMPED to `[0, orig_amp_width]`** — so ALL OUTPUTS stay in
the caller's original bp frame with **no negative coordinates** (segment start/end AND the NUC dyad
are clamped ≥0; an off-edge +1 nuc reports e.g. `NUC(0,65,dyad=0)`). Because it's inside
`viterbi_decode`, every caller (main decode, EM re-decode, 3-pass discovery re-decode,
`inspect_molecule_v5.py`, `discover_footprint_vocabulary.py`) gets it transparently — no other code
changed. Caveat: clamped edge dyads pile at 0, so nucleosome PHASING off dyads should treat dyad==0
(abutting bp 0) as an edge marker, not a measured position (the true center is off-array/unmeasurable).

**Fix 2 — NUC evidence floor `nuc_min_prot_gpcs` (default 3).** The generous pad reintroduces a
risk: a NUC can hide its unsupported length off-screen, so 1–2 protected edge GpCs could call a
spurious +1 (min real-DNA extent = `nuc_min − left_pad` < 0 at `left_pad=115`). The pad is a
fragile lever for this (the bound `left_pad < nuc_min − g₂` is per-amplicon, GpC-spacing-dependent),
so instead a NUC segment must contain **≥ `nuc_min_prot_gpcs` observed-protected GpCs** (per
molecule), enforced vectorized in the NUC branch via a leading-zero cumsum of `obs1`. Robust and
amplicon-agnostic; keeps the pad symmetric at 115. Applies to ALL nucs but on these GpC-dense
amplicons essentially only bites the edge. CLI `--nuc_min_prot_gpcs` (0 disables).

**Validation (H4C5 0x, 3 samples; scripts in the session scratchpad).**
- 8/8 synthetic tests still pass.
- Edge-UNID 30.5% → **1.4%**; truncated-edge molecules now decode NUC; 65% of edge NUCs anchor at
  start==0 (dyad off-edge, the intended ghost behavior).
- NUC/TF-decontaminated bulk at the low-edge GpCs (26/35/44/48) drops from raw **0.45–0.63** to
  **0.00–0.11** (< the 0.25 discovery floor) → the spurious `32-47` footprint will **no longer be
  discovered**.
- Evidence floor: short edge blips (next GpCs accessible) with 1 protected GpC went 33%→**0%** NUC,
  2 GpCs 68%→**0%**; ≥3-GpC calls unchanged (92%/100%). Interior mean nucs/molecule 3.64 → 3.55.

**Why the floor is 3, not 4 (threshold discussion 2026-07-19, user's call = 3).** Full data +
scripts: `data/binding_model_code_updating/260719_edge_nuc_padding/` (`RESULTS.md`).
- The INTERIOR naturally never calls a nuc on <4 protected GpCs (0 of 42,620 across 4 promoters) —
  the emission-vs-cost balance enforces it, since a full ~147 bp span at ~1 GpC/10 bp is ~14 GpCs.
  I first proposed matching that (N=4), but that comparison is INVALID: edge nucs are boundary-
  TRUNCATED and physically can't show as many GpCs, so holding them to the fully-observed interior
  bar is wrong. Truncation + Occam + "nucs live in the gene body" ⇒ a short contiguous protected run
  abutting the boundary is probably a truncated nuc and deserves a LOWER bar. A GLOBAL floor of 3
  gives position-aware behavior for free: interior self-selects ≥4, edges need 3.
- **3 is independently endorsed by the data.** The PRE-pad low-edge UNIDs ARE the truncated +1 nucs;
  they span 47–80 bp (median 65) and carry ≥3 protected GpCs EVERY time (min 3, median 4, none <3).
  Mechanistic reason: UNID start cost 6.0 / ~2.3 nats-per-protected-GpC ≈ 2.6 ⇒ the UNID channel
  already needed ~3 protected GpCs to fire. So a second, independent part of the model had already
  picked 3. Before/after is consistent: pre-pad 1–2 protected edge GpCs → OPEN, ≥3 → UNID (spurious);
  post-pad floor=3: 1–2 → OPEN (unchanged), ≥3 → NUC (correct, decontaminated out of discovery).
- Fraction-protected does NOT discriminate (interior 0.96 vs low-edge 0.97) ⇒ raw COUNT is the lever.
- The gene-body/downstream region is CONSTANT across promoters, so the edge GpC coords are conserved
  and "≥3 protected GpCs" ≡ a fixed geometric "sufficiently over the amplicon" overlap — i.e. the
  count floor IS the "pad penalty" idea in evidence units; no separate penalty mechanism needed. The
  exact `left_pad` is therefore not critical (kept symmetric at 115). CAVEAT: "interior min ≈4 /
  low-edge min =3" is calibrated to this panel's GpC density (~1/10 bp); re-derive for a much sparser
  panel. A graded pad penalty was considered and rejected (hard count already separates cleanly).

**Not yet done / open:** (a) regenerate the pooled vocabulary + full decode on the real panels
(`260713_modelrunning_v5_attempt1/log.sh` STAGE 1+2 — heavy, run via sbatch, not login node) to
confirm the spurious low-edge entries vanish across all promoters end-to-end; (b) the separate
"treat the +1 nuc as functional" idea is a DOWNSTREAM label, not a decode change — identify the NUC
segment abutting bp 0 in the planned v5-native aggregator as a `plus1_nuc` track. `score_segmentation`
(diagnostic-only, used by `inspect_molecule_v5.py --compare`) was left un-padded — it scores
user-supplied original-frame segments; only matters for hand-scoring an explicitly off-edge NUC.

## Updates 2026-07-20 — allow NEGATIVE NUC dyads (off-low-edge +1 nucleosomes) + 16-promoter check

**16-promoter cross-check of the 2026-07-19 fix (COMPLETED).** Ran the three 2026-07-19 checks
(RESULTS.md) panel-wide (all 16 opoBD9 promoters, 0xTetO, 3 samples each) via
`data/binding_model_code_updating/260719_edge_nuc_padding/check_16promoters.py`
(summary: `check_16promoters.summary.txt`). Findings:
- **The UNID→NUC flip works everywhere.** Across all 16 promoters the mostly-protected low edge now
  decodes as NUC, not UNID (edge cover@40 e.g. 3435 NUC / 26 UNID; edge-NUC `start==0` frac 0.53–0.89
  = anchored off the low edge as intended).
- **Interior floor impact negligible panel-wide.** `frac<3 protected GpCs` = 0.000–0.004 across all
  promoters — floor=3 barely touches the interior (a handful of rare 2–3-GpC interior nucs exist, so
  "interior min ≥4" is *almost* but not strictly universal; substance holds).
- **Pre-pad low-edge UNID min protected GpCs = 3 or 4 in every promoter** (never <3) — the N=3 floor
  choice is endorsed panel-wide.
- **Spurious low-edge footprint:** my checker's decontam-vs-0.25 proxy said 7/9 clearly gone, 2/9
  (LYL1 0.25, PGK 0.28) borderline — BUT the checker was deliberately conservative (applied the OLD
  spurious vocab during decode, biasing toward UNID, and used only 0x). The definitive test is a
  fresh STAGE-1 regeneration (decodes with NO vocab, pools across all copy-number variants) — see
  `footprint_vocabulary.postfix.txt` diff in the results dir for the real answer.

**NEGATIVE NUC dyads now allowed (was: clamped ≥0 → all +1s piled at dyad==0).** In `viterbi_decode`
the reported NUC dyad is the true modeled center `0.5*(start_bp+end_bp)` in the ORIGINAL frame; for an
off-low-edge +1 nucleosome this is < 0 and is now KEPT negative (only the high side is capped, a
defensive no-op since `orig_amp_width` already contains the high ghost pad). Rationale: clamping all
+1 dyads to 0 discarded the real ordering their observed (interior) protection boundary carries;
different +1s with different interior extents now get different (more/less negative) dyads. **Safe** —
the dyad is a reported VALUE, never an array/column index (only segment start/end and column coords
must stay ≥0, and still are), so no negative-index wraparound. **Caveat unchanged in spirit:** an
off-edge dyad is a *prior-regularized estimate* (pinned by the interior boundary + duration prior, NOT
measured on the off-edge side) — good for relative ordering / "how buried is this +1," NOT for absolute
nucleosome phasing. `paths_to_tidy` was also fixed (the old `anchor >= 0` guard would have written a
negative NUC dyad as `''`; now NUC always emits its dyad, TF keeps the motif-index guard). **+1
identification is `start == 0`** (unchanged, robust); `dyad < 0` / `dyad < tss` are secondary; `dyad
== 0` is no longer special. 8/8 synthetic tests + a new negative-dyad assertion pass.

## Relationship to the v4 inter-TF work

The uncommitted v4 working tree (branch `v4-nuc-model-updates`) has the TF-boundary penalty
(OUTER flanks only) + run-length Occam cap + Option-A collapse. It does NOT have an inter-TF-gap
penalty — inter-TF gaps are exempt (`classify_single_molecule_binding_v4.py:90,372`). In v5 the
inter-TF "usually-but-not-always accessible" behavior falls out for free: adjacent bound TFs
decode as `TF → OPEN(short linker) → TF`, and the linker OPEN's emission naturally prefers an
accessible linker GpC but tolerates a protected one (test `two_adjacent_tfs` confirms).

## Updates 2026-08-03 — stable site names in the wide output

`tfbs_{k}` is POSITIONAL and not stable: the positions file sorts promoter sites BEFORE the array
(promoter = LOW cols, TetO = HIGH), so inserting one promoter site renumbers every TetO and silently
invalidates anything assuming `tfbs_1..6 == TetO1..6`. Added `tfbs_site_names()` (sanitize +
ordinal-suffix duplicates → `TetO1..TetO6`, with a collision guard), `is_teto_site()`,
`write_tfbs_index()`. `paths_to_wide` now also emits `site_{name}` bools — deliberately a DIFFERENT
prefix so existing `df.filter(like='tfbs_')` calls don't double-count — plus `n_tf_teto` /
`n_tf_other`. `paths_to_tidy` gains `motif_name`. New CLI: `--tfbs_index_output`,
`--teto_name_prefix`; sidecar `<output>.tfbs_index.txt` maps index → (name, lo, hi, is_teto).

**⚠️ Use `n_tf_teto`, not `n_tf`, for "how many TetOs are bound".** `n_tf` counts every named site and
grows as promoter footprints (TATA/BRE/Inr/…) get added to the positions file. Same trap for
`df.filter(like='tfbs_').sum(axis=1)`. Downstream code should key on the NAME, never the index.

## Updates 2026-08-04 — region-anchored per-molecule columns (`--regions`)

The per-molecule output was coordinate-free (counts + raw intervals); nothing in it knew what a TSS or
promoter was. Added region-anchored columns, computed by intersecting the Viterbi segments with the
named windows in a `build_promoter_regions.py` region file.

- **Headline:** `tss_nuc` (a NUC segment covers the TSS column itself — chosen over "overlaps
  `TSS_Inr`" because it needs no window-width convention), `tss_footprint` (UNID/TF over the TSS),
  `tss_open` (neither). Mutually exclusive and exhaustive — asserted on real data. The three-way split
  matters: a TSS under a PIC/Inr footprint is *protected* but mechanistically the opposite of closed,
  so a single `tss_accessible` bool would lump it with either nucleosomal or naked molecules.
- **`promoter_nuc_gt50`** = >50% of the promoter INSERT covered by NUC, + `nuc_bases_promoter`,
  `frac_nuc_promoter`, `promoter_len`, and an optional absolute-bp column (`--promoter_nuc_min_bp`).
- **Generic:** `{region}_bases_{type}` / `_frac_` / `_any_` for every region × {NUC,TF,UNID,OPEN}, so
  `plus1_nuc_frac_NUC`, `TATA_frac_UNID` etc. come for free.
- **Two entry points, one implementation.** `--regions` on the classifier (inline, when re-running
  anyway) and the standalone `annotate_molecule_regions.py` (post-hoc, against segments already on
  disk — no re-decode). Both call `annotate_one()` on the same tidy segments table; verified identical
  across all 82 shared columns × 681 molecules.
- **`build_promoter_regions.py` now emits a `promoter` region** = the insert block
  `[prom_start, prom_start+Length-1]`, plus a `# prom_start=… prom_end=… prom_len=…` comment.
  `prom_end` is the INCLUSIVE last column, matching the region line and Layer B's inclusive
  `min(b,hi)-max(a,lo)+1` arithmetic — one convention in the file. Regenerating is purely additive
  (diffed against `260720_opoBD9_promoter_regions.txt`: no drift in any existing region or tss_col).
  New file: `260804_opoBD9_promoter_regions.txt`.
- **Copy-number variants:** regions are matched by PROMOTER name, so the single `_opJS4_6xTetO` entry
  covers the whole TetO series. Verified valid — the array starts at a FIXED low coord (JUNB 481,
  RPS9 417; the 64 bp difference = the 264 vs 200 bp promoter length, with an identical 81 bp gap) and
  grows upward, so everything below it is coordinate-identical. The annotator warns if any region
  extends past a variant's array start.

**⚠️ CAVEAT — `promoter_nuc_gt50` is NOT comparable across promoters of different length.** Insert
lengths are not uniform: 264 bp for most, RPS9 200, **minCMV 59**. At 59 bp the threshold is 30 bp
(nearly any grazing nucleosome clears it); at 264 bp it's 132 bp, about a whole nucleosome. Observed
on one plusDox sample: minCMV `gt50`=0.53 with mean NUC 30.0 bp, vs LYL1 `gt50`=0.074 with mean
31.5 bp — near-identical nucleosome coverage, 7× different boolean, purely from the denominator.
Use `--promoter_nuc_min_bp` for cross-promoter statements.

Synthetic suite 8/8 after both changes.

## Updates 2026-08-04 — `plus1_nuc` redefined to ABSOLUTE construct coords; negative-dyad note corrected

**The "+1 nucleosome" is construct-positioned, not TSS-positioned.** Measured on 498,302 NUC segments
(2 samples × 16 promoters × all copy-number variants). Using the nucleosome's interior (upper)
boundary — the MEASURED quantity, not the prior-regularized dyad, so this is not a padding artifact:

| frame | mean | sd | range |
|---|---|---|---|
| absolute column | 108.4 | **10.1** | 80–115 |
| relative to TSS | −59.9 | **23.6** | −121 → −31 |

2.3× tighter in absolute coordinates; 12/16 promoters have a median boundary of exactly 115, and
`end − prom_start` = −21 for nearly all (H4C5 is the outlier at −56). Peak dyad column is 46–51 for
EVERY promoter while `tss_col` ranges 136–236, so the apparent "+1 distance" (90 bp minCMV → 190 bp
BAX) was purely an artifact of promoter insert length. No separate peak at the canonical +120 for
high-`tss_col` promoters ⇒ **no distinct TSS-positioned +1 nucleosome in this data.** Mechanistically,
cols [0, prom_start) are the constant downstream backbone shared by all constructs; the nucleosome
fills it and is bounded ~21 bp below the insert start.

**Change:** `build_promoter_regions.py` now emits `plus1_nuc,0,{prom_start-1},NUC` — identical for
every promoter — instead of the TSS-relative +20..+170. **Name deliberately kept** for continuity with
the P110 JUNB analyses (`junb_cluster_states.py` etc. key on the literal string). Defined on segment
start/end, NOT on the dyad (see below).

**Impact (old vs new, same molecules):** the permissive `plus1_nuc_any_NUC` boolean is nearly
unchanged — 98.9% molecule-level agreement, overall rate 0.741 → 0.741 — so existing JUNB numbers are
safe (JUNB r=0.996, Δfrac +0.021). The *fraction* moves where the old window was displaced, and the
correlation degrades monotonically with `tss_col`: r = 1.000 at tss≈160, 0.915 LYL1 (189), 0.884 RPS9
(201), **0.621 BAX (236, Δfrac +0.235, +26 bp)**. Cross-promoter SD of the mean fraction drops
0.068 → 0.055, i.e. the new definition measures the same thing at every promoter.

**CORRECTION to the 2026-07-20 negative-dyad entry.** Negative dyads are implemented and real
(verified in code at the `anchor = 0.5*(start_bp+end_bp)` / `min(anchor, orig_amp_width)` site — only
the HIGH side is capped, and `anchor` is computed BEFORE the span clamp; segment start/end are still
clamped to [0, amp_width] since those index arrays). In one sample, 2.4% of NUC segments have a
negative dyad, min −35. **But they are CENSORED, not distributed**: quantiles of the negatives are
0%/25%/50%/75% = −35 exactly, max −2, i.e. ≥75% sit precisely on the `left_pad` floor. The stated
rationale for un-clamping — "different +1s with different interior extents get different, more/less
negative dyads, good for relative ordering" — largely does NOT hold in practice. Do not treat the
negative dyad tail as a usable continuous coordinate. This is why `plus1_nuc` is defined on
start/end overlap rather than on the dyad.
