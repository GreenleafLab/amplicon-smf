# Notes: v4 nucleosome-model changes (Option A) — design + rollback

Design record for the 2026-07-08 discussion on v4 binding-classifier failure modes and the
chosen fix ("Option A"). Read alongside [`NOTES_v4_status.md`](./NOTES_v4_status.md) (overall v4
status) and [`NOTES_classify_binding_scripts.md`](./NOTES_classify_binding_scripts.md) (bug
catalog). **This file exists so the change is easy to understand and roll back later.**

## Motivating observations (three real molecules, opoBD9-style 6x-TetO amplicon)

- **Good molecule** `M00653:302:000000000-LC8HW:1:1102:21757:18698`: 3 nucleosomes + 4 TFs, with
  the unbound 3' TetO pair correctly absorbed into the downstream nucleosome. Core signal→call
  mapping works.
- **Bad molecule 1 (missed TF)** `...1119:19194:19277`: model tiled ~4 nucleosomes across the
  array and called NO TF, but there is punctate protected-at-motif / accessible-linker structure
  in the middle = a bound TF that got swallowed by nucleosome tiling.
- **Bad molecule 2 (extra TF)** `...1109:4975:13836`: one continuous ~160 bp protection block was
  split into nucleosome + an adjacent TF at its 3' edge, with no accessible linker between them —
  evidence better explained by a single longer nucleosome.

## Diagnosis (confirmed by code read)

- v4 nucleosomes are **fixed-footprint** via a distance-from-dyad logistic
  `p_U = expit((d_edge - dist)/softness)` (`build_nuc_protection_matrix`), NOT variable-length.
  This is correct/desirable (the sigmoid encodes "very bad to methylate at the dyad, fine on the
  flanks") — the nuc channel is *legitimately* different from the flat TF channel (`prob_unmeth_given_tf`).
- **Dead zone (drives molecule 2):** one nuc covers ~130 bp (d_edge=65 → half-max at ±65). Two
  nucs must be ≥140 bp apart (min spacing = `nuc_length // bin_size`) and then leave a protection
  dip between them + overshoot to ~270 bp. So a continuous **130–270 bp** interior streak can be
  neither one nuc (too short) nor two (dip+overshoot) → the model bridges the overhang with a
  spurious TF. Ghost-position extension already handles this at read *edges*, not the interior.
- **Sliding degeneracy (drives molecule 1):** soft flanks let a nucleosome be positioned so
  accessible linkers fall on the cheap low-`p_U` flank while protected motifs sit near the dyad,
  so an all-nucleosome tiling can out-score TFs+linkers.
- **Likelihood already weights each observed GpC equally** (`classify_all_molecules:316` is an
  unweighted sum) — there is no per-base upweighting to remove; that was a v2 idea, not in v4.
- **EM will not fix these:** `run_em_mstep` estimates params conditioned on the *current*
  assignments, so a nuc-biased assignment self-reinforces (`prob_unmeth_given_tf` left at default
  when `n_tf < min_obs`). No parsimony/identifiability logic in EM.

## Decision: Option A (chosen 2026-07-08)

Keep the **pure-dyad parametrization** (a nucleosome = a dyad position; "size" is entirely the
protection curve — there is NO length state). Chosen over Option B (per-nuc discrete length menu
{~147, ~180} as a new state dimension), which resolves the dead-zone tension more completely but
reintroduces length as state and multiplies the nuc state count. **Revisit B only if the residual
dead-zone TF-bridging is common after A.**

Known limitation accepted for A: the footprint (`d_edge`) is a **single global parameter**, so
"reach up to ~180 bp" applies to *all* nucs — it cannot be tight for typical ~147 bp nucs and wide
for the occasional 180–200 bp streak simultaneously. Generous `softness` mitigates (forgiving
flanks so a slightly-too-long footprint doesn't miss normal nucs) but does not fully resolve it.

### Concrete changes for A
1. **Decouple dyad spacing from footprint.** New config `nuc_min_spacing = 130` (was implicitly
   `nuc_length = 140`), used by `recursive_1_placer` + `check_nuc_state_sanity`. Footprint is
   `d_edge`/`softness` only.
2. **Softer, wider logistic (option a, not a length penalty).** Larger default `d_edge` (toward a
   ~180 bp span, i.e. d_edge ≈ 90) and larger `softness` so edges decay gradually and
   self-penalize. Let EM fit both; make EM bounds consistent (current cap `d_edge ≤ 90` already
   permits a 180 bp span; softness cap 30 already fine).
3. **Collapse observationally-equivalent states → keep most parsimonious.** EXACT `p_U` equality
   does NOT work (continuous sigmoid → floating-point-distinct vectors, collapses ~nothing), so
   collapse on a **quantized** signature (round per-GpC `p_U`, or per-GpC channel category
   nuc/TF/open + coarse protected/accessible bin). Tie-break: **fewest TFs first, then fewest
   nucs.** Quantization granularity exposed as a config knob; over-merge biases toward parsimony
   (acceptable — aligned with the goal).
4. **Move nuc params to config** (`nuc_min_spacing`, `d_edge`/`softness` defaults, collapse
   granularity), EM bounds kept consistent.

### Explicitly deferred (do NOT change now)
- Overlapping di-nucleosomes (octasome/hexasome): would require relaxing the non-overlap/spacing
  constraint; the 130 bp spacing intentionally still forbids it. A wide+soft footprint models ONE
  nuc's breathing/linker reach, NOT two overlapping cores.
- Per-nuc length menu (Option B).
- Parsimony *prior* on feature count (user lukewarm; the identifiability-collapse tie-break gives
  a parsimony effect without an explicit prior for now).

## Update 2026-07-08b: TF-boundary penalty (+ Option-A footprint reverted)

After the first Option-A run (dir `260708_v4updates_fromclaude`, sample L650 opJS4_6x), results:
- GOOD `1102`: 4 TF / 3 nuc preserved (fit worse, −13.0→−20.8, from the wide footprint).
- EXTRA-TF `1109`: **fixed** — spurious TF removed (1→0 TF), fit improved (−13.1→−8.9).
- MISSED-TF `1119`: 0→3 TF (right direction; possibly over-called — the accessible dips make 3
  punctate TFs genuinely higher-likelihood than a nucleosome).

So Option A worked but (a) the wide global footprint hurt normal-nucleosome fit and (b) nothing
stopped TFs being placed against nucleosomes. Added a **TF-boundary penalty** to encode the user's
likelihood intuition and **reverted the footprint** (`nuc_d_edge` 90→65, `nuc_softness` 8→5).

### The principle (user's, validated)
A real TF is a *punctate protected blob bounded by accessible GpCs on both OUTER sides of its
contiguous run*. One accessible flank is ambiguous with a nucleosome edge; the TF↔nucleosome
junction is the hard case ("to place a TF against a nuc AND get the intervening flank wrong, the
model must be damn sure"). The base likelihood cannot express this: it is site-independent Bernoulli
with `max`-combined channels, so an adjacent nucleosome *provides* the TF's junction protection for
free — there is no representation of runs/boundaries. The boundary term injects that structure.

### Implementation (`build_tf_boundary_penalty` + applied in `classify_all_molecules`)
- Per contiguous TF run, penalize PROTECTED observations in each **outer** flank zone, distance-
  ramped: 0 within `tf_flank_breath` (10 bp; tolerates TF breathing / a too-close flanking GpC that
  reads protected — the user's distance-dependence point), ramping to `tf_flank_lambda` (3.0) by
  `tf_flank_dist` (35 bp). Vectorized as `edgeL @ (wL @ sub) + edgeR @ (wR @ sub)`.
- **Run-length Occam cap**: runs longer than `max_tf_run` (3) pay `tf_run_lambda` (3.0) per extra
  site → 4–5 contiguous protected sites default to a nucleosome.
- All five params are config/CLI knobs. **λ and the window (`tf_flank_breath`/`dist`) WILL need
  tuning** — calibrate `tf_flank_lambda` via the ±TF nucleosome-length analysis (below).

### KNOWN SIMPLIFICATION to revisit (user clarification 2026-07-08)
TFs *can* bind directly adjacent to each other. Between adjacent TFs the intervening base is
**usually methylated (accessible) but not always** — a soft, probabilistic expectation, NOT a hard
rule. The current code **fully exempts inter-TF gaps** (only outer run edges are checked), which is
too permissive: it never expects the between-TF accessibility. Options: (a) add a *soft* (smaller-λ)
accessibility expectation on inter-TF gaps, or (b) leave it and let a future HMM handle it. This
"usually-but-not-always" is precisely the kind of thing a generative model expresses as an emission
probability rather than a hand-coded rule.

### Considering an HMM / semi-Markov rewrite (open, likely preferable long-term)
Hand-coding boundary/adjacency/parsimony as penalty terms is accreting hacks. A hidden semi-Markov
segmentation (states: open / nucleosome / TF-bound; explicit durations; emissions = the existing
sigmoid-nuc / flat-TF / low-open protection; transitions encode adjacency + parsimony) would express
ALL of this natively and replace enumeration + pruning + collapse + boundary-penalty + run-cap with
one principled model decoded by Viterbi (linear per molecule, no state explosion). Estimated effort:
~400–700 lines + modeling/tuning; reuses the emission model and I/O. It is a genuine rewrite, not a
tweak. Decision pending: ship the boundary-term patch for now vs. invest in the HMM.

## Rollback plan — commit hashes and exact commands (2026-07-08)

All work is on branch **`v4-nuc-model-updates`** (off master `9dd6970`). Commits so far:

| Hash | What it is |
|---|---|
| `9dd6970` | master tip when we branched (pre-everything) |
| `6b400fb` | **CpG feature** — per-molecule clean-CpG side output (validated, keep) |
| `cfc8f67` | **v4 baseline** — classifier exactly as-is BEFORE Option A. ← **ROLLBACK ANCHOR** |
| *(pending)* | Option A nucleosome-model changes (committed only after test-run looks good) |

### How to undo (copy/paste; run from the `amplicon-smf/` directory)

**See where you are:** `git log --oneline -5` and `git status`

**Undo the Option-A edits BEFORE they're committed** (just discard working-tree changes to the
one file, restoring the committed baseline):
```
git checkout -- workflow/scripts/classify_single_molecule_binding_v4.py
```

**Undo Option A AFTER it's committed** (throw away the Option-A commit, keep CpG + baseline):
```
git reset --hard cfc8f67
```
`cfc8f67` is the v4 baseline; the CpG commit `6b400fb` is its parent, so it survives this reset.

**Keep the CpG feature on master later** (when ready — this is a push-free local merge):
```
git checkout master
git merge v4-nuc-model-updates      # or: git cherry-pick 6b400fb   to take just CpG
```

**Throw the whole experiment away** (back to clean master, lose the branch):
```
git checkout master
git branch -D v4-nuc-model-updates
```
(Nothing has been pushed to any remote, so all of this is local and safe.)
