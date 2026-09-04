# v5 HSMM calibration log

Running record of edge-case molecules, what the model got wrong, the diagnosis, the parameter
change, and its validation. Read [`EXPLAINER_hsmm_v5.md`](./EXPLAINER_hsmm_v5.md) for what the
parameters mean and [`NOTES_v5_hsmm_status.md`](./NOTES_v5_hsmm_status.md) for the full flag list.

**Purpose:** every default change to the model should be traceable to a real molecule + a measured
justification, not vibes. Add a new dated entry each time we tune something.

---

## The workflow (how to work an edge case)

1. **Find a bad molecule.** Note its `read_id`, sample, amplicon (from a v5 PDF, the bulk sanity
   panel, or the wide table).
2. **Inspect it** with `inspect_molecule_v5.py` — prints the raw GpC data, the decoded Viterbi
   path with a **per-segment score breakdown** (emission / start / duration / transition), and
   optionally scores a hypothesized "correct" segmentation so you can see the exact score delta and
   which term drives the wrong call:

   ```bash
   # run on a compute node (srun/sh_dev), from workflow/scripts/
   python inspect_molecule_v5.py \
     --input  <sample>.<amplicon>.dedup.full_unclustered.matrix \
     --positions <positions.long.txt> --amplicon_name <amplicon> \
     --read_id <read_id> --region 250,560 \
     --compare "TF@0;TF@1;TF@2;TF@3;TF@4;TF@5;TF@6"
   ```
   `--compare` accepts `TF@k` (TF at motif k) and/or `TYPE:start-end` (OPEN gaps auto-filled).
   Sweep a knob inline, e.g. add `--start_cost_tf 1.0`, and re-run to see if the call flips.

3. **Diagnose** which term is responsible (is the wrong path winning on emission, or on structural
   cost?). Emission wins => the p_U / conversion / sigmoid params. Structural cost wins => the
   start/transition/duration knobs.
4. **Pick the minimal knob** that flips the call. Prefer the one that maps to the real cause.
5. **Validate it doesn't break the bulk** — re-decode a few hundred molecules and check the change
   is targeted (recovers the intended calls without wholesale conversion of one state into another;
   watch mean n_nuc / n_tf / n_unid and per-motif occupancy).
6. **Record it here**, change the default in `classify_single_molecule_binding_v5_hsmm.py` (with a
   code comment pointing back to the dated entry), and update the param reference in
   `NOTES_v5_hsmm_status.md`.

---

## Edge cases

### 2026-07-10a — nucleosome swallows a run of bound TFs → lower `start_cost_tf` (2.0 → 1.0)

**Molecule:** `M00653:302:000000000-LC8HW:1:1109:15970:16317`, sample `L650`, amplicon
`opJS4_7x_TetO_21bp_no_CG`.

**Symptom (user):** "calls nucleosomes over many consecutive bound TFs that happen to have the
spots methylated between." The data over the array is textbook punctate binding — protected at each
TetO, accessible in every linker (…281P 302P **307a 313a** 321P 342P 345P **352a** 361P 382P…) — i.e.
all 7 TFs bound. The model (default params) called TFs 3–7 but fused **TetO1+TetO2 into one
nucleosome** (`NUC 235-350`, dyad ~292), swallowing the accessible linker at 307/313.

**Diagnosis (from `inspect_molecule_v5.py` breakdown):** there is genuine nucleosome-level
protection *immediately upstream* of the array (159–277), with **no accessible gap** before TetO1
(277P 281P contiguous). The nucleosome uses that upstream protection as a free anchor and extends
right to absorb TetO1+TetO2; the only cost is ~−9 for mis-predicting the two accessible linker GpCs
(307/313). Calling TF1+TF2 separately was *more* expensive only because it forced the upstream
protection to become its own extra footprint segment — i.e. the wrong call won on **start-cost
parsimony**, not on emission. TFs 3–7 were called correctly because they are isolated (accessible
linkers on both sides), so a nucleosome over them would mis-predict many linker GpCs.

**Fix:** `start_cost_tf` 2.0 → **1.0**. TF binding at these TetO arrays is common and shouldn't be
penalized when the linker evidence supports it. At 1.0 the molecule decodes correctly:
`NUC 150-278` (the real upstream nuc, now stopping before TetO1) → TF at all 7 motifs →
`NUC 545-680`.

**Validation (600-read subsample, L650 opJS4_7x):**

| start_cost_tf | mean n_tf | mean n_nuc | mean n_unid | per-motif TF occupancy |
|---|---|---|---|---|
| 2.0 (old) | 1.24 | 3.26 | 0.21 | 0.11 0.14 0.22 0.20 0.29 **0.10** 0.18 |
| **1.0 (new)** | **1.71** | **3.20** | **0.21** | 0.16 0.19 0.28 0.24 0.34 **0.29** 0.20 |

Key safety signals: mean **n_nuc barely moves** (3.26→3.20) and **n_unid is flat** — we are *not*
converting nucleosomes into TF calls wholesale, only flipping the swallowed cases. Motif 5, which
was anomalously low at 2.0 (0.10, systematically swallowed), rises to 0.29, in line with its
neighbors. If 1.0 were over-calling, n_nuc would collapse and occupancy would explode; neither
happens.

**Open follow-ups this raises:**
- Watch for the *opposite* failure now (spurious TF at a motif that's genuinely nucleosomal) as more
  molecules come in — if it appears, the more targeted lever is a stronger nucleosome-core
  accessibility penalty rather than raising `start_cost_tf` back up.
- The upstream-nuc/TetO1 boundary is intrinsically ambiguous (contiguous protection, no gap); the
  model now places the nuc edge at the motif boundary, which is reasonable but not data-determined
  there.

### 2026-07-10b — nuc-tail vs first-TF is a genuine coin-flip; do NOT hand-tune d_edge/trans_nuc_tf

**Molecule:** `M00653:302:000000000-LC8HW:1:2113:18379:6876`, `L650`, `opJS4_4x_TetO_21bp_no_CG`.
A contiguous ~143 bp protected block (159–302, no internal accessible gap) was split into
`NUC 150-278` + a spurious `TF 278-304` — the model bolted a TF onto the nucleosome's 3′ tail.

**Diagnosis:** the nuc+TF path beats the one-nucleosome reading by only **0.22 nats** — a near-tie,
not a confident error. Root cause: `d_edge=65` gives a ~130 bp footprint, too narrow to cover a full
~143 bp nucleosome from one dyad, so the tail (281/302) needs an extra segment.

**Why NOT to hand-fix:** cross-checked both candidate levers against the 2026-07-10a molecule (7x),
where TetO1 is a *real* bound TF that also abuts an upstream nucleosome. Both `trans_nuc_tf=2` and
`nuc_d_edge=73` **remove that real TF** — the two molecules are locally near-identical
(protected block ending at TetO1, right after a nucleosome), so no single hand-tuned cost separates
"nucleosome tail" from "first TF of an array." **Left defaults unchanged** (`d_edge=65`,
`trans_nuc_tf=1`); accepted as a coin-flip per the user. `d_edge=65` was always a placeholder to be
learned, not a confident value.

**0x nucleosome-length calibration + EM findings (same session):**
- The no-TF sample (`opJS4_0x`, 2344 molecules) is a clean nucleosome-length source. Model-free
  protection streaks: median **122 bp** GpC-to-GpC → true footprint **~140 bp**; real tail to
  200–400 bp (di-nucleosomes). So `nuc_mode ~140` is data-supported (default is 147; minor).
- Built Viterbi (hard) EM (`--do_em`). **Hard EM degenerates on the nucleosome shape/spread
  params** and must not fit them: on 0x it drove `nuc_d_edge` 65→103, `nuc_softness` 5→1.3 (a wide
  razor-sharp top-hat), inflated `prob_unmeth_given_open` 0.05→0.22, and collapsed `nuc_sigma`
  25→5. Cause: hard assignment gives the soft-edge transition zone entirely to OPEN (so the sigmoid
  fit never sees the edge), and decoded lengths cluster at the prior mode (so measured variance
  shrinks). This confirms the v4 "EM self-reinforces" warning.
- **Resolution:** `EM_FIT_SAFE` fits only `prob_unmeth_given_tf`, `prob_unmeth_given_unid`,
  `nuc_mode` by default. `d_edge`/`softness`/`p_open`/`nuc_sigma` are excluded (shown as `est_*`
  diagnostics in the EM log but not applied). Properly fitting the footprint shape needs
  **soft/forward-backward EM** (future) or external 0x-streak calibration. `--em_fit` can override
  the set if you know why.
