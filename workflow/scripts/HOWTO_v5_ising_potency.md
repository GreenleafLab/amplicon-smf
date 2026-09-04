# HOWTO: run the v5 → Ising → potency chain on new samples

End-to-end recipe for taking raw amplicon-SMF data to (a) per-molecule footprint calls,
(b) equilibrium array parameters `h, J, mu, delta`, and (c) potency `k_pot` with its
instantaneous/integrated split `w`.

Validated 2026-09-04 by reproducing **Nature 2024 Fig. 3** from the original 2022 opJS45 FASTQs.
That reference run — including its samplesheet, config, sbatch scripts and a `NOTES.md` recording
every decision — lives at
`data/260903_v5_hsmm_partition_function_test_old_data/` in the parent project.
**Copy that directory's layout when starting a new one.** What it reproduced:

| quantity | published | this chain |
|---|---|---|
| potency `k_pot` | ~0.090 | 0.085–0.103 |
| basal `p0` | ~0.03 | 0.027–0.035 |
| nucleosome-displacement term | "the important one" | `delta` = −1.09, biggest single AIC gain |

New: `w` ≈ 0.33–0.54 (fraction of potency that is *integrated* rather than instantaneous),
consistent across two independent axes.

---

## Step 0 — set up a contained run directory

```bash
RUN=/path/to/data/YYMMDD_my_experiment
mkdir -p $RUN/logs && cd $RUN
ln -s /path/to/amplicon-smf amplicon-smf     # rules invoke scripts by RELATIVE path
```

The `amplicon-smf` symlink is required: `workflow/rules/other.smk` calls
`python amplicon-smf/workflow/scripts/...`, resolved against your working directory.

Write `<name>.tsv` (samplesheet) and `<name>.yaml` (config). Copy both from the reference run and
edit. Also write a **sample metadata table** (`series`, `dox_ng_ml`, `biorep`, …) and group
downstream analyses on *that*, never on parsed sample names.

> **Trap:** when pulling rows from an older project directory, use its **top-level
> `<experiment>.tsv`**, not `amplicon-smf/config/samples.tsv` — the latter are often stale
> single-row leftovers that silently process one sample and look like they worked.

## Step 1 — core pipeline (FASTQ → matrices)

```bash
snakemake -s amplicon-smf/workflow/Snakefile --configfile <name>.yaml \
  -c 8 -k -w 15 --use-conda --conda-frontend conda \
  --conda-prefix <an existing run>/.snakemake/conda --rerun-incomplete
```

- `--conda-frontend conda` is required (no `mamba` on Sherlock). Reusing a `--conda-prefix` whose
  `envs/python3_v6.yaml` is byte-identical skips a ~2.3 GB rebuild.
- Set `threads: 1` in the yaml and get parallelism from `-c`. In `other.smk` the aligner's thread
  count is a `params`, not Snakemake's `threads:`, so Snakemake cannot schedule around it and
  `threads: N` with `-c N` oversubscribes.
- Snakemake exits **non-zero** if any job fails under `-k`, so Slurm will report FAILED even when
  the run is fine. Check *which* rules failed before worrying.

**Then check the QC before going further:**

```bash
cat results/*/*/stats/*.background_cpg_methylation.sample_level.txt
cat results/*/*/stats/*.duplication_rate.stats.txt
```

`no_endog_meth` **must** be decided from this, not assumed:

| `cpg_methylation_pct` | meaning | set |
|---|---|---|
| low (few %) | endogenous CpG methylation negligible | `no_endog_meth=TRUE` (keeps GCG-ambiguous GpCs) |
| high (>~10%) | GCG sites are contaminated | `no_endog_meth=FALSE` (drops them) |

In the 2022 opJS45 data this was **14–17%**, so `FALSE` was correct — even though `TRUE` matches
what v5 was calibrated on. If you need both, run the pipeline twice into two `experiment` names and
use the `FALSE` tree for anything promoter-flavoured. For opJS4 the only ambiguous positions are
GpC **54 and 130**; 130 is inside the TSS window, and **neither is in the TetO array (282–581)**,
which is why the Ising/array fit was unaffected while promoter/potency was not. **Re-derive those
positions for a new amplicon** — count `GCG` in the post-revcomp reference.

## Step 2 — v5 HSMM decode

Copy `run_v5_decode.sbatch` from the reference run (array over samples, `$EXPERIMENT` selects the
tree, resumable, ~53 s per 9.3k reads).

```bash
sbatch --export=ALL,EXPERIMENT=<experiment> run_v5_decode.sbatch
```

- **Positions file must use bare-operator coordinates**
  (`convert_fa_to_positions_for_script.py --l_offset 0 --r_offset 0`). The old `+2/-1`-asymmetric
  opJS45 file mis-calls TetOs.
- **Decode the 0x / zero-motif amplicons too.** They carry no TF signal but are the cleanest
  constraint on nucleosome fugacity and they supply `--mu_ambient` (Step 3). v5 handles K=0 and
  simply emits no `tfbs_`/`site_` columns.
- Chain dependent jobs with `--dependency=afterany`, **not `afterok`** (see the non-zero-exit note
  above).
- Expect fewer outputs than samples × amplicons if any samples are single-library.

Output schema (differs from `260715_v5_downstream_design.md`, which predicted `seg_loglik`):

```
main:     read_id, tfbs_1..K, site_<name>1..K, n_tf, n_tf_teto, n_tf_other,
          n_nuc, n_unid, nucs, unids, log_likelihood
segments: read_id, seg_index, type, start, end, dyad_or_motif, motif_name
sidecar:  <sample>.<amplicon>.tfbs_index.txt   index -> name -> lo/hi -> is_teto
```

`tfbs_{k}` is ordered by ascending operator coordinate (the sidecar proves it), which is what
`n_pairs` depends on. Segment coordinates share the positions-file frame. Use **`n_tf_teto`**, not
`n_tf`, when off-target motifs are also named.

## Step 3 — Ising / partition-function fit

```bash
python amplicon-smf/workflow/scripts/fit_ising_model.py selftest      # 49 assertions, run first

# measure the ambient nucleosome coverage from the ZERO-OPERATOR amplicons,
# excluding the promoter, then convert it to a fugacity:
python .../fit_ising_model.py mu_ambient --coverage 0.7826 --nuc_len 147   # -> -0.1784

python .../fit_ising_model.py fit \
  --main_files binding_v5/<exp>/<SAMPLE>/<SAMPLE>.<amp glob>.single_molecule_classification.txt \
  --matrices   results/<exp>/<SAMPLE>/matrices/<SAMPLE>.<amp glob>.dedup.full_unclustered.matrix \
  --positions  <bare positions file> \
  --mu_ambient -0.1784 --mu_ambient_scan 0.75 0.82 5 \
  --specs 2param 3param_tfcoop 3param_nuc 4param --out_prefix YYMMDD_ising
```

`mu_ambient` is a **log-fugacity in kT, not a log-probability** — get it from measured coverage and
invert; never guess. `python fit_ising_model.py mu_ambient` with no `--coverage` prints the table.
Re-measure it for every new construct.

**Read the emitted `.ising_model_checks.pdf`, and trust it over AIC.** Each panel is labelled
FIT MOMENT (near-tautological) vs PREDICTION (evidence). We hit a case where AIC preferred a
variant (`3param_nucadj`) whose out-of-sample predictions were ~4× worse, because the likelihood
only ever sees the mean moments.

Known limits, all diagnosed and none blocking:

- `J` is **not robustly determined** — it ranged +1.06 → +0.02 → −0.28 depending on which
  nucleosome coupling was in the model. Do not quote it as cooperativity without the spec.
- `<n_nuc>`'s range is **compressed** vs data, because `delta`'s switch form
  `1[n_tf>=1]` saturates as soon as P(n_tf≥1) does. Candidate fix: the graded dosage form
  `delta*n_nuc*n_tf` (needs joint rod-count × TF-count DP state, ~60× bigger). Not built.
- Operators are idealized as **single bins at the motif midpoint** (±10 bp slop vs a 20 bp TetO).
- `--tie_mu_ambient` exists but its **gradient is inconsistent and it does not converge**. You
  should not need it: including the 0x amplicon already pulls fitted `mu` to within 0.013 kT of
  `mu_ambient`.

## Step 4 — potency and the two-timescale split

```bash
python .../fit_two_timescale_potency.py \
  --binding_dir binding_v5/<the no_endog=FALSE tree> \
  --samples <S1> <S2> ... \
  --amplicon_glob 'opJS4_?x_TetO_21bp_no_CG' \
  --promoter_lo 136 --promoter_hi 195 \
  --out_prefix YYMMDD_potency
```

Activity is v5-native: active = **no NUC segment overlapping the promoter window**. `(136,195)` is
the opJS4 minimal-promoter window from `assign_promoter_state_from_model_new.py` — **you must
change it for a different promoter.** Sanity-check it against a known bulk number before trusting
anything downstream (it gave 0.295 at 6x vs ~0.25–0.30 published).

Outputs `k_pot`, `w`, `tau_ratio` = `w/(1-w)` = τ_int/τ_c, plus the two-panel PDF (bulk Fig. 3d on
the left, panel-h coordinates on the right).

### The two axes, and why you need both

```bash
# axis 1 -- copy number: many amplicons, one sample
--samples <one sample> --amplicon_glob 'opJS4_?x_TetO_21bp_no_CG'

# axis 2 -- dox at FIXED geometry: one amplicon, many dox levels, pooled
--samples <all dox levels> --amplicon_glob 'opJS4_6x_TetO_21bp_no_CG' --pool_samples
```

They carry **different confounds** — copy number changes which sites are bound; dox changes global
TF level and so possibly global cell state. Agreement between them is the actual evidence for
`w > 0`; neither alone is decisive. In the reference run: 0.34–0.54 (copy number) vs 0.33–0.44
(dox), with `k_slow` significant at fixed geometry (z = 2.7 and 4.1).

**Check `k_interact_z` (parallelism).** |z| ≳ 3 means the panel-h lines fan out and the
two-timescale form is wrong. Do **not** use `identity_residual` as a test — `k_fast + k_slow` equals
the bulk slope by algebra, on any data.

---

## Porting to a new TF or promoter — the checklist

1. **New positions file**, bare coordinates, with the new motif's windows. Name TetO-equivalent
   sites so `--teto_name_prefix` / `n_tf_teto` picks them up.
2. **Re-derive the GCG-ambiguous positions** for the new reference and re-decide `no_endog_meth`
   from that sample's own background-CpG QC.
3. **Re-measure `mu_ambient`** from a zero-motif amplicon of the *new* construct. Do not reuse
   −0.1784.
4. **Change the promoter window** `--promoter_lo/--promoter_hi` and validate the resulting
   activity against an independent bulk measurement.
5. **Include a zero-motif amplicon** in the panel if at all possible. It pins `mu`, gives
   `mu_ambient`, and doubles as the no-binding control.
6. **Get a genuine no-TF control** (no TF expressed at all) rather than relying on zero-inducer:
   inducible TFs are leaky, so "0 dox" signal is a mixture of real leaky binding and decoder false
   positives and only bounds the latter. In the reference data that combined signal was ~6.6% per
   site, rising with copy number — the rising pattern is more consistent with cooperative leaky
   binding than with a flat error rate.
7. **Expect batch-to-batch variation in the absolute numbers.** In the published work this was
   handled by making claims only against same-day controls. Design comparisons that way.
8. **Run both potency axes** and require agreement before believing `w`.
