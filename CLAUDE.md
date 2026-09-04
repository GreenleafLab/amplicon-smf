# CLAUDE.md — amplicon-smf pipeline

Context for the Snakemake pipeline in this directory. This is the FASTQ → single-molecule
matrix → binding/promoter/potency processing code referenced by the parent project's
`/oak/stanford/groups/wjg/bgrd/papers/ad_smf/CLAUDE.md`. That file covers the paper/analysis
side (notebooks, `data/`, figures); this one covers how the raw data actually gets processed.

---

## What this pipeline does, in one paragraph

Amplicon SMF measures chromatin accessibility by treating cells/nuclei with a GpC
methyltransferase (M.CviPI) — or, for some samples, a cytosine deaminase — before bisulfite-type
conversion. Accessible GpCs get methylated (deaminase: converted) and read out as "protected"
(unconverted) after conversion+sequencing; DNA occupied by protein (nucleosome or TF) stays
unmethylated and reads as "accessible" (converted) — i.e. **the protection signal is inverted
relative to intuition**: methylated/unconverted = accessible-to-enzyme = NOT protein-bound;
converted = inaccessible = protein-bound. Reads are amplicon-targeted (PCR against a known small
set of loci, not genome-wide). The pipeline aligns reads with a bisulfite-aware aligner
(bwa-meth), scores conversion state at GpC (and optionally CpG/all-C) positions per read, and
assembles a single-molecule matrix per amplicon: rows = molecules, columns = genomic positions,
values `1` = protected/bound, `0` = accessible, `-1` = no info. Downstream scripts turn that
matrix into TF/nucleosome occupancy calls, promoter chromatin states, a thermodynamic
partition-function fit, and a scalar "potency" per sample.

## Two halves of this repo

1. **Core Snakemake pipeline** (`workflow/Snakefile`, `workflow/rules/*.smk`) — FASTQ in,
   single-molecule matrix + bulk QC plots out. Fully wired, runs end-to-end via `snakemake`.
2. **Downstream analysis scripts** (`workflow/scripts/`) — binding classification →
   partition-function fit → promoter-state calling → potency → consolidation. Runnable
   standalone (full `argparse` CLIs), but **not yet wired into any Snakemake rule** — you invoke
   them by hand per the README's "Downstream analyses" section. Several are recent/untracked
   (`git status` shows them as new files) and are actively being rewritten/consolidated — see
   [Downstream scripts: current state](#downstream-scripts-current-state) below.

---

## Core pipeline: rule graph (`workflow/rules/other.smk`)

```
fastq_R1/R2  ─┐
              ├─ reverse_complement_fastq (×2)          [swaps R1/R2! see gotcha below]
amplicon.fa  ─┴─ reverse_complement_fasta                [RCs the amplicon if bottom_strand]
                       │
                       ├─ index_fasta (bwameth.py index)
                       │
              align_bwameth  or  align_bwameth_all       [choice driven by filter_contigs]
                       │
        correct_mismatched_amplicons (filter_bam_by_matching_contigs2.py)
                       │
                  sam_to_bam → sort_index_bam
                       │
              filter_uncoverted (mark-nonconverted-reads-and-plot.py)
                       │
        ┌──────────────┴───────────────┐
   run_methyldackel               join_reads_and_first_cluster
        │                          (dSMF_footprints_clustering_py3.py)
   plot_bulk_methylation           │
   (from bedgraph)          ├─ matrices/*.full_unclustered.matrix (+ .dedup., + .clustered.)
                             ├─ amplicon_stats.txt
                             ├─ plot_bulk_methylation2 (from matrices)
                             └─ plot_nuc_qc (nucleosome length GMM)
```

Everything lands under `results/{experiment}/{sample}/...`; final plots under
`results/{experiment}/plots/`.

### Key gotchas in the core pipeline (worth knowing before debugging)

- **R1/R2 are swapped throughout.** `get_fastq()` in `other.smk:32` maps wildcard `read` to
  `fastq_R{3-read}` — i.e. Snakemake's "read1" is actually the samplesheet's `fastq_R2` and vice
  versa. The README's `read1_length`/`read2_length` docs call this out explicitly: "read1_length
  should actually be length of read2 from the sequencer." If you're chasing a length-threshold or
  wrong-primer bug, check you have this backwards before anything else.
- **Amplicon FASTA orientation matters and is non-obvious.** `bottom_strand` (default `TRUE`)
  reverse-complements the amplicon FASTA before indexing. README: amplicon FASTA should be on
  the *read1-primer strand*, i.e. typically the "top" strand on Benchling — but because of how
  primers are designed, only one actual strand is captured, hence the RC step. If plots look
  reversed, this is why (see README FAQ "Why is everything backwards?").
- **`filter_contigs` picks between two bwameth variants.** `bwameth.py` (best-alignment-only) vs.
  `bwameth_all_alignments.py` (adds `-a` to emit all alignments, needed so
  `filter_bam_by_matching_contigs2.py` can find pairs where R1/R2 map to the *same* amplicon
  contig — necessary because the amplicon FASTA can contain many similar/overlapping constructs).
  Should basically always be `TRUE` per the samplesheet doc ("necessary for backwards
  compatibility reasons").
- **Conversion filtering (`mark-nonconverted-reads-and-plot.py`) is skipped for `deaminase`
  samples** (`c_frac` forced to `0.001` in `other.smk:279`) since deaminase readout doesn't have
  a separate non-GpC-C conversion-rate QC the same way M.CviPI does.
- **`c_type` for matrix-building is derived, not a direct samplesheet column**
  (`other.smk:332-339`, `get_c_type`): `deaminase→allC`, else `include_cpg→both_dimers`, else
  `GC`. Endogenous CpG-overlapping GpCs (GCG trinucleotides) are excluded by default unless
  `no_endog_meth`/`-noEndogenousMethylation` is set.
- **`dedup_on` (recent addition, see `git log`)** controls whether duplicate-molecule detection
  in `dSMF_footprints_clustering_py3.py` uses only the scored context (`c_type`, default) or all
  Cs (`allC`). All-Cs dedup is more permissive (more reads called duplicates) and per the commit
  message may be "too lenient" due to sequencing error — this is an open judgment call per
  experiment, not a settled default.
- **Config/sample schema validation is effectively disabled.** `workflow/schemas/config.schema.yaml`
  and `samples.schema.yaml` are still the cookiecutter-template stubs (`samples`/`condition`
  columns) — they don't describe the real samplesheet columns at all, and the `validate(...)`
  calls in `common.smk` are commented out. Don't trust the schema files as documentation of valid
  samplesheet columns; the README's "Input Files" section is the actual source of truth.
- **`bwameth.py` is vendored, unmodified** (its own header says "THIS VERSION HAS NO EDITS",
  standard bwa-meth v0.2.9). `bwameth_all_alignments.py` is the one intentional fork, adding `-a`
  only (header: "THIS VERSION HAS THE -A EDIT!"). There are also `*_OLDIGNORE.py` variants in
  `workflow/scripts/` that are dead/superseded — don't confuse them with the live versions.
- **Off-by-one in matrix coordinates is a known, deliberately-deferred bug.** Per the README's
  "To Do": the first base in the reference is chopped off in matrix files, likely from a 1- vs
  0-indexing mismatch in `convert_amplicon_fa_to_peaklist.py`. Left alone because fixing it would
  ripple into the binding-model code that has patched around the offset.

### Matrix/data conventions

- Single-molecule matrix values: `1` = protected (bound/nucleosome), `0` = accessible
  (methylated/converted), `-1` = no information at that position for that molecule.
- Matrix columns for `-` strand regions are reversed to canonical orientation
  (`dSMF_footprints_clustering_py3.py`).
- `common.load_single_molecule_matrix()` treats `-1` as NaN and drops uninformative columns.
- Three matrix flavors per sample/amplicon: `.full_unclustered.matrix` (all reads),
  `.dedup.full_unclustered.matrix` (unique molecules only — this is what almost everything
  downstream actually consumes), `.clustered.matrix` (k-means-clustered, only if `-cluster`
  passed).

### Conda envs

`workflow/rules/envs/` has many versioned python envs (`python3_v2` through `_v6`,
`smf_py3_v7.yaml`); per the README, the pipeline has converged on a single py3 env
(`smf_py3_v7.yaml` looks to be the current one — check which one `*.smk` rules actually
reference if you need to add a dependency, most currently point at `python3_v6.yaml`).

### Tests

`tests/` has unit tests (`unittest`, run via `tests/Snakefile` →
`python -m unittest discover`) for the alignment-scoring internals of
`dSMF_footprints_clustering_py3.py`: `test_score_reads.py` (per-read/per-region C-context
scoring), `test_merge_alignments.py` (multi-segment alignment reconstruction),
`test_reconstruct_alignment.py`. No tests exist yet for the downstream binding/promoter/potency
scripts.

---

## GpC / CpG / GCG ambiguity handling (investigated 2026-07-01)

The GpC MTase methylates accessible GpC dinucleotides (stays unconverted → matrix `0`); protein
occupancy blocks methylation (converts → matrix `1`). A **GCG trinucleotide** has a C that is
simultaneously the C of a GpC *and* the C of a CpG — if that C reads unconverted, you can't tell
whether that's real GpC-driven accessibility signal or incidental/background CpG methylation
unrelated to chromatin state. Only the "unconverted" reading is confusable this way; a converted
call at a GCG is unambiguous (no enzyme methylated it, full stop).

### Two independent mechanisms — different robustness properties

1. **Automatic, always-recomputed exclusion** — `score_read_against_region()` in
   `dSMF_footprints_clustering_py3.py:225-236`. For `c_type=GC` (or the GC half of
   `both_dimers`), a GpC at `(i, i+1)` is excluded from the matrix (scored `-1`) if
   `ref_arr[i+2]=='G'` (i.e. it's a GCG) — *unless* the samplesheet sets `no_endog_meth=TRUE`,
   which flows through `other.smk:353` to the CLI flag `-noEndogenousMethylation` and disables the
   exclusion. **This is robust and strand-correct by construction**: the check and the position it
   acts on (`c_pos = i+1`) are computed from the same `i`, in the same function call, fresh from
   whatever reference file is actually being used that run. There's no hardcoded number to go
   stale here, and it automatically adapts to any amplicon design.
2. **Manual, hardcoded imputation** — `adjust_gcgs(mat, gcg_pos, fix_pos)` in `common.py:316-325`.
   For molecules that read `0` (ambiguous) at a *specific, pre-identified* ambiguous GCG column,
   overwrites that value with a *specific, pre-identified* nearby trustworthy GpC column's value.
   Used with hardcoded position pairs in `assign_promoter_state_from_model_new.py:469-470`
   (`(54,50)` and `(130,123)`) and in the older/non-canonical `classify_single_molecule_binding_v2.py`/`_v3.py`/`_v4.py`
   (`(53,49)` and `(129,122)` — commented out in the canonical `_v2_carosversion.py`). **This is
   inherently fragile** — the numbers are hand-derived by inspecting a specific matrix/amplicon at
   a point in time and don't get recomputed. This is intentional and should stay that way: whether
   a specific nearby GpC is "close enough" to trust as a stand-in for an ambiguous site is a
   judgment call (proximity, whether it's on the correct side of a footprint boundary, etc.) that a
   string-scan can't make — it should remain a deliberate, documented, opt-in correction per
   amplicon, not something derived/applied automatically. Don't automate this away.

### Strand direction (`bottom_strand`) — why checking literal "GCG" post-RC is correct

GC and CG are each self-reverse-complementary as 2-mers (revcomp("GC")="GC", revcomp("CG")="CG"),
so plain GpC/CpG identification is strand-symmetric — a GpC on one strand is "the same" dinucleotide
on the other strand, just at a mirrored coordinate. **GCG is not** (revcomp("GCG")="CGC"), so the
ambiguity check is inherently strand-specific and this needed to be verified carefully rather than
assumed.

Resolution (confirmed by the user, whose lab designs these libraries): after bisulfite/enzymatic
conversion, the top and bottom strands of the amplicon are no longer reverse complements of each
other, so primers are designed to selectively amplify one specific strand — for these experiments,
**the bottom strand**. The `bottom_strand` samplesheet flag (default `TRUE`) drives the
`reverse_complement_fasta` rule, which reverse-complements the user-supplied (top-strand-oriented)
amplicon FASTA before it's used as the alignment/scoring reference. Since
RC(top strand) = the bottom strand written in its own 5'→3' frame, **the post-RC reference file
already directly *is* the bottom strand's sequence** — so checking literal "GCG" (not "CGC") in
that file, which is exactly what `score_read_against_region` does, is correct. No additional
strand translation is needed; the one RC step already IS the translation.

### Verification methodology (redo this if the amplicon design or matrix format changes)

Checked against the actual production amplicon file —
`/oak/stanford/groups/wjg/mhinks/projects/smf/221106_P035_rTetR_only_opJS45/amplicon-smf/amplicon-info/opJS45.amplicon.long.fa`
(**not** `.../220829_P026_opJS45/amplicon-info/opJS4.amplicon.fa`, a similarly-named file from a
different project directory with genuinely different sequence — this is a real footgun given how
many near-identically-named amplicon files exist across project directories under
`/oak/stanford/groups/wjg/*/projects/smf/`; always confirm which one a given analysis actually
used before trusting a hardcoded default path like `aggregate_binding_model.py`'s `--pos_dicts`).

For `opJS4_6x_TetO_21bp_no_CG`, reverse-complemented exactly as `bottom_strand=TRUE` would:
`rc[53:56]=='GCG'`, i.e. the ambiguous C is at index **54** (the middle base — a regex/string
match's start index is the *first* base of "GCG", the G, not the scored C; this tripped me up
once mid-investigation, worth double-checking anyone's arithmetic here including your own).
Confirmed `rc[54]=='C'`. The neighboring imputation targets are also correctly C-indexed:
`rc[50]=='C'` (with `rc[49]=='G'`) and `rc[123]=='C'` (with `rc[122]=='G'`). So
`adjust_gcgs(mat, 54, 50)` and `adjust_gcgs(mat, 130, 123)` in
`assign_promoter_state_from_model_new.py` are **correct as written** for this amplicon.

Caveats on generalizing this:
- Only verified for the `opJS4_6x_TetO` variant. It should hold for the other `opJS4_Nx_TetO`
  copy-number variants too, structurally: the TATA/promoter/reporter-start block is a
  constant-sequence region that sits *after* the variable-length TetO array in the original
  (pre-RC) orientation, so after RC it lands at a fixed low-column position regardless of TetO
  copy number — consistent with one hardcoded pair being applied uniformly across all of them in
  `assign_promoter_state_from_model_new.py`.
- Does **not** extend to structurally different amplicons in the same fasta (e.g. `CTCF`, `BD24`)
  — if those are ever passed in `--amplicons` to the same script, the same hardcoded correction
  would be applied at a meaningless position. No check currently guards against this.
- The older/non-canonical classifiers' hardcoded pairs (`53,49` / `129,122`) are the G-indexed
  (off-by-one) version of the same positions — consistent with the "off by 1" comment in the new
  script being a deliberate, correct fix when porting to the current one-column-per-GpC format,
  not a leftover bug. Only matters if v2/v3/v4 are ever run with default
  `--convert_ambiguous_gcgs`/`--fix_ambiguous_gcg` against current-format matrices; low priority
  since the canonical classifier doesn't call `adjust_gcgs` at all.

### Net assessment

Default GCG exclusion (governed by `c_type`/`no_endog_meth`) is strand-correct and self-maintaining
— safe to trust without re-verification when the amplicon design changes. The hand-curated
`adjust_gcgs` imputation is correct for the one amplicon/position pair checked here, but by design
stays a manual, opt-in, documented judgment call rather than something to auto-derive — re-verify
by hand (using the method above) whenever it's applied to a new amplicon.

### Background CpG methylation QC (added 2026-07-01)

`plot_bulk_methylation_signal.py` now always runs a second, unconditional check (independent of
`include_cpg`/`no_endog_meth`/`deaminase`) via `plot_background_methylation_qc()`, so background/
endogenous CpG methylation can be empirically confirmed negligible (or not) regardless of which
enzyme mode a sample was run in, rather than assumed. Per-amplicon it computes, from the same
`_CHG`/`_CHH`/`_CpG` bedgraphs `run_methyldackel` already produces for every sample:
- **"clean" CpG-context methylation**: C followed by G, *excluding* ambiguous GCG sites (C also
  preceded by G) — deliberately excluded so this reads as background methylation, not GpC MTase
  activity leaking into the CpG channel (see the ambiguity discussion above).
- **"other C" methylation**: not preceded or followed by G — same definition as
  `filter_snps()` in `mark-nonconverted-reads-and-plot.py`, so this number is directly comparable
  to that script's `c_frac` conversion-efficiency threshold.

Outputs (wired into the existing `plot_bulk_methylation` rule in `other.smk`, no new rule needed
since inputs are identical): a trace PDF (`{sample}.background_cpg_methylation.pdf`, own
`PdfPages`, two pages per amplicon — clean-CpG trace and other-C trace), a per-sample summary plot
across amplicons (`{sample}.background_cpg_methylation.summary.pdf`, grouped bars per amplicon),
and two stats tables: `{sample}.background_cpg_methylation.amplicon_level.txt` (one row per
amplicon) and `.sample_level.txt` (one row — **mean of the per-amplicon means**, i.e. every
amplicon weighted equally regardless of how many qualifying positions it has; reconsider this
weighting if amplicons in a panel end up with very different CpG-site counts). The sample-level
number is the one meant to get pulled into a future cross-experiment tracking table, following the
same `_amplicon_level`/`_sample_level` naming convention as `aggregate_binding_model.py` and
`assign_promoter_state_from_model_new.py` so it can plug into that aggregation machinery later
without renaming anything.

### Per-molecule clean-CpG methylation side output (added 2026-07-08)

`dSMF_footprints_clustering_py3.py` gained an optional `--cpg_meth_stats <path>` that writes a
per-read **GpC-uncontaminated** CpG-methylation table (columns: `amplicon, read_id,
n_clean_cpg_covered, n_clean_cpg_methylated, frac_cpg_meth`) — intended as a per-molecule
endogenous-CpG-methylation QC/filter (drop molecules above some `frac_cpg_meth`). "Clean CpG" uses
the same definition as the bulk background QC above (`clean_cpg_positions()`: C followed by G, not
preceded by G → excludes GCG), verified to match `cpgs - gpcs` across all 160 amplicons of a real
panel. Methylated = unconverted cytosine, same convention as the accessibility matrix. Rows are
emitted for **all retained (pre-dedup) reads** — a superset of the dedup matrix's read IDs, so it
joins by `read_id` to either the `.full_unclustered` or `.dedup.full_unclustered` matrix.

**Side output only — does not touch any existing matrix/output.** Independent of `--include_cpg`
(which instead folds CpGs *into* the accessibility matrix via `both_dimers`). Gated to M.CviPI
samples: the file is always created (fixed Snakemake output) but rows are only computed when
`c_type != 'allC'` (deaminase ⇒ `allC`, where the stat isn't meaningful). Wired into the
`join_reads_and_first_cluster` rule (`other.smk`) as output `cpg_per_read` and listed in
`common.smk` `wanted_input`. Note this is per-molecule (single-read granularity), distinct from the
bulk per-amplicon background-CpG QC above.

### Cross-experiment QC stats tracking (added 2026-07-01)

Distinct from `consolidate_data.py`/`aggregate_binding_model.py` (which build the paper's
scientific tables for a specific set of samples you already know about), there's now a
lightweight, generic mechanism for pulling simple QC numbers (duplication rate, mapped/unmapped
counts, background CpG methylation, nucleosome length QC, ...) across *every* amplicon-smf project
you've ever run, even ones living in entirely separate directories under
`/oak/stanford/groups/wjg/*/projects/smf/*`.

**Convention**: every "simple" per-sample QC number goes in `results/{experiment}/{sample}/stats/`
as a headerless, two-column `metric\tvalue` TSV (one file per check, one or more metric rows per
file) - e.g. `{sample}.bwameth.contig_filtered.stats.txt`, `{sample}.nuc_len_qc.stats.txt`,
`{sample}.duplication_rate.stats.txt`, `{sample}.background_cpg_methylation.sample_level.txt`.
Per-amplicon/per-construct breakdowns that you want to plot directly stay **wide and headered**
instead (e.g. `{sample}.amplicon_stats.txt`, `{sample}.background_cpg_methylation.amplicon_level.txt`)
and live wherever's convenient - `amplicon_stats.txt` deliberately stays at the sample's top level,
not in `stats/`, with `compute_duplication_rate.py`/`{sample}.duplication_rate.stats.txt` as its
tidy per-sample summary sibling in `stats/`.

**Collector**: `workflow/scripts/collect_sample_stats.py` - not a Snakemake rule, a standalone
tool since it deliberately spans multiple independent project directories, each with its own
`results/` tree. Takes `--dirs <dir1> <dir2> ...`, walks each `results/*/*/stats/*.txt`, and
auto-detects which files are tidy (`is_tidy_metric_value_file`: every line has exactly 2 tab
fields, second field parses as a number) vs. wide/headered (skipped, since a header row's second
field is text, not a number, and/or there are more than 2 columns). No hardcoded list of check
names or metrics - any future step that writes a tidy metric/value file into `stats/` is picked
up automatically next time this runs. Output is one long table:
`project_dir, experiment, sample, check, metric, value` - pivot however's needed downstream
(e.g. `df.pivot_table(index=['experiment','sample'], columns=['check','metric'], values='value')`).

Judgment call worth revisiting: `compute_duplication_rate.py` reduces `amplicon_stats.txt`'s
per-amplicon `reads_per_state` into a single **pooled** (sum of totals / sum of uniques) per-sample
rate, not an unweighted mean of each amplicon's ratio - chosen so low-coverage amplicons don't get
equal say to well-covered ones. This differs from the background-CpG sample-level number, which
*is* an unweighted mean across amplicons (a rate, not a raw-count ratio, so equal weighting made
more sense there) - two different reduction strategies for two different kinds of metric, not an
inconsistency to "fix."

## Downstream scripts: current state

These consume the core pipeline's `matrices/*.dedup.full_unclustered.matrix` output and are run
by hand (see README "Downstream analyses"). As of 2026-07-01 they are mid-rewrite; treat the
below as a snapshot, not a settled architecture.

### Binding classification (five generations; v2_carosversion canonical for the paper, v5 in active dev)

| Script | Status |
|---|---|
| `classify_single_molecule_binding_v2_carosversion.py` | **Believed canonical for the paper** — simplest, binary protection model. Referenced by `data/binding_model_standardization/` in the parent project. |
| `classify_single_molecule_binding_v2.py` | Experimental 4-state extension; the "important open" state and TF/nuc weighting are both dead code due to bugs (see `NOTES_classify_binding_scripts.md`). Not believed to be used for real results. |
| `classify_single_molecule_binding_v3.py` | Redesign with fixed-length (140bp) bin-based nucleosome placement + partial EM. EM is broken (hardcoded overwrite of the fitted parameter). Not wired anywhere. |
| `classify_single_molecule_binding_v4.py` | In-progress clean-room consolidation of the above three, fixes several cataloged bugs, adds a real EM step. **Not wired into Snakemake, not validated at real amplicon scale** — the state-space-explosion concern from the old notes is still open. See `workflow/scripts/NOTES_v4_status.md` before touching this. |
| `classify_single_molecule_binding_v5_hsmm.py` | **Current active development.** Segmental hidden-semi-Markov (per-molecule Viterbi over OPEN/NUC/TF/UNID), replacing v4's enumerate-microstates approach; `UNID` calls off-motif "unidentified" footprints. Nucleosome footprint width is tied to segment length (no `nuc_d_edge`). Has an optional empirical-Bayes **footprint-discovery** pass (`--discover_footprints 3pass`): decode → nucleosome-decontaminated bulk → discover recurrent off-motif footprints (local-background-relative, strength+count-graded) → re-decode with a position-specific reduced UNID start cost; always emits a per-bp NUC/TF/UNID-occupancy + raw/decontam-bulk diagnostic page. For robustness, `discover_footprint_vocabulary.py` builds a **POOLED per-promoter** vocabulary (pools decontaminated bulk across samples + copy-number variants, discovers once per promoter over the coord-conserved region upstream of the array) that the classifier applies via `--discovered_footprints_file` — avoids noisy small-N per-amplicon bulks. See `workflow/scripts/NOTES_v5_hsmm_status.md` (params) + `EXPLAINER_hsmm_v5.md` (concepts) + `NOTES_v5_calibration_log.md`. **Positions files must use bare-operator coords** (`convert_fa_to_positions_for_script.py --l_offset 0 --r_offset 0`); the old `+2/-1`-asymmetric opJS45 file causes TetO mis-calls — corrected file at `data/binding_model_code_updating/260713_opJS45.positions.long.bare.txt`. Optional cosmetic motif overlay via `--motif_track_file` (per-amplicon `>name`/`lo,hi,label[,color]` in the model bp frame; strip on diagnostic + per-read plots; opt-in, never affects decoding; producer `export_motif_track.py` in the P110 project builds the opoBD9 curated track from `motif_annotations.py`). Densely-bound-TetO undercalling was investigated 2026-07-14 (TF-penalty sweep + P110 ±dox swallow/obvious-FN analysis) and found ~negligible (model 99.4% right on unambiguous cases) — see `NOTES_v5_hsmm_status.md`. |

Full bug catalog for v2/v2_carosversion/v3: `workflow/scripts/NOTES_classify_binding_scripts.md`.
v4 status/next-steps: `workflow/scripts/NOTES_v4_status.md`.
v5 status/params: `workflow/scripts/NOTES_v5_hsmm_status.md`.

All variants **except v5**: enumerate physically-valid chromatin microstates (nuc + TF
combinations) from a `positions.txt` file (built via `convert_fa_to_positions_for_script.py` from
the amplicon FASTA + a motif FASTA), then assign each molecule its max-likelihood microstate under
a Bernoulli protection model. Output: `{sample}.{amplicon}.single_molecule_classification.txt`
(per-molecule state) + `{sample}.{amplicon}.valid_states_table.txt` (the enumerated state space).

### v5 output schema + downstream design (discussion 2026-07-15)

Full design writeup (with the partition-function math typeset):
`workflow/scripts/260715_v5_downstream_design.md` / `.pdf` (source + PDF; regenerate the PDF with
`ml system texlive/2019 && pandoc <md> -o <pdf> --pdf-engine=pdflatex`).

**v5 emits a deliberately new, non-legacy schema** (does NOT enumerate a global microstate space,
so **no `idx` column and no `valid_states_table.txt`**):
- Main per-molecule file (same `single_molecule_classification.txt` filename slot, indexed by
  `read_id`): `tfbs_1..K` (bool), `n_tf`/`n_nuc`/`n_unid` (counts), `nucs` (`start:end:dyad;...`),
  `unids` (`start:end;...`), `log_likelihood`.
- `.segments.txt`: tidy long, one row per Viterbi segment (`read_id, seg_index, type, start, end,
  dyad_or_motif, seg_loglik`) — the preferred substrate for nucleosome/occupancy analyses.

**Compatibility:** `aggregate_binding_model.py` runs on v5 output unchanged (it only does
`df.filter(like='tfbs_').sum()`; none of v5's extra columns contain `tfbs_`; peripheral
accessibility comes from the matrix, classifier-independent). Nucleosome analyses must move off the
old ragged `nuc{N}_start/end` columns to the `nucs` string / segments file. **The partition
function fit is conceptually incompatible** — `fit_partition_function_model_v3.py` needs
`valid_states_table.txt` + `idx`, which v5 doesn't produce. Decision: do NOT coerce v5 back to the
old schema; write new non-legacy downstream code.

**Planned (not yet built) downstream work:**
- **v5-native aggregator** (new script, reads main file + `.segments.txt`): tidy 3-level output
  losing no columns, adding `nuc_bases_covered`/`frac_nuc`, `unid_bases_covered`/`frac_unid`,
  `open_bases`/`frac_open`, and a `footprint_bases_covered`/`frac_footprint` union track
  (TF∪UNID∪NUC) = the "is there a footprint here" query. `frac_*` denominator = decoded span
  (first-to-last GpC), not raw amplicon length.
- **Nucleosome meta-plot**: occupancy(b) = (#molecules with a NUC segment spanning b) / (#total
  molecules). Denominator is clean — Viterbi assigns a state to every base of every molecule, so
  no `-1`/no-info correction; every `read_id` appears in the segments file.
- **UNID footprints**: single recurrent sites currently under-called because UNID entry cost
  (`start_cost_unid=6.0`) needs ~4 protected GpCs (each worth ≈+1.5 LLR) to clear, vs `start_cost_tf=1.0`.
  Two levers, same idea: (1) `--discovered_footprints_file` lowers UNID entry cost position-specifically;
  (2) **name the motif → it calls as TF** (entry cost 1.0, exact motif-window width, not bounded by
  `unid_min=15bp`). Endgame: name every recurrent footprint into the TF file; UNID persists as the
  residual/discovery channel, never forced to zero.

### ⭐ Running the v5 → Ising → potency chain on new samples

**Start here: `workflow/scripts/HOWTO_v5_ising_potency.md`** — the end-to-end recipe (pipeline →
v5 decode → Ising fit → potency + `w`), with every trap and every "change this for a new
construct" item. Reference run, with samplesheet/config/sbatch/NOTES.md to copy:
`../data/260903_v5_hsmm_partition_function_test_old_data/`.

**Validated 2026-09-04 by reproducing Nature 2024 Fig. 3 from the original 2022 opJS45 FASTQs**:
potency `k_pot` 0.085–0.103 (published ~0.090), basal 0.027–0.035 (~0.03), and nucleosome
displacement `delta` = −1.09 recovered as the dominant term — the paper's own conclusion, reached
here by a completely different route (moment matching + transfer matrix vs enumerate-microstates
ML). New result: `w` ≈ 0.33–0.54, the fraction of potency that is *integrated* rather than
instantaneous, agreeing across two independent axes.

Three things that bite hardest, all documented in the HOWTO:
- **`no_endog_meth` must be decided from each run's own background-CpG QC**, not assumed. It was
  14–17% in the 2022 data, so GCG-ambiguous GpCs had to be dropped. The ambiguous positions
  (opJS4: GpC 54 and 130) miss the TetO array entirely but 130 sits in the TSS window — so the
  Ising/array fit was unaffected while promoter/potency was not. Re-derive per amplicon.
- **`mu_ambient` is a log-fugacity in kT, not a log-probability.** Measure ambient nucleosome
  coverage on a zero-motif amplicon and invert with
  `fit_ising_model.py mu_ambient --coverage <x>`. Re-measure per construct.
- **Trust the `.ising_model_checks.pdf` over AIC.** The likelihood only sees the mean moments, so
  AIC cannot see a variant wrecking the out-of-sample predictions — we hit exactly that case.

New scripts from this work: `fit_ising_model.py` (49 self-tests; `selftest` subcommand),
`fit_two_timescale_potency.py`, `260902_w_regimes_cartoon.py`.

### Partition function fit — REPARAMETRIZATION PLANNED (Ising/MaxEnt)

Docs, in reading order: `260715_v5_downstream_design.md` §7 (formal derivation) →
`260717_ising_intuition.md` (concepts; **§7 is superseded, see the potency note below**) →
`260717_fit_lattice_gas_ising_sketch.py` (runnable architectural sketch; `python <script>` with no
args runs a synthetic self-check).

**Old approach (`fit_partition_function_model_v3.py`) is being replaced for v5.** The old enumerate-
microstates + fit-Boltzmann-weights approach IS already an Ising/MaxEnt model, just parametrized by
enumeration. Plan: keep energy = Σ_k θ_k φ_k(σ) linear in discrete interaction features, but **fit
by moment-matching** (⟨φ_k⟩_model = ⟨φ_k⟩_data) — data enters only through empirical feature
averages read directly off v5 per-molecule calls (no `idx`, no enumeration). Fit is convex (log Z is
log-sum-exp), and the Hessian (= N·Cov_model(φ)) gives free error bars. Z computed exactly via
**transfer matrix** for the 1D TetO chain and a **hard-rod lattice-gas DP** for nucleosomes.

**Decisions settled (do not re-open):**
- **Hard-rod lattice gas is the nucleosome model.** Coarse region indicators are a *disposable*
  one-off cross-check to reproduce the old fit, explicitly not a stepping stone — they can't
  represent packing or a reservoir.
- **Grand-canonical boundary, never hard walls at the amplicon edge.** Fugacity + pad ≥ one
  footprint each side + bulk-occupancy seed. Correctness check: a barrier-free uniform stretch must
  give *flat* occupancy. Hard walls manufacture fake edge phasing.
- **Shared operator field h — no free per-site h_i.** TetOs are sequence-identical; the middle>edge
  gradient must *emerge* from excluded volume as a prediction.
- **Model flavors = a boolean mask over one superset feature vector**, not separate
  `assign_energy_*` functions. Makes LRT nestedness structural and tells you which posterior checks
  are tautological (fit moments) vs informative (everything else).

**Three different δ's — don't conflate them** (the old `--model` names are ambiguous):
- **δ_steric**: TF under a nucleosome. **Hard exclusion, not a fitted parameter** — v5's Viterbi
  segmentation makes the configuration unrepresentable, so its observed count is structurally zero.
- **δ_soft-steric**: relaxing the above. Unidentifiable here; ignore.
- **δ_remodel**: what the old `3param_nuc` actually fits — `num_nuc * (nuc_e - delta_e *
  (num_tfs>0))`, i.e. nucleosome fugacity shifted when *any* TF is bound (recruited remodeler).
  **This is the one worth fitting.** Non-local, but exact two ways: (a) three DP calls via
  Z = Z_noTF(μ) + Z_all(μ+δ) − Z_noTF(μ+δ) (needs a log-difference, loses precision when few TFs
  bind), or (b) augment the DP state with `(has_tf, n_nuc)` (~3.5k states, numerically robust).
  Build (b), keep (a) as a fast path, and unit-test them against each other.

**⚠️ Potency is NOT an Ising coupling** — `260717_ising_intuition.md` §7 ("promoter as one more
node, potency = J_prom in k_BT") is **superseded** by `260902_two_timescale_potency.md` (+ `.pdf`,
12pp) and the cartoon `260902_w_regimes_cartoon.py`. Two independent reasons, both from published
Fig. 3: the occupancy→activity link is line-like, not convex (so not a Boltzmann coupling — no
promoter parameter is a free energy), and promoter activity depends on the *available* site count
as well as the currently-bound count (time integration, which no single-snapshot equilibrium model
can represent). Current model: P = p_0 + k_pot[(1−w)·n + w·⟨n⟩], where the w terms **cancel exactly**
in the bulk average — so the published potency = k_pot stands unrevised, while any same-molecule
coupling only ever measures k_pot(1−w). The equilibrium *array* model (h, J, μ, δ_remodel) is
unaffected; only the promoter node leaves equilibrium.

### Partition function fit (old / v2-era)

`fit_partition_function_model_v3.py` — fits a thermodynamic (Boltzmann) model over the enumerated
microstates via `scipy.optimize.minimize`, in three flavors (`2param`, `3param_nuc` [nucleosome
disruption by bound TF], `3param_tfcoop` [TF cooperativity]). Consumes the binding classifier's
`valid_states_table.txt` + `single_molecule_classification.txt`; outputs a 2-column
`{sample}.fit.{model_lib}.independent_fit.txt` (param name, fitted energy value) per sample.

### New consolidation layer (untracked, being actively built out)

Five new scripts implement a **three-level aggregation pattern** (single-molecule → per
sample×amplicon → per sample) that's consistent across the pipeline:

- **`aggregate_binding_model.py`** — collates binding classifications across samples/amplicons;
  computes `n_bound`, `tf_bound`, `peripheral_accessibility` (GpC accessibility downstream of the
  last TFBS, a chromatin-remodeling proxy). Has hardcoded amplicon names
  (`opJS4_6x_TetO_21bp_no_CG`, `opJS4_4x_TetO_21bp_no_CG`) baked in for sample-level extraction —
  will silently produce missing rows if your amplicon set doesn't include those exact names.
- **`aggregate_partition_function.py`** — collates `fit_partition_function_model_v3.py` outputs
  across samples into one params-by-sample table. Silently skips samples with no fit file; no
  check that all samples were fit with the same model type (columns would be ragged if not).
- **`assign_promoter_state_from_model_new.py`** — the most complex of the five. Detects and
  strips nucleosome protection over the promoter (`get_all_protection_streaks`, configurable
  `--min_nuc_len`/`--length_method`), then runs the nucleosome-free promoter pattern through
  *two* k-means models (a 40-cluster legacy model + a 20-cluster new model), maps clusters to
  named states (`pic`/`tbp`/`pause`/etc.) via a lookup table, and calls TATA/TSS/pause footprints
  by mean occupancy in fixed windows (`tata=(155,170)`, `tss=(120,135)`, `pause=(80,100)`, hardcoded
  at file top). Depends on three pickled/text model files with hardcoded Stanford absolute paths
  as CLI defaults — check these still exist/are current before relying on the defaults. Supersedes
  the older, much simpler `assign_promoter_state_from_model.py` (single k-means pass only, no
  nucleosome handling).
- **`compute_potency.py`** — fits either a mechanistic saturation curve (`frac_on = (k_basal +
  k_tf·n_bound) / (1 + k_basal + k_tf·n_bound)`, via `curve_fit`) or a linear model (weighted
  least squares) of promoter-open-fraction vs. TF-bound-count, per sample. This is what produces
  the "potency" scalar referenced throughout the parent project's `CLAUDE.md`.
- **`consolidate_data.py`** — outer-joins all of the above three-level outputs into final
  `*.single_molecule.txt.gz` / `*.amplicon_level.txt.gz` / `*.sample_level.txt.gz` tables. Missing
  input tables are silently dropped (resulting columns are just NaN) rather than erroring — worth
  double-checking which inputs were actually supplied if a consolidated table looks sparse.

None of the five have Snakemake rules yet; they're invoked by hand, presumably from a notebook or
shell script that isn't in this repo (check `activation_domain_smf_paper/notebooks/` in the
parent project).

### `common.py` (shared utilities, currently mid-edit — see `git diff`)

Recently gained three functions supporting the new promoter-state clustering:
`assign_cluster()` (predict from a pickled k-means model), `annotate_40_clusters()` (numpy
broadcasting version, replacing a commented-out pandas `.isin()` version — rewritten for memory
efficiency on large (~5M row) molecule tables), `assign_promoter_footprints()` (mean-occupancy
footprint calling used by `assign_promoter_state_from_model_new.py`). Also has known technical
debt flagged in its own comments: `cluster_single_reads()` calls seaborn's `clustermap` just to
get a linkage/dendrogram, which the code's own comment calls "super clunky."

---

## Known bugs (glaring, found 2026-07-01)

Verified by direct code read (not just skimmed) — line numbers are current as of this writing.

### 1. Data-corrupting: `every_other=args.no_dedup` in `assign_promoter_state_from_model_new.py:461`

```python
mat = load_single_molecule_matrix(p, every_other=args.no_dedup)
```

`load_single_molecule_matrix()`'s `every_other` flag exists to drop what used to be duplicate
GpC columns from the *old* matrix format ("Georgi's script has 2 columns per GpC" — see
`common.py:198`). The matrix format was migrated to one column per GpC (`git log`: "changing the
load single molecule every other to be False since now we only have one column per GpC"), and
every other caller in the repo now correctly hardcodes `every_other=False`
(`aggregate_binding_model.py:69`, `plot_bulk_methylation_from_matrices.py:32`). This one call site
was missed: whenever `--no_dedup` is passed, `every_other` becomes `True`, and under the current
one-column-per-GpC format this **silently drops every other real GpC position** from the
promoter matrix — not a duplicate-removal, an actual data-loss bug. Anyone who has ever run
`assign_promoter_state_from_model_new.py --no_dedup` has gotten a promoter matrix missing ~half
its real GpCs, with no error or warning. Fix: hardcode `every_other=False` here like the other
callers.

### 2. Crash risk: unguarded division by zero in `mark-nonconverted-reads-and-plot.py:299`

```python
converted_fraction = 1 - (my_unconverted / my_total_c)
```

`my_total_c` (count of non-GpC/non-CpG Cs seen on the read) has no floor check. Any read where
this count is 0 — short reads, heavily soft-clipped reads (see `softclip()`), or reads that
happen to land entirely in a GpC/CpG-dense stretch — raises an unhandled `ZeroDivisionError` and
kills the whole `filter_uncoverted` rule for that sample. This runs unconditionally in the core
pipeline (`other.smk` rule `filter_uncoverted`), including for `deaminase` samples (only the
threshold `c_frac` changes for those, not whether this code path runs).

### Minor / worth a look

- **Dead code, not a correctness bug:** `mark-nonconverted-reads-and-plot.py:134,159` — the
  original `--c_count`-based pre-filter (`if read...count("C") ... >= args.c_count:`) was replaced
  with a bare `if True:`, so `filter_snps()` now runs on every single read. `--c_count` is a
  vestigial, unused CLI argument.
- **Silent sample loss:** `aggregate_partition_function.py:29-34` skips any sample whose
  `{sample}.fit.{model_lib}.independent_fit.txt` doesn't exist, with no warning printed — a sample
  can silently vanish from the aggregated table.
- **Bare crash on missing required args:** `consolidate_data.py` does `pd.read_table(args.promoter_single_molecule, ...)`
  (line 35) with no check that the arg was actually provided — passing an incomplete set of
  `--*` flags fails with a raw pandas/Python exception instead of a usable error message.

### Already-catalogued, still live (see `NOTES_classify_binding_scripts.md` for full detail)

The single-molecule binding classifiers have several confirmed, more severe bugs from an earlier
review pass — flagging the two worst here since they silently change scientific results rather
than just crashing:
- **v2 `classify_single_molecule_binding_v2.py` B4:** the 4-state model's "important open" state
  (code `3`) can never be assigned — `enhancer_upper`/`enhancer_lower` are swapped, so the
  triggering condition is always false. The 4-state model silently degrades to 3 states.
- **v3 `classify_single_molecule_binding_v3.py` B11:** `prob_unmeth_given_open` is correctly
  estimated from data, then immediately overwritten by a hardcoded `0.05` on the next line — the
  EM loop for this parameter never actually learns anything, despite appearing to.

Per N2 in that file, `classify_single_molecule_binding_v2_carosversion.py` (believed canonical
for the paper) doesn't have the EM/4-state code at all, so B4/B11 likely don't affect published
results — but worth confirming which classifier version any given downstream table actually came
from before trusting it.

---

## Where to look next

- Parent project context (paper, notebooks, data conventions): `../CLAUDE.md`.
- Binding-classifier bug catalog: `workflow/scripts/NOTES_classify_binding_scripts.md`.
- v4 classifier rewrite status: `workflow/scripts/NOTES_v4_status.md`.
- README.md in this directory has the authoritative samplesheet/config column docs and the
  downstream-analysis CLI invocation examples — schema files under `workflow/schemas/` are stale,
  do not rely on them.
