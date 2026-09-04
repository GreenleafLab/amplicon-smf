# Notes: `classify_single_molecule_binding_v4.py` — status as of 2026-06-30

Follow-up to [`NOTES_classify_binding_scripts.md`](./NOTES_classify_binding_scripts.md) (the
bug catalog for carosversion/v2/v3). This file tracks where the v4 rewrite stands.

## What v4 is

A clean-room consolidation of the three older classifiers — not a patch on any one of them.
Takes v3's continuous distance-based nucleosome idea and reworks it with a cleaner state
representation, a real EM step, and explicit fixes for several bugs cataloged in the notes file
above.

**Not yet wired into the Snakemake pipeline** — no rule in `workflow/rules/` or the Snakefile
references it yet. Still a standalone script.

## Bugs fixed vs. old scripts (confirmed by reading the code)

| Bug | Old problem | v4 fix | Where |
|---|---|---|---|
| B1 | `--reads_to_use` parsed but never applied | Actually subsamples | [`classify_single_molecule_binding_v4.py:769`](./classify_single_molecule_binding_v4.py#L769) |
| B6 | `list == tuple` comparison for `(0,0)` always False | Parses to list, compares `promoter_lo < promoter_hi` | [`classify_single_molecule_binding_v4.py:782`](./classify_single_molecule_binding_v4.py#L782) |
| B8 | `expand_state` loop vars shadowed the `idx` param | Renamed to `bin_k`/`tf_k`, `state_idx` preserved | [`classify_single_molecule_binding_v4.py:482`](./classify_single_molecule_binding_v4.py#L482) |
| B11 | EM param estimate immediately overwritten by hardcoded stub | Real inversion of observed T-fraction through conversion model (`invert_t_fraction`), with damping + `min_obs` guard | `run_em_mstep` |
| B12 | EM debug files dumped to cwd every iteration | Only written if `em_log_path` given; defaults to `<output>.em_log.tsv` | `compute_classifications` |

## Design changes vs. all three predecessors

- Protection model = three additive "channels" (nuc / TF / open) combined via `max`. This is
  safe *because* `prune_states` already guarantees bound TFs and nucleosomes never physically
  overlap within a surviving state — so at most one non-open channel is ever active at a given
  GpC.
- Nucleosome protection is a continuous logistic-distance function (`nuc_d_edge`, `nuc_softness`)
  fit by `scipy.optimize.minimize` (L-BFGS-B) during EM — this is the actual thing EM tunes,
  unlike v3 where it was faked.
- `amp_width` auto-detects from the data (`auto_detect_amp_width`) instead of being a hardcoded
  module-level global — fixes the old "non-reusable across amplicons" complaint.
- Two-pass argparse: optional `--config` YAML supplies defaults, CLI flags override.

## Real-amplicon test run (2026-06-30) — N1 resolved, EM still untested

First real-scale run lives at
`/oak/stanford/groups/wjg/bgrd/papers/ad_smf/data/binding_model_code_updating/260507_v4firstrun`
(driver: `log.sh` there). 5 samples (L647–L651) × ~29 amplicons-with-reads (opJS4 0–8x, CTCF,
opJS5 0–9x b1/b2), ~55 min/sample, ran **pre-EM** (no `--do_em`, so no `em_log` — expected, not a
bug; empty BD24/opJS4_3x outputs are just no-reads). Invocation used `--promoter_positions 75,175`,
`--p_t_given_unmeth 0.99 --p_t_given_meth 0.05`.

- **N1 (state-space explosion) is NOT a problem for these amplicons.** Enumerated state counts
  (rows in `valid_states_table`) are the same order of magnitude as the old canonical
  carosversion classifier on the same matrices — *smaller* for low TetO (opJS4_0x: 105k→41k), ~1x
  at opJS4_8x (331k→371k), ~2x only at the extreme opJS5_9x (263k→516k). Since carosversion
  already ran 331k states in production, v4 is not a new scaling risk here. Did not hang/OOM.
- The alarming ~114 MB `valid_states_table.v4.txt` sizes are **column width, not row count**: v4
  writes a wide fixed-slot format (`nuc0..nuc67` × present/start/end ≈ 205 cols). On-disk size
  overstates state count ~5x. Keep in mind when loading these tables.
- `recursive_1_placer` is still an uncapped recursive combinatorial placer pruned against TF
  states only after generation — fine at observed scale, but no guardrail if a much larger amplicon
  is ever passed.
- **EM path (B11/B12 fixes) is coded but still UNEXERCISED on real data** — this run was pre-EM.
  Running with `--do_em` on one real sample is now the top open item.
- `test_v4_smoke.py` docstring says "Do not commit this file" but it is **not** in `.gitignore`
  — currently shows as untracked in `git status`. Decide whether to gitignore it, promote it to a
  real test, or delete it before it accidentally gets committed.
- Not yet integrated into any Snakemake rule — need a rule + config wiring if/when this is meant
  to replace v2/v3/carosversion in the actual pipeline.

## Other untracked files in this dir (broader consolidation effort, not yet reviewed in depth)

These match the pipeline script table in the project's `CLAUDE.md` and appear to be part of the
same in-progress cleanup pass as v4, but haven't been read/reviewed yet:

- `aggregate_binding_model.py`
- `aggregate_partition_function.py`
- `assign_promoter_state_from_model_new.py`
- `compute_potency.py`
- `consolidate_data.py`

## Next steps (pick one to start with next session)

1. ~~Stress-test v4 on a real single-molecule matrix~~ — DONE 2026-06-30 (see section above); N1
   fine, runtime ~55 min/sample pre-EM.
2. **Exercise the EM path on real data** (`--do_em` on one real sample) — coded but never run at
   scale; verify B11/B12 fixes actually converge and that `em_log.tsv` is sane. Top open item.
3. Decide `test_v4_smoke.py`'s fate (gitignore / promote / delete).
4. If v4 is meant to become canonical, add a Snakemake rule for it and decide whether it replaces
   or runs alongside `classify_single_molecule_binding_v2_carosversion.py` (currently believed to
   be the paper's canonical binding model per N2 in the other notes file).
