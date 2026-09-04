# Notes: `classify_single_molecule_binding` Script Comparison

Written: 2026-05-06

These notes compare the three versions of the single-molecule binding classification script and document bugs, edge cases, and performance issues found. **Nothing has been changed** — this is documentation only.

---

## Overview: What the scripts do

All three scripts implement the same high-level algorithm:

1. Load a single-molecule methylation matrix (GpCs × molecules)
2. Enumerate all physically valid chromatin microstates (combinations of nucleosome and rTetR occupancy)
3. For each molecule, compute the maximum-likelihood microstate given its observed methylation pattern (Bernoulli likelihood)
4. Write per-molecule state assignments to a TSV

The core output is: for each molecule, which rTetR binding sites are occupied, and which nucleosomes are present.

---

## Script-by-Script Summary

### `classify_single_molecule_binding_v2_carosversion.py` ("Caro's version")

The simplest of the three. Uses the original binary protection model: each GpC is either protected (1) or accessible (0). Nucleosome states are enumerated by finding all valid GpC-to-GpC spans (in reference coordinates) that fall within a size window.

- **Nuc span:** `(110, 140)` bp
- **Likelihood model:** Binary: bound → high p(protected); unbound → low p(protected). Promoter region gets a special "don't-care" treatment for accessible positions.
- **Promoter positions:** Hardcoded as `(75, 175)`, cannot be changed via CLI
- **Weights/EM:** None. Simpler than v2 in this respect.
- **CLI additions vs. v2:** Missing `--individual_reads_to_plot`, `--fix_ambiguous_gcg`, `--promoter_positions`

### `classify_single_molecule_binding_v2.py`

Extends Caro's version with a 4-state protection code to differentiate types of protected/accessible bases, allowing more nuanced likelihood assignments. Also adds flanking-GpC upweighting to penalize states that get the TetO-adjacent bases wrong.

- **Nuc span:** `(110, 147)` bp — slightly wider than Caro's version
- **Likelihood model:** 4-state codes: `0`=open, `1`=nuc, `2`=TF, `3`=important-open (between TetOs). Each gets a different Bernoulli probability. The TF-bound probability is softer (`p_t_given_bound * 0.85`) and the "important open" sites are penalized for being protected.
- **Promoter positions:** Passed as CLI arg `--promoter_positions`, default `"75,175"`
- **Weights/EM:** Computes `tf_weights`/`nuc_weights` diagonal matrices, but the weighted matmul line in `classify_all_molecules` is **commented out** — weights are computed but never applied.

### `classify_single_molecule_binding_v3.py`

A substantial redesign. Nucleosomes are now parametrized by their midpoint position (in bins of 10 bp) rather than GpC spans, giving fixed-length 140 bp nucleosomes. Introduces a logistic regression model for position-dependent methylation probability within a nucleosome. Has a (partially implemented) EM loop to estimate hyperparameters from data.

- **Nuc model:** Fixed-length (140 bp), bin-based midpoints on a 600 bp amplicon
- **Likelihood model:** Explicitly separates p(methylated | state) from p(T | methylated) and p(T | unmethylated). Nucleosome protection probability is a logistic function of distance to nuc center. TF protection is a flat probability within TFBS bounds.
- **EM:** Optional `--do_em` flag. The EM loop iterates MLE state assignments → parameter re-estimation → convergence check. Partially implemented (see bugs below).
- **Global state:** Uses module-level globals `amp_width=600, nuc_length=140, bin_size=10` defined at import time, which makes the code non-reusable across different amplicons without code changes.

---

## Key Differences Between Scripts

| Feature | carosversion | v2 | v3 |
|---|---|---|---|
| Nuc parametrization | GpC-span (variable length) | GpC-span (variable length) | Bin-midpoint (fixed 140 bp) |
| Nuc size | 110–140 bp | 110–147 bp | Fixed 140 bp |
| Protection codes | Binary (0/1) | 4-state (0/1/2/3) | Continuous LR probability |
| TF likelihood | Same as nuc | Softened (0.85×) | Separate p(unmeth|TF) param |
| Promoter treatment | Hardcoded `(75,175)` | CLI arg, default `(75,175)` | Hardcoded default `(0,250)` |
| Weights | None | Computed but not applied | N/A |
| EM parameter estimation | No | No | Yes (partially) |
| Chunked matmul | Yes | Yes | Yes |
| Precomputed states | Yes (pickle) | Yes (pickle) | Yes (pickle) |

---

## Bugs and Issues

### All three scripts

**B1: `--reads_to_use` argument is silently ignored.**
All three scripts accept a `--reads_to_use` argument promising to subsample the matrix, but none of them ever actually apply it. The parsed value is never used after `args = parser.parse_args()`. If someone passes `--reads_to_use 5000` expecting a downsampled run, they will silently get the full matrix.

**B2: `find_first_matching_sublist` does not early-exit when span is already too large.**
(Applies to carosversion and v2 which share this function.)
```python
while end < len(lst):
    temp_span = lst[end] - lst[start] + 1
    if temp_span >= min_span and temp_span <= max_span:
        return end
    else:
        end += 1  # keeps going even after temp_span > max_span
```
Since `lst` is sorted, once `temp_span > max_span` it can only increase — continuing the loop is wasted work. The fix is to add `elif temp_span > max_span: return -1` to the else branch. For a long list this can cause unnecessary O(N) scans that could be O(1) exits.

**B3: Promoter probability hack leaves the likelihood un-normalized.**
(Applies to `create_bernoulli_logprob_matrices` in carosversion and v2.)

The promoter "hack" sets `prob_states_unmeth[idx,i]` = `(1 - p_t_given_unbound)` for accessible promoter positions, but leaves `prob_states_meth[idx,i]` unchanged. Before the hack, `prob_states_unmeth + prob_states_meth = 1` elementwise. After the hack, for an accessible promoter GpC:
- `prob_states_unmeth = (1 - p_t_given_unbound)` ≈ 0.94
- `prob_states_meth = (1 - p_t_given_unbound)` ≈ 0.94 (unchanged, computed with p_m=0)
- Sum ≈ 1.88, not 1.

For argmax classification this may not change outcomes (the same offset applies to all states at that position), but it's conceptually wrong and could affect any downstream likelihood-based comparisons.

---

### `classify_single_molecule_binding_v2.py` specific

**B4: "Important open" (code 3) is never assigned — core 4-state feature is broken.**

In `generate_predicted_protection_more_states`:
```python
enhancer_upper, enhancer_lower = flattened_flanks[0]-1, flattened_flanks[-1]+1
for idx, p in enumerate(methyl_positions):
    if idx >= enhancer_lower and idx <= enhancer_upper:  # ← impossible condition
```

`enhancer_upper` gets the *smaller* index (`flattened_flanks[0]-1`) and `enhancer_lower` gets the *larger* index (`flattened_flanks[-1]+1`). So `enhancer_lower > enhancer_upper` always, making the condition `idx >= enhancer_lower and idx <= enhancer_upper` **always false**. No GpC ever receives protection code `3`.

The intended logic (based on context) was to mark all bases *between* the TetO-flanking GpCs as "important open" — i.e., accessible bases flanked by TetOs that should strongly constrain what states are plausible. As written, the 4-state model silently collapses to 3 states (0, 1, 2 only), and the entire `np.where(p_m_array == 3, ...)` branch in `create_bernoulli_logprob_matrices_more_states` is dead code.

Fix: swap the variable assignment or flip the comparison:
```python
# Option A: fix assignment
enhancer_lower, enhancer_upper = flattened_flanks[0]-1, flattened_flanks[-1]+1
# Option B: fix check
if idx >= enhancer_upper and idx <= enhancer_lower:
```

**B5: `compute_weights` mixes index and position — weights always equal 1.**

```python
enhancer_upper, enhancer_lower = flattened_flanks[0], flattened_flanks[-1]
enhancer_upper_idx, enhancer_lower_idx = methyl_positions.index(enhancer_upper), methyl_positions.index(enhancer_lower)
...
if idx >= enhancer_upper_idx - 1 and p <= enhancer_upper_idx + 1:
```

Here `idx` is a loop index (small integer, e.g. 0–30) and `p` is a genomic position (potentially hundreds to thousands). `enhancer_upper_idx` is also a small integer (array index). The condition compares `p` (a position) to `enhancer_upper_idx + 1` (an index + 1) — since genomic positions >> array indices, `p <= enhancer_upper_idx + 1` is almost certainly never true for any real data. All positions get weight 1.

This is moot anyway because the weighted matmul in `classify_all_molecules` is commented out:
```python
res = np.argmax(np.exp(log_prob_states_unmeth @ sub_mat + ...), axis=0)
# res = np.argmax(np.exp(log_prob_states_unmeth @ nuc_weights @ ...), axis=0)  # ← never runs
```
So `compute_weights` and both weight matrices are **computed but never used** — dead code.

**B6: `promoter_positions == (0,0)` check compares list to tuple — always False.**

```python
promoter_positions = list(map(int, args.promoter_positions.split(',')))
if promoter_positions == (0,0):   # ← list [0,0] != tuple (0,0) in Python
    promoter_positions = None
```

`list(map(int, ...))` returns `[0, 0]` (a list). Comparing to the tuple `(0,0)` is always `False` in Python. Passing `--promoter_positions 0,0` as documented (to disable promoter handling) **silently does nothing** — the promoter hack remains active. Fix: compare to `[0, 0]` instead.

**B7: `get_flanking_gcs` assert will crash for non-standard amplicons.**

```python
assert len(pos) == 2
```

This assertion hard-codes the assumption that every TFBS has exactly 2 GpC positions. If the amplicon lacks a GpC flanking one side of a TetO (or has more than 2), this will crash with an AssertionError and no helpful message.

---

### `classify_single_molecule_binding_v3.py` specific

**B8: `expand_state` shadows its own `idx` parameter with loop variables.**

```python
def expand_state(tf_state, nuc_state, idx, tfbs_positions, ...):
    ...
    for idx in nuc_state.index:   # ← overwrites `idx` parameter
        expanded_state['nuc{}_present'.format(idx)] = ...
    for idx in tf_state.index:    # ← overwrites again
        expanded_state['tfbs_{}'.format(idx)] = ...
    expanded_state['idx'] = idx   # ← this is now the LAST TF index, not the state index!
```

The `idx` parameter (intended to store the state's index for downstream lookup) is shadowed by both loop variables. The final `expanded_state['idx']` will always be the last value from `tf_state.index` (e.g., the index of the last TF site), not the state's position in the valid-states list. Any downstream code that uses `assignments['idx']` to look up states is silently getting wrong values.

**B9: `write_valid_states` calls `expand_state` with the v2 signature — would crash if ever called.**

```python
def write_valid_states(states, out_path, tfbs_positions):
    for idx, state in enumerate(states):
        expanded_assignments.append(expand_state(state, idx, tfbs_positions))  # ← wrong signature
```

v3's `expand_state` expects `(tf_state, nuc_state, idx, tfbs_positions)` but this calls it as `(state, idx, tfbs_positions)`. This is a copy-paste from v2 that was never updated. It would raise a `TypeError` at runtime. However, in `compute_classifications`, the `write_valid_states` call is **commented out** (`# if valid_states_path: write_valid_states(...)`), so this is a latent crash that would surface if the `--all_states_output` functionality is ever re-enabled.

**B10: `prune_states` uses `tfbs_positions` as an undeclared implicit global.**

```python
def prune_states(tf_states, nuc_states):
    ...
    tf_dict = construct_tf_positions_array(tf_states, tfbs_positions, ...)  # ← not a parameter!
```

`tfbs_positions` is not in the function signature. It works when run as `__main__` because it's assigned as a local in the `if __name__ == '__main__':` block before `compute_classifications` is called. But if this module is ever imported and `compute_classifications` is called directly, this will raise a `NameError`. The function should take `tfbs_positions` as an explicit parameter.

**B11: EM estimation of `prob_unmeth_given_open` is immediately overwritten with hardcoded value.**

In `estimate_latent_params_from_state_assignments`:
```python
prob_unmeth_given_open = (obs_prob_t - p_t_given_meth) / (p_t_given_unmeth - p_t_given_meth)
prob_unmeth_given_open = 0.05   # ← overwrites the line above immediately
```

The first line correctly estimates `prob_unmeth_given_open` from the observed data. The second line discards it with a hardcoded value. The EM loop for this parameter is therefore broken — it never actually updates. The comment above the second line says "now need to convert..." suggesting this was left as a stub.

**B12: EM debug output writes files to cwd unconditionally.**

`estimate_latent_params_from_state_assignments` writes several debug files (`nuc_states.{rep}.txt`, `tf_states.{rep}.txt`, `tf_data.{rep}.txt`, `open_data.{rep}.txt`, `ml_assignments.{rep}.txt`) to the **current working directory** every EM iteration. These are presumably development artifacts. In production use on a cluster, this would scatter files wherever the job is launched from.

**B13: v3 promoter window default is different from v2 and undocumented.**

v3's `create_bernoulli_logprob_matrices` defaults to `promoter_positions=(0,250)`. v2 defaults to `(75,175)`. The difference is undocumented and the v3 value is not exposed as a CLI argument, so users cannot change it without editing source.

---

## Design / Conceptual Notes

**N1: State-space explosion in carosversion/v2.**
The `sublist_finder` recursion can generate very large numbers of states for long amplicons (exponential in the number of possible nuc placements). For a 6xTetO construct with a long amplicon, this is manageable but could become expensive. v3's fixed-midpoint binning gives a more predictable and typically smaller state count.

**N2: Caro's version is likely the "canonical" one for the paper's binding model.**
Based on the CLAUDE.md reference to `data/binding_model_standardization/` and "Caro version," this is probably the version actually used for the analyses in the paper. The v2 and v3 scripts contain partially-implemented experimental extensions (4-state codes, EM) that were apparently never finished or validated.

**N3: Nucleosome size parameter differences between scripts are meaningful.**
Caro's version uses `(110, 140)`, v2 uses `(110, 147)`. These are the span between the first and last GpC covered by the nucleosome, not the actual nucleosome length. v3 uses an actual fixed length of 140 bp and a bin size of 10 bp, which is a different parametrization entirely. Results would differ between scripts, especially at nucleosome boundaries.

**N4: v3's logistic regression for nucleosome protection is not trained from data in practice.**
`get_initial_em_hyperparams` initializes the LR from hardcoded fake training data (`X_tmp`, `y_tmp`). Without a completed EM loop, this never gets updated, so v3 uses a fixed LR with arbitrary parameters. The EM step could in principle learn these from the data, but per B11 it is broken.

**N5: Chunked matmul is consistent across all three scripts but the chunk size formula is conservative.**
`determine_chunk_size` returns `int((n_mols * n_states * data_size) // (memory_threshold * 1e9)) + 1`. This always produces at least 2 chunks (due to `+ 1`), even when a single chunk would fit in memory. This is slightly wasteful but not harmful.

**N6: `enumerate_states` in v3 returns `valid_tf_states, valid_nuc_states` as paired DataFrames.**
After `prune_states`, both DataFrames have the same number of rows, where row `i` jointly describes the TF configuration and nuc configuration of valid microstate `i`. This pairing is load-bearing — `create_bernoulli_logprob_matrices` in v3 relies on both DataFrames being aligned row-wise. This is a non-obvious invariant with no assertion to enforce it.
