#!/oak/stanford/groups/wjg/bgrd/bin/miniconda3/bin/python3.8
"""
classify_single_molecule_binding_v4.py

Single-molecule chromatin footprinting classification pipeline.

Data flow:
  1. Load & preprocess  -- load_single_molecule_matrix, filter, fix_missing_data, subsample
  2. State enumeration  -- enumerate_tf_states, enumerate_nucleosomal_states, prune_states
  3. Likelihood construction -- build_log_prob_matrices (nuc/TF/open protection channels)
  4. MLE classification -- classify_all_molecules (chunked matmul, argmax)
  5. EM (optional)      -- run_em_mstep, parameter damping, convergence check, re-classify
  6. Output             -- summarize_assignments, write_valid_states, plots
"""

import pandas as pd
import numpy as np
import argparse
import math
import os
import os.path
import pickle
import warnings
import copy
from copy import copy as _copy_shallow
from itertools import product
from dataclasses import dataclass

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')

from scipy.special import expit
from scipy.optimize import minimize

try:
    import yaml
except ImportError:
    yaml = None  # handled lazily in load_yaml_config

warnings.filterwarnings('ignore')

from common import (
    load_single_molecule_matrix,
    get_methyl_positions,
    load_tfbs_positions,
    filter_all_converted_reads,
    decorate_single_read_plot2,
    plot_single_read,
    adjust_gcgs,
    fix_missing_data,
)


# ---------------------------------------------------------------------------
# ModelParams dataclass
# ---------------------------------------------------------------------------

@dataclass
class ModelParams:
    amp_width: int
    nuc_length: int = 140
    bin_size: int = 10
    p_t_given_unmeth: float = 0.95
    p_t_given_meth: float = 0.15
    nuc_d_edge: float = 65.0
    nuc_softness: float = 5.0
    prob_unmeth_given_tf: float = 0.9
    prob_unmeth_given_open: float = 0.05
    promoter_lo: int = 0
    promoter_hi: int = 0
    prob_unmeth_given_open_promoter: float = 0.5


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def auto_detect_amp_width(methyl_positions, nuc_length):
    return int(max(methyl_positions)) + nuc_length // 2 + 10


def get_nuc_start_and_end_from_midpoint(bin_idx, bin_size, nuc_length):
    center = bin_idx * bin_size + bin_size // 2
    return center - nuc_length // 2, center + nuc_length // 2


# ---------------------------------------------------------------------------
# State enumeration
# ---------------------------------------------------------------------------

def recursive_1_placer(lst, spacing, start_idx):
    """
    Places 1s into an empty array of 0s subject to spacing constraints between adjacent 1s.
    Helper function for enumerate_nucleosomal_states; recursive.
    """
    to_return = [_copy_shallow(lst)]
    if start_idx < len(lst):
        for idx in range(start_idx, len(lst)):
            new_lst = _copy_shallow(lst)
            assert new_lst[idx] == 0
            new_lst[idx] = 1
            to_return += recursive_1_placer(_copy_shallow(new_lst), spacing, idx + spacing)
    return to_return


def enumerate_nucleosomal_states(amp_width, nuc_length, bin_size):
    """
    Puts nucleosome midpoints (1s) into an array of bins.
    Returns a DataFrame where each row is a valid nucleosome placement state.
    """
    nuc_bin_width = nuc_length // bin_size
    input_lst = [0 for _ in range(amp_width // bin_size)]
    nucs = recursive_1_placer(input_lst, nuc_bin_width, 0)
    return pd.DataFrame(nucs, columns=list(range(len(input_lst))))


def powerset_generator(n):
    """Returns a list of lists of all possible binary strings of length n."""
    if n == 0:
        return [[]]
    else:
        return [[0] + x for x in powerset_generator(n - 1)] + \
               [[1] + x for x in powerset_generator(n - 1)]


def enumerate_tf_states(tfbs_positions):
    """Returns a DataFrame of all possible TF binding configurations."""
    return pd.DataFrame(
        powerset_generator(len(tfbs_positions)),
        columns=list(range(len(tfbs_positions)))
    )


def construct_nuc_occupancy_matrix(nuc_states, amp_width, nuc_length, bin_size):
    """
    (n_states, amp_width) int8 matrix: 1 where any bound nuc covers that base.
    Processed bin-by-bin to stay within O(n_states * amp_width) memory.
    """
    n_states = len(nuc_states)
    mat = np.zeros((n_states, amp_width), dtype=np.int8)
    for k in range(len(nuc_states.columns)):
        start, end = get_nuc_start_and_end_from_midpoint(k, bin_size, nuc_length)
        start, end = max(0, start), min(amp_width, end)
        if end > start:
            occupied = nuc_states.values[:, k].astype(np.int8)
            mat[:, start:end] = np.maximum(mat[:, start:end], occupied[:, np.newaxis])
    return mat


def construct_tf_occupancy_matrix(tf_states, tfbs_positions, amp_width):
    """
    (n_states, amp_width) int8 matrix: 1 where any bound TF covers that base.
    """
    n_states = len(tf_states)
    mat = np.zeros((n_states, amp_width), dtype=np.int8)
    for k, tfbs in enumerate(tfbs_positions):
        start, end = max(0, int(tfbs[0])), min(amp_width, int(tfbs[1]))
        if end > start:
            occupied = tf_states.values[:, k].astype(np.int8)
            mat[:, start:end] = np.maximum(mat[:, start:end], occupied[:, np.newaxis])
    return mat


def prune_states(tf_states, nuc_states, tfbs_positions, amp_width, nuc_length, bin_size):
    """
    Returns paired (tf_states, nuc_states) DataFrames of valid microstates,
    where row i of both DataFrames jointly describes microstate i.
    Uses vectorized matmul overlap check.
    """
    nuc_mat = construct_nuc_occupancy_matrix(nuc_states, amp_width, nuc_length, bin_size)
    tf_mat = construct_tf_occupancy_matrix(tf_states, tfbs_positions, amp_width)
    # overlap[i,j] > 0 means tf state i and nuc state j conflict
    overlap = (tf_mat.astype(np.int32) @ nuc_mat.astype(np.int32).T) > 0
    valid = np.argwhere(~overlap)  # (n_valid, 2): columns are [tf_idx, nuc_idx]
    valid_tf = pd.DataFrame(
        tf_states.values[valid[:, 0]],
        columns=tf_states.columns,
        index=range(len(valid))
    )
    valid_nuc = pd.DataFrame(
        nuc_states.values[valid[:, 1]],
        columns=nuc_states.columns,
        index=range(len(valid))
    )
    return valid_tf, valid_nuc


def check_nuc_state_sanity(nuc_states, bin_size, nuc_length, n_sample=200):
    """Assert no two occupied nuc bins overlap in any sampled state."""
    import random
    indices = list(range(min(100, len(nuc_states))))
    if len(nuc_states) > 100:
        indices += random.sample(
            range(100, len(nuc_states)),
            min(100, len(nuc_states) - 100)
        )
    for idx in indices:
        state = nuc_states.iloc[idx]
        occupied = [k for k in state.index if state[k] > 0]
        windows = [get_nuc_start_and_end_from_midpoint(k, bin_size, nuc_length) for k in occupied]
        for i in range(len(windows)):
            for j in range(i + 1, len(windows)):
                s1, e1 = windows[i]
                s2, e2 = windows[j]
                assert e1 <= s2 or e2 <= s1, (
                    f"Nuc overlap in sampled state {idx}: bins {occupied[i]} and {occupied[j]} "
                    f"produce overlapping windows {windows[i]} and {windows[j]}"
                )


# ---------------------------------------------------------------------------
# Protection probability matrices
# ---------------------------------------------------------------------------

def build_nuc_protection_matrix(nuc_states, methyl_positions, params):
    """
    (n_states, n_gpcs) matrix of p(U=1 | nuc channel, state, j).
    p_nuc(j, bin_k) = sigmoid((d_edge - |pos_j - center_k|) / softness).
    Combined across occupied bins via max, so unoccupied bins contribute 0.
    Processed bin-by-bin to keep memory at O(n_states * n_gpcs).
    """
    methyl_arr = np.array(methyl_positions)
    n_states = len(nuc_states)
    n_gpcs = len(methyl_positions)
    protection = np.zeros((n_states, n_gpcs))
    for k in range(len(nuc_states.columns)):
        center = k * params.bin_size + params.bin_size // 2
        dist = np.abs(methyl_arr - center)
        prot_bin = expit((params.nuc_d_edge - dist) / params.nuc_softness)  # (n_gpcs,)
        occupied = nuc_states.values[:, k]  # (n_states,) 0 or 1
        contribution = occupied[:, np.newaxis] * prot_bin[np.newaxis, :]   # (n_states, n_gpcs)
        protection = np.maximum(protection, contribution)
    return protection


def build_tf_protection_matrix(tf_states, methyl_positions, tfbs_positions, prob_unmeth_given_tf):
    """
    (n_states, n_gpcs) matrix of p(U=1 | TF channel, state, j).
    Flat probability prob_unmeth_given_tf inside [tf_start-2, tf_end+2] per bound TF,
    combined across bound TFs via max.
    """
    methyl_arr = np.array(methyl_positions)
    n_states = len(tf_states)
    n_gpcs = len(methyl_positions)
    protection = np.zeros((n_states, n_gpcs))
    for k, tfbs in enumerate(tfbs_positions):
        start, end = int(tfbs[0]) - 2, int(tfbs[1]) + 2
        in_window = ((methyl_arr >= start) & (methyl_arr <= end)).astype(float)  # (n_gpcs,)
        prot_tf = in_window * prob_unmeth_given_tf
        occupied = tf_states.values[:, k]  # (n_states,)
        contribution = occupied[:, np.newaxis] * prot_tf[np.newaxis, :]
        protection = np.maximum(protection, contribution)
    return protection


def build_open_protection_array(methyl_positions, params):
    """
    (n_gpcs,) vector of p(U=1 | open channel, j).
    All positions get prob_unmeth_given_open; promoter window (if enabled) gets
    prob_unmeth_given_open_promoter. Promoter is disabled when promoter_lo == promoter_hi == 0.
    """
    arr = np.full(len(methyl_positions), params.prob_unmeth_given_open)
    if params.promoter_lo < params.promoter_hi:
        for i, p in enumerate(methyl_positions):
            if params.promoter_lo <= p <= params.promoter_hi:
                arr[i] = params.prob_unmeth_given_open_promoter
    return arr


def build_log_prob_matrices(nuc_states, tf_states, methyl_positions, tfbs_positions, params):
    """
    Build (n_states, n_gpcs) log p(T|state,j) and log p(C|state,j) matrices.
    p(T|state,j) = p_U(state,j) * p_T_given_unmeth + (1-p_U(state,j)) * p_T_given_meth
    where p_U = max(nuc_channel, tf_channel, open_channel).
    """
    nuc_prot = build_nuc_protection_matrix(nuc_states, methyl_positions, params)
    tf_prot = build_tf_protection_matrix(
        tf_states, methyl_positions, tfbs_positions, params.prob_unmeth_given_tf
    )
    open_prot = build_open_protection_array(methyl_positions, params)  # (n_gpcs,)

    p_unmeth = np.maximum(np.maximum(nuc_prot, tf_prot), open_prot[np.newaxis, :])
    p_t = p_unmeth * params.p_t_given_unmeth + (1.0 - p_unmeth) * params.p_t_given_meth
    p_c = 1.0 - p_t
    log_prob_t = np.log(np.clip(p_t, 1e-15, 1.0))
    log_prob_c = np.log(np.clip(p_c, 1e-15, 1.0))
    return log_prob_t, log_prob_c


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def determine_chunk_size(n_mols, n_states, memory_threshold_gb):
    """Minimum chunks so each chunk's (n_states x chunk_size x 8 bytes) fits in memory_threshold_gb."""
    return max(1, math.ceil(n_mols * n_states * 8 / (memory_threshold_gb * 1e9)))


def classify_all_molecules(single_molecules, log_prob_t, log_prob_c, n_chunks):
    """
    Chunked matmul: log_prob_t @ sub_mat + log_prob_c @ (1 - sub_mat) gives score
    matrix (n_states x chunk_size). Argmax per column = MLE state. Also return
    per-molecule log-likelihood = score of the assigned state.
    single_molecules: (n_gpcs, n_mols).
    """
    n_mols = single_molecules.shape[1]
    chunk_size = (n_mols + n_chunks - 1) // n_chunks
    assignments = []
    log_likelihoods = []
    for i in range(n_chunks):
        sub = single_molecules[:, i * chunk_size:(i + 1) * chunk_size]
        if sub.shape[1] == 0:
            continue
        scores = log_prob_t @ sub + log_prob_c @ (1.0 - sub)  # (n_states, chunk)
        best = np.argmax(scores, axis=0)
        lls = scores[best, np.arange(sub.shape[1])]
        assignments.extend(best.tolist())
        log_likelihoods.extend(lls.tolist())
    return np.array(assignments, dtype=np.int64), np.array(log_likelihoods)


# ---------------------------------------------------------------------------
# EM helper functions
# ---------------------------------------------------------------------------

def invert_t_fraction(t_frac, p_t_given_unmeth, p_t_given_meth):
    """
    Given observed mean T-fraction, recover p(unmeth) by inverting:
      E[T] = p_U * p_T_given_unmeth + (1-p_U) * p_T_given_meth
    Clamps result to [0, 1].
    """
    denom = p_t_given_unmeth - p_t_given_meth
    if abs(denom) < 1e-10:
        return None
    return float(np.clip((t_frac - p_t_given_meth) / denom, 0.0, 1.0))


def _build_tf_gpc_mask(methyl_positions, tfbs_positions):
    """(n_tfs, n_gpcs) bool array: True where GpC j is within [tf_start-2, tf_end+2]."""
    methyl_arr = np.array(methyl_positions)
    n_tfs = len(tfbs_positions)
    mask = np.zeros((n_tfs, len(methyl_positions)), dtype=bool)
    for k, tfbs in enumerate(tfbs_positions):
        mask[k] = (methyl_arr >= int(tfbs[0]) - 2) & (methyl_arr <= int(tfbs[1]) + 2)
    return mask


def _build_nuc_binary_coverage(nuc_states, methyl_positions, bin_size, nuc_length):
    """
    (n_states, n_gpcs) bool array: True where any bound nuc in that state covers GpC j
    (within the nuc's [center-nuc_length//2, center+nuc_length//2) window).
    """
    methyl_arr = np.array(methyl_positions)
    n_states, n_bins = nuc_states.values.shape
    n_gpcs = len(methyl_positions)
    coverage = np.zeros((n_states, n_gpcs), dtype=bool)
    for k in range(n_bins):
        start, end = get_nuc_start_and_end_from_midpoint(k, bin_size, nuc_length)
        in_window = (methyl_arr >= start) & (methyl_arr < end)   # (n_gpcs,)
        occupied = nuc_states.values[:, k].astype(bool)           # (n_states,)
        coverage |= occupied[:, np.newaxis] & in_window[np.newaxis, :]
    return coverage


def _build_distance_to_nearest_nuc(nuc_states, methyl_positions, bin_size):
    """
    (n_states, n_gpcs) float array: distance from GpC j to the closest occupied nuc
    center in each state. inf where no nuc is bound.
    Uses BenP's broadcasting trick: set unoccupied entries to inf, then take row-wise min.
    """
    methyl_arr = np.array(methyl_positions)
    n_states, n_bins = nuc_states.values.shape
    dist_to_nearest = np.full((n_states, len(methyl_positions)), np.inf)
    for k in range(n_bins):
        center = k * bin_size + bin_size // 2
        dist = np.abs(methyl_arr - center)  # (n_gpcs,)
        occupied = nuc_states.values[:, k].astype(bool)
        # contrib: dist where occupied, inf where not
        contrib = np.where(occupied[:, np.newaxis], dist[np.newaxis, :], np.inf)
        dist_to_nearest = np.minimum(dist_to_nearest, contrib)
    return dist_to_nearest


def run_em_mstep(single_molecules, assignments, tf_states, nuc_states,
                 methyl_positions, tfbs_positions, params,
                 nuc_distance_threshold=90, min_obs=100):
    """
    M-step: estimate prob_unmeth_given_tf, prob_unmeth_given_open, nuc_d_edge,
    nuc_softness from observed data at positions assigned to each protection class.
    All methylation probabilities are recovered by inverting through conversion params
    (NOT using raw T-fractions, which conflate conversion noise with methylation).
    Returns (estimates_dict, n_obs_dict) where estimates_dict may have None for any
    param where n_obs < min_obs.
    """
    methyl_arr = np.array(methyl_positions)
    n_mols = single_molecules.shape[1]  # single_molecules is (n_gpcs, n_mols)
    obs = single_molecules.T  # (n_mols, n_gpcs)

    tf_gpc_mask = _build_tf_gpc_mask(methyl_positions, tfbs_positions)       # (n_tfs, n_gpcs)
    nuc_coverage = _build_nuc_binary_coverage(
        nuc_states, methyl_positions, params.bin_size, params.nuc_length
    )   # (n_states, n_gpcs)
    dist_to_nuc = _build_distance_to_nearest_nuc(
        nuc_states, methyl_positions, params.bin_size
    )   # (n_states, n_gpcs)

    # Per-molecule masks
    tf_state_per_mol = tf_states.values[assignments, :]          # (n_mols, n_tfs)
    nuc_cov_per_mol = nuc_coverage[assignments, :]               # (n_mols, n_gpcs)
    dist_per_mol = dist_to_nuc[assignments, :]                   # (n_mols, n_gpcs)
    tf_cov_per_mol = (tf_state_per_mol @ tf_gpc_mask) > 0       # (n_mols, n_gpcs)

    # Promoter mask
    if params.promoter_lo < params.promoter_hi:
        in_promoter = (methyl_arr >= params.promoter_lo) & (methyl_arr <= params.promoter_hi)
    else:
        in_promoter = np.zeros(len(methyl_positions), dtype=bool)

    # --- estimate prob_unmeth_given_tf ---
    tf_mask = tf_cov_per_mol                                     # (n_mols, n_gpcs)
    tf_obs = obs[tf_mask]
    n_tf = len(tf_obs)
    est_tf = None
    if n_tf >= min_obs:
        est_tf = invert_t_fraction(tf_obs.mean(), params.p_t_given_unmeth, params.p_t_given_meth)

    # --- estimate prob_unmeth_given_open ---
    open_mask = (~nuc_cov_per_mol & ~tf_cov_per_mol
                 & ~in_promoter[np.newaxis, :])                  # (n_mols, n_gpcs)
    open_obs = obs[open_mask]
    n_open = len(open_obs)
    est_open = None
    if n_open >= min_obs:
        est_open = invert_t_fraction(
            open_obs.mean(), params.p_t_given_unmeth, params.p_t_given_meth
        )

    # --- fit nuc_d_edge and nuc_softness ---
    nuc_fit_mask = (dist_per_mol < nuc_distance_threshold) & ~tf_cov_per_mol
    dist_flat = dist_per_mol[nuc_fit_mask]
    nuc_obs_flat = obs[nuc_fit_mask]
    n_nuc = len(nuc_obs_flat)
    est_d_edge, est_softness = None, None
    if n_nuc >= min_obs:
        p_tu = params.p_t_given_unmeth
        p_tm = params.p_t_given_meth

        def nll(pv):
            d, s = pv
            p_u = expit((d - dist_flat) / s)
            p_t = p_u * p_tu + (1.0 - p_u) * p_tm
            p_t = np.clip(p_t, 1e-10, 1.0 - 1e-10)
            return -np.sum(
                nuc_obs_flat * np.log(p_t) + (1.0 - nuc_obs_flat) * np.log(1.0 - p_t)
            )

        res = minimize(
            nll,
            [params.nuc_d_edge, params.nuc_softness],
            method='L-BFGS-B',
            bounds=[(20.0, 90.0), (1.0, 30.0)]
        )
        if res.success or res.fun < nll([params.nuc_d_edge, params.nuc_softness]):
            est_d_edge, est_softness = float(res.x[0]), float(res.x[1])

    estimates = {
        'prob_unmeth_given_tf': est_tf,
        'prob_unmeth_given_open': est_open,
        'nuc_d_edge': est_d_edge,
        'nuc_softness': est_softness,
    }
    n_obs = {'n_obs_tf': n_tf, 'n_obs_open': n_open, 'n_obs_nuc': n_nuc}
    return estimates, n_obs


# ---------------------------------------------------------------------------
# Output functions
# ---------------------------------------------------------------------------

def expand_state(tf_state, nuc_state, state_idx, tfbs_positions, bin_size, nuc_length):
    """
    Convert a (tf_state, nuc_state) pair into a flat dict for DataFrame construction.
    state_idx is preserved as 'idx'. Loop variables use bin_k / tf_k so they do NOT
    shadow the state_idx parameter.
    """
    expanded = {}
    for bin_k in nuc_state.index:
        expanded[f'nuc{bin_k}_present'] = bool(nuc_state[bin_k] == 1)
        bounds = get_nuc_start_and_end_from_midpoint(bin_k, bin_size, nuc_length)
        expanded[f'nuc{bin_k}_start'] = bounds[0]
        expanded[f'nuc{bin_k}_end'] = bounds[1]
    for tf_k in tf_state.index:
        expanded[f'tfbs_{tf_k + 1}'] = bool(tf_state[tf_k] == 1)
    expanded['idx'] = state_idx
    return expanded


def summarize_assignments(tf_states, nuc_states, assignments, log_likelihoods,
                          tfbs_positions, bin_size, nuc_length):
    rows = []
    for i, a in enumerate(assignments):
        row = expand_state(
            tf_states.iloc[a], nuc_states.iloc[a], int(a),
            tfbs_positions, bin_size, nuc_length
        )
        row['log_likelihood'] = float(log_likelihoods[i])
        rows.append(row)
    return pd.DataFrame(rows)


def write_valid_states(tf_states, nuc_states, out_path, tfbs_positions, bin_size, nuc_length):
    """Write valid state enumeration TSV."""
    rows = []
    for i in range(len(tf_states)):
        rows.append(
            expand_state(tf_states.iloc[i], nuc_states.iloc[i], i,
                         tfbs_positions, bin_size, nuc_length)
        )
    pd.DataFrame(rows).to_csv(out_path, sep='\t', header=True, index=False)


# ---------------------------------------------------------------------------
# Top-level compute_classifications
# ---------------------------------------------------------------------------

def compute_classifications(single_molecules, methyl_positions, tfbs_positions, params,
                            valid_states_path=None, precomputed_states_file=None,
                            memory_threshold=1.0, do_em=False, em_max_iters=5,
                            em_log_path=None, nuc_distance_threshold=90, min_obs=100):
    """
    Main classification entry point.
    single_molecules: (n_gpcs, n_mols) numpy array with values in {0, 1} (0.5 for missing).
    methyl_positions: list of int GpC positions, length n_gpcs.
    tfbs_positions: list of [start, end, name] from load_tfbs_positions.
    params: ModelParams instance.
    Returns DataFrame with one row per molecule: nuc/TF state columns + log_likelihood.
    """
    # --- State enumeration ---
    if precomputed_states_file and os.path.exists(precomputed_states_file):
        print('loading precomputed states')
        with open(precomputed_states_file, 'rb') as fh:
            tf_states, nuc_states = pickle.load(fh)
    else:
        print('enumerating states')
        raw_tf = enumerate_tf_states(tfbs_positions)
        raw_nuc = enumerate_nucleosomal_states(params.amp_width, params.nuc_length, params.bin_size)
        print(f'  raw tf states: {len(raw_tf)}, raw nuc states: {len(raw_nuc)}')
        tf_states, nuc_states = prune_states(
            raw_tf, raw_nuc, tfbs_positions,
            params.amp_width, params.nuc_length, params.bin_size
        )
        if precomputed_states_file:
            with open(precomputed_states_file, 'wb') as fh:
                pickle.dump((tf_states, nuc_states), fh)
            print(f'  saved states to {precomputed_states_file}')

    print(f'valid states: {len(tf_states)}')
    check_nuc_state_sanity(nuc_states, params.bin_size, params.nuc_length)

    n_states = len(tf_states)
    n_mols = single_molecules.shape[1]
    n_chunks = determine_chunk_size(n_mols, n_states, memory_threshold)
    print(f'n_chunks: {n_chunks}')

    # --- EM or single-pass ---
    em_log_rows = []

    if do_em:
        print('starting EM')
        current_params = params  # will be replaced each iter via dataclass copy
        converged = False
        for em_iter in range(em_max_iters):
            log_prob_t, log_prob_c = build_log_prob_matrices(
                nuc_states, tf_states, methyl_positions, tfbs_positions, current_params
            )
            assignments, log_likelihoods = classify_all_molecules(
                single_molecules, log_prob_t, log_prob_c, n_chunks
            )

            estimates, n_obs = run_em_mstep(
                single_molecules, assignments, tf_states, nuc_states,
                methyl_positions, tfbs_positions, current_params,
                nuc_distance_threshold=nuc_distance_threshold, min_obs=min_obs
            )

            # Apply damping and compute fractional changes; warn on insufficient data
            param_names = ['prob_unmeth_given_tf', 'prob_unmeth_given_open',
                           'nuc_d_edge', 'nuc_softness']
            old_vals = {p: getattr(current_params, p) for p in param_names}
            new_vals = {}
            frac_changes = []
            for pname in param_names:
                est = estimates[pname]
                old = old_vals[pname]
                if est is None:
                    obs_key = {
                        'prob_unmeth_given_tf': 'n_obs_tf',
                        'prob_unmeth_given_open': 'n_obs_open',
                        'nuc_d_edge': 'n_obs_nuc',
                        'nuc_softness': 'n_obs_nuc',
                    }[pname]
                    print(f'  WARNING: {pname} not updated ({n_obs[obs_key]} obs < {min_obs})')
                    new_vals[pname] = old
                    frac_changes.append(0.0)
                else:
                    damped = 0.5 * old + 0.5 * est
                    new_vals[pname] = damped
                    frac_changes.append(abs(damped - old) / max(abs(old), 1e-10))

            max_frac = max(frac_changes)
            converged = max_frac < 0.05

            em_log_rows.append({
                'iter': em_iter,
                'prob_unmeth_tf': new_vals['prob_unmeth_given_tf'],
                'prob_unmeth_open': new_vals['prob_unmeth_given_open'],
                'nuc_d_edge': new_vals['nuc_d_edge'],
                'nuc_softness': new_vals['nuc_softness'],
                'n_obs_tf': n_obs['n_obs_tf'],
                'n_obs_open': n_obs['n_obs_open'],
                'n_obs_nuc': n_obs['n_obs_nuc'],
                'max_frac_change': max_frac,
                'converged': converged,
            })
            print(f'  iter {em_iter}: max_frac_change={max_frac:.4f} converged={converged}')

            # Rebuild params with updated values (dataclass copy)
            import copy as _copy
            current_params = _copy.copy(current_params)
            for pname, val in new_vals.items():
                setattr(current_params, pname, val)

            if converged:
                print('  EM converged')
                break

        # Final pass with converged params
        log_prob_t, log_prob_c = build_log_prob_matrices(
            nuc_states, tf_states, methyl_positions, tfbs_positions, current_params
        )
        assignments, log_likelihoods = classify_all_molecules(
            single_molecules, log_prob_t, log_prob_c, n_chunks
        )

    else:
        print('classification (no EM)')
        log_prob_t, log_prob_c = build_log_prob_matrices(
            nuc_states, tf_states, methyl_positions, tfbs_positions, params
        )
        assignments, log_likelihoods = classify_all_molecules(
            single_molecules, log_prob_t, log_prob_c, n_chunks
        )

    # Write EM log (never to cwd; only if path provided)
    if em_log_rows and em_log_path:
        pd.DataFrame(em_log_rows).to_csv(em_log_path, sep='\t', index=False)

    # Write valid states if requested
    if valid_states_path:
        write_valid_states(
            tf_states, nuc_states, valid_states_path,
            tfbs_positions, params.bin_size, params.nuc_length
        )

    print('summarizing output')
    return summarize_assignments(
        tf_states, nuc_states, assignments, log_likelihoods,
        tfbs_positions, params.bin_size, params.nuc_length
    )


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_yaml_config(path):
    try:
        import yaml as _yaml
    except ImportError:
        raise ImportError(
            "PyYAML is required to use --config. Install with: pip install pyyaml"
        )
    with open(path) as fh:
        return _yaml.safe_load(fh) or {}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    # Two-pass parse: first extract --config, load YAML defaults, then full parse
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument('--config', default=None)
    pre_args, _ = pre.parse_known_args()
    config_defaults = load_yaml_config(pre_args.config) if pre_args.config else {}

    parser = argparse.ArgumentParser(
        description='Single-molecule binding classification v4')
    parser.set_defaults(**config_defaults)

    # I/O
    parser.add_argument('--input', dest='input', type=str, required=True)
    parser.add_argument('--output', dest='output', type=str, required=True)
    parser.add_argument('--all_states_output', default=None)
    parser.add_argument('--precomputed_states_file', default=None)
    parser.add_argument('--plot', dest='plot_file', type=str, required=True)
    parser.add_argument('--positions', dest='positions_file', type=str, required=True)
    parser.add_argument('--amplicon_name', type=str, required=True)
    parser.add_argument('--config', default=None, help='Path to YAML config file')

    # Preprocessing
    parser.add_argument('--reads_to_use', type=int, default=0,
                        help='Subsample to this many reads after loading (0=use all)')
    parser.add_argument('--filter_threshold', type=float, default=1.0)
    parser.add_argument('--convert_ambiguous_gcgs', type=str, default='',
                        help='Comma-separated pairs of positions for ambiguous GCG correction')

    # Amplicon constants
    parser.add_argument('--amp_width', type=int, default=None,
                        help='Amplicon width in bp (default: auto-detect from methyl_positions)')
    parser.add_argument('--nuc_length', type=int, default=140)
    parser.add_argument('--bin_size', type=int, default=10)

    # Conversion params
    parser.add_argument('--p_t_given_unmeth', type=float, default=0.95)
    parser.add_argument('--p_t_given_meth', type=float, default=0.15)

    # Protection params
    parser.add_argument('--nuc_d_edge', type=float, default=65.0)
    parser.add_argument('--nuc_softness', type=float, default=5.0)
    parser.add_argument('--prob_unmeth_given_tf', type=float, default=0.9)
    parser.add_argument('--prob_unmeth_given_open', type=float, default=0.05)
    parser.add_argument('--promoter_positions', type=str, default='0,0',
                        help='lo,hi for promoter window (0,0 = disabled)')
    parser.add_argument('--prob_unmeth_given_open_promoter', type=float, default=0.5)

    # Plotting
    parser.add_argument('--reads_to_plot', type=int, default=10)
    parser.add_argument('--individual_reads_to_plot', type=str, default='')

    # Memory
    parser.add_argument('--memory_threshold', type=float, default=1.0)

    # EM
    parser.add_argument('--do_em', action='store_true', default=False)
    parser.add_argument('--em_max_iters', type=int, default=5)
    parser.add_argument('--em_log', dest='em_log', type=str, default=None,
                        help='Path for EM convergence log TSV (default: <output>.em_log.tsv)')
    parser.add_argument('--nuc_distance_threshold', type=float, default=90.0)
    parser.add_argument('--em_min_obs', type=int, default=100)

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        print('loading data')
        mat = load_single_molecule_matrix(args.input)
        methyl_positions = get_methyl_positions(mat)

        pos_file = load_tfbs_positions(args.positions_file)
        tfbs_positions = pos_file[args.amplicon_name]

        mat = filter_all_converted_reads(mat, args.filter_threshold)
        mat = fix_missing_data(mat)

        # B1 fix: actually apply --reads_to_use
        if args.reads_to_use > 0 and len(mat) > args.reads_to_use:
            mat = mat.sample(args.reads_to_use, random_state=42)
            print(f'subsampled to {args.reads_to_use} reads')

        if args.convert_ambiguous_gcgs:
            positions_to_fix = list(map(int, args.convert_ambiguous_gcgs.split(',')))
            assert len(positions_to_fix) % 2 == 0
            for i in range(len(positions_to_fix) // 2):
                mat = adjust_gcgs(mat, positions_to_fix[2 * i], positions_to_fix[2 * i + 1])

        print(f'num molecules: {len(mat)}')

        # Parse promoter positions (B6 fix: compare [0,0] list not (0,0) tuple)
        prom = list(map(int, args.promoter_positions.split(',')))
        promoter_lo, promoter_hi = prom[0], prom[1]

        # Determine amp_width
        amp_width = args.amp_width
        if amp_width is None:
            amp_width = auto_detect_amp_width(methyl_positions, args.nuc_length)
            print(f'auto-detected amp_width: {amp_width}')

        params = ModelParams(
            amp_width=amp_width,
            nuc_length=args.nuc_length,
            bin_size=args.bin_size,
            p_t_given_unmeth=args.p_t_given_unmeth,
            p_t_given_meth=args.p_t_given_meth,
            nuc_d_edge=args.nuc_d_edge,
            nuc_softness=args.nuc_softness,
            prob_unmeth_given_tf=args.prob_unmeth_given_tf,
            prob_unmeth_given_open=args.prob_unmeth_given_open,
            promoter_lo=promoter_lo,
            promoter_hi=promoter_hi,
            prob_unmeth_given_open_promoter=args.prob_unmeth_given_open_promoter,
        )

        single_molecules = mat.T.values

        em_log_path = args.em_log if args.em_log else (args.output + '.em_log.tsv')

        print('starting classification')
        assignments = compute_classifications(
            single_molecules, methyl_positions, tfbs_positions, params,
            valid_states_path=args.all_states_output,
            precomputed_states_file=args.precomputed_states_file,
            memory_threshold=args.memory_threshold,
            do_em=args.do_em,
            em_max_iters=args.em_max_iters,
            em_log_path=em_log_path if args.do_em else None,
            nuc_distance_threshold=args.nuc_distance_threshold,
            min_obs=args.em_min_obs,
        )
        assignments.index = mat.index
        assignments.to_csv(args.output, sep='\t', header=True, index=True)

        if args.reads_to_plot:
            print('plotting')
            n_to_plot = min(args.reads_to_plot, len(mat))
            to_plot = assignments.sample(n_to_plot)
            for idx, assignment in to_plot.iterrows():
                fig, ax = plt.subplots(figsize=(12, 8))
                plot_single_read(mat, idx, ax, fillbetween=tfbs_positions)
                decorate_single_read_plot2(assignment, ax, tfbs_positions)
                plots.savefig()
                plt.close()

            if args.individual_reads_to_plot:
                for idx in args.individual_reads_to_plot.split(','):
                    if idx in mat.index.tolist():
                        fig, ax = plt.subplots(figsize=(12, 8))
                        plot_single_read(mat, idx, ax, fillbetween=tfbs_positions)
                        decorate_single_read_plot2(assignments.loc[idx], ax, tfbs_positions)
                        plots.savefig()
                        plt.close()
