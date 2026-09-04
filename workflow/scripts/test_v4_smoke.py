#!/usr/bin/env python3
"""
test_v4_smoke.py — smoke test for classify_single_molecule_binding_v4.py
Do not commit this file.

Constructs synthetic data with two clearly separable molecule classes,
runs compute_classifications both with and without EM, and asserts
basic correctness without validating exact classification accuracy.
"""

import sys
import os
import tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import classify_single_molecule_binding_v4 as v4


def make_synthetic_data():
    """
    30 molecules x 25 GpCs on a ~330bp amplicon.

    Class A (molecules 0-14, "TF"): all 4 TFs bound, no nuc.
      Protected (1) at the 8 GpC positions that fall inside the 4 TFBS windows.
      Accessible (0) everywhere else.

    Class B (molecules 15-29, "nuc"): one nuc centered around position 215,
      covering GpCs 150-280. Open at TF region.
      Protected (1) at GpCs in [145, 285).
      Accessible (0) everywhere else.

    Returns: (single_molecules: (n_gpcs, n_mols), methyl_positions, tfbs_positions)
    """
    methyl_positions = list(range(10, 260, 10))   # [10, 20, ..., 250], 25 positions
    assert len(methyl_positions) == 25

    # tfbs_positions format from common.load_tfbs_positions: [[start, end, name], ...]
    tfbs_positions = [
        [50,  60,  'TetO1'],
        [80,  90,  'TetO2'],
        [110, 120, 'TetO3'],
        [140, 150, 'TetO4'],
    ]

    # GpCs covered by TFs (within [start-2, end+2])
    tf_covered = set()
    for tfbs in tfbs_positions:
        for p in methyl_positions:
            if tfbs[0] - 2 <= p <= tfbs[1] + 2:
                tf_covered.add(p)

    # GpCs covered by nuc (center=215, nuc_length=140 -> [145, 285))
    nuc_center = 215
    nuc_half = 70
    nuc_covered = {p for p in methyl_positions
                   if nuc_center - nuc_half <= p < nuc_center + nuc_half}

    n_gpcs = len(methyl_positions)
    mat = np.zeros((30, n_gpcs), dtype=float)

    for i in range(15):           # TF molecules
        for j, p in enumerate(methyl_positions):
            mat[i, j] = 1.0 if p in tf_covered else 0.0

    for i in range(15, 30):       # nuc molecules
        for j, p in enumerate(methyl_positions):
            mat[i, j] = 1.0 if p in nuc_covered else 0.0

    return mat.T, methyl_positions, tfbs_positions   # (n_gpcs, n_mols)


def test_no_em():
    print("=== Test 1: no EM ===")
    single_molecules, methyl_positions, tfbs_positions = make_synthetic_data()

    amp_width = v4.auto_detect_amp_width(methyl_positions, nuc_length=140)
    print(f"  amp_width: {amp_width}")
    params = v4.ModelParams(amp_width=amp_width)

    result = v4.compute_classifications(
        single_molecules, methyl_positions, tfbs_positions, params,
        memory_threshold=1.0, do_em=False,
    )

    assert len(result) == 30, f"Expected 30 rows, got {len(result)}"
    assert 'log_likelihood' in result.columns, "Missing log_likelihood column"
    assert result['log_likelihood'].notna().all(), "NaN log_likelihoods found"

    tf_cols  = [f'tfbs_{k}' for k in range(4) if f'tfbs_{k}' in result.columns]
    nuc_cols = [c for c in result.columns if c.endswith('_present')]

    tf_correct = 0
    for i in range(15):
        n_bound = sum(int(result.iloc[i][c]) for c in tf_cols)
        if n_bound >= 3:     # at least 3 of 4 TFs bound
            tf_correct += 1

    nuc_correct = 0
    for i in range(15, 30):
        n_nucs = sum(int(result.iloc[i][c]) for c in nuc_cols)
        if n_nucs >= 1:      # at least one nuc present
            nuc_correct += 1

    total_correct = tf_correct + nuc_correct
    print(f"  TF molecules correctly classified:  {tf_correct}/15")
    print(f"  Nuc molecules correctly classified: {nuc_correct}/15")
    print(f"  Total: {total_correct}/30 (allowing up to 5 miscalls)")

    assert total_correct >= 25, (
        f"Too many miscalls: {30 - total_correct}/30 wrong (>5 allowed). "
        "Check protection parameters or synthetic data construction."
    )
    print("  PASS")
    return result


def test_em():
    print("=== Test 2: with EM ===")
    single_molecules, methyl_positions, tfbs_positions = make_synthetic_data()

    amp_width = v4.auto_detect_amp_width(methyl_positions, nuc_length=140)
    params = v4.ModelParams(amp_width=amp_width)

    em_log_fd, em_log_path = tempfile.mkstemp(suffix='.em_log.tsv')
    os.close(em_log_fd)

    try:
        result = v4.compute_classifications(
            single_molecules, methyl_positions, tfbs_positions, params,
            memory_threshold=1.0,
            do_em=True,
            em_max_iters=3,
            em_log_path=em_log_path,
            min_obs=5,    # lower than default 100 so tiny test data isn't all warnings
        )

        assert len(result) == 30, f"Expected 30 rows, got {len(result)}"
        assert 'log_likelihood' in result.columns

        log_df = pd.read_csv(em_log_path, sep='\t')
        assert len(log_df) >= 1, f"EM log has {len(log_df)} rows, expected >= 1"
        expected_cols = {'iter', 'prob_unmeth_tf', 'prob_unmeth_open',
                         'nuc_d_edge', 'nuc_softness',
                         'n_obs_tf', 'n_obs_open', 'n_obs_nuc',
                         'max_frac_change', 'converged'}
        missing = expected_cols - set(log_df.columns)
        assert not missing, f"EM log missing columns: {missing}"

        print(f"  EM ran {len(log_df)} iteration(s)")
        print(f"  EM log columns present: OK")
        print("  PASS")
    finally:
        if os.path.exists(em_log_path):
            os.unlink(em_log_path)

    return result


if __name__ == '__main__':
    test_no_em()
    test_em()
    print("\nAll smoke tests PASSED")
