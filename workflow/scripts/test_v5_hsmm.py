#!/usr/bin/env python
"""
Synthetic unit tests for classify_single_molecule_binding_v5_hsmm.py (SPEC section 7).

Each molecule has a KNOWN, clean (noise-free) footprint pattern; we assert the Viterbi
segmentation recovers the intended state sequence. Run from workflow/scripts/:

    python test_v5_hsmm.py
"""
import numpy as np

import classify_single_molecule_binding_v5_hsmm as v5

# amplicon: GpCs every 10 bp, motifs mimic opJS4 6x TetO (280-302, 320-342, ...)
GPC = np.arange(0, 601, 10)
MOTIFS = [[280, 302, 'TetO'], [320, 342, 'TetO'], [360, 382, 'TetO'],
          [400, 422, 'TetO'], [440, 462, 'TetO'], [480, 502, 'TetO']]


def make_obs(protected_intervals, missing=()):
    """Build obs1/obs0 (n_gpc, 1) for a single molecule; protected inside intervals, else accessible."""
    obs = np.zeros(len(GPC), dtype=int)          # 0 accessible
    for (lo, hi) in protected_intervals:
        obs[(GPC >= lo) & (GPC < hi)] = 1        # 1 protected
    for (lo, hi) in missing:
        obs[(GPC >= lo) & (GPC < hi)] = -1
    obs1 = (obs == 1).astype(float)[:, None]
    obs0 = (obs == 0).astype(float)[:, None]
    return obs1, obs0


def decode_one(protected_intervals, missing=(), **overrides):
    params = v5.ModelParams(amp_width=610, grid_step=5)
    for k, val in overrides.items():
        setattr(params, k, val)
    obs1, obs0 = make_obs(protected_intervals, missing)
    paths, ll = v5.viterbi_decode(GPC, obs1, obs0, MOTIFS, params)
    return paths[0], ll[0]


def counts(path):
    c = {'OPEN': 0, 'NUC': 0, 'TF': 0, 'UNID': 0}
    tf_motifs = []
    for (s, a, b, anchor) in path:
        c[v5.STATE_NAMES[s]] += 1
        if s == v5.TF:
            tf_motifs.append(int(anchor))
    return c, sorted(tf_motifs)


def show(name, path):
    segs = ' '.join('{}[{}-{}{}]'.format(
        v5.STATE_NAMES[s], a, b, ('@dyad{:.0f}'.format(anchor) if s == v5.NUC else
                                  ('@motif{}'.format(int(anchor)) if s == v5.TF else '')))
        for (s, a, b, anchor) in path)
    print('  {}: {}'.format(name, segs))


def test_lone_tf():
    path, _ = decode_one([(278, 304)])              # single motif-0 blob
    show('lone_tf', path)
    c, motifs = counts(path)
    assert c['TF'] == 1, c
    assert motifs == [0], motifs
    assert c['NUC'] == 0 and c['UNID'] == 0, c


def test_two_adjacent_tfs():
    path, _ = decode_one([(278, 304), (318, 344)])  # motif 0 and 1 protected, linker accessible
    show('two_adjacent_tfs', path)
    c, motifs = counts(path)
    assert motifs == [0, 1], motifs
    assert c['NUC'] == 0 and c['UNID'] == 0, c


def test_single_nucleosome():
    path, _ = decode_one([(85, 216)])               # ~130 bp protected, off-motif region
    show('single_nucleosome', path)
    c, motifs = counts(path)
    assert c['NUC'] == 1, c
    assert c['TF'] == 0 and c['UNID'] == 0, c
    # dyad should land near the center (~150)
    dyad = [anchor for (s, a, b, anchor) in path if s == v5.NUC][0]
    assert 120 <= dyad <= 180, dyad


def test_nucleosome_plus_tf():
    path, _ = decode_one([(85, 216), (278, 304)])   # nuc then a TF at motif 0
    show('nucleosome_plus_tf', path)
    c, motifs = counts(path)
    assert c['NUC'] == 1, c
    assert motifs == [0], motifs


def test_unid_offmotif_short():
    path, _ = decode_one([(125, 176)])              # ~50 bp protected, off-motif, too short for a nuc
    show('unid_offmotif_short', path)
    c, motifs = counts(path)
    assert c['UNID'] == 1, c
    assert c['TF'] == 0 and c['NUC'] == 0, c


def test_di_nucleosome():
    path, _ = decode_one([(60, 340)])               # ~280 bp continuous protection => two nucs
    show('di_nucleosome', path)
    c, motifs = counts(path)
    # should be explained by two nucleosomes (di-nuc), not a single over-long nuc or spurious TFs
    assert c['NUC'] == 2, c


def test_all_open():
    path, _ = decode_one([])                         # nothing protected
    show('all_open', path)
    c, motifs = counts(path)
    assert c['NUC'] == 0 and c['TF'] == 0 and c['UNID'] == 0, c


def test_missing_data_skipped():
    # a lone TF at motif 0 (flanks still accessible so the call is well-posed), plus a missing
    # patch off in an open region. Missing GpCs must be skipped => no phantom footprint there,
    # and the TF is still recovered.
    path, _ = decode_one([(278, 304)], missing=[(120, 200)])
    show('missing_data_skipped', path)
    c, motifs = counts(path)
    assert motifs == [0], motifs
    assert c['NUC'] == 0 and c['UNID'] == 0, c


if __name__ == '__main__':
    tests = [test_lone_tf, test_two_adjacent_tfs, test_single_nucleosome,
             test_nucleosome_plus_tf, test_unid_offmotif_short, test_di_nucleosome,
             test_all_open, test_missing_data_skipped]
    failed = 0
    for t in tests:
        print('== {} =='.format(t.__name__))
        try:
            t()
            print('  PASS')
        except AssertionError as e:
            failed += 1
            print('  FAIL: {}'.format(e))
    print('\n{}/{} passed'.format(len(tests) - failed, len(tests)))
    raise SystemExit(1 if failed else 0)
