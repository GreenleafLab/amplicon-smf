#!/usr/bin/env python
"""
inspect_molecule_v5.py -- single-molecule diagnostic / calibration tool for the v5 HSMM classifier.

Decode one (or a few) molecules, print the raw GpC data and the Viterbi path WITH a per-segment
score breakdown, and optionally score a hypothesized "correct" segmentation to see exactly which
term (emission vs start/duration/transition cost) drives the model's choice. This is the tool used
to work through edge cases in NOTES_v5_calibration_log.md.

Examples (run from workflow/scripts/, or it adds itself to sys.path):

  # inspect one molecule with default params
  inspect_molecule_v5.py --input <matrix> --positions <pos.txt> --amplicon_name opJS4_7x_TetO_21bp_no_CG \
      --read_id M00653:...:1109:15970:16317

  # compare the model's call to "all TFs bound" (TF@k = TF at motif k), and try a lower TF cost
  inspect_molecule_v5.py ... --read_id <id> --start_cost_tf 1.0 \
      --compare "TF@0;TF@1;TF@2;TF@3;TF@4;TF@5;TF@6"

  # compare to an explicit segmentation (OPEN gaps are auto-filled)
  inspect_molecule_v5.py ... --compare "NUC:135-277;TF:278-304;TF:318-344"
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_single_molecule_binding_v5_hsmm as v5
from common import load_single_molecule_matrix, get_methyl_positions, load_tfbs_positions


def parse_compare(spec, tfbs_positions, tf_margin):
    """Parse a footprint spec into [(state, a, b), ...]. Accepts 'TYPE:start-end' or 'TF@k'."""
    fps = []
    for tok in spec.split(';'):
        tok = tok.strip()
        if not tok:
            continue
        if tok.upper().startswith('TF@'):
            k = int(tok[3:])
            a, b = int(tfbs_positions[k][0]) - tf_margin, int(tfbs_positions[k][1]) + tf_margin
            fps.append((v5.TF, a, b))
            continue
        typ, span = tok.split(':')
        a, b = span.split('-')
        state = {'OPEN': v5.OPEN, 'NUC': v5.NUC, 'TF': v5.TF, 'UNID': v5.UNID}[typ.strip().upper()]
        fps.append((state, int(a), int(b)))
    return fps


def print_path(title, total, rows):
    print('\n=== {}  (total log-lik = {:.2f}) ==='.format(title, total))
    print('  {:5s} {:>4s}-{:<4s} {:>8s} {:>6s} {:>7s} {:>6s} {:>8s}  {}'.format(
        'type', 'a', 'b', 'emission', 'start', 'dur', 'trans', 'subtotal', 'note'))
    for r in rows:
        print('  {:5s} {:>4d}-{:<4d} {:8.2f} {:6.2f} {:7.2f} {:6.2f} {:8.2f}  {}'.format(
            r['type'], r['start'], r['end'], r['emission'], r['start_cost'],
            r['dur_logprob'], r['trans_from_prev'], r['subtotal'], r['note']))


def main():
    p = argparse.ArgumentParser(description='v5 HSMM single-molecule diagnostic')
    p.add_argument('--input', required=True)
    p.add_argument('--positions', required=True)
    p.add_argument('--amplicon_name', required=True)
    p.add_argument('--read_id', required=True, help='comma-separated read_ids')
    p.add_argument('--region', default='', help='lo,hi bp range to print data for (default: all)')
    p.add_argument('--compare', default='', help='hypothesized footprints; see module docstring')
    p.add_argument('--config', default=None, help='YAML of model params (defaults; flags override)')
    # the knobs most worth sweeping; everything else uses ModelParams defaults / --config
    p.add_argument('--grid_step', type=int, default=5)
    p.add_argument('--p_t_given_unmeth', type=float, default=0.99)
    p.add_argument('--p_t_given_meth', type=float, default=0.05)
    p.add_argument('--start_cost_tf', type=float, default=None)
    p.add_argument('--start_cost_nuc', type=float, default=None)
    p.add_argument('--start_cost_unid', type=float, default=None)
    p.add_argument('--trans_nuc_tf', type=float, default=None)
    p.add_argument('--nuc_softness', type=float, default=None)
    p.add_argument('--promoter_positions', default='0,0')
    args = p.parse_args()

    cfg = {}
    if args.config:
        import yaml
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh) or {}
    prom = list(map(int, args.promoter_positions.split(',')))
    overrides = dict(grid_step=args.grid_step, p_t_given_unmeth=args.p_t_given_unmeth,
                     p_t_given_meth=args.p_t_given_meth, promoter_lo=prom[0], promoter_hi=prom[1])
    for k in ('start_cost_tf', 'start_cost_nuc', 'start_cost_unid', 'trans_nuc_tf',
              'nuc_softness'):
        v = getattr(args, k)
        if v is not None:
            overrides[k] = v
    cfg.update(overrides)

    mat = load_single_molecule_matrix(args.input)
    gpc = np.array(get_methyl_positions(mat), dtype=np.int64)
    tfbs = load_tfbs_positions(args.positions).get(args.amplicon_name, [])
    amp_width = int(gpc.max()) + 115
    params = v5.ModelParams(amp_width=amp_width, **cfg)

    region = None
    if args.region:
        lo, hi = map(int, args.region.split(','))
        region = (lo, hi)

    for rid in args.read_id.split(','):
        rid = rid.strip()
        if rid not in mat.index:
            print('!! read_id not found: {}'.format(rid))
            continue
        obs = mat.loc[[rid]].values.T.astype(float)
        obs1 = (obs == 1).astype(float)
        obs0 = (obs == 0).astype(float)

        print('\n' + '=' * 70)
        print('read: {}'.format(rid))
        data = []
        for pos, o in zip(gpc, obs[:, 0]):
            if region and not (region[0] <= pos <= region[1]):
                continue
            data.append('{}{}'.format(pos, {1: 'P', 0: 'a', -1: '.'}.get(int(o), '?')))
        print('data{}: {}'.format('' if region is None else ' [{}-{}]'.format(*region), ' '.join(data)))

        paths, ll = v5.viterbi_decode(gpc, obs1, obs0, tfbs, params)
        best = paths[0]
        tot, rows = v5.score_segmentation(best, gpc, obs1[:, 0], obs0[:, 0], tfbs, params)
        tf_called = sorted(int(a) for (s, _, _, a) in best if s == v5.TF)
        print_path('DECODED (model best)', tot, rows)
        print('  TF motifs called: {}   n_nuc={}   n_unid={}'.format(
            tf_called, sum(1 for (s, *_r) in best if s == v5.NUC),
            sum(1 for (s, *_r) in best if s == v5.UNID)))

        if args.compare:
            fps = parse_compare(args.compare, tfbs, params.tf_margin)
            alt = v5.fill_open_gaps(fps, amp_width, gpc, tfbs, params)
            tot2, rows2 = v5.score_segmentation(alt, gpc, obs1[:, 0], obs0[:, 0], tfbs, params)
            print_path('COMPARE (your hypothesis)', tot2, rows2)
            print('\n  DELTA (compare - decoded) = {:.2f}  =>  {}'.format(
                tot2 - tot, 'your hypothesis WINS' if tot2 > tot else 'model best wins'))


if __name__ == '__main__':
    main()
