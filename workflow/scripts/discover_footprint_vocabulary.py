#!/usr/bin/env python
"""
Build a POOLED per-promoter footprint VOCABULARY for the v5 HSMM classifier.

Footprint discovery in classify_single_molecule_binding_v5_hsmm.py normally runs per
(sample x amplicon) matrix, so low-molecule amplicons have noisy bulks and give wacko UNID calls.
This script instead POOLS the nucleosome(+TF)-decontaminated bulk across BOTH samples AND
copy-number variants that share a promoter (e.g. all BAX_opJS4_{0..8}xTetO), discovers the
off-motif footprint vocabulary ONCE per promoter, and writes a `lo,hi,cost` block per amplicon.
Feed the output to the classifier via --discovered_footprints_file (it then just applies the
prior, no per-amplicon discovery).

COORDINATE SAFETY: the upstream promoter (< the promoter's first-TetO coord, which is fixed across
its copy-number variants thanks to the RC encoding + array-grows-downstream design) is byte-identical
at the same bp coords across variants, so we pool by bp position with NO alignment -- but ONLY over
positions < the (per-promoter) first-TetO boundary. Verify per promoter by sequence before trusting.

Run on a compute node (decodes every grouped matrix once).
"""
import argparse, os, re, sys
from collections import defaultdict
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_single_molecule_binding_v5_hsmm as v5
from common import load_single_molecule_matrix, get_methyl_positions, load_tfbs_positions


def promoter_of(amp):
    return amp.split('_opJS4_')[0] if '_opJS4_' in amp else None   # None => not a TetO-promoter amplicon


def main():
    p = argparse.ArgumentParser(description='pooled per-promoter footprint vocabulary for v5')
    p.add_argument('--results_dir', required=True, help='.../results')
    p.add_argument('--experiment', required=True)
    p.add_argument('--samples', nargs='+', required=True)
    p.add_argument('--positions', required=True, help='positions.long.txt (for TetO coords)')
    p.add_argument('--output_vocab', required=True)
    p.add_argument('--plot', default=None, help='pooled diagnostic PDF (one page/promoter)')
    p.add_argument('--matrix_suffix', default='.dedup.full_unclustered.matrix')
    p.add_argument('--amplicons', nargs='*', default=None, help='default = all TetO amplicons in positions')
    # model params (must match the decode run)
    p.add_argument('--grid_step', type=int, default=5)
    p.add_argument('--p_t_given_unmeth', type=float, default=0.99)
    p.add_argument('--p_t_given_meth', type=float, default=0.05)
    # discovery + grading params (same meaning as the classifier)
    p.add_argument('--footprint_abs_floor', type=float, default=0.25)
    p.add_argument('--footprint_rel_delta', type=float, default=0.20)
    p.add_argument('--footprint_bg_window', type=int, default=80)
    p.add_argument('--footprint_score_hi', type=float, default=1.0)
    p.add_argument('--start_cost_unid', type=float, default=6.0)
    p.add_argument('--unid_discovered_start_cost', type=float, default=2.0)
    p.add_argument('--min_pooled_molecules', type=int, default=50, help='skip a promoter with fewer pooled reads')
    args = p.parse_args()

    pos = load_tfbs_positions(args.positions)
    amplicons = args.amplicons or [a for a in pos if promoter_of(a) is not None]
    groups = defaultdict(list)
    for a in amplicons:
        pr = promoter_of(a)
        if pr is not None:
            groups[pr].append(a)

    def first_teto(amp):
        ts = [int(t[0]) for t in pos.get(amp, [])]
        return min(ts) if ts else None

    def matrix_path(samp, amp):
        return os.path.join(args.results_dir, args.experiment, samp, 'matrices',
                            '{}.{}{}'.format(samp, amp, args.matrix_suffix))

    plots = PdfPages(args.plot) if args.plot else None
    vocab = {}                                   # amplicon -> [(lo,hi,cost),...]
    for pr in sorted(groups):
        amps = groups[pr]
        boundary = min([ft for ft in (first_teto(a) for a in amps) if ft is not None], default=None)
        if boundary is None:
            print('promoter {}: no TetO in any variant, skipping'.format(pr)); continue
        prot = defaultdict(float); cov = defaultdict(float)      # bp pos -> pooled protected / covered
        n_reads = 0
        for amp in amps:
            for samp in args.samples:
                mp = matrix_path(samp, amp)
                if not (os.path.exists(mp) and os.path.getsize(mp) > 0):
                    continue
                try:
                    mat = load_single_molecule_matrix(mp)
                except Exception as e:
                    print('  skip {} {}: {}'.format(samp, amp, e)); continue
                gpc = np.asarray(get_methyl_positions(mat), dtype=np.int64)
                amp_width = int(max(gpc)) + int(210 // 2) + 10
                params = v5.ModelParams(amp_width=amp_width, grid_step=args.grid_step,
                                        p_t_given_unmeth=args.p_t_given_unmeth,
                                        p_t_given_meth=args.p_t_given_meth)
                obs = mat.values.T; obs1 = (obs == 1).astype(float); obs0 = (obs == 0).astype(float)
                base_paths, _ = v5.viterbi_decode(gpc, obs1, obs0, pos.get(amp, []), params)
                frac, nused = v5.per_gpc_protection(mat, gpc, paths=base_paths,
                                                    exclude_states=(v5.NUC, v5.TF))
                n_reads += mat.shape[0]
                for j, gp in enumerate(gpc):
                    if gp < boundary and nused[j] > 0 and np.isfinite(frac[j]):
                        prot[int(gp)] += frac[j] * nused[j]; cov[int(gp)] += nused[j]
        if n_reads < args.min_pooled_molecules or not cov:
            print('promoter {}: too few pooled reads ({}), skipping'.format(pr, n_reads)); continue
        positions = np.array(sorted(cov))
        pooled_frac = np.array([prot[q] / cov[q] for q in positions])
        # reuse the classifier's relative-background discovery on the pooled bulk (upstream => no motifs)
        params0 = v5.ModelParams(amp_width=int(boundary), grid_step=args.grid_step)
        intervals, _ = v5.discover_footprints_from_bulk(
            None, positions, [], params0, frac=pooled_frac,
            abs_floor=args.footprint_abs_floor, rel_delta=args.footprint_rel_delta,
            bg_window=args.footprint_bg_window)
        graded = []
        for (lo, hi, score) in intervals:
            cost = float(np.interp(score, [args.footprint_rel_delta, args.footprint_score_hi],
                                   [args.start_cost_unid, args.unid_discovered_start_cost]))
            graded.append((int(lo), int(hi), round(cost, 2)))
        for amp in amps:
            vocab[amp] = graded
        print('promoter {}: pooled {} reads over {} amplicons; boundary<{}; footprints={}'.format(
            pr, n_reads, len(amps), boundary, graded))
        if plots is not None:
            fig, ax = plt.subplots(figsize=(12, 3.6))
            ax.plot(positions, pooled_frac, '-o', ms=3, color='darkorange',
                    label='pooled decontaminated bulk (n_reads={})'.format(n_reads))
            for thr, cc in [(0.3, '0.75'), (0.4, '0.6'), (0.5, '0.45')]:
                ax.axhline(thr, color=cc, lw=0.8, ls='--')
            for (lo, hi, c) in graded:
                ax.axvspan(lo, hi, color='blueviolet', alpha=0.15)
            ax.axvline(boundary, color='red', lw=1, ls=':', label='first TetO (pool boundary)')
            ax.set_xlim(0, boundary + 40); ax.set_ylim(0, 1)
            ax.set_title('{}  pooled footprint vocabulary: {}'.format(pr, graded), fontsize=9)
            ax.set_xlabel('bp'); ax.set_ylabel('protected frac'); ax.legend(fontsize=8, loc='upper right')
            fig.tight_layout(); plots.savefig(fig); plt.close(fig)

    if plots is not None:
        plots.close()
    with open(args.output_vocab, 'w') as fh:
        for amp in amplicons:
            fh.write('>{}\n'.format(amp))
            for (lo, hi, cost) in vocab.get(amp, []):
                fh.write('{},{},{}\n'.format(lo, hi, cost))
    print('\nwrote vocabulary {} ({} amplicons, {} with footprints)'.format(
        args.output_vocab, len(amplicons), sum(1 for a in amplicons if vocab.get(a))))


if __name__ == '__main__':
    main()
