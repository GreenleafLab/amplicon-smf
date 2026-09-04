#!/usr/bin/env python3
"""
Protection-streak histogram: +/- TF control for the nuc-TF-adjacency question.

Model-free test of whether nucleosome-TF adjacency is REAL (bound TF sits contiguous with a
nucleosome, giving longer protected stretches) vs. a classifier artifact. Takes the single-molecule
matrices from two samples at the SAME amplicon -- one that CAN bind the TF and one that ABSOLUTELY
CANNOT (e.g. 0xTetO or no-dox) -- computes contiguous protection streaks per molecule from the raw
data, and overlays the two streak-length distributions.

Interpretation: if the TF-capable sample has a longer-streak tail than the TF-incapable one, the
long protected stretches are TF-dependent -> nuc-TF adjacency is real, and penalizing it in the
classifier would remove real signal. If the distributions match, the long streaks are not
TF-dependent (pure nucleosomes), and an adjacency penalty is justified.

Self-contained (no repo imports) so it can be run from any directory. Matrix format is the pipeline's
`*.matrix`: header `#<chrom>\\t<pos>\\t<pos>...` (one column per bp of the region), then one row per
molecule `read_id\\t<val>...` with val in {1=protected, 0=accessible, -1=no info}.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from scipy.stats import ks_2samp
except ImportError:
    ks_2samp = None


def load_matrix(path):
    """Return (positions: int array (bp), values: float array (n_mols, n_pos) with NaN for -1)."""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    positions = np.array([int(x) for x in header[1:]])   # first field is '#<chrom>'
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            vals = np.array([float(v) for v in parts[1:]])
            rows.append(vals)
    mat = np.vstack(rows)
    mat[mat < 0] = np.nan                                  # -1 -> NaN (no info)
    return positions, mat


def streaks_for_molecule(positions, row, max_gap):
    """Maximal runs of consecutive INFORMATIVE positions all == protected(1). A run of accessible(0)
    breaks a streak; NaN positions are skipped; a gap > max_gap bp between consecutive informative
    protected positions also breaks the streak. Returns list of streak lengths in bp."""
    inf = ~np.isnan(row)
    p = positions[inf]
    v = row[inf]
    streaks = []
    start = None
    prev = None
    for pos, val in zip(p, v):
        if val == 1 and (start is None or (pos - prev) <= max_gap):
            if start is None:
                start = pos
            prev = pos
        else:
            if start is not None and prev is not None:
                streaks.append(prev - start)
            start = pos if val == 1 else None
            prev = pos if val == 1 else None
    if start is not None and prev is not None:
        streaks.append(prev - start)
    return streaks


def collect(positions, mat, mode, min_streak, max_gap, region):
    if region is not None:
        keep = (positions >= region[0]) & (positions <= region[1])
        positions, mat = positions[keep], mat[:, keep]
    out = []
    for i in range(mat.shape[0]):
        s = streaks_for_molecule(positions, mat[i], max_gap)
        if not s:
            continue
        if mode == 'longest':
            out.append(max(s))
        else:  # 'all'
            out.extend([x for x in s if x >= min_streak])
    return np.array(out, dtype=float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix_can_bind', required=True, help='matrix of the sample that CAN bind TF')
    ap.add_argument('--matrix_cannot_bind', required=True, help='matrix of the sample that CANNOT bind TF')
    ap.add_argument('--label_can', default='TF-capable')
    ap.add_argument('--label_cannot', default='TF-incapable')
    ap.add_argument('--mode', choices=['longest', 'all'], default='longest',
                    help="'longest' = longest streak per molecule; 'all' = every streak >= --min_streak")
    ap.add_argument('--min_streak', type=float, default=100.0,
                    help='bp; in "all" mode keep streaks >= this; also used for the long-streak fraction stat')
    ap.add_argument('--max_gap', type=float, default=30.0,
                    help='bp; gap between consecutive informative protected GpCs above which a streak breaks')
    ap.add_argument('--region', type=str, default=None,
                    help='optional "start,end" bp window to restrict to (e.g. the nuc/TetO region)')
    ap.add_argument('--bins', type=int, default=40)
    ap.add_argument('--output', required=True, help='output histogram PDF')
    ap.add_argument('--stats_output', default=None, help='optional tidy stats TSV')
    args = ap.parse_args()

    region = tuple(int(x) for x in args.region.split(',')) if args.region else None

    pos_can, mat_can = load_matrix(args.matrix_can_bind)
    pos_cannot, mat_cannot = load_matrix(args.matrix_cannot_bind)

    s_can = collect(pos_can, mat_can, args.mode, args.min_streak, args.max_gap, region)
    s_cannot = collect(pos_cannot, mat_cannot, args.mode, args.min_streak, args.max_gap, region)

    fig, ax = plt.subplots(figsize=(8, 5))
    lo = 0
    hi = max(s_can.max() if len(s_can) else 0, s_cannot.max() if len(s_cannot) else 0, 1)
    bins = np.linspace(lo, hi, args.bins + 1)
    ax.hist(s_cannot, bins=bins, density=True, alpha=0.5, label=f'{args.label_cannot} (n={len(s_cannot)})')
    ax.hist(s_can, bins=bins, density=True, alpha=0.5, label=f'{args.label_can} (n={len(s_can)})')
    ax.axvline(147, color='k', ls='--', lw=1, label='147 bp (mono-nuc)')
    ax.set_xlabel('protection streak length (bp)' + (f'  [{args.mode}]'))
    ax.set_ylabel('density')
    ttl = 'Protection streaks: TF-capable vs TF-incapable'
    if region:
        ttl += f'  (region {region[0]}-{region[1]})'
    ax.set_title(ttl)
    ax.legend()
    plt.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    print(f'wrote {args.output}')

    def summ(a):
        if not len(a):
            return dict(n=0, mean=float('nan'), median=float('nan'), p90=float('nan'), frac_long=float('nan'))
        return dict(n=len(a), mean=float(np.mean(a)), median=float(np.median(a)),
                    p90=float(np.percentile(a, 90)), frac_long=float(np.mean(a >= args.min_streak)))
    sc, sn = summ(s_can), summ(s_cannot)
    ks = ks_2samp(s_can, s_cannot) if (ks_2samp and len(s_can) and len(s_cannot)) else None
    print(f'{args.label_can:16s} median={sc["median"]:.0f} p90={sc["p90"]:.0f} '
          f'frac>= {args.min_streak:.0f}bp={sc["frac_long"]:.3f} (n={sc["n"]})')
    print(f'{args.label_cannot:16s} median={sn["median"]:.0f} p90={sn["p90"]:.0f} '
          f'frac>= {args.min_streak:.0f}bp={sn["frac_long"]:.3f} (n={sn["n"]})')
    if ks:
        print(f'KS test: D={ks.statistic:.3f} p={ks.pvalue:.2e}  '
              f'(rightward shift in TF-capable => nuc-TF adjacency likely REAL)')

    if args.stats_output:
        with open(args.stats_output, 'w') as f:
            for lab, d in [(args.label_can, sc), (args.label_cannot, sn)]:
                for k, v in d.items():
                    f.write(f'{lab}\t{k}\t{v}\n')
            if ks:
                f.write(f'ks\tstatistic\t{ks.statistic}\n')
                f.write(f'ks\tpvalue\t{ks.pvalue}\n')
        print(f'wrote {args.stats_output}')


if __name__ == '__main__':
    main()
