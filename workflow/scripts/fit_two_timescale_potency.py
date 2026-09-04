#!/usr/bin/env python
"""
fit_two_timescale_potency.py -- potency from v5 output, with the instantaneous
vs integrated decomposition.

Model (see 260902_two_timescale_potency.md; identity link, because Fig. 3d is
straighter than any saturating form):

    P(active | n, <n>) = p0 + k_pot * [ (1-w)*n + w*<n> ]
                       = p0 + k_fast*n + k_slow*<n>

    n    = TFs bound on THIS molecule in THIS snapshot          (v5 n_tf_teto)
    <n>  = mean occupancy of this molecule's CONDITION           (per amplicon)
    k_pot = k_fast + k_slow  = POTENCY (the published Fig. 3d slope)
    w     = k_slow / k_pot   = integrated fraction = tau_int/(tau_int+tau_c)

Because E[n] = <n> by construction, the w terms CANCEL EXACTLY on averaging, so
the bulk Fig. 3d line is p0 + k_pot*<n> for any w. That is why bulk potency is
blind to w while panel h is not: k_pot is identified by BETWEEN-condition
variation and w by WITHIN-condition variation.

    !! k_fast + k_slow == bulk slope is an ALGEBRAIC IDENTITY, not a test.
    Within a condition the N-weighted mean of n IS <n>, so n = <n> + d makes d
    within-condition centred and therefore N-orthogonal to every
    condition-level function; the design splits into orthogonal
    between/within blocks. Verified to hold to machine precision on random
    data. `identity_residual` is only a numerical self-check.

The two tests that DO have content:
  * PARALLELISM -- the panel-h lines must share one slope k_pot*w. Reported as
    `k_interact` (coefficient on n x <n>) and `k_interact_z`; |z| >~ 3 falsifies
    the two-timescale form.
  * CROSS-AXIS AGREEMENT -- w from the copy-number axis and w from the dox axis
    (fixed geometry, `--amplicon_glob <one amplicon> --pool_samples`) should
    agree. They have DIFFERENT confounds (copy number changes which sites are
    bound; dox changes global rTetR and so possibly global cell state), so
    agreement is the real evidence and neither axis alone is decisive.

Promoter activity is taken v5-natively: a molecule is "active" (promoter
nucleosome-free) if NO NUC segment overlaps the promoter window, default
(136,195) -- the same window `assign_promoter_state_from_model_new.py` used, in
the same coordinate frame (verified: v5 segment coords match the positions
file).

USE THE no_endog_meth=FALSE TREE for promoter work: GpC 130 is GCG-ambiguous and
falls inside the TSS window (120,135), and endogenous CpG methylation in these
samples is 14-17%.
"""

import argparse
import glob
import os
from os import path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.switch_backend('agg')

MPLSTYLE = ('/oak/stanford/groups/wjg/bgrd/papers/ad_smf/'
            'activation_domain_smf_paper/ad_smf.mplstyle')


def copy_number(amplicon):
    for tok in amplicon.replace('_', ' ').replace('x', 'x ').split():
        if tok.endswith('x') and tok[:-1].isdigit():
            return int(tok[:-1])
    return np.nan


def load_molecules(binding_dir, sample, amplicon_glob, prom_lo, prom_hi):
    """One row per molecule: sample, amplicon, n_tf, active."""
    rows = []
    pat = path.join(binding_dir, sample,
                    f'{sample}.{amplicon_glob}.single_molecule_classification.txt')
    for f in sorted(glob.glob(pat)):
        amp = path.basename(f).split('.')[-3]
        seg_f = f.replace('.single_molecule_classification.txt', '.segments.txt')
        if not path.exists(seg_f):
            continue
        main = pd.read_table(f)
        if 'n_tf_teto' in main.columns:
            n_tf = main['n_tf_teto'].to_numpy()
        else:                                   # older schema fallback
            tf = [c for c in main.columns
                  if c.startswith('tfbs_') and c.split('_')[-1].isdigit()]
            n_tf = main[tf].to_numpy().sum(1) if tf else np.zeros(len(main))
        seg = pd.read_table(seg_f)
        nuc = seg[seg['type'].str.upper() == 'NUC']
        over = nuc[(nuc.start < prom_hi) & (nuc.end > prom_lo)].read_id.unique()
        active = ~np.isin(main['read_id'].to_numpy(), over)
        rows.append(pd.DataFrame({'sample': sample, 'amplicon': amp,
                                  'copy_number': copy_number(amp),
                                  'n_tf': n_tf, 'active': active}))
    return pd.concat(rows, ignore_index=True) if rows else None


def wls(x, y, wt):
    """Weighted least squares for y ~ 1 + x (columns of x allowed)."""
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float)
                                             for c in np.atleast_2d(x)])
    W = np.asarray(wt, float)
    XtW = X.T * W
    beta = np.linalg.solve(XtW @ X, XtW @ np.asarray(y, float))
    resid = y - X @ beta
    dof = max(len(y) - X.shape[1], 1)
    s2 = float((W * resid ** 2).sum() / dof)
    cov = s2 * np.linalg.inv(XtW @ X)
    return beta, np.sqrt(np.clip(np.diag(cov), 0, None))


def analyse(mol, min_cell=30):
    """Returns (per-amplicon bulk table, per-(amplicon,n_tf) cell table,
    fitted parameters)."""
    bulk = (mol.groupby(['sample', 'amplicon', 'copy_number'])
            .agg(mean_n_tf=('n_tf', 'mean'), frac_active=('active', 'mean'),
                 N=('active', 'size')).reset_index())

    cell = (mol.groupby(['sample', 'amplicon', 'n_tf'])
            .agg(frac_active=('active', 'mean'), N=('active', 'size'))
            .reset_index())
    cell = cell.merge(bulk[['sample', 'amplicon', 'mean_n_tf']],
                      on=['sample', 'amplicon'])
    cell = cell[cell.N >= min_cell].copy()

    out = {}
    # --- bulk (Fig. 3d): potency = the slope, and it is w-independent
    b, se = wls(bulk.mean_n_tf.values, bulk.frac_active.values, bulk.N.values)
    out['p0_bulk'], out['k_pot_bulk'] = b[0], b[1]
    out['p0_bulk_se'], out['k_pot_bulk_se'] = se[0], se[1]

    # --- two-timescale: P = p0 + k_fast*n + k_slow*<n>, molecule-level via
    #     the cell table weighted by cell size (equivalent for a linear
    #     probability model, and far cheaper)
    b2, se2 = wls([cell.n_tf.values, cell.mean_n_tf.values],
                  cell.frac_active.values, cell.N.values)
    p0, k_fast, k_slow = b2
    out.update({'p0': p0, 'k_fast': k_fast, 'k_slow': k_slow,
                'p0_se': se2[0], 'k_fast_se': se2[1], 'k_slow_se': se2[2]})
    k_pot = k_fast + k_slow
    out['k_pot'] = k_pot
    out['w'] = k_slow / k_pot if k_pot else np.nan
    out['tau_ratio'] = (out['w'] / (1 - out['w'])
                        if 0 <= out['w'] < 1 else np.nan)
    # --- NOT a test: k_fast + k_slow == bulk slope is an ALGEBRAIC IDENTITY.
    # Within each amplicon the N-weighted mean of n_tf IS <n>, so
    # n_tf = <n> + d with d within-amplicon centred; d is then N-orthogonal to
    # every function of the amplicon (including 1 and <n>), so the design splits
    # into orthogonal between/within blocks. (p0, k_fast+k_slow) is exactly the
    # bulk fit and k_fast comes only from within-amplicon spread -- they cannot
    # disagree. Verified to hold on random data. Kept only as a numerical
    # self-check that the two fits are consistent.
    out['identity_residual'] = k_pot - out['k_pot_bulk']

    # --- the REAL test: PARALLELISM. The model says every panel-h line has the
    # same slope in <n>, k_pot*w. Allow an interaction n*<n>; if it is non-zero
    # the lines fan out and the two-timescale form is wrong.
    b3, se3 = wls([cell.n_tf.values, cell.mean_n_tf.values,
                   cell.n_tf.values * cell.mean_n_tf.values],
                  cell.frac_active.values, cell.N.values)
    out['k_interact'] = b3[3]
    out['k_interact_se'] = se3[3]
    out['k_interact_z'] = b3[3] / se3[3] if se3[3] > 0 else np.nan
    return bulk, cell, out


def plot(bulk, cell, par, out_pdf, title=''):
    if path.exists(MPLSTYLE):
        plt.style.use(MPLSTYLE)
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)

    ax = axs[0]
    sc = ax.scatter(bulk.mean_n_tf, bulk.frac_active, c=bulk.copy_number,
                    s=np.clip(bulk.N / bulk.N.max() * 120 + 20, 20, 160),
                    cmap='viridis', edgecolors='k', linewidth=0.4, zorder=3)
    xs = np.linspace(0, bulk.mean_n_tf.max() * 1.05, 50)
    ax.plot(xs, par['p0_bulk'] + par['k_pot_bulk'] * xs, 'r-', lw=1.5,
            label=f"k_pot = {par['k_pot_bulk']:.4f} $\\pm$ {par['k_pot_bulk_se']:.4f}")
    ax.set_xlabel('average TF occupancy  $\\langle n\\rangle$')
    ax.set_ylabel('fraction of promoters active')
    ax.set_title('BULK (paper Fig. 3d)\nslope = potency, blind to $w$', fontsize=9)
    ax.legend(fontsize=7, loc='upper left')
    plt.colorbar(sc, ax=ax, label='TetO copy number')

    ax = axs[1]
    ns = sorted(cell.n_tf.unique())
    greens = plt.cm.Greens(np.linspace(0.3, 0.95, len(ns)))
    for c, n in zip(greens, ns):
        sub = cell[cell.n_tf == n].sort_values('mean_n_tf')
        ax.plot(sub.mean_n_tf, sub.frac_active, 'o-', ms=4, lw=1.1, color=c,
                label=str(int(n)))
    ax.set_xlabel('average TF occupancy  $\\langle n\\rangle$')
    ax.set_ylabel('fraction active')
    ax.set_title('PANEL h coordinates: $P$ vs $\\langle n\\rangle$ at fixed $n$\n'
                 'slope $=k_{pot}w$, spacing $=k_{pot}(1-w)$', fontsize=9)
    ax.legend(title='TFs bound', fontsize=6, title_fontsize=6, ncol=2,
              loc='upper left')
    ax.text(0.97, 0.04,
            f"slope   $k_{{pot}}w$ = {par['k_slow']:.4f}\n"
            f"spacing $k_{{pot}}(1-w)$ = {par['k_fast']:.4f}\n"
            f"sum = {par['k_pot']:.4f}  vs bulk {par['k_pot_bulk']:.4f}\n"
            f"$w$ = {par['w']:.3f}",
            transform=ax.transAxes, fontsize=7, ha='right', va='bottom',
            bbox=dict(boxstyle='round', fc='wheat', alpha=0.7))

    fig.suptitle(title, fontsize=10)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f'wrote {out_pdf}')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binding_dir', required=True,
                    help='binding_v5/<experiment> (use the no_endog=FALSE tree)')
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--amplicon_glob', default='opJS4_?x_TetO_21bp_no_CG')
    ap.add_argument('--promoter_lo', type=int, default=136)
    ap.add_argument('--promoter_hi', type=int, default=195)
    ap.add_argument('--min_cell', type=int, default=30,
                    help='min molecules per (amplicon, n_tf) cell')
    ap.add_argument('--pool_samples', action='store_true',
                    help='one joint fit over all samples instead of per-sample')
    ap.add_argument('--output_dir', default='.')
    ap.add_argument('--out_prefix', required=True)
    args = ap.parse_args()

    mols = [m for m in (load_molecules(args.binding_dir, s, args.amplicon_glob,
                                       args.promoter_lo, args.promoter_hi)
                        for s in args.samples) if m is not None]
    if not mols:
        raise SystemExit('no molecules loaded')
    mol = pd.concat(mols, ignore_index=True)
    os.makedirs(args.output_dir, exist_ok=True)
    p = path.join(args.output_dir, args.out_prefix)

    groups = [('ALL', mol)] if args.pool_samples else \
             [(s, g) for s, g in mol.groupby('sample')]
    rows = []
    for name, g in groups:
        bulk, cell, par = analyse(g, min_cell=args.min_cell)
        par['group'] = name
        par['n_molecules'] = len(g)
        par['n_amplicons'] = bulk.shape[0]
        rows.append(par)
        plot(bulk, cell, par, f'{p}.{name}.potency.pdf',
             title=f'{name}  (n={len(g)} molecules)')
        bulk.to_csv(f'{p}.{name}.bulk.txt', sep='\t', index=False)
        cell.to_csv(f'{p}.{name}.cells.txt', sep='\t', index=False)

    res = pd.DataFrame(rows)
    cols = ['group', 'n_molecules', 'n_amplicons', 'k_pot_bulk',
            'k_pot_bulk_se', 'p0_bulk', 'k_fast', 'k_fast_se', 'k_slow',
            'k_slow_se', 'k_pot', 'w', 'tau_ratio', 'k_interact',
            'k_interact_se', 'k_interact_z', 'identity_residual']
    res = res[[c for c in cols if c in res.columns]]
    res.to_csv(f'{p}.two_timescale_potency.txt', sep='\t', index=False)
    print('\n=== two-timescale potency ===')
    print(res.to_string(index=False))
    print('\nNOTE: k_fast + k_slow == k_pot_bulk is an ALGEBRAIC IDENTITY, not a')
    print('test (`identity_residual` is only a numerical self-check).')
    print('The REAL test is PARALLELISM: k_interact should be 0. |z| >~ 3 means')
    print('the panel-h lines fan out and the two-timescale form is wrong.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
