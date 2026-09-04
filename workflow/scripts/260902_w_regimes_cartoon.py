#!/usr/bin/env python
"""
260902_w_regimes_cartoon.py

CARTOON SIMULATION of the two-timescale potency model (see
`260902_two_timescale_potency.md`), demonstrating the central claim:

    The bulk potency line (paper Fig. 3d) is IDENTICAL for every w.
    Panel h (activity stratified by instantaneous n) is the ONLY place
    the w regimes differ.

Model (identity link -- Fig. 3d is straighter than the hyperbolic additive
model, so we use a linear probability model over the observed range):

    P(active | n, <n>) = p0 + k_pot * [ (1-w)*n + w*<n> ]

    n     = TFs bound on THIS molecule in THIS snapshot
    <n>   = mean occupancy of this molecule's CONDITION (K, dox, AD)
    w     = integrated fraction of potency = tau_int / (tau_int + tau_c)

Because E[n] = <n> by construction, averaging over molecules gives

    Pbar = p0 + k_pot*<n>        EXACTLY, for any w

so k_pot (= potency) is the conserved quantity and w redistributes it:

    slope of each panel-h curve in <n>   =  k_pot * w
    spacing between consecutive n curves =  k_pot * (1-w)
    sum                                  =  k_pot

Three regimes:
    w = 0    curves FLAT in <n>, maximally separated by n   (pure instantaneous)
    w = 1    curves COLLAPSE onto one sloping line          (pure history;
                                                             real potency,
                                                             ZERO single-molecule
                                                             signal)
    0<w<1    parallel sloping lines, partially separated

Run with no args for the default figure.
"""

import argparse
from os import path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.switch_backend('agg')

MPLSTYLE = ('/oak/stanford/groups/wjg/bgrd/papers/ad_smf/'
            'activation_domain_smf_paper/ad_smf.mplstyle')

# Mean occupancy per TetO copy number, eyeballed off Fig. 3c (background 0/2).
# Sigmoidal in K -- that shape comes from cooperative binding and is taken as
# given here; this cartoon is about the promoter layer, not the array layer.
MEAN_OCC = {0: 0.0, 1: 0.08, 2: 0.15, 3: 0.45, 4: 1.05, 5: 1.95,
            6: 2.85, 7: 3.55, 8: 4.10}

# Parameters chosen to reproduce Fig. 3d: Pbar runs 0.03 -> ~0.40 as <n> -> 4.1
P0_DEFAULT = 0.03
KPOT_DEFAULT = 0.09


def simulate(w, k_pot=KPOT_DEFAULT, p0=P0_DEFAULT, n_mols=20000,
             mean_occ=None, overdispersion=0.0, rng=None):
    """Generate synthetic molecules under one w regime.

    Occupancy is drawn per molecule as Binomial(K, <n>_K / K) so that the
    condition mean matches `mean_occ` by construction (that identity is what
    makes the bulk cancellation exact). `overdispersion` > 0 switches to a
    beta-binomial with that concentration, making binding more all-or-none
    (closer to real cooperative binding) without changing the mean.
    """
    rng = np.random.default_rng(0 if rng is None else rng)
    mean_occ = MEAN_OCC if mean_occ is None else mean_occ

    rows = []
    for K, mu in mean_occ.items():
        if K == 0:
            n = np.zeros(n_mols, dtype=int)
        else:
            p = np.clip(mu / K, 0.0, 1.0)
            if overdispersion > 0:
                # beta-binomial: same mean p, extra molecule-to-molecule spread
                c = overdispersion
                p_mol = rng.beta(max(p * c, 1e-6), max((1 - p) * c, 1e-6),
                                 size=n_mols)
                n = rng.binomial(K, p_mol)
            else:
                n = rng.binomial(K, p, size=n_mols)

        # the two-timescale drive, then a Bernoulli activity draw
        drive = (1.0 - w) * n + w * mu
        prob = np.clip(p0 + k_pot * drive, 0.0, 1.0)
        active = rng.random(n_mols) < prob

        rows.append(pd.DataFrame({'K': K, 'mean_occ': mu, 'n': n,
                                  'active': active}))

    return pd.concat(rows, ignore_index=True)


def fit_back(df):
    """Recover (p0, k_fast, k_slow) -> (k_pot, w) by OLS on the linear
    probability model: active ~ 1 + n + mean_occ.

    This is the estimator recommended in the design doc: linear in the
    parameters, so no convexity worries, and it is the same weighted-least-
    squares family `compute_potency.py` already uses.
    """
    X = np.column_stack([np.ones(len(df)), df['n'].to_numpy(float),
                         df['mean_occ'].to_numpy(float)])
    y = df['active'].to_numpy(float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    p0_hat, k_fast, k_slow = coef
    k_pot = k_fast + k_slow
    w_hat = k_slow / k_pot if k_pot != 0 else np.nan
    return dict(p0=p0_hat, k_fast=k_fast, k_slow=k_slow,
                k_pot=k_pot, w=w_hat)


def panel_h_table(df, mol_thresh=10, max_n=5):
    """Fraction active per (condition, instantaneous n) -- i.e. paper panel h,
    but with <n> on x (the coordinate the model is linear in) rather than K.
    `mol_thresh` mirrors the real pipeline's per-point molecule floor."""
    g = (df[df['n'] <= max_n]
         .groupby(['K', 'mean_occ', 'n'])['active']
         .agg(['mean', 'size'])
         .reset_index())
    return g[g['size'] > mol_thresh]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--w', type=float, nargs='+', default=[0.0, 0.25, 1.0],
                    help='w regimes to simulate (default: 0 0.25 1)')
    ap.add_argument('--k_pot', type=float, default=KPOT_DEFAULT)
    ap.add_argument('--p0', type=float, default=P0_DEFAULT)
    ap.add_argument('--n_mols', type=int, default=20000,
                    help='molecules per condition')
    ap.add_argument('--overdispersion', type=float, default=0.0,
                    help='beta-binomial concentration; 0 = plain binomial')
    ap.add_argument('--output', type=str,
                    default=path.join(path.dirname(path.abspath(__file__)),
                                      '260902_w_regimes_cartoon.pdf'))
    args = ap.parse_args()

    if path.exists(MPLSTYLE):
        plt.style.use(MPLSTYLE)

    sims = {w: simulate(w, k_pot=args.k_pot, p0=args.p0, n_mols=args.n_mols,
                        overdispersion=args.overdispersion, rng=1)
            for w in args.w}

    nw = len(args.w)
    # the repo mplstyle sets constrained_layout=True; keep it and skip tight_layout
    fig, axs = plt.subplots(2, nw, figsize=(3.6 * nw, 7.0), squeeze=False,
                            constrained_layout=True)

    # ---- row 0: the BULK potency plot (paper Fig. 3d) -- identical for all w
    for j, w in enumerate(args.w):
        ax = axs[0][j]
        bulk = (sims[w].groupby('mean_occ')['active'].mean().reset_index())
        ax.plot(bulk['mean_occ'], bulk['active'], 'o', ms=6, color='k',
                label='simulated')
        xs = np.linspace(0, max(MEAN_OCC.values()), 100)
        ax.plot(xs, args.p0 + args.k_pot * xs, '-', color='crimson', lw=1.5,
                label=r'$p_0 + k_{pot}\langle n\rangle$')
        ax.set_title(f'$w$ = {w:g}', fontweight='bold')
        ax.set_xlabel(r'average TF occupancy $\langle n\rangle$')
        ax.set_ylabel('fraction of promoters active')
        ax.set_ylim(0, 0.55)
        if j == 0:
            ax.legend(fontsize=7, loc='upper left')
            ax.text(0.03, 0.62, 'BULK (Fig. 3d)\nsame line for every $w$',
                    transform=ax.transAxes, fontsize=8, style='italic')

    # ---- row 1: panel h -- this is where the regimes separate
    for j, w in enumerate(args.w):
        ax = axs[1][j]
        tab = panel_h_table(sims[w])
        ns = sorted(tab['n'].unique())
        greens = plt.cm.Greens(np.linspace(0.35, 0.95, len(ns)))
        for c, n in zip(greens, ns):
            sub = tab[tab['n'] == n].sort_values('mean_occ')
            ax.plot(sub['mean_occ'], sub['mean'], 'o-', ms=4, lw=1.2,
                    color=c, label=f'{n}')
        ax.set_xlabel(r'average TF occupancy $\langle n\rangle$')
        ax.set_ylabel('fraction active')
        ax.set_ylim(0, 0.85)
        ax.legend(title='TFs bound', fontsize=6, title_fontsize=6, ncol=2,
                  loc='upper left')

        # annotate the decomposition
        ax.text(0.97, 0.04,
                f'slope $= k_{{pot}}w = {args.k_pot * w:.3f}$\n'
                f'spacing $= k_{{pot}}(1-w) = {args.k_pot * (1 - w):.3f}$',
                transform=ax.transAxes, fontsize=7, ha='right', va='bottom',
                bbox=dict(boxstyle='round', fc='wheat', alpha=0.6))
        if j == 0:
            ax.text(0.40, 0.97, 'PANEL h\n$w$ is visible here',
                    transform=ax.transAxes, fontsize=8, style='italic',
                    va='top')

    fig.savefig(args.output)
    plt.close(fig)
    print(f'wrote {args.output}')

    # ---- recovery check: can we get w back out?
    print('\nrecovery by OLS on  active ~ 1 + n + <n>   '
          f'(true k_pot = {args.k_pot:g}, p0 = {args.p0:g})')
    print(f"{'w_true':>8} {'w_hat':>8} {'k_pot_hat':>10} "
          f"{'k_fast':>8} {'k_slow':>8} {'p0_hat':>8}")
    for w in args.w:
        f = fit_back(sims[w])
        print(f"{w:8.2f} {f['w']:8.3f} {f['k_pot']:10.4f} "
              f"{f['k_fast']:8.4f} {f['k_slow']:8.4f} {f['p0']:8.4f}")

    # ---- and the point of the whole exercise: bulk cannot tell them apart
    print('\nbulk potency slope fitted per regime (should be ~identical):')
    for w in args.w:
        bulk = sims[w].groupby('mean_occ')['active'].mean().reset_index()
        A = np.column_stack([np.ones(len(bulk)),
                             bulk['mean_occ'].to_numpy(float)])
        c, *_ = np.linalg.lstsq(A, bulk['active'].to_numpy(float), rcond=None)
        print(f'   w = {w:4.2f}   intercept = {c[0]:.4f}   slope = {c[1]:.4f}')


if __name__ == '__main__':
    main()
