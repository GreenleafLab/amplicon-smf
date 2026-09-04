#!/usr/bin/env python
"""
fit_ising_model.py -- Ising / MaxEnt partition-function fit for v5 HSMM output.

Replaces `fit_partition_function_model_v3.py` for v5. Design docs, in order:
    260715_v5_downstream_design.md   §7  (formal derivation)
    260717_ising_intuition.md            (concepts; §7 SUPERSEDED, see below)
    260902_two_timescale_potency.md      (the promoter/potency layer, which is
                                          deliberately NOT part of this model)

WHAT THIS FITS
--------------
An equilibrium lattice model of the TetO array + nucleosomes on one amplicon.
Energy is linear in four integer features of a configuration sigma:

    -E(sigma) = h*n_tf + J*n_pairs + mu*n_nuc + delta*n_nuc*1[n_tf>=1]
    P(sigma)  = exp(-E) / Z

    h     : shared TetO field (intrinsic affinity; ONE value for all operators,
            never per-site -- the middle>edge occupancy gradient must EMERGE
            from excluded volume, see intuition doc §5)
    J     : nearest-neighbour TetO cooperativity (adjacent bound pairs, in
            operator index order)
    mu    : nucleosome fugacity (grand-canonical; NOT a fixed count)
    delta : delta_REMODEL -- nucleosome fugacity shifts when ANY TF is bound
            (recruited remodeler). This is what the old `3param_nuc` fit:
            `num_nuc * (nuc_e - delta_e * (num_tfs>0))`.

delta_STERIC (a TF cannot sit under a nucleosome) is NOT a parameter: it is hard
exclusion, built into the DP. v5's Viterbi segmentation makes the violating
configuration unrepresentable, so its observed count is structurally zero and
moment-matching would drive any soft penalty to infinity regardless.

The promoter is deliberately absent. Published Fig. 3 shows the occupancy->
activity link is line-like rather than convex, so it is not a Boltzmann
coupling, and activity depends on the *available* site count as well as the
bound count (time integration). Potency therefore is NOT an Ising coupling --
it lives in 260902_two_timescale_potency.md. This script fits the array only.

HOW IT FITS
-----------
The likelihood depends on the data ONLY through the mean feature counts (they
are the exponential family's sufficient statistics), so the MLE is moment
matching, <phi_k>_model = <phi_k>_data, and no enumerated microstate table or
per-molecule `idx` join is needed -- which is exactly what lets it consume v5.

    NLL(theta) = sum_a N_a [ log Z_a(theta) - theta . <phi>_data,a ]
    dNLL/dtheta_k = sum_a N_a ( <phi_k>_model,a - <phi_k>_data,a )

Z is computed EXACTLY by an O(T) forward DP (transfer matrix for the operator
chain x hard-rod lattice gas for nucleosomes), never by enumeration. The
boundary is grand-canonical: the lattice is padded by >= one footprint on each
side and rods may overhang inward, because the amplicon edge is a measurement
window, not a physical wall -- hard walls manufacture fake edge phasing.

A "model flavor" is a boolean MASK over the one superset parameter vector, not
a separate energy function. That makes nestedness structural (so the LRT is
valid by construction) and tells you which posterior checks are tautological
(a fitted moment) versus informative (everything else).

USAGE
-----
    # verify the implementation with no data at all (recommended first):
    python fit_ising_model.py selftest

    # inspect what a v5 output file actually contains before trusting it:
    python fit_ising_model.py columns --main_files <f> [<f> ...]

    # fit, jointly across amplicons, and write checks:
    python fit_ising_model.py fit \
        --main_files  <sample>.<amp>.single_molecule_classification.txt ... \
        --matrices    <sample>.<amp>.dedup.full_unclustered.matrix ... \
        --positions   260713_opJS45.positions.long.bare.txt \
        --specs 2param 3param_tfcoop 3param_nuc 4param \
        --out_prefix  260903_ising

Positions files MUST use bare-operator coordinates
(`convert_fa_to_positions_for_script.py --l_offset 0 --r_offset 0`); the old
+2/-1-asymmetric opJS45 file mis-calls TetOs.
"""

import argparse
import itertools
import os
import sys
from os import path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp

# ---------------------------------------------------------------------------
# Parameter vector and model flavors
# ---------------------------------------------------------------------------
# Canonical order. theta is ALWAYS length 4; a spec pins the unused entries to 0.
FEATURES = ['h', 'J', 'mu', 'delta', 'delta_adj']
# The empirical feature conjugate to each parameter, in the same order.
FEATURE_STATS = ['n_tf', 'n_pairs', 'n_nuc', 'n_nuc_and_tf', 'n_nuc_adj']

SPECS = {
    '2param':        ['h', 'mu'],
    '3param_tfcoop': ['h', 'mu', 'J'],
    # delta = delta_remodel, SWITCH form: shifts every nucleosome's fugacity if
    # ANY TF is bound. Saturates immediately -- one bound TF gives the same
    # effect as eight -- so it cannot produce a graded <n_nuc> decline with
    # copy number, which is what the data show.
    '3param_nuc':    ['h', 'mu', 'delta'],
    '4param':        ['h', 'mu', 'J', 'delta'],
    # delta_adj = ADJACENCY-LOCAL remodelling: a nucleosome pays the shift only
    # if the nearest operator 5' of it is bound. Naturally GRADED (more bound
    # operators => more affected rod positions) and LOCAL, so it needs no extra
    # DP state -- it reuses the cooperativity bit `s`. Asymmetric by
    # construction (5'/promoter-side only); symmetric adjacency would cost a
    # further state bit.
    '3param_nucadj': ['h', 'mu', 'delta_adj'],
    '4param_adj':    ['h', 'mu', 'J', 'delta_adj'],
    '5param':        ['h', 'mu', 'J', 'delta', 'delta_adj'],
}

NEG = -np.inf

# Fixed ambient nucleosome fugacity of the EXTERIOR. Not fitted: it is the
# one number summarising the infinite bulk array outside the measurement
# window (intuition doc §5.4, "the bath is one number"). Keeping it separate
# from the fitted `mu` is what makes d logZ/d mu equal the observable count.
MU_AMBIENT_DEFAULT = -1.0


def spec_mask(spec):
    """Boolean mask over FEATURES for a named spec."""
    if spec not in SPECS:
        raise ValueError(f'unknown spec {spec!r}; choose from {sorted(SPECS)}')
    active = set(SPECS[spec])
    return np.array([f in active for f in FEATURES], dtype=bool)


def expand(theta_active, mask):
    """Scatter the free parameters into the full length-4 vector (rest = 0)."""
    theta = np.zeros(len(FEATURES))
    theta[mask] = theta_active
    return theta


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
class Geometry:
    """1D bp lattice for one amplicon.

    operator_positions : bp coordinate of each operator, in the model frame.
        Each operator is idealized as a SINGLE bin at the motif midpoint. The
        real TetO is ~20 bp, so steric exclusion against a 147 bp rod has about
        +/-10 bp of slop at the boundary. That is deliberate: representing the
        full TF span would need another coverage-countdown dimension in the DP
        state, and 10 bp is small against a 147 bp footprint. Revisit only if a
        claim turns on the exact exclusion boundary.
    scored_lo/scored_hi : first/last scored GpC (the measurement window). Taken
        from the matrix columns, which is exactly v5's decoded span.
    L    : nucleosome rod length.
    pad  : bins added on EACH side beyond the scored window, so rods may
        overhang inward from the reservoir. Must be >= L. NEVER a hard wall at
        the amplicon edge (intuition doc §5.4).
    """

    def __init__(self, operator_positions, scored_lo, scored_hi, L=147, pad=None,
                 name=None):
        self.name = name
        self.L = int(L)
        self.pad = int(pad) if pad is not None else self.L
        if self.pad < self.L:
            raise ValueError(f'pad ({self.pad}) must be >= L ({self.L}): a pad '
                             'shorter than one footprint reintroduces a hard wall')
        self.scored_lo, self.scored_hi = int(scored_lo), int(scored_hi)
        self.offset = self.scored_lo - self.pad
        self.T = (self.scored_hi + self.pad) - self.offset + 1
        ops = sorted(int(round(p)) - self.offset for p in operator_positions)
        self.op_bin = np.array([p for p in ops if 0 <= p < self.T], dtype=int)
        if len(self.op_bin) != len(ops):
            dropped = len(ops) - len(self.op_bin)
            print(f'  [warn] {name}: {dropped} operator(s) fall outside the '
                  'padded lattice and were dropped', file=sys.stderr)
        self.K = len(self.op_bin)
        self.is_teto = np.zeros(self.T, dtype=bool)
        self.is_teto[self.op_bin] = True
        self.op_of_bin = -np.ones(self.T, dtype=int)
        for k, b in enumerate(self.op_bin):
            self.op_of_bin[b] = k
        # Which bins are inside the SCORED window (as opposed to the pad). The
        # fitted `mu` is conjugate to rods starting in the window, because that
        # is what v5 can actually observe; the pad carries the fixed ambient
        # fugacity instead. Keeping those separate is what makes
        # d logZ / d mu == <n_nuc observed in the window> hold exactly.
        self.in_window = np.zeros(self.T, dtype=bool)
        self.in_window[self.pad:self.T - self.pad] = True
        # max rods that fit; sets the DP's nucleosome-count dimension
        self.max_nuc = self.T // self.L + 1

    def __repr__(self):
        return (f'Geometry({self.name}, K={self.K}, T={self.T}, L={self.L}, '
                f'pad={self.pad}, window=[{self.scored_lo},{self.scored_hi}])')


def load_positions(positions_file, label_filter='TetO'):
    """Parse a bare-operator positions file -> {amplicon: [(lo, hi, label), ...]}.

    Only rows whose label contains `label_filter` are kept (pass None to keep
    all). This matters: the same file can carry non-operator annotations.
    """
    out, amp = {}, None
    with open(positions_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                amp = line[1:]
                out[amp] = []
            else:
                parts = line.split(',')
                lo, hi = int(parts[0]), int(parts[1])
                lab = parts[2] if len(parts) > 2 else ''
                if label_filter is None or label_filter in lab:
                    out[amp].append((lo, hi, lab))
    return out


def scored_window_from_matrix(matrix_file):
    """(first, last) scored GpC coordinate, read from a matrix file's header."""
    with open(matrix_file) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    cols = [int(c) for c in header[1:]]
    if not cols:
        raise ValueError(f'{matrix_file}: no position columns in header')
    return min(cols), max(cols)


# ---------------------------------------------------------------------------
# log Z  --  the whole model
# ---------------------------------------------------------------------------
def _bulk_log_lambda(mu, L):
    """log of the dominant eigenvalue of the bulk hard-rod transfer matrix.

    Solves  lambda^(L-1) * (lambda - 1) = z = e^mu  for lambda > 1, in log
    space. Verified against np.linalg.eigvals of the explicit matrix.
    """
    from scipy.optimize import brentq
    f = lambda x: (L - 1) * x + np.log(np.expm1(x)) - mu
    lo, hi = 1e-13, 1.0
    while f(hi) < 0.0:                      # push the upper end above the root
        hi *= 2.0
        if hi > 1e6:
            raise RuntimeError(f'no bulk lambda for mu={mu}, L={L}')
    while f(lo) > 0.0:                      # ...and the lower end below it.
        lo *= 1e-3                          # log(expm1(x)) -> -inf as x -> 0+,
        if lo < 1e-300:                     # so this always terminates.
            raise RuntimeError(f'no bulk lambda for mu={mu}, L={L}')
    return brentq(f, lo, hi, xtol=1e-15, rtol=8.9e-16)


def _bulk_boundary(mu, L):
    """(log r, log l) -- the stationary boundary vectors for coverage state c.

    r is the dominant RIGHT eigenvector: the distribution of coverage state
    entering the window from an infinite bulk array on the left.
    l is the dominant LEFT eigenvector: the weight of continuing a given
    coverage state into the infinite bulk on the right.

        r[0] = 1,  r[j] = lambda^(j-1) * (lambda - 1)
        l[c] = lambda^(-c)

    This is the whole of "the exterior is bulk chromatin at ambient density":
    the infinite bath enters as two vectors and one number (mu). It replaces
    the hard walls that a plain padded lattice still has -- padding alone does
    NOT remove them, and the resulting edge depletion oscillates inward with
    period ~L for far longer than any affordable pad.
    """
    x = _bulk_log_lambda(mu, L)
    j = np.arange(1, L)
    log_r = np.concatenate([[0.0], (j - 1) * x + np.log(np.expm1(x))])
    log_l = -np.arange(L) * x
    return log_r, log_l


def _rod_weights(geom, mu, mu_ambient, bump=None):
    """Per-bin log-weight for starting a rod at that bin.

    Window bins get the fitted `mu`; pad bins get the fixed `mu_ambient`.
    `bump` (length T) is added on top, used to probe occupancy by finite
    difference without touching the boundary condition.
    """
    w = np.where(geom.in_window, float(mu), float(mu_ambient))
    if bump is not None:
        w = w + np.asarray(bump, dtype=float)
    return w


def _log_Z_core(h, J, rod_w, geom, mu_ambient, extra_field=None,
                allow_binding=True, boundary='stationary', seed_bump=None,
                delta_adj=0.0):
    """log Z with delta = 0. DP state V[c, s]:
         c in 0..L-1 : rod coverage still owed after this bin (0 = free)
         s in {0,1}  : was the most recent operator bound (carries J)

    `rod_w[i]` is the log-weight of starting a rod at bin i (see _rod_weights).

    boundary='stationary' : grand-canonical. Seed with the bulk right
        eigenvector, let rods start anywhere (they may overhang either edge into
        the reservoir), and contract with the bulk left eigenvector. The window
        is then an exact marginal of an infinite array.
    boundary='wall' : the naive padded-but-walled lattice. Retained ONLY so the
        brute-force enumeration test has something finite to compare against --
        it produces ~36x edge depletion oscillating inward with period ~L, so
        never use it for a fit.

    Steric exclusion is automatic: an operator reached via the coverage
    countdown cannot bind, and a rod starting at i only covers bins >= i, so
    every operator it occludes is still in the future. No lookahead needed.
    """
    L, T = geom.L, geom.T
    op_field = np.full(max(geom.K, 1), float(h))
    if extra_field is not None:
        op_field = op_field + np.asarray(extra_field, dtype=float)

    V = np.full((L, 2), NEG)
    if boundary == 'stationary':
        log_r, log_l = _bulk_boundary(mu_ambient, L)
        if seed_bump is not None:
            # Rods that began OUTSIDE the lattice are folded into the seed
            # vector rather than appearing as explicit start events, so
            # probing their fugacity means bumping seed components directly.
            log_r = log_r + np.asarray(seed_bump, dtype=float)
        V[:, 0] = log_r
    else:
        V[0, 0] = 0.0

    for i in range(T):
        teto = geom.is_teto[i]
        Vn = np.full((L, 2), NEG)
        # (1) continuation: c -> c-1. A covered operator is unbound, which
        #     collapses the cooperativity state.
        if teto:
            Vn[0:L - 1, 0] = np.logaddexp(V[1:L, 0], V[1:L, 1])
        else:
            Vn[0:L - 1, :] = V[1:L, :]
        v0 = V[0, :]
        # (2a) start a rod at i. delta_adj is paid when the most recent
        #      operator was bound (s == 1) -- that is the whole adjacency-local
        #      remodelling term, and it is free because `s` is already carried.
        if boundary == 'stationary' or i + L - 1 < T:
            w = rod_w[i]
            if teto:
                # bin i is an operator and the new rod covers it, so that
                # operator is forced unbound and s collapses to 0
                Vn[L - 1, 0] = np.logaddexp(
                    Vn[L - 1, 0],
                    np.logaddexp(v0[0] + w, v0[1] + w + delta_adj))
            else:
                Vn[L - 1, 0] = np.logaddexp(Vn[L - 1, 0], v0[0] + w)
                Vn[L - 1, 1] = np.logaddexp(Vn[L - 1, 1],
                                            v0[1] + w + delta_adj)
        # (2b) no rod starts at i
        if teto:
            if allow_binding:
                kf = op_field[geom.op_of_bin[i]]
                Vn[0, 1] = np.logaddexp(Vn[0, 1],
                                        np.logaddexp(v0[0] + kf, v0[1] + kf + J))
            Vn[0, 0] = np.logaddexp(Vn[0, 0], np.logaddexp(v0[0], v0[1]))
        else:
            Vn[0, :] = np.logaddexp(Vn[0, :], v0)
        V = Vn

    if boundary == 'stationary':
        return logsumexp(V + log_l[:, None])
    return logsumexp(V[0, :])


def _log_Z_delta_exact(theta, geom, mu_ambient, extra_field=None,
                       boundary='stationary', bump=None):
    """log Z including delta_remodel, by augmenting the DP state.

    delta multiplies n_nuc * 1[n_tf>=1], which is NON-LOCAL: a rod placed at the
    left end cannot know whether a TF will bind at the right end. Fix: carry
    both the running rod count and a "has any TF bound yet" bit, then pay
    delta*m*f once at termination.

    State V[c, s, f, m]. `m` counts rods started INSIDE THE SCORED WINDOW only,
    matching what v5 can observe (and hence what the empirical n_nuc counts).
    ~L*2*2*(max_nuc+1) states, all-positive log-sum-exp, numerically robust.
    """
    h, J, mu, delta, delta_adj = theta
    L, T, M = geom.L, geom.T, geom.max_nuc
    rod_w = _rod_weights(geom, mu, mu_ambient, bump)
    op_field = np.full(max(geom.K, 1), float(h))
    if extra_field is not None:
        op_field = op_field + np.asarray(extra_field, dtype=float)

    V = np.full((L, 2, 2, M + 1), NEG)
    if boundary == 'stationary':
        log_r, log_l = _bulk_boundary(mu_ambient, L)
        V[:, 0, 0, 0] = log_r
    else:
        V[0, 0, 0, 0] = 0.0

    for i in range(T):
        teto = geom.is_teto[i]
        counts = geom.in_window[i]        # does a rod starting here count?
        Vn = np.full((L, 2, 2, M + 1), NEG)
        if teto:
            Vn[0:L - 1, 0] = np.logaddexp(V[1:L, 0], V[1:L, 1])
        else:
            Vn[0:L - 1] = V[1:L]
        v0 = V[0]                                   # [s, f, m]
        if boundary == 'stationary' or i + L - 1 < T:
            w = rod_w[i]
            if teto:
                # s collapses to 0; s=1 sources pay the adjacency term
                src = np.logaddexp(v0[0] + w, v0[1] + w + delta_adj)  # [f, m]
                if counts:
                    Vn[L - 1, 0, :, 1:] = np.logaddexp(Vn[L - 1, 0, :, 1:],
                                                       src[:, :-1])
                else:
                    Vn[L - 1, 0] = np.logaddexp(Vn[L - 1, 0], src)
            else:
                add = np.stack([np.full_like(v0[0], w),
                                np.full_like(v0[1], w + delta_adj)])
                if counts:
                    Vn[L - 1, :, :, 1:] = np.logaddexp(
                        Vn[L - 1, :, :, 1:], (v0 + add)[:, :, :-1])
                else:
                    Vn[L - 1] = np.logaddexp(Vn[L - 1], v0 + add)
        if teto:
            kf = op_field[geom.op_of_bin[i]]
            bound = np.logaddexp(v0[0] + kf, v0[1] + kf + J)   # [f, m]
            # binding sets f = 1 regardless of its previous value
            Vn[0, 1, 1, :] = np.logaddexp(Vn[0, 1, 1, :],
                                          logsumexp(bound, axis=0))
            Vn[0, 0] = np.logaddexp(Vn[0, 0], np.logaddexp(v0[0], v0[1]))
        else:
            Vn[0] = np.logaddexp(Vn[0], v0)
        V = Vn

    m_idx = np.arange(M + 1)[None, None, :]
    f_idx = np.array([0.0, 1.0])[None, :, None]
    payoff = delta * m_idx * f_idx
    if boundary == 'stationary':
        return logsumexp(V + payoff[None, :, :, :] + log_l[:, None, None, None])
    return logsumexp(V[0] + payoff)


def _log_Z_delta_identity(theta, geom, mu_ambient, extra_field=None,
                          boundary='stationary', bump=None):
    """log Z including delta_remodel, by conditional decomposition -- the same
    field-shift trick that handles any binary molecule-level variable:

        Z = Z_noTF(mu) + [ Z_all(mu + delta) - Z_noTF(mu + delta) ]

    Three DP calls, no extra state, so much cheaper than the exact version. But
    it needs a log-DIFFERENCE, which loses precision when few TFs bind
    (Z_all -> Z_noTF). Kept as a fast path and, mainly, as an independent
    implementation to cross-test the augmented DP against.
    """
    h, J, mu, delta, delta_adj = theta
    w_base = _rod_weights(geom, mu, mu_ambient, bump)
    w_shift = _rod_weights(geom, mu + delta, mu_ambient, bump)
    lz_no_base = _log_Z_core(h, J, w_base, geom, mu_ambient, extra_field,
                             allow_binding=False, boundary=boundary,
                             delta_adj=delta_adj)
    lz_all_shift = _log_Z_core(h, J, w_shift, geom, mu_ambient, extra_field,
                               allow_binding=True, boundary=boundary,
                               delta_adj=delta_adj)
    lz_no_shift = _log_Z_core(h, J, w_shift, geom, mu_ambient, extra_field,
                              allow_binding=False, boundary=boundary,
                              delta_adj=delta_adj)
    d = np.minimum(lz_no_shift - lz_all_shift, -1e-16)
    diff = lz_all_shift + np.log1p(-np.exp(d))
    return np.logaddexp(lz_no_base, diff)


def log_Z(theta, geom, mu_ambient, extra_field=None, method='exact',
          boundary='stationary', bump=None, seed_bump=None):
    """Dispatch. When delta == 0 the compact (c, s) DP is exact and fastest."""
    theta = np.asarray(theta, dtype=float)
    h, J, mu, delta, delta_adj = theta
    if delta == 0.0:
        # delta_adj is local, so the compact (c, s) DP stays exact
        rod_w = _rod_weights(geom, mu, mu_ambient, bump)
        return _log_Z_core(h, J, rod_w, geom, mu_ambient, extra_field,
                           boundary=boundary, seed_bump=seed_bump,
                           delta_adj=delta_adj)
    if method == 'identity':
        return _log_Z_delta_identity(theta, geom, mu_ambient, extra_field,
                                     boundary, bump)
    return _log_Z_delta_exact(theta, geom, mu_ambient, extra_field,
                              boundary, bump)


# ---------------------------------------------------------------------------
# Moments, occupancy, and the n_tf distribution
# ---------------------------------------------------------------------------
def model_moments(theta, geom, mu_ambient, eps=1e-4, method='exact',
                  boundary='stationary'):
    """<phi>_model = d log Z / d theta  (exponential-family identity),
    by central finite difference. Returns length-4, in FEATURES order.

    Because `mu` weights only rods started inside the scored window while the
    boundary is set by the fixed `mu_ambient`, d logZ / d mu is exactly
    <n_nuc observed in the window> -- the statistic v5 actually reports. Tying
    mu to the boundary instead would add an exterior term and silently break
    moment matching for the nucleosome parameters.
    """
    theta = np.asarray(theta, dtype=float)
    g = np.zeros(len(FEATURES))
    for k in range(len(FEATURES)):
        tp, tm = theta.copy(), theta.copy()
        tp[k] += eps
        tm[k] -= eps
        g[k] = (log_Z(tp, geom, mu_ambient, method=method, boundary=boundary)
                - log_Z(tm, geom, mu_ambient, method=method,
                        boundary=boundary)) / (2 * eps)
    return g


def _h_nodes(sigma_h, n_nodes=5):
    """Gauss-Hermite nodes/weights for h ~ N(0, sigma_h^2), as OFFSETS from h0.

    Each molecule is taken to draw its own h and then equilibrate, so the model
    is a MIXTURE of Boltzmann distributions and every model quantity is the
    weighted average of the per-node quantity. sigma_h = 0 collapses to a single
    node, recovering the homogeneous model exactly.
    """
    if sigma_h <= 0:
        return np.array([0.0]), np.array([1.0])
    x, w = np.polynomial.hermite.hermgauss(int(n_nodes))
    return np.sqrt(2.0) * float(sigma_h) * x, w / np.sqrt(np.pi)


def model_moments_het(theta, geom, mu_ambient, sigma_h=0.0, n_nodes=5,
                      eps=1e-4, method='exact', boundary='stationary'):
    """<phi>_model marginalised over per-molecule heterogeneity in h."""
    offs, wts = _h_nodes(sigma_h, n_nodes)
    acc = np.zeros(len(FEATURES))
    for off, w in zip(offs, wts):
        th = np.asarray(theta, dtype=float).copy()
        th[FEATURES.index('h')] += off
        acc += w * model_moments(th, geom, mu_ambient, eps=eps, method=method,
                                 boundary=boundary)
    return acc


def n_tf_distribution_het(theta, geom, mu_ambient, sigma_h=0.0, n_nodes=5,
                          boundary='stationary'):
    """P(n_tf) marginalised over h heterogeneity. This is the statistic that
    actually IDENTIFIES sigma_h -- the mean moments are only weakly sensitive
    to it (via Jensen), whereas the distribution's spread is strongly so."""
    offs, wts = _h_nodes(sigma_h, n_nodes)
    acc = None
    for off, w in zip(offs, wts):
        th = np.asarray(theta, dtype=float).copy()
        th[FEATURES.index('h')] += off
        pj = n_tf_distribution(th, geom, mu_ambient, boundary=boundary)
        acc = w * pj if acc is None else acc + w * pj
    return acc


def fit_het(datasets, spec, sigma_h, mu_ambient, n_nodes=5, method='exact',
            boundary='stationary', theta0=None):
    """Moment-match at FIXED sigma_h by least squares on the moment residuals.

    Not maximum likelihood: with a mixture the likelihood no longer reduces to
    the mean moments, so this is a method-of-moments (GMM) estimator. That is
    fine for the intended use -- asking whether h, J, mu, delta MOVE when
    plausible per-molecule h heterogeneity is allowed.
    """
    from scipy.optimize import least_squares
    mask = spec_mask(spec)
    if theta0 is None:
        start = {'h': 1.0, 'J': 0.0, 'mu': -0.2, 'delta': -1.0,
                 'delta_adj': 0.0}
        theta0 = np.array([start[f] for f in FEATURES])[mask]
    ntot = sum(d.N for d in datasets)

    def resid(x):
        theta = expand(x, mask)
        r = []
        for d in datasets:
            mm = model_moments_het(theta, d.geom, mu_ambient, sigma_h,
                                   n_nodes, method=method, boundary=boundary)
            r.append(np.sqrt(d.N / ntot) * (mm - d.m)[mask])
        return np.concatenate(r)

    res = least_squares(resid, np.asarray(theta0, float), method='lm',
                        xtol=1e-10, ftol=1e-10)
    return expand(res.x, mask), res


def per_site_occupancy(theta, geom, mu_ambient, eps=1e-4, method='exact',
                       boundary='stationary'):
    """P(operator k bound) = d log Z / d h_k, via a per-site field bump.

    NOT a fitted quantity -- only the TOTAL <n_tf> enters the likelihood -- so
    the shape of this profile is a genuine out-of-sample prediction, and the
    middle>edge gradient emerging from a single shared h is the headline claim.
    """
    occ = np.zeros(geom.K)
    for k in range(geom.K):
        bump = np.zeros(geom.K)
        bump[k] = eps
        lp = log_Z(theta, geom, mu_ambient, extra_field=bump, method=method,
                   boundary=boundary)
        lm = log_Z(theta, geom, mu_ambient, extra_field=-bump, method=method,
                   boundary=boundary)
        occ[k] = (lp - lm) / (2 * eps)
    return occ


def nuc_occupancy_profile(theta, geom, mu_ambient, eps=1e-3,
                          boundary='stationary'):
    """P(bin b covered by a nucleosome), by bumping the fugacity of every rod
    that covers b.

    Used for the boundary-condition correctness check: a barrier-free uniform
    stretch must give FLAT interior occupancy. Any phasing there means the
    lattice edge is acting as a spurious barrier.
    """
    T, L = geom.T, geom.L
    occ = np.zeros(T)
    for b in range(T):
        lo = max(0, b - L + 1)
        hi = b + 1                       # rods starting in [lo, hi) cover b
        bump = np.zeros(T)
        bump[lo:hi] = eps
        # Bin b is also covered by an exterior-initiated rod whenever the
        # entering coverage state exceeds b; those live in the seed vector.
        sb = np.where(np.arange(L) > b, eps, 0.0) \
            if boundary == 'stationary' else None
        lp = log_Z(theta, geom, mu_ambient, boundary=boundary, bump=bump,
                   seed_bump=sb)
        lm = log_Z(theta, geom, mu_ambient, boundary=boundary, bump=-bump,
                   seed_bump=None if sb is None else -sb)
        occ[b] = (lp - lm) / (2 * eps)
    return occ


def n_tf_distribution(theta, geom, mu_ambient, boundary='stationary'):
    """P(n_tf = j) for j = 0..K, by carrying the bound count in the DP state.

    Only the MEAN of n_tf is fit, so the whole shape is a prediction. This is
    the check that detects population heterogeneity: a mixture of "on" and
    "off" cells shows up as an over-dispersed / bimodal n_tf distribution that
    no single parameter set can reproduce.
    """
    h, J, mu, delta, delta_adj = np.asarray(theta, dtype=float)
    L, T, K = geom.L, geom.T, geom.K
    M = geom.max_nuc
    rod_w = _rod_weights(geom, mu, mu_ambient)
    nm = (M + 1) if delta != 0.0 else 1

    V = np.full((L, 2, K + 1, nm), NEG)
    if boundary == 'stationary':
        log_r, log_l = _bulk_boundary(mu_ambient, L)
        V[:, 0, 0, 0] = log_r
    else:
        V[0, 0, 0, 0] = 0.0

    for i in range(T):
        teto = geom.is_teto[i]
        counts = geom.in_window[i] and delta != 0.0
        Vn = np.full((L, 2, K + 1, nm), NEG)
        if teto:
            Vn[0:L - 1, 0] = np.logaddexp(V[1:L, 0], V[1:L, 1])
        else:
            Vn[0:L - 1] = V[1:L]
        v0 = V[0]                                    # [s, j, m]
        if boundary == 'stationary' or i + L - 1 < T:
            w = rod_w[i]
            if teto:
                src = np.logaddexp(v0[0] + w, v0[1] + w + delta_adj)  # [j, m]
                if counts:
                    Vn[L - 1, 0, :, 1:] = np.logaddexp(Vn[L - 1, 0, :, 1:],
                                                       src[:, :-1])
                else:
                    Vn[L - 1, 0] = np.logaddexp(Vn[L - 1, 0], src)
            else:
                add = np.stack([np.full_like(v0[0], w),
                                np.full_like(v0[1], w + delta_adj)])
                if counts:
                    Vn[L - 1, :, :, 1:] = np.logaddexp(
                        Vn[L - 1, :, :, 1:], (v0 + add)[:, :, :-1])
                else:
                    Vn[L - 1] = np.logaddexp(Vn[L - 1], v0 + add)
        if teto:
            bound = np.logaddexp(v0[0] + h, v0[1] + h + J)     # [j, m]
            Vn[0, 1, 1:] = np.logaddexp(Vn[0, 1, 1:], bound[:-1])
            Vn[0, 0] = np.logaddexp(Vn[0, 0], np.logaddexp(v0[0], v0[1]))
        else:
            Vn[0] = np.logaddexp(Vn[0], v0)
        V = Vn

    if delta != 0.0:
        m_idx = np.arange(M + 1)[None, None, :]
        f_idx = (np.arange(K + 1) >= 1).astype(float)[None, :, None]
        payoff = delta * m_idx * f_idx
    else:
        payoff = 0.0
    if boundary == 'stationary':
        tot = V + payoff + log_l[:, None, None, None]
        lz_j = logsumexp(tot, axis=(0, 1, 3))
    else:
        lz_j = logsumexp(V[0] + payoff, axis=(0, 2))
    return np.exp(lz_j - logsumexp(lz_j))


# ---------------------------------------------------------------------------
# Empirical moments from v5 output
# ---------------------------------------------------------------------------
def _resolve_site_matrix(df, expect_K=None):
    """Extract the per-operator boolean occupancy matrix, in operator order.

    v5's schema is in flux (see NOTES_v5_hsmm_status.md and the stable-site
    naming work), so try the documented conventions in order and fail LOUDLY
    with the available columns rather than guessing. Order matters: n_pairs is
    defined on adjacent operators, so columns must be sorted by operator index.
    """
    # 1) documented v5 schema: tfbs_1..K
    tf_cols = [c for c in df.columns
               if c.startswith('tfbs_') and c.split('_')[-1].isdigit()]
    if tf_cols:
        tf_cols.sort(key=lambda c: int(c.split('_')[-1]))
        return df[tf_cols].to_numpy().astype(bool), tf_cols
    # 2) stable-site naming: site_<name>, TetO sites only
    site_cols = [c for c in df.columns if c.startswith('site_')]
    teto_cols = [c for c in site_cols if 'teto' in c.lower()]
    if teto_cols:
        def keyf(c):
            digits = ''.join(ch for ch in c if ch.isdigit())
            return int(digits) if digits else 0
        teto_cols.sort(key=keyf)
        print(f'  [warn] no tfbs_N columns; falling back to site_* columns '
              f'{teto_cols} -- CHECK that this ordering is by operator '
              f'position, since n_pairs depends on it', file=sys.stderr)
        return df[teto_cols].to_numpy().astype(bool), teto_cols
    if expect_K == 0:
        # A legitimately operator-free amplicon (opJS4_0x / opJS5_0x): v5 emits
        # no tfbs_/site_ columns at all. Keep it -- with zero TF competition it
        # is the cleanest constraint on the nucleosome fugacity mu, and it
        # contributes n_tf = n_pairs = 0 by construction.
        return np.zeros((len(df), 0), dtype=bool), []
    raise ValueError(
        'could not find per-operator occupancy columns (looked for tfbs_<N> '
        'then site_*teto*). Available columns:\n  ' + '\n  '.join(map(str, df.columns)))


def _n_nuc_adj_from_segments(segments_file, sites, read_ids, op_mid):
    """Per-molecule count of nucleosomes whose nearest PRECEDING operator is
    bound -- the empirical statistic conjugate to delta_adj.

    Uses operator MIDPOINTS, because the DP idealizes each operator as a single
    bin at its midpoint and sets `s` when it passes that bin. v5's segment
    coordinates share the positions-file frame (verified: a TetO1 TF segment is
    280-303 against a motif window of 282-301, the difference being
    tf_margin=2), so the two are directly comparable.
    """
    seg = pd.read_table(segments_file)
    nuc = seg[seg['type'].str.upper() == 'NUC']
    counts = pd.Series(0, index=read_ids, dtype=float)
    if not len(nuc) or not len(op_mid):
        return counts.to_numpy()
    row_of = {r: i for i, r in enumerate(read_ids)}
    rows = nuc.read_id.map(row_of)
    keep = rows.notna().to_numpy()
    if not keep.any():
        return counts.to_numpy()
    rows = rows.to_numpy()[keep].astype(int)
    starts = nuc.start.to_numpy()[keep]
    op = np.sort(np.asarray(op_mid, dtype=float))
    # index of the nearest operator strictly before the rod start
    idx = np.searchsorted(op, starts, side='left') - 1
    ok = idx >= 0
    hit = np.zeros(len(rows), dtype=bool)
    hit[ok] = sites[rows[ok], idx[ok]]
    per = pd.Series(hit.astype(float)).groupby(rows).sum()
    counts.iloc[per.index.to_numpy()] = per.to_numpy()
    return counts.to_numpy()


def empirical_moments(main_file, segments_file=None, expect_K=None,
                      op_mid=None):
    """Mean feature counts from one v5 per-molecule file.

    Returns dict with the four FEATURE_STATS, N, K, and the observed per-site
    occupancy profile + n_tf histogram (for the posterior checks).
    """
    df = pd.read_table(main_file)
    sites, used_cols = _resolve_site_matrix(df, expect_K=expect_K)
    K = sites.shape[1]
    if expect_K is not None and K != expect_K:
        raise ValueError(
            f'{path.basename(main_file)}: found {K} operator columns '
            f'({used_cols}) but the positions file gives {expect_K} operators '
            'for this amplicon -- these must agree or n_pairs and the '
            'per-site occupancy check are meaningless')

    n_tf = sites.sum(axis=1) if K else np.zeros(len(df), dtype=int)
    n_pairs = ((sites[:, :-1] & sites[:, 1:]).sum(axis=1) if K > 1
               else np.zeros(len(df), dtype=int))

    if 'n_nuc' in df.columns:
        n_nuc = df['n_nuc'].to_numpy(dtype=float)
    elif segments_file is not None:
        seg = pd.read_table(segments_file)
        counts = (seg[seg['type'].str.upper() == 'NUC']
                  .groupby('read_id').size())
        n_nuc = counts.reindex(df.index if df.index.name == 'read_id'
                               else df.iloc[:, 0]).fillna(0).to_numpy(dtype=float)
    else:
        raise ValueError(f'{main_file}: no `n_nuc` column and no --segments '
                         'file given, so the nucleosome moment is unavailable')

    # delta_adj's conjugate statistic needs the segments file
    if segments_file is None:
        cand = main_file.replace('.single_molecule_classification.txt',
                                 '.segments.txt')
        segments_file = cand if path.exists(cand) else None
    if segments_file is not None and K and op_mid is not None:
        rid = (df['read_id'].tolist() if 'read_id' in df.columns
               else df.index.tolist())
        n_nuc_adj = _n_nuc_adj_from_segments(segments_file, sites, rid, op_mid)
    else:
        n_nuc_adj = np.zeros(len(df))

    return {
        'n_tf': float(n_tf.mean()),
        'n_pairs': float(n_pairs.mean()),
        'n_nuc': float(n_nuc.mean()),
        'n_nuc_and_tf': float((n_nuc * (n_tf >= 1)).mean()),
        'n_nuc_adj': float(np.mean(n_nuc_adj)),
        'N': int(len(df)),
        'K': K,
        'obs_site_occ': sites.mean(axis=0),
        'obs_n_tf_hist': np.bincount(n_tf, minlength=K + 1) / len(df),
    }


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
class Dataset:
    """One (sample, amplicon) unit: its geometry and its empirical moments."""

    def __init__(self, key, geom, moments):
        self.key = key
        self.geom = geom
        self.m = np.array([moments[s] for s in FEATURE_STATS])
        self.N = moments['N']
        self.moments = moments


def _nll_grad(theta_active, mask, datasets, method, mu_ambient, boundary,
              tie_mu=False):
    theta = expand(theta_active, mask)
    if tie_mu:
        mu_ambient = theta[FEATURES.index('mu')]
    ntot = sum(d.N for d in datasets)
    nll, grad = 0.0, np.zeros(int(mask.sum()))
    for d in datasets:
        lz = log_Z(theta, d.geom, mu_ambient, method=method, boundary=boundary)
        nll += d.N * (lz - float(np.dot(theta, d.m)))
        mm = model_moments(theta, d.geom, mu_ambient, method=method,
                           boundary=boundary)
        grad += d.N * (mm - d.m)[mask]
    return nll / ntot, grad / ntot


def fit(datasets, spec='4param', theta0=None, method='exact', verbose=True,
        mu_ambient=MU_AMBIENT_DEFAULT, boundary='stationary', tie_mu=False):
    """Joint fit across all datasets with SHARED parameters.

    Convex in theta (log Z is log-sum-exp, and a sum of convex functions is
    convex), so the optimum is unique and L-BFGS is appropriate. The
    copy-number series is what makes h, J, mu identifiable: different K give
    different <n_tf>/<n_pairs> that one shared parameter set must reproduce
    simultaneously.
    """
    mask = spec_mask(spec)
    if theta0 is None:
        start = {'h': 0.0, 'J': 0.0, 'mu': -1.0, 'delta': 0.0,
                 'delta_adj': 0.0}
        theta0 = np.array([start[f] for f in FEATURES])[mask]
    res = minimize(_nll_grad, np.asarray(theta0, dtype=float),
                   args=(mask, datasets, method, mu_ambient, boundary, tie_mu),
                   jac=True, method='L-BFGS-B')

    theta = expand(res.x, mask)
    if tie_mu:
        # One nucleosome fugacity inside and outside the window. Physically the
        # more natural model, but it costs the exact moment-matching identity:
        # d logZ/d mu now also picks up the exterior's response, so the fit is a
        # method-of-moments estimator rather than exact ML. Use the untied fit
        # as the reference and this as the parsimony check.
        mu_ambient = theta[FEATURES.index('mu')]
    ntot = sum(d.N for d in datasets)
    # total log-likelihood = sum_a N_a [ theta.<phi>_a - log Z_a ]
    loglik = sum(d.N * (float(np.dot(theta, d.m))
                        - log_Z(theta, d.geom, mu_ambient, method=method,
                                boundary=boundary))
                 for d in datasets)

    # Hessian of the TOTAL NLL = sum_a N_a Cov_model,a(phi); its inverse is the
    # covariance of theta_hat, so the standard errors come for free.
    k = int(mask.sum())
    H = np.zeros((k, k))
    e = 1e-3
    for j in range(k):
        xp, xm = res.x.copy(), res.x.copy()
        xp[j] += e
        xm[j] -= e
        _, gp = _nll_grad(xp, mask, datasets, method, mu_ambient, boundary, tie_mu)
        _, gm = _nll_grad(xm, mask, datasets, method, mu_ambient, boundary, tie_mu)
        H[:, j] = ntot * (gp - gm) / (2 * e)
    H = 0.5 * (H + H.T)
    try:
        cov = np.linalg.inv(H)
        se_active = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se_active = np.full(k, np.nan)
    se = np.full(len(FEATURES), np.nan)
    se[mask] = se_active

    out = {'spec': spec, 'theta': theta, 'se': se, 'mask': mask,
           'mu_ambient': mu_ambient, 'tie_mu': tie_mu,
           'loglik': loglik, 'n_params': k, 'n_molecules': ntot,
           'aic': -2 * loglik + 2 * k,
           'bic': -2 * loglik + k * np.log(ntot),
           'converged': bool(res.success), 'message': res.message,
           'nll_per_mol': float(res.fun)}
    if verbose:
        print(f'  {spec:16s} loglik={loglik:14.2f}  k={k}  '
              f'AIC={out["aic"]:14.2f}  {"ok" if res.success else "FAILED"}')
        for f, v, s in zip(FEATURES, theta, se):
            if f in SPECS[spec]:
                print(f'      {f:6s} {v:+8.4f} +/- {s:.4f}')
    return out


def compare(results):
    """Model-comparison table + LRT for nested pairs.

    Nestedness is structural here (a spec is a mask over one parameter vector),
    so a pair is nested iff one mask is a subset of the other -- no special
    casing, and no risk of running an LRT on a non-nested pair.
    """
    rows = []
    for r in results:
        rows.append({'spec': r['spec'], 'n_params': r['n_params'],
                     'loglik': r['loglik'], 'aic': r['aic'], 'bic': r['bic'],
                     'converged': r['converged'],
                     **{f'{f}': r['theta'][i] for i, f in enumerate(FEATURES)},
                     **{f'{f}_se': r['se'][i] for i, f in enumerate(FEATURES)}})
    tab = pd.DataFrame(rows).sort_values('aic').reset_index(drop=True)

    from scipy.stats import chi2
    lrt = []
    for a, b in itertools.combinations(results, 2):
        sa, sb = set(SPECS[a['spec']]), set(SPECS[b['spec']])
        if sa < sb:
            small, big = a, b
        elif sb < sa:
            small, big = b, a
        else:
            continue                      # not nested -> LRT invalid, skip
        dk = big['n_params'] - small['n_params']
        stat = 2 * (big['loglik'] - small['loglik'])
        lrt.append({'reduced': small['spec'], 'full': big['spec'], 'df': dk,
                    'chi2': stat,
                    'p_value': float(chi2.sf(stat, dk)) if stat > 0 else 1.0})
    return tab, pd.DataFrame(lrt)


def posterior_checks(theta, datasets, method='exact',
                     mu_ambient=MU_AMBIENT_DEFAULT, boundary='stationary'):
    """Observed vs model-predicted statistics, flagged tautological or not.

    Under a JOINT fit with shared parameters the gradient only forces the
    N-weighted SUM of each fitted moment to match, so the PER-AMPLICON fitted
    moments need not agree -- which makes them a real check. Statistics of
    features that the spec pinned to zero, and every statistic that is not a
    fitted moment at all (per-site profile, n_tf distribution), are fully
    out-of-sample.
    """
    rows = []
    active = set()
    for i, f in enumerate(FEATURES):
        if theta[i] != 0.0:
            active.add(FEATURE_STATS[i])
    for d in datasets:
        mm = model_moments(theta, d.geom, mu_ambient, method=method,
                           boundary=boundary)
        for stat, pred in zip(FEATURE_STATS, mm):
            rows.append({'check': 'moment', 'key': d.key, 'statistic': stat,
                         'observed': d.moments[stat], 'predicted': pred,
                         'n_molecules': d.N,
                         'is_fit_moment': stat in active})
        occ = per_site_occupancy(theta, d.geom, mu_ambient, method=method,
                                 boundary=boundary)
        obs = d.moments['obs_site_occ']
        for k in range(min(len(occ), len(obs))):
            rows.append({'check': 'site_occupancy', 'key': d.key,
                         'statistic': f'operator_{k + 1}',
                         'observed': float(obs[k]), 'predicted': float(occ[k]),
                         'n_molecules': d.N, 'is_fit_moment': False})
        pj = n_tf_distribution(theta, d.geom, mu_ambient,
                               boundary=boundary)
        obs_h = d.moments['obs_n_tf_hist']
        for j in range(min(len(pj), len(obs_h))):
            rows.append({'check': 'n_tf_distribution', 'key': d.key,
                         'statistic': f'n_tf_{j}',
                         'observed': float(obs_h[j]), 'predicted': float(pj[j]),
                         'n_molecules': d.N, 'is_fit_moment': False})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
def _enumerate_configs(geom):
    """Yield (rod_start_list, bound_operator_tuple) for every valid config on a
    WALLED lattice. Only tractable for tiny lattices -- which is the point: it
    gives ground truth to validate the DP's indexing against."""
    L, T = geom.L, geom.T
    starts = [i for i in range(T) if i + L - 1 < T]

    def rod_sets(avail, cur):
        yield list(cur)
        for idx, i in enumerate(avail):
            nxt = [j for j in avail[idx + 1:] if j >= i + L]
            yield from rod_sets(nxt, cur + [i])

    for rods in rod_sets(starts, []):
        covered = np.zeros(T, dtype=bool)
        for i in rods:
            covered[i:i + L] = True
        free_ops = [k for k, b in enumerate(geom.op_bin) if not covered[b]]
        for r in range(len(free_ops) + 1):
            for chosen in itertools.combinations(free_ops, r):
                yield rods, chosen


def _brute_force_features(rods, chosen, geom):
    """The features of one enumerated configuration. n_nuc counts rods started
    INSIDE the scored window, matching the DP and v5's observable.

    n_nuc_adj counts rods whose nearest PRECEDING operator is bound -- which is
    exactly what the DP's `s` bit encodes at a rod start, so this is the
    independent check that the adjacency term is wired up correctly."""
    occ = np.zeros(max(geom.K, 1), dtype=bool)
    occ[list(chosen)] = True
    n_tf = int(len(chosen))
    n_pairs = int((occ[:-1] & occ[1:]).sum()) if geom.K > 1 else 0
    n_nuc = int(sum(1 for i in rods if geom.in_window[i]))
    n_adj = 0
    for i in rods:
        prev = [k for k, b in enumerate(geom.op_bin) if b < i]
        if prev and occ[prev[-1]]:
            n_adj += 1
    return np.array([n_tf, n_pairs, n_nuc, n_nuc * (1 if n_tf >= 1 else 0),
                     n_adj], dtype=float)


def _brute_force_log_Z(theta, geom, mu_ambient):
    """Ground truth log Z by explicit enumeration (walled lattice)."""
    h, J, mu, delta, delta_adj = theta
    terms = []
    for rods, chosen in _enumerate_configs(geom):
        f = _brute_force_features(rods, chosen, geom)
        n_pad = sum(1 for i in rods if not geom.in_window[i])
        terms.append(h * f[0] + J * f[1] + mu * f[2] + delta * f[3]
                     + delta_adj * f[4] + mu_ambient * n_pad)
    return logsumexp(terms)


def _brute_force_moments(theta, geom, mu_ambient):
    """Ground truth <phi> by explicit enumeration (walled lattice)."""
    h, J, mu, delta, delta_adj = theta
    lz = _brute_force_log_Z(theta, geom, mu_ambient)
    acc = np.zeros(len(FEATURES))
    for rods, chosen in _enumerate_configs(geom):
        f = _brute_force_features(rods, chosen, geom)
        n_pad = sum(1 for i in rods if not geom.in_window[i])
        w = np.exp(h * f[0] + J * f[1] + mu * f[2] + delta * f[3]
                   + delta_adj * f[4] + mu_ambient * n_pad - lz)
        acc += w * f
    return acc


def _bulk_coverage_analytic(mu, L):
    """Analytic hard-rod coverage fraction in an infinite 1D lattice gas.

    Rod-start density is rho = d log(lambda) / d mu, so coverage = L * rho.
    Completely independent of the DP, which makes it a real check on the
    stationary boundary rather than a self-consistency test.
    """
    x = _bulk_log_lambda(mu, L)
    # Written with expm1 to avoid catastrophic cancellation: for strongly
    # negative mu, x is ~1e-18 and exp(x) rounds to exactly 1.0 in float64, so
    # (lambda - 1) underflows to 0 and 1/(lambda-1) blows up. Algebraically
    #   coverage = L / [ (L-1) + lambda/(lambda-1) ]
    # and lambda/(lambda-1) = exp(x)/expm1(x) is stable all the way down.
    return L / ((L - 1) + np.exp(x) / np.expm1(x))


def mu_from_coverage(coverage, L):
    """Invert _bulk_coverage_analytic: the mu whose bulk hard-rod coverage is
    `coverage`. This is how --mu_ambient should be chosen -- pick the number
    from a MEASURED nucleosome coverage fraction, not by intuition about mu.
    """
    if not 0.0 < coverage < 1.0:
        raise ValueError('coverage must be strictly between 0 and 1')
    from scipy.optimize import brentq
    # Solve on LOG coverage: coverage falls exponentially with mu, so the
    # linear-scale residual is hopelessly ill-conditioned at low density.
    target = np.log(coverage)
    f = lambda mu: np.log(_bulk_coverage_analytic(mu, L)) - target
    lo, hi = -60.0, 60.0
    if f(lo) > 0 or f(hi) < 0:
        raise ValueError(f'coverage {coverage} unreachable for L={L} '
                         f'(range is {_bulk_coverage_analytic(lo, L):.3e}..'
                         f'{_bulk_coverage_analytic(hi, L):.6f})')
    return brentq(f, lo, hi, xtol=1e-10, rtol=8.9e-16)


def selftest():
    ok = True

    def check(name, cond, detail=''):
        nonlocal ok
        ok = ok and bool(cond)
        print(f'  [{"PASS" if cond else "FAIL"}] {name}' +
              (f'   {detail}' if detail else ''))

    print('\n== 1. DP vs brute-force enumeration (tiny lattice, WALLED) ==')
    # boundary='wall' is the only finitely-enumerable case, so this test targets
    # the DP recursion itself; test 5 targets the stationary boundary.
    g = Geometry([12, 17, 22], scored_lo=10, scored_hi=24, L=5, pad=5, name='tiny')
    for theta in [(0.0, 0.0, 0.0, 0.0, 0.0), (1.5, 0.0, -1.0, 0.0, 0.0),
                  (1.5, 0.7, -1.0, 0.0, 0.0), (0.5, -0.3, -0.5, 0.9, 0.0),
                  (2.0, 0.5, -1.5, -1.2, 0.0),
                  # delta_adj alone, and both delta terms together
                  (1.5, 0.0, -1.0, 0.0, -0.8), (1.5, 0.4, -1.0, 0.0, 1.1),
                  (0.7, -0.2, -0.6, 0.5, -0.9), (2.0, 0.3, -1.2, -0.7, 0.6)]:
        ma = -0.7                                  # distinct from mu on purpose
        bf = _brute_force_log_Z(theta, g, ma)
        dp = log_Z(theta, g, ma, boundary='wall')
        check(f'log Z  theta={theta}', abs(bf - dp) < 1e-9,
              f'brute={bf:.10f} dp={dp:.10f}')

    print('\n== 2. moments vs brute force (incl. window-restricted n_nuc) ==')
    for theta in [(1.2, 0.4, -0.8, 0.0, 0.0), (1.2, 0.4, -0.8, 0.6, 0.0),
                  (0.3, -0.5, -0.2, -0.9, 0.0),
                  (1.2, 0.4, -0.8, 0.0, -0.7), (1.0, 0.2, -0.5, 0.4, 0.9)]:
        ma = -0.7
        bf = _brute_force_moments(theta, g, ma)
        dp = model_moments(theta, g, ma, eps=1e-5, boundary='wall')
        err = np.abs(bf - dp).max()
        check(f'<phi>  theta={theta}', err < 1e-4,
              'brute=' + np.array2string(bf, precision=6) +
              ' dp=' + np.array2string(dp, precision=6))

    print('\n== 3. the two delta implementations agree (stationary) ==')
    g2 = Geometry([300, 340, 380, 420], scored_lo=120, scored_hi=560,
                  L=147, pad=147, name='delta_xcheck')
    for theta in [(2.0, 0.5, -1.0, 0.8, 0.0), (1.0, 0.0, -0.5, -1.5, 0.0),
                  (3.0, 0.2, -2.0, 2.0, 0.0)]:
        a = _log_Z_delta_exact(np.array(theta), g2, -1.0)
        b = _log_Z_delta_identity(np.array(theta), g2, -1.0)
        check(f'exact vs identity  theta={theta}', abs(a - b) < 1e-6,
              f'exact={a:.8f} identity={b:.8f}')

    print('\n== 4. n_tf distribution is consistent with <n_tf> ==')
    for theta in [(1.2, 0.4, -0.8, 0.0, 0.0), (1.2, 0.4, -0.8, 0.6, 0.0)]:
        pj = n_tf_distribution(theta, g, -0.7)
        mean_dist = float(np.dot(np.arange(len(pj)), pj))
        mean_mom = model_moments(theta, g, -0.7, eps=1e-5)[0]
        check(f'sum P(n_tf)=1  theta={theta}', abs(pj.sum() - 1) < 1e-9)
        check(f'mean agrees    theta={theta}', abs(mean_dist - mean_mom) < 1e-4,
              f'dist={mean_dist:.6f} moment={mean_mom:.6f}')

    print('\n== 5. stationary boundary: barrier-free lattice is EXACTLY flat ==')
    # No operators, and mu == mu_ambient so the whole line is uniform. The
    # window is then an exact marginal of an infinite uniform array, so
    # occupancy must be flat across EVERY bin -- edges included. This is the
    # test the naive padded-but-walled lattice failed by a factor of 36.
    for mu in (-2.0, -1.0, 0.0):
        gf = Geometry([], scored_lo=100, scored_hi=200, L=20, pad=40, name='flat')
        occ = nuc_occupancy_profile((0.0, 0.0, mu, 0.0, 0.0), gf, mu)
        spread = float(occ.max() - occ.min())
        analytic = _bulk_coverage_analytic(mu, gf.L)
        check(f'flat  mu={mu:+.1f}', spread < 1e-6,
              f'max-min={spread:.2e} over all {gf.T} bins')
        check(f'matches analytic hard-rod coverage  mu={mu:+.1f}',
              abs(occ.mean() - analytic) < 1e-6,
              f'dp={occ.mean():.8f} analytic={analytic:.8f}')

    print('\n== 5b. mu <-> coverage map is invertible and robust ==')
    for L in (20, 147):
        for mu in (-20.0, -10.0, -6.0, -1.0, 0.0, 5.0):
            cov = _bulk_coverage_analytic(mu, L)
            back = mu_from_coverage(cov, L)
            check(f'roundtrip L={L:3d} mu={mu:+.1f}', abs(back - mu) < 1e-6,
                  f'coverage={cov:.6f} -> mu={back:+.8f}')

    print('\n== 6. middle > edge TF occupancy from a SINGLE shared h ==')
    g6 = Geometry([282, 322, 362, 402, 442, 482, 522, 562],
                  scored_lo=120, scored_hi=600, L=147, pad=147, name='opJS4_8x')
    occ = per_site_occupancy((2.0, 0.0, -1.0, 0.0, 0.0), g6, -1.0)
    check('interior operators more occupied than edges',
          occ[len(occ) // 2] > occ[0] and occ[len(occ) // 2] > occ[-1],
          'profile=' + np.array2string(occ, precision=4))

    print('\n== 7. fit round-trip (model moments in -> theta out) ==')
    for spec, true in [('2param', (1.5, 0.0, -1.0, 0.0, 0.0)),
                       ('3param_tfcoop', (1.5, 0.6, -1.0, 0.0, 0.0)),
                       ('3param_nuc', (1.5, 0.0, -1.0, 0.7, 0.0)),
                       ('4param', (1.5, 0.6, -1.0, 0.7, 0.0)),
                       ('3param_nucadj', (1.5, 0.0, -1.0, 0.0, -0.8)),
                       ('4param_adj', (1.5, 0.5, -1.0, 0.0, -0.8)),
                       ('5param', (1.5, 0.5, -1.0, 0.6, -0.7))]:
        ds = []
        for name, ops in [('4x', [282, 322, 362, 402]),
                          ('8x', [282, 322, 362, 402, 442, 482, 522, 562])]:
            gg = Geometry(ops, scored_lo=120, scored_hi=600, L=147, pad=147,
                          name=name)
            mm = model_moments(true, gg, -1.0)
            mom = {st: mm[i] for i, st in enumerate(FEATURE_STATS)}
            mom.update({'N': 5000, 'K': gg.K,
                        'obs_site_occ': per_site_occupancy(true, gg, -1.0),
                        'obs_n_tf_hist': n_tf_distribution(true, gg, -1.0)})
            ds.append(Dataset(name, gg, mom))
        r = fit(ds, spec=spec, verbose=False, mu_ambient=-1.0)
        err = max(abs(r['theta'][i] - true[i]) for i in range(4))
        check(f'{spec:16s} recovers theta', err < 5e-3,
              'got ' + np.array2string(r['theta'], precision=4) +
              ' want ' + str(true))

    print('\n== 8. posterior checks: fitted moments tautological, others not ==')
    ds = []
    true = (1.8, 0.4, -1.0, 0.0, 0.0)
    for name, ops in [('4x', [282, 322, 362, 402]),
                      ('8x', [282, 322, 362, 402, 442, 482, 522, 562])]:
        gg = Geometry(ops, scored_lo=120, scored_hi=600, L=147, pad=147, name=name)
        mm = model_moments(true, gg, -1.0)
        mom = {st: mm[i] for i, st in enumerate(FEATURE_STATS)}
        mom.update({'N': 5000, 'K': gg.K,
                    'obs_site_occ': per_site_occupancy(true, gg, -1.0),
                    'obs_n_tf_hist': n_tf_distribution(true, gg, -1.0)})
        ds.append(Dataset(name, gg, mom))
    r = fit(ds, spec='3param_tfcoop', verbose=False, mu_ambient=-1.0)
    chk = posterior_checks(r['theta'], ds, mu_ambient=-1.0)
    fitm = chk[(chk.check == 'moment') & chk.is_fit_moment]
    resid = (fitm.observed - fitm.predicted).abs().max()
    check('fitted moments reproduced (data generated from the model)',
          resid < 1e-3, f'max|obs-pred|={resid:.2e}')
    check('per-site rows flagged out-of-sample',
          not chk[chk.check == 'site_occupancy'].is_fit_moment.any())

    print('\n' + ('ALL SELF-TESTS PASSED' if ok else 'SELF-TESTS FAILED'))
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _amplicon_of(fname):
    """`{sample}.{amplicon}.single_molecule_classification.txt` -> amplicon."""
    base = path.basename(fname)
    for suf in ('.single_molecule_classification.txt',
                '.dedup.full_unclustered.matrix'):
        if base.endswith(suf):
            base = base[:-len(suf)]
            break
    parts = base.split('.')
    return parts[-1] if len(parts) > 1 else base


def cmd_columns(args):
    for f in args.main_files:
        print(f'\n=== {f}')
        if not path.exists(f):
            print('  MISSING')
            continue
        df = pd.read_table(f, nrows=5)
        print(f'  {len(df.columns)} columns: {list(df.columns)}')
        try:
            sites, used = _resolve_site_matrix(df, expect_K=args.expect_K)
            print(f'  -> resolved {sites.shape[1]} operator columns: {used}')
        except ValueError as e:
            print(f'  -> COULD NOT RESOLVE: {e}')
    return 0


def plot_checks(checks, out_pdf, title=''):
    """Render a model_checks table as a multipage PDF.

    Kept separate from the fit so an existing `.ising_model_checks.txt` can be
    replotted without refitting. Panels are labelled FIT MOMENT vs PREDICTION,
    because a fitted moment agreeing is near-tautological while a prediction
    agreeing is evidence (see the `posterior_checks` docstring).
    """
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    style = ('/oak/stanford/groups/wjg/bgrd/papers/ad_smf/'
             'activation_domain_smf_paper/ad_smf.mplstyle')
    if path.exists(style):
        plt.style.use(style)

    def copyno(key):
        amp = key.split('|')[-1]
        for tok in amp.replace('_', ' ').split():
            if tok.endswith('x') and tok[:-1].isdigit():
                return int(tok[:-1])
        return 99

    keys = sorted(checks.key.unique(), key=copyno)

    with PdfPages(out_pdf) as pdf:
        # ---- page 1: moments, observed vs predicted
        mo = checks[checks.check == 'moment']
        stats = [st for st in FEATURE_STATS if (mo.statistic == st).any()]
        fig, axs = plt.subplots(1, len(stats), figsize=(3.2 * len(stats), 3.4),
                                squeeze=False)
        for ax, st in zip(axs[0], stats):
            sub = mo[mo.statistic == st]
            fit_m = bool(sub.is_fit_moment.iloc[0])
            ax.scatter(sub.predicted, sub.observed,
                       c=[copyno(k) for k in sub.key], cmap='viridis',
                       s=40, edgecolors='k', linewidth=0.4, zorder=3)
            both = np.concatenate([sub.predicted.values, sub.observed.values])
            lim = (both.min() - 0.05 * np.ptp(both), both.max() + 0.05 * np.ptp(both))
            ax.plot(lim, lim, 'r--', lw=1, zorder=1)
            ax.set_xlim(lim); ax.set_ylim(lim); ax.set_box_aspect(1)
            ax.set_xlabel('model'); ax.set_ylabel('observed')
            ax.set_title(f'{st}\n' + ('FIT MOMENT' if fit_m else 'PREDICTION'),
                         fontsize=8,
                         color=('dimgray' if fit_m else 'crimson'))
        fig.suptitle(f'{title} -- per-amplicon moments (colour = TetO copy number)',
                     fontsize=9)
        pdf.savefig(fig); plt.close(fig)

        # ---- page 2: per-site occupancy (never fitted)
        so = checks[checks.check == 'site_occupancy']
        if len(so):
            ks = [k for k in keys if (so.key == k).any()]
            nc = min(4, len(ks)); nr = int(np.ceil(len(ks) / nc))
            fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.7 * nr),
                                    squeeze=False)
            for ax, k in zip(axs.ravel(), ks):
                sub = so[so.key == k]
                x = np.arange(1, len(sub) + 1)
                ax.plot(x, sub.observed.values, 'o-', color='k', ms=4,
                        label='observed')
                ax.plot(x, sub.predicted.values, 's--', color='crimson', ms=4,
                        label='model')
                ax.set_title(k.split('|')[-1], fontsize=7)
                ax.set_xlabel('operator (promoter-proximal = 1)', fontsize=7)
                ax.set_ylabel('P(bound)', fontsize=7)
                ax.set_ylim(0, 1)
            for ax in axs.ravel()[len(ks):]:
                ax.axis('off')
            axs[0][0].legend(fontsize=6)
            fig.suptitle(f'{title} -- per-site occupancy: OUT-OF-SAMPLE '
                         '(only total <n_tf> was fitted)', fontsize=9)
            pdf.savefig(fig); plt.close(fig)

        # ---- page 3: n_tf distribution (only the mean was fitted)
        nd = checks[checks.check == 'n_tf_distribution']
        if len(nd):
            ks = [k for k in keys if (nd.key == k).any()]
            nc = min(4, len(ks)); nr = int(np.ceil(len(ks) / nc))
            fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.7 * nr),
                                    squeeze=False)
            for ax, k in zip(axs.ravel(), ks):
                sub = nd[nd.key == k]
                x = np.arange(len(sub))
                ax.bar(x - 0.2, sub.observed.values, width=0.4, color='k',
                       label='observed')
                ax.bar(x + 0.2, sub.predicted.values, width=0.4,
                       color='crimson', label='model')
                ax.set_title(k.split('|')[-1], fontsize=7)
                ax.set_xlabel('n_tf', fontsize=7)
                ax.set_ylabel('P(n_tf)', fontsize=7)
            for ax in axs.ravel()[len(ks):]:
                ax.axis('off')
            axs[0][0].legend(fontsize=6)
            fig.suptitle(f'{title} -- n_tf distribution: OUT-OF-SAMPLE shape '
                         '(only the MEAN was fitted). Over-dispersion here '
                         '= population heterogeneity.', fontsize=9)
            pdf.savefig(fig); plt.close(fig)
    print(f'wrote {out_pdf}')


def cmd_plot(args):
    plot_checks(pd.read_table(args.checks_file), args.out,
                title=args.title or path.basename(args.checks_file))
    return 0


def cmd_mu_ambient(args):
    L = args.nuc_len
    if args.coverage is not None:
        mu = mu_from_coverage(args.coverage, L)
        print(f'coverage {args.coverage:.4f}  (L={L})  ->  --mu_ambient {mu:.4f}')
        return 0
    print(f'mu <-> bulk nucleosome coverage, L={L} (mu is a log-fugacity in kT,'
          f' NOT a log-probability)\n')
    print(f'  {"mu (kT)":>9}  {"coverage":>9}   {"isolated-site e^mu/(1+e^mu)":>28}')
    for mu in (-8, -6, -5, -4, -3, -2, -1, 0, 1, 2):
        cov = _bulk_coverage_analytic(float(mu), L)
        iso = 1.0 / (1.0 + np.exp(-mu))
        print(f'  {mu:+9.1f}  {cov:9.4f}   {iso:28.4f}')
    print('\nThe two right-hand columns differ because ~L placements can cover '
          'any given bp;\nthe DP sums that positional degeneracy, which is '
          'exactly what a naive\nper-site logistic would miss.')
    print('\nPick --mu_ambient by measuring coverage, then inverting:')
    print('  python fit_ising_model.py mu_ambient --coverage 0.80 '
          f'--nuc_len {L}')
    return 0


def cmd_fit(args):
    pos = load_positions(args.positions)
    mats = {_amplicon_of(m): m for m in (args.matrices or [])}

    datasets = []
    for f in args.main_files:
        amp = _amplicon_of(f)
        if amp not in pos:
            print(f'  [skip] {path.basename(f)}: amplicon {amp!r} not in the '
                  'positions file', file=sys.stderr)
            continue
        ops = [(lo + hi) / 2.0 for lo, hi, _ in pos[amp]]
        # NB: do NOT skip when ops is empty. The 0x amplicons have no operators
        # and so carry no TF signal, but with zero TF competition they are the
        # cleanest constraint on the nucleosome fugacity mu. K=0 is handled
        # throughout (Geometry, log_Z, per_site_occupancy, _resolve_site_matrix).
        if amp in mats:
            lo, hi = scored_window_from_matrix(mats[amp])
        elif args.scored_lo is not None and args.scored_hi is not None:
            lo, hi = args.scored_lo, args.scored_hi
        else:
            print(f'  [skip] {path.basename(f)}: no matrix for {amp} and no '
                  '--scored_lo/--scored_hi given', file=sys.stderr)
            continue
        geom = Geometry(ops, lo, hi, L=args.nuc_len, pad=args.pad, name=amp)
        try:
            mom = empirical_moments(f, args.segments_file,
                                    expect_K=geom.K, op_mid=ops)
        except ValueError as e:
            print(f'  [skip] {path.basename(f)}: {e}', file=sys.stderr)
            continue
        key = path.basename(f).split('.')[0] + '|' + amp
        datasets.append(Dataset(key, geom, mom))
        print(f'  loaded {key}: N={mom["N"]} K={geom.K} '
              f'<n_tf>={mom["n_tf"]:.3f} <n_nuc>={mom["n_nuc"]:.3f}')

    if not datasets:
        print('no usable datasets', file=sys.stderr)
        return 1
    print(f'\nfitting {len(datasets)} dataset(s), '
          f'{sum(d.N for d in datasets)} molecules\n')

    results = [fit(datasets, spec=s, method=args.method,
                   mu_ambient=args.mu_ambient, boundary=args.boundary,
                   tie_mu=args.tie_mu_ambient)
               for s in args.specs]
    tab, lrt = compare(results)
    print('\n=== model comparison (sorted by AIC) ===')
    print(tab.to_string(index=False))
    if len(lrt):
        print('\n=== likelihood-ratio tests (nested pairs only) ===')
        print(lrt.to_string(index=False))

    # --- mu_ambient sensitivity band -------------------------------------
    # mu_ambient is FIXED, not fitted, and the fit is strongly sensitive to it:
    # a 5-percentage-point error in the assumed ambient coverage swings J by
    # ~35% and mu by ~60% (they trade off -- less outside competition is
    # absorbed as more intrinsic cooperativity). So never report a point
    # estimate without this band.
    if args.mu_ambient_scan:
        lo_c, hi_c, n_c = args.mu_ambient_scan
        band = []
        for c in np.linspace(lo_c, hi_c, int(n_c)):
            ma = mu_from_coverage(c, args.nuc_len)
            for spec in args.specs:
                rr = fit(datasets, spec=spec, method=args.method,
                         mu_ambient=ma, boundary=args.boundary, verbose=False)
                # Per-site occupancy is NOT a fitted moment (only the total
                # <n_tf> is), so its discrepancy is a valid out-of-sample
                # criterion for mu_ambient -- and unlike the log-likelihood it
                # needs no partition-function normalisation, so it IS
                # comparable across mu_ambient. The signal lives almost
                # entirely in the EDGE operators, which is where competition
                # from the exterior bath bites.
                sq, npt = 0.0, 0
                for d in datasets:
                    pred = per_site_occupancy(rr['theta'], d.geom, ma,
                                              method=args.method,
                                              boundary=args.boundary)
                    obs = np.asarray(d.moments['obs_site_occ'], dtype=float)
                    k = min(len(pred), len(obs))
                    sq += float(((pred[:k] - obs[:k]) ** 2).sum())
                    npt += k
                band.append({'assumed_coverage': c, 'mu_ambient': ma,
                             'spec': spec,
                             'site_occ_rmse': np.sqrt(sq / max(npt, 1)),
                             **{f: rr['theta'][i]
                                for i, f in enumerate(FEATURES)}})
        band = pd.DataFrame(band)
        print('\n=== mu_ambient scan ===')
        print('site_occ_rmse is the out-of-sample criterion -- its MINIMUM is the')
        print('data-driven estimate of mu_ambient. (log-likelihood is deliberately')
        print('NOT reported here: logZ carries the exterior normalisation, so it is')
        print('not comparable across mu_ambient, only across specs at fixed mu.)')
        print(band.to_string(index=False))
        for spec in args.specs:
            sub = band[band.spec == spec]
            if len(sub) > 1:
                b = sub.loc[sub.site_occ_rmse.idxmin()]
                print(f'  best mu_ambient for {spec}: coverage '
                      f'{b.assumed_coverage:.4f} (mu={b.mu_ambient:+.3f}), '
                      f'rmse={b.site_occ_rmse:.5f}')
    else:
        band = None

    # --- h-heterogeneity robustness scan --------------------------------
    # The point is NOT to report sigma_h. It is to ask whether the headline
    # parameters MOVE when per-molecule heterogeneity in h is allowed -- the
    # n_tf distribution is visibly too all-or-none without it. If they do not
    # move, that is a robustness result and the homogeneous fit stands.
    if args.sigma_h_scan:
        lo_s, hi_s, n_s = args.sigma_h_scan
        spec = args.sigma_h_spec
        het = []
        for sg in np.linspace(lo_s, hi_s, int(n_s)):
            th, res = fit_het(datasets, spec, sg, args.mu_ambient,
                              n_nodes=args.sigma_h_nodes, method=args.method,
                              boundary=args.boundary)
            sq, npt = 0.0, 0
            for d in datasets:
                pj = n_tf_distribution_het(th, d.geom, args.mu_ambient, sg,
                                           args.sigma_h_nodes,
                                           boundary=args.boundary)
                obs = np.asarray(d.moments['obs_n_tf_hist'], dtype=float)
                k = min(len(pj), len(obs))
                sq += float(((pj[:k] - obs[:k]) ** 2).sum()); npt += k
            het.append({'sigma_h': sg, 'spec': spec,
                        'n_tf_dist_rms': np.sqrt(sq / max(npt, 1)),
                        'converged': bool(res.success),
                        **{f: th[i] for i, f in enumerate(FEATURES)}})
        het = pd.DataFrame(het)
        print(f'\n=== h-heterogeneity scan (spec={spec}, GMM at fixed sigma_h) ===')
        print('n_tf_dist_rms is out-of-sample; its minimum estimates sigma_h.')
        print(het.to_string(index=False))
    else:
        het = None

    best = min(results, key=lambda r: r['aic'])
    print(f'\nbest by AIC: {best["spec"]}; writing posterior checks for it')
    chk = posterior_checks(best['theta'], datasets, method=args.method,
                           mu_ambient=best['mu_ambient'],
                           boundary=args.boundary)

    os.makedirs(args.output_dir, exist_ok=True)
    p = path.join(args.output_dir, args.out_prefix)
    tab.to_csv(f'{p}.ising_fits.txt', sep='\t', index=False)
    if len(lrt):
        lrt.to_csv(f'{p}.ising_lrt.txt', sep='\t', index=False)
    chk.to_csv(f'{p}.ising_model_checks.txt', sep='\t', index=False)
    if band is not None:
        band.to_csv(f'{p}.ising_mu_ambient_band.txt', sep='\t', index=False)
    if het is not None:
        het.to_csv(f'{p}.ising_sigma_h_scan.txt', sep='\t', index=False)
    plot_checks(chk, f'{p}.ising_model_checks.pdf',
                title=f'{args.out_prefix} [{best["spec"]}]')
    print(f'wrote {p}.ising_fits.txt / .ising_model_checks.txt / '
          '.ising_model_checks.pdf')
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    sub.add_parser('selftest', help='verify the implementation, no data needed')

    pc = sub.add_parser('columns', help='show what a v5 output file contains')
    pc.add_argument('--main_files', nargs='+', required=True)
    pc.add_argument('--expect_K', type=int, default=None,
                    help='pass 0 to allow operator-free amplicons')

    pp = sub.add_parser('plot', help='render an existing model_checks table')
    pp.add_argument('--checks_file', required=True)
    pp.add_argument('--out', required=True)
    pp.add_argument('--title', default=None)

    pm = sub.add_parser('mu_ambient',
                        help='convert a measured nucleosome coverage fraction '
                             'into the --mu_ambient value to pass to `fit`')
    pm.add_argument('--coverage', type=float, default=None,
                    help='fraction of bp covered by a nucleosome in bulk '
                         '(0-1), measured away from the array and promoter')
    pm.add_argument('--nuc_len', type=int, default=147)

    pf = sub.add_parser('fit', help='fit the model')
    pf.add_argument('--main_files', nargs='+', required=True,
                    help='v5 {sample}.{amplicon}.single_molecule_classification.txt')
    pf.add_argument('--matrices', nargs='+', default=None,
                    help='matching .dedup.full_unclustered.matrix files, used '
                         'only to read the scored window from their headers')
    pf.add_argument('--segments_file', default=None,
                    help='v5 .segments.txt, only needed if the main file has '
                         'no n_nuc column')
    pf.add_argument('--positions', required=True,
                    help='BARE-operator positions file (--l_offset 0 --r_offset 0)')
    pf.add_argument('--specs', nargs='+', default=list(SPECS),
                    choices=list(SPECS))
    pf.add_argument('--nuc_len', type=int, default=147)
    pf.add_argument('--pad', type=int, default=None,
                    help='lattice pad per side; default = nuc_len (must be >=)')
    pf.add_argument('--scored_lo', type=int, default=None)
    pf.add_argument('--scored_hi', type=int, default=None)
    pf.add_argument('--method', default='exact', choices=['exact', 'identity'])
    pf.add_argument('--mu_ambient', type=float, default=MU_AMBIENT_DEFAULT,
                    help='FIXED nucleosome fugacity of the exterior bulk array. '
                         'Not fitted -- it is the single number summarising the '
                         'reservoir. Set it from the measured bulk nucleosome '
                         'occupancy; the fitted mu is the WINDOW fugacity.')
    pf.add_argument('--boundary', default='stationary',
                    choices=['stationary', 'wall'],
                    help="'wall' exists only for the enumeration test; it gives "
                         'huge spurious edge depletion. Do not fit with it.')
    pf.add_argument('--sigma_h_scan', type=float, nargs=3, default=None,
                    metavar=('LO', 'HI', 'N'),
                    help='scan per-molecule h heterogeneity (SD in kT) and '
                         'refit at each value. Robustness check on h/J/mu/'
                         'delta, e.g. 0.0 1.0 5')
    pf.add_argument('--sigma_h_spec', default='3param_nuc', choices=list(SPECS))
    pf.add_argument('--sigma_h_nodes', type=int, default=5,
                    help='Gauss-Hermite nodes (cost scales linearly)')
    pf.add_argument('--tie_mu_ambient', action='store_true',
                    help='use ONE nucleosome fugacity inside and outside the '
                         'window (mu_ambient follows the fitted mu). More '
                         'parsimonious and physically natural, but forfeits '
                         'the exact moment-matching identity -- treat as a '
                         'method-of-moments check against the untied fit.')
    pf.add_argument('--mu_ambient_scan', type=float, nargs=3, default=None,
                    metavar=('LO_COV', 'HI_COV', 'N'),
                    help='refit across this range of ASSUMED ambient coverage '
                         'and emit the resulting parameter band. Strongly '
                         'recommended: the fit is sensitive to mu_ambient '
                         '(5 points of coverage swings J ~35%%), so a point '
                         'estimate alone is misleading. e.g. 0.75 0.82 5')
    pf.add_argument('--output_dir', default='.')
    pf.add_argument('--out_prefix', required=True)

    args = ap.parse_args()
    if args.cmd == 'selftest':
        return selftest()
    if args.cmd == 'columns':
        return cmd_columns(args)
    if args.cmd == 'plot':
        return cmd_plot(args)
    if args.cmd == 'mu_ambient':
        return cmd_mu_ambient(args)
    return cmd_fit(args)


if __name__ == '__main__':
    sys.exit(main())
