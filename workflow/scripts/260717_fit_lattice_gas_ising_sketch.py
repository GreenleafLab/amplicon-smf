#!/usr/bin/env python
"""
260717_fit_lattice_gas_ising_sketch.py

SKETCH (not wired into Snakemake, not production) of the reparametrized
partition-function fit described in:
    260715_v5_downstream_design.md   (formal derivation)
    260717_ising_intuition.md        (intuition; see esp. sections 4-5)

Goal of the sketch: show the *architecture* of the new fit end-to-end ---
  (1) read empirical feature moments off v5's HARD per-molecule calls,
  (2) compute log Z of the TetO transfer-matrix x nucleosome hard-rod
      lattice-gas by an O(length) forward DP with a padded, grand-canonical
      boundary (NO hard walls at the amplicon edge),
  (3) fit the four Ising parameters by moment-matching (convex NLL), and
  (4) demonstrate the headline claim: with a SINGLE SHARED operator field h,
      the middle>edge TF-occupancy gradient EMERGES from excluded volume ---
      no per-site h_i.

Everything here is in kT units. Energy of a config sigma:
    E(sigma) = -( h * n_tf + J * n_pairs + mu * n_nuc )     [+ delta term, see below]
    weight(sigma) = exp(-E) ,  P = weight / Z
so h, J, mu are log-weights (fields/fugacity). n_tf = # bound operators,
n_pairs = # adjacent bound-operator pairs (cooperativity), n_nuc = # nucleosomes.

TF-nucleosome antagonism (delta): implemented here as HARD exclusion (a base
covered by a nucleosome cannot host a bound TF), i.e. delta -> +inf. The place
to relax it to a soft 4th parameter is flagged inline (search "SOFT delta").

Run with no args for the synthetic self-check (prints the emergent gradient).
Point --main_file / --segments_file at real v5 outputs to fit real data.
"""

import argparse
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------
class Geometry:
    """1D bin lattice for one amplicon, in bp (bin size = 1 bp for correctness).

    operator_positions: bin coords of each TetO (single-bin TF footprint here;
        widen to a span if you want the TF to exclude more than 1 bp).
    scored_lo/scored_hi: first/last scored GpC bin (the *measurement window*).
    L: nucleosome rod length (147; make a small distribution to match v5's
       duration prior if a nuc-energetics claim needs it -- see design doc 5.6).
    pad: bins added on EACH side beyond the scored window so nucleosomes may
         overhang inward from the reservoir. Must be >= L (grand-canonical
         boundary; NEVER a hard wall at the amplicon edge -- see intuition 5.4).
    """

    def __init__(self, operator_positions, scored_lo, scored_hi, L=147, pad=None):
        self.L = int(L)
        self.pad = int(pad) if pad is not None else self.L
        # shift everything so the padded lattice starts at 0
        self.offset = scored_lo - self.pad
        self.T = (scored_hi + self.pad) - self.offset + 1           # padded length
        self.operators = sorted(int(p) - self.offset for p in operator_positions)
        self.is_teto = np.zeros(self.T, dtype=bool)
        for p in self.operators:
            if 0 <= p < self.T:
                self.is_teto[p] = True
        # operator index at each teto bin, for reading per-site occupancy back out
        self.op_bin = np.array(self.operators, dtype=int)
        self.K = len(self.operators)


# ----------------------------------------------------------------------------
# log Z by forward DP  (this is the whole model)
# ----------------------------------------------------------------------------
# State = (c, s):
#   c in 0..L-1  = nucleosome coverage still owed after this bin (0 => free)
#   s in {0,1}   = whether the most recent operator was TF-bound (for coop J)
# 2*L states. Transfer is applied bin by bin, left to right.

def log_Z(theta, geom, extra_field=None):
    """theta = (h, J, mu). extra_field: optional per-operator additive field
    bump (len K) used to read off occupancies by finite difference.

    Vectorized forward DP. State V has shape (L, 2): V[c, s] = log-weight of
    configs with `c` bins of nucleosome coverage still owed after this position
    and last-operator TF-state `s` (for cooperativity J). The coverage countdown
    (c -> c-1) is a single array shift, so each position is O(L) numpy, not a
    Python double loop -- ~1000x faster, which matters because the fit calls this
    hundreds of times."""
    h, J, mu = theta
    L, T = geom.L, geom.T
    NEG = -np.inf

    op_field = np.full(geom.K, float(h))
    if extra_field is not None:
        op_field = op_field + extra_field
    op_of_bin = -np.ones(T, dtype=int)
    for k, b in enumerate(geom.op_bin):
        if 0 <= b < T:
            op_of_bin[b] = k

    V = np.full((L, 2), NEG)
    V[0, 0] = 0.0                                # start free (c=0, s=0)

    for i in range(T):
        teto = geom.is_teto[i]
        can_start = (i + L - 1 < T)
        Vn = np.full((L, 2), NEG)

        # (1) continuation: coverage c=1..L-1 -> c-1 (array shift). Covered
        #     operator is unbound => collapses the coop state s to 0 at a teto.
        if teto:
            cont = np.logaddexp(V[1:L, 0], V[1:L, 1])
            Vn[0:L - 1, 0] = np.logaddexp(Vn[0:L - 1, 0], cont)
        else:
            Vn[0:L - 1, :] = V[1:L, :]

        v0 = V[0, :]                             # free state (c==0), s in {0,1}

        # (2a) start a nucleosome here (pay fugacity mu); SOFT delta would attach
        #      here as a penalty + an extra "TF-under-rod" state bit.
        if can_start:
            if teto:
                Vn[L - 1, 0] = np.logaddexp(Vn[L - 1, 0],
                                            np.logaddexp(v0[0], v0[1]) + mu)
            else:
                Vn[L - 1, :] = np.logaddexp(Vn[L - 1, :], v0 + mu)

        # (2b) no nucleosome start at i
        if teto:
            kf = op_field[op_of_bin[i]]
            bound = np.logaddexp(v0[0] + kf, v0[1] + kf + J)   # -> (0, s=1)
            Vn[0, 1] = np.logaddexp(Vn[0, 1], bound)
            Vn[0, 0] = np.logaddexp(Vn[0, 0], np.logaddexp(v0[0], v0[1]))  # unbound
        else:
            Vn[0, :] = np.logaddexp(Vn[0, :], v0)

        V = Vn

    # nucleosome starts were only allowed if they fit, so terminal c==0 always.
    return logsumexp(V[0, :])


# ----------------------------------------------------------------------------
# Model moments and per-site occupancy  (by finite difference of log Z)
# ----------------------------------------------------------------------------
# <phi_k>_model = d logZ / d theta_k   (exponential-family identity).
def model_moments(theta, geom, eps=1e-4):
    g = np.zeros(3)
    for k in range(3):
        tp = np.array(theta, float); tp[k] += eps
        tm = np.array(theta, float); tm[k] -= eps
        g[k] = (log_Z(tp, geom) - log_Z(tm, geom)) / (2 * eps)
    return {"n_tf": g[0], "n_pairs": g[1], "n_nuc": g[2]}


def per_site_occupancy(theta, geom, eps=1e-4):
    """P(operator k bound) = d logZ / d h_k, via a per-site field bump."""
    occ = np.zeros(geom.K)
    for k in range(geom.K):
        bump = np.zeros(geom.K); bump[k] = eps
        lp = log_Z(theta, geom, extra_field=bump)
        lm = log_Z(theta, geom, extra_field=-bump)
        occ[k] = (lp - lm) / (2 * eps)
    return occ


# ----------------------------------------------------------------------------
# Promoter node (optional extension -- see intuition doc Section 7)
# ----------------------------------------------------------------------------
# !! SUPERSEDED for the potency question. Published Fig. 3 shows (a) the
# !! occupancy->activity link is line-like, not convex, so it is NOT a Boltzmann
# !! coupling, and (b) activity depends on the AVAILABLE site count as well as
# !! the bound count (time integration), which no single-snapshot equilibrium
# !! model can represent. See 260902_two_timescale_potency.md. Kept here because
# !! the field-shift identity below is still the right machinery for any genuine
# !! molecule-level binary coupling (e.g. delta_remodel), and because the w=0
# !! limit of the two-timescale model is exactly this.
#
# One binary sigma_prom in {closed, open}, field h_prom, coupled to the TF
# configuration by J_prom * f(sigma_TF). The coupling enters as a per-operator
# field bump (the "shape" below); conditioning the promoter open == boosting
# each operator's field by J_prom * shape_i (the field-shift identity, 7.4):
#
#   Z = Z_array(theta)  +  exp(h_prom) * Z_array(theta; h_i -> h_i + J_prom*shape_i)
#
# so no new DP -- two calls to the same array log_Z, reusing extra_field.
def coupling_shape(geom, mode="count", lam=20.0, side="low"):
    """Fixed per-operator shape (linear in J_prom).

    `side` says which end of the amplicon the promoter is on, and therefore
    which operator is promoter-PROXIMAL:
      "low"  -- promoter at low coordinates  => operators[0]  is proximal
      "high" -- promoter at high coordinates => operators[-1] is proximal

    Default is "low", which is correct for opJS4/opJS5: the minimal promoter
    sits at ~136-195 while the TetO array runs 282 -> 581, and copy number adds
    DISTAL sites without moving existing ones (all variants are 630 bp with
    shared absolute operator coordinates). Getting this backwards silently
    inverts `nearest` and `decay` -- it puts all the coupling weight on the
    most distal operator.
    """
    K = geom.K
    if mode == "count":                                   # dosage: all equal
        return np.ones(K)
    if side not in ("low", "high"):
        raise ValueError(f"side must be 'low' or 'high', got {side!r}")
    if mode == "nearest":                                 # only proximal operator
        v = np.zeros(K)
        v[0 if side == "low" else -1] = 1.0
        return v
    if mode == "decay":                                   # range lambda (bp)
        # distance from the proximal operator; 0 at the proximal op
        d = (geom.op_bin - geom.op_bin.min()) if side == "low" \
            else (geom.op_bin.max() - geom.op_bin)
        return np.exp(-d / lam)
    raise ValueError(mode)


def log_Z_promoter(params5, geom, mode="count", lam=20.0, side="low"):
    """params5 = (h, J, mu, h_prom, J_prom). `side` per coupling_shape()."""
    h, J, mu, h_prom, J_prom = params5
    theta = (h, J, mu)
    shape = coupling_shape(geom, mode, lam, side)
    lz_closed = log_Z(theta, geom)
    lz_open = h_prom + log_Z(theta, geom, extra_field=J_prom * shape)
    return np.logaddexp(lz_closed, lz_open)


def promoter_occupancy(params5, geom, mode="count", lam=20.0, side="low"):
    """P(promoter open)."""
    h, J, mu, h_prom, J_prom = params5
    shape = coupling_shape(geom, mode, lam, side)
    lz_open = h_prom + log_Z((h, J, mu), geom, extra_field=J_prom * shape)
    return float(np.exp(lz_open - log_Z_promoter(params5, geom, mode, lam, side)))


def model_moments_promoter(params5, geom, mode="count", lam=20.0, eps=1e-4,
                           side="low"):
    """Five conjugate moments, by finite difference of log_Z_promoter:
       d/dh->n_tf, d/dJ->n_pairs, d/dmu->n_nuc, d/dh_prom-><sigma_prom>,
       d/dJ_prom-><f(sigma_TF)*sigma_prom> (the coupling feature = 'potency' obs)."""
    names = ["n_tf", "n_pairs", "n_nuc", "sigma_prom", "coupling"]
    out = {}
    for k, nm in enumerate(names):
        pp = np.array(params5, float); pp[k] += eps
        pm = np.array(params5, float); pm[k] -= eps
        out[nm] = (log_Z_promoter(pp, geom, mode, lam, side)
                   - log_Z_promoter(pm, geom, mode, lam, side)) / (2 * eps)
    return out


def fit_promoter(data_moments, geom, mode="count", lam=20.0,
                 theta0=(0.0, 0.0, -1.0, 0.0, 0.0), side="low"):
    """Fit 5 params by moment matching. data_moments: dict with keys
    n_tf, n_pairs, n_nuc, sigma_prom, coupling (all per-molecule means)."""
    m = np.array([data_moments[k] for k in
                  ["n_tf", "n_pairs", "n_nuc", "sigma_prom", "coupling"]])

    def nll(p):
        return log_Z_promoter(p, geom, mode, lam, side) - float(np.dot(p, m))

    def grad(p):
        mm = model_moments_promoter(p, geom, mode, lam, side=side)
        return np.array([mm[k] for k in
                         ["n_tf", "n_pairs", "n_nuc", "sigma_prom", "coupling"]]) - m

    return minimize(nll, np.array(theta0, float), jac=grad, method="L-BFGS-B")


# ----------------------------------------------------------------------------
# Empirical moments from v5 HARD calls
# ----------------------------------------------------------------------------
def empirical_moments_from_v5(main_file):
    """Read v5 single_molecule_classification.txt: tfbs_1..K bool + n_nuc.
    Returns (<n_tf>, <n_pairs>, <n_nuc>, N). n_pairs uses operator column order."""
    import pandas as pd
    df = pd.read_table(main_file)
    tf_cols = sorted([c for c in df.columns if c.startswith("tfbs_")],
                     key=lambda c: int(c.split("_")[1]))
    T = df[tf_cols].to_numpy().astype(bool)                 # (N, K)
    n_tf = T.sum(1)
    n_pairs = (T[:, :-1] & T[:, 1:]).sum(1)                 # adjacent co-bound
    n_nuc = df["n_nuc"].to_numpy() if "n_nuc" in df else np.zeros(len(df))
    return dict(n_tf=n_tf.mean(), n_pairs=n_pairs.mean(),
                n_nuc=float(np.mean(n_nuc)), N=len(df))


# ----------------------------------------------------------------------------
# Fit: minimize convex NLL   NLL(theta)/N = logZ - sum_k theta_k * <phi_k>_data
# gradient = <phi>_model - <phi>_data  (moment matching at the optimum)
# ----------------------------------------------------------------------------
def fit(data_moments, geom, theta0=(0.0, 0.0, -1.0)):
    m = np.array([data_moments["n_tf"], data_moments["n_pairs"],
                  data_moments["n_nuc"]])

    def nll(theta):
        return log_Z(theta, geom) - float(np.dot(theta, m))

    def grad(theta):
        mm = model_moments(theta, geom)
        return np.array([mm["n_tf"], mm["n_pairs"], mm["n_nuc"]]) - m

    res = minimize(nll, np.array(theta0, float), jac=grad, method="L-BFGS-B")
    # Hessian = Cov_model(phi) -> error bars (free, per intuition 3.3). Here we
    # approximate it by finite-differencing the gradient at the optimum.
    H = np.zeros((3, 3)); e = 1e-3
    for k in range(3):
        tp = res.x.copy(); tp[k] += e
        tm = res.x.copy(); tm[k] -= e
        H[:, k] = (grad(tp) - grad(tm)) / (2 * e)
    cov = np.linalg.pinv(0.5 * (H + H.T))
    se = np.sqrt(np.clip(np.diag(cov), 0, None)) / max(data_moments.get("N", 1), 1) ** 0.5
    return res.x, se, res


# ----------------------------------------------------------------------------
# Synthetic self-check: uniform h, does the middle>edge gradient emerge?
# ----------------------------------------------------------------------------
def _synthetic_demo():
    # 6 operators, evenly spaced, inside a scored window; nucleosomes L=147.
    ops = [300, 321, 342, 363, 384, 405]          # 6x TetO, 21 bp spacing
    geom = Geometry(ops, scored_lo=120, scored_hi=560, L=147, pad=147)
    theta = (2.0, 0.5, -1.0)                       # h, J, mu (illustrative)

    occ = per_site_occupancy(theta, geom)
    mm = model_moments(theta, geom)
    print("== synthetic self-check (single shared h = %.2f) ==" % theta[0])
    print("per-operator TF occupancy (edge..middle..edge):")
    for k, o in enumerate(occ):
        print("   op %d: %.3f" % (k + 1, o))
    print("--> interior > edge from excluded volume alone, "
          "with NO per-site h_i.\n")
    print("model moments:", {k: round(v, 3) for k, v in mm.items()})

    # round-trip: treat model moments as 'data' and recover theta
    data = dict(mm); data["N"] = 5000
    hat, se, _ = fit(data, geom, theta0=(0.0, 0.0, -1.5))
    print("\nrecovered theta (should match %.2f, %.2f, %.2f):" % theta)
    print("   h=%.3f  J=%.3f  mu=%.3f" % tuple(hat))

    # -- promoter node (Section 7): potency = the coupling J_prom --
    print("\n== promoter node: potency as a coupling (count model) ==")
    print("J_prom  P(open)  <n_tf*sigma_prom>   <n_tf|open>-<n_tf|closed>")
    for Jp in (0.0, 0.5, 1.0, 2.0):
        p5 = (2.0, 0.5, -1.0, -1.0, Jp)           # h,J,mu,h_prom,J_prom
        popen = promoter_occupancy(p5, geom)
        mmp = model_moments_promoter(p5, geom)
        # mutual stabilization: mean TFs conditional on promoter state
        ntf_open = model_moments(  # <n_tf> with all fields boosted (open branch)
            (2.0 + Jp, 0.5, -1.0), geom)["n_tf"]
        ntf_closed = model_moments((2.0, 0.5, -1.0), geom)["n_tf"]
        print("  %.1f     %.3f      %.3f              +%.3f"
              % (Jp, popen, mmp["coupling"], ntf_open - ntf_closed))
    print("--> larger J_prom => promoter opens more AND draws more TF binding "
          "(one symmetric coupling).")

    # round-trip the 5-param promoter fit
    true5 = (2.0, 0.5, -1.0, -1.0, 1.0)
    dmom = model_moments_promoter(true5, geom)
    res5 = fit_promoter(dmom, geom, theta0=(0.0, 0.0, -1.5, 0.0, 0.0))
    print("\nrecovered 5 params (match 2.0,0.5,-1.0,-1.0,1.0):")
    print("   h=%.3f J=%.3f mu=%.3f h_prom=%.3f J_prom=%.3f" % tuple(res5.x))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main_file", help="v5 single_molecule_classification.txt")
    ap.add_argument("--operators", type=int, nargs="+",
                    help="operator bin positions (model frame)")
    ap.add_argument("--scored_lo", type=int)
    ap.add_argument("--scored_hi", type=int)
    ap.add_argument("--L", type=int, default=147)
    args = ap.parse_args()

    if not args.main_file:
        _synthetic_demo()
        return

    geom = Geometry(args.operators, args.scored_lo, args.scored_hi, L=args.L)
    data = empirical_moments_from_v5(args.main_file)
    print("empirical moments:", {k: round(v, 4) for k, v in data.items()})
    hat, se, res = fit(data, geom)
    names = ["h (TetO field)", "J (cooperativity)", "mu (nuc fugacity)"]
    print("\nfit (kT units):")
    for n, v, s in zip(names, hat, se):
        print("   %-20s %+.3f +/- %.3f" % (n, v, s))
    print("\npredicted per-site occupancy:", np.round(per_site_occupancy(hat, geom), 3))


if __name__ == "__main__":
    main()
