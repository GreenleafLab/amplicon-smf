#!/usr/bin/env python
"""
classify_single_molecule_binding_v5_hsmm.py

Hidden semi-Markov model (HSMM) segmentation of single-molecule SMF footprints.

This is the clean-room v5 alternative to classify_single_molecule_binding_v4.py, built from
SPEC_hmm_footprint_model.md. Instead of enumerating every global microstate (nuc x TF combos),
pruning overlaps, collapsing, and bolting on boundary/run-length penalties, it decodes each
molecule independently with a segmental (semi-Markov) Viterbi pass over base-pair space. Every
rule v4 hand-coded as a penalty (TF boundaries, run-length Occam cap, parsimony, di-nucleosome
handling) instead falls out of emissions + duration priors + transition costs.

Latent segment types (states):
  OPEN  - accessible linker. Anywhere. Geometric-ish (flat) duration. Low protection prob
          (prob_unmeth_given_open), with an optional promoter window override.
  NUC   - nucleosome. Dyad at the segment center; POSITION-DEPENDENT emission via a
          logistic-distance sigmoid whose half-max is TIED to the segment edge (footprint
          FWHM == segment duration), with nuc_softness a fixed edge-taper ("breathing") width.
          Extent/size variability is carried ENTIRELY by the duration prior ~ Gaussian peaked at
          nuc_mode (~147bp), allowing ~nuc_min..nuc_max. Adjacent nucleosomes = two NUC
          segments back-to-back (NUC->NUC transition), NOT one long NUC; whether an OPEN linker
          sits between them is decided by the linker GpCs' emissions.
  TF    - bound TF. Only allowed spanning an annotated TFBS motif (motif +/- tf_margin). Fixed
          duration = that motif window. Flat high protection (prob_unmeth_given_tf).
  UNID  - unidentified footprint (NEW capability). A protected segment that is NOT at a motif
          AND NOT nucleosome-length (duration in [unid_min,unid_max], and forbidden to overlap
          any motif). Flat high protection. Carries a higher start cost so it is the explanation
          of last resort, only winning when neither NUC nor TF fits.

Emissions pass through the same conversion model as v4:
  p(obs=protected | state,j) = p_U(state,j) * p_t_given_unmeth + (1-p_U) * p_t_given_meth
NaN / no-info GpCs (matrix value -1) contribute nothing to any segment's emission.

Decoding is O(G^2 * S) grid-segments per molecule (G = boundary grid points, S = 4), vectorized
across all molecules at once. No global state enumeration, no pruning, no collapse.

Outputs (new clean schema -- deliberately NOT coerced back to v4's bin-based columns):
  1. <output>                : one row per molecule -- read_id, tfbs_{k} binary calls (legacy,
                               positional), site_{name} binary calls (stable, from the positions
                               file's name column), n_tf / n_tf_teto / n_tf_other / n_nuc / n_unid,
                               semicolon-joined nuc and unid interval lists, log_likelihood.
                               ** Use n_tf_teto, not n_tf, for "how many TetOs are bound" --
                               n_tf counts every named site and so grows as promoter footprints
                               (TATA/BRE/Inr/...) are added to the positions file. Same trap
                               applies to `df.filter(like='tfbs_').sum(axis=1)`. **
  2. <output>.segments.txt   : tidy long-format Viterbi path -- one row per segment
                               (read_id, seg_index, type, start, end, dyad_or_motif, motif_name).
  3. <output>.tfbs_index.txt : sidecar index -> (name, lo, hi, is_teto) map, so an old file's
                               positional tfbs_{k} columns stay interpretable after the
                               positions file is edited.
  4. <plot>                  : per-read traces with NUC (gray) / TF (red) / UNID (purple) drawn.

With --regions <build_promoter_regions.py file>, the main file also gains REGION-ANCHORED
per-molecule columns -- the only ones that know what a TSS or promoter is:
  tss_nuc / tss_footprint / tss_open    state of the TSS column itself (NUC / UNID|TF / neither)
  promoter_nuc_gt50                     >50% of the promoter INSERT covered by NUC
  nuc_bases_promoter, frac_nuc_promoter, promoter_len   (carry these -- see caveat)
  {region}_bases_{type}, _frac_, _any_  for every region x segment type (plus1_nuc, TATA, ...)
** Promoter insert lengths are NOT uniform (264 bp for most opoBD9 promoters, RPS9 200,
minCMV 59), so `promoter_nuc_gt50` is NOT comparable across promoters of different length --
at 59 bp the threshold is 30 bp, which almost any grazing nucleosome clears. Use
--promoter_nuc_min_bp for an absolute-bp column when comparing across promoters. **
Computed via annotate_molecule_regions.annotate_one() -- the same function the post-hoc
annotator uses, on the same tidy segments table, so inline and post-hoc agree exactly.

EM parameter fitting is available behind --do_em (Viterbi / hard EM). By default it fits only the
safe set (prob_unmeth_given_tf, prob_unmeth_given_unid, nuc_mode); structural costs, conversion
rates and nuc_softness stay fixed. prob_unmeth_given_open and nuc_sigma are unstable under hard EM
and must be opted into explicitly via --em_fit. EM is OFF by default -- the shipped defaults are
hand-calibrated and work as-is.
"""

import argparse
import os
import os.path
import re
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace as dc_replace
from typing import List, Tuple

import numpy as np
import pandas as pd

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Patch, Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.transforms import blended_transform_factory
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')

from scipy.special import expit
from scipy.optimize import minimize

warnings.filterwarnings('ignore')

from common import (
    load_single_molecule_matrix,
    get_methyl_positions,
    load_tfbs_positions,
    filter_all_converted_reads,
    plot_single_read,
    adjust_gcgs,
)

# State indices
OPEN, NUC, TF, UNID = 0, 1, 2, 3
STATE_NAMES = {OPEN: 'OPEN', NUC: 'NUC', TF: 'TF', UNID: 'UNID'}
N_STATES = 4
NEG_INF = -1e18


# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------

@dataclass
class ModelParams:
    amp_width: int

    # --- decoding grid ---
    # Segment boundaries are placed on a bp grid of this spacing (plus every motif edge is added
    # exactly). Finer = more precise dyad/boundary placement, slower. Dyad resolution ~ grid_step/2.
    grid_step: int = 5

    # --- conversion model (shared by all emissions; same meaning as v4) ---
    p_t_given_unmeth: float = 0.95   # P(observe protected | truly protected)
    p_t_given_meth: float = 0.15     # P(observe protected | truly accessible)

    # --- per-state "truly protected" probabilities (p_U) ---
    prob_unmeth_given_open: float = 0.05
    prob_unmeth_given_tf: float = 0.9
    prob_unmeth_given_unid: float = 0.9
    # promoter window gets a different open prob (chromatin at the min promoter is more open)
    promoter_lo: int = 0
    promoter_hi: int = 0
    prob_unmeth_given_open_promoter: float = 0.5

    # --- nucleosome protection footprint (logistic-in-distance taper) ---
    # The footprint half-max is TIED to the segment half-length (see nuc_emission): the sigmoid
    # crosses p_U=0.5 at the segment edge, so FWHM == segment duration. Extent/size variability is
    # therefore carried ENTIRELY by the duration prior (nuc_mode/nuc_sigma); nuc_softness is a fixed
    # mechanistic "breathing" constant (edge transition width ~= 2*softness bp), NOT a free width
    # parameter and NOT EM-fit. This removes the old d_edge<->nuc_sigma degeneracy (two knobs both
    # encoding single-nucleosome extent) and the "footprint narrower than segment -> spurious tail
    # segment" bug. See NOTES_v5_calibration_log.md 2026-07-13.
    nuc_softness: float = 5.0

    # --- nucleosome duration prior (Gaussian, log-unnormalized) ---
    nuc_mode: float = 147.0
    nuc_sigma: float = 25.0
    nuc_min: int = 100
    nuc_max: int = 210
    # Minimum number of OBSERVED-protected GpCs a NUC segment must contain to be callable (per
    # molecule). Evidence floor -- prevents a nucleosome being called on just 1-2 protected GpCs.
    # Matters most at the amplicon edges: the ghost padding lets a NUC hide its unsupported length
    # off-screen, so without this floor a +1 nuc could be called from a single protected edge GpC
    # (a spurious call). Applies to ALL nucs; on these GpC-dense amplicons it essentially only ever
    # bites the edge case. Set 0 to disable. Mirrors the "needs a few protected GpCs" intuition on
    # the TF/TetO side.
    nuc_min_prot_gpcs: int = 3

    # --- TF geometry ---
    # bp added on each side of a motif (matches v4's -2/+2). Besides the intended flanking-GpC
    # capture, this ALSO silently absorbs a known 1bp right-flank off-by-one in the positions file
    # (the operator's right-flank GpC lands on the half-open, exclusive `end` and is dropped by the
    # file; tf_margin>=2 re-includes it). Keep tf_margin>=2 unless positions files are regenerated
    # with bare-operator coords AND re-verified. See NOTES_v5_hsmm_status.md (2026-07-13).
    tf_margin: int = 2

    # --- UNID geometry (off-motif, non-nucleosome length) ---
    unid_min: int = 15
    unid_max: int = 90               # <= nuc_min so a UNID can never be nucleosome-length

    # --- UNID / motif exclusivity (2026-08-12) ---
    # Historically a UNID could not overlap ANY positions-file motif (hard `continue` in the DP).
    # That ban is REDUNDANT for its apparent purpose: prob_unmeth_given_tf == prob_unmeth_given_unid
    # (both 0.9), so a UNID whose span coincides with a motif has an IDENTICAL emission to the TF
    # segment and costs start_cost_unid - start_cost_tf = 5.0 nats more. TF therefore wins every tie
    # by construction; no constraint is needed to stop UNID "stealing" a named site.
    # What the ban actually forbids is EXTENT REVISION -- a UNID spanning the motif *plus more* --
    # i.e. it hard-codes "the annotation is exactly right about the footprint's width". And because
    # unid_min=15, a modest extension cannot be expressed as a separate abutting UNID either (a 10 bp
    # flank is below the minimum), so with the ban on, the model literally cannot say "the footprint
    # is this motif and then a bit more".
    # Setting unid_may_overlap_motifs=True converts that hard constraint into a soft ~5-nat prior
    # (~2.25 protected GpCs at 2.22 nats/GpC): a UNID only beats the named site by demonstrating
    # that much extra protection outside the motif. Well-posed likelihood-ratio test on extent.
    # ⚠ Scope it. TetO array sites stay exclusive by default (unid_overlap_teto=False) -- a wide UNID
    # straddling two operators would muddy n_tf_teto, the "how many TetOs bound" readout.
    # DEFAULT FLIPPED TO True 2026-08-12 (bgrd). Rounds before this date, including attempt2 and the
    # attempt3 JUNB r3 / RPS9 r2 decodes, ran with the ban ON -- do not compare across the flip
    # without re-running. Pass --no_unid_may_overlap_motifs to restore the old behavior.
    unid_may_overlap_motifs: bool = True
    unid_overlap_teto: bool = False
    teto_name_prefix: str = 'TetO'   # sites named with this prefix stay UNID-exclusive

    # --- structural costs (all POSITIVE; subtracted as log-penalties) ---
    # per-segment "start" costs => parsimony (higher = fewer of that segment type)
    start_cost_open: float = 0.5
    start_cost_nuc: float = 2.0
    # 1.0 (was 2.0): at 2.0 a nucleosome anchored to adjacent protection would "swallow" bound TFs
    # (the linker-accessible spots cost less than the extra footprint segments needed to call the
    # TFs separately). 1.0 recovers those TFs without cannibalizing real nucleosomes -- see
    # NOTES_v5_calibration_log.md, edge case 2026-07-10a.
    start_cost_tf: float = 1.0
    start_cost_unid: float = 6.0     # last resort
    # reduced UNID start cost applied ONLY where a bulk-discovered footprint says one is expected
    # (empirical-Bayes / position-specific prior; see discover_footprints_from_bulk + the
    # `discovered_footprints` arg to viterbi_decode). Borrows strength across molecules: a site with
    # recurrent bulk protection becomes a data-derived "soft annotation", sitting between a known
    # motif (start_cost_tf=1.0) and a truly-unexpected footprint (start_cost_unid=6.0). MOCKUP
    # 2026-07-13 -- default equals start_cost_unid so behavior is unchanged unless intervals passed.
    unid_discovered_start_cost: float = 6.0
    # transition costs between adjacent segment types
    trans_nuc_nuc: float = 1.0       # di-nucleosome (adjacent nucs, no linker)
    trans_tf_tf: float = 1.0         # directly adjacent TFs (no linker)
    trans_nuc_tf: float = 1.0        # nuc directly against a TF (calibrate via streak histogram)
    trans_unid_adj: float = 1.0      # UNID directly against a nuc or TF

    def __post_init__(self):
        if self.unid_max > self.nuc_min:
            # keep UNID strictly non-nucleosome-length
            self.unid_max = self.nuc_min


# ---------------------------------------------------------------------------
# Grid + geometry
# ---------------------------------------------------------------------------

def build_boundary_grid(amp_width, grid_step, tfbs_positions, tf_margin):
    """
    Boundary grid = regular bp lattice UNION every motif's [start-margin, end+margin] edges, so TF
    segments align to motifs exactly. Returns a sorted, unique int array of bp positions in [0, amp_width].
    """
    pts = set(range(0, amp_width + 1, grid_step))
    pts.add(0)
    pts.add(amp_width)
    for tfbs in tfbs_positions:
        s = int(tfbs[0]) - tf_margin
        e = int(tfbs[1]) + tf_margin
        for p in (s, e):
            if 0 <= p <= amp_width:
                pts.add(p)
    return np.array(sorted(pts), dtype=np.int64)


def motif_segments_on_grid(bnd, tfbs_positions, tf_margin):
    """
    For each motif, the (start_boundary_idx, end_boundary_idx, motif_k) of its TF segment.
    Motif edges were added to the grid so these align exactly.
    """
    segs = []
    for k, tfbs in enumerate(tfbs_positions):
        s = int(tfbs[0]) - tf_margin
        e = int(tfbs[1]) + tf_margin
        ai = int(np.searchsorted(bnd, s))
        bi = int(np.searchsorted(bnd, e))
        if ai < bi and bnd[ai] == s and bnd[bi] == e:
            segs.append((ai, bi, k))
    return segs


def motif_intervals(tfbs_positions, tf_margin):
    """All motif (lo, hi) bp intervals, margin-padded. Defines where a TF segment is LEGAL."""
    return [(int(t[0]) - tf_margin, int(t[1]) + tf_margin) for t in tfbs_positions]


def unid_blocked_intervals(tfbs_positions, tf_margin, params):
    """
    Subset of motif_intervals() that a UNID segment may NOT overlap. See the
    unid_may_overlap_motifs block in ModelParams for why this is a subset and not all of them.

      unid_may_overlap_motifs=False (default)  -> every motif blocks (historical behavior)
      True, unid_overlap_teto=False            -> only TetO-array sites block; promoter sites are
                                                  open to extent revision
      True, unid_overlap_teto=True             -> nothing blocks

    Name matching uses the same case-insensitive prefix rule as is_teto_site(), so an unnamed
    positions-file entry (no name column) is treated as a promoter site, not an array site.
    """
    if not params.unid_may_overlap_motifs:
        return motif_intervals(tfbs_positions, tf_margin)
    if params.unid_overlap_teto:
        return []
    out = []
    for t in tfbs_positions:
        name = str(t[2]).strip() if (len(t) > 2 and t[2] is not None) else ''
        if is_teto_site(name, params.teto_name_prefix):
            out.append((int(t[0]) - tf_margin, int(t[1]) + tf_margin))
    return out


# ---------------------------------------------------------------------------
# Emission precompute
# ---------------------------------------------------------------------------

def _logw(p_u, p_tu, p_tm):
    """Return (log P(obs=1|p_u), log P(obs=0|p_u)) for a scalar or array p_u."""
    p_t = np.clip(p_u * p_tu + (1.0 - p_u) * p_tm, 1e-12, 1.0 - 1e-12)
    return np.log(p_t), np.log(1.0 - p_t)


def flat_state_pu(state, gpc_pos, params):
    """Per-GpC 'truly protected' probability p_U for a FLAT state (OPEN/TF/UNID)."""
    n_gpc = len(gpc_pos)
    if state == OPEN:
        pu = np.full(n_gpc, params.prob_unmeth_given_open)
        if params.promoter_lo < params.promoter_hi:
            in_prom = (gpc_pos >= params.promoter_lo) & (gpc_pos <= params.promoter_hi)
            pu[in_prom] = params.prob_unmeth_given_open_promoter
        return pu
    if state == TF:
        return np.full(n_gpc, params.prob_unmeth_given_tf)
    if state == UNID:
        return np.full(n_gpc, params.prob_unmeth_given_unid)
    raise ValueError('flat_state_pu called on non-flat state {}'.format(state))


def build_flat_cumulatives(gpc_pos, obs1, obs0, params):
    """
    For the flat-emission states (OPEN with promoter override, TF, UNID) precompute per-state
    cumulative-over-GpC contribution arrays so any window's emission is an O(1) difference.

    contribution_s[j, m] = logw1_s[j]*obs1[j,m] + logw0_s[j]*obs0[j,m]   (0 where obs is missing)
    C_s = cumsum over j (leading zero row) -> (n_gpc+1, n_mol)
    emission over GpC index window [jlo, jhi) = C_s[jhi] - C_s[jlo]

    Missing GpCs (obs1==obs0==0) contribute nothing, which is correct per SPEC.
    Returns dict {OPEN|TF|UNID: C_s}.
    """
    p_tu, p_tm = params.p_t_given_unmeth, params.p_t_given_meth
    n_gpc = len(gpc_pos)

    # OPEN p_U per GpC (promoter override)
    open_pu = np.full(n_gpc, params.prob_unmeth_given_open)
    if params.promoter_lo < params.promoter_hi:
        in_prom = (gpc_pos >= params.promoter_lo) & (gpc_pos <= params.promoter_hi)
        open_pu[in_prom] = params.prob_unmeth_given_open_promoter

    cumulatives = {}
    for state, pu in (
        (OPEN, open_pu),
        (TF, np.full(n_gpc, params.prob_unmeth_given_tf)),
        (UNID, np.full(n_gpc, params.prob_unmeth_given_unid)),
    ):
        lw1, lw0 = _logw(pu, p_tu, p_tm)                       # (n_gpc,)
        contrib = lw1[:, None] * obs1 + lw0[:, None] * obs0    # (n_gpc, n_mol)
        C = np.zeros((n_gpc + 1, contrib.shape[1]))
        np.cumsum(contrib, axis=0, out=C[1:])
        cumulatives[state] = C
    return cumulatives


def nuc_emission(jlo, jhi, dyad, half_len, gpc_pos, obs1, obs0, params):
    """
    On-the-fly NUC emission for the GpC index window [jlo, jhi) given a dyad bp position and the
    segment half-length half_len (= (end-start)/2).
    p_U(j) = sigmoid((half_len - |pos_j - dyad|) / softness): position-dependent, with the half-max
    TIED to the segment edge (dist == half_len -> p_U = 0.5). So the footprint's full-width-at-
    half-max equals the segment duration; softness is the fixed edge-taper ("breathing") width.
    Returns an (n_mol,) vector. Missing GpCs contribute nothing (obs1==obs0==0 there).
    """
    n_mol = obs1.shape[1]
    if jhi <= jlo:
        return np.zeros(n_mol)
    dist = np.abs(gpc_pos[jlo:jhi].astype(float) - dyad)
    p_u = expit((half_len - dist) / params.nuc_softness)
    lw1, lw0 = _logw(p_u, params.p_t_given_unmeth, params.p_t_given_meth)   # (w,)
    return lw1 @ obs1[jlo:jhi] + lw0 @ obs0[jlo:jhi]                        # (n_mol,)


# ---------------------------------------------------------------------------
# Structural (transition / duration / start) costs
# ---------------------------------------------------------------------------

def build_transition_matrix(params):
    """
    (S, S) additive log-score matrix trans[s_prev, s_new]. Disallowed transitions = NEG_INF.
    Same-type self transitions are disallowed EXCEPT NUC->NUC (di-nucleosome) and TF->TF
    (adjacent TFs); OPEN->OPEN and UNID->UNID are disallowed (a segment is already maximal).
    """
    t = np.full((N_STATES, N_STATES), NEG_INF)

    # returning to OPEN from a real footprint is always free
    t[NUC, OPEN] = 0.0
    t[TF, OPEN] = 0.0
    t[UNID, OPEN] = 0.0
    # OPEN -> OPEN disallowed (stays NEG_INF)

    # entering NUC
    t[OPEN, NUC] = 0.0
    t[NUC, NUC] = -params.trans_nuc_nuc
    t[TF, NUC] = -params.trans_nuc_tf
    t[UNID, NUC] = -params.trans_unid_adj

    # entering TF
    t[OPEN, TF] = 0.0
    t[TF, TF] = -params.trans_tf_tf
    t[NUC, TF] = -params.trans_nuc_tf
    t[UNID, TF] = -params.trans_unid_adj

    # entering UNID
    t[OPEN, UNID] = 0.0
    t[NUC, UNID] = -params.trans_unid_adj
    t[TF, UNID] = -params.trans_unid_adj
    # UNID -> UNID disallowed

    return t


def start_costs(params):
    """(S,) additive log-score for beginning a segment of each type (parsimony)."""
    return np.array([
        -params.start_cost_open,
        -params.start_cost_nuc,
        -params.start_cost_tf,
        -params.start_cost_unid,
    ])


def nuc_duration_logprob(dur, params):
    """Gaussian (log-unnormalized) duration prior for a nucleosome, -inf outside [min,max]."""
    if dur < params.nuc_min or dur > params.nuc_max:
        return NEG_INF
    z = (dur - params.nuc_mode) / params.nuc_sigma
    return -0.5 * z * z


# ---------------------------------------------------------------------------
# Segmental Viterbi (vectorized across molecules)
# ---------------------------------------------------------------------------

def viterbi_decode(gpc_pos, obs1, obs0, tfbs_positions, params, discovered_footprints=None):
    """
    Decode every molecule at once. Returns:
      paths     : list (len n_mol) of lists of segments (state, start_bp, end_bp, dyad_or_motifk)
      total_ll  : (n_mol,) total Viterbi log-likelihood

    discovered_footprints: optional list of (lo, hi) bp intervals (from bulk discovery or a manual
      annotation). A UNID segment whose CENTER falls in one gets the reduced
      `unid_discovered_start_cost` instead of the default `start_cost_unid` -- the position-specific
      empirical-Bayes prior. None (default) => scalar start cost everywhere (unchanged behavior).
    """
    n_mol = obs1.shape[1]

    # --- symmetric ghost-nucleosome padding on the LOW edge --------------------------------------
    # The auto amp_width only pads the HIGH side (beyond the last GpC), so a nucleosome whose dyad
    # sits at/off the LOW edge -- e.g. the +1 nucleosome downstream of the TSS -- cannot be placed
    # as a full NUC segment (its truncated span falls below nuc_min) and gets misassigned to UNID
    # (which then pollutes the discovered-footprint vocabulary). We decode in an internal frame
    # shifted RIGHT by left_pad so a NUC can anchor its dyad in the [0,left_pad) ghost zone -- there
    # are no GpCs there, so the ghost half contributes no emission, exactly mirroring the existing
    # high-side ghost room. Every GpC / motif / promoter / discovered coord is shifted into this
    # internal frame; on return, segment coords are shifted back and CLAMPED to [0, orig_amp_width],
    # so all OUTPUTS stay in the caller's original bp frame with no negative coordinates.
    left_pad = int(params.nuc_max // 2) + 10
    orig_amp_width = int(params.amp_width)
    gpc_pos = np.asarray(gpc_pos) + left_pad
    tfbs_positions = [tuple([int(t[0]) + left_pad, int(t[1]) + left_pad] + list(t[2:]))
                      for t in tfbs_positions]
    _plo, _phi = params.promoter_lo, params.promoter_hi
    if _plo < _phi:                                   # promoter window is on (0,0 => off)
        _plo, _phi = _plo + left_pad, _phi + left_pad
    params = dc_replace(params, amp_width=orig_amp_width + left_pad,
                        promoter_lo=_plo, promoter_hi=_phi)

    # each discovered footprint is (lo, hi) [uses params.unid_discovered_start_cost] or
    # (lo, hi, cost) [per-interval cost, e.g. frequency-graded]. Reduced cost applies to a UNID
    # segment whose CENTER falls in the interval. Shifted into the internal (padded) frame.
    disc = []
    for t in (discovered_footprints or []):
        disc.append((float(t[0]) + left_pad, float(t[1]) + left_pad,
                     float(t[2]) if len(t) == 3 else params.unid_discovered_start_cost))
    bnd = build_boundary_grid(params.amp_width, params.grid_step, tfbs_positions, params.tf_margin)
    G = len(bnd) - 1                                    # number of grid intervals
    bidx = np.searchsorted(gpc_pos, bnd, side='left')   # GpC index at each boundary; (G+1,)

    flat_cum = build_flat_cumulatives(gpc_pos, obs1, obs0, params)
    trans = build_transition_matrix(params)
    scost = start_costs(params)

    # leading-zero cumsum of observed-protected calls over the GpC axis, so the count of protected
    # GpCs in a boundary window [ai, j) is cum_obs1[bidx[j]] - cum_obs1[bidx[ai]] (per molecule).
    # Used to enforce params.nuc_min_prot_gpcs (NUC evidence floor).
    cum_obs1 = np.zeros((obs1.shape[0] + 1, n_mol))
    np.cumsum(obs1, axis=0, out=cum_obs1[1:])

    tf_segs = motif_segments_on_grid(bnd, tfbs_positions, params.tf_margin)
    tf_segs_by_end = {}
    for (ai, bi, k) in tf_segs:
        tf_segs_by_end.setdefault(bi, []).append((ai, k))
    mot_iv = unid_blocked_intervals(tfbs_positions, params.tf_margin, params)

    # V[j, s, m]  = best log-score of a segmentation of [0, bnd[j]) whose last segment is type s
    # bpA[j, s, m] = start boundary index of that last segment
    # bpP[j, s, m] = state of the segment preceding it (or -1 for START)
    V = np.full((G + 1, N_STATES, n_mol), NEG_INF)
    bpA = np.full((G + 1, N_STATES, n_mol), -1, dtype=np.int32)
    bpP = np.full((G + 1, N_STATES, n_mol), -1, dtype=np.int8)

    # incoming[a, s_new, m] = max over s_prev (V[a, s_prev, m] + trans[s_prev, s_new]); + argmax
    inc = np.full((G + 1, N_STATES, n_mol), NEG_INF)
    inc_arg = np.full((G + 1, N_STATES, n_mol), -1, dtype=np.int8)
    inc[0, :, :] = 0.0            # START: any state may begin the molecule, no transition cost
    inc_arg[0, :, :] = -1

    def flat_emission(state, ai, bi):
        C = flat_cum[state]
        return C[bidx[bi]] - C[bidx[ai]]

    for j in range(1, G + 1):
        end_bp = bnd[j]

        for s_new in (OPEN, NUC, TF, UNID):
            best = np.full(n_mol, NEG_INF)
            best_a = np.full(n_mol, -1, dtype=np.int32)
            best_p = np.full(n_mol, -1, dtype=np.int8)

            if s_new == TF:
                candidates = tf_segs_by_end.get(j, [])       # [(ai, k), ...]
                for (ai, _k) in candidates:
                    dur = end_bp - bnd[ai]
                    em = flat_emission(TF, ai, j)
                    cand = inc[ai, TF, :] + scost[TF] + em
                    upd = cand > best
                    best = np.where(upd, cand, best)
                    best_a = np.where(upd, ai, best_a)
                    best_p = np.where(upd, inc_arg[ai, TF, :], best_p)

            elif s_new == OPEN:
                for ai in range(j):
                    dur = end_bp - bnd[ai]
                    if dur < 1:
                        continue
                    em = flat_emission(OPEN, ai, j)
                    cand = inc[ai, OPEN, :] + scost[OPEN] + em
                    upd = cand > best
                    best = np.where(upd, cand, best)
                    best_a = np.where(upd, ai, best_a)
                    best_p = np.where(upd, inc_arg[ai, OPEN, :], best_p)

            elif s_new == NUC:
                for ai in range(j):
                    dur = end_bp - bnd[ai]
                    dlp = nuc_duration_logprob(dur, params)
                    if dlp <= NEG_INF / 2:
                        continue
                    dyad = 0.5 * (bnd[ai] + end_bp)
                    half_len = 0.5 * (end_bp - bnd[ai])
                    em = nuc_emission(bidx[ai], bidx[j], dyad, half_len,
                                      gpc_pos, obs1, obs0, params)
                    cand = inc[ai, NUC, :] + scost[NUC] + dlp + em
                    if params.nuc_min_prot_gpcs > 0:                # evidence floor (per molecule)
                        n_prot = cum_obs1[bidx[j], :] - cum_obs1[bidx[ai], :]
                        cand = np.where(n_prot >= params.nuc_min_prot_gpcs, cand, NEG_INF)
                    upd = cand > best
                    best = np.where(upd, cand, best)
                    best_a = np.where(upd, ai, best_a)
                    best_p = np.where(upd, inc_arg[ai, NUC, :], best_p)

            elif s_new == UNID:
                for ai in range(j):
                    dur = end_bp - bnd[ai]
                    if dur < params.unid_min or dur > params.unid_max:
                        continue
                    lo, hi = bnd[ai], end_bp
                    if any(not (hi <= mlo or lo >= mhi) for (mlo, mhi) in mot_iv):
                        continue          # UNID may not overlap a BLOCKING motif (see mot_iv above)
                    em = flat_emission(UNID, ai, j)
                    sc_unid = scost[UNID]
                    if disc:
                        c = 0.5 * (lo + hi)
                        for (dlo, dhi, dcost) in disc:                     # position-specific prior
                            if dlo <= c <= dhi:
                                sc_unid = -dcost
                                break
                    cand = inc[ai, UNID, :] + sc_unid + em
                    upd = cand > best
                    best = np.where(upd, cand, best)
                    best_a = np.where(upd, ai, best_a)
                    best_p = np.where(upd, inc_arg[ai, UNID, :], best_p)

            V[j, s_new, :] = best
            bpA[j, s_new, :] = best_a
            bpP[j, s_new, :] = best_p

        # finalize incoming[j] from V[j]
        # cand_prev[s_prev, s_new, m] = V[j, s_prev, m] + trans[s_prev, s_new]
        cand_prev = V[j][:, None, :] + trans[:, :, None]      # (S_prev, S_new, n_mol)
        inc[j] = cand_prev.max(axis=0)                        # (S_new, n_mol)
        inc_arg[j] = cand_prev.argmax(axis=0).astype(np.int8)

    # --- backtrace per molecule ---
    total_ll = V[G].max(axis=0)
    end_state = V[G].argmax(axis=0)

    paths = []
    for m in range(n_mol):
        segs = []
        j = G
        s = int(end_state[m])
        while j > 0:
            a = int(bpA[j, s, m])
            p = int(bpP[j, s, m])
            if a < 0:
                break
            # shift back to the caller's original frame (undo left_pad)
            start_bp, end_bp = int(bnd[a]) - left_pad, int(bnd[j]) - left_pad
            if s == NUC:
                # dyad = true modeled center of the footprint. For an off-LOW-edge +1 nucleosome this
                # is < 0; KEEP it negative (2026-07-20) -- clamping all +1s to 0 discards the real
                # ordering their observed (interior) protection boundary carries. Safe: the dyad is a
                # reported VALUE, never an array/column index (only start/end/columns must stay >=0,
                # and they are, below), so no negative-index wraparound. Treat an off-edge dyad as a
                # prior-regularized estimate (pinned by the interior boundary + duration prior, NOT a
                # measured position) -- do not use edge dyads for absolute nucleosome phasing.
                # High side is capped as a defensive no-op: orig_amp_width already contains the high
                # ghost pad, so an interior/high nuc center never exceeds it.
                anchor = 0.5 * (start_bp + end_bp)
                anchor = float(min(anchor, orig_amp_width))
            elif s == TF:
                anchor = _motif_for_segment(a, j, tf_segs)    # motif index (0-based) or -1
            else:
                anchor = -1
            start_bp = max(0, start_bp)                       # clamp span to [0, orig_amp_width]
            end_bp = min(orig_amp_width, end_bp)
            if end_bp > start_bp:                             # drop any all-ghost (zero-length) seg
                segs.append((s, start_bp, end_bp, anchor))
            j = a
            s = p
            if s < 0:
                break
        segs.reverse()
        paths.append(segs)

    return paths, total_ll


def _motif_for_segment(ai, bi, tf_segs):
    for (a, b, k) in tf_segs:
        if a == ai and b == bi:
            return k
    return -1


# ---------------------------------------------------------------------------
# Scoring an arbitrary segmentation (diagnostics / calibration)
# ---------------------------------------------------------------------------

def _seg_emission(state, a, b, dyad, gpc_pos, obs1_col, obs0_col, params):
    """Emission of one segment [a,b) for a single molecule column (obs*_col are (n_gpc,))."""
    jlo = int(np.searchsorted(gpc_pos, a, side='left'))
    jhi = int(np.searchsorted(gpc_pos, b, side='left'))
    if jhi <= jlo:
        return 0.0
    if state == NUC:
        half_len = 0.5 * (b - a)
        return float(nuc_emission(jlo, jhi, dyad, half_len, gpc_pos,
                                  obs1_col[:, None], obs0_col[:, None], params)[0])
    lw1, lw0 = _logw(flat_state_pu(state, gpc_pos, params), params.p_t_given_unmeth,
                     params.p_t_given_meth)
    return float(lw1[jlo:jhi] @ obs1_col[jlo:jhi] + lw0[jlo:jhi] @ obs0_col[jlo:jhi])


def per_gpc_protection(mat, gpc_pos, paths=None, exclude_states=(NUC,)):
    """
    Per-GpC bulk protected fraction (ignoring missing). If `paths` (a prior decode) is given,
    DROP every (molecule, position) observation the decode assigned to a segment in `exclude_states`
    -- i.e. nucleosome-decontaminate the bulk, so recurrent NUCLEOSOME coverage doesn't masquerade as
    a footprint prior. Returns (frac[n_gpc], n_used[n_gpc]). This is the same conditioning idea used
    for the TetO pseudobulk (condition out NUC-explained signal). Excluding TF too (exclude_states=
    (NUC,TF)) additionally removes annotated-motif protection.
    """
    gp = np.asarray(gpc_pos)
    V = mat.values
    keep = (V >= 0)
    if paths is not None:
        for m, segs in enumerate(paths):
            for (s, a, b, _anc) in segs:
                if s in exclude_states:
                    keep[m, (gp >= a) & (gp < b)] = False
    frac = np.full(V.shape[1], np.nan)
    n_used = keep.sum(axis=0)
    for j in range(V.shape[1]):
        if n_used[j] > 0:
            frac[j] = V[keep[:, j], j].mean()
    return frac, n_used


def discover_footprints_from_bulk(mat, gpc_pos, tfbs_positions, params, frac=None,
                                  abs_floor=0.25, rel_delta=0.20, bg_window=80, bg_pct=25,
                                  min_gpcs=1, merge_gap=25):
    """
    Empirical-Bayes footprint discovery. Find OFF-MOTIF sites with recurrent bulk protection and
    return (lo, hi) bp intervals for viterbi_decode(discovered_footprints). Cheap approximation to a
    joint cross-molecule fit: a site's recurrent bulk protection becomes a per-position prior that
    lowers the UNID start cost there, so a footprint the per-molecule decode would miss gets called
    because the POPULATION says one lives there.

    LOCAL-BACKGROUND-RELATIVE criterion (an absolute 0.5 is too strict and ignores local baseline):
    a GpC is flagged if its (decontaminated) protected fraction >= abs_floor AND exceeds the local
    background by >= rel_delta, where local background = the bg_pct-th percentile of frac within
    +-bg_window bp. So a peak on a low local baseline (e.g. BAX 171/188 ~0.40 over a ~0.10 floor)
    is caught even though it's below 0.5. Off-motif only. Adjacent flags (<= merge_gap) merge;
    intervals with >= min_gpcs flagged GpCs and span <= unid_max are kept, padded to >= unid_min.

    Each returned interval is (lo, hi, score) where score = SUM over the interval's flagged GpCs of
    (frac - local_background) -- i.e. it folds in BOTH per-GpC strength (bigger excess) AND number of
    supporting GpCs (more terms). A lone but very-strong GpC and a multi-GpC moderate site can both
    score high; the caller maps score -> UNID start cost smoothly. Returns (intervals, per_gpc_frac).
    """
    gpc = np.asarray(gpc_pos)
    if frac is None:                     # raw bulk (nucleosome-CONTAMINATED); prefer a decontaminated
        frac, _ = per_gpc_protection(mat, gpc_pos, paths=None)   # frac (per_gpc_protection(paths=...))
    fr = np.asarray(frac, dtype=float)
    # DELIBERATELY uses motif_intervals(), not unid_blocked_intervals(): vocabulary discovery is
    # stage 1 and is FROZEN from attempt2 (README_attempt3.md sec.3). Making it follow
    # unid_may_overlap_motifs would make the vocabulary itself a second decode-affecting variable
    # per round. Keep discovery blind to named sites regardless of the DP's exclusivity setting.
    mot_iv = motif_intervals(tfbs_positions, params.tf_margin)
    off_motif = lambda p: not any(mlo <= p <= mhi for (mlo, mhi) in mot_iv)
    flagged = []                                          # (pos, excess_over_background)
    for j in range(len(gpc)):
        if not np.isfinite(fr[j]) or fr[j] < abs_floor or not off_motif(int(gpc[j])):
            continue
        win = (gpc >= gpc[j] - bg_window) & (gpc <= gpc[j] + bg_window) & np.isfinite(fr)
        bg = np.percentile(fr[win], bg_pct) if win.sum() >= 3 else 0.0
        if fr[j] - bg >= rel_delta:
            flagged.append((int(gpc[j]), float(fr[j] - bg)))
    groups, cur = [], []
    for (p, ex) in flagged:
        if cur and p - cur[-1][0] > merge_gap:
            groups.append(cur); cur = []
        cur.append((p, ex))
    if cur:
        groups.append(cur)
    out = []
    for grp in groups:
        ps = [p for (p, ex) in grp]
        lo, hi, score = min(ps), max(ps), float(sum(ex for (p, ex) in grp))
        if len(grp) < min_gpcs or (hi - lo) > params.unid_max:
            continue
        if hi - lo < params.unid_min:        # pad narrow/single-GpC peaks so a UNID segment can match
            c = 0.5 * (lo + hi)
            lo, hi = int(c - params.unid_min / 2.0), int(c + params.unid_min / 2.0)
        out.append((lo, hi, score))
    return out, fr


def load_discovered_footprints_file(path, amplicon_name):
    """
    Load a precomputed footprint VOCABULARY (from discover_footprint_vocabulary.py). Same block
    layout as positions.long.txt but data lines are `lo,hi,cost` (cost = reduced UNID start cost).
    Returns [(lo, hi, cost), ...] for `amplicon_name` (empty if absent).
    """
    out, cur = [], None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                cur = s[1:]
                continue
            if cur == amplicon_name:
                parts = s.split(',')
                out.append((int(parts[0]), int(parts[1]), float(parts[2])))
    return out


def score_segmentation(segs, gpc_pos, obs1_col, obs0_col, tfbs_positions, params):
    """
    Score an arbitrary segmentation of ONE molecule exactly as the Viterbi DP would, returning
    (total, breakdown_rows). Each row: type, start, end, emission, start_cost, dur_logprob,
    trans_from_prev, subtotal. Used to compare the model's chosen path against a hypothesized
    'correct' one and see which term (emission vs structural cost) drives the difference.
    Reports NEG_INF contributions for structurally-disallowed transitions/placements.
    """
    trans = build_transition_matrix(params)
    scost = start_costs(params)
    mot_iv = motif_intervals(tfbs_positions, params.tf_margin)              # TF legality
    unid_iv = unid_blocked_intervals(tfbs_positions, params.tf_margin, params)  # UNID legality
    gpc_pos = np.asarray(gpc_pos, dtype=np.int64)

    rows = []
    total = 0.0
    prev = None
    for (s, a, b, anchor) in segs:
        dyad = 0.5 * (a + b) if s == NUC else anchor
        em = _seg_emission(s, a, b, dyad, gpc_pos, obs1_col, obs0_col, params)
        sc = float(scost[s])
        dlp = nuc_duration_logprob(b - a, params) if s == NUC else 0.0
        tr = 0.0 if prev is None else float(trans[prev, s])
        # structural validity flags (don't silently pass an illegal path)
        note = ''
        if s == TF and not any(a >= mlo and b <= mhi for (mlo, mhi) in mot_iv):
            note = 'TF-not-at-motif(illegal)'
        if s == UNID and any(not (b <= mlo or a >= mhi) for (mlo, mhi) in unid_iv):
            note = 'UNID-overlaps-motif(illegal)'
        if s == UNID and not (params.unid_min <= (b - a) <= params.unid_max):
            note = (note + ' ' if note else '') + 'UNID-len-out-of-band'
        sub = em + sc + dlp + tr
        total += sub
        rows.append({'type': STATE_NAMES[s], 'start': a, 'end': b,
                     'emission': em, 'start_cost': sc, 'dur_logprob': dlp,
                     'trans_from_prev': tr, 'subtotal': sub, 'note': note})
        prev = s
    return total, rows


def fill_open_gaps(footprints, amp_width, gpc_pos, tfbs_positions, params):
    """
    Given only the footprint segments (list of (state, a, b)), return a full tiling with OPEN
    segments filling every gap (incl. before the first and after the last), and anchors filled
    (NUC dyad = center; TF motif index = the motif it spans, else -1).
    """
    fps = sorted(footprints, key=lambda x: x[1])
    full = []
    cur = 0
    for (s, a, b) in fps:
        if a > cur:
            full.append((OPEN, cur, a, -1))
        if s == NUC:
            anchor = 0.5 * (a + b)
        elif s == TF:
            anchor = _motif_for_segment_bp(a, b, tfbs_positions, params.tf_margin)
        else:
            anchor = -1
        full.append((s, a, b, anchor))
        cur = b
    if cur < amp_width:
        full.append((OPEN, cur, amp_width, -1))
    return full


def _motif_for_segment_bp(a, b, tfbs_positions, tf_margin):
    for k, t in enumerate(tfbs_positions):
        if a == int(t[0]) - tf_margin and b == int(t[1]) + tf_margin:
            return k
    # fall back: nearest motif center
    for k, t in enumerate(tfbs_positions):
        if a <= 0.5 * (int(t[0]) + int(t[1])) <= b:
            return k
    return -1


# ---------------------------------------------------------------------------
# EM (Viterbi / hard EM) parameter fitting
# ---------------------------------------------------------------------------

def invert_t_fraction(t_frac, p_t_given_unmeth, p_t_given_meth):
    """Recover p_U from observed protected-fraction by inverting the conversion model. v4's math."""
    denom = p_t_given_unmeth - p_t_given_meth
    if abs(denom) < 1e-10:
        return None
    return float(np.clip((t_frac - p_t_given_meth) / denom, 0.0, 1.0))


def em_mstep(paths, gpc_pos, obs, params, min_obs=100):
    """
    Hard-EM M-step. From the current Viterbi segmentation of every molecule, re-estimate emission
    params (prob_unmeth_given_open/tf/unid) and the nucleosome duration prior (nuc_mode, nuc_sigma).
    Conversion rates, nuc_softness (fixed breathing constant), and structural (start/transition)
    costs are held fixed -- and there is no d_edge to fit anymore: the footprint width is tied to
    the segment length, so extent is fit only via the duration prior. Returns (estimates_dict with
    None where n_obs<min_obs, n_obs_dict). obs: (n_gpc, n_mol) with values in {1, 0, -1}.
    """
    gp = np.asarray(gpc_pos, dtype=np.int64)
    open_o, tf_o, unid_o = [], [], []
    nuc_lens = []
    for m, segs in enumerate(paths):
        col = obs[:, m]
        for (s, a, b, anchor) in segs:
            jlo = int(np.searchsorted(gp, a, side='left'))
            jhi = int(np.searchsorted(gp, b, side='left'))
            if jhi > jlo:
                vv = col[jlo:jhi]
                vals = vv[vv >= 0].astype(float)
                if s == OPEN:
                    open_o.append(vals)
                elif s == TF:
                    tf_o.append(vals)
                elif s == UNID:
                    unid_o.append(vals)
            if s == NUC:
                nuc_lens.append(b - a)

    def cat(lst):
        return np.concatenate(lst) if lst else np.array([])

    open_o, tf_o, unid_o = cat(open_o), cat(tf_o), cat(unid_o)
    p_tu, p_tm = params.p_t_given_unmeth, params.p_t_given_meth

    est = {'prob_unmeth_given_open': None, 'prob_unmeth_given_tf': None,
           'prob_unmeth_given_unid': None, 'nuc_mode': None, 'nuc_sigma': None}
    if len(open_o) >= min_obs:
        est['prob_unmeth_given_open'] = invert_t_fraction(open_o.mean(), p_tu, p_tm)
    if len(tf_o) >= min_obs:
        est['prob_unmeth_given_tf'] = invert_t_fraction(tf_o.mean(), p_tu, p_tm)
    if len(unid_o) >= min_obs:
        est['prob_unmeth_given_unid'] = invert_t_fraction(unid_o.mean(), p_tu, p_tm)
    if len(nuc_lens) >= max(10, min_obs // 10):
        est['nuc_mode'] = float(np.mean(nuc_lens))
        est['nuc_sigma'] = float(max(np.std(nuc_lens), 5.0))

    n_obs = {'open': len(open_o), 'tf': len(tf_o), 'unid': len(unid_o),
             'nuc_segs': len(nuc_lens)}
    return est, n_obs


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def tfbs_site_names(tfbs_positions):
    """
    Stable, unique, machine-filterable name for each positions-file site.

    The positions file is `lo,hi,name[,strand]`. The name column is parsed by
    load_tfbs_positions but was historically never used for output, so molecules were
    reported only as `tfbs_{k}` where k is ORDER OF APPEARANCE in the positions file.
    That index is not stable: inserting a promoter site (which sorts BEFORE the array,
    since the promoter is at LOW matrix columns and TetO at HIGH) renumbers every TetO
    and silently invalidates any analysis keyed on `tfbs_1..6 == TetO1..6`.

    Names fix that. Every opoBD9 entry is literally `TetO`, so names are not unique on
    their own -- duplicates get a 1-based ordinal suffix IN FILE ORDER
    (TetO x6 -> TetO1..TetO6). Sites with no name column fall back to `site{k+1}`.

    Downstream code should key on the NAME, never the index.
    """
    raw = []
    for k, t in enumerate(tfbs_positions):
        nm = str(t[2]).strip() if (len(t) > 2 and str(t[2]).strip()) else ''
        nm = re.sub(r'[^0-9A-Za-z]+', '_', nm).strip('_')
        raw.append(nm or 'site{}'.format(k + 1))

    counts = Counter(raw)
    seen = defaultdict(int)
    out = []
    for nm in raw:
        if counts[nm] > 1:
            seen[nm] += 1
            out.append('{}{}'.format(nm, seen[nm]))
        else:
            out.append(nm)

    # Guard the pathological case where suffixing collides with a literal name
    # (e.g. a file holding both `TetO` x6 and a site actually called `TetO1`).
    if len(set(out)) != len(out):
        used = set()
        for i, nm in enumerate(out):
            cand, j = nm, 1
            while cand in used:
                cand, j = '{}_{}'.format(nm, j), j + 1
            used.add(cand)
            out[i] = cand
    return out


def is_teto_site(name, teto_prefix='TetO'):
    """A site is 'privileged' (part of the synthetic array) iff its NAME starts with the
    prefix, case-insensitively. Everything else -- TATA, BRE, Inr, ATF1, pause -- is a
    discovered/annotated promoter footprint and is counted separately."""
    return bool(teto_prefix) and name.lower().startswith(teto_prefix.lower())


def paths_to_wide(paths, total_ll, read_ids, n_tfbs, site_names=None, teto_prefix='TetO'):
    """
    One row per molecule.

    Columns, in order:
      tfbs_{k}        legacy positional booleans, 1-based, UNCHANGED. Kept because
                      common.py and the partition-function/lattice-gas fitters read them
                      positionally or via `df.filter(like='tfbs_')`.
      site_{name}     the same booleans keyed by stable site name (see tfbs_site_names).
                      Deliberately a DIFFERENT prefix from `tfbs_` so it does not get
                      swept up by existing `filter(like='tfbs_')` calls, which would
                      otherwise double-count every site.
      n_tf            all TF segments (legacy meaning: grows as the site inventory grows)
      n_tf_teto       TF segments at array sites only  <-- use this for "num TetO bound"
      n_tf_other      TF segments at named promoter sites
      n_nuc / n_unid / nucs / unids / log_likelihood
    """
    if site_names is None:
        site_names = ['site{}'.format(k + 1) for k in range(n_tfbs)]
    teto_idx = {k for k, nm in enumerate(site_names) if is_teto_site(nm, teto_prefix)}

    rows = []
    for m, segs in enumerate(paths):
        tf_bound = set()
        nucs = []
        unids = []
        for (s, a, b, anchor) in segs:
            if s == TF and anchor >= 0:
                tf_bound.add(int(anchor))
            elif s == NUC:
                nucs.append((a, b, anchor))
            elif s == UNID:
                unids.append((a, b))
        row = {'read_id': read_ids[m]}
        for k in range(n_tfbs):
            row['tfbs_{}'.format(k + 1)] = (k in tf_bound)
        for k in range(n_tfbs):
            row['site_{}'.format(site_names[k])] = (k in tf_bound)
        row['n_tf'] = len(tf_bound)
        row['n_tf_teto'] = len(tf_bound & teto_idx)
        row['n_tf_other'] = len(tf_bound - teto_idx)
        row['n_nuc'] = len(nucs)
        row['n_unid'] = len(unids)
        row['nucs'] = ';'.join('{}:{}:{:.0f}'.format(a, b, d) for (a, b, d) in nucs)
        row['unids'] = ';'.join('{}:{}'.format(a, b) for (a, b) in unids)
        row['log_likelihood'] = float(total_ll[m])
        rows.append(row)
    return pd.DataFrame(rows).set_index('read_id')


def paths_to_tidy(paths, read_ids, site_names=None):
    """Tidy long-format Viterbi path: one row per segment.

    `motif_name` resolves the TF anchor index to its stable site name (blank for
    non-TF segments and for anchor-less TF), so the segments table can be filtered
    on TetO-vs-promoter without joining back to the positions file.
    """
    rows = []
    for m, segs in enumerate(paths):
        for i, (s, a, b, anchor) in enumerate(segs):
            is_named_tf = (s == TF and anchor >= 0)
            name = ''
            if is_named_tf and site_names is not None and int(anchor) < len(site_names):
                name = site_names[int(anchor)]
            rows.append({
                'read_id': read_ids[m],
                'seg_index': i,
                'type': STATE_NAMES[s],
                'start': a,
                'end': b,
                # NUC always writes its dyad (may be <0 for an off-low-edge +1 nuc; see viterbi_decode
                # 2026-07-20). TF writes its motif index only when it maps to one (anchor -1 = no motif).
                'dyad_or_motif': (int(anchor) if (s == NUC or is_named_tf) else ''),
                'motif_name': name,
            })
    return pd.DataFrame(rows)


def write_tfbs_index(path, amplicon_name, tfbs_positions, site_names, teto_prefix='TetO'):
    """Sidecar mapping index -> (name, lo, hi, is_teto) so an existing `tfbs_{k}` file
    stays interpretable after the positions file changes."""
    with open(path, 'w') as fh:
        fh.write('# tfbs index map for {}\n'.format(amplicon_name))
        fh.write('# `index` is 1-based and matches the tfbs_{k} columns; it is NOT stable\n'
                 '# across positions-file edits. `name` is. Key downstream code on `name`.\n')
        fh.write('index\tname\tlo\thi\tis_teto\n')
        for k, (t, nm) in enumerate(zip(tfbs_positions, site_names)):
            fh.write('{}\t{}\t{}\t{}\t{}\n'.format(
                k + 1, nm, int(t[0]), int(t[1]), is_teto_site(nm, teto_prefix)))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def load_motif_track(path, amplicon_name):
    """
    Load a per-amplicon motif annotation track (from export_motif_track.py) for overlaying known
    motifs on the diagnostic/per-read plots. Same block layout as positions.long.txt but data lines
    are `lo,hi,label[,color]` in the SAME matrix-column/bp frame the model uses (colors optional;
    default gray). Returns [(lo, hi, label, color), ...] for `amplicon_name` (empty if absent/None).
    Purely cosmetic -- never touches decoding.
    """
    if not path:
        return []
    out, cur = [], None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            if s.startswith('>'):
                cur = s[1:].strip()
                continue
            if cur == amplicon_name:
                parts = [x.strip() for x in s.split(',')]
                if len(parts) >= 3:
                    lo, hi, label = int(parts[0]), int(parts[1]), parts[2]
                    color = parts[3] if len(parts) >= 4 and parts[3] else '0.4'
                    out.append((lo, hi, label, color))
    return out


def draw_motif_strip(ax, motifs, y0=1.01, h=0.045, fontsize=6, label=True):
    """
    Draw a motif annotation strip just ABOVE the axis (blended transform: x=data bp, y=axes frac),
    so it never overlaps the data regardless of ylim. Each motif = a colored bar + rotated label.
    """
    if not motifs:
        return
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for (lo, hi, name, color) in motifs:
        ax.add_patch(Rectangle((lo, y0), max(hi - lo, 1), h, transform=trans,
                               color=color, lw=0, clip_on=False, zorder=5))
        ax.axvspan(lo, hi, color=color, alpha=0.07, lw=0, zorder=0)   # faint in-plot guide band
        if label:
            ax.text(0.5 * (lo + hi), y0 + h + 0.01, name, transform=trans, rotation=90,
                    fontsize=fontsize, ha='center', va='bottom', color=color, clip_on=False)


def decorate_hsmm_read(ax, segs, tfbs_positions):
    """Draw the Viterbi segmentation for one molecule: NUC gray, TF red, UNID purple."""
    for (s, a, b, anchor) in segs:
        if s == NUC:
            ax.add_patch(Ellipse((0.5 * (a + b), 0.1), (b - a), 0.1, color='0.4'))
        elif s == TF and anchor >= 0:
            ts, te = int(tfbs_positions[int(anchor)][0]), int(tfbs_positions[int(anchor)][1])
            ax.add_patch(Ellipse((0.5 * (ts + te) + 3, 0.1), 10, 0.1, color='r'))
        elif s == UNID:
            ax.add_patch(Ellipse((0.5 * (a + b), 0.1), (b - a), 0.1, color='blueviolet'))


# state colormap shared by the bulk plot (index = state code: OPEN, NUC, TF, UNID)
STATE_CMAP = ListedColormap(['white', '0.5', 'red', 'blueviolet'])
STATE_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], STATE_CMAP.N)


def build_state_matrix_bp(paths_subset, amp_width):
    """(n_reads, amp_width) int8 of the Viterbi state at each bp (OPEN=0 baseline)."""
    mat = np.zeros((len(paths_subset), amp_width), dtype=np.int8)
    for i, segs in enumerate(paths_subset):
        for (s, a, b, _anchor) in segs:
            a2, b2 = max(0, int(a)), min(amp_width, int(b))
            if b2 > a2:
                mat[i, a2:b2] = s
    return mat


def build_data_matrix_bp(vals, gpc_pos, amp_width, halfwidth=2):
    """
    (n_reads, amp_width) grayscale display of the raw data on the bp axis: each GpC widened to
    +/- halfwidth bp so it is visible. protected -> 0.0 (black), accessible -> 0.75 (gray),
    missing/no-data -> NaN (white). vals is (n_reads, n_gpc).
    """
    D = np.full((vals.shape[0], amp_width), np.nan)
    for j, p in enumerate(gpc_pos):
        lo, hi = max(0, int(p) - halfwidth), min(amp_width, int(p) + halfwidth + 1)
        if hi <= lo:
            continue
        col = vals[:, j].astype(float)
        disp = np.where(col == 1, 0.0, np.where(col == 0, 0.75, np.nan))
        D[:, lo:hi] = disp[:, None]
    return D


def plot_bulk_sanity(mat, gpc_pos, paths, params, plots, tfbs_positions,
                     n_reads=1000, seed=42):
    """
    Two stacked panels sharing the bp axis: observed data (top) and the HSMM state prediction
    (bottom, 4 colors), for the same reads, sorted by predicted-state pattern so coherent calls
    form bands. Lets you eyeball whether protected stretches map to sensible states in bulk.
    """
    n = len(mat)
    if n == 0:
        return
    k = min(n_reads, n)
    sel = np.random.RandomState(seed).choice(n, k, replace=False)

    paths_sub = [paths[i] for i in sel]
    state_bp = build_state_matrix_bp(paths_sub, params.amp_width)
    order = np.lexsort(state_bp[:, ::-1].T)          # group similar predictions (primary = col 0)
    state_bp = state_bp[order]
    vals_sub = mat.values[sel][order]
    data_bp = build_data_matrix_bp(vals_sub, gpc_pos, params.amp_width)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12), sharex=True)

    ax1.imshow(data_bp, aspect='auto', cmap='gray', vmin=0, vmax=1, interpolation='none')
    ax1.set_title('observed data (n={} reads, sorted by prediction): '
                  'black=protected, gray=accessible, white=no data'.format(k))
    ax1.set_ylabel('molecules')

    ax2.imshow(state_bp, aspect='auto', cmap=STATE_CMAP, norm=STATE_NORM, interpolation='none')
    ax2.set_title('HSMM prediction')
    ax2.set_xlabel('position (bp)')
    ax2.set_ylabel('molecules')
    ax2.legend(handles=[Patch(facecolor='white', edgecolor='k', label='OPEN'),
                        Patch(facecolor='0.5', label='NUC'),
                        Patch(facecolor='red', label='TF'),
                        Patch(facecolor='blueviolet', label='UNID')],
               loc='upper right', ncol=4, fontsize=8, framealpha=0.9)

    # mark motif positions on both panels for orientation
    for t in tfbs_positions:
        for ax in (ax1, ax2):
            ax.axvline(int(t[0]), color='orange', lw=0.4, alpha=0.6)
            ax.axvline(int(t[1]), color='orange', lw=0.4, alpha=0.6)

    fig.tight_layout()
    plots.savefig(fig)
    plt.close(fig)


def plot_footprint_diagnostic(mat, gpc_pos, paths, params, plots, tfbs_positions, diag,
                              motifs=None):
    """
    Diagnostic pseudobulk page (always emitted): TOP = per-bp occupancy of NUC/TF/UNID across all
    molecules (where each footprint type is called); BOTTOM = raw bulk protection and, if a
    discovery pass ran, the nucleosome(+TF)-decontaminated bulk, with discovered footprint intervals
    shaded. Lets you eyeball whether calls match the observed signal and judge the discovery cutoff.
    """
    gpc_pos = np.asarray(gpc_pos, dtype=np.int64)
    W = params.amp_width
    sm = build_state_matrix_bp(paths, W)
    x = np.arange(W)
    fig, (axT, axB) = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
    axT.plot(x, (sm == NUC).mean(0), color='0.5', label='NUC')
    axT.plot(x, (sm == TF).mean(0), color='crimson', label='TF')
    axT.plot(x, (sm == UNID).mean(0), color='blueviolet', label='UNID')
    for t in tfbs_positions:
        axT.axvspan(int(t[0]), int(t[1]), color='red', alpha=0.06)
    for (lo, hi, c) in diag.get('discovered', []):
        axT.axvspan(lo, hi, color='blueviolet', alpha=0.12)
    draw_motif_strip(axT, motifs)
    axT.set_ylim(0, 1); axT.set_ylabel('occupancy'); axT.legend(fontsize=8, loc='upper left')
    axT.set_title('per-bp state occupancy (n={})   discovered footprints: {}'.format(
        len(mat), [(l, h) for l, h, _ in diag.get('discovered', [])]), fontsize=9,
        pad=(34 if motifs else None))
    raw = diag.get('raw_frac'); decon = diag.get('decon_frac')
    if raw is not None:
        axB.plot(gpc_pos, raw, '-o', ms=2.5, color='black', label='raw bulk protection')
    if decon is not None:
        axB.plot(gpc_pos, decon, '-o', ms=2.5, color='darkorange',
                 label='nucleosome(+TF)-decontaminated bulk')
    for thr, cc in [(0.3, '0.75'), (0.4, '0.6'), (0.5, '0.45')]:
        axB.axhline(thr, color=cc, lw=0.8, ls='--')
    for (lo, hi, c) in diag.get('discovered', []):
        axB.axvspan(lo, hi, color='blueviolet', alpha=0.12)
    axB.set_ylim(0, 1); axB.set_ylabel('protected fraction'); axB.set_xlabel('bp')
    axB.legend(fontsize=8, loc='upper left')
    fig.tight_layout()
    plots.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------------

EM_FIT_ALL = ['prob_unmeth_given_open', 'prob_unmeth_given_tf', 'prob_unmeth_given_unid',
              'nuc_mode', 'nuc_sigma']
# Params SAFE to fit with hard (Viterbi) EM. EXCLUDED and why (verified on the 0x sample
# 2026-07-10, see NOTES_v5_calibration_log.md):
#   prob_unmeth_given_open -- hard assignment hands the nucleosome's soft-edge transition zone
#     entirely to OPEN, so p_open is biased upward by residual protected edge GpCs.
#   nuc_sigma -- decoded segment lengths cluster at the duration-prior mode, so the measured
#     variance shrinks each iteration, tightening the prior further (collapses 25->~5). The real
#     spread (0x model-free protection streaks: median 122, p90 188) is much larger.
# Those need soft/forward-backward EM or external (0x streak) calibration; set them by hand for now.
# nuc_mode is kept: its EM estimate (~133-140) agrees with the 0x model-free footprint (~140).
# (There is no longer a d_edge/softness footprint-shape M-step: the footprint width is tied to the
#  segment length and nuc_softness is a fixed constant -- see NOTES_v5_calibration_log.md 2026-07-13.)
EM_FIT_SAFE = ['prob_unmeth_given_tf', 'prob_unmeth_given_unid', 'nuc_mode']


def run_em(gpc_pos, obs, obs1, obs0, tfbs_positions, params, max_iters=8, min_obs=100,
           damping=0.5, tol=0.02, em_log_path=None, fit_params=None):
    """
    Viterbi (hard) EM. Alternates decode (E) and em_mstep (M), damping each update, until the
    max fractional parameter change < tol or max_iters. Structural costs + conversion rates fixed.
    Only params in `fit_params` (default EM_FIT_SAFE) are applied; the M-step still ESTIMATES the
    excluded ones and logs them (as est_<name>) for diagnostics, but does not apply them.
    Returns (final_params, em_log_rows).
    """
    import copy as _copy
    if fit_params is None:
        fit_params = EM_FIT_SAFE
    cur = _copy.copy(params)
    log_rows = []
    for it in range(max_iters):
        paths, _ = viterbi_decode(gpc_pos, obs1, obs0, tfbs_positions, cur)
        est, n_obs = em_mstep(paths, gpc_pos, obs, cur, min_obs=min_obs)
        fracs = []
        newp = _copy.copy(cur)
        for name in fit_params:
            e = est[name]
            old = getattr(cur, name)
            if e is None:
                continue
            damped = (1 - damping) * old + damping * e
            setattr(newp, name, damped)
            fracs.append(abs(damped - old) / max(abs(old), 1e-9))
        newp.__post_init__()   # re-apply unid_max<=nuc_min clamp
        maxfrac = max(fracs) if fracs else 0.0
        row = {'iter': it, 'max_frac_change': maxfrac,
               **{k: getattr(newp, k) for k in EM_FIT_ALL},                 # applied values
               **{'est_' + k: est[k] for k in EM_FIT_ALL},                 # raw M-step estimates
               **{'n_obs_' + k: v for k, v in n_obs.items()}}
        log_rows.append(row)
        print('  EM iter {}: max_frac_change={:.4f}  nuc_mode={:.0f} nuc_sigma={:.0f} '
              'p_tf={:.2f} p_unid={:.2f}   [est-only p_open={}]'.format(
                  it, maxfrac, newp.nuc_mode, newp.nuc_sigma,
                  newp.prob_unmeth_given_tf, newp.prob_unmeth_given_unid,
                  None if est['prob_unmeth_given_open'] is None else round(est['prob_unmeth_given_open'], 2)))
        cur = newp
        if maxfrac < tol:
            print('  EM converged')
            break
    if em_log_path and log_rows:
        pd.DataFrame(log_rows).to_csv(em_log_path, sep='\t', index=False)
        print('  wrote EM log {}'.format(em_log_path))
    return cur, log_rows


def compute_classifications(mat, gpc_pos, tfbs_positions, params,
                            do_em=False, em_max_iters=8, em_min_obs=100, em_log_path=None,
                            em_fit_params=None, discover_mode='off',
                            fp_abs_floor=0.25, fp_rel_delta=0.20, fp_bg_window=80, fp_score_hi=1.0,
                            precomputed_discovered=None, teto_prefix='TetO',
                            region_entry=None, promoter_thresh=0.5, promoter_min_bp=None):
    """
    mat: DataFrame (n_mol x n_gpc) with values in {1 protected, 0 accessible, -1 missing}.
    discover_mode: 'off'  = single decode, scalar UNID start cost (unchanged behavior);
                   '2pass'= discover footprints from RAW bulk, position-specific UNID prior, decode;
                   '3pass'= decode once -> nucleosome(+TF)-DECONTAMINATED bulk -> discover -> decode.
    Discovery uses a LOCAL-BACKGROUND-RELATIVE criterion (fp_abs_floor / fp_rel_delta / fp_bg_window)
    so a footprint on a low local baseline is caught even below an absolute 0.5. The discovered UNID
    start cost is graded by (decontaminated) strength between start_cost_unid (weak) and
    unid_discovered_start_cost (strong).
    Returns (wide_df, tidy_df, paths, total_ll, fitted_params, diag) where diag has
    raw_frac / decon_frac / discovered for the diagnostic plot.
    """
    gpc_pos = np.asarray(gpc_pos, dtype=np.int64)
    obs = mat.values.T                                  # (n_gpc, n_mol)
    obs1 = (obs == 1).astype(float)
    obs0 = (obs == 0).astype(float)

    if do_em:
        print('running Viterbi-EM (fitting: {})'.format(em_fit_params or EM_FIT_SAFE))
        params, _ = run_em(gpc_pos, obs, obs1, obs0, tfbs_positions, params,
                           max_iters=em_max_iters, min_obs=em_min_obs, em_log_path=em_log_path,
                           fit_params=em_fit_params)
        print('EM done; final params: softness={:.1f} nuc_mode={:.0f} '
              'nuc_sigma={:.0f} p_tf={:.3f} p_open={:.3f} p_unid={:.3f}'.format(
                  params.nuc_softness, params.nuc_mode, params.nuc_sigma,
                  params.prob_unmeth_given_tf, params.prob_unmeth_given_open,
                  params.prob_unmeth_given_unid))

    # always compute raw bulk (for the diagnostic); decon only when we did a baseline pass
    raw_frac, _ = per_gpc_protection(mat, gpc_pos, paths=None)
    decon_frac, discovered = None, []
    if precomputed_discovered is not None:
        # use a shared vocabulary (e.g. pooled per-promoter from discover_footprint_vocabulary.py);
        # no per-amplicon discovery / baseline decode.
        discovered = [(int(l), int(h), float(c)) for (l, h, c) in precomputed_discovered]
        print('using {} precomputed footprint(s): {}'.format(
            len(discovered), [(l, h, round(c, 1)) for l, h, c in discovered]))
    elif discover_mode in ('2pass', '3pass'):
        if discover_mode == '3pass':
            base_paths, _ = viterbi_decode(gpc_pos, obs1, obs0, tfbs_positions, params)
            decon_frac, _ = per_gpc_protection(mat, gpc_pos, paths=base_paths,
                                               exclude_states=(NUC, TF))
            disc_frac = decon_frac
        else:
            disc_frac = raw_frac
        intervals, _ = discover_footprints_from_bulk(
            mat, gpc_pos, tfbs_positions, params, frac=disc_frac,
            abs_floor=fp_abs_floor, rel_delta=fp_rel_delta, bg_window=fp_bg_window)
        # map evidence score (strength x count) -> UNID start cost, smoothly: a barely-flagged site
        # (score ~ rel_delta) keeps ~full cost; score >= fp_score_hi gets the reduced cost.
        for (lo, hi, score) in intervals:
            cost = float(np.interp(score, [fp_rel_delta, fp_score_hi],
                                   [params.start_cost_unid, params.unid_discovered_start_cost]))
            discovered.append((lo, hi, cost))
        print('discovered {} footprint(s) ({}): {}'.format(
            len(discovered), discover_mode, [(l, h, round(c, 1)) for l, h, c in discovered]))

    paths, total_ll = viterbi_decode(gpc_pos, obs1, obs0, tfbs_positions, params,
                                     discovered_footprints=discovered)

    read_ids = list(mat.index)
    site_names = tfbs_site_names(tfbs_positions)
    wide = paths_to_wide(paths, total_ll, read_ids, len(tfbs_positions),
                         site_names=site_names, teto_prefix=teto_prefix)
    tidy = paths_to_tidy(paths, read_ids, site_names=site_names)

    # Region-anchored per-molecule columns (tss_nuc, promoter_nuc_gt50, {region}_frac_{type}...).
    # Computed by the SAME function the post-hoc annotator uses, fed the same tidy segments
    # table, so inline and post-hoc numbers are identical by construction.
    if region_entry is not None:
        from annotate_molecule_regions import annotate_one
        reg = annotate_one(tidy, region_entry, promoter_thresh, promoter_min_bp)
        wide = wide.join(reg, how='left')

    diag = {'raw_frac': raw_frac, 'decon_frac': decon_frac, 'discovered': discovered,
            'site_names': site_names}
    return wide, tidy, paths, total_ll, params, diag


def load_yaml_config(path):
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def main():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument('--config', default=None)
    pre_args, _ = pre.parse_known_args()
    config_defaults = load_yaml_config(pre_args.config) if pre_args.config else {}

    p = argparse.ArgumentParser(description='HSMM single-molecule footprint segmentation (v5)')
    p.set_defaults(**config_defaults)

    # I/O
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--segments_output', default=None,
                   help='Tidy long-format per-segment table (default: <output>.segments.txt)')
    p.add_argument('--tfbs_index_output', default=None,
                   help='Sidecar index->name/lo/hi map (default: <output>.tfbs_index.txt)')
    p.add_argument('--teto_name_prefix', default='TetO',
                   help='Positions-file sites whose NAME starts with this (case-insensitive) are '
                        'counted in n_tf_teto; all others in n_tf_other. Set to "" to disable the '
                        'split. Default: TetO')
    p.add_argument('--regions', default=None,
                   help='Optional region file from build_promoter_regions.py. If given, adds '
                        'region-anchored per-molecule columns to the main output: tss_nuc, '
                        'tss_footprint, tss_open, promoter_nuc_gt50 (+ nuc_bases_promoter / '
                        'frac_nuc_promoter / promoter_len), and {region}_bases/frac/any_{type} '
                        'for every region x segment type. Regions are matched by PROMOTER name '
                        '(text before --amplicon_sep), so one 6xTetO entry covers every '
                        'copy-number variant -- valid because the promoter block sits below the '
                        'array start, which is fixed across the series.')
    p.add_argument('--amplicon_sep', default='_opJS4_',
                   help='Separator splitting promoter from amplicon name for --regions lookup')
    p.add_argument('--promoter_nuc_thresh', type=float, default=0.5,
                   help='Fraction of the promoter insert covered by NUC for promoter_nuc_gt50')
    p.add_argument('--promoter_nuc_min_bp', type=int, default=None,
                   help='Also emit an absolute-bp NUC threshold column. RECOMMENDED for '
                        'cross-promoter comparisons: insert lengths are NOT uniform (264 bp for '
                        'most opoBD9 promoters, RPS9 200, minCMV 59), so a fixed FRACTION means '
                        'a different number of bases per promoter.')
    p.add_argument('--plot', required=True)
    p.add_argument('--positions', dest='positions_file', required=True)
    p.add_argument('--amplicon_name', required=True)
    p.add_argument('--motif_track_file', default=None,
                   help='OPTIONAL cosmetic motif annotation track (from export_motif_track.py): '
                        'per-amplicon `>name` blocks of `lo,hi,label[,color]` in the model bp frame, '
                        'overlaid as a strip on the diagnostic + per-read plots. Amplicons absent '
                        'from the file (or --motif_track_file omitted) draw no strip. Never affects '
                        'decoding.')
    p.add_argument('--config', default=None, help='YAML config (defaults; CLI overrides)')

    # preprocessing
    p.add_argument('--reads_to_use', type=int, default=0)
    p.add_argument('--filter_threshold', type=float, default=1.0)
    p.add_argument('--convert_ambiguous_gcgs', default='',
                   help='comma-separated pairs of positions for ambiguous GCG imputation')

    # geometry / grid
    p.add_argument('--amp_width', type=int, default=None)
    p.add_argument('--grid_step', type=int, default=5)

    # conversion
    p.add_argument('--p_t_given_unmeth', type=float, default=0.95)
    p.add_argument('--p_t_given_meth', type=float, default=0.15)

    # protection probs
    p.add_argument('--prob_unmeth_given_open', type=float, default=0.05)
    p.add_argument('--prob_unmeth_given_tf', type=float, default=0.9)
    p.add_argument('--prob_unmeth_given_unid', type=float, default=0.9)
    p.add_argument('--promoter_positions', default='0,0')
    p.add_argument('--prob_unmeth_given_open_promoter', type=float, default=0.5)

    # nucleosome footprint + duration
    # (no --nuc_d_edge: the footprint half-max is tied to the segment half-length; nuc_softness is
    #  a fixed edge-taper/breathing width, see NOTES_v5_calibration_log.md 2026-07-13)
    p.add_argument('--nuc_softness', type=float, default=5.0)
    p.add_argument('--nuc_mode', type=float, default=147.0)
    p.add_argument('--nuc_sigma', type=float, default=25.0)
    p.add_argument('--nuc_min', type=int, default=100)
    p.add_argument('--nuc_max', type=int, default=210)
    p.add_argument('--nuc_min_prot_gpcs', type=int, default=3,
                   help='min observed-protected GpCs for a NUC segment (evidence floor; prevents '
                        'spurious edge/+1 nucleosomes from 1-2 protected GpCs). 0 disables.')

    # TF / UNID geometry
    p.add_argument('--tf_margin', type=int, default=2)
    p.add_argument('--unid_min', type=int, default=15)
    p.add_argument('--unid_max', type=int, default=90)
    p.add_argument('--unid_may_overlap_motifs', action='store_true', default=True,
                   help='allow a UNID segment to span a named positions-file motif. ON by default '
                        'since 2026-08-12 (was OFF; --no_unid_may_overlap_motifs restores that). '
                        'The named site still wins every tie automatically -- '
                        'TF and UNID share prob_unmeth_given_* = 0.9, so a coincident UNID has an '
                        'identical emission and costs 5 nats more -- so this does NOT let UNID steal '
                        'named sites. What it enables is EXTENT REVISION: a footprint wider than the '
                        'annotation can be called, at a cost of ~2.25 protected GpCs of extra '
                        'evidence. TetO array sites stay exclusive unless --unid_overlap_teto.')
    p.add_argument('--no_unid_may_overlap_motifs', dest='unid_may_overlap_motifs',
                   action='store_false',
                   help='restore the pre-2026-08-12 behavior: a UNID may not overlap ANY named motif. '
                        'Needed to reproduce attempt2 / attempt3 JUNB-r3 / RPS9-r2 exactly.')
    p.add_argument('--unid_overlap_teto', action='store_true',
                   help='with --unid_may_overlap_motifs, also let UNID span TetO array sites. NOT '
                        'recommended: a UNID straddling two operators muddies n_tf_teto.')

    # structural costs
    p.add_argument('--start_cost_open', type=float, default=0.5)
    p.add_argument('--start_cost_nuc', type=float, default=2.0)
    p.add_argument('--start_cost_tf', type=float, default=1.0)
    p.add_argument('--start_cost_unid', type=float, default=6.0)
    p.add_argument('--trans_nuc_nuc', type=float, default=1.0)
    p.add_argument('--trans_tf_tf', type=float, default=1.0)
    p.add_argument('--trans_nuc_tf', type=float, default=1.0)
    p.add_argument('--trans_unid_adj', type=float, default=1.0)

    # footprint discovery / position-specific UNID prior (empirical-Bayes across molecules)
    p.add_argument('--discover_footprints', choices=['off', '2pass', '3pass'], default='off',
                   help="off=1-pass, scalar UNID cost (default). 2pass=discover from raw bulk. "
                        "3pass=decode -> nucleosome(+TF)-decontaminated bulk -> discover -> re-decode "
                        "(recommended: removes nucleosome contamination from the footprint prior).")
    p.add_argument('--discovered_footprints_file', default=None,
                   help="use a PRECOMPUTED footprint vocabulary (from discover_footprint_vocabulary.py; "
                        "blocks of `lo,hi,cost` per amplicon). Overrides --discover_footprints: no "
                        "per-amplicon discovery, just applies these intervals' UNID prior. Use this "
                        "for the pooled per-promoter vocabulary (avoids noisy small-N per-amplicon bulks).")
    p.add_argument('--unid_discovered_start_cost', type=float, default=2.0,
                   help='reduced UNID start cost at a STRONG discovered footprint (graded up toward '
                        '--start_cost_unid for weaker sites)')
    p.add_argument('--footprint_abs_floor', type=float, default=0.25,
                   help='discovery: min (decontaminated) bulk protected-fraction to flag a GpC')
    p.add_argument('--footprint_rel_delta', type=float, default=0.20,
                   help='discovery: min excess of bulk protection over local background to flag')
    p.add_argument('--footprint_bg_window', type=int, default=80,
                   help='discovery: +/- bp window for the local-background percentile')
    p.add_argument('--footprint_score_hi', type=float, default=1.0,
                   help='discovery: evidence score (sum of per-GpC excess-over-bg) at/above which a '
                        'site gets the full --unid_discovered_start_cost discount')

    # EM (Viterbi / hard EM) -- fit emission + nucleosome footprint/duration params from data
    p.add_argument('--do_em', action='store_true', default=False,
                   help='fit prob_unmeth_given_{open,tf,unid}, nuc_mode/sigma from the data by '
                        'Viterbi-EM (structural costs, conversion rates, and nuc_softness stay fixed)')
    p.add_argument('--em_max_iters', type=int, default=8)
    p.add_argument('--em_min_obs', type=int, default=100)
    p.add_argument('--em_log', default=None, help='EM convergence TSV (default <output>.em_log.tsv)')
    p.add_argument('--em_fit', default=None,
                   help='comma-separated params EM may update (default: the safe set '
                        'prob_unmeth_given_tf,prob_unmeth_given_unid,nuc_mode). '
                        'prob_unmeth_given_open/nuc_sigma are UNSTABLE under hard EM; only add them '
                        'if you know why. (nuc_softness is fixed; the footprint width is tied to '
                        'the segment length so there is no d_edge to fit.)')

    # plotting
    p.add_argument('--reads_to_plot', type=int, default=10)
    p.add_argument('--individual_reads_to_plot', default='')
    p.add_argument('--bulk_reads', type=int, default=1000,
                   help='reads in the bulk data-vs-prediction sanity panel (0 disables)')

    args = p.parse_args()

    with PdfPages(args.plot) as plots:
        print('loading data')
        mat = load_single_molecule_matrix(args.input)
        gpc_pos = get_methyl_positions(mat)

        pos_file = load_tfbs_positions(args.positions_file)
        tfbs_positions = pos_file.get(args.amplicon_name, [])
        site_names = tfbs_site_names(tfbs_positions)
        idx_out = args.tfbs_index_output or (args.output + '.tfbs_index.txt')
        write_tfbs_index(idx_out, args.amplicon_name, tfbs_positions, site_names,
                         args.teto_name_prefix)
        n_teto = sum(is_teto_site(nm, args.teto_name_prefix) for nm in site_names)
        print('{} site(s): {} TetO, {} other -> {}'.format(
            len(site_names), n_teto, len(site_names) - n_teto, idx_out))

        # Resolve the region entry once, up front, so a typo'd/absent promoter fails loudly
        # here rather than silently producing an output with no region columns.
        region_entry = None
        if args.regions:
            from annotate_molecule_regions import load_regions, build_region_lookup, promoter_of
            all_regions = load_regions(args.regions)
            region_entry = all_regions.get(args.amplicon_name)
            if region_entry is None:
                prom = promoter_of(args.amplicon_name, args.amplicon_sep)
                region_entry = build_region_lookup(all_regions, args.amplicon_sep).get(prom)
            if region_entry is None:
                raise SystemExit(
                    'ERROR: --regions {} has no entry for amplicon {} (promoter {}). Available '
                    'promoters: {}'.format(args.regions, args.amplicon_name,
                                           promoter_of(args.amplicon_name, args.amplicon_sep),
                                           sorted({promoter_of(a, args.amplicon_sep)
                                                   for a in all_regions})))
            print('regions: {} window(s), tss_col={}, promoter=[{},{}]'.format(
                len(region_entry['regions']), region_entry['tss'],
                region_entry['prom_start'], region_entry['prom_end']))

        mat = filter_all_converted_reads(mat, args.filter_threshold)

        if args.reads_to_use > 0 and len(mat) > args.reads_to_use:
            mat = mat.sample(args.reads_to_use, random_state=42)
            print('subsampled to {} reads'.format(args.reads_to_use))

        if args.convert_ambiguous_gcgs:
            fixpos = list(map(int, args.convert_ambiguous_gcgs.split(',')))
            assert len(fixpos) % 2 == 0
            for i in range(len(fixpos) // 2):
                mat = adjust_gcgs(mat, fixpos[2 * i], fixpos[2 * i + 1])

        print('num molecules: {}'.format(len(mat)))
        if len(mat) == 0:
            # write empty outputs and exit cleanly (matches v4 "no reads" behavior)
            pd.DataFrame().to_csv(args.output, sep='\t')
            seg_out = args.segments_output or (args.output + '.segments.txt')
            pd.DataFrame(columns=['read_id', 'seg_index', 'type', 'start', 'end',
                                  'dyad_or_motif', 'motif_name']).to_csv(seg_out, sep='\t',
                                                                         index=False)
            return

        prom = list(map(int, args.promoter_positions.split(',')))

        amp_width = args.amp_width
        if amp_width is None:
            pad = int(args.nuc_max // 2) + 10
            amp_width = int(max(gpc_pos)) + pad
            print('auto-detected amp_width: {}'.format(amp_width))

        params = ModelParams(
            amp_width=amp_width,
            grid_step=args.grid_step,
            p_t_given_unmeth=args.p_t_given_unmeth,
            p_t_given_meth=args.p_t_given_meth,
            prob_unmeth_given_open=args.prob_unmeth_given_open,
            prob_unmeth_given_tf=args.prob_unmeth_given_tf,
            prob_unmeth_given_unid=args.prob_unmeth_given_unid,
            promoter_lo=prom[0], promoter_hi=prom[1],
            prob_unmeth_given_open_promoter=args.prob_unmeth_given_open_promoter,
            nuc_softness=args.nuc_softness,
            nuc_mode=args.nuc_mode, nuc_sigma=args.nuc_sigma,
            nuc_min=args.nuc_min, nuc_max=args.nuc_max,
            nuc_min_prot_gpcs=args.nuc_min_prot_gpcs,
            tf_margin=args.tf_margin,
            unid_min=args.unid_min, unid_max=args.unid_max,
            unid_may_overlap_motifs=args.unid_may_overlap_motifs,
            unid_overlap_teto=args.unid_overlap_teto,
            teto_name_prefix=args.teto_name_prefix,
            start_cost_open=args.start_cost_open, start_cost_nuc=args.start_cost_nuc,
            start_cost_tf=args.start_cost_tf, start_cost_unid=args.start_cost_unid,
            unid_discovered_start_cost=args.unid_discovered_start_cost,
            trans_nuc_nuc=args.trans_nuc_nuc, trans_tf_tf=args.trans_tf_tf,
            trans_nuc_tf=args.trans_nuc_tf, trans_unid_adj=args.trans_unid_adj,
        )

        print('decoding {} molecules ({} TFBS motifs, grid_step={})'.format(
            len(mat), len(tfbs_positions), params.grid_step))
        # decode-affecting and easy to forget which way a round was run -- say it in the log
        _nblock = len(unid_blocked_intervals(tfbs_positions, params.tf_margin, params))
        print('UNID/motif exclusivity: may_overlap={} overlap_teto={} -> {}/{} motifs block UNID'
              .format(params.unid_may_overlap_motifs, params.unid_overlap_teto,
                      _nblock, len(tfbs_positions)))
        em_log_path = args.em_log or (args.output + '.em_log.tsv')
        em_fit_params = [s.strip() for s in args.em_fit.split(',')] if args.em_fit else None
        precomputed = None
        if args.discovered_footprints_file:
            precomputed = load_discovered_footprints_file(args.discovered_footprints_file,
                                                          args.amplicon_name)
            print('loaded {} precomputed footprints for {} from {}'.format(
                len(precomputed), args.amplicon_name, args.discovered_footprints_file))
        wide, tidy, paths, total_ll, params, diag = compute_classifications(
            mat, gpc_pos, tfbs_positions, params,
            do_em=args.do_em, em_max_iters=args.em_max_iters, em_min_obs=args.em_min_obs,
            em_log_path=(em_log_path if args.do_em else None), em_fit_params=em_fit_params,
            discover_mode=args.discover_footprints, fp_abs_floor=args.footprint_abs_floor,
            fp_rel_delta=args.footprint_rel_delta, fp_bg_window=args.footprint_bg_window,
            fp_score_hi=args.footprint_score_hi, precomputed_discovered=precomputed,
            teto_prefix=args.teto_name_prefix, region_entry=region_entry,
            promoter_thresh=args.promoter_nuc_thresh, promoter_min_bp=args.promoter_nuc_min_bp)

        wide.to_csv(args.output, sep='\t', header=True, index=True)
        seg_out = args.segments_output or (args.output + '.segments.txt')
        tidy.to_csv(seg_out, sep='\t', index=False)
        print('wrote {} and {}'.format(args.output, seg_out))

        # plotting
        idx_list = list(mat.index)
        motifs = load_motif_track(args.motif_track_file, args.amplicon_name)
        if args.motif_track_file:
            print('loaded {} motif annotations for {} (overlay)'.format(len(motifs), args.amplicon_name))
        # diagnostic pseudobulk page (always): per-bp NUC/TF/UNID occupancy + raw/decontam bulk
        print('plotting footprint diagnostic panel')
        plot_footprint_diagnostic(mat, np.asarray(gpc_pos, dtype=np.int64), paths, params, plots,
                                  tfbs_positions, diag, motifs=motifs)
        if args.bulk_reads:
            print('plotting bulk sanity panel')
            plot_bulk_sanity(mat, np.asarray(gpc_pos, dtype=np.int64), paths, params, plots,
                             tfbs_positions, n_reads=args.bulk_reads)
        if args.reads_to_plot:
            n_to_plot = min(args.reads_to_plot, len(mat))
            rng = np.random.RandomState(42)
            chosen = rng.choice(len(mat), n_to_plot, replace=False)
            for mi in chosen:
                fig, ax = plt.subplots(figsize=(12, 8))
                plot_single_read(mat, idx_list[mi], ax, fillbetween=tfbs_positions)
                decorate_hsmm_read(ax, paths[mi], tfbs_positions)
                draw_motif_strip(ax, motifs)
                ax.set_title('{}  (ll={:.1f})'.format(idx_list[mi], total_ll[mi]),
                             pad=(46 if motifs else None))
                plots.savefig()
                plt.close()

        if args.individual_reads_to_plot:
            id_to_mi = {rid: mi for mi, rid in enumerate(idx_list)}
            for rid in args.individual_reads_to_plot.split(','):
                if rid in id_to_mi:
                    mi = id_to_mi[rid]
                    fig, ax = plt.subplots(figsize=(12, 8))
                    plot_single_read(mat, rid, ax, fillbetween=tfbs_positions)
                    decorate_hsmm_read(ax, paths[mi], tfbs_positions)
                    draw_motif_strip(ax, motifs)
                    ax.set_title('{}  (ll={:.1f})'.format(rid, total_ll[mi]),
                                 pad=(46 if motifs else None))
                    plots.savefig()
                    plt.close()


if __name__ == '__main__':
    main()
