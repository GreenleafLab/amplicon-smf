#!/usr/bin/env python
"""
PER-MOLECULE region-anchored annotations from v5 HSMM segment calls.

The v5 classifier emits coordinate-free per-molecule summaries (counts `n_nuc`/`n_unid`,
raw intervals `nucs`/`unids`). Neither knows what a TSS or a promoter is. This script
intersects the Viterbi segments against the NAMED windows in a region file from
`build_promoter_regions.py` and emits one row per molecule with region-anchored columns
-- the per-molecule layer that `analyze_footprint_regions.py` (Layer B) computes and then
immediately aggregates away.

Deliberately a POST-HOC script, not a classifier change: the segments files already exist
on disk (thousands of them), so this runs in minutes instead of requiring a full re-decode
of the panel. It is also promoter-agnostic input-wise -- point it at a different region
file and every column follows.

Headline columns (what this was written for):
  tss_nuc            bool  -- a NUC segment covers the TSS column itself
  promoter_nuc_gt50  bool  -- >THRESH of the promoter insert block is covered by NUC

...plus the supporting numbers, because a bare boolean hides a real comparability trap:
promoter insert LENGTHS ARE NOT UNIFORM (264 bp for most opoBD9 promoters, RPS9 200,
**minCMV 59**). At 59 bp, "50% covered" is 30 bp -- which almost any nucleosome grazing
the block satisfies -- whereas for a 264 bp promoter it is 132 bp, about a full
nucleosome. So `promoter_nuc_gt50` is NOT comparable between minCMV and the rest. Always
carry `promoter_len` / `nuc_bases_promoter` alongside it, and prefer an absolute-bp
threshold (`--promoter_nuc_min_bp`) if you need a cross-promoter statement.

Generic columns, for every region R in the region file x every segment type T:
  {R}_bases_{T}   bp of R covered by segments of type T
  {R}_frac_{T}    that over len(R)
  {R}_any_{T}     bool, >=1 bp overlap
This means `plus1_nuc_frac_NUC`, `TATA_frac_UNID`, etc. all come for free.

Coordinate frame: model bp = matrix column, same as the region file and the segments file.
Segment intervals are treated as INCLUSIVE [start, end] (matching `analyze_footprint_regions.py`
and `junb_cluster_states.py`, which both do `min(b,hi) - max(a,lo) + 1`).

Copy-number variants: region files are emitted for `_opJS4_6xTetO` only, but the promoter
block is coordinate-conserved across the TetO series (the array starts at a fixed low coord
-- JUNB 481, RPS9 417 -- and grows upward, so everything below it is identical). Amplicons
are therefore matched to regions by PROMOTER NAME, and every region is asserted to lie below
the amplicon's array start. Use --strict_amplicon to require an exact amplicon-name match.

Usage:
  python annotate_molecule_regions.py \
      --regions 260804_opoBD9_promoter_regions.txt \
      --model_dir 260720_modelrunning_v5_attempt2 \
      --output molecule_regions.tidy.txt.gz
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

SEG_TYPES = ('NUC', 'TF', 'UNID', 'OPEN')


# ---- region file ----------------------------------------------------------------
def load_regions(path):
    """-> {amplicon: {'tss': int|None, 'prom_start':int|None, 'prom_end':int|None,
                      'regions': [(name, lo, hi, type)]}}

    Same format/parser semantics as analyze_footprint_regions.load_regions, plus the
    prom_start/prom_end comment added by build_promoter_regions.py (2026-08-04).
    """
    out = {}
    cur = None
    for line in open(path):
        line = line.rstrip('\n')
        if not line or (line.startswith('#') and cur is None):
            continue
        if line.startswith('>'):
            cur = line[1:].strip()
            out[cur] = {'tss': None, 'tata': None, 'conf': None,
                        'prom_start': None, 'prom_end': None, 'regions': []}
        elif line.startswith('#'):
            for key, field in (('tss_col', 'tss'), ('tata_col', 'tata'),
                               ('prom_start', 'prom_start'), ('prom_end', 'prom_end')):
                m = re.search(r'{}=(-?\d+)'.format(key), line)
                if m:
                    out[cur][field] = int(m.group(1))
            m = re.search(r'confidence=(\w+)', line)
            if m:
                out[cur]['conf'] = m.group(1)
        else:
            parts = [p.strip() for p in line.split(',')]
            name, lo, hi, styp = parts[0], parts[1], parts[2], parts[3]
            rule, rparams = parse_rule(parts[4] if len(parts) > 4 else '',
                                       parts[5] if len(parts) > 5 else '',
                                       int(lo), int(hi), name)
            out[cur]['regions'].append((name, int(lo), int(hi), styp, rule, rparams))
    return out


# ---- region call rules -----------------------------------------------------------
# A region line is  name,lo,hi,segment_type[,RULE[,params]]  -- RULE defaults to MIDPOINT so
# every pre-2026-08-12 region file keeps working untouched.
#
#   MIDPOINT  segment CENTER falls in [lo,hi].            params: none
#   EDGES     segment START in [s_lo,s_hi] AND END in [e_lo,e_hi].
#             params: `tol=N` (symmetric +/-N around lo and hi; default 8), or
#                     `start=a:b;end=c:d` for explicit asymmetric windows.
#   GPCS      segment SPAN covers every `require` position and no `exclude` position.
#             params: `require=145+159+164+171;exclude=130`   ('+' separates; ',' is taken)
#
# WHY GPCS EXISTS: segment boundaries live on a 5 bp grid (grid_step) and, inside a GpC desert,
# are unidentifiable -- the decode has no evidence preferring one grid point over another
# between two GpCs. An EDGES rule in bp therefore partly scores quantization noise wherever
# GpCs are sparse (at FTH1 the spacing either side of col 145 is 14-15 bp). GPCS states the
# same distinction in the units the data actually has, and it fails loudly: if two hypotheses
# differ only across a desert you cannot write the rule, which is the correct answer rather
# than a coin flip.
RULES = ('MIDPOINT', 'EDGES', 'GPCS')
EDGE_TOL_DEFAULT = 8


def parse_rule(rule, params, lo, hi, name):
    """-> (RULE, params_dict). Empty/absent rule => MIDPOINT (backward compatible)."""
    rule = (rule or 'MIDPOINT').strip().upper()
    if rule not in RULES:
        raise ValueError('region {!r}: unknown rule {!r} (expected one of {})'
                         .format(name, rule, '/'.join(RULES)))
    kv = {}
    for tok in (params or '').split(';'):
        tok = tok.strip()
        if not tok:
            continue
        if '=' not in tok:
            raise ValueError('region {!r}: bad rule param {!r} (want key=value)'.format(name, tok))
        k, v = tok.split('=', 1)
        kv[k.strip().lower()] = v.strip()

    out = {}
    if rule == 'EDGES':
        if 'start' in kv or 'end' in kv:
            if not ('start' in kv and 'end' in kv):
                raise ValueError('region {!r}: EDGES needs BOTH start= and end= (or just tol=)'
                                 .format(name))
            for key in ('start', 'end'):
                a, b = kv[key].split(':')
                out[key] = (int(a), int(b))
        else:
            tol = int(kv.get('tol', EDGE_TOL_DEFAULT))
            out['start'] = (lo - tol, lo + tol)
            out['end'] = (hi - tol, hi + tol)
    elif rule == 'GPCS':
        req = [int(x) for x in kv.get('require', '').split('+') if x != '']
        exc = [int(x) for x in kv.get('exclude', '').split('+') if x != '']
        if not req:
            raise ValueError('region {!r}: GPCS needs require=p1+p2+...'.format(name))
        out['require'], out['exclude'] = req, exc
    return rule, out


def rule_suffix(rule):
    """Column-name fragment, so all three rules can coexist across a region file."""
    return {'MIDPOINT': 'mid', 'EDGES': 'edges', 'GPCS': 'gpcs'}[rule]


def segment_satisfies(rule, rparams, lo, hi, a, b):
    """Does segment [a,b] (inclusive) satisfy `rule` for region [lo,hi]?"""
    if rule == 'MIDPOINT':
        return lo <= 0.5 * (a + b) <= hi
    if rule == 'EDGES':
        s_lo, s_hi = rparams['start']
        e_lo, e_hi = rparams['end']
        return (s_lo <= a <= s_hi) and (e_lo <= b <= e_hi)
    if rule == 'GPCS':
        return (all(a <= p <= b for p in rparams['require']) and
                not any(a <= p <= b for p in rparams['exclude']))
    raise ValueError(rule)


def promoter_of(amplicon, sep='_opJS4_'):
    """`JUNB_opJS4_2xTetO` -> `JUNB`. Falls back to the whole name if the separator is absent."""
    return amplicon.split(sep)[0] if sep in amplicon else amplicon


def build_region_lookup(regions, sep='_opJS4_'):
    """promoter name -> region entry, so copy-number variants share the 6x definitions."""
    lut = {}
    for amp, entry in regions.items():
        prom = promoter_of(amp, sep)
        if prom in lut:
            print('WARNING: two region entries map to promoter {} ({} and the earlier one); '
                  'keeping the first'.format(prom, amp), file=sys.stderr)
            continue
        lut[prom] = entry
    return lut


# ---- per-molecule intersection --------------------------------------------------
def overlap_bp(a, b, lo, hi):
    """Inclusive-interval overlap in bp between segment [a,b] and region [lo,hi]."""
    return max(0, min(b, hi) - max(a, lo) + 1)


def annotate_one(seg, entry, promoter_thresh=0.5, promoter_min_bp=None,
                 region_thresh=0.5):
    """
    seg: DataFrame for ONE (sample, amplicon) -- columns read_id, type, start, end.
    entry: the region-file entry for this promoter.
    -> DataFrame indexed by read_id.
    """
    reg_list = entry['regions']
    tss = entry['tss']

    # Region name -> (lo, hi, length). Length uses inclusive bounds to match overlap_bp.
    # reg_list entries are (name, lo, hi, seg_type, rule, rule_params); the last two default
    # to ('MIDPOINT', {}) for region files written before rules existed.
    reg_span = {r[0]: (r[1], r[2], r[2] - r[1] + 1) for r in reg_list}
    reg_rule = {r[0]: (r[4] if len(r) > 4 else 'MIDPOINT',
                       r[5] if len(r) > 5 else {}) for r in reg_list}

    types = seg['type'].astype(str).values
    starts = seg['start'].astype(int).values
    ends = seg['end'].astype(int).values
    reads = seg['read_id'].astype(str).values

    uniq_reads, read_idx = np.unique(reads, return_inverse=True)
    n = len(uniq_reads)

    # bases[region][type] accumulated per molecule
    acc = {name: {t: np.zeros(n, dtype=np.int32) for t in SEG_TYPES} for name in reg_span}
    # rule_hit[region][type] -- does ANY segment of that type satisfy the region's call rule?
    rule_hit = {name: {t: np.zeros(n, dtype=bool) for t in SEG_TYPES} for name in reg_span}
    # ...and the EXTENT of the segment that satisfied it. This is the point of a require-only
    # GPCS rule: write one region for "any footprint covering these GpCs", then read the width
    # distribution post-hoc instead of pre-enumerating a state per extent. Without these the
    # rule is only a boolean and that workflow is impossible.
    # A segmentation TILES, so for a GPCS rule with >=1 required position at most ONE segment of
    # a given type can satisfy it -- the match is unique. For MIDPOINT/EDGES two segments of the
    # same type could in principle both qualify; keep the WIDEST and count the collisions.
    rule_lo = {name: {t: np.full(n, -1, dtype=np.int32) for t in SEG_TYPES} for name in reg_span}
    rule_hi = {name: {t: np.full(n, -1, dtype=np.int32) for t in SEG_TYPES} for name in reg_span}
    rule_multi = {name: {t: np.zeros(n, dtype=np.int16) for t in SEG_TYPES} for name in reg_span}
    tss_cov = {t: np.zeros(n, dtype=bool) for t in SEG_TYPES}

    for i in range(len(reads)):
        t = types[i]
        if t not in SEG_TYPES:
            continue
        a, b, m = starts[i], ends[i], read_idx[i]
        for name, (lo, hi, _L) in reg_span.items():
            ov = overlap_bp(a, b, lo, hi)
            if ov:
                acc[name][t][m] += ov
            rule, rparams = reg_rule[name]
            # NOT gated on `ov`: an EDGES/GPCS rule can legitimately be evaluated against a
            # segment that grazes or misses [lo,hi] -- gating on overlap would silently make
            # every rule an `_any_` rule underneath.
            if segment_satisfies(rule, rparams, lo, hi, a, b):
                rule_hit[name][t][m] = True
                rule_multi[name][t][m] += 1
                if (b - a) > (rule_hi[name][t][m] - rule_lo[name][t][m]):
                    rule_lo[name][t][m], rule_hi[name][t][m] = a, b
        if tss is not None and a <= tss <= b:
            tss_cov[t][m] = True

    df = pd.DataFrame(index=pd.Index(uniq_reads, name='read_id'))
    for name, (lo, hi, L) in reg_span.items():
        for t in SEG_TYPES:
            bases = acc[name][t]
            frac = bases / float(L)
            df['{}_bases_{}'.format(name, t)] = bases
            df['{}_frac_{}'.format(name, t)] = frac
            df['{}_any_{}'.format(name, t)] = bases > 0
            # Majority-coverage boolean (--region_thresh). Prefer this over `_any_` for
            # abundance: `_any_` is the >=1bp overlap rule, which OVER-COUNTS when a wide
            # adjacent footprint spills across the edge -- the +TRP PIC clipping the pause
            # window, which is exactly why analyze_footprint_regions.py switched its
            # aggregate statistic to center-in-window. Like `promoter_nuc_gt50` this is a
            # fraction-OF-REGION threshold, so it is not comparable between regions of very
            # different length.
            df['{}_gt{:g}_{}'.format(name, region_thresh * 100, t)] = frac > region_thresh
            # The region's own call rule (MIDPOINT default / EDGES / GPCS). This is the
            # column to prefer for abundance: `_mid_` is the center-in-window rule the
            # aggregate path already uses, and `_edges_`/`_gpcs_` express extent hypotheses
            # that no coverage fraction can distinguish.
            rs = rule_suffix(reg_rule[name][0])
            df['{}_{}_{}'.format(name, rs, t)] = rule_hit[name][t]
            # extent of the satisfying segment; -1 / 0 where the rule did not fire
            width = np.where(rule_hit[name][t],
                             rule_hi[name][t] - rule_lo[name][t] + 1, 0).astype(np.int32)
            df['{}_{}_{}_lo'.format(name, rs, t)] = rule_lo[name][t]
            df['{}_{}_{}_hi'.format(name, rs, t)] = rule_hi[name][t]
            df['{}_{}_{}_width'.format(name, rs, t)] = width
            if rule_multi[name][t].max() > 1:
                print('NOTE: region {!r} rule {} matched >1 {} segment on {} molecule(s); '
                      'kept the widest'.format(name, reg_rule[name][0], t,
                                               int((rule_multi[name][t] > 1).sum())),
                      file=sys.stderr)

    # ---- headline columns -------------------------------------------------------
    # TSS covered by a nucleosome: a NUC segment spans the TSS column itself. Chosen over
    # "overlaps the TSS_Inr window" because it needs no window-width convention; the
    # window version is still available as TSS_Inr_frac_NUC / TSS_Inr_any_NUC.
    df['tss_nuc'] = tss_cov['NUC'] if tss is not None else False
    # Companion state of the TSS base, so "not nucleosomal" isn't silently read as "naked".
    # A TSS under a UNID/TF footprint is protected but mechanistically the opposite of closed.
    df['tss_footprint'] = tss_cov['UNID'] | tss_cov['TF']
    df['tss_open'] = ~(df['tss_nuc'].values | df['tss_footprint'].values)

    if 'promoter' in reg_span:
        L = reg_span['promoter'][2]
        nuc_bp = acc['promoter']['NUC']
        df['promoter_len'] = L
        df['nuc_bases_promoter'] = nuc_bp
        df['frac_nuc_promoter'] = nuc_bp / float(L)
        df['promoter_nuc_gt50'] = df['frac_nuc_promoter'].values > promoter_thresh
        if promoter_min_bp is not None:
            df['promoter_nuc_min_bp'] = nuc_bp > promoter_min_bp
    else:
        print('WARNING: no `promoter` region in the region file -- promoter columns skipped. '
              'Regenerate it with the 2026-08-04 build_promoter_regions.py.', file=sys.stderr)

    if tss is not None:
        df['tss_col'] = tss
    df['conf'] = entry.get('conf')
    return df


def check_regions_below_array(entry, amplicon, array_lo):
    """Regions come from the 6x amplicon; verify they sit below THIS amplicon's array."""
    if array_lo is None:
        return
    bad = [r[0] for r in entry['regions'] if r[2] >= array_lo]
    if bad:
        print('WARNING: {} -- region(s) {} extend into/past the TetO array start ({}); '
              'coordinate conservation across copy-number variants does not hold there.'
              .format(amplicon, ','.join(bad), array_lo), file=sys.stderr)


def load_array_starts(positions_file):
    """amplicon -> lowest TFBS start (the array start), or {} if no file given."""
    if not positions_file:
        return {}
    out, cur = {}, None
    for line in open(positions_file):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            cur = line[1:].strip()
        else:
            lo = int(line.split(',')[0])
            out[cur] = min(out.get(cur, lo), lo)
    return out


# ---- driver ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description='Per-molecule region-anchored annotations from v5 segments.')
    ap.add_argument('--regions', required=True,
                    help='region file from build_promoter_regions.py')
    ap.add_argument('--model_dir', required=True,
                    help='dir containing {sample}/{sample}.{amplicon}.segments*.txt')
    ap.add_argument('--segments_glob', default='*/*.segments*.txt',
                    help='glob under --model_dir (default: %(default)s)')
    ap.add_argument('--positions', default=None,
                    help='optional positions file; used to verify regions sit below the '
                         'TetO array start for each copy-number variant')
    ap.add_argument('--samples', nargs='+', default=None,
                    help='restrict to these sample names (default: all found)')
    ap.add_argument('--amplicons', nargs='+', default=None,
                    help='restrict to these amplicon names (default: all found)')
    ap.add_argument('--amplicon_sep', default='_opJS4_',
                    help='promoter/amplicon name separator (default: %(default)s)')
    ap.add_argument('--strict_amplicon', action='store_true',
                    help='require an exact amplicon match in the region file instead of '
                         'matching by promoter name across copy-number variants')
    ap.add_argument('--promoter_thresh', type=float, default=0.5,
                    help='fraction of the promoter insert covered by NUC for '
                         'promoter_nuc_gt50 (default: %(default)s)')
    ap.add_argument('--promoter_nuc_min_bp', type=int, default=None,
                    help='ALSO emit an absolute-bp threshold column (recommended for '
                         'cross-promoter comparisons -- insert lengths are not uniform)')
    ap.add_argument('--region_thresh', type=float, default=0.5,
                    help='fraction-of-region coverage for the generic {R}_gt{N}_{T} '
                         'majority booleans (default: %(default)s)')
    ap.add_argument('--output', required=True, help='tidy per-molecule table (.txt/.txt.gz)')
    args = ap.parse_args()

    regions = load_regions(args.regions)
    lut = build_region_lookup(regions, args.amplicon_sep)
    array_lo = load_array_starts(args.positions)
    print('Loaded {} region entries -> {} promoters'.format(len(regions), len(lut)))

    paths = sorted(glob.glob(os.path.join(args.model_dir, args.segments_glob)))
    if not paths:
        sys.exit('No segments files matched {}/{}'.format(args.model_dir, args.segments_glob))
    print('Found {} segments files'.format(len(paths)))

    frames = []
    n_skip_noregion, n_skip_empty = 0, 0
    for p in paths:
        base = os.path.basename(p)
        # {sample}.{amplicon}.segments[.v5].txt
        m = re.match(r'^(.+?)\.(.+?)\.segments.*\.txt$', base)
        if not m:
            continue
        sample, amplicon = m.group(1), m.group(2)
        if args.samples and sample not in args.samples:
            continue
        if args.amplicons and amplicon not in args.amplicons:
            continue

        if args.strict_amplicon:
            entry = regions.get(amplicon)
        else:
            entry = lut.get(promoter_of(amplicon, args.amplicon_sep))
        if entry is None:
            n_skip_noregion += 1
            continue
        check_regions_below_array(entry, amplicon, array_lo.get(amplicon))

        seg = pd.read_csv(p, sep='\t')
        if seg.empty:
            n_skip_empty += 1
            continue

        df = annotate_one(seg, entry, args.promoter_thresh, args.promoter_nuc_min_bp,
                          args.region_thresh)
        df.insert(0, 'amplicon', amplicon)
        df.insert(0, 'sample', sample)
        df.insert(0, 'promoter', promoter_of(amplicon, args.amplicon_sep))
        frames.append(df.reset_index())

    if not frames:
        sys.exit('No (sample, amplicon) pairs produced rows.')

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(args.output, sep='\t', index=False,
               compression=('gzip' if args.output.endswith('.gz') else None))

    print('\nWrote {}  ({} molecules, {} sample x amplicon pairs)'
          .format(args.output, len(out), len(frames)))
    if n_skip_noregion:
        print('  skipped {} file(s): no region entry for that promoter'.format(n_skip_noregion))
    if n_skip_empty:
        print('  skipped {} file(s): empty segments table'.format(n_skip_empty))
    print('\nHeadline rates (all molecules pooled -- per-promoter breakdown below):')
    for c in ('tss_nuc', 'tss_footprint', 'tss_open', 'promoter_nuc_gt50'):
        if c in out:
            print('  {:20s} {:.3f}'.format(c, out[c].mean()))
    if 'promoter_nuc_gt50' in out:
        print('\nPer promoter (NOTE: promoter_len differs -- gt50 is not comparable across '
              'promoters of different length):')
        summ = out.groupby('promoter').agg(
            n=('read_id', 'size'), promoter_len=('promoter_len', 'first'),
            tss_nuc=('tss_nuc', 'mean'), nuc_bp=('nuc_bases_promoter', 'mean'),
            gt50=('promoter_nuc_gt50', 'mean'))
        print(summ.round(3).to_string())


if __name__ == '__main__':
    main()
