#!/usr/bin/env python
"""
Build a per-amplicon PROMOTER REGION file for the v5 footprint-region analysis
(`analyze_footprint_regions.py`).

Layer A of the footprint-region tool (see the project's footprint-region-tool
plan). Regions are defined two ways and materialized to model-bp columns so the
downstream aggregator never has to know about annotations:

  1. TSS-relative TEMPLATE regions (apply to every promoter given its TSS):
       core_PIC, TSS_Inr, pause   (see REGION_TEMPLATE)
  2. ANNOTATION-anchored regions (only where the feature exists):
       TATA      (from the promoter's TATA bounds / motif-track TATA)
       promoter  (the promoter INSERT block, [prom_start, prom_start+Length-1];
                  always emitted -- the TSS sits inside it, near the low end)
       plus1_nuc (the constant downstream backbone, [0, prom_start-1]) -- ABSOLUTE
                  construct coords, NOT TSS-relative as of 2026-08-04; the data say
                  this nucleosome is construct-positioned. See the emission site.

TSS is chosen per promoter with a transparent priority and a CONFIDENCE flag,
NOT silently guessed:
  motif-track Inr  ->  derived TSS_Jason  ->  derived TSS_CAGE
where derived col = prom_end - annotation, prom_end = PROM_START + length.
Confidence is 'high' when the chosen TSS agrees with an independent source
within TSS_AGREE_BP, else 'low' (dispersed / discordant promoter -- the
PIC/pause region concept is weaker there anyway). All alternatives and the
per-promoter provenance are written as comments so they are reviewable.

Output format (mirrors the positions / motif_track files):
    >AMPLICON
    # tss_col=.. tata_col=.. confidence=.. tss_source=.. alts=..
    region_name,lo_col,hi_col,segment_type
    ...

`segment_type` is the footprint type that region is about (UNID/NUC/TF), so the
aggregator knows what to count; it still computes every type in every region.

Coordinate frame: model bp = matrix column. dist_to_TSS = tss_col - col
(upstream/TetO negative, gene body/downstream positive). A template region given
as (dist_lo, dist_hi) maps to cols [tss_col - dist_hi, tss_col - dist_lo].
"""
import argparse
import re
import sys
import pandas as pd

PROM_START_DEFAULT = 136
TSS_AGREE_BP = 6          # |a-b| <= this  => sources 'agree'
TATA_AGREE_BP = 12        # derived vs motif-track TATA agreement tolerance

# TSS-relative template. (dist_lo, dist_hi, segment_type). dist upstream negative.
# pause low edge shifted to +30 so a wide core/PIC footprint spilling past the TSS
# (median ~50 bp under initiation stall) no longer clips the pause window edge.
REGION_TEMPLATE = [
    ("core_PIC",  -40,  15, "UNID"),   # PIC footprint over core promoter, incl. spillover past TSS
    ("TSS_Inr",    -8,   8, "UNID"),   # right at the initiator
    ("pause",      30,  75, "UNID"),   # promoter-proximal paused Pol II
]
# plus1_nuc is NOT TSS-relative -- see the block comment where it is emitted below.
# TATA is anchored (not template): tata_col +/- TATA_HALF
TATA_HALF = 8


def num(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def parse_motif_track(path):
    """-> {promoter: [(lo,hi,label), ...]} keyed by the promoter part of the 6xTetO header."""
    mt = {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            m = re.match(r">(\w+)_opJS4_6xTetO$", line)
            cur = m.group(1) if m else None
            if cur:
                mt.setdefault(cur, [])
        elif cur and line.strip():
            p = line.split(",")
            mt[cur].append((int(p[0]), int(p[1]), p[2]))
    return mt


def anchor(mt, prom, keys):
    """center col of the first motif whose label contains any key (case-insensitive)."""
    for lo, hi, lab in mt.get(prom, []):
        if any(k.lower() in lab.lower() for k in keys):
            return (lo + hi) // 2
    return None


def amplicons_for(positions_path):
    """-> {promoter: [amplicon headers]} from the positions file."""
    out = {}
    for line in open(positions_path):
        if line.startswith(">"):
            amp = line[1:].strip()
            m = re.match(r"(\w+?)_opJS4_", amp)
            if m:
                out.setdefault(m.group(1), []).append(amp)
    return out


def choose_tss(inr, tj, tc):
    """priority TSS_Jason -> TSS_CAGE -> motif_Inr; confidence high if the chosen value
    agrees with an independent source.

    PRIORITY CHANGED 2026-08-05. It was motif_Inr FIRST, on the reasoning that an Inr is a
    direct bp-resolution positional element while TSS_Jason is a derived annotation. That
    was backwards in practice, for two reasons:

      1. TSS_Jason is not an annotation guess -- it is Jason's empirical 5'-read-end TSS
         mapping, i.e. measured initiation. Checked against the histogram mode: BAX
         Jason=230 vs mode 229, LYL1 Jason=235 vs mode 234 (1 bp each).
      2. The "motif Inr" comes from the ENCODE BED track, projected through genomic
         coordinates, and that source is demonstrably unreliable on this panel -- it
         shifted every RPS9 motif by +30 bp, and on LYL1 it labels two minus-strand GATA
         sites as "TA-Inr".

    Silently overriding measured initiation with that produced three materially wrong
    anchors: BAX 236 vs 170 (+66), FTH1GABPA 162 vs 189 (-27), LYL1 189 vs 165 (+24).
    Now Jason wins and a disagreeing Inr is FLAGGED rather than used.
    """
    cands = [("TSS_Jason", tj), ("TSS_CAGE", tc), ("motif_Inr", inr)]
    cands = [(s, v) for s, v in cands if v is not None]
    if not cands:
        return None, "none", "no_tss_annotation"
    src, val = cands[0]
    others = [v for s, v in cands[1:]]
    conf = "high" if any(abs(val - o) <= TSS_AGREE_BP for o in others) else "low"
    alts = " ".join(f"{s}={v}" for s, v in cands)
    note = f"{src}(conf={conf}); alts: {alts}"
    if inr is not None and abs(val - inr) > TSS_AGREE_BP:
        note += (f" | INR-DISCORDANT-review: motif_Inr={inr} disagrees by {inr - val:+d} bp "
                 f"and was NOT used (encode-derived Inr calls are unreliable on this panel)")
    return val, conf, note


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations", required=True, help="promoter annotations TSV")
    ap.add_argument("--motif_track", required=True, help="model-frame curated motif track")
    ap.add_argument("--positions", required=True, help="positions file (source of amplicon list)")
    ap.add_argument("--output", required=True, help="output region file")
    ap.add_argument("--prom_start", type=int, default=PROM_START_DEFAULT)
    ap.add_argument("--amplicon_suffix", default="_opJS4_6xTetO",
                    help="only emit regions for {promoter}{suffix} amplicons")
    args = ap.parse_args()

    ann = pd.read_csv(args.annotations, sep="\t", dtype=str)
    mt = parse_motif_track(args.motif_track)
    amp_map = amplicons_for(args.positions)

    diag = []
    with open(args.output, "w") as out:
        out.write("# per-amplicon promoter regions for analyze_footprint_regions.py\n")
        out.write(f"# generated by build_promoter_regions.py; prom_start={args.prom_start}\n")
        out.write("# region line: name,lo_col,hi_col,segment_type   (model-bp cols)\n")
        for _, r in ann.iterrows():
            prom = r["Promoter"]
            L = num(r["Length"])
            if L is None:
                continue
            amp = f"{prom}{args.amplicon_suffix}"
            if amp not in amp_map.get(prom, []):
                continue
            pe = args.prom_start + L
            tj = num(r.get("TSS_Jason"))
            tc = num(r.get("TSS_CAGE"))
            tb = num(r.get("TATA bounds"))
            inr = anchor(mt, prom, ["inr", "ca-", "ta-"])
            tss_col, conf, prov = choose_tss(inr,
                                             pe - tj if tj is not None else None,
                                             pe - tc if tc is not None else None)
            if tss_col is None:
                diag.append((prom, "SKIP: no TSS"))
                continue

            # TATA: prefer derived (TATA bounds) but cross-check against motif-track TATA
            tata_der = pe - tb if (tb is not None and tb > 0) else None
            tata_mt = anchor(mt, prom, ["tbp", "tata"])
            tata_col = None
            tata_note = "none"
            if tata_der is not None:
                tata_col = tata_der
                if tata_mt is not None:
                    d = tata_der - tata_mt
                    tata_note = f"derived={tata_der} motif={tata_mt} diff={d}"
                    if abs(d) > TATA_AGREE_BP:
                        tata_note += " DISCORDANT-review"
                else:
                    tata_note = f"derived={tata_der} (no motif TATA)"
            elif tata_mt is not None:
                tata_col = tata_mt
                tata_note = f"motif={tata_mt} (no derived)"

            out.write(f">{amp}\n")
            out.write(f"# tss_col={tss_col} tata_col={tata_col} confidence={conf} "
                      f"init={r.get('Initiation')} | tss_src: {prov} | tata: {tata_note}\n")
            # prom_start/prom_end/prom_len are the promoter INSERT block (the cloned promoter
            # sequence itself), not a TSS-relative window -- downstream needs them to normalize
            # "fraction of the promoter covered by X" without re-reading the annotations TSV.
            # prom_end is the INCLUSIVE last column (pe-1), matching the `promoter` region line
            # and Layer B's inclusive overlap arithmetic. Deliberately not the exclusive
            # `prom_start+L` used internally above, so there is only one convention in the file.
            out.write(f"# prom_start={args.prom_start} prom_end={pe - 1} prom_len={L}\n")
            regions = []
            for name, dlo, dhi, styp in REGION_TEMPLATE:
                clo = tss_col - dhi
                chi = tss_col - dlo
                regions.append((name, clo, chi, styp))
            if tata_col is not None:
                regions.append(("TATA", tata_col - TATA_HALF, tata_col + TATA_HALF, "UNID"))
            # The promoter insert block, as an ordinary region so Layer B treats it like any
            # other window. Anchored (not TSS-relative): the TSS sits INSIDE it, typically near
            # the low/downstream end (upstream = HIGH cols), e.g. JUNB tss=168 in [136,400].
            # NOTE lengths are NOT uniform (264 bp for most, RPS9 200, minCMV 59) -- see the
            # caveat in annotate_molecule_regions.py about thresholding a *fraction* of this.
            # `pe - 1`: region bounds are INCLUSIVE downstream (Layer B counts
            # min(b,hi)-max(a,lo)+1), so [prom_start, pe] would span L+1 bp. The insert is
            # exactly L bp, and here L is a real sequence length we can check against, so
            # the usual template sloppiness (a (dlo,dhi) window spans dhi-dlo+1) isn't
            # acceptable -- it would inflate every promoter fraction by 1/L.
            regions.append(("promoter", args.prom_start, pe - 1, "NUC"))
            # plus1_nuc -- ABSOLUTE construct coordinates [0, prom_start-1], NOT TSS-relative.
            #
            # It was TSS-relative (+20..+170) until 2026-08-04. Measured across 16 promoters x all
            # copy-number variants (498,302 NUC segments), the nucleosome's interior boundary --
            # the MEASURED quantity, not the prior-regularized dyad -- is pinned to absolute
            # column ~115 (mean 108.4, sd 10.1, 12/16 promoters median exactly 115), whereas in
            # TSS coordinates it scatters (mean -59.9, sd 23.6, range -121..-31). 2.3x tighter in
            # absolute coords. Peak dyad column is 46-51 for EVERY promoter while tss_col ranges
            # 136-236, so the apparent "+1 distance" (90 bp for minCMV, 190 for BAX) was purely an
            # artifact of promoter insert length. There is no separate peak at the canonical +120
            # for high-tss_col promoters => no distinct TSS-positioned +1 in this data.
            #
            # Mechanistically: cols [0, prom_start) are the CONSTANT downstream backbone shared by
            # every construct (the insert is [prom_start, prom_start+L]). The nucleosome fills that
            # backbone and is bounded ~21 bp below the insert start. Construct-positioned, not
            # TSS-positioned. The name is kept for continuity with existing P110 analyses.
            #
            # Defined on segment start/end (always clamped to [0, amp_width]), NOT on the dyad --
            # deliberately: off-low-edge dyads are negative AND censored, piling up on the
            # left_pad floor at -35 (>=75% of negatives sit exactly there), so they are not a
            # usable continuous coordinate.
            regions.append(("plus1_nuc", 0, args.prom_start - 1, "NUC"))
            for name, clo, chi, styp in regions:
                out.write(f"{name},{clo},{chi},{styp}\n")
            diag.append((prom, f"tss={tss_col} conf={conf} tata={tata_col} [{tata_note}]"))

    print(f"Wrote {args.output}\n")
    print("Per-promoter provenance (review 'low'/'DISCORDANT' before trusting):")
    for prom, note in diag:
        flag = "  <-- REVIEW" if ("low" in note or "DISCORDANT" in note or "SKIP" in note) else ""
        print(f"  {prom:12s} {note}{flag}")


if __name__ == "__main__":
    main()
