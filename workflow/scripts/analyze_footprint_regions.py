#!/usr/bin/env python
"""
Per-region v5 footprint frequencies across promoters x samples (Layer B).

Consumes the tidy v5 `.segments.v5.txt` (one row per segment: read_id, type in
{OPEN,NUC,TF,UNID}, start, end) for each (sample, amplicon), plus a per-amplicon
REGION file from `build_promoter_regions.py`, and a set of samples whose
condition/batch is parsed from the directory name. For each
(promoter, amplicon, condition, batch|combined, region, segment_type) it reports:

  - abundance, two DEFINITIONS x two DENOMINATORS:
      overlap : molecule has >=1 segment of the region's type overlapping the region (>=1bp)
      center  : ... whose CENTER falls in the region  (robust to a wide adjacent
                footprint spilling across the region edge -- the generalization of
                hand-shifting a window bound)
      denom 'all' = all decoded molecules ; denom 'inf' = molecules whose decoded
      span (first..last segment) overlaps the region.
  - width: mean / median bp of the region's overlapping segments.
  - per-bp occupancy vectors (for the profile plots).

Everything is model-bp / matrix-column frame; plots are reindexed to distance-to-TSS
(tss_col from the region-file comment; dist = tss_col - col, upstream negative).

Outputs:
  {outdir}/footprint_regions.tidy.csv         -- the long table (the reusable artifact)
  {outdir}/by_promoter/{promoter}.pdf          -- per-promoter profile + region bars
"""
import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# ---- condition parsing / display ----------------------------------------------
DRUGS = ["A485", "BRM", "TRP"]
DRUG_COLOR = {"": "#333333", "A485": "#4daf4a", "BRM": "#984ea3", "TRP": "#e41a1c"}

def parse_sample(dirname):
    """L####_..._{minusDox|plusDox}[_DRUG]_b{1,2} -> (condition, dox, drug, batch)."""
    dox = "plus" if "plusDox" in dirname else ("minus" if "minusDox" in dirname else "?")
    drug = ""
    for d in DRUGS:
        if f"_{d}_" in dirname or dirname.endswith(f"_{d}") or f"_{d}_b" in dirname:
            drug = d
            break
    m = re.search(r"_b(\d)\b", dirname) or re.search(r"_b(\d)$", dirname)
    batch = f"b{m.group(1)}" if m else "b?"
    cond = f"{dox}Dox" + (f"_{drug}" if drug else "")
    return cond, dox, drug, batch

def cond_style(cond):
    dox = "plus" if cond.startswith("plus") else "minus"
    drug = ""
    for d in DRUGS:
        if cond.endswith(d):
            drug = d
    color = DRUG_COLOR[drug]
    ls = "-" if dox == "plus" else "--"
    alpha = 1.0 if dox == "plus" else 0.75
    return color, ls, alpha

# ---- region file ---------------------------------------------------------------
def load_regions(path):
    """-> {amplicon: {'tss': int|None, 'regions': [(name, lo, hi, type)]}}"""
    out = {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if not line or (line.startswith("#") and cur is None):
            continue
        if line.startswith(">"):
            cur = line[1:].strip()
            out[cur] = {"tss": None, "tata": None, "conf": None, "init": None,
                        "regions": []}
        elif line.startswith("#"):
            m = re.search(r"tss_col=(\d+)", line)
            if m:
                out[cur]["tss"] = int(m.group(1))
            m = re.search(r"tata_col=(\d+)", line)
            if m:
                out[cur]["tata"] = int(m.group(1))
            m = re.search(r"confidence=(\w+)", line)
            if m:
                out[cur]["conf"] = m.group(1)
            m = re.search(r"init=(\w+)", line)
            if m:
                out[cur]["init"] = m.group(1)
        else:
            # Region lines gained optional RULE / rule-param fields (2026-08-12; see
            # annotate_molecule_regions.parse_rule). This is the AGGREGATE path and always
            # uses center-in-window, so a per-region rule is not applicable here -- drop the
            # extra fields rather than crash on a rule-bearing region file. The per-molecule
            # rule booleans come from annotate_molecule_regions.py.
            parts = [p.strip() for p in line.split(",")]
            name, lo, hi, styp = parts[0], parts[1], parts[2], parts[3]
            out[cur]["regions"].append((name, int(lo), int(hi), styp))
    return out

# ---- segments loading -----------------------------------------------------------
def load_segments(seg_path):
    """-> (segs_by_type dict{type: list[(read,start,end)]}, spans dict{read:(lo,hi)}, n_mol)."""
    seg = pd.read_csv(seg_path, sep="\t")
    spans = (seg.groupby("read_id").agg(lo=("start", "min"), hi=("end", "max")))
    spans = {r: (row.lo, row.hi) for r, row in spans.iterrows()}
    by_type = {}
    for t, sub in seg.groupby("type"):
        by_type[t] = list(zip(sub["read_id"], sub["start"], sub["end"]))
    return by_type, spans, len(spans)

def ov(a, b, lo, hi):
    return a <= hi and b >= lo

# ---- aggregation ----------------------------------------------------------------
def aggregate_amplicon(seg_by_type, spans, n_mol, regions):
    """Return list of per-region stat dicts for one (sample, amplicon)."""
    rows = []
    for name, lo, hi, styp in regions:
        segs = seg_by_type.get(styp, [])
        # molecules (read ids) hitting the region, by overlap and by center
        ov_reads, ctr_reads, widths = set(), set(), []
        for read, a, b in segs:
            if ov(a, b, lo, hi):
                ov_reads.add(read)
                widths.append(b - a)
                c = (a + b) / 2.0
                if lo <= c <= hi:
                    ctr_reads.add(read)
        # region-informative molecules: decoded span overlaps region
        n_inf = sum(1 for (slo, shi) in spans.values() if ov(slo, shi, lo, hi))
        w = np.array(widths) if widths else np.array([])
        rows.append({
            "region": name, "segment_type": styp, "lo_col": lo, "hi_col": hi,
            "n_molecules": n_mol, "n_informative": n_inf,
            "n_overlap": len(ov_reads), "n_center": len(ctr_reads),
            "frac_overlap_all": len(ov_reads) / n_mol if n_mol else np.nan,
            "frac_overlap_inf": len(ov_reads) / n_inf if n_inf else np.nan,
            "frac_center_all": len(ctr_reads) / n_mol if n_mol else np.nan,
            "frac_center_inf": len(ctr_reads) / n_inf if n_inf else np.nan,
            "mean_width": float(np.mean(w)) if len(w) else np.nan,
            "median_width": float(np.median(w)) if len(w) else np.nan,
        })
    return rows

def occupancy(seg_by_type, n_mol, styp, col_min, col_max):
    cols = np.arange(col_min, col_max + 1)
    cov = np.zeros(len(cols))
    for read, a, b in seg_by_type.get(styp, []):
        lo = max(a, col_min); hi = min(b, col_max)
        if hi >= lo:
            cov[lo - col_min: hi - col_min + 1] += 1
    return cols, (cov / n_mol if n_mol else cov)

# ---- main -----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--segments_dir", required=True,
                    help="dir with {sample}/ subdirs holding {sample}.{amplicon}.segments.v5.txt")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", nargs="*", default=None,
                    help="sample dir names; default = all L#### subdirs of --segments_dir")
    ap.add_argument("--promoters", nargs="*", default=None,
                    help="restrict to these promoters (default all in region file)")
    ap.add_argument("--sample_map", default=None,
                    help="TSV 'dir<TAB>condition<TAB>batch' mapping sample dirs to "
                         "condition labels (for projects whose dir names don't parse "
                         "as {minusDox|plusDox}[_DRUG]_b#, e.g. P116). Overrides parse_sample.")
    ap.add_argument("--cross_conditions", nargs="*", default=None,
                    help="condition labels to show in the cross-promoter figures "
                         "(default: minusDox plusDox plusDox_TRP).")
    ap.add_argument("--cross_reference", default=None,
                    help="reference condition for the cross-promoter heatmap deltas "
                         "(default: first of --cross_conditions, else minusDox).")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "by_promoter"), exist_ok=True)
    regions = load_regions(args.regions)

    # sample dir -> (condition, batch): explicit map (P116) or parsed (P110).
    sample_meta = {}
    if args.sample_map:
        for line in open(args.sample_map):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            sample_meta[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else "b?")

    def resolve(sd):
        if sd in sample_meta:
            cond, batch = sample_meta[sd]
            return cond, None, None, batch
        return parse_sample(sd)

    # cross-promoter condition set + colors + heatmap reference (config-overridable).
    global CROSS_CONDS, CROSS_COND_COLOR, CROSS_REFERENCE
    if args.cross_conditions:
        CROSS_CONDS, CROSS_COND_COLOR = _build_cross(args.cross_conditions)
    CROSS_REFERENCE = (args.cross_reference or
                       (args.cross_conditions[0] if args.cross_conditions else "minusDox"))
    if args.samples:
        samples = args.samples
    else:
        samples = sorted(d for d in os.listdir(args.segments_dir)
                         if re.match(r"^L\d+_", d) and
                         os.path.isdir(os.path.join(args.segments_dir, d)))

    # amplicon -> promoter
    amp_prom = {amp: amp.split("_opJS4_")[0] for amp in regions}
    if args.promoters:
        keep = set(args.promoters)
        amp_prom = {a: p for a, p in amp_prom.items() if p in keep}

    tidy = []
    # occupancy cache for plotting: {(amp, cond, batch, type): (cols, occ, n)}
    occ_cache = {}
    plot_ctx = {}  # prom -> dict(amp, reg, regs, col_min, col_max) for cross-promoter figs
    for amp, prom in amp_prom.items():
        reg = regions[amp]
        regs = reg["regions"]
        types_needed = sorted({styp for _, _, _, styp in regs})
        # plotting window from region extents (stays inside the promoter, excludes TetO)
        lo_all = min(r[1] for r in regs); hi_all = max(r[2] for r in regs)
        col_min = max(0, lo_all - 15); col_max = hi_all + 15
        # gather per (cond,batch) segment data, plus pooled 'combined'
        pooled = {}  # (cond) -> merged seg_by_type/spans across batches
        for sd in samples:
            seg_path = os.path.join(args.segments_dir, sd,
                                    f"{sd}.{amp}.segments.v5.txt")
            if not os.path.exists(seg_path) or os.path.getsize(seg_path) == 0:
                continue
            cond, dox, drug, batch = resolve(sd)
            try:
                by_type, spans, n_mol = load_segments(seg_path)
            except Exception as e:
                print(f"  WARN {sd} {amp}: {e}", file=sys.stderr)
                continue
            # per-batch rows
            for row in aggregate_amplicon(by_type, spans, n_mol, regs):
                tidy.append({"promoter": prom, "amplicon": amp, "condition": cond,
                             "batch": batch, "sample": sd, **row})
            for t in types_needed:
                occ_cache[(amp, cond, batch, t)] = occupancy(by_type, n_mol, t, col_min, col_max)
            # accumulate for combined (prefix read ids by sample to keep unique)
            pc = pooled.setdefault(cond, {"by_type": {}, "spans": {}})
            for t, lst in by_type.items():
                pc["by_type"].setdefault(t, []).extend(
                    [(f"{sd}::{r}", a, b) for r, a, b in lst])
            for r, sp in spans.items():
                pc["spans"][f"{sd}::{r}"] = sp
        # combined rows + occupancy
        for cond, pc in pooled.items():
            n_mol = len(pc["spans"])
            for row in aggregate_amplicon(pc["by_type"], pc["spans"], n_mol, regs):
                tidy.append({"promoter": prom, "amplicon": amp, "condition": cond,
                             "batch": "combined", "sample": "combined", **row})
            for t in types_needed:
                occ_cache[(amp, cond, "combined", t)] = occupancy(pc["by_type"], n_mol, t, col_min, col_max)
        plot_ctx[prom] = dict(amp=amp, reg=reg, regs=regs,
                              col_min=col_min, col_max=col_max)
        # ---- per-promoter figure ----
        make_promoter_pdf(prom, amp, reg, regs, occ_cache, tidy,
                          col_min, col_max, os.path.join(args.outdir, "by_promoter", f"{prom}.pdf"))
        print(f"  {prom} ({amp}) done")

    df = pd.DataFrame(tidy)
    csv = os.path.join(args.outdir, "footprint_regions.tidy.csv")
    df.to_csv(csv, index=False)
    print(f"\nWrote {csv}  ({len(df)} rows)")

    # ---- cross-promoter summary figures ----
    make_cross_abundance_grid(df, plot_ctx,
                              os.path.join(args.outdir, "cross_promoter_abundance.pdf"))
    make_cross_profile_grid(occ_cache, plot_ctx,
                            os.path.join(args.outdir, "cross_promoter_profiles.pdf"))
    make_cross_heatmap(df, plot_ctx,
                       os.path.join(args.outdir, "cross_promoter_heatmap.pdf"))

def make_promoter_pdf(prom, amp, reg, regs, occ_cache, tidy, col_min, col_max, path):
    tss = reg["tss"]
    def d2t(col):
        return (tss - col) if tss is not None else col
    conds = [c for c in ["minusDox", "plusDox", "minusDox_A485", "plusDox_A485",
                          "minusDox_BRM", "plusDox_BRM", "minusDox_TRP", "plusDox_TRP"]
             if (amp, c, "combined", regs[0][3]) in occ_cache]
    types_needed = sorted({r[3] for r in regs})
    with PdfPages(path) as pdf:
        # occupancy profiles (one panel per segment type present: UNID, NUC, ...)
        for t in types_needed:
            fig, ax = plt.subplots(figsize=(11, 5.2))
            nmax = 0
            for c in conds:
                key = (amp, c, "combined", t)
                if key not in occ_cache:
                    continue
                cols, occ = occ_cache[key]
                col, ls, al = cond_style(c)
                xs = d2t(cols); o = np.argsort(xs)
                ax.plot(xs[o], occ[o], color=col, ls=ls, alpha=al, lw=1.6, label=c)
            for name, lo, hi, styp in regs:
                if styp != t:
                    continue
                ax.axvspan(d2t(hi), d2t(lo), color="grey", alpha=0.12)
                ax.text((d2t(lo) + d2t(hi)) / 2, ax.get_ylim()[1] * 0.97, name,
                        ha="center", va="top", fontsize=7)
            if tss is not None:
                ax.axvline(0, color="k", ls=":", lw=0.8)
                if reg["tata"] is not None:
                    ax.axvline(d2t(reg["tata"]), color="green", ls=":", lw=0.7)
            ax.set_xlabel("distance to TSS (bp)   [upstream/TetO <-- 0 --> gene body]"
                          if tss is not None else "matrix column")
            ax.set_ylabel(f"frac molecules with {t} footprint")
            conf = "?"
            ax.set_title(f"{prom} {amp}: {t} occupancy (combined b1+b2)")
            ax.legend(fontsize=7, ncol=2)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # region abundance bars (center-based, informative denom) + width
        sub = pd.DataFrame([r for r in tidy if r["promoter"] == prom and r["batch"] == "combined"])
        for metric, ylab in [("frac_center_inf", "abundance (center, region-informative denom)"),
                             ("mean_width", "mean footprint width (bp)")]:
            rnames = [r[0] for r in regs]
            fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(rnames)), 5))
            x = np.arange(len(rnames)); w = 0.8 / max(1, len(conds))
            for i, c in enumerate(conds):
                vals = []
                for rn in rnames:
                    m = sub[(sub.region == rn) & (sub.condition == c)]
                    vals.append(m[metric].values[0] if len(m) else np.nan)
                col, ls, al = cond_style(c)
                ax.bar(x + (i - (len(conds) - 1) / 2) * w, vals, w, color=col, alpha=al,
                       edgecolor="k" if ls == "-" else "none", linewidth=0.4, label=c)
            ax.set_xticks(x); ax.set_xticklabels(rnames, fontsize=8, rotation=20)
            ax.set_ylabel(ylab); ax.set_title(f"{prom} {amp}: {ylab}")
            ax.legend(fontsize=7, ncol=2)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

# ---- cross-promoter summary figures ---------------------------------------------
# The three conditions requested for the cross-promoter story, and the four
# regions.  Abundance is always the center-in-window / region-informative metric
# (frac_center_inf): a molecule counts only if its footprint CENTER lands in the
# region, so a wide adjacent footprint clipping the region edge does not inflate
# the count (the boundary-overlap confound).
CROSS_CONDS = [("minusDox", "−Dox"), ("plusDox", "+Dox"),
               ("plusDox_TRP", "+Dox+TRP")]
CROSS_COND_COLOR = {"minusDox": "#bdbdbd", "plusDox": "#3182bd",
                    "plusDox_TRP": "#de2d26"}
CROSS_REFERENCE = "minusDox"   # heatmap deltas are each condition minus this

# Known P110 display names/colors; anything else gets a tab10 colour + raw label.
_CROSS_KNOWN_COLOR = {"minusDox": "#bdbdbd", "plusDox": "#3182bd",
                      "plusDox_TRP": "#de2d26"}
_CROSS_KNOWN_DISP = {"minusDox": "−Dox", "plusDox": "+Dox",
                     "plusDox_TRP": "+Dox+TRP"}

def _cross_disp(cond):
    return _CROSS_KNOWN_DISP.get(cond, cond.replace("_plusDox", "").replace(
        "_minusDox", " −dox").replace("plusDox", "+Dox").replace("minusDox", "−Dox"))

def _build_cross(conds):
    """(CROSS_CONDS list, CROSS_COND_COLOR dict) for an arbitrary condition set."""
    cmap = plt.cm.tab10.colors
    extra = [c for c in conds if c not in _CROSS_KNOWN_COLOR]
    cross_conds = [(c, _cross_disp(c)) for c in conds]
    colors = {c: (_CROSS_KNOWN_COLOR.get(c) or cmap[extra.index(c) % len(cmap)])
              for c in conds}
    return cross_conds, colors
CROSS_REGIONS = ["TATA", "core_PIC", "TSS_Inr", "pause"]
REGION_LABEL = {"TATA": "TATA", "core_PIC": "PIC", "TSS_Inr": "TSS/Inr",
                "pause": "Pause"}
CROSS_METRIC = "frac_center_inf"

def _grid_shape(n, ncol=4):
    return ncol, int(np.ceil(n / ncol))

def _is_low_conf(reg):
    return (reg.get("conf") or "").lower() != "high"

def make_cross_abundance_grid(df, plot_ctx, path):
    """One panel per promoter: grouped bars (region x 3 conditions), abundance."""
    from matplotlib.patches import Patch
    proms = list(plot_ctx)
    ncol, nrow = _grid_shape(len(proms))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.7 * nrow),
                             squeeze=False)
    for idx, prom in enumerate(proms):
        ax = axes[idx // ncol][idx % ncol]
        sub = df[(df.promoter == prom) & (df.batch == "combined")]
        present = [r for r in CROSS_REGIONS if r in set(sub.region)]
        nc = len(CROSS_CONDS)
        x = np.arange(len(present)); w = 0.8 / max(1, nc)
        for i, (cond, _) in enumerate(CROSS_CONDS):
            vals = []
            for rn in present:
                m = sub[(sub.region == rn) & (sub.condition == cond)]
                vals.append(m[CROSS_METRIC].values[0] if len(m) else np.nan)
            ax.bar(x + (i - (nc - 1) / 2) * w, vals, w,
                   color=CROSS_COND_COLOR[cond], edgecolor="k", linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels([REGION_LABEL.get(r, r) for r in present],
                           fontsize=7, rotation=20)
        ax.tick_params(axis="y", labelsize=7)
        low = _is_low_conf(plot_ctx[prom]["reg"])
        ax.set_title(prom + (" *" if low else ""), fontsize=8,
                     color="0.45" if low else "k",
                     style="italic" if low else "normal")
        if idx % ncol == 0:
            ax.set_ylabel("frac (center)", fontsize=7)
    for j in range(len(proms), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    handles = [Patch(facecolor=CROSS_COND_COLOR[c], edgecolor="k", label=l)
               for c, l in CROSS_CONDS]
    fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Region footprint abundance across promoters  "
                 "(center-in-window; * = low-confidence TSS)",
                 y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(path); plt.close(fig)
    print(f"Wrote {path}")

def make_cross_profile_grid(occ_cache, plot_ctx, path, styp="UNID"):
    """One panel per promoter: UNID occupancy profile, 3 conditions overlaid."""
    from matplotlib.lines import Line2D
    proms = list(plot_ctx)
    ncol, nrow = _grid_shape(len(proms))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.7 * nrow),
                             squeeze=False)
    for idx, prom in enumerate(proms):
        ax = axes[idx // ncol][idx % ncol]
        ctx = plot_ctx[prom]; amp = ctx["amp"]; reg = ctx["reg"]
        tss = reg["tss"]
        d2t = (lambda col: (tss - col) if tss is not None else col)
        for cond, _ in CROSS_CONDS:
            key = (amp, cond, "combined", styp)
            if key not in occ_cache:
                continue
            cols, occ = occ_cache[key]
            xs = d2t(cols); o = np.argsort(xs)
            ax.plot(xs[o], occ[o], color=CROSS_COND_COLOR[cond], lw=1.3)
        for name, lo, hi, t in ctx["regs"]:
            if t != styp or name not in CROSS_REGIONS:
                continue
            ax.axvspan(d2t(hi), d2t(lo), color="grey", alpha=0.10)
        if tss is not None:
            ax.axvline(0, color="k", ls=":", lw=0.7)
            if reg["tata"] is not None:
                ax.axvline(d2t(reg["tata"]), color="green", ls=":", lw=0.6)
        ax.tick_params(labelsize=7)
        low = _is_low_conf(reg)
        ax.set_title(prom + (" *" if low else ""), fontsize=8,
                     color="0.45" if low else "k",
                     style="italic" if low else "normal")
        if idx % ncol == 0:
            ax.set_ylabel(f"frac {styp}", fontsize=7)
        if idx // ncol == nrow - 1:
            ax.set_xlabel("dist to TSS (bp)", fontsize=7)
    for j in range(len(proms), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    handles = [Line2D([0], [0], color=CROSS_COND_COLOR[c], lw=1.6, label=l)
               for c, l in CROSS_CONDS]
    fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(f"{styp} occupancy profiles across promoters  "
                 "(upstream/TetO <-- 0 --> gene body; * = low-conf TSS)",
                 y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(path); plt.close(fig)
    print(f"Wrote {path}")

def make_cross_heatmap(df, plot_ctx, path):
    """Promoter x region abundance heatmaps: one per condition + deltas.

    Rows = promoters (region-file order, low-conf flagged), cols = the four
    regions.  Panels: one sequential heatmap per condition (shared vmax) then
    two diverging delta panels (+Dox - -Dox  and  +Dox+TRP - +Dox).  NaN cells
    (region absent, e.g. TATA in dispersed promoters) are drawn grey.
    """
    proms = list(plot_ctx)
    regions = CROSS_REGIONS
    reglab = [REGION_LABEL.get(r, r) for r in regions]

    # group rows: focused (high-conf TSS) on top, dispersed/low-conf below
    focused = [p for p in proms if not _is_low_conf(plot_ctx[p]["reg"])]
    dispersed = [p for p in proms if _is_low_conf(plot_ctx[p]["reg"])]
    row_proms = focused + dispersed
    split = len(focused)  # divider drawn between rows split-1 and split

    def mat(cond):
        M = np.full((len(row_proms), len(regions)), np.nan)
        for i, prom in enumerate(row_proms):
            sub = df[(df.promoter == prom) & (df.batch == "combined")
                     & (df.condition == cond)]
            for j, rn in enumerate(regions):
                m = sub[sub.region == rn]
                if len(m):
                    M[i, j] = m[CROSS_METRIC].values[0]
        return M

    mats = {c: mat(c) for c, _ in CROSS_CONDS}
    vmax = np.nanmax([np.nanmax(m) for m in mats.values()])
    # absolute panel per condition, then a diverging delta panel for every
    # non-reference condition (condition − reference).
    ref = CROSS_REFERENCE if CROSS_REFERENCE in mats else CROSS_CONDS[0][0]
    ref_disp = dict(CROSS_CONDS).get(ref, ref)
    deltas = [(c, disp, mats[c] - mats[ref]) for c, disp in CROSS_CONDS if c != ref]
    dvals = np.concatenate([d.ravel() for _, _, d in deltas]) if deltas else np.array([0.0])
    dmax = np.nanmax(np.abs(dvals)) or 1.0

    panels = [(disp, mats[c], "viridis", 0, vmax) for c, disp in CROSS_CONDS]
    panels += [(f"{disp} − {ref_disp}", d, "RdBu_r", -dmax, dmax)
               for _, disp, d in deltas]

    ylab = [p + (" *" if _is_low_conf(plot_ctx[p]["reg"]) else "")
            for p in row_proms]
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(2.6 * len(panels), 0.34 * len(row_proms) + 2),
                             squeeze=False)
    for k, (title, M, cmap_name, vmin, vmx) in enumerate(panels):
        ax = axes[0][k]
        cmap = plt.get_cmap(cmap_name).copy(); cmap.set_bad("0.85")
        im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmx)
        ax.set_xticks(range(len(regions)))
        ax.set_xticklabels(reglab, fontsize=7, rotation=30, ha="right")
        ax.set_title(title, fontsize=9)
        if 0 < split < len(row_proms):
            ax.axhline(split - 0.5, color="k", lw=1.5)
        if k == 0:
            ax.set_yticks(range(len(row_proms)))
            ax.set_yticklabels(ylab, fontsize=6.5)
        else:
            ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).ax.tick_params(labelsize=6)
    fig.suptitle("Region footprint abundance (center-in-window)  "
                 "top block = focused TSS, bottom = dispersed/low-conf (*)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path); plt.close(fig)
    print(f"Wrote {path}")

if __name__ == "__main__":
    main()
