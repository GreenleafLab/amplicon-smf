#!/usr/bin/env python3
"""
Python3 drop-in replacement for Georgi's original protection/cluster
script. Uses array masking logic for C contexts instead of dicts.

Outputs:
  - <prefix>.<label>.full_unclustered.matrix
  - <prefix>.<label>.clustered.matrix (when -cluster is used)
  - QC summary to out_done_file
"""
import argparse
import sys
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def parse_args():
    p = argparse.ArgumentParser(
        description=("Compute SMF cytosine protection matrices from a BAM.")
    )
    # Keep positional args order/meaning identical to the legacy script
    p.add_argument("bam",           help="BAM file")
    p.add_argument("fa",            help="Reference FASTA for amplicons")
    p.add_argument("c_type",        choices=["CG", "GC", "both_dimers", "allC"], help="C context(s) to score")
    p.add_argument("peaks",         help="BED-like regions file")
    p.add_argument("chr_col",       type=int)
    p.add_argument("start_col",     type=int)
    p.add_argument("end_col",       type=int)
    p.add_argument("strand_col",    type=int)
    p.add_argument("out_prefix",    help="Prefix for matrix outputs")
    p.add_argument("out_qc_file",   help="Per-region QC summary output")

    # Legacy flags (same spelling/behavior)
    p.add_argument("-subset", type=int, default=None,
                   help="Keep only N best-covered fragments per region")
    p.add_argument("-label", dest="label_col", type=int, default=None,
                   help="Column index to label regions (instead of coordinates)")
    p.add_argument("-noEndogenousMethylation", action="store_true",
                   help="In GC mode, do NOT exclude GCG sites")
    p.add_argument("-minCov", type=float, default=None,
                   help="Minimum fractional region coverage required for a read")
    p.add_argument("-unstranded", action="store_true",
                   help="Ignore strand; treat all regions as '+'")
    p.add_argument("-heatmap", action="store_true",
                   help="Generate heatmap PNG")
    p.add_argument("-cluster", action="store_true",
                   help="Cluster rows with Ward linkage (distance threshold=0)")
    p.add_argument("-dedup", action="store_true", default=False,
                   help="Deduplicate rows with identical C protection patterns before clustering and plotting")
    return p.parse_args()

def reconstruct_aligned_query(aln: pysam.AlignedSegment) -> Tuple[int, str]:
    """
    Reconstruct the query sequence aligned to the reference:
      M/=/X  -> copy bases
      I      -> advance query (skip on reference)
      D      -> pad with 'N'
      S      -> skip soft-clipped
    Returns (ref_start, aligned_query_seq).
    """
    q = aln.query_sequence or ""
    out = []
    qpos = 0
    rstart = aln.reference_start  # 0-based
    for op, bp in (aln.cigartuples or []):
        if op == 4:            # S
            qpos += bp
        elif op in (0, 7, 8):  # M, =, X
            out.append(q[qpos:qpos+bp])
            qpos += bp
        elif op == 1:          # I
            qpos += bp
        elif op == 2:          # D
            out.append("N" * bp)
        # other ops (H,N,P) are ignored for these amplicons
    return rstart, "".join(out)

def merge_alignments(
    segs: List[Tuple[int, str]],
    region_start: int,
    ref: str,
    prefer_ref_on_conflict: bool = True,
) -> Tuple[int, str, np.ndarray]:
    """
    Merge multiple aligned query segments (same orientation) into one.
    segs: list of (rstart, seq) where 'seq' is already aligned to reference
          (your reconstruct_aligned_query output; may include 'N' for deletions).
    region_start: genomic start of the region we're scoring (0-based).
    ref: reference substring for the whole region (len = region_len).
    prefer_ref_on_conflict: if two non-'N' bases disagree, choose the one that equals ref,
                            else write 'N'.

    Returns:
      rstart_merged: genomic start of merged sequence
      merged_seq: merged sequence string
      cov: int coverage vector (per base in merged_seq)
    """
    if not segs:
        return region_start, "", np.zeros(0, dtype=int)

    # sort by genomic start
    segs = sorted(segs, key=lambda x: x[0])
    min_r = segs[0][0]
    max_end = max(r + len(s) for r, s in segs)

    # workspace spanning [min_r, max_end)
    L = max_end - min_r
    buf = np.full(L, "N", dtype="<U1")
    cov = np.zeros(L, dtype=int)

    # place segments
    for r0, s in segs:
        off = r0 - min_r
        for i, b in enumerate(s):
            if b == "":  # just in case
                continue
            j = off + i
            if j < 0 or j >= L:
                continue
            if b == "N":
                continue
            if buf[j] == "N":
                buf[j] = b
            elif buf[j] != b:
                # conflict between two informative bases
                if prefer_ref_on_conflict:
                    # choose base matching reference if possible
                    ref_base = ref[(min_r - region_start) + j] if 0 <= (min_r - region_start) + j < len(ref) else None
                    if ref_base == b:
                        buf[j] = b
                    elif ref_base == buf[j]:
                        pass  # keep existing
                    else:
                        buf[j] = "N"
                else:
                    buf[j] = "N"
            cov[j] += 1

    # trim leading/trailing Ns to produce a compact merged chunk
    nonN = np.where(buf != "N")[0]
    if nonN.size == 0:
        return region_start, "", np.zeros(0, dtype=int)

    left = int(nonN.min())
    right = int(nonN.max()) + 1
    merged_seq = "".join(buf[left:right].tolist())
    merged_rstart = min_r + left
    merged_cov = cov[left:right]

    return merged_rstart, merged_seq, merged_cov

def score_read_against_region(
    rstart: int,
    rseq: str,
    region_start: int,
    c_type: str,
    no_endog: bool,
    ref: str,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Build scores for a single read over [region_start, region_start+len(ref)).

    scores: int8 array of length len(ref) where:
        1 = protected / unmethylated
        0 = methylated / unprotected
       -1 = not a scored site (non-target context or uncovered)

    Contexts:

      - GC + no_endog=True  → GCN  (include GCG)
      - GC + no_endog=False → GCH  (exclude GCG)
      - CG                  → NCG (C of a CG)
      - allC                → all cytosines
      - both_dimers         → union of GC(+/-GCG filter) and CG

    The *C* is the base scored in the motif.
    """
    # Length of reference region
    L = len(ref)
    if L == 0:
        return None, None

    scores = np.full(L, -1, dtype=np.int8)
    covered = np.zeros(L, dtype=bool)

    # Overlap of this read's aligned segment with the region
    ovl_start = max(0, rstart - region_start)
    ovl_end   = min(L, (rstart - region_start) + len(rseq))
    if ovl_end <= ovl_start:
        return None, None

    covered[ovl_start:ovl_end] = True

    # reference as bytes
    ref_arr = np.frombuffer(ref.encode("ascii"), dtype="S1")

    # Decide which positions to score
    score_pos = np.zeros(L, dtype=bool)

    # allCs mode: just every cytosine in the region
    if c_type == "allC":
        score_pos |= (ref_arr == b"C")

    # CG / both_dimers: C of NCG (the C in "CG")
    if c_type in ("CG", "both_dimers"):
        for i in range(0, L - 1):
            if ref_arr[i] == b"C" and ref_arr[i + 1] == b"G":
                score_pos[i] = True  # score the C

    # GC / both_dimers: C of GCN or GCH depending on no_endog
    if c_type in ("GC", "both_dimers"):
        for i in range(0, L - 1):
            if ref_arr[i] == b"G" and ref_arr[i + 1] == b"C":
                include = True
                if not no_endog:
                    # exclude GCG when we care about endogenous methylation
                    if i + 2 < L and ref_arr[i + 2] == b"G":
                        include = False
                if include:
                    c_pos = i + 1
                    score_pos[c_pos] = True

    # Now actually score the read bases at those positions
    for i in range(ovl_start, ovl_end):
        if not score_pos[i]:
            continue
        # read-relative index of this genomic base
        rr = (region_start + i) - rstart
        if not (0 <= rr < len(rseq)):
            continue
        base = rseq[rr]
        # C = unconverted (protected), anything else = converted
        scores[i] = 0 if base == "C" else 1

    return scores, covered

def cluster_rows_kmeans(X: np.ndarray, n_clusters: int = 4):
    """
    K-means on ternary {-1,0,1} matrices.

    - Enforces column rule: keep only columns with no -1 anywhere.
      (Columns with any -1 → dropped; columns with all -1 → dropped.)
    - Returns (labels, order):
        labels : 1..k cluster labels (np.int32)
        order  : row permutation, sorted by (cluster id, distance to its center)

    If nothing informative remains (no kept columns) or n<=1,
    returns labels=all ones and identity order.
    """
    n = X.shape[0]
    if n <= 1 or X.size == 0:
        return np.ones(n, dtype=np.int32), np.arange(n, dtype=np.int64)

    # Keep only columns with no -1 anywhere (fully observed)
    keep = (X != -1).all(axis=0)
    Xf = X[:, keep]
    if Xf.shape[1] == 0:
        return np.ones(n, dtype=np.int32), np.arange(n, dtype=np.int64)

    # Float for k-means
    Xf = Xf.astype(float, copy=False)

    # Clamp k to number of rows
    k = max(1, min(n_clusters, n))

    km = KMeans(n_clusters=k, init="k-means++", n_init="auto",
                max_iter=300, random_state=0)
    km.fit(Xf)
    labels0 = km.labels_                  # 0..k-1
    # distance of each point to its assigned center
    d = km.transform(Xf)[np.arange(n), labels0]

    # Order: by cluster id, then by increasing distance
    order = np.lexsort((d, labels0)).astype(int)

    # 1-based labels
    labels = (labels0 + 1).astype(int)
    return labels, order

def write_matrix(path: str, chrom: str, start: int, end: int,
                 X: np.ndarray, ids: List[str], strand: str):
    """Write a tab-separated matrix"""
    header = "#" + chrom + "".join(f"\t{i}" for i in range(start, end))
    with open(path, "w") as fh:
        fh.write(header + "\n")
        mat = X[:, ::-1] if strand == "-" else X
        for rid, row in zip(ids, mat):
            fh.write(str(rid) + "\t" + "\t".join(map(str, row.tolist())) + "\n")

def qc_dedup(arr: np.ndarray, do_dedup: bool=False):
    """
    Returns: total, unique_states, kept_idx, matrix_out
    """
    if arr.ndim != 2 or arr.size == 0:
        return 0, 0, float("nan"), np.array([], dtype=np.int64), arr

    total = arr.shape[0]
    row_view = arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))
    _, kept_idx = np.unique(row_view, return_index=True)
    kept_idx.sort()
    dedup = arr[kept_idx]
    unique_states = dedup.shape[0]

    if do_dedup:
        return total, unique_states, kept_idx, dedup
    
    return total, unique_states, np.arange(total, dtype=np.int64), arr

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    fa = pysam.FastaFile(args.fa)

    with open(args.out_qc_file, "w") as out_qc, open(args.peaks) as pl:
        out_qc.write("amplicon\ttotal_reads\tobserved_states\treads_per_state\n")
        for raw in pl:
            # Skip empty lines and lines starting with '#'
            if (not raw.strip()) or (raw.startswith("#")):
                continue
            
            # Try to unpack each line of FA, skip bad lines
            fields = raw.rstrip("\n").split("\t")
            try:
                chrom = fields[args.chr_col]
                start = max(0, int(fields[args.start_col]))
                end = int(fields[args.end_col])
                strand = "+" if args.unstranded else fields[args.strand_col]
            except Exception as e:
                print(f"Skipping bad line: {raw.strip()} ({e})", file=sys.stderr)
                continue
            
            # Decide how to label
            if args.label_col is not None and args.label_col < len(fields):
                label = fields[args.label_col]
            else:
                label = f"{chrom}_{start}-{end}_{'for' if strand == '+' else 'rev'}"

            # Check for non-empty region
            region_len = end - start
            if region_len <= 0:
                print(f"Skipping empty region: {label}", file=sys.stderr)
                continue

            # Reference & masks
            ref = fa.fetch(chrom, start, end).upper()
            
            # Collect alignments overlapping region and put in DF
            rows = []
            for aln in bam.fetch(chrom, start, end):
                if aln.is_unmapped:
                    continue
                r0, rseq = reconstruct_aligned_query(aln)
                if not rseq:
                    continue
                rows.append({"qname": aln.query_name, "rstart": r0, "seq": rseq})
            
            # Put in DF
            if not rows:
                print(f"No reads for: {label}")
                continue
            df = pd.DataFrame(rows)

            # Group and score
            scores_list = []
            cover_counts = []
            ids = []
            # Group by read id; each group may have 2 forward segments with slight overlap
            for qname, g in df.groupby("qname", sort=False):
                segs = list(zip(g["rstart"].to_numpy(), g["seq"].to_numpy()))
                rstart_m, seq_m, _ = merge_alignments(segs,
                                                      region_start=start,
                                                      ref=ref,
                                                      prefer_ref_on_conflict=True)
                if not seq_m:
                    continue
                sc, cov = score_read_against_region(rstart_m,   # read start
                                                    seq_m,      # merged seq
                                                    start,      # region start      
                                                    args.c_type,
                                                    args.noEndogenousMethylation)
                # TODO: Check CG/GC/GCG ambiguous when different args passed in

                if sc is None:
                    continue
                if (args.minCov is not None) and (cov.mean() < args.minCov):
                    continue

                scores_list.append(sc)
                cover_counts.append(int(cov.sum()))
                ids.append(qname)

            if not scores_list:
                print(f"No reads retained for: {label}")
                continue
            
            # TODO: test this block more thoroughly
            # Subset best-covered reads if requested (stable, tie-aware)
            if args.subset is not None and len(scores_list) > args.subset:
                cov = np.asarray(cover_counts)
                order = np.argsort(-cov)  # descending by coverage
                cutoff_cov = cov[order[args.subset-1]]
                top_mask = cov >= cutoff_cov
                top_idx = np.where(top_mask)[0]
                # if too many due to ties, randomly sample among the tied boundary
                if top_idx.size > args.subset:
                    # keep strictly > cutoff, and sample the rest
                    strict_idx = np.where(cov > cutoff_cov)[0]
                    need = args.subset - strict_idx.size
                    tied_idx = np.where(cov == cutoff_cov)[0]
                    pick = np.random.choice(tied_idx, size=need, replace=False)
                    final_idx = np.concatenate([strict_idx, pick])
                else:
                    final_idx = top_idx
                final_idx.sort()
                scores_list = [scores_list[i] for i in final_idx]
                ids         = [ids[i] for i in final_idx]
            
            # Make array, write some QC metrics, dedup if requested
            all_scores = np.vstack(scores_list).astype(int, copy=False)
            total, uniq, kept_idx, all_scores = qc_dedup(all_scores, do_dedup=args.dedup)
            ids = [ids[i] for i in kept_idx]  # keep IDs in sync
            out_qc.write(f"{label}\t{total}\t{uniq}\t{total/uniq:.2f}\n")
            out_qc.flush()

            # Write full (unclustered) matrix
            full_path = f"{args.out_prefix}.{label}.full_unclustered.matrix"
            write_matrix(full_path, chrom, start, end, all_scores, ids, strand)
            print(f"Wrote full matrix: {full_path}")

            # Optional clustered matrix
            if args.cluster and len(all_scores) > 1:
                labels, order = cluster_rows_kmeans(all_scores, n_clusters=4)
                all_scores = all_scores[order]
                ids        = [ids[i] for i in order]
                labels     = labels[order]

                clust_path = f"{args.out_prefix}.{label}.clustered.matrix"
                header = "#" + chrom + "".join(f"\t{i}" for i in range(start, end))
                with open(clust_path, "w") as fh:
                    fh.write(header + "\n")
                    mat = all_scores[:, ::-1] if strand == "-" else all_scores
                    for lab, rid, row in zip(labels, ids, mat):
                        fh.write(f"{lab}\t" + "\t".join(map(str, row.tolist())) + "\n")
                print(f"wrote clustered matrix: {clust_path}")

                # Optional heatmap
                if args.heatmap:
                    xps = 10.0      # x pixel size
                    yps = 3.0       # y pixel size
                    inches = 10.0   # figure width in inches

                    # Convert to float and mask -1 as np.nan
                    plot_arr = all_scores.astype(float)
                    plot_arr[plot_arr < 0] = np.nan
                    plot_masked = np.ma.masked_invalid(plot_arr)

                    # Check if there's anything to plot
                    n_rows, n_cols = plot_masked.shape
                    if n_rows == 0 or n_cols == 0:
                        continue

                    # Match old script: Height = n_rows*yps, Width = n_cols*xps,
                    # fig size = (inches, inches * Height/Width)
                    height_px = n_rows * yps
                    width_px = n_cols * xps
                    height_inches = inches * (height_px / width_px)

                    # 'binary' colormap with NaNs transparent
                    cmap = plt.get_cmap("binary").copy()
                    cmap.set_bad((1.0, 1.0, 1.0, 0.0))  # RGBA, alpha=0 for NaNs

                    fig, ax = plt.subplots(figsize=(inches, height_inches), dpi=600)
                    ax.imshow(
                        plot_masked,
                        cmap=cmap,
                        vmin=0.0,
                        vmax=1.0,
                        aspect="auto",
                        interpolation="nearest"
                    )
                    ax.set_axis_off()
                    fig.patch.set_alpha(0.0)  # transparent figure background

                    out_png = f"{args.out_prefix}.{label}.matrix.png"
                    plt.savefig(out_png, dpi=600, transparent=True, 
                                bbox_inches="tight", pad_inches=0)
                    plt.close(fig)

    bam.close()
    fa.close()

if __name__ == "__main__":
    main()


