#!/usr/bin/env python3
"""
Python3 drop-in replacement for Georgi's original protection/cluster
script. Uses array masking logic for C contexts instead of dicts.

Outputs:
  - <prefix>.<label>.full_unclustered.matrix
  - <prefix>.<label>.dedup.full_unclustered.matrix
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
                   help="Generate heatmap PNG, requires -cluster")
    p.add_argument("-cluster", action="store_true",
                   help="Cluster rows with sklearn k-means on Cs")
    p.add_argument(
        "--dedup_on",
        choices=["c_type", "allC"],
        default="c_type",
        help=(
            "Which positions to consider when deduplicating rows: "
            "'c_type' = only Cs in the requested context (default); "
            "'allC' = all cytosines in the region."
        ),
    )

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
                    ref_idx = (min_r - region_start) + j
                    ref_base = ref[ref_idx] if 0 <= ref_idx < len(ref) else None  

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
        1 = converted cytosine
        0 = unconverted cytosine
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

    # allCs mode: every cytosine in the region
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

    # Actually score the read bases at those positions
    for i in range(ovl_start, ovl_end):
        if not score_pos[i]:
            continue
        rr = (region_start + i) - rstart    # read-relative index of this base
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
    - Returns (labels, order)
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
        return 0, 0, np.array([], dtype=np.int64), arr

    total = arr.shape[0]
    row_view = arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))
    _, kept_idx = np.unique(row_view, return_index=True)
    kept_idx.sort()
    dedup = arr[kept_idx]
    unique_states = dedup.shape[0]

    if do_dedup:
        return total, unique_states, kept_idx, dedup
    
    return total, unique_states, np.arange(total, dtype=np.int64), arr

def subset_and_cluster(
    all_scores: np.ndarray,
    ids: List[str],
    cover_counts: List[int],
    subset: Optional[int],
    n_clusters: int = 4,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Apply optional coverage-based subsetting and k-means clustering.

    Returns:
      clust_scores : (n_rows_subset, n_cols) int array
      labels       : (n_rows_subset,) int cluster labels (1..k)
      ids_out      : list of IDs in clustered order
    """
    if all_scores.size == 0 or len(ids) == 0:
        return all_scores, np.ones(0, dtype=int), ids

    scores = all_scores
    ids_out = list(ids)
    cov = np.asarray(cover_counts, dtype=float)

    # Optional subsetting by coverage
    if subset is not None and scores.shape[0] > subset:
        order = np.argsort(-cov)  # descending
        cutoff_cov = cov[order[subset - 1]]
        top_mask = cov >= cutoff_cov
        top_idx = np.where(top_mask)[0]
        if top_idx.size > subset:
            strict_idx = np.where(cov > cutoff_cov)[0]
            need = subset - strict_idx.size
            tied_idx = np.where(cov == cutoff_cov)[0]
            pick = np.random.choice(tied_idx, size=need, replace=False)
            final_idx = np.concatenate([strict_idx, pick])
        else:
            final_idx = top_idx
        final_idx.sort()
        scores = scores[final_idx]
        ids_out = [ids_out[i] for i in final_idx]

    # Cluster
    if scores.shape[0] <= 1:
        labels = np.ones(scores.shape[0], dtype=int)
        return scores, labels, ids_out

    labels, order = cluster_rows_kmeans(scores, n_clusters=n_clusters)
    scores = scores[order]
    labels = labels[order]
    ids_out = [ids_out[i] for i in order]

    return scores, labels, ids_out


def make_heatmap(all_scores: np.ndarray, out_png: str):
    """
    Draw a binary heatmap (0/1 with -1 as transparent) roughly matching
    Georgi's heatmap.py layout logic.
    """
    if all_scores.size == 0:
        return

    xps = 10.0      # x pixel size
    yps = 3.0       # y pixel size
    inches = 10.0   # figure width in inches

    plot_arr = all_scores.astype(float)
    plot_arr[plot_arr < 0] = np.nan
    plot_masked = np.ma.masked_invalid(plot_arr)

    n_rows, n_cols = plot_masked.shape
    if n_rows == 0 or n_cols == 0:
        return

    height_px = n_rows * yps
    width_px = n_cols * xps
    height_inches = inches * (height_px / width_px)

    cmap = plt.get_cmap("binary").copy()
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))  # NaNs transparent

    fig, ax = plt.subplots(figsize=(inches, height_inches), dpi=600)
    ax.imshow(
        plot_masked,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_axis_off()
    fig.patch.set_alpha(0.0)

    plt.savefig(out_png, dpi=600, transparent=True,
                bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    fa = pysam.FastaFile(args.fa)

    with open(args.out_qc_file, "w") as out_qc, open(args.peaks) as pl:
        out_qc.write("amplicon\ttotal_reads\tobserved_states\treads_per_state\n")

        for raw in pl:
            # Skip empty lines and comments
            if (not raw.strip()) or raw.startswith("#"):
                continue

            fields = raw.rstrip("\n").split("\t")
            try:
                chrom = fields[args.chr_col]
                start = max(0, int(fields[args.start_col]))
                end = int(fields[args.end_col])
                strand = "+" if args.unstranded else fields[args.strand_col]
            except Exception as e:
                print(f"Skipping bad line: {raw.strip()} ({e})", file=sys.stderr)
                continue

            if args.label_col is not None and args.label_col < len(fields):
                label = fields[args.label_col]
            else:
                label = f"{chrom}_{start}-{end}_{'for' if strand == '+' else 'rev'}"

            if (end - start) <= 0:
                print(f"Skipping empty region: {label}", file=sys.stderr)
                continue

            # Reference
            ref = fa.fetch(chrom, start, end).upper()

            # Collect alignments overlapping region
            rows = []
            for aln in bam.fetch(chrom, start, end):
                if aln.is_unmapped:
                    continue
                r0, rseq = reconstruct_aligned_query(aln)
                if not rseq:
                    continue
                rows.append({"qname": aln.query_name, "rstart": r0, "seq": rseq})

            if not rows:
                print(f"No reads for: {label}")
                continue

            df = pd.DataFrame(rows)

            # Group and score per read ID
            scores_ctx_list: List[np.ndarray] = []      # requested context (CG/GC/both/allC)
            allC_scores_list: List[np.ndarray] = []     # all Cs for deduplication (only used if dedup_on == 'allC')
            cover_counts: List[int] = []
            ids: List[str] = []

            for qname, g in df.groupby("qname", sort=False):
                segs = list(zip(g["rstart"].to_numpy(), g["seq"].to_numpy()))
                rstart_m, seq_m, cov_m = merge_alignments(
                    segs,
                    region_start=start,
                    ref=ref,
                    prefer_ref_on_conflict=True,
                )
                if not seq_m:
                    continue

                # First: score in the requested context (this is what we'll output)
                sc_ctx, cov = score_read_against_region(
                    rstart_m,
                    seq_m,
                    start,
                    args.c_type,
                    args.noEndogenousMethylation,
                    ref,
                )
                if sc_ctx is None:
                    continue
                if (args.minCov is not None) and (cov.mean() < args.minCov):
                    continue

                scores_ctx_list.append(sc_ctx)
                cover_counts.append(int(cov.sum()))
                ids.append(qname)

                # Optionally: score all cytosines for deduplication basis
                if args.dedup_on == "allC":
                    if args.c_type == "allC":
                        # already all Cs
                        sc_allC = sc_ctx
                    else:
                        sc_allC, _ = score_read_against_region(
                            rstart_m,
                            seq_m,
                            start,
                            "allC",   # score every C
                            False,    # no_endog irrelevant for allC
                            ref,
                        )
                    allC_scores_list.append(sc_allC)

            if not scores_ctx_list:
                print(f"No reads retained for: {label}")
                continue

            # 1) Build matrix for requested context
            scores_ctx_arr = np.vstack(scores_ctx_list).astype(int, copy=False)

            # Decide which matrix to use for deduplication
            if args.dedup_on == "allC":
                if args.c_type == "allC":
                    dedup_basis_arr = scores_ctx_arr
                else:
                    dedup_basis_arr = np.vstack(allC_scores_list).astype(int, copy=False)
            else:
                # dedup only on the scored context (c_type)
                dedup_basis_arr = scores_ctx_arr

            # Deduplicate patterns based on the chosen basis
            total, uniq, kept_idx, _ = qc_dedup(dedup_basis_arr, do_dedup=True)

            reads_per_state = (total / uniq) if uniq > 0 else float("nan")
            out_qc.write(f"{label}\t{total}\t{uniq}\t{reads_per_state:.2f}\n")
            out_qc.flush()

            # Keep IDs and coverage in sync with deduped matrix
            ids_dedup = [ids[i] for i in kept_idx]
            cov_dedup = [cover_counts[i] for i in kept_idx]

            # Context-specific deduped matrix (same columns as requested context)
            scores_ctx_dedup = scores_ctx_arr[kept_idx]

            # 2) Write non-dedup (all reads) matrix for requested context
            full_path = f"{args.out_prefix}.{label}.full_unclustered.matrix"
            write_matrix(full_path, chrom, start, end, scores_ctx_arr, ids, strand)
            print(f"Wrote full matrix (all reads): {full_path}")

            # 3) Write dedup matrix for requested context (fixed suffix)
            dedup_path = f"{args.out_prefix}.{label}.dedup.full_unclustered.matrix"
            write_matrix(dedup_path, chrom, start, end, scores_ctx_dedup, ids_dedup, strand)
            print(f"Wrote dedup matrix (basis={args.dedup_on}): {dedup_path}")

            # 4) Subset + cluster on the deduped matrix if requested
            if args.cluster and scores_ctx_dedup.shape[0] > 1:
                clust_scores, labels, ids_clust = subset_and_cluster(
                    scores_ctx_dedup,
                    ids_dedup,
                    cov_dedup,
                    args.subset,
                    n_clusters=4,
                )

                clust_path = f"{args.out_prefix}.{label}.clustered.matrix"
                header = "#" + chrom + "".join(f"\t{i}" for i in range(start, end))
                with open(clust_path, "w") as fh:
                    fh.write(header + "\n")
                    mat = clust_scores[:, ::-1] if strand == "-" else clust_scores
                    for lab, row in zip(labels, mat):
                        fh.write(f"{lab}\t" + "\t".join(map(str, row.tolist())) + "\n")
                print(f"Wrote clustered matrix (dedup basis={args.dedup_on}): {clust_path}")

                # 5) Heatmap on clustered matrix if requested
                if args.heatmap:
                    out_png = f"{args.out_prefix}.{label}.matrix.png"
                    make_heatmap(clust_scores, out_png)

    bam.close()
    fa.close()


if __name__ == "__main__":
    main()