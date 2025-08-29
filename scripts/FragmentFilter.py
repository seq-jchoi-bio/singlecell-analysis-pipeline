# Filtering: Integrated scATAC-seq cell filtering using fragments.tsv.gz only (rice-ready).
# Outputs (under Filtering_results/):
#  - filtered_samples/<sample>_filtered.csv   # Per-sample filtered cells
#  - summary_pre_post.csv                     # Pre/post counts per sample with retention
#  - summary.txt                              # Human-readable summary (no raw logs)
#  - (optional) h5ad/                         # Cached per-sample h5ad if 

import os
import sys
import gzip
import pandas as pd
import numpy as np
import snapatac2 as snap
from datetime import datetime
from typing import Dict, List, Tuple

__version__ = 1.1

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq Cell Filter"
HEADER += "\n      Version {0}".format(__version__)
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

desc_txt = """
This tool aggregates per-sample QC, attaches TSSE scores computed from fragments,
and filters cells with (1) high-confidence barcodes, (2) TSSE >= cutoff,
and (3) UMI within the specified min/max range. Results are saved per sample,
along with a concise summary across all samples.
"""

# User variables (EDIT ME)
base_dir         = "/Users/jchoi/Desktop/Test"   # <sample>/outs/fragments.tsv.gz
output_root      = os.path.join(base_dir, "Filtering_results")
chrom_sizes_path = "/Users/jchoi/Desktop/Test/genome.chrom.sizes"  # rice chrom.sizes
annot_path       = "/Users/jchoi/Desktop/Test/genes.gtf"           # rice GTF (unzipped)

# Core thresholds
TSS_CUTOFF = 6.0
UMI_MIN    = 5000
UMI_MAX    = 50000

# Behavior toggles
ENABLE_H5AD_CACHE    = True    # Save h5ad for reuse
STRIP_BARCODE_SUFFIX = False   # If True, strip "-1" from barcodes (not typical for scATAC)

# Utilities
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def load_chrom_sizes(path: str) -> Dict[str, int]:
    """Load chrom.sizes into dict usable by SnapATAC2."""
    sizes: Dict[str, int] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom, size = parts[0], int(parts[1])
            sizes[chrom] = size
    if not sizes:
        raise ValueError(f"No entries parsed from chrom sizes: {path}")
    return sizes

def find_samples(root: str) -> List[Tuple[str, str]]:
    """
    Discover samples directly under base_dir (non-recursive).
    Returns list of (sample_name, fragments_path).
    """
    found: List[Tuple[str, str]] = []
    for name in sorted(os.listdir(root)):
        outs = os.path.join(root, name, "outs")
        frag = os.path.join(outs, "fragments.tsv.gz")
        if os.path.isfile(frag):
            found.append((name, frag))
    return found

def compute_tsse(adata: "snap.AnnData", gtf_path: str) -> None:
    """Robust TSSE call across SnapATAC2 versions (species-agnostic)."""
    try:
        snap.metrics.tsse(adata, gtf_path)   # preferred API
        return
    except TypeError:
        if hasattr(snap.genome, "read_gtf"):
            try:
                ann = snap.genome.read_gtf(gtf_path)
                snap.metrics.tsse(adata, ann)
                return
            except Exception:
                pass
        for kw in ("gtf_file", "annotation_file", "gtf"):
            try:
                snap.metrics.tsse(adata, **{kw: gtf_path})
                return
            except TypeError:
                continue
    raise RuntimeError("TSSE call failed: incompatible SnapATAC2 API for GTF input.")

def compute_tsse_map(fragment_path: str,
                     chrom_sizes: Dict[str, int],
                     gtf_path: str) -> Tuple[Dict[str, float], "snap.AnnData"]:
    """
    Import fragments with rice chrom sizes and compute TSSE using rice GTF.
    Returns ({barcode: tsse}, AnnData)
    """
    adata = snap.pp.import_fragments(
        fragment_path,
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=False,
    )
    compute_tsse(adata, gtf_path)
    tsse_map = dict(zip(adata.obs_names, adata.obs["tsse"]))
    return tsse_map, adata

def quick_chrom_name_sanity(fragment_path: str, chrom_sizes: Dict[str, int], n_lines: int = 1000) -> None:
    """Ensure fragment chromosome names exist in chrom_sizes dict."""
    with gzip.open(fragment_path, "rt") as f:
        for _ in range(n_lines):
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in chrom_sizes:
                raise ValueError(
                    f"Chrom '{chrom}' not found in chrom_sizes_dict — "
                    f"check rice chrom.sizes vs fragments' chromosome naming."
                )
            return  # stop after first valid record

def build_qc_from_fragments(adata: "snap.AnnData",
                            tsse_map: Dict[str, float]) -> pd.DataFrame:
    """
    Build QC dataframe using fragments only:
    - index: barcodes (adata.obs_names)
    - is__cell_barcode: 1 (treated as high-confidence by default)
    - passed_filters: n_fragment (as UMI-like proxy)
    - TSS_enrichment: from tsse_map
    """
    df = pd.DataFrame(index=adata.obs_names)
    if STRIP_BARCODE_SUFFIX:
        df.index = pd.Index(df.index).str.replace(r"-1$", "", regex=True)

    # All barcodes as cells (no singlecell.csv available)
    df["is__cell_barcode"] = 1

    # Use n_fragment as UMI-like metric
    if "n_fragment" in adata.obs:
        nfrag = adata.obs.reindex(df.index)["n_fragment"].fillna(0).astype(int)
    else:
        nfrag = pd.Series(0, index=df.index, dtype=int)
    df["passed_filters"] = nfrag

    # Attach TSSE
    df["TSS_enrichment"] = pd.Series(tsse_map).reindex(df.index).astype(float)

    return df

def filter_cells(df: pd.DataFrame,
                 tss_cutoff: float,
                 umi_min: int,
                 umi_max: int) -> pd.DataFrame:
    """
    Keep rows by TSSE and UMI range; '+1' follows prior logic.
    """
    required = ["is__cell_barcode", "passed_filters", "TSS_enrichment"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns in QC dataframe: {missing}")

    cells_only = df.query("is__cell_barcode == 1")
    filtered = cells_only[
        (cells_only["TSS_enrichment"] >= tss_cutoff) &
        (cells_only["passed_filters"] + 1 >= umi_min) &
        (cells_only["passed_filters"] + 1 <= umi_max)
    ]
    return filtered

# Main
def main():
    print(HEADER)
    print(desc_txt)

    # Pre-flight
    if annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip: use an uncompressed .gtf for 'annot_path'.")

    ensure_dir(output_root)
    out_filtered_dir = os.path.join(output_root, "filtered_samples")
    ensure_dir(out_filtered_dir)
    out_h5_dir = os.path.join(output_root, "h5ad")
    if ENABLE_H5AD_CACHE:
        ensure_dir(out_h5_dir)

    # Load rice chrom sizes
    chrom_sizes_dict = load_chrom_sizes(chrom_sizes_path)

    # Discover samples by fragments only
    samples = find_samples(base_dir)
    if not samples:
        print(f"[{ts()}] No samples found under: {base_dir}")
        # Optional hint
        try:
            print(" - First-level entries:", sorted(os.listdir(base_dir))[:10])
        except Exception:
            pass
        sys.exit(1)

    # Optional: quick chromosome naming sanity
    try:
        first_frag = samples[0][1]
        quick_chrom_name_sanity(first_frag, chrom_sizes_dict, n_lines=1000)
    except Exception as e:
        print(f"[{ts()}] Chromosome naming check failed: {e}")
        sys.exit(1)

    # Accumulators
    qc_df_dic: Dict[str, pd.DataFrame] = {}
    pre_counts: Dict[str, int] = {}
    post_counts: Dict[str, int] = {}
    skipped_samples: List[str] = []

    # Step 1: TSSE + QC from fragments; optional h5ad cache
    for sample, frag in samples:
        print(f"[{ts()}] Processing sample: {sample}")
        try:
            tsse_map, adata = compute_tsse_map(frag, chrom_sizes_dict, annot_path)

            if ENABLE_H5AD_CACHE:
                ad_out = os.path.join(out_h5_dir, f"{sample}.h5ad")
                if os.path.exists(ad_out):
                    os.remove(ad_out)
                adata.write_h5ad(ad_out)

            qc_df = build_qc_from_fragments(adata, tsse_map)
            qc_df_dic[sample] = qc_df

        except Exception as e:
            print(f"[{ts()}]  ! Skipped {sample} due to error: {e}")
            skipped_samples.append(sample)

    samples_kept = [s for s in samples if s[0] not in skipped_samples]

    # Step 2: filter per sample and save
    for sample, _ in samples_kept:
        df = qc_df_dic[sample]
        pre_counts[sample] = len(df)

        filtered = filter_cells(df, TSS_CUTOFF, UMI_MIN, UMI_MAX)
        post_counts[sample] = len(filtered)

        out_csv = os.path.join(out_filtered_dir, f"{sample}_filtered.csv")
        filtered.to_csv(out_csv, index=True)

        print(f"[{ts()}]  - Saved filtered cells: {sample} -> {out_csv}  "
              f"(kept {post_counts[sample]} / {pre_counts[sample]})")

    # Step 3: summary tables
    if pre_counts:
        summary = (
            pd.DataFrame({"pre_cells": pre_counts, "post_cells": post_counts})
              .assign(retention=lambda d: d["post_cells"] / d["pre_cells"])
              .sort_index()
        )
    else:
        summary = pd.DataFrame(columns=["pre_cells", "post_cells", "retention"])

    out_summary_csv = os.path.join(output_root, "summary_pre_post.csv")
    summary.to_csv(out_summary_csv)

    total_pre = int(summary["pre_cells"].sum()) if not summary.empty else 0
    total_post = int(summary["post_cells"].sum()) if not summary.empty else 0
    total_ret = (total_post / total_pre) if total_pre > 0 else float("nan")

    # Build human-readable TXT (data only; no raw logs)
    txt_lines: List[str] = []
    txt_lines.append(HEADER)
    txt_lines.append("Filtering Summary")
    txt_lines.append(f"Generated at: {ts()}")
    txt_lines.append("")
    txt_lines.append("Parameters")
    txt_lines.append(f"- Base dir           : {base_dir}")
    txt_lines.append(f"- chrom.sizes (rice) : {chrom_sizes_path}")
    txt_lines.append(f"- GTF (rice)         : {annot_path}")
    txt_lines.append(f"- TSS cutoff         : {TSS_CUTOFF}")
    txt_lines.append(f"- UMI min/max        : {UMI_MIN} / {UMI_MAX}")
    txt_lines.append(f"- Barcode suffix -1  : {STRIP_BARCODE_SUFFIX}")
    txt_lines.append("")
    if not summary.empty:
        txt_lines.append("Per-sample summary (pre/post/retention):")
        txt_lines.append(summary.to_string())
        txt_lines.append("")
        txt_lines.append(f"TOTAL: pre={total_pre}  post={total_post}  retention={total_ret:.4f}")
    else:
        txt_lines.append("No valid samples were processed.")
    if skipped_samples:
        txt_lines.append("")
        txt_lines.append(f"Skipped samples: {', '.join(skipped_samples)}")
    txt_lines.append("")
    txt_lines.append(f"(CSV saved at) {out_summary_csv}")

    out_summary_txt = os.path.join(output_root, "summary.txt")
    with open(out_summary_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(txt_lines))

    print(f"[{ts()}] Summary CSV  : {out_summary_csv}")
    print(f"[{ts()}] Summary TXT  : {out_summary_txt}")
    print(f"[{ts()}] Done. All outputs saved under: {output_root}")

if __name__ == "__main__":
    main()
