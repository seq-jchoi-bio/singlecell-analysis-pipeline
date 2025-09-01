#!/usr/bin/env python3
# Filtering: Integrated scATAC-seq cell filtering using fragments.tsv.gz only (rice-ready).
# Outputs (under Filter_results/):
#  - <sample>_filtered.csv    # Per-sample filtered cells
#  - summary.txt              # Human-readable summary (no raw logs)

import os
import sys
import gzip
import pandas as pd
import snapatac2 as snap
from datetime import datetime
from typing import Dict, List, Tuple

__version__ = 1.2

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq Cell Filter"
HEADER += "\n      Version {0}".format(__version__)
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

desc_txt = """This tool aggregates per-sample QC, attaches TSSE scores computed from fragments,
and filters cells with (1) high-confidence barcodes, (2) TSSE >= cutoff,
and (3) UMI within the specified min/max range. Results are saved per sample,
along with a concise summary across all samples."""

# =========================
# User-configurable section
# =========================
base_dir         = "examples"   # <sample>/outs/fragments.tsv.gz
output_root      = os.path.join("results", "Filter_results")
chrom_sizes_path = "work/genome.chrom.sizes"
annot_path       = "work/genes.gtf"

# Core thresholds
TSS_CUTOFF = 4.69
UMI_MIN    = 5000
UMI_MAX    = 60000

# Behavior toggles
STRIP_BARCODE_SUFFIX = False   # If True, strip "-1" from barcodes (not typical for scATAC)

# =========================
# Utilities
# =========================
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def load_chrom_sizes(path: str) -> Dict[str, int]:
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
    found: List[Tuple[str, str]] = []
    for name in sorted(os.listdir(root)):
        outs = os.path.join(root, name, "outs")
        frag = os.path.join(outs, "fragments.tsv.gz")
        if os.path.isfile(frag):
            found.append((name, frag))
    return found

def compute_tsse(adata: "snap.AnnData", gtf_path: str) -> None:
    try:
        snap.metrics.tsse(adata, gtf_path)
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
    adata = snap.pp.import_fragments(
        fragment_path,
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=False,
    )
    compute_tsse(adata, gtf_path)
    tsse_map = dict(zip(adata.obs_names, adata.obs["tsse"]))
    return tsse_map, adata

def quick_chrom_name_sanity(fragment_path: str, chrom_sizes: Dict[str, int], n_lines: int = 1000) -> None:
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
            return

def build_qc_from_fragments(adata: "snap.AnnData",
                            tsse_map: Dict[str, float]) -> pd.DataFrame:
    df = pd.DataFrame(index=adata.obs_names)
    if STRIP_BARCODE_SUFFIX:
        df.index = pd.Index(df.index).str.replace(r"-1$", "", regex=True)
    df["is__cell_barcode"] = 1
    if "n_fragment" in adata.obs:
        nfrag = adata.obs.reindex(df.index)["n_fragment"].fillna(0).astype(int)
    else:
        nfrag = pd.Series(0, index=df.index, dtype=int)
    df["passed_filters"] = nfrag
    df["TSS_enrichment"] = pd.Series(tsse_map).reindex(df.index).astype(float)
    return df

def filter_cells(df: pd.DataFrame,
                 tss_cutoff: float,
                 umi_min: int,
                 umi_max: int) -> pd.DataFrame:
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

# =========================
# Main
# =========================
def main():
    print(HEADER)
    print(desc_txt)

    if annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip: use an uncompressed .gtf for 'annot_path'.")

    ensure_dir(output_root)

    chrom_sizes_dict = load_chrom_sizes(chrom_sizes_path)

    samples = find_samples(base_dir)
    if not samples:
        print(f"[{ts()}] No samples found under: {base_dir}")
        try:
            print(" - First-level entries:", sorted(os.listdir(base_dir))[:10])
        except Exception:
            pass
        sys.exit(1)

    try:
        first_frag = samples[0][1]
        quick_chrom_name_sanity(first_frag, chrom_sizes_dict, n_lines=1000)
    except Exception as e:
        print(f"[{ts()}] Chromosome naming check failed: {e}")
        sys.exit(1)

    pre_counts: Dict[str, int] = {}
    post_counts: Dict[str, int] = {}
    skipped_samples: List[str] = []

    for sample, frag in samples:
        print(f"[{ts()}] Processing sample: {sample}")
        try:
            tsse_map, adata = compute_tsse_map(frag, chrom_sizes_dict, annot_path)
            qc_df = build_qc_from_fragments(adata, tsse_map)
            pre_counts[sample] = len(qc_df)
            filtered = filter_cells(qc_df, TSS_CUTOFF, UMI_MIN, UMI_MAX)
            post_counts[sample] = len(filtered)
            out_csv = os.path.join(output_root, f"{sample}_filtered.csv")
            filtered.to_csv(out_csv, index=True)
            print(f"[{ts()}]  - Saved filtered cells: {sample} -> {out_csv}  "
                  f"(kept {post_counts[sample]} / {pre_counts[sample]})")
        except Exception as e:
            print(f"[{ts()}]  ! Skipped {sample} due to error: {e}")
            skipped_samples.append(sample)

    total_pre = sum(pre_counts.values()) if pre_counts else 0
    total_post = sum(post_counts.values()) if post_counts else 0
    total_ret = (total_post / total_pre) if total_pre > 0 else float("nan")

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
    if pre_counts:
        txt_lines.append("Per-sample summary (pre/post/retention):")
        for sample in pre_counts:
            pre = pre_counts[sample]
            post = post_counts.get(sample, 0)
            ret = (post / pre) if pre > 0 else float("nan")
            txt_lines.append(f"{sample}: pre={pre}  post={post}  retention={ret:.4f}")
        txt_lines.append("")
        txt_lines.append(f"TOTAL: pre={total_pre}  post={total_post}  retention={total_ret:.4f}")
    else:
        txt_lines.append("No valid samples were processed.")
    if skipped_samples:
        txt_lines.append("")
        txt_lines.append(f"Skipped samples: {', '.join(skipped_samples)}")

    out_summary_txt = os.path.join(output_root, "summary.txt")
    with open(out_summary_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(txt_lines))

    print(f"[{ts()}] Summary TXT  : {out_summary_txt}")
    print(f"[{ts()}] Done. All outputs saved under: {output_root}")

if __name__ == "__main__":
    main()
