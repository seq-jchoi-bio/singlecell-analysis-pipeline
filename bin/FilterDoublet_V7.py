#!/usr/bin/env python3

import os
import sys
import gzip
import math
import warnings

warnings.simplefilter("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from scipy import sparse as sp
import snapatac2 as snap
from datetime import datetime
from typing import Dict, List, Tuple

try:
    import scanpy as sc
    _HAS_SCANPY = True
except Exception:
    _HAS_SCANPY = False

__version__ = "2.2"

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq: Filter + Doublet Removal"
HEADER += f"\n      Version {__version__}"
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

desc_txt = """This script:
  (1) computes TSSE from fragments + GTF,
  (2) builds a minimal QC table (is__cell_barcode, passed_filters=n_fragment, TSSE),
  (3) filters cells by TSSE >= cutoff and UMI within [min, max],
  (4) applies Scrublet-based doublet removal on filtered cells,
  and saves one consolidated summary under Filter_results/summary.txt."""

# =========================
# User-configurable section (wrapper-friendly)
# =========================
chrom_sizes_path = None  # wrapper may inject
annot_path        = None  # wrapper may inject
base_dir          = None  # wrapper may inject
out_dir           = None  # wrapper sets to CWD/Filter_results by default

# Thresholds / knobs
TSS_CUTOFF        = 1.0   # default; wrapper -tss can override
UMI_MIN           = 500   # default; wrapper -min can override
UMI_MAX           = 100000 # default; wrapper -max can override
N_FEATURES        = 250000 # default; wrapper -feat can override
STRIP_BARCODE_SUFFIX = False # wrapper --strip can set True

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
            sizes[parts[0]] = int(parts[1])
    if not sizes:
        raise ValueError(f"No entries parsed from chrom sizes: {path}")
    return sizes

def find_samples(root: str) -> List[Tuple[str, str]]:
    found: List[Tuple[str, str]] = []
    for name in sorted(os.listdir(root)):
        frag = os.path.join(root, name, "outs", "fragments.tsv.gz")
        if os.path.isfile(frag):
            found.append((name, frag))
    return found

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
                    f"Chrom '{chrom}' not found in chrom.sizes — "
                    f"check rice chrom.sizes vs fragments' chromosome naming."
                )
            return

def _ensure_str_indices(adata: "snap.AnnData") -> None:
    try:
        adata.obs_names = adata.obs_names.astype(str)
    except Exception:
        pass
    try:
        adata.var_names = adata.var_names.astype(str)
    except Exception:
        pass

def safe_import_fragments(fragment_path: str, chrom_sizes: Dict[str, int]) -> "snap.AnnData":
    try:
        adata = snap.pp.import_fragments(
            fragment_path, chrom_sizes=chrom_sizes, sorted_by_barcode=False,
        )
    except Exception:
        adata = snap.pp.import_data(  # older API
            fragment_path, chrom_sizes=chrom_sizes, sorted_by_barcode=False,
        )
    _ensure_str_indices(adata)
    return adata

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

def compute_tsse_map(adata: "snap.AnnData", gtf_path: str) -> Dict[str, float]:
    compute_tsse(adata, gtf_path)
    return dict(zip(adata.obs_names, adata.obs["tsse"]))

def build_qc_from_fragments(adata: "snap.AnnData",
                            tsse_map: Dict[str, float]) -> pd.DataFrame:
    df = pd.DataFrame(index=adata.obs_names)
    if STRIP_BARCODE_SUFFIX:
        df.index = pd.Index(df.index).str.replace(r"-\d+$", "", regex=True)
    df["is__cell_barcode"] = 1
    if "n_fragment" in adata.obs:
        nfrag = adata.obs.reindex(adata.obs_names)["n_fragment"].fillna(0)
        try:
            nfrag = nfrag.astype(int)
        except Exception:
            nfrag = nfrag.astype(float).round().astype(int)
            nfrag = pd.Series(0, index=adata.obs_names, dtype=int)
    df["passed_filters"]   = nfrag
    df["TSS_enrichment"]   = pd.Series(tsse_map).reindex(df.index).astype(float)
    return df

def filter_cells(df: pd.DataFrame,
                 tss_cutoff: float,
                 umi_min: int,
                 umi_max: int) -> pd.Index:
    required = ["is__cell_barcode", "passed_filters", "TSS_enrichment"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns in QC dataframe: {missing}")
    cells_only = df.query("is__cell_barcode == 1")
    mask = (
        (cells_only["TSS_enrichment"] >= tss_cutoff) &
        (cells_only["passed_filters"] + 1 >= umi_min) &
        (cells_only["passed_filters"] + 1 <= umi_max)
    )
    return cells_only.index[mask]

def normalize_barcodes(idx: pd.Index) -> pd.Index:
    if STRIP_BARCODE_SUFFIX:
        return pd.Index(idx).str.replace(r"-\d+$", "", regex=True)
    return idx

def _ensure_csr_matrix(adata: "snap.AnnData") -> None:
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.X = adata.X.tocsr()

def run_doublet_pipeline(adata: "snap.AnnData") -> "snap.AnnData":
    """Doublet pipeline (no feature slicing; keep full grid for merge safety)."""
    _ensure_str_indices(adata)

    # 1) Full grid
    snap.pp.add_tile_matrix(adata)
    _ensure_str_indices(adata)
    _ensure_csr_matrix(adata)

    try:
        snap.pp.select_features(adata, n_features=N_FEATURES)
    except Exception:
        pass

    # 2) Scrublet on full grid
    snap.pp.scrublet(adata, features=None)
    snap.pp.filter_doublets(adata)
    _ensure_str_indices(adata)
    _ensure_csr_matrix(adata)
    return adata

def write_summary(lines: List[str], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

# =========================
# Main
# =========================
# Wrapper injection
def _apply_wrapper_injections_and_defaults() -> None:
    global base_dir, out_dir, chrom_sizes_path, annot_path
    global TSS_CUTOFF, UMI_MIN, UMI_MAX, N_FEATURES, STRIP_BARCODE_SUFFIX
    g = globals()

    # Accept both lowercase and uppercase names from wrappers
    if g.get("BASE_DIR") is not None:
        base_dir = g.get("BASE_DIR")
    if g.get("OUT_DIR") is not None:
        out_dir = g.get("OUT_DIR")
    if g.get("CHROM_SIZES_PATH") is not None and (g.get("chrom_sizes_path") is None):
        chrom_sizes_path = g.get("CHROM_SIZES_PATH")
    if g.get("GTF_PATH") is not None and (g.get("annot_path") is None):
        annot_path = g.get("GTF_PATH")
    # Lowercase overrides (already common in this module)
    if g.get("chrom_sizes_path") is not None:
        chrom_sizes_path = g.get("chrom_sizes_path")
    if g.get("annot_path") is not None:
        annot_path = g.get("annot_path")
    if g.get("base_dir") is not None:
        base_dir = g.get("base_dir")
    if g.get("out_dir") is not None:
        out_dir = g.get("out_dir")

    # Ensure out_dir default
    if out_dir is None:
        out_dir = os.path.join(os.getcwd(), "Filter_results")
    os.makedirs(out_dir, exist_ok=True)

def main():
    _apply_wrapper_injections_and_defaults()
    print(HEADER)
    print(desc_txt)

    if annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip the GTF: provide a plain .gtf for 'annot_path'.")

    ensure_dir(out_dir)
    chrom_sizes = load_chrom_sizes(chrom_sizes_path)

    samples = find_samples(base_dir)
    if not samples:
        print(f"[{ts()}] No samples found under: {base_dir}")
        try:
            print(" - First-level entries:", sorted(os.listdir(base_dir))[:10])
        except Exception:
            pass
        sys.exit(1)

    try:
        quick_chrom_name_sanity(samples[0][1], chrom_sizes, n_lines=1000)
    except Exception as e:
        print(f"[{ts()}] Chromosome naming check failed: {e}")
        sys.exit(1)

    filt_pre: Dict[str, int]  = {}   # before TSSE/UMI filtering
    filt_post: Dict[str, int] = {}   # after TSSE/UMI filtering
    dbl_pre: Dict[str, int]   = {}   # input to doublet removal (= filt_post)
    dbl_post: Dict[str, int]  = {}   # after doublet removal
    skipped: List[str]        = []

    for sample, frag in samples:
        print(f"[{ts()}] Processing sample: {sample}")
        try:
            # Import fragments (one AnnData per sample)
            adata_all = safe_import_fragments(frag, chrom_sizes)

            # TSSE
            tsse_map = compute_tsse_map(adata_all, annot_path)

            # Build QC+filter
            qc_df = build_qc_from_fragments(adata_all, tsse_map)
            filt_pre[sample] = qc_df.shape[0]

            keep_barcodes = filter_cells(qc_df, TSS_CUTOFF, UMI_MIN, UMI_MAX)
            filt_post[sample] = len(keep_barcodes)

            out_csv = os.path.join(out_dir, f"{sample}_filtered.csv")
            qc_df.loc[keep_barcodes, ["TSS_enrichment", "passed_filters"]].to_csv(out_csv, index=True)
            print(f"[{ts()}]  - Saved filtered cells: {sample} -> {out_csv} "
                  f"(kept {filt_post[sample]} / {filt_pre[sample]})")

            if filt_post[sample] == 0:
                print(f"[{ts()}]  - No cells passed thresholds; skip doublet removal.")
                dbl_pre[sample]  = 0
                dbl_post[sample] = 0
                continue

            # Subset to filtered barcodes on the same in-memory AnnData
            obs_idx_norm   = normalize_barcodes(pd.Index(adata_all.obs_names.astype(str)))
            keep_idx_norm  = normalize_barcodes(pd.Index(keep_barcodes))
            mask           = obs_idx_norm.isin(set(keep_idx_norm))
            adata_filt     = adata_all[np.asarray(mask, dtype=bool), :].copy()
            _ensure_str_indices(adata_filt)

            dbl_pre[sample] = adata_filt.shape[0]

            # Doublet pipeline (mode-dependent)
            adata_nodbl = run_doublet_pipeline(adata_filt)
            dbl_post[sample] = adata_nodbl.shape[0]

            # Save h5ad to Filter_results
            out_h5 = os.path.join(out_dir, f"{sample}_doublets.h5ad")
            adata_nodbl.write_h5ad(out_h5)
            print(f"[{ts()}]  - Saved: {out_h5}  (kept {dbl_post[sample]} / {dbl_pre[sample]})")

        except Exception as e:
            print(f"[{ts()}]  ! Skipped {sample} due to error: {e}")
            skipped.append(sample)
            filt_pre.setdefault(sample, 0)
            filt_post.setdefault(sample, 0)
            dbl_pre.setdefault(sample, 0)
            dbl_post.setdefault(sample, 0)
            continue

    # Merge step (strict simple concatenation)
    try:
        import re as _re
        # gather per-sample outputs
        _h5s = []
        for _name in sorted(os.listdir(out_dir)):
            if _name.endswith("_doublets.h5ad"):
                _path = os.path.join(out_dir, _name)
                _m = _re.match(r"(.+?)_doublets\.h5ad$", _name)
                _sample = _m.group(1) if _m else os.path.splitext(_name)[0]
                _h5s.append((_sample, _path))
        if len(_h5s) >= 2:
            merge_output_path = os.path.join(out_dir, "merged_doublets.h5ads")
            print(f"[{ts()}] Merging {len(_h5s)} files -> {merge_output_path}")
            data = snap.AnnDataSet(adatas=_h5s, filename=merge_output_path)
            print(f"[{ts()}]  - Merged data saved to: {merge_output_path}")
        else:
            print(f"[{ts()}] Merge skipped (found {len(_h5s)} *_doublets.h5ad under {out_dir})")
    except Exception as _e:
        print(f"[{ts()}] Merge failed (strict concat): {_e}")
        print(f"[{ts()}] Hint: ensure all inputs share identical var_names (features & order).")

    lines: List[str] = []
    lines.append(HEADER)
    lines.append("")
    lines.append("===Consolidated Summary (Filter + Doublet Removal)===")
    lines.append(f"Generated at: {ts()}")
    lines.append("")
    lines.append("Parameters")
    lines.append(f"- Base dir             : {base_dir}")
    lines.append(f"- Output root          : {out_dir}")
    lines.append(f"- chrom.sizes (rice)   : {chrom_sizes_path}")
    lines.append(f"- GTF (rice)           : {annot_path}")
    lines.append(f"- TSS cutoff           : {TSS_CUTOFF}")
    lines.append(f"- UMI min/max          : {UMI_MIN} / {UMI_MAX}")
    lines.append(f"- ")
    lines.append(f"- N_FEATURES (per-sample): {N_FEATURES}")
    lines.append(f"- Strip barcode suffix : {STRIP_BARCODE_SUFFIX}")
    lines.append(f"- Scanpy fallback      : {_HAS_SCANPY}")
    lines.append("")

    if filt_pre:
        lines.append("A) Filtering (TSSE/UMI): pre → post → retention")
        tot_f_pre  = sum(filt_pre.values())
        tot_f_post = sum(filt_post.get(s, 0) for s in filt_pre)
        for s in filt_pre:
            pre = filt_pre[s]
            post = filt_post.get(s, 0)
            ret = (post / pre) if pre > 0 else float("nan")
            lines.append(f"- {s}: {pre} → {post}  (retention={ret:.4f})")
        lines.append(f"TOTAL: {tot_f_pre} → {tot_f_post}  (retention={(tot_f_post/tot_f_pre if tot_f_pre>0 else float('nan')):.4f})")
        lines.append("")

        lines.append("B) Doublet removal (Scrublet): pre → post → retention")
        tot_d_pre  = sum(dbl_pre.values())
        tot_d_post = sum(dbl_post.get(s, 0) for s in dbl_pre)
        for s in dbl_pre:
            pre = dbl_pre[s]
            post = dbl_post.get(s, 0)
            ret = (post / pre) if pre > 0 else float("nan")
            lines.append(f"- {s}: {pre} → {post}  (retention={ret:.4f})")
        lines.append(f"TOTAL: {tot_d_pre} → {tot_d_post}  (retention={(tot_d_post/tot_d_pre if tot_d_pre>0 else float('nan')):.4f})")
        lines.append("No valid samples were processed.")

    if skipped:
        lines.append("")
        lines.append("Skipped samples:")
        lines.append(", ".join(skipped))

    out_summary = os.path.join(out_dir, "summary.txt")
    write_summary(lines, out_summary)
    print(f"[{ts()}] Summary written: {out_summary}")
    print(f"[{ts()}] Done. All outputs saved under: {out_dir}")

if __name__ == "__main__":
    main()