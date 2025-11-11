#!/usr/bin/env python3

import os
import sys
import gzip
import math
import warnings

import os as _os
_os.environ.setdefault('HDF5_USE_FILE_LOCKING','FALSE')

import time
import gc
warnings.simplefilter("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from scipy import sparse as sp
import snapatac2 as snap
from datetime import datetime
from joblib import Parallel, delayed
from typing import Dict, List, Tuple

try:
    import scanpy as sc
    _HAS_SCANPY = True
except Exception:
    _HAS_SCANPY = False

__version__ = "3.0"

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq: Filter + Doublet Removal"
HEADER += f"\n      Version {__version__}"
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

# ---- I/O defaults ----
chrom_sizes_path = None
annot_path        = None
base_dir          = None
out_dir           = None

n_jobs            = 1
TSS_CUTOFF        = 1.0
UMI_MIN           = 500
UMI_MAX           = 100000
N_FEATURES        = 250000
STRIP_BARCODE_SUFFIX = False

# ---- Utills ----
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

def _best_barcode_mapping(merged_idx, sample_idx, sample_name):
    mset = set(merged_idx)
    sset = set(sample_idx)
    inter = mset & sset
    if inter:
        return {k: k for k in inter}

    cand_maps = []
    def strip_prefix(v, pre):
        return v[len(pre):] if v.startswith(pre) else v
    def strip_suffix(v, suf):
        return v[:-len(suf)] if v.endswith(suf) else v

    prefixes = [f"{sample_name}_", f"{sample_name}-", f"{sample_name}:", f"{sample_name}#", f"{sample_name}."]
    suffixes = [f"_{sample_name}", f"-{sample_name}", f":{sample_name}", f"#{sample_name}", f".{sample_name}"]

    for pre in prefixes:
        mapped = {m: strip_prefix(m, pre) for m in merged_idx}
        inter2 = set(mapped.values()) & sset
        if inter2:
            d = {m: v for m, v in mapped.items() if v in sset}
            cand_maps.append(d)
    for suf in suffixes:
        mapped = {m: strip_suffix(m, suf) for m in merged_idx}
        inter2 = set(mapped.values()) & sset
        if inter2:
            d = {m: v for m, v in mapped.items() if v in sset}
            cand_maps.append(d)

    import re as _re
    def split_tokens(v):
        parts = _re.split(r'[:#._\-|]', v)
        return [p for p in parts if p]
    mapped_last = {m: (split_tokens(m)[-1] if split_tokens(m) else m) for m in merged_idx}
    inter3 = set(mapped_last.values()) & sset
    if inter3:
        d = {m: v for m, v in mapped_last.items() if v in sset}
        cand_maps.append(d)
    mapped_first = {m: (split_tokens(m)[0] if split_tokens(m) else m) for m in merged_idx}
    inter4 = set(mapped_first.values()) & sset
    if inter4:
        d = {m: v for m, v in mapped_first.items() if v in sset}
        cand_maps.append(d)

    if not cand_maps:
        return {}
    best = max(cand_maps, key=lambda d: len(d))
    return best

def augment_merged_with_nfragment(merged_path, out_dir):
    try:
        import anndata as ad
        import numpy as np
        import os, glob
    except Exception as e:
        print(f"[{ts()}] Warning: cannot import ad dependencies for augment step: {e}")
        return

    if not os.path.exists(merged_path):
        print(f"[{ts()}] Info: merged h5ads not found: {merged_path}")
        return

    try:
        merged = ad.read_h5ad(merged_path)
    except Exception as e:
        print(f"[{ts()}] Warning: failed to read merged h5ads: {e}")
        return

    if 'sample' not in merged.obs:
        print(f"[{ts()}] Info: 'sample' column absent in merged; cannot augment n_fragment.")
        return

    if 'n_fragment' not in merged.obs:
        merged.obs['n_fragment'] = np.nan

    samples = sorted(merged.obs['sample'].astype(str).unique().tolist())

    # expect pattern: <sample>_doublets.h5ad
    for s in samples:
        cand = os.path.join(out_dir, f"{s}_doublets.h5ad")
        if not os.path.exists(cand):
            gl = glob.glob(os.path.join(out_dir, f"*{s}*doublets*.h5ad"))
            cand = gl[0] if gl else None
        if not cand or not os.path.exists(cand):
            print(f"[{ts()}] Info: per-sample file not found for sample={s}; skip.")
            continue
        try:
            per = ad.read_h5ad(cand)
        except Exception as e:
            print(f"[{ts()}] Warning: read fail {cand}: {e}")
            continue
        if 'n_fragment' not in per.obs:
            print(f"[{ts()}] Info: n_fragment missing in {cand}; skip.")
            continue

        m_mask = (merged.obs['sample'].astype(str) == str(s))
        m_idx = merged.obs.index[m_mask]

        mapping = _best_barcode_mapping(m_idx, per.obs.index, s)
        if not mapping:
            print(f"[{ts()}] Info: could not map barcodes for sample={s}; skip.")
            continue
        per_vals = per.obs['n_fragment']
        count = 0
        for mkey, skey in mapping.items():
            try:
                merged.obs.at[mkey, 'n_fragment'] = per_vals.get(skey, np.nan)
                count += 1
            except Exception:
                pass
        print(f"[{ts()}] Augmented n_fragment: sample={s}, mapped={count}/{len(m_idx)}")

    try:
        merged.write_h5ad(merged_path)
        print(f"[{ts()}] Saved augmented merged h5ads with n_fragment: {merged_path}")
    except Exception as e:
        print(f"[{ts()}] Warning: failed to write augmented merged: {e}")

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
    _ensure_str_indices(adata)

    snap.pp.add_tile_matrix(adata)
    _ensure_str_indices(adata)
    _ensure_csr_matrix(adata)

    try:
        snap.pp.select_features(adata, n_features=N_FEATURES)
    except Exception:
        pass

    snap.pp.scrublet(adata, features=None)
    snap.pp.filter_doublets(adata)
    _ensure_str_indices(adata)
    _ensure_csr_matrix(adata)
    return adata

def write_summary(lines: List[str], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

def _safe_read_h5ads(path: str, tries: int = 3, delay: float = 0.35):
    try:
        import anndata as ad
    except Exception:
        return None
    last_err = None
    for _ in range(max(1, tries)):
        try:
            return ad.read_h5ad(path)
        except Exception as e:
            last_err = e
            time.sleep(delay)
    # keep silent (caller will log)
    return None

# ---------------- Main ----------------
def _apply_wrapper_injections_and_defaults() -> None:
    global base_dir, out_dir, chrom_sizes_path, annot_path
    global TSS_CUTOFF, UMI_MIN, UMI_MAX, N_FEATURES, STRIP_BARCODE_SUFFIX
    g = globals()

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

    if out_dir is None:
        out_dir = os.path.join(os.getcwd(), "Filter_results")
    os.makedirs(out_dir, exist_ok=True)

def main():
    _apply_wrapper_injections_and_defaults()

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

    filt_pre: Dict[str, int]  = {}
    filt_post: Dict[str, int] = {}
    dbl_pre: Dict[str, int]   = {}
    dbl_post: Dict[str, int]  = {}
    skipped: List[str]        = []

    
    def _process_one_sample(_sample: str, _frag: str) -> tuple:
        try:
            print(f"[{ts()}] Processing sample: {_sample}")
            # Import fragments
            adata_all = safe_import_fragments(_frag, chrom_sizes)

            # TSSE
            tsse_map = compute_tsse_map(adata_all, annot_path)

            # Build QC+filter
            qc_df = build_qc_from_fragments(adata_all, tsse_map)
            _filt_pre = qc_df.shape[0]

            keep_barcodes = filter_cells(qc_df, TSS_CUTOFF, UMI_MIN, UMI_MAX)
            _filt_post = len(keep_barcodes)

            out_csv = os.path.join(out_dir, f"{_sample}_filtered.csv")
            qc_df.loc[keep_barcodes, ["TSS_enrichment", "passed_filters"]].to_csv(out_csv, index=True)
            print(f"[{ts()}]  - Saved filtered cells: {_sample} -> {out_csv} (kept {_filt_post} / {_filt_pre})")

            if _filt_post == 0:
                return (_sample, _filt_pre, _filt_post, 0, 0, False, None)

            obs_idx_norm  = normalize_barcodes(pd.Index(adata_all.obs_names.astype(str)))
            keep_idx_norm = normalize_barcodes(pd.Index(keep_barcodes))
            mask = obs_idx_norm.isin(set(keep_idx_norm))
            adata_filt = adata_all[np.asarray(mask, dtype=bool), :].copy()
            _ensure_str_indices(adata_filt)

            _dbl_pre = adata_filt.shape[0]
            adata_nodbl = run_doublet_pipeline(adata_filt)
            _dbl_post = adata_nodbl.shape[0]

            out_h5 = os.path.join(out_dir, f"{_sample}_doublets.h5ad")
            adata_nodbl.write_h5ad(out_h5)
            print(f"[{ts()}]  - Saved: {out_h5}  (kept {_dbl_post} / {_dbl_pre})")
            return (_sample, _filt_pre, _filt_post, _dbl_pre, _dbl_post, False, None)
        except Exception as _e:
            print(f"[{ts()}]  - ERROR in sample {_sample}: {_e}")
            return (_sample, 0, 0, 0, 0, True, str(_e))

    if int(n_jobs) and int(n_jobs) > 1:
        _results = Parallel(n_jobs=int(n_jobs))(delayed(_process_one_sample)(s, f) for s, f in samples)
    else:
        _results = [_process_one_sample(s, f) for s, f in samples]

    for (_s, _fpre, _fpost, _dpre, _dpost, _skipped, _err) in _results:
        if _skipped:
            skipped.append(_s)
            continue
        filt_pre[_s] = _fpre
        filt_post[_s] = _fpost
        dbl_pre[_s] = _dpre
        dbl_post[_s] = _dpost

    if False:
        for sample, frag in samples:
            print(f"[{ts()}] Processing sample: {sample}")
            try:
                adata_all = safe_import_fragments(frag, chrom_sizes)
    
                tsse_map = compute_tsse_map(adata_all, annot_path)
    
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
    
                obs_idx_norm   = normalize_barcodes(pd.Index(adata_all.obs_names.astype(str)))
                keep_idx_norm  = normalize_barcodes(pd.Index(keep_barcodes))
                mask           = obs_idx_norm.isin(set(keep_idx_norm))
                adata_filt     = adata_all[np.asarray(mask, dtype=bool), :].copy()
                _ensure_str_indices(adata_filt)
    
                dbl_pre[sample] = adata_filt.shape[0]
    
                # Doublet pipeline
                adata_nodbl = run_doublet_pipeline(adata_filt)
                dbl_post[sample] = adata_nodbl.shape[0]
    
                # Save
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

    try:
        import re as _re
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
            del data; gc.collect()
            time.sleep(0.8)
        else:
            print(f"[{ts()}] Merge skipped (found {len(_h5s)} *_doublets.h5ad under {out_dir})")
    except Exception as _e:
        print(f"[{ts()}] Merge failed (strict concat): {_e}")

    try:
        time.sleep(0.35)
        _chk = _safe_read_h5ads(merge_output_path, tries=3, delay=0.4)
        if _chk is None:
            print(f"[{ts()}] Warning: failed to read merged h5ads immediately; augment anyway.")
        try:
            augment_merged_with_nfragment(merge_output_path, out_dir)
            print(f"[{ts()}]  - Augmented merged h5ads with n_fragment.")
        except Exception as _e_aug:
            print(f"[{ts()}] Warning: augment_merged_with_nfragment failed: {_e_aug}")
        _verify = _safe_read_h5ads(merge_output_path, tries=3, delay=0.4)
        if _verify is not None:
            try:
                _cols = list(getattr(_verify, 'obs', {}).columns)
            except Exception:
                _cols = []
            print(f"[{ts()}]  - merged.obs columns: {_cols}")
    except Exception as _e_post:
        print(f"[{ts()}] Warning: post-merge verification failed: {_e_post}")
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

    try:
        merged_path = os.path.join(out_dir, "merged_doublets.h5ads")
        augment_merged_with_nfragment(merged_path, out_dir)
    except Exception as e:
        print(f"[{ts()}] Warning: augment step failed: {e}")

    out_summary = os.path.join(out_dir, "summary.txt")
    write_summary(lines, out_summary)
    print(f"[{ts()}] Summary written: {out_summary}")
    print(f"[{ts()}] Done. All outputs saved under: {out_dir}")

if __name__ == "__main__":
    main()
