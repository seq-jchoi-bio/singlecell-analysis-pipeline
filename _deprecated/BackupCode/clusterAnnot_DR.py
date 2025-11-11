#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Annot_V5 (no-BED) with NumPy>=2.0 compatibility shim for snapatac2.leiden
- Adds np.compat.unicode = str and legacy dtype aliases before leiden()
"""

import os
import re
from datetime import datetime
from typing import Dict, Optional, List, Tuple
from collections import defaultdict

# Thread caps
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"]      = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"
os.environ["RAYON_NUM_THREADS"]    = "1"
os.environ["MALLOC_ARENA_MAX"]     = "2"

# Plot backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Core libs
import numpy as np
import pandas as pd
import anndata as an
import snapatac2 as snap
from scipy import sparse

import warnings
warnings.filterwarnings("ignore", message=r"`_import_from_c` is deprecated", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden", category=UserWarning, module=r"umap\.umap_")

HEADER = (
    "=================================================================================\n"
    "  Integrated scATAC-seq: Clustering + Annotation (no-BED)\n"
    "  Version 2.9-nobed (np2fix)\n"
    "================================================================================="
)
desc_txt = ""

# ---------------- User paths (relative OK) ----------------
BASE_DIR        = "Test"
H5ADS_INPUT     = "Test/Filter_results/merged_doublets.h5ads"
AGGR_DIR        = "Test/Annot_results"
FIG_DIR         = "Test/Annot_results/Plots"

H5ADS_OUT       = "Test/Annot_results/annot_merged.h5ads"
ANNOT_H5AD      = "Test/Annot_results/annot_merged_cells.h5ad"
ANNOT_GENE_H5AD = "Test/Annot_results/annot_gene_activity.h5ad"
SUMMARY_TXT     = "Test/Annot_results/summary.txt"
CSV_OUT         = "Test/Annot_results/annotation_results.csv"

CHROM_SIZES_PATH = os.path.expanduser("~/singlecell-analysis-pipeline/refGenome/rice/genome.chrom.sizes")
GTF_PATH         = os.path.expanduser("~/singlecell-analysis-pipeline/refGenome/rice/genes.gtf")
REF_CSV_PATH     = os.path.expanduser("~/singlecell-analysis-pipeline/refGenome/rice/PCMDB_Marker_ENDO.csv")

# Defaults
REF_LABEL_KEY       = "celltype_id"
MARKER_MIN_LOG2FC   = 0.5
MARKER_MAX_PADJ     = 0.05
MIN_MARKERS_PER_SET = 5
NO_SCORE_MODE       = "strict"
NO_SCORE_UBI_K      = 3
NO_SCORE_CAP        = None
NO_SCORE_PRIORITY   = None
N_FEATURES          = 250_000
N_NEIGHBORS         = 10
USE_DIMS            = None
MAKE_GENE_ACTIVITY  = True
TRANSFER_MODE       = "auc"
UNKNOWN_CONF_THRESH = 0.5
SAVE_UMAP           = True
LABEL_MIN_COUNT     = 20
LABEL_OUTLINE       = True

# ---------------- Utils ----------------
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dirs() -> None:
    os.makedirs(AGGR_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

def write_summary(lines: List[str]) -> None:
    with open(SUMMARY_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[{ts()}] Summary written: {SUMMARY_TXT}")

def _patch_numpy_for_leiden():
    """Shim NumPy>=2 removals used by dependencies (leiden/harmonypy/sklearn)."""
    try:
        # Legacy dtype aliases
        for name, val in {"float": float, "int": int, "bool": bool, "object": object}.items():
            if not hasattr(np, name):
                setattr(np, name, val)
        # compat namespace
        compat = getattr(np, "compat", None)
        if compat is None:
            class _Compat: ...
            np.compat = _Compat()
            compat = np.compat
        if not hasattr(compat, "unicode"):
            setattr(compat, "unicode", str)
    except Exception:
        pass

# Plot helpers
def _obs_as_str_array(adata, key: str) -> np.ndarray:
    try:
        vals = adata.obs[key]
    except Exception:
        return np.array(["NA"] * adata.n_obs, dtype=object)
    try:
        arr = np.asarray(vals, dtype=object)
    except Exception:
        if hasattr(vals, "to_numpy"):
            arr = vals.to_numpy(dtype=object)
        else:
            arr = np.array(list(vals), dtype=object)
    return np.asarray(arr, dtype=object)

def _scatter_by_category(X: np.ndarray, col: np.ndarray, embed_name: str, basename: str, figdir: str) -> None:
    cats = np.array(sorted(np.unique(col), key=str))
    cmap = plt.get_cmap("tab20", len(cats))
    colors = cmap(range(len(cats)))

    fig, ax = plt.subplots(figsize=(6, 6))
    for i, cat in enumerate(cats):
        m = (col == cat)
        ax.scatter(X[m, 0], X[m, 1], s=5, alpha=0.5, color=colors[i], linewidths=0, label=str(cat))

    ax.set_xlabel(f"{embed_name} 1")
    ax.set_ylabel(f"{embed_name} 2")
    ax.set_title(f"{embed_name} by category")
    ax.set_aspect("equal", adjustable="box")

    if len(cats) <= 60:
        ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left",
                  fontsize=7, frameon=False, markerscale=1.6, handlelength=1.0, handletextpad=0.4)

    os.makedirs(figdir, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 0.80, 1.0])
    for ext in ("png", "svg"):
        out = os.path.join(figdir, f"{basename}.{ext}")
        fig.savefig(out, dpi=300, facecolor="white")
        print(f"[{ts()}] Saved figure: {out}")
    plt.close(fig)

def plot_umap_matplotlib(adata, color: str, basename: str, figdir: str) -> None:
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    col = _obs_as_str_array(adata, color)
    _scatter_by_category(X, col, "UMAP", basename, figdir)

# Chromosome normalization
def _normalize_chrom(s: str) -> str:
    s0 = str(s).strip()
    s0 = re.sub(r'^(chr|Chr)', '', s0)
    s0 = re.sub(r'^NC_[0-9.]+', '', s0)
    return s0.lower()

# GTF handling
def _parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_str):
        out[m.group(1)] = m.group(2)
    return out

def detect_has_gene_rows(gtf_in: str, max_check: int = 100000) -> bool:
    cnt = 0
    with open(gtf_in, "r") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2].lower() == "gene":
                return True
            cnt += 1
            if cnt >= max_check:
                break
    return False

def build_gene_only_gtf_from_any(gtf_in: str, gtf_out: str) -> str:
    pat_gene_id = re.compile(r'gene_id\s+"([^"]+)"')
    pat_gene_nm = re.compile(r'gene_name\s+"([^"]+)"')

    has_gene = detect_has_gene_rows(gtf_in)
    if has_gene:
        with open(gtf_in, "r") as fin, open(gtf_out, "w") as fout:
            for line in fin:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2].lower() != "gene":
                    continue
                attrs = cols[8]
                if not pat_gene_id.search(attrs):
                    continue
                if not pat_gene_nm.search(attrs):
                    if attrs and attrs[-1] != ";":
                        attrs += ";"
                    gid = pat_gene_id.search(attrs).group(1)
                    attrs += f' gene_name "{gid}";'
                    cols[8] = attrs
                fout.write("\t".join(cols) + "\n")
        return gtf_out

    agg: Dict[str, Dict[str, object]] = {}
    with open(gtf_in, "r") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2].lower()
            if feature not in ("transcript", "exon", "cds", "mrna"):
                continue
            seqname, source, start, end, score, strand, frame, attrs = cols[0], cols[1], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
            try:
                start_i = int(start); end_i = int(end)
            except Exception:
                continue
            a = _parse_gtf_attrs(attrs)
            gid = a.get("gene_id", None)
            if gid is None:
                continue
            rec = agg.get(gid)
            if rec is None:
                agg[gid] = {
                    "seqname": seqname,
                    "source": source if source else "collapsed",
                    "start": start_i,
                    "end": end_i,
                    "score": ".",
                    "strand": strand if strand in ("+", "-") else ".",
                    "frame": ".",
                }
            else:
                rec["start"] = min(rec["start"], start_i)
                rec["end"]   = max(rec["end"], end_i)
                if rec["strand"] not in ("+", "-") and strand in ("+", "-"):
                    rec["strand"] = strand

    if not agg:
        raise RuntimeError("Failed to build gene-only GTF: no usable gene_id found.")

    with open(gtf_out, "w") as fout:
        for gid, rec in agg.items():
            attrs = f'gene_id "{gid}"; gene_name "{gid}";'
            cols = [
                str(rec["seqname"]),
                str(rec["source"]),
                "gene",
                str(rec["start"]),
                str(rec["end"]),
                str(rec["score"]),
                str(rec["strand"]),
                str(rec["frame"]),
                attrs
            ]
            fout.write("\t".join(cols) + "\n")
    return gtf_out

def make_gene_matrix_geneid_only(data, gtf_path: str, chrom_sizes: Optional[Dict[str, int]] = None):
    norm_gtf = os.path.join(os.path.dirname(gtf_path), "genes.gene_only.auto.gtf")
    norm_gtf = build_gene_only_gtf_from_any(gtf_path, norm_gtf)
    try:
        return snap.pp.make_gene_matrix(data, annotation=norm_gtf, id_type="gene_id")
    except Exception:
        pass
    try:
        return snap.pp.make_gene_matrix(data, annotation=norm_gtf)
    except Exception:
        pass
    try:
        if chrom_sizes is not None:
            return snap.pp.make_gene_matrix(data, gtf=norm_gtf, chrom_sizes=chrom_sizes)
        else:
            return snap.pp.make_gene_matrix(data, gtf=norm_gtf)
    except Exception:
        pass
    try:
        if hasattr(snap.genome, "read_gtf"):
            ann = snap.genome.read_gtf(norm_gtf)
            try:
                return snap.pp.make_gene_matrix(data, annotation=ann)
            except Exception:
                return snap.pp.make_gene_matrix(data, ann)
    except Exception:
        pass
    raise RuntimeError("make_gene_matrix failed with auto gene-only GTF.")

# Fallback gene activity
def _parse_peak_name(name: str):
    s = str(name)
    m = re.match(r'^([^:]+):(\d+)-(\d+)$', s)
    if not m:
        m = re.match(r'^([^_]+)_(\d+)_(\d+)$', s)
    if not m:
        return None
    chrom = _normalize_chrom(m.group(1))
    start = int(m.group(2)); end = int(m.group(3))
    if end < start:
        start, end = end, start
    return chrom, start, end

def _load_gene_intervals(gtf_path: str, flank_up=2000, flank_down=0) -> Dict[str, List[Tuple[int,int,str]]]:
    genes_by_chr: Dict[str, List[Tuple[int,int,str]]] = defaultdict(list)
    pat_gid = re.compile(r'gene_id\s+"([^"]+)"')
    with open(gtf_path, "r") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2].lower() != "gene":
                continue
            chrom = _normalize_chrom(cols[0])
            try:
                start = int(cols[3]); end = int(cols[4])
            except Exception:
                continue
            m = pat_gid.search(cols[8])
            if not m:
                continue
            gid = m.group(1)
            s = max(0, start - int(flank_up))
            e = end + int(flank_down)
            if e <= s:
                e = s + 1
            genes_by_chr[chrom].append((s, e, gid))
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda x: x[0])
    return genes_by_chr

def _map_peaks_to_genes(peak_coords_by_chr, gene_intervals_by_chr):
    g2p: Dict[str, set] = defaultdict(set)
    for chrom, peaks in peak_coords_by_chr.items():
        genes = gene_intervals_by_chr.get(chrom, [])
        if not genes or not peaks:
            continue
        gi = 0
        for (p_start, p_end, p_idx) in peaks:
            center = (p_start + p_end) // 2
            while gi < len(genes) and genes[gi][1] < center:
                gi += 1
            for k in (gi-2, gi-1, gi, gi+1, gi+2):
                if 0 <= k < len(genes):
                    g_start, g_end, gid = genes[k]
                    if g_start <= center <= g_end:
                        g2p[gid].add(p_idx)
    return g2p

def build_gene_activity_fallback(data: snap.AnnDataSet, gene_gtf: str, flank_up=2000, flank_down=0) -> an.AnnData:
    ad = data.to_adata()
    try:
        ad.obs_names_make_unique()
    except Exception:
        pass
    peak_coords_by_chr: Dict[str, List[Tuple[int,int,int]]] = defaultdict(list)
    for j, name in enumerate(map(str, ad.var_names)):
        parsed = _parse_peak_name(name)
        if parsed is None:
            continue
        chrom, s, e = parsed
        peak_coords_by_chr[chrom].append((s, e, j))
    if not peak_coords_by_chr:
        raise RuntimeError("Fallback gene activity: var_names are not interval-like.")
    for chrom in peak_coords_by_chr:
        peak_coords_by_chr[chrom].sort(key=lambda x: x[0])
    genes_by_chr = _load_gene_intervals(gene_gtf, flank_up=flank_up, flank_down=flank_down)
    if not genes_by_chr:
        raise RuntimeError("Fallback gene activity: no gene intervals loaded.")
    g2p = _map_peaks_to_genes(peak_coords_by_chr, genes_by_chr)
    if not g2p:
        raise RuntimeError("Fallback gene activity: no overlapping peaks and genes.")
    gene_ids = sorted(g2p.keys())
    gene_pos = {g:i for i,g in enumerate(gene_ids)}
    rows, cols, vals = [], [], []
    for g, peaks in g2p.items():
        gi = gene_pos[g]
        for pj in peaks:
            rows.append(gi); cols.append(pj); vals.append(1.0)
    from scipy.sparse import csr_matrix
    P = csr_matrix((vals, (rows, cols)), shape=(len(gene_ids), ad.n_vars), dtype=ad.X.dtype)
    X = ad.X
    X_gene = X.dot(P.T) if sparse.issparse(X) else np.asarray(X).dot(P.T.toarray())
    ga = an.AnnData(
        X=X_gene,
        obs=ad.obs.copy(),
        var=pd.DataFrame(index=pd.Index(gene_ids, name="gene_id")),
        obsm=ad.obsm.copy() if hasattr(ad, "obsm") else None,
        uns=ad.uns.copy() if isinstance(ad.uns, dict) else {},
    )
    return ga

def _obs_numeric_depth_series(adata):
    cands = [
        "log10_Fragments", "log10_fragments", "log10_depth",
        "n_fragment", "n_fragments", "nFrags", "n_frag", "Fragments", "fragments", "depth"
    ]
    for nm in cands:
        if nm in adata.obs:
            vals = adata.obs[nm]
            try:
                arr = np.asarray(vals, dtype=float)
            except Exception:
                try:
                    arr = np.asarray(vals.astype(float))
                except Exception:
                    continue
            if "log10" not in nm.lower():
                arr = np.log10(np.clip(arr, 1.0, None))
                lbl = "log10 Fragments"
            else:
                lbl = "log10 Fragments"
            return arr, lbl
    return None, None

def plot_umap_by_depth(adata, basename: str, figdir: str) -> None:
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    vals, label = _obs_numeric_depth_series(adata)
    if vals is None:
        raise RuntimeError("No depth-like column found in .obs (e.g., n_fragment).")
    fig, ax = plt.subplots(figsize=(6, 6))
    sc = ax.scatter(X[:, 0], X[:, 1], c=vals, s=5, alpha=0.65, linewidths=0)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Cell Clustering (By Depth)")
    ax.set_aspect("equal", adjustable="box")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(label)
    os.makedirs(figdir, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 1])
    for ext in ("png", "svg"):
        out = os.path.join(figdir, f"{basename}.{ext}")
        fig.savefig(out, dpi=300, facecolor="white")
        print(f"[{ts()}] Saved figure: {out}")
    plt.close(fig)

def _apply_wrapper_injections_and_defaults() -> None:
    """No wrapper expected; keep as-is to allow external injection if present."""
    # Nothing to adjust since user provided direct paths above.
    pass

def main():
    _apply_wrapper_injections_and_defaults()
    print(HEADER)
    print(desc_txt)
    print(f"[{ts()}] Start: clustering + annotation")
    ensure_dirs()

    # Load dataset
    if not H5ADS_INPUT or not os.path.exists(H5ADS_INPUT):
        raise FileNotFoundError(f"Pre-merged dataset not found: {H5ADS_INPUT}")
    data = snap.read_dataset(H5ADS_INPUT)
    print(f"[{ts()}] Using pre-merged dataset: {H5ADS_INPUT}")

    # Pipeline
    snap.pp.select_features(data, n_features=min(N_FEATURES, data.n_vars))
    snap.tl.spectral(data)
    snap.pp.harmony(data, batch="sample", use_dims=USE_DIMS, max_iter_harmony=20)
    n_neighbors_eff = min(N_NEIGHBORS, max(2, data.n_obs - 1))
    snap.pp.knn(data, n_neighbors=n_neighbors_eff, use_dims=USE_DIMS, use_rep="X_spectral_harmony")

    # ---- NumPy 2.0 compat patch BEFORE leiden ----
    _patch_numpy_for_leiden()
    snap.tl.leiden(data)
    print(f"[{ts()}] Leiden done. (n_neighbors={n_neighbors_eff})")

    # UMAP + plots
    if SAVE_UMAP:
        snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=USE_DIMS)
        try:
            plot_umap_matplotlib(data, color="sample", basename="umap_by_sample", figdir=FIG_DIR)
        except Exception as e:
            print(f"[{ts()}] Warning: UMAP by sample failed: {e}")
        try:
            plot_umap_matplotlib(data, color="leiden", basename="umap_by_leiden", figdir=FIG_DIR)
        except Exception as e:
            print(f"[{ts()}] Warning: UMAP by leiden failed: {e}")
        try:
            plot_umap_by_depth(data, basename="umap_by_depth", figdir=FIG_DIR)
        except Exception as e:
            print(f"[{ts()}] Warning: UMAP by depth failed: {e}")

    # Gene activity (optional)
    gene_matrix = None
    if MAKE_GENE_ACTIVITY:
        try:
            try:
                chrom_sizes_dict = None
                if CHROM_SIZES_PATH and os.path.exists(CHROM_SIZES_PATH):
                    chrom_sizes_dict = {}
                    with open(CHROM_SIZES_PATH, "r") as f:
                        for line in f:
                            if not line.strip() or line.startswith(("#","track")):
                                continue
                            t = line.split()
                            if len(t) >= 2:
                                chrom = re.sub(r'^(chr|Chr)','',t[0]).lower()
                                chrom = re.sub(r'^NC_[0-9.]+','',chrom)
                                chrom_sizes_dict[chrom] = int(t[1])
            except Exception:
                chrom_sizes_dict = None
            try:
                gene_matrix = make_gene_matrix_geneid_only(data, gtf_path=GTF_PATH, chrom_sizes=chrom_sizes_dict)
                print(f"[{ts()}] Gene activity via snapatac2.make_gene_matrix OK")
            except Exception as e1:
                print(f"[{ts()}] make_gene_matrix failed; switching to fallback: {e1}")
                auto_gtf = os.path.join(os.path.dirname(GTF_PATH), "genes.gene_only.auto.gtf")
                if not os.path.exists(auto_gtf):
                    auto_gtf = build_gene_only_gtf_from_any(GTF_PATH, auto_gtf)
                gene_matrix = build_gene_activity_fallback(data, auto_gtf, flank_up=2000, flank_down=0)
                print(f"[{ts()}] Gene activity via fallback overlap OK")
            try:
                if "X_umap" in data.obsm:
                    gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
            except Exception:
                pass
            try:
                gene_matrix.write(ANNOT_GENE_H5AD)
                print(f"[{ts()}] Saved gene activity AnnData: {ANNOT_GENE_H5AD}")
            except Exception as e:
                print(f"[{ts()}] Warning: failed to save gene activity h5ad: {e}")
        except Exception as e:
            print(f"[{ts()}] Gene activity skipped (failed): {e}")
            gene_matrix = None

    # Summary
    lines = [
        HEADER, "", "===Annotation Summary===",
        f"Generated at     : {ts()}", "",
        f"Input (pre-merged) AnnDataSet: {H5ADS_INPUT}",
        f"Annotated h5ad   : {ANNOT_H5AD}",
        f"Gene activity    : {ANNOT_GENE_H5AD}",
        f"CSV output       : {CSV_OUT}",
        f"Figures dir      : {FIG_DIR}", "", "Genome:",
        f"  - chrom.sizes  : {CHROM_SIZES_PATH}",
        f"  - GTF          : {GTF_PATH}", "", "Parameters:",
        f"  - USE_DIMS     : {'all' if USE_DIMS is None else USE_DIMS}",
        f"  - N_NEIGHBORS  : {n_neighbors_eff}",
        f"  - N_FEATURES   : {min(N_FEATURES, data.n_vars)}",
        f"  - REF_CSV      : {REF_CSV_PATH}",
        f"  - LABEL_KEY    : {REF_LABEL_KEY}", "",
        f"Cells (n_obs)    : {data.n_obs}",
        f"Features (n_vars): {data.n_vars}",
        f"Clusters (Leiden): {len(set(map(str, data.obs['leiden']))) if 'leiden' in data.obs else 0}", "",
    ]
    try:
        cn_counts = data.obs["clusterName"].value_counts().sort_index()
        for s, cnt in cn_counts.items():
            lines.append(f"  - {s}: {int(cnt)}")
    except Exception:
        lines.append("  (unavailable)")
    write_summary(lines)
    print(f"[{ts()}] All done.")

if __name__ == "__main__":
    main()
