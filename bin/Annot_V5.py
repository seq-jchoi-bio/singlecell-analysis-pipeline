#!/usr/bin/env python3

import os
import re
import glob
from datetime import datetime
from typing import Dict, Optional, List, Tuple
from collections import defaultdict

# Environment (thread caps)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"]      = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"
os.environ["RAYON_NUM_THREADS"]    = "1"
os.environ["MALLOC_ARENA_MAX"]     = "2"

# Plot backends
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
mpl.rcParams["savefig.facecolor"] = "white"
import matplotlib.pyplot as plt

# Core libs
import numpy as np
import pandas as pd
import anndata as an
import snapatac2 as snap
from snapatac2 import AnnDataSet as _SnapAnnDataSet
from scipy import sparse

import warnings
warnings.filterwarnings("ignore", message=r"`_import_from_c` is deprecated", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden", category=UserWarning, module=r"umap\.umap_")

__version__ = "2.9"

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq: Clustering + Annotation Pipeline"
HEADER += f"\n      Version {__version__}"
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

desc_txt = """This script:
  (1) loads a pre-merged .h5ads dataset
  (2) builds a gene activity matrix (gene_id-based) with robust GTF handling
  (3) parses a reference CSV of marker genes -> per-cell-type signatures
  (4) transfers labels via AUC or weighted scoring
  (5) saves annotated .h5ad, CSV, UMAP (Leiden-majority by clusterName), and
      BED files per (Leiden × sample) under AGGR_DIR/bed as "<Lxx>.<Sample>.bed"
"""

# =========================
# User-configurable section (wrapper-friendly)
# =========================
BASE_DIR            = None   # injected by wrapper (-i)
H5ADS_INPUT         = None   # default: BASE_DIR/Filter_results/merged_doublets.h5ads
AGGR_DIR            = None   # default: CWD/Annot_results
FIG_DIR             = None   # default: AGGR_DIR/Plots

H5ADS_OUT           = None   # default: AGGR_DIR/annot_merged.h5ads
ANNOT_H5AD          = None   # default: AGGR_DIR/annot_merged_cells.h5ad
ANNOT_GENE_H5AD     = None   # default: AGGR_DIR/annot_gene_activity.h5ad
SUMMARY_TXT         = None   # default: AGGR_DIR/summary.txt
CSV_OUT             = None   # default: AGGR_DIR/annotation_results.csv

CHROM_SIZES_PATH    = None   # injected by wrapper (-s/-cs)
GTF_PATH            = None   # injected by wrapper (-s/-gtf)
REF_CSV_PATH        = None   # auto-detected from refGenome or overridden by --ref-csv

# Defaults (added)
REF_LABEL_KEY       = "celltype_id"
MARKER_MIN_LOG2FC   = 0.5
MARKER_MAX_PADJ     = 0.05
MIN_MARKERS_PER_SET = 5
N_FEATURES          = 250_000
N_NEIGHBORS         = 10
USE_DIMS            = None            # None => use all spectral dims; or set int/list
MAKE_GENE_ACTIVITY  = True            # required for CSV marker-based transfer
TRANSFER_MODE       = "auc"           # "auc" or "weighted"
UNKNOWN_CONF_THRESH = 0.5
SAVE_UMAP           = True
LABEL_MIN_COUNT     = 20
LABEL_OUTLINE       = True

# =========================
# Utilities
# =========================
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dirs() -> None:
    os.makedirs(AGGR_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

def write_summary(lines: List[str]) -> None:
    with open(SUMMARY_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[{ts()}] Summary written: {SUMMARY_TXT}")

def _patch_numpy_for_leiden() -> None:
    if not hasattr(np, "compat"):
        class _Compat: pass
        np.compat = _Compat()
    if not hasattr(np.compat, "unicode"):
        np.compat.unicode = str

# Plot
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
        ax.scatter(X[m, 0], X[m, 1],
                   s=5, alpha=0.5, color=colors[i], linewidths=0, label=str(cat))

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

def _compute_label_centers(X: np.ndarray, labels: np.ndarray, min_count: int = LABEL_MIN_COUNT):
    centers = {}
    for lab in np.unique(labels):
        mask = (labels == lab)
        n = int(mask.sum())
        if n < max(1, min_count):
            continue
        x_med = float(np.median(X[mask, 0]))
        y_med = float(np.median(X[mask, 1]))
        centers[str(lab)] = (x_med, y_med, n)
    return centers

def _draw_text_labels(ax, centers: Dict[str, tuple]):
    for lab, (x, y, n) in centers.items():
        if LABEL_OUTLINE:
            for lw in (3.5, 1.5):
                ax.text(x, y, lab, ha="center", va="center",
                        fontsize=10, weight="bold",
                        color="white", zorder=5, path_effects=None)
        ax.text(x, y, lab, ha="center", va="center",
                fontsize=10, weight="bold", color="black", zorder=6)

def plot_umap_matplotlib(adata, color: str, basename: str, figdir: str) -> None:
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    col = _obs_as_str_array(adata, color)
    _scatter_by_category(X, col, "UMAP", basename, figdir)

def plot_umap_with_text(adata, label_key: str, basename: str, figdir: str) -> None:
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    labels = _obs_as_str_array(adata, label_key)
    cats = np.array(sorted(np.unique(labels), key=str))
    cmap = plt.get_cmap("tab20", len(cats))
    colors = cmap(range(len(cats)))

    fig, ax = plt.subplots(figsize=(6, 6))
    for i, cat in enumerate(cats):
        m = (labels == cat)
        ax.scatter(X[m, 0], X[m, 1],
                   s=6, alpha=0.6, color=colors[i], linewidths=0, label=str(cat))

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP by {label_key} (text labels)")
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

# Chromosome name normalization
def _normalize_chrom(s: str) -> str:
    s0 = str(s).strip()
    s0 = re.sub(r'^(chr|Chr)', '', s0)
    s0 = re.sub(r'^NC_[0-9.]+', '', s0)
    return s0.lower()

# GTF processing
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
                    if attrs and attrs[-1] != ";": attrs += ";"
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
        raise RuntimeError("Failed to build gene-only GTF: no usable gene_id found in transcript/exon lines.")

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

# Gene activity matrix
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

# Reference CSV -> signatures/weights
def load_reference_markers(csv_path: str,
                           label_key: str = REF_LABEL_KEY,
                           min_log2fc: float = MARKER_MIN_LOG2FC,
                           max_padj: float = MARKER_MAX_PADJ,
                           min_markers: int = MIN_MARKERS_PER_SET) -> Dict[str, List[str]]:
    df = pd.read_csv(csv_path, sep=None, engine="python")
    key = label_key if label_key in df.columns else ("clusterName" if "clusterName" in df.columns else None)
    if key is None:
        raise ValueError("Neither 'celltype_id' nor 'clusterName' is in the reference CSV.")
    if "avg_log2FC" not in df.columns: df["avg_log2FC"] = 0.0
    if "p_val_adj" not in df.columns:  df["p_val_adj"]  = 0.0
    m = (df["avg_log2FC"] >= float(min_log2fc)) & (df["p_val_adj"] <= float(max_padj))
    df_f = df.loc[m].copy()
    sigs: Dict[str, List[str]] = {}
    for lab, sub in df_f.groupby(key):
        genes = list(map(str, sub["gene"].astype(str)))
        genes = [g for g in genes if len(g) > 0]
        if len(genes) >= int(min_markers):
            sigs[str(lab)] = sorted(set(genes))
    if not sigs:
        raise ValueError("No signature built after filtering; relax thresholds or check CSV.")
    return sigs

def load_reference_weights(csv_path: str,
                           label_key: str = REF_LABEL_KEY,
                           min_log2fc: float = MARKER_MIN_LOG2FC,
                           max_padj: float = MARKER_MAX_PADJ,
                           min_markers: int = MIN_MARKERS_PER_SET) -> Dict[str, Dict[str, float]]:
    df = pd.read_csv(csv_path, sep=None, engine="python")
    key = label_key if label_key in df.columns else ("clusterName" if "clusterName" in df.columns else None)
    if key is None:
        raise ValueError("Neither 'celltype_id' nor 'clusterName' is in the reference CSV.")
    if "avg_log2FC" not in df.columns: df["avg_log2FC"] = 0.0
    if "p_val_adj" not in df.columns:  df["p_val_adj"]  = 1.0
    m = (df["avg_log2FC"] >= float(min_log2fc)) & (df["p_val_adj"] <= float(max_padj))
    df_f = df.loc[m].copy()
    df_f["weight"] = df_f["avg_log2FC"] * (-np.log10(df_f["p_val_adj"].astype(float) + 1e-300))
    weights: Dict[str, Dict[str, float]] = {}
    for lab, sub in df_f.groupby(key):
        w = {str(r["gene"]): float(r["weight"]) for _, r in sub.iterrows()}
        w = {g: wt for g, wt in w.items() if wt > 0}
        if len(w) >= int(min_markers):
            weights[str(lab)] = w
    if not weights:
        raise ValueError("No weighted signature built; relax thresholds or check CSV.")
    return weights

# Scoring logic
def _ensure_dense(X) -> np.ndarray:
    if sparse.issparse(X): return X.toarray()
    return np.asarray(X)

def score_auc_rank(ga: an.AnnData, signatures: Dict[str, List[str]]) -> pd.DataFrame:
    gene_index = {g: i for i, g in enumerate(map(str, ga.var_names))}
    sig_idx: Dict[str, np.ndarray] = {}
    for lab, gene_list in signatures.items():
        idx = [gene_index[g] for g in gene_list if g in gene_index]
        if len(idx) >= 1:
            sig_idx[lab] = np.array(idx, dtype=int)
    if not sig_idx:
        raise RuntimeError("No signature genes overlap with gene_activity var_names.")
    X = ga.X
    if sparse.issparse(X): X = X.tocsr()
    n_cells = ga.n_obs
    labels = list(sig_idx.keys())
    out = np.zeros((n_cells, len(labels)), dtype=float)
    for i in range(n_cells):
        row = X.getrow(i).toarray().ravel() if sparse.issparse(X) else np.asarray(X[i, :]).ravel()
        order = np.argsort(row)
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, row.size + 1, dtype=float)
        ranks_norm = ranks / float(row.size)
        for j, lab in enumerate(labels):
            idx = sig_idx[lab]
            r_mean = ranks_norm[idx].mean() if idx.size > 0 else 1.0
            out[i, j] = 1.0 - r_mean
    out_min = out.min(axis=0, keepdims=True)
    out_max = out.max(axis=0, keepdims=True)
    denom = np.clip(out_max - out_min, 1e-9, None)
    out = (out - out_min) / denom
    return pd.DataFrame(out, index=ga.obs_names, columns=labels)

def score_weighted_sum(ga: an.AnnData, weights: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    gene_index = {g: i for i, g in enumerate(map(str, ga.var_names))}
    label_genes = {lab: [gene_index[g] for g in wg if g in gene_index] for lab, wg in weights.items()}
    label_w     = {lab: np.array([weights[lab][g] for g in weights[lab] if g in gene_index], dtype=float)
                   for lab in weights.keys()}
    label_genes = {lab: idx for lab, idx in label_genes.items() if len(idx) > 0}
    label_w     = {lab: label_w[lab] for lab in label_genes.keys()}
    if not label_genes:
        raise RuntimeError("No overlap between weighted markers and gene_activity var_names.")
    X = _ensure_dense(ga.X)
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True)
    sd = np.where(sd < 1e-9, 1.0, sd)
    Z = (X - mu) / sd
    cells = ga.n_obs
    labels = list(label_genes.keys())
    out = np.zeros((cells, len(labels)), dtype=float)
    for j, lab in enumerate(labels):
        idx = np.array(label_genes[lab], dtype=int)
        w   = label_w[lab]
        z_sub = Z[:, idx]
        if w.size != z_sub.shape[1]:
            genes_lab = [ga.var_names[i] for i in idx]
            w_lookup = weights[lab]
            w = np.array([w_lookup[str(g)] for g in genes_lab], dtype=float)
        out[:, j] = np.maximum(0.0, (z_sub * w).sum(axis=1))
    out_min = out.min(axis=0, keepdims=True)
    out_max = out.max(axis=0, keepdims=True)
    denom = np.clip(out_max - out_min, 1e-9, None)
    out = (out - out_min) / denom
    return pd.DataFrame(out, index=ga.obs_names, columns=labels)

def pick_labels(score_df: pd.DataFrame, unknown_thresh: float = UNKNOWN_CONF_THRESH) -> Tuple[pd.Series, pd.Series]:
    max_lab = score_df.idxmax(axis=1)
    max_val = score_df.max(axis=1)
    pred = max_lab.copy()
    pred[max_val < float(unknown_thresh)] = "unknown"
    return pred, max_val

# Leiden majority label
def add_leiden_majority_label(ad_like, leiden_key: str = "leiden", label_key: str = "pred_cell_type",
                              out_key: str = "leiden_majority_celltype") -> None:
    try:
        leiden = np.asarray(ad_like.obs[leiden_key]).astype(str)
        label  = np.asarray(ad_like.obs[label_key]).astype(str)
    except Exception:
        return
    if leiden.size != ad_like.n_obs or label.size != ad_like.n_obs:
        return
    df = pd.DataFrame({"leiden": leiden, "label": label})
    maj = df.groupby("leiden")["label"].agg(lambda s: s.value_counts().idxmax()).to_dict()
    mapped = np.array([maj.get(l, "unknown") for l in leiden], dtype=object)
    ad_like.obs[out_key] = mapped

# Mapping: predicted label -> clusterName
def build_clustername_mapping(csv_path: str, ref_label_key: str) -> Dict[str, str]:
    df = pd.read_csv(csv_path, sep=None, engine="python")
    if "clusterName" not in df.columns:
        return {}
    key = ref_label_key if ref_label_key in df.columns else ("clusterName" if "clusterName" in df.columns else None)
    if key is None:
        return {}
    mapping: Dict[str, str] = {}
    for k, sub in df.groupby(key):
        cn = sub["clusterName"].dropna().astype(str)
        if len(cn) == 0:
            continue
        mapping[str(k)] = cn.value_counts().idxmax()
    return mapping

# Write BED per leiden and sample
def _format_leiden_tag(val) -> str:
    try:
        ival = int(str(val))
        return f"L{ival:02d}"
    except Exception:
        return f"L{str(val)}"

def _parse_all_peaks(var_names):
    parsed = []
    for j, v in enumerate(map(str, var_names)):
        p = _parse_peak_name(v)
        if p is None:
            continue
        chrom, s, e = p
        parsed.append((chrom, int(s), int(e), v))
    return parsed

def _cluster_purity(labels: np.ndarray, majority: str) -> float:
    if labels.size == 0:
        return 0.0
    return float((labels == majority).sum()) / float(labels.size)

def write_bed_per_leiden_and_sample(ad_cells: an.AnnData,
                                    out_dir: str,
                                    leiden_key: str = "leiden",
                                    majority_key: str = "leiden_majority_clusterName",
                                    sample_key: str = "sample",
                                    max_files: int = 10000):

    os.makedirs(out_dir, exist_ok=True)
    if leiden_key not in ad_cells.obs or sample_key not in ad_cells.obs:
        print(f"[{ts()}] Skip BED: missing obs keys ('{leiden_key}' or '{sample_key}')")
        return

    X = ad_cells.X
    is_sparse = sparse.issparse(X)
    if is_sparse:
        X = X.tocsr()

    peak_meta = _parse_all_peaks(ad_cells.var_names)
    if len(peak_meta) == 0:
        print(f"[{ts()}] Skip BED: var_names not interval-like.")
        return
    name_to_parsed = {pname: (chrom, start, end) for (chrom, start, end, pname) in peak_meta}
    parsed_mask = np.array([vn in name_to_parsed for vn in map(str, ad_cells.var_names)], dtype=bool)
    keep_idx = np.where(parsed_mask)[0]
    if keep_idx.size == 0:
        print(f"[{ts()}] Skip BED: no parseable peaks.")
        return

    leiden_vals = np.asarray(ad_cells.obs[leiden_key]).astype(str)
    sample_vals = np.asarray(ad_cells.obs[sample_key]).astype(str)

    if majority_key in ad_cells.obs:
        majority_vals = np.asarray(ad_cells.obs[majority_key]).astype(str)
    else:
        if "clusterName" in ad_cells.obs:
            tmp_df = pd.DataFrame({"leiden": leiden_vals,
                                   "clusterName": np.asarray(ad_cells.obs["clusterName"]).astype(str)})
            maj_map = tmp_df.groupby("leiden")["clusterName"].agg(lambda s: s.value_counts().idxmax()).to_dict()
        else:
            maj_map = {}
        majority_vals = np.array([maj_map.get(l, "unknown") for l in leiden_vals], dtype=object)

    u_leiden = np.unique(leiden_vals)
    u_sample = np.unique(sample_vals)

    files_written = 0
    for lval in u_leiden:
        Ltag = _format_leiden_tag(lval)
        m_in_cluster = (leiden_vals == lval)
        if m_in_cluster.any():
            maj_name = pd.Series(majority_vals[m_in_cluster]).mode().iloc[0]
            purity = _cluster_purity(majority_vals[m_in_cluster], maj_name)
        else:
            maj_name = "unknown"
            purity = 0.0

        for sval in u_sample:
            mask = (leiden_vals == lval) & (sample_vals == sval)
            n_cells = int(mask.sum())
            if n_cells == 0:
                continue

            if is_sparse:
                X_mean = X[mask, :][:, keep_idx].mean(axis=0).A1
            else:
                X_mean = X[mask, :][:, keep_idx].mean(axis=0)

            if X_mean.size == 0:
                continue
            lo = np.percentile(X_mean, 5.0)
            hi = np.percentile(X_mean, 95.0)
            denom = max(hi - lo, 1e-9)
            score_vec = np.clip((X_mean - lo) / denom, 0.0, 1.0)
            score_vec = np.round(score_vec * 1000).astype(int)

            rows = []
            for j_local, j in enumerate(keep_idx):
                raw_name = str(ad_cells.var_names[j])
                chrom, start, end = name_to_parsed[raw_name]
                name_field = f"{Ltag}|{maj_name}|{raw_name}"
                rows.append((chrom, int(start), int(end), name_field,
                             sval, maj_name, f"{purity:.3f}", n_cells, int(score_vec[j_local])))

            fname = f"{Ltag}_{_sanitize_for_fs(maj_name)}_{_sanitize_for_fs(sval)}.bed"
            out_path = os.path.join(out_dir, fname)

            with open(out_path, "w") as fout:
                for r in rows:
                    # BED4 + extras: chrom, start, end, name, sample, majority, purity, n_cells, mean_accessibility_0_1000
                    fout.write("\t".join(map(str, r)) + "\n")

            print(f"[{ts()}] BED written: {out_path} (peaks={len(rows)}, cells={n_cells}, purity={purity:.3f})")
            files_written += 1
            if files_written >= max_files:
                print(f"[{ts()}] Reached max_files={max_files}; stopping BED export.")
                return

def _sanitize_for_fs(name: str) -> str:
    """Make a string safe for filenames: spaces -> underscores, remove unsafe chars."""
    s = re.sub(r'\s+', '_', str(name).strip())
    s = re.sub(r'[^A-Za-z0-9._-]+', '_', s)
    s = re.sub(r'_{2,}', '_', s).strip('_')
    return s or "unknown"

# =========================
# Main
# =========================
# Wrapper injection
def _apply_wrapper_injections_and_defaults() -> None:
    """Resolve wrapper-injected globals and set safe defaults for outputs/inputs."""
    global BASE_DIR, H5ADS_INPUT, AGGR_DIR, FIG_DIR, SUMMARY_TXT, CSV_OUT
    global CHROM_SIZES_PATH, GTF_PATH, REF_CSV_PATH, REF_H5AD_PATH

    g = globals()
    # Accept lowercase aliases as well
    if g.get("base_dir"): BASE_DIR = g.get("base_dir")
    if g.get("out_dir"): AGGR_DIR = g.get("out_dir")
    if g.get("chrom_sizes_path") and not g.get("CHROM_SIZES_PATH"):
        CHROM_SIZES_PATH = g.get("chrom_sizes_path")
    if g.get("annot_path") and not g.get("GTF_PATH"):
        GTF_PATH = g.get("annot_path")

    # Default outputs to CWD/Annot_results
    if not AGGR_DIR:
        AGGR_DIR = os.path.join(os.getcwd(), "Annot_results")
    os.makedirs(AGGR_DIR, exist_ok=True)

    # Set standard derived outputs if variables exist in logic
    if not FIG_DIR:
        FIG_DIR = os.path.join(AGGR_DIR, "Plots")
    os.makedirs(FIG_DIR, exist_ok=True)

    if not SUMMARY_TXT:
        SUMMARY_TXT = os.path.join(AGGR_DIR, "summary.txt")
    if not CSV_OUT:
        CSV_OUT = os.path.join(AGGR_DIR, "annotation_results.csv")

    # H5ADS_INPUT default if not provided
    if not H5ADS_INPUT and BASE_DIR:
        candidate = os.path.join(BASE_DIR, "Filter_results", "merged_doublets.h5ads")
        H5ADS_INPUT = candidate

def _apply_wrapper_injections_and_defaults() -> None:
    global BASE_DIR, H5ADS_INPUT, AGGR_DIR, FIG_DIR, H5ADS_OUT, ANNOT_H5AD, ANNOT_GENE_H5AD
    global SUMMARY_TXT, CSV_OUT, CHROM_SIZES_PATH, GTF_PATH, REF_CSV_PATH

    g = globals()
    # Accept lowercase aliases from wrappers
    if g.get("base_dir"):
        BASE_DIR = g.get("base_dir")
    if g.get("chrom_sizes_path") and not g.get("CHROM_SIZES_PATH"):
        CHROM_SIZES_PATH = g.get("chrom_sizes_path")
    if g.get("annot_path") and not g.get("GTF_PATH"):
        GTF_PATH = g.get("annot_path")
    if g.get("out_dir"):
        AGGR_DIR = g.get("out_dir")

    # Outputs under CWD/Annot_results by default
    if not AGGR_DIR:
        AGGR_DIR = os.path.join(os.getcwd(), "Annot_results")
    os.makedirs(AGGR_DIR, exist_ok=True)

    if not FIG_DIR:
        FIG_DIR = os.path.join(AGGR_DIR, "Plots")
    os.makedirs(FIG_DIR, exist_ok=True)

    # Standard output files
    if not H5ADS_OUT:
        H5ADS_OUT = os.path.join(AGGR_DIR, "annot_merged.h5ads")
    if not ANNOT_H5AD:
        ANNOT_H5AD = os.path.join(AGGR_DIR, "annot_merged_cells.h5ad")
    if not ANNOT_GENE_H5AD:
        ANNOT_GENE_H5AD = os.path.join(AGGR_DIR, "annot_gene_activity.h5ad")
    if not SUMMARY_TXT:
        SUMMARY_TXT = os.path.join(AGGR_DIR, "summary.txt")
    if not CSV_OUT:
        CSV_OUT = os.path.join(AGGR_DIR, "annotation_results.csv")

    # Input default
    if not H5ADS_INPUT and BASE_DIR:
        H5ADS_INPUT = os.path.join(BASE_DIR, "Filter_results", "merged_doublets.h5ads")

def main():
    _apply_wrapper_injections_and_defaults()
    _apply_wrapper_injections_and_defaults()
    print(HEADER)
    print(desc_txt)
    print(f"[{ts()}] Start: clustering + annotation")
    ensure_dirs()
    # 1) Load pre-merged AnnDataSet (*.h5ads)
    if not os.path.exists(H5ADS_INPUT):
        raise FileNotFoundError(f"Pre-merged dataset not found: {H5ADS_INPUT}")
    try:
        data = snap.read_dataset(H5ADS_INPUT)
        print(f"[{ts()}] Using pre-merged dataset: {H5ADS_INPUT}")
    except Exception as e:
        raise RuntimeError(f"Failed to open merged dataset: {H5ADS_INPUT} ({e})")

    # 2) Feature selection -> spectral -> Harmony -> KNN -> Leiden
    snap.pp.select_features(data, n_features=min(N_FEATURES, data.n_vars))
    snap.tl.spectral(data)
    snap.pp.harmony(data, batch="sample", use_dims=USE_DIMS, max_iter_harmony=20)
    n_neighbors_eff = min(N_NEIGHBORS, max(2, data.n_obs - 1))
    snap.pp.knn(data, n_neighbors=n_neighbors_eff, use_dims=USE_DIMS, use_rep="X_spectral_harmony")
    _patch_numpy_for_leiden()
    snap.tl.leiden(data)
    print(f"[{ts()}] Leiden done. (n_neighbors={n_neighbors_eff})")

    # 3) UMAP
    if SAVE_UMAP:
        snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=USE_DIMS)
        try:
            plot_umap_matplotlib(data, color="sample", basename="umap_by_sample", figdir=FIG_DIR)
        except Exception as e: print(f"[{ts()}] Warning: UMAP by sample failed: {e}")
        try:
            plot_umap_matplotlib(data, color="leiden", basename="umap_by_leiden", figdir=FIG_DIR)
        except Exception as e: print(f"[{ts()}] Warning: UMAP by leiden failed: {e}")

    # 4) Gene activity (robust+fallback)
    gene_matrix = None
    if MAKE_GENE_ACTIVITY:
        try:
            try:
                chrom_sizes_dict = load_chrom_sizes(CHROM_SIZES_PATH)
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
                if "X_umap" in data.obsm: gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
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
    else:
        print(f"[{ts()}] Skipped gene activity by configuration.")

    # 5) Load reference CSV -> signatures -> scoring
    if gene_matrix is None:
        raise RuntimeError("Gene activity matrix is required for CSV marker-based annotation.")
    if TRANSFER_MODE == "auc":
        signatures = load_reference_markers(REF_CSV_PATH, REF_LABEL_KEY,
                                            MARKER_MIN_LOG2FC, MARKER_MAX_PADJ, MIN_MARKERS_PER_SET)
        score_df = score_auc_rank(gene_matrix, signatures)
    elif TRANSFER_MODE == "weighted":
        weights = load_reference_weights(REF_CSV_PATH, REF_LABEL_KEY,
                                         MARKER_MIN_LOG2FC, MARKER_MAX_PADJ, MIN_MARKERS_PER_SET)
        score_df = score_weighted_sum(gene_matrix, weights)
    else:
        raise ValueError("TRANSFER_MODE must be 'auc' or 'weighted'.")

    pred_label, conf = pick_labels(score_df, UNKNOWN_CONF_THRESH)

    # 6) Attach annotation and clusterName mapping
    labels_aligned = (
        pd.Series(pred_label, index=gene_matrix.obs_names)
          .reindex(data.obs_names)
          .fillna("unknown")
          .astype(str)
          .to_numpy()
    )
    conf_aligned = (
        pd.Series(conf, index=gene_matrix.obs_names)
          .reindex(data.obs_names)
          .fillna(0.0)
          .astype(float)
          .to_numpy()
    )
    data.obs["pred_cell_type"]  = labels_aligned
    data.obs["pred_confidence"] = conf_aligned

    lbl2cn = build_clustername_mapping(REF_CSV_PATH, REF_LABEL_KEY)
    if not lbl2cn:
        lbl2cn = {s: s for s in np.unique(labels_aligned)}
    cluster_arr = np.array([lbl2cn.get(str(s), str(s)) if str(s) != "unknown" else "unknown"
                            for s in labels_aligned], dtype=object)
    data.obs["clusterName"] = cluster_arr

    if "leiden" in data.obs:
        try:
            data.obs["leiden"] = [str(x) for x in data.obs["leiden"]]
        except Exception:
            pass
        add_leiden_majority_label(data, "leiden", "pred_cell_type", "leiden_majority_celltype")
        add_leiden_majority_label(data, "leiden", "clusterName",    "leiden_majority_clusterName")

    # 7) Save annotated per-cell h5ad
    try:
        ad_cells = data.to_adata()
        ad_cells.write(ANNOT_H5AD)
        print(f"[{ts()}] Saved annotated cells h5ad: {ANNOT_H5AD}")
    except Exception as e:
        print(f"[{ts()}] Warning: failed to save annotated per-cell h5ad: {e}")
        ad_cells = data.to_adata()

    # 8) Plots
    if SAVE_UMAP:
        if "leiden_majority_clusterName" in data.obs:
            try:
                plot_umap_with_text(data, label_key="leiden_majority_clusterName",
                                    basename="umap_by_leiden_majority_clusterName_labeled", figdir=FIG_DIR)
            except Exception as e:
                print(f"[{ts()}] Warning: UMAP with text labels (leiden_majority_clusterName) failed: {e}")

    # 9) Export CSV
    df_obs = ad_cells.obs.copy()
    if "X_umap" in ad_cells.obsm:
        umap = np.asarray(ad_cells.obsm["X_umap"])
        df_obs["umap1"] = umap[:, 0]
        df_obs["umap2"] = umap[:, 1]
    cols = [c for c in [
        "sample", "leiden",
        "pred_cell_type", "pred_confidence",
        "clusterName", "leiden_majority_clusterName",
        "leiden_majority_celltype",
        "umap1", "umap2"
    ] if c in df_obs.columns]
    out_df = df_obs[cols].copy()
    out_df.insert(0, "cell_id", df_obs.index.astype(str))
    out_df.to_csv(CSV_OUT, index=False)
    print(f"[{ts()}] Saved CSV: {CSV_OUT}")

    # 10) BED export
    try:
        bed_out_dir = os.path.join(AGGR_DIR, "bed")
        os.makedirs(bed_out_dir, exist_ok=True)
        write_bed_per_leiden_and_sample(ad_cells, bed_out_dir,
                                        leiden_key="leiden",
                                        majority_key="leiden_majority_clusterName",
                                        sample_key="sample")
    except Exception as e:
        print(f"[{ts()}] Warning: BED export failed: {e}")

    # 11) Summary
    lines = [
        HEADER,
        "",
        "===Annotation Summary===",
        f"Generated at     : {ts()}",
        "",
        f"Input (pre-merged) AnnDataSet: {H5ADS_INPUT}",
        f"Annotated h5ad   : {ANNOT_H5AD}",
        f"Gene activity    : {ANNOT_GENE_H5AD}",
        f"CSV output       : {CSV_OUT}",
        f"Figures dir      : {FIG_DIR}",
        f"BED out dir      : {os.path.join(AGGR_DIR, 'bed')}",
        "",
        "Genome:",
        f"  - chrom.sizes  : {CHROM_SIZES_PATH}",
        f"  - GTF          : {GTF_PATH}",
        "",
        "Parameters:",
        f"",
        f"  - USE_DIMS     : {'all' if USE_DIMS is None else USE_DIMS}",
        f"  - N_NEIGHBORS  : {n_neighbors_eff}",
        f"  - N_FEATURES   : {min(N_FEATURES, data.n_vars)}",
        f"  - TRANSFER     : {TRANSFER_MODE} (unknown_thresh={UNKNOWN_CONF_THRESH})",
        f"  - REF_CSV      : {REF_CSV_PATH}",
        f"  - LABEL_KEY    : {REF_LABEL_KEY}",
        f"  - Marker filt  : log2FC≥{MARKER_MIN_LOG2FC}, padj≤{MARKER_MAX_PADJ}, min_n={MIN_MARKERS_PER_SET}",
        "",
        f"Cells (n_obs)    : {data.n_obs}",
        f"Features (n_vars): {data.n_vars}",
        f"Clusters (Leiden): {len(set(map(str, data.obs['leiden']))) if 'leiden' in data.obs else 0}",
        "",
        "Predicted cell types (counts):",
    ]
    try:
        ct_counts = ad_cells.obs["pred_cell_type"].value_counts().sort_index()
        for s, cnt in ct_counts.items():
            lines.append(f"  - {s}: {int(cnt)}")
        lines.append("")
        lines.append("clusterName (counts):")
        cn_counts = ad_cells.obs["clusterName"].value_counts().sort_index()
        for s, cnt in cn_counts.items():
            lines.append(f"  - {s}: {int(cnt)}")
    except Exception:
        lines.append("  (unavailable)")
    write_summary(lines)
    print(f"[{ts()}] All done.")

if __name__ == "__main__":
    main()

