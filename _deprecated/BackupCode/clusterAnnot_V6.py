

# Default parallel jobs (figures excluded from parallelization)
NJOB = 8
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
from joblib import Parallel, delayed

# Core libs
import numpy as np

def _obs_has(_adata, _key):
    _obs = getattr(_adata, "obs", None)
    if _obs is None:
        return False
    # dict-like access
    try:
        _ = _obs[_key]
        return True
    except Exception:
        pass
    # attribute-style
    try:
        _ = getattr(_obs, _key)
        return True
    except Exception:
        pass
    # pandas-like columns attribute if present (best-effort only)
    _cols = getattr(_obs, "columns", None)
    if _cols is not None:
        try:
            return _key in _cols
        except Exception:
            return False
    return False

def _obs_get(_adata, _key):
    _obs = getattr(_adata, "obs", None)
    if _obs is None:
        raise KeyError("AnnData has no .obs")
    try:
        _col = _obs[_key]
    except Exception:
        _col = getattr(_obs, _key)  # may still raise
    try:
        import numpy as _np
        return _np.asarray(_col)
    except Exception:
        return [v for v in _col]
import pandas as pd
import anndata as an
import snapatac2 as snap
from scipy import sparse

from matplotlib import MatplotlibDeprecationWarning
import warnings
warnings.filterwarnings("ignore", message=r"`_import_from_c` is deprecated", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden", category=UserWarning, module=r"umap\.umap_")
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning, message=r".*get_cmap function was deprecated.*")
warnings.filterwarnings("ignore",message=r".*`np\.object`.*", category=FutureWarning)

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
REF_LABEL_KEY       = "clusterName"
N_FEATURES          = 250_000
N_NEIGHBORS         = 10
USE_DIMS            = None
MAKE_GENE_ACTIVITY  = True
SAVE_UMAP           = True


# --- Threshold Defaults for label assignment ---
MIN_SCORE  = 0.50   # clusters with best_score below this become 'unknown'
MIN_MARGIN = 0.25   # if USE_MARGIN is True: clusters with margin below this become 'unknown'
USE_MARGIN = False   # set False to ignore margin rule and use score-only

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
    # Parallel over chromosomes when NJOB > 1
    if int(NJOB) <= 1:
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

    def _map_one_chrom(task):
        chrom, peaks = task
        genes = gene_intervals_by_chr.get(chrom, [])
        local: Dict[str, set] = defaultdict(set)
        if not genes or not peaks:
            return local
        gi = 0
        for (p_start, p_end, p_idx) in peaks:
            center = (p_start + p_end) // 2
            while gi < len(genes) and genes[gi][1] < center:
                gi += 1
            for k in (gi-2, gi-1, gi, gi+1, gi+2):
                if 0 <= k < len(genes):
                    g_start, g_end, gid = genes[k]
                    if g_start <= center <= g_end:
                        local[gid].add(p_idx)
        return local

    parts = Parallel(n_jobs=int(NJOB), backend="loky", prefer="processes")(
        delayed(_map_one_chrom)(item) for item in peak_coords_by_chr.items()
    )
    merged: Dict[str, set] = defaultdict(set)
    for d in parts:
        for gid, idxs in d.items():
            merged[gid].update(idxs)
    return merged

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



# ---------------- Heatmap & Scoring (reference markers × Leiden) ----------------
# (Dual exports: robust-z cluster-median and pseudobulk mean logCPM)

from typing import List, Tuple

def _load_reference_markers_with_labels(csv_path: str, label_col: str) -> List[Tuple[str, str]]:
    import csv, os
    out: List[Tuple[str,str]] = []
    if not csv_path or not os.path.exists(csv_path):
        return out
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return out
        lower = [c.lower() for c in reader.fieldnames]
        gene_cand = None
        for name in ["gene", "symbol", "gene_symbol", "gene_name", "genes", "marker", "markers"]:
            if name in lower:
                gene_cand = reader.fieldnames[lower.index(name)]
                break
        if gene_cand is None:
            gene_cand = reader.fieldnames[0]
        lab_col = label_col if label_col in reader.fieldnames else None
        if lab_col is None:
            for nm in ["clusterName","cluster","celltype","celltype_id","label","group"]:
                if nm in reader.fieldnames:
                    lab_col = nm
                    break
        if lab_col is None:
            lab_col = reader.fieldnames[-1]
        seen = set()
        for row in reader:
            g = str(row.get(gene_cand, "")).strip()
            lab = str(row.get(lab_col, "")).strip()
            if not g:
                continue
            key = (g, lab)
            if key in seen:
                continue
            seen.add(key)
            out.append((g, lab))
    return out

def _robust_z_per_gene_dense(X):
    import numpy as np
    med = np.median(X, axis=0)
    mad = np.median(np.abs(X - med), axis=0)
    mad = 1.4826 * mad
    mad[mad == 0] = 1.0
    Z = (X - med) / mad
    return Z

def _to_dense_matrix(mat):
    try:
        import scipy.sparse as sp
        if sp.issparse(mat):
            return mat.tocsr().astype(float).toarray()
    except Exception:
        pass
    import numpy as np
    return np.asarray(mat, dtype=float)

def _align_obs(an_a, an_b):
    A = an_a.to_adata() if hasattr(an_a, "to_adata") else an_a
    B = an_b.to_adata() if hasattr(an_b, "to_adata") else an_b
    try:
        A.obs_names_make_unique(); B.obs_names_make_unique()
    except Exception:
        pass
    common = A.obs_names.intersection(B.obs_names)
    if len(common) == 0:
        raise RuntimeError("No overlapping cell barcodes between objects.")
    return A[common].copy(), B[common].copy()

def _auroc_pos_vs_neg(pos_vals, neg_vals):
    import numpy as np
    from scipy.stats import rankdata
    x = np.concatenate([pos_vals, neg_vals])
    r = rankdata(x, method="average")
    n1 = len(pos_vals); n2 = len(neg_vals)
    R1 = r[:n1].sum()
    auc = (R1 - n1*(n1+1)/2) / (n1*n2) if n1>0 and n2>0 else float("nan")
    return float(auc)

def _cluster_median_by_leiden(z_mat, obs_leiden, var_names):
    import pandas as pd
    labs = list(map(str, obs_leiden))
    df = pd.DataFrame(z_mat, columns=list(map(str, var_names)))
    df["__lab__"] = labs
    med = df.groupby("__lab__", sort=True).median().sort_index(key=lambda idx: [int(x) if str(x).isdigit() else 9_999_999 for x in idx])
    med = med.drop(columns=["__lab__"], errors="ignore")
    return med

def _build_marker_groups(marker_pairs: List[Tuple[str,str]]):
    """
    Build ordered gene list and (label, start, end) groups
    with labels sorted alphabetically (e.g., A, B, C),
    preserving per-label gene order as in marker_pairs.
    """
    from collections import defaultdict
    # Group genes by clusterName (label)
    buckets = defaultdict(list)
    for g, lab in marker_pairs:
        buckets[str(lab)].append(g)
    # Sort labels alphabetically / natural string order
    ordered_labels = sorted(buckets.keys(), key=lambda s: str(s))
    ordered = []
    groups = []
    for lab in ordered_labels:
        start = len(ordered)
        ordered.extend(buckets[lab])
        groups.append((lab, start, len(ordered)))
    return ordered, groups


def _pseudobulk_mean_logcpm(G_an, A_obs_leiden, var_names, libsize_add=1.0):
    import numpy as np, pandas as pd
    X = None
    try:
        if "counts" in getattr(G_an, "layers", {}):
            X = G_an.layers["counts"]
    except Exception:
        X = None
    if X is None:
        X = G_an.X
    X = _to_dense_matrix(X).astype(float)
    labs = list(map(str, A_obs_leiden))
    df = pd.DataFrame(X, columns=list(map(str, var_names)))
    df["__lab__"] = labs
    sums = df.groupby("__lab__", sort=True).sum()
    if "__lab__" in sums.columns:
        sums = sums.drop(columns=["__lab__"], errors="ignore")
    lib = sums.sum(axis=1).to_numpy().reshape(-1, 1) + libsize_add
    cpm = (sums.to_numpy() / lib) * 1e6
    logcpm = np.log1p(cpm)
    out = pd.DataFrame(logcpm, index=sums.index, columns=sums.columns)
    return out

def plot_marker_heatmap_and_export(gene_ad, ad_with_obs, marker_pairs: List[Tuple[str,str]], figdir: str, basename="heatmap_reference_markers",
                                   label_key="leiden", z_std_thresh=0.05, logcpm_std_thresh=0.02,
                                   clamp_robustz=None, clamp_logcpm=None):
    import numpy as np, pandas as pd, matplotlib.pyplot as plt, anndata as an, os

    if not isinstance(gene_ad, an.AnnData):
        gene_ad = gene_ad.to_adata()
    A, G = _align_obs(ad_with_obs, gene_ad)

    # --- NEW: keep full reference order/groups for display ---
    ordered_all_genes, groups_all = _build_marker_groups(marker_pairs)

    # present-only for scoring
    present = set(map(str, G.var_names))
    filtered_pairs = [(g, lab) for (g, lab) in marker_pairs if g in present]
    if len(filtered_pairs) == 0:
        print(f"[{ts()}] None of the markers were found in gene activity; skip heatmap/scoring.")
        return
    # ordered_genes used previously only for plotting; we will use full ref groups
    ordered_genes, groups = _build_marker_groups(filtered_pairs)

    X_log1p = np.log1p(_to_dense_matrix(G.X).astype(float))
    Z = _robust_z_per_gene_dense(X_log1p)
    medZ = _cluster_median_by_leiden(Z, A.obs[label_key], G.var_names)
    meanLogCPM = _pseudobulk_mean_logcpm(G, A.obs[label_key], G.var_names)

    # --- NEW: zero-pad missing genes to show all reference markers ---
    for g in ordered_all_genes:
        if g not in medZ.columns:       medZ[g] = 0.0
        if g not in meanLogCPM.columns: meanLogCPM[g] = 0.0
    # reindex to full reference order for CSV/plotting
    medZ = medZ.reindex(columns=ordered_all_genes)
    meanLogCPM = meanLogCPM.reindex(columns=ordered_all_genes)

    os.makedirs(figdir, exist_ok=True)
    csv_z = os.path.join(figdir, f"{basename}_robustz_cluster_median.csv")
    csv_cpm = os.path.join(figdir, f"{basename}_meanLogCPM.csv")
    medZ.to_csv(csv_z, index=True); print(f"[{ts()}] Saved CSV: {csv_z}")
    meanLogCPM.to_csv(csv_cpm, index=True); print(f"[{ts()}] Saved CSV: {csv_cpm}")

    def _filter_low_var(df, thr):
        if thr is None or thr <= 0: 
            return df
        sd = df.std(axis=0)
        keep = sd[sd >= thr].index.tolist()
        if len(keep) == 0:
            return df
        return df[keep]

    medZ_f = _filter_low_var(medZ, z_std_thresh)
    meanLogCPM_f = _filter_low_var(meanLogCPM, logcpm_std_thresh)

    # --- NEW: union of columns across both metrics and fillna(0) ---
    cols_union = [g for g in ordered_all_genes if (g in set(medZ_f.columns)) or (g in set(meanLogCPM_f.columns))]
    if len(cols_union) == 0:
        cols_union = ordered_all_genes
    medZ_f = medZ_f.reindex(columns=cols_union).fillna(0.0)
    meanLogCPM_f = meanLogCPM_f.reindex(columns=cols_union).fillna(0.0)

    from collections import OrderedDict
    var_index = {v:i for i,v in enumerate(map(str, G.var_names))}
    label2genes = OrderedDict()
    for g, lab in filtered_pairs:
        if g in medZ.columns:
            label2genes.setdefault(lab, []).append(g)

    set_scores_cell = {}
    for lab, genes in label2genes.items():
        idxs = [var_index[g] for g in genes if g in var_index]
        if not idxs: 
            continue
        vals = Z[:, idxs].mean(axis=1)
        set_scores_cell[lab] = vals

    labs_leiden = list(map(str, A.obs[label_key]))
    unique_leiden = sorted(set(labs_leiden), key=lambda x: (int(x) if x.isdigit() else 9_999_999, x))
    rec = []
    for c in unique_leiden:
        mask_pos = np.array([lab==c for lab in labs_leiden], dtype=bool)
        mask_neg = ~mask_pos
        for lab, vals in set_scores_cell.items():
            pos_vals = vals[mask_pos]; neg_vals = vals[mask_neg]
            auc = _auroc_pos_vs_neg(pos_vals, neg_vals)
            gene_list = label2genes[lab]
            mean_z = float(medZ.loc[c, [g for g in gene_list if g in medZ.columns]].mean())
            rec.append({"leiden": c, "ref_label": lab, "mean_z": mean_z, "auroc": float(auc)})
    score_df = pd.DataFrame(rec)
    score_df["combined"] = (score_df["auroc"] + (score_df["mean_z"] - score_df["mean_z"].min()) / (score_df["mean_z"].max() - score_df["mean_z"].min() + 1e-9)) / 2.0
    csv2 = os.path.join(figdir, "marker_set_scores.csv"); score_df.to_csv(csv2, index=False); print(f"[{ts()}] Saved CSV: {csv2}")

    assigns = []
    for c, sub in score_df.groupby("leiden"):
        sub = sub.sort_values("combined", ascending=False).reset_index(drop=True)
        best = sub.iloc[0]
        margin = float("nan")
        if len(sub) > 1:
            margin = float(best["combined"] - sub.iloc[1]["combined"])
        assigns.append({
            "leiden": c, "best_label": best["ref_label"], "best_score": float(best["combined"]), "margin": float(margin),
            "best_mean_z": float(best["mean_z"]), "best_auroc": float(best["auroc"])
        })
    
    # Build DataFrame and apply unknown rule once using configured thresholds
    assign_df = pd.DataFrame(assigns)
    # ensure consistent ordering by numeric leiden if possible
    try:
        assign_df['leiden_order'] = pd.to_numeric(assign_df['leiden'], errors='coerce')
        assign_df = assign_df.sort_values(['leiden_order','leiden']).drop(columns=['leiden_order'])
    except Exception:
        assign_df = assign_df.sort_values('leiden', ascending=True)
    # Apply unknown rule
    try:
        _score  = pd.to_numeric(assign_df['best_score'], errors='coerce')
        _margin = pd.to_numeric(assign_df['margin'], errors='coerce')
        if USE_MARGIN:
            _mask_unknown = (_score < MIN_SCORE) | (_margin < MIN_MARGIN)
        else:
            _mask_unknown = (_score < MIN_SCORE)
        assign_df.loc[_mask_unknown, 'best_label'] = 'unknown'
        print(f"[{ts()}] Apply unknown rule: MIN_SCORE={MIN_SCORE}, USE_MARGIN={USE_MARGIN}, MIN_MARGIN={MIN_MARGIN}; unknown={_mask_unknown.sum()} clusters")
    except Exception as _err_rule:
        print(f"[{ts()}] Warning: failed to apply unknown rule: {{_err_rule}}")
    # Save cluster assignments
    csv3 = os.path.join(figdir, "cluster_assignments.csv")
    assign_df.to_csv(csv3, index=False)
    print(f"[{ts()}] Saved CSV: {csv3}")

    def _draw(df, title, out_base, clamp=None):
        import matplotlib.pyplot as plt

        cols = list(df.columns)
        idx_map = {g:i for i,g in enumerate(cols)}
        # --- NEW: use full reference groups for separators/labels ---
        group_spans = []
        for lab, s0, e0 in groups_all:
            sub = [g for g in cols if g in ordered_all_genes[s0:e0]]
            if not sub: 
                continue
            s = idx_map[sub[0]]; e = idx_map[sub[-1]]+1
            group_spans.append((lab, s, e))
        fig_w = max(8.0, min(28.0, 0.45 * len(cols) + 4.0))
        fig_h = max(4.0, min(20.0, 0.5 * max(2, len(df.index)) + 2.0))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
        vmin = vmax = None
        if clamp is not None:
            vmin, vmax = clamp
        im = ax.imshow(df.values, aspect="auto", origin="upper", vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
        ax.set_yticks(range(len(df.index))); ax.set_yticklabels(list(map(str, df.index)), fontsize=9)
        ax.set_xlabel("Reference marker genes (grouped by clusterName)"); ax.set_ylabel("Leiden cluster (0..N)")
        for _, s, e in group_spans:
            ax.axvline(s-0.5, linestyle=":", linewidth=0.6); ax.axvline(e-0.5, linestyle=":", linewidth=0.6)
        y_top = -0.7
        for lab, s, e in group_spans:
            mid = (s + e - 1)/2.0
            ax.text(mid, y_top, str(lab), ha="center", va="bottom", rotation=90, fontsize=8)
        # === Top color bar per clusterName ===
        try:
            import matplotlib.pyplot as _plt
            import matplotlib as _mpl
            # create a slim inset axis above the heatmap
            ax_top = ax.inset_axes([0, 1.02, 1, 0.05], transform=ax.transAxes)
            ax_top.set_xlim(0, len(cols)); ax_top.set_ylim(0, 1)
            ax_top.axis("off")
            # deterministic palette across up to 20 labels
            cmap = _mpl.cm.get_cmap("tab20")
            lab_list = [lab for lab, _, _ in group_spans]
            lab_to_color = {lab: cmap(i % 20) for i, lab in enumerate(lab_list)}
            for lab, s, e in group_spans:
                ax_top.add_patch(_mpl.patches.Rectangle((s, 0), e - s, 1, transform=ax_top.transData, linewidth=0, facecolor=lab_to_color[lab]))
        except Exception as _e:
            # fail-safe: if anything goes wrong, just skip the color bar
            pass
        cbar = fig.colorbar(im, ax=ax); cbar.set_label(title)
        fig.tight_layout()
        for ext in ("png","svg"):
            out = os.path.join(figdir, f"{out_base}.{ext}")
            fig.savefig(out, dpi=300, facecolor="white"); print(f"[{ts()}] Saved figure: {out}")
        plt.close(fig)

    _draw(medZ_f, "Robust z (cluster median)", f"{basename}_robustz", clamp=clamp_robustz)
    # --- NEW: clearer fragCPM label ---
    _draw(meanLogCPM_f, "Pseudobulk accessibility (log1p, fragCPM)", f"{basename}_meanlogcpm", clamp=clamp_logcpm)
    print(f"[{ts()}] Heatmaps + CSV export completed.")
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



    # Reference-marker heatmap + scoring + CSV export
    try:
        marker_pairs = _load_reference_markers_with_labels(REF_CSV_PATH, REF_LABEL_KEY)
    except Exception as e:
        print(f"[{ts()}] Marker CSV load failed: {e}")
        marker_pairs = []
    gene_ad_for_heatmap = None
    try:
        if gene_matrix is not None:
            gene_ad_for_heatmap = gene_matrix
    except Exception:
        pass
    if gene_ad_for_heatmap is not None and len(marker_pairs) > 0:
        plot_marker_heatmap_and_export(gene_ad_for_heatmap, data, marker_pairs, figdir=FIG_DIR, basename="heatmap_reference_markers", label_key="leiden")
    else:
        print(f"[{ts()}] Skip heatmap (gene_activity or markers unavailable).")
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
    # === UMAP by majority (pure-CSV; no pandas; robust obs access) ===
    try:
        _csv_assign = os.path.join(FIG_DIR, "cluster_assignments.csv")
        if os.path.exists(_csv_assign):
            import csv as _csv
            # Sniff delimiter
            with open(_csv_assign, "r", encoding="utf-8", newline="") as _fh:
                _sample = _fh.read(4096)
                _fh.seek(0)
                try:
                    _dialect = _csv.Sniffer().sniff(_sample, delimiters=[",", "\t", ";", "|"])
                except Exception:
                    class _D: delimiter = ","
                    _dialect = _D()
                _reader = _csv.DictReader(_fh, delimiter=_dialect.delimiter)
                _rows = []
                for _row in _reader:
                    if not _row:
                        continue
                    _norm = {}
                    for _k, _v in _row.items():
                        _kk = ("" if _k is None else str(_k)).strip().lower().replace(" ", "_")
                        _vv = "" if _v is None else str(_v).strip()
                        _norm[_kk] = _vv
                    _rows.append(_norm)

            # Accept variants of headers
            _leiden_keys = ["leiden", "cluster", "cluster_id", "clusterid"]
            _best_keys = ["best_label", "bestlabel", "best", "best_label_name"]
            _have_keys = list(_rows[0].keys()) if _rows else []

            _lk = next((k for k in _leiden_keys if k in _have_keys), None)
            _bk = next((k for k in _best_keys if k in _have_keys), None)
            if _lk is None or _bk is None:
                raise KeyError(f"cluster_assignments.csv missing required headers; have={_have_keys}")

            _map = {}
            for _r in _rows:
                _lv = str(_r.get(_lk, "")).strip()
                _bv = str(_r.get(_bk, "")).strip()
                if _lv == "":
                    continue
                _map[_lv] = _bv if _bv != "" else "unknown"

            # Robust obs access (no .columns / no .astype)
            def _obs_get_str_labels(_adata, _cands):
                _obs = getattr(_adata, "obs", None)
                if _obs is None:
                    raise KeyError("AnnData has no .obs")
                # dict-like first
                for _k in _cands:
                    try:
                        _col = _obs[_k]
                        try:
                            import numpy as _np
                            return list(map(str, _np.asarray(_col).tolist()))
                        except Exception:
                            return [str(x) for x in _col]
                    except Exception:
                        continue
                # attribute-style
                for _k in _cands:
                    try:
                        _col = getattr(_obs, _k)
                        try:
                            import numpy as _np
                            return list(map(str, _np.asarray(_col).tolist()))
                        except Exception:
                            return [str(x) for x in _col]
                    except Exception:
                        continue
                raise KeyError("No leiden/cluster-like column found in AnnData.obs")

            _labs = _obs_get_str_labels(data, ["leiden", "cluster", "clusters", "clusterName", "cluster_name"])
            data.obs["majority"] = np.array([_map.get(x, "unknown") for x in _labs], dtype=object)

            # Ensure UMAP exists
            if "X_umap" not in data.obsm:
                try:
                    snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=USE_DIMS)
                except Exception:
                    pass

            plot_umap_matplotlib(data, color="majority", basename="umap_by_majority", figdir=FIG_DIR)
            # ---- HOOK: post-majority update/export (auto-inserted) ----
            try:
                from fragmentBED import update_and_export

                h5ad_dir    = os.path.join(BASE_DIR, "Filter_results")          # per-sample *.h5ad folder
                out_update  = os.path.join(AGGR_DIR, "update_h5ad")   # *_updated.h5ad output
                out_bed     = os.path.join(AGGR_DIR, "bed")           # BED/BED4 output
                summary_txt = os.path.join(AGGR_DIR, "summary.txt")

                update_and_export(
        h5ad_dir=os.path.join(BASE_DIR, 'Filter_results'),
        out_update_dir=os.path.join(AGGR_DIR, 'update_h5ad'),
        out_bed_dir=os.path.join(AGGR_DIR, 'bed'),
        njob=NJOB,
        summary_path=os.path.join(AGGR_DIR, 'summary.txt'),
        majority_key='majority',
        leiden_key='leiden',
        sample_key='sample',
        cluster_csv=os.path.join(AGGR_DIR, 'Plots', 'cluster_assignments.csv')
    )
                print(f"[{ts()}] Post-majority update/export completed.")
            except Exception as _e:
                print(f"[{ts()}] WARNING: post-majority step skipped: {_e}")
            # ---- END HOOK ----

        else:
            print(f"[{ts()}] Warning: cluster_assignments.csv not found; skipping umap_by_majority.")
    except Exception as e:
        print(f"[{ts()}] Warning: UMAP by majority failed: {e}")
if __name__ == "__main__":
    main()