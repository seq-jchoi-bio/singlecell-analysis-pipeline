#!/usr/bin/env python3

"""
Post-clustering utility:
- Load clustered .h5ad (AnnData)
- Auto-detect genome resources from species (-s): ~/singlecell-analysis-pipeline/refGenome/<species>/{genes.gtf, genome.chrom.sizes}
- Build gene activity matrix (snapATAC2 native -> robust fallback overlap)
- Copy UMAP from clustered object to gene matrix
- Save gene activity as Annot_results/annot_gene_activity.h5ad
- Plot UMAPs by sample/leiden/depth under Annot_results/Plots
- Export UMAP CSV (full obs + coords)
- Add gene_name column into .var alongside gene_id-index
CLI:
  python postCLS.py -i <.h5ad or project_dir> [-o <file_or_dir>] -s <marker_path.csv>
"""

import os
import re
import sys
import argparse
from datetime import datetime
from typing import Dict, Optional, List, Tuple
from collections import defaultdict

# ---- Thread caps ----
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"]      = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"
os.environ["RAYON_NUM_THREADS"]    = "1"
os.environ["MALLOC_ARENA_MAX"]     = "2"

# ---- Plot backend and core lib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import anndata as an
import snapatac2 as snap
from scipy import sparse

# ---- Silence warnings and logs ----
import warnings, logging
from matplotlib import MatplotlibDeprecationWarning

warnings.filterwarnings("ignore", message=r"`_import_from_c` is deprecated", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden", category=UserWarning, module=r"umap\.umap_")
warnings.filterwarnings("ignore", message=r".*`np\.object`.*", category=FutureWarning)
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

warnings.filterwarnings("ignore")

logging.getLogger().setLevel(logging.ERROR)
for name in (
    "harmonypy",
    "scanpy",
    "anndata",
    "matplotlib",
    "numba",
    "pyo3_polars",
    "umap",
):
    logging.getLogger(name).setLevel(logging.ERROR)

__version__ = "1.3"

# ---- I/O defaults ----
BASE_DIR        = "."
ANNOT_DIR_DEF   = os.path.join(BASE_DIR, "Annot_results")
FIG_DIR_NAME    = "plots"
CLUSTERED_DEF   = os.path.join(ANNOT_DIR_DEF, "clustered.h5ad")
GENE_H5AD_NAME  = "annot_gene_activity.h5ad"

# ---- Utils ----
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log(msg: str) -> None:
    print(f"[{ts()}] {msg}")

def _abs_path(base: str, p: str) -> str:
    if not p:
        return p
    return p if os.path.isabs(p) else os.path.abspath(os.path.join(base, p))

def _ensure_dirs(root: str) -> Tuple[str, str]:
    annot = os.path.join(root, "Annot_results")
    plots = os.path.join(annot, FIG_DIR_NAME)
    os.makedirs(annot, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    return annot, plots

def _normalize_outroot_and_paths(out_arg: Optional[str], annotdir_arg: Optional[str], cwd: str) -> Tuple[str, str, str]:
    if out_arg:
        out_abs = _abs_path(cwd, out_arg)
        if out_abs.lower().endswith((".h5ad", ".h5ads")):
            root = os.path.dirname(os.path.dirname(out_abs)) if os.path.basename(out_abs).startswith("clustered") else os.path.dirname(out_abs)
        else:
            root = out_abs
        annot, plots = _ensure_dirs(root)
        return root, annot, plots

    if annotdir_arg:
        annot = _abs_path(cwd, annotdir_arg)
        root  = os.path.dirname(annot) if os.path.basename(annot) == "Annot_results" else os.path.abspath(os.path.join(annot, ".."))
        os.makedirs(annot, exist_ok=True)
        plots = os.path.join(annot, FIG_DIR_NAME)
        os.makedirs(plots, exist_ok=True)
        return root, annot, plots

    root = os.path.abspath(cwd)
    annot, plots = _ensure_dirs(root)
    return root, annot, plots

def _detect_genome_paths_from_species(species: str) -> Tuple[Optional[str], Optional[str]]:
    home = os.path.expanduser("~")
    base = os.path.join(home, "singlecell-analysis-pipeline", "refGenome", species)
    gtf  = os.path.join(base, "genes.gtf")
    chrom = os.path.join(base, "genome.chrom.sizes")
    if not os.path.exists(gtf):
        gtf = None
    if not os.path.exists(chrom):
        chrom = None
    return gtf, chrom

def _read_chrom_sizes(path: Optional[str]) -> Optional[Dict[str, int]]:
    if not path or not os.path.exists(path):
        return None
    d: Dict[str, int] = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            t = line.split()
            if len(t) >= 2:
                chrom = re.sub(r'^(chr|Chr)', '', t[0]).lower()
                chrom = re.sub(r'^NC_[0-9.]+', '', chrom)
                try:
                    d[chrom] = int(t[1])
                except Exception:
                    pass
    return d or None

def _normalize_chrom(s: str) -> str:
    s0 = str(s).strip()
    s0 = re.sub(r'^(chr|Chr)', '', s0)
    s0 = re.sub(r'^NC_[0-9.]+', '', s0)
    return s0.lower()

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

def _parse_peak_varname(name: str):
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

def _load_gene_intervals(gtf_path: str, flank_up=3000, flank_down=3000) -> Dict[str, List[Tuple[int,int,str]]]:
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

def load_gene_id_name_mapping(gtf_gene_only: str) -> Dict[str, str]:
    pat_gid = re.compile(r'gene_id\s+"([^"]+)"')
    pat_gnm = re.compile(r'gene_name\s+"([^"]+)"')
    mapping: Dict[str, str] = {}
    with open(gtf_gene_only, "r") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2].lower() != "gene":
                continue
            attrs = cols[8]
            mid = pat_gid.search(attrs)
            mname = pat_gnm.search(attrs)
            if not mid:
                continue
            gid = mid.group(1)
            gname = mname.group(1) if mname else gid
            mapping[gid] = gname
    return mapping

def attach_gene_names(adata: an.AnnData, gid_to_name: Dict[str, str]) -> None:
    """IMPORTANT!! ensure adata.var has a 'gene_name' column aligned to var index (gene_id)."""
    if adata is None or adata.var is None:
        return
    idx = adata.var.index.astype(str)
    mapped = pd.Index([gid_to_name.get(g, g) for g in idx], name="gene_name")
    adata.var["gene_name"] = pd.Series(mapped.values, index=adata.var.index)

def _as_csr_float32(X):
    """Return X as CSR float32 -> required for making matrix!! ."""
    try:
        if sparse.issparse(X):
            return X.tocsr().astype(np.float32, copy=False)
        arr = np.asarray(X, dtype=np.float32)
        return sparse.csr_matrix(arr)
    except Exception:
        return sparse.csr_matrix(np.array(X, dtype=np.float32, copy=False))

# ----- Gene activity (matrix) builders -----
def make_gene_matrix_geneid_only(dataset: snap.AnnDataSet, gtf_path: str, chrom_sizes: Optional[Dict[str, int]] = None) -> an.AnnData:
    norm_gtf = os.path.join(os.path.dirname(gtf_path), "genes.gene_only.auto.gtf")
    norm_gtf = build_gene_only_gtf_from_any(gtf_path, norm_gtf)
    gid_to_name = load_gene_id_name_mapping(norm_gtf)

    try:
        ga = snap.pp.make_gene_matrix(dataset, annotation=norm_gtf, id_type="gene_id")
        attach_gene_names(ga, gid_to_name)
        return ga
    except Exception:
        pass
    try:
        ga = snap.pp.make_gene_matrix(dataset, annotation=norm_gtf)
        attach_gene_names(ga, gid_to_name)
        return ga
    except Exception:
        pass
    try:
        if chrom_sizes is not None:
            ga = snap.pp.make_gene_matrix(dataset, gtf=norm_gtf, chrom_sizes=chrom_sizes)
        else:
            ga = snap.pp.make_gene_matrix(dataset, gtf=norm_gtf)
        attach_gene_names(ga, gid_to_name)
        return ga
    except Exception:
        pass
    try:
        if hasattr(snap.genome, "read_gtf"):
            ann = snap.genome.read_gtf(norm_gtf)
            try:
                ga = snap.pp.make_gene_matrix(dataset, annotation=ann)
            except Exception:
                ga = snap.pp.make_gene_matrix(dataset, ann)
            attach_gene_names(ga, gid_to_name)
            return ga
    except Exception:
        pass

    ga = _fallback_gene_activity_via_overlap(dataset, norm_gtf, flank_up=3000, flank_down=3000)
    attach_gene_names(ga, gid_to_name)
    return ga

def _fallback_gene_activity_via_overlap(dataset, gene_gtf: str, flank_up=3000, flank_down=3000) -> an.AnnData:
    ad = dataset.to_adata() if hasattr(dataset, "to_adata") else dataset
    try:
        ad.obs_names_make_unique()
    except Exception:
        pass

    peak_coords_by_chr: Dict[str, List[Tuple[int,int,int]]] = defaultdict(list)
    for j, name in enumerate(map(str, ad.var_names)):
        parsed = _parse_peak_varname(name)
        if parsed is None:
            continue
        chrom, s, e = parsed
        peak_coords_by_chr[chrom].append((s, e, j))
    if not peak_coords_by_chr:
        raise RuntimeError("Fallback gene activity: var_names are not interval-like (no chrom:start-end).")
    for chrom in peak_coords_by_chr:
        peak_coords_by_chr[chrom].sort(key=lambda x: x[0])

    genes_by_chr = _load_gene_intervals(gene_gtf, flank_up=flank_up, flank_down=flank_down)
    if not genes_by_chr:
        raise RuntimeError("Fallback gene activity: no gene intervals loaded from GTF.")

    g2p = _map_peaks_to_genes(peak_coords_by_chr, genes_by_chr)
    if not g2p:
        raise RuntimeError("Fallback gene activity: no overlapping peaks and genes.")

    gene_ids = sorted(g2p.keys())
    gene_pos = {g:i for i,g in enumerate(gene_ids)}
    rows, cols = [], []
    for g, peaks in g2p.items():
        gi = gene_pos[g]
        for pj in peaks:
            rows.append(gi); cols.append(pj)
    vals = np.ones(len(rows), dtype=np.float32)
    P = sparse.csr_matrix((vals, (rows, cols)), shape=(len(gene_ids), ad.n_vars), dtype=np.float32)

    X_csr = _as_csr_float32(ad.X)
    G = (X_csr @ P.T.tocsc()).astype(np.float32)

    ga = an.AnnData(
        X=G,
        obs=ad.obs.copy(),
        var=pd.DataFrame(index=pd.Index(gene_ids, name="gene_id")),
        obsm=ad.obsm.copy() if hasattr(ad, "obsm") else None,
        uns=ad.uns.copy() if isinstance(ad.uns, dict) else {},
    )
    return ga

# ---- Plot ----
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
        ax.scatter(X[m, 0], X[m, 1], s=5, alpha=0.6, color=colors[i], linewidths=0, label=str(cat))

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
    
        if not any(key in basename for key in ("umap_by_sample", "umap_by_leiden")):
            log(f"Saved figure: {out}")
    
    plt.close(fig)

def plot_umap_matplotlib(adata, color: str, basename: str, figdir: str) -> None:
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    col = _obs_as_str_array(adata, color)
    _scatter_by_category(X, col, "UMAP", basename, figdir)

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
            label = "log10 Fragments"
            return arr, label
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
    
        if "umap_by_depth" not in basename:
            log(f"Saved figure: {out}")
    
    plt.close(fig)

def _obs_to_full_dataframe(adata):
    n = getattr(adata, "n_obs", None)
    obs = getattr(adata, "obs", None)
    if obs is None:
        raise RuntimeError("AnnData has no .obs")
    try:
        if hasattr(obs, "to_pandas"):
            df = obs.to_pandas()
        elif hasattr(obs, "to_df"):
            df = obs.to_df()
        else:
            df = pd.DataFrame(obs.copy())
        if n is not None and len(df) != n:
            df = df.reset_index(drop=True)
        try:
            cells = list(map(str, adata.obs_names))
        except Exception:
            try:
                cells = list(map(str, adata.obs.index))
            except Exception:
                cells = [str(i) for i in range(n or len(df))]
        if len(df) != len(cells):
            df = df.iloc[:len(cells)].reset_index(drop=True)
        df.insert(0, "cell", cells[:len(df)])
        return df
    except Exception:
        keys = []
        try:
            cols = getattr(obs, "columns", None)
            if cols is not None:
                keys = list(cols)
        except Exception:
            pass
        d = {}
        for k in keys:
            try:
                col = obs[k]
                try:
                    d[k] = col.tolist()
                except Exception:
                    d[k] = list(col)
            except Exception:
                continue
        m = None
        for v in d.values():
            if isinstance(v, (list, tuple)):
                m = len(v); break
        if m is None:
            m = int(n or 0)
        try:
            cells = list(map(str, adata.obs_names))
        except Exception:
            try:
                cells = list(map(str, adata.obs.index))
            except Exception:
                cells = [str(i) for i in range(m)]
        df = pd.DataFrame({"cell": cells})
        for k, v in d.items():
            df[k] = v
        return df

def export_umap_full_any(adata, out_path: str):
    try:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        df = _obs_to_full_dataframe(adata)
        if "X_umap" in adata.obsm:
            um = adata.obsm["X_umap"]
            if hasattr(um, "toarray"):
                um = um.toarray()
            um = np.asarray(um)
            if um.shape[1] >= 2:
                df["umap_x"] = um[:, 0]
                df["umap_y"] = um[:, 1]
        df.to_csv(out_path, index=False)
        log(f"Saved CSV: {out_path}")
        return True
    except Exception as e:
        log(f"Warning: failed to write {out_path}: {e}")
        return False

# ---- Main ----
def main():
    ap = argparse.ArgumentParser(description="Post-clustering: make gene-matrix and plots")
    ap.add_argument("-i", "--input", dest="inp", default=CLUSTERED_DEF, help="Input clustered gene_annot_activity.h5ad path")
    ap.add_argument("-s", "--species", dest="species", required=True, help="Species tag (e.g., hg38, mm10, rice)")
    ap.add_argument("-o", "--output", dest="out_root", default=None, help="Output root directory (creates Annot_results/ and plots/ under it)")
    ap.add_argument("--annotdir", default=None, help="Explicit Annot_results directory (legacy option)")
    args = ap.parse_args()

    print("\n====== Step 2. Post-clustering: make gene-matrix (bin/postCLS.py) ======\n")
    cwd = os.getcwd()

    root_dir, annot_dir, fig_dir = _normalize_outroot_and_paths(args.out_root, args.annotdir, cwd)

    # Resolve clustered input
    inp_path = _abs_path(cwd, args.inp)
    if os.path.isdir(inp_path):
        inp_path = os.path.join(inp_path, "Annot_results", "clustered.h5ad")
    if not os.path.exists(inp_path):
        # fallback to normalized annot_dir
        cand = os.path.join(annot_dir, "clustered.h5ad")
        if os.path.exists(cand):
            inp_path = cand
    if not os.path.exists(inp_path):
        raise FileNotFoundError(f"Clustered h5ad not found: {inp_path}")

    # Genome paths
    species = args.species
    gtf_path, chrom_sizes_path = _detect_genome_paths_from_species(species)
    if not gtf_path:
        raise FileNotFoundError(f"GTF not found for species '{species}'. Expected: ~/singlecell-analysis-pipeline/refGenome/{species}/genes.gtf")
    if not chrom_sizes_path:
        log(f"Warning: chrom.sizes not found for species '{species}'. Proceeding without chrom_sizes.")

    log(f"Input clustered: {inp_path}")
    log(f"GTF: {gtf_path}")
    if chrom_sizes_path:
        log(f"chrom.sizes: {chrom_sizes_path}")

    # Load clustered AnnData and depict plots
    log("Loading object...")
    adata = an.read_h5ad(inp_path)

    if "X_umap" not in adata.obsm:
        log("UMAP missing in clustered; computing...")
        try:
            snap.tl.umap(adata, use_rep="X_spectral_harmony")
        except Exception as e:
            log(f"Warning: UMAP computation failed: {e}")

    try:
        plot_umap_matplotlib(adata, color="sample", basename="umap_by_sample", figdir=fig_dir)
    except Exception as e:
        log(f"Warning: UMAP by sample failed: {e}")
    try:
        plot_umap_matplotlib(adata, color="leiden", basename="umap_by_leiden", figdir=fig_dir)
    except Exception as e:
        log(f"Warning: UMAP by leiden failed: {e}")
    try:
        plot_umap_by_depth(adata, basename="umap_by_depth", figdir=fig_dir)
    except Exception as e:
        log(f"Warning: UMAP by depth failed: {e}")

    export_umap_full_any(adata, os.path.join(fig_dir, "umap.csv"))

    # Make gene matrix
    log("Making gene matrix...")
    chrom_sizes_dict = _read_chrom_sizes(chrom_sizes_path)
    try:
        gene_ad = make_gene_matrix_geneid_only(adata, gtf_path=gtf_path, chrom_sizes=chrom_sizes_dict)
        log("Gene activity via snapatac2 / normalized GTF OK")
    except Exception as e1:
        log(f"make_gene_matrix failed; trying strict fallback: {e1}")
        try:
            auto_gtf = os.path.join(os.path.dirname(gtf_path), "genes.gene_only.auto.gtf")
            if not os.path.exists(auto_gtf):
                auto_gtf = build_gene_only_gtf_from_any(gtf_path, auto_gtf)
            gene_ad = _fallback_gene_activity_via_overlap(adata, auto_gtf, flank_up=2000, flank_down=0)
            gid_to_name = load_gene_id_name_mapping(auto_gtf)
            attach_gene_names(gene_ad, gid_to_name)
            log("Gene activity via strict fallback overlap OK")
        except Exception as e2:
            raise RuntimeError(f"Gene activity failed (native+fallback): {e2}")

    # Copy UMAP
    try:
        if "X_umap" in adata.obsm:
            gene_ad.obsm["X_umap"] = adata.obsm["X_umap"]
    except Exception:
        pass

    # Save gene matrix h5ad
    gene_out = os.path.join(annot_dir, GENE_H5AD_NAME)
    try:
        gene_ad.write(gene_out)
        log(f"Saved gene activity AnnData: {gene_out}")
    except Exception as e:
        log(f"Warning: failed to save gene activity h5ad: {e}")

    # Additional UMAPs from gene_ad
    try:
        plot_umap_matplotlib(gene_ad, color="leiden", basename="umap_gene_by_leiden", figdir=fig_dir)
    except Exception as e:
        log(f"Warning: gene_umap by leiden failed: {e}")
    try:
        plot_umap_matplotlib(gene_ad, color="sample", basename="umap_gene_by_sample", figdir=fig_dir)
    except Exception as e:
        log(f"Warning: gene_umap by sample failed: {e}")

    log("Done.")

if __name__ == "__main__":
    main()
