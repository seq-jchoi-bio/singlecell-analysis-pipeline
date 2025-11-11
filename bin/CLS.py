#!/usr/bin/env python3

"""
- Input accepts a .h5ads file OR a project directory (e.g., '.')
- Pipeline: select_features → spectral → harmony(batch='sample') → knn → leiden → umap
- Output rules:
  * If -o is a FILE path (.h5ads/.h5ad): save clustered there (fallback .h5ad if needed)
  * If -o is a DIR path: save to <DIR>/Annot_results/clustered.h5ads (and figures to <DIR>/Annot_results/plots)
  * If -o omitted: save to ./Annot_results/clustered.h5ads (figures to ./Annot_results/plots)
- Figures: UMAP(by sample / by leiden / by depth), plus umap.csv and summary.txt

CLI:
  python CLS.py -i <.h5ads or project_dir> [-o <file_or_dir>] [--neighbors 10] [--dims None]
"""

import os
import sys
import argparse
from datetime import datetime

# ---- Thread caps (keep small and predictable) ----
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

import numpy as np
import pandas as pd
import snapatac2 as snap

# ---- Silence Python warnings & noisy library logs ----
import warnings, logging
from matplotlib import MatplotlibDeprecationWarning

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden", module=r"umap\.umap_")
warnings.filterwarnings("ignore", message=r"converged after", module=r"harmonypy")

logging.getLogger().setLevel(logging.ERROR)
for name in (
    "harmonypy",
    "scanpy",
    "anndata",
    "matplotlib",
    "numba",
    "pyo3_polars",
    "umap",
    "h5py",
):
    logging.getLogger(name).setLevel(logging.ERROR)

__version__ = "1.3"

# ---- Utils ----
def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def _patch_numpy_for_leiden():
    try:
        for name, val in {"float": float, "int": int, "bool": bool, "object": object}.items():
            if not hasattr(np, name):
                setattr(np, name, val)
        compat = getattr(np, "compat", None)
        if compat is None:
            class _Compat: ...
            np.compat = _Compat()
            compat = np.compat
        if not hasattr(compat, "unicode"):
            setattr(compat, "unicode", str)
    except Exception:
        pass

def _obs_as_str_array(adata, key: str):
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

def _scatter_by_category(X: np.ndarray, col: np.ndarray, embed_name: str, basename: str, figdir: str):
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

def plot_umap_matplotlib(adata, color: str, basename: str, figdir: str):
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    col = _obs_as_str_array(adata, color)
    _scatter_by_category(X, col, "UMAP", basename, figdir)

def _first_key(adata, candidates):
    for k in candidates:
        try:
            if k in adata.obs.columns:
                return k
        except Exception:
            pass
        try:
            _ = adata.obs[k]
            return k
        except Exception:
            pass
    return None

def export_umap_csv(adata, out_csv: str, include_nfrag: bool = True):
    """This defenition will be removed (redundancy!!!!)"""
    try:
        if "X_umap" not in adata.obsm:
            return False
        um = adata.obsm["X_umap"]
        if hasattr(um, "toarray"):
            um = um.toarray()
        um = np.asarray(um)
        if um.shape[1] < 2:
            return False

        sk = _first_key(adata, ["sample","Sample","sample_name","sampleid","sample_id"])
        lk = _first_key(adata, ["leiden","Leiden","cluster","Cluster"])
        nk = _first_key(adata, ["n_fragment","n_fragments","passed_filters","nFragments","depth"]) if include_nfrag else None

        sample_vals = adata.obs[sk].astype(str).tolist() if sk else [""] * adata.n_obs
        leiden_vals = adata.obs[lk].astype(str).tolist() if lk else [""] * adata.n_obs
        nfrag_vals  = adata.obs[nk].astype(float).tolist() if nk else ["" for _ in range(adata.n_obs)]

        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        df = pd.DataFrame({
            "sample": sample_vals,
            "leiden": leiden_vals,
            "n_fragment": nfrag_vals,
            "umap_x": um[:, 0],
            "umap_y": um[:, 1],
        })
        df.to_csv(out_csv, index=False)
        print(f"[{ts()}] Saved CSV: {out_csv}")
        return True
    except Exception as e:
        print(f"[{ts()}] Warning: failed to write {out_csv}: {e}")
        return False

def write_summary(out_txt: str, data, n_neighbors: int, use_dims):
    os.makedirs(os.path.dirname(out_txt), exist_ok=True)
    lines = [
        f"=== Clustering Summary (v{__version__}) ===",
        f"Generated at : {ts()}",
        "",
        f"Cells (n_obs)    : {data.n_obs}",
        f"Features (n_vars): {data.n_vars}",
        f"n_neighbors      : {n_neighbors}",
        f"use_dims         : {'all' if use_dims is None else use_dims}",
        "",
        f"Leiden clusters  : {len(set(map(str, data.obs['leiden']))) if 'leiden' in data.obs else 0}",
    ]
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[{ts()}] Summary written: {out_txt}")

def _resolve_input_path(user_inp: str) -> str:
    p = os.path.abspath(user_inp)
    if os.path.isfile(p) and p.lower().endswith(".h5ads"):
        return p
    if os.path.isdir(p):
        cands = [
            os.path.join(p, "Filter_results", "merged_doublets.h5ads"),
            os.path.join(p, "merged_doublets.h5ads"),
        ]
        for c in cands:
            if os.path.isfile(c):
                return c
        fr = os.path.join(p, "Filter_results")
        if os.path.isdir(fr):
            h5ads = [os.path.join(fr, f) for f in os.listdir(fr) if f.lower().endswith(".h5ads")]
            if h5ads:
                h5ads.sort()
                return h5ads[0]
        h5ads2 = [os.path.join(p, f) for f in os.listdir(p) if f.lower().endswith(".h5ads")]
        if h5ads2:
            h5ads2.sort()
            return h5ads2[0]
    return user_inp

def _normalize_output_paths(out_arg: str | None, fig_arg: str | None, cwd: str, project_hint: str | None) -> tuple[str, str]:
    base = os.path.abspath(project_hint or cwd)
    if out_arg is None:
        file_out = os.path.join(base, "Annot_results", "clustered.h5ads")
        figdir   = os.path.join(base, "Annot_results", "plots")
        return file_out, figdir

    out_abs = os.path.abspath(out_arg)
    if out_abs.lower().endswith(".h5ads") or out_abs.lower().endswith(".h5ad"):
        file_out = out_abs if out_abs.lower().endswith(".h5ads") else out_abs  # keep suffix; we convert to .h5ad later if needed
        # Figure dir: default next to file under <root>/Annot_results/plots
        root_dir = os.path.dirname(os.path.dirname(file_out)) if os.path.basename(file_out).lower().startswith("clustered") else os.path.dirname(file_out)
        if fig_arg:
            figdir = os.path.abspath(fig_arg) if os.path.isabs(fig_arg) else os.path.join(root_dir, fig_arg)
        else:
            figdir = os.path.join(root_dir, "plots")
        return file_out, figdir

    # Treat as directory root
    root = out_abs
    file_out = os.path.join(root, "Annot_results", "clustered.h5ads")
    if fig_arg:
        figdir = os.path.abspath(fig_arg) if os.path.isabs(fig_arg) else os.path.join(root, fig_arg)
    else:
        figdir = os.path.join(root, "Annot_results", "plots")
    return file_out, figdir

def _save_dataset_or_adata(data, out_path: str):
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    saved = False
    try:
        if hasattr(snap, "write_dataset"):
            snap.write_dataset(data, out_path)
            print(f"[{ts()}] Saved clustered dataset (SnapATAC2): {out_path}")
            saved = True
        else:
            raise AttributeError("snap.write_dataset not available")
    except Exception as e:
        print(f"[{ts()}] snap.write_dataset failed ({e}); fallback to AnnData .h5ad")

    if not saved:
        ad = data.to_adata()
        # If user used .h5ads, switch to .h5ad to avoid confusion
        if out_path.endswith(".h5ads"):
            out_path = out_path[:-1]
        ad.write(out_path)
        print(f"[{ts()}] Saved clustered AnnData: {out_path}")
    return True

def _obs_numeric_depth_series(adata):
    cands = [
        "log10_Fragments", "log10_fragments", "log10_depth",
        "n_fragment", "n_fragments", "passed_filters",
        "nFrags", "n_frag", "Fragments", "fragments", "depth"
    ]
    for nm in cands:
        try:
            if nm in adata.obs:
                vals = adata.obs[nm]
            else:
                continue
            try:
                arr = np.asarray(vals, dtype=float)
            except Exception:
                arr = np.asarray(vals.astype(float))
            if "log10" in nm.lower():
                return arr, "log10 Fragments"
            arr = np.log10(np.clip(arr, 1.0, None))
            return arr, "log10 Fragments"
        except Exception:
            continue
    return None, None

def plot_umap_by_depth(adata, basename: str, figdir: str):
    if "X_umap" not in adata.obsm:
        raise RuntimeError("X_umap embedding not found.")
    X = np.asarray(adata.obsm["X_umap"])
    vals, label = _obs_numeric_depth_series(adata)
    if vals is None:
        raise RuntimeError("No depth-like column found in .obs (e.g., n_fragment).")
    fig, ax = plt.subplots(figsize=(6, 6))
    sc = ax.scatter(X[:, 0], X[:, 1], c=vals, s=5, alpha=0.65, linewidth=0)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("UMAP by depth")
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

# ----- main -----
def main():
    ap = argparse.ArgumentParser(description="scATAC clustering (peak-matrix based)")
    ap.add_argument("-i", dest="inp", required=True, help="Input .h5ads OR a project dir (e.g., '.') to auto-resolve")
    ap.add_argument("-o", dest="outp", required=False, help="omitted -> ./Annot_results/clustered.h5ad")
    ap.add_argument("--figdir", default=None, help="Figure dir; default depends on -o normalization (…/Annot_results/plots)")
    ap.add_argument("--neighbors", type=int, default=10, help="k for KNN graph (default: 10)")
    ap.add_argument("--dims", default="None", help="Use dims (e.g., 0:50). 'None' = use all")
    ap.add_argument("--n_features", type=int, default=250_000, help="Feature cap for select_features")
    args = ap.parse_args()

    print("\n====== Step 1. Clustering (bin/CLS.py) ======\n")
    
    inp_resolved = _resolve_input_path(args.inp)
    if not os.path.exists(inp_resolved):
        raise FileNotFoundError(f"Input not found: {inp_resolved}")
    print(f"[{ts()}] Load dataset: {inp_resolved}")
    data = snap.read_dataset(inp_resolved)

    norm_out, norm_figdir = _normalize_output_paths(args.outp, args.figdir, cwd=os.getcwd(),
                                                    project_hint=os.path.dirname(inp_resolved) if os.path.isfile(inp_resolved) else os.path.abspath(args.inp))
    print(f"[{ts()}] select_features (cap={min(int(args.n_features), data.n_vars)})")
    snap.pp.select_features(data, n_features=min(int(args.n_features), data.n_vars))

    print(f"[{ts()}] spectral")
    snap.tl.spectral(data)

    print(f"[{ts()}] harmony (batch='sample')")
    use_dims = None
    if args.dims and str(args.dims).lower() != "none":
        if ":" in args.dims:
            s, e = args.dims.split(":", 1)
            try:
                s = int(s) if s != "" else None
                e = int(e) if e != "" else None
                use_dims = slice(s, e)
            except Exception:
                use_dims = None
        else:
            try:
                e = int(args.dims)
                use_dims = slice(0, e)
            except Exception:
                use_dims = None
    snap.pp.harmony(data, batch="sample", use_dims=use_dims, max_iter_harmony=20)

    n_neighbors_eff = max(2, min(args.neighbors, max(2, data.n_obs - 1)))
    print(f"[{ts()}] knn (n_neighbors={n_neighbors_eff})")
    snap.pp.knn(data, n_neighbors=n_neighbors_eff, use_dims=use_dims, use_rep="X_spectral_harmony")

    _patch_numpy_for_leiden()
    print(f"[{ts()}] leiden")
    snap.tl.leiden(data)

    print(f"[{ts()}] umap")
    snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=use_dims)

    _save_dataset_or_adata(data, norm_out)

    # Figures: CSV logic will be removed!!!
    try:
        plot_umap_matplotlib(data, color="sample", basename="umap_by_sample", figdir=norm_figdir)
    except Exception as e:
        print(f"[{ts()}] Warning: UMAP by sample failed: {e}")

    try:
        plot_umap_matplotlib(data, color="leiden", basename="umap_by_leiden", figdir=norm_figdir)
    except Exception as e:
        print(f"[{ts()}] Warning: UMAP by leiden failed: {e}")

    try:
        plot_umap_by_depth(data, basename="umap_by_depth", figdir=norm_figdir)
    except Exception as e:
        print(f"[{ts()}] Warning: UMAP by depth failed: {e}")

    try:
        export_umap_csv(data, os.path.join(norm_figdir, "umap.csv"))
    except Exception:
        pass

    try:
        write_summary(os.path.join(norm_figdir, "summary.txt"), data, n_neighbors_eff, use_dims)
    except Exception as e:
        print(f"[{ts()}] Warning: summary write failed: {e}")

    print(f"[{ts()}] Done.")

if __name__ == "__main__":
    main()
