#!/usr/bin/env python3

"""
Scanpy-based dotplot from annot_gene_activity.h5ad + marker CSV.
CSV schema (required):
  - 'type' (marker group name)
  - one of {'gene_id', 'gene_name'}
Behavior:
  - If matching by gene_id: x-axis shows gene_id.
  - If matching by gene_name: x-axis shows gene_name.
Outputs (default):
  - <OUT>/dot_mean_<groupby>.csv, <OUT>/dot_frac_<groupby>.csv
  - <OUT>/dot_mean_<groupby>_dendro.csv, <OUT>/dot_frac_<groupby>_dendro.csv
  - <OUT>/dotplot_dendro.png / .svg
  - <OUT>/dotplot_dendro_auto.png / .svg
  - <OUT>/label_scores.csv, <OUT>/label_map.csv (headers: leiden,type)

CLI:
  python dotplot.py -i <annot_gene_activity.h5ad> [-o <file_or_dir>] -r <marker.csv> [other params refer to the arg section]
"""

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Tuple

# ---- Silence all warnings and logs ----
import warnings, logging
from matplotlib import MatplotlibDeprecationWarning

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

warnings.filterwarnings("ignore", message=r".*dendrogram data not found.*")
warnings.filterwarnings("ignore", message=r".*agg function failed.*")
warnings.filterwarnings("ignore", message=r".*Running `sc\.tl\.dendrogram`.*")
warnings.filterwarnings("ignore", message=r".*Retrying without dendrogram.*")

logging.getLogger().setLevel(logging.ERROR)
for name in ("scanpy", "anndata", "umap", "matplotlib", "harmonypy", "numba"):
    logging.getLogger(name).setLevel(logging.ERROR)

def log(msg: str):
    # suppress all warnings lines
    if "Warning:" in msg or "WARNING:" in msg:
        return
    print(f"[scanpy-dotplot] {msg}")

# ---- Argparse ----
def parse_args():
    ap = argparse.ArgumentParser(description="Draw Scanpy dotplot from .h5ad and marker CSV")
    ap.add_argument("-i", "--input", required=True,
                    help="Path to .h5ad OR a project dir containing Annot_results/annot_gene_activity.h5ad")
    ap.add_argument("-r", "--reference", required=True,
                    help="Marker CSV with headers: type + (gene_id | gene_name)")
    ap.add_argument("-o", "--output", default=None,
                    help="Output directory (default: Annot_results/dotplots if detectable, else ./dotplots)")
    ap.add_argument("--groupby", default=None,
                    help="obs column to group by (default: auto: leiden → clusterName → sample)")
    ap.add_argument("--standard-scale", default="var", choices=["var", "group", "none", "None"],
                    help='Dot color scaling (default: "var")')
    ap.add_argument("--dot-min", type=float, default=0.1, help="Minimum dot size (default: 0.1)")
    ap.add_argument("--dot-max", type=float, default=1.0, help="Maximum dot size (default: 1.0)")
    ap.add_argument("--figsize", default="12,10", help="Figure size W,H (default: 12,10)")
    ap.add_argument("--cmap", default="Blues", help='Matplotlib colormap (default: "Blues")')
    ap.add_argument("--case-insensitive", action="store_true",
                    help="Case-insensitive match for gene_name")
    return ap.parse_args()

# ---- Paths ----
def resolve_h5ad_path(inp: str) -> Path:
    p = Path(inp).resolve()
    if p.is_file() and p.suffix.lower() == ".h5ad":
        return p
    if p.is_dir():
        cand = p / "Annot_results" / "annot_gene_activity.h5ad"
        if cand.exists():
            return cand
        cand2 = p / "annot_gene_activity.h5ad"
        if cand2.exists():
            return cand2
        h5s = sorted(p.glob("*.h5ad"), key=lambda x: x.stat().st_mtime, reverse=True)
        if h5s:
            return h5s[0]
    raise FileNotFoundError(f"Could not resolve .h5ad from: {inp}")

def default_outdir_for_input(inp: Path) -> Path:
    p = inp.resolve()
    if p.is_file():
        base = p.parent
        if base.name != "Annot_results":
            proj = base
            for _ in range(4):
                if (proj / "Annot_results").exists():
                    base = proj / "Annot_results"
                    break
                proj = proj.parent
            else:
                base = p.parent / "Annot_results"
                base.mkdir(parents=True, exist_ok=True)
    else:
        base = (p / "Annot_results") if p.name != "Annot_results" else p
        base.mkdir(parents=True, exist_ok=True)
    outdir = base / "dotplots"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

# ---- Marker loading ----
def load_markers_csv(csv_path: Path, case_insensitive: bool):
    import pandas as pd
    df = pd.read_csv(csv_path)
    cols_lower = {c.lower(): c for c in df.columns}
    if "type" not in cols_lower:
        raise ValueError("Marker CSV must contain 'type' header.")
    type_col = cols_lower["type"]

    gene_col = None
    if "gene_id" in cols_lower:
        gene_col = cols_lower["gene_id"]
        id_mode = "gene_id"
    elif "gene_name" in cols_lower:
        gene_col = cols_lower["gene_name"]
        id_mode = "gene_name"
    else:
        raise ValueError("Marker CSV must contain one of: 'gene_id' or 'gene_name'.")

    df = df[[type_col, gene_col]].dropna()
    df[gene_col] = df[gene_col].astype(str)
    df["_orig_label_"] = df[gene_col].astype(str)

    if id_mode == "gene_name" and case_insensitive:
        df[gene_col] = df[gene_col].str.upper()

    return df, type_col, gene_col, id_mode

# ---- Plots ----
def choose_groupby(adata, user_key: str | None) -> str:
    if user_key and user_key in adata.obs.columns:
        return user_key
    if "leiden" in adata.obs.columns:
        return "leiden"
    if "clusterName" in adata.obs.columns:
        return "clusterName"
    if "sample" in adata.obs.columns:
        return "sample"
    for c in adata.obs.columns:
        try:
            s = adata.obs[c].astype(str)
            if s.nunique() <= max(50, int(adata.n_obs * 0.5)):
                return c
        except Exception:
            pass
    raise KeyError("No suitable groupby column found. Please specify --groupby.")

def _ensure_categories_ordered(adata, groupby: str):
    import pandas as pd
    try:
        s = adata.obs[groupby].astype(str)
        s = pd.Categorical(s).remove_unused_categories()
        cats = list(pd.Series(s.astype(str)).unique())
        adata.obs[groupby] = pd.Categorical(s.astype(str), categories=cats, ordered=False)
    except Exception:
        pass

def build_selection_and_labels(adata, markers_df, type_col, gene_col, id_mode, case_insensitive: bool):
    selected_ids = {}
    id_to_label = {}
    missing = 0

    if id_mode == "gene_id":
        var_set = set(map(str, adata.var_names))
        for t, sub in markers_df.groupby(type_col, sort=False):
            genes_in = []
            for g, lab in zip(sub[gene_col].tolist(), sub["_orig_label_"].tolist()):
                if g in var_set:
                    genes_in.append(g)
                    if g not in id_to_label:
                        id_to_label[g] = g
                else:
                    missing += 1
            if genes_in:
                selected_ids[str(t)] = genes_in
    else:
        if "gene_name" not in adata.var.columns:
            raise KeyError("AnnData .var does not contain 'gene_name' column.")
        series = adata.var["gene_name"].astype(str)
        if case_insensitive:
            series_upper = series.str.upper()
            name_to_var = {}
            for var_id, nmU in zip(adata.var_names.astype(str), series_upper.values):
                if nmU and nmU not in name_to_var:
                    name_to_var[nmU] = var_id
        else:
            name_to_var = {}
            for var_id, nm in zip(adata.var_names.astype(str), series.values):
                if nm and nm not in name_to_var:
                    name_to_var[nm] = var_id

        for t, sub in markers_df.groupby(type_col, sort=False):
            genes_in = []
            for g, lab in zip(sub[gene_col].tolist(), sub["_orig_label_"].tolist()):
                key = g.upper() if case_insensitive else g
                var_id = name_to_var.get(key)
                if var_id is not None:
                    genes_in.append(var_id)
                    if var_id not in id_to_label:
                        id_to_label[var_id] = lab
                else:
                    missing += 1
            if genes_in:
                selected_ids[str(t)] = genes_in

    mode_msg = f"[scanpy-dotplot] Markers loaded using '{id_mode}'"
    if id_mode == "gene_name" and case_insensitive:
        mode_msg += " (case-insensitive)"
    print(mode_msg)
    if missing > 0:
        print(f"[scanpy-dotplot] Note: {missing} marker entries were not found and skipped.")

    if not selected_ids:
        raise ValueError("No markers matched AnnData var_names.")
    return selected_ids, id_to_label

def _make_relabeled_subset(adata, selected_ids: Dict[str, List[str]], id_to_label: Dict[str, str]):
    import anndata as ad
    var_ids_ordered = []
    for genes in selected_ids.values():
        var_ids_ordered.extend(genes)
    var_ids_ordered = list(dict.fromkeys(var_ids_ordered))

    A_sub = ad.AnnData(
        X=adata[:, var_ids_ordered].X.copy(),
        obs=adata.obs.copy(),
        var=adata[:, var_ids_ordered].var.copy(),
        obsm=adata.obsm.copy() if hasattr(adata, "obsm") else None,
        uns=adata.uns.copy() if isinstance(adata.uns, dict) else {},
    )

    orig_ids = [str(x) for x in var_ids_ordered]
    labels = [id_to_label[v] for v in orig_ids]

    seen = {}
    uniq_labels = []
    for x in labels:
        k = str(x)
        c = seen.get(k, 0)
        uniq_labels.append(k if c == 0 else f"{k}#{c+1}")
        seen[k] = c + 1

    A_sub.var_names = uniq_labels

    varid_to_label = {var_id: lbl for var_id, lbl in zip(orig_ids, uniq_labels)}
    selected_labels = {t: [varid_to_label[v] for v in genes if v in varid_to_label] for t, genes in selected_ids.items()}
    return A_sub, selected_labels

def _to_dense(X):
    try:
        if hasattr(X, "toarray"):
            return X.toarray()
        return X
    except Exception:
        return X

def _ensure_numeric_X(adata):
    import numpy as np
    from scipy import sparse
    X = adata.X
    if sparse.issparse(X):
        adata.X = X.tocsr().astype(np.float32)
    else:
        adata.X = np.asarray(X, dtype=np.float32)

def _group_categories(adata, groupby: str) -> List[str]:
    s = adata.obs[groupby]
    if not hasattr(s, "cat"):
        s = s.astype("category")
    return list(s.cat.categories)

class Out:
    def __init__(self, outdir: Path):
        self.outdir = outdir.resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)

    def save_csv(self, df, name: str):
        p = self.outdir / name
        df.to_csv(p, index=False)
        print(f"[scanpy-dotplot] Saved: {p}")
        return p

    def save_fig(self, stem: str, plt):
        for ext in ("png", "svg"):
            p = self.outdir / f"{stem}.{ext}"
            plt.savefig(p, bbox_inches="tight", dpi=300)
            print(f"[scanpy-dotplot] Saved: {p}")

def _export_dot_values(A_sub, groupby: str, var_labels: List[str], out: Out, suffix: str = ""):
    import numpy as np
    import pandas as pd
    cats = _group_categories(A_sub, groupby)
    means = []
    fracs = []
    for g in cats:
        mask = (A_sub.obs[groupby].astype(str) == g).values
        Xg = _to_dense(A_sub.X[mask, :])
        Xg = np.asarray(Xg, dtype=float)
        m = Xg.mean(axis=0)
        f = (Xg > 0).sum(axis=0) / max(1, Xg.shape[0])
        means.append(m)
        fracs.append(f)
    df_mean = pd.DataFrame(means, index=cats, columns=var_labels).reset_index(names=[groupby])
    df_frac = pd.DataFrame(fracs, index=cats, columns=var_labels).reset_index(names=[groupby])
    suf = f"_{suffix}" if suffix else ""
    out.save_csv(df_mean, f"dot_mean_{groupby}{suf}.csv")
    out.save_csv(df_frac, f"dot_frac_{groupby}{suf}.csv")

def _run_dotplot(A_sub, selected_labels, groupby: str, std_scale, dot_min, dot_max, cmap, figsize, out: Out, out_stem: str, dendro: bool = False):
    import matplotlib.pyplot as plt
    import scanpy as sc
    def _plot(_dendro: bool):
        sc.pl.dotplot(
            A_sub,
            var_names=selected_labels,
            groupby=groupby,
            standard_scale=std_scale,
            dot_min=dot_min,
            dot_max=dot_max,
            cmap=cmap,
            figsize=figsize,
            dendrogram=_dendro,
            use_raw=False,
            layer=None,
            show=False,
        )
    try:
        _plot(dendro)
    except Exception as e:
        sys.stderr.write(f"[scanpy-dotplot] Warning: plot failed with dendrogram={dendro} ({e}); retrying without dendrogram.\n")
        plt.close()
        _plot(False)
    out.save_fig(out_stem, plt)
    plt.close()

def _apply_category_order(adata, groupby: str, ordered: List[str]):
    import pandas as pd
    s = adata.obs[groupby].astype(str)
    adata.obs[groupby] = pd.Categorical(s, categories=ordered, ordered=True)

def _get_dendrogram_order(adata, groupby: str) -> List[str]:
    key = f"dendrogram_{groupby}"
    info = adata.uns.get(key, None)
    if not info:
        return _group_categories(adata, groupby)
    ordered = info.get("categories_ordered", None)
    if not ordered:
        return _group_categories(adata, groupby)
    return list(ordered)

def _compute_signature_scores(A_sub, groupby: str, selected_labels: Dict[str, List[str]]) -> Tuple[List[str], "pd.DataFrame"]:
    import numpy as np
    import pandas as pd
    groups = _group_categories(A_sub, groupby)
    label_to_idx = {lbl: i for i, lbl in enumerate(A_sub.var_names.astype(str).tolist())}
    type_to_cols = {}
    for t, lbls in selected_labels.items():
        idxs = [label_to_idx[lbl] for lbl in lbls if lbl in label_to_idx]
        if idxs:
            type_to_cols[t] = idxs
    scores = []
    for g in groups:
        mask = (A_sub.obs[groupby].astype(str) == g).values
        Xg = _to_dense(A_sub.X[mask, :])
        Xg = np.asarray(Xg, dtype=float)
        row = {}
        for t, cols in type_to_cols.items():
            sub = Xg[:, cols]
            row[t] = float(sub.mean()) if sub.size else 0.0
        scores.append(row)
    import pandas as pd
    scores_df = pd.DataFrame(scores, index=groups, columns=list(type_to_cols.keys()))
    return groups, scores_df

def _assign_auto_labels(groups: List[str], scores_df, prefer_prefix: str = "") -> Dict[str, str]:
    mapping = {}
    for g in groups:
        row = scores_df.loc[g]
        if row.isna().all() or (row.fillna(0).sum() == 0):
            mapping[g] = ("Unassigned", f"{g} · Unassigned")
            continue
        top = row.idxmax()
        label_text = f"{prefer_prefix}{g} · {top}" if prefer_prefix else f"{g} · {top}"
        mapping[g] = (top, label_text)
    return mapping  # group -> (type, display_text)

def _manual_dendro_order(A_sub, groupby):
    import numpy as np
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    cats = _group_categories(A_sub, groupby)
    means = []
    for g in cats:
        mask = (A_sub.obs[groupby].astype(str) == g).values
        Xg = _to_dense(A_sub.X[mask, :])
        Xg = np.asarray(Xg, dtype=float)
        means.append(Xg.mean(axis=0))
    if not means:
        return cats
    M = np.vstack(means)
    Z = linkage(pdist(M, metric="euclidean"), method="ward")
    order = leaves_list(Z).tolist()
    return [cats[i] for i in order]

# ---- Main ----
def main():
    args = parse_args()
    print("\n====== Step 3. Dotplot (bin/dotplot.py) ======\n")
    
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc
    import anndata as ad
    import pandas as pd

    h5ad_path = resolve_h5ad_path(args.input)
    markers_path = Path(args.reference).resolve()

    outdir = default_outdir_for_input(h5ad_path)
    out = Out(outdir)
    print(f"[scanpy-dotplot] Output dir: {outdir}")

    A: ad.AnnData = sc.read_h5ad(h5ad_path)

    groupby = choose_groupby(A, args.groupby)
    _ensure_categories_ordered(A, groupby)

    markers_df, type_col, gene_col, id_mode = load_markers_csv(
        markers_path, case_insensitive=args.case_insensitive
    )
    selected_ids, id_to_label = build_selection_and_labels(
        A, markers_df, type_col, gene_col, id_mode, case_insensitive=args.case_insensitive
    )

    try:
        w, h = (float(x) for x in args.figsize.split(","))
        figsize = (w, h)
    except Exception:
        figsize = (12, 10)

    std_scale = args.standard_scale
    if isinstance(std_scale, str) and std_scale.lower() in ("none", "null", "false"):
        std_scale = None

    A_sub, selected_labels = _make_relabeled_subset(A, selected_ids, id_to_label)

    # 1) Base tables
    base_var_labels = A_sub.var_names.astype(str).tolist()
    _export_dot_values(A_sub, groupby, base_var_labels, out=out, suffix="")  # dot_mean, dot_frac

    # 2) Dendrogram-ordered dotplot + CSV (dendro order)
    try:
        _ensure_numeric_X(A_sub)
        import scanpy as sc
        sc.tl.dendrogram(A_sub, groupby=groupby)
        ordered = _get_dendrogram_order(A_sub, groupby)
        _apply_category_order(A_sub, groupby, ordered)
        _export_dot_values(A_sub, groupby, A_sub.var_names.astype(str).tolist(), out=out, suffix="dendro")
        _run_dotplot(
            A_sub, selected_labels, groupby, std_scale,
            args.dot_min, args.dot_max, args.cmap, figsize,
            out=out, out_stem="dotplot_dendro", dendro=True
        )
    except Exception as e:
        sys.stderr.write(f"[scanpy-dotplot] Warning: dendrogram failed ({e}). Using manual order.\n")
        try:
            ordered = _manual_dendro_order(A_sub, groupby)
            _apply_category_order(A_sub, groupby, ordered)
            _export_dot_values(A_sub, groupby, A_sub.var_names.astype(str).tolist(), out=out, suffix="dendro")
            _run_dotplot(
                A_sub, selected_labels, groupby, std_scale,
                args.dot_min, args.dot_max, args.cmap, figsize,
                out=out, out_stem="dotplot_dendro", dendro=False
            )
        except Exception as e2:
            sys.stderr.write(f"[scanpy-dotplot] Warning: manual dendrogram failed ({e2}). Skipping dendrogram plot.\n")

    # 3) Auto-labeled plot + label tables
    try:
        groups, scores_df = _compute_signature_scores(A_sub, groupby, selected_labels)
        scores_df.index.name = "leiden"
        out.save_csv(scores_df.reset_index(), "label_scores.csv")

        raw_mapping = _assign_auto_labels(groups, scores_df, prefer_prefix="")
        rows = [(g, raw_mapping[g][0]) for g in groups]
        label_map_df = pd.DataFrame(rows, columns=["leiden", "type"])
        out.save_csv(label_map_df, "label_map.csv")

        auto_col = f"{groupby}_auto"
        display_map = {g: raw_mapping[g][1] for g in groups}
        A_sub.obs[auto_col] = A_sub.obs[groupby].astype(str).map(display_map).astype("category")

        cur_groups = _group_categories(A_sub, groupby)
        auto_levels = [display_map[g] for g in cur_groups if g in display_map]
        A_sub.obs[auto_col] = pd.Categorical(A_sub.obs[auto_col].astype(str), categories=auto_levels, ordered=True)

        _ensure_numeric_X(A_sub)
        try:
            sc.tl.dendrogram(A_sub, groupby=auto_col)
        except Exception as e:
            sys.stderr.write(f"[scanpy-dotplot] Warning: auto dendrogram precompute failed ({e}).\n")

        _run_dotplot(
            A_sub, selected_labels, auto_col, std_scale,
            args.dot_min, args.dot_max, args.cmap, figsize,
            out=out, out_stem="dotplot_dendro_auto", dendro=True
        )
    except Exception as e:
        sys.stderr.write(f"[scanpy-dotplot] Warning: auto-labeling failed ({e}). Skipping auto-labeled plot.\n")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[scanpy-dotplot] Error: {e}\n")
        sys.exit(1)
