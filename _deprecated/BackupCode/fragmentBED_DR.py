# fragmentBED.py (project-tuned, robust): uses annot_gene_activity.h5ad by default to fill 'leiden'
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module=r"^anndata(\.|$)")
warnings.filterwarnings("ignore", category=FutureWarning, message=r".*Importing read_.* from `anndata` is deprecated.*")
import os
import glob
import time
from datetime import datetime
from typing import List, Tuple, Optional, Dict

import pandas as pd

try:
    import scanpy as sc
except Exception:
    sc = None

try:
    import anndata as ad
except Exception:
    ad = None

try:
    import snapatac2 as snap
except Exception:
    snap = None

# ----------------- small utils -----------------

def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _append_summary(summary_path: str, lines: List[str]) -> None:
    _ensure_dir(os.path.dirname(summary_path) or ".")
    with open(summary_path, "a", encoding="utf-8") as f:
        for line in lines:
            f.write(line.rstrip("\n") + "\n")


def _describe_dataframe(df: pd.DataFrame, max_rows: int = 30) -> List[str]:
    out = []
    df_small = df.copy()
    if df_small.shape[0] > max_rows:
        df_small = df_small.iloc[:max_rows, :]
        out.append(f"(showing top {max_rows} rows)")
    out.append(df_small.to_string())
    return out


def _list_h5ad(input_dir: str) -> List[str]:
    paths = sorted(glob.glob(os.path.join(input_dir, "*.h5ad")))
    # skip already-updated ones
    return [p for p in paths if not os.path.basename(p).startswith("updated_") and not os.path.basename(p).endswith("_updated.h5ad")]


def _read_h5ad(path: str):
    if sc is None:
        raise ImportError("scanpy is required. Please install scanpy.")
    return sc.read_h5ad(path)


def _write_h5ad(a, path: str) -> None:
    a.write(path)


# ----------------- CSV & mapping -----------------

_LEIDEN_ALIASES = ["leiden", "cluster", "clusters", "cluster_id", "louvain", "leiden_label"]
_LABEL_ALIASES  = ["best_label", "majority", "annotation", "celltype", "cell_type", "type", "label", "name"]
_BC_ALIASES     = ["barcode", "cell", "cell_id", "obs_name", "obs_names", "cellid"]


def _normalize_col(df: pd.DataFrame, target: str, aliases: List[str]) -> pd.DataFrame:
    if target in df.columns:
        return df
    lower = {c.lower(): c for c in df.columns}
    for a in aliases:
        if a in lower:
            return df.rename(columns={lower[a]: target})
    if target == "best_label":
        for c in df.columns:
            if c.lower().startswith("major"):
                return df.rename(columns={c: target})
    return df


def _load_cluster_csv(cluster_csv: Optional[str]) -> Optional[Dict[str, str]]:
    if not cluster_csv:
        return None
    if not os.path.exists(cluster_csv):
        print(f"[{_ts()}] Warning: cluster CSV not found: {cluster_csv}")
        return None
    try:
        df = pd.read_csv(cluster_csv, dtype=str)
        df = _normalize_col(df, "leiden", _LEIDEN_ALIASES)
        df = _normalize_col(df, "best_label", _LABEL_ALIASES)
        if "leiden" not in df.columns or "best_label" not in df.columns:
            print(f"[{_ts()}] Warning: CSV must contain columns 'leiden' and 'best_label'; got columns={list(df.columns)}")
            return None
        mapping = dict(zip(df["leiden"].astype(str), df["best_label"].astype(str)))
        print(f"[{_ts()}] cluster CSV loaded: {len(mapping)} leiden -> best_label mappings")
        return mapping
    except Exception as e:
        print(f"[{_ts()}] Warning: failed to read cluster CSV: {e}")
        return None


def _load_cells_csv(cells_csv: Optional[str]) -> Optional[pd.DataFrame]:
    if not cells_csv:
        return None
    if not os.path.exists(cells_csv):
        print(f"[{_ts()}] Warning: cells CSV not found: {cells_csv}")
        return None
    try:
        df = pd.read_csv(cells_csv, dtype=str)
        for c in df.columns:
            if c.lower() in _BC_ALIASES:
                df = df.set_index(c)
                break
        df = _normalize_col(df, "leiden", _LEIDEN_ALIASES)
        df.index = df.index.astype(str)
        df["leiden"] = df["leiden"].astype(str)
        return df[["leiden"]]
    except Exception as e:
        print(f"[{_ts()}] Warning: failed to read cells CSV: {e}")
        return None


def _load_aggr_h5ad(aggr_h5ad: Optional[str]) -> Optional[pd.DataFrame]:
    if not aggr_h5ad:
        return None
    if not os.path.exists(aggr_h5ad):
        print(f"[{_ts()}] Warning: aggr h5ad not found: {aggr_h5ad}")
        return None
    try:
        A = _read_h5ad(aggr_h5ad)
        if "leiden" not in A.obs.columns:
            print(f"[{_ts()}] Warning: aggr h5ad has no 'leiden' in obs")
            return None
        df = A.obs[["leiden"]].copy()
        df.index = df.index.astype(str)
        df["leiden"] = df["leiden"].astype(str)
        return df
    except Exception as e:
        print(f"[{_ts()}] Warning: failed to read aggr h5ad: {e}")
        return None


# ----------------- barcode match helpers -----------------

def _norm(bc: str) -> str:
    bc = str(bc).strip()
    if bc.endswith("-1"):
        return bc[:-2]
    return bc


def _variants(bc: str, tags: List[str]) -> List[str]:
    seps = [":", "_", ".", "-", "#", "/"]
    base = [bc, _norm(bc)]
    if not bc.endswith("-1"):
        base.append(bc + "-1")
    else:
        base.append(bc[:-2])
    out = list(dict.fromkeys(base))
    for tag in tags:
        for s in seps:
            out.append(f"{tag}{s}{bc}")
            out.append(f"{bc}{s}{tag}")
            nbc = _norm(bc)
            out.append(f"{tag}{s}{nbc}")
            out.append(f"{nbc}{s}{tag}")
    return list(dict.fromkeys(out))


def _fill_leiden_from_sources(a_obs: pd.DataFrame, stem: str, sample_value: str,
                              cells_df: Optional[pd.DataFrame], aggr_df: Optional[pd.DataFrame]) -> int:
    filled = 0
    need_mask = a_obs["leiden"].astype(str).isin(["", "NA", "nan"])
    if not need_mask.any():
        return 0
    idx_need = a_obs.index[need_mask].astype(str)
    tags = list(dict.fromkeys([sample_value, stem]))
    def try_map(src_df: Optional[pd.DataFrame]) -> int:
        nonlocal filled
        if src_df is None or src_df.empty:
            return 0
        src_index = set(src_df.index.astype(str))
        mapping = {}
        for bc in idx_need:
            for v in _variants(bc, tags):
                if v in src_index:
                    mapping[bc] = v
                    break
        if mapping:
            src = pd.Index(mapping.keys())
            dst = pd.Index([mapping[k] for k in mapping.keys()])
            a_obs.loc[src, "leiden"] = src_df.loc[dst, "leiden"].astype(str).values
            filled += len(mapping)
        return filled
    try_map(cells_df)
    try_map(aggr_df)
    return filled


# ----------------- main -----------------

def update_and_export(
    h5ad_dir: str,
    out_update_dir: str,
    out_bed_dir: str,
    # Back-compat (accepted but ignored):
    majority_key: str = "majority",
    leiden_key: str = "leiden",
    sample_key: str = "sample",
    summary_path: Optional[str] = None,
    dataset_name: str = "merged_updated.h5ads",
    bed_prefix: str = "group",
    bed_suffix: str = ".bed",
    convert_to_bed4: bool = False,   # unused; BED only
    # Defaults tuned for this project:
    cluster_csv: Optional[str] = "Test/Annot_results/Plots/cluster_assignments.csv",
    label_key: str = "best_label",
    cells_assignments_csv: Optional[str] = None,
    aggr_h5ad: Optional[str] = "Test/Annot_results/annot_gene_activity.h5ad",
) -> Dict[str, object]:

    t0 = time.time()
    if snap is None:
        raise ImportError("snapatac2 is required. Please install snapatac2.")
    if sc is None or ad is None:
        raise ImportError("scanpy/anndata are required. Please install scanpy and anndata.")

    _ensure_dir(out_update_dir)
    _ensure_dir(out_bed_dir)

    input_h5ads = _list_h5ad(h5ad_dir)
    if not input_h5ads:
        raise FileNotFoundError(f"No *.h5ad found under: {h5ad_dir}")

    # Load mappings/sources
    l2label = _load_cluster_csv(cluster_csv)
    cells_df = _load_cells_csv(cells_assignments_csv)
    aggr_df  = _load_aggr_h5ad(aggr_h5ad)

    updated: List[Tuple[str, str]] = []
    n_cells_total = 0
    print(f"[{_ts()}] Found {len(input_h5ads)} h5ad files under {h5ad_dir}")
    for p in input_h5ads:
        stem = os.path.splitext(os.path.basename(p))[0]
        try:
            a = _read_h5ad(p)
            a = a[~a.obs.index.duplicated(keep='first')].copy()
            n_cells_total += a.n_obs

            # Ensure required columns
            if sample_key not in a.obs:
                a.obs[sample_key] = stem
            a.obs[sample_key] = a.obs[sample_key].astype(str)

            if leiden_key not in a.obs:
                a.obs[leiden_key] = "NA"
                print(f"[{_ts()}] Info: {stem} has no '{leiden_key}' in obs; set to 'NA'")

            a.obs[leiden_key] = a.obs[leiden_key].astype(str)

            # Fill leiden if missing using aggr_h5ad or cells_assignments_csv
            filled = _fill_leiden_from_sources(
                a.obs, stem, a.obs[sample_key].iloc[0],
                cells_df=cells_df, aggr_df=aggr_df
            )
            if filled > 0:
                print(f"[{_ts()}] {stem}: filled leiden for {filled} cells from external assignments")

            # Map label by leiden (cluster-level mapping)
            if label_key not in a.obs:
                a.obs[label_key] = "Unknown"
            if l2label is not None:
                a.obs[label_key] = a.obs[leiden_key].map(l2label).fillna(a.obs[label_key]).astype(str)

            # sample_cluster for reference
            a.obs["sample_cluster"] = a.obs[label_key].astype(str) + "." + a.obs[sample_key].astype(str)

            # Save updated h5ad
            out_path = os.path.join(out_update_dir, f"updated_{stem}.h5ad")
            _write_h5ad(a, out_path)
            updated.append((stem, out_path))
            print(f"[{_ts()}] Updated: {out_path} (cells={a.n_obs})")
        except Exception as e:
            print(f"[{_ts()}] Error updating {p}: {e}")

    if not updated:
        raise RuntimeError("No files were updated; aborting.")

    # Build AnnDataSet
    dataset_path = os.path.join(out_update_dir, dataset_name)
    ds = None
    try:
        ds = snap.AnnDataSet(adatas=updated, filename=dataset_path)
        print(f"[{_ts()}] AnnDataSet created: {dataset_path}")

        # Group labels: use label_key, fallback to 'leidenX'
        labels = []
        groups_seen = set()
        for name, path in updated:
            a = _read_h5ad(path)
            lbl = a.obs[label_key].astype(str)
            lei = a.obs[leiden_key].astype(str).radd("leiden")
            label = lbl.mask(lbl.isin(['', 'Unknown', 'NA', 'nan']), lei)
            grp = label + "_" + a.obs[sample_key].astype(str)
            groups_seen.update(grp.unique().tolist())
            labels.append(grp.astype(str))
        groupby_series = pd.concat(labels)

        print(f"[{_ts()}] Exporting fragments by ({label_key}|leiden, sample) to: {out_bed_dir}")
        snap.ex.export_fragments(
            ds,
            groupby=groupby_series,
            out_dir=out_bed_dir,
            prefix="group",
            suffix=".bed",
        )
    finally:
        try:
            if ds is not None:
                ds.close()
        except Exception:
            pass

    # Normalize bed filenames to '<label>.bed'
    produced = sorted(glob.glob(os.path.join(out_bed_dir, "*.bed")))
    renamed = []
    for bp in produced:
        if not os.path.isfile(bp):
            continue
        base = os.path.basename(bp)
        root, _ = os.path.splitext(base)
        match = None
        for g in sorted(groups_seen, key=len, reverse=True):
            if g in root:
                match = g
                break
        if match:
            desired = f"{match}.bed"
            dst = os.path.join(out_bed_dir, desired)
            if bp != dst:
                try:
                    os.replace(bp, dst)
                    renamed.append(dst)
                except Exception:
                    pass

    # Summary
    print(f"[{_ts()}] Computing summary tables...")
    all_list, keys = [], []
    for name, path in updated:
        a = _read_h5ad(path)
        for col in [leiden_key, label_key, sample_key]:
            if col not in a.obs:
                raise KeyError(f"{os.path.basename(path)} missing '{col}'")
            a.obs[col] = a.obs[col].astype(str)
        all_list.append(a)
        keys.append(name)
    ad_all = ad.concat(all_list, join="outer", label="sample", keys=keys, index_unique=None)
    obs = ad_all.obs

    cluster_counts = obs[leiden_key].value_counts().sort_index()
    label_counts = obs[label_key].value_counts()
    cluster_label = pd.crosstab(obs[leiden_key], obs[label_key])
    sample_label = pd.crosstab(obs[sample_key], obs[label_key])

    if summary_path:
        lines = []
        lines.append(f"[{_ts()}] Post-label update summary")
        lines.append(f"Updated h5ad dir: {out_update_dir}")
        lines.append(f"BED dir: {out_bed_dir}")
        lines.append(f"AnnDataSet: {dataset_path}")
        lines.append(f"Files updated: {len(updated)}")
        lines.append(f"Total cells (approx): {n_cells_total}")
        lines.append("— Cluster counts (leiden) —")
        lines += _describe_dataframe(cluster_counts.to_frame("cells"))
        lines.append("— Label counts ({label_key}) —".format(label_key=label_key))
        lines += _describe_dataframe(label_counts.to_frame("cells"))
        lines.append("— Cluster × Label —")
        lines += _describe_dataframe(cluster_label)
        lines.append("— Sample × Label —")
        lines += _describe_dataframe(sample_label)
        elapsed = time.time() - t0
        lines.append(f"Elapsed: {elapsed:.1f}s")
        lines.append("")
        _append_summary(summary_path, lines)
        print(f"[{_ts()}] Summary appended to: {summary_path}")

    return {
        "updated_h5ads": [p for _, p in updated],
        "dataset_path": dataset_path,
        "bed_dir": out_bed_dir,
        "n_cells_total": n_cells_total,
    }


def cli():
    import argparse
    p = argparse.ArgumentParser(description="Update per-sample h5ad: fill 'leiden' (aggr_h5ad/cells CSV), map to label via cluster CSV, export BED per (label,sample).")
    p.add_argument("--h5ad_dir", required=True, help="Input directory containing per-sample *.h5ad")
    p.add_argument("--out_update_dir", default="Test/Annot_results/update_h5ad", help="Directory to write updated_*.h5ad")
    p.add_argument("--out_bed_dir", default="Test/Annot_results/bed", help="Directory to write BED files")
    p.add_argument("--summary", default="Test/Annot_results/summary.txt", help="Path to append summary (summary.txt)")
    p.add_argument("--leiden_key", default="leiden")
    p.add_argument("--sample_key", default="sample")
    p.add_argument("--cluster_csv", default="Test/Annot_results/Plots/cluster_assignments.csv", help="CSV with columns [leiden,best_label]")
    p.add_argument("--label_key", default="best_label")
    p.add_argument("--cells_assignments_csv", default=None, help="Optional: CSV [barcode,leiden]")
    p.add_argument("--aggr_h5ad", default="Test/Annot_results/annot_gene_activity.h5ad", help="Aggregated h5ad with obs['leiden']")
    # Back-compat flag; ignored
    p.add_argument("--majority_key", default="majority")
    args = p.parse_args()

    update_and_export(
        h5ad_dir=args.h5ad_dir,
        out_update_dir=args.out_update_dir,
        out_bed_dir=args.out_bed_dir,
        summary_path=args.summary,
        leiden_key=args.leiden_key,
        sample_key=args.sample_key,
        cluster_csv=args.cluster_csv,
        label_key=args.label_key,
        cells_assignments_csv=args.cells_assignments_csv,
        aggr_h5ad=args.aggr_h5ad,
        majority_key=args.majority_key,  # accepted but ignored
    )


if __name__ == "__main__":
    cli()
