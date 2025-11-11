#!/usr/bin/env python3

"""
Update per-sample h5ad with sample/leiden/celltype, build an AnnDataSet (.h5ads),
and export fragments as BED files grouped by celltype.sample (per-sample).

CLI:
  python annotate.py -b {/path/to/BASE_DIR}
"""

import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import scanpy as sc
import snapatac2 as snap


# ---- I/O helpers ----
def read_label_map(label_csv: Path) -> Dict[str, str]:
    df = pd.read_csv(label_csv)
    if "leiden" not in df.columns or "type" not in df.columns:
        raise ValueError(f"'leiden' and 'type' columns are required in: {label_csv}")
    df = df[["leiden", "type"]].copy()
    df["leiden"] = df["leiden"].astype(str)
    df["type"] = df["type"].astype(str)
    df = df.drop_duplicates(subset=["leiden"], keep="first")
    return dict(zip(df["leiden"], df["type"]))


def load_annot_meta(annot_h5ad: Path) -> pd.DataFrame:
    ad = sc.read_h5ad(annot_h5ad)
    missing = [c for c in ("sample", "leiden") if c not in ad.obs.columns]
    if missing:
        raise KeyError(f"Missing columns in annot h5ad: {', '.join(missing)}")

    meta = ad.obs[["sample", "leiden"]].copy()
    meta = meta[~meta.index.duplicated(keep="first")].copy()
    meta.index = meta.index.astype(str)
    meta["sample"] = meta["sample"].astype(str)
    meta["leiden"] = meta["leiden"].astype(str)

    try:
        ad.file.close()
    except Exception:
        pass
    del ad
    return meta

def infer_sample_name(h5ad_path: Path) -> str:
    stem = h5ad_path.stem
    return stem.split("_")[0] if "_" in stem else stem

# ---- Core lib ----
def update_one_h5ad(
    in_path: Path,
    sample_name: str,
    annot_meta: pd.DataFrame,
    lmap: Dict[str, str],
    out_dir: Path,
) -> Tuple[str, Path, Dict[str, int]]:
    ad = sc.read_h5ad(in_path)

    # Ensure unique barcodes and string index
    ad = ad[~ad.obs.index.duplicated(keep="first")].copy()
    ad.obs.index = ad.obs.index.astype(str)

    # Align annot meta to current sample barcodes
    meta_s = annot_meta[annot_meta["sample"] == sample_name]
    aligned = meta_s.reindex(ad.obs.index)

    # Leiden string; NA -> "NA"
    leiden_series = aligned["leiden"].astype(str).fillna("NA")

    # Celltype via mapping
    celltype_series = leiden_series.map(lambda x: lmap.get(x, "Unknown")).astype(str)

    if "sample" in ad.obs.columns:
        sample_series = ad.obs["sample"].astype(str).fillna(sample_name)
    else:
        sample_series = pd.Series(sample_name, index=ad.obs.index, dtype="string").astype(str)

    # Write obs
    ad.obs["sample"] = pd.Categorical(sample_series.values)
    ad.obs["leiden"] = pd.Categorical(leiden_series.values)
    ad.obs["celltype"] = pd.Categorical(celltype_series.values)
    ad.obs["sample_cluster"] = pd.Categorical(
        (ad.obs["celltype"].astype(str) + "." + ad.obs["sample"].astype(str)).values
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"update_{sample_name}.h5ad"
    ad.write(out_path)

    total = int(ad.shape[0])
    matched = int(aligned["leiden"].notna().sum())
    unknown = int((ad.obs["celltype"].astype(str) == "Unknown").sum())

    try:
        ad.file.close()
    except Exception:
        pass
    del ad

    return sample_name, out_path, {
        "cells_total": total,
        "cells_matched_leiden": matched,
        "cells_celltype_unknown": unknown,
    }

# ---- Make AnnDataSet and exporting BED ----
def build_h5ads_and_export_bed(
    updated_file_paths: List[Tuple[str, Path]],
    bed_dir: Path,
    h5ads_path: Path,
):
    dataset = snap.AnnDataSet(adatas=updated_file_paths, filename=str(h5ads_path))

    labels = []
    for sample_name, h5ad_path in updated_file_paths:
        ad = sc.read_h5ad(h5ad_path)
        for col in ("celltype", "sample", "sample_cluster"):
            if col not in ad.obs:
                raise KeyError(f"'{col}' not found in {h5ad_path}")
        s = ad.obs["sample_cluster"].astype(str).fillna("Unknown")
        labels.append(s)
        try:
            ad.file.close()
        except Exception:
            pass
        del ad

    groupby_series = pd.concat(labels)  # index: concatenated barcodes

    # Export BED grouped by celltype.sample
    bed_dir.mkdir(parents=True, exist_ok=True)
    # Try to get filenames as <type>.<sample>.bed
    snap.ex.export_fragments(
        dataset,
        groupby=groupby_series,
        out_dir=str(bed_dir),
        prefix="",
        suffix=".bed",
    )

    groups = sorted(groupby_series.unique())
    for g in groups:
        candidates = [
            bed_dir / f"{g}.bed",
            bed_dir / f"_{g}.bed",
            bed_dir / f"{g}_fragments.bed",
            bed_dir / f"_{g}_fragments.bed",
            bed_dir / f"ctype_{g}.bed",
            bed_dir / f"ctype_{g}_fragments.bed",
            bed_dir / f"mcluster_{g}.bed",
            bed_dir / f"mcluster_{g}_fragments.bed",
        ]
        target = bed_dir / f"{g}.bed"
        for c in candidates:
            if c.exists():
                if c != target:
                    try:
                        c.replace(target)
                    except Exception:
                        # Fallback: copy then remove
                        import shutil
                        shutil.copy2(c, target)
                        try:
                            c.unlink()
                        except Exception:
                            pass
                break

    counts = groupby_series.value_counts().sort_index()
    print("[Export] Groups (celltype.sample) and cell counts:")
    print(counts.to_string())

# ---- Main ----
def main():
    parser = argparse.ArgumentParser(description="Update per-sample h5ad and export BED by celltype.sample.")
    parser.add_argument("-b", "--base-dir", required=True, help="Base directory path")
    args = parser.parse_args()

    base = Path(args.base_dir).resolve()
    filter_dir = base / "Filter_results"
    annot_h5ad = base / "Annot_results" / "annot_gene_activity.h5ad"
    label_csv = base / "Annot_results" / "dotplots" / "label_map.csv"

    out_dir_h5ad = base / "Annot_results" / "update_h5ad"
    out_h5ads = out_dir_h5ad / "merged_updated.h5ads"
    out_dir_bed = base / "Annot_results" / "bed"

    # Sanity checks
    if not filter_dir.is_dir():
        raise FileNotFoundError(f"Not found: {filter_dir}")
    if not annot_h5ad.is_file():
        raise FileNotFoundError(f"Not found: {annot_h5ad}")
    if not label_csv.is_file():
        raise FileNotFoundError(f"Not found: {label_csv}")

    # Load resources
    annot_meta = load_annot_meta(annot_h5ad)
    lmap = read_label_map(label_csv)

    # Collect *.h5ad under Filter_results
    h5ads = sorted([p for p in filter_dir.glob("*.h5ad") if p.is_file()])
    if not h5ads:
        raise FileNotFoundError(f"No .h5ad files under: {filter_dir}")

    # Update each per-sample file
    summaries = []
    updated_file_paths: List[Tuple[str, Path]] = []
    for in_path in h5ads:
        sample_name = infer_sample_name(in_path)
        try:
            sample, out_path, summ = update_one_h5ad(
                in_path, sample_name, annot_meta, lmap, out_dir_h5ad
            )
            print(
                f"[OK] {sample}: total={summ['cells_total']}, "
                f"matched_leiden={summ['cells_matched_leiden']}, "
                f"celltype_unknown={summ['cells_celltype_unknown']} → {out_path}"
            )
            summaries.append(
                {
                    "input": str(in_path),
                    "output": str(out_path),
                    "sample": sample,
                    **summ,
                }
            )
            updated_file_paths.append((sample, out_path))
        except Exception as e:
            print(f"[ERR] {sample_name}: {in_path} → {e}")

    # Save summary CSV
    if summaries:
        df = pd.DataFrame(summaries)
        out_dir_h5ad.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir_h5ad / "update_summary.csv", index=False)
        print(f"Summary saved: {out_dir_h5ad/'update_summary.csv'}")
    else:
        print("No files updated.")
        return

    # Build .h5ads and export BED grouped by celltype.sample
    print(f"[Build] AnnDataSet → {out_h5ads}")
    build_h5ads_and_export_bed(
        updated_file_paths=updated_file_paths,
        bed_dir=out_dir_bed,
        h5ads_path=out_h5ads,
    )
    print(f"[Done] BED files are saved in: {out_dir_bed}")


if __name__ == "__main__":
    main()
