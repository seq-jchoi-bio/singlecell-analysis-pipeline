#!/usr/bin/env python3

import os
import sys
import math
import gzip
import numpy as np
import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import snapatac2 as snap
from anndata import AnnData
from typing import Optional, List, Dict, Tuple
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from datetime import datetime

__version__ = "2.2"

HEADER = "\n ================================================================================="
HEADER += "\n     Integrated scATAC-seq: Quality Control Pipeline"
HEADER += "\n     Version {0}".format(__version__)
HEADER += "\n     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================"

desc_txt = """
This tool provides an integrated quality-control workflow for single-cell ATAC-seq data.
It automatically processes fragment files to generate fragment-size distribution plots,
TSSE (TSS enrichment) metrics, per-sample QC figures (as a single combined grid),
and summary outputs. All results are saved as PNG/SVG files, plus a concise TXT summary
with per-sample percentile statistics (p5/p10/median/p90/p95).
"""

warnings.filterwarnings("ignore", category=FutureWarning, module="snapatac2")

# =========================
# User-configurable section (wrapper-friendly)
# =========================
chrom_sizes_path = None  # wrapper may inject CHROM_SIZES_PATH or chrom_sizes_path
annot_path        = None  # wrapper may inject GTF_PATH or annot_path
base_dir          = None  # wrapper may inject BASE_DIR
out_dir           = None  # wrapper or runtime will set to CWD/QC_results

# Fragment-size grid
max_recorded_size = 1000          # upper bound for fragment-size hist
subsample_lines   = None           # e.g., 2_000_000 for fast preview in fallback
grid_rows         = 3              # rows for both FSD grid and TSSE grid
dpi_save          = 300            # figure DPI

# =========================
# Utilities
# =========================
def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def save_fig(fig: Figure, out_base: str, dpi: int = dpi_save) -> None:
    png = f"{out_base}.png"
    svg = f"{out_base}.svg"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    Saved: {png}")
    print(f"    Saved: {svg}")

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
            chrom, size = parts[0], int(parts[1])
            sizes[chrom] = size
    if not sizes:
        raise ValueError(f"No entries parsed from chrom sizes: {path}")
    return sizes

def collect_fragment_files(base_dir: str) -> List[str]:
    files: List[str] = []
    for name in sorted(os.listdir(base_dir)):
        p = os.path.join(base_dir, name, "outs", "fragments.tsv.gz")
        if os.path.isfile(p):
            files.append(p)
    return files

def compute_size_hist_from_fragments(fragment_file: str,
                                     max_size: int,
                                     subsample: Optional[int] = None) -> np.ndarray:
    hist = np.zeros(max_size, dtype=np.int64)
    n = 0
    with gzip.open(fragment_file, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            size = end - start
            if 1 <= size < max_size:
                hist[size] += 1
            n += 1
            if subsample and n >= subsample:
                break
    return hist

def frag_size_distr_matplotlib(adata: AnnData,
                               fragment_file: str,
                               use_rep: str = "frag_size_distr",
                               max_recorded_size: int = 1000,
                               ax: Optional[plt.Axes] = None) -> None:
    dist = None
    try:
        if use_rep not in adata.uns or len(adata.uns[use_rep]) < max_recorded_size:
            snap.metrics.frag_size_distr(
                adata, add_key=use_rep, max_recorded_size=max_recorded_size
            )
        dist = adata.uns[use_rep]
    except Exception:
        hist = compute_size_hist_from_fragments(fragment_file, max_recorded_size,
                                                subsample=subsample_lines)
        dist = hist.tolist()

    x, y = zip(*enumerate(dist))
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x[1:], y[1:])
    ax.set_xlabel("Fragment size")
    ax.set_ylabel("Count")

def compute_tsse(adata: AnnData, gtf_path: str) -> None:
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

def _pct(arr: np.ndarray, p: float) -> float:
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.percentile(arr, p))

def _mean(arr: np.ndarray) -> float:
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.mean(arr))

def _sd(arr: np.ndarray) -> float:
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size <= 1:
        return float("nan")
    return float(np.std(arr, ddof=1))

def _tsse_plot_to_image(adata: AnnData, cutoff_line: float) -> Optional[np.ndarray]:
    try:
        import plotly.io as pio
        fig = None
        try:
            fig = snap.pl.tsse(adata, interactive=False, width=600, height=480, show=False)
        except TypeError:
            fig = snap.pl.tsse(adata, interactive=False, width=600, height=480)
        if hasattr(fig, "add_shape") and np.isfinite(cutoff_line):
            fig.add_shape(
                type="line", xref="paper", yref="y", x0=0, x1=1,
                y0=cutoff_line, y1=cutoff_line, line=dict(dash="dash", color="red")
            )
        png_bytes = pio.to_image(fig, format="png", scale=2)  # requires kaleido
        from io import BytesIO
        img = mpimg.imread(BytesIO(png_bytes), format="png")
        return img
    except Exception:
        return None

# =========================
# Main
# =========================
# Wrapper injection
def _apply_wrapper_injections_and_defaults() -> None:
    """Map uppercase injections (if any) and set safe defaults."""
    global base_dir, out_dir, chrom_sizes_path, annot_path, max_recorded_size, subsample_lines, grid_rows, dpi_save
    g = globals()
    # Map uppercase-injected names to lowercase vars used in this module
    if g.get("BASE_DIR") is not None:
        base_dir = g.get("BASE_DIR")
    if g.get("OUT_DIR") is not None:
        out_dir = g.get("OUT_DIR")
    if g.get("CHROM_SIZES_PATH") is not None and chrom_sizes_path is None:
        chrom_sizes_path = g.get("CHROM_SIZES_PATH")
    if g.get("GTF_PATH") is not None and annot_path is None:
        annot_path = g.get("GTF_PATH")
    # Also accept alternative attribute names
    if g.get("chrom_sizes_path") is not None:
        chrom_sizes_path = g.get("chrom_sizes_path")
    if g.get("annot_path") is not None:
        annot_path = g.get("annot_path")
    # Defaults
    if out_dir is None:
        out_dir = os.path.join(os.getcwd(), "QC_results")
    os.makedirs(out_dir, exist_ok=True)

def main():
    _apply_wrapper_injections_and_defaults()
    import matplotlib as mpl
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["figure.dpi"]   = dpi_save

    print(HEADER)
    print(desc_txt)
    # Pre-flight checks
    if annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip: use an uncompressed .gtf for 'annot_path'.")

    ensure_dir(out_dir)

    print("[1/5] Collecting fragment files ...", flush=True)
    fragment_files = collect_fragment_files(base_dir)
    if len(fragment_files) == 0:
        raise FileNotFoundError(f"No fragments.tsv.gz found under: {base_dir}")
    print(f"    Found {len(fragment_files)} samples.", flush=True)
    print("[2/5] Loading chrom.sizes ...", flush=True)
    chrom_sizes_dict = load_chrom_sizes(chrom_sizes_path)
    print(f"    Loaded {len(chrom_sizes_dict)} contigs.", flush=True)
    n = len(fragment_files)

    # (1) Fragment-size distribution grid -> saved directly under out_dir
    print("[3/5] Computing fragment-size distributions (FSD, grid) ...", flush=True)
    ncols = math.ceil(math.sqrt(n)) if n > 0 else 1
    nrows = math.ceil(n / ncols) if n > 0 else 1
    if grid_rows and grid_rows > 0:
        nrows = min(grid_rows, n)
        ncols = math.ceil(n / nrows)

    fig_fsd, axes_fsd = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes_list = np.atleast_1d(axes_fsd).ravel().tolist()
    iterable_fsd = fragment_files
    for i, fragment_file in enumerate(iterable_fsd):
        try:
            adata = snap.pp.import_fragments(
                fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
            )
            ax = axes_list[i]
            frag_size_distr_matplotlib(
                adata, fragment_file=fragment_file,
                max_recorded_size=max_recorded_size, ax=ax
            )
            sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
            ax.set_title(sample_name)
        except Exception:
            ax = axes_list[i]
            ax.text(0.5, 0.5, "(error)", ha="center", va="center"); ax.axis("off")
    for j in range(len(fragment_files), len(axes_list)):
        fig_fsd.delaxes(axes_list[j])

    out_base = os.path.join(out_dir, f"fragment_size")
    save_fig(fig_fsd, out_base)
    # (2) TSSE per sample + violin (compute GLOBAL MEDIAN TSSE for cutoff line)
    print("[4/5] Computing TSSE per sample + violin plot ...", flush=True)
    all_tsse_scores_raw: List[List[float]] = []
    sample_labels: List[str] = []
    kept_counts = 0
    skipped = []
    iterable_tsse = fragment_files
    for i, fragment_file in enumerate(iterable_tsse):
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
        try:
            adata = snap.pp.import_fragments(
                fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
            )
            compute_tsse(adata, annot_path)
            tsse_vals = np.asarray(adata.obs.get("tsse"), dtype=float)
            tsse_vals = tsse_vals[np.isfinite(tsse_vals)]
            tsse_vals = tsse_vals[tsse_vals > 0]

            if tsse_vals.size == 0:
                skipped.append(sample_name)
            else:
                all_tsse_scores_raw.append(tsse_vals.tolist())
                sample_labels.append(sample_name)
                kept_counts += 1
        except Exception:
            skipped.append(sample_name)
    if kept_counts > 0:
        all_vec = np.concatenate([np.asarray(v, dtype=float) for v in all_tsse_scores_raw])
        global_median_tsse = float(np.median(all_vec)) if all_vec.size > 0 else 6.0
    else:
        global_median_tsse = 6.0

    if kept_counts > 0:
        x_values = np.arange(kept_counts)
        fig_violin = plt.figure(figsize=(12, 8))
        plt.violinplot(all_tsse_scores_raw, positions=x_values, showmeans=False, showmedians=True)
        plt.axhline(y=global_median_tsse, linestyle="--",
                    label=f"Global median TSSE = {global_median_tsse:.2f}", color="red")
        plt.xlabel("Sample")
        plt.ylabel("TSS Enrichment Score")
        plt.title("Violin Plot of TSSE Across Samples")
        plt.xticks(ticks=x_values, labels=sample_labels, rotation=30, ha="right")
        plt.legend()
        plt.tight_layout()

        out_base = os.path.join(out_dir, f"TSSE_violin")
        save_fig(fig_violin, out_base)
    else:
        print("    [warn] No sample had valid TSSE values. Violin plot skipped.")
    if skipped:
        print(f"    Skipped samples (no/invalid TSSE): {', '.join(skipped)}")

    # (3) Combined per-sample TSSE grid
    print("[5/5] Rendering combined per-sample TSSE grid (single figure) ...", flush=True)
    ncols = math.ceil(math.sqrt(n)) if n > 0 else 1
    nrows = math.ceil(n / ncols) if n > 0 else 1
    if grid_rows and grid_rows > 0:
        nrows = min(grid_rows, n)
        ncols = math.ceil(n / nrows)

    fig_grid, axes_grid = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows), squeeze=False)
    iterable_ps = fragment_files
    for i, fragment_file in enumerate(iterable_ps):
        r, c = divmod(i, ncols)
        ax = axes_grid[r, c]
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
        try:
            adata = snap.pp.import_fragments(
                fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
            )
            compute_tsse(adata, annot_path)

            # Prefer Plotly path for the classic KDE
            img = _tsse_plot_to_image(adata, cutoff_line=global_median_tsse)
            if img is not None:
                ax.imshow(img)
                ax.axis("off")
                ax.set_title(sample_name)
            else:
                # Fallback: minimal scatter to ensure something renders
                ts = np.asarray(adata.obs.get("tsse"), dtype=float)
                nf = np.asarray(adata.obs.get("n_fragment"), dtype=float)
                m = np.isfinite(ts) & np.isfinite(nf) & (ts > 0) & (nf > 0)
                if m.sum() == 0:
                    ax.text(0.5, 0.5, "(no valid TSSE)", ha="center", va="center")
                    ax.axis("off")
                else:
                    ax.scatter(nf[m], ts[m], s=6, alpha=0.25)
                    if np.isfinite(global_median_tsse):
                        ax.axhline(global_median_tsse, ls="--", color="red")
                    ax.set_xlabel("n_fragment"); ax.set_ylabel("TSSE")
                    ax.set_title(sample_name)
        except Exception:
            ax.text(0.5, 0.5, "(error)", ha="center", va="center")
            ax.axis("off")
    # Hide any unused axes
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        ax = axes_grid[r, c]
        ax.axis("off")

    plt.tight_layout()
    out_base = os.path.join(out_dir, f"TSSE_grid")
    save_fig(fig_grid, out_base)
    # (4) Concise QC summary (TXT/CSV)
    print("[*] Writing QC summary (TXT/CSV) ...", flush=True)
    import pandas as pd

    summary_lines = []
    summary_lines.append(HEADER)
    summary_lines.append("")
    summary_lines.append("===== QC SUMMARY =====")
    summary_lines.append(f"Generated at: {timestamp()}")
    summary_lines.append("")
    summary_lines.append(f"Base dir     : {base_dir}")
    summary_lines.append(f"chrom.sizes  : {chrom_sizes_path}")
    summary_lines.append(f"GTF          : {annot_path}")
    summary_lines.append("")

    # Build per-sample table
    rows = []
    headers = [
        "sample",
        "cells_total",
        "nFrag_min","nFrag_p5","nFrag_p10","nFrag_med","nFrag_p90","nFrag_p95","nFrag_max","nFrag_mean","nFrag_sd",
        "TSSE_min","TSSE_p5","TSSE_p10","TSSE_med","TSSE_p90","TSSE_p95","TSSE_max","TSSE_mean","TSSE_sd"
    ]

    for fragment_file in fragment_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
        try:
            adata_tmp = snap.pp.import_fragments(
                fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
            )
            compute_tsse(adata_tmp, annot_path)
            tsse_vals = np.asarray(adata_tmp.obs.get("tsse"), dtype=float)
            tsse_vals = tsse_vals[np.isfinite(tsse_vals) & (tsse_vals > 0)]
            nfrag_vals = np.asarray(adata_tmp.obs.get("n_fragment"), dtype=float)
            nfrag_vals = nfrag_vals[np.isfinite(nfrag_vals) & (nfrag_vals > 0)]
            cells_total = int(min(len(tsse_vals), len(nfrag_vals))) if (tsse_vals.size>0 and nfrag_vals.size>0) else int(max(len(tsse_vals), len(nfrag_vals)))

            row = {
                "sample": sample_name,
                "cells_total": cells_total,
                "nFrag_min": float(np.min(nfrag_vals)) if nfrag_vals.size>0 else float("nan"),
                "nFrag_p5":  _pct(nfrag_vals, 5),
                "nFrag_p10": _pct(nfrag_vals, 10),
                "nFrag_med": _pct(nfrag_vals, 50),
                "nFrag_p90": _pct(nfrag_vals, 90),
                "nFrag_p95": _pct(nfrag_vals, 95),
                "nFrag_max": float(np.max(nfrag_vals)) if nfrag_vals.size>0 else float("nan"),
                "nFrag_mean": _mean(nfrag_vals),
                "nFrag_sd":   _sd(nfrag_vals),
                "TSSE_min": float(np.min(tsse_vals)) if tsse_vals.size>0 else float("nan"),
                "TSSE_p5":  _pct(tsse_vals, 5),
                "TSSE_p10": _pct(tsse_vals, 10),
                "TSSE_med": _pct(tsse_vals, 50),
                "TSSE_p90": _pct(tsse_vals, 90),
                "TSSE_p95": _pct(tsse_vals, 95),
                "TSSE_max": float(np.max(tsse_vals)) if tsse_vals.size>0 else float("nan"),
                "TSSE_mean": _mean(tsse_vals),
                "TSSE_sd":   _sd(tsse_vals),
            }
            rows.append(row)
        except Exception as e:
            summary_lines.append(f"[WARN] Skipped {sample_name} ({e})")

    if rows:
        import pandas as pd
        df = pd.DataFrame(rows, columns=headers).set_index("sample").sort_index()
        summary_lines.append("Per-sample QC table (TSV columns):")
        summary_lines.append("\t".join(["sample"] + headers[1:]))
        for s, r in df.sort_index().iterrows():
            vals = [s] + [r[c] for c in headers[1:]]
            summary_lines.append("\t".join(str(v) for v in vals))
        summary_lines.append("")
        summary_lines.append(f"TOTAL cells (sum of cells_total) : {int(df['cells_total'].sum())}")
        summary_lines.append(f"Global median TSSE               : {global_median_tsse:.2f}")
        try:
            summary_lines.append(f"Global median nFrag              : {float(np.median(df['nFrag_med'].values[np.isfinite(df['nFrag_med'].values)])):.0f}")
        except Exception:
            pass
    else:
        summary_lines.append("No valid samples for summary.")

    out_summary_txt = os.path.join(out_dir, "qc_summary.txt")
    with open(out_summary_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))
    print(f"    Saved summary TXT: {out_summary_txt}")

    try:
        if rows:
            df = pd.DataFrame(rows, columns=headers).set_index("sample").sort_index()
            out_csv = os.path.join(out_dir, "qc_summary_table.csv")
            df.to_csv(out_csv, index=True)
            print(f"    Saved summary CSV: {out_csv}")
    except Exception as e:
        print(f"    [warn] CSV export failed: {e}")
    print("Done. All figures and summaries saved under:", out_dir, flush=True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[Interrupted] Exiting.", file=sys.stderr)
        sys.exit(130)
