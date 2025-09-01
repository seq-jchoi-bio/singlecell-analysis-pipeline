#!/usr/bin/env python3
# Integrated scATAC-seq QC Pipeline (v1.10)
#  - Fragment-size distribution grid (PNG/SVG)
#  - TSSE violin across samples (PNG/SVG) with GLOBAL MEDIAN TSSE line
#  - Per-sample TSSE plots using snap.pl.tsse (classic KDE look) with robust detection,
#    PNG/SVG parity, and multiple fallbacks to avoid blank images
#  - TSSE mosaic (PNG/SVG)
#  - Concise QC summary (TXT)
#  - Percentile-based cutoff suggestions (CSV, optional)
# Layout:
#  - All outputs directly under QC_results/ (per-sample images under QC_results/per_sample_tsse/)

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

try:
    from tqdm.auto import tqdm
except Exception:
    tqdm = None  # fallback if tqdm is unavailable

__version__ = "1.10"

HEADER = "\n ================================================================================="
HEADER += "\n     Integrated scATAC-seq QC Pipeline"
HEADER += "\n     Version {0}".format(__version__)
HEADER += "\n     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================"

desc_txt = """
This tool provides an integrated quality-control workflow for single-cell ATAC-seq data.
It automatically processes fragment files to generate fragment-size distribution plots,
TSSE (TSS enrichment) metrics, per-sample QC figures, and summary mosaics.
All results are saved as publication-ready PNG and SVG files, plus a concise TXT summary
and optional CSV cutoff suggestions.
"""

# Silence FutureWarnings from plotting internals
warnings.filterwarnings("ignore", category=FutureWarning, module="snapatac2")

# =========================
# User-configurable section
# (Temp) genome관련 내용은 config로 빼서 human/mouse는 직접 다운받아 사용하고, custom genome은 경로를 지정하는 것으로 수정.
# (Temp) base_dir의 하위 폴더를 모두 인식.
# =========================
chrom_sizes_path = "work/genome.chrom.sizes"  # path to chrom.sizes
annot_path        = "work/genes.gtf"          # path to unzipped GTF
base_dir          = "examples"                    # contains <sample>/outs/fragments.tsv.gz
out_dir           = "results/QC_results"         # all outputs directly here

# Fragment-size grid
max_recorded_size = 1000          # upper bound for fragment-size hist
subsample_lines   = None          # e.g., 2_000_000 for fast preview in fallback
mosaic_rows       = 3
dpi_save          = 600
USE_TQDM          = False         # True to keep tqdm bars (progress will be duplicated)

# Optional: write percentile-based cutoff suggestions CSV
WRITE_SUGGESTIONS = True

# =========================
# Utilities
# Kernel density plot에서 레스터라이즈가 아닌 벡터로 표현하고자 했는데, 레스터라이즈한 이유가 있을 것으로 보임(소형 코멘트), 그래서 그냥 레스터라이즈 kde사용했음
# =========================
def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def save_fig(fig: Figure, out_base: str, dpi: int = dpi_save) -> None:
    """Save both PNG and SVG consistently."""
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
    """Collect <base>/<sample>/outs/fragments.tsv.gz (non-recursive)."""
    files: List[str] = []
    for name in sorted(os.listdir(base_dir)):
        p = os.path.join(base_dir, name, "outs", "fragments.tsv.gz")
        if os.path.isfile(p):
            files.append(p)
    return files

def compute_size_hist_from_fragments(fragment_file: str,
                                     max_size: int,
                                     subsample: Optional[int] = None) -> np.ndarray:
    """Directly parse fragments to make a size histogram (fallback)."""
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
    """Try SnapATAC2 metrics; fallback to direct parsing if needed."""
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
    """Robust TSSE call across SnapATAC2 versions."""
    try:
        # Preferred: pass GTF path directly
        snap.metrics.tsse(adata, gtf_path)
        return
    except TypeError:
        # Older builds: accept an annotation object
        if hasattr(snap.genome, "read_gtf"):
            try:
                ann = snap.genome.read_gtf(gtf_path)
                snap.metrics.tsse(adata, ann)
                return
            except Exception:
                pass
        # Try alternative kw names
        for kw in ("gtf_file", "annotation_file", "gtf"):
            try:
                snap.metrics.tsse(adata, **{kw: gtf_path})
                return
            except TypeError:
                continue
    raise RuntimeError("TSSE call failed: incompatible SnapATAC2 API for GTF input.")

# Robust KDE plotter using snap.pl.tsse
def _axes_artist_counts(ax: Axes) -> int:
    """Heuristic: count data-bearing artists to detect if drawing occurred."""
    return len(ax.images) + len(ax.collections) + len(ax.lines)

def _figure_has_content(fig: Figure) -> bool:
    for a in fig.get_axes():
        if _axes_artist_counts(a) > 0:
            return True
    return False

def tsse_plot_snap(adata: AnnData, cutoff_line: float) -> Figure:
    """
    Render per-sample TSSE with the classic KDE look via snap.pl.tsse.
    Robustly detects whether drawing happened on our Axes or on a new Figure.
    Fallbacks to PNG-bytes or scatter to avoid blanks.
    """
    # 1) Validate TSSE quickly
    y = np.asarray(adata.obs.get("tsse"), dtype=float)
    x = np.asarray(adata.obs.get("n_fragment"), dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if m.sum() == 0:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(1, 1, 1)
        ax.text(0.5, 0.5, "No valid TSSE (check GTF/chrom names)", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        return fig

    # 2) Try to draw on our Axes
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(1, 1, 1)
    pre = _axes_artist_counts(ax)
    ret = None
    drew_here = False
    try:
        try:
            ret = snap.pl.tsse(adata, ax=ax, interactive=False)
        except TypeError:
            plt.sca(ax)
            ret = snap.pl.tsse(adata, interactive=False)
    except Exception:
        ret = None

    post = _axes_artist_counts(ax)
    if post > pre:
        drew_here = True

    # If it drew on our Axes, add cutoff and return
    if drew_here:
        if np.isfinite(cutoff_line):
            ax.axhline(cutoff_line, ls="--", color="red", label=f"TSSE = {cutoff_line:.2f}")
            ax.legend()
        fig.tight_layout()
        return fig

    # 3) If snap created a NEW figure (ret is None commonly), capture it
    try:
        fig2 = plt.gcf()
        if fig2 is not fig and isinstance(fig2, Figure) and _figure_has_content(fig2):
            # Add cutoff to the first content axis if possible
            axs = fig2.get_axes()
            if axs and np.isfinite(cutoff_line):
                a0 = axs[0]
                a0.axhline(cutoff_line, ls="--", color="red", label=f"TSSE = {cutoff_line:.2f}")
                try:
                    a0.legend()
                except Exception:
                    pass
            fig2.tight_layout()
            return fig2
    except Exception:
        pass

    # 4) Handle direct return types
    try:
        if isinstance(ret, Figure) and _figure_has_content(ret):
            if ret.axes and np.isfinite(cutoff_line):
                a0 = ret.axes[0]
                a0.axhline(cutoff_line, ls="--", color="red", label=f"TSSE = {cutoff_line:.2f}")
                try:
                    a0.legend()
                except Exception:
                    pass
            ret.tight_layout()
            return ret

        if isinstance(ret, Axes) and _axes_artist_counts(ret) > 0:
            fig_tsse = ret.figure
            if np.isfinite(cutoff_line):
                ret.axhline(cutoff_line, ls="--", color="red", label=f"TSSE = {cutoff_line:.2f}")
                try:
                    ret.legend()
                except Exception:
                    pass
            fig_tsse.tight_layout()
            return fig_tsse
    except Exception:
        pass

    # 5) IPython Image (PNG bytes) fallback -> show bytes in a Figure
    from io import BytesIO
    png_bytes = None
    try:
        from IPython.display import Image as IPyImage
        if isinstance(ret, IPyImage) and getattr(ret, "data", None):
            png_bytes = ret.data
    except Exception:
        png_bytes = None

    fig_fb = plt.figure(figsize=(6, 5))
    ax_fb = fig_fb.add_subplot(1, 1, 1)
    if png_bytes:
        try:
            import PIL.Image as PILImage
            img = PILImage.open(BytesIO(png_bytes))
            ax_fb.imshow(img); ax_fb.axis("off")
            fig_fb.tight_layout()
            return fig_fb
        except Exception:
            try:
                img = mpimg.imread(BytesIO(png_bytes), format="png")
                ax_fb.imshow(img); ax_fb.axis("off")
                fig_fb.tight_layout()
                return fig_fb
            except Exception:
                pass

    # 6) Last resort: minimal scatter so that something shows
    xx, yy = x[m], y[m]
    ax_fb.scatter(xx, yy, s=6, alpha=0.25)
    if np.isfinite(cutoff_line):
        ax_fb.axhline(cutoff_line, ls="--", color="red", label=f"TSSE = {cutoff_line:.2f}")
        try:
            ax_fb.legend()
        except Exception:
            pass
    ax_fb.set_xlabel("n_fragment"); ax_fb.set_ylabel("TSSE"); ax_fb.set_title("TSSE QC (fallback)")
    fig_fb.tight_layout()
    return fig_fb

def print_progress(pct: float, msg: str = ""):
    pct = max(0, min(100, pct))
    print(f"Progress: {pct:6.2f}%  {msg[:80]:<80}")

class GlobalProgress:
    def __init__(self):
        self.total_units = 0
        self.done_units  = 0
    def add_units(self, n: int):
        self.total_units += int(n)
    def inc(self, n: int = 1, msg: str = ""):
        self.done_units += int(n)
        if self.total_units > 0:
            pct = (self.done_units / self.total_units) * 100.0
            print_progress(pct, msg=msg)

# =========================
# Main
# =========================
def main():
    # Configure SVG text/vector consistency
    import matplotlib as mpl
    mpl.rcParams["svg.fonttype"] = "none"   # keep text as vectors
    mpl.rcParams["figure.dpi"]   = dpi_save

    print(HEADER)
    print(desc_txt)

    GP = GlobalProgress()
    print_progress(0.0, "Initializing")

    # Pre-flight checks
    if annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip: use an uncompressed .gtf for 'annot_path'.")

    ensure_dir(out_dir)
    per_sample_dir = os.path.join(out_dir, "per_sample_tsse"); ensure_dir(per_sample_dir)

    print("[1/5] Collecting fragment files ...", flush=True)
    fragment_files = collect_fragment_files(base_dir)
    if len(fragment_files) == 0:
        raise FileNotFoundError(f"No fragments.tsv.gz found under: {base_dir}")
    print(f"    Found {len(fragment_files)} samples.", flush=True)
    GP.inc(1, "Collected fragment files")

    print("[2/5] Loading chrom.sizes ...", flush=True)
    chrom_sizes_dict = load_chrom_sizes(chrom_sizes_path)
    print(f"    Loaded {len(chrom_sizes_dict)} contigs.", flush=True)
    GP.inc(1, "Loaded chrom sizes")

    # Define total workload after discovering n
    n = len(fragment_files)
    GP.add_units(n)          # FSD loop
    GP.add_units(2)          # FSD figure save + overhead
    GP.add_units(n)          # TSSE per-sample gather for violin
    GP.add_units(2)          # Violin build + save
    GP.add_units(n)          # Per-sample TSSE plots
    GP.add_units(2)          # Mosaic build + save
    GP.add_units(2)          # Summary + suggestions

    # (1) Fragment-size distribution grid -> saved directly under out_dir
    print("[3/5] Computing fragment-size distributions (grid) ...", flush=True)
    ncols = math.ceil(math.sqrt(n)) if n > 0 else 1
    nrows = math.ceil(n / ncols) if n > 0 else 1

    fig_fsd, axes_fsd = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes_list = np.atleast_1d(axes_fsd).ravel().tolist()

    iterable_fsd = fragment_files
    if USE_TQDM and tqdm is not None:
        iterable_fsd = tqdm(fragment_files, desc="  FSD per sample")

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
        GP.inc(1, f"FSD {i+1}/{n}")

    for j in range(len(fragment_files), len(axes_list)):
        fig_fsd.delaxes(axes_list[j])

    out_base = os.path.join(out_dir, f"fragment_size_grid_{timestamp()}")
    save_fig(fig_fsd, out_base)
    GP.inc(2, "Saved FSD grid")

    # (2) TSSE per sample + violin (compute GLOBAL MEDIAN TSSE for cutoff line)
    print("[4/5] Computing TSSE per sample + violin plot ...", flush=True)
    all_tsse_scores_raw: List[List[float]] = []
    sample_labels: List[str] = []
    kept_counts = 0
    skipped = []

    iterable_tsse = fragment_files
    if USE_TQDM and tqdm is not None:
        iterable_tsse = tqdm(fragment_files, desc="  TSSE per sample")

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
        GP.inc(1, f"TSSE {i+1}/{n}")

    # Compute GLOBAL MEDIAN TSSE (fallback to 6.0 if none)
    if kept_counts > 0:
        all_vec = np.concatenate([np.asarray(v, dtype=float) for v in all_tsse_scores_raw])
        global_median_tsse = float(np.median(all_vec)) if all_vec.size > 0 else 6.0
    else:
        global_median_tsse = 6.0

    if kept_counts > 0:
        x_values = np.arange(kept_counts)
        fig_violin = plt.figure(figsize=(12, 8))
        plt.violinplot(all_tsse_scores_raw, positions=x_values, showmeans=False, showmedians=True)
        plt.axhline(y=global_median_tsse, linestyle="--", label=f"Global median TSSE = {global_median_tsse:.2f}", color="red")
        plt.xlabel("Sample")
        plt.ylabel("TSS Enrichment Score")
        plt.title("Violin Plot of TSSE Across Samples")
        plt.xticks(ticks=x_values, labels=sample_labels, rotation=30, ha="right")
        plt.legend()
        plt.tight_layout()

        out_base = os.path.join(out_dir, f"tsse_violin_{timestamp()}")
        save_fig(fig_violin, out_base)
        GP.inc(2, "Saved TSSE violin")
    else:
        print("    [warn] No sample had valid TSSE values. Violin plot skipped.")
        GP.inc(2, "Violin skipped (no valid TSSE)")

    if skipped:
        print(f"    Skipped samples (no/invalid TSSE): {', '.join(skipped)}")

    # (3) Per-sample TSSE plots via snap.pl.tsse (KDE-like), PNG/SVG parity
    print("[5/5] Rendering per-sample TSSE plots (PNG & SVG) ...", flush=True)
    saved_imgs: List[Optional[str]] = []

    iterable_ps = fragment_files
    if USE_TQDM and tqdm is not None:
        iterable_ps = tqdm(fragment_files, desc="  TSSE plot per sample")

    for i, fragment_file in enumerate(iterable_ps):
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
        if not os.path.exists(fragment_file):
            saved_imgs.append(None)
            GP.inc(1, f"Per-sample TSSE {i+1}/{n}")
            continue
        try:
            adata = snap.pp.import_fragments(
                fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
            )
            compute_tsse(adata, annot_path)

            ts = np.asarray(adata.obs.get("tsse"), dtype=float)
            ts = ts[np.isfinite(ts)]; ts = ts[ts > 0]
            if ts.size == 0:
                saved_imgs.append(None)
                GP.inc(1, f"Per-sample TSSE {i+1}/{n}")
                continue

            base = os.path.join(per_sample_dir, f"{sample_name}_tsse_{timestamp()}")
            # Do NOT close all figures here; snap may draw on current/new figure
            fig_tsse = tsse_plot_snap(adata, cutoff_line=global_median_tsse)
            save_fig(fig_tsse, base)  # saves .png and .svg
            saved_imgs.append(base + ".png")

        except Exception:
            saved_imgs.append(None)
        GP.inc(1, f"Per-sample TSSE {i+1}/{n}")

    # Mosaic (only if any) -> saved directly under out_dir
    valid_imgs = [p for p in saved_imgs if p and os.path.exists(p)]
    if valid_imgs:
        num = len(valid_imgs)
        rows = min(mosaic_rows, num)
        cols = math.ceil(num / rows)
        fig_mosaic, axes_mosaic = plt.subplots(rows, cols, figsize=(5*cols, 5*rows), squeeze=False)

        for i in range(rows * cols):
            r, c = divmod(i, cols)
            ax = axes_mosaic[r, c]
            if i < num:
                img = mpimg.imread(valid_imgs[i])
                ax.imshow(img)
                title = os.path.basename(valid_imgs[i]).replace(".png", "")
                ax.set_title(title)
            else:
                ax.text(0.5, 0.5, "(missing)", ha="center", va="center")
            ax.axis("off")

        plt.tight_layout()
        out_base = os.path.join(out_dir, f"tsse_mosaic_{timestamp()}")
        save_fig(fig_mosaic, out_base)
        GP.inc(2, "Saved TSSE mosaic")
    else:
        print("    No valid per-sample TSSE images. Mosaic skipped.")
        GP.inc(2, "Mosaic skipped")

    # (4) Concise QC summary
    print("[*] Writing QC summary (TXT) ...", flush=True)
    import pandas as pd

    summary_lines = []
    summary_lines.append("===== QC SUMMARY =====")
    summary_lines.append(f"Generated at: {timestamp()}")
    summary_lines.append(f"Base dir     : {base_dir}")
    summary_lines.append(f"chrom.sizes  : {chrom_sizes_path}")
    summary_lines.append(f"GTF          : {annot_path}")
    summary_lines.append(f"Global median TSSE (used for lines) : {global_median_tsse:.2f}")
    summary_lines.append("")

    all_stats = []
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

            if tsse_vals.size > 0 and nfrag_vals.size > 0:
                all_stats.append({
                    "sample": sample_name,
                    "cells": int(tsse_vals.size),
                    "median_TSSE": float(np.median(tsse_vals)),
                    "median_nFrag": float(np.median(nfrag_vals))
                })
        except Exception as e:
            summary_lines.append(f"[WARN] Skipped {sample_name} ({e})")

    if all_stats:
        df_stats = pd.DataFrame(all_stats).set_index("sample").sort_index()
        summary_lines.append("Per-sample QC (cells / median TSSE / median nFrag):")
        summary_lines.append(df_stats.to_string())
        summary_lines.append("")
        summary_lines.append(f"TOTAL cells          : {df_stats['cells'].sum()}")
        summary_lines.append(f"Global median TSSE   : {global_median_tsse:.2f}")
        summary_lines.append(f"Global median nFrag  : {df_stats['median_nFrag'].median():.0f}")
    else:
        summary_lines.append("No valid samples for summary.")

    out_summary_txt = os.path.join(out_dir, f"qc_summary_{timestamp()}.txt")
    with open(out_summary_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))
    print(f"    Saved summary TXT: {out_summary_txt}")

    # (5) Optional: percentile-based cutoff suggestions
    if WRITE_SUGGESTIONS:
        print("[*] Writing percentile-based cutoff suggestions (CSV) ...", flush=True)

        def q(arr, p):
            arr = np.asarray(arr, dtype=float)
            arr = arr[np.isfinite(arr)]
            return float(np.percentile(arr, p)) if arr.size > 0 else float("nan")

        all_tsse_vec = []
        all_nfrag_vec = []
        for fragment_file in fragment_files:
            try:
                ad = snap.pp.import_fragments(
                    fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
                )
                compute_tsse(ad, annot_path)
                tv = np.asarray(ad.obs.get("tsse"), dtype=float)
                tv = tv[np.isfinite(tv) & (tv > 0)]
                nv = np.asarray(ad.obs.get("n_fragment"), dtype=float)
                nv = nv[np.isfinite(nv) & (nv > 0)]
                if tv.size: all_tsse_vec.append(tv)
                if nv.size: all_nfrag_vec.append(nv)
            except Exception:
                pass

        all_tsse_vec = np.concatenate(all_tsse_vec) if len(all_tsse_vec) else np.array([])
        all_nfrag_vec = np.concatenate(all_nfrag_vec) if len(all_nfrag_vec) else np.array([])

        suggest = {
            "TSSE_q10": q(all_tsse_vec, 10),
            "TSSE_q50": q(all_tsse_vec, 50),
            "TSSE_q90": q(all_tsse_vec, 90),
            "UMI_q05":  q(all_nfrag_vec, 5),
            "UMI_q95":  q(all_nfrag_vec, 95),
            "UMI_q99":  q(all_nfrag_vec, 99),
        }

        heuristic = {
            "TSSE_LINE_used":        float(global_median_tsse),
            "UMI_MIN_suggest":       int(max(3000,  suggest["UMI_q05"]  if np.isfinite(suggest["UMI_q05"])  else 5000)),
            "UMI_MAX_suggest":       int(min(100000, suggest["UMI_q99"]  if np.isfinite(suggest["UMI_q99"])  else 50000)),
        }

        import pandas as pd
        df_suggest = pd.DataFrame([{**suggest, **heuristic}])
        out_suggest = os.path.join(out_dir, f"qc_threshold_{timestamp()}.csv")
        df_suggest.to_csv(out_suggest, index=False)
        print(f"    Saved suggestions: {out_suggest}")

    print_progress(100.0, "Done")
    print("Done. All figures and summaries saved under:", out_dir, flush=True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[Interrupted] Exiting.", file=sys.stderr)
        sys.exit(130)
