# Full pipeline (headless, robust) with Global 0→100% progress:
#  - Fragment-size distribution grid  -> PNG & SVG
#  - TSSE violin across samples       -> PNG & SVG (skip empty samples safely)
#  - Per-sample TSSE plots            -> PNG & SVG (robust save: Figure/Image/fallback)
#  - TSSE mosaic                      -> PNG & SVG (only valid images)
# Notes (EN):
# 1) Use Cell Ranger refdata: genome.fa(.fai)->chrom.sizes, genes.gtf (unzipped)
# 2) import_fragments requires a dict for chrom_sizes
# 3) TSSE wrapper: pass GTF path directly; fallback to older APIs if needed
# 4) No plt.show; progress via single-line print + optional tqdm
# 5) GlobalProgress ensures the overall progress bar spans exactly 0~100% regardless of sample count

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
from typing import Optional, List, Dict
from matplotlib.figure import Figure
from datetime import datetime
from tqdm.auto import tqdm

__version__ = 1.4

HEADER = "\n ================================================================================="
HEADER += "\n     Integrated scATAC-seq QC Pipeline"
HEADER += "\n     Version {0}".format(__version__)
HEADER += "\n     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================"

desc_txt = """
This tool provides an integrated quality-control workflow for
single-cell ATAC-seq data. It automatically processes fragment
files to generate fragment-size distribution plots, TSSE (TSS
enrichment) metrics, per-sample QC figures, and summary mosaics.
All results are saved as publication-ready PNG and SVG files,
ensuring reproducible and standardized QC across samples.
"""

# Optional: silence FutureWarnings from plotting internals
warnings.filterwarnings("ignore", category=FutureWarning, module="snapatac2.plotting")

# User variables (EDIT ME) -> genome관련 내용은 config로 빼서 human/mouse는 직접 다운받아 사용하고, custom genome은 경로를 지정하는 것으로 수정.
chrom_sizes_path = "/Users/jchoi/Desktop/Test/genome.chrom.sizes"
annot_path        = "/Users/jchoi/Desktop/Test/genes.gtf"
base_dir          = "/Users/jchoi/Desktop/Test"            # contains <sample>/outs/fragments.tsv.gz, 특정 폴더의 모든 하위 폴더를 인식.
out_dir           = "/Users/jchoi/Desktop/Test/qc_outputs"
tsse_cutoff       = 6.0
max_recorded_size = 1000
subsample_lines   = None           # e.g., 2_000_000 for fast preview (FSD fallback)
mosaic_rows       = 3
dpi_save          = 300
USE_TQDM          = False          # Set True to keep per-loop tqdm bars (may duplicate progress)

# Utilities
def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")

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
    sizes = {}
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
    files = []
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

# TSSE wrapper: pass path directly; fallback to older APIs if needed
def compute_tsse(adata: AnnData, gtf_path: str) -> None:
    """Robust TSSE call across SnapATAC2 versions."""
    try:
        # Preferred: path directly (gene_anno accepts GTF/GFF path)
        snap.metrics.tsse(adata, gtf_path)
        return
    except TypeError:
        # Older/alt builds: accept an annotation object or keyword args
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
    # If all fail, re-raise clear error
    raise RuntimeError("TSSE call failed: incompatible SnapATAC2 API for GTF input.")

# Matplotlib-only TSSE plot (no IPython/Plotly/Kaleido)
def tsse_plot_matplotlib(adata: AnnData,
                         tsse_cutoff: float = 6.0,
                         max_points: int = 50000) -> Figure:
    """Scatter: TSSE vs n_fragment + cutoff line."""
    x = np.asarray(adata.obs.get("n_fragment"), dtype=float)
    y = np.asarray(adata.obs.get("tsse"), dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[m]; y = y[m]
    fig = plt.figure(figsize=(6, 5))
    ax = plt.gca()

    if x.size == 0:
        ax.text(0.5, 0.5, "No valid TSSE values", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        return fig

    if x.size > max_points:
        idx = np.random.default_rng(0).choice(x.size, size=max_points, replace=False)
        x = x[idx]; y = y[idx]

    ax.scatter(x, y, s=6, alpha=0.25)
    ax.axhline(tsse_cutoff, linestyle="--", color="red", label=f"TSSE = {tsse_cutoff}")
    ax.set_xlabel("n_fragment")
    ax.set_ylabel("TSSE")
    ax.set_title("TSSE QC (matplotlib)")
    ax.legend()
    fig.tight_layout()
    return fig

# Global progress (0~100%)
def print_progress(pct: float, msg: str = ""):
    """Progress: print on its own line instead of overwriting."""
    pct = max(0, min(100, pct))
    print(f"Progress: {pct:6.2f}%  {msg[:80]:<80}")

class GlobalProgress:
    """
    Tracks global progress as a sum of work units.
    Each unit contributes equally to final 100%.
    """
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

print(HEADER)
print(desc_txt)

GP = GlobalProgress()
print_progress(0.0, "Initializing")

# (0) Pre-flight checks
if annot_path.endswith(".gtf.gz"):
    raise ValueError("Please gunzip: use an uncompressed .gtf for 'annot_path'.")

ensure_dir(out_dir)
fig_dir = os.path.join(out_dir, "figures"); ensure_dir(fig_dir)
per_sample_dir = os.path.join(fig_dir, "per_sample_tsse"); ensure_dir(per_sample_dir)

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
# Weights: tune as you like (only relative sizes matter)
GP.add_units(n)          # FSD loop
GP.add_units(2)          # FSD figure save + overhead
GP.add_units(n)          # TSSE per-sample loop
GP.add_units(2)          # Violin build + save
GP.add_units(n)          # Per-sample TSSE plots
GP.add_units(2)          # Mosaic build + save
GP.add_units(2)          # Pre/Post housekeeping

# (1) Fragment-size distribution grid (save only)
print("[3/5] Computing fragment-size distributions (grid) ...", flush=True)
ncols = math.ceil(math.sqrt(n)) if n > 0 else 1
nrows = math.ceil(n / ncols) if n > 0 else 1

fig_fsd, axes_fsd = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
axes_list = np.atleast_1d(axes_fsd).ravel().tolist()

iterable_fsd = fragment_files if not USE_TQDM else tqdm(fragment_files, desc="  FSD per sample")
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

out_base = os.path.join(fig_dir, f"fragment_size_grid_{timestamp()}")
save_fig(fig_fsd, out_base)
GP.inc(2, "Saved FSD grid")

# (2) TSSE per sample + violin (save only, robust filtering)
print("[4/5] Computing TSSE per sample + violin plot ...", flush=True)

all_tsse_scores_raw: List[List[float]] = []
sample_labels: List[str] = []
kept_counts = 0
skipped = []

iterable_tsse = fragment_files if not USE_TQDM else tqdm(fragment_files, desc="  TSSE per sample")
for i, fragment_file in enumerate(iterable_tsse):
    sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
    try:
        adata = snap.pp.import_fragments(
            fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False
        )
        compute_tsse(adata, annot_path)
        tsse_vals = np.asarray(adata.obs["tsse"], dtype=float)
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

if kept_counts > 0:
    x_values = np.arange(kept_counts)
    fig_violin = plt.figure(figsize=(12, 8))
    plt.violinplot(all_tsse_scores_raw, positions=x_values, showmeans=False, showmedians=True)
    plt.axhline(y=tsse_cutoff, linestyle="--", label=f"TSS Enrichment = {tsse_cutoff}", color="red")
    plt.xlabel("Sample")
    plt.ylabel("TSS Enrichment Score")
    plt.title("Violin Plot of TSSE Across Samples")
    plt.xticks(ticks=x_values, labels=sample_labels, rotation=30, ha="right")
    plt.legend()
    plt.tight_layout()

    out_base = os.path.join(fig_dir, f"tsse_violin_{timestamp()}")
    save_fig(fig_violin, out_base)
    GP.inc(2, "Saved TSSE violin")
else:
    print("    [warn] No sample had valid TSSE values. Violin plot skipped.")
    GP.inc(2, "Violin skipped (no valid TSSE)")

if skipped:
    print(f"    Skipped samples (no/invalid TSSE): {', '.join(skipped)}")

# (3) Per-sample TSSE plots (robust save: Figure/Image/fallback)
print("[5/5] Rendering per-sample TSSE plots (PNG & SVG) ...", flush=True)
saved_imgs = []

iterable_ps = fragment_files if not USE_TQDM else tqdm(fragment_files, desc="  TSSE plot per sample")
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

        ts = np.asarray(adata.obs["tsse"], dtype=float)
        ts = ts[np.isfinite(ts)]; ts = ts[ts > 0]
        if ts.size == 0:
            saved_imgs.append(None)
            GP.inc(1, f"Per-sample TSSE {i+1}/{n}")
            continue

        base = os.path.join(per_sample_dir, f"{sample_name}_tsse_{timestamp()}")
        plt.close('all')

        # Try native plotter
        ret = None
        try:
            ret = snap.pl.tsse(adata, interactive=False)
        except Exception:
            ret = None

        # Case A: matplotlib Figure
        if isinstance(ret, Figure):
            fig_tsse = ret
            save_fig(fig_tsse, base)
            saved_imgs.append(base + ".png")
            GP.inc(1, f"Per-sample TSSE {i+1}/{n}")
            continue

        # Case B: IPython.display.Image (PNG bytes)
        try:
            from IPython.display import Image as IPyImage
            if isinstance(ret, IPyImage) and getattr(ret, "data", None):
                with open(base + ".png", "wb") as f:
                    f.write(ret.data)
                print(f"    Saved: {base}.png")
                # SVG via matplotlib fallback
                fig_svg = tsse_plot_matplotlib(adata, tsse_cutoff=tsse_cutoff)
                fig_svg.savefig(base + ".svg", dpi=dpi_save, bbox_inches="tight", facecolor="white")
                plt.close(fig_svg)
                print(f"    Saved: {base}.svg")
                saved_imgs.append(base + ".png")
                GP.inc(1, f"Per-sample TSSE {i+1}/{n}")
                continue
        except Exception:
            pass

        # Case C: pure-matplotlib fallback
        fig_fb = tsse_plot_matplotlib(adata, tsse_cutoff=tsse_cutoff)
        save_fig(fig_fb, base)
        saved_imgs.append(base + ".png")
    except Exception:
        saved_imgs.append(None)
    GP.inc(1, f"Per-sample TSSE {i+1}/{n}")

# Mosaic (only if any)
valid_imgs = [p for p in saved_imgs if p and os.path.exists(p)]
if valid_imgs:
    num = len(valid_imgs)

    # rows: limited by mosaic_rows, cols computed automatically
    rows = min(mosaic_rows, num)
    cols = math.ceil(num / rows)

    # always return 2D array -> 매우중요, 1일때 디멘션 desc, 없으면 1일때 r,c에서 인덱싱 에러발생했음.
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
    out_base = os.path.join(fig_dir, f"tsse_mosaic_{timestamp()}")
    save_fig(fig_mosaic, out_base)
    GP.inc(2, "Saved TSSE mosaic")
else:
    print("    No valid per-sample TSSE images. Mosaic skipped.")
    GP.inc(2, "Mosaic skipped")

GP.inc(2, "Finalizing")
print_progress(100.0, "Done")
print("Done. All figures saved under:", fig_dir, flush=True)
