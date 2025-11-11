
#!/usr/bin/env python3

import os, sys, math, gzip, tempfile, shutil, uuid, json
from datetime import datetime
from typing import Dict, List, Tuple, Optional

import numpy as np
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.figure import Figure

from joblib import Parallel, delayed

try:
    import snapatac2 as snap
    from anndata import AnnData
except Exception as e:
    print("[ERROR] snapatac2 is required:", e, file=sys.stderr)
    raise

warnings.filterwarnings("ignore", category=FutureWarning, module="snapatac2")

HEADER = """
 =============================================================================
     Integrated scATAC-seq: Quality Control Pipeline
     Version 2.2
     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*
 ============================================================================="""

# -------- I/O defaults --------
chrom_sizes_path: Optional[str] = None
annot_path: Optional[str]        = None
base_dir: Optional[str]          = None
out_dir: Optional[str]           = None

max_recorded_size: int = 1000
subsample_lines: Optional[int]   = None

n_jobs: int = 1
grid_rows: int = 3
dpi_save: int = 300

backend_default: str = "hdf5"
chunk_size_default: int = 3000
tmp_root: Optional[str] = None

# -------- Utils --------
def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def save_fig(fig: Figure, out_base: str, dpi: int = None) -> None:
    if dpi is None:
        dpi = dpi_save
    figures_dir = os.path.join(out_dir, "Figures")
    os.makedirs(figures_dir, exist_ok=True)
    
    base_name = os.path.basename(out_base)
    png = os.path.join(figures_dir, f"{base_name}.png")
    svg = os.path.join(figures_dir, f"{base_name}.svg")

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
            if not line or line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 2: continue
            sizes[parts[0]] = int(parts[1])
    if not sizes:
        raise ValueError(f"No entries parsed from chrom sizes: {path}")
    return sizes

def collect_fragment_files(base: str) -> List[str]:
    files: List[str] = []
    for name in sorted(os.listdir(base)):
        p = os.path.join(base, name, "outs", "fragments.tsv.gz")
        if os.path.isfile(p): files.append(p)
    return files

def _available_cores() -> int:
    try:
        return max(1, os.cpu_count() or 1)
    except Exception:
        return 1

def plan_jobs(requested: int, sample_count: int) -> Tuple[int, int]:
    C = max(1, min(requested, _available_cores()))
    S = max(1, int(sample_count))
    if S >= max(1, C // 2):
        return (min(S, C), 1)  # outer dominates
    outer = min(S, C)
    inner = max(1, C // outer)
    while outer * inner > C and inner > 1:
        inner -= 1
    return (outer, max(1, inner))

def _mk_worker_tempdir(root: Optional[str]) -> str:
    base = root or tempfile.gettempdir()
    run_root = os.path.join(base, "qc_run", str(os.getpid()), str(uuid.uuid4()))
    os.makedirs(run_root, exist_ok=True)
    return tempfile.mkdtemp(dir=run_root, prefix="wkr_")

# -------- Robust FSD --------
def compute_size_hist_from_fragments(fragment_file: str, max_recorded_size: int, subsample: Optional[int] = None):
    hist = np.zeros(max_recorded_size, dtype=np.int64)
    n = 0
    opener = gzip.open if fragment_file.endswith(".gz") else open
    with opener(fragment_file, "rt", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if subsample is not None and i >= subsample:
                break
            if not line or line[0] == "#": continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            try:
                start = int(parts[1]); end = int(parts[2])
                L = end - start
                if 1 <= L <= max_recorded_size:
                    hist[L-1] += 1; n += 1
            except Exception:
                continue
    return hist if n > 0 else None

# -------- Worker functions --------
def _fsd_image_worker(fragment_file: str, chrom_sizes_dict: Dict[str, int],
                      max_recorded_size: int, dpi_save_: int,
                      inner_jobs: int, tmp_root_use: Optional[str]):
    import io
    hist = compute_size_hist_from_fragments(fragment_file, max_recorded_size, subsample=subsample_lines)
    sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
    fig, ax = plt.subplots(figsize=(5, 4))
    if hist is None or int(np.sum(hist)) == 0:
        ax.text(0.5, 0.5, "(no fragments)", ha="center", va="center"); ax.axis("off")
    else:
        xs = np.arange(1, len(hist)+1)
        ax.plot(xs, hist, linewidth=1.2)
        ax.set_xlim(0, max_recorded_size)
        ax.set_xlabel("Fragment length (bp)")
        ax.set_ylabel("Count")
        ax.set_title(sample_name)
    fig.tight_layout()
    bio = io.BytesIO()
    fig.savefig(bio, format="png", dpi=dpi_save_)
    plt.close(fig); bio.seek(0)
    return (sample_name, bio.read(), hist.tolist() if hist is not None else None, None)

def _tsse_worker(fragment_file: str, chrom_sizes_dict: Dict[str, int],
                 annot_path: str, inner_jobs: int, tmp_root_use: Optional[str]) -> Tuple[str, List[float], List[float]]:
    tmpdir = _mk_worker_tempdir(tmp_root_use)
    try:
        ad = snap.pp.import_fragments(
            fragment_file, chrom_sizes=chrom_sizes_dict, sorted_by_barcode=False,
            backend=backend_default, chunk_size=chunk_size_default,
            tempdir=tmpdir, n_jobs=max(1, int(inner_jobs))
        )
        try:
            snap.pp.compute_tsse(ad, annot_path)
        except Exception:
            import snapatac2.metrics as sm
            sm.tsse(ad, annot_path)

        ts = np.asarray(ad.obs.get("tsse"), dtype=float)
        nf = np.asarray(ad.obs.get("n_fragment"), dtype=float)
        m = np.isfinite(ts) & np.isfinite(nf) & (ts > 0) & (nf > 0)
        ts_pair = ts[m]
        nf_pair = nf[m]
        return (os.path.basename(os.path.dirname(os.path.dirname(fragment_file))), ts_pair.tolist(), nf_pair.tolist())
    except Exception:
        return (os.path.basename(os.path.dirname(os.path.dirname(fragment_file))), [], [])
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

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

# -------- Main --------
def _apply_wrapper_injections_and_defaults() -> None:
    global out_dir, base_dir, chrom_sizes_path, annot_path
    g = globals()
    base_dir  = g.get("BASE_DIR", base_dir) or os.getcwd()
    out_dir   = g.get("OUT_DIR", out_dir) or os.path.join(os.getcwd(), "QC_results")
    chrom_sizes_path = g.get("CHROM_SIZES_PATH", chrom_sizes_path)
    annot_path       = g.get("GTF_PATH", annot_path)
    os.makedirs(out_dir, exist_ok=True)

def main() -> None:
    _apply_wrapper_injections_and_defaults()
    import matplotlib as mpl
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["figure.dpi"]   = dpi_save

    if annot_path and annot_path.endswith(".gtf.gz"):
        raise ValueError("Please gunzip your GTF: use uncompressed .gtf for annot_path.")

    print("\n====== Step 1. Collecting fragment files ======", flush=True)
    fragment_files = collect_fragment_files(base_dir)
    print(f"    Found {len(fragment_files)} samples.", flush=True)

    if not chrom_sizes_path or not os.path.exists(chrom_sizes_path):
        raise FileNotFoundError("chrom.sizes path not set or not found.")
    print("\n====== Step 2. Loading chrom.sizes ======", flush=True)
    chrom_sizes = load_chrom_sizes(chrom_sizes_path)
    print(f"    Loaded {len(chrom_sizes)} contigs.", flush=True)

    # (3) FSD grid
    print("\n====== Step 3. Computing fragment-size distributions (FSD) ======", flush=True)
    n = len(fragment_files)
    ncols = 3 if n >= 3 else max(1, n)
    nrows = math.ceil(n / ncols) if n > 0 else 1
    if grid_rows and grid_rows > 0:
        nrows = min(grid_rows, n if n>0 else 1)
        ncols = math.ceil(n / max(1,nrows))
    fig_fsd, axes_fsd = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes_list = np.atleast_1d(axes_fsd).ravel().tolist()

    outer_jobs, inner_jobs = plan_jobs(int(n_jobs), len(fragment_files))
    tmp_root_use = tmp_root

    results_fsd = Parallel(n_jobs=max(1, int(outer_jobs)))(
        delayed(_fsd_image_worker)(f, chrom_sizes, max_recorded_size, dpi_save, inner_jobs, tmp_root_use)
        for f in fragment_files
    )

    for i, (sample_name, png_bytes, hist, err) in enumerate(results_fsd):
        ax = axes_list[i] if i < len(axes_list) else axes_list[-1]
        if png_bytes is not None:
            from io import BytesIO
            img = mpimg.imread(BytesIO(png_bytes), format="png")
            ax.imshow(img); ax.axis("off"); ax.set_title(sample_name)
        else:
            ax.text(0.5,0.5, err or "(error)", ha="center", va="center"); ax.axis("off")

    for j in range(len(results_fsd), len(axes_list)):
        fig_fsd.delaxes(axes_list[j])
    save_fig(fig_fsd, os.path.join(out_dir, "fragment_size"))

    # (4) TSSE violin
    print("\n====== Step 4. Computing TSSE per sample and depicting plots ======", flush=True)
    results = Parallel(n_jobs=max(1, int(outer_jobs)))(
        delayed(_tsse_worker)(f, chrom_sizes, annot_path, inner_jobs, tmp_root_use)
        for f in fragment_files
    )

    all_tsse_scores_raw: List[List[float]] = []
    sample_labels: List[str] = []
    for name, vec, _nf in results:
        if len(vec) == 0:
            continue
        all_tsse_scores_raw.append(vec)
        sample_labels.append(name)

    if len(all_tsse_scores_raw) > 0:
        all_vals = np.concatenate([np.asarray(v, dtype=float) for v in all_tsse_scores_raw])
        all_vals = all_vals[np.isfinite(all_vals) & (all_vals>0)]
        global_median_tsse = float(np.median(all_vals)) if all_vals.size>0 else float("nan")
        fig, ax = plt.subplots(figsize=(max(6, 1.5*len(sample_labels)), 4))
        ax.violinplot(all_tsse_scores_raw, showmeans=False, showmedians=True, showextrema=False)
        ax.set_xticks(np.arange(1, len(sample_labels)+1))
        ax.set_xticklabels(sample_labels, rotation=45, ha="right")
        ax.set_ylabel("TSSE")
        ax.set_title("TSSE (per-sample) Violin")
        if np.isfinite(global_median_tsse):
            ax.axhline(global_median_tsse, linestyle="--", color="red", linewidth=1)
        save_fig(fig, os.path.join(out_dir, "TSSE_violin"))
    # (4.5) Combined per-sample TSSE grid (classic style)
    try:
        print("\n====== Step 5. Rendering combined per-sample TSSE grid ======", flush=True)
        n = len(fragment_files)
        ncols = math.ceil(math.sqrt(n)) if n > 0 else 1
        nrows = math.ceil(n / ncols) if n > 0 else 1
        if grid_rows and grid_rows > 0:
            nrows = min(grid_rows, n if n>0 else 1)
            ncols = math.ceil(n / nrows) if n>0 else 1
        fig_grid, axes_grid = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows), squeeze=False)
        for i, fragment_file in enumerate(fragment_files):
            r, c = divmod(i, ncols)
            ax = axes_grid[r, c]
            sample_name = os.path.basename(os.path.dirname(os.path.dirname(fragment_file)))
            try:
                adata_grid = snap.pp.import_fragments(
                    fragment_file, chrom_sizes=chrom_sizes, sorted_by_barcode=False
                )
                # robust TSSE computation
                try:
                    compute_tsse(adata_grid, annot_path)
                except Exception:
                    try:
                        snap.pp.compute_tsse(adata_grid, annot_path)
                    except Exception:
                        import snapatac2.metrics as sm
                        sm.tsse(adata_grid, annot_path)
                img = _tsse_plot_to_image(adata_grid, cutoff_line=float(global_median_tsse) if np.isfinite(global_median_tsse) else float("nan"))
                if img is not None:
                    ax.imshow(img); ax.axis("off"); ax.set_title(sample_name)
                else:
                    ts = np.asarray(adata_grid.obs.get("tsse"), dtype=float)
                    nf = np.asarray(adata_grid.obs.get("n_fragment"), dtype=float)
                    m = np.isfinite(ts) & np.isfinite(nf) & (ts > 0) & (nf > 0)
                    if m.sum() > 0:
                        ax.scatter(nf[m], ts[m], s=1, alpha=0.3)
                        ax.set_xscale("log")
                        ax.set_xlabel("n_fragment (log)"); ax.set_ylabel("TSSE")
                        ax.set_title(sample_name)
                        if np.isfinite(global_median_tsse):
                            ax.axhline(global_median_tsse, linestyle="--", color="red", linewidth=1)
                    else:
                        ax.text(0.5, 0.5, "(no data)", ha="center", va="center")
                        ax.axis("off")
            except Exception:
                ax.text(0.5, 0.5, "(error)", ha="center", va="center")
                ax.axis("off")
        for i in range(n, nrows*ncols):
            r, c = divmod(i, ncols)
            axes_grid[r, c].axis("off")
        plt.tight_layout()
        save_fig(fig_grid, os.path.join(out_dir, "TSSE_grid"))
    except Exception as e:
        print(f"    [warn] TSSE_grid rendering failed: {e}")

    else:
        global_median_tsse = float("nan")
        print("    [warn] No valid TSSE vectors; TSSE violin skipped.")

    # (5) Extra plots using n_fragment from results
    nfrag_lists = []
    nfrag_labels = []
    for name, _vec, nfv in results:
        if nfv:
            nfrag_lists.append(np.asarray(nfv, dtype=float))
            nfrag_labels.append(name)

    if len(nfrag_lists) > 0:
        # Fragment violin
        fig_fv = plt.figure(figsize=(max(10, 1.8*len(nfrag_lists)), 6))
        xpos = np.arange(len(nfrag_lists))
        plt.violinplot(nfrag_lists, positions=xpos, showmeans=False, showmedians=True, showextrema=False)
        plt.xticks(ticks=xpos, labels=nfrag_labels, rotation=30, ha="right")
        plt.xlabel("Sample"); plt.ylabel("n_fragment")
        plt.yscale("log")
        plt.title("Violin Plot of Fragment Counts Across Samples")
        save_fig(fig_fv, os.path.join(out_dir, "fragment_violin"))

        # Knee plot (rank–count) log–log
        fig_knee = plt.figure(figsize=(10, 8))
        for lab, vec in zip(nfrag_labels, nfrag_lists):
            v = np.sort(np.asarray(vec, dtype=float))[::-1]
            if v.size == 0: continue
            ranks = np.arange(1, v.size+1)
            plt.plot(ranks, v, label=str(lab))
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Rank (cells)"); plt.ylabel("n_fragment")
        plt.title("Knee Plot (rank–count)")
        if len(nfrag_labels) <= 20: plt.legend(loc="best", fontsize=8, ncol=2)
        save_fig(fig_knee, os.path.join(out_dir, "knee_plot"))
    else:
        print("    [warn] No valid n_fragment values; fragment violin/knee skipped.")

    # (6) Summary TXT/CSV
    print("\n====== Writing summary ======", flush=True)
    headers = ["sample","cells_total","nFrag_min","nFrag_p5","nFrag_p10","nFrag_med","nFrag_p90","nFrag_p95","TSSE_min","TSSE_p5","TSSE_p10","TSSE_med","TSSE_p90","TSSE_p95"]
    rows: List[List[object]] = []

    for name, vec, nfv in results:
        tsse = np.asarray(vec, dtype=float)
        tsse = tsse[np.isfinite(tsse) & (tsse>0)]
        def pct(a,p): 
            a = np.asarray(a, dtype=float); a = a[np.isfinite(a)]
            return float(np.percentile(a, p)) if a.size>0 else float("nan")
        row = [
            name, int(len(tsse)),
            pct(nfv,0) if len(nfv)>0 else float("nan"), pct(nfv,5), pct(nfv,10), pct(nfv,50), pct(nfv,90), pct(nfv,95),
            np.min(tsse) if tsse.size>0 else float("nan"),
            pct(tsse,5), pct(tsse,10), pct(tsse,50), pct(tsse,90), pct(tsse,95)
        ]
        rows.append(row)

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
    if len(rows)>0:
        summary_lines.append("\t".join(headers))
        for vals in rows:
            summary_lines.append("\t".join(str(v) for v in vals))
        all_tsse_all = []
        for _, vec, _ in results:
            v = np.asarray(vec, dtype=float)
            v = v[np.isfinite(v) & (v>0)]
            if v.size>0: all_tsse_all.append(v)
        if len(all_tsse_all)>0:
            gmed = float(np.median(np.concatenate(all_tsse_all)))
            summary_lines.append(f"Global median TSSE: {gmed:.2f}")
    else:
        summary_lines.append("No valid samples for summary.")

    ensure_dir(out_dir)
    out_txt = os.path.join(out_dir, "qc_summary.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))
    print(f"    Saved summary TXT: {out_txt}")

    try:
        import pandas as pd
        if len(rows)>0:
            df = pd.DataFrame(rows, columns=headers).set_index("sample").sort_index()
            out_csv = os.path.join(out_dir, "qc_summary_table.csv")
            df.to_csv(out_csv, index=True)
            print(f"    Saved summary CSV: {out_csv}")
    except Exception as e:
        print(f"    [warn] CSV export failed: {e}")

    # (7) JSON sidecars for interactive HTML
    try:
        assets = os.path.join(out_dir, "QC_reports", "data")
        os.makedirs(assets, exist_ok=True)
        # fragment_size histograms
        fsd_payload = {"max_recorded_size": int(max_recorded_size), "samples": []}
        for name, png_bytes, hist, err in results_fsd:
            if hist is not None:
                fsd_payload["samples"].append({"name": name, "hist": hist})
        with open(os.path.join(assets, "fragment_size.json"), "w", encoding="utf-8") as fjson:
            json.dump(fsd_payload, fjson)

        # TSSE + n_fragment
        tsse_labels = []; tsse_data = []; nfrag_data = []
        for name, vec, nfv in results:
            tsse_labels.append(name); tsse_data.append(vec); nfrag_data.append(nfv)
        payload_tv = {"labels": tsse_labels, "tsse": tsse_data, "nfrag": nfrag_data,
                      "global_median_tsse": float(global_median_tsse) if np.isfinite(global_median_tsse) else None}
        with open(os.path.join(assets, "tsse_violin.json"), "w", encoding="utf-8") as fjson:
            json.dump(payload_tv, fjson)

        # TSSE grid density per sample
        grid_bins_x = 60; grid_bins_y = 60
        tg = {"samples": [], "bins_x": grid_bins_x, "bins_y": grid_bins_y}
        for name, vec, nfv in results:
            x = np.asarray(nfv, dtype=float); y = np.asarray(vec, dtype=float)
            # data are already pair-filtered in _tsse_worker; keep finite guard
            k = np.isfinite(x) & np.isfinite(y)
            x = x[k]; y = y[k]
            if x.size == 0:
                continue
            try:
                x_edges = np.geomspace(max(1, x.min()), x.max(), grid_bins_x+1)
            except Exception:
                x_edges = np.linspace(x.min(), x.max(), grid_bins_x+1)
            y_edges = np.linspace(y.min(), y.max(), grid_bins_y+1)
            H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges])
            tg["samples"].append({"name": name, "z": H.T.tolist(), "x_edges": xe.tolist(), "y_edges": ye.tolist()})
        with open(os.path.join(assets, "tsse_grid.json"), "w", encoding="utf-8") as fjson:
            json.dump(tg, fjson)
    except Exception as e:
        print(f"    [warn] JSON export failed: {e}")

    print("Done. All figures and summaries saved under:", out_dir, flush=True)

if __name__ == "__main__":
    main()
