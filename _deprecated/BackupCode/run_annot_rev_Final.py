#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime
from typing import Optional, List

__version__ = "2.0"

# Helpers
def _as_bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {v}")

def _int_or_list(v):
    s = str(v).strip()
    if s.lower() == "none":
        return None
    if "," in s:
        return [int(x) for x in s.split(",") if x.strip() != ""]
    return int(s)

def _ensure_dirs(*paths):
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)

def _tee_setup(log_path: Optional[str]):
    if not log_path:
        return

    class Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, data):
            for f in self.files:
                f.write(data)
        def flush(self):
            for f in self.files:
                f.flush()

    fh = open(log_path, "a", encoding="utf-8")
    sys.stdout = Tee(sys.stdout, fh)
    sys.stderr = Tee(sys.stderr, fh)

def _import_logic(project_root: Path):
    bin_dir = project_root / "bin"
    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(bin_dir))
    try:
        import bin.clusterAnnot as logic
        return logic
    except Exception:
        try:
            import clusterAnnot as logic  # fallback
            return logic
        except Exception as e2:
            raise ImportError(
                "Import failed: expected bin/clusterAnnot.py (or clusterAnnot.py) in the project root."
            ) from e2

def _find_default_h5ads(base_dir: Path) -> Path:
    filt = base_dir / "Filter_results"
    candidate = filt / "merged_doublets.h5ads"
    if candidate.exists():
        return candidate
    if filt.exists():
        h5ads = sorted(filt.rglob("*.h5ads"), key=lambda p: p.stat().st_mtime, reverse=True)
        if h5ads:
            return h5ads[0]
    raise FileNotFoundError(f"No .h5ads found under: {filt}")

def _suffix_name(name: str, suffix: Optional[str]) -> str:
    if not suffix:
        return name
    stem, ext = os.path.splitext(name)
    return f"{stem}_{suffix}{ext}"

def _apply_overrides_module(logic, args):
    mapping = {
        # IO roots
        "NJOB": args.njob,
        "BASE_DIR": str(args.base_dir),
        "H5ADS_INPUT": str(args.h5ads_input),
        "AGGR_DIR": str(args.aggr_dir),
        "FIG_DIR": str(args.fig_dir),

        # Outputs
        "H5ADS_OUT": str(args.h5ads_out),
        "ANNOT_H5AD": str(args.annot_h5ad),
        "ANNOT_GENE_H5AD": str(args.annot_gene_h5ad),
        "SUMMARY_TXT": str(args.summary_txt),
        "CSV_OUT": str(args.csv_out),

        # References
        "CHROM_SIZES_PATH": str(args.chrom_sizes),
        "GTF_PATH": str(args.gtf),
        "REF_CSV_PATH": str(args.ref_csv),

        # Algo
        "REF_LABEL_KEY": args.ref_label_key,
        "N_FEATURES": args.n_features,
        "N_NEIGHBORS": args.n_neighbors,
        "USE_DIMS": args.use_dims,
        "MAKE_GENE_ACTIVITY": args.make_gene_activity,
        "SAVE_UMAP": args.save_umap,

        # Thresholds
        "MIN_SCORE": args.min_score,
        "MIN_MARGIN": args.min_margin,
        "USE_MARGIN": args.use_margin,
    }
    for k, v in mapping.items():
        if hasattr(logic, k):
            setattr(logic, k, v)

# ---------- Marker CSV auto-detection ----------
def _auto_pick_marker_csv(ref_base: Path) -> str:
    # collect candidates
    all_csv = [p for p in ref_base.glob("*.csv") if p.is_file() and not p.name.startswith(".")]
    if not all_csv:
        raise FileNotFoundError(f"No CSV found under: {ref_base}")

    # prefer names that contain 'marker'
    marker_like = [p for p in all_csv if "marker" in p.name.lower()]
    pool: List[Path] = marker_like or all_csv

    if len(pool) == 1:
        return str(pool[0])

    def _rank(p: Path):
        name = p.name.lower()
        return (
            0 if "sc" in name else 1,
            0 if "marker" in name else 1,
            -p.stat().st_mtime,
        )

    pool.sort(key=_rank)
    return str(pool[0])

def _maybe_set_refs_from_species(args, project_root: Path):
    if not args.species:
        return

    local_base = project_root / "refGenome" / args.species
    home_base = Path(os.path.expanduser(f"~/singlecell-analysis-pipeline/refGenome/{args.species}"))
    ref_base = local_base if local_base.exists() else home_base

    # Only fill when user didn't provide manual overrides
    if args.chrom_sizes is None:
        args.chrom_sizes = str(ref_base / "genome.chrom.sizes")
    if args.gtf is None:
        args.gtf = str(ref_base / "genes.gtf")
    if args.ref_csv is None:
        args.ref_csv = _auto_pick_marker_csv(ref_base)
        print(f"[info] Auto-detected marker CSV: {args.ref_csv}")

    # Existence checks for clearer errors
    missing = [p for p in [args.chrom_sizes, args.gtf, args.ref_csv] if not Path(p).exists()]
    if missing:
        m = "\n  - ".join(missing)
        raise FileNotFoundError(
            f"[refs missing] One or more reference files do not exist for species '{args.species}':\n  - {m}"
        )

# CLI
def build_parser():
    p = argparse.ArgumentParser(
        description="Wrapper for bin/clusterAnnot.py. Project root is the folder containing this script."
    )
    # Input root (relative to project root; required)
    p.add_argument(
        "-i", "--input",
        dest="input_root",
        required=True,
        help="Root folder containing data (relative to project root). Example: 'Test'"
    )

    # Species shortcut for references
    p.add_argument(
        "-s", "--species",
        default=None,
        help="Species key under ./refGenome/ (e.g., rice). Falls back to ~/singlecell-analysis-pipeline if not found."
    )

    # Optional suffix/prefix/out-root
    p.add_argument(
        "--prefix",
        default=None,
        help="Optional suffix appended to output filenames (e.g., v2 or 20251010)."
    )
    p.add_argument(
        "--out-root",
        default=None,
        help="Optional override for output root (default: project root)."
    )

    # System & logging
    p.add_argument("--njob", type=int, default=8)
    p.add_argument("--seed", type=int, default=777)
    p.add_argument("--log-file", default=None)

    # Algorithmic params
    p.add_argument("--ref-label-key", default="clusterName")
    p.add_argument("--n-features", type=int, default=250_000)
    p.add_argument("--n-neighbors", type=int, default=10)
    p.add_argument("--use-dims", type=_int_or_list, default=None)  # "None" | "30" | "1,2,3"
    p.add_argument("--make-gene-activity", type=_as_bool, default=True)
    p.add_argument("--save-umap", type=_as_bool, default=True)

    # Thresholds
    p.add_argument("--min-score", type=float, default=0.50)
    p.add_argument("--min-margin", type=float, default=0.25)
    p.add_argument("--use-margin", type=_as_bool, default=False)

    # References (manual overrides; if omitted and --species is set, they will be filled from species)
    p.add_argument("--chrom-sizes", default=None, help="Override chrom sizes path (default: from -s).")
    p.add_argument("--gtf", default=None, help="Override GTF path (default: from -s).")
    p.add_argument("--ref-csv", default=None, help="Override marker CSV path (default: auto-detected from -s).")

    # (Advanced) Manual IO overrides (optional; normally auto)
    p.add_argument("--h5ads-input", default=None, help="Manual override for input .h5ads path.")
    p.add_argument("--aggr-dir", default=None, help="Manual override for output Annot_results dir.")
    p.add_argument("--fig-dir", default=None, help="Manual override for output Plots dir.")
    p.add_argument("--h5ads-out", default=None)
    p.add_argument("--annot-h5ad", default=None)
    p.add_argument("--annot-gene-h5ad", default=None)
    p.add_argument("--summary-txt", default=None)
    p.add_argument("--csv-out", default=None)
    return p

# Main
def main():
    args = build_parser().parse_args()

    PROJECT_ROOT = Path(__file__).resolve().parent

    base = (PROJECT_ROOT / args.input_root).resolve()
    out_root = (PROJECT_ROOT / args.out_root).resolve() if args.out_root else PROJECT_ROOT

    _maybe_set_refs_from_species(args, PROJECT_ROOT)

    if args.h5ads_input is None:
        args.h5ads_input = str(_find_default_h5ads(base))

    if args.aggr_dir is None:
        args.aggr_dir = str(out_root / "Annot_results")
    if args.fig_dir is None:
        args.fig_dir = str(out_root / "Annot_results" / "Plots")

    suffix = args.prefix
    if suffix is None:
        # Keep stable names by default; uncomment below for timestamp suffix
        # suffix = datetime.now().strftime("%Y%m%d_%H%M%S")
        suffix = None

    def _p(name):
        return str(out_root / "Annot_results" / name)

    if args.h5ads_out is None:
        args.h5ads_out = _p(_suffix_name("annot_merged.h5ads", suffix))
    if args.annot_h5ad is None:
        args.annot_h5ad = _p(_suffix_name("annot_merged_cells.h5ad", suffix))
    if args.annot_gene_h5ad is None:
        args.annot_gene_h5ad = _p(_suffix_name("annot_gene_activity.h5ad", suffix))
    if args.summary_txt is None:
        args.summary_txt = _p(_suffix_name("summary.txt", suffix))
    if args.csv_out is None:
        args.csv_out = _p(_suffix_name("annotation_results.csv", suffix))

    args.base_dir = str(base)

    _ensure_dirs(args.aggr_dir, args.fig_dir)
    if args.log_file:
        _ensure_dirs(Path(args.log_file).parent)

    _tee_setup(args.log_file)

    logic = _import_logic(PROJECT_ROOT)
    _apply_overrides_module(logic, args)

    if hasattr(logic, "main") and callable(getattr(logic, "main")):
        return logic.main(args)   # recommended
    elif hasattr(logic, "run") and callable(getattr(logic, "run")):
        return logic.run()        # legacy
    else:
        raise RuntimeError("No entrypoint found in bin/clusterAnnot.py (define main(args) or run()).")

if __name__ == "__main__":
    sys.exit(main())
