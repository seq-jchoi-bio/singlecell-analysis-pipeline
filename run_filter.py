#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Optional, Tuple

__version__ = "3.0"

HEADER = "\n ================================================================================="
HEADER += "\n      Integrated scATAC-seq: Filter + Doublet Removal"
HEADER += f"\n      Version {__version__}"
HEADER += "\n      (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*"
HEADER += "\n ================================================================================="

desc_txt = """This script:
  (1) computes TSSE from fragments + GTF,
  (2) builds a minimal QC table (is__cell_barcode, passed_filters=n_fragment, TSSE),
  (3) filters cells by TSSE >= cutoff and UMI within [min, max],
  (4) applies Scrublet-based doublet removal on filtered cells,
and saves one consolidated summary under Filter_results/summary.txt.
"""

# ----- Strict rule of genome detection -------

def _strict_species_paths(project_root: Path, species: str) -> Tuple[Optional[str], Optional[str]]:
    """Look only for refGenome/<species>/genes.gtf and genome.chrom.sizes (exact names)."""
    base = project_root / "refGenome" / species
    if not base.exists() or not base.is_dir():
        return None, None
    gtf = base / "genes.gtf"
    chrom = base / "genome.chrom.sizes"
    if gtf.exists() and gtf.is_file() and chrom.exists() and chrom.is_file():
        return str(chrom), str(gtf)
    return None, None

def _strict_auto_detect(project_root: Path) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Auto-detect a species folder under refGenome/ that contains BOTH:
      - genes.gtf
      - genome.chrom.sizes
    Priority: 'rice' if valid → otherwise first valid folder in alphabetical order.
    """
    ref_root = project_root / "refGenome"
    if not ref_root.exists() or not ref_root.is_dir():
        return None, None, None

    # 1) Prefer 'rice' if valid
    c, g = _strict_species_paths(project_root, "rice")
    if c and g:
        return c, g, "rice"

    # 2) Else pick first valid folder by name
    for d in sorted([p for p in ref_root.iterdir() if p.is_dir()], key=lambda p: p.name.lower()):
        c2, g2 = _strict_species_paths(project_root, d.name)
        if c2 and g2:
            return c2, g2, d.name

    return None, None, None

# ---------- CLI ----------

def str2bool(v):
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ("yes", "true", "t", "1"):
        return True
    if v in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("Expected True or False")

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_filter",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument(
        "-i", "--input-dir", type=Path, required=True, metavar="{INPUT_DIR}",
        help="Input base directory; assigned to logic base_dir."
    )
    p.add_argument(
        "-s", "--species", type=str, required=False, default=None, metavar="{species}",
        help="Species folder name under refGenome/. Must contain genes.gtf and genome.chrom.sizes."
    )
    # Simplified thresholds/knobs
    p.add_argument("-tss", type=float, metavar="{FLOAT}", default=None,
                   help="TSSE cutoff (default inside logic if unset).")
    p.add_argument("-min", type=int, metavar="{INT}", default=None,
                   help="Minimum UMI (default inside logic if unset).")
    p.add_argument("-max", type=int, metavar="{INT}", default=None,
                   help="Maximum UMI (default inside logic if unset).")
    p.add_argument("-feat", type=int, metavar="{INT}", default=None,
                   help="Number of features for scrublet/kNN graph (default inside logic if unset).")
    p.add_argument("--strip", type=str2bool, default=False, metavar="{False/True}",
                   help="Strip trailing '-<digits>' from barcodes.")
    p.add_argument("-j", "--n-jobs", type=int, default=None, metavar="{int}",
                   help="Parallel workers for per-sample filtering (like QC).")
    return p

# ---------- Main ----------

def main(argv: Optional[list] = None) -> int:
    print(HEADER)
    print(desc_txt)
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(__file__).resolve().parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    try:
        logic = importlib.import_module("bin.FilterDoublet")
    except Exception as e:
        print("[ERROR] Could not import 'bin.FilterDoublet'. Ensure bin/ is a package and file exists.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        return 2

    # Inject base_dir and out_dir (fixed to CWD/Filter_results) as strings
    base_dir_str = str(Path(args.input_dir).resolve())
    outdir_str = str((Path.cwd() / "Filter_results").resolve())
    try:
        if hasattr(logic, "base_dir"):
            setattr(logic, "base_dir", base_dir_str)
        if hasattr(logic, "out_dir"):
            setattr(logic, "out_dir", outdir_str)
    except Exception as e:
        print(f"[WARN] Failed to inject base/out dirs: {e}", file=sys.stderr)

    # -------- Strict genome loading --------
    chrom_str: Optional[str] = None
    gtf_str: Optional[str] = None
    chosen_species: Optional[str] = None

    if args.species:
        chrom_str, gtf_str = _strict_species_paths(project_root, args.species)
        chosen_species = args.species if (chrom_str and gtf_str) else None
        if not (chrom_str and gtf_str):
            print(f"[ERROR] Missing genome files under refGenome/{args.species}: require genes.gtf and genome.chrom.sizes", file=sys.stderr)
            return 3
    else:
        c3, g3, sp3 = _strict_auto_detect(project_root)
        if c3 and g3 and sp3:
            chrom_str, gtf_str, chosen_species = c3, g3, sp3
            print(f"[INFO] Auto-detected species: {sp3}")
        else:
            print("[ERROR] Auto-detect failed in refGenome/: no folder containing both genes.gtf and genome.chrom.sizes", file=sys.stderr)
            return 3

    # Inject genome paths into logic
    try:
        if hasattr(logic, "chrom_sizes_path"):
            setattr(logic, "chrom_sizes_path", chrom_str)
        if hasattr(logic, "annot_path"):
            setattr(logic, "annot_path", gtf_str)
        # Also set alternative variable names if present
        if hasattr(logic, "CHROM_SIZES_PATH"):
            setattr(logic, "CHROM_SIZES_PATH", chrom_str)
        if hasattr(logic, "GTF_PATH"):
            setattr(logic, "GTF_PATH", gtf_str)
    except Exception as e:
        print(f"[WARN] Failed to inject genome paths: {e}", file=sys.stderr)

    # Optional threshold overrides
    if hasattr(logic, "TSS_CUTOFF") and args.tss is not None:
        logic.TSS_CUTOFF = float(args.tss)
    if hasattr(logic, "UMI_MIN") and args.min is not None:
        logic.UMI_MIN = int(args.min)
    if hasattr(logic, "UMI_MAX") and args.max is not None:
        logic.UMI_MAX = int(args.max)
    if hasattr(logic, "N_FEATURES") and args.feat is not None:
        logic.N_FEATURES = int(args.feat)
    if hasattr(logic, "STRIP_BARCODE_SUFFIX"):
        logic.STRIP_BARCODE_SUFFIX = bool(getattr(args, "strip", False))
    if getattr(args, 'n_jobs', None) is not None:
        try:
            setattr(logic, 'n_jobs', int(args.n_jobs))
        except Exception:
            pass

    # Run logic main()
    ran = False
    if hasattr(logic, "main") and callable(getattr(logic, "main")):
        try:
            logic.main()
            ran = True
        except Exception as e:
            print(f"[WARN] logic.main() failed: {e}", file=sys.stderr)

    if not ran:
        print("[INFO] No callable entry found in bin/FilterDoublet. "
              "Variables have been set; please invoke the appropriate function inside the module.", file=sys.stderr)

    print(f"[OK] Finished. (species={chosen_species})")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
