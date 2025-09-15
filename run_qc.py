#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Optional, List

__version__ = "1.0"

def _is_readable_file(p: Path) -> bool:
    try:
        return p.exists() and p.is_file()
    except Exception:
        return False

def _first_existing(candidates: List[Path]):
    for p in candidates:
        if _is_readable_file(p):
            return p
    return None

def detect_rice_paths(project_root: Path):
    base = project_root / "refGenome" / "rice"
    if not base.exists():
        return None, None
    chrom = _first_existing([
        base / "chrom.sizes",
        base / "rice.chrom.sizes",
        *list(base.glob("*.chrom.sizes"))
    ])
    gtf = _first_existing([
        base / "genes.gtf",
        base / "rice.gtf",
        *list(base.glob("*.gtf")),
        *list(base.glob("*.gtf.gz"))
    ])
    return (str(chrom) if chrom else None, str(gtf) if gtf else None)

# CLI
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_qc",
        description="Wrapper for bin/IntegQC_V8.py with project-root refGenome detection and standard defaults.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("-i", "--input-dir", type=Path, required=True,
                   default=argparse.SUPPRESS, metavar="{INPUT_DIR}",
                   help="Input base directory; assigned to logic base_dir/BASE_DIR.")
    p.add_argument("-s", "--species", type=str, choices=["rice"], required=True,
                   default=argparse.SUPPRESS,
                   help="Species key for auto-detect. Currently supports 'rice'.")
    p.add_argument("-mrs", "--max-recorded-size", type=int, default=1000,
                   metavar="{NUMBER}",
                   help="Upper bound for fragment-size histogram.")
    p.add_argument("-sl", "--subsample-lines", type=int, default=None,
                   metavar="{NUMBER}",
                   help="Subsample number of lines for fast preview.")
    p.add_argument("-grid", "--grid-rows", type=int, default=3,
                   metavar="{NUMBER}",
                   help="Grid rows for plots.")
    p.add_argument("-dpi", "--dpi", type=int, default=300,
                   metavar="{NUMBER}",
                   help="Figure DPI.")
    return p

# main
def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(__file__).resolve().parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    try:
        logic = importlib.import_module("bin.IntegQC_V8")
    except Exception as e:
        print("[ERROR] Could not import 'bin.IntegQC_V8'. Ensure bin/ is a package and file exists.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        return 2

    # Inject base_dir and out_dir as strings
    base_dir_str = str(Path(args.input_dir).resolve())
    outdir_str = str((Path.cwd() / "QC_results").resolve())
    try:
        if hasattr(logic, "BASE_DIR"):
            setattr(logic, "BASE_DIR", base_dir_str)
        if hasattr(logic, "base_dir"):
            setattr(logic, "base_dir", base_dir_str)
        if hasattr(logic, "OUT_DIR"):
            setattr(logic, "OUT_DIR", outdir_str)
        if hasattr(logic, "out_dir"):
            setattr(logic, "out_dir", outdir_str)
    except Exception as e:
        print(f"[WARN] Failed to inject base/out dirs: {e}", file=sys.stderr)

    # Species detection (strings)
    if args.species == "rice":
        chrom_str, gtf_str = detect_rice_paths(project_root)
        if chrom_str:
            for name in ["CHROM_SIZES_PATH", "chrom_sizes_path"]:
                if hasattr(logic, name):
                    setattr(logic, name, chrom_str)
        else:
            print("[WARN] refGenome/rice chrom.sizes not found under project root.", file=sys.stderr)
        if gtf_str:
            for name in ["GTF_PATH", "annot_path"]:
                if hasattr(logic, name):
                    setattr(logic, name, gtf_str)
        else:
            print("[WARN] refGenome/rice GTF not found under project root.", file=sys.stderr)

    # Inject knobs if present
    if hasattr(logic, "max_recorded_size") and getattr(args, "max_recorded_size", None) is not None:
        logic.max_recorded_size = int(args.max_recorded_size)
    if hasattr(logic, "subsample_lines") and getattr(args, "subsample_lines", None) is not None:
        logic.subsample_lines = int(args.subsample_lines)
    if hasattr(logic, "grid_rows") and getattr(args, "grid_rows", None) is not None:
        logic.grid_rows = int(args.grid_rows)
    if hasattr(logic, "dpi_save") and getattr(args, "dpi", None) is not None:
        logic.dpi_save = int(args.dpi)

    # Prefer logic.main(); fallback to run_pipeline()
    ran = False
    if hasattr(logic, "main") and callable(getattr(logic, "main")):
        try:
            logic.main()
            ran = True
        except Exception as e:
            print(f"[WARN] logic.main() failed: {e}", file=sys.stderr)

    if not ran and hasattr(logic, "run_pipeline") and callable(getattr(logic, "run_pipeline")):
        try:
            logic.run_pipeline()
            ran = True
        except Exception as e:
            print(f"[WARN] run_pipeline invocation failed: {e}", file=sys.stderr)

    if not ran:
        print("[INFO] No callable entry found in bin/IntegQC_V8.py. "
              "Variables have been set; please invoke the appropriate function inside the module.", file=sys.stderr)

    print(f"[OK] Finished.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
