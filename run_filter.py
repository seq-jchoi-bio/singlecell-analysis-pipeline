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

def _first_existing(candidates: List[Path]) -> Optional[Path]:
    for p in candidates:
        if _is_readable_file(p):
            return p
    return None

def detect_rice_paths(project_root: Path) -> tuple[Optional[str], Optional[str]]:
    base = project_root / "refGenome" / "rice"
    if not base.exists():
        return None, None
    chrom = _first_existing([
        base / "chrom.sizes",
        base / "rice.chrom.sizes",
        *list(base.glob("*.chrom.sizes")),
    ])
    gtf = _first_existing([
        base / "genes.gtf",
        base / "rice.gtf",
        *list(base.glob("*.gtf")),
        *list(base.glob("*.gtf.gz")),
    ])
    return (str(chrom) if chrom else None, str(gtf) if gtf else None)

def auto_detect_genome(project_root: Path) -> tuple[Optional[str], Optional[str], Optional[str]]:
    base = project_root / "refGenome"
    if not base.exists() or not base.is_dir():
        return None, None, None
    candidates: list[tuple[str, Path, Path]] = []
    for d in base.iterdir():
        if not d.is_dir():
            continue
        chrom = _first_existing([
            d / "chrom.sizes",
            d / f"{d.name}.chrom.sizes",
            *list(d.glob("*.chrom.sizes")),
        ])
        gtf = _first_existing([
            d / "genes.gtf",
            d / f"{d.name}.gtf",
            *list(d.glob("*.gtf")),
            *list(d.glob("*.gtf.gz")),
        ])
        if chrom and gtf:
            candidates.append((d.name, chrom, gtf))
    if not candidates:
        return None, None, None
    # Prefer rice
    candidates.sort(key=lambda x: (0 if x[0].lower() == "rice" else 1, x[0].lower()))
    sp, chrom, gtf = candidates[0]
    return str(chrom), str(gtf), sp

# CLI
def str2bool(v):
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ("yes", "true", "t", "1"):
        return True
    elif v in ("no", "false", "f", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Expected True or False")
        
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_filter",
        description="Wrapper for Filter+Doublet pipeline with genome auto-detection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("-i", "--input-dir", type=Path, required=True,
                   default=argparse.SUPPRESS, metavar="{INPUT_DIR}",
                   help="Input base directory; assigned to logic base_dir.")
    p.add_argument("-s", "--species", type=str, choices=["rice"], required=True,
                   default=argparse.SUPPRESS,
                   help="Species key for auto-detect. Currently supports 'rice'.")
    p.add_argument("-cs", "--chrom-sizes", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Path to chrom.sizes (overrides species/auto-detect).")
    p.add_argument("-gtf", "--gtf", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Path to gene annotation GTF or GTF.gz (overrides species/auto-detect).")
    # Simplified thresholds/knobs
    p.add_argument("-tss", type=float, metavar="{FLOAT}", default=argparse.SUPPRESS,
                   help="TSSE cutoff (default: 1.0).")
    p.add_argument("-min", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Minimum UMI (default: 500).")
    p.add_argument("-max", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Maximum UMI (default: 100000).")
    p.add_argument("-feat", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Number of features for scrublet/kNN graph (default: 250000).")
    p.add_argument("--strip", type=str2bool, default=False, metavar="{False/True}",
                   help="Strip trailing '-<digits>' from barcodes.")
    return p

# main
def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(__file__).resolve().parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    try:
        logic = importlib.import_module("bin.FilterDoublet_V7")
    except Exception as e:
        print("[ERROR] Could not import 'bin.FilterDoublet_V7'. Ensure bin/ is a package and file exists.", file=sys.stderr)
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

    # Genome references
    chrom_str = str(args.chrom_sizes.resolve()) if args.chrom_sizes else None
    gtf_str = str(args.gtf.resolve()) if args.gtf else None

    if not chrom_str or not gtf_str:
        if args.species == "rice":
            c2, g2 = detect_rice_paths(project_root)
            chrom_str = chrom_str or c2
            gtf_str = gtf_str or g2
        else:
            c3, g3, sp3 = auto_detect_genome(project_root)
            if c3 and g3:
                chrom_str = chrom_str or c3
                gtf_str = gtf_str or g3
                if args.species is None and sp3:
                    print(f"[INFO] Auto-detected species: {sp3}")

    # Early validation
    if not chrom_str or not gtf_str:
        print("[ERROR] Genome files not found. Provide -s rice or -cs <chrom.sizes> and -gtf <genes.gtf(.gz)>.", file=sys.stderr)
        return 3

    # Inject genome paths into logic
    if hasattr(logic, "chrom_sizes_path"):
        setattr(logic, "chrom_sizes_path", chrom_str)
    if hasattr(logic, "annot_path"):
        setattr(logic, "annot_path", gtf_str)

    # Optional threshold overrides
    if hasattr(logic, "TSS_CUTOFF") and getattr(args, "tss", None) is not None:
        logic.TSS_CUTOFF = float(args.tss)
    if hasattr(logic, "UMI_MIN") and getattr(args, "min", None) is not None:
        logic.UMI_MIN = int(args.min)
    if hasattr(logic, "UMI_MAX") and getattr(args, "max", None) is not None:
        logic.UMI_MAX = int(args.max)
    if hasattr(logic, "N_FEATURES") and getattr(args, "feat", None) is not None:
        logic.N_FEATURES = int(args.feat)
    if hasattr(logic, "STRIP_BARCODE_SUFFIX"):
        logic.STRIP_BARCODE_SUFFIX = bool(getattr(args, "strip", False))

    # Run logic main()
    ran = False
    if hasattr(logic, "main") and callable(getattr(logic, "main")):
        try:
            logic.main()
            ran = True
        except Exception as e:
            print(f"[WARN] logic.main() failed: {e}", file=sys.stderr)

    if not ran:
        print("[INFO] No callable entry found in bin/FilterDoublet_V7.py. "
              "Variables have been set; please invoke the appropriate function inside the module.", file=sys.stderr)

    print(f"[OK] Finished.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
