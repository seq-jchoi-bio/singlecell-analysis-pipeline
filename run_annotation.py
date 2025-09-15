#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Optional, List, Any

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

def detect_species_genome(project_root: Path, species: Optional[str]) -> tuple[Optional[str], Optional[str], Optional[str]]:
    base = project_root / "refGenome"
    if not base.exists() or not base.is_dir():
        return None, None, None

    def find_for_dir(d: Path) -> tuple[Optional[Path], Optional[Path]]:
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
        return chrom, gtf

    if species:
        d = base / species
        if d.is_dir():
            chrom, gtf = find_for_dir(d)
            return (str(chrom) if chrom else None,
                    str(gtf) if gtf else None,
                    species)
        return None, None, species

    candidates: list[tuple[str, Path, Path]] = []
    for d in base.iterdir():
        if not d.is_dir():
            continue
        chrom, gtf = find_for_dir(d)
        if chrom and gtf:
            candidates.append((d.name, chrom, gtf))
    if not candidates:
        return None, None, None
    candidates.sort(key=lambda x: (0 if x[0].lower() == "rice" else 1, x[0].lower()))
    sp, chrom, gtf = candidates[0]
    return str(chrom), str(gtf), sp

def detect_reference_assets(project_root: Path, species: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    base = project_root / "refGenome"
    if not base.exists() or not base.is_dir():
        return None, None

    def scan_dir(d: Path):
        csvs = list(d.glob("*.csv"))
        csvs.sort(key=lambda p: (0 if "marker" in p.name.lower() else 1, p.name.lower()))
        h5ads = list(d.glob("*.h5ad")) + list(d.glob("*.h5ads"))
        csv = csvs[0] if csvs else None
        h5 = h5ads[0] if h5ads else None
        return csv, h5

    if species:
        d = base / species
        if d.is_dir():
            csv, h5 = scan_dir(d)
            if csv or h5:
                return (str(csv) if csv else None, str(h5) if h5 else None)

    for d in base.iterdir():
        if not d.is_dir():
            continue
        csv, h5 = scan_dir(d)
        if csv or h5:
            return (str(csv) if csv else None, str(h5) if h5 else None)
    return None, None

def parse_bool(val: Optional[str]) -> Optional[bool]:
    if val is None:
        return None
    s = str(val).strip().lower()
    if s in {"1","true","t","yes","y","on"}:
        return True
    if s in {"0","false","f","no","n","off"}:
        return False
    raise argparse.ArgumentTypeError("Expected TRUE/FALSE (or 1/0, yes/no).")

def parse_dims_arg(val: Optional[str]) -> Optional[Any]:
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    if "," in s:
        out: list[int] = []
        for tok in s.split(","):
            tok = tok.strip()
            if tok:
                out.append(int(tok))
        return out if out else None
    try:
        return int(s)
    except ValueError:
        raise argparse.ArgumentTypeError("USE_DIMS expects int or comma-separated ints (e.g., 30 or 1,2,3).")

# CLI
def build_parser(logic_defaults: dict) -> argparse.ArgumentParser:
    def d(key, fallback):
        val = logic_defaults.get(key, fallback)
        return val

    p = argparse.ArgumentParser(
        prog="run_annotation",
        description="Wrapper for Annot_V5 with genome/reference auto-detection and flexible options.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("-i", "--input-dir", type=Path, required=True,
                   default=argparse.SUPPRESS, metavar="{INPUT_DIR}",
                   help="Project base directory; logic will auto-locate H5ADS_INPUT from here.")
    p.add_argument("-s", "--species", type=str, choices=["rice"], required=True,
                   default=argparse.SUPPRESS,
                   help="Species key for genome/reference auto-detect. Currently supports 'rice'.")
    p.add_argument("-cs", "--chrom-sizes", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Path to chrom.sizes.")
    p.add_argument("-gtf", "--gtf", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Path to gene annotation GTF.")

    # Reference overrides
    p.add_argument("--ref-csv", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Reference marker CSV.")
    p.add_argument("--ref-h5ad", type=Path, default=None, metavar="{PATH}",
                   help="(Optional) Reference atlas H5AD.")

    # Mapping flags
    p.add_argument("--label", type=str, metavar="{KEY}", default=argparse.SUPPRESS,
                   help="Map to REF_LABEL_KEY (default: celltype_id).")
    p.add_argument("--log2fc", type=float, metavar="{FLOAT}", default=argparse.SUPPRESS,
                   help="Map to MARKER_MIN_LOG2FC (default: 0.5).")
    p.add_argument("--pval", type=float, metavar="{FLOAT}", default=argparse.SUPPRESS,
                   help="Map to MARKER_MAX_PADJ (default: 0.05).")
    p.add_argument("--markerSet", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Map to MIN_MARKERS_PER_SET (default: 5).")
    p.add_argument("--feat", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Map to N_FEATURES (default: 250000).")
    p.add_argument("--neigh", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help="Map to N_NEIGHBORS (default: 10).")
    p.add_argument("--dim", type=parse_dims_arg, metavar="{INT/None}", default=argparse.SUPPRESS,
                   help="Map to USE_DIMS. Accepts '30' or '1,2,3' (default: None).")
    p.add_argument("--matrix", type=parse_bool, metavar="{TRUE/FALSE}", default=argparse.SUPPRESS,
                   help="Set MAKE_GENE_ACTIVITY (default: True).")
    p.add_argument("--mode", type=str, choices=["auc","weighted"], default=argparse.SUPPRESS,
                   help="Map to TRANSFER_MODE (default: auc).")
    p.add_argument("--coreTh", type=float, metavar="{FLOAT}", default=argparse.SUPPRESS,
                   help="Map to UNKNOWN_CONF_THRESH (default: 0.5).")
    p.add_argument("--umap", type=parse_bool, metavar="{TRUE/FALSE}", default=argparse.SUPPRESS,
                   help="Set SAVE_UMAP (default: True).")
    p.add_argument("--umapCount", type=int, metavar="{INT}", default=argparse.SUPPRESS,
                   help=f"Map to LABEL_MIN_COUNT (default: 10).")
    p.add_argument("--umapLabel", type=parse_bool, metavar="{TRUE/FALSE}", default=argparse.SUPPRESS,
                   help="Set LABEL_OUTLINE (default: True).")

    return p

# main
def main(argv: Optional[list] = None) -> int:
    project_root = Path(__file__).resolve().parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    # Import logic early to read defaults for help text
    logic_modname = None
    logic = None
    for mod in ("bin.Annot_V5", "bin.Annot_V5-1"):
        try:
            logic = importlib.import_module(mod)
            logic_modname = mod
            break
        except Exception:
            continue
    if logic is None:
        print("[ERROR] Could not import 'bin/Annot_V5'. Ensure bin/ is a package and file exists.", file=sys.stderr)
        return 2

    # Collect defaults from logic to show in help
    keys = [
        "REF_LABEL_KEY","MARKER_MIN_LOG2FC","MARKER_MAX_PADJ","MIN_MARKERS_PER_SET",
        "N_FEATURES","N_NEIGHBORS","USE_DIMS","MAKE_GENE_ACTIVITY",
        "TRANSFER_MODE","UNKNOWN_CONF_THRESH","SAVE_UMAP","LABEL_MIN_COUNT","LABEL_OUTLINE",
        "REF_CSV_PATH","REF_H5AD_PATH"
    ]
    logic_defaults = {k: getattr(logic, k, None) for k in keys}

    parser = build_parser(logic_defaults)
    args = parser.parse_args(argv)

    # Re-import logic in case argparse modified sys.path etc.
    logic = importlib.import_module(logic_modname) if logic_modname else importlib.import_module("bin.Annot_V5")

    # Inject base_dir / outputs
    base_dir_str = str(Path(args.input_dir).resolve())
    outdir_str = str((Path.cwd() / "Annot_results").resolve())
    try:
        if hasattr(logic, "BASE_DIR"):
            setattr(logic, "BASE_DIR", base_dir_str)
        if hasattr(logic, "base_dir"):
            setattr(logic, "base_dir", base_dir_str)
        if hasattr(logic, "AGGR_DIR"):
            setattr(logic, "AGGR_DIR", outdir_str)
        if hasattr(logic, "FIG_DIR"):
            setattr(logic, "FIG_DIR", str(Path(outdir_str) / "Plots"))
        if hasattr(logic, "SUMMARY_TXT"):
            setattr(logic, "SUMMARY_TXT", str(Path(outdir_str) / "summary.txt"))
        if hasattr(logic, "CSV_OUT"):
            setattr(logic, "CSV_OUT", str(Path(outdir_str) / "annotation_results.csv"))
        if hasattr(logic, "H5ADS_INPUT"):
            setattr(logic, "H5ADS_INPUT", str(Path(base_dir_str) / "Filter_results" / "merged_doublets.h5ads"))
    except Exception as e:
        print(f"[WARN] Failed to inject base/out vars: {e}", file=sys.stderr)

    # Genome: flags -> species -> auto-detect
    chrom_str = str(args.chrom_sizes.resolve()) if args.chrom_sizes else None
    gtf_str = str(args.gtf.resolve()) if args.gtf else None
    if not chrom_str or not gtf_str:
        c2, g2, sp2 = detect_species_genome(project_root, args.species)
        if c2 and g2:
            chrom_str = chrom_str or c2
            gtf_str = gtf_str or g2
            if args.species is None and sp2:
                print(f"[INFO] Auto-detected species: {sp2}")

    if not chrom_str or not gtf_str:
        print("[ERROR] Genome files not found. Provide -s rice or -cs <chrom.sizes> and -gtf <genes.gtf(.gz)>.", file=sys.stderr)
        return 3

    if hasattr(logic, "CHROM_SIZES_PATH"):
        setattr(logic, "CHROM_SIZES_PATH", chrom_str)
    if hasattr(logic, "GTF_PATH"):
        setattr(logic, "GTF_PATH", gtf_str)

    # References: flags -> auto-detect
    ref_csv = str(args.ref_csv.resolve()) if args.ref_csv else None
    ref_h5ad = str(args.ref_h5ad.resolve()) if args.ref_h5ad else None
    if not ref_csv or not ref_h5ad:
        c_csv, c_h5 = detect_reference_assets(project_root, args.species)
        ref_csv = ref_csv or (str(c_csv) if c_csv else None)
        ref_h5ad = ref_h5ad or (str(c_h5) if c_h5 else None)

    if ref_csv and hasattr(logic, "REF_CSV_PATH"):
        setattr(logic, "REF_CSV_PATH", ref_csv)
    if ref_h5ad and hasattr(logic, "REF_H5AD_PATH"):
        setattr(logic, "REF_H5AD_PATH", ref_h5ad)

    # Overrides only if user provided
    if hasattr(logic, "REF_LABEL_KEY") and getattr(args, "label", None) is not None:
        logic.REF_LABEL_KEY = str(args.label)
    if hasattr(logic, "MARKER_MIN_LOG2FC") and getattr(args, "log2fc", None) is not None:
        logic.MARKER_MIN_LOG2FC = float(args.log2fc)
    if hasattr(logic, "MARKER_MAX_PADJ") and getattr(args, "pval", None) is not None:
        logic.MARKER_MAX_PADJ = float(args.pval)
    if hasattr(logic, "MIN_MARKERS_PER_SET") and getattr(args, "markerSet", None) is not None:
        logic.MIN_MARKERS_PER_SET = int(args.markerSet)
    if hasattr(logic, "N_FEATURES") and getattr(args, "feat", None) is not None:
        logic.N_FEATURES = int(args.feat)
    if hasattr(logic, "N_NEIGHBORS") and getattr(args, "neigh", None) is not None:
        logic.N_NEIGHBORS = int(args.neigh)
    if hasattr(logic, "USE_DIMS") and getattr(args, "dim", None) is not None:
        logic.USE_DIMS = args.dim  # already parsed to int or list[int]
    if hasattr(logic, "MAKE_GENE_ACTIVITY") and getattr(args, "matrix", None) is not None:
        logic.MAKE_GENE_ACTIVITY = bool(args.matrix)
    if hasattr(logic, "TRANSFER_MODE") and getattr(args, "mode", None) is not None:
        logic.TRANSFER_MODE = str(args.mode)
    if hasattr(logic, "UNKNOWN_CONF_THRESH") and getattr(args, "coreTh", None) is not None:
        logic.UNKNOWN_CONF_THRESH = float(args.coreTh)
    if hasattr(logic, "SAVE_UMAP") and getattr(args, "umap", None) is not None:
        logic.SAVE_UMAP = bool(args.umap)
    if hasattr(logic, "LABEL_MIN_COUNT") and getattr(args, "umapCount", None) is not None:
        logic.LABEL_MIN_COUNT = int(args.umapCount)
    if hasattr(logic, "LABEL_OUTLINE") and getattr(args, "umapLabel", None) is not None:
        logic.LABEL_OUTLINE = bool(args.umapLabel)

    # Ensure output root exists
    Path(outdir_str).mkdir(parents=True, exist_ok=True)

    # Prefer logic.main(); fallback to logic.run_pipeline()
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
            print(f"[WARN] run_pipeline() failed: {e}", file=sys.stderr)

    if not ran:
        print("[INFO] No callable entry found in bin/Annot_V5.py. "
              "Variables have been set; please invoke the appropriate function inside the module.", file=sys.stderr)

    print(f"[OK] Finished.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
