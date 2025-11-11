#!/usr/bin/env python3

import argparse, sys, subprocess
from pathlib import Path

HEADER = """
 =============================================================================
     Integrated scATAC-seq: Annotation and Reporting
     Version 4.0
     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*
 =============================================================================
"""

desc_txt = """This script sequentially executes two integrated modules for scATAC-seq downstream analysis:
  1) bin/annotate.py     → Updates per-sample h5ad files with Leiden and cell-type information,
                           merges them into a unified AnnDataSet, and exports BED files for each group.
  2) bin/FinalReport.py  → Generates a fully embedded HTML report summarizing clustering,
                           annotation results, and marker-based dotplots.
"""

def _auto_base_dir(p: Path) -> Path:
    p = p.resolve()
    cur = p
    for _ in range(6):
        if (cur / "Filter_results").is_dir() and (cur / "Annot_results").is_dir():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return p

def main():
    print(HEADER)
    print(desc_txt)
    
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True, help="Project root (or any child path)")
    ap.add_argument("-o","--output-dir", default=None, help="Report output directory (default: project root)")
    ap.add_argument("-n","--name", default="FinalReport.html", help="Report file name (default: FinalReport.html)")
    ap.add_argument("-r","--marker", default=None, help="Marker CSV (columns: type + gene_name/gene_id/gene)")
    args = ap.parse_args()

    # Paths relative to this script
    script_dir = Path(__file__).resolve().parent
    annotate_py = script_dir / "bin" / "annotate.py"
    finalreport_py = script_dir / "bin" / "FinalReport.py"

    # Sanity checks
    if not annotate_py.is_file():
        raise FileNotFoundError(f"annotate.py not found: {annotate_py}")
    if not finalreport_py.is_file():
        raise FileNotFoundError(f"FinalReport.py not found: {finalreport_py}")

    base = _auto_base_dir(Path(args.input))
    outdir = Path(args.output_dir) if args.output_dir else base

    print(f"[Preparing...] Detection base DIR: {base}")
    print("====== Step 1. Annotate and update h5ad (bin/annotate.py)====== ")
    subprocess.check_call([sys.executable, str(annotate_py), "-b", str(base)])

    print("====== Step 2. Make report (bin/FinalReport.py)======")
    cmd2 = [sys.executable, str(finalreport_py), "-b", str(base), "-o", str(outdir), "-n", args.name]
    if args.marker:
        cmd2 += ["-r", str(Path(args.marker).resolve())]
    subprocess.check_call(cmd2)

    print("All steps completed.")

if __name__ == "__main__":
    main()
