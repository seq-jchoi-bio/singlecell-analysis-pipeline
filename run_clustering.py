#!/usr/bin/env python3

import argparse, os, sys, subprocess
from pathlib import Path

HEADER = """
 =============================================================================
     Integrated scATAC-seq: Clustering
     Version 1.3
     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*
 ============================================================================="""

desc_txt = """This code sequentially runs three modules for integrated scATAC-seq clustering:
  1) bin/CLS.py       → Performs feature selection, spectral embedding, and Leiden clustering.
  2) bin/postCLS.py   → Generates gene activity matrices using the given genome reference.
  3) bin/dotplot.py   → Produces marker-based dotplots for annotated clusters.
All intermediate and final results are stored under Annot_results/, with plots saved
in Annot_results/plots and dotplots in Annot_results/dotplots.
"""

def main():
    print(HEADER)
    print(desc_txt)

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input project directory (must contain Filter_results/*.h5ads)")
    ap.add_argument("-o", "--output", default=None, help="Output directory (default: ./Annot_results)")
    ap.add_argument("-s", "--species", required=True, help="Reference genome (e.g., hg38, mm10, or rice)")
    ap.add_argument("-r", "--markers", required=True, help="Marker CSV file for dotplot")
    args = ap.parse_args()
    
    proj = Path(args.input).resolve()
    out_root = Path(args.output).resolve() if args.output else proj

    rc = subprocess.call([sys.executable, str(Path(__file__).parent / "bin" / "CLS.py"),
                          "-i", str(proj), "-o", str(out_root)], cwd=str(proj))
    if rc != 0: sys.exit(1)
    rc = subprocess.call([sys.executable, str(Path(__file__).parent / "bin" / "postCLS.py"),
                          "-i", str(proj), "-o", str(out_root), "-s", args.species], cwd=str(proj))
    if rc != 0: sys.exit(1)
    rc = subprocess.call([sys.executable, str(Path(__file__).parent / "bin" / "dotplot.py"),
                          "-i", str(proj), "-o", str(out_root), "-r", str(Path(args.markers).resolve())], cwd=str(proj))
    if rc != 0: sys.exit(1)

    print("All steps completed.")

if __name__ == "__main__":
    main()