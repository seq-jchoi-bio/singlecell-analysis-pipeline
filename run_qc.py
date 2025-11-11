#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Optional

HEADER = """
 =============================================================================
     Integrated scATAC-seq: Quality Control Pipeline
     Version 2.2
     (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*
 ============================================================================="""
desc_txt = """
This tool provides an integrated quality-control workflow for single-cell ATAC-seq data.
It automatically processes fragment files to generate fragment-size distribution plots,
TSSE (TSS enrichment) metrics, per-sample QC figures (as a single combined grid),
and summary outputs. All results are saved as PNG/SVG files, plus a concise TXT summary
with per-sample percentile statistics (p5/p10/median/p90/p95).
"""

def build_parser() -> argparse.ArgumentParser:
    import argparse
    p = argparse.ArgumentParser()

    p.add_argument("-i", "--input", type=str, required=False, default=None,
                   help="Input base directory; assigned to logic base_dir/BASE_DIR.")
    p.add_argument("-o", "--out", type=str, required=False, default=None,
                   help="Output directory; assigned to logic out_dir/OUT_DIR. Default: ./QC_results")
    p.add_argument("-s", "--species", type=str, required=False, default=None,
                   help="Species folder name under refGenome/. Must contain genes.gtf and genome.chrom.sizes.")
    p.add_argument("-j", "--n-jobs", type=int, default=None,
                   help="Parallel workers for per-sample tasks.")

    # hidden
    p.add_argument("-mrs", "--max-recorded-size", type=int, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("-sl", "--subsample-lines", type=int, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("-grid", "--grid-rows", type=int, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("-dpi", "--dpi", type=int, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("--make-html", action="store_true", default=True,
                   help=argparse.SUPPRESS)
    p.add_argument("--html-title", type=str, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("--embed-assets", action="store_true",
                   help=argparse.SUPPRESS)

    return p

def main(argv: Optional[list] = None) -> int:
    print(HEADER)
    print(desc_txt)
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(__file__).resolve().parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    try:
        logic = importlib.import_module("bin.IntegQC")
    except Exception as e:
        print(f"[ERROR] Could not import 'bin.IntegQC': {e}", file=sys.stderr)
        return 1

    if args.input:
        for name in ["base_dir", "BASE_DIR"]:
            setattr(logic, name, args.input)
    out_dir = args.out or None
    if out_dir:
        for name in ["out_dir", "OUT_DIR"]:
            setattr(logic, name, out_dir)

    # Species-specific: this logic will be removed!!!
    if args.species:
        base = project_root / "refGenome" / args.species
        chrom_path = base / "genome.chrom.sizes"
        gtf_path = base / "genes.gtf"
        if chrom_path.exists():
            for name in ["CHROM_SIZES_PATH", "chrom_sizes_path"]:
                setattr(logic, name, str(chrom_path))
        else:
            print(f"[WARN] Not found: {chrom_path}", file=sys.stderr)
        if gtf_path.exists():
            for name in ["GTF_PATH", "annot_path"]:
                setattr(logic, name, str(gtf_path))
        else:
            print(f"[WARN] Not found: {gtf_path}", file=sys.stderr)

    # Optional params
    if args.max_recorded_size is not None:
        logic.max_recorded_size = int(args.max_recorded_size)
    if args.subsample_lines is not None:
        logic.subsample_lines = int(args.subsample_lines)
    if args.grid_rows is not None:
        logic.grid_rows = int(args.grid_rows)
    if args.dpi is not None:
        logic.dpi_save = int(args.dpi)
    if args.n_jobs is not None:
        logic.n_jobs = int(args.n_jobs)

    ran = False
    if hasattr(logic, "main") and callable(logic.main):
        try:
            logic.main()
            ran = True
        except Exception as e:
            print(f"[WARN] logic.main() failed: {e}", file=sys.stderr)

    if not ran and hasattr(logic, "run_pipeline") and callable(logic.run_pipeline):
        try:
            logic.run_pipeline()
            ran = True
        except Exception as e:
            print(f"[WARN] run_pipeline invocation failed: {e}", file=sys.stderr)

    # HTML report
    if args.make_html:
        try:
            from bin import makeQCReport as report
            used_out = getattr(logic, "out_dir", None) or getattr(logic, "OUT_DIR", None)
            if used_out is None:
                used_out = str(Path.cwd() / "QC_results")
            report.build_html_report(
                out_dir=used_out,
                html_title=args.html_title or "QC Summary",
                embed_assets=args.embed_assets
            )
            print(f"[OK] HTML report generated at: {Path(used_out) / 'qc_summary.html'}")
        except Exception as e:
            print(f"[WARN] HTML report generation failed: {e}", file=sys.stderr)

    print("[OK] Finished.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
