#!/usr/bin/env python3

import sys
import argparse
import anndata as an

__version__ = "1.1"
__authors__ = "Janghyun Choi"

BANNER = f"""
========================================================================
  compare_h5ad.py — Compare var(features) between two .h5ad files
  Version {__version__} | Authors: {__authors__}
  Usage  : python compare_h5ad.py A.h5ad B.h5ad [--examples 10 --sort]
========================================================================
""".strip()

# Core
def load_vars(path):
    ad = an.read_h5ad(path)
    return ad.n_obs, ad.n_vars, set(map(str, ad.var_names))

def print_examples(title, items, k):
    print(f"\n--- {title} ---")
    if k <= 0 or not items:
        print("[]")
        return
    lst = list(items)[:k]
    print(lst)

# CLI
def build_argparser():
    desc = "Compare feature sets (var_names) between two .h5ad files."
    epilog = (
        "Examples:\n"
        "  python compare_h5ad.py A.h5ad B.h5ad\n"
        "  python compare_h5ad.py A.h5ad B.h5ad --examples 20 --sort\n"
    )
    ap = argparse.ArgumentParser(
        prog="compare_h5ad.py",
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog,
    )
    ap.add_argument("fileA", help="Path to SampleA .h5ad")
    ap.add_argument("fileB", help="Path to SampleB .h5ad")
    ap.add_argument("--examples", type=int, default=10, help="Number of example features to print (default: 10)")
    ap.add_argument("--sort", action="store_true", help="Sort features before showing examples")
    ap.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    return ap

def maybe_print_banner():
    argv = set(sys.argv[1:])
    if {"-h", "--help", "-v", "--version"} & argv:
        return
    print(BANNER)

def main():
    maybe_print_banner()
    ap = build_argparser()
    args = ap.parse_args()

    print(f"Loading {args.fileA} ...")
    nA_cells, nA_vars, setA = load_vars(args.fileA)
    print(f"  - cells: {nA_cells}, features: {nA_vars}")

    print(f"Loading {args.fileB} ...")
    nB_cells, nB_vars, setB = load_vars(args.fileB)
    print(f"  - cells: {nB_cells}, features: {nB_vars}")

    common = setA & setB
    onlyA  = setA - setB
    onlyB  = setB - setA

    if args.sort:
        setA = set(sorted(setA))
        setB = set(sorted(setB))
        common = set(sorted(common))
        onlyA  = set(sorted(onlyA))
        onlyB  = set(sorted(onlyB))

    print("\n=== Feature comparison ===")
    print(f"SampleA features: {len(setA)}")
    print(f"SampleB features: {len(setB)}")
    print(f"Common features : {len(common)}")
    print(f"A only features : {len(onlyA)}")
    print(f"B only features : {len(onlyB)}")

    k = max(0, int(args.examples))
    print_examples("Example features from SampleA", setA, k)
    print_examples("Example features from SampleB", setB, k)
    print_examples("Example common features", common, k)
    print_examples("Example A-only features", onlyA, k)
    print_examples("Example B-only features", onlyB, k)

if __name__ == "__main__":
    main()
