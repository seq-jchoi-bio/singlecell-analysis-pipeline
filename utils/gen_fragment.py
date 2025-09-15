#!/usr/bin/env python3

import os, argparse, random, gzip, math, sys
from collections import defaultdict

__version__ = "1.5"
__authors__ = "Janghyun Choi"

HEADER = f"""
=================================================================================
  gen_fragment.py — Synthetic scATAC fragments generator (A/B)
  Version {__version__} | Authors: {__authors__}
  Summary : Generate realistic fragments.tsv(.gz) for two samples (A/B).
  Usage   : python gen_fragment.py --chrom-sizes genome.chrom.sizes --outdir out \\
             --mode mixed --overlap 0.4 --cells 100 --umi-per-cell 6000 [--gtf genes.gtf]
=================================================================================
""".strip()

# I/O
def open_w(path, gz):
    return gzip.open(path, "wt") if gz else open(path, "w")

def rand_barcode():
    return "".join(random.choice("ACGT") for _ in range(16)) + "-1"

# Genome / annotation utils
def read_chrom_sizes(path, min_len=5000, bin_size=500):
    chroms = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c, ln = line.split()[:2]
            ln = int(ln)
            if ln >= max(min_len, bin_size * 5):
                chroms.append((c, ln))
    if not chroms:
        raise SystemExit("No usable chromosomes in chrom.sizes.")
    return chroms

def read_tss_from_gtf(gtf_path):
    tss = defaultdict(list)
    try:
        with open(gtf_path) as f:
            for ln in f:
                if not ln.strip() or ln.startswith("#"):
                    continue
                cols = ln.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                chrom, _, feat, start, end, _, strand = cols[:7]
                feat = feat.lower()
                if feat not in ("gene", "transcript"):
                    continue
                try:
                    s = int(start); e = int(end)
                except Exception:
                    continue
                pos = s if strand == "+" else e
                tss[chrom].append(pos)
    except Exception:
        return {}
    return tss

def choose_center_len_residue(chrom_len, bin_size, residue, margin=2000):
    if chrom_len <= 2 * margin + bin_size:
        mid = max(margin + 1, min(chrom_len - margin - 1, chrom_len // 2))
    else:
        mid = random.randint(margin + 1, chrom_len - margin - 1)
    r = mid % bin_size
    center = mid - ((r - residue) % bin_size)
    if center < margin + 1:
        step = ((margin + 1 - center + bin_size - 1) // bin_size)
        center += bin_size * step
    if center > chrom_len - margin - 1:
        step = ((center - (chrom_len - margin - 1) + bin_size - 1) // bin_size)
        center -= bin_size * step
    return max(margin + 1, min(chrom_len - margin - 1, center))

def choose_center_near_tss(chrom, tss_map, chrom_len, jitter=120, margin=1000):
    if (chrom not in tss_map) or (not tss_map[chrom]):
        return None
    base = random.choice(tss_map[chrom])
    center = base + random.randint(-jitter, jitter)
    center = max(margin + 1, min(chrom_len - margin - 1, center))
    return center

# Fragment-length samplers
def sample_len_uniform(Lmin, Lmax):
    return random.randint(Lmin, Lmax)

def sample_len_mixture(nfr_mu=50, nfr_sd=12,
                       mono_mu=180, mono_sd=20,
                       di_mu=360, di_sd=30,
                       tail_mean=120, tail_base=180,
                       weights=(0.55, 0.30, 0.11, 0.04),
                       nuc_wiggle=2, Lmin=20, Lmax=1000):
    r = random.random()
    w1, w2, w3, w4 = weights
    if r < w1:
        L = int(random.gauss(nfr_mu, nfr_sd))
    elif r < w1 + w2:
        L = int(random.gauss(mono_mu, mono_sd))
    elif r < w1 + w2 + w3:
        L = int(random.gauss(di_mu, di_sd))
    else:
        L = tail_base + int(random.expovariate(1 / float(tail_mean)))
    if nuc_wiggle:
        L += int(nuc_wiggle * math.sin(L * (2 * math.pi / 10.4)))
    return max(Lmin, min(Lmax, L))

# Core writer
def write_sample(
    path, chroms, cells, umi_per_cell,
    mode, bin_size, overlap, residue_a=0, residue_b_shift=1,
    tss_map=None, tss_fraction=0.0, gzip_out=False, seed=7,
    shared_pool=None, collect_shared=False,
    umi_jitter=0.0, count_min=1, count_max=5,
    len_model="mixture", len_min=60, len_max=120,
    nuc_wiggle=2
):

    random.seed(seed)
    with open_w(path, gzip_out) as f:
        barcodes = [rand_barcode() for _ in range(cells)]
        p_shared = 1.0 if mode == "common" else (overlap if mode == "mixed" else 0.0)

        for bc in barcodes:
            if umi_jitter > 0:
                scale = random.lognormavariate(0.0, umi_jitter) if hasattr(random, "lognormavariate") else random.lognormvariate(0.0, umi_jitter)  # Python compat
                target_umi = max(200, int(round(umi_per_cell * scale)))
            else:
                target_umi = int(umi_per_cell)

            remain = target_umi
            chrom_idx = 0

            while remain > 0:
                chrom, L = chroms[chrom_idx % len(chroms)]
                batch_lines = max(1, min(50, (remain + count_max - 1) // max(1, count_max)))

                for _ in range(batch_lines):
                    if remain <= 0:
                        break

                    use_shared = (random.random() < p_shared)
                    center = None

                    if (mode == "common") or (mode == "mixed" and use_shared):
                        if tss_map and (random.random() < tss_fraction):
                            center = choose_center_near_tss(chrom, tss_map, L)
                        if center is None:
                            center = choose_center_len_residue(L, bin_size, residue=residue_a, margin=2000)
                        if collect_shared:
                            shared_pool[chrom].append(center)
                        elif shared_pool is not None and shared_pool[chrom]:
                            center = random.choice(shared_pool[chrom])
                    else:
                        if tss_map and (random.random() < tss_fraction):
                            center = choose_center_near_tss(chrom, tss_map, L)
                        if center is None:
                            alt_res = (residue_a + residue_b_shift) % bin_size
                            center = choose_center_len_residue(L, bin_size, residue=alt_res, margin=2000)

                    if len_model == "uniform":
                        frag_len = sample_len_uniform(len_min, len_max)
                    else:
                        frag_len = sample_len_mixture(nuc_wiggle=nuc_wiggle, Lmin=20, Lmax=1000)

                    s = max(0, center - frag_len // 2)
                    e = min(L, s + frag_len)
                    if e <= s:
                        e = s + 1
                    cnt = min(remain, random.randint(max(1, count_min), max(1, count_max)))
                    strand = random.choice(["+", "-"])
                    f.write(f"{chrom}\t{s}\t{e}\t{bc}\t{cnt}\t{strand}\n")
                    remain -= cnt

                chrom_idx += 1

# CLI
def build_argparser():
    desc = (
        "Generate realistic fragments.tsv(.gz) for two samples (A/B), "
        "supporting overlap modes: common, disjoint, mixed."
    )
    epilog = (
        "Examples:\n"
        "  1) Common bins (identical tiles)\n"
        "     python gen_fragment.py --chrom-sizes genome.chrom.sizes --outdir toy --mode common\n\n"
        "  2) Disjoint bins (no overlap)\n"
        "     python gen_fragment.py --chrom-sizes genome.chrom.sizes --mode disjoint --cells 200\n\n"
        "  3) Mixed overlap with TSS bias\n"
        "     python gen_fragment.py --chrom-sizes genome.chrom.sizes --mode mixed --overlap 0.5 \\\n"
        "                            --gtf genes.gtf --tss-fraction 0.2 --gzip"
    )
    ap = argparse.ArgumentParser(
        prog="gen_fragment.py",
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog
    )
    ap.add_argument("--chrom-sizes", required=True, help="chrom.sizes file (two columns: chrom length)")
    ap.add_argument("--outdir", default="toy_cells", help="output directory")
    ap.add_argument("--bin", type=int, default=500, help="tile bin size controlling overlap logic")
    ap.add_argument("--cells", type=int, default=100, help="cells per sample (unique barcodes)")
    ap.add_argument("--umi-per-cell", type=int, default=6000, help="target fragments per cell before jitter")
    ap.add_argument("--umi-jitter", type=float, default=0.3, help="lognormal sigma for per-cell UMI variation (0 disables)")
    ap.add_argument("--count-min", type=int, default=1, help="min count per line (count column)")
    ap.add_argument("--count-max", type=int, default=5, help="max count per line (count column)")
    ap.add_argument("--mode", choices=["common", "disjoint", "mixed"], default="common", help="overlap mode")
    ap.add_argument("--overlap", type=float, default=0.4, help="mixed: fraction of A's bins reused by B [0..1]")
    ap.add_argument("--gtf", help="optional GTF; if set, some fragments are placed near TSS")
    ap.add_argument("--tss-fraction", type=float, default=0.0, help="fraction of fragments near TSS (0..1)")
    ap.add_argument("--gzip", action="store_true", help="write .tsv.gz instead of .tsv")
    ap.add_argument("--seed", type=int, default=7, help="base random seed")
    ap.add_argument("--len-model", choices=["mixture", "uniform"], default="mixture", help="fragment-length model")
    ap.add_argument("--len-min", type=int, default=60, help="uniform model: minimum length")
    ap.add_argument("--len-max", type=int, default=120, help="uniform model: maximum length")
    ap.add_argument("--nuc-wiggle", type=int, default=2, help="mixture model: 10.4-bp wiggle amplitude (0 disables)")
    ap.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    return ap

def print_header_if_needed():
    argv = set(sys.argv[1:])
    if {"-h", "--help", "-v", "--version"} & argv:
        return
    print(HEADER)

def main():
    print_header_if_needed()
    ap = build_argparser()
    args = ap.parse_args()

    if args.count_min < 1:
        args.count_min = 1
    if args.count_max < args.count_min:
        args.count_max = args.count_min

    random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    chroms = read_chrom_sizes(args.chrom_sizes, min_len=5000, bin_size=args.bin)
    tss_map = read_tss_from_gtf(args.gtf) if args.gtf else None

    outA = os.path.join(args.outdir, "fragments_A.tsv" + (".gz" if args.gzip else ""))
    outB = os.path.join(args.outdir, "fragments_B.tsv" + (".gz" if args.gzip else ""))

    shared_pool = defaultdict(list)

    write_sample(
        outA, chroms,
        cells=args.cells, umi_per_cell=args.umi_per_cell,
        mode=args.mode, bin_size=args.bin, overlap=args.overlap,
        residue_a=0, residue_b_shift=1,
        tss_map=tss_map, tss_fraction=(args.tss_fraction if tss_map else 0.0),
        gzip_out=args.gzip, seed=args.seed,
        shared_pool=shared_pool, collect_shared=(args.mode in ("common", "mixed")),
        umi_jitter=max(0.0, args.umi_jitter),
        count_min=args.count_min, count_max=args.count_max,
        len_model=args.len_model, len_min=args.len_min, len_max=args.len_max,
        nuc_wiggle=args.nuc_wiggle
    )

    write_sample(
        outB, chroms,
        cells=args.cells, umi_per_cell=args.umi_per_cell,
        mode=args.mode, bin_size=args.bin, overlap=args.overlap,
        residue_a=0, residue_b_shift=1,
        tss_map=tss_map, tss_fraction=(args.tss_fraction if tss_map else 0.0),
        gzip_out=args.gzip, seed=args.seed + 1,
        shared_pool=shared_pool, collect_shared=False,
        umi_jitter=max(0.0, args.umi_jitter),
        count_min=args.count_min, count_max=args.count_max,
        len_model=args.len_model, len_min=args.len_min, len_max=args.len_max,
        nuc_wiggle=args.nuc_wiggle
    )

    print(f"[OK] wrote:\n  {outA}\n  {outB}")
    print(f"  mode={args.mode}, bin={args.bin}, cells={args.cells}, umi/cell={args.umi_per_cell}, "
          f"umi_jitter={args.umi_jitter}, count=[{args.count_min},{args.count_max}], "
          f"len_model={args.len_model}, nuc_wiggle={args.nuc_wiggle}")
    if args.gtf:
        print(f"  TSS bias: {args.tss_fraction:.2f} near TSS (from {os.path.basename(args.gtf)})")

if __name__ == "__main__":
    main()
