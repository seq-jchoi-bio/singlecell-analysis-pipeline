#!/bin/bash
#SBATCH --job-name=fetch_fastq
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --partition=cpu
#SBATCH --mem=4G

#set -o errexit
#set -o nounset
#set -o pipefail

srid="$1"    # e.g., SRR12345678
odir="$2"    # output directory
sname="$3"   # desired sample name used in FASTQ file names

# Constants for 10x-style file naming
lane="L001"
sidx="S1"
chunk="001"

# Load environment for SRA Toolkit (BioGrids)
source /programs/biogrids.shrc >/dev/null

# Download and convert SRA to FASTQ
prefetch --max-size unlimited "$srid" -O "$odir" # Unlimited check!!!!
fasterq-dump "$odir/$srid/${srid}.sra" \
  --outdir "$odir/$srid" \
  --split-files \
  --include-technical

# Gzip compress
for i in 1 2 3 4; do
  f="$odir/$srid/${srid}_$i.fastq"
  if [[ -f "$f" ]]; then
    gzip -f "$f"
  fi
done

# Rename to 10X pattern:
#    _1 -> I1, _2 -> R1, _3 -> R2, _4 -> R3
declare -A MAP
MAP["${srid}_1.fastq.gz"]="${sname}_${sidx}_${lane}_I1_${chunk}.fastq.gz"
MAP["${srid}_2.fastq.gz"]="${sname}_${sidx}_${lane}_R1_${chunk}.fastq.gz"
MAP["${srid}_3.fastq.gz"]="${sname}_${sidx}_${lane}_R2_${chunk}.fastq.gz"
MAP["${srid}_4.fastq.gz"]="${sname}_${sidx}_${lane}_R3_${chunk}.fastq.gz"

for src in "${!MAP[@]}"; do
  src_path="$odir/$srid/$src"
  dst_path="$odir/$srid/${MAP[$src]}"
  if [[ -f "$src_path" ]]; then
    mv -f "$src_path" "$dst_path"
  fi
done

echo "Renamed files in: $odir/$srid"
ls -1 "$odir/$srid" | grep -E "${sname}_${sidx}_${lane}_(I1|R1|R2|R3)_${chunk}\.fastq\.gz" || true

