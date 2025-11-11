#!/bin/bash
# peakcall for scATAC via macs2 (BED-only, per-file subfolders)

#SBATCH --job-name=macs2_scATAC
#SBATCH --time=23:59:59
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -eo pipefail

usage() {
  cat <<'EOF'
Usage:
  sbatch macs2_scATAC.sh <GENOME> <INPUT_DIR> <OUTDIR>

Arguments:
  GENOME     Genome key (e.g., hg19/hg38 → hs, rice/os → 373245519)
  INPUT_DIR  Directory containing *.bed files (non-recursive)
  OUTDIR     Output root directory (must not exist)

Examples:
  sbatch macs2_scATAC.sh rice /data/bed_inputs   /data/macs2_out

Notes:
  - Only *.bed files are processed (no *.bed.gz).
  - OUTDIR must not exist; the script refuses to overwrite.
  - Each BED produces results under <OUTDIR>/<basename>/
Options:
  -h, --help, help   Show this help and exit
EOF
}

# help switch
case "${1:-}" in
  -h|--help|help) usage; exit 0 ;;
esac

# load conda environment
source /data/shared/init_conda.sh
conda activate /data/pipelines/atac_bulk/conda-env

# variables
GENOME="${1:-}"
INPUT_DIR="${2:-}"
OUTDIR="${3:-}"

# missing check
if [[ -z "${GENOME}" || -z "${INPUT_DIR}" || -z "${OUTDIR}" ]]; then
  usage
  exit 1
fi
if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "[ERROR] INPUT_DIR not found: ${INPUT_DIR}" >&2
  exit 1
fi
if [[ -e "${OUTDIR}" ]]; then
  echo "[ERROR] OUTDIR already exists: ${OUTDIR} (refuse to overwrite)" >&2
  exit 1
fi
mkdir -p "${OUTDIR}"

# species params -> need to improve....
gsize="${GENOME:0:2}"
if [[ "$gsize" == "hg" ]]; then
  gsize="hs"
elif [[ "$gsize" == "rn" ]]; then
  gsize="mm"
elif [[ "$gsize" == "ri" ]]; then
  gsize="373245519"
fi
echo "[INFO] macs2 -g: ${gsize}"

# auto detection BED
mapfile -t BED_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -iname "*.bed" | sort)
if [[ "${#BED_FILES[@]}" -eq 0 ]]; then
  echo "[ERROR] No .bed files found under: ${INPUT_DIR}" >&2
  exit 1
fi
echo "[INFO] Found ${#BED_FILES[@]} BED file(s)."

# Running
for bed in "${BED_FILES[@]}"; do
  bn="$(basename "${bed}")"
  sample="${bn%.bed}"
  sub_out="${OUTDIR}/${sample}"
  mkdir -p "${sub_out}"

  echo "[RUN] macs2 for ${bn} -> ${sub_out}"
  macs2 callpeak -n "${sample}" \
    -g "${gsize}" --nomodel --nolambda \
    --shift -75 --extsize 150 \
    --keep-dup all -f BED -t "${bed}" \
    -B --SPMR --call-summits \
    --buffer-size 1000000 \
    --outdir "${sub_out}"

  [[ -f "${sub_out}/${sample}_treat_pileup.bdg" ]] && gzip -f "${sub_out}/${sample}_treat_pileup.bdg"
  [[ -f "${sub_out}/${sample}_control_lambda.bdg" ]] && gzip -f "${sub_out}/${sample}_control_lambda.bdg"

  echo "[OK] ${sample}"
done

echo "[DONE] All BED files processed into: ${OUTDIR}"
