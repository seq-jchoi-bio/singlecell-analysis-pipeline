#!/bin/bash
#SBATCH --job-name=macs2_callpeak
#SBATCH --account=lab             
#SBATCH --partition=normal        
#SBATCH --cpus-per-task=4        
#SBATCH --mem=16G                
#SBATCH --time=24:00:00           
#SBATCH --output=logs/macs2_%j.out
#SBATCH --error=logs/macs2_%j.err

set -euo pipefail

# ===== Directories =====
INPUT_DIR="/home/sohyeong/projects/better_atac_pipeline/output_so/fragments"
OUTPUT_BASE="/home/sohyeong/projects/better_atac_pipeline/output_so/peakcalling_re"
LOG_BASE="${OUTPUT_BASE}/_logs_local"
mkdir -p "$OUTPUT_BASE" "$LOG_BASE"

# ===== Conda activate =====
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)" || true
fi
source /data/shared/init_conda.sh
conda activate /data/pipelines/atac_bulk/conda-env

if ! command -v macs2 >/dev/null 2>&1; then
  echo "ERROR: macs2 not found in conda env" >&2
  exit 1
fi
echo "macs2 version: $(macs2 --version 2>&1 || true)"

# ===== Run one sample =====
run_one_sample () {
  local CELLTYPE="$1"
  local i="$2"
  local SAMPLE_NAME="mcluster${CELLTYPE}_sample${i}"
  local TN5_BED="${INPUT_DIR}/mcluster${CELLTYPE}.sample${i}_fragments.bed"
  local OUTDIR="${OUTPUT_BASE}/${SAMPLE_NAME}/macs2_${SAMPLE_NAME}"
  local LOG_OUT="${LOG_BASE}/${SAMPLE_NAME}.out"
  local LOG_ERR="${LOG_BASE}/${SAMPLE_NAME}.err"

  mkdir -p "$OUTDIR"

  {
    echo "==== [$(date)] START ${SAMPLE_NAME} ===="
    echo "Input : $TN5_BED"
    echo "Output: $OUTDIR"

    if [ ! -s "$TN5_BED" ]; then
      echo "ERROR: Input file missing or empty: $TN5_BED" >&2
      return 2
    fi

    macs2 callpeak -n "${SAMPLE_NAME}" \
      -g hs \
      -p 0.01 \
      --nomodel \
      --shift -75 --extsize 150 \
      --keep-dup all \
      --nolambda \
      -f BED \
      -t "$TN5_BED" \
      -B --SPMR --call-summits \
      --buffer-size 1000000 \
      --outdir "$OUTDIR"

    for f in "${SAMPLE_NAME}_treat_pileup.bdg" "${SAMPLE_NAME}_control_lambda.bdg"; do
      [ -f "$OUTDIR/$f" ] && gzip -f "$OUTDIR/$f" || true
    done

    echo "==== [$(date)] DONE  ${SAMPLE_NAME} ===="
  } >"$LOG_OUT" 2>"$LOG_ERR"
}

# ===== Main loop =====
declare -a SUCCESS=()
declare -a FAILED=()

CELLTYPES=("Nervous_Cells" "Endothelial" "Fibroblasts" "Atrial_Cardiomyocytes" "Ventricular_Cardiomyocytes" "Primitive_Endoderm" "Myofibroblasts" "Smooth_Muscle" "Macrophages")

for CELLTYPE in "${CELLTYPES[@]}"; do
  for i in {1..9}; do
    NAME="mcluster${CELLTYPE}_sample${i}"
    echo "[INFO] Running ${NAME}..."
    if run_one_sample "$CELLTYPE" "$i"; then
      echo "[OK] ${NAME} succeeded."
      SUCCESS+=("${NAME}")
    else
      rc=$?
      echo "[FAIL] ${NAME} failed with rc=${rc}. See logs:"
      echo "       ${LOG_BASE}/${NAME}.out"
      echo "       ${LOG_BASE}/${NAME}.err"
      FAILED+=("${NAME}")
    fi
  done
done

# ===== Summary =====
echo
echo "===== SUMMARY ====="
echo "SUCCESS (${#SUCCESS[@]}): ${SUCCESS[*]:-none}"
echo "FAILED  (${#FAILED[@]}): ${FAILED[*]:-none}"
echo "Logs at: ${LOG_BASE}"
