#!/bin/bash
#SBATCH --job-name=macs2_callpeak
#SBATCH --account=lab
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/macs2_%j.out
#SBATCH --error=logs/macs2_%j.err

set -eo pipefail

# ===== Paths =====
INPUT_DIR="$(pwd)/Annot_results/bed"
OUTPUT_BASE="$(pwd)/Peakcall_results"
LOG_BASE="${OUTPUT_BASE}/_logs_local"
mkdir -p "$OUTPUT_BASE" "$LOG_BASE"

# ===== Conda activation =====
source /data/shared/init_conda.sh
conda activate /data/pipelines/atac_bulk/conda-env

if ! command -v macs2 >/dev/null 2>&1; then
  echo "ERROR: macs2 not found in conda environment" >&2
  exit 1
fi
echo "[INFO] macs2 version: $(macs2 --version 2>&1 || true)"

# ===== Genome size mapping =====
GENOME=${GENOME:-os}   # default = rice
gsize=${GENOME:0:2}
if [[ "$gsize" == "hg" ]]; then
    gsize="hs"
elif [[ "$gsize" == "rn" ]]; then
    gsize="mm"
elif [[ "$gsize" == "os" ]]; then
    gsize="373245519"
else
    echo "[WARN] Unknown GENOME=${GENOME}, fallback hs"
    gsize="hs"
fi
echo "[INFO] GENOME=${GENOME}, using gsize=${gsize}"

# ===== Function: run peak calling =====
run_one_sample () {
  local FN="$1"
  local BASENAME
  BASENAME=$(basename "$FN")
  BASENAME=${BASENAME%.bed}
  BASENAME=${BASENAME%.bed.gz}
  local SAMPLE_NAME="$BASENAME"
  local OUTDIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
  local LOG_OUT="${LOG_BASE}/${SAMPLE_NAME}.out"
  local LOG_ERR="${LOG_BASE}/${SAMPLE_NAME}.err"

  mkdir -p "$OUTDIR"

  # 임시 BED3 파일 (seekable, 안전하게 /tmp 에 저장)
  local TMPBED
  TMPBED=$(mktemp "/tmp/${SAMPLE_NAME}.XXXXXX.bed3")

  {
    echo "==== [$(date)] START ${SAMPLE_NAME} ===="
    echo "INPUT  : $FN"
    echo "OUTDIR : $OUTDIR"

    if [ ! -s "$FN" ]; then
      echo "ERROR: Input file missing or empty: $FN" >&2
      return 2
    fi

    # 원본에서 앞 3컬럼만 추출
    cut -f1-3 "$FN" > "$TMPBED"

    local N=$(wc -l < "$TMPBED" || echo 0)
    echo "[INFO] BED3 written to $TMPBED, lines=$N"

    if [ "$N" -eq 0 ]; then
      echo "ERROR: Temp BED3 is empty: $TMPBED" >&2
      return 3
    fi

    echo "CMD    : macs2 callpeak -n ${SAMPLE_NAME} -g ${gsize} -t ${TMPBED} ..."

    macs2 callpeak \
      -n "${SAMPLE_NAME}" \
      -g "$gsize" \
      -p 0.05 \
      --nomodel \
      --shift -75 --extsize 150 \
      --keep-dup all \
      --nolambda \
      -f BED \
      -t "$TMPBED" \
      -B --SPMR --call-summits \
      --buffer-size 1000000 \
      --outdir "$OUTDIR"

    # compress auxiliary bdg files
    [ -f "$OUTDIR/${SAMPLE_NAME}_treat_pileup.bdg" ] && gzip -f "$OUTDIR/${SAMPLE_NAME}_treat_pileup.bdg"
    [ -f "$OUTDIR/${SAMPLE_NAME}_control_lambda.bdg" ] && gzip -f "$OUTDIR/${SAMPLE_NAME}_control_lambda.bdg"

    echo "==== [$(date)] DONE  ${SAMPLE_NAME} ===="
  } >"$LOG_OUT" 2>"$LOG_ERR"

  # 임시 파일 정리
  rm -f "$TMPBED" || true
}

# ===== Main loop =====
declare -a SUCCESS=()
declare -a FAILED=()

shopt -s nullglob
files=( "${INPUT_DIR}"/*.bed "${INPUT_DIR}"/*.bed.gz )
if ((${#files[@]} == 0)); then
  echo "ERROR: No BED files found in $INPUT_DIR" >&2
  exit 2
fi

for FN in "${files[@]}"; do
  NAME=$(basename "$FN")
  NAME=${NAME%.bed}
  NAME=${NAME%.bed.gz}
  echo "[INFO] Running ${NAME}..."
  if run_one_sample "$FN"; then
    SUCCESS+=("${NAME}")
  else
    rc=$?
    echo "[FAIL] ${NAME} failed (rc=${rc}). See logs:"
    echo "       ${LOG_BASE}/${NAME}.out"
    echo "       ${LOG_BASE}/${NAME}.err"
    FAILED+=("${NAME}")
  fi
done

# ===== Summary =====
echo
echo "===== SUMMARY ====="
echo "SUCCESS (${#SUCCESS[@]}): ${SUCCESS[*]:-none}"
echo "FAILED  (${#FAILED[@]}): ${FAILED[*]:-none}"
echo "Logs at: ${LOG_BASE}"
