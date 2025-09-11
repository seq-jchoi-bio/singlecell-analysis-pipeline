#!/bin/bash
# run_macs2_sequential.sh
# 지정된 세포타입 × sample1~9 조합에 대해 macs2 peak calling을 순차적으로 실행하는 스크립트
# 실행 순서: Nervous_Cells → Endothelial → Fibroblasts → … (리스트 순서대로)
# 실패한 경우에도 전체 루프는 계속 진행하며, 결과는 SUCCESS / FAILED 로 요약 출력

set -euo pipefail
# -e : 에러 발생 시 즉시 종료 (단, 함수 실행 후 루프 내에서 에러를 잡아 처리하므로 전체는 계속 진행)
# -u : 정의되지 않은 변수를 사용하면 에러
# -o pipefail : 파이프라인 중 하나라도 실패 시 에러로 간주

# ===== 경로 설정 =====
INPUT_DIR="/home/sohyeong/projects/better_atac_pipeline/output_so/tn5_ins"   # 입력 BED 파일 경로
OUTPUT_BASE="/home/sohyeong/projects/better_atac_pipeline/output_so/peakcalling" # 결과 저장 기본 디렉토리
LOG_BASE="${OUTPUT_BASE}/_logs_local"                                       # 로그 저장 디렉토리
mkdir -p "$OUTPUT_BASE" "$LOG_BASE"                                         # 필요시 디렉토리 생성

# ===== Conda 활성화 =====
# set -u 상태에서는 정의되지 않은 변수 사용 시 에러가 나므로, 일시적으로 +u로 해제
if command -v conda >/dev/null 2>&1; then
  set +u
  eval "$(conda shell.bash hook)" || true   # conda 초기화
  set -u
fi
set +u
source /data/shared/init_conda.sh          # 공용 초기화 스크립트
conda activate /data/pipelines/atac_bulk/conda-env   # ATAC 전용 conda 환경 활성화
set -u

# macs2 실행 가능 여부 확인
if ! command -v macs2 >/dev/null 2>&1; then
  echo "ERROR: macs2 not found in conda env" >&2
  exit 1
fi
echo "macs2 version: $(macs2 --version 2>&1 || true)"  # macs2 버전 출력

# ===== 실행 함수 =====
run_one_sample () {
  local CELLTYPE="$1"        # 첫 번째 인자: 세포 타입
  local i="$2"               # 두 번째 인자: 샘플 번호 (1~9)
  local SAMPLE="sample${i}"  # 샘플 이름 (예: sample1)
  local SAMPLE_NAME="mcluster${CELLTYPE}_${SAMPLE}"  # mcluster 접두 + 세포타입 + 샘플
  local TN5_BED="${INPUT_DIR}/mcluster${CELLTYPE}.${SAMPLE}_tn5_ins.bed" # 입력 tn5 insertion BED
  local OUTDIR="${OUTPUT_BASE}/${SAMPLE_NAME}/macs2_${SAMPLE_NAME}"      # macs2 결과 디렉토리
  local LOG_OUT="${LOG_BASE}/${SAMPLE_NAME}.out"                         # 표준출력 로그
  local LOG_ERR="${LOG_BASE}/${SAMPLE_NAME}.err"                         # 표준에러 로그

  mkdir -p "$OUTDIR" "$LOG_BASE"   # 출력 및 로그 디렉토리 생성

  {
    echo "==== [$(date)] START ${SAMPLE_NAME} ===="
    echo "INPUT  : $TN5_BED"
    echo "OUTDIR : $OUTDIR"

    # 입력 파일이 존재하지 않거나 비어 있으면 에러 처리
    if [ ! -s "$TN5_BED" ]; then
      echo "ERROR: 입력 파일이 없거나 비어있음: $TN5_BED" >&2
      return 2
    fi

    # MACS2 실행 (ATAC-seq 권장 파라미터 사용)
    macs2 callpeak -n "${SAMPLE_NAME}" \
      -g hs \                                # reference genome: human
      -p 0.01 \                              # p-value cutoff
      --nomodel \                            # model 기반 추정 비활성화 (ATAC 전용)
      --shift -75 --extsize 150 \            # ATAC 권장 shift/extsize
      --keep-dup all \                       # 모든 중복 리드 유지
      --nolambda \                           # local lambda 추정 비활성화
      -f BED \                               # 입력 포맷: BED
      -t "$TN5_BED" \                        # 입력 파일
      -B --SPMR --call-summits \             # bedGraph 출력, signal per million reads, summit 호출
      --buffer-size 1000000 \                # 메모리 버퍼 크기 조정
      --outdir "$OUTDIR"                     # 결과 출력 디렉토리

    # coverage bedGraph 파일은 용량 절약을 위해 압축 (존재할 경우만)
    [ -f "$OUTDIR/${SAMPLE_NAME}_treat_pileup.bdg" ] && gzip -f "$OUTDIR/${SAMPLE_NAME}_treat_pileup.bdg" || true
    [ -f "$OUTDIR/${SAMPLE_NAME}_control_lambda.bdg" ] && gzip -f "$OUTDIR/${SAMPLE_NAME}_control_lambda.bdg" || true

    echo "==== [$(date)] DONE  ${SAMPLE_NAME} ===="
  } >"$LOG_OUT" 2>"$LOG_ERR"   # 표준출력과 표준에러를 각각 로그 파일로 저장
}

# ===== 메인 루프 =====
declare -a SUCCESS=()   # 성공한 샘플 리스트
declare -a FAILED=()    # 실패한 샘플 리스트

# 실행할 세포타입 순서 지정 (원하는 순서대로 실행됨)
CELLTYPES=("Nervous_Cells" "Endothelial" "Fibroblasts" "Atrial_Cardiomyocytes" "Ventricular_Cardiomyocytes" "Primitive_Endoderm" "Myofibroblasts" "Smooth_Muscle" "Macrophages")

# 셀타입별, 샘플별 순차 실행
for CELLTYPE in "${CELLTYPES[@]}"; do
  for i in {1..9}; do
    SAMPLE="sample${i}"
    NAME="mcluster${CELLTYPE}_${SAMPLE}"
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

# ===== 요약 출력 =====
echo
echo "===== SUMMARY ====="
echo "SUCCESS (${#SUCCESS[@]}): ${SUCCESS[*]:-none}"   # 성공한 샘플 이름 목록
echo "FAILED  (${#FAILED[@]}): ${FAILED[*]:-none}"     # 실패한 샘플 이름 목록
echo "Logs at: ${LOG_BASE}"                            # 로그 경로 안내