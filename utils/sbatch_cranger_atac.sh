#!/bin/bash
#SBATCH --job-name=cranger_atac
#SBATCH --time=2-0
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=cpu
##SBATCH --account=lab -> by JChoi
#SBATCH --mem=80G

set -e

##
# usage: sbatch /data/programs/sbatch_cranger_atac.sh [SAMPLE_ID] [REFERENCE_DIR] [FASTQ_DIR] [ATAC_VER]
##

SAMPLE_ID=$1
REF_PATH=/data/programs/cellranger/refdata/$2
FASTQS_PATH=$3
ATAC_VER=$4

PROG_DIR=/data/programs/cellranger/atac-${ATAC_VER}

source $PROG_DIR/sourceme.bash
$PROG_DIR/cellranger-atac count \
  --id=$SAMPLE_ID \
  --reference=$REF_PATH \
  --fastqs=$FASTQS_PATH \
  --sample=$SAMPLE_ID \
  --localcores=12 --localmem=80

