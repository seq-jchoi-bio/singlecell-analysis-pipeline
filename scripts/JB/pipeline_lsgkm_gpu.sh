#!/bin/bash

#SBATCH --job-name=get_fasta
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=gpu
#SBATCH --account=lab
#SBATCH --mem=4G

#set -o errexit
#set -o nounset

BASEDIR=/data/pipelines/dsvm

source /data/programs/conda/bin/activate $BASEDIR/data/conda-env

inputf=$1   # bed file
extlen=$2   # PEAKLEN/2 (ie. 150 for 300bp, 300 for 600bp)
column=$3
topn=$4     # number of top peaks for training
genome=$5   # genome version : 나는 IRGSP
rseed=$6

T=4
L=11
K=7
D=3
NCV=5

expid=${inputf%.bed}.t${topn}
posbed=${expid}.pos${rseed}.bed
negbed=${expid}.neg${rseed}.bed
posfa=${expid}.pos${rseed}.fa
negfa=${expid}.neg${rseed}.fa
svmid=${expid}.r${rseed}.lsgkm.${T}.${L}.${K}.${D}

# 1. sort out positive set
sort -gr -k${column},${column} ${inputf} | head -n ${topn} | sortBed > ${posbed}

conda deactivate

# 2. generate negative set and retrieve fasta file
gqcpath=/data/pipelines/gkmqc

conda activate ${gqcpath}/conda-env/
python ${gqcpath}/scripts/seqs_nullgen.py -p ${posbed} -n ${negbed} -g ${genome} -t $((extlen * 2)) -s ${rseed} -@ 16
conda deactivate

# 3. train 5-CV for just one set
if [ $rseed == 1 ]; then
	restxt=`sbatch ${BASEDIR}/scripts/sbatch_gkmtrain.sh -t $T -l $L -k $K -d $D -x $NCV -i 1 $posfa $negfa $svmid`
	restxtarr=($restxt)
	jobid1=${restxtarr[3]}
fi
# 4. train model using all data
sbatch --export=NONE --dependency=afterany:${jobid1} ${BASEDIR}/scripts/sbatch_gkmtrain.sh -t $T -l $L -k $K -d $D $posfa $negfa $svmid
