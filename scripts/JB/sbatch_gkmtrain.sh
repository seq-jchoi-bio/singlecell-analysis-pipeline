#!/bin/bash

#Wrapper script for sbatch with 16-threads
# For slurm comptability
#SBATCH --job-name=gkmt8_100G_gpu
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=lab
#SBATCH --partition=gpu
#SBATCH --mem=100G

#set -o errexit
#set -o nounset

plineDir="/data/pipelines/dsvm"
$plineDir/bin/gkmtrain2 -T 16 -m 200000 "$@"
