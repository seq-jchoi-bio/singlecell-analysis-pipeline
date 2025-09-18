#!/bin/bash
#Wrapper script for sbatch with 4-threads

#SBATCH --job-name=gkmpred4
#SBATCH --time=20:59:59
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=sampson
#SBATCH --partition=sampson-compute
#SBATCH --mem=20G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dwlee@jhu.edu

set -o errexit
set -o nounset

plineDir="/data/pipelines/dsvm"
$plineDir/bin/gkmpredict2 -T 4 "$@"
