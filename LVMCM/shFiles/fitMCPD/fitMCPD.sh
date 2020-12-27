#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l h_rt=1:00:00
#$ -l h_vmem=2G # memory request
#$ -t 1-900
#$ -tc 40

# load required modules
echo "Estimating mixing parameter from LVMCM"
module load R
PARFILE=/data/SBCS-RossbergLab/jack/post_doc/removal/bFileRemoval.txt

Rscript /data/home/btx188/LVMCM_src/LVMCM/RCode/FitMCPD_exec.R $(sed -n "${SGE_TASK_ID}p" $PARFILE) 0 1 A ml

exit 0
