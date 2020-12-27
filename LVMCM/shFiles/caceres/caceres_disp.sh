#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=30:00:00
#$ -l h_vmem=6G # memory request
#$ -t 1-90

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing Caaceres et al model fitting"
let "SMO=$SGE_TASK_ID-1"
let "PARLINE=1+($SMO-$SMO%10)/10"
echo $PARLINE
let "REP=($SGE_TASK_ID%10)-1"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/caceres/caceres_disp_pars.txt

/data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${PARLINE}p" $PARFILE) $REP

exit 0
