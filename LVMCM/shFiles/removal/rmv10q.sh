#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=1:00:00
#$ -l h_vmem=6G # memory request
#$ -t 1-200
#$ -tc 20

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Running removal experiment (N=10 with environmental response function)"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/removal/bFileRemovalq10.txt

/data/home/btx188/LVMCM_src/LVMCM/build/LVMCM -b $(sed -n "${SGE_TASK_ID}p" $PARFILE) -K 0 -t 5000

exit 0
