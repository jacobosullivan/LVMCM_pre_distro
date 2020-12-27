#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=1:00:00
#$ -l h_vmem=6G # memory request
#$ -t 1-20

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing Caaceres et al removal experiment"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/caceres/bFileList.txt

/data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE) -K 0 -t 5000 -sc "/data/home/btx188/LVMCM_src/LVMCM/fitNetwork/SC_Caceres2017_med.mat"

exit 0
