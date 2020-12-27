#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=3:00:00
#$ -l h_vmem=6G # memory request
#$ -t 1-2

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing Caaceres et al model fitting"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/caceres/caceres_pars.txt

for i in {0..9}
do
  /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE) $i
done

exit 0
