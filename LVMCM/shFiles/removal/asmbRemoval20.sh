#!/bin/sh
#$ -cwd
#$ -pe smp 4 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=3:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-10
#$ -t 10

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing node removal assembly"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/removal/removal_pars_20.txt

for i in {0..9}
do
  /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE) $i
done

exit 0
