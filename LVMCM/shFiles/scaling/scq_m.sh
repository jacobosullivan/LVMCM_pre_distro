#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=1:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-9
#$ -tc 10

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing interaction coefficient scaling assembly"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/scaling/scaling_quad_pars_missing.txt

/data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)


exit 0
