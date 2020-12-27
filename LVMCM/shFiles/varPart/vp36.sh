#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=12:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-40

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing variance partitioning assembly"
PARFILE=/data/home/btx188/LVMCM_src_quad/LVMCM/parFiles/varPar/varPar_pars_N=36.txt
/data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
