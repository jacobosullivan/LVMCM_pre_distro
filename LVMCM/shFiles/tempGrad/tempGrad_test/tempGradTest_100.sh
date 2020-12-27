#!/bin/sh
#$ -cwd
#$ -pe smp 8 # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=1:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-3

# parameterization, edited by LVMCM_SUB_A.sh
PARFILE="/data/home/btx188/LVMCM_src/LVMCM/Params/tempGrad/tempGrad_pars.csv"

# load required modules
module load singularity openmpi

echo "Initializing assembly"
SINGULARITYENV_OMP_NUM_THREADS=1 mpirun -np $NSLOTS singularity exec /data/home/btx188/containers/LVMCM_mpi.img \
  /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
