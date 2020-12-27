#!/bin/sh
#$ -cwd
#$ -pe smp 2 # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=0:30:00
#$ -l h_vmem=1G # memory request
#$ -t 31

# parameterization, edited by LVMCM_SUB_A.sh
PARFILE="/data/home/btx188/LVMCM_src_mm/LVMCM/parFiles/perfTest/perfTest_pars_gcc.csv"

# load required modules
module load singularity openmpi

echo "Initializing performance test"
SINGULARITYENV_OMP_NUM_THREADS=1 mpirun -np $NSLOTS singularity exec /data/home/btx188/containers/LVMCM_mpi.img \
  /data/home/btx188/LVMCM_src_mm/LVMCM/build_gcc/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
