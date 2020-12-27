#!/bin/sh
#$ -cwd
#$ -pe smp CPU # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=0:10:00
#$ -l h_vmem=2G # memory request
#$ -t START-END

# parameterization, edited by LVMCM_SUB_A.sh
PARFILE=FILE

# load required modules
module load singularity openmpi

echo "Initializing performance test"
SINGULARITYENV_OMP_NUM_THREADS=1 mpirun -np $NSLOTS singularity exec /data/home/btx188/containers/LVMCM_mpi.img \
  /data/home/btx188/LVMCM_src_mm/LVMCM/build_gcc/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
