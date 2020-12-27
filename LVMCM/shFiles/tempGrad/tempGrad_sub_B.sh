#!/bin/sh
#$ -cwd
#$ -pe smp CPU # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=H_RT:00:00
#$ -l h_vmem=H_VMEMG # memory request
#$ -t START-END

# parameterization, edited by tempGrad_sub_A.sh
PARFILE=FILE

# load required modules
module load singularity openmpi

echo "Initializing temperature gradient assembly"

if [ "$NSLOTS" -gt 1 ]; then
  SINGULARITYENV_OMP_NUM_THREADS=1 mpirun -np $NSLOTS singularity exec /data/home/btx188/containers/LVMCM_mpi.img \
    /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)
else
  SINGULARITYENV_OMP_NUM_THREADS=1 singularity exec /data/home/btx188/containers/LVMCM_mpi.img \
    /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)
fi

exit 0
