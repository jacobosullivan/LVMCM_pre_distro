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
module load use.own gsl cmake btx188
module load intel/2020.0 intelmpi/2020.0

echo "Initializing performance test"
set OMP_NUM_THREADS=1 
mpirun -np $NSLOTS /data/home/aax010/git/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
