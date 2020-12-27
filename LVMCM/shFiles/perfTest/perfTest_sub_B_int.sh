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
module load use.own
module load intel intelmpi dependencies

echo "Initializing performance test"
export I_MPI_HYDRA_TOPOLIB=ipl
mpirun -np $NSLOTS /data/home/btx188/LVMCM_src_mm/LVMCM/build_int/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
