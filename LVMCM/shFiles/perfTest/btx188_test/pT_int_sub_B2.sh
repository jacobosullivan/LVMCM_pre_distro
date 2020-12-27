#!/bin/sh
#$ -cwd
#$ -pe smp 2 # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=0:10:00
#$ -l h_vmem=1G # memory request
#$ -t 31

# parameterization, edited by LVMCM_SUB_A.sh
PARFILE="/data/home/btx188/LVMCM_src_mm/LVMCM/parFiles/perfTest/perfTest_pars_int.csv"

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing performance test"
export I_MPI_HYDRA_TOPOLIB=ipl
mpirun -np $NSLOTS /data/home/btx188/LVMCM_src_mm/LVMCM/build_intel/LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

exit 0
