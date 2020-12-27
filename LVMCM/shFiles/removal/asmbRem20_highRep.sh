#!/bin/sh
#$ -cwd
#$ -pe smp 4 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=3:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-100
#$ -t 20

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Initializing node removal assembly"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/parFiles/removal/removal_pars_20.txt

let "SMO=$SGE_TASK_ID-1"
let "PARLINE=1+($SMO-$SMO%10)/10"
echo $PARLINE
for ITER in {0..9}
  do
    let "REP=$ITER + 10*($SMO%10)"
    echo $REP
    /data/home/btx188/LVMCM_src/LVMCM/build/LVMCM $(sed -n "${PARLINE}p" $PARFILE) $REP
  done

exit 0
