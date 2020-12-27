#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G # memory request
#$ -t 1-10

# load required modules
module load use.own
module load intel intelmpi dependencies

echo "Running harvesting experiment"

PARFILE=/data/home/btx188/LVMCM_src/LVMCM/shFiles/scaling/sc16/bFileList16.txt

for i in {1..10}
do
  let "LINE=($SGE_TASK_ID-1)*10+$i"
  /data/home/btx188/LVMCM_src/LVMCM/build_int/LVMCM -b $(sed -n "${LINE}p" $PARFILE) -H T
done

exit 0
