#!/bin/sh
#$ -cwd
#$ -pe smp 1 # select number of CPUs
#$ -l node_type=nxv|sdx|sdv
#$ -l h_rt=1:00:00
#$ -l h_vmem=8G # memory request

echo "Analysing scaling"

module load R

Rscript /data/home/btx188/LVMCM_src/LVMCM/RCode/scaling_C_ij.R "/data/home/btx188/LVMCM_src/LVMCM/shFiles/scaling/bFileScaling.txt" "/data/SBCS-RossbergLab/jack/post_doc/scaling/scaling_quad_2-14.csv"

exit 0
