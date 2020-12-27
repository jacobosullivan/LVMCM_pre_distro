#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l node_type=nxv # compile on nxv
#$ -l h_rt=00:10:00
#$ -l h_vmem=1G

module load singularity openmpi

# clean/navigate to build folder
cd /data/home/btx188/LVMCM_src_mm/LVMCM/

# select intel CMake project
cp CMakeLists_gcc.txt CMakeLists.txt

rm -rf build_gcc
mkdir build_gcc; cd build_gcc

# generate make files
singularity exec /data/home/btx188/containers/LVMCM_mpi.img cmake ..

# compile
singularity exec /data/home/btx188/containers/LVMCM_mpi.img make VERBOSE=1

exit 0
