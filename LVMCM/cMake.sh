#!/bin/bash

rm -rf build
mkdir build
cd build
#if module list 2>&1 >/dev/null | grep -q 'singularity'; then
#echo container modules loaded
#else
#echo loading container modules...
#module load singularity
#fi
#singularity exec /data/home/btx188/containers/LVMCM_mpi.img cmake ..
#singularity exec /data/home/btx188/containers/LVMCM_mpi.img make
#cd ..
cmake ..
make
