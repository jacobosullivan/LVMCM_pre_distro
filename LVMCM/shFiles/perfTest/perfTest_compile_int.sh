#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l node_type=nxv # compile on nxv
#$ -l h_rt=00:10:00
#$ -l h_vmem=1G

module load use.own
module load intel intelmpi dependencies cmake

# clean/navigate to build folder
cd /data/home/btx188/LVMCM_src_mm/LVMCM/

# select intel CMake project
cp CMakeLists_int.txt CMakeLists.txt

rm -rf build_int
mkdir build_int; cd build_int

# generate make files
cmake .. \
-DARMADILLO_LIBRARY=$HOME/install/armadillo/lib64/libarmadillo.so \
-DSUNDIALS_NVECSERIAL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_nvecserial.so \
-DSUNDIALS_KINSOL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_kinsol.so \
-DSUNDIALS_CVODE_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_cvode.so \
-DCMAKE_BUILD_TYPE=Release

# compile
make VERBOSE=1

exit 0
