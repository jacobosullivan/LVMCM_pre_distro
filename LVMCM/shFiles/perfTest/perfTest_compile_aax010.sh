#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l node_type=nxv # compile on nxv
#$ -l h_rt=00:10:00
#$ -l h_vmem=1G

module load use.own gsl cmake btx188
module load intel/2020.0 intelmpi/2020.0

# clean/navigate to build folder
cd /data/home/aax010/git/LVMCM_src/LVMCM/
rm -rf build
mkdir build; cd build

# generate make files
cmake .. \
-DARMADILLO_LIBRARY=$HOME/usr/lib64/libarmadillo.so \
-DSUNDIALS_NVECSERIAL_LIB:FILEPATH=$HOME/usr/lib64/libsundials_nvecserial.so \
-DSUNDIALS_KINSOL_LIB:FILEPATH=$HOME/usr/lib64/libsundials_kinsol.so \
-DSUNDIALS_CVODE_LIB:FILEPATH=$HOME/usr/lib64/libsundials_cvode.so

# compile
make VERBOSE=1

exit 0
