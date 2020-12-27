#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=00:20:00
#$ -l h_vmem=0.5G
#$ -t 1-5 # select lines of parameter file

# copy and edit template shell script
## set parameter file
PARFILE="/data/home/aax010/git/LVMCM_src/LVMCM/parFiles/perfTest/perfTest_pars_aax010.csv"
## format for regular expression
PARFILE="$(echo $PARFILE | sed 's/\//\\\//g')"
PARFILE="$(echo $PARFILE | sed 's/\./\\\./g')"

## set simulation script template file
TEMPLATE="perfTest_sub_B_aax010.sh"

## extract job/task ID from environmental variables...
JB_TSK_ID="$JOB_ID.$SGE_TASK_ID"
## use to generate file name of similation script
SHFILE="_$JOB_ID.$SGE_TASK_ID.sh"

# make and navigate into output file directory
## name by parent job/task ID
mkdir $JB_TSK_ID
cd $JB_TSK_ID
# copy template
## name by parent job/task ID
cp ../$TEMPLATE $SHFILE
## move output/error files current job
mv ../*[a-z]$JB_TSK_ID .

# edit submission script
## set parameter file
sed -i "s/^PARFILE.*/PARFILE\=$PARFILE/" $SHFILE
## set line number
START=$((1 +  30 * ( ${SGE_TASK_ID} -1 )))
END=$((30 * ${SGE_TASK_ID}))
sed -i "s/START/$START/" $SHFILE
sed -i "s/END/$END/" $SHFILE
# set CPU number
CPU=$((2 ** ( ${SGE_TASK_ID} -1)))
sed -i "s/CPU/$CPU/" $SHFILE

# submit concurrent simulation job
echo "Submitting Performance test array job"
qsub -notify $SHFILE # -notify flag required for signal handler

exit 0
