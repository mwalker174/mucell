
#PBS -S /bin/bash
#PBS -j oe
#PBS -o output.txt
#PBS -m ae
#PBS -M mwalke49@jhmi.edu
#PBS -q rlwq
#PBS -l nodes=1:ppn=8
#PBS -l walltime=6:00:00
#PBS -N test

#Your model directory here
BASE_DIR=/home/mwalker/models/Stoch3D
BUILD_DIR=$BASE_DIR
EXEC=$BASE_DIR/stoch3d

#Place params.txt in WORK_DIR before running; output will also be saved to WORK_DIR
WORK_DIR=$BASE_DIR/output
#If you want to save 3D VTK files, create a "stochout" subdirectory in WORK_DIR
mkdir $WORK_DIR/stochout
cd $WORK_DIR

# load the modules
. /etc/profile.modules
module load openmpi/openmpi-1.6.5-intel
echo Modules loaded

#Do not recompile if using Intel because it is slow (see Notes.txt)
#cd $BUILD_DIR
#make clean
#make

echo pbs nodefile:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
cat $PBS_NODEFILE > host.list

#Try this line if you get errors about limited virtual memory
#ulimit -l unlimited

echo Running job...
export OMP_NUM_THREADS=$NPROCS
mpirun --mca btl openib,self -np 1 -hostfile host.list $EXEC

#Try ^openib flag if you get errors about openfabrics:
#mpirun --mca btl ^openib -x OMP_NUM_THREADS=8 -np 1 -hostfile host.list $EXEC
