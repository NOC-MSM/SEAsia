#!/bin/bash
#SBATCH --job-name=nemo_unforced
#SBATCH --time=00:20:00
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --account=n01-ACCORD
#SBATCH --partition=standard
#SBATCH --qos=standard

# Ensure the correct modules are loaded for ARCHER2
# Load modules listed in /work/n01/shared/nemo/setup
# Tested 11Jan22
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1
export UCX_IB_REG_METHODS=direct
ulimit -s unlimited
ulimit -c unlimited

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

# Launch the parallel job
#   srun picks up the distribution from the sbatch options
srun --het-group=0 --nodes=1 --tasks-per-node=2 ${SHARED_ARGS} ./xios_server.exe : \
     --het-group=1 --nodes=7 --tasks-per-node=120 ${SHARED_ARGS} ./nemo.exe


