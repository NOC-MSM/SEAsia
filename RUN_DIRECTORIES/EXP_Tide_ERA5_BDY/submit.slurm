#!/bin/bash
#SBATCH --job-name=nemo_test
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --account=n01-ACCORD
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the job environment (this module needs to be loaded before any other modules)
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
#   Using 512 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options
srun --het-group=0 --nodes=1 --tasks-per-node=2 ${SHARED_ARGS} /work/n01/n01/$USER/SEVERN-SWOT/BUILD_EXE/XIOS/xios-2.5/bin/xios_server.exe : \
     --het-group=1 --nodes=3 --tasks-per-node=126 ${SHARED_ARGS} /work/n01/n01/$USER/SEVERN-SWOT/BUILD_EXE/NEMO/4.0.6/cfgs/NEMOconstTS/BLD/bin/nemo.exe

