#!/bin/bash
#SBATCH --job-name=nemo_dom
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --account=n01-ACCORD
#SBATCH --partition=serial
#SBATCH --qos=serial

# Define memory required for this jobs. By default, you would 
# get just under 2 GB, but you can ask for up to 125 GB.
#SBATCH --mem=16G

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

# Ensure the correct modules are loaded for ARCHER2
# Load modules listed in /work/n01/shared/nemo/setup
# Tested 10Jan22
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

# Launch the serial job
XXX_TDIR_XXX/DOMAINcfg/BLD/bin/make_domain_cfg.exe

# After create copy it and store it for further use
cp XXX_TDIR_XXX/DOMAINcfg/domain_cfg.nc XXX_DOMAIN_XXX/domain_cfg_XXX_REPO_XXX.nc
