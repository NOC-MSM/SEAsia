#!/bin/bash
#SBATCH --job-name=nemo_dom
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --account=n01-ACCORD
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH hetjob

#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=120

# Setup the job environment (this module needs to be loaded before any other modules)
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
#export UCX_IB_REG_METHODS=direct

# Launch the parallel job
#   Using 512 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options
srun --distribution=block:block --hint=nomultithread --het-group=0 XXX_TDIR_XXX/DOMAINcfg/BLD/bin/make_domain_cfg.exe


### NOTE - THESE COMMANDS WORK ON THE COMMANDLINE BUT NOT IN THIS SCRIPT, BUT THEY NEED TO HAPPEN AFTER THE make_domain_cfg.exe IS COMPLETE.

  # Rebuild the files. Here there are 8 tiles (and rebuilding on a single thread)
srun --distribution=block:block --hint=nomultithread --het-group=0  XXX_TDIR_XXX/REBUILD_NEMO/rebuild_nemo -t 1 domain_cfg 8

  # After create copy it and store it for further use
srun --distribution=block:block --hint=nomultithread --het-group=0  cp XXX_TDIR_XXX/DOMAINcfg/domain_cfg.nc XXX_DOMAIN_XXX/domain_cfg_SEAsia.nc

