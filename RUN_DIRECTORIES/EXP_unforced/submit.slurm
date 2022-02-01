#!/bin/bash
#SBATCH --job-name=nemo_test
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=n01-ACCORD
#SBATCH --partition=standard
#SBATCH --qos=standard
##BATCH --reservation=shortqos
##BATCH --qos=short
# Created by: mkslurm -S 1 -s 16 -m 1 -C  2 -c  2 -t 00:10:00 -a n01-ACCORD -j nemo_test

# Setup the job environment (this module needs to be loaded before any other modules)
#module -s restore /work/n01/shared/acc/n01_modules/ucx_env
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

export OMP_NUM_THREADS=1
#
cat > myscript_wrapper2.sh << EOFB
#!/bin/ksh
#
set -A map /work/n01/n01/micdom/SEVERN-SWOT/BUILD_EXE/XIOS/xios-2.5/bin/xios_server.exe /work/n01/n01/micdom/SEVERN-SWOT/BUILD_EXE/NEMO/4.0.6/cfgs/NEMOconstTS/BLD/bin/nemo.exe
exec_map=( 0 1 1 )
#
exec \${map[\${exec_map[\$SLURM_PROCID]}]} 
##
EOFB
chmod u+x ./myscript_wrapper2.sh
#
srun --mem-bind=local --cpu-bind=v,map_cpu:00,0x10,0x12, ./myscript_wrapper2.sh
