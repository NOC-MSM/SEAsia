#!/bin/bash
#PBS -N rebuild
#PBS -l select=serial=true:ncpus=1
#PBS -l walltime=00:20:00
#PBS -A n01-ACCORD
#PBS -q short

## LOAD MODULES
module load cray-libsci/16.11.1
module load cray-mpich/7.5.5
module load cray-netcdf-hdf5parallel/4.4.1.1
module load cray-hdf5-parallel/1.10.0.1
module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf-hdf5parallel cray-hdf5-parallel

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
ulimit -c unlimited
ulimit -s unlimited

cd $PBS_O_WORKDIR
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 25hourm_grid_V 12
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 25hourm_grid_U 12
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 25hourm_grid_T 12
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 shelftmb_grid_T 12
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 output.abort 1716
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 mesh_mask 168
#$TDIR/REBUILD_NEMO/rebuild_nemo -t 1 mesh_mask 1716  
/work/n01/n01/nibrun/NEMO/AMM15_v3_6_STABLE_package/NEMOGCM/TOOLS/REBUILD_NEMO/rebuild_nemo -t 1 output.init 88
