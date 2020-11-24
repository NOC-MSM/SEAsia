#!/bin/bash

:'

*******************************
make_domain_cfg.sh
*******************************

This script is generates hybrid z-sigma vertical coordinates with
the provided coordinates and bathymetry netCDF files.
'
#::


  cd $TDIR

  # Load modules
  module unload nco cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  ## Obtain the appropriate namelist (modify it if necessary)
  # Hybrid z-sigma vertical coordinates
  cp $GITCLONE/DOMAIN/hyb-z-s_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg
  # Stretched-sigma vertical coordinates
  #cp $GITCLONE/DOMAIN/s-sig_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg
  # z-partial-step vertical coordinates
  #cp $GITCLONE/DOMAIN/z=ps_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg

  # link coordinates and bathymetry files, previously generated
  ln -s $DOMAIN/coordinates.nc $TDIR/DOMAINcfg/.
  ln -s $DOMAIN/bathy_meter_ORCA12.nc $TDIR/DOMAINcfg/bathy_meter.nc

  # Make an adjustment to the source code for the hybrid coords
  cp $GITCLONE/DOMAIN/domzgr_jelt_changes.f90 $TDIR/DOMAINcfg/src/domzgr.f90

  cd $TDIR

  # Recompile DOMAINcfg tool with code fix
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg clean
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg

  cd $TDIR/DOMAINcfg

  # Submit the domain creation as a job, so copy the jobs script and
  # modify email notification is desired.
  cp $GITCLONE/DOMAIN/job_create.sh $TDIR/DOMAINcfg/job_create.sh
  qsub -q short job_create.sh

  # After create copy it and store it for durther use
  cp $TDIR/DOMAINcfg/domain_cfg.nc $DOMAIN/domain_cfg_ORCA12.nc

  cd $WORK
