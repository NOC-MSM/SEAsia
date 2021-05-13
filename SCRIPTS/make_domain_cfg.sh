#!/bin/bash

:'

*******************************
make_domain_cfg.sh
*******************************

This script is generates s-sigma vertical coordinates with
the provided coordinates and bathymetry netCDF files.
'
#::

  ## Obtain the appropriate namelist (modify it if necessary)
  # Hybrid z-sigma vertical coordinates
  #cp $DOMAIN/hyb-z-s_namelist_cfg $DOMAIN/namelist_cfg
  # Stretched-sigma vertical coordinates
  cp $DOMAIN/s-sig_namelist_cfg $DOMAIN/namelist_cfg
  # z-partial-step vertical coordinates
  #cp $DOMAIN/z=ps_namelist_cfg $DOMAIN/namelist_cfg

  # Ensure the coordinates and bathymetry files, previously generated, are in place.
  ls $DOMAIN/coordinates.nc 
  ls $DOMAIN/bathy_meter.nc

  # Make an adjustment to the source code for the hybrid coords
  #cp $GITCLONE/DOMAIN/domzgr_jelt_changes.f90 $TDIR/DOMAINcfg/src/domzgr.f90

  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?" job_create_template.slurm > job_create.slurm
  # Submit the domain creation as a job,
  sbatch job_create.slurm

  # After create copy it and store it for durther use
  cp $DOMAIN/domain_cfg.nc $DOMAIN/domain_cfg_SEVERN.nc

  cd $WORK
