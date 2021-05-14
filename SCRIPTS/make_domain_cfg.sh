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
  #cp $DOMAIN/hyb-z-s_DOMAINcfg_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg
  # Stretched-sigma vertical coordinates
  cp $DOMAIN/s-sig_DOMAINcfg_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg
  # z-partial-step vertical coordinates
  #cp $DOMAIN/z-ps_DOMAINcfg_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg

  # Ensure that the namelist_cfg has the appropriate parameters and number of lat,lon,depth levels set

  # Ensure the coordinates and bathymetry files, previously generated, are in place.
  ln -s $DOMAIN/coordinates.nc $TDIR/DOMAINcfg/.
  ln -s $DOMAIN/bathy_meter.nc $TDIR/DOMAINcfg/.

  # Make an adjustment to the source code for the hybrid coords
  #cp $GITCLONE/DOMAIN/domzgr_jelt_changes.f90 $TDIR/DOMAINcfg/src/domzgr.f90

  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?" $DOMAIN/job_create_template.slurm > $TDIR/DOMAINcfg/job_create.slurm
  
  # Submit the domain creation as a job,
  cd $TDIR/DOMAINcfg
  sbatch job_create.slurm

  # Rebuild the files. Here there are 8 tiles (and rebuilding on a single thread) 
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 1 domain_cfg 8

  # After create copy it and store it for further use
  cp $TDIR/DOMAINcfg/domain_cfg.nc $DOMAIN/domain_cfg_SEVERN.nc

  cd $WORK
