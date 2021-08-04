#!/bin/bash

#:'
#
#*******************************
#make_domain_cfg.sh
#*******************************
#
# This script is generates s-sigma vertical coordinates with
# the provided coordinates and bathymetry netCDF files.
#'
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

  ## Make an adjustment to the DOMAINcfg source code to accomodate more varied vertical coords.
  ## Done in make_tools.sh
  #cp $DOMAIN/domzgr.f90.melange $TDIR/DOMAINcfg/src/domzgr.f90

  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?g" $DOMAIN/job_create_domain_template.slurm > $TDIR/DOMAINcfg/job_create_domain.slurm
  sed -i "s?XXX_DOMAIN_XXX?$DOMAIN?g" $TDIR/DOMAINcfg/job_create_domain.slurm
  
  # Submit job script to build domain_cfg tiles
  cd $TDIR/DOMAINcfg
  sbatch job_create_domain.slurm
  
  #wait for domain creation job to finish
  for i in {0..7}; do #8 tiles
  while [ ! -f domain_cfg_000$i.nc ] ;
  do
      echo  "wait for domain creation job to finish"
      sleep 60
  done
  done
  
  # Rebuild the files. Here there are 8 tiles (and rebuilding on a single thread) 
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 1 domain_cfg 8

  # After create copy it and store it for further use
  cp $TDIR/DOMAINcfg/domain_cfg.nc $DOMAIN/domain_cfg_SEAsia.nc
  rm $TDIR/DOMAINcfg/domain_cfg_000*.nc #remove tiles
  cd $WDIR/SCRIPTS