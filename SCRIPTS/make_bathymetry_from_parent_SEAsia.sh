#!/bin/bash

#:'
#
#*************************************
#make_bathymetry_from_parent_SEAsia.sh
#*************************************
#
# In the following the bathmetry file is constructed from the gridded bathymetry
# of the parent model. This dataset is available on the JASMIN compute service. E.g.
# at ``/gws/nopw/j04/nemo_vol6/acc/eORCA12-N512-ay652/domain/eORCA12_bathymetry_v2.4.nc``
#
# The process is similar to ``make_bathymetry_from_gebco.sh``, which creates the
# bathymetry from a GEBCO bathymetric data. Though here preprocessing steps are
# not required before the bathymetry is mapped to the child grid using the SCRIP
# tools.
# '
#::

  cd $DOMAIN

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env

  # Obtain the bathymetry from the ORCA12 global model from JASMIN. The e(xtended)ORCA domain has a higher latitude reach. 
  # Or on ARCHER you can take it from an existing directory E.g. /gws/nopw/j04/nemo_vol6/acc/eORCA12-N512-ay652/domain/eORCA12_bathymetry_v2.4.nc 
  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/bathymetry_ORCA12_V3.3.nc -O $DOWNLOADS/bathymetry_ORCA12_V3.3.nc

  ln -s $DOWNLOADS/bathymetry_ORCA12_V3.3.nc $DOMAIN/bathymetry_ORCA12_V3.3.nc

  # Execute first SCRIP process::
  ## Generate the new coordinates file with the namelist.input settings using AGRIF tool
  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?g" $DOMAIN/job_create_SCRIP_template.slurm > $DOMAIN/job_create_SCRIP.slurm

  # Submit the coordinates creation as a job (uses namelist_reshape_bilin_ORCA12 for settings)
  cd $DOMAIN
  sbatch   job_create_SCRIP.slurm 

# This batch script uses the SCRIP tool to perform a number of steps:
# Firstly, $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_ORCA12
#  Output files::
#    remap_nemo_grid_gebco.nc
#    remap_data_grid_gebco.nc
#'
# Then, execute 2nd SCRIP process:
# $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_ORCA12
#  Output files::
#    data_nemo_bilin_gebco.nc

# Finally, execute third SCRIP process:
# $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_ORCA12
#  Output files::
#    bathy_meter.nc
#

cd $WDIR/SCRIPTS
