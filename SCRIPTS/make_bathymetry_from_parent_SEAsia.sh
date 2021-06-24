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

  # Obtain the bathymetry from the ORCA12 global model from JASMIN.
  # Or on ARCHER you can take it from an existing directory E.g.
  cp /work/n01/n01/annkat/EXTRA_TOOLS/BATH/eORCA12_bathymetry_v2.4.nc $DOMAIN/.
  /gws/nopw/j04/nemo_vol6/acc/eORCA12-N512-ay652/domain/eORCA12_bathymetry_v2.4.nc 
  #wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/bathymetry_ORCA12_V3.3.nc -O $DOMAIN/bathymetry_ORCA12_V3.3.nc

  # Execute first SCRIP process::
  ## Generate the new coordinates file with the namelist.input settings using AGRIF tool
  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?g" $DOMAIN/job_create_SCRIP_template.slurm > $DOMAIN/job_create_SCRIP.slurm

  # Submit the coordinates creation as a job (uses namelist_reshape_bilin_ORCA12 for settings)
  cd $DOMAIN
  sbatch   job_create_SCRIP.slurm 

# This batch script uses the SCRIP tool to perform a number of steps:
# Firstly, $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_eORCA12

#:'
#
#  Output files::
#
#    remap_nemo_grid_gebco.nc
#    remap_data_grid_gebco.nc
#'
#::

# Then, execute 2nd SCRIP process:
# $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_eORCA12

#:'
#
#  Output files::
#
#    data_nemo_bilin_gebco.nc
#'
#::

# Finally, execute third SCRIP process:
# $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_eORCA12

#:'

#  Output files::
#
#    bathy_meter.nc
#
# Finally, load ``nco`` tools to tidy some grid scale issues
#'
#::

  # load nco modules. Modules do not currently exist on ARCHER2 so do elsewhere.
  module load nco

  # Remove weirdness with negative bathymetry and make minimum bathymetry
  #equal to 10 m (resolve any possible wet-drying problems)
  ncap2 -s 'where(Bathymetry < 0) Bathymetry=0' bathy_meter.nc tmp1.nc
  ncap2 -s 'where(Bathymetry < 10 && Bathymetry > 0) Bathymetry=10' tmp1.nc -O bathy_meter.nc
  rm tmp1.nc

  # Copy it if you want for safe keeping
  cp bathy_meter.nc bathy_meter_ORCA12.nc

  # Fix bathymetry to deal with instabilities (opening some straights that
  #have only 2 grid points)
  ncap2 -s 'Bathymetry(140,464)=200' bathy_meter_ORCA12.nc bathy_meter_ORCA12.nc -O
  ncap2 -s 'Bathymetry(141,464)=200' bathy_meter_ORCA12.nc bathy_meter_ORCA12.nc -O
  ncap2 -s 'Bathymetry(145,563)=400' bathy_meter_ORCA12.nc bathy_meter_ORCA12.nc -O
  ncap2 -s 'Bathymetry(145,564)=400' bathy_meter_ORCA12.nc bathy_meter_ORCA12.nc -O
  ncap2 -s 'Bathymetry(140,467)=80' bathy_meter_ORCA12.nc bathy_meter_ORCA12.nc -O

  cd $WORK
