#!/bin/bash

:'

*******************************
make_bathymetry_from_gebco.sh
*******************************

In the following I use the 1-minute 2008 data. The download was called:
``GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc``

The process is similar t0 ``make_bathymetry_from_parent.sh``, which creates the
bathymetry from a parent configuration bathymetry file. Here additionally some data is
pre-processing is required to 1) flatten out land elevations, 2) make depths
positive. Then the bathymetry is mapped from the parent to the child grid using
the SCRIP tools.

'
#::

  cd $DOMAIN

  # On ARCHER modules had to be unloaded for NCO tools to work
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  # Flatten the land
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc tmp.nc
  # Make negative depths positive
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

  # Restore the original parallel modules::
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

:'

Then the bathymetry is interpolated onto the coordinates grid using SCRIP tools
'
#::

  # Copy namelist for reshaping GEBCO data
  # If necessary, edit namelist to point to correct input file and variable names
  cp $GITCLONE/DOMAIN/namelist_reshape_bilin_gebco $DOMAIN/.

  # Execute first SCRIP process::
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

:'

  Output files::

    remap_nemo_grid_gebco.nc
    remap_data_grid_gebco.nc
'
#::

  #Execute second SCRIP process:
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

:'

  Output files::

    data_nemo_bilin_gebco.nc
'
#::

  # Execute third SCRIP process:
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

:'

  Output files::

    bathy_meter.nc

  Finally, load ``nco`` tools to tidy some grid scale issues
'
#::

  # load nco modules
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

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
