#!/bin/bash

:'

*******************************
make_bathymetry_from_parent.sh
*******************************

In the following the bathmetry file is constructed from the gridded bathymetry
of the parent model. This dataset is available on the JASMIN compute service. E.g.
at ``/gws/nopw/j04/nemo_vol6/acc/eORCA12-N512-ay652/domain/eORCA12_bathymetry_v2.4.nc``

The process is similar to ``make_bathymetry_from_gebco.sh``, which creates the
bathymetry from a GEBCO bathymetric data. Though here preprocessing steps are
not required before the bathymetry is mapped to the child grid using the SCRIP
tools.
'
#::

  cd $DOMAIN

  #load modules
  module unload nco cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  # Obtain the bathymetry from the ORCA12 global model from JASMIN.
  # Or on ARCHER you can take it from an existing directory
  cp /work/n01/n01/annkat/EXTRA_TOOLS/BATH/eORCA12_bathymetry_v2.4.nc $DOMAIN

  # Copy namelist for reshaping the parent data.
  # If necessary, edit namelist to point to correct input file and variable names
  cp $GITCLONE/DOMAIN/namelist_reshape_bilin_eORCA12 $DOMAIN

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
