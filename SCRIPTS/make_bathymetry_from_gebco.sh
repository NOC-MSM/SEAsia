#!/bin/bash

#:'
#
#*******************************
#make_bathymetry_from_gebco.sh
#*******************************
#
#In the following we use gridded data from GEBCO https://www.gebco.net/data_and_products/gridded_bathymetry_data/.
# Specifically the 2020 at 15 arc-second intervals.
#
# For a Severn demonstration: Bounds N 51.8 W -6.9 S 50.1 E -2
# Grid dimensions: W 1176 H 408
# File name: gebco_2020_n51.8_s50.1_w-6.9_e-2.0.nc
#
#The process is similar t0 ``make_bathymetry_from_parent.sh``, which creates the
#bathymetry from a parent configuration bathymetry file. Here additionally some data is
#pre-processing is required to 1) flatten out land elevations, 2) make depths
#positive. Then the bathymetry is mapped from the parent to the child grid using
#the SCRIP tools.

#'
#::

  export BATHYFILE=gebco_in.nc
  #export BATHYFILE=gebco_2020_n51.8_s50.1_w-6.9_e-2.0.nc
  #export BATHYFILE=GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc


  ## On a machine process the bathymetry. Modify the 'elevation' variable
  ## such that water depths as positive; land points are zero
  #
  # module load nco
  #
  ## Flatten the land
  # ncap2 -s 'where(elevation > 0) elevation=0' $BATHYFILE tmp.nc
  ## Make negative depths positive
  # ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  # rm tmp.nc


## Copy BATHYFILE into $DOMAIN
## E.g. copy the gebco_in.nc file to ARCHER2
## E.g. scp gebco_in.nc  jelt@archer2:/work/n01/n01/$USER/NEMO-REGION/BUILD_CFG/DOMAIN/.

#:'
#
#Then the bathymetry is interpolated onto the coordinates grid using SCRIP tools
#'
#::

  # Check namelist for reshaping GEBCO data
  # If necessary, edit namelist to point to correct input file and variable names
  cd $DOMAIN

  # Link in the new coordinates file (required but already there)
  # ln -s $TDIR/NESTING/1_coordinates_AMM15.nc coordinates.nc

  # Execute first SCRIP process::
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

#  Output files::
#
#    remap_nemo_grid_gebco.nc
#    remap_data_grid_gebco.nc

  #Execute second SCRIP process:
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

#  Output files::
#
#    data_nemo_bilin_gebco.nc


  # Execute third SCRIP process:
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

#  Output files::
#
#    bathy_meter.nc
#
#  Finally, load ``nco`` tools to tidy some grid scale issues

#
# Finally make minor modifications to the bathymetry file
#
#  # load nco modules (this was done locally)
#  module load nco
#
#  # Remove weirdness with negative bathymetry and make minimum bathymetry
#  #equal to 10 m (resolve any possible wet-drying problems)
#  ncap2 -s 'where(Bathymetry < 0) Bathymetry=0' bathy_meter.nc tmp1.nc
#  ncap2 -s 'where(Bathymetry < 10 && Bathymetry > 0) Bathymetry=10' tmp1.nc -O bathy_meter.nc
#  rm tmp1.nc
#
#  # Copy it if you want for safe keeping
#  cp bathy_meter.nc bathy_meter_10m.nc
#
#  # Fix bathymetry to deal with instabilities (opening some straights that
#  #have only 2 grid points)
#  # E.g. ncap2 -s 'Bathymetry(0,0)=0' bathy_meter_10m.nc bathy_meter_10m.nc -O



  cd $WORK
