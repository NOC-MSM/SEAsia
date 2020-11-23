#!/bin/bash

:'

*******************************
make_coordinates_from_parent.sh
*******************************

Make a netCDF coordinates file holds that holds the horizontal grid coordinate
position and spacings information. This is used to construct the 3D version:
the domain_cfg.nc file.

This script makes a coordinates.nc file using a coordinate.nc file for a parent
configuration.

'
#::

  cd $TDIR/NESTING

  #load modules
  module unload nco cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  # subdomain of ORCA global
  # you can download the ORCA R12 coordinates
  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/coordinates.nc -O $TDIR/NESTING/coordinates_ORCA_R12.nc

  #  or in ARCHER you can take it from my directory
  # cp /work/n01/n01/jelt/EXTRA_TOOLS/GRIDS/coordinates_ORCA_R12.nc $TDIR/NESTING/.

  # Get the namelist file from the cloned NEMO-RELOC repository
  cp $GITCLONE/DOMAIN/namelist.input $TDIR/NESTING/

:'

The namelist.input file controls the create process. Edit the bounding
coordinates (imax, ..., jmax) and scale factor (rho, rhot) to suit the child
coordinates. rho=rhot=7 with increase the resolution by a factor of 7.
'

#::

  cd$TDIR/NESTING/

  # Generate the new coordinates file with the namelist.input settings
  ./agrif_create_coordinates.exe

  # copy it to your DOMAIN folder where it will be used to create the domin_cfg.nc file
  cp 1_coordinates_ORCA_R12.nc $DOMAIN/coordinates.nc

  cd $WORK
