#!/bin/bash

:'

*******************************
make_coordinates_from_parent.sh
*******************************

Make a netCDF coordinates file holds that holds the horizontal grid coordinate
position and spacings information. This is used to construct the 3D version:
the domain_cfg.nc file.

This script generates a coordinates.nc file using a parent NEMO coordinates file.

'
#::

  cd $TDIR/NESTING

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env

  # subdomain of ORCA global
  # you can download the ORCA R12 coordinates
  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/coordinates.nc -O $TDIR/NESTING/coordinates_ORCA_R12.nc

  # Get the namelist file from the cloned NEMO-RELOC repository
  cp $DOMAIN/namelist.input $TDIR/NESTING/

:'

The namelist.input file controls the create process. Edit the bounding
coordinates (imax, ..., jmax) and scale factor (rho, rhot) to suit the child
coordinates. rho=rhot=7 with increase the resolution by a factor of 7.
'

#::

  # Generate the new coordinates file with the namelist.input settings
  ./agrif_create_coordinates.exe

  # copy it to your DOMAIN folder where it will be used to create the domin_cfg.nc file
  cp 1_coordinates_ORCA_R12.nc $DOMAIN/coordinates.nc

  cd $WORK
