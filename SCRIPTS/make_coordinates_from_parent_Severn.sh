#!/bin/bash

#:'
#
#*******************************
#make_coordinates_from_parent.sh
#*******************************
#
# Make a netCDF coordinates file holds that holds the horizontal grid coordinate
# position and spacings information. This is used to construct the 3D version:
# the domain_cfg.nc file.



# NB There was a problem with the GRIDGEN create_coordinates.exe tool. So here we
#   use the AGRIF version. The problem was that the spacing variables were not being
#  scaled by the zoom factor.

# Assuming that the NESTING tools have been build in make_tools.sh



#  Link in a parent coordinates file (from somewhere)::

  export coordinates_parent_file=coordinates_AMM15.nc
  #export coordinates_parent_file=coordinates_ORca_R12.nc

    ln -s $DOWNLOADS/$coordinates_parent_file $TDIR/NESTING/.


  # Get the namelist file from the cloned repository
  cp $DOMAIN/namelist.input_SEVERN $TDIR/NESTING/namelist.input

  # The namelist.input file controls the create process. Edit the bounding
  # coordinates (imax, ..., jmax) and scale factor (rho, rhot) to suit the child
  # coordinates. rho=rhot=3 with increase the resolution by a factor of 3.




  # Execute tool::

    ./agrif_create_coordinates.exe

  # This creates a coordinate file::
  # 1_coordinates_AMM15.nc


  # Copy it to the $INPUTS directory::

    cp 1_$coordinates_parent_file $INPUTS/coordinates.nc
