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

  Write a namelist file::

    vi namelist.input

    &input_output
        iom_activated = true
    /
    &coarse_grid_files
        parent_coordinate_file = 'coordinates_AMM15.nc'
    /
    &bathymetry
    /
    &nesting
        imin = 694
        imax = 807
        jmin = 400
        jmax = 490
        rho  = 3
        rhot = 3
        bathy_update = false
    /
    &vertical_grid
    /
    &partial_cells
    /
    &nemo_coarse_grid
    /
    &forcing_files
    /
    &interp
    /
    &restart
    /
    &restart_trc
    /

  ---

  Execute tool::

    ./agrif_create_coordinates.exe

  This creates a coordinate file::

   1_coordinates_AMM15.nc



  Copy it to the $INPUTS directory::

    cp 1_coordinates_AMM15.nc $INPUTS/coordinates.nc
