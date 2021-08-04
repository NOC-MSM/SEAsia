#!/bin/bash

#**************************************
#make_coordinates_from_parent_SEAsia.sh
#**************************************
#
# Make a netCDF coordinates file holds that holds the horizontal grid coordinate
# position and spacings information. This is used to construct the 3D version:
# the domain_cfg.nc file.
#
# This script generates a coordinates.nc file using a parent NEMO coordinates file.
#
# To demonstrate the prinicple, here a straight forward cut-out is done from a parent grid.
# For more complex relocable domain, with for example resolution refinement, the AGRIF tools can be used.
# (outlined below).

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env
  module load nco

  # subdomain of ORCA global
  # you can download the ORCA R12 coordinates
  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/coordinates.nc -O $DOMAIN/coordinates_ORCA_R12.nc

  
  cd $DOMAIN
  
  ncks -d x,50,730 -d y,1250,1800 coordinates_ORCA_R12.nc -O cropped.nc 
  ncwa -a time cropped.nc coordinates_ORCA_R12_cropped.nc
  ln -s coordinates_ORCA_R12_cropped.nc coordinates.nc
  rm cropped.nc
  
  cd $WDIR/SCRIPTS
  

### AN OUTLINE OF AGRIF TOOL USEAGE.
#
#  cd $TDIR/NESTING
#
#  #load modules
#  module -s restore /work/n01/shared/acc/n01_modules/ucx_env
#
#  # subdomain of ORCA global
#  # you can download the ORCA R12 coordinates
#  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/coordinates.nc -O $TDIR/NESTING/coordinates_ORCA_R12.nc
#
#  # Get the namelist file from the cloned NEMO-RELOC repository
#  cp $DOMAIN/namelist.input_SEAsia $TDIR/NESTING/namelist.input
#
## The namelist.input file controls the create process. Edit the bounding
## coordinates (imax, ..., jmax) and scale factor (rho, rhot) to suit the child
## coordinates. rho=rhot=7 with increase the resolution by a factor of 7.
#
#  ## Generate the new coordinates file with the namelist.input settings using AGRIF tool
#  # Edit job script
#  sed "s?XXX_TDIR_XXX?$TDIR?" $DOMAIN/job_create_coords_template.slurm > $TDIR/NESTING/job_create_coords.slurm
#
#  # Submit the coordinates creation as a job (uses namelist.input for settings)
#  cd $TDIR/NESTING
#  sbatch job_create_coords.slurm
#
#  # copy it to your DOMAIN folder where it will be used to create the domin_cfg.nc file
#  cp 1_coordinates_ORCA_R12.nc $DOMAIN/coordinates.nc
#
#  cd $WDIR/SCRIPTS





