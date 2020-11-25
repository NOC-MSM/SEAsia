#!/bin/bash

:'

***********************
make_paths.sh
***********************
'

# This script should be editted and sets all the directory paths that subsequently
# used.

#::

  export CONFIG=SEAsia_test
  export WORK=/work/n01/n01/$USER
  export WDIR=$WORK/$CONFIG
  export NEMO=$WORK/$CONFIG/NEMO_ERSEM_CMEMS
  export DOMAIN=$WDIR/DOMAIN
  export SBC=$WDIR/SURFACE_FORCING
  export TIDE=$WDIR/TIDAL_FORCING
  export OBC=$WDIR/BOUNDARY_FORCING
  export ICS=$WDIR/INITIAL_CONDITIONS
  export CDIR=$NEMO/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$NEMO/trunk_NEMOGCM_r8395/TOOLS
  export RUND=$WDIR/RUN_DIRECTORY
  export EXP=$RUND/EXP_TEST
  export XIOS_DIR=$WDIR/xios-2.5
  export XIOS1_DIR=$WDIR/xios-1.0
  export GITCLONE=$WDIR/START_FILES
  export SCRIPTS=$WDIR/SCRIPTS
  export RIVER=$WDIR/RIVERS
  export FABM=$WDIR/NEMO_fabm
