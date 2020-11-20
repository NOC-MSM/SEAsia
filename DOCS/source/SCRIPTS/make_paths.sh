#!/bin/bash

:'
***********************
make_paths.sh
***********************
'

# This script should be editted and sets all the directory paths that subsequently
# used.

#::

  # [quote-snippet]
  export CONFIG=SEAsia_test
  export WORK=/work/n01/n01/$USER
  # [quote-snippet]
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
  export GITCLONE=$WDIR/FILES_START
  export RIVER=$WDIR/RIVERS
  export FABM=$WDIR/NEMO_fabm
