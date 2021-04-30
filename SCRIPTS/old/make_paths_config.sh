#!/bin/bash

:'

***********************
make_paths_config.sh
***********************
'

# This script should be editted and sets all the directory paths that subsequently
# used.

#::

  # Paths required for initial setup
  export CONFIG=NEMO_FABM_ERSEM
  export WORK=/work/n01/n01/$USER
  export WDIR=$WORK/$CONFIG
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export XIOS_DIR=$WORK/xios-2.5
  export XIOS1_DIR=$WORK/xios-1.0
  export FABM=$WDIR/NEMO_fabm
  export GITCLONE=$WDIR/BUILD_CFG

  # Paths also required for configuration build
  export DOMAIN=$WDIR/DOMAIN
  export SBC=$WDIR/SURFACE_FORCING
  export TIDE=$WDIR/TIDAL_FORCING
  export OBC=$WDIR/BOUNDARY_FORCING
  export ICS=$WDIR/INITIAL_CONDITIONS
  export RUND=$WDIR/RUN_DIRECTORY
  export EXP=$RUND/EXP_TEST
  export SCRIPTS=$WDIR/SCRIPTS
  export RIVER=$WDIR/RIVERS
