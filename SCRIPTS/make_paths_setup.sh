#!/bin/bash

#:'
#
#***********************
#make_paths_setup.sh
#***********************
#'

# This script should be editted and sets all the directory paths that subsequently
# used.

#::

# Paths required for initial setup
export CONFIG=SEAsia
export WORK=/work/n01/n01/$USER
export WDIR=$WORK/$CONFIG
export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
export XIOS_DIR=$WORK/xios-2.5
export XIOS1_DIR=$WORK/xios-1.0
export FABM=$WDIR/NEMO_fabm
export GITCLONE=$WORK/NEMO-RELOC # Location of the git clone repo
