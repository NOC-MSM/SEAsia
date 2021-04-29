#!/bin/bash

#:'
#
#***********************
#make_paths.sh
#***********************
#'

# This script should be editted and sets all the directory paths that subsequently
# used to compile NEMO/OPA and XIOS

# Assumed:
# git clone https://github.com/NOC-MSM/NEMO-RELOC.git NEMO-REGION


#::

# Paths required to compile (and run) NEMO/OPA and XIOS
export CONFIG=TESTCOMPILE
export WORK=/work/n01/n01/$USER
export WDIR=$WORK/NEMO-REGION # Also the git clone directory
export NEMO_VER=4.0.6
export NEMO=$WDIR/BUILD_EXE/NEMO/$NEMO_VER
#export CDIR=$NEMO/cfgs
export EXP=$WDIR/RUN_DIRECTORIES/EXP_TEST
export XIOS_DIR=$WDIR/BUILD_EXE/XIOS/xios-2.5
