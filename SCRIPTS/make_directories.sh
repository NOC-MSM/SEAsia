#!/bin/bash

#:'
#
#***********************
#make_directories.sh
#***********************
#
#Create the expected directory structure.
#'
#::

  # Choose an appropriate directory for your XIOS installation
  if [ ! -d "$XIOS_DIR" ]; then
    mkdir $XIOS_DIR
  fi

  # Choose an appropriate directory for your EXPeriment location
  if [ ! -d "$EXP" ]; then
    mkdir $EXP
  fi

  # Choose an appropriate directory for your NEMO installation
  if [ ! -d "$NEMO" ]; then
    mkdir $NEMO
  fi
