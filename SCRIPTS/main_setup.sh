#!/bin/bash

#***********************
#main_setup.sh
#***********************
#
# This script builds the exectuables needed when setting up a new configuration.

  # Generic setup steps
  echo "Making Paths"
  . ./SCRIPTS/make_paths.sh                         > main_output_1.txt 2>&1
  echo "Making Directories"
  . ./SCRIPTS/make_directories.sh                  >> main_output_1.txt 2>&1
  echo "Installing XIOS_2.5 - this will take 5-10 mins"
  . ./SCRIPTS/make_xios.sh                               >> main_output_1.txt 2>&1
  echo "Installing NEMO - this will take a good 10/15 mins!"
  echo "WARNING - this automatically chooses OPA_SRC only"
  . ./SCRIPTS/make_nemo.sh                         >> main_output_1.txt 2>&1
  echo "Compiling various grid tools"
  . ./make_tools.sh                                >> main_output_1.txt 2>&1
  
