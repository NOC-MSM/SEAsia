#!/bin/bash

:'

***********************
main_setup.sh
***********************

Execute this script for end-to-end setup (assuming all the paths and modules
are correct).
The first part complete generic setup steps. The second part completes
configuration specific steps.
'
#::

  # Generic setup steps
  echo "Making Paths"
  . ./make_paths_setup.sh                         > main_output_1.txt 2>&1
  echo "Making Directories"
  . ./make_directories_setup.sh                  >> main_output_1.txt 2>&1
  echo "Installing XIOS_2.5 - this will take 5-10 mins"
  . ./make_xios.sh                               >> main_output_1.txt 2>&1
  echo "Compiling various grid tools"
  . ./make_tools.sh                              >> main_output_1.txt 2>&1
  #echo "Installing NEMO-FABM-ERSEM - this will take a good 10/15 mins"
  #echo "WARNING - this automatically chooses OPA_SRC and TOP_SRC"
  #. ./make_nemo_fabm_ersem.sh                      >> main_output.txt 2>&1
  echo "Installing NEMO - this will take a good 10/15 mins!"
  echo "WARNING - this automatically chooses OPA_SRC only"
  echo "If you want to choose anything else e.g. LIM_SRC, remove the printf in make_nemo.sh"
  . ./make_nemo.sh                               >> main_output_1.txt 2>&1

  # Configuration specific steps
  echo "Making Paths"
  . ./make_paths_config.sh                        > main_output_2.txt 2>&1
  echo "Making Directories"
  . ./make_directories_config.sh                 >> main_output_2.txt 2>&1
  echo "Creating coordinate file"
  . ./make_coordinates_from_parent.sh            >> main_output_2.txt 2>&1
  echo "Creating bathymetry for our domain"
  . ./make_bathymetry.sh                         >> main_output_2.txt 2>&1
  echo "Creating domain configuration file"
  . ./make_domain_cfg.sh                         >> main_output_2.txt 2>&1
  echo "Creating rivers"
  . ./create_rivers.sh                           >> main_output_2.txt 2>&1
  echo "Creating IC"
  . ./make_CMEMS_IC.sh                           >> main_output_2.txt 2>&1
  echo "Atmospheric forcing weights"
  . ./create_WA.sh                               >> main_output_2.txt 2>&1
  echo "Prepare experiment"
  . ./prepare_EXP.sh                             >> main_output_2.txt 2>&1
  #echo "DONE"
