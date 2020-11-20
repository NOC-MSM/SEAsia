#!/bin/bash

:'

***********************
make_xios.sh
***********************

Checkout and compile the XIOS2.5 executable for I/O management
You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register
'
#::

  cd $WDIR
  #load modules
  module load cdt/15.11
  module unload PrgEnv-cray PrgEnv-gnu
  module load PrgEnv-intel/5.2.82
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.3.3.1
  module load cray-hdf5-parallel/1.8.14

  # download xios
  svn checkout -r 1566 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_DIR
  cd $XIOS_DIR

  # copy the arch
  cp $GITCLONE/XIOS/arch-XC30_ARCHER* arch/.

  # compile xios
  ./make_xios  --arch XC30_ARCHER --full --job 4
  #######################################
  cd $WORK
