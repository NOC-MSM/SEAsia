#!/bin/bash

:'

***********************
make_xios.sh
***********************

Checkout and compile the XIOS2.5 executable for I/O management
You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register

Check the path in variable ``%XIOS_HOME`` in ``XC_ARCHER_INTEL_XIOS1`` and
``arch-XC30_ARCHER`` are consistent with your settings

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

#::

:'

At the time of writing XIOS1 was needed to build the DOMAINcfg, REBUILD_NEMO and
WEIGHTS tools. This is done as follows::
'

#::

  cd $WDIR

  # download xios
  svn checkout -r 703 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0 $XIOS1_DIR
  cd $XIOS1_DIR

  # copy the arch
  cp $GITCLONE/XIOS/arch-XC_ARCHER_INTEL_XIOS1* arch/.

  # compile xios
  ./make_xios  --arch XC_ARCHER_INTEL_XIOS1 --full --job 4
  #######################################
  cd $WORK
