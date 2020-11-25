#!/bin/bash

:'

***********************
make_tools.sh
***********************

At the time of writing the tools used here required an older version of XIOS to
work (XIOS1). This is compiled in make_xios.sh. Refer to the README in the
 $TDIR/DOMAINcfg for more details.

Check the path in variable ``%XIOS_HOME`` in ``XC_ARCHER_INTEL_XIOS1`` are
consistent with your settings

Also make the SOSIE tool for flood filling land in parent data when interpolating
 onto child grid
'
#::

  # Start in the TOOLS directory
  cd $TDIR
  # Copy the compiler option files from the NEMO-RELOC repository clone
  #  into the architecture directory
  cp $GITCLONE/NEMO-FABM-ERSM/arch-XC_ARCHER_INTEL_NOXIOS.fcm ../ARCH/.
  cp $GITCLONE/NEMO-FABM-ERSM//arch-XC_ARCHER_INTEL_XIOS1.fcm  ../ARCH/.

  # Apply patches for the weight file code
  cd $TDIR/WEIGHTS/src
  patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripinterp_mod.patch
  patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripinterp.patch
  patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scrip.patch
  patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripshape.patch
  patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripgrid.patch

  # Load modules
  module unload nco cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  # compile tools
  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_NOXIOS -n NESTING -j 6
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n WEIGHTS


  # Make SOSIE tool.
  $SCRIPTS/install_sosie.sh

  cd $WORK
