#!/bin/bash

:'

************************
make_nemo_fabm_ersem.sh
************************

Checkout and compile the NEMO executable with FABM ERSEM for physics with
biogeochemistry

You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register

Bare in mind that NEMO is being compiled for use with a particular version of
XIOS. Ensure edits to ``%XIOS_HOME`` in the ``arch-*`` file are consistent with
the XIOS build.

This revision is referred to as NEMOvp4 (proceeds v4), as it contains some
significant structural changes that were introduced in version 4, though much
of the code is still at v3.6.
'
#::

  cd $WDIR
  #################################################################
  # First get/download NEMO and FABM ERSEM
  #################################################################
  # checkout NEMO from the paris repository
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

  # Replace TOP_SRC with a cloned version from the NEMO-RELOC repository clone
  mv trunk_NEMOGCM_r8395/NEMO/TOP_SRC trunk_NEMOGCM_r8395/NEMO/TOP_SRC_old
  # The NEMO-FABM coupler ia maintained here: https://github.com/NOC-MSM/NEMO_ERSEM/tree/master/TOP_SRC_r8395_FABM
  cp -r $GITCLONE/BUILD_CFG/TOP_SRC_r8395_FABM trunk_NEMOGCM_r8395/NEMO/TOP_SRC

  cd $WDIR
  # get ERSEM from the NEMO-RELOC repository clone
  cp -r $GITCLONE/BUILD_CFG/ERSEM-master ./

  cd $FABM
  # clone FABM from its GitHub repository
  git clone https://github.com/fabm-model/fabm.git

  #######################################################################
  # Compile FABM
  #######################################################################
  # load modules
  module load cdt/15.11
  module unload PrgEnv-cray PrgEnv-gnu
  module load PrgEnv-intel/5.2.82
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.3.3.1
  module load cray-hdf5-parallel/1.8.14

  module load cmake

  cd $FABM

  cmake $FABM/fabm/src -DFABM_HOST=nemo -DCMAKE_Fortran_COMPILER=ifort \
        -DFABM_ERSEM_BASE=$WDIR/ERSEM-master -DFABM_EMBED_VERSION=ON

  make install

  # Move the compiled code to the $FABM directory
  cp -r /home/n01/n01/$USER/local/fabm/nemo/lib ./
  cp -r /home/n01/n01/$USER/local/fabm/nemo/include ./

  #################################################################
  # Compile nemo and update TOP source
  #################################################################
  # get arch
  # ATTENTION modify the following file to have the correct paths
  cp $WDIR/trunk_NEMOGCM_r8395/NEMO/TOP_SRC/arch-XC_ARCHER_INTEL_FABM.fcm $WDIR/trunk_NEMOGCM_r8395/ARCH/arch-XC_ARCHER_INTEL_FABM.fcm

  cd $CDIR
  # make your configuration (**ATTENTION**: here we choose to include in our configuration only OPA and TOP (no ice etc.))
  printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 3

  # changes the keys and copy MY_SRC to your configurations
  cd $CDIR/$CONFIG
  cp $GITCLONE/BUILD_CFG/cpp_SEAsia_FABM.fcm cpp_$CONFIG.fcm
  cp -r -f $GITCLONE/BUILD_CFG/MY_SRC ./

  # Add fabm and ERSEM options in compiler (you can add or just copy the file)
  #in bldxag.cfg add
  #bld::excl_dep        use::fabm
  #bld::excl_dep        use::fabm_config
  #bld::excl_dep        use::fabm_types
  #bld::excl_dep        use::fabm_driver
  #bld::excl_dep        use::fabm_version
  #OR instead take the ready file
  cp $GITCLONE/BUILD_CFG/bldxag_FABM.cfg $WDIR/trunk_NEMOGCM_r8395/TOOLS/COMPILE/bldxag.cfg

  # Make configuration with all the above updates included
  cd $CDIR
  ./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 4 clean
  ./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 4
  #################################################################
  cd $WORK


  #  A successful compilation will generate a ``nemo.exe`` executable in
  #   ``$NEMO/trunk_NEMOGCM_r8395/$CONFIG/BLD/bin/``

  #  Further information on the NEMO-ERSEM coupling system can be found here:
  #   `NEMO-ERSEM <https://github.com/NOC-MSM/NEMO_ERSEM>`_
