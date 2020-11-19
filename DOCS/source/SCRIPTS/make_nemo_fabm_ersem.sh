#!/bin/bash

:'

************************
make_nemo_fabm_ersem.sh
************************

Checkout and compile the NEMO executable with FABM ERSEM for physics with
biogeochemistry
'
#::

  cd $NEMO
  #################################################################
  # First get/download NEMO and FABM ERSEM
  #################################################################
  # checkout NEMO from the paris repository
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

  # Replace TOP_SRC with a cloned version from the NEMO-RELOC repository clone
  mv trunk_NEMOGCM_r8395/NEMO/TOP_SRC trunk_NEMOGCM_r8395/NEMO/TOP_SRC_old
  # The NEMO-FABM coupler ia maintained here: https://github.com/NOC-MSM/NEMO_ERSEM/tree/master/TOP_SRC_r8395_FABM
  cp -r $GITCLONE/NEMO-FABM-ERSEM/TOP_SRC_r8395_FABM trunk_NEMOGCM_r8395/NEMO/TOP_SRC

  cd $WDIR
  # get ERSEM from the NEMO-RELOC repository clone
  cp -r $GITCLONE/NEMO-FABM-ERSEM/ERSEM-master ./

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
  # Compile nemo
  #################################################################
  # get arch
  #ATTENTION modify the following file to have the correct paths
  cp $NEMO/trunk_NEMOGCM_r8395/NEMO/TOP_SRC/arch-XC_ARCHER_INTEL_FABM.fcm $NEMO/trunk_NEMOGCM_r8395/ARCH/arch-XC_ARCHER_INTEL_FABM.fcm

  cd $CDIR
  # make configuration first
  printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 0

  # changes the keys and copy MY_SRC to your configurations
  cd $CDIR/$CONFIG
  cp $GITCLONE/NEMO-FABM-ERSEM/cpp_SEAsia_FABM.fcm cpp_$CONFIG.fcm
  cp -r -f $GITCLONE/NEMO-FABM-ERSEM/MY_SRC ./

  # Add fabm and ERSEM options in compiler (you can add or just copy the file)
  #in bldxag.cfg add
  #bld::excl_dep        use::fabm
  #bld::excl_dep        use::fabm_config
  #bld::excl_dep        use::fabm_types
  #bld::excl_dep        use::fabm_driver
  #bld::excl_dep        use::fabm_version
  #OR instead take the ready file
  cp $GITCLONE/NEMO-FABM-ERSEM/bldxag_FABM.cfg $NEMO/trunk_NEMOGCM_r8395/TOOLS/COMPILE/bldxag.cfg

  # Make configuration with updates included
  cd $CDIR
  ./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 4 clean
  ./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 4
  #################################################################
  cd $WORK
