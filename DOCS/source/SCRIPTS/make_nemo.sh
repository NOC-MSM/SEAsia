#!/bin/bash

:'

***********************
make_nemo.sh
***********************

Checkout and compile the NEMO executable for physics only simulations
'
#::

  cd $NEMO
  # Checkout the code from the paris repository
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395
  # copy the appropriate acrhitecture file into place
  cp $GFILE/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/


  cd $CDIR
  # Start building...
  printf 'y\nn\nn\nn\nn\nn\nn\nn\n' | ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10
  # The previous command will fail but it has established some structures
  # continue with...
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

  ## Copy cloned files into their appropriate locations
  # FORTRAN modifications
  cp $GFILE/f_files/* $CDIR/$CONFIG/MY_SRC/.
  cp $WORK/MY_SRC/* $CDIR/$CONFIG/MY_SRC/.
  # list of compiler flags
  cp $GFILE/cpp_SANH.fcm $CONFIG/cpp_$CONFIG.fcm


  # Then we are ready to compile NEMO
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

  # Finally move files to the required locations
  cp  $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe
  cp $CDIR/$CONFIG/EXP00/* $EXP/

  cd $SCRIPTS
