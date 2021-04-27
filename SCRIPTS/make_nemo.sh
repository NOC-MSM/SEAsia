#!/bin/bash

:'

***********************
make_nemo.sh
***********************

Checkout and compile the NEMO executable for physics only simulations

You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register

Bare in mind that NEMO is being compiled for use with a particular version of
XIOS. Ensure edits to ``%XIOS_HOME`` in the ``arch-*`` file are consistent with
the XIOS build.

'
#::

  cd $WDIR
  # Checkout the code from the paris repository
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

  # Checkout the NEMO code from the SVN Paris repository
  echo "Checking out NEMO repository"
  case "${NEMO_VER}"
    in
    4.0.6)   echo "NEMO Verion 4.0.6 will be checked out"
             ;;
    *)       echo "NEMO Version not recognised"
             echo "Versions available at present: 4.0.6"
             exit 1
  esac
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER --depth empty $NEMO
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/src --depth infinity $NEMO/src
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/SHARED $NEMO/cfgs/SHARED
  svn export http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/ref_cfgs.txt $NEMO/cfgs/ref_cfgs.txt


  exit;
  

  # copy the appropriate architecture file into place
  cp $GFILE/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/

  cd $CDIR

:'
Start building...
The prescribed options make a new config directory structure and say say YES to
OPA_SRC only. Edit if you have other plans)
'#::

  printf 'y\nn\nn\nn\nn\nn\nn\nn\n' | ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10
  # The previous command will fail but it has established some structures
  # continue with...
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

  ## Copy cloned files into their appropriate locations
  # FORTRAN modifications
  cp $GFILE/f_files/* $CDIR/$CONFIG/MY_SRC/.
  cp $WORK/MY_SRC/* $CDIR/$CONFIG/MY_SRC/.
  # copy the list of compiler flags
  cp $GFILE/cpp_SANH.fcm $CONFIG/cpp_$CONFIG.fcm

  # Then we are ready to compile NEMO
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

  # Finally move files to the required locations
  cp $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe
  cp $CDIR/$CONFIG/EXP00/* $EXP/

  cd $SCRIPTS

  #::
