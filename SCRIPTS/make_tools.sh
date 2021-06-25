#!/bin/bash

#:'
#
#***********************
#make_tools.sh
#***********************
#

  cd $NEMO
  for ext_name in tools
    do
    ext=`svn propget svn:externals | grep $ext_name | cut -c2-`
    echo $ext
    svn co http://forge.ipsl.jussieu.fr/nemo/svn/$ext
  done


  # Make an adjustment to the DOMAINcfg source code to accomodate more varied vertical coords
  cp $DOMAIN/domzgr.f90.melange $TDIR/DOMAINcfg/src/domzgr.f90

  # Apply patches for the weight file code
  cd $NEMO/tools/WEIGHTS/src
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp_mod.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scrip.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripshape.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripgrid.patch

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env

  # compile tools
  cd $NEMO/tools
  ./maketools -m X86_ARCHER2-Cray -n NESTING
  ./maketools -m X86_ARCHER2-Cray -n REBUILD_NEMO
  ./maketools -m X86_ARCHER2-Cray -n WEIGHTS
  ./maketools -m X86_ARCHER2-Cray -n DOMAINcfg

  # Make SOSIE tool.
  #$SCRIPTS/install_sosie.sh

  # Install PyNEMO
  #$SCRIPTS/install_pynemo.sh

  cd $WORK
