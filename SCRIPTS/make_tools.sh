#!/bin/bash

#:'
#
#***********************
#make_tools.sh
#***********************
#

  # Checkout the NEMO code from the SVN Paris repository. The tools have not
  # been updated for a while. But will be soon... 
  
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/utils/tools_r4.0-HEAD/ $WDIR/BUILD_EXE/TOOLS_HEAD

  cd $WDIR/BUILD_EXE/TOOLS_HEAD
  
  # Assuming that make_nemo.sh has already been done then all the external can be reused.
  cp $NEMO/arch/ .
  ln -s $NEMO/mk .
  ln -s $NEMO/ext .

  # Apply patches for the weight file code
  cd $WDIR/BUILD_EXE/TOOLS_HEAD/WEIGHTS/src
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp_mod.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scrip.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripshape.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripgrid.patch

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env

  # compile tools
  cd $WDIR/BUILD_EXE/TOOLS_HEAD
  ./maketools -m X86_ARCHER2-Cray -n NESTING
  ./maketools -m X86_ARCHER2-Cray -n REBUILD_NEMO
  ./maketools -m X86_ARCHER2-Cray -n WEIGHTS
  ./maketools -m X86_ARCHER2-Cray -n DOMAINcfg


  # Make SOSIE tool.
  #$SCRIPTS/install_sosie.sh

  # Install PyNEMO
  #$SCRIPTS/install_pynemo.sh

  cd $WORK
