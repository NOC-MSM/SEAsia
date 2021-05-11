#!/bin/bash

#:'
#
#***********************
#make_tools.sh
#***********************
#
# At the time of writing the tools used here required an older version of XIOS to
#work (XIOS1). This is compiled in make_xios.sh. Refer to the README in the
# $TDIR/TOOLS/DOMAINcfg for more details.
#
#Also make the SOSIE tool for flood filling land in parent data when interpolating
# onto child grid
#'
#::

  # Checkout the NEMO code from the SVN Paris repository. The tools have not
  # been updated for a while. Consequently there is a rather clumsy nested
  # TOOLS/TOOLS directory because a parent directory (to the child TOOLS) was
  # required for the compatible fcm tools.
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 $WDIR/BUILD_EXE/TOOLS_r8395



  # Copy the compiler option files from the NEMO-RELOC repository clone
  #  into the architecture directory
  cp $WDIR/HPC_ARCH_FILES/NEMO/arch-X86_ARCHER2-Cray.fcm $WDIR/BUILD_EXE/TOOLS_r8395/ARCH/arch-X86_ARCHER2-Cray.fcm

  # Edit ARCH file
  # Dirty fix to hard wire path otherwise user will have to set XIOS_DIR in every new shell session
  sed "s?XXX_XIOS_DIR_XXX?$XIOS1_DIR?" $WDIR/BUILD_EXE/TOOLS_r8395/ARCH/arch-X86_ARCHER2-Cray.fcm > tmp_arch
  mv tmp_arch $WDIR/BUILD_EXE/TOOLS_r8395/ARCH/arch-X86_ARCHER2-Cray_XIOS1.fcm


  # Apply patches for the weight file code
  cd $WDIR/BUILD_EXE/TOOLS_r8395/TOOLS/WEIGHTS/src
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp_mod.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripinterp.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scrip.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripshape.patch
  patch -b < $WDIR/BUILD_EXE/patch_files/scripgrid.patch

  #load modules
  module -s restore /work/n01/shared/acc/n01_modules/ucx_env

  # compile tools
  cd $WDIR/BUILD_EXE/TOOLS_r8395/TOOLS
  ./maketools -m X86_ARCHER2-Cray_XIOS1 -n NESTING
  ./maketools -m X86_ARCHER2-Cray_XIOS1 -n REBUILD_NEMO
  ./maketools -m X86_ARCHER2-Cray_XIOS1 -n WEIGHTS
  ./maketools -m X86_ARCHER2-Cray_XIOS1 -n DOMAINcfg



  # Make SOSIE tool.
  #$SCRIPTS/install_sosie.sh

  # Install PyNEMO
  #$SCRIPTS/install_pynemo.sh

  cd $WORK
