#!/bin/bash

#:'
#
#***********************
#make_xios1.sh
#***********************
#
#Checkout and compile the XIOS1 executable for I/O management
#You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register
#
#Check the path in variable ``%XIOS_HOME`` in ``XC_ARCHER_INTEL_XIOS1`` and
#``arch-X86_ARCHER2-Cray*`` are consistent with your settings
#
# THIS IS USED TO BUILD NEMO TOOLS WHICH ARE OLD AND REQUIRE AN OLD XIOS VERSION.
#'
#::


cd $WDIR
# Ensure the correct modules are loaded for ARHCER2
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

#download xios
#svn checkout -r 1964 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_DIR
svn checkout -r 703 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0 $XIOS1_DIR

cd $XIOS1_DIR

#copy the arch files to build location
cp $WDIR/HPC_ARCH_FILES/XIOS/arch-X86_ARCHER2-Cray_XIOS1.* $XIOS1_DIR/arch/.


#compile xios
./make_xios --prod --arch X86_ARCHER2-Cray_XIOS1 --netcdf_lib netcdf4_par --job 16 --full

# First time compile will fail
# got to $XIOS_DIR/tools/FCM/lib/Fcm/Config.pm and change
# FC_MODSEARCH => '',             # FC flag, specify "module" path
#to
#FC_MODSEARCH => '-J',           # FC flag, specify "module" path
sed -i "s/FC_MODSEARCH => ''/FC_MODSEARCH => '-J'/g" tools/FCM/lib/Fcm/Config.pm

#recompile xios
./make_xios --prod --arch X86_ARCHER2-Cray_XIOS1 --netcdf_lib netcdf4_par --job 16 --full

echo "Executable is $XIOS1_DIR/bin/xios_server.exe"

#######################################
cd $WDIR
