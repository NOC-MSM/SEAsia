#!/bin/bash

#:'
#
#***********************
#make_xios.sh
#***********************
#
#Checkout and compile the XIOS2.5 executable for I/O management
#You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register
#
#Check the path in variable ``%XIOS_HOME`` in ``XC_ARCHER_INTEL_XIOS1`` and
#``arch-X86_ARCHER2-Cray*`` are consistent with your settings
#
#'
#::


cd $WDIR
# Ensure the correct modules are loaded for ARHCER2
module restore $WDIR/HPC_ARCH_FILES/envs/ucx_env_${HPC_TARG}_${COMPILER}.fcm

#download xios
svn checkout -r 1964 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_DIR
cd $XIOS_DIR

#copy the arch files to build location
cp $WDIR/HPC_ARCH_FILES/XIOS/arch-${HPC_TARG}_${COMPILER}.* $XIOS_DIR/arch/.

#compile xios
./make_xios --prod --arch ${HPC_TARG}_${COMPILER} --netcdf_lib netcdf4_par --job 16 --full

# First time compile will fail
# got to $XIOS_DIR/tools/FCM/lib/Fcm/Config.pm and change
# FC_MODSEARCH => '',             # FC flag, specify "module" path
#to
#FC_MODSEARCH => '-J',           # FC flag, specify "module" path
sed -i "s/FC_MODSEARCH => ''/FC_MODSEARCH => '-J'/g" tools/FCM/lib/Fcm/Config.pm

#recompile xios
./make_xios --prod --arch ${HPC_TARG}_${COMPILER} --netcdf_lib netcdf4_par --job 16 --full

echo "Executable is $XIOS_DIR/bin/xios_server.exe"

#######################################
cd $WDIR/SCRIPTS
