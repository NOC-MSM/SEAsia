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
#module -s restore /work/n01/shared/acc/n01_modules/ucx_env
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

#download xios
#svn checkout -r 1964 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_DIR
cd $XIOS_DIR
cp -r /work/n01/shared/nemo/xios-2.5/* .

#copy the arch files to build location
cp $WDIR/HPC_ARCH_FILES/XIOS/arch-X86_ARCHER2-Cray.* $XIOS_DIR/arch/.

#compile xios
./make_xios --prod --arch X86_ARCHER2-Cray --netcdf_lib netcdf4_par --job 16 --full

# First time compile will fail
# got to $XIOS_DIR/tools/FCM/lib/Fcm/Config.pm and change
# FC_MODSEARCH => '',             # FC flag, specify "module" path
#to
#FC_MODSEARCH => '-J',           # FC flag, specify "module" path
sed -i "s/FC_MODSEARCH => ''/FC_MODSEARCH => '-J'/g" tools/FCM/lib/Fcm/Config.pm

#recompile xios
./make_xios --prod --arch X86_ARCHER2-Cray --netcdf_lib netcdf4_par --job 16 --full

echo "Executable is $XIOS_DIR/bin/xios_server.exe"

#######################################
cd $WDIR/SCRIPTS
