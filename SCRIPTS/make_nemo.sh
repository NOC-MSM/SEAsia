#!/bin/bash

#:'
#
#***********************
#make_nemo.sh
#***********************
#
# Checkout and compile the NEMO executable for physics only simulations
#
# You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register
#
#Bare in mind that NEMO is being compiled for use with a particular version of
#XIOS. Ensure edits to ``%XIOS_HOME`` in the ``arch-*`` file are consistent with
#the XIOS build.
# Also any code hardwiring (E.g. user defined initial conditions) must be put
# into MY_SRC before makenemo is executed.
#
#'
#::

cd $WDIR

# Ensure the correct modules are loaded for ARCHER2
# Load modules listed in /work/n01/shared/nemo/setup
# Tested 10Jan22
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

# Checkout the code from the paris repository
#svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

# Checkout the NEMO code from the SVN Paris repository
echo "Checking out NEMO repository"
case "${NEMO_VER}"
  in
  4.0.6)   echo "NEMO Verion 4.0.6 will be checked out"
           ;;
  4.0.7)   echo "NEMO Verion 4.0.7 will be checked out"
           ;;
  *)       echo "NEMO Version not recognised"
           echo "Versions available at present: 4.0.6, 4.0.7"
           exit 1
esac
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER --depth empty $NEMO
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/src --depth infinity $NEMO/src
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/SHARED $NEMO/cfgs/SHARED
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/AMM12 $NEMO/cfgs/AMM12
svn export http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/ref_cfgs.txt $NEMO/cfgs/ref_cfgs.txt

cd $NEMO
# Now check EXTERNALS revision number before checking out the rest
for ext_name in mk FCM IOIPSL
  do
  ext=`svn propget svn:externals | grep $ext_name | cut -c2-`
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/$ext
done

ext=`svn propget svn:externals | grep makenemo | cut -c2-`
svn export http://forge.ipsl.jussieu.fr/nemo/svn/$ext

cd $NEMO/ext/FCM/lib/Fcm
sed -e "s/FC_MODSEARCH => '',  /FC_MODSEARCH => '-J',/" Config.pm > tmp_file
mv tmp_file Config.pm

# copy the appropriate architecture file into place
mkdir $NEMO/arch
cp $WDIR/HPC_ARCH_FILES/NEMO/arch-${HPC_TARG}_${COMPILER}.fcm $NEMO/arch/arch-${HPC_TARG}_${COMPILER}.fcm

# Edit ARCH file
# Dirty fix to hard wire path otherwise user will have to set XIOS_DIR in every new shell session
sed "s?XXX_XIOS_DIR_XXX?$XIOS_DIR?" $NEMO/arch/arch-${HPC_TARG}_${COMPILER}.fcm > tmp_arch
mv tmp_arch $NEMO/arch/arch-${HPC_TARG}_${COMPILER}.fcm
# Add which modules are required for NEMO build
#echo $CONFIG 'OCE' >> $NEMO/cfgs/work_cfgs.txt

# Build from a reference configuration. This only uses OCE module
cd $NEMO
./makenemo -n $CONFIG -r AMM12  -m ${HPC_TARG}_${COMPILER} -j 16

# Change the keys and copy MY_SRC to your configurations
cp $NEMO/../cpp_file.fcm $NEMO/cfgs/$CONFIG/cpp_$CONFIG.fcm
cp -rf $NEMO/../MY_SRC $NEMO/cfgs/$CONFIG/.

# At this point make changes in MY_SRC if initial conditions are to be
# hardwired in at compile time. This can be done manually or triggered (as
# here) by CONFIG name
case "${CONFIG}"
  in
  NEMOconstTS)   cp $NEMO/cfgs/$CONFIG/MY_SRC/usrdef_istate.F90_constTS $NEMO/cfgs/$CONFIG/MY_SRC/usrdef_istate.F90
           ;;    # constant T,S initial condition:
  NEMOhorizTS)   cp $NEMO/cfgs/$CONFIG/MY_SRC/usrdef_istate.F90_horizTS $NEMO/cfgs/$CONFIG/MY_SRC/usrdef_istate.F90
           ;;    # horiz const, depth varying T,S initial condition
  *)       echo "Use default usrdef_istate.F90"
           exit 1
esac

# Make a refence experiment directory with various dummy files to allow it to compile
mkdir $NEMO/cfgs/$CONFIG/EXPREF
cp $NEMO/cfgs/SHARED/*namelist* $NEMO/cfgs/$CONFIG/EXPREF/.
cp $NEMO/cfgs/SHARED/*.xml $NEMO/cfgs/$CONFIG/EXPREF/.

#make configuration with updates included
./makenemo -r $CONFIG -m ${HPC_TARG}_${COMPILER} -j 16 clean
./makenemo -r $CONFIG -m ${HPC_TARG}_${COMPILER} -j 16

echo "Executable is $NEMO/cfgs/$CONFIG/BLD/bin/nemo.exe"

cd $WDIR/SCRIPTS
  #::
