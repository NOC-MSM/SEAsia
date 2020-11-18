cd $NEMO 
#################################################################
#first get/download NEMO and FABM ERSM
#################################################################
#get nemo
svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

#replace TOP_SRC_old
mv trunk_NEMOGCM_r8395/NEMO/TOP_SRC trunk_NEMOGCM_r8395/NEMO/TOP_SRC_old
#original NEMO-FABM coupler can be found form here: https://github.com/NOC-MSM/NEMO_ERSEM/tree/master/TOP_SRC_r8395_FABM 
cp -r $GITCLONE/NEMO-FABM-ERSM/TOP_SRC_r8395_FABM trunk_NEMOGCM_r8395/NEMO/TOP_SRC

cd $WDIR
# get ERSM
cp -r $GITCLONE/NEMO-FABM-ERSM/ERSEM-master ./

cd $FABM
# get FABM
git clone https://github.com/fabm-model/fabm.git
#######################################################################

#compile FABM
#######################################################################
#load modules
module load cdt/15.11
module unload PrgEnv-cray PrgEnv-gnu
module load PrgEnv-intel/5.2.82
module unload cray-netcdf cray-hdf5
module load cray-netcdf-hdf5parallel/4.3.3.1
module load cray-hdf5-parallel/1.8.14

module load cmake

cd $FABM

cmake $FABM/fabm/src -DFABM_HOST=nemo -DCMAKE_Fortran_COMPILER=ifort -DFABM_ERSEM_BASE=$WDIR/ERSEM-master -DFABM_EMBED_VERSION=ON

make install
#change the directory you installed fabm to your directory
cp -r /home/n01/n01/annkat/local/fabm/nemo/lib ./
cp -r /home/n01/n01/annkat/local/fabm/nemo/include ./
#################################################################

#compile nemo
#################################################################
# get arch 
#ATTENTION modify the following file to have the correct paths
cp $NEMO/trunk_NEMOGCM_r8395/NEMO/TOP_SRC/arch-XC_ARCHER_INTEL_FABM.fcm $NEMO/trunk_NEMOGCM_r8395/ARCH/arch-XC_ARCHER_INTEL_FABM.fcm

cd $CDIR
#make configuration first
printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 0

#changes the keys and copy MY_SRC to your configurations
cd $CDIR/$CONFIG
cp $GITCLONE/NEMO-FABM-ERSM/cpp_SEAsia_FABM.fcm cpp_$CONFIG.fcm
cp -r -f $GITCLONE/NEMO-FABM-ERSM/MY_SRC ./

#add fabm and ERSEM options in compiler (you can add or just copy the file)
#in bldxag.cfg add
#bld::excl_dep        use::fabm
#bld::excl_dep        use::fabm_config
#bld::excl_dep        use::fabm_types
#bld::excl_dep        use::fabm_driver
#bld::excl_dep        use::fabm_version
#OR instead take the ready file
cp $GITCLONE/NEMO-FABM-ERSM/bldxag_FABM.cfg $NEMO/trunk_NEMOGCM_r8395/TOOLS/COMPILE/bldxag.cfg

#make configuration with updates included 
cd $CDIR
./makenemo -n $CONFIG -m XC_ARCHER_INTEL_FABM  -j 4 clean
./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 4
#################################################################
cd $WORK
