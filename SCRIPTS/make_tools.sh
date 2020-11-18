cd $TDIR

#copy ARCH files that need old version of xios (tools sometimes do not work with new version)
#but you can try. For the XIOS1 to work you need access to jelt xios.1 or else you will have to
#download and compile this xios version
cp $GITCLONE/NEMO-FABM-ERSM/arch-XC_ARCHER_INTEL_NOXIOS.fcm ../ARCH/.
cp $GITCLONE/NEMO-FABM-ERSM//arch-XC_ARCHER_INTEL_XIOS1.fcm  ../ARCH/.

#patch for the weight files
cd $TDIR/WEIGHTS/src
patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripinterp_mod.patch
patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripinterp.patch
patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scrip.patch
patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripshape.patch
patch -b < $GITCLONE/NEMO-FABM-ERSM/patch_files/scripgrid.patch

#load modules
module unload nco cray-netcdf cray-hdf5
module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf-hdf5parallel
module load cray-hdf5-parallel
~
~
# compile tools
cd $TDIR
./maketools -m XC_ARCHER_INTEL_NOXIOS -n NESTING -j 6
./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO
./maketools -m XC_ARCHER_INTEL_XIOS1 -n WEIGHTS

cd $WORK
