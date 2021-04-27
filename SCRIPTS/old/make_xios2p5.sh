cd $WDIR
#####################################
#download install xios
####################################
#load modules
#module unload cray-mpich
#module load craype-network-ucx
#module load cray-mpich-ucx
#module load libfabric
#module load cray-hdf5-parallel
#module load cray-netcdf-hdf5parallel
#module load gcc
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

#download xios
svn checkout -r 1566 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_DIR
cd $XIOS_DIR

#copy the arch
cp /work/n01/shared/acc/xios-2.5/arch/arch-X86_ARCHER2-Cray.* arch/.
# change %BASE_CFLAGS    -w -D_GLIBCXX_USE_CXX11_ABI=0 into arch-X86_ARCHER2-Cray.fcm

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

#######################################
cd $WORK

