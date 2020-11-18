cd $TDIR/NESTING

#load modules
module unload nco cray-netcdf cray-hdf5
module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf-hdf5parallel
module load cray-hdf5-parallel

# subdomain of ORCA global
#you can download the ORCA R12 coordinates or in ARCHER you can take it from my directory
cp /work/n01/n01/annkat/EXTRA_TOOLS/GRIDS/coordinates_ORCA_R12.nc $TDIR/NESTING/.
cp $GITCLONE/DOMAIN/namelist.input $TDIR/NESTING/

./agrif_create_coordinates.exe

#copy it to you DOMAIN to safe keep it
cp 1_coordinates_ORCA_R12.nc $DOMAIN/coordinates.nc

cd $WORK
