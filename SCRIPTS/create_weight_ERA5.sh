#ERA5 files have been downloaded in livlojobs and transferred to $SBC in ARCHER2

cd $SBC

#update namelists to reflect paths and names (only if doing something new)
# namelist_reshape_bicubic_atmos
# namelist_reshape_bilin_atmos

#link coordinate file
ln -s $DOMAIN/coordinates.nc $SBC/.

#load modules
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7


#generate weights
$TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
$TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

cd $WDIR
