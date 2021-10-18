#ERA5 files have been downloaded in livlojobs and transferred to $SBC in ARCHER2

cd $SBC

#update namelists to reflect paths and names (only if doing something new)
# namelist_reshape_bicubic_atmos
# namelist_reshape_bilin_atmos

#link coordinate file
ln -s $DOMAIN/coordinates.nc $SBC/.

#load modules
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

#generate weights
$TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
$TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
$TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

cd $WDIR
