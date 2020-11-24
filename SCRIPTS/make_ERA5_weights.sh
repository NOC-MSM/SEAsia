#!/bin/bash

:'

*******************************
make_ERA5_weights.sh
*******************************

Create WEIGHTS for the atmospheric forcing
^^^^^^^^^^^^^^^^^^^^

To use the ERA5 forcing fields weights file need to be generated that map the
source grid onto the target grid. This is done using the SCRIP tool, that was
previously built.
'

#::

  cd $SBC

  #The ERA5 fields have been downloaded in livljobs
  ln -s $DOMAIN/coordinates.nc $SBC/.

  #load modules
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

  #obtain namelists
  cp $GITCLONE/ATMOS_FORCING/namelist_reshape_bicubic_atmos $SBC/.
  cp $GITCLONE/ATMOS_FORCING/namelist_reshape_bilin_atmos $SBC/.

  #changes in the namelists manually if you have too
  #sed -i 's,ERA5_MSL_y2017.nc,/work/n01/n01/annkat/SEAsia_ERSEM_CMEMS/SURFACE_FORCING/ERA5_MSDWLWRF_y2017.nc,g' namelist_reshape_bilin_atmos
  #sed -i 's,/work/n01/n01/acc/DFS/DFS5.2/1960/drowned_precip_DFS5.2_y1960.nc,/work/n01/n01/annkat/SEAsia_ERSEM_CMEMS/SURFACE_FORCING/ERA5_MSDWLWRF_y2017.nc,g' namelist_reshape_bicubic_atmos

  #generate weights. Requiring make_tools.sh to have previously been run
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

  cd $WORK

:`

Alternatively makes weights files for the DFS forcing fields, that forced the
parent model. These, therefore need to pre-processing but weights files are needed
for the reduced domain.

These could be manually generated as follows.

Obtain namelist files and data file::

  cp $START_FILES/namelist_reshape_bilin_atmos $INPUTS/.
  cp $START_FILES/namelist_reshape_bicubic_atmos $INPUTS/.

Edit namelist to reflect source filenames (just a year change)::

  vi $WDIR/INPUTS/namelist_reshape_bilin_atmos
  ...
  &grid_inputs
    input_file = '/work/n01/n01/acc/DFS/DFS5.2/1960/drowned_precip_DFS5.2_y1960.nc'

  vi $WDIR/INPUTS/namelist_reshape_bicubic_atmos
  ...
  &grid_inputs
  input_file = '/work/n01/n01/acc/DFS/DFS5.2/1960/drowned_precip_DFS5.2_y1960.nc'


Setup weights files for the atmospheric forcing. Use the pre-compiled tools::

  export OLD_TDIR=$WORK/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

Setup weights files for the atmospheric forcing::

  cd $INPUTS
  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos

Generate  remap files ``remap_nemo_grid_atmos.nc`` and ``remap_data_grid_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos

Generates ``data_nemo_bilin_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos

Generates ``weights_bilinear_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos

Generates ``data_nemo_bicubic_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

Generates ``weights_bicubic_atmos.nc``.

'
