#!/bin/bash

:'

*******************************
make_coordinates_from_gebco.sh
*******************************

Make a netCDF coordinates file holds that holds the horizontal grid coordinate
position and spacings information. This is used to construct the 3D version:
the domain_cfg.nc file.

This script generates a coordinates.nc file using the gebco bathymetry data.

'





Hmm try and jump straight to DOMAINcfg with only lat and lon points


Extract from the parent bathymetry the size, start lat and lon, and resolution::


  y = 1363
        x = 2410

::

  vi namelist_cfg

  jpidta      =    2410   !  1st lateral dimension ( >= jpi )
  jpjdta      =    1363   !  2nd    "         "    ( >= jpj )
  ...
  jpiglo      =    2410   !  1st dimension of global domain --> i =jpidta
  jpjglo      =    1363   !  2nd    -                  -    --> j  =jpjdta


  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
    jphgr_msh   =       1               !  type of horizontal mesh
    ppglam0     =  -1.6638888888888899             !  longitude of first raw and column T-point (jphgr_msh = 1)
    ppgphi0     =  50.53833333333690803             ! latitude  of first raw and column T-point (jphgr_msh = 1)
    ppe1_deg    =  0.0002777777778             !  zonal      grid-spacing (degrees)
    ppe2_deg    =  0.0002777777778             !  meridional grid-spacing (degrees)
    ppdzmin     =  999999.              !  Minimum vertical spacing
    pphmax      =  999999.              !  Maximum depth
    ldbletanh   =  .FALSE.              !  Use/do not use double tanf function for vertical coordinates
    ppa2        =  999999.              !  Double tanh function parameters
    ppkth2      =  999999.              !
    ppacr2      =  999999.              !
  /



This seems to work, in populating e1 and e2 fields.

NExt steps:
Do rigourously.


Then build into a coordinate.nc file::

  $archer
  module load nco

  # copy the desired variables
  ncks -O -C -a -v nav_lon,nav_lat,glamt,glamu,glamv,glamf,gphit,gphiu,gphiv,gphif,e1t,e1u,e1v,e1f,e2t,e2u,e2v,e2f,ff_f,ff_t domain_cfg.nc coordinates_Solent_R3600.nc


This creates a new coordinates file with contents, which is now copied to
  INPUTS::

    cp coordinates_Solent_R3600.nc $INPUTS/coordinates.nc



Now we need to generate a bathymetry on this new grid.





#::



  cd $TDIR/NESTING

  #load modules
  module unload nco cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  # subdomain of ORCA global
  # you can download the ORCA R12 coordinates
  wget http://gws-access.jasmin.ac.uk/public/nemo/runs/ORCA0083-N06/domain/coordinates.nc -O $TDIR/NESTING/coordinates_ORCA_R12.nc

  #  or in ARCHER you can take it from my directory
  # cp /work/n01/n01/jelt/EXTRA_TOOLS/GRIDS/coordinates_ORCA_R12.nc $TDIR/NESTING/.

  # Get the namelist file from the cloned NEMO-RELOC repository
  cp $GITCLONE/DOMAIN/namelist.input $TDIR/NESTING/

:'

The namelist.input file controls the create process. Edit the bounding
coordinates (imax, ..., jmax) and scale factor (rho, rhot) to suit the child
coordinates. rho=rhot=7 with increase the resolution by a factor of 7.
'

#::

  cd$TDIR/NESTING/

  # Generate the new coordinates file with the namelist.input settings
  ./agrif_create_coordinates.exe

  # copy it to your DOMAIN folder where it will be used to create the domin_cfg.nc file
  cp 1_coordinates_ORCA_R12.nc $DOMAIN/coordinates.nc

  cd $WORK
