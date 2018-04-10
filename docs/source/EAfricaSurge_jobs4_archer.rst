.. _EAfrica_surge:

**************************************
Setting up surge model for East Africa
**************************************

To link to me use::
 
  :ref:`EAfrica_surge`

**NOTE** This receipe is still being updated and hasn't had a second run through yet so there might be some mistakes that need tweaking to get going. 
  
Summary of the following:

1. Get and setup the Met Office NEMO v.4 surge code on Archer. **Note** two features of this code are: first, a namelist_cfg flag (ln_2d) which sets a constant temp
   and salinity, and second a namelist_cfg flag (ln_usr, ln_use_sbc) which allows the user to set all surface boundary conditions to zero or to read in wind and pressure.
   This is trunk revision 8814. For reference, the flag ln_2d appears in istate.f90 and ln_usr appears in sbcmod.f90 and sbcice_lim.F90.

2. Create 12th degree coordinates file for East Africa region on Archer 

3. Generate bathymetry file for this region on Archer

4. Generate a domain configuration file suitable for surge model on Archer

5. Generate tide boundary conditions on livljobs4

6. Generate wind and pressure boundary conditions   

7. Run with tides only forcing on Archer 
  
  
1a) Get and build Met Surge Config
=================================

Login to Archer ::

  ssh -l $USER login.archer.ac.uk

Create file 'temporary_path_names_for_NEMO_build_surge' and add the following :: 
  
  export CONFIG=EAFRICA_SURGE
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/CONFIG
  export TDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

execute e.g. ::

  . ~/temporary_path_names_for_NEMO_build_surge

Get some required files ::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES
  cp $WORK/jelt/LBay/START_FILES/dommsk.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/bdyini.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WORK/jelt/LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  
Load modules ::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

1b) Build XIOS2 @ r1080
======================

Follow instructions at :ref:`build_XIOS2`
(Note the final instruction to link the xios_server.exe may not work if the file structure has not been set
up, leave it, we do it here anyway)

1c) Build NEMO
=============

Get NEMO branch ::

  mkdir $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM/dev_r8814_surge_modelling_Nemo4

Get the correct archer compiler options file ::

  cp /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Create a link to xios_server.exe ::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

Make NEMO ::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

If build finished then jump to next section. If build failed try :: 

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Check compile flags ::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys  key_nosignedzero key_diainstant key_mpp_mpi key_iomput

Build ::

 ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10


2) Generate coordinates file
============================

First obtain the parent coordinate file coordinates_ORCA_R12.nc ::

  cp $START_FILES/coordinates_ORCA_R12.nc $INPUTS/.

Using this NetCDF file first decide upon which indices will define your domain. Use ferret or python, for example, to explore the domain ::

  livljobs4$ cd /work/anwise/tpxo7_2
  ferret
  use coordinates_ORCA_R12.nc 
  shade/X=3745:4225/Y=1000:1434 E2T, nav_lon, nav_lat ; go fland

To find specific longitudes and latitude you can for example ::

  ipython
  from netCDF4 import Dataset
  import numpy as np
  fn = 'coordinates_ORCA_R12.nc'
  nc_fid = Dataset(fn,'r')
  lat = nc_fid.variables['nav_lat'][:]
  lon = nc_fid.variables['nav_lon'][:]
  np.abs(lon[1000,:] - 20.0).argmin()
  np.abs(lon[1000,:] - 65.0).argmin()
  np.abs(lat[:,4225] - -5.0).argmin()
  np.abs(lat[:,4225] - -38.0).argmin()

It could also be useful to look at the TPXO harmonic amplitudes to find good cut off locations for boundaries. For example try ::

  livljobs4$ cd /work/anwise/tpxo7_2
  ferret
  go plot_EAfrica_harmonics.jnl

In this case we are using the interval i=3685:4225, j=1000:1434 which is approximately 20E-65E and 38S-5S. To obtain coordinates for
this domain create a namelist ::

  cd $TDIR/NESTING
  vim namelist.input

  &input_output
      iom_activated = true
  /
  &coarse_grid_files
      parent_coordinate_file = 'coordinates_ORCA_R12.nc'
  /
  &bathymetry
  /
  &nesting
      imin = 3685
      imax = 4225
      jmin = 1000
      jmax = 1434
      rho  = 1
      rhot = 1
      bathy_update = false
  /
  &vertical_grid
  /
  &partial_cells
  /
  &nemo_coarse_grid
  /
  &forcing_files
  /
  &interp
  /
  &restart
  /
  &restart_trc
  /

To build coordinates file see :ref:`build_and_create_coordinates`

Now copy to INPUTS ::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

3) Generate bathymetry file
===========================

For GEBCO bathymetry data head to BODC and download desired domain. Here we use 1-minute 2D dataset (2008) for 19E:66E, 39S:4S (we want the dataset to be spatially larger than the desired domain).
Copy NetCDF file to inputs e.g. ::

  scp GRIDONE_2D_19.0_-39.0_66.0_-4.0.nc $USER@login.archer.ac.uk:$INPUTS/.

Copy over namelist for reshaping bathymetry ::

  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Check that the lat and lon variable names are the same as in the data nc file. Now we need to flatten out the land elevations
and make the depths positive ::

  cd $INPUTS

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2D_19.0_-39.0_66.0_-4.0.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

Restore original modules ::
  
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

Execute script to map bathymetry to grid and generate remap_nemo_grid_gebco.nc and remap_data_grid_gebco.nc files ::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Execute script to generate data_nemo_bilin_gebco.nc file ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Execute script to generate bath_meter.nc file ::

 $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

4) Generate a domain configuration file
=======================================

Now generate a domain_cfg.nc file describing the vertial grid of the model.
In previous NEMO versions this would have been part of the main namelist_cfg.

Copy required files into DOMAINcfg directory ::

  cp $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
  cp $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Now edit the namelist_cfg file in the DOMAINcfg dirctory by following the instructions in :ref:`build_domain_cfg_file.rst`
for your desired domain setup. Here we use a 3 level s-coordinate set up ::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
    nn_no       =       0   !  job number (no more used...)
    cn_exp      =  "domaincfg"  !  experience name
    nn_it000    =       1   !  first time step
    nn_itend    =      75   !  last  time step (std 5475)
  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
    !
    ln_e3_dep   = .false.   ! =T : e3=dk[depth] in discret sens.
    !                       !      ===>>> will become the only possibility in v4.0
    !                       ! =F : e3 analytical derivative of depth function
    !                       !      only there for backward compatibility test with v3.6
    !                       !
    cp_cfg      =  "orca"   !  name of the configuration
    jp_cfg      =      12   !  resolution of the configuration
    jpidta      =     544   !  1st lateral dimension ( >= jpi )
    jpjdta      =     438   !  2nd    "         "    ( >= jpj )
    jpkdta      =      3    !  number of levels      ( >= jpk )
    jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
    jpjglo      =     438   !  2nd    -                  -    --> j  =jpjdta
    jpizoom     =       1   !  left bottom (i,j) indices of the zoom
    jpjzoom     =       1   !  in data domain indices
    jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
    ln_zco      = .false.   !  z-coordinate - full    steps
    ln_zps      = .false.   !  z-coordinate - partial steps
    ln_sco      = .true.   !  s- or hybrid z-s-coordinate
    ln_isfcav   = .false.   !  ice shelf cavity
    ln_linssh   = .false.   !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
    ln_s_sh94   = .true.    !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
    ln_s_sf12   = .false.   !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
    ln_sigcrit  = .false.   !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                            !  stretching coefficients for all functions
    rn_sbot_min =   6.0     !  minimum depth of s-bottom surface (>0) (m)
    rn_sbot_max =   7000.0  !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
    rn_hc       =   0.0     !  critical depth for transition to stretched coordinates
           !!!!!!!  Envelop bathymetry
    rn_rmax     =   0.3     !  maximum cut-off r-value allowed (0<r_max<1)
           !!!!!!!  SH94 stretching coefficients  (ln_s_sh94 = .true.)
    rn_theta    =   20.0    !  surface control parameter (0<=theta<=20)
    rn_bb       =   0.8     !  stretching with SH94 s-sigma
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
    nn_msh      =    0      !  create (=1) a mesh file or not (=0)
    rn_rdt      =   360.     !  time step for the dynamics (and tracer if nn_acc=0)
    ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
    ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
    ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
    ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
    ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
    ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
    ppsur       =  999999.0             !  ORCA r4, r2 and r05 coefficients
    ppa0        =  999999.0             ! (default coefficients)
    ppa1        =  999999.0             !
    ppkth       =      23.563           !
    ppacr       =       9.0             !
    ppdzmin     =       6.0             !  Minimum vertical spacing
    pphmax      =    5720.              !  Maximum depth
    ldbletanh   =  .FALSE.              !  Use/do not use double tanf function for vertical coordinates
    ppa2        =  999999.              !  Double tanh function parameters
    ppkth2      =  999999.              !
    ppacr2      =  999999.
  /
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
    ln_teos10   = .true.         !  = Use TEOS-10 equation of state
  /

Build a script to run the executable (change the email) ::

  #!/bin/bash
  #PBS -N domain_cfg
  #PBS -l walltime=00:20:00
  #PBS -l select=1
  #PBS -j oe
  #PBS -A n01-NOCL
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M anwise@noc.ac.uk
  #! -----------------------------------------------------------------------------

  # Change to the directory that the job was submitted from
  cd $PBS_O_WORKDIR

  # Set the number of threads to 1
  #   This prevents any system libraries from automatically
  #   using threading.
  export OMP_NUM_THREADS=1
  # Change to the directory that the job was submitted from
  ulimit -s unlimited

  #===============================================================
  # LAUNCH JOB
  #===============================================================
  echo `date` : Launch Job
  aprun -n 1 -N 1 ./make_domain_cfg.exe >&  stdouterr_cfg

  exit

Run it ::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Copy to EXP directory and also change permissions to ensure readable to others ::

  chmod a+rx $TDIR/DOMAINcfg/domain_cfg.nc
  rsync -uvt $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.      

5) Generate boundary conditions
===============================

First install pyNEMO/NRCT :ref:`install_nrct` (on livljobs4 currently)

Now set up the directory structure on livljobs4 ::

  cat > ~/temporary_path_names_for_NEMO_build << EOL
  export CONFIG=EAFRICA_SURGE
  export WORK=/work
  export WDIR=$WORK/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  EOL

Execute ::

  . ~/temporary_path_names_for_NEMO_build

Copy files from ARCHER ::

  cd /work/$USER
  mkdir NEMO
  mkdir $CONFIG
  mkdir $INPUTS
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/bathy_meter.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/domain_cfg.nc  $INPUTS/.

We require 6 files ::

  namelist.bdy
  inputs_src.ncml
  inputs_dst.ncml
  mask_src.nc
  mesh_hgr_src.nc
  mesh_zgr_src.nc

The last 3 of these files can be copied locally, try ::

  cd $INPUTS
  cp /work/anwise/NEMO/ACCORD/INPUTS/mask_src.nc .
  cp /work/anwise/NEMO/ACCORD/INPUTS/mesh_hgr_src.nc .
  cp /work/anwise/NEMO/ACCORD/INPUTS/mesh_zgr_src.nc .

Create ncml source file for ORCA12 source data ::

  vim inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>

Create ncml file for mapping variables in destination grid to what pyNEMO expects ::

  vim inputs_dst.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf location="file:domain_cfg.nc">
      <ns0:variable name="mbathy" orgName="top_level" />
      <ns0:variable name="gdept" orgName="gdept_0" />
      <ns0:variable name="gdepw" orgName="gdepw_0" />
      <ns0:variable name="e3u" orgName="e3u_0" />
      <ns0:variable name="e3v" orgName="e3v_0" />
      <ns0:variable name="e3t" orgName="e3t_0" />
      <ns0:variable name="e3w" orgName="e3w_0" />
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>

Create the namelist.bdy file. For tides only boundary forcing we only have 2D variables so vertical coordinates choice is arbitrary ::

  vim namelist.bdy

  !-----------------------------------------------------------------------
  !   vertical coordinate
  !-----------------------------------------------------------------------
  ln_zco      = .true.   !  z-coordinate - full    steps   (T/F)
  ln_zps      = .false.    !  z-coordinate - partial steps   (T/F)
  ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
  rn_hmin     =   -5     !  min depth of the ocean (>0) or
  !  min number of ocean level (<0)
  
  !-----------------------------------------------------------------------
  !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
  rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
  rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
  !  (= ocean depth) (>0) (m)
  ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
  rn_hc       =  50.0    !  critical depth with s-sigma
  
  !-----------------------------------------------------------------------
  !  grid information
  !-----------------------------------------------------------------------
  sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
  sn_src_zgr = './mesh_zgr_src.nc'   !  parent
  sn_dst_hgr = './domain_cfg.nc'
  sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
  sn_src_msk = './mask_src.nc'       ! parent
  sn_bathy   = './bathy_meter.nc'
  
  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
  sn_src_dir = './inputs_src.ncml'       ! src_files/'
  sn_dst_dir = '/work/anwise/NEMO/EAFRICA_SURGE/INPUTS/'
  sn_fn      = 'ACCORD'                 ! prefix for output files
  nn_fv      = -1e20                     !  set fill value for output files
  nn_src_time_adj = 0                                    ! src time adjustment
  sn_dst_metainfo = 'metadata info: jelt'
  
  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
  ln_coords_file = .true.               !  =T : produce bdy coordinates files
  cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
  ln_mask_file   = .false.              !  =T : read mask from file
  cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
  ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
  ln_tra         = .false.               !  boundary conditions for T and S
  ln_ice         = .false.               !  ice boundary condition
  nn_rimwidth    = 1                    !  width of the relaxation zone
  
  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
  ln_tide        = .true.               !  =T : produce bdy tidal conditions
  clname(1) ='M2'
  clname(2)='S2'
  clname(3)='K2'
  ln_trans       = .false.
  sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
  sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'
  
  !-----------------------------------------------------------------------
  !  Time information
  !-----------------------------------------------------------------------
  nn_year_000     = 1979        !  year start
  nn_year_end     = 1979        !  year end
  nn_month_000    = 11          !  month start (default = 1 is years>1)
  nn_month_end    = 11          !  month end (default = 12 is years>1)
  sn_dst_calendar = 'gregorian' !  output calendar format
  nn_base_year    = 1978        !  base year for time counter
  sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
  
  !-----------------------------------------------------------------------
  !  Additional parameters
  !-----------------------------------------------------------------------
  nn_wei  = 1                   !  smoothing filter weights
  rn_r0   = 0.041666666         !  decorrelation distance use in gauss
  !  smoothing onto dst points. Need to
  !  make this a funct. of dlon
  sn_history  = 'bdy files produced by jelt from ORCA0083-N01'
  !  history for netcdf file
  ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
  nn_alpha    = 0               !  Euler rotation angle
  nn_beta     = 0               !  Euler rotation angle
  nn_gamma    = 0               !  Euler rotation angle
  rn_mask_max_depth = 7000.0    !  Maximum depth to be ignored for the mask
  rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break
 
Generate the boundary condition files with pyNEMO ::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

Output required (for M2,S2,K2) ::

  ls -lh $INPUTS

  coordinates.bdy.nc
  ACCORD_bdytide_rotT_M2_grid_T.nc
  ACCORD_bdytide_rotT_K2_grid_T.nc
  ACCORD_bdytide_rotT_S2_grid_T.nc
  ACCORD_bdytide_rotT_M2_grid_U.nc
  ACCORD_bdytide_rotT_K2_grid_U.nc
  ACCORD_bdytide_rotT_S2_grid_U.nc
  ACCORD_bdytide_rotT_M2_grid_V.nc
  ACCORD_bdytide_rotT_K2_grid_V.nc
  ACCORD_bdytide_rotT_S2_grid_V.nc

Copy files back to ARCHER ::

  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.bdy.nc
  for file in $CONFIG*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done

**Note** that for tide only boundary conditions, it appears that changing the domain_cfg.nc file does not alter these files.

6a) Get wind and pressure data for surface boundary conditions on livljobs4
===========================================================================

First get data to feed into model. As an example we use some 2010 data from the ERA5 dataset from ECMWF that Nico downloaded.
The files are large so extract a subset using the latitude and longitude of our domain (slightly larger in fact)::

  cd $INPUTS
  module load nco/gcc/4.4.2.ncwa 
  
  ncea -d latitude,-39.0,-4.0 -d longitude,19.0,66.0 /projectsa/NEMO/Forcing/ERA5/INST/msl/msl_2010.nc EAfrica_msl_2010.nc
  
  ncea -d latitude,-39.0,-4.0 -d longitude,19.0,66.0 /projectsa/NEMO/Forcing/ERA5/INST/u10/u10_2010.nc EAfrica_u10_2010.nc

  ncea -d latitude,-39.0,-4.0 -d longitude,19.0,66.0 /projectsa/NEMO/Forcing/ERA5/INST/v10/v10_2010.nc EAfrica_v10_2010.nc

Now extract a forcing time period. Extract 2 weekly periods which coincide with a tropical depression in the 
vicinity. For example 21 Jan 10 - 3 Feb 10. To find the appropriate indices ::

  ipython

  import netCDF4 as nc
  import numpy as np
  from datetime import datetime

  ncfile = nc.Dataset('/work/anwise/DATA/ERA5/EAfrica_msl_2010.nc', 'r')
  time = ncfile.variables['time']
  dates = nc.num2date(time[:], time.units, time.calendar)

  # print index of date (year,month,day,hour,min) in time array
  idx = np.argwhere(dates==datetime(2010,1,21,0,0))
  idx
  time[idx]
  dates[idx].strftime("%T %A %d. %B %Y")
  ncfile.close()

Now extract the 2 weeks from the dataset (the out files are in required NEMO format) ::

  ncea -d time,480,647 EAfrica_msl_2010.nc EAfrica_msl_y2010m01d21.nc
  ncea -d time,648,815 EAfrica_msl_2010.nc EAfrica_msl_y2010m01d28.nc
  
  ncea -d time,480,647 EAfrica_u10_2010.nc EAfrica_u10_y2010m01d21.nc
  ncea -d time,648,815 EAfrica_u10_2010.nc EAfrica_u10_y2010m01d28.nc
  
  ncea -d time,480,647 EAfrica_v10_2010.nc EAfrica_v10_y2010m01d21.nc
  ncea -d time,648,815 EAfrica_v10_2010.nc EAfrica_v10_y2010m01d28.nc


6b) Generate surface boundary conditions (wind and pressure)
============================================================

Copy files over to Archer ::

  for file in EAfrica_*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done 

Copy over file required to create weight files ::

  ssh -l $USER login.archer.ac.uk
  
  cp $START_FILES/namelist_reshape_bilin_atmos $INPUTS/.
  cp $START_FILES/namelist_reshape_bicubic_atmos $INPUTS/.

Create a script file ::

  vim make_weights_bicubic.sh

  #!/bin/bash
  # Within the bilinear and bicubic namelists ensure that:
  # 1. the name of the data input_file is correct
  # 2. the name of the NEMO coordinate file for nemo_file is correct
  # 3. the input_lon, input_lat, nemo_lon, nemo_lat variables
  #    correspond to the approriate fields in your input_file and nemo_file

  /work/n01/n01/anwise/ACCORD/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
  /work/n01/n01/anwise/ACCORD/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
  /work/n01/n01/anwise/ACCORD/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
  /work/n01/n01/anwise/ACCORD/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
  /work/n01/n01/anwise/ACCORD/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

  rm remap_nemo_grid_atmos.nc
  rm remap_data_grid_atmos.nc
  rm data_nemo_bilin_atmos.nc
  rm weights_bilinear_atmos.nc
  rm data_nemo_bicubic_atmos.nc

We will run this script file 3 times, once for msl, once for u10 and once for v10.
Each time before running the script change some things in the 2 namelist files: 
1. input_file must be set to one of the 6 data files e.g. EAfrica_msl_y2010m01d21.nc
2. nemo_file must be set to the NEMO coordinate file we created eariler
3. input_lon, input_lat, nemo_lon, nemo_lat must be set to match what is used in the input_file and nemo_file respectively
4. In namelist_reshape_bicubic_atmos only, under &shape_inputs section, rename output_file to e.g. EAfrica_msl_weights_bicubic.nc

Run the script ::

  ./make_weights_bicubic.sh

This will have created a weight file. Change namelist as above for msl,u10,v10 and repeat. 
Should now have 3 files EAfrica_msl_weights_bicubic.nc, EAfrica_u10_weights_bicubic.nc and EAfrica_v10_weights_bicubic.nc.
Copy the data and weight files to experiment directory ::

  mkdir $EXP/fluxes
  for file in EAfrica_*nc; do rsync -utv $file $EXP/fluxes/$file ; done

7) Running model with tidal forcing at the boundaries on ARCHER
===============================================================

Copy files to EXP directory ::

  cd $EXP
  rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.

Link to the tide data ::

  ln -s $INPUTS $EXP/bdydta

Edit the namelist_cfg file ::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
    cn_exp      =  "AMMSURGE"  !  experience name
    nn_it000    = 1   !  first time step
    nn_itend    =  43200    !  last  time step (for dt = 6 min, 240*dt = 1 day)
    nn_date0    =  20130101 !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
    nn_time0    =       0   !  initial time of day in hhmm
    nn_leapy    =       1   !  Leap year calendar (1) or not (0)
    ln_rstart   =  .false.  !  start from rest (F) or from a restart file (T)
      nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
      nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
      !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
      !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
      !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
      cn_ocerst_in    = "AMMSURGE_00043920_restart"   !  suffix of ocean restart name (input)
      cn_ocerst_indir = "./Restart_files"         !  directory from which to read input ocean restarts
      cn_ocerst_out   = "restart"   !  suffix of ocean restart name (output)
      cn_ocerst_outdir= "./Restart_files"         !  directory in which to write output ocean restarts
    nn_istate   =       0   !  output the initial state (1) or not (0)
    nn_stock    =   43200    !  frequency of creation of a restart file (modulo referenced to 1)
    nn_write    =   43200    !  frequency of write in the output file   (modulo referenced to nit000)
  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
                            !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
     cn_domcfg = "domain_cfg"         ! domain configuration filename
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     ln_2d        = .true.  !  (=T) run in 2D barotropic mode (no tracer processes or vertical diffusion)
     rn_rdt      =   360.    !  time step for the dynamics (and tracer if nn_acc=0)
  /

  !-----------------------------------------------------------------------
  &namtsd    !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
     ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)
  /
  !-----------------------------------------------------------------------
  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
     nn_fsbc     = 1         !  frequency of surface boundary condition computation
                             !     (also = the frequency of sea-ice model call)
     ln_usr = .true.
     ln_blk =  .false.
     ln_apr_dyn  = .false.    !  Patm gradient added in ocean & ice Eqs.   (T => fill namsbc_apr )
     nn_ice      = 0         !  =0 no ice boundary condition   ,
     ln_rnf      = .false.   !  Runoffs                                   (T => fill namsbc_rnf)
     ln_ssr      = .false.   !  Sea Surface Restoring on T and/or S       (T => fill namsbc_ssr)
     ln_traqsr   = .false.   !  Light penetration in the ocean            (T => fill namtra_qsr)
     nn_fwb      = 0         !  FreshWater Budget: =0 unchecked
  /
  !-----------------------------------------------------------------------
  &namsbc_usr  !   namsbc_surge   surge model fluxes
  !-----------------------------------------------------------------------
     ln_use_sbc  = .false.    ! (T) to turn on surge fluxes (wind and pressure only)
                              ! (F) for no fluxes (ie tide only case)

  !
  !              !  file name                    ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation !
  !              !                               !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  !
     sn_wndi     = 'windspd_u_amm7'              ,       1           ,'x_wind',   .true.     , .false. , 'daily'  ,'' , ''
     sn_wndj     = 'windspd_v_amm7'              ,       1           ,'y_wind',   .true.     , .false. , 'daily'  ,'' , ''
     cn_dir      = './fluxes/'          !  root directory for the location of the bulk files
     rn_vfac     = 1.                   !  multiplicative factor for ocean/ice velocity
                                        !  in the calculation of the wind stress (0.=absolute winds or 1.=relative winds)
     rn_charn_const = 0.0275
  /
  !-----------------------------------------------------------------------
  &namtra_qsr    !   penetrative solar radiation
  !-----------------------------------------------------------------------
     ln_traqsr   = .false.   !  Light penetration (T) or not (F)
     nn_chldta   =      0    !  RGB : Chl data (=1) or cst value (=0)
  /
  !-----------------------------------------------------------------------
  &namsbc_apr    !   Atmospheric pressure used as ocean forcing or in bulk
  !-----------------------------------------------------------------------
  !          !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !          !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_apr= 'pressure_amm7',        1         ,   'air_pressure_at_sea_level' ,    .true.    , .false., 'daily'   ,  ''      ,   ''     ,  ''
     cn_dir      = './fluxes/'!  root directory for the location of the bulk files
     rn_pref     = 101200.    !  reference atmospheric pressure   [N/m2]/
     ln_ref_apr  = .false.    !  ref. pressure: global mean Patm (T) or a constant (F)
     ln_apr_obc  = .true.     !  inverse barometer added to OBC ssh data
  /
  !-----------------------------------------------------------------------
  &namlbc        !   lateral momentum boundary condition
  !-----------------------------------------------------------------------
  !   rn_shlat    =     0     !  shlat = 0  !  0 < shlat < 2  !  shlat = 2  !  2 < shlat
                             !  free slip  !   partial slip  !   no slip   ! strong slip
  /

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     rdttideramp =    1.
     clname(1)     =   'M2'   !  name of constituent
     clname(2)     =   'S2'
     clname(3)     =   'K2'
  /
  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
     ln_bdy     = .true.
     nb_bdy         = 1                    !  number of open boundary sets
     cn_coords_file = 'bdydta/coordinates.bdy.nc' !  bdy coordinates files
     cn_dyn2d       = 'flather'            !
     nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                           !  = 1, bdy data are read in 'bdydata   .nc' files
                                           !  = 2, use tidal harmonic forcing data from files
                                           !  = 3, use external data AND tidal harmonic forcing
     cn_tra        =  'frs'                !
     nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state
                                           !  = 1, bdy data are read in 'bdydata   .nc' files
     nn_rimwidth   = 1                    !  width of the relaxation zone
  /
  !-----------------------------------------------------------------------
  &nambdy_tide     ! tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/ACCORD_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.
  /
  !-----------------------------------------------------------------------
  &nambfr        !   bottom friction
  !-----------------------------------------------------------------------
     nn_bfr      =    2      !  type of bottom friction :   = 0 : free slip,  = 1 : linear friction
                             !                              = 2 : nonlinear friction
     rn_bfri2    =    2.4e-3 !  bottom drag coefficient (non linear case)
     rn_bfeb2    =    0.0e0  !  bottom turbulent kinetic energy background  (m2/s2)
     ln_loglayer =   .false. !  loglayer bottom friction (only effect when nn_bfr = 2)
     rn_bfrz0    =    0.003  !  bottom roughness (only effect when ln_loglayer = .true.)
  /
  !-----------------------------------------------------------------------
  &nambbc        !   bottom temperature boundary condition
  !-----------------------------------------------------------------------
     ln_trabbc   = .false.   !  Apply a geothermal heating at the ocean bottom
  /
  !-----------------------------------------------------------------------
  &nambbl        !   bottom boundary layer scheme
  !-----------------------------------------------------------------------
     nn_bbl_ldf  =  0      !  diffusive bbl (=1)   or not (=0)
  /
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .true.         !  = Use TEOS-10 equation of state
  /
  !-----------------------------------------------------------------------
  &namdyn_vor    !   option of physics/algorithm (not control by CPP keys)
  !-----------------------------------------------------------------------
     ln_dynvor_een = .true.  !  energy & enstrophy scheme
  /
  !-----------------------------------------------------------------------
  &namdyn_hpg    !   Hydrostatic pressure gradient option
  !-----------------------------------------------------------------------
     ln_hpg_zps  = .false.   !  z-coordinate - partial steps (interpolation)
     ln_hpg_sco  = .true.    !  s-coordinate (Standard Jacobian scheme)
  /
  !-----------------------------------------------------------------------
  &namdyn_spg    !   surface pressure gradient   (CPP key only)
  !-----------------------------------------------------------------------
     ln_dynspg_ts = .true.    ! split-explicit free surface
     ln_bt_auto =    .true.           !  Set nn_baro automatically to be just below
                                         !  a user defined maximum courant number (rn_bt_cmax)
  /
  !-----------------------------------------------------------------------
  &namdyn_ldf    !   lateral diffusion on momentum
  !-----------------------------------------------------------------------
     !                       !  Type of the operator :
     ln_dynldf_blp  =  .true.   !  bilaplacian operator
     ln_dynldf_lap    =  .false.  !  bilaplacian operator
     !                       !  Direction of action  :
     ln_dynldf_lev  =  .true.   !  iso-level
                             !  Coefficient
     rn_ahm_0     = 60.0      !  horizontal laplacian eddy viscosity   [m2/s]
     rn_bhm_0     = -1.0e+9   !  horizontal bilaplacian eddy viscosity [m4/s]
  /
  !-----------------------------------------------------------------------
  &namzdf        !   vertical physics
  !-----------------------------------------------------------------------
     rn_avm0     =   0.1e-6  !  vertical eddy viscosity   [m2/s]          (background Kz if not "key_zdfcst")
     rn_avt0     =   0.1e-6  !  vertical eddy diffusivity [m2/s]          (background Kz if not "key_zdfcst")
     ln_zdfevd   = .false.   !  enhanced vertical diffusion (evd) (T) or not (F)
     nn_evdm     =    1      !  evd apply on tracer (=0) or on tracer and momentum (=1)
  /
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents ('key_diaharm')
  !-----------------------------------------------------------------------
      nit000_han = 1     ! First time step used for harmonic analysis
      nitend_han = 43200     ! Last time step used for harmonic analysis
      nstep_han  = 5         ! Time step frequency for harmonic analysis
      tname(1)   = 'M2'      ! Name of tidal constituents
      tname(2)   = 'S2'
      tname(3)   = 'K2'
  /
  !-----------------------------------------------------------------------
  &namwad       !   Wetting and Drying namelist
  !-----------------------------------------------------------------------
     ln_wd = .false.   !: key to turn on/off wetting/drying (T: on, F: off)
     rn_wdmin1=0.1     !: minimum water depth on dried cells
     rn_wdmin2 = 0.01  !: tolerrance of minimum water depth on dried cells
     rn_wdld   = 20.0  !: land elevation below which wetting/drying will be considered
     nn_wdit   =   10  !: maximum number of iteration for W/D limiter
  /

Edit to have 1 hr SSH output ::

  vi file_def_nemo.xml
  ...
  <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
   <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"   />
   </file>

   <file id="file20" name_suffix="_Tides" description="tidal harmonics" >
     <field field_ref="K2x"          name="K2x"      long_name="K2 Elevation harmonic real part"                       />
     <field field_ref="K2y"          name="K2y"      long_name="K2 Elevation harmonic imaginary part"                  />
     <field field_ref="M2x"          name="M2x"      long_name="M2 Elevation harmonic real part"                       />
     <field field_ref="M2y"          name="M2y"      long_name="M2 Elevation harmonic imaginary part"                  />
     <field field_ref="S2x"          name="S2x"      long_name="S2 Elevation harmonic real part"                       />
     <field field_ref="S2y"          name="S2y"      long_name="S2 Elevation harmonic imaginary part"                  />
    </file> 
   </file_group>

Ensure that file_def_nemo.xml is pointed to ::

  vim context_nemo.xml
  ...
  <!--
  ============================================================================================== 
      NEMO context
  ============================================================================================== 
  -->
  <context id="nemo">
  <!-- $id$ -->
  <!-- Fields definition -->
      <field_definition src="./field_def_nemo-opa.xml"/>   <!--  NEMO ocean dynamics                     -->

  <!-- Files definition -->
      <file_definition src="./file_def_nemo.xml"/>     <!--  NEMO ocean dynamics                     -->
      <!-- 
  ...
  </context>

Create short queue runscript (Change the email address) ::

  #!/bin/bash
  # ---------------------------
  #===============================================================
  # CLUSTER BITS
  #===============================================================
  #PBS -N EA_R12
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL
  #PBS -j oe
  #PBS -r n
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M anwise@noc.ac.uk

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  # Change to the direcotry that the job was submitted from
  cd $PBS_O_WORKDIR


  # Set the number of threads to 1
  #   This prevents any system libraries from automatically
  #   using threading.
  export OMP_NUM_THREADS=1
  # Change to the directory that the job was submitted from
  ulimit -s unlimited
  ulimit -c unlimited

  export NEMOproc=96 #550
  export XIOSproc=1

  #===============================================================
  # LAUNCH JOB
  #===============================================================
  echo `date` : Launch Job
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $NEMOproc -N 24 ./opa
  exit

Submit the job ::

  cd $EXP
  qsub -q short runscript
