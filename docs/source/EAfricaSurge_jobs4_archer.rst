.. _AMM7_SURGE_archer:

**************************************
Setting up surge model for East Africa
**************************************

To link to me use::
 
  :ref:`EAfrica_surge`

1. Get and build Met Surge Config
=================================

Login to Archer ::

  ssh -l $USER login.archer.ac.uk

Create file 'temporary_path_names_for_NEMO_build_surge' and add the following :: 
  
  export CONFIG=AMM7_SURGE
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/CONFIG
  export TDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

execute e.g. ::

  . ~/temporary_path_names_for_NEMO_build_surge

Load modules ::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

2. Build XIOS2 @ r1080
======================

Follow instructions at :ref:`build_XIOS2`
(Note the final instruction to link the xios_server.exe may not work if the file structure has not been set
up, leave it, we do it here anyway)

3. Build NEMO
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

To change the variable name in the bdy_mask.nc file ::

  ncrename -v mask,bdy_msk bdy_mask.nc
  for i in 'SWPacific_bdytide*'; do scp $i anwise@login.archer.ac.uk:/work/n01/n01/anwise/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00_surge/bdydta/; done
  scp coordinates.bdy.nc anwise@login.archer.ac.uk:/work/n01/n01/anwise/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00_surge/bdydta/coordinates.bdy.nc

Setup for East Africa
=====================

Create useful shortcuts in temporary_path_names_for_NEMO_build_ACCORD ::

  export CONFIG=ACCORD
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_EAFRICA
  export WORK=/work/n01/n01
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
 
Make directories ::
  
  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES

Copy over essential files ::

  cp $WORK/jelt/LBay/START_FILES/dommsk.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/bdyini.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WORK/jelt/LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  cp $WORK/jelt/SEAsia/START_FILES/usrdef_istate.F90 $START_FILES/.
  cp $WORK/jelt/SEAsia/START_FILES/usrdef_sbc.F90    $START_FILES/.

3. Generate coordinates file
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

Generate bathymetry file
========================

For GEBCO bathymetry data head to BODC and download desired domain. Here we use 1-minute 2D dataset (2008) for 19E:66E, 39S:4S (we want the dataset to be spatially larger than the desired domain).
Copy NetCDF file to inputs ::

  scp GRIDONE_2D_19.0_-39.0_66.0_-4.0.nc anwise@login.archer.ac.uk:/work/n01/n01/anwise/ACCORD/INPUTS/

Copy over namelist for reshaping bathymetry ::

  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Check that the lat and lon variable names are the same as in the data nc file. Now we need to flatten out the land elevations
and make the depths positive ::

  cd $INPUTS

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

Restore original modules ::
  
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

Execute script to map bathymetry with grid generating remap_nemo_grid_gebco.nc and remap_data_grid_gebco.nc files ::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Execute script to generate data_nemo_bilin_gebco.nc file ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Execute script to generate bath_meter.nc file ::

 $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

Generate a domain configuration file
====================================

The idea is to copy the namelist_cfg file into the DOMAINcfg directory along with the input files and generate a domain_cfg.nc.
This contains information about the domain. In previous NEMO version this would have been part of the main namelist_cfg.

Copy required files into DOMAINcfg directory ::

  cp $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
  cp $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Now edit the namelist_cfg file in the DOMAINcfg dirctory by following the instructions in :ref:`build_domain_cfg_file.rst`
for your desired domain setup. Here we use the 2 level hybrid z-s coordinate set up ::

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
     jpkdta      =       2   !  number of levels      ( >= jpk )
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
  !  cn_coord_hgr = 'coordinates.nc'  ! File containing gphit (latitude) coordinate for use if ln_eq_taper=.true.
    rn_sbot_min =   6.0     !  minimum depth of s-bottom surface (>0) (m)
    rn_sbot_max = 7000.0    !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
    rn_hc       =   0.0     !  critical depth for transition to stretched coordinates
                         !!!!!!!  Envelop bathymetry
    rn_rmax     =    1.0    !  maximum cut-off r-value allowed (0<r_max<1)
                         !!!!!!!  SH94 stretching coefficients  (ln_s_sh94 = .true.)
    rn_theta    =    0.0    !  surface control parameter (0<=theta<=20)
    rn_bb       =    0.8    !  stretching with SH94 s-sigma
                         !!!!!!!  SF12 stretching coefficient  (ln_s_sf12 = .true.)
    rn_alpha    =    4.4    !  stretching with SF12 s-sigma
    rn_efold    =    0.0    !  efold length scale for transition to stretched coord
    rn_zs       =    1.0    !  depth of surface grid box
                            !  bottom cell depth (Zb) is a linear function of water depth Zb = H*a + b
    rn_zb_a     =    0.024  !  bathymetry scaling factor for calculating Zb
    rn_zb_b     =   -0.2    !  offset for calculating Zb
                         !!!!!!!! Other stretching (not SH94 or SF12) [also uses rn_theta above]
    rn_thetb    =    1.0    !  bottom control parameter  (0<=thetb<= 1)
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     rn_rdt      =   360.    !  time step for the dynamics (and tracer if nn_acc=0)
  !    rn_rdtmin   =   600.          !  minimum time step on tracers (used if nn_acc=1)
  !   rn_rdtmax   =   600.          !  maximum time step on tracers (used if nn_acc=1)
  !   rn_rdth     =   600.          !  depth variation of tracer time step  (used if nn_acc=1)
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
  &namlbc        !   lateral momentum boundary condition
  !-----------------------------------------------------------------------
     rn_shlat    =     0     !  shlat = 0  !  0 < shlat < 2  !  shlat = 2  !  2 < shlat
                             !  free slip  !   partial slip  !   no slip   ! strong slip
  /
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .true.         !  = Use TEOS-10 equation of state
  /

Build a script to run the executable ::

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

Run is ::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Copy to EXP directory and also change permissions to ensure readable to others ::

  chmod a+rx $TDIR/DOMAINcfg/domain_cfg.nc
  rsync -uvt $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.      

Generate boundary conditions
============================

First install pyNEMO/NRCT :ref:`install_nrct` (on livljobs4 currently)

Now set up the directory structure on livljobs4 ::

  cat > ~/temporary_path_names_for_NEMO_build << EOL
  export CONFIG=ACCORD
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
  cd NEMO
  mkdir $CONFIG
  cd $CONFIG
  mkdir INPUTS
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

The last 3 of these files can be copied locally try ::

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

Create the namelist.bdy file. For tides only we have only 2D variables so vertical coordinates choice is arbitrary ::

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
  sn_dst_dir = '/work/anwise/NEMO/ACCORD/INPUTS_2lvl/'
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

  ls -1 $INPUTS

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

Note that for tide only boundary conditions, it appears that changing the domain_cfg.nc file does not alter the bdytide
grid files.

Running model with tidal forcing at the boundaries on ARCHER
============================================================

Copy files to EXP directory ::

  cd $EXP
  rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.

Link to the tide data ::

  ln -s $INPUTS $EXP/bdydta

Edit to have 1 hr SSH output ::

  vi file_def_nemo.xml
  ...
  <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
   <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"   />
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
