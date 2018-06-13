.. _Solent_surge:

**************************************
Setting up surge model for East Africa
**************************************

To link to me use::

  :ref:`Solent_surge`

**NOTE** This receipe is a copy from :ref:`EAfrica_Surge` and using a working
version of the full 3D (constant T,S) in `<Solent_archer_livljobs_TCarta.rst>_`

Summary of the following:

1. Get and setup the Met Office NEMO v.4 surge code on Archer. **Note** two features of this code are: first, a namelist_cfg flag (ln_2d) which sets a constant temp
   and salinity, and second a namelist_cfg flag (ln_usr, ln_use_sbc) which allows the user to set all surface boundary conditions to zero or to read in wind and pressure.
   This is trunk revision 8814. For reference, the flag ln_2d appears in istate.f90 and ln_usr appears in sbcmod.f90 and sbcice_lim.F90.

2. copy coordinates file

3. copy bathymetry file

4. Generate a domain configuration file suitable for surge model on Archer

5. Generate tide boundary conditions on livljobs4

6. Run with tides only forcing on Archer

7. Output xml and running on Archer


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
  cp $WORK/jelt/SEAsia/git_repo/SEAsia/MY_SRC/stpctl.F90 $START_FILES/.
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

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM dev_r8814_surge_modelling_Nemo4

Get the correct archer compiler options file ::

  cp /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make NEMO ::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

If build finished then jump to next section. If build failed try ::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Check compile flags ::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys  key_nosignedzero key_diainstant key_mpp_mpi key_iomput

Minor edit to solver.stat output::

  cp $START_FILES/stpctl.F90  $CDIR/$CONFIG/MY_SRC/.

Build ::

 ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Create a link to xios_server.exe ::

 ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe


2) Generate coordinates file
============================

Copy the coordinates file from Solent config::

  cp /work/n01/n01/jelt/Solent/INPUTS/coordinates.nc $INPUTS/coordinates.nc


3) Generate bathymetry file
===========================

Copy the bathymetry file from Solent config::

  cp /work/n01/n01/jelt/Solent/INPUTS/bathy_meter.nc $INPUTS/bathy_meter.nc



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
    ln_e3_dep   = .true.   ! =T : e3=dk[depth] in discret sens.
    !                       !      ===>>> will become the only possibility in v4.0
    !                       ! =F : e3 analytical derivative of depth function
    !                       !      only there for backward compatibility test with v3.6
    !                       !
    cp_cfg      =  "orca"   !  name of the configuration
    jp_cfg      =    3600   !  resolution of the configuration
    jpidta      =    2410   !  1st lateral dimension ( >= jpi )
    jpjdta      =    1363   !  2nd    "         "    ( >= jpj )
    jpkdta      =       3   !  number of levels      ( >= jpk )
    jpiglo      =    2410   !  1st dimension of global domain --> i =jpidta
    jpjglo      =    1363   !  2nd    -                  -    --> j  =jpjdta
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
    rn_sbot_max =   100.0  !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
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
    rn_rdt      =   1.     !  time step for the dynamics (and tracer if nn_acc=0)
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
  #PBS -A n01-ACCORD
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk
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

Check the executable is there (or add it e.g.)::

  ln -s /work/n01/n01/jelt/Solent/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/make_domain_cfg.exe $TDIR/DOMAINcfg/.

Run it ::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Copy to EXP directory and also change permissions to ensure readable to others ::

  chmod a+rx $TDIR/DOMAINcfg/domain_cfg.nc
  rsync -uvt $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.

5) Generate boundary conditions
===============================

This is done for TPXO and FES data.

These were generated using TPXO data previously
::


  livljobs4
  cd /work/jelt/NEMO/Solent/INPUTS
  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/.

  cd TPXO
  for file in Solent_*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/TPXO/$file ; done

  cd ../FES
  for file in Solent_*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/FES/$file ; done





6) Running model with tidal forcing at the boundaries on ARCHER
===============================================================

Copy files to EXP directory ::

  cd $EXP
  rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.

.. note : Hmm I'm sure I don't need to copy bathy_meter.nc to EXP

Link to the tide data ::

  ln -s $INPUTS $EXP/bdydta

Edit the namelist_cfg file.
(chanage the lateral diffusion to laplacian = 25) ::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
    cn_exp      =  "Solent_surge"  !  experience name
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
      cn_ocerst_in    = "Solent_surge_00043920_restart"   !  suffix of ocean restart name (input)
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
     rn_rdt      =   1.    !  time step for the dynamics (and tracer if nn_acc=0)
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
     ln_tide_ramp = .true.
     rdttideramp =    0.166 # 4 hours
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
     filtide      = 'bdydta/Solent_bdytide_rotT_'         !  file name root of tidal forcing files
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
     ln_dynldf_blp  =  .false.   !  bilaplacian operator
     ln_dynldf_lap    = .true.  !  bilaplacian operator
     !                       !  Direction of action  :
     ln_dynldf_lev  =  .true.   !  iso-level
                             !  Coefficient
     rn_ahm_0     = 25.0      !  horizontal laplacian eddy viscosity   [m2/s]
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


9) Output xml and running
=========================

Edit to have 1 hr SSH output ::

  vi file_def_nemo.xml
  ...
  <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
   <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"      operation="instant"   />
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
  #PBS -N SolentSurg
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-ACCORD
  #PBS -j oe
  #PBS -r n
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk

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


---

progress
+++++++++

rn_rdt=1
ahm=25
kt = 252

stp_ctl : the ssh is larger than 10m
=======
kt=   252 max ssh:   10.42    , i j:  2400 1086

rn_rdt=1
ahm=60

stpctl: the zonal velocity is larger than 20 m/s
======
kt=    48 max abs(U):   23.62    , i j k:  2404 1083    2


rn_rdt=1
ahm=10

Runs 1hr( (3600steps) in about 10mins.
Run again with 1hrly output for 6 hours nt=21600 (1hr walltime)

stp_ctl : the ssh is larger than 10m
=======
kt=  9507 max ssh:   10.02    , i j:   190  502


Got two hours out with SSH fields.
Rebuild on 96 processors - something numerical happens west of Hurst in open water.

**Plan** generate a restart and fiddle with viscosity parameters.

Run for 9500 to generate a restart 25mins walltime.

*(11 June)*
restart
biharm + laplacian? - No. Use one or none.

rn_ahm_0 = 1.

SSH2 seems to go down. Run for 2hours and see if the currents are better. (Was
getting 3.9m/s downstream of Hurst point, to west, with ahm_0 = 10

stp_ctl : the ssh is larger than 10m
=======
kt=  3969 max ssh:   10.66    , i j:  1273  120


Try restoring laplacian diffusion ramping up the tides instead.
Run for 4 hours. Ramp up over 4 hours

ln_tide_ramp = .true.
rdttideramp =    0.1666
rn_ahm_0 = 10.
rn_rdt = 1.

This worked and ran out to 4hrs (14400 steps).

Restart and run for 14400 more hours. (no tidal ramp).

stp_ctl : the ssh is larger than 10m
=======
kt= 27929 max ssh:   10.27    , i j:  1308  121

*2 June 2018*
Ooo just missed finishing.
Try increasing the bottom friction (from 2.4e-3 to 3.e-3)

rn_bfri2    =    3.e-3 ! 2.4e-3 !  bottom drag coefficient (non linear case)

restart from 21600. Two hours on short queue.

stp_ctl : the ssh is larger than 10m
=======
kt= 28641 max ssh:   10.15    , i j:  1308  120


Try increasing bottom friction some more:
rn_bfri2    =    4.e-3 !  bottom drag coefficient (non linear case)

restart from 21600. Two hours on short queue.

Will blow up.

Increase friction by 10
2.3e-2
Cold start.

Worked ran to 7200.
Retarted. Ran and past peak in vel2 and sum(ssh2).

Restart for a long run and leave.

13 hours in 6 x 20mins + 10mins = 2 hours 10mins

 14401 --> 14400 + 46,800 = 61200

 **PENDING**



* Should inspect domain_cfg.nc. What are the e3t units? cm? ulikely...
* Should have restarting tides.
