==========================================
Setting up AMM60 in NEMO v4
==========================================

Machines: livljobs4, ARCHER

URL:: *to add*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds AMM60 in NEMOv4

Build on a combination of livljobs4 and ARCHER.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: *skip: Generate initial conditions and atmospheric forcings*
#. ARCHER: Run simulation

This will borrow forcings from Guihou et al 2018 AMM60 run
``/work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/``

Issues that arose
=================

* ...

.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on ARCHER.
Starting on ARCHER::

  ssh login.archer.ac.uk

  export CONFIG=AMM60
  export WORK=/work/n01/n01
  export WDIR=$WORK/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/
  export EXP=$CDIR/$CONFIG/EXP00

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel


  export AMM60DIR=/work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/

---

Collect essential files
=======================

Note you might have to mkdir the odd directory or two...::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES

  cp $WDIR/../LBay/START_FILES/dommsk.F90 $START_FILES/.
  cp $WDIR/../LBay/START_FILES/bdyini.F90 $START_FILES/.
  #cp $WDIR/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  #cp $WDIR/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  #cp $WDIR/../SWPacific/START_FILES/usrdef_istate.F90 $START_FILES/.
  #cp $WDIR/../SWPacific/START_FILES/usrdef_sbc.F90    $START_FILES/.


Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build (if it is already downloaded). Note here we use user defined
 functions for the initial state (constant T and S) and surface forcing (zero forcing)::

  cd $CDIR
  #cp $START_FILES/usrdef_istate.F90 $CDIR/$CONFIG/MY_SRC/.
  #cp $START_FILES/usrdef_sbc.F90    $CDIR/$CONFIG/MY_SRC/.

  cp $START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dommsk.F90 $CDIR/$CONFIG/MY_SRC/.

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

---

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

Or just link XIOS executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---


Build TOOLS
===========

To generate domain coords and rebuild tools we first need
to compile some of the NEMO TOOLS.

.. note: These are compiled with XIOS2. However DOMAINcfg has to be compiled
  with XIOS1. There is a README in the $TDIR/DOMAINcfg on what to do.

First build DOMAINcfg (which is relatively new and in NEMOv4). Use my XIOS1 file
(see userid and path in variable ``%XIOS_HOME``). Copy from ARCH *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $CDIR/../ARCH/.
  cd $TDIR

  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO

For the generation of bathymetry and met forcing weights files we need to patch
the code (to allow direct passing of arguments. NB this code has not been
updated in 7 years.)::

  cd $TDIR/WEIGHTS/src
  patch -b < $START_FILES/scripinterp_mod.patch
  patch -b < $START_FILES/scripinterp.patch
  patch -b < $START_FILES/scrip.patch
  patch -b < $START_FILES/scripshape.patch
  patch -b < $START_FILES/scripgrid.patch

  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n WEIGHTS




1. Generate new coordinates file
++++++++++++++++++++++++++++++++

Or just copy it::

  cp $AMM60DIR/EXP_harmIT2/WDIR/coordinates.nc $INPUTS/.



2. Generate bathymetry file
+++++++++++++++++++++++++++

Or just copy it::

  cp $AMM60DIR/EXP_harmIT2/WDIR/bathy_meter.nc $INPUTS/.



3. Generate initial conditions
++++++++++++++++++++++++++++++

SKIP

Skip this first time round. First test for stability with constant T and S.
Then try with tides.
Then try with initial conditions.

For constant T and S use the user defined functions in ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``.


.. note: Skip this for now.

    Copy ``make.macro`` file and edit the path if necessary::
    **FIX** to the notes (copied from jdha instead): ``cp $WDIR/INPUTS/make.macro ./``::

      cp /home/n01/n01/jdha/sosie/make.macro /home/n01/n01/jelt/sosie/.

      vi /home/n01/n01/jelt/sosie/make.macro
      # Directory to install binaries:
      INSTALL_DIR = /home/n01/n01/jelt/local

    Proceed with Step 6 (of Lighhouse Reef Readthedocs)::

      cd ~
      mkdir local
      svn co svn://svn.code.sf.net/p/sosie/code/trunk sosie
      cd sosie

      make
      make install
      export PATH=~/local/bin:$PATH
      cd $WDIR/INPUTS


    Obtain the fields to interpolate. Interpolate AMM60
    data. Get the namelists::

      cp $START_FILES/initcd_votemper.namelist $INPUTS/.
      cp $START_FILES/initcd_vosaline.namelist $INPUTS/.

    Generate the actual files. Cut them out of something bigger. Use the same indices
    as used in coordinates.nc (note that the nco tools don't like the
    parallel modules)::

    ----

    *(3 March 2017)*
    Insert new method to use AMM60 data for initial conditions.
    /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT
    AMM60_5d_20131013_20131129_grid_T.nc

    Find the AMM60 indices using FERRET on the bathy_meter.nc file: ``shade log(Bathymetry[I=540:750, J=520:820])``

    Note that the temperature and salinity variables are ``thetao`` and ``so``

    ::

      module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
      module load cray-netcdf cray-hdf5
      module load nco/4.5.0
      cd $INPUTS

      ncks -d x,560,620 -d y,720,800 /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT/AMM60_5d_20131013_20131129_grid_T.nc $INPUTS/cut_down_20131013_LBay_grid_T.nc

    Average over time and restore the parallel modules::

      ncwa -a time_counter $START_FILES/cut_down_20131013_LBay_grid_T.nc  $INPUTS/cut_down_201310_LBay_grid_T.nc

      module unload nco cray-netcdf cray-hdf5
      module load cray-netcdf-hdf5parallel cray-hdf5-parallel



    Edit namelists::

      vi initcd_votemper.namelist
      cf_in     = 'cut_down_201310_LBay_grid_T.nc'
      cv_in     = 'thetao'
      cf_x_in   = 'cut_down_201310_LBay_grid_T.nc'
      cv_out   = 'thetao'
      csource  = 'AMM60'
      ctarget  = 'LBay'

      vi initcd_vosaline.namelist
      ...
      cv_out   = 'so'
      ...



    Do stuff. I think the intention was for SOSIE to flood fill the land::

      sosie.x -f initcd_votemper.namelist

    Creates::

      thetao_AMM60-LBay_2013.nc4
      sosie_mapping_AMM60-LBay.nc

    Repeat for salinity::

      sosie.x -f initcd_vosaline.namelist

    Creates::

      so_AMM60-LBay_2013.nc4


    Now do interpolation as before. First copy the namelists::

      cp $START_FILES/namelist_reshape_bilin_initcd_votemper $INPUTS/.
      cp $START_FILES/namelist_reshape_bilin_initcd_vosaline $INPUTS/.

    Edit the input files::

      vi $INPUTS/namelist_reshape_bilin_initcd_votemper
      &grid_inputs
        input_file = 'thetao_AMM60-LBay_2013.nc4'
      ...

      &interp_inputs
        input_file = "thetao_AMM60-LBay_2013.nc4"
      ...

    Simiarly for the *vosaline.nc file::

      vi $INPUTS/namelist_reshape_bilin_initcd_vosaline
      &grid_inputs
        input_file = 'so_AMM60-LBay_2013.nc4'
      ...

      &interp_inputs
        input_file = "so_AMM60-LBay_2013.nc4"
      ...


    Produce the remap files::

      $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

    Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

      $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

    Creates ``data_nemo_bilin_R12.nc``. Then::

      $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

    Creates ``initcd_votemper.nc``. Then::

      $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

    Creates ``initcd_vosaline.nc``.



4. Generate a domain configuration file
=======================================

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg``
directory along with all the inputs files that would have previously been needed
get v3.6 running. The reason being that all the non-time stepping stuff, like
grid generating, has been abstracted from the core OPA code and is now done as
a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory::

    ln -s $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
    ln -s $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template ``namelist_cfg`` with only the essenetial domain building stuff.
Get the size of the new domain from ``ncdump -h bathy_meter.nc``.

Follow recipe of hybrid z-s coordinates in `build_domain_cfg_file.rst`_

NB Copy namelist_cfg settings from

``$AMM60DIR/EXP_harmIT2/namelist_cfg``::


  cd $TDIR/DOMAINcfg
  vi namelist_cfg

  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ...
     jpidta      =     1120               !  1st lateral dimension ( >= jpi )
     jpjdta      =     1440               !  2nd    "         "    ( >= jpj )
     jpkdta      =      51               !  number of levels      ( >= jpk )
     jpiglo      =     1120               !  1st dimension of global domain --> i =jpidta
     jpjglo      =     1440              !  2nd    -                  -    --> j  =jpjdta

     ...



When the ``domain_cfg.nc`` is built, copy it to the EXP directory (also copy it
 to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.

.. mote :  should check the difference between the homemade sco version the AMM60
  verison did:      ``diff namelist_cfg_sco_WIP namelist_cfg_AMM60``

.. note : alternativly should check the difference between the AMM60 and local
  output.namelist.dyn: ``diff output.namelist.dyn /work/n01/n01/jelt/NEMO/NEMOGCM/CONFIG/AMM60smago/EXP_NSea/output.namelist.dyn``
  I notice that rmax is different.



5. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++



6. Generate boundary conditions with NRCT/PyNEMO: Create netcdf abstraction wrapper
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


7. Edit the namelist_cfg
++++++++++++++++++++++++

::
  cd $AMM60DIR
  cp -r EXP_harmIT2 EXP_v4

Hack the namcfg fields::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
  ln_read_cfg = .true.   !  (=T) read the domain configuration file
     !                   !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
     cn_domcfg = "domain_cfg"         ! domain configuration filename
     !
  ln_write_cfg= .false.   !  (=T) create the domain configuration file
     cn_domcfg_out = "domain_cfg_out" ! newly created domain configuration filename
     !
  ln_use_jattr = .false.  !  use (T) the file attribute: open_ocean_jstart, if present
  !                       !  in netcdf input files, as the start j-row for reading
  /



  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .false.         !  = Use TEOS-10 equation of state
     ln_eos80    = .true.         !  = Use EOS80 equation of state
     ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)
                                   !
     !                     ! S-EOS coefficients (ln_seos=T):
     !                             !  rd(T,S,Z)*rau0 = -a0*(1+.5*lambda*dT+mu*Z+nu*dS)*dT+b0*dS
     rn_a0       =  1.6550e-1      !  thermal expension coefficient
     rn_b0       =  7.6554e-1      !  saline  expension coefficient
     rn_lambda1  =  5.9520e-2      !  cabbeling coeff in T^2  (=0 for linear eos)
     rn_lambda2  =  7.4914e-4      !  cabbeling coeff in S^2  (=0 for linear eos)
     rn_mu1      =  1.4970e-4      !  thermobaric coeff. in T (=0 for linear eos)
     rn_mu2      =  1.1090e-5      !  thermobaric coeff. in S (=0 for linear eos)
     rn_nu       =  2.4341e-3      !  cabbeling coeff in T*S  (=0 for linear eos)
  /

Edit namelist_ref::

  vi namelist_ref
  ...
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .false.         !  = Use TEOS-10 equation of state
     ln_eos80    = .false.         !  = Use EOS80 equation of state
     ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)
                                   !
     !                     ! S-EOS coefficients (ln_seos=T):
     !                             !  rd(T,S,Z)*rau0 = -a0*(1+.5*lambda*dT+mu*Z+nu*dS)*dT+b0*dS
     rn_a0       =  1.6550e-1      !  thermal expension coefficient
     rn_b0       =  7.6554e-1      !  saline  expension coefficient
     rn_lambda1  =  5.9520e-2      !  cabbeling coeff in T^2  (=0 for linear eos)
     rn_lambda2  =  7.4914e-4      !  cabbeling coeff in S^2  (=0 for linear eos)
     rn_mu1      =  1.4970e-4      !  thermobaric coeff. in T (=0 for linear eos)
     rn_mu2      =  1.1090e-5      !  thermobaric coeff. in S (=0 for linear eos)
     rn_nu       =  2.4341e-3      !  cabbeling coeff in T*S  (=0 for linear eos)
  /


SWitch in new namrun::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     nn_no       =       0   !  job number (no more used...)
     cn_exp      =    "AMM60_v4"  !  experience name
     nn_it000    =      1   !  first time step
     nn_itend    =      1440   !  last  time step (std 5475)
     nn_date0    =  20010101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       0   !  Leap year calendar (1) or not (0)
     ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "restart"   !  suffix of ocean restart name (input)
        cn_ocerst_indir = "."         !  directory from which to read input ocean restarts
        cn_ocerst_out   = "restart"   !  suffix of ocean restart name (output)
        cn_ocerst_outdir= "."         !  directory in which to write output ocean restarts
     ln_iscpl    = .false.   !  cavity evolution forcing or coupling to ice sheet model
     nn_istate   =       0   !  output the initial state (1) or not (0)
     ln_rst_list = .false.   !  output restarts at list of times using nn_stocklist (T) or at set frequency with nn_stock (F)
     nn_stock    =    1440   !  frequency of creation of a restart file (modulo referenced to 1)
     nn_stocklist = 0,0,0,0,0,0,0,0,0,0 ! List of timesteps when a restart file is to be written
     nn_write    =    1440   !  frequency of write in the output file   (modulo referenced to nn_it000)
     ln_mskland  = .false.   !  mask land points in NetCDF outputs (costly: + ~15%)
     ln_cfmeta   = .false.   !  output additional data to netCDF files required for compliance with the CF metadata standard
     ln_clobber  = .true.    !  clobber (overwrite) an existing file
     nn_chunksz  =       0   !  chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)
  /

Remove namzgr and namzgr_sco

New namdom::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     ln_linssh   = .false.   !  =T  linear free surface  ==>>  model level are fixed in time
     nn_closea   =    0      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
     !
     nn_msh      =    0      !  create (>0) a mesh file or not (=0)
     rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice
     !
     rn_rdt      =  360.     !  time step for the dynamics (and tracer if nn_acc=0)
     rn_atfp     =    0.1    !  asselin time filter parameter
     !
     ln_crs      = .false.   !  Logical switch for coarsening module
  /


Change namcrs::

  !-----------------------------------------------------------------------
  &namcrs        !   coarsened grid (for outputs and/or TOP)              ("key_crs")
  !-----------------------------------------------------------------------
     nn_factx    = 3         !  Reduction factor of x-direction
     nn_facty    = 3         !  Reduction factor of y-direction
     nn_binref   = 0         !  Bin centering preference: NORTH or EQUAT
                             !  0, coarse grid is binned with preferential treatment of the north fold
                             !  1, coarse grid is binned with centering at the equator
                             !    Symmetry with nn_facty being odd-numbered. Asymmetry with even-numbered nn_facty.
     nn_msh_crs  = 1         !  create (=1) a mesh file or not (=0)
     nn_crs_kz   = 0         ! 0, MEAN of volume boxes
                             ! 1, MAX of boxes
                             ! 2, MIN of boxes
     ln_crs_wn   = .true.    ! wn coarsened (T) or computed using horizontal divergence ( F )
  /


Add WnD::

  !-----------------------------------------------------------------------
  &namwad        !   Wetting and drying                                   (default F)
  !-----------------------------------------------------------------------
     ln_wd       = .false.   !  T/F activation of wetting and drying
     rn_wdmin1   =  0.1      !  Minimum wet depth on dried cells
     rn_wdmin2   =  0.01     !  Tolerance of min wet depth on dried cells
     rn_wdld     =  20.0     !  Land elevation below which wetting/drying is allowed
     nn_wdit     =  10       !  Max iterations for W/D limiter
  /

Add in 1D config options::

  !-----------------------------------------------------------------------
  &namc1d        !   1D configuration options                             ("key_c1d")
  !-----------------------------------------------------------------------
     rn_lat1d    =      50   !  Column latitude (default at PAPA station)
     rn_lon1d    =    -145   !  Column longitude (default at PAPA station)
     ln_c1d_locpt=  .true.   ! Localization of 1D config in a grid (T) or independant point (F)
  /
  !-----------------------------------------------------------------------
  &namc1d_dyndmp !   U & V newtonian damping                              ("key_c1d")
  !-----------------------------------------------------------------------
     ln_dyndmp   =  .false.  !  add a damping term (T) or not (F)
  /
  !-----------------------------------------------------------------------
  &namc1d_uvd    !   data: U & V currents                                 ("key_c1d")
  !-----------------------------------------------------------------------
  !              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_ucur     = 'ucurrent'  ,         -1        ,'u_current',   .false.    , .true. , 'monthly' ,  ''      ,  'Ume'   , ''
     sn_vcur     = 'vcurrent'  ,         -1        ,'v_current',   .false.    , .true. , 'monthly' ,  ''      ,  'Vme'   , ''
  !
     cn_dir        = './'    !  root directory for the location of the files
     ln_uvd_init   = .false. !  Initialisation of ocean U & V with U & V input data (T) or not (F)
     ln_uvd_dyndmp = .false. !  damping of ocean U & V toward U & V input data (T) or not (F)
  /

  !!======================================================================
  !!            ***  Surface Boundary Condition namelists  ***
  !!======================================================================
  !!   namsbc          surface boundary condition
  !!   namsbc_flx      flux               formulation                     (ln_flx     =T)
  !!   namsbc_blk      Bulk formulae formulation                          (ln_blk     =T)
  !!   namsbc_cpl      CouPLed            formulation                     ("key_oasis3" )
  !!   namsbc_sas      Stand-Alone Surface module
  !!   namtra_qsr      penetrative solar radiation                        (ln_traqsr  =T)
  !!   namsbc_rnf      river runoffs                                      (ln_rnf     =T)
  !!   namsbc_isf      ice shelf melting/freezing                         (nn_isf     >0)
  !!   namsbc_iscpl    coupling option between land ice model and ocean
  !!   namsbc_apr      Atmospheric Pressure                               (ln_apr_dyn =T)
  !!   namsbc_ssr      sea surface restoring term (for T and/or S)        (ln_ssr     =T)
  !!   namsbc_alb      albedo parameters
  !!   namsbc_wave     external fields from wave model                    (ln_wave    =T)
  !!   namberg         iceberg floats                                     (ln_icebergs=T)
  !!======================================================================
  !

There is a few changes with teh structure to how the Surface Boundry Condition is applied

Need blk apr and rnf

Remove ana, clio

core --> namsbc_blk
Add in land/sea mask column

Missing two variables::

sn_slp      = 'slp.15JUNE2009_fill'        ,         6         , 'SLP'     ,   .false.    , .true. , 'yearly'  , 'weights_core_orca2_bilinear_noc.nc'  , ''       , ''
   sn_tdif     = 'taudif_core'                ,        24         , 'taudif'  ,   .false.    , .true. , 'yearly'  , 'weights_core_orca2_bilinear_noc.nc'  , ''       , ''

Added in slp, below::

  sn_slp      = 'ggas'             ,       3          , 'MSL'     ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''

Updated T,Q, reference height to 2m (from 10m). Choose COARE_3p0

Add extra variables::

  !-----------------------------------------------------------------------
  &namsbc_blk   !   namsbc_blk  generic Bulk formula                      (ln_blk = T)
  !-----------------------------------------------------------------------
  !              !  file name                    ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                               !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_wndi     = 'ggas'             ,       6          , 'U10'     ,   .true.     , .false. , 'yearly'  ,'weights_bicub.nc' , 'Uwnd'      , ''
     sn_wndj     = 'ggas'             ,       6          , 'V10'     ,   .true.     , .false. , 'yearly'  ,'weights_bicub.nc' , 'Vwnd'      , ''
     sn_qsr      = 'gafs'             ,       3          , 'SSRD'    ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     sn_qlw      = 'gafs'             ,       3          , 'STRD'    ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     sn_tair     = 'ggas'             ,       6          , 'T2'      ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     sn_humi     = 'ggas'             ,       6          , 'Q'      ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''            , ''
     sn_prec     = 'gafs'             ,       3          , 'TP'      ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     sn_snow     = 'gafs'             ,       3          , 'SF'      ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     sn_slp      = 'ggas'             ,       3          , 'MSL'     ,   .true.     , .false. , 'yearly'  ,'weights_bilin.nc', ''           , ''
     !                    !  bulk algorithm :
     ln_NCAR     = .false.   ! "NCAR"      algorithm   (Large and Yeager 2008)
     ln_COARE_3p0= .true.   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
     ln_COARE_3p5= .false.   ! "COARE 3.5" algorithm   (Edson et al. 2013)
     ln_ECMWF    = .false.   ! "ECMWF"     algorithm   (IFS cycle 31)
     !
     cn_dir      = '/work/n01/n01/kariho40/NEMO/FORCINGS/ATM/ERAint/'      !  root directory for the location of the bulk files
     ln_taudif   = .false.   !  HF tau contribution: use "mean of stress module - module of the mean stress" data
     rn_zqt      = 2.       !  Air temperature and humidity reference height (m)
     rn_zu       = 10.       !  Wind vector reference height (m)
     rn_pfac     = 1.        !  multiplicative factor for precipitation (total & snow)
     rn_efac     = 1.        !  multiplicative factor for evaporation (0. or 1.)
     rn_vfac     = 0.        !  multiplicative factor for ocean/ice velocity
                             !  in the calculation of the wind stress (0.=absolute winds or 1.=relative winds)
     ln_Cd_L12   = .false.   !  Modify the drag ice-atm and oce-atm depending on ice concentration
                             !  This parameterization is from Lupkes et al. (JGR 2012)
  /

Delelte namsbc_mfs

Try leaving out::

  !-----------------------------------------------------------------------
  &namsbc_sas    !   Stand Alone Surface boundary condition
  !-----------------------------------------------------------------------


Add empty land/sea mask column to run off VARIABLES::

  !-----------------------------------------------------------------------
  &namsbc_rnf    !   runoffs namelist surface boundary condition
  !-----------------------------------------------------------------------
  !              !  file name           ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                      !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_rnf      = 'rivers'   ,        24         , 'sorunoff',   .true.    , .true. , 'yearly'  , ''       , ''                  , ''
     sn_cnf      = 'rivers'   ,        0          , 'socoefr',   .true.    , .true. , 'yearly'  , ''       , ''                   , ''
 ...


Add extra runoff options (not used)::


   ln_rnf_depth_ini = .false. ! compute depth at initialisation from runoff file
   rn_rnf_max  = 5.735e-4  !  max value of the runoff climatologie over global domain ( ln_rnf_depth_ini = .true )
   rn_dep_max  = 150.      !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )
   nn_rnf_depth_file = 0   !  create (=1) a runoff depth file or not (=0)


Edit nambdy::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  3                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'none'               !
      nn_dyn3d_dta  =  1                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  1                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_ice_lim      =  'none'             !
      nn_ice_lim_dta  =  0                  !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      rn_ice_tem      = 270.                !  lim3 only: arbitrary temperature of incoming sea ice
      rn_ice_sal      = 10.                 !  lim3 only:      --   salinity           --
      rn_ice_age      = 30.                 !  lim3 only:      --   age                --

      ln_tra_dmp    =.false.                !  open boudaries conditions for tracers
      ln_dyn3d_dmp  =.false.                !  open boundary condition for baroclinic velocities
      rn_time_dmp   =  1.                   ! Damping time scale in days
      rn_time_dmp_out =  1.                 ! Outflow damping time scale
      nn_rimwidth   = 10                    !  width of the relaxation zone
      ln_vol        = .false.               !  total volume correction (see nn_volctl parameter)
      nn_volctl     = 1                     !  = 0, the total water flux across open boundaries is zero
      nb_jpk_bdy    = 75                    ! number of levels in the bdy data (set < 0 if consistent with planned run)
  /

Though not sure about ln_tra_dmp, which was a sponge in AMM60::
..
  ln_sponge     = .true.                 ! Sponge added by enda
  rn_sponge     = 10                     ! Sponge diffusion multiplier


Add ln_tide::

    !-----------------------------------------------------------------------
    &nam_tide      !   tide parameters
    !-----------------------------------------------------------------------
       ln_tide     = .true.
       ln_tide_pot = .true.    !  use tidal potential forcing
       ln_tide_ramp= .false.   !

Copy appropriate namelist_ref::

  cp  $EXP/namelist_ref namelist_ref



Copy necessary file across::

  rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc .


Add in restarts::

  cd EXP_v4/RESTART
  ln -s /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXPD376/RESTART/01264320 .






#Edit run_counter.txt::
#
#  vi run_counter.txt
#  1 1 7200 20100105
#  2 1264321 1265761

Edit submit script::

  vi submit_nemo.pbs
  ...
  #PBS -N AMM60_v4
  #PBS -l select=92
  #PBS -l walltime=00:20:00

..
  Edit runscript (EXP dir + XIOS version)::

    vi run_nemo
    ...
    export RUNNAME=EXP_v4
    ...
    export XIOSDIR=/work/n01/n01/jelt/xios-2.0_r1080/


  Submit and see what it throws up::

    . run_nemo

::

  qsub submit_nemo.pbs



Issue in nambdy
cd  /work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4
less namelist_cfg


REference against
/work/n01/n01/jelt/ACCORD/trunk_NEMOGCM_r8395/CONFIG/ACCORD/EXP_EAFRICA> vi namelist_cfg






Hmmm not sure the restart WORKS



Try hacking submit_nemo.pbs  ::

  aprun -b -n $NEMOproc -N 24 ./nemo.exe : -N 5 -n $XIOSproc ./xios_server.exe >&stdouterr

Other bits::


  cd /work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4
  cp /work/n01/n01/jelt/AMM60/trunk_NEMOGCM_r8395/CONFIG/AMM60/EXP00/*xml .

  qsub submit_nemo.pbs


Missing init.nc file for cold start

Turn on restart flag and get nn_date0 from restart (it didn't work taking it
 from the restart).
Symbolic link restarts into EXPeriment directory::

  ln -s RESTART/00*/restart*nc .


Missing coordinates.bdy.nc file
Check what other files might be missing. Less run_nemo. (NB the
year might not be consistent with the restart - 20010101. Note also that the xml
files already copied in)::

  export YEARrun='2012'
  export GRIDDIR=/work/n01/n01/kariho40/NEMO/GRID         # Where to get forcings
  export DATADIR=/work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013           # Where to get forcings
  export WDIR=/work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4                        # Working directory


  #===============================================================
  # INPUT FILES
  #===============================================================
  #---------------------------------------------------------------
  # Coordinates
  #---------------------------------------------------------------
  echo `date`: Link coordinates
  ln -s $GRIDDIR/coordinates_AMM60.nc        ./coordinates.nc

  #---------------------------------------------------------------
  # Bathymetry
  #---------------------------------------------------------------
  echo `date`: Link Bathymetry
  ln -s $GRIDDIR/bathyfile_AMM60_nosmooth.nc ./bathy_meter.nc


  #===============================================================
  # INPUT FILES
  #===============================================================
  #---------------------------------------------------------------
  # BDY
  #---------------------------------------------------------------
  echo `date`: Link bdy data
  TIDEDIR=$DATADIR
  BDYDIR=/work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013

  for fic in `ls $TIDEDIR/NNA_AMM60bdy__bdytide*nc`; do
      ficdest=`basename $fic`
      ln -s $fic $WDIR/$ficdest
  done

  for yyyy in $YEARrun; do
    for fic in `ls $BDYDIR/AMM60bdy_NNA_R12_*${yyyy}*nc`; do
      ficdest=`basename $fic`
      ln -s $fic $WDIR/$ficdest
    done
  done

  ln -s $BDYDIR/coordinates.bdy.nc ./coordinates.bdy.nc
  #ln -s $DATADIR/runoff_AMM60_allindex_bathynosmooth.nc ./rivers.nc


Changed tracer ldf parameters
::

  !-----------------------------------------------------------------------
  &namtra_ldf    !   lateral diffusion scheme for tracers                 (default: NO diffusion)
  !-----------------------------------------------------------------------
     !                       !  Operator type:
     !                           !  no diffusion: set ln_traldf_lap=..._blp=F
     ln_traldf_lap   =  .false.  !    laplacian operator
     ln_traldf_blp   =  .false.  !  bilaplacian operator
     !
     !                       !  Direction of action:
     ln_traldf_lev   =  .false.  !  iso-level
     ln_traldf_hor   =  .false.  !  horizontal (geopotential)
     ln_traldf_iso   =  .true.  !  iso-neutral (standard operator)
     ln_traldf_triad =  .false.  !  iso-neutral (triad    operator)
     !
     !                             !  iso-neutral options:
     ln_traldf_msc   =  .false.  !  Method of Stabilizing Correction (both operators)
     rn_slpmax       =   0.01    !  slope limit                      (both operators)
     ln_triad_iso    =  .false.  !  pure horizontal mixing in ML              (triad only)
     rn_sw_triad     =  1        !  =1 switching triad ; =0 all 4 triads used (triad only)
     ln_botmix_triad =  .false.  !  lateral mixing on bottom                  (triad only)
     !
     !                       !  Coefficients:
     nn_aht_ijk_t    = 0         !  space/time variation of eddy coef
     !                                !   =-20 (=-30)    read in eddy_diffusivity_2D.nc (..._3D.nc) file
     !                                !   =  0           constant
     !                                !   = 10 F(k)      =ldf_c1d
     !                                !   = 20 F(i,j)    =ldf_c2d
     !                                !   = 21 F(i,j,t)  =Treguier et al. JPO 1997 formulation
     !                                !   = 30 F(i,j,k)  =ldf_c2d * ldf_c1d
     !                                !   = 31 F(i,j,k,t)=F(local velocity and grid-spacing)
     rn_aht_0        = 125.     !  lateral eddy diffusivity   (lap. operator) [m2/s]
     rn_bht_0        = 1.e+12    !  lateral eddy diffusivity (bilap. operator) [m4/s]
  /



  !-----------------------------------------------------------------------
  &namdyn_ldf    !   lateral diffusion on momentum
  !-----------------------------------------------------------------------
     !                       !  Type of the operator :
     !                           !  no diffusion: set ln_dynldf_lap=..._blp=F
     ln_dynldf_lap =  .true.    !    laplacian operator
     ln_dynldf_blp =  .false.    !  bilaplacian operator
     !                       !  Direction of action  :
     ln_dynldf_lev =  .false.    !  iso-level
     ln_dynldf_hor =  .true.    !  horizontal (geopotential)
     ln_dynldf_iso =  .false.    !  iso-neutral
     !                       !  Coefficient
     nn_ahm_ijk_t  = 0           !  space/time variation of eddy coef
     !                                !  =-30  read in eddy_viscosity_3D.nc file
     !                                !  =-20  read in eddy_viscosity_2D.nc file
     !                                !  =  0  constant
     !                                !  = 10  F(k)=c1d
     !                                !  = 20  F(i,j)=F(grid spacing)=c2d
     !                                !  = 30  F(i,j,k)=c2d*c1d
     !                                !  = 31  F(i,j,k)=F(grid spacing and local velocity)
     !                                !  = 32  F(i,j,k)=F(local gridscale and deformation rate)
     ! Caution in 20 and 30 cases the coefficient have to be given for a 1 degree grid (~111km)
     rn_ahm_0      =  50.     !  horizontal laplacian eddy viscosity   [m2/s]
     rn_ahm_b      =      0.     !  background eddy viscosity for ldf_iso [m2/s]
     rn_bhm_0      = -1.e+10      !  horizontal bilaplacian eddy viscosity [m4/s]
     !                       !  Smagorinsky settings (nn_ahm_ijk_t  = 32) :
     rn_csmc       = 3.5         !  Smagorinsky constant of proportionality
     rn_minfac     = 1.0         !  multiplier of theorectical lower limit
     rn_maxfac     = 1.0         !  multiplier of theorectical upper limit
  /

NOTE THERE WAS SOME CONFUISON OVER THE LAP COEFF 50 or 4000?


Add tracer advection. AMM60 used TVD, but that is not an option now. Use "FCT",
which with the boxS (seen in SWPacific)
::

  !-----------------------------------------------------------------------
  &namtra_adv    !   advection scheme for tracer
  !-----------------------------------------------------------------------
     ln_traadv_cen = .false. !  2nd order centered scheme
        nn_cen_h   =  4            !  =2/4, horizontal 2nd order CEN / 4th order CEN
        nn_cen_v   =  4            !  =2/4, vertical   2nd order CEN / 4th order COMPACT
     ln_traadv_fct = .true. !  FCT scheme
        nn_fct_h   =  2            !  =2/4, horizontal 2nd / 4th order
        nn_fct_v   =  2            !  =2/4, vertical   2nd / COMPACT 4th order
        nn_fct_zts =  0            !  >=1,  2nd order FCT scheme with vertical sub-timestepping
        !                          !        (number of sub-timestep = nn_fct_zts)
     ln_traadv_mus = .false. !  MUSCL scheme
        ln_mus_ups = .false.       !  use upstream scheme near river mouths
     ln_traadv_ubs = .false. !  UBS scheme
        nn_ubs_v   =  2            !  =2  , vertical 2nd order FCT / COMPACT 4th order
     ln_traadv_qck = .false. !  QUICKEST scheme
  /


Add more vvl lines::

  !-----------------------------------------------------------------------
  &nam_vvl  !
  !-----------------------------------------------------------------------
     ln_vvl_zstar  = .true.           !  zstar vertical coordinate
     ln_vvl_ztilde = .false.          !  ztilde vertical coordinate: only high frequency variations
     ln_vvl_layer  = .false.          !  full layer vertical coordinate
     ln_vvl_ztilde_as_zstar = .false. !  ztilde vertical coordinate emulating zstar
     ln_vvl_zstar_at_eqtor  = .false. !  ztilde near the equator
     rn_ahe3       = 0.0e0            !  thickness diffusion coefficient
     rn_rst_e3t    = 30.e0            !  ztilde to zstar restoration timescale [days]
     rn_lf_cutoff  = 5.0e0            !  cutoff frequency for low-pass filter  [days]
     rn_zdef_max   = 0.9e0            !  maximum fractional e3t deformation
     ln_vvl_dbg    = .true.           !  debug prints    (T/F)
  /


Add dyn advection and dynamics vorticity::


  !-----------------------------------------------------------------------
  &namdyn_adv    !   formulation of the momentum advection
  !-----------------------------------------------------------------------
     ln_dynadv_vec = .true.  !  vector form (T) or flux form (F)
     nn_dynkeg     = 0       ! scheme for grad(KE): =0   C2  ;  =1   Hollingsworth correction
     ln_dynadv_cen2= .false. !  flux form - 2nd order centered scheme
     ln_dynadv_ubs = .false. !  flux form - 3rd order UBS      scheme
     ln_dynzad_zts = .false. !  Use (T) sub timestepping for vertical momentum advection
  /
  !-----------------------------------------------------------------------
  &namdyn_vor    !   option of physics/algorithm (not control by CPP keys)
  !-----------------------------------------------------------------------
     ln_dynvor_ene = .false. !  enstrophy conserving scheme
     ln_dynvor_ens = .false. !  energy conserving scheme
     ln_dynvor_mix = .false. !  mixed scheme
     ln_dynvor_een = .true. !  energy & enstrophy scheme
        nn_een_e3f = 1          ! e3f = masked averaging of e3t divided by 4 (=0) or by the sum of mask (=1)
     ln_dynvor_msk = .false. !  vorticity multiplied by fmask (=T) or not (=F) (all vorticity schemes)  ! PLEASE DO NOT ACTIVATE
  /



Unhelpful error: It looks like something to do with XIOS. Try Copying
XML files from SWPacific instead of a freshly build EXP directory (as above)::

  cp /work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/*xml
  /work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4/.


Read date from namelist::

  nn_date0    =  20100101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
  nn_time0    =       0   !  initial time of day in hhmm
  nn_leapy    =       0   !  Leap year calendar (1) or not (0)
  ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)


Resubmit.

Try again::

  qsub submit_nemo.PBS


Initial conditions from AMM60 run. Cold starts
++++++++++++++++++++++++++++++++++++++++++++++

Insert new method to use AMM60 data for initial conditions.
/work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT
AMM60_5d_20131013_20131129_grid_T.nc

Note that the temperature and salinity variables are ``thetao`` and ``so``

::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0
  cd $INPUTS

  ncks /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT/AMM60_5d_20120620_20120818_grid_T.nc $INPUTS/AMM60_ave_20120620_20120818_grid_T.nc
  cp $INPUTS/AMM60_ave_20120620_20120818_grid_T.nc $START_FILES/AMM60_ave_20120620_20120818_grid_T.nc

Average over time and restore the parallel modules (this went into a serial queue
in the script ``$INPUTS/ncwa_AMM60_5d``. Took 11 minutes.)::

  ncwa -a time_counter $START_FILES/AMM60_ave_20120620_20120818_grid_T.nc $INPUTS/AMM60_ave_20120620_20120818_grid_T.nc
  cp $INPUTS/AMM60_ave_20120620_20120818_grid_T.nc $START_FILES/AMM60_ave_20120620_20120818_grid_T.nc

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel


----

Cold start
++++++++++

Edit namelist::

  vi namelist_cfg
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "AMM60"  !  experience name
     nn_it000    =  1   !  first time step
     nn_itend    =  1440 ! 30day=7200   !  last  time step (std 5475)
     nn_date0    =  20120701   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       0   !  Leap year calendar (1) or not (0)
     ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)


init.nc --> AMM60_ave_20120620_20120818_grid_T.nc with path settings ::

  !-----------------------------------------------------------------------
  &namtsd    !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  ! file name ! frequency (hours)    ! variable ! time interp. ! clim  !'yearly' or ! weights  ! rotation !
  !           !  (if <0  months)     !   name   !  (logical)   ! (T/F) ! 'monthly'  ! filename ! pairing  !
     sn_tem  = 'AMM60_ave_20120620_20120818_grid_T.nc', -1,'thetao',  .true.  , .true., 'yearly'   , ' '      , ' '
     sn_sal  = 'AMM60_ave_20120620_20120818_grid_T.nc', -1,'so',  .true.  , .true., 'yearly'   , ''       , ' '
  !
     cn_dir        = '/work/n01/n01/jelt/AMM60/INPUTS/'     !  root directory for the location of the runoff files
     ln_tsd_init   = .true.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)
  /

---
Injected *(9 Feb 2018)*. To try and do cold starts from AMM60 T,S output. It works.
----

Dead end reading in baroclinic boundary conditions. Run terminates without an error::

  ...
  Number of levels for vobtcrtx is            1
          read vobtcrtx (rec:      1) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2010m01.nc ok
  fld_init : time-interpolation for vobtcrtx read previous record =      1 at time =   -0.50 days

  ===>>> : W A R N I N G
          ===============

  previous year/month/week/day file: /work/n01/n01/kariho40/NEMO/FORCINGS/2010_20
  13/AMM60bdy_NNA_R12_bdyV not present -> back to current year/month/week/day
                      iom_nf90_open ~~~ open existing file: /work/n01/n01/kariho4
  0/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyV_y2010m01.nc in READ mode
                     ---> /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy
  _NNA_R12_bdyV_y2010m01.nc OK
  Dim size for vobtcrty is         2642
  Number of levels for vobtcrty is            1
            read vobtcrty (rec:      1) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyV_y2010m01.nc ok
    fld_init : time-interpolation for vobtcrty read previous record =      1 at time =   -0.50 days

  ===>>> : W A R N I N G
          ===============

  previous year/month/week/day file: /work/n01/n01/kariho40/NEMO/FORCINGS/2010_20
  13/AMM60bdy_NNA_R12_bdyT not present -> back to current year/month/week/day
                      iom_nf90_open ~~~ open existing file: /work/n01/n01/kariho4
  0/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyT_y2010m01.nc in READ mode
                     ---> /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy
  _NNA_R12_bdyT_y2010m01.nc OK
  Dim size for votemper is        26342
  Number of levels for votemper is           51
            read votemper (rec:      1) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyT_y2010m01.nc ok


Try neumann boundary conditions (as James does in his ORCHESTRA run)::

  cn_dyn3d      =  'neumann'               !
  nn_dyn3d_dta  =  1

I thought that I had some trouble in fldread, but I think it was actually something to do with a missing met
variable...

.. old notes (now to remove)
  Make some edits in fldread. It looks like perhaps the new version does not read in bdy data with different numbers of levels to the target grid as well as it did.
  Try copying out a chunk of code from the old fldread.F90 and see what happens.
  Code modifcation are all in MY_SRC/fldread.F90
  Compile and resubmit::

    cd $CDIR
    ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

  cd $AMM60DIR/EXP_v4
  qsub submit_nemo.pbs


  **OPTIONS**
  **
  SWitch from false
  ln_full_vel = .true.
  This is belt and braces for my config. It removes the mean for the 3d vel.
  **

  Could try
  bdy_dta + write statements
  Or wrap statement in step.F90:  call bdy_dta

  Try -g compile option.


**PENDING**
Try
ln_bdy = .false.
**PENDING**
This works. It blows up with fast velocities at fisrt time step.


Try again with ln_bdy=TRUE
**PENDING**
Crashes::

  jelt:I think we got stuck here in fldread.F90. jpk_bdy          75
            read vozocrtx (rec:      1) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2010m01.nc ok
  jelt" pre interp


Note to James::

  Update before I jet off to AGU

  ln_bdy = F —> It ran to the first time step, when it blew up - which is kind of OK as the restart field was botched.
  ln_bdy = T —> It fell over within the sophisticated ERROR trapping:

  tail ocean.output
  ...
  previous year/month/week/day file: /work/n01/n01/kariho40/NEMO/FORCINGS/2010_20
   13/AMM60bdy_NNA_R12_bdyU not present -> back to current year/month/week/day
                       iom_nf90_open ~~~ open existing file: /work/n01/n01/kariho4
   0/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2010m01.nc in READ mode
                      ---> /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy
   _NNA_R12_bdyU_y2010m01.nc OK
   Dim size for vozocrtx is         2642
   Number of levels for vozocrtx is           51
   jelt:I think we got stuck here in fldread.F90. jpk_bdy          75
             read vozocrtx (rec:      1) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2010m01.nc ok
   jelt" pre interp


  Recall fldread.F90:

           WRITE(numout,*)'jelt" pre interp'
           IF ( ln_bdy ) THEN
             CALL fld_bdy_interp(dta_read, dta_read_z, dta_read_dz, map, jpk_bdy, igrd, ibdy, fv, dta, fvl, ilendta)
             CALL flush(numout)
             WRITE(numout,*)'jelt" post interp'
           ENDIF


Check to see if there is a problem with gdept not existing...
Inserted write statemtsn to see if arrays are empty.

Fails without a clear error. Adding MAXVAL write end of ocean.output::

  Number of levels for sossheig is            1
          read sossheig (rec:     30) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bt_bdyT_y2012m06.nc ok
  fld_init : time-interpolation for sossheig read previous record =     30 at time =  180.50 days
                    iom_nf90_open ~~~ open existing file: /work/n01/n01/kariho4
  0/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2012m06.nc in READ mode
                     ---> /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy
  _NNA_R12_bdyU_y2012m06.nc OK
  Dim size for vozocrtx is         2642
  Number of levels for vozocrtx is           51
  jelt:I think we got stuck here in fldread.F90. jpk_bdy          75
            read vozocrtx (rec:     30) in /work/n01/n01/kariho40/NEMO/FORCINGS/2010_2013/AMM60bdy_NNA_R12_bdyU_y2012m06.nc ok
  jelt: pre interp. Maxvel dta_read_z:   0.000000000000000E+000
   Maxval dta_read_dz:   0.000000000000000E+000

Suggests that dz and z are not getting written. James spotted that this is because igrd=0 when
he expect igrd=2

**HERE**


Could debug on short queue? 2 XIOS servers. remaining 168 on opa (NEMOPROC?)
--> short_submit_nemo.pbs (1 XIOS, 96 NEMOPROC. Also edit namelist_cfg in nammpp)
Cold start, so I do't have to rebuild restart file.

**Short queue should work now**


I will look into fld_bdy_interp when I come back to this fun problem.

James has a fix in::

  /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/

Copy files across::

  export jdhaMY_SRC=/work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/AMM60/MY_SRC/
  cp $jdhaMY_SRC/bdy_oce.F90 $CDIR/$CONFIG/MY_SRC/bdy_oce.F90
  cp $jdhaMY_SRC/bdydta.F90 $CDIR/$CONFIG/MY_SRC/bdydta.F90

Recompile::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Resubmit::

  qsub short_submit_nemo.pbs

Core dump error. tail ocean.output::

  dyn_keg : kinetic energy gradient trend, scheme number=           0
  ~~~~~~~

  dyn_zad : arakawa advection scheme

  dyn:vor_een : vorticity term: energy and enstrophy conserving scheme


Karen's AMM60
/work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_harmIT2/WDIR/output.namelist.dyn
::

  &NAMDYN_VOR
  LN_DYNVOR_ENS   = F,
  LN_DYNVOR_ENE   = F,
  LN_DYNVOR_MIX   = F,
  LN_DYNVOR_EEN   = T,
  LN_DYNVOR_EEN_OLD       = F
  /

This run
/work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4/output.namelist.dyn
::

  &NAMDYN_VOR
  LN_DYNVOR_ENS   = F,
  LN_DYNVOR_ENE   = F,
  LN_DYNVOR_MIX   = F,
  LN_DYNVOR_EEN   = T,
  NN_EEN_E3F      =           1,
  LN_DYNVOR_MSK   = F
/


SWPacific
/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/output.namelist.dyn
::

  &NAMDYN_VOR
  LN_DYNVOR_ENS   = F,
  LN_DYNVOR_ENE   = F,
  LN_DYNVOR_MIX   = F,
  LN_DYNVOR_EEN   = T,
  NN_EEN_E3F      =           1,
  LN_DYNVOR_MSK   = F
  /


Karen's tail ocean.output ::

  less /work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_harmIT2/LOGS/restart/ocean.output_EXP_harmIT2
  ...
  dyn_zad : arakawa advection scheme

  dyn:vor_een : vorticity term: energy and enstrophy conserving scheme
  ~~~~~~~~~~~

  dyn_ldf_iso : iso-neutral laplacian diffusive operator or
  ~~~~~~~~~~~   s-coordinate horizontal diffusive operator

  dyn:hpg_prj : hydrostatic pressure gradient trend
  ~~~~~~~~~~~   s-coordinate case, cubic spline pressure Jacobian
  ...



This run::

  tail -100 /work/n01/n01/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_v4/ocean.output

  ...

  dyn_zad : arakawa advection scheme

  dyn:vor_een : vorticity term: energy and enstrophy conserving scheme
  ~~~~~~~~~~~


So perhaps it is an issue with dyn_ldf_iso?

Comparing with SEAsia (LN_DYNLDF_BLP=T) instead of Karen's (LN_DYNLDF_LAP=T)
Try switching to BLP=T
::

  vi namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &namdyn_ldf    !   lateral diffusion on momentum
  !-----------------------------------------------------------------------
     !                       !  Type of the operator :
     !                           !  no diffusion: set ln_dynldf_lap=..._blp=F
     ln_dynldf_lap =  .false.    !    laplacian operator
     ln_dynldf_blp =  .true.    !  bilaplacian operator
     !                       !  Direction of action  :
     ln_dynldf_lev =  .true.    !  iso-level
     ln_dynldf_hor =  .false.    !  horizontal (geopotential)
     ln_dynldf_iso =  .false.    !  iso-neutral
     !                       !  Coefficient
     nn_ahm_ijk_t  = 0           !  space/time variation of eddy coef
     !                                !  =-30  read in eddy_viscosity_3D.nc file
     !                                !  =-20  read in eddy_viscosity_2D.nc file
     !                                !  =  0  constant
     !                                !  = 10  F(k)=c1d
     !                                !  = 20  F(i,j)=F(grid spacing)=c2d
     !                                !  = 30  F(i,j,k)=c2d*c1d
     !                                !  = 31  F(i,j,k)=F(grid spacing and local velocity)
     !                                !  = 32  F(i,j,k)=F(local gridscale and deformation rate)
     ! Caution in 20 and 30 cases the coefficient have to be given for a 1 degree grid (~111km)
     rn_ahm_0      =  0.     !  horizontal laplacian eddy viscosity   [m2/s]
     rn_ahm_b      =      0.     !  background eddy viscosity for ldf_iso [m2/s]
     rn_bhm_0      = -1.25e+10      !  horizontal bilaplacian eddy viscosity [m4/s]
     !                       !  Smagorinsky settings (nn_ahm_ijk_t  = 32) :
     rn_csmc       = 3.5         !  Smagorinsky constant of proportionality
     rn_minfac     = 1.0         !  multiplier of theorectical lower limit
     rn_maxfac     = 1.0         !  multiplier of theorectical upper limit
  /

Submit::

  qsub short_submit_nemo

Still get a core dump with the last ocean.output at::

  dyn:vor_een : vorticity term: energy and enstrophy conserving scheme



















































Next steps
++++++++++



Update tides code with Nico's version.
++++++++++++++++++++++++++++++++++++++


MPP decomposition for land suppression
++++++++++++++++++++++++++++++++++++++


Backup to repo key files
++++++++++++++++++++++++


Rebuild the output and inspect `rebuild_and_inspect_NEMO_output.rst`_
++++++++++++++++++++++++++++++

---
