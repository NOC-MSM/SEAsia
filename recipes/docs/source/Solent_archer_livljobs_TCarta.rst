==========================================
Setting up a Solent NEMO v4 configuration
==========================================

Machines: livljobs4, ARCHER

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a Solent regional model using TCarta bathymetry, FES tidal boundaries.

Build on a combination of livljobs4 and ARCHER.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, domain_cfg.nc
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO. Generate bathymetry
#. ARCHER: Run simulation


Issues that arose
=================

* SCRIP tools didn't work on ARCHER. Walltime too long for interactive mode,
and couldn't get serial job to work with intel compiler.
* ...

----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO and SCRIP bathymetry
generation bits are done on livljobs4.
Starting on ARCHER::

  ssh login.archer.ac.uk

  export CONFIG=Solent
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

---

Collect essential files
=======================

Note you might have to mkdir the odd directory or two... Some of these will need
to be done on ARCHER and LIVLJOBS4::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES

  # SCRIP tools
  rsync -uvt $WDIR/../SEAsia/START_FILES/scrip.patch $START_FILES/.
  rsync -uvt $WDIR/../SEAsia/START_FILES/scripgrid.patch $START_FILES/.
  rsync -uvt $WDIR/../SEAsia/START_FILES/scripinterp.patch $START_FILES/.
  rsync -uvt $WDIR/../SEAsia/START_FILES/scripinterp_mod.patch $START_FILES/.
  rsync -uvt $WDIR/../SEAsia/START_FILES/scripshape.patch $START_FILES/.

  rsync -uvt $WDIR/../LBay/START_FILES/dommsk.F90 $START_FILES/.
  rsync -uvt $WDIR/../LBay/START_FILES/bdyini.F90 $START_FILES/.
  rsync -uvt $WDIR/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  rsync -uvt $WDIR/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  rsync -uvt $WDIR/../SWPacific/START_FILES/usrdef_istate.F90 $START_FILES/.
  rsync -uvt $WDIR/../SWPacific/START_FILES/usrdef_sbc.F90    $START_FILES/.

  # Interpolate z-coordinates for initial conditions on the fly
  rsync -uvt /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/par_oce.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/dtatsd.F90 $START_FILES/.

  # FES tides and POLCOMS harmonic analysis `<FES2014_NEMO.rst>`_
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/bdyini.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm_fast.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/sbctide.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step_oce.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_FES14.h90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tideini.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_mod.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/ $START_FILES/.


Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build (if it is already downloaded). Note here we use user defined
 functions for the initial state (constant T and S) and surface forcing (zero forcing)::

  cd $CDIR
  cp $START_FILES/usrdef_istate.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/usrdef_sbc.F90    $CDIR/$CONFIG/MY_SRC/.
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

---

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

Or just link XIOS executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---


Build TOOLS
===========

To generate domain coords and rebuild tools we first need
to compile some of the NEMO TOOLS. These are split across ARCHER and LIVLJOBS4.
In particular the SCRIP tools didn't work on ARCHER.

.. note: These are compiled with XIOS2. However DOMAINcfg has to be compiled
  with XIOS1. There is a README in the $TDIR/DOMAINcfg on what to do.

First build DOMAINcfg (which is relatively new and in NEMOv4). Use my XIOS1 file
(see userid and path in variable ``%XIOS_HOME``). Copy from ARCH *store*::

  # ARCHER
  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $TDIR/../ARCH/.
  cd $TDIR

  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO


---

To generate bathymetry, initial conditions and grid information we first need
to compile some of the NEMO TOOLS (after a small bugfix - and to allow direct
passing of arguments). **For some reason GRIDGEN doesn’t like INTEL.**
**Do this on livljobs4**::

  ssh livljobs4

Apply patches::

  cd $TDIR/WEIGHTS/src
  patch -b < $START_FILES/scripinterp_mod.patch
  patch -b < $START_FILES/scripinterp.patch
  patch -b < $START_FILES/scrip.patch
  patch -b < $START_FILES/scripshape.patch
  patch -b < $START_FILES/scripgrid.patch

Setup for PGI modules and compile::

  cd $TDIR
  cp $START_FILES/arch-pgf90_linux_jb.fcm $TDIR/../ARCH/arch-pgf90_linux_jb.fcm

  module add netcdf/gcc/4.1.3
  module add pgi/15.4

  ./maketools -n WEIGHTS -m pgf90_linux_jb
  ./maketools -n REBUILD_NEMO -m pgf90_linux_jb
  ./maketools -n GRIDGEN -m pgf90_linux_jb
  ./maketools -n NESTING -m pgf90_linux_jb

.. note: Originally I tried to do this on ARCHER but hit a limit with the walltime
 for the interactive mode. Then I couldn't get it to work on the serial queue. So
 just ported to livljobs4. Don't necessarily need to compile all these on livljobs


1. Generate new coordinates file
++++++++++++++++++++++++++++++++

...

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


2. Generate bathymetry file
+++++++++++++++++++++++++++

.. note:
  This was first done on livljobs4. Here I do it on ARCHER.

Take parent bathymetry from TCarta. Put the data on ARCHER in $INPUTS::

 scp TCarta_30m_dbm.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/SEAsia/INPUTS/.

and::

  cp $START_FILES/namelist_reshape_bilin_tcarta $INPUTS/.


Do some things to 0) set the elevation variable name, 1) flatten out land
 elevations, 2) make depths positive.
Have to swap around with the modules to get nco working *(James
noted a problem with the default nco module)*
Here force a minimum depth of 10m by lowering the data by 12m::

  cd $INPUTS

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncrename -v TCarta_30m_dbm.tif,elevation   TCarta_30m_dbm.nc
  ncap2 -s 'elevation=elevation-12' TCarta_30m_dbm.nc tmp.nc

  # Fix the missing values: reassign the values and zero. (Check what the values are first)
  ncap2 -s 'where(elevation < -32000) elevation=0' tmp.nc tmp2.nc

  ncflint --fix_rec_crd -w -1.0,0.0 tmp2.nc tmp2.nc TCarta_Solent_10m.nc
  rm tmp*.nc

Restore the original parallel modules::

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel


Port (TCarta_Solent_10m.nc, coordinates.nc) to livljobs4 for processing
 (convert to v3 and run SCRIP routines)::

  # Convert to netCDF3_classic
  module load nco/gcc/4.4.2.ncwa
  ncks --fl_fmt=classic TCarta_Solent_10m.nc TCarta_Solent_10m_v3.nc


Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h TCarta_Solent_10m_v3.nc`` to get input variable names)::

  vi $INPUTS/namelist_reshape_bilin_tcarta
  ...
  &grid_inputs
    input_file = 'TCarta_Solent_10m_v3.nc'
    nemo_file = 'coordinates.nc'
    ...
    input_lon = 'lon'
    input_lat = 'lat'
    nemo_lon = 'glamt'
    nemo_lat = 'gphit'
    ...

    &interp_inputs
    input_file = "TCarta_Solent_10m_v3.nc"
    ...
    input_name = "elevation"


Restore the original parallel modules (used livljobs4)::

  # livljobs4
  module purge
  module add netcdf/gcc/4.1.3
  module add pgi/15.4

Execute first scrip thing::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_tcarta

Output files::

  remap_nemo_grid_tcarta.nc
  remap_data_grid_tcarta.nc

Execute second scip thing (~6 hrs)::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_tcarta

Output files::

  data_nemo_bilin_tcarta.nc

Execute third scip thing::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_tcarta

Output files::

  bathy_meter.nc

Copy to Archer::

  rsync -urtv bathy_meter.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.

3. Generate initial conditions
++++++++++++++++++++++++++++++

Skip this for a tide-only run. Using user defined constant T and S.

For constant T and S use the user defined functions in ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``.



4. Generate a domain configuration file
=======================================

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg``
directory along with all the inputs files that would have previously been needed
get v3.6 running. The reason being that all the non-time stepping stuff, like
grid generating, has been abstracted from the core OPA code and is now done as
a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory::

    #ln -s $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.  # Already there
    ln -s $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template ``namelist_cfg`` with only the essenetial domain building stuff.
Get the size of the new domain from ``ncdump -h bathy_meter.nc``.

Follow recipe of hybrid z-s coordinates in `build_domain_cfg_file.rst`_

.. note : couldn't compile make_domain_cfg.exe (don't know what...) so linked to SEAsia executable...

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.


5. Generate tidal boundary conditions
=====================================

livljobs4: get all the necessary files onto this machine::

  cd $INPUTS
  rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/domain_cfg.nc .
  rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/coordinates.nc .
  rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/jelt/Solent/INPUTS/bathy_meter.nc .

Need to generate 3 more files: A ``namelist.bdy`` which drives PyNEMO and which
has two input files: ``inputs_src.ncml`` which points to the data source and
``inputs_dst.ncml`` which remaps some variable names in the destination files::

  cp ../../SWPacific/INPUTS/namelist.bdy .


First install PyNEMO `install_nrct`_ if not already done so. Use branch ``Generalise-tide-input``::

  cd /work/$USER/nrct
  git checkout Generalise-tide-input

Copy across some parent mesh files and a mask file (even though they are not
used. This is because this old version of PyNEMO didn't anticipate tide-only usage)::

  cp ../../SEAsia/INPUTS/mesh_?gr_src.nc $INPUTS/.
  cp ../../SEAsia/mask_src.nc $INPUTS/.


If I don't make a boundary mask then it doesn't work... This can also be done with
the PyNEMO GUI. The mask variable takes values (-1 mask, 1 wet, 0 land). Get a
template from domain_cfg.nc and then modify as desired around the boundary::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython::

  import netCDF4
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  dset.variables['mask'][0,:]  = -1     # Southern boundary
  dset.variables['mask'][-1,:] = -1    # Northern boundary
  dset.variables['mask'][:,-1] = -1    # Eastern boundary
  dset.variables['mask'][:,0] = -1        # Western boundary
  dset.close()



Couldn't find the FES data (they have moved from Tom's work dir). Tide data source
is clumsily set in ``nemo_bdy_tide3.py``

::
  vi namelist.bdy

  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  : namelist for BDY generation tool
  !!
  !!             User inputs for generating open boundary conditions
  !!             employed by the BDY module in NEMO. Boundary data
  !!             can be set up for v3.2 NEMO and above.
  !!
  !!             More info here.....
  !!
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
     sn_src_dir = './cut_inputs_src.ncml'       ! src_files/'
     sn_dst_dir = '/work/jelt/NEMO/Solent/INPUTS/'
     sn_fn      = 'Solent'                 ! prefix for output files
     nn_fv      = -1e20                     !  set fill value for output files
     nn_src_time_adj = 0                                    ! src time adjustment
     sn_dst_metainfo = 'metadata info: jelt'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_coords_file = .false.               !  =T : produce bdy coordinates files
      cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
      ln_mask_file   = .true.              !  =T : read mask from file
      cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
      ln_tra         = .false.               !  boundary conditions for T and S
      ln_ice         = .false.               !  ice boundary condition
      nn_rimwidth    = 1                    !  width of the relaxation zone

  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
      ln_tide        = .true.               !  =T : produce bdy tidal conditions
  !TESTING
  !               clname(1) ='M2'
  !               clname(2)='S2'
  !               clname(3)='K2'
  !TPXO
      clname(1)='m2'
      clname(2)='s2'
      clname(3)='n2'
      clname(4)='k2'
      clname(5)='k1'
      clname(6)='o1'
      clname(7)='p1'
      clname(8)='q1'
      clname(9)='mf'
      clname(10)='mm'
      clname(11)='m4'
      clname(12)='ms4'
      clname(13)='mn4'
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



Generate the boundary conditions with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

This creates::

  Solent_bdytide_rotT_M4_grid_T.nc
  Solent_bdytide_rotT_MM_grid_T.nc
  Solent_bdytide_rotT_MN4_grid_T.nc
  Solent_bdytide_rotT_MS4_grid_T.nc
  Solent_bdytide_rotT_M2_grid_T.nc
  Solent_bdytide_rotT_N2_grid_T.nc
  Solent_bdytide_rotT_S2_grid_T.nc
  Solent_bdytide_rotT_K1_grid_T.nc
  Solent_bdytide_rotT_K2_grid_T.nc
  Solent_bdytide_rotT_P1_grid_T.nc
  Solent_bdytide_rotT_O1_grid_T.nc
  Solent_bdytide_rotT_MF_grid_T.nc
  Solent_bdytide_rotT_Q1_grid_T.nc
  Solent_bdytide_rotT_M4_grid_U.nc
  Solent_bdytide_rotT_MM_grid_U.nc
  Solent_bdytide_rotT_MN4_grid_U.nc
  Solent_bdytide_rotT_MS4_grid_U.nc
  Solent_bdytide_rotT_M2_grid_U.nc
  Solent_bdytide_rotT_N2_grid_U.nc
  Solent_bdytide_rotT_S2_grid_U.nc
  Solent_bdytide_rotT_K1_grid_U.nc
  Solent_bdytide_rotT_K2_grid_U.nc
  Solent_bdytide_rotT_P1_grid_U.nc
  Solent_bdytide_rotT_O1_grid_U.nc
  Solent_bdytide_rotT_MF_grid_U.nc
  Solent_bdytide_rotT_Q1_grid_U.nc
  Solent_bdytide_rotT_M4_grid_V.nc
  Solent_bdytide_rotT_MM_grid_V.nc
  Solent_bdytide_rotT_MN4_grid_V.nc
  Solent_bdytide_rotT_MS4_grid_V.nc
  Solent_bdytide_rotT_M2_grid_V.nc
  Solent_bdytide_rotT_N2_grid_V.nc
  Solent_bdytide_rotT_S2_grid_V.nc
  Solent_bdytide_rotT_K1_grid_V.nc
  Solent_bdytide_rotT_K2_grid_V.nc
  Solent_bdytide_rotT_P1_grid_V.nc
  Solent_bdytide_rotT_O1_grid_V.nc
  Solent_bdytide_rotT_MF_grid_V.nc
  Solent_bdytide_rotT_Q1_grid_V.nc
  coordinates.bdy.nc

Copy the new files back onto ARCHER::

  livljobs4$
  cd $INPUTS
  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.bdy.nc
  for file in $CONFIG*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done


6. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

Get set up::

  ssh archer
  . ~/temporary_path_names_for_NEMO_build

Get important files into EXP directory. Should already have ``domain_cfg.nc``::


  cd $EXP
  rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.
  rsync -tuv $START_FILES/namelist_cfg $EXP/.

Create symbolic links from EXP directory::

  ln -s $INPUTS $EXP/bdydta

Edit the output to have 1hrly SSH::

 vi file_def_nemo.xml
 ...
 <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
  <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
    <field field_ref="ssh"          name="zos"   />
  </file>
 </file_group>

---

Create a short queue runscript. (Note: PBS -N jobname, PBS -m email)::

  vi runscript

  #!/bin/bash
  #PBS -N Solent
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  #  echo $(readlink -f $PBS_O_WORKDIR)
  # export OMP_NUM_THREADS=1

  cd $PBS_O_WORKDIR
  #
    echo " ";
    OCEANCORES=74
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

  exit

Change the notification email to your own address::

  sed -i "s/xxx@noc/$USER@noc/g" runscript



Follow `<MMP_decomp_land_suppression.rst>`_ to set the number of processors etc.
Edit namelist_cfg key parameters::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
  !-----------------------------------------------------------------------
   ...
   jpni        =  12    !  jpni   number of processors following i (set automatically if < 1)
   jpnj        =  8     !  jpnj   number of processors following j (set automatically if < 1)
   jpnij       =  74    !  jpnij  number of local domains (set automatically if < 1)

  rn_rdt = 1.
  ln_bdy=.false.
  ln_tide=.true.
  ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
  ln_tsd_interp = .false.    !  Interpolation of T & S in the verticalinput data (T) or not (F)


Adjust runscript parameters::

  vi runscript
  ...
  #PBS -l select=5
  ...
  OCEANCORES=74
  ...
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 19 ./opa




27 May 2018
**BLOWS AFTER A FEW STEPS. PROBABLY SOMETHING WRONG WITH THE VERTICAL GRID IN
THE DEEP CHANNEL**

7 Jun 2018
Regenerate domain_cfg.nc with stretched sigma and resubmit.
Still blows after a few steps. Didn't get a chance to look why..















### SCRATCH
###########

export CONFIG=Solent
export WORK=/work
export WDIR=$WORK/jelt/NEMO/$CONFIG
export INPUTS=$WDIR/INPUTS
export START_FILES=$WDIR/START_FILES


## Make a Solent directory

livljobs4:
mkdir $WDIR
mkdir $INPUTS
mkdir $START_FILES


cd $START_FILES

## Ckeckout the TCarta data for the Solent. 1/3600 deg

import netCDF4
import matplotlib.pyplot as plt
dset = netCDF4.Dataset('TCarta_30m_dbm.nc','a')

H = dset.variables['TCarta_30m_dbm.tif'][:,:]
lon = dset.variables['lon'][:]
lat = dset.variables['lat'][:]

plt.pcolormesh(lon,lat,H)
plt.colorbar()
plt.show()

lat:
50.9166666666666714, 50.53833333333690803
lon:
-1.6638888888888899, 0.99472222222211859

## Find the coordinates for the Solent.


import netCDF4
import matplotlib.pyplot as plt
import numpy as np
dset = netCDF4.Dataset('coordinates_ORCA_R12.nc','a')

nav_lon = dset.variables['nav_lon'][:]
nav_lat = dset.variables['nav_lat'][:]

plt.pcolormesh(lon,lat,H)
plt.colorbar()
plt.show()


I = 2203,3411 # SW Corner.

In [46]: print 'lat {}, lon {}'.format(nav_lat[I], nav_lon[I])
lat 50.5043754578, lon -1.67069804668

---

livljobs6:


module load nco/gcc/4.4.2.ncwa
# create new file with correct dimensions within
ncks -O -C -a -v TCarta_30m_dbm.tif TCarta_30m_dbm.nc coordinates_Solent_R3600.nc

# rename variable and dimensions
ncrename -v TCarta_30m_dbm.tif,dummy   coordinates_Solent_R3600.nc
ncrename -d lat,y -d lon,x  coordinates_Solent_R3600.nc

# Add a depth dimension
ncap2 -s 'defdim("z",1); dum3[$z,$y,$x]=dummy' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -s 'defdim("time",1); dum4[$time,$z,$y,$x]=dum3' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

# add a blank variables
ncap2 -O -s 'nav_lat[y,x]=1f' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'nav_lon[y,x]=1f' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

ncap2 -O -s 'gphit[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'gphif[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'gphiu[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'gphiv[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

ncap2 -O -s 'glamt[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'glamf[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'glamu[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'glamv[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

ncap2 -O -s 'e1t[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e1u[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e1v[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e1f[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

ncap2 -O -s 'e2t[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e2u[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e2v[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc
ncap2 -O -s 'e2f[time,z,y,x]=1d' coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc


# add needed attributes
ncatted -a units,nav_lon,c,c,'degrees_east'      \
        -a units,nav_lon,c,c,'degrees_north'  coordinates_Solent_R3600.nc

# Remove now unnecessary bathymetry variable
ncks -x -v dummy,dum3,dum4 coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

# Extend the time dimension to UNLIMITED
ncks -O --mk_rec time coordinates_Solent_R3600.nc coordinates_Solent_R3600.nc

---
# Replace lat lon variable with bathy data
ipython

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
dout = netCDF4.Dataset('coordinates_Solent_R3600.nc','a')

din = netCDF4.Dataset('TCarta_30m_dbm.nc','a')

TC_lon = din.variables['lon'][:]
TC_lat = din.variables['lat'][:]

yv, xv = np.meshgrid(np.flipud(TC_lat), TC_lon)

dout.variables['glamt'][0,0,:,:] = xv.T
dout.variables['gphit'][0,0,:,:] = yv.T
dout.variables['nav_lon'][:,:] = xv.T
dout.variables['nav_lat'][:,:] = yv.T

dout.close()
din.close()



# copy lat and lon variables


----
Target variables
++++++++++++++++

variables::

  float nav_lon(y, x) ;
    nav_lon:units = "degrees_east" ;
    nav_lon:valid_min = -3.540191f ;
    nav_lon:valid_max = -2.716762f ;
    nav_lon:long_name = "Longitude" ;
  float nav_lat(y, x) ;
    nav_lat:units = "degrees_north" ;
    nav_lat:valid_min = 53.11438f ;
    nav_lat:valid_max = 53.96194f ;
    nav_lat:long_name = "Latitude" ;
  double glamt(y, x) ;
    glamt:missing_value = 1.e+20f ;
  double glamu(y, x) ;
    glamu:missing_value = 1.e+20f ;
  double glamv(y, x) ;
    glamv:missing_value = 1.e+20f ;
  double glamf(y, x) ;
    glamf:missing_value = 1.e+20f ;
  double gphit(y, x) ;
    gphit:missing_value = 1.e+20f ;
  double gphiu(y, x) ;
    gphiu:missing_value = 1.e+20f ;
  double gphiv(y, x) ;
    gphiv:missing_value = 1.e+20f ;
  double gphif(y, x) ;
    gphif:missing_value = 1.e+20f ;
  double e1t(y, x) ;
    e1t:missing_value = 1.e+20f ;
  double e1u(y, x) ;
    e1u:missing_value = 1.e+20f ;
  double e1v(y, x) ;
    e1v:missing_value = 1.e+20f ;
  double e1f(y, x) ;
    e1f:missing_value = 1.e+20f ;
  double e2t(y, x) ;
    e2t:missing_value = 1.e+20f ;
  double e2u(y, x) ;
    e2u:missing_value = 1.e+20f ;
  double e2v(y, x) ;
    e2v:missing_value = 1.e+20f ;
  double e2f(y, x) ;
    e2f:missing_value = 1.e+20f ;
