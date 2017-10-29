=============================================
Setting up a SW Pacific NEMO v4 configuration
=============================================

Machines: livljobs4, ARCHER

URL:: *to add*
*Actually I've taken this off readthedocs as it is internal. However writing in*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a SW Pacific regional tide-only model using GEBCO bathymetry, FES tidal
boundaries.

Build on a combination of livljobs4 and ARCHER.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. ARCHER: Run simulation

It is a tide only run.

Issues that arose
=================

* I had real trouble getting this to be stable is s-ccordinates, thought it worked
for the SEAsia doman. S-coordinates affects the domain_cfg.nc generation and PyNEMO
process. So I have made "protected" copies of the ``ARCHER:$TIDR/DOMAINcfg/SCO``
``LIVLJOBS:$INPUTS/SCO`` directories, and tried again with z-coords (following \
James' examples).


.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO bits are done on livljobs4.
(There was a problem with some java bits working on ARCHER)
Starting on ARCHER::

  ssh login.archer.ac.uk

  export CONFIG=SWPacific
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

Note you might have to mkdir the odd directory or two...::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES
  cp $WDIR/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WDIR/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.


Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build (if it is already downloaded)::

  cd $CDIR
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



.. note :
  Skip 1 and 2 - COPY bathy_meter.nc and coordinates.nc from Tom
    /work/thopri/NEMO/SWPacific_ver3.6/INPUTS/bathy_meter.nc
    /work/thopri/NEMO/SWPacific_ver3.6/INPUTS/coordinates.nc

    Otherwise I'd follow the following instuctions.



1. Generate new coordinates file
++++++++++++++++++++++++++++++++

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements. Originally we tried to interpolate to a finer grid with
``rho=rhot=3``. Now keeping resolution at R12 incase it sorts out a stability issue.

.. warning:
  Using the GRIDGEN/create_coordinates.exe tool runs into a problem for zoom factor
  >1, since the horizontal spacing metric e.g. e[12]t always match
  the parent grid. I think that this is a bug. The agrif version works.

Edit namelist::

  cd $TDIR/NESTING
  vi namelist.input

  &input_output
      iom_activated = true
  /
  &coarse_grid_files
      parent_coordinate_file = 'coordinates_ORCA_R12.nc'
  /
  &bathymetry
  /
  &nesting
      imin = 865
      imax = 1405
      jmin = 1116
      jmax = 1494
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

Build and execute agrif version of create_coordinates.exe.
See `build_and_create_coordinates.rst`_

This creates a new coordinatesfile with contents, which is now copied to
INPUTS::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

Now we need to generate a bathymetry on this new grid.



2. Generate bathymetry file
+++++++++++++++++++++++++++

thopri downloaded some and merged some 1 minute GEBCO data for (-30N :0N , -170E : 145E ).
Method in ``/work/thopri/NEMO/SWPacific/START_FILES/gebco_lon_convertor.py``
Copy to ARCHER::

  livljobs4$ scp /work/thopri/NEMO/SWPacific_ver3.6/START_FILES/GRIDONE_2D_140_-35.0_-165.0_5.0.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/.

Copy namelist for reshaping GEBCO data::

  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h GRIDONE_2D_140_-35.0_-165.0_5.0.nc`` to get input
variable names)::

  vi $INPUTS/namelist_reshape_bilin_gebco
  ...
  &grid_inputs
    input_file = 'gebco_in.nc'
    nemo_file = 'coordinates.nc'
    ...
    input_lon = 'lon'
    input_lat = 'lat'
    nemo_lon = 'glamt'
    nemo_lat = 'gphit'
    ...

    &interp_inputs
    input_file = "gebco_in.nc"
    ...
    input_name = "elevation"


Do some things to 1) flatten out land elevations, 2) make depths positive.
Have to swap around with the modules to get nco working *(James
noted a problem with the default nco module)*::

  cd $INPUTS

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2D_140_-35.0_-165.0_5.0.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

Restore the original parallel modules::

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

Execute first scrip thing::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Output files::

  remap_nemo_grid_gebco.nc
  remap_data_grid_gebco.nc

Execute second scip thing::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Output files::

  data_nemo_bilin_gebco.nc

Execute third scip thing::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

Output files::

  bathy_meter.nc


3. Generate initial conditions
++++++++++++++++++++++++++++++

.. note : Skip this

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

    cp $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
    cp $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template namelist_cfg with only the essenetial domain building stuff.
Get the indices from ``ncdump -h coordinates.nc``.

When the increase in resolution was x3, the size of the coordinate were::

  jpidta      =    1624   !  1st lateral dimension ( >= jpi )
  jpjdta      =    1138   !  2nd    "         "    ( >= jpj )

When the increase in resolution is x1, the size of the coordinates were::

  jpidta      =    544   !  1st lateral dimension ( >= jpi )
  jpjdta      =    382   !  2nd    "         "    ( >= jpj )


Somewhat arbitrarily I am going to use **31** ``(jpkdta=31; rn_jpk=31)`` levels.


s-coordinates `SWPacific_DOMAINcfg_namelist_cfg`_ (didn't work)
-------------

 (If I tried 5 levels then the function to compute the depth range breaks).
 Added in s-coordinate settings from AMM60. (There is a function that tapers dz
 near the equation. This is activated in ``domzgr.F90`` and not by a logical flag)::

  cd $TDIR/DOMAINcfg
  vi namelist_cfg

  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  Configuration namelist : used to overwrite defaults values defined in SHARED/namelist_ref
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !
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
     jpjdta      =     382   !  2nd    "         "    ( >= jpj )
     jpkdta      =      31   !  number of levels      ( >= jpk )
     jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     382   !  2nd    -                  -    --> j  =jpjdta
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
    ln_s_sh94   = .false.    !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
    ln_s_sf12   = .true.   !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
    ln_sigcrit  = .true.    !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                            !  stretching coefficients for all functions
    rn_jpk      =   31       ! Number of S levels
    !ln_eq_taper = .false.   !  Tapering of S coords near equator
    cn_coord_hgr = 'coordinates.nc'  ! File containing gphit (latitude) coordinate for use if ln_eq_taper=.true.
    rn_sbot_min =   10.0    !  minimum depth of s-bottom surface (>0) (m)
    rn_sbot_max = 7000.0    !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
    rn_hc       =   50.0    !  critical depth for transition to stretched coordinates
                         !!!!!!!  Envelop bathymetry
    rn_rmax     =    0.3    !  maximum cut-off r-value allowed (0<r_max<1)
                         !!!!!!!  SH94 stretching coefficients  (ln_s_sh94 = .true.)
    rn_theta    =    6.0    !  surface control parameter (0<=theta<=20)
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


z-coordinates
-------------

::
  cd $TDIR/DOMAINcfg
  vi namelist_cfg

  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  Configuration namelist : used to overwrite defaults values defined in SHARED/namelist_ref
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !
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
     jpjdta      =     382   !  2nd    "         "    ( >= jpj )
     jpkdta      =      31   !  number of levels      ( >= jpk )
     jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     382   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
    ln_zco      = .false.   !  z-coordinate - full    steps
    ln_zps      = .true.   !  z-coordinate - partial steps
    ln_sco      = .false.   !  s- or hybrid z-s-coordinate
    ln_isfcav   = .false.   !  ice shelf cavity
    ln_linssh   = .false.   !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     jphgr_msh   =       0               !  type of horizontal mesh
     ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
     ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
     ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
     ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
     ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
     ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
     ppsur       =   -4762.96143546300   !  ORCA r4, r2 and r05 coefficients
     ppa0        =     255.58049070440   ! (default coefficients)
     ppa1        =     245.58132232490   !
     ppkth       =      21.43336197938   !
     ppacr       =       3.0             !
     ppdzmin     =  999999.              !  Minimum vertical spacing
     pphmax      =  999999.              !  Maximum depth
     ldbletanh   =  .FALSE.              !  Use/do not use double tanf function for vertical coordinates
     ppa2        =  999999.              !  Double tanh function parameters
     ppkth2      =  999999.              !
     ppacr2      =  999999.              !
  /

----



.. note:

  No gdept output in the offical v4 release. Though it was acheived here setting
  ln_e3_dep = F. This is needed for PyNEMO, though could be constructed from e3[tw].

Build a script to run the executable::

  vi $TDIR/DOMAINcfg/rs

  #!/bin/bash
  #PBS -N domain_cfg
  #PBS -l walltime=00:20:00
  #PBS -l select=1
  #PBS -j oe
  #PBS -A n01-NOCL
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


Try running it::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Put a copy of the namelist_cfg in $INPUTS for safe keeping::

  cp $TDIR/DOMAINcfg/namelist_cfg $INPUTS/namelist_cfg_generateDOMAINcfg

#Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

  cp $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  cp $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.




5. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

.. note : skip this for now

  Generate cut down drowned precip file (note that the nco tools don't like the
  parallel modules). **HEALTH WARNING** *Cut out files with only one index in that lat direction broke NEMO*::

    module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
    module load cray-netcdf cray-hdf5
    module load nco/4.5.0
    ncks -d lon,355.,360. -d lat,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_precip_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_precip_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_u10_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_u10_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_v10_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_v10_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_radsw_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_radsw_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_radlw_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_radlw_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_t2_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_t2_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_q2_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_q2_DFS5.1.1_y2000.nc
    ncks -d lon0,355.,360. -d lat0,48.,55. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_snow_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_snow_DFS5.1.1_y2000.nc

    module unload nco/4.5.0
    module unload cray-netcdf cray-hdf5
    module load cray-netcdf-hdf5parallel cray-hdf5-parallel

  Obtain namelist files and data file::

    cp $START_FILES/namelist_reshape_bilin_atmos $INPUTS/.
    cp $START_FILES/namelist_reshape_bicubic_atmos $INPUTS/.

  Edit namelist to reflect source filenames (just a year change)::

    vi $WDIR/INPUTS/namelist_reshape_bilin_atmos
    ...
    &grid_inputs
        input_file = 'cutdown_drowned_precip_DFS5.1.1_y2000.nc'

    vi $WDIR/INPUTS/namelist_reshape_bicubic_atmos
    ...
    &grid_inputs
      input_file = 'cutdown_drowned_precip_DFS5.1.1_y2000.nc'


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



THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO *(20 Sept 2017)*
If all the files are ready to go jump straight to `7. Generate boundary conditions with PyNEMO: Run PyNEMO`_

----

Follow SEAsia notes `SEAsia_archer_livljobs4.rst`_ for PyNEMO bit. Changing instances of `SEAsia` for `SWPacific`
Note a lot of the files are actually used for a tide-only simulation but it is not
trivial to prevent PyNEMO looking for them.

----

Set up directory structure and files in livljobs4
=================================================

Define paths::

cat > ~/temporary_path_names_for_NEMO_build << EOL
export CONFIG=SWPacific
export WORK=/work
export WDIR=\$WORK/$USER/NEMO/\$CONFIG
export INPUTS=\$WDIR/INPUTS
export START_FILES=\$WDIR/START_FILES
#export CDIR=\$WDIR/trunk_NEMOGCM_r8395/CONFIG
#export TDIR=\$WDIR/trunk_NEMOGCM_r8395/TOOLS
#export EXP=\$CDIR/\$CONFIG/EXP00

EOL

Execute path settings::

  . ~/temporary_path_names_for_NEMO_build

Synchronise with key files from ARCHER::

  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/bathy_meter.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/domain_cfg.nc  $INPUTS/.


Statement about external forcing and files
==========================================

Use ORCA 1/12 data via a thredds server. Also use ORCA 1/12 mesh and mask files
 via a thredds server.
 Copy necessary files into INPUTS. (Be careful of symbolic links in PyNEMO).::

   ls -lh $INPUTS/bathy_meter.nc
   ls -lh $INPUTS/coordinates.nc
   ls -lh $INPUTS/domain_cfg.nc

Need to generate 6 more files: A ``namelist.bdy`` (archived as `SWPacific_namelist.bdy`_)
which drives PyNEMO and which has two input files: ``inputs_src.ncml``
which points to the remote ORCA data source and ``inputs_dst.ncml`` which
remaps some variable names in the destination files. Then generate ncml files
to get the mesh and mask files (``mask_src.ncml  mesh_hgr_src.ncml  mesh_zgr_src.ncml``).

Create ncml file for ORCA12 source data::

  vi  $INPUTS/inputs_src.ncml

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


Create NCML file for mapping variables in DESTINATION grid to what PyNEMO expects::

  vi inputs_dst.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf location="file:domain_cfg.nc">
      <ns0:variable name="mbathy" orgName="top_level" />
      <ns0:variable name="gdept" orgName="gdept_0" />
      <ns0:variable name="gdepw" orgName="gdepw_0" />
      <ns0:variable name="e3u" orgName="e3u_0" />
      <ns0:variable name="e3v" orgName="e3v_0" />
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>

Create a PyNEMO namelist.bdy file. Note this this is tides only and so z-coordinates
are OK as it is only 2D variables::

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
     sn_src_dir = './inputs_src.ncml'       ! src_files/'
     sn_dst_dir = '/work/jelt/NEMO/SWPacific/INPUTS/'
     sn_fn      = 'SWPacific'                 ! prefix for output files
     nn_fv      = -1e20                     !  set fill value for output files
     nn_src_time_adj = 0                                    ! src time adjustment
     sn_dst_metainfo = 'metadata info: jelt'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
  ln_coords_file = .true.               !  =T : produce bdy coordinates files
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
  clname(1)      = 'M2'                 ! constituent name
  clname(2)      = 'S2'
  clname(3)      = 'K2'
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


Tried to put the mesh and mask files as remote access via ncml files. Failed, so
downloaded them instead...

.. note:

    Finally, create separate ncml files for each mesh mask file (because I don't know
    how to join them...) and they are on a thredds server::

      vi mesh_hgr_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mesh_hgr.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>


      vi mesh_zgr_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mesh_zgr.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>


      vi mask_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mask.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>

Finally I am going to create a boundary mask file (I don't think PyNEMO likes it if you
don't have one but have a rimwidth>1)::

  module load nco/gcc/4.4.2.ncwa
  ncks -v top_level domain_cfg.nc tmp.nc
  ncrename -h -v top_level,mask tmp.nc bdy_mask.nc
  rm tmp.nc

In ipython::

  import netCDF4
  dset = netCDF4.Dataset('bdy_mask.nc')
  dset['mask'][:][dset['mask'][:,0,:] ] = 0 # zero the rim
  dset['mask'][:][0,0:4,:] = 0   # Extra rim to remove blow up pt
  dset['mask'][:][dset['mask'][:,1,0] ] = 0 # zero the rim
  dset.close() # if you want to write the variable back to disk

.. note: Could copy all these files to START_FILES or the repo...


Run PyNEMO / NRCT to generate boundary conditions
+++++++++++++++++++++++++++++++++++++++++++++++++

First install PyNEMO `install_nrct`_ if not already done so.

Generate the boundary conditions with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

.. note : Can use the ``-g`` option if you want the GUI. Sometime problems are
  circumvented using the "-g" option.

.. note : log revealed:
    ...
    INFO:pynemo.profile:horizontal grid info
    ERROR:pynemo.reader.ncml:Cannot find the requested variable glamt
    ERROR:pynemo.reader.ncml:Cannot find the requested variable gphit
    ...
    Hmm. Something to look at is it doesn't work.

  .. note : Error
  File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/tide/nemo_bdy_tide3.py", line 50, in nemo_bdy_tpx7p2_rot
    dst_lon = DC.bdy_lonlat[Grid_U.grid_type]['lon'][Grid_U.bdy_r == 0]
    As before fixed but commenting out the [Grid_U.bdy_r == 0] subsetting. x 4.

This is because I have tried to get PyNEMO to use s-coordinates. These are not
needed for a tide only run

.. note PyNEMO does not work for s-coordinates. But to just use it for tidal
 boundary conditions it doesn't matter what the vertical grid is doing.

.. note: crashed in the tra routine. A problem with source coordinates zt not working as expected

27 Oct. z-ccords. Try rimwidth=1 because I can't get PyNemo to work with wider widths
(I think that if I'm imposing the same tide over a rim with sloping bathymetry I get problems)
Then expand the bdy_mask.nc file to inclide the blow-up point.

28 Oct. Fiddle with PyNEMO gui to get a 10pt border. Output to bdy_mask.nc variables:
nav_lat, nav_lon, mask.
---


This generates::
  ls -1 $INPUTS

  coordinates.bdy.nc
  SWPacific_bdytide_rotT_M2_grid_T.nc
  SWPacific_bdytide_rotT_K2_grid_T.nc
  SWPacific_bdytide_rotT_S2_grid_T.nc
  SWPacific_bdytide_rotT_M2_grid_U.nc
  SWPacific_bdytide_rotT_K2_grid_U.nc
  SWPacific_bdytide_rotT_S2_grid_U.nc
  SWPacific_bdytide_rotT_M2_grid_V.nc
  SWPacific_bdytide_rotT_K2_grid_V.nc
  SWPacific_bdytide_rotT_S2_grid_V.nc




Copy the new files back onto ARCHER
::

  livljobs4$
  cd $INPUTS
  for file in SWPacific*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/. ; done
  scp coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/.
  scp bdy_mask.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/. # variable mask - for pynemo
  #scp bdy_msk.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/.  # variable bdy_msk - for nemo


8. Run the configuration ON ARCHER. Turn on the tides
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

Create a short queue runscript::

  vi runscript

  #!/bin/bash
  #PBS -N SWPacific
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
    OCEANCORES=96
    XIOCORES=1
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core

  #aprun -n $OCEANCORES -N 24 ./opa
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa
  #aprun -b -n $XIOCORES -N 1 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

  exit

---

Edit ``namelist_cfg`` to make sure it is OK

Synchronise the namelist_cfg with the GitLab repo::

  e.g. rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/namelist_cfg SWPacific_EXP_namelist_cfg


Blow up at 79 ts
Use tide_ramp over 1 day.

Blow up after 197 steps (=197 minutes).

Will try lateral diffusion to stabilize. makes it worse...

Sea level blows up at an island in the rim region. Try decreasing rimwidth to
accomodate rimwidth=5.
---

Submit::

 cd $EXP
 qsub -q short runscript

27 Oct. z-coords. Need to archive namelists (especially if it works)...
Blows up on the southern boundary again....

Plan. 1) Go back to rimwidth=1 (perhaps blow up comes from same tide imposed
over varying bathymetry across the rimwidth).
Blowup at  kt=   100 max ssh:   10.32    , i j:   172    4
Then try and edit the mask to remove the bad points. This blew up in the same place
Check the output.abort files. Mesh hasn't worked as expected...

.. note : I notice that the mask file I created does not have nav_lat or nav_lon.
  Also there is a time axis in the mask variable. The GUI generated mask file is
  cleaner: 	float nav_lat(y, x) ;
	float nav_lon(y, x) ;
	float mask(y, x) ;

28 Oct. Tried using PyNEMO to generate a 10pt mask. Turned off the mask in the EXP namelist_cfg
Rimwidth=1.
**PENDING** 4877760::
  stp_ctl : the ssh is larger than 10m
  =======
  kt=   410 max ssh:   10.08    , i j:   446   11

Chop out SSH also create mask file ::


  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncks -v sossheig output.abort.nc output.abort_ssh.nc

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

It blew up just inside the maks. Apply mask in the EXP run. REsubmit.
James' examply `bdy_msk.nc` has a variable with the same name.

**livljobs4**::

  cd $INPUTS
  module load nco/gcc/4.4.2.ncwa
  ncrename -h -v mask,bdy_msk bdy_mask.nc bdy_msk.nc
  scp bdy_msk.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/.

Resubmit: 4877863.sdb **PENDING** Does the mask work? Does it blow up somewhere else?
stp_ctl : the ssh is larger than 10m
=======
kt=    33 max ssh:   10.59    , i j:   535   11

Blows up in the corner. But the grids work.
Try the tideramp=1. Turn off tide_pot too

Does the tideramp slow the blow up? Yes try again with 10 day ramp
Blows up with checker boarding. Turn on bilaplacian diffusion.
stp_ctl : the ssh is larger than 10m
=======
kt=   103 max ssh:   11.45    , i j:   179   71
checker board blow up - away from the domain.

Copy James' bilaplacian coefficients.
stp_ctl : the ssh is larger than 10m
=======
kt=    62 max ssh:   12.66    , i j:   507    8
It still blows up within the rimwidth. Something is not working as I would expect.

The rimwidth=1. If I make it larger then the PyNEMO doesn't seem to work. This
seems to be the best line of pursuit.
Perhaps PyNEMO needed a mask file with the same rimwidth as the rimwidth variable.
Could quickly try that.

livljobs4: namelist.bdy
rimwidth=10
This breaks. I do not understand why it doesn't work but it seems to be a sticking point.


---

Backup to repo key files
========================

::

  cd ~/GitLab/NEMO-RELOC/docs/source
  # DOMANcfg namelist_cfg for domain_cfg.nc (for s-coordinates)
  rsync -utv jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/namelist_cfg SWPacific_DOMAINcfg_namelist_cfg

  # EXP namelist_cfg (for s-coordinates)
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/namelist_cfg SWPacific_EXP_namelist_cfg

  # PyNEMO namelist.bdy (for s-coordinates)
  rsync -utv jelt@livljobs4:/work/jelt/NEMO/SWPacific/INPUTS/namelist.bdy SWPacific_namelist.bdy

  # Python quick plot of SSH in the output.abort.nc file
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/quickplotNEMO.py quickplotNEMO.py
