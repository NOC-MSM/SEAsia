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

* ...

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
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL_XIOS1

For the generation of bathymetry we need to patch the code (to allow direct
passing of arguments. NB this code has not been updated in 7 years.)::

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
resolution grid elements.

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
      rho  = 5
      rhot = 5
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

Execute first scrip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Output files::

  remap_nemo_grid_gebco.nc
  remap_data_grid_gebco.nc

Execute second scip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Output files::

  data_nemo_bilin_gebco.nc

Execute third scip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

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

Somewhat arbitrarily I am going to use **5** ``(jpkdta=5)`` levels::

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
     ln_e3_dep   = .false.    ! =T : e3=dk[depth] in discret sens.
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !
     cp_cfg      =  "orca"   !  name of the configuration
     jp_cfg      =       0   !  resolution of the configuration
     jpidta      =    1624   !  1st lateral dimension ( >= jpi )
     jpjdta      =    1138   !  2nd    "         "    ( >= jpj )
     jpkdta      =       5   !  number of levels      ( >= jpk )
     jpiglo      =    1624   !  1st dimension of global domain --> i =jpidta
     jpjglo      =    1138   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_sco      = .true.    !  z-coordinate - partial steps
     ln_linssh   = .false.    !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     jphgr_msh   =       0               !  type of horizontal mesh
    ...


Also add in s-coordinate namelist entries. These were copied from AMM60::

  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
  ln_s_sh94   = .false.   !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
  ln_s_sf12   = .true.    !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
  ln_sigcrit  = .true.    !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                          !  stretching coefficients for all functions
  rn_hc       =   50.0    !  critical depth for transition to stretched coordinates
  rn_rmax     =    0.1    !  maximum cut-off r-value allowed (0<r_max<1)
                       !!!!!!!  SF12 stretching coefficient  (ln_s_sf12 = .true.)
  rn_alpha    =    4.4    !  stretching with SF12 s-sigma
  rn_efold    =    0.0    !  efold length scale for transition to stretched coord
  rn_zs       =    1.0    !  depth of surface grid box
  rn_sbot_min =   10.0    !  minimum depth of s-bottom surface (>0) (m)
  rn_sbot_max = 7000.0    !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
                      !  bottom cell depth (Zb) is a linear function of water depth Zb = H*a + b
  rn_zb_a     =    0.024  !  bathymetry scaling factor for calculating Zb
  rn_zb_b     =   -0.2    !  offset for calculating Zb
/



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

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
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

Also need the parent grid mask and mesh files, from ORCA12,  in $START_FILES. (e.g. for AMM60)::

  #rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/START_FILES/mask_AMM60.nc  $START_FILES/.
  #rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/START_FILES/mesh_zgr_AMM60.nc  $START_FILES/.
  #rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/START_FILES/mesh_hgr_AMM60.nc  $START_FILES/.


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

.. note : Can use the ``-g`` option if you want the GUI.

.. note : Previously I actually had to do a fix to get PyNEMO to behave See: https://bitbucket.org/jdha/nrct/issues/22/issue-with-applying-bc-over-the-rim

Now I get the error message::

  ...
  [1621 1134]
  [1621 1135]]
  INFO:pynemo.profile:done  bdy t,u,v,f
  INFO:pynemo.profile:1.82
  INFO:pynemo.profile:bdy_ind f (40089, 2) (40089,)
  INFO:pynemo.profile:bdy_ind u (40125, 2) (40125,)
  INFO:pynemo.profile:bdy_ind t (40150, 2) (40150,)
  INFO:pynemo.profile:bdy_ind v (40116, 2) (40116,)
  INFO:pynemo.profile:done coord gen
  INFO:pynemo.profile:0.0
  INFO:pynemo.profile:./domain_cfg.nc
  INFO:pynemo.profile:done coord pop
  INFO:pynemo.profile:0.4
  INFO:pynemo.profile:gather grid info
  INFO:pynemo.profile:0.01
  INFO:pynemo.profile:generating depth info
  INFO:pynemo.reader.ncml:(1, 1138, 1624)
  INFO:pynemo.reader.ncml:(5,)
  ERROR:pynemo.reader.ncml:Cannot find the requested variable hbatt
  Traceback (most recent call last):
   File "/login/jelt/.conda/envs/nrct_env/bin/pynemo", line 11, in <module>
     load_entry_point('pynemo==0.2', 'console_scripts', 'pynemo')()
   File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_exe.py", line 44, in main
     profile.process_bdy(setup_file, mask_gui)
   File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 175, in process_bdy
     z = zgrv.Depth(grid_t.bdy_i, grid_u.bdy_i, grid_v.bdy_i, settings)
   File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/nemo_bdy_zgrv2.py", line 59, in __init__
     hbatt[mbathy == 0] = np.NaN
  TypeError: 'NoneType' object does not support item assignment


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

.. note : There is a problem with the processing the other input data, which I am not using*
  and which results in a pynemo crash::

  File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/nemo_bdy_extr_tm3.py", line 140,
  dejavu_sorted_index:  [] is emtpy



Copy the new files back onto ARCHER
::

  livljobs4$
  cd $INPUTS
  for file in SWPacific*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/. ; done
  scp coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/INPUTS/.


8. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

Get important files into EXP directory. Should already have ``domain_cfg.nc``::


  cd $EXP
  cp $INPUTS/bathy_meter.nc $EXP/.
  cp $INPUTS/coordinates.nc $EXP/.
  cp $INPUTS/coordinates.bdy.nc $EXP/.
  cp $START_FILES/namelist_cfg $EXP/.

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

  e.g. rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/namelist_cfg SWPacific_v4_namelist_cfg


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
