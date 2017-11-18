==========================================
Setting up a SE Asia NEMO v4 configuration
==========================================

Machines: livljobs4, ARCHER

URL:: *to add*
*Actually I've taken this off readthedocs as it is internal. However writing in*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a SE Asia regional tide-only model using GEBCO bathymetry, TPXO(later FES?)
 tidal boundaries.

Build on a combination of livljobs4 and ARCHER.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: *skip: Generate initial conditions and atmospheric forcings*
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

  export CONFIG=SEAsia
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

  cp $WDIR/../LBay/START_FILES/dommsk.F90 $START_FILES/.
  cp $WDIR/../LBay/START_FILES/bdyini.F90 $START_FILES/.
  cp $WDIR/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WDIR/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  cp $WDIR/../SWPacific/START_FILES/usrdef_istate.F90 $START_FILES/.
  cp $WDIR/../SWPacific/START_FILES/usrdef_sbc.F90    $START_FILES/.


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

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements.

.. warning:
  Using the GRIDGEN/create_coordinates.exe tool runs into a problem for zoom factor
  >1, since the horizontal spacing metric e.g. e[12]t always match
  the parent grid. I think that this is a bug. The agrif version works.

1a. Choose domain
=================

First we need to figure out the indices for the new domain, from the parent grid.
Move parent grid into INPUTS::

  cp $START_FILES/coordinates_ORCA_R12.nc $INPUTS/.

Inspect this parent coordinates file to define the boundary indices for the new config.

Note, I used FERRET locally::

  $livljobs2$ scp $USER@login.archer.ac.uk:/work/n01/n01/$USER/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.
  ferret etc
  shade/i=3385:3392/j=2251:2266 NAV_LAT
  shade/i=3385:3392/j=2251:2266 NAV_LON

**Use indices  i=50:730 j=1250:1800**

---

**Longer version**

Inspect this parent coordinates file to define the boundary indices for the new config.

Note, I used FERRET on livljobs4.

*(27 Sept 2017)*

Decide coordinates for new SE Asia configuration at 1/12 degree, R12
====================================================================

Inspect TPXO harmonic amplitudes to find a good cut off location for boundaries::

  livljobs4$ cd /work/jelt/tpxo7.2
  ferret
  go plot_SEAsia_harmonics.jnl

... note::

  ! plot_SEAsia_harmonics.jnl
  ! Plot tpxo harmonics for the SE Asia region.
  ! Want to build a NEMO config without significant amphidromes on the boundary

  use h_tpxo7.2.nc

  set win 1
  set viewport ul
  shade/k=1/j=300:700/i=250:500/levels=(0,1,0.1)/title="M2" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=2/j=300:700/i=250:500/levels=(0,1,0.1)/title="S2" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=3/j=300:700/i=250:500/levels=(0,1,0.1)/title="N2" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=4/j=300:700/i=250:500/levels=(0,1,0.1)/title="K2" HA, lon_z, lat_z; go fland

  set win 2
  set viewport ul
  shade/k=5/j=300:700/i=250:500/levels=(0,1,0.1)/title="K1" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=6/j=300:700/i=250:500/levels=(0,1,0.1)/title="O1" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=7/j=300:700/i=250:500/levels=(0,1,0.1)/title="P1" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=8/j=300:700/i=250:500/levels=(0,1,0.1)/title="Q1" HA, lon_z, lat_z; go fland


Conclusion. Plot the proposed domain::

  $livljobs2$ scp $USER@login.archer.ac.uk:/work/n01/n01/$USER/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.

  ferret
  use coordinates_ORCA_R12.nc
  set win 1; shade/X=50:730/Y=1250:1800 E2T, nav_lon, nav_lat ; go fland
  set win 2; set viewport upper; shade/i=50:730/j=1250:1800 NAV_LAT
  set win 2; set viewport lower; shade/i=50:730/j=1250:1800 NAV_LON

Use indices  **i=50:730 j=1250:1800**

---


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
      imin = 50
      imax = 730
      jmin = 1250
      jmax = 1800
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

Or just execute tool::

  ./agrif_create_coordinates.exe

This creates a new coordinates file with contents, which is now copied to
INPUTS::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

Now we need to generate a bathymetry on this new grid.


2. Generate bathymetry file
+++++++++++++++++++++++++++

.. note:
  This was first done on livljobs4. Here I do it on ARCHER.

Take GEBCO bathymetry. For a domain as large as SE Asia, the 30-minute GEBCO
data is too large to process and needs some spatial filtering. BODC also host a
 1-minute data set (2008) which should work without pre-processing but is not
  updated.

.. warning:

  A 30-second GEBCO cutout is too large to process for the SE Asia domain
  (7081 x 5521 pts). The older 1-minute data is fine.

Download some GEBCO 2014 and 2008 data (75E,-21N,134E,25N) and copy to $INPUTS::

 scp GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/SEAsia/INPUTS/.

.. note: Copying to livljobs4

  livmaf$
  scp ~/Downloads/RN-9621_1506544326915/GEBCO_2014_2D_75.0_-21.0_134.0_25.0.nc jelt@livljobs4.nerc-liv.ac.uk:$INPUTS/GEBCO_2014_2D5.0_-21.0_134.0_25.0.nc
  scp ~/Downloads/RN-6060_1506606001516/GRIDONE_2D_74.0_-21.0_134.0_25.0.nc jelt@livljobs4.nerc-liv.ac.uk:$INPUTS/GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc

**In the following I use the 2008 data**
Copy namelist for reshaping GEBCO data::

  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc`` to get input
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
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc tmp.nc
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

.. note: ferret

 use bathy_meter.nc
 shade log(BATHYMETRY), nav_lon, nav_lat; go land


3. Generate initial conditions
++++++++++++++++++++++++++++++

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

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
   rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.

.. mote :  should check the difference between the homemade sco version the AMM60
  verison did:      ``diff namelist_cfg_sco_WIP namelist_cfg_AMM60``

.. note : alternativly should check the difference between the AMM60 and local
  output.namelist.dyn: ``diff output.namelist.dyn /work/n01/n01/jelt/NEMO/NEMOGCM/CONFIG/AMM60smago/EXP_NSea/output.namelist.dyn``
  I notice that rmax is different.

.. note : There are a lot of unknown parameters in these settings. And I don't
  want to find i made some naive error in six months. Looking at the domain there
  are some serious trenches near land. S-coords will not work well there. Conversely,
  ODA is all abount the near-coastal environment. There is a strong case for using
  hybrid s-z coordinates a la NNA...

  James noted that high resolution neat the bed caused significant difficulty in
  deep water stability. Whereas you want it on the shelf. Hence regular stretched
  s-coords wont really work.

  PyNEMO outputs boundary conditions on the parent z-grid. This can be interpolated
  at run-time to the child grid.

5. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

Initially use zero atm forcing. Specified in usr defined functions in MY_SRC.

.. note: Comment out weight for atm forcing

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

Statement about external forcing
================================

Uses ORCA 1/12 via a thredds server.
I have the mesh and mask files ``mask_src.nc  mesh_hgr_src.nc  mesh_zgr_src.nc``
 stored locally (from the lighthouse reef experiment).

 Copy necessary files into INPUTS::

   cp $START_FILES/mask_src.nc     $INPUTS/.
   cp $START_FILES/mesh_hgr_src.nc $INPUTS/.
   cp $START_FILES/mesh_zgr_src.nc $INPUTS/.

   ls -lh $INPUTS/bathy_meter.nc
   ls -lh $INPUTS/coordinates.nc
   ls -lh $INPUTS/domain_cfg.nc

Need to generate 3 more files: A ``thredds_namelist.bdy`` which drives PyNEMO and which
has two input files: ``thredds_inputs_src.ncml`` which points to the data source and
``inputs_dst.ncml`` which remaps some variable names in the destination files.

6. Generate boundary conditions with NRCT/PyNEMO: Create netcdf abstraction wrapper
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

First install PyNEMO `install_nrct`_ if not already done so.


6a. Generate ncml file that points to the external data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

This can be done with the automatic generator (*pynemo_ncml_generator*) or manually

Here the object is to generate a ncml file that is read in by PyNEMO as the ``sn_src_dir``
(in the ``namelist.bdy`` file)

.. note: If using the generator, fill in the Tracer and Dynamics for T,S,U,V,Z
 tabs: using T,T & U,V,T in the reg expressions e.g. .*T\.nc$. To generate an
  e.g. ``inputs_src.ncml`` file click  **generate**. Defining the filename seems
   to work better with the file selector rather than direct typing.


Note need to set the time variables and new ``sn_src_dir`` in namelist.bdy.
 (Time variables correspond to simulation window and the time_origin for the time
axis of these data). Actually upated the following with all the Nov 1979 files::

  cd $INPUTS
  vi thredds_inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791201d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791126d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791121d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791116d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791111d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05U.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791201d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791126d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791121d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791116d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791111d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05V.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>




6b. Generate the namelist.bdy file for PyNEMO / NRCT
+++++++++++++++++++++++++++++++++++++++++++++++++++

Copy the NRCT template namelist.bdy from the START_FILES::

  cd $INPUTS
  cp $START_FILES/thredds_namelist.bdy $INPUTS/.

Edit namelist.bdy to for the configuration name and ``ncml`` file name::

  vi thredds_namelist.bdy
  sn_src_dir = './thredds_inputs_src.ncml'       ! src_files/'
  sn_dst_dir = '/work/n01/n01/jelt/SEAsia/INPUTS/'
  sn_fn      = 'SEAsia'                 ! prefix for output files
  ...

Make sure the timestamps correspond to the input data in ``*_inputs_src.ncml``.
Turn off as many things as possible to help it along.
Turned off ``ln_mask_file``. James said it was for outputting a new mask file
but it might have given me trouble. *Actually I also turn off all the ORCA inputs*.

Point to the correct source and destination mesh and mask files/variables.
 ::

   vi thredds_namelist.bdy

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
      ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)
      ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
      ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
      rn_hmin     =   -10     !  min depth of the ocean (>0) or
                              !  min number of ocean level (<0)

   !-----------------------------------------------------------------------
   !   s-coordinate or hybrid z-s-coordinate
   !-----------------------------------------------------------------------
      rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
      rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
                              !  (= ocean depth) (>0) (m)
      ln_s_sigma  = .true.   !  hybrid s-sigma coordinates
      rn_hc       =  150.0    !  critical depth with s-sigma

   !-----------------------------------------------------------------------
   !  grid information
   !-----------------------------------------------------------------------
      sn_src_hgr = './mesh_hgr_src.nc'   !  /grid/
      sn_src_zgr = './mesh_zgr_src.nc'
      sn_dst_hgr = './domain_cfg.nc'
      sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
      sn_src_msk = './mask_src.nc'
      sn_bathy   = './bathy_meter.nc'

   !-----------------------------------------------------------------------
   !  I/O
   !-----------------------------------------------------------------------
      sn_src_dir = './thredds_inputs_src.ncml'       ! src_files/'
      sn_dst_dir = '/work/jelt/NEMO/SEAsia/INPUTS/'
      sn_fn      = 'SEAsia'                 ! prefix for output files
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
        ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
        ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
        ln_tra         = .false.               !  boundary conditions for T and S
        ln_ice         = .false.               !  ice boundary condition
        nn_rimwidth    = 1                     !  width of the relaxation zone

   !-----------------------------------------------------------------------
   !  unstructured open boundaries tidal parameters
   !-----------------------------------------------------------------------
       ln_tide        = .true.               !  =T : produce bdy tidal conditions
       clname(1) =  'M2'
       clname(2) =  'S2'
       clname(3) =  'N2'
       clname(4) =  'K2'
       clname(5) =  'K1'
       clname(6) =  'O1'
       clname(7) =  'P1'
       clname(8) =  'Q1'
       clname(9) =  'M4'
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
       rn_mask_max_depth = 300.0     !  Maximum depth to be ignored for the mask
       rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break

.. warning:

  It doesn't quite work with ``ln_tra = .false.``

.. note :

  I thought that I needed to create a bdy_mask.nc file so I did this from doman_cfg.nc
  though it turns out not to have been needed. Nevertheless I did the following::

      module load nco/gcc/4.4.2.ncwa
      ncks -v top_level domain_cfg.nc tmp.nc
      ncrename -h -v top_level,mask tmp.nc bdy_mask.nc
      rm tmp.nc


Also had to check/create ``inputs_dst.ncml``, that it has the correct file name within:
 *Now domain_cfg.nc, formerly mesh_zgr.nc*. Note also that some variables in
  domain_cfg.nc have different names e.g. ``mbathy`` --> ``bottom_level``. Check the mapping
  in ``inputs_dst.ncml``::

   vi inputs_dst.ncml

    <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf location="file:domain_cfg.nc">
        <ns0:variable name="mbathy" orgName="bottom_level" />
        <ns0:variable name="e3u" orgName="e3u_0" />
        <ns0:variable name="e3v" orgName="e3v_0" />
        <ns0:variable name="e3t" orgName="e3t_0" />
        <ns0:variable name="e3w" orgName="e3w_0" />
        </ns0:netcdf>
      </ns0:aggregation>
    </ns0:netcdf>


  .. warning:
    In the actual v4 release domain_cfg.nc  will not have gdept or gdepw. These
    will need to be reconstructed from e3[tw].

  .. note : 18 Nov.  comment out the gdept and gdepw lines and
     inserted e3t and e3w. Previouly the inputs_dst.ncml looked like::

    <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf location="file:domain_cfg.nc">
        <ns0:variable name="mbathy" orgName="bottom_level" />
        <ns0:variable name="gdept" orgName="gdept_0" />
        <ns0:variable name="gdepw" orgName="gdepw_0" />
        <ns0:variable name="e3u" orgName="e3u_0" />
        <ns0:variable name="e3v" orgName="e3v_0" />
        </ns0:netcdf>
      </ns0:aggregation>
    </ns0:netcdf>





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

.. note : I have a PyNEMO mod to use FES tides instead of TPXO tides for these boundary
  forcing. It is currently a hardwire fix in ``tide/nemo_bdy_tide3.py``

This generates::
  ls -1 $INPUTS

  coordinates.bdy.nc
  SEAsia_bdytide_rotT_M2_grid_T.nc
  SEAsia_bdytide_rotT_N2_grid_T.nc
  SEAsia_bdytide_rotT_S2_grid_T.nc
  SEAsia_bdytide_rotT_K1_grid_T.nc
  SEAsia_bdytide_rotT_K2_grid_T.nc
  SEAsia_bdytide_rotT_P1_grid_T.nc
  SEAsia_bdytide_rotT_O1_grid_T.nc
  SEAsia_bdytide_rotT_M4_grid_T.nc
  SEAsia_bdytide_rotT_Q1_grid_T.nc
  SEAsia_bdytide_rotT_M2_grid_U.nc
  SEAsia_bdytide_rotT_N2_grid_U.nc
  SEAsia_bdytide_rotT_S2_grid_U.nc
  SEAsia_bdytide_rotT_K1_grid_U.nc
  SEAsia_bdytide_rotT_K2_grid_U.nc
  SEAsia_bdytide_rotT_P1_grid_U.nc
  SEAsia_bdytide_rotT_O1_grid_U.nc
  SEAsia_bdytide_rotT_M4_grid_U.nc
  SEAsia_bdytide_rotT_Q1_grid_U.nc
  SEAsia_bdytide_rotT_M2_grid_V.nc
  SEAsia_bdytide_rotT_N2_grid_V.nc
  SEAsia_bdytide_rotT_S2_grid_V.nc
  SEAsia_bdytide_rotT_K1_grid_V.nc
  SEAsia_bdytide_rotT_K2_grid_V.nc
  SEAsia_bdytide_rotT_P1_grid_V.nc
  SEAsia_bdytide_rotT_O1_grid_V.nc
  SEAsia_bdytide_rotT_M4_grid_V.nc
  SEAsia_bdytide_rotT_Q1_grid_V.nc


Copy the new files back onto ARCHER
::

  livljobs4$
  cd $INPUTS
  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.bdy.nc
  for file in $CONFIG*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done




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

Edit the output to have 1hrly SSH, and harmonic output::

 vi file_def_nemo.xml
 ...
 <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
  <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
    <field field_ref="ssh"          name="zos"   />
  </file>
 </file_group>
 ...
 <file_group id="5d" output_freq="5d"  output_level="10" enabled=".TRUE.">  <!-- 5d files -->
   <file id="file8" name_suffix="_D2_Tides" description="tidal harmonics" >
     <field field_ref="M2x"          name="M2x"      long_name="M2 Elevation harmonic real part"                       />
     <field field_ref="M2y"          name="M2y"      long_name="M2 Elevation harmonic imaginary part"                  />
     <field field_ref="M2x_u"        name="M2x_u"    long_name="M2 current barotropic along i-axis harmonic real part "       />
     <field field_ref="M2y_u"        name="M2y_u"    long_name="M2 current barotropic along i-axis harmonic imaginary part "  />
     <field field_ref="M2x_v"        name="M2x_v"    long_name="M2 current barotropic along j-axis harmonic real part "       />
     <field field_ref="M2y_v"        name="M2y_v"    long_name="M2 current barotropic along j-axis harmonic imaginary part "  />
     ...
   </file>
 </file_group>

---

Create a short queue runscript (Note: PBS -N jobname, PBS -m email)::

  vi runscript

  #!/bin/bash
  #PBS -N SEAsia
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
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

  exit

Change the notification email to your own address::

  sed -i "s/xxx@noc/$USER@noc/g" runscript

Might also want to change the account name. E.g. ``n01-ACCORD``

---

Edit ``namelist_cfg`` to make sure it is OK

---
*IT WORKS. IF IT WORKS, ARCHIVE namelist_cfg too**

---
*(17 Nov 17)* build new 75 level hybrid z-s coordinates. Submitted cold start
 20 min job.
DID IT WORK? Yes. (Just completed the 1440 steps in 20mins)

*(18 Nov 17)* Add in lots of TPXO harmonics. Run again with 40mins. Completed in 21mins.
With rn_rdt=60 this is only 1 day of simulation.
Try increasing the timestep.

rn_rdt = 360 and resubmit. **PENDING** 4922962.sdb

----

No met (missing slp) ``ln_usr=T``. rn_rdt=60s. Output more harmonics (20-30days).
Run for 30 days::

 cd $EXP
 qsub -q short runscript

**IT WORKS!** Hit wall time of 20mins after ~2 hours

Submit a big job on 2k processors to get through the spin up *(Need to do this efficiently)*::

  vi runscript
  #PBS -N SEAsia
  #PBS -l select=92
  #PBS -l walltime=00:20:00
  …
    echo " ";
    OCEANCORES=200
    XIOSCORES=40
  …
  aprun -b -n $XIOSCORES -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa


Ran for 60 hrs before hitting 20 min wall time. (NB 51 levels)
use python script to plot SSH animation (NB need to put the python script somewhere better!)::

  % python SEAsia_SSH_anim.py

Creates an animation of hours 35 - 60 in SSH.


---

MPP decomposition for land suppression
++++++++++++++++++++++++++++++++++++++

`MPP_decomp_lanf_suppression.rst`_

Before doing long runs it is important to optimise MPP decompositoin by invoking
 land supression to save redundant ocean processors.
Resulting decomposition::

   vi namelist_cfg
   ...
   !-----------------------------------------------------------------------
   &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
   !-----------------------------------------------------------------------
      ...
      jpni        =  12       !  jpni   number of processors following i (set automatically if < 1)
      jpnj        =  8    !  jpnj   number of processors following j (set automatically if < 1)
      jpnij       =  92    !  jpnij  number of local domains (set automatically if < 1)

Inspect ``ocean_output`` to find ``jpnij``. In my simulation ``jpni=12, jpnj=8 --> jpnij = 92``
Update OCEANCORES in runscript (make sure the ``aprun`` statement is as expected too)::

  vi runscript
  ...
  OCEANCORES=92

And submit again.

----

2074 timesteps ~ 34 hours with dt=1 min

Not much point using more processors as the tiles are already quite small. Instead
need more walltime.

With dt=60, 1 day = 1440 steps. Run one day on the short queue to see what is
happening with SSH etc.

Do I need to a run with constant T, S? **YES**





Note about timestep
===================

This domain has some very deep water ~7km and large areas with 6km bathymetry.
This should limit the timestep quite considerably:
A deep water wave would travel across a 1/12 degree model cell in
``110000/12 / sqrt(7000*9.8) ~ 35 seconds``...
Consequently a timestep of 1 minute might be as far as I can push it.





Backup to repo key files
========================

::

  cd ~/GitLab/NEMO-RELOC/docs/source
  # DOMANcfg namelist_cfg for domain_cfg.nc (for hybrid z-s coordinates)
  rsync -utv jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/namelist_cfg SEAsia_DOMAINcfg_namelist_cfg

  # EXP namelist_cfg
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP00/namelist_cfg SEAsia_EXP_namelist_cfg

  # PyNEMO namelist.bdy
  rsync -utv jelt@livljobs4:/work/jelt/NEMO/SEAsia/INPUTS/namelist.bdy SEAsia_namelist.bdy

  # Python quick plot of SSH in the output.abort.nc file
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP00/quickplotNEMO.py quickplotNEMO.py


---










Rebuild the output and inspect `rebuild_and_inspect_NEMO_output.rst`_
++++++++++++++++++++++++++++++

---









Next steps
++++++++++

# Tidy up recipe
# Freeze it





---

.. note::

  **TO DO** another time / for Solent config

  * Add in met forcing
  * Turn on T,S at boundares
  * Change namelist to include tidal harmonic analysis::

  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents ('key_diaharm')
  !-----------------------------------------------------------------------
       nit000_han = 1440         ! First time step used for harmonic analysis
       nitend_han = 14400        ! Last time step used for harmonic analysis
       nstep_han  = 15        ! Time step frequency for harmonic analysis
       tname(1)     =   'O1'  !  name of constituent
       tname(2)     =   'P1'
       tname(3)     =   'K1'
       tname(4)     =   'N2'
       tname(5)     =   'M2'
       tname(6)     =   'S2'
       tname(7)     =   'K2'
       tname(8)     =   'Q1'
       tname(9)     =   'M4'

  * Harmonise all wet forcing to use AMM60 data.
