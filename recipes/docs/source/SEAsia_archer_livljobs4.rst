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
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. ARCHER: Run simulation

It is a tide only run.

Issues that arose
=================

* The THREDDS server to access ORCA0083 data is broken and the JASMIN server
 hosting it is not a thredds server. Therefore PyNEMO needs to run on local
  copies of the parent data.

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

  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/par_oce.F90 $START_FILES/.
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/dtatsd.F90 $START_FILES/.

  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharmana.F90 $START_FILES/.
  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step_oce.F90 $CSTART_FILES/.
  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step.F90 $START_FILES/.
  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/bdytides.F90 $START_FILES/.

.. note : jelt 10 May 2018: I think that the link to nico's harmonic analysis files above are out of date.
   I think that Nico updated it, but in the following I am still using his first version which I
   stored in START_FILES.

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

2.5 Generate River forcing
++++++++++++++++++++++++++

`Matlab script to generate river forcing <Generate_river_forcing.rst>`_

Link the river forcing output into ``$INPUTS`` e.g.::

  scp river_test.nc $USER@login.archer.ac.uk:$INPUTS/$CONFIG_rivers.nc



3. Generate a domain configuration file
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



4. Generate initial conditions
==============================


For a new configuration you probably want to start with idealised, or homogenous
initial conditions. This is done with user defined initial conditions ``ln_usr=T``
with the expression being compiled into the executable. (In ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``).

To use initial conditions from an existing T,S field you might need to do a bit
of interpolation. It is advisable to let NEMO do the heavy lifting for vertical
interpolation (rquiring some FORTRAN modifictions), though SOSIE tools can be user
to do simple horizontal interpolation. See notes in

`<generate_initial_conditions.rst>`_

This will generates files typically called ``initcd_votemper.nc``
and ``initcd_vosaline.nc``, and corresponding mask and depth variables files
 ``initcd_mask.nc`` and ``initcd_depth.nc``.




5. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

Initially use zero atm forcing. Specified in usr defined functions in MY_SRC.

Second time around add in met forcing.

Generate cut down drowned precip file (note that the nco tools don't like the
parallel modules). **HEALTH WARNING** *Cut out files with only one index in that lat direction broke NEMO*

**NOTE THAT I'VE LABELLED THE CUTDOWN FILES AS y1979 WHEN THEY ARE 2000. THIS IS TO GET THINGS MOVING AS BCS ARE 1979**::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0
  ncks -d lon,70.,140. -d lat,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_precip_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_precip_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_u10_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_u10_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_v10_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_v10_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_radsw_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_radsw_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_radlw_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_radlw_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_t2_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_t2_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_q2_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_q2_DFS5.1.1_y1979.nc
  ncks -d lon0,70.,140. -d lat0,-21.,25. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_snow_DFS5.1.1_y2000.nc $WDIR/INPUTS/cutdown_drowned_snow_DFS5.1.1_y1979.nc

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
      input_file = 'cutdown_drowned_precip_DFS5.1.1_y1979.nc'

  vi $WDIR/INPUTS/namelist_reshape_bicubic_atmos
  ...
  &grid_inputs
    input_file = 'cutdown_drowned_precip_DFS5.1.1_y1979.nc'


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


*(27 Apr 2018)*
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

Originally, for barotropic forcing::

        ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
        ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
        ln_tra         = .false.               !  boundary conditions for T and S
        ln_ice         = .false.               !  ice boundary condition
        nn_rimwidth    = 1                     !  width of the relaxation zone

Change for baroclinic forcing::

  ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
  ln_tra         = .true.               !  boundary conditions for T and S
  ln_ice         = .false.               !  ice boundary condition
  nn_rimwidth    = 9                    !  width of the relaxation zone

Continuing::

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

  vi runscript_short

  #!/bin/bash
  #PBS -N SEAsia
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk
  #PBS -q short

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  #  echo $(readlink -f $PBS_O_WORKDIR)
  # export OMP_NUM_THREADS=1

  cd $PBS_O_WORKDIR
  #
    echo " ";
    OCEANCORES=88
    XIOSCORES=4
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core
  aprun -b -n $XIOSCORES -N 4 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

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

rn_rdt = 360 and resubmit. Completes in 20min 30days.
**But fills with NaNs from NE Boundary**

---
*(23 Jan 18)* rn_rdt = 120. 7200 steps. 20 mins. Ran 1278 steps in 20mins (~42hrs). STABLE.


Update tides code with Nico's version.
++++++++++++++++++++++++++++++++++++++

Add the POLCOMS harmonic analysis to the executable (as in now in)
`<build_opa_orchestra.rst>`_
This requires some changes to the standard ``namelist_cfg``

Add the final (extra) three variables in your namelist_cfg / nambdy_tide ::

  vi $EXP/namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                    !
                                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.          ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                 ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.false.                    ! Read haronic analyisis from a restart
  /


Edit xml files to output harmonics as amplitudes and phases (e.g.)::

  vi file_def_nemo.xml
  ...
  <file_group id="tidal_harmonics" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1d files -->
    <file id="tidalanalysis.grid_T" name="harmonic_grid_T" description="ocean T grid variables"  enabled=".TRUE.">

      <field field_ref="O1amp"         name="O1amp"       operation="instant" enabled=".TRUE." />
      <field field_ref="O1phase"       name="O1phase"     operation="instant" enabled=".TRUE." />


  vi field_def_nemo-opa.xml
  ...
      <field_group id="Tides_T" grid_ref="grid_T_2D" operation="once" >
      <!-- tidal composante -->
      ...
      <field id="Q1amp"        long_name="Q1 Elevation harmonic Amplitude"                              unit="m"        />
      <field id="Q1phase"      long_name="Q1 Elevation harmonic Phase"                                  unit="degree"   />

*Recall there are elevation, u-vel and v-vel harmonics*. Also editted suffixes
 in velocity fields, adding ``_2D``.


* As before the constituents you want to analyse are set-up in ``nam_diaharm``
 namelist.

* The harmonic analysis is done at the end only as well as the restart dumping
so you can only restart from the last time step so make sure you output the full
 restart at the end. To restart, you just need to turn on the ``ln_harmana_read``
  and to map the files to something like ``restart_harm_ana_*``  as this bit as
   not been developed with a prefix to load the files. You can look at this
    python script if needed:
  ``/work/n01/n01/nibrun/RUNS/SWPacific/SIMU/01_harm_links.py``



Resubmit::

  cd $EXP
  qsub -q short runscript_short


Ran for 20mins. Simulated 45hrs (though I guess it hit the wall limit before
doing the harmonic analysis)



----

*(This short section was done with the old tides routine)*

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
Resulting decomposition (with notes for both standard or short queue
 configurations)::

   vi namelist_cfg
   ...
   !-----------------------------------------------------------------------
   &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
   !-----------------------------------------------------------------------
      cn_mpi_send =  'I'      !  mpi send/recieve type   ='S', 'B', or 'I' for standard send,
                              !  buffer blocking send or immediate non-blocking sends, resp.
      nn_buffer   =   0       !  size in bytes of exported buffer ('B' case), 0 no exportation
      ln_nnogather=  .false.  !  activate code to avoid mpi_allgather use at the northfold
      jpni        =  12    ! standardqueue:12 ! shortqueue:11      !  jpni   number of processors following i (set automatically if < 1)
      jpnj        =  8     !  jpnj   number of processors following j (set automatically if < 1)
      jpnij       =  92    ! standardqueue:92 ! shortqueue:88 !  jpnij  number of local domains (set automatically if < 1)

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



Horizontally constant T and S initial condition
===============================================

Emulate James ORCHESTRA method (first moved usrdef_istate.F90 to usrdef_istate.F90_constTS for safe keeping)::

  cd $CDIR/$CONFIG/MY_SRC
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/usrdef_istate.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/lapack.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/splines.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/types.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/utils.F90 .

Edit usrdef_istate.F90 to put in profile
Data from Hamburg WOCE Climatology Live Access Server at (-2N, 95E).
With constant variable padding below 4000m to make it up to 75 levels::

      zdep(:) = (/     0,    10,    20,    30,    40,    50,   75,    100,   125,   &
                &    150,   175,   200,   250,   300,   350,  400,    500,   600,   &
                &    700,   800,   900,  1000,  1100,  1200,  1300,  1400,  1500,   &
                &   1750,  2000,  2250,  2500,  2750,  3000,  3250,  3500,  3750,   &
                &   4000,  4100,  4200,  4300,  4400,  4500,  4550,  4600,  4700,   &
                &   4800,  4900,  5000,  5100,  5200,  5300,  5400,  5500,  5600,   &
                &   5700,  5800,  5900,  6000,  6100,  6200,  6300,  6400,  6500,   &
                &   6600,  6700,  6800,  6900,  7000,  7100,  7200,  7300,  7400,   &
                &   7500,  7600,  7700 /)

      zsal(:) = (/ 34.05, 34.05, 34.10, 34.13, 34.25, 34.42, 34.88, 35.08, 35.13,   &
                &  35.08, 35.07, 35.06, 35.06, 35.03, 35.01, 34.99, 34.96, 34.97,   &
                &  34.97, 34.95, 34.92, 34.91, 34.88, 34.87, 34.85, 34.83, 34.82,   &
                &  34.80, 34.77, 34.76, 34.75, 34.74, 34.73, 34.73, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72 /)

      ztmp(:) = (/ 28.87, 28.87, 28.87, 28.74, 28.33, 28.01, 25.21, 21.99, 18.51,   &
                &  15.55, 14.39, 13.43, 12.27, 11.48, 11.02, 10.51,  9.58,  8.95,   &
                &   8.35,  7.78,  7.16,  6.52,  5.88,  5.44,  5.02,  4.57,  4.14,   &
                &   3.34,  2.64,  2.31,  2.05,  1.86,  1.69,  1.58,  1.41,  1.23,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15  /)

I fixed the top two layers to be constant salinity and temperature to avoid possible cabling.
The raw data had inversions at 10m

NB avoding mistakes: namely a duplicate depth in the input depth funciton --> spline failure




Turn off tidal forcing
======================

Need to turn off tidal forcing. See ``&nam_tide``::

  ln_tide     = .false.
  ln_tide_pot = .false.    !  use tidal potential forcing


Comment out tidal analysis line ::

  vi ../cpp_SEAsia.fcm
    bld::tool::fppkeys key_zdfgls        \
  //REMOVE//                 key_harm_ana       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Will have to recompile.


Set boudaries to initial condition
++++++++++++++++++++++++++++++++++

(nn_dyn2d_dta was 2)::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .false.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state



Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move the executable to a special name::

  mv $CONFIG/BLD/bin/nemo.exe $CONFIG/BLD/bin/nemo_notide_TSprofile.exe

Move to experiment dir and link executable file in::

  cd $EXP/../EXP-hpg_err
  ln -s $CDIR/$CONFIG/BLD/bin/nemo_notide_TSprofile.exe opa



James spotted that I didn't have lateral diffusion of momentum. Made some changes (following ORCHESTRA namelist_cfg)
Submitted run (EXP01) to test timestep. rn_rdt=360 ran 1304 in 20mins ==> 5.4 days

There an SSH instability in the NE corner when using::

 cn_dyn3d      =  'neumann'

Switch to::

 cn_dyn3d      =  'zerograd'

NB Tried masking, using mask_bdy.nc, but this didn't work.

Cold start run::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  01   !  first time step
     nn_itend    =  7200 ! 30day=7200   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)


Run::

  qsub runscript


Yes. That works.

Though the SSS min decreases by 0.017 over 30 days. A bit odd? Perhaps it is
the non conservative nature of the advective schemes...

NB The effect is coastal. There is no problem around the boundaries.
*(26 March 2018)*

---


Performance note::

  With tides turned off. Initial conditions: T(z),S(z) profiles and u=v=w=0.

  Submit for 10 days dt=60s. nt=14400
  Try 20mins queue. --> 1128 steps.

  Try 30mins with 5day mean output.

  Ran 10days simulation in 2hrs 24mins (rn_rdt=60s, nt=14400, on 92 ocean cores, 120 cpus(total)).


---



HPG errors
++++++++++

Submit a 30 day simulation, from rest, with depth varying spatially homogeneous
temperature and salinity profiles, with no forcing, boundary conditions off:
``ln_bdy = F``


Edit runscript: 2hrs walltime. It took 1h 50mins

Edit namelist_cfg
360s, nt=7200 --> 30 days::

  ln_tide     = .false.
  ln_tide_pot = .false.    !  use tidal potential forcing

cd $EXP/../EXP_hpg_err

Scrape ``umax`` from ``solver.state`` using plot_solver_stat.py

Some along rim currents started but these are small compared to interior currents.
Restart for another 30 days.
After 30 days umax is still growing. Restart run and continue::

  mv solver.stat solver.stat_part1

Check progress with::

   hpg_error_plotNEMO.py
   plot_solver_stat.py

Edit namelist_cfg for restarting::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  7201   !  first time step
     nn_itend    =  14400 ! 30day=7200   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        cn_ocerst_in    = "SEAsia_00007200_restart"   !  suffix of ocean restart name (input)


Resubmit::

  qsub runscript

Ran in 1hr 46

Joing the solver.stat files together::

cp solver.stat solver.stat_part2
cp solver.stat_part1 solver.stat
cat solver.stat_part2 >> solver.stat

module load anaconda
python plot_solver_stat.py


..note::

  In the end I didn't do the restart for the last run I did.



Internal tide with idealised stratification
+++++++++++++++++++++++++++++++++++++++++++
Try with tides turned on.
(Recompile a tide and no tide version. Save in $CONFIG/BLD/bin as nemo_tide.exe
and nemo_notide.exe, then link as appropriate)::

  cd EXP00
  ls -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tide.exe opa

  <    ln_tide     = .true.
  <    ln_tide_pot = .true.    !  use tidal potential forcing
  ...
  <     nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state

  qsub -q short runscript


---



Tide only simulation
++++++++++++++++++++

directory ``EXP_tideonly``

Only tidal forcing. constant T and S
Include: key_harm_ana
EXEC: nemo_tideonyl_TSconst.exe


Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  cp $CONFIG/MY_SRC/usrdef_istate.F90_horizTS $CONFIG/MY_SRC/usrdef_istate.F90

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move executable to something permanent::

  cd $CDIR/$CONFIG/BLD/bin
  mv nemo.exe nemo_tideonly_TSconst.exe

  cd $CDIR/$CONFIG/EXP_tideonly
  ln -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tideonly_TSconst.exe opa

Edit namelist_cfg::

  nn_itend = 14400
  ln_rstart = .false.
  ...
  ln_tide     = .true.
  ln_tide_pot = .true.    !  use tidal potential forcing
  ...
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
      nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state

Didn't bother with the tidal harmonics. It will run but I am spinning up.

Run for 4 hours. nt = 14400, dt =360, 60 days. Completed in 3hr 31.

Edit XML output to produce 5d output.
Resubmit::

  qsub runscript

*(23 Mar 2018)*
Turn on 19 harmonics using the POLCOMS harmonic analysis (Nico's instructions and edits)
Run for another 60 days with harmonic analysis restarting capbability.

Works. Did tidal analysis plots ::

  ~/GitLab/JMMP_tools
  python Tidal_analysis_amplitude.py --verbose
  python Tidal_analysis_plot.py --verbose


Need to continue run

::

  namelist_cfg

  ...
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  28801   !  first time step
     nn_itend    =  43200 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00028800_restart"   !  suffix of ocean restart name (input)

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                   !
                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.                   ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                     ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.true.                      ! Read haronic analyisis from a restart



   !-----------------------------------------------------------------------
   &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
   !-----------------------------------------------------------------------
       nit000_han = 28801         ! First time step used for harmonic analysis
       nitend_han = 43200 ! 1440 !      ! Last time step used for harmonic analysis



Submit::

  qsub runscript

Took 3:36 mins. The tidal analysis (https://www.evernote.com/shard/s652/nl/85824826/674115d9-29be-480a-ba71-6814de98df4b/) doesn't show significant improvement
::

  ~/GitLab/JMMP_tools
  python Tidal_analysis_amplitude.py --verbose
  python Tidal_analysis_plot.py --verbose

But 4 months of simulation might still be on the short side. Run for another two months.
::

  vi namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  43201   !  first time step
     nn_itend    =  57600 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00043200_restart"   !  suffix of ocean restart name (input)




   !-----------------------------------------------------------------------
   &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
   !-----------------------------------------------------------------------
       nit000_han = 43201         ! First time step used for harmonic analysis
       nitend_han = 57600 ! 1440 !      ! Last time step used for harmonic analysis



Submit::

  qsub runscript

Ran with minor change to statistcs.



Tide with idealised stratification and rivers
+++++++++++++++++++++++++++++++++++++++++++++

directory ``EXP_tide_TSprofile_river``

Only tidal forcing. profile of T and S
Include: key_harm_ana
EXEC: nemo_tide_TSprofile_nomet.exe

::

  cp EXP_tideonly/* EXP_tide_TSprofile_river/.
  ln -s $INPUTS bdydta

Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  cp $CONFIG/MY_SRC/usrdef_istate.F90_horizTS $CONFIG/MY_SRC/usrdef_istate.F90

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move executable to something permanent::

  cd $CDIR/$CONFIG/BLD/bin
  mv nemo.exe nemo_tide_TSprofile_nomet.exe

  cd $CDIR/$CONFIG/EXP_tide_TSprofile_river
  ln -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tide_TSprofile_nomet.exe opa

Edit namelist_cfg::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  57601   !  first time step
     nn_itend    =  58560 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00057600_restart"   !  suffix of ocean restart name (input)

  ...

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing
  ...

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
      nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state
      ...

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                   !
                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.                   ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                     ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.true.                      ! Read haronic analyisis from a restart

  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 57601         ! First time step used for harmonic analysis
      nitend_han = 58560 ! 1440 !      ! Last time step used for harmonic analysis




Finally turn rivers on::

  !-----------------------------------------------------------------------
  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
  ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)
  ...
  !-----------------------------------------------------------------------
  &namsbc_rnf    !   runoffs namelist surface boundary condition          (ln_rnf=T)
  !-----------------------------------------------------------------------
  !              !  file name           ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                      !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_rnf      = 'SEAsia_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_cnf      = 'SEAsia_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''
     ...
     cn_dir      = 'bdydta/'      !  root directory for the location of the runoff files
     ln_rnf_mouth= .false.    !  specific treatment at rivers mouths


Edit runscript::

  vi runscript
  #PBS -l walltime=00:30:00

Resubmit::

  qsub runscript

**PENDING**
*(27 Apr 2018)*

/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_tide_TSprofile_river

Run for 30 mins. nt = 960, dt =360, 4 days. Completed in  15 mins 06s.

CHECK OUTPUT - This works, but 4 days is not very exciting.



Initial conditions
++++++++++++++++++

directory: EXP_hpg_initcd

Switch in initial conditions from the existing (working) tide only homogeneous run.
::
  mkdir EXP_hpg_initcd

Copy files from ??? to EXP_hpg_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_tideonly/* EXP_hpg_initcd/


Try from rest.


Rerun the last successful tidal only simulation but with initial conditions.
Switch off bcs and tides run for 60 days::


  vi namelist_cfg

  ...
  cn_exp      =    "SEAsia"  !  experience name
  nn_it000    =  1   !  first time step
  nn_itend    =  7200 ! 30day=7200; dt=360  !  last  time step (std 5475)

  ln_restart = .false.

Turn off open boundaries::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .false.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state


Turn off tides::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .false.
     ln_tide_pot = .false.    !  use tidal potential forcing

Make sure the harmonic restart is also off (not so important yet since ln_tide=F)::

  ln_harmana_read=.false.                      ! Read haronic analyisis from a restart
  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 1         ! First time step used for harmonic analysis
      nitend_han = 2400 ! 1440 !      ! Last time step used for harmonic analysis

Turn on intial conditions::

  !-----------------------------------------------------------------------
  &namtsd        !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  !              !  file name                 ! frequency (hours) ! variable ! time interp.!  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                            !  (if <0  months)  !   name   !  (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
  sn_tem  = 'votemper_ORCA0083-N01-SEAsia_1978.nc',         -12        ,'votemper' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_sal  = 'vosaline_ORCA0083-N01-SEAsia_1978.nc',         -12        ,'vosaline' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_dep  = 'initcd_depth.nc'   ,         -12        ,'gdept_4D',   .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''
  sn_msk  = 'initcd_mask.nc'    ,         -12        ,'mask',       .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''

    !
     cn_dir        = '../../../../INPUTS/'     !  root directory for the location of the runoff files
     ln_tsd_init   = .true.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_interp = .true.    !  Interpolation of T & S in the verticalinput data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)


Add some daily output::

  vi file_def_nemo.xml

    <file_group id="1d" output_freq="1d"  output_level="10" enabled=".TRUE.">  <!-- 1d files -->
    <file id="file19" name_suffix="_SS" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"   />
     <field field_ref="sss"          name="sss"   />
    </file>
  </file_group>


Link in missing sybolic link files::

  rm opa
  ln -s ../BLD/bin/nemo_tide_nomet.exe opa

I also fudged the dates on the boundary conditions files::

  cd $INPUTS
  ln -s SEAsia_bdyT_y1979m11.nc SEAsia_bdyT_y2000m01.nc
  ln -s SEAsia_bdyU_y1979m11.nc SEAsia_bdyU_y2000m01.nc
  ln -s SEAsia_bdyV_y1979m11.nc SEAsia_bdyV_y2000m01.nc
  ln -s SEAsia_bt_bdyT_y1979m11.nc SEAsia_bt_bdyT_y2000m01.nc

Run on short queue::

  cd SEAsia/EXP_tide_initcd
  qsub runscript_short

Later in evening run on standard queue (2hrs)::

  cd SEAsia/EXP_tide_initcd
  qsub runscript


**THIS WORKS** 13 May 2018



Initial conditions + tides
--------------------------

directory: EXP_tide_initcd

EXEC: ``nemo_tide_nomet.exe``

Switch in initial conditions from the existing (working) tide only homogeneous run.
::
  mkdir EXP_tide_initcd

Copy files from ??? to EXP_tide_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_hpg_initcd/* EXP_tide_initcd/


Try from rest.::

  vi namelist_cfg

  ...
  cn_exp      =    "SEAsia"  !  experience name
  nn_it000    =  1   !  first time step
  nn_itend    =  7200 ! 30day=7200; dt=360  !  last  time step (std 5475)

  ln_restart = .false.

Turn off open boundaries::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .false.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state


Turn on tides::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing

Make sure the harmonic restart is also off ::

  ln_harmana_read=.false.                      ! Read haronic analyisis from a restart
  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 1         ! First time step used for harmonic analysis
      nitend_han = 2400 ! 1440 !      ! Last time step used for harmonic analysis



cd /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_tide_initcd

**PENDING** 5324836.sdb
13 May 2018

completed 1hr 39




Try lateral boundary conditions T,S,u
++++++++++++++++++++++++++++++++++++++


**Next steps**

* lateral boundary conditions (T,S,u)
* tides + lateral boundary conditions


---

*(24 Apr 2018)*
Build the 3D boundary conditions from ORCA0083-N01
==================================================

**The following gets a bit messy towards the end of this section and needs tidying up. But it seems to work**

Use PyNEMO to generate 3D bcs on livljobs4::

  ssh livljobs4

  . ~/temporary_path_names_for_NEMO_build
  cd $INPUTS

Edit the namelist.bdy for 3D fields. There are a few key things to note:

* Turn off the tides ``ln_tide =.false.`` and change ``rimwidth``: 1 --> 9
(tides don't work with rimwidth != 1 and baroclinc bcs don't work with rimwidth=1)

* For the ORCA0083-N001 run there was a glitch with the model timestamp. It is
 out by 36 days and 16 hours. This can be adjusted with a new namelist variable
 ``nn_src_time_adj``.

* The parent files could not be accessed on the thredd server because it is broken.
 For now I've copied Nov 1978 into ``/projectsa/accord/ORCA0083``. This is reflected
 in local_inputs_src.ncml

* Make sure the namelist time variables (years and months) match the timestamps in the files you
 load otherwise not much happens...

* You might need to have the time data spanning a whole month, so that a month
 can be extracted.

* Finally, there may have been an issue with a pynemo datetime function that
 determined the datetime for end of the month in profile.py
I  made change to the master branch to use datetime and timedelta instead of
 netcdfdatime, or something. This only matters is pynemo breaks with an
associated datetime error. Then this is the fix.

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
     sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
     sn_src_zgr = './mesh_zgr_src.nc'   !  parent
     sn_dst_hgr = './domain_cfg.nc'
     sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
     sn_src_msk = './mask_src.nc'       ! parent
     sn_bathy   = './bathy_meter.nc'

  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
     sn_src_dir = './local_inputs_src.ncml'       ! src_files/'
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
      ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
      ln_tra         = .true.               !  boundary conditions for T and S
      ln_ice         = .false.               !  ice boundary condition
      nn_rimwidth    = 9                    !  width of the relaxation zone

  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
      ln_tide        = .false.               !  =T : produce bdy tidal conditions
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
      nn_month_end    = 11         !  month end (default = 12 is years>1)
      sn_dst_calendar = 'gregorian' !  output calendar format
      nn_base_year    = 1978        !  base year for time counter
      sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
      nn_src_time_adj    = -3254400  ! -3168000 - 86400 ! fix to align model time stamp
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


The ORCA0083-N01 parent files. (This could be streamlined with some wildcards)::

  vi local_inputs_src.ncml::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>


Run pynemo (on branch ORCA0083)::


  cd $INPUTS

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  pynemo -s namelist.bdy


*(25 Apr 2018)*

* ``nemo_bdy_extr_tm3.py`` Editted ind_NAN variable and caught exception for
empty index list. Resubmitted.

Outputs::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc

This took 2 hours! To produce a month of boundary conditions.


Copy files from SAN to ARCHER::

  livljobs4
  cd $INPUTS
  for file in  SEAsia_b*nc; do rsync -uvt $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done





Tide + 2D forcing simulation
++++++++++++++++++++++++++++

directory ``EXP_openbcs``

EXEC: ``nemo_tide_nomet.exe``


Switch in initial conditions from the existing (working) tide only homogeneous
 run::

  mkdir EXP_openbcs
  cd EXP_openbcs

Copy files from ??? to EXP_tide_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=$CONFIG*.nc EXP_tide_initcd/* EXP_openbcs/

Link in missing sybolic link files::

  rm opa
  ln -s ../BLD/bin/nemo_tide_nomet.exe opa

Try from rest.
First try tides only with ln_bdy=T. ``nn_dyn2d_dta   =  2``. Ran to 400 steps (40hrs) OK
Second try with 2D bcs not tides: ``nn_dyn2d_dta   =  1``.

===>>> : E R R O R
        ===========

 stpctl: the speed is larger than 20 m/s
 ======
kt=     1 max abs(U):  1.0000E+20, i j k:   288    2    1

          output of last fields in numwso

--- **PROBLEM WITH 2D BCS**

Third try with tides + tracers::

  cn_dyn2d       = 'flather'               !
  nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
  cn_dyn3d      =  'zerograd'               !
  nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
  cn_tra        =  'frs'               !
  nn_tra_dta    =  1                    !  = 0, bdy data are equal to the initial state

Also blows immediately.

Revisit PyNEMO. Turn off tidal output and run.

----

trouble
+++++++

19 May 2018

Make a new mask file to only process the southern boundary.
Turn off vertical interpolation in PyNEMO so that it creates output on the parent
vertical grid.

livljobs4/6


The mask variable takes values (-1 mask, 1 wet/active domain, 0 land). Need to
only mask a single point around the edge since the rimwidth is considered to be
part of the active domain.

PyNEMO looks for the interface between -1 and 1 to generate boundary forcings. Get a
template from domain_cfg.nc and then modify as desired around the boundary.

Maked the entire boundary "land=0", except for the wet bits along the southern boundary
which are "mask=-1"::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython::

import netCDF4
import numpy as np
dset = netCDF4.Dataset('bdy_mask2.nc','a')
[ny,nx] = np.shape(dset.variables['mask'][:])
for i in range(nx):
  if dset.variables['mask'][1,i] == 1:
    dset.variables['mask'][0,i] = -1
  else:
    dset.variables['mask'][0,i] = 0
dset.close()

Rebuld pynemo as necessary `<install_nrct.rst>`_ with any fixes.

Rebuild boundary conditions::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

----

**21 May 2018**

Get James' to generate open bcs and try them out.
Also created a new domain_cfg.nc file.
From /work/jdha/jelt copy to $INPUTS::

  cd /work/jdha/jelt/
  rsync -uvat  SEAsia_full_bdyT_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bt_bdyT_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdyU_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdyV_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdytide*nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  domain_cfg_jdha.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.

Also has a bdy_mask.nc and a domain_cfg.nc that were the same as mine.

Edit the namlist for accomdate new files::

  cd $EXP/../EXP_openbcs
  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
        !                   !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
        cn_domcfg = "domain_cfg_jdha"         ! domain configuration filename

  ...

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_full_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_full_bdyU',  24        , 'vobtcrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_full_bdyV',  24        , 'vobtcrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_full_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_full_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_full_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_full_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''


 !-----------------------------------------------------------------------
 &nambdy_tide   !  tidal forcing at open boundaries
 !-----------------------------------------------------------------------
    filtide      = 'bdydta/SEAsia_full_bdytide_'         !  file name root of tidal forcing files





Catch up with James' edits
++++++++++++++++++++++++++
*4 June 2018*

Add Git repo of MY_SRC and do a "sparse checkout" of the MY_SRC folder::

  mv $CDIR/$CONFIG/MY_SRC $CDIR/$CONFIG/MY_SRC/MY_SRC_safe
  cd $WORK/$USER/$CONFIG
  mkdir git_repo
  cd git_repo
  git init
  git remote add -f origin https://github.com/NOC-MSM/NEMO_cfgs.git  # Creates and empty repo without checking out the features
  git config core.sparseCheckout true

Specify the folders I want::

  echo $CONFIG/MY_SRC >> .git/info/sparse-checkout

Last but not least, update your empty repo with the state from the remote. (And
set the upstream remote master)::

  git pull origin master
  git push --set-upstream origin master

Symbolic link to where MY_SRC actually sits::

  ln -s $WORK/$USER/$CONFIG/git_repo/$CONFIG/MY_SRC $CDIR/$CONFIG/MY_SRC


.. note : the above git integration needs to be abstracted to the setup phase.

Rebuild opa (I removed the tidal analysis key from the cpp file, as we aren't
using it yet) At this time also removed Nico's harmonic analysis changes::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Link in executable::

  cd $EXP/../EXP_openbcs
  rm opa
  ln -s ../BLD/bin/nemo.exe opa

Edit namelist_cfg for boundary conditions. Use full velocities with FRS for tracers
and specified for 3d velcities::


  vi namelist_cfg
  ...
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
      nn_dyn2d_dta   =  3                  !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'specified'               !
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
      ln_dyn3d_dmp  =.true.                !  open boundary condition for baroclinic velocities
      rn_time_dmp   =  1.                   ! Damping time scale in days
      rn_time_dmp_out =  150.                 ! Outflow damping time scale
      nn_rimwidth   = 9                    !  width of the relaxation zone
      ln_vol        = .false.               !  total volume correction (see nn_volctl parameter)
      nn_volctl     = 1                     !  = 0, the total water flux across open boundaries is zero
      nb_jpk_bdy    =  75                   ! number of levels in the bdy data (set < 0 if consistent with planned run)
  /

Presently need to use James' forcing files::

  ...
  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_full_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_full2_bdyU',  24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_full2_bdyV',  24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_full2_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_full2_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_full_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_full_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  ! for lim2
  !   bn_frld    = 'amm12_bdyT_ice',         24        , 'ileadfra',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_hicif   = 'amm12_bdyT_ice',         24        , 'iicethic',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_hsnif   = 'amm12_bdyT_ice',         24        , 'isnowthi',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  ! for lim3
  !   bn_a_i     = 'amm12_bdyT_ice',         24        , 'ileadfra',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_ht_i    = 'amm12_bdyT_ice',         24        , 'iicethic',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_ht_s    = 'amm12_bdyT_ice',         24        , 'isnowthi',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''

     cn_dir      = 'bdydta/' !  root directory for the location of the bulk files
     ln_full_vel = .true.   !
  /

Submit for 3 hours (nt=7200; 30 days) and see what happens::

  qsub runscript
  5388562.sdb

*11 June 18* Got a lead from Andrew C about XIOS2. Sumbmitted and waiting
**PENDING**
It should output a sub domain...
Hmm it doesn't seem to output anything.
Perhaps start again with fresh XML files.


---

*19 June 2018*

Check that my PyNEMO implementation matches James' output.
::

  livljobs6:
  cd $INPUTS
  vi namelist.bdy
  ...
  ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
  ln_tra         = .true.               !  boundary conditions for T and S
  ...
  ln_tide = .false.


Rebuild PyNEMO for 3d fields (branch ORCA0083).
run time 1729.83568501 = 28 mins

Output::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc


The rimwidth doesn't get picked up for the u2d and u3d files. Though it does
 for the T,S. Try again using the master branch.
 Execution Time: 1601.68916202::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc

Not sure about the size of the velocity arrays. They don't use rimwidth != 1 ..


*(3 Oct 18)*
Rebuild PyNEMO for 3d fields (branch ORCA0083). And execute (namelist.bdy as above)
Build in environment ```nrct_env2``

Try building tides separately ``ln_tides=.False.``. Otherwise problem with the rimwidth variable.
Had to edit ``inputs_dst.ncml``  to remove e3t_0 and e3w_0 renaming.

(30mins)::

  pynemo -s namelist.bdy

Generates::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc



.. also had to regenerate the tidal files and copy.
 Caution the coordinates.bdy.nc for tides only should not be used as one needs
 to set rimwidth=1 to get it to work

Copy files from SAN to ARCHER::

  livljobs4
  cd $INPUTS
  for file in  SEAsia_b*nc; do rsync -uvt $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done
  rsync -uvt coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.

Also copy coordinates.bdy.nc to EXP directory

Move to ARCHER::

  cd $EXP/../EXP_openbcs

Remove links to *_full* in forcing files

vi namelist_cfg::

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_bdyU',  24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_bdyV',  24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''


     ...

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files

This means that NEMO will take the 3d velocities and construct its own depth averaged
and perturbation components.

Restore missing namelist_cfg.nc file.
Rsync and make a symbolic link from INPUTS::

  livljobs6 INPUTS $ rsync -uvt domain_cfg.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  ln -s /work/n01/n01/jelt/SEAsia/INPUTS/domain_cfg.nc domain_cfg.nc

PENDING

Submit (nt=240) and see what happens::

 qsub runscript_short
 5649089.sdb


**The above needs tidying up. But it seems to work**


The full ocean: Initial conditions, tides, rivers, open boundaries - no met
===========================================================================

directory: EXP_fullocean

EXEC: ``nemo.exe``

::
  mkdir EXP_fullocean

Copy files from EXP_tide_initcd to EXP_fullocean **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_tide_initcd/* EXP_fullocean/

Need to add in:

* open boundary conditions
* rivers

Get the executable in place (I think that the standard nemo.exe is OK)::

  cd EXP_fullocean
  rm opa
  ln -s ../BLD/bin/nemo.exe opa


Add rivers
++++++++++

::

  cp ../EXP_openbcs/namelist_cfg namelist_cfg
  vi namelist_cfg

  ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)

  sn_rnf      = 'SEAsia_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
  sn_cnf      = 'SEAsia_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''

  cn_dir      = 'bdydta/'      !  root directory for the location of the runoff files
  ln_rnf_mouth= .false.    !  specific treatment at rivers mouths

Run on the standard queue::

  jpni        =  12    ! standardqueue:12 ! shortqueue:11      !  jpni   number of processors following i (set automatically if < 1)
  jpnj        =  8     !  jpnj   number of processors following j (set automatically if < 1)
  jpnij       =  92    ! standardqueue:92 ! shortqueue:88 !  jpnij  number of local domains (set automatically if < 1)

Submit::

  qsub runscript

NB I don't have Nico's tide mods in. They need to be restored.
30 days ran in 1h44.


Add in FES tides. Tidal Potential and harmonic analyisis
========================================================

Following `<FES2014_Nemo.rst>`_.
Changes to use FES2014 tides::

  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_FES14.h90 $CDIR/$CONFIG/MY_SRC/.

Compile with flag *key_FES14_tides*

Long period tidal potential and variable Love number adjustments (latter including a new namelist variable)::

  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/sbctide.F90 $CDIR/$CONFIG/MY_SRC/.
  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_mod.F90 $CDIR/$CONFIG/MY_SRC/.
  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tideini.F90 $CDIR/$CONFIG/MY_SRC/.

For restarts, 3d and fast harmonic analysis::

  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm_fast.F90 $CDIR/$CONFIG/MY_SRC/.

with new compilation key *key_diaharm_fast*

The resulting compilation file is now::

  vi $CDIR/$CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                key_mpp_mpi       \
                key_iomput        \
                key_nosignedzero  \
                key_FES14_tides   \
                key_diaharm_fast


Namelist changes::

  !-----------------------------------------------------------------------
  &nam_diaharm_fast   !   Harmonic analysis of tidal constituents               ("key_diaharm_fast")
  !-----------------------------------------------------------------------
  ln_diaharm_store = .true.
  ln_diaharm_compute = .true.
  ln_diaharm_read_restart = .false.
  ln_ana_ssh   = .true.
  ln_ana_uvbar = .false.
  ln_ana_bfric = .false.
  ln_ana_rho  = .false.
  ln_ana_uv3d = .false.
  ln_ana_w3d  = .false.
  tname(1) = 'O1',
  tname(2) = 'M2',
  /
  ...
  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing
     ln_tide_ramp= .true.   !
     rdttideramp =    1.     !
     dn_love_number = 0.69!   clname(1)   =  'M2'   !  name of constituent - all tidal components must be set in namelist_cfg
   ...

Changes to the XML file: ``file_def_nemo.xml`` and ``field_def_nemo-opa.xml`` which I already had


Move the working executable *BLD/bin/nemo.exe --> nemo_8Oct18.exe*. Build::
./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10


Submitted run to check if it works. PENDING.





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





Note about timestep
===================

This domain has some very deep water ~7km and large areas with 6km bathymetry.
This should limit the timestep quite considerably:
A deep water wave would travel across a 1/12 degree model cell in
``110000/12 / sqrt(7000*9.8) ~ 35 seconds``...
Consequently a timestep of 1 minute might be as far as I can push it.
---




Rebuild the output and inspect `rebuild_and_inspect_NEMO_output.rst`_
++++++++++++++++++++++++++++++

::
  qsub -q short rebuild.pbs

---
