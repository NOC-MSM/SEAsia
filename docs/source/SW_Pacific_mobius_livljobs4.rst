=============================================
Setting up a SW Pacific NEMO v4 configuration
=============================================

Machines: livljobs4, MOBIUS

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a SW Pacific regional tide-only model using GEBCO bathymetry, FES tidal
boundaries.

Build on a combination of livljobs4 and MOBIUS.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. MOBIUS: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. MOBIUS: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. MOBIUS: Run simulation

It is a tide only run.

Issues that arose
=================

* ...

.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on MOBIUS but the PyNEMO bits are done on livljobs4.

Starting on MOBIUS::

  ssh mobius

  export CONFIG=SWPacific
  export USER=thopri
  export WORK=/work/thopri/NEMO
  export WDIR=/work/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/
  export EXP=$CDIR/$CONFIG/EXP00

  module purge
  module load shared intel/compiler/64/14.0/2013_sp1.3.174 mvapich2/intel/64/2.0b slurm/14.03.0 cluster-tools/7.0


  Collect essential files
=======================

Note you might have to mkdir the odd directory or two...::
  mkdir $WDIR
  cd $WDIR
  mkdir $START_FILES
  cp $WDIR/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WDIR/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

  cd $WORK
  svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1080
  cd xios-2.0_r1080
  cp $WORK/Mobius/arch* arch/.

Build XIOS::

  ./make_xios --full --prod --arch mobius_intel  --netcdf_lib netcdf4_par --jobs 6

Link the xios-2.0_r1080 to a generic XIOS directory name::

  ln -s  $WORK/xios-2.0_r108  $WORK/XIOS

ARCH file needs to be modifed to point to correct XIOS instance::

  vi $CDIR/../ARCH/arch-mobius_intel.fcm

Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_::
  
  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395
  cp $WORK/Mobius/1arch-mobius_intel.fcm $CDIR/../ARCH/arch-mobius_intel.fcm

Then build::

  cd $CDIR
  ./makenemo -n $CONFIG -m mobius_intel -j 10
  say yes to OPA_SRC no to everything else
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

For the generation of bathymetry: I actually use some old WEIGHTS tools that I
patched and have previously compiled.
I have not reproduced the compilation here (need to keep if the source needs patching
again). If it didn't need patching::

  #./maketools -n WEIGHTS -m XC_ARCHER_INTEL_XIOS1

Otherwise we will use the weights tool in::

  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/






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
