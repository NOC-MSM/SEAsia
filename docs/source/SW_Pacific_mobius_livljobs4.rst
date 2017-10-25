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
  #export USER=thopri
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
  cp $WORK/../LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WORK/../LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

  cd $WORK
  svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1080
  cd xios-2.0_r1080
  cp $WORK/Mobius/arch* arch/.

Build XIOS::

  ./make_xios --full --prod --arch mobius_intel  --netcdf_lib netcdf4_par

Commented out the ``--job 8`` command, just incase.
Results in error with building::

  ...
  touch /work/jelt/NEMO/xios-2.0_r1080/flags/LDFLAGS__test.flags
  touch /work/jelt/NEMO/xios-2.0_r1080/flags/LDFLAGS__test__test_remap.flags
  touch /work/jelt/NEMO/xios-2.0_r1080/done/mod_wait.done
  fcm_internal load:F test test_remap.o test_remap.exe
  ar: creating /work/jelt/NEMO/xios-2.0_r1080/tmp/lib__fcm__test_remap.a
  mpif90 -nofor-main -o test_remap.exe /work/jelt/NEMO/xios-2.0_r1080/obj/test_remap.o -L/work/jelt/NEMO/xios-2.0_r1080/lib -l__fcm__test_remap -Wl,"--allow-multiple-definition" -L/login/jdha/utils/netcdf_mob_intel/lib -L/login/jdha/utils/hdf5_mob_intel/lib   -lnetcdf -lnetcdff  -lhdf5_hl -lhdf5 -lz -lstdc++
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_genatt.o): In function `nf_copy_att':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_genatt.f90:235: undefined reference to `nc_copy_att'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_genatt.o): In function `nf_rename_att':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_genatt.f90:270: undefined reference to `nc_rename_att'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_genatt.o): In function `nf_del_att':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_genatt.f90:300: undefined reference to `nc_del_att'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_genvar.o): In function `nf_copy_var':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_genvar.f90:380: undefined reference to `nc_copy_var'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_def_enum':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1161: undefined reference to `nc_def_enum'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_insert_enum':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1200: undefined reference to `nc_insert_enum'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_inq_enum':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1232: undefined reference to `nc_inq_enum'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_inq_enum_member':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1274: undefined reference to `nc_inq_enum_member'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_inq_enum_ident':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1309: undefined reference to `nc_inq_enum_ident'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_def_opaque':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1344: undefined reference to `nc_def_opaque'
  /login/jdha/utils/netcdf_mob_intel/lib/libnetcdff.a(nf_nc4.o): In function `nf_inq_opaque':
  /work/jdha/TOSORT/utils/netcdf-fortran/fortran/nf_nc4.f90:1377: undefined reference to `nc_inq_opaque'
  fcm_internal load failed (256)
  gmake: *** [test_remap.exe] Error 1
  gmake -f /work/jelt/NEMO/xios-2.0_r1080/Makefile -j 1 all failed (2) at /work/jelt/NEMO/xios-2.0_r1080/tools/FCM/bin/../lib/Fcm/Build.pm line 597
  ->Make: 854 seconds
  ->TOTAL: 878 seconds
  Build failed on Mon Oct 23 10:16:37 2017.

*(END OF JEFF'S INVESTIGATION)*
... have looked at Jeff's arch files for ARCHER and tried modifiying some bits of my configuration file that could be the issue with no luck. However I have found a webpage on the NEMO site that has identified a similar problem::

  https://forge.ipsl.jussieu.fr/orchidee/wiki/DevelopmentActivities/ORCHIDEE-MICT-IMBALANCE-P/knownissues

This seems to be similar so have tried the fix by modifying the make_xios file and the arch*.path file.

I changed a line in the makenemo::

  From: XIOS_LIB="$XIOS_LIB $NETCDF_LIBDIR $HDF5_LIBDIR $MPI_LIBDIR $NETCDF_LIB $HDF5_LIB $MPI_LIB"
  To:   XIOS_LIB="$XIOS_LIB $NETCDF_LIBDIR $HDF5_LIBDIR $MPI_LIBDIR $NETCDF_LIB $HDF5_LIB -L$MPI_LIB -lnetcdff"

And changes the HDF5_LIB line in the .path file too::

  HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz -lcurl"

However this produced an error with lcurl so was removed resulting in ::

  HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz"

Still doesn't work. Hmmmmmmm sure its something simple but what? Current error message is as follows::

  fcm_internal load failed (256)
  gmake: *** [test_remap.exe] Error 1
  gmake -f /work/thopri/NEMO/xios-2.0_r1080/Makefile -j 1 all failed (2) at /work/thopri/NEMO/xios-2.0_r1080/tools/FCM/bin/../lib/Fcm/Build.pm line 597
  ->Make: 867 seconds
  ->TOTAL: 891 seconds
  Build failed on Fri Oct 20 15:14:40 2017.

Have found James Harle documentation on how to run NEMO on mobius which is helpful if I need to remake arch files. but am trying to build XIOS again, have downloaded NEMO code as well as it may be required to build XIOS (downloading NEMO comes first in many guides).

I think the issue is with the libraries that James has put together, they work for XIOS1 but the errors relate to undefined references in these libraries so I don't think they are compatable with XIOS2.







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
  ./makenemo -n $CONFIG -m mobius_intel -j 10 clean
  say yes to OPA_SRC no to everything else
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
