================================================================================
LBay practise prior to Setting up a Solomon Islands and Fiji NEMO configuration
================================================================================

* Build notes on creating 2D/3D tide model covering the Solomon Islands and Fiji.
* Geographical extent is currently 140 degrees W to 185 degrees W and 0 degrees to 30 degrees S
* Following Jeff's Recipe notes for Liverpool Bay and Examples on pynemo website.
* Aim is to run this on Mobius

Recipe Notes
============

22/09/2017
++++++++++

Working within Mobius

Step One:-

Define working and other directory, load relevent modules::

	module purge
	module load shared intel/compiler/64/14.0/2013_sp1.3.174 mvapich2/intel/64/2.0b slurm/14.03.0 cluster-tools/7.0
  export WDIR=/work/$USER/NEMO/Solo
	export INPUTS=/work/$USER/NEMO/INPUTS
	export MOBIUS=/work/thopri/NEMO/Mobius
	export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

CDIR, TDIR and INPUTS do not currently exist. Lets make them!::

.. Tom::

  cd $WDIR
  cp $INPUTS/INPUTS.tar.gz $WDIR
  tar xvfz INPUTS.tar.gz
  rm INPUTS.tar.gz

.. Jeff::
  ln -s /work/thopri/NEMO/INPUTS $INPUTS

Checkout NEMO and XIOS paynote to revision number::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP@5709
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@629

Need to get arch files from Ash's files::

  cd $WDIR/xios-1.0
  cp $MOBIUS/arch* ./arch

Compile XIOS::

 	./make_xios --full --prod --arch mobius_intel  --netcdf_lib netcdf4_par --jobs 6


.. comment::

  I notice that you have two versions of XIOS:
  livljobs4 ARCH $ ls -l  /work/thopri/NEMO/Solo/xios-1.0/bin/xios_server.exe
-rwxr-xr-x 1 thopri pol 11818229 Sep 26 14:38 /work/thopri/NEMO/Solo/xios-1.0/bin/xios_server.exe
livljobs4 ARCH $ ls -l  /work/thopri/NEMO/Mobius/xios-1.0/bin/xios_server.exe
-rwxr-xr-x 1 thopri pol 11852488 Sep 26 09:35 /work/thopri/NEMO/Mobius/xios-1.0/bin/xios_server.exe

  And that your $CDIR/../ARCH/arch-mobius_intel.fcm files points to the older one (in NEMO/Mobius)
  I am just building one. In Solo/.

Step two. Obtain and apply patches::

	export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS
	cd $CDIR/../NEMO/OPA_SRC/SBC
	patch -b < $INPUTS/fldread.patch
	cd ../DOM
	patch -b < $INPUTS/dommsk.patch
	cd ../BDY
	patch -b < $INPUTS/bdyini.patch
	cd $CDIR
	rm $CDIR/../NEMO/OPA_SRC/TRD/trdmod.F90
  #cp $INPUTS/arch-* ../ARCH   # You don't want this. You want the mobius file from Ash::
  cp /scratch/ashbre/NEMO_xios/1arch-mobius_intel.fcm $CDIR/../ARCH/arch-mobius_intel.fcm

25/09/2017
+++++++++++

Ok have got input.gz from Jeff. Will copy into INPUTS directory then untar (is that a word?) and remove tar file. Have ammended recipe from last week to account for new INPUTS file. So will now try and apply patches.

**Patches applied successfully!**

Copy some input files to new configuration path::

  ./makenemo -n Solo -m mobius_intel -j 10 clean

  cp $INPUTS/cpp_LH_REEF.fcm $CDIR/Solo/cpp_Solo.fcm
  cp $INPUTS/dtatsd.F90 $CDIR/Solo/MY_SRC/

Compile NEMO::

	./makenemo -n Solo -m mobius_intel -j 10


26/09/2017
+++++++++++

Not a day I want to relive


27/09/2017
++++++++++

New dawn new day and I am carrying on with Jeff's recipe.

To generate bathymetry, initial conditions and grid information we first need
to compile some of the NEMO TOOLS (after a small bugfix - and to allow direct
passing of arguments). For some reason GRIDGEN doesn’t like INTEL::

  ssh livljobs4

Copy PATHS again::

	export WDIR=/work/$USER/NEMO/Solo
	export INPUTS=/work/$USER/NEMO/INPUTS
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

Apply patches::

	cd $TDIR/WEIGHTS/src
  patch -b < $INPUTS/scripinterp_mod.patch
  patch -b < $INPUTS/scripinterp.patch
  patch -b < $INPUTS/scrip.patch
  patch -b < $INPUTS/scripshape.patch
  patch -b < $INPUTS/scripgrid.patch

Setup for PGI modules and compile::

  cd $TDIR
	cp /work/jelt/NEMO/SEAsia/INPUTS/arch-pgf90_linux_jb.fcm $TDIR/../ARCH/arch-pgf90_linux_jb.fcm

	module add netcdf/gcc/4.1.3
	module add pgi/15.4

  ./maketools -n WEIGHTS -m pgf90_linux_jb
  ./maketools -n REBUILD_NEMO -m pgf90_linux_jb
  ./maketools -n GRIDGEN -m pgf90_linux_jb

Next we use these tools.

1. Generate new coordinates file
++++++++++++++++++++++++++++++++

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements.

First we need to figure out the indices for the new domain, from the parent grid.
Move parent grid into INPUTS::

  #cp $INPUTS/coordinates_ORCA_R12.nc $WDIR/INPUTS/. # Doesn't work for me. As same directory

Inspect this parent coordinates file to define the boundary indices for the new config.

Note, I used FERRET on Livljobs4::

  module load ferret
  FERRET
  use coordinates_ORCA_R12.nc
  shade/i=3385:3392/j=2251:2266 NAV_LAT
  shade/i=3385:3392/j=2251:2266 NAV_LON

Copy namelist file from LH_reef and edit with new indices, retaining use of
ORCA_R12 as course parent grid. (I changed a path somewhere so had to add .. to
``cn_parent_coordinate_file`` path)::

  cd $TDIR/GRIDGEN
  cp $INPUTS/namelist_R12 ./
  vi namelist_R12
  ...
  cn_parent_coordinate_file = '../../../../../INPUTS/coordinates_ORCA_R12.nc'
  ...
  nn_imin = 3385
  nn_imax = 3392
  nn_jmin = 2251
  nn_jmax = 2266
  nn_rhox  = 7
  nn_rhoy = 7

  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe

This generates ``1_coordinates_ORCA_R12.nc``
