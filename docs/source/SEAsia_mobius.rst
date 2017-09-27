=============================================================================
Grid generation and inputs files using Liverpool machine. Application SE Asia
=============================================================================

* Port of LBay.rst that was written for ARCHER
* Concerned here with generating grids, bathymetry and forcing files for a new config
* Aim is to run this anywhere (ARCHER), but generate setup locally.

Recipe Notes
============

27/09/2017
++++++++++

Working on Mobius::

  ssh mobius

Step One:-

Define working and other directory, load relevent modules::

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
	export INPUTS=/work/$USER/NEMO/INPUTS
	#export MOBIUS=/work/thopri/NEMO/Mobius
	export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  module purge
	module load shared intel/compiler/64/14.0/2013_sp1.3.174 mvapich2/intel/64/2.0b slurm/14.03.0 cluster-tools/7.0

CDIR, TDIR and INPUTS do not currently exist. Lets make them!::

  mkdir $WDIR
  mkdir $INPUTS

..
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

Need to get arch files. NB These point to jdha utils paths::

  cd $WDIR/xios-1.0
  cp $WDIR/../ARCH/arch* ./arch


.. Do I need this?

      Compile XIOS::

       	./make_xios --full --prod --arch mobius_intel  --netcdf_lib netcdf4_par --jobs 6




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
  cp $WDIR/../ARCH/1arch-mobius_intel.fcm $CDIR/../ARCH/arch-mobius_intel.fcm

Edit path in arch file::

  vi $CDIR/../ARCH/arch-mobius_intel.fcm
  ...
  %XIOS_HOME           /work/jelt/NEMO/SEAsia/xios-1.0
  ...


25/09/2017
+++++++++++

Ok have got input.gz from Jeff. Will copy into INPUTS directory then untar (is that a word?) and remove tar file. Have ammended recipe from last week to account for new INPUTS file. So will now try and apply patches.

**Patches applied successfully!**

Create new config. Select only OPA_SRC option::

  ./makenemo -n $CONFIG -m mobius_intel -j 10 clean

Create / Edit new cpp keys file::

  echo "bld::tool::fppkeys   key_dynspg_ts key_ldfslp key_zdfgls key_vvl key_mpp_mpi key_netcdf4 key_nosignedzero key_iomput key_gen_IC key_bdy" > $CDIR/$CONFIG/cpp_$CONFIG.fcm


**GOT HERE**

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

  cd $WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/WEIGHTS/src
  patch -b < $INPUTS/scripinterp_mod.patch
  patch -b < $INPUTS/scripinterp.patch
  patch -b < $INPUTS/scrip.patch
  patch -b < $INPUTS/scripshape.patch
  patch -b < $INPUTS/scripgrid.patch

  cd ../../
  ./maketools -n WEIGHTS -m mobius_intel
  ./maketools -n REBUILD_NEMO -m mobius_intel

  module load netcdf hdf5
  ./maketools -n GRIDGEN -m mobius_intel

Need to take a more structured approach to setting up this new configuration

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
