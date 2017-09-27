================================================================================
LBAY practise prior to Setting up a Solomon Islands and Fiji NEMO configuration
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
	export WDIR=/work/thopri/NEMO/Solo
	export INPUTS=/work/thopri/NEMO/INPUTS
	export MOBIUS=/work/thopri/NEMO/Mobius
	export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

CDIR, TDIR and INPUTS do not currently exist. Lets make them!::

  cd $WDIR
  cp $INPUTS/INPUTS.tar.gz $WDIRls
  tar xvfz INPUTS.tar.gz
  rm INPUTS.tar.gz

Checkout NEMO and XIOS paynote to revision number::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP@5709
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@629

Need to get arch files from Ash's files::

  cd $WDIR/xios-1.0
  cp $MOBIUS/arch* ./arch

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
	cp $INPUTS/arch-* ../ARCH

25/09/2017
+++++++++++

Ok have got input.gz from Jeff. Will copy into INPUTS directory then untar (is that a word?) and remove tar file. Have ammended recipe from last week to account for new INPUTS file. So will now try and apply patches.

**Patches applied successfully!**

Copy some input files to new configuration path::

  cp $INPUTS/cpp_LH_REEF.fcm ./Solo/cpp_Solo.fcm
  cp $INPUTS/dtatsd.F90 LBay/MY_SRC/

Compile NEMO::

	./makenemo -n Solo -m mobius_intel -j 10 del_key 'key_lim2'

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

  cp $INPUTS/coordinates_ORCA_R12.nc $WDIR/INPUTS/.

Inspect this parent coordinates file to define the boundary indices for the new config.

Note, I used FERRET on Livljobs4::

  module load ferret
  FERRET
  use coordinates_ORCA_R12.nc
  shade/i=3385:3392/j=2251:2266 NAV_LAT
  shade/i=3385:3392/j=2251:2266 NAV_LON

Copy namelist file from LH_reef and edit with new indices, retaining use of
ORCA_R12 as course
parent grid::

  cd $TDIR/GRIDGEN
  cp $INPUTS/namelist_R12 ./
  vi namelist_R12
  ...
  cn_parent_coordinate_file = '../../../../INPUTS/coordinates_ORCA_R12.nc'
  ...
  nn_imin = 3385
  nn_imax = 3392
  nn_jmin = 2251
  nn_jmax = 2266
  nn_rhox  = 7
  nn_rhoy = 7

  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe

This executable hangs have put terminal output below::

	 *** Reading coordinates file: ../../../../INPUTS/coordinates_ORCA_R12.nc


	 Size of input matrix:
	 (        4322 ;        3059 )

	 Domain:

	                min(1,1:nsizey)             max(nsizex,1:nsizey)
	 longitude:    72.9166658980293       -->    73.0416672469239
	 latitude:    -77.0104751586914       -->    89.9376449584961

	 Domain defined by user:

	                      min                          max
	 longitude:   -3.56828775817968       -->   -2.73814083428448
	 latitude:     53.1041603088379       -->    53.9494680943349

	 Size of domain: REGIONAL

	  ### SUBROUTINE define_domain ###

	  ******************************
	  *** WITHOUT NORTH BOUNDARY ***
	  ******************************

	  ### END SUBROUTINE define_domain ###


	  ### SUBROUTINE define_mixed_grid ###

	  *** CHECKING SIZE OF COARSE DOMAIN ***
	          10 x          18


	 *** SIZE OF MIXED GRID ***
	         140  x          252


	  ### SUBROUTINE write_mixed_grid ###


	  ### END SUBROUTINE write_mixed_grid ###


	  ### END SUBROUTINE define_mixed_grid ###


	 *** SIZE OF FINE GRID ***
	          57  x          113


	 ### SUBROUTINE interp_grid ###


	 *** FUNCTION pol_coef ***


	 *** CHECK LAGRANGE COEFFICIENTS: ***

	 point #           1
	 dlcoef(:)=  -3.790087463556851E-002  0.909620991253644
	  0.151603498542274      -2.332361516034985E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           2
	 dlcoef(:)=  -5.830903790087463E-002  0.787172011661808
	  0.314868804664723      -4.373177842565597E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           3
	 dlcoef(:)=  -6.413994169096210E-002  0.641399416909621
	  0.481049562682216      -5.830903790087463E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           4
	 dlcoef(:)=  -5.830903790087463E-002  0.481049562682216
	  0.641399416909621      -6.413994169096210E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           5
	 dlcoef(:)=  -4.373177842565597E-002  0.314868804664723
	  0.787172011661808      -5.830903790087463E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           6
	 dlcoef(:)=  -2.332361516034985E-002  0.151603498542274
	  0.909620991253644      -3.790087463556852E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 Interpolation along longitude

	 *** FUNCTION pol_coef ***


	 *** CHECK LAGRANGE COEFFICIENTS: ***

	 point #           1
	 dlcoef(:)=  -3.790087463556851E-002  0.909620991253644
	  0.151603498542274      -2.332361516034985E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           2
	 dlcoef(:)=  -5.830903790087463E-002  0.787172011661808
	  0.314868804664723      -4.373177842565597E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           3
	 dlcoef(:)=  -6.413994169096210E-002  0.641399416909621
	  0.481049562682216      -5.830903790087463E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           4
	 dlcoef(:)=  -5.830903790087463E-002  0.481049562682216
	  0.641399416909621      -6.413994169096210E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           5
	 dlcoef(:)=  -4.373177842565597E-002  0.314868804664723
	  0.787172011661808      -5.830903790087463E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 point #           6
	 dlcoef(:)=  -2.332361516034985E-002  0.151603498542274
	  0.909620991253644      -3.790087463556852E-002
	 SUM(dlcoef(:)) =   1.00000000000000

	 Interpolation along latitude
	           0

	 ### END SUBROUTINE interp_grid ###

Only way to end is to close terminal.
