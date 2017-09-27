=============================================================================
Grid generation and inputs files using Liverpool machine. Application SE Asia
=============================================================================

* Port of LBay.rst, which was written for ARCHER. This is for liverpool machines
* Concerned here with generating grids, bathymetry and forcing files for a new config
* Aim is to run this anywhere (ARCHER), but generate setup locally.

Recipe Notes
============

*(27/09/2017)*

Build NEMO on Mobius::

  ssh mobius

Step One:-

Define working and other directory, load relevent modules::

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
  export START_FILES=$WDIR/START_FILES # generic stuff for making more stuff. Mostly code.
	export INPUTS=$WDIR/INPUTS         # config specific stuff that gets made and is for running NEMO
  export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  module purge
	module load shared intel/compiler/64/14.0/2013_sp1.3.174 mvapich2/intel/64/2.0b slurm/14.03.0 cluster-tools/7.0

CDIR, TDIR and INPUTS do not currently exist. Lets make them!::

  mkdir $WDIR
  mkdir $START_FILES # Untar start_files.tar if it is complete

..
      .. Tom::

        cd $WDIR
        cp $START_FILES/INPUTS.tar.gz $WDIR
        tar xvfz INPUTS.tar.gz
        rm INPUTS.tar.gz

      .. Jeff::
        ln -s /work/thopri/NEMO/INPUTS $START_FILES

        cp /work/thopri/NEMO/J_INPUTS/*patch $START_FILES/.

Checkout NEMO and XIOS paynote to revision number::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP@5709
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@629

Need to get arch files. NB These point to jdha utils paths::

  cd $WDIR/xios-1.0
  cp $START_FILES/arch* ./arch


I think XIOS is needed to make NEMO run, which I need to generate mesh files.
Compile XIOS::

  ./make_xios --full --prod --arch mobius_intel  --netcdf_lib netcdf4_par --jobs 6




Step two. Obtain and apply patches::

	export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS
	cd $CDIR/../NEMO/OPA_SRC/SBC
	patch -b < $START_FILES/fldread.patch
	cd ../DOM
	patch -b < $START_FILES/dommsk.patch
	cd ../BDY
	patch -b < $START_FILES/bdyini.patch
	cd $CDIR
	rm $CDIR/../NEMO/OPA_SRC/TRD/trdmod.F90
  cp $START_FILES/1arch-mobius_intel.fcm $CDIR/../ARCH/arch-mobius_intel.fcm

Edit XIOS path in arch file::

  vi $CDIR/../ARCH/arch-mobius_intel.fcm
  ...
  %XIOS_HOME           /work/jelt/NEMO/SEAsia/xios-1.0
  ...


Blind bake a fresh NEMO config
++++++++++++++++++++++++++++++

Create new config. Select only OPA_SRC option::

  ./makenemo -n $CONFIG -m mobius_intel -j 10 clean

Create / Edit new cpp keys file::

  echo "bld::tool::fppkeys   key_dynspg_ts key_ldfslp key_zdfgls key_vvl key_mpp_mpi key_netcdf4 key_nosignedzero key_iomput key_gen_IC key_bdy" > $CDIR/$CONFIG/cpp_$CONFIG.fcm


Add a F90 file that handles initial conditions to MY_SRC::

  cp $START_FILES/dtatsd.F90 $CDIR/$CONFIG/MY_SRC/

Compile NEMO::

	./makenemo -n $CONFIG -m mobius_intel -j 10


Build tools on livljobs4
========================

To generate bathymetry, initial conditions and grid information we first need
to compile some of the NEMO TOOLS (after a small bugfix - and to allow direct
passing of arguments). **For some reason GRIDGEN doesn’t like INTEL.**
**Do this on livljobs4**::

  ssh livljobs4

Copy PATHS again::

	export WDIR=/work/$USER/NEMO/Solo
	export INPUTS=/work/$USER/NEMO/INPUTS
	export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

Apply patches::

  cd $TDIR/WEIGHTS/src
  patch -b < $START_FILES/scripinterp_mod.patch
  patch -b < $START_FILES/scripinterp.patch
  patch -b < $START_FILES/scrip.patch
  patch -b < $START_FILES/scripshape.patch
  patch -b < $START_FILES/scripgrid.patch

Setup for PGI modules and compile::

  cd $TDIR
  cp $START_FILES/arch-pgf90_linux_jb.fcm $CDIR/../ARCH/arch-pgf90_linux_jb.fcm

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
It is from global NEMO 1/12, and in INPUTS::

  ls -lh $START_FILES/coordinates_ORCA_R12.nc

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

  $livljobs2$ scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.

  ferret
  use coordinates_ORCA_R12.nc
  set win 1; shade/X=50:730/Y=1250:1800 E2T, nav_lon, nav_lat ; go fland
  set win 2; set viewport upper; shade/i=50:730/j=1250:1800 NAV_LAT
  set win 2; set viewport lower; shade/i=50:730/j=1250:1800 NAV_LON

Use indices  **i=50:730 j=1250:1800**



---


Copy namelist file from INPUTS and edit with new indices, retaining use of
ORCA_R12 as course parent grid. (Still on livljobs4)::

  cd $TDIR/GRIDGEN
  cp $START_FILES/namelist_R12 ./
  vi namelist_R12
  ...
  cn_parent_coordinate_file = '../../../../INPUTS/coordinates_ORCA_R12.nc'
  ...
  nn_imin = 50
  nn_imax = 730
  nn_jmin = 1250
  nn_jmax = 1800
  nn_rhox  = 7
  nn_rhoy = 7

  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe

This generates ``1_coordinates_ORCA_R12.nc``,

Collect built items specific to the new configuration in INPUTS.
Move this coords file there as ``coordinates.nc::

  mv 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

File summary::

  ncdump -h $WDIR/INPUTS/coordinates.nc
  netcdf coordinates {
  dimensions:
   x = 4768 ;
   y = 3858 ;
   z = 1 ;
   time = UNLIMITED ; // (1 currently)
  variables:
   float nav_lon(y, x) ;
     nav_lon:units = "degrees_east" ;
     nav_lon:valid_min = 76.9642f ;
     nav_lon:valid_max = 133.7143f ;
     nav_lon:long_name = "Longitude" ;
   float nav_lat(y, x) ;
     nav_lat:units = "degrees_north" ;
     nav_lat:valid_min = -20.03138f ;
     nav_lat:valid_max = 24.65656f ;
     nav_lat:long_name = "Latitude" ;
   float nav_lev(z) ;
   float time(time) ;
   int time_steps(time) ;
   double glamt(z, y, x) ;
   double glamu(z, y, x) ;
   double glamv(z, y, x) ;
   double glamf(z, y, x) ;
   double gphit(z, y, x) ;
   double gphiu(z, y, x) ;
   double gphiv(z, y, x) ;
   double gphif(z, y, x) ;
   double e1t(z, y, x) ;
   double e1u(z, y, x) ;
   double e1v(z, y, x) ;
   double e1f(z, y, x) ;
   double e2t(z, y, x) ;
  }


Now we need to generate a bathymetry on this new grid.

----

2. Generate bathymetry file
+++++++++++++++++++++++++++

Download some GEBCO 2014 data (75E,-21N,134E,25N) and copy to $INPUTS::

..  scp ~/Downloads/RN-5922_1488296787410/GEBCO_2014_2D_-4.7361_53.0299_-2.5941_54.4256.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/.

Copy namelist for reshaping GEBCO data::

  cp $START_FILES/namelist_reshape_bilin_gebco $WDIR/INPUTS/.

Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h GEBCO_2014_2D_-4.7361_53.0299_-2.5941_54.4256.nc`` to get input
variable names)::

  vi $WDIR/INPUTS/namelist_reshape_bilin_gebco
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


Do some things to 1) flatten out land elevations, 2) make depths positive. *(James
noted a problem with the default nco module)*::

  cd $WDIR/INPUTS
  module load nco/4.5.0
  ncap2 -s 'where(elevation > 0) elevation=0' GEBCO_2014_2D_-4.7361_53.0299_-2.5941_54.4256.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc


Restore the original parallel modules, which were removed to fix tool building issue::

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
