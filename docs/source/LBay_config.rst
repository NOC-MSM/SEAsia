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

  export CONFIG=LBay
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
      .. Tom:: Create Start File Folder from Jeff's INPUTS tar amd Workspace and Ash's ARCH files

        cd /work/thopri/NEMO/
        tar xvfz INPUTS.tar.gz
        cp INPUTS/*.patch $START_FILES
        cp INPUTS/coordinates_ORCA_R12.nc $START_FILES
        cp Mobius/arch* $START_FILES
        cp Mobius/1arch* $START_FILES
        cp /work/jelt/NEMO/SEAsia/START_FILES/arch-pgf90_linux_jb.fcm $START_FILES
        cp INPUTS/dtatsd.F90 $START_FILES
        cp INPUTS/namelist_R12 $START_FILES
        cp INPUTS/namelist_reshape_bilin_gebco $START_FILES
        rm INPUTS

      .. Jeff::
        ln -s /work/thopri/NEMO/INPUTS $START_FILES

        cp /work/thopri/NEMO/J_INPUTS/*patch $START_FILES/. #I have removed this directory to reduce duplication so will need changing in future (currently /NEMO/SEAsia/START_FILES) tar file is in my NEMO directory

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
  %XIOS_HOME           /work/thopri/NEMO/LBay/xios-1.0
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

Copy PATHS again:: #I added some paths here and changed some to match the ones used in MOBIUS.

  export CONFIG=LBay
  export WDIR=/work/$USER/NEMO/$CONFIG
  export START_FILES=$WDIR/START_FILES # generic stuff for making more stuff. Mostly code.
  export INPUTS=$WDIR/INPUTS         # config specific stuff that gets made and is for running NEMO
  export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
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
  #get arch file from Jeff's workspace first
  cp $START_FILES/arch-pgf90_linux_jb.fcm $TDIR/../ARCH/arch-pgf90_linux_jb.fcm

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

Use indices  **i=3385:3392 j=2251:2266**

---

Copy namelist file from INPUTS and edit with new indices, retaining use of
ORCA_R12 as course parent grid. Keep same grid ie. 1/12 degree, so scale factors are unitary. (Still on livljobs4)::

  cd $TDIR/GRIDGEN
  cp $START_FILES/namelist_R12 ./
  vi namelist_R12
  ...

  JEFF::
  cn_parent_coordinate_file = '../../../../INPUTS/coordinates_ORCA_R12.nc'

  TOM::
  cn_parent_coordinate_file = '../../../../START_FILES/coordinates_ORCA_R12.nc'
  cn_parent_coordinate_file = '../../../../START_FILES/coordinates_ORCA_R12.nc'

  ...
  nn_imin = 3385
  nn_imax = 3392
  nn_jmin = 2251
  nn_jmax = 2266
  nn_rhox  = 1
  nn_rhoy = 1

  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe

This generates ``1_coordinates_ORCA_R12.nc``,

Collect built items specific to the new configuration in INPUTS.
Move this coords file there as ``coordinates.nc::

TOM::
  cd $WDIR
  mkdir INPUTS 
  mv $WDIR/GRIDGEN/1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

File summary::

  ncdump -h $WDIR/INPUTS/coordinates.nc
netcdf coordinates {
dimensions:
        x = 10 ;
        y = 18 ;
        z = 1 ;
        time = UNLIMITED ; // (1 currently)
variables:
        float nav_lon(y, x) ;
                nav_lon:units = "degrees_east" ;
                nav_lon:valid_min = -3.666563f ;
                nav_lon:valid_max = -2.68825f ;
                nav_lon:long_name = "Longitude" ;
        float nav_lat(y, x) ;
                nav_lat:units = "degrees_north" ;
                nav_lat:valid_min = 53.04568f ;
                nav_lat:valid_max = 53.98522f ;
                nav_lat:long_name = "Latitude" ;
        float nav_lev(z) ;
                nav_lev:units = "model_levels" ;
                nav_lev:valid_min = 0.f ;
                nav_lev:valid_max = 0.f ;
                nav_lev:long_name = "Model levels" ;
        float time(time) ;
                time:units = "seconds since 0000-01-01 00:00:00" ;
                time:long_name = "Time axis" ;
                time:calendar = "gregorian" ;
                time:title = "Time" ;
                time:time_origin = " 0000-JAN-01 00:00:00" ;
        int time_steps(time) ;
                time_steps:units = "timesteps since 0000-01-01 00:00:00" ;
                time_steps:long_name = "Time step axis" ;
                time_steps:title = "Time steps" ;
                time_steps:time_origin = " 0000-JAN-01 00:00:00" ;
                time_steps:tstep_sec = 0.f ;
        double glamt(z, y, x) ;
                glamt:missing_value = 1.e+20f ;
        double glamu(z, y, x) ;
                glamu:missing_value = 1.e+20f ;
        double glamv(z, y, x) ;
                glamv:missing_value = 1.e+20f ;
        double glamf(z, y, x) ;
                glamf:missing_value = 1.e+20f ;
        double gphit(z, y, x) ;
                gphit:missing_value = 1.e+20f ;
        double gphiu(z, y, x) ;
                gphiu:missing_value = 1.e+20f ;
        double gphiv(z, y, x) ;
                gphiv:missing_value = 1.e+20f ;
        double gphif(z, y, x) ;
                gphif:missing_value = 1.e+20f ;
        double e1t(z, y, x) ;
                e1t:missing_value = 1.e+20f ;
        double e1u(z, y, x) ;
                e1u:missing_value = 1.e+20f ;
        double e1v(z, y, x) ;
                e1v:missing_value = 1.e+20f ;
        double e1f(z, y, x) ;
                e1f:missing_value = 1.e+20f ;
        double e2t(z, y, x) ;
                e2t:missing_value = 1.e+20f ;
        double e2u(z, y, x) ;
                e2u:missing_value = 1.e+20f ;
        double e2v(z, y, x) ;
                e2v:missing_value = 1.e+20f ;
        double e2f(z, y, x) ;
                e2f:missing_value = 1.e+20f ;



Now we need to generate a bathymetry on this new grid.

----

2. Generate bathymetry file
+++++++++++++++++++++++++++

Download some GEBCO 2014 data (-4E,53N,-2E,54N) and copy to $INPUTS::

Tom::

  File from BODC: GEBCO_2014_2D_-4.0_53.0_-3.0_54.0.nc

Jeff::

  livmaf$ scp ~/Downloads/RN-9621_1506544326915/GEBCO_2014_2D_75.0_-21.0_134.0_25.0.nc jelt@livljobs4.nerc-liv.ac.uk:$INPUTS/GEBCO_2014_2D5.0_-21.0_134.0_25.0.nc

Copy namelist for reshaping GEBCO data::

  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h GEBCO_2014_2D5.0_-21.0_134.0_25.0.nc`` to get input
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


Do some things to 1) flatten out land elevations, 2) make depths positive. *(James
noted a problem with the default nco module)*::

  cd $INPUTS
  module load nco/gcc/4.4.2.ncwa
  ncap2 -s 'where(elevation > 0) elevation=0' GEBCO_2014_2D_-4.0_53.0_-3.0_54.0.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

Restore the original modules for building tools, which were tampered with to fix a bathy building issue::

  module purge
  module add netcdf/gcc/4.1.3
  module add pgi/15.4

Execute first scrip thing::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Output files::

  remap_nemo_grid_gebco.nc
  remap_data_grid_gebco.nc

*(28 Sept 2017)*

Execute second scrip thing::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Output files::

  data_nemo_bilin_gebco.nc

.. note::
  This is fast on ARCHER but does not work for me on livljobs4. It just hangs. Running on ARCHER the output is::

    Using latitude bins to restrict search.
     Computing remappings between:
    Remapped regular grid for SCRIP

                          and
    Remapped NEMO grid for SCRIP

    jelt@archer$>


I need to get things moving so I will use ARCHER to build these files. (Not very helpful for some I know)
Perhaps it is a problem with the config not the tools...

.. note::

  Make a temporary directory on ARCHER and copy in ``gebco_in.nc``, ``cordinates.nc`` and ``namelist_reshape_bilin_gebco``

  ssh archer
  mkdir /work/n01/n01/jelt/tmp/
  scp $INPUTS/gebco

Execute third scip thing::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

Output files::

  bathy_meter.nc

3. Generate initial conditions
++++++++++++++++++++++++++++++








**AT END OF PROCESS NEED TO BUILD A start_files.tar BALL**
