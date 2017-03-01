=============================================
Setting up a Liverpool Bay NEMO configuration
=============================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/LBay.html

Issues that arose
=================

*


Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Follow PyNEMO recipe for Lighthouse Reef: ``http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html``

----

Recipe Notes
============

Define working directory and other useful shortcuts::

  export WDIR=/work/n01/n01/jelt/LBay/
  export INPUTS=/work/n01/n01/jelt/lighthousereef/INPUTS
  export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

Follow recipe. Step 1 inlcuded getting INPUT files. For LHReef these were all
prepared. Now they are not so make them as and when they are required::

  cd $WDIR
  mkdir INPUTS

Old code::

  cd INPUTS
  wget ftp.nerc-liv.ac.uk:/pub/general/jdha/inputs.tar.gz
  tar xvfz inputs.tar.gz
  rm inputs.tar.gz

*Now* code::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP@5709
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@629

Need to get arch files from INPUTS::

  cd $WDIR/xios-1.0
  cp $INPUTS/arch-XC30_ARCHER.* ./arch

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par


Step 2. Obtain and apply patches::

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

On first make only choose OPA_SRC::

  ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10

It breaks. Remove key_lim2 from cpp*fcm file and remake::

  vi LBay/cpp_LBay.fcm
  ...

  ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10

Copy some input files to new configuration path::

  cp $INPUTS/cpp_LH_REEF.fcm ./LBay
  cp $INPUTS/dtatsd.F90 LBay/MY_SRC/

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
  ./maketools -n WEIGHTS -m XC_ARCHER_INTEL
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module swap PrgEnv-intel PrgEnv-cray
  module load cray-netcdf cray-hdf5
  ./maketools -n GRIDGEN -m XC_ARCHER

  module swap PrgEnv-cray PrgEnv-intel

*(1 March 2017)*

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

Note, I used FERRET locally::

  $livljobs2$ scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.
  ferret etc
  shade/i=3369:3392/j=2251:2266 NAV_LAT
  shade/i=3369:3392/j=2251:2266 NAV_LON


Copy namelist file from LH_reef and edit with new indices, retaining use of
ORCA_R12 as course
parent grid::

  cd $TDIR/GRIDGEN
  cp $INPUTS/namelist_R12 ./
  vi namelist_R12
  ...
  cn_parent_coordinate_file = '../../../../INPUTS/coordinates_ORCA_R12.nc'
  ...
  nn_imin = 3369
  nn_imax = 3392
  nn_jmin = 2251
  nn_jmax = 2266
  nn_rhox  = 7
  nn_rhoy = 7

  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe
  cp 1_coordinates_ORCA_R12.nc $WDIR/INPUTS/coordinates.nc

This creates a coordinates.nc file with contents, which are now copied to
INPUTS::

  dimensions:
  	x = 176 ;
  	y = 113 ;
  	z = 1 ;
  	time = UNLIMITED ; // (1 currently)
  variables:
    float nav_lon(y, x) ;
    float nav_lat(y, x) ;
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
    double e2u(z, y, x) ;
    double e2v(z, y, x) ;
    double e2f(z, y, x) ;

Now we need to generate a bathymetry on this new grid.



2. Generate bathymetry file
+++++++++++++++++++++++++++

Download some GEBCO data and copy to ARCHER::

  scp ~/Downloads/RN-5922_1488296787410/GEBCO_2014_2D_-4.7361_53.0299_-2.5941_54.4256.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/.

Copy namelist for reshaping GEBCO data::

  cp $INPUTS/namelist_reshape_bilin_gebco $WDIR/INPUTS/.

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


Peeking at this output using FERRET suggests that something went a bit wrong at
the outer western limits of the region. I think that perhaps the domain extends
beyond where I got GEBCO data so there is some strange extrapolation between Angelsey. However, it works!

**ACTION:** Need to fix this. Make the coordinates.nc have reduced westward extent.


3. Generate initial conditions
++++++++++++++++++++++++++++++

4. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

5. Generate mesh and mask files for open boundary conditions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

6. Generate boundary conditions with PyNEMO
+++++++++++++++++++++++++++++++++++++++++++

----



Need to make a new cutdown GEBCO file. Should try and make it match James' variables names::

  jelt@eslogin007:/work/n01/n01/jelt/LBay/INPUTS> ncdump -h gebco_1_cutdown.nc
  netcdf gebco_1_cutdown {
  dimensions:
  	latitude = 361 ;
  	longitude = 601 ;
  variables:
  	float latitude(latitude) ;
  		latitude:units = "degrees N" ;
  		latitude:standard_name = "Latitude" ;
  		latitude:long_name = "Latitude degrees N" ;
  	float longitude(longitude) ;
  		longitude:units = "degrees E" ;
  		longitude:standard_name = "Longitude" ;
  		longitude:long_name = "Longitude degrees E" ;
  	float topo(latitude, longitude) ;
  		topo:units = "metres" ;
  		topo:standard_name = "topography" ;
  		topo:long_name = "topography" ;
  		topo:_FillValue = -1.e+34f ;

  // global attributes:
  		:description = "Gebco 1 minute topography from Gebco Atlas CD." ;
  		:author = "James Harle" ;
  		:date = "24/05/2005" ;
  		:history = "Thu Oct 22 08:13:23 2015: ncks -d longitude,-90.,-80. -d latitude,15.,21. gebco_1.nc gebco_1_cd_v2.nc" ;
  		:NCO = "4.4.2" ;

E.g.::

  ncea -d lat,50,54 -d lon,350,360 /projectsa/FASTNEt/kariho40/AMM60/BATHY/bathy_AMM60.nc  /scratch/jelt/tmp/GEBCO_cutdown.nc

Hmm got this far. Downloaded some GEBCO data from BODC. Have not made a cutdown
bathymetry file yet.


When it is ready proceed with the following

To do:
++++++


To do: "To create the bathymetry we use the gebco dataset. On ARCHER I had to use a
non-default nco module for netcdf operations to work. I also had to cut down
the gebco data as the SCRIP routines failed for some unknown reason"::

  cp $INPUTS/gebco_1_cutdown.nc $WDIR/INPUTS/.
  cp $INPUTS/namelist_reshape_bilin_gebco $WDIR/INPUTS/.
  cd $WDIR/INPUTS
  module load nco/4.5.0
  ncap2 -s 'where(topo > 0) topo=0' gebco_1_cutdown.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco


Hang on. Isn't PyNEMO supposed to do all this hard work in defining the domain?

----

Try again from the start
++++++++++++++++++++++++

Following ``http://pynemo.readthedocs.io/en/latest/examples.html``::

  export WDIR=/work/n01/n01/jelt/LBay/

  cd pynemo/trunk/Python
  cp data/namelist.bdy $WDIR
  cd $WDIR
  vi namelist.bdy
  ...
  sn_src_dir = ‘/work/n01/n01/jelt/LBay/test.ncml’
  sn_dst_dir = ‘/work/n01/n01/jelt/LBay/OUTPUT’
  cn_mask_file = ‘/work/n01/n01/jelt/LBay/mask.nc’)



  ssh -Y espp1
  module load anaconda
  source activate pynemo_env
  export WDIR=/work/n01/n01/jelt/LBay/

  cd $WDIR
  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/pynemo/lib/python2.7/site-packages/:$PYTHONPATH

  ~/.conda/envs/pynemo/bin/pynemo_ncml_generator

On generation got an error saying **Not all the variables have been defined**
Perhaps I could get a look at a proper ncml file so I can see what all these variables are?
Had a look in lighthousereef/INPUTS. There are two different files there.
