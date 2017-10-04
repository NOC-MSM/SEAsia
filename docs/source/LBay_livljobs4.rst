============================================================
Setting up a Liverpool Bay NEMO configuration, inc. livljob4
============================================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/LBay_livljobs4.html

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Originally tried building everything on ARCHER but ran into a java problem with PyNEMO that they needed to fix. Also MDP can't use ARCHER.

The summary procedure:
1) start with running NEMO on a new domain to generate grid files
2) generate boundary forcings (e.g. tides, T,S,U,V) from a parent grid / model
3) run the new configuration

Step 2 was done on livljobs4. The other steps were on ARCHER.

Issues that arose
=================

* PyNEMO doesn't yet deal with sigma parent grids.

Other related recipes:
Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Follow PyNEMO recipe for Lighthouse Reef: ``http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html``
Follow PyNEMO recipe for LBay on ARCHER (not complete because PyNEMO java issue): ``http://nemo-reloc.readthedocs.io/en/latest/LBay.html``.
(Perhaps this other note is superceded by this present one, though the other
note includes information of using PyNEMO with different data sources).

.. note::

  It is very easy to break the code with bad edits to the iodef.xml file. Don't
  change the iodef.xml file at the same time as something else.

----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO bits are done on livljobs4.
(There was a problem with some java bits working on ARCHER)
Starting on ARCHER::

  ssh login.archer.ac.uk

  export WDIR=/work/n01/n01/jelt/LBay/
  export INPUTS=$WDIR/INPUTS
  export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS


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

Copy some input files to new configuration path::

  cp $INPUTS/cpp_LH_REEF.fcm ./LBay/cpp_LBay.fcm
  cp $INPUTS/dtatsd.F90 LBay/MY_SRC/

On first make only choose OPA_SRC::

  ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10

It breaks. Remove key_lim2 from cpp*fcm file and remake::

  vi LBay/cpp_LBay.fcm
  ...

  ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10



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
  shade/i=3385:3392/j=2251:2266 NAV_LAT
  shade/i=3385:3392/j=2251:2266 NAV_LON


Copy namelist file from LH_reef and edit with new indices, retaining use of
ORCA_R12 as course parent grid::

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
  cp 1_coordinates_ORCA_R12.nc $WDIR/INPUTS/coordinates.nc

This creates a coordinates.nc file with contents, which are now copied to
INPUTS::

  dimensions:
  	x = 57 ;
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



3. Generate initial conditions
++++++++++++++++++++++++++++++


Copy ``make.macro`` file and edit the path if necessary::
**FIX** to the notes (copied from jdha instead): ``cp $WDIR/INPUTS/make.macro ./``::

  cp /home/n01/n01/jdha/sosie/make.macro /home/n01/n01/jelt/sosie/.

  vi /home/n01/n01/jelt/sosie/make.macro
  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/jelt/local

Proceed with Step 6::

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

  cp $INPUTS/initcd_votemper.namelist .
  cp $INPUTS/initcd_vosaline.namelist .

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
  cd $WDIR/INPUTS

  ncks -d x,560,620 -d y,720,800 /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT/AMM60_5d_20131013_20131129_grid_T.nc $WDIR/INPUTS/cut_down_20131013_LBay_grid_T.nc

Average over time and restore the parallel modules::

  ncwa -a time_counter $WDIR/INPUTS/cut_down_20131013_LBay_grid_T.nc  $WDIR/INPUTS/cut_down_201310_LBay_grid_T.nc

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

  cp $INPUTS/namelist_reshape_bilin_initcd_votemper $WDIR/INPUTS/.
  cp $INPUTS/namelist_reshape_bilin_initcd_vosaline $WDIR/INPUTS/.

Edit the input files::

  vi $WDIR/INPUTS/namelist_reshape_bilin_initcd_votemper
  &grid_inputs
    input_file = 'thetao_AMM60-LBay_2013.nc4'
  ...

  &interp_inputs
    input_file = "thetao_AMM60-LBay_2013.nc4"
  ...

Simiarly for the *vosaline.nc file::

  vi $WDIR/INPUTS/namelist_reshape_bilin_initcd_vosaline
  &grid_inputs
    input_file = 'so_AMM60-LBay_2013.nc4'
  ...

  &interp_inputs
    input_file = "so_AMM60-LBay_2013.nc4"
  ...


Produce the remap files::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates ``data_nemo_bilin_R12.nc``. Then::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates ``initcd_votemper.nc``. Then::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates ``initcd_vosaline.nc``.


4. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

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

  cp $INPUTS/namelist_reshape_bilin_atmos $WDIR/INPUTS/.
  cp $INPUTS/namelist_reshape_bicubic_atmos $WDIR/INPUTS/.

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

  cd $WDIR/INPUTS
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos

Generate  remap files ``remap_nemo_grid_atmos.nc`` and ``remap_data_grid_atmos.nc``. Then::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos

Generates ``data_nemo_bilin_atmos.nc``. Then::

  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos

Generates ``weights_bilinear_atmos.nc``. Then::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos

Generates ``data_nemo_bicubic_atmos.nc``. Then::

  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

Generates ``weights_bicubic_atmos.nc``.


5. Generate mesh and mask files for open boundary conditions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Run the model to generate the mesh and mask files::

  cd $CDIR
  cp $INPUTS/cpp_LH_REEF.fcm LBay/cpp_LBay.fcm
  ln -s $WDIR/INPUTS/bathy_meter.nc $CDIR/LBay/EXP00/bathy_meter.nc
  ln -s $WDIR/INPUTS/coordinates.nc $CDIR/LBay/EXP00/coordinates.nc
  cp $INPUTS/runscript $CDIR/LBay/EXP00
  cp $INPUTS/namelist_cfg $CDIR/LBay/EXP00/namelist_cfg
  cp $INPUTS/namelist_ref $CDIR/LBay/EXP00/namelist_ref
  ./makenemo clean
  ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10
  cd LBay/EXP00
  ln -s $WDIR/xios-1.0/bin/xios_server.exe xios_server.exe

Edit the namelist files for this configuration::

  ncdump -h coordinates.nc
  x = 57 ;
  y = 113 ;

  vi namelist.cfg
  ...
  cn_exp      =   "LBay"  !  experience name
  ...
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     cp_cfg      =  "lbay"                !  name of the configuration
     jp_cfg      =     084               !  resolution of the configuration
     jpidta      =      57               !  1st lateral dimension ( >= jpi )
     jpjdta      =     113               !  2nd    "         "    ( >= jpj )
     jpkdta      =      51               !  number of levels      ( >= jpk )
     jpiglo      =      57               !  1st dimension of global domain --> i =jpidta
     jpjglo      =     113               !  2nd    -                  -    --> j  =jpjdta

**ACTION: There are further edits to be made for when the model is actually run**
**E.g. other filename instances of Lbay**

Note that the old LH_REEF has the following
| jpidta      =     358               !  1st lateral dimension ( >= jpi )
| jpjdta      =     428               !  2nd    "         "    ( >= jpj )

with the dimensions in the LH_REFF coordinates file as
| ncdump -h coordinates.nc
| x = 358 ;
| y = 428 ;

Edit the runscript to include modules and the Account name (n01-NOCL)::

  vi runscript

  #!/bin/bash
  #PBS -N LBay
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel
  ...

Submit::

  qsub -q short runscript


*(6 March 2017)*

If that works, we then need to rebuild the mesh and mask files in to single files for the next step::

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_zgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_hgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mask 96
  mv mesh_zgr.nc mesh_hgr.nc mask.nc $WDIR/INPUTS
  rm mesh_* mask_* LBay_0000*
  cd $WDIR/INPUTS


THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO *(20 Sept 2017)*



6. Generate boundary conditions with PyNEMO: Create netcdf abstraction wrapper
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this section there are two stages.
* generate a ncml file which describes the files needed to create boundary conditions
* generate a namelist.bdy file which controls the actual boundary condition generation.

For each parent data set a new pair of (``*.ncml``, ``namelist.bdy``) are needed.
Here I attempt to use parent data from:
* AMM60 local data (doesn't yet work because of the sigma levels)
* thredds server (as in the LH_REEF example)
* NNA local data (easiest ?)

First install PyNEMO if not already done so. Full description::

  ssh -Y livljobs4
  cd /work/jelt
  mkdir /work/jelt/NEMO/LBay
  export WDIR=/work/jelt/NEMO/LBay/
  module load anaconda/2.1.0  # Want python2
  conda create --name nrct_env scipy=0.16.0 numpy matplotlib=1.5.1 basemap netcdf4 libgfortran=1.0.0
  source activate nrct_env
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 # Note had to add https path
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

Find java object by doing a which java and then following the trail
find  /usr/lib/jvm/jre-1.7.0-openjdk.x86_64/ -name libjvm.so -print
::

  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS # Didn't need this on ARCHER...
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd nrct/Python
  python setup.py build
  export PYTHONPATH=/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $WDIR/INPUTS


I suggest managing the namelist.bdy file after the ``ncml`` file is generated.
 A fresh ``ncml`` file can be generated automatically or an existing one can be edited.


6a. Generate ncml files (Not tested on livljobs4)
+++++++++++++++++++++++


Activate generator:

Start up pynemo and generate boundary conditions. First we need to create a
few ncml files to gather input data and map variable names. Then using pynemo
we define the area we want to model.
Redefine ``WDIR``. Launch from WDIR::

  ssh -Y espp1
  module load anaconda
  source activate pynemo_env
  #  export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  #  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/:$PYTHONPATH
  cd $WDIR/INPUTS
  pynemo_ncml_generator


Here the object is to generate a ncml file that is read in by PyNEMO as the ``sn_src_dir``
(in the ``namelist.bdy`` file)

Fill in the Tracer and Dynamics for T,S,U,V,Z tabs: using T,T & U,V,T in the reg
expressions e.g. .*T\.nc$
To generate a e.g. ``inputs_src.ncml`` file click  **generate**. Defining the
filename seems to work better with the file selector rather than direct typing.

In the following I have three ncml files.
* One for using the thredds server to get remote ORCA12 data.
* One for using local AMM60 data, with ackward s-sigma levels
* One for using local NNA data

NNA_inputs_src.ncml
+++++++++++++++++++

Note need to set the time variables and new ``sn_src_dir`` in namelist.bdy.
Actually upated the following with all the Jan 2000 files::

  cd $WDIR/INPUTS
  vi NNA_inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="votemper" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/jdha/LBay/INPUTS/NNA" regExp=".*T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="vosaline" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/jdha/LBay/INPUTS/NNA" regExp=".*T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="vozocrtx" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/jdha/LBay/INPUTS/NNA" regExp=".*U\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="vomecrty" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/jdha/LBay/INPUTS/NNA" regExp=".*V\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="sossheig" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/jdha/LBay/INPUTS/NNA" regExp=".*T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>




6b. Generate the namelist.bdy file for PyNEMO
+++++++++++++++++++++++++++++++++++++++++++++


Copy the PyNEMO template namelist.bdy from the lighthouse project::

  cd $WDIR/INPUTS
  cp $INPUTS/namelist.bdy $WDIR/INPUTS/.

Edit namelist.bdy to for the configuration name and ``ncml`` file name. **Note
need the slash following OUTPUT**::

  vi namelist.bdy
  sn_src_dir = './inputs_src.ncml'       ! src_files/'
  sn_dst_dir = '/work/n01/n01/jelt/LBay/OUTPUT/'
  sn_fn      = 'LBay'                 ! prefix for output files
  ...
  cn_mask_file   = './mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)

Now edit the pynemo namelist file. Add location of grid information. Note had to
 hunt for a mesh.nc file. Incase this doesn't work, there were a couple of
 options. (Tried both) Note also that mesh_zgr includes gdept_0, gdepw_0, e3t_0, e3u_0,
 e3v_0, e3w_0, so use ncml to convert to variables without *_0. (Also didn't convert e3w_0).

 Make sure the timestamps correspond to the input data.
Turn off as many things as possible to help it along.
Turned off ``ln_mask_file``. James said it was for outputting a new mask file
but it might have given me trouble.

I have a namelist.bdy file for each ncml configuration
* namelist.bdy_AMM60
* namelist.bdy_thredds (uses global 1/12 degree data)
* namelist.bdy_NNA




7. Generate boundary conditions with PyNEMO: Run PyNEMO
+++++++++++++++++++++++++++++++++++++++++++++++++++++++



Using livljobs4

*(20/21 Sept 2017)*

**Start the process again on livljobs4: LBay_livljobs4.rst**

Need to grab some INPUT files from ARCHER because I am not building them with livljobs4::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/namelist.bdy $WDIR/INPUTS/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/bathy_meter.nc $WDIR/INPUTS/.
  rsync -uartv jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/ $WDIR/INPUTS
  cd $WDIR/INPUTS

Make sure the NNA data is available::

  mkdir $WDIR/INPUTS/NNA
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mask.nc $WDIR/INPUTS/NNA/.
  for file in NNA_*200001*nc ; do scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/$file $WDIR/INPUTS/NNA/. ; done

Make sure the destination meshes exist. These were generated by running the new config for a timestep::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/mesh_zgr.nc $WDIR/INPUTS/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/mesh_hgr.nc $WDIR/INPUTS/.



namelist.bdy_NNA
++++++++++++++++

Edit namelist.bdy to reflect locally stored mesh and mask files. Also
NNA_inputs_src.ncml. Set the date info back to (Nov?) 1979.

 ::

   vi namelist.bdy_NNA

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
      sn_src_hgr = '/work/jelt/NEMO/LBay/INPUTS/NNA/mesh_hgr.nc'   !  /grid/
      sn_src_zgr = '/work/jelt/NEMO/LBay/INPUTS/NNA/mesh_zgr.nc'
      sn_dst_hgr = './mesh_hgr.nc'
      sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
      sn_src_msk = '/work/jelt/NEMO/LBay/INPUTS/NNA/mask.nc'
      sn_bathy   = './bathy_meter.nc'

   !-----------------------------------------------------------------------
   !  I/O
   !-----------------------------------------------------------------------
      sn_src_dir = './NNA_inputs_src.ncml'       ! src_files/'
      sn_dst_dir = '/work/jelt/NEMO/LBay/INPUTS/'
      sn_fn      = 'LBay'                 ! prefix for output files
      nn_fv      = -1e20                     !  set fill value for output files
      nn_src_time_adj = 0                                    ! src time adjustment
      sn_dst_metainfo = 'metadata info: jelt'

   !-----------------------------------------------------------------------
   !  unstructured open boundaries
   !-----------------------------------------------------------------------
       ln_coords_file = .true.               !  =T : produce bdy coordinates files
       cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
       ln_mask_file   = .false.              !  =T : read mask from file
       cn_mask_file   = './mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
       ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
       ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
       ln_tra         = .true.               !  boundary conditions for T and S
       ln_ice         = .false.               !  ice boundary condition
       nn_rimwidth    = 9                    !  width of the relaxation zone

   !-----------------------------------------------------------------------
   !  unstructured open boundaries tidal parameters
   !-----------------------------------------------------------------------
       ln_tide        = .true.               !  =T : produce bdy tidal conditions
       clname(1)      = 'M2'                 ! constituent name
       clname(2)      = 'S2'
       clname(3)      = 'K2'
       ln_trans       = .false.
       sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
       sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'


   !-----------------------------------------------------------------------
   !  Time information
   !-----------------------------------------------------------------------
       nn_year_000     = 2000        !  year start
       nn_year_end     = 2000        !  year end
       nn_month_000    = 01          !  month start (default = 1 is years>1)
       nn_month_end    = 01          !  month end (default = 12 is years>1)
       sn_dst_calendar = 'gregorian' !  output calendar format
       nn_base_year    = 1979        !  base year for time counter
       sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'

   !-----------------------------------------------------------------------
   !  Additional parameters
   !-----------------------------------------------------------------------
       nn_wei  = 1                   !  smoothing filter weights
       rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                     !  smoothing onto dst points. Need to
                                     !  make this a funct. of dlon
       sn_history  = 'bdy files produced by jelt from AMM60 for testing'
                                     !  history for netcdf file
       ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
       nn_alpha    = 0               !  Euler rotation angle
       nn_beta     = 0               !  Euler rotation angle
       nn_gamma    = 0               !  Euler rotation angle
       rn_mask_max_depth = 300.0     !  Maximum depth to be ignored for the mask
       rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break



Generate the boundary conditions again, with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  export WDIR=/work/jelt/NEMO/LBay/
  cd $WDIR/INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -g -s namelist.bdy_NNA

Once the area of interest is selected and the close button is clicked, open
boundary data should be generated in $WDIR/OUTPUT

The SAVE button only updates the ``namelist.bdy`` file. The CLOSE button activates the process.

This generates::
  ls -1 /work/n01/n01/jelt/LBay/OUTPUT

  coordinates.bdy.nc
  LBay_bdytide_rotT_M2_grid_T.nc
  LBay_bdytide_rotT_K2_grid_T.nc
  LBay_bdytide_rotT_S2_grid_T.nc
  LBay_bdytide_rotT_M2_grid_U.nc
  LBay_bdytide_rotT_K2_grid_U.nc
  LBay_bdytide_rotT_S2_grid_U.nc
  LBay_bdytide_rotT_M2_grid_V.nc
  LBay_bdytide_rotT_K2_grid_V.nc
  LBay_bdytide_rotT_S2_grid_V.nc
  LBay_bdyT_y2000m01.nc
  LBay_bt_bdyT_y2000m01.nc
  LBay_bdyU_y2000m01.nc
  LBay_bdyV_y2000m01.nc

.. Warning::

  Though it doesn't quite work with ``ln_tra = .false.``

  This wont work because variable ``ft`` which deals with the number of time steps
  if only defined using the T fields, but needed for the velocity bcs. Wont work
  with ln_dyn3d = .true. either
  See e.g.::

    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 435, in process_bdy
     ft, num_bdy, time_counter, unit_origin)
     UnboundLocalError: local variable 'ft' referenced before assignment




8. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

*(21 Sept 2017)*

When I've got all the bdy files need to fix some variable names.
Copy the new files back onto ARCHER
::

  livljobs4$
  cd /work/jelt/NEMO/LBay/INPUTS
  for file in LBay*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/. ; done


Open a terminal on ARCHER.

Compile code with new flags (not sure if I need the module commands). Note in particular,
``key_gen_IC`` and ``key_tide``::

 module swap PrgEnv-cray PrgEnv-intel
 module load cray-netcdf-hdf5parallel
 module load cray-hdf5-parallel
 export WDIR=/work/n01/n01/jelt/LBay/
 export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG

 vi $CDIR/LBay/cpp_LBay.fcm
 bld::tool::fppkeys   key_dynspg_ts key_ldfslp  key_zdfgls  key_vvl key_mpp_mpi key_netcdf4 key_nosignedzero  key_iomput key_gen_IC key_bdy key_tide


 cd $CDIR
 ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10 clean
 ./makenemo -n LBay -m XC_ARCHER_INTEL -j 10

Prepare the boundary files (need to fix some variable names)::

 export WDIR=/work/n01/n01/jelt/LBay/
 export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG

 cd $WDIR/INPUTS

 module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
 module load nco/4.5.0
 ncrename -v deptht,gdept LBay_bdyT_y2000m01.nc
 ncrename -v depthu,gdepu LBay_bdyU_y2000m01.nc
 ncrename -v depthv,gdepv LBay_bdyV_y2000m01.nc

 ncrename -v deptht,gdept initcd_votemper.nc
 module unload nco
 module load cray-netcdf-hdf5parallel cray-hdf5-parallel


Link the boundary data to the EXP direcory and update the namelist_cfg for
running, not mesh generation::

 cd $CDIR/LBay/EXP00
 ln -s $WDIR/INPUTS/coordinates.bdy.nc $CDIR/LBay/EXP00/coordinates.bdy.nc
 ln -s $WDIR/INPUTS/LBay_bdyT_y2000m01.nc $CDIR/LBay/EXP00/LBay_bdyT_y2000m01.nc
 ln -s $WDIR/INPUTS/LBay_bdyU_y2000m01.nc    $CDIR/LBay/EXP00/LBay_bdyU_y2000m01.nc
 ln -s $WDIR/INPUTS/LBay_bdyV_y2000m01.nc    $CDIR/LBay/EXP00/LBay_bdyV_y2000m01.nc
 ln -s $WDIR/INPUTS/LBay_bt_bdyT_y2000m01.nc  $CDIR/LBay/EXP00/LBay_bt_bdyT_y2000m01.nc
 ln -s $WDIR/INPUTS $CDIR/LBay/EXP00/bdydta
 sed -e 's/nn_msh      =    3/nn_msh      =    0/' namelist_cfg > tmp
 sed -e 's/nn_itend    =      1/nn_itend    =       1440 /' tmp > namelist_cfg

Edit the link to the tidal boundary files to fix file names::

  vi namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &nambdy_tide     ! tidal forcing at open boundaries
  !-----------------------------------------------------------------------
    filtide      = 'bdydta/LBay_bdytide_rotT_'         !  file name root of tidal forcing files

Change flag so that boundary data is read in. Spotted ``nn_dyn2d_dta = 3`` in AMM60 run.
Switch from 1 to 3. I.e.::

   &nambdy        !  unstructured open boundaries                          ("key_bdy")
       nb_bdy         = 1                    !  number of open boundary sets
       ln_coords_file = .true.               !  =T : read bdy coordinates from file
       cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
       ln_mask_file   = .true.              !  =T : read mask from file
       cn_mask_file   = '../../../../../INPUTS/LBay_bdyT_y2000m01.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
       cn_dyn2d       = 'flather'               !
       nn_dyn2d_dta   =  3                   !  = 0, bdy data are equal to the initial state
                                             !  = 1, bdy data are read in 'bdydata   .nc' files
                                             !  = 2, use tidal harmonic forcing data from files
                                             !  = 3, use external data AND tidal harmonic forcing

Also need to make sure the harmonic tidal boundary files are consistent with the
 harmonics expected e.g.::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters (#ifdef key_tide)
  !-----------------------------------------------------------------------
      clname(1) = 'K2'
      clname(2) = 'M2'
      clname(3) = 'S2'



.. Note:

  I had a problem with initial T,S conditions because the generated netCDF files
   only had vector fields for the z-coordinate. However, Using ``key_gen_IC``
   generates the vertical velocity on the fly.

  Completes. Works as a restart or from initial conditions::

    ln_rstart   =  .false.  !  start from rest (F) or from a restart file (T)
    ln_tsd_init   = .true.   !  Initialisation of ocean T & S with T &S input
     data (T) or not (F)

  OR as::

    ln_rstart   =  .true.  !  start from rest (F) or from a restart file (T)
    ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input
     data (T) or not (F)





Should also check the xml files. There was something **fishy** with the
``field_def.xml`` and ``iodef.xml`` files where variables were not defined in
the NEMO checkout. Copy these files from lighthouse reef::

 cp /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/EXP00/iodef.xml  $CDIR/LBay/EXP00/iodef.xml
 cp /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/SHARED/field_def.xml  $CDIR/SHARED/field_def.xml

Edit the output to have 1hrly SSH::

 vi $CDIR/LBay/EXP00/iodef.xml
 ...



Edit runscript: Increase the number of XIOS cores used::

 vi runscript
 ...
 aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa
 #aprun -b -n $XIOCORES -N 1 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

Then submit::

 cd $CDIR/LBay/EXP00
 qsub -q short runscript

 4806706.sdb



*(26 Sept 2017)*


Rebuild the files::

  export WDIR=/work/n01/n01/jelt/LBay/
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_T 5
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_U 5
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_V 5
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_W 5

Should remove individual processor files once the build is verified::

  rm LBay_1h_20000102_20000106_grid_?_*nc

Inspect locally e.g.::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LBay/EXP00/LBay_1h_20000102_20000106_grid_T.nc .

  ferret
  use LBay_1h_20000102_20000106_grid_T.nc
  plot /i=25/j=70 SOSSHEIG

---


**Whoo hoo! It works! With a semi-diurnal SSH signal**

At this point I'll move onto the SE Asia configuration.


----


*(4 Oct 2017)*

Compile and test in ORCHESTRA code base
+++++++++++++++++++++++++++++++++++++++

This may not work straight off as the new code needs different namelist, and other configuration management files
.. note:

  I made a *store* of architecure files in work, for compiling NEMO

  mkdir /work/n01/n01/jelt/ARCH
  cp /work/n01/n01/jelt/ACCORD/trunk_NEMOGCM_r8395/ARCH/arch-XC_ARCHER_INTEL.fcm /work/n01/n01/jelt/ARCH/.
  echo "arch-XC_ARCHER_INTEL.fcm is for compiling OPA on ARCHER with my XIOS2" > /work/n01/n01/jelt/ARCH/README_arch-XC_ARCHER_INTEL.fcm
  vi /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL.fcm
  ...
  %XIOS_HOME           /work/n01/n01/jelt/XIOS

Load paths and modules::

  export CONFIG=LBay
  export WDIR=/work/n01/n01/jelt/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel


Build NEMO ORCHESTRA trunk @ r8395::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

Use my XIOS file (see ``%XIOS_HOME``). Copy from *store*::

  cp /work/n01/n01/$USER/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Edit the CPP flags::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_diaharm       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

It compiles.

Link XIOS2 executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/XIOS/bin/xios_server.exe $EXP/xios_server.exe

----

Copy all the EXP files from the working experiment into this EXP directory::

  cp  /work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LBay/EXP00/* $EXP/.
  ln -s /work/n01/n01/jelt/LBay/INPUTS $EXP/bdydta

Fix the links with the xios and opa exectutables (incase they were broken by the cp)::
  ln -s  /work/n01/n01/$USER/XIOS/bin/xios_server.exe $EXP/xios_server.exe
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa


Submit::

   qsub -q short runscript

**Does it work?**


---

Rebuild the output and inspect
++++++++++++++++++++++++++++++

Rebuild the SSH files using old tools::

 export WDIR=/work/n01/n01/jelt/LBay/
 export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

 cd $EXP
 $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_T 5


Should remove individual processor files once the build is verified::

 rm LBay_1h_20000102_20000106_grid_?_*nc

Inspect locally e.g.::

 scp jelt@login.archer.ac.uk:$EXP/LBay_1h_20000102_20000106_grid_T.nc .

 ferret
 use LBay_1h_20000102_20000106_grid_T.nc
 plot /i=25/j=70 SOSSHEIG


**Does it work?**


---

.. note::

  **TO DO** another time / for Solent config

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
