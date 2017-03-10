=============================================
Setting up a Liverpool Bay NEMO configuration
=============================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/LBay.html

Issues that arose
=================

* Notes the ``cpp_*.fcm`` files are different::

  diff /work/n01/n01/jelt/Solent/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/Solent/cpp_Solent.fcm /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/cpp_LH_REEF.fcm
  1c1
  <  bld::tool::fppkeys key_trabbl key_dynspg_flt key_diaeiv key_ldfslp key_traldf_c2d key_traldf_eiv key_dynldf_c3d key_zdftke key_zdfddm key_zdftmx key_iomput key_mpp_mpi
  ---
  > bld::tool::fppkeys   key_dynspg_ts key_ldfslp  key_zdfgls  key_vvl key_mpp_mpi key_netcdf4 key_nosignedzero  key_iomput key_gen_IC key_bdy

* I didn't copy some stuff into the NEMO compilation. This may be a problem later::

  cp $WDIR/INPUTS/cpp_LH_REEF.fcm ./LH_REEF
  cp $WDIR/INPUTS/dtatsd.F90 LH_REEF/MY_SRC/


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

*(3 March )*
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

Obtain namelist files and data file::

  cp $INPUTS/namelist_reshape_bilin_atmos $WDIR/INPUTS/.
  cp $INPUTS/namelist_reshape_bicubic_atmos $WDIR/INPUTS/.

Generate cut down drowned precip file (note that the nco tools don't like the
parallel modules)::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  ncks -d lon,355.,358. -d lat,53.,54. /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING/drowned_precip_DFS5.1.1_y1979.nc $WDIR/INPUTS/cutdown_drowned_precip_DFS5.1.1_y1979.nc

Setup weights files for the atmospheric forcing::

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


LH_REEF
jpidta      =     358               !  1st lateral dimension ( >= jpi )
jpjdta      =     428               !  2nd    "         "    ( >= jpj )

ncdump -h coordinates.nc
x = 358 ;
y = 428 ;

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

6. Generate boundary conditions with PyNEMO: Create netcdf abstraction wrapper
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*(Reinstall pyNEMO, with updates, 10 March 2017)*

Install: Full description::

  cd ~
  module load anaconda
  conda create --name pynemo_env scipy=0.16.0 numpy matplotlib=1.5.1 basemap netcdf4 libgfortran=1.0.0
  source activate pynemo_env
  conda install -c conda-forge seawater=3.3.4
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius
  export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  svn checkout https://ccpforge.cse.rl.ac.uk/svn/pynemo
  cd pynemo/trunk/Python
  python setup.py build
  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/pynemo/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/pynemo_env
  #cp data/namelist.bdy $WDIR
  cd $WDIR

But since I have done the installation I will copy the template namelist.bdy
from the lighthouse project. (There are less things to edit this way)::

  cd $WDIR/INPUTS
  cp $INPUTS/namelist.bdy $WDIR/INPUTS/.


Edit namelist.bdy::

  vi namelist.bdy
  sn_src_dir = './inputs_src.ncml'       ! src_files/'
  sn_dst_dir = '/work/n01/n01/jelt/LBay/OUTPUT'
  sn_fn      = 'LBay'                 ! prefix for output files
  ...
  cn_mask_file   = './mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)


Activate generator:

Start up pynemo and generate boundary conditions. First we need to create a
few ncml files to gather input data and map variable names. Then using pynemo
we define the area we want to mode.
Redefine ``WDIR``. Launch from WDIR::

  ssh -Y espp1
  module load anaconda
  source activate pynemo_env
  export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/:$PYTHONPATH
  cd $WDIR/INPUTS
  pynemo_ncml_generator

*Old*

| Note the file path for output filename is
 ``/work/n01/n01/jelt/LBay/INPUTS/inputs_src.ncml`` for the work dir. Has to
  match the ``sn_src_dir``
| Source directory: ``http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data``
| For the Grid tab, added source directory only - no Reg expressions
|For the Ice tab, added source directory and a duff reg exp (latter prob not
required)
| Filled in the Tracer and Dynamics for T,S,U,V,Z tabs: using T,T & U,V,T in the reg
expressions e.g. .*T\.nc$
| Click to generate ``inputs_src.ncml`` file.

Manually editted input_src.ncml. Copied it to inputs_src.ncml_SAFE.

*New*

Decided to manually edit and create a *ncml file for local AMM60 data::

  cd $WDIR/INPUTS
  cp inputs_src.ncml AMM60_inputs_src.ncml

Manually edit to remove all the bits I don't want. Use a specific AMM60_1d*
file set::

  vi AMM60_inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT" regExp="AMM60_1d_20120221_20120420_grid_T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT" regExp="AMM60_1d_20120221_20120420_grid_T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT" regExp="AMM60_1d_20120221_20120420_grid_U\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT" regExp="AMM60_1d_20120221_20120420_grid_V\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
      <ns0:netcdf>
        <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:scan location="file://work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT" regExp="AMM60_1d_20120221_20120420_grid_T\.nc$" />
        </ns0:aggregation>
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>



Now edit the pynemo namelist file. Add location of grid information. Note had to
 hunt for a mesh.nc file. Incase this doesn't work, there were a couple of
 options. (Tried both) Note also that mesh_zgr includes gdept_0, gdepw_0, e3t_0, e3u_0,
 e3v_0, e3w_0, so use ncml to convert to variables without *_0. (Also didn't convert e3w_0).

 Make sure the timestamps correspond to the input data.
Turn off as many things as possible to help it along.
Turned off ``ln_mask_file``. James said it was for outputting a new mask file
but it might have given me trouble.

 ::

   vi namelist.bdy

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
      sn_src_hgr = '/work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/mesh_hgr.nc'   !  /grid/
      sn_src_zgr = '/work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/mesh_zgr.nc'
      sn_dst_hgr = './mesh_hgr.nc'
      sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
      sn_src_msk = '/work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP00/OUTPUT/00007200/mask.nc'
      sn_bathy   = './bathy_meter.nc'

   !-----------------------------------------------------------------------
   !  I/O
   !-----------------------------------------------------------------------
      sn_src_dir = './AMM60_inputs_src.ncml'       ! src_files/'
      sn_dst_dir = '/work/n01/n01/jelt/LBay/OUTPUT'
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
       ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
       ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
       ln_tra         = .true.               !  boundary conditions for T and S
       ln_ice         = .false.               !  ice boundary condition
       nn_rimwidth    = 9                    !  width of the relaxation zone

   !-----------------------------------------------------------------------
   !  unstructured open boundaries tidal parameters
   !-----------------------------------------------------------------------
       ln_tide        = .false.               !  =T : produce bdy tidal conditions
       clname(1)      = 'M2'                 ! constituent name
       clname(2)      = 'S2'
       clname(3)      = 'K2'
       ln_trans       = .false.
       sn_tide_h     = '/Users/jdha/Projects/pynemo_data/DATA/h_tpxo7.2.nc'
       sn_tide_u     = '/Users/jdha/Projects/pynemo_data/DATA/u_tpxo7.2.nc'

   !-----------------------------------------------------------------------
   !  Time information
   !-----------------------------------------------------------------------
       nn_year_000     = 2012        !  year start
       nn_year_end     = 2012        !  year end
       nn_month_000    = 3           !  month start (default = 1 is years>1)
       nn_month_end    = 3           !  month end (default = 12 is years>1)
       sn_dst_calendar = 'gregorian' !  output calendar format
       nn_base_year    = 1950        !  base year for time counter
       sn_tide_grid    = '/Users/jdha/Projects/pynemo_data/DATA/grid_tpxo7.2.nc'

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





7. Generate boundary conditions with PyNEMO: Run PyNMEO
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

::

  ssh -Y espp1
  module load anaconda
  source activate pynemo_env

  export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  export PYTHONPATH=~/.conda/envs/pynemo_env/lib/python2.7/site-packages:$PYTHONPATH
  cd $WDIR/INPUTS
  pynemo -g -s namelist.bdy

Once the area of interest is selected and the close button is clicked, open
boundary data should be generated in $WDIR/OUTPUT

The SAVE button only updates the ``namelist.bdy`` file. The CLOSE button activates the process.

8. Run the configuration
++++++++++++++++++++++++
