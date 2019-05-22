Generate Initial conditions from CMEMS data
+++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++++++++++++++++++++++++

*(17 May 2019)*

This is based on a `more general recipe <generate_initial_conditions.rst>`_

To use initial conditions from an existing CMEMS T,S field you might need to do a bit
of interpolation. It is advisable to let NEMO do the heavy lifting for vertical
interpolation (requiring some FORTRAN modifications in MY_SRC), though SOSIE tools can be user
to do simple horizontal interpolation and also to make sure the initial conditions
have a matching number of vertical levels (or the NEMO interpolation doesn't work).


Building T,S field initial conditions from existing fields
==========================================================

Build 3D initial conditions for T/S. Velocities will start from rest.


Outline:

* Cut out rough domain T and S.
* Use SCRIP tools to remap onto configurations horizontal coords
* Use SOSIE to remove land by extrapolating water laterally.
* Interpolate on the fly in NEMO to convert z-level to hybrid coords.


Set some ARCHER paths
---------------------

::

  export CONFIG=BoBEAS
  export WORK=/work/n01/n01
  export WDIR=$WORK/jelt/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/NAMELISTS_AND_FORTRAN_FILES
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00


Git clone the repositoty
------------------------

Clone the SCRIPTS and NAMELISTS_AND_FORTRAN_FILES folders::

  cd $WORK/$USER
  git clone https://github.com/NOC-MSM/BoBEAS.git



Rough cut some initial conditions from parent (global) dataset
--------------------------------------------------------------

Make/get cut down parent file using CMEMS Mercator 1/12 ocean data.
Define a rough cut box::

  longitude: 60, 110
  latitude: 0, 30

NB Aim: Get the data (1st Apr 2019 - 10 May 2019).
Get 1st Apr 2019 for initial conditions. In the following a file CMEMS_01042019_T.nc
is generated with 3D T, S for the cut out (wish written it as YYMMDD, Doh)::

  livljobs4:
  /projectsa/accord/BoBEAS/INPUTS
  python -m pip install motuclient

  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-04-01 12:00:00" --date-max "2019-04-01 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable thetao --variable so --out-name CMEMS_01042019_T_download.nc --user jpolton --pwd JeffPCMEMS2018

Unpack the netcdf offsets and compression (SCRIP tools don't like it). Copy
 parent file to ARCHER INPUTS::

  livljobs4
  cd /projectsa/accord/BoBEAS/INPUTS/
  ncpdq -U CMEMS_01042019_T_download.nc CMEMS_01042019_T.nc
  scp /projectsa/accord/BoBEAS/INPUTS/CMEMS_01042019_T.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/BoBEAS/INPUTS/.



Use SCRIP tools to remap to the new grid
----------------------------------------

Now do interpolation onto child lateral grid.  The ``scrip`` tools are built in ``TDIR``
with a few modifications see
 `REPO:SCRIPTS/make_tools.sh <https://github.com/NOC-MSM/BoBEAS/blob/master/SCRIPTS/make_tools.sh>`_.
Then
::

  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/

First copy the namelists::

  cp $START_FILES/INITIAL_CONDITION/namelist_reshape_bilin_initcd_votemper $INPUTS/.
  cp $START_FILES/INITIAL_CONDITION/namelist_reshape_bilin_initcd_vosaline $INPUTS/.

Edit the input files::

  vi $INPUTS/namelist_reshape_bilin_initcd_votemper
  &grid_inputs
    input_file = 'CMEMS_01042019_T.nc'
  ...
    input_name = "thetao"

  &interp_inputs
    input_file = "CMEMS_01042019_T.nc"
    ...
    input_vars = "depth", "time"


Similarly for the *vosaline.nc file::

  vi $INPUTS/namelist_reshape_bilin_initcd_vosaline
  &grid_inputs
    input_file = 'CMEMS_01042019_T.nc'
    ...
    input_name = "vosaline"
  ...

  &interp_inputs
    input_file = 'CMEMS_01042019_T.nc'
    ...
    input_vars = "depth", "time"




Produce the remap files::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates ``data_nemo_bilin_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates ``initcd_origz_votemper.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates ``initcd_origz_vosaline.nc``.

These files have the same vertical grid as the parent data.
---



Use SOSIE tools to flood fill the parent initial conditions
-----------------------------------------------------------

Interpolating the T,S on z-levels onto hybrid levels can create water where
there was previously only land. Convert all the land in the parent initial conditions
to water by "flooding" the domain. This can be done with the SOSIE tool.

Before building and using the tool first make a land mask file to tell the SOSIE
what needs flooding. Use the salinity field to do this since we know the
salinity field is zero on land. Using NCO tools (mask out the fresh coastal water
as it makes a mess of the flood filling and subsequent z-interpolation)::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  ncks -d time_counter,0,0,1 -v vosaline initcd_origz_vosaline.nc initcd_origz_mask.nc
  ncap2 -O -s 'where(vosaline <=30.) vosaline=0' initcd_origz_mask.nc initcd_origz_mask.nc
  ncap2 -O -s 'where(vosaline >0.) vosaline=1' initcd_origz_mask.nc initcd_origz_mask.nc
  ncrename -v vosaline,mask initcd_origz_mask.nc

Restore modules::

  module unload nco/4.5.0
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

This has created a file ``initcd_origz_mask.nc`` with a variable ``mask``.

Now build the SOSIE tool.
Copy ``make.macro`` file and edit the path if necessary::

  cp $START_FILES/make.macro /home/n01/n01/jelt/sosie/.

  vi /home/n01/n01/jelt/sosie/make.macro
  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/jelt/local

Install. This might be best done in a clean terminal::

  cd ~
  mkdir local
  git clone https://github.com/brodeau/sosie.git
  cd sosie

  make
  make install
  export PATH=~/local/bin:$PATH
  cd $WDIR/INPUTS

Obtain the fields to interpolate. E.g interpolate CMEMS, AMM60 or ORCA
data. Get the namelists::

  cp $START_FILES/INITIAL_CONDITION/initcd_votemper.namelist $INPUTS/.
  cp $START_FILES/INITIAL_CONDITION/initcd_vosaline.namelist $INPUTS/.

The sosie routine is VERY slow first time round (1hr). This is when it
makes a ``sosie_mapping`` file that can be reused for other variables.

It is advisable to let NEMO do the details of vertical interpolation. Use SOSIE
 tools for the flood filling and for getting the *same number of levels* as
 appear in the child grid.

 Edit namelists to the variables you want::

  vi initcd_vosaline.namelist
  &ninput
  ivect     = 0
  lregin    = T
  cf_in     = 'initcd_origz_vosaline.nc'
  cv_in     = 'vosaline'
  cv_t_in   = 'time_counter'
  jt1       = 0
  jt2       = 0
  jplev     = 0
  cf_x_in   = 'initcd_origz_vosaline.nc'
  cv_lon_in = 'x'
  cv_lat_in = 'y'
  cf_lsm_in = 'initcd_origz_mask.nc'
  cv_lsm_in = 'mask'
  ldrown    = T
  ...

  &n3d
  cf_z_in  = 'initcd_origz_vosaline.nc'
  cv_z_in  = 'gdept'
  cf_z_out = 'domain_cfg.nc'
  cv_z_out = 'nav_lev'
  cv_z_out_name = 'gdept'
  ctype_z_in = 'z'
  ctype_z_out = 'z'
  /


  &nhtarget
  lregout    = F
  cf_x_out   = 'initcd_origz_vosaline.nc'
  cv_lon_out = 'x'
  cv_lat_out = 'y'
  cf_lsm_out = ''
  cv_lsm_out = ''
  lmout      = F

  &noutput
  cmethod  = 'bilin'
  cv_t_out = 'time_counter'
  cv_out   = 'vosaline'
  cu_out   = 'PSU'
  cln_out  = 'Salinity'
  cd_out   = '.'
  !!
  csource  = 'CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024'
  ctarget  = 'BoBEAS'
  /

Similarly for ``initcd_votemper.namelist``::

  vi initcd_votemper.namelist

  vosaline --> votemper
  ...
  cu_out   = 'C'
  cln_out  = 'Temperature'


Executing SOSIE tools is fine in interactive mode if you already have generated
the sosie_mapping file. (I.e. run it once before). For the first run I had to submit
it as a serial job  **IT TOOK 1hrs 1m**

PBS submission script::

  cd $INPUTS
  vi sosie_initcd_T

  #!/bin/bash
  #PBS -N init_T
  #PBS -l select=serial=true:ncpus=1
  #PBS -l walltime=06:00:00
  #PBS -o init_T.log
  #PBS -e init_T.err
  #PBS -A n01-ACCORD
  ###################################################

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel


  cd /home/n01/n01/jelt/sosie
  make clean
  make
  make install

  #set up paths
  cd /work/n01/n01/jelt/BoBEAS/INPUTS

  /home/n01/n01/jelt/local/bin/sosie.x -f initcd_votemper.namelist
  /home/n01/n01/jelt/local/bin/sosie.x -f initcd_vosaline.namelist

  # qsub -q serial <filename>
  ###################################################


Launch job::

  qsub -q serial sosie_initcd_T

Subsequent jobs could be in interactive mode::

  #sosie.x -f initcd_vosaline.namelist
  #sosie.x -f initcd_votemper.namelist

Whether as a serial job or from the command line, the temperature process creates::

  sosie_mapping_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-BoBEAS.nc
  votemper_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-BoBEAS_Apr2019.nc
  vosaline_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-BoBEAS_Apr2019.nc

Where the ``sosie_mapping*.nc`` file is created by the first sosie iteration only.

Check these fields are OK.

---
By this stage should have initial conditions T and S files ``votemper_*_Apr2019.nc``
and ``vosaline_*_Apr2019.nc`` on the configurations horizontal grid
and on a z-level grid with the same number of levels as the target. The z-levels are
not on the target's hybrid vertical coordinates. These will be the initial conditions used.
NEMO can do on the fly vertical interpolation.
It might be convenient to sym link them to::

   ln -s vosaline_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-BoBEAS_Apr2019.nc initcd_vosaline.nc
   ln -s votemper_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-BoBEAS_Apr2019.nc initcd_votemper.nc

NB These two files will have to be linked into the ICS dir in the EXP dir.


Interpolate in z on the fly
===========================

For vertical interpolation we let NEMO do the heavy lifting. This requires some changes
to the FORTRAN using ``par_oce.F90`` and ``dtatsd.F90`` in ``MY_SRC``. See
`<build_opa_orchestra.rst>`_

Maybe move the executable to something memorable e.g.::

  cd $CDIR
  mv $CONFIG/BLD/bin/nemo.exe $CONFIG/BLD/bin/nemo_tide_nomet.exe

To interpolate the initial conditions on-the-fly need to pass information to
NEMO about the parent vertical grid and parent mask file. Appropriate variables
are created in external files that are read into the namelist.

These mask and depth variables need to be 4D variables, where length(t)=1.
They can be created with NCO tools by manipulating a parent initial condition file.
On archer, load the appropriate modules::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

If the depth (gdept) variable is 1D and the file has dimensions
[time,z,y,x] then first we make it 3D and call it something like gdept_3D::

  cd $INPUTS/
  ncap2 -O -s 'gdept_3D[z,y,x]=gdept' initcd_votemper.nc tmp.nc

Then add a time dimension::

  ncap2 -O -s 'gdept_4D[time_counter,z,y,x]=gdept_3D' tmp.nc initcd_depth.nc
  rm tmp.nc

For the mask variable use one of the tracer variables (in this case salinity
 and we know the land values are set to zero). NB if following progressively,
 a similar mask file (with a different limit salinity and potentially different
 number of vertical levels) was created just before the SOSIE step
 ``initcd_origz_mask.nc``. A mask with the correct number of vertical levels is needed
 of the number of levels changes between parent and child::

  ncks -d time_counter,0,0,1 -v vosaline initcd_vosaline.nc initcd_mask.nc
  #ncap2 -O -s 'where(vosaline <=0.) vosaline=0' initcd_mask.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline <=0.) vosaline=1' initcd_mask.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline >0.) vosaline=1' initcd_mask.nc initcd_mask.nc
  ncrename -v vosaline,mask initcd_mask.nc

.. note: Changed the above so that ALL the values are 1. (I.e. a rubbish mask).
  The problem was that in the child bathymetry some of the sea mounts have moved
  and so using a mask from the old grid meant sst where pulled down to deep water.
  Since the parent is flood filled there is no need for a mask anyway. Easiest to
  set all values to one instead of disactivating.

Restore modules::

  module unload nco/4.5.0
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

The resulting files are ``initcd_mask.nc`` and ``initcd_depth.nc`` which are read
into the namelist.

NB These two files will have to be linked into the ICS dir in the EXP dir.

Edit, or add, new **mask** and **depth** variables to the namelist_cfg. Also
add the logical switch to do vertical interpolation ``ln_tsd_interp=T``::

  cd $EXP/../EXP_Apr19
  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namtsd        !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  !              !  file name                 ! frequency (hours) ! variable ! time interp.!  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                            !  (if <0  months)  !   name   !  (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
  sn_tem  = 'initcd_votemper.nc',         -12        ,'votemper' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_sal  = 'initcd_votemper.nc',         -12        ,'vosaline' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_dep  = 'initcd_depth.nc'   ,         -12        ,'gdept_4D',   .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''
  sn_msk  = 'initcd_mask.nc',       -12        ,'mask',       .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''

    !
     cn_dir        = './ICS/'     !  root directory for the location of the runoff files
     ln_tsd_init   = .true.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_interp = .true.    !  Interpolation of T & S in the verticalinput data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)

.. Note: Can interpolate the initcd_fields in time if that is appropriate. Can in
 principle apply a weightings file so that the initcd_field file are uncut parent grid
 data at some other resolution.

 However, do not do use the weights files to perform horizontal interpolation combined
 with  ln_tsd_interp = .true. to perform vertical interpolation as the mask file
 will be rendered useless! If you are going to take this approach flood-fill all
 the land and then set the mask array to equal 1 everywhere. That way it won’t be
 corrupted when using the weights files to interpolate onto the child grid.


For some reason the initial condition files generated are statically unstable...
This can be resolved using NEMO machinery with non penetrative convective
adjustment. E.g.::

  &namzdf
  ln_zdfevd   = .true. --> .false.
  ln_zdfnpc   = .false. --> .true.
