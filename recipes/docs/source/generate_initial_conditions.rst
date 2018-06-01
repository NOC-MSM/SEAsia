Generate Initial conditions
+++++++++++++++++++++++++++
+++++++++++++++++++++++++++

**10 May 2018: THIS IS NOT PROPERLY GENERALISED BECAUSE I'VE ONLY DONE IT ONCE!**

For a new configuration you probably want to start with idealised, or homogenous
initial conditions. This is done with user defined initial conditions ``ln_usr``
with the expression being compiled into the executable

To use initial conditions from an existing T,S field you might need to do a bit
of interpolation. It is advisable to let NEMO do the heavy lifting for vertical
interpolation (rquiring some FORTRAN modifictions), though SOSIE tools can be user
to do simple horizontal interpolation.


User defined initial initial conditions
=======================================

For constant T and S use the user defined functions in ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``. Compile and save executable with
  telegraphic names that point to compile options. e.g.::

    nemo_notide_TSprofile.exe
    nemo_tideonly_TSconst.exe


Building T,S field initial conditions from existing fields
==========================================================

Second time around we build 3D initial conditions
*(27 Apr 2018)*

*Since my parent and child are on the same grid I'm not sure I need all these steps.
Indeed there is almost certainly a more efficient method with hindsight. However
 I am marching onwards*

Outline:

* Cut out rough domain T and S.
* Use SCRIP tools to remap onto configurations horizontal coords
* Use SOSIE to remove land by extrapolating water laterally.
* Interpolate on the fly in NEMO to convert z-level to hybrid coords.

Rough cut some initial conditions from parent (global) dataset
--------------------------------------------------------------

Make cut down parent file using ORCA0083-N01.
Copy parent file to ARCHER INPUTS (need to generalise / improve)::

  livljobs4
  scp /projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.

Cut down based on coordintaes from *create coordinates* namelist. (Add a bit of
a buffer)::

    module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
    module load cray-netcdf cray-hdf5
    module load nco/4.5.0
    cd $WDIR/INPUTS

    ncks -d x,45,735 -d y,1245,1810 ORCA0083-N01_19791101d05T.nc $WDIR/INPUTS/cut_down_19791101d05_SEAsia_grid_T.nc

Average over time and restore the parallel modules (Not necessary for this data with 1 time point)::

    #e.g. ncwa -a time_counter $WDIR/INPUTS/cut_down_20131013_LBay_grid_T.nc  $WDIR/INPUTS/cut_down_201310_LBay_grid_T.nc

    module unload nco cray-netcdf cray-hdf5
    module load cray-netcdf-hdf5parallel cray-hdf5-parallel



Use SCRIP tools to remap to the new grid
----------------------------------------

Now do interpolation onto child lateral grid.  The ``scrip`` tools are build in ``TDIR``
e.g. in `Build Tools<SEAsia_archer_livljobs4.rst>`_
::

  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/

First copy the namelists::

  cp $START_FILES/namelist_reshape_bilin_initcd_votemper $INPUTS/.
  cp $START_FILES/namelist_reshape_bilin_initcd_vosaline $INPUTS/.

Edit the input files::

  vi $INPUTS/namelist_reshape_bilin_initcd_votemper
  &grid_inputs
    input_file = 'cut_down_19791101d05_SEAsia_grid_T.nc'
  ...
    input_name = "votemper"

  &interp_inputs
    input_file = "cut_down_19791101d05_SEAsia_grid_T.nc"
    ...
    input_vars = "deptht", "time_counter"


Similarly for the *vosaline.nc file::

  vi $INPUTS/namelist_reshape_bilin_initcd_vosaline
  &grid_inputs
    input_file = 'cut_down_19791101d05_SEAsia_grid_T.nc'
    ...
    input_name = "vosaline"
  ...

  &interp_inputs
    input_file = 'cut_down_19791101d05_SEAsia_grid_T.nc'
    ...
    input_vars = "deptht", "time_counter"




Produce the remap files::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates ``data_nemo_bilin_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates ``initcd_votemper.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates ``initcd_vosaline.nc``.

---



Use SOSIE tools to flood fill the parent initial conditions
-----------------------------------------------------------

Interpolating the T,S on z-levels onto hybrid levels can create water where
there was previously only land. Convert all the land in the parent initial conditions
to water by `flooding` the domain. This can be done with the SOSIE tool.

Before building and using the tool first make a land mask file to tell the SOSIE
what needs flooding. Use the salinity field to do this since we know the
salinity field is zero on land. Using NCO tools (mask out the fresh coastal water
as it makes a mess of the flood filling and subsequent z-interpolation)::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  ncks -d time_counter,0,0,1 -v vosaline initcd_vosaline.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline <=30.) vosaline=0' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline >0.) vosaline=1' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncrename -v vosaline,mask sosie_initcd_mask.nc

Restore modules::

  module unload nco/4.5.0
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

This has created a file ``initcd_mask`` with a variable ``mask``.

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

Obtain the fields to interpolate. E.g interpolate AMM60 or ORCA
data. Get the namelists::

  cp $START_FILES/initcd_votemper.namelist $INPUTS/.
  cp $START_FILES/initcd_vosaline.namelist $INPUTS/.

The sosie routine is VERY slow first time round (4hr 25 mins). This is when it
makes a ``sosie_mapping`` file that can be reused for other variables.

It is advisable to let NEMO do the vertical interpolation so only use SOSIE
 tools for the flood filling. Though it can do other things.

 Edit namelists to the variables you want::

  vi initcd_vosaline.namelist
  &ninput
  ivect     = 0
  lregin    = F
  cf_in     = 'initcd_vosaline.nc'
  cv_in     = 'vosaline'
  cv_t_in   = 'time_counter'
  jt1       = 0
  jt2       = 0
  jplev     = 0
  cf_x_in   = 'initcd_vosaline.nc'
  cv_lon_in = 'x'
  cv_lat_in = 'y'
  cf_lsm_in = 'sosie_initcd_mask.nc'
  cv_lsm_in = 'mask'
  ldrown    = T
  ...

  &n3d
  cf_z_in  = 'initcd_vosaline.nc'
  cv_z_in  = 'gdept'
  cf_z_out = 'initcd_vosaline.nc'
  cv_z_out = 'gdept'
  cv_z_out_name = 'gdept'
  ctype_z_in = 'z'
  ctype_z_out = 'z'
  /


  &nhtarget
  lregout    = F
  cf_x_out   = 'initcd_vosaline.nc'
  cv_lon_out = 'x'
  cv_lat_out = 'y'
  cf_lsm_out = ''
  cv_lsm_out = ''
  lmout      = F

  &noutput
  cmethod  = 'bilin'
  cv_t_out = 'time_counter'
  cv_out   = 'vosaline'
  cu_out   = 'psu'
  cln_out  = 'Salinity'
  cd_out   = '.'
  !!
  csource  = 'ORCA0083-N01'
  ctarget  = 'SEAsia'
  cextra   = '1978'
  /

Similarly for ``initcd_votemper.namelist``::

  vi initcd_votemper.namelist

  vosaline --> votemper
  ...
  cu_out   = 'C'
  cln_out  = 'Temperature'


Executing SOSIE tools is fine in interactive mode if you already have generated
the sosie_mapping file. (I.e. run it once before). For the first run I had to submit
it as a serial job **IT TOOK 4hrs 25m TO DO 3D**

PBS submission script::

  cd $INPUTS
  vi $INPUTS/sosie_initcd_T

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
  cd /work/n01/n01/jelt/SEAsia/INPUTS

  /home/n01/n01/jelt/local/bin/sosie.x -f initcd_votemper.namelist

  # qsub -q serial <filename>
  ###################################################


Subsequent jobs can be in interactive mode::

  qsub -q serial sosie_initcd_T
  sosie.x -f initcd_vosaline.namelist
  #sosie.x -f initcd_votemper.namelist

Whether as a serial job or from the commandline, the temperature process creates::

  sosie_mapping_ORCA0083-N01-SEAsia.nc
  votemper_ORCA0083-N01-SEAsia_1978.nc

And the salinity process creates::

  vosaline_ORCA0083-N01-SEAsia_1978.nc

Check these fields are OK.

---
By this stage should have initial conditions T and S files ``votemper_ORCA0083-N01-SEAsia_1978.nc``
and ``vosaline_ORCA0083-N01-SEAsia_1978.nc`` on the configurations horizontal grid
and on the ORCA parent z-level grid.


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

  ncap2 -O -s 'gdept_4D[time_counter,z,y,x]=gdept_4D' tmp.nc initcd_depth.nc
  rm tmp.nc

For the mask variable use one of the tracer variables (in this case salinity
 and we know the land values are set to zero). NB if following progressively,
 a similar mask file (with a different limit salinity) was created in the SOSIE step
 **DONT NEED THIS**::

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

Edit, or add, new **mask** and **depth** variables to the namelist_cfg. Also
add the logical switch to do vertical interpolation ``ln_tsd_interp=T``::

  cd $EXP/../EXP_tide_initcd
  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namtsd        !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  !              !  file name                 ! frequency (hours) ! variable ! time interp.!  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                            !  (if <0  months)  !   name   !  (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
  sn_tem  = 'votemper_ORCA0083-N01-SEAsia_1978.nc',         -12        ,'votemper' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_sal  = 'vosaline_ORCA0083-N01-SEAsia_1978.nc',         -12        ,'vosaline' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_dep  = 'initcd_depth.nc'   ,         -12        ,'gdept_4D',   .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''
  sn_msk  = 'initcd_mask.nc'    ,         -12        ,'mask',       .false.   , .true. , 'yearly'   , ''  ,    ''    ,      ''

    !
     cn_dir        = '../../../../INPUTS/'     !  root directory for the location of the runoff files
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
