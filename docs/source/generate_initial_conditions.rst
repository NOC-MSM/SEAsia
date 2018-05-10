Generate Initial conditions
+++++++++++++++++++++++++++
+++++++++++++++++++++++++++

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

*Since my parent and child are on the same grid I'm not sure I need the SOSIE
step. Try skipping it for now. You would want to use SOSIE to flood fill the land
if you are horizontally interpolating*

Jump straight to using SCRIP tools to remap parent fields onto child grid


Use SOSIE tools to interpolate onto a finer grid
------------------------------------------------

Copy ``make.macro`` file and edit the path if necessary::
**FIX** to the notes (copied from jdha instead): ``cp $WDIR/INPUTS/make.macro ./``::

  cp /home/n01/n01/jdha/sosie/make.macro /home/n01/n01/jelt/sosie/.

  vi /home/n01/n01/jelt/sosie/make.macro
  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/jelt/local

Proceed with Step 6 (of Lighhouse Reef Readthedocs). This is best done in a clean terminal
::

  cd ~
  mkdir local
  svn co svn://svn.code.sf.net/p/sosie/code/trunk sosie
  cd sosie

  make
  make install
  export PATH=~/local/bin:$PATH
  cd $WDIR/INPUTS


Obtain the fields to interpolate. E.g interpolate AMM60 or ORCA
data. Get the namelists::

  cp $START_FILES/initcd_votemper.namelist $INPUTS/.
  cp $START_FILES/initcd_vosaline.namelist $INPUTS/.

The sosie routine is VERY slow. (2.5 hrs)
*UPDATE 30 Apr 18: possible because I was doing unnecessary 3D interpolation*.
So, make a cut down parent file using ORCA0083-N01.
Cut down based on coordintaes from create coordinates namelist. (Add a buffer as
I'm not sure how the sosie extraction works)::

    module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
    module load cray-netcdf cray-hdf5
    module load nco/4.5.0
    cd $WDIR/INPUTS

    ncks -d x,45,735 -d y,1245,1810 ORCA0083-N01_19791101d05T.nc $WDIR/INPUTS/cut_down_19791101d05_SEAsia_grid_T.nc

Average over time and restore the parallel modules (Not necessary for this data with 1 time point)::

    #ncwa -a time_counter $WDIR/INPUTS/cut_down_20131013_LBay_grid_T.nc  $WDIR/INPUTS/cut_down_201310_LBay_grid_T.nc

    module unload nco cray-netcdf cray-hdf5
    module load cray-netcdf-hdf5parallel cray-hdf5-parallel



Edit namelists to the variables you want.
It is advisable to let NEMO do the vertical interpolation so only use SOSIE tools for 2D
interpolation.

In this exmaple I am now pointlessly interpolating something onto itself. This is a template and would
be useful when I have a child grid at a finer resolution than the parent.

Set jplev = 1
To interpolate only on the first level (and output level..). **NOT SURE HOW YOU FIX IT TO DO 2D INTERPOLATION**
::

  vi initcd_votemper.namelist
  cf_in     = 'cut_down_19791101d05_SEAsia_grid_T.nc'
  cv_in     = 'votemper'
  cf_x_in   = 'cut_down_19791101d05_SEAsia_grid_T.nc'
  cv_out   = 'votemper'
  csource  = 'ORCA0083-N01'
  ctarget  = 'SEAsia'

  vi initcd_vosaline.namelist
  ...
  cv_out   = 'vosaline'
  ...

Copy parent file to ARCHER INPUTS (need to generalise / improve)::

  livljobs4
  scp /projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.


Do stuff (on ARCHER). I think the intention was for SOSIE to flood fill the land::

  cd $INPUTS
  sosie.x -f initcd_votemper.namelist
  sosie.x -f initcd_vosaline.namelist


I had trouble getting ARCHER to run this. (Though with the cut down parent the
 commandline is fine, though it runs out of walltime in Interactive Mode)
Initially, running in the commandline the job failed with insufficient memory,
 because I didn't cut the ORCA data down first.
 In the end I submitted it as a pair of serial jobs. **IT TOOK 4hrs 25m TO DO 3D**::

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
  #/home/n01/n01/jelt/local/bin/sosie.x -f initcd_vosaline.namelist


  # qsub -q serial <filename>
  ###################################################


Similarly for ``sosie_initcd_S``. Then::

  qsub -q serial sosie_initcd_T
  qsub -q serial sosie_initcd_S

3 hours not enough - resubmit with 6 hrs! *It took 4h 25min*

Whether as a serial job or from the commandline, the temperature process creates::

  sosie_mapping_ORCA0083-N01-SEAsia.nc
  votemper_ORCA0083-N01-SEAsia_1978.nc4

And the salinity process creates::

  vosaline_ORCA0083-N01-SEAsia_1978.nc4


Use SCRIP tools to remap to the new grid
----------------------------------------

Now do interpolation onto child grid.  The ``scrip`` tools are build in ``TDIR``
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

For the mask variable use one of the tracer variables (in this case salinity and we know the land values are set to zero)::

  ncks -d time_counter,0,0,1 -v vosaline initcd_vosaline.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline <=0.) vosaline=0' initcd_mask.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline >0.) vosaline=1' initcd_mask.nc initcd_mask.nc
  ncrename -v vosaline,mask initcd_mask.nc

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
  sn_tem  = 'initcd_votemper.nc',         -12        ,'votemper' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
  sn_sal  = 'initcd_vosaline.nc',         -12        ,'vosaline' ,  .false.   , .true. , 'yearly'   , ''   ,   ''    ,    ''
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
