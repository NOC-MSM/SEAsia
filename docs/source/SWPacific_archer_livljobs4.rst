=============================================
Setting up a SW Pacific NEMO v4 configuration
=============================================

Machines: livljobs4, ARCHER

URL:: *to add*
*Actually I've taken this off readthedocs as it is internal. However writing in*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a SW Pacific regional tide-only model using GEBCO bathymetry, FES tidal
boundaries.

Build on a combination of livljobs4 and ARCHER.

Uses a prerelease of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. ARCHER: Run simulation

It is a tide only run.

Issues that arose
=================

* ...

.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO bits are done on livljobs4.
(There was a problem with some java bits working on ARCHER)
Starting on ARCHER::

  ssh login.archer.ac.uk

  export CONFIG=SWPacific
  export WORK=/work/n01/n01
  export WDIR=$WORK/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export OLD_TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/
  export EXP=$CDIR/$CONFIG/EXP00

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel


---

.. note: you might have to mkdir the odd directory or two.

Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

---

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

Or just link XIOS executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---


Build TOOLS
===========

To generate bathymetry, initial conditions and grid information we first need
to compile some of the NEMO TOOLS. NB for some reason GRIDGEN doesn’t like INTEL

.. note: These are compiled with XIOS2. However DOMAINcfg has to be compiled
  with XIOS1. There is a README in the $TDIR/DOMAINcfg on what to do.

First build DOMAINcfg (which is relatively new and in NEMOv4). Use my XIOS1 file
(see userid and path in variable ``%XIOS_HOME``). Copy from ARCH *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $CDIR/../ARCH/.
  cd $TDIR

  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg

.. note: Check which arch file this is. Surely should be consistent.
   Though I don't attempt to change the GRIDGEN build
::

  ./maketools -n WEIGHTS -m XC_ARCHER_INTEL_XIOS1
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL_XIOS1

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module swap PrgEnv-intel PrgEnv-cray
  module load cray-netcdf cray-hdf5
  ./maketools -n GRIDGEN -m XC_ARCHER    # This uses acc's XIOS_r474

  module unload cray-netcdf cray-hdf5
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel





1. Generate new coordinates file
++++++++++++++++++++++++++++++++

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements.

.. warning:
  Using the GRIDGEN/create_coordinates.exe tool runs into a problem for zoom factor
  >1, since the horizontal spacing metric e.g. e[12]t always match
  the parent grid. I think that this is a bug. The agrif version works.

Build and execute agrif version of create_coordinates.exe.
See `build_and_create_coordinates.rst`_

This creates a new coordinatesfile with contents, which is now copied to
INPUTS::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

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

Execute first scrip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Output files::

  remap_nemo_grid_gebco.nc
  remap_data_grid_gebco.nc

Execute second scip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Output files::

  data_nemo_bilin_gebco.nc

Execute third scip thing (use old tools - already compiled)::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

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

Proceed with Step 6 (of Lighhouse Reef Readthedocs)::

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

  cp $START_FILES/initcd_votemper.namelist $INPUTS/.
  cp $START_FILES/initcd_vosaline.namelist $INPUTS/.

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
  cd $INPUTS

  ncks -d x,560,620 -d y,720,800 /work/n01/n01/kariho40/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/AMM60smago/EXP_notradiff/OUTPUT/AMM60_5d_20131013_20131129_grid_T.nc $INPUTS/cut_down_20131013_LBay_grid_T.nc

Average over time and restore the parallel modules::

  ncwa -a time_counter $START_FILES/cut_down_20131013_LBay_grid_T.nc  $INPUTS/cut_down_201310_LBay_grid_T.nc

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

  cp $START_FILES/namelist_reshape_bilin_initcd_votemper $INPUTS/.
  cp $START_FILES/namelist_reshape_bilin_initcd_vosaline $INPUTS/.

Edit the input files::

  vi $INPUTS/namelist_reshape_bilin_initcd_votemper
  &grid_inputs
    input_file = 'thetao_AMM60-LBay_2013.nc4'
  ...

  &interp_inputs
    input_file = "thetao_AMM60-LBay_2013.nc4"
  ...

Simiarly for the *vosaline.nc file::

  vi $INPUTS/namelist_reshape_bilin_initcd_vosaline
  &grid_inputs
    input_file = 'so_AMM60-LBay_2013.nc4'
  ...

  &interp_inputs
    input_file = "so_AMM60-LBay_2013.nc4"
  ...


Produce the remap files::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates ``data_nemo_bilin_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates ``initcd_votemper.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates ``initcd_vosaline.nc``.



4. Generate a domain configuration file
=======================================

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg``
directory along with all the inputs files that would have previously been needed
get v3.6 running. The reason being that all the non-time stepping stuff, like
grid generating, has been abstracted from the core OPA code and is now done as
a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory::

    cp $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
    cp $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template namelist_cfg with only the essenetial domain building stuff::

  cd $TDIR/DOMAINcfg
  vi namelist_cfg

  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  Configuration namelist : used to overwrite defaults values defined in SHARED/namelist_ref
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     nn_no       =       0   !  job number (no more used...)
     cn_exp      =  "domaincfg"  !  experience name
     nn_it000    =       1   !  first time step
     nn_itend    =      75   !  last  time step (std 5475)
  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     !
     ln_e3_dep   = .false.    ! =T : e3=dk[depth] in discret sens.
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !
     cp_cfg      =  "orca"   !  name of the configuration
     jp_cfg      =       0   !  resolution of the configuration
     jpidta      =      53   !  1st lateral dimension ( >= jpi )
     jpjdta      =     109   !  2nd    "         "    ( >= jpj )
     jpkdta      =      51   !  number of levels      ( >= jpk )
     jpiglo      =      53   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     109   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_sco      = .true.    !  z-coordinate - partial steps
     ln_linssh   = .false.    !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     jphgr_msh   =       0               !  type of horizontal mesh
  ...

.. note:

  No gdept output in the offical v4 release. Though it was acheived here setting
  ln_e3_dep = F. This is needed for PyNEMO, though could be constructed from e3[tw].

Build a script to run the executable::

  vi $TDIR/DOMAINcdf/rs

  #!/bin/bash
  #PBS -N domain_cfg
  #PBS -l walltime=00:20:00
  #PBS -l select=1
  #PBS -j oe
  #PBS -A n01-NOCL
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk
  #! -----------------------------------------------------------------------------

  # Change to the directory that the job was submitted from
  cd $PBS_O_WORKDIR

  # Set the number of threads to 1
  #   This prevents any system libraries from automatically
  #   using threading.
  export OMP_NUM_THREADS=1
  # Change to the directory that the job was submitted from
  ulimit -s unlimited

  #===============================================================
  # LAUNCH JOB
  #===============================================================
  echo `date` : Launch Job
  aprun -n 1 -N 1 ./make_domain_cfg.exe >&  stdouterr_cfg

  exit


Try running it::

  cd $TDIR/DOMAINcfg
  qsub -q short rs


Put a copy in $INPUTS for safe keeping::

    cp $TDIR/DOMAINcfg/namelist_cfg $INPUTS/namelist_cfg_generateDOMAINcfg

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

  cp $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  cp $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.




5. Generate weights for atm forcing
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

  cp $START_FILES/namelist_reshape_bilin_atmos $INPUTS/.
  cp $START_FILES/namelist_reshape_bicubic_atmos $INPUTS/.

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

  cd $INPUTS
  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos

Generate  remap files ``remap_nemo_grid_atmos.nc`` and ``remap_data_grid_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos

Generates ``data_nemo_bilin_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos

Generates ``weights_bilinear_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos

Generates ``data_nemo_bicubic_atmos.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

Generates ``weights_bicubic_atmos.nc``.



THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO *(20 Sept 2017)*
If all the files are ready to go jump straight to `7. Generate boundary conditions with PyNEMO: Run PyNEMO`_


6. Generate boundary conditions with NRCT/PyNEMO: Create netcdf abstraction wrapper
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this section there are two stages.
* generate a ncml file which describes the files needed to create boundary conditions
* generate a namelist.bdy file which controls the actual boundary condition generation.

For each parent data set a new pair of (``*.ncml``, ``namelist.bdy``) are needed.
Here I attempt to use parent data from NNA. I could use data from:
* AMM60 local data (might not yet work because of the sigma levels)
* thredds server (as in the LH_REEF example, though this is turned off!)
* NNA local data (easiest ?)

First install PyNEMO if not already done so. Full description (If this is already
installed then follow through anyway but skip the mkdir / create / install / clone
 and build commands)::

  ssh -Y livljobs4

  export CONFIG=LBay
  export WORK=/work
  export WDIR=$WORK/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  #export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  #export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  #export EXP=$CDIR/$CONFIG/EXP00

  cd $WORK/$USER
  mkdir $WDIR
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
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $INPUTS

.. note: It might be best to abstract the above into a separate recipe that deals
 with either livljobs4 or archer

I suggest managing the namelist.bdy file after the ``ncml`` file is generated.
 A fresh ``ncml`` file can be generated automatically or an existing one can be edited.


6a. Generate ncml files
+++++++++++++++++++++++

.. note: Not tested on livljobs4. It used to work on ARCHER, before switching to
 domain_cfg.nc files)

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

In ``$INPUTS`` I have three ncml files.
* One for using the thredds server to get remote ORCA12 data.
* One for using local AMM60 data, with ackward s-sigma levels
* One for using local NNA data

The first two are work in progress / templates. The latter is used here.

NNA_inputs_src.ncml
+++++++++++++++++++

Note need to set the time variables and new ``sn_src_dir`` in namelist.bdy.
Actually upated the following with all the Jan 2000 files::

  cd $INPUTS
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




6b. Generate the namelist.bdy file for PyNEMO / NRCT
+++++++++++++++++++++++++++++++++++++++++++++++++++

Copy the NRCT template namelist.bdy from the START_FILES::

  cd $INPUTS
  cp $START_FILES/namelist.bdy $INPUTS/.

.. note: There is an old namelist.bdy files namelist.bdy_old_mesh_files that does
 not rely on the new domain_cfg.nc file

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

Need to grab some INPUT files. (File bathy_meter.nc and domain_cfg.nc should be
 there already)::

   scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/domain_cfg.nc $INPUTS/.
   scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/bathy_meter.nc $INPUTS/.


Copy the pynemo namelist and ncml files::

  cp $START_FILES/namelist.bdy_NNA    $INPUTS/.
  cp $START_FILES/NNA_inputs_src.ncml $INPUTS/.
  cp $START_FILES/inputs_dst.ncml     $INPUTS/.
  cd $INPUTS

Make sure the NNA data is available::

  mkdir $WDIR/INPUTS/NNA
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mask.nc $WDIR/INPUTS/NNA/.
  for file in NNA_*200001*nc ; do scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/$file $WDIR/INPUTS/NNA/. ; done

.. note: I have not done this as a clean build with the new domain_cfg.nc files

namelist.bdy_NNA
++++++++++++++++

Edit namelist.bdy_NNA to reflect locally stored mesh and mask files. Also
inputs_dst.ncml. Set the date info back to (Nov?) 1979.

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
      sn_dst_hgr = './domain_cfg.nc'
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
       clname(1)      = 'M2'                  ! constituent name
       clname(2)      = 'S2'
       clname(3)      = 'N2'
       clname(4)      = 'K2'
       clname(5)      = 'K1'
       clname(6)      = 'O1'
       clname(7)      = 'P1'
       clname(8)      = 'Q1'
       clname(9)      = 'MF'
       clname(10)     = 'MM'
       clname(11)     = 'M4'
       clname(12)     = 'MS4'
       clname(13)     = 'MN4'
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

Also had to check that ``inputs_dst.ncml`` has the correct file name within:
 *Now domain_cfg.nc, formerly mesh_zgr.nc*. Note also that some variables in
  domain_cfg.nc have different names e.g. ``mbathy`` --> ``bottom_level``. Check the mapping
  in ``inputs_dst.ncml``::

   vi inputs_dst.ncml

   <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
     <ns0:aggregation type="union">
       <ns0:netcdf location="file:domain_cfg.nc">
       <ns0:variable name="mbathy" orgName="bottom_level" />
       <ns0:variable name="gdept" orgName="gdept_0" />
       <ns0:variable name="gdepw" orgName="gdepw_0" />
       <ns0:variable name="e3u" orgName="e3u_0" />
       <ns0:variable name="e3v" orgName="e3v_0" />
       </ns0:netcdf>
     </ns0:aggregation>
   </ns0:netcdf>

.. warning:
  In the actual v4 release domain_cfg.nc  will not have gdept or gdepw. These
  will need to be reconstructed from e3[tw].

Generate the boundary conditions again, with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  #pynemo -g -s namelist.bdy_NNA
  pynemo -s namelist.bdy_NNA

Once the area of interest is selected and the close button is clicked, open
boundary data should be generated in the current directory (NB I dont fiddle
with the GUI; I just click CLOSE to activiate, if everything is already sorted
in the input files).

The SAVE button only updates the ``namelist.bdy`` file. The CLOSE button activates the process.

Alternatively don't invoke the `-g` option and there is no GUI.

This generates::
  ls -1 $INPUTS

  coordinates.bdy.nc
  LBay_bdytide_rotT_M4_grid_T.nc
  LBay_bdytide_rotT_MM_grid_T.nc
  LBay_bdytide_rotT_MN4_grid_T.nc
  LBay_bdytide_rotT_MS4_grid_T.nc
  LBay_bdytide_rotT_M2_grid_T.nc
  LBay_bdytide_rotT_N2_grid_T.nc
  LBay_bdytide_rotT_S2_grid_T.nc
  LBay_bdytide_rotT_K1_grid_T.nc
  LBay_bdytide_rotT_K2_grid_T.nc
  LBay_bdytide_rotT_P1_grid_T.nc
  LBay_bdytide_rotT_O1_grid_T.nc
  LBay_bdytide_rotT_MF_grid_T.nc
  LBay_bdytide_rotT_Q1_grid_T.nc
  LBay_bdytide_rotT_M4_grid_U.nc
  LBay_bdytide_rotT_MM_grid_U.nc
  LBay_bdytide_rotT_MN4_grid_U.nc
  LBay_bdytide_rotT_MS4_grid_U.nc
  LBay_bdytide_rotT_M2_grid_U.nc
  LBay_bdytide_rotT_N2_grid_U.nc
  LBay_bdytide_rotT_S2_grid_U.nc
  LBay_bdytide_rotT_K1_grid_U.nc
  LBay_bdytide_rotT_K2_grid_U.nc
  LBay_bdytide_rotT_P1_grid_U.nc
  LBay_bdytide_rotT_O1_grid_U.nc
  LBay_bdytide_rotT_MF_grid_U.nc
  LBay_bdytide_rotT_Q1_grid_U.nc
  LBay_bdytide_rotT_M4_grid_V.nc
  LBay_bdytide_rotT_MM_grid_V.nc
  LBay_bdytide_rotT_MN4_grid_V.nc
  LBay_bdytide_rotT_MS4_grid_V.nc
  LBay_bdytide_rotT_M2_grid_V.nc
  LBay_bdytide_rotT_N2_grid_V.nc
  LBay_bdytide_rotT_S2_grid_V.nc
  LBay_bdytide_rotT_K1_grid_V.nc
  LBay_bdytide_rotT_K2_grid_V.nc
  LBay_bdytide_rotT_P1_grid_V.nc
  LBay_bdytide_rotT_O1_grid_V.nc
  LBay_bdytide_rotT_MF_grid_V.nc
  LBay_bdytide_rotT_Q1_grid_V.nc
  LBay_bdyT_y2000m01.nc
  LBay_bt_bdyT_y2000m01.nc
  LBay_bdyU_y2000m01.nc
  LBay_bdyV_y2000m01.nc

.. Warning:

   It doesn't quite work with ``ln_tra = .false.``

  This wont work because variable ``ft`` which deals with the number of time steps
  if only defined using the T fields, but needed for the velocity bcs. Wont work
  with ln_dyn3d = .true. either
  See e.g.::

    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 435, in process_bdy
     ft, num_bdy, time_counter, unit_origin)
     UnboundLocalError: local variable 'ft' referenced before assignment


Prepare the boundary files (need to fix some variable names)::

  cd $INPUTS

  module load nco/gcc/4.4.2.ncwa

  ncrename -v depthu,gdepu LBay_bdyU_y2000m01.nc
  ncrename -v depthv,gdepv LBay_bdyV_y2000m01.nc
  module unload nco

.. note:
    Previously the initial T and S fields also needed fixing. But these seem to
    be OK now as they already contain gdept:
      ncrename -v deptht,gdept initcd_votemper.nc
      ncrename -v deptht,gdept initcd_vosaline.nc

Copy the new files back onto ARCHER
::

  livljobs4$
  cd $INPUTS
  for file in LBay*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/. ; done
  scp coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/.

I previously also had the following. But now I comment it out as these files are
generated on ARCHER::

  #for file in initcd_vo*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/. ; done


8. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

*(16 OCt 17)*


Edit the output to have 1hrly SSH::

 vi file_def_nemo.xml
 ...
 <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
  <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
    <field field_ref="ssh"          name="zos"   />
  </file>
 </file_group>

---

Create a short queue runscript::

  vi runscript
  #!/bin/bash
  #PBS -N LBay
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M jelt@noc.ac.uk


  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  #  echo $(readlink -f $PBS_O_WORKDIR)
  # export OMP_NUM_THREADS=1

  cd $PBS_O_WORKDIR
  #
    echo " ";
    OCEANCORES=96
    XIOCORES=1
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core

  #aprun -n $OCEANCORES -N 24 ./opa
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa
  #aprun -b -n $XIOCORES -N 1 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

  exit

---

No met (missing slp) ``ln_usr=T``. rn_rdt=60s. Output more harmonics (20-30days).
Run for 30 days::

 cd $EXP
 qsub -q short runscript

**IT WORKS!**

* Tidy up this note. Freeze it. Add met?

 ---




Note about changes to the namelist_cfg
======================================

In order to make a simulation with the new v4 namelist_cfg,  a number of
changes appear with its format. Some of these are noted below. But since they are
a cut down from a debug log they may not be exhaustive.

I might have added some tidal harmonic constituent names to the template::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters (#ifdef key_tide)
  !-----------------------------------------------------------------------
  clname(1)    = 'Q1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(2)    = 'O1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(3)    = 'P1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(4)    = 'K1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(5)    = 'N2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(6)   =  'M2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(7)   = 'S2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(8)   = 'K2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(9)   = 'M4'   !  name of constituent - all tidal components must be set in namelist_cfg


Also noted that the sbc namelist variables have changed. Now use ``ln_blk`` and
 ``ln_COARE_3p5= .true.`` instead of ``ln_blk_core``. *(Except that I use*
 *ln_blk=F because I don't have a SLP variable for the COARE algorithm to work...)*

There was an extra variable in namdom. Comment out ldbletanh (which was essential in DOMAINcfg)::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
  ...
  !    ldbletanh   =    .false.             !  Use/do not use double tanf function for vertical coordinates



The configuration namelist variable is simplified with the new domain_cfg.nc file::

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
        !                   !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
        cn_domcfg = "domain_cfg"         ! domain configuration filename
        !
     ln_write_cfg= .false.   !  (=T) create the domain configuration file
        cn_domcfg_out = "domain_cfg_out" ! newly created domain configuration filename
        !
     ln_use_jattr = .false.  !  use (T) the file attribute: open_ocean_jstart, if present
     !                       !  in netcdf input files, as the start j-row for reading
  /


Equation of state needed fixing::

  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .false.         !  = Use TEOS-10 equation of state
     ln_eos80    = .true.         !  = Use EOS80 equation of state
     ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)


``&namdom`` has a new format (I think)::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
   ln_linssh   = .false.   !  =T  linear free surface  ==>>  model level are fixed in time
   nn_closea   =    0      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
   !
   nn_msh      =    0      !  create (>0) a mesh file or not (=0)
   rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice
   !
   rn_rdt      =  60.     !  time step for the dynamics (and tracer if nn_acc=0)
   rn_atfp     =    0.1    !  asselin time filter parameter
   !
   ln_crs      = .false.   !  Logical switch for coarsening module


Note ``ln_linssh`` seems to have moved to ``&namdom``. (I think it was in ``&namzgr_sco``
 before)

I added in s-coordinate option and which seemed to have vanished from the default
namelist_cfg template::

  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate                                  (default: NO selection)
  !-----------------------------------------------------------------------
     ln_zco      = .false.   !  z-coordinate - full    steps
     ln_zps      = .false.   !  z-coordinate - partial steps
     ln_sco      = .true.   !  s- or hybrid z-s-coordinate
     ln_isfcav   = .false.   !  ice shelf cavity
  /
  !-----------------------------------------------------------------------
  &namzgr_sco
  !-----------------------------------------------------------------------
     ln_s_sf12=.true.,
     ln_s_sh94=.false.,
     ln_sigcrit=.true.,
     rn_alpha=4.4,
     rn_bb=0.8,
     rn_efold=0.0,
     rn_hc=10000.0,
     rn_rmax=1.0,
     rn_sbot_max=10000.0,
     rn_sbot_min=1.0e-5,
     rn_theta=6.0,
     rn_thetb=1.0,
     rn_zb_a=0.024,
     rn_zb_b=-0.2,
     rn_zs=1.0,
  /


Though the fields are set up (except for SLP field), no met forcing is used::

  !-----------------------------------------------------------------------
  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
     ln_usr      = .true.    !  user defined formulation                  (T => check usrdef_sbc)
     ln_flx      = .false.   !  flux formulation                          (T => fill namsbc_flx )
     ln_blk      = .false.    !  Bulk formulation                          (T => fill namsbc_blk )

Though the fields are ready to go I have not used inital T,S conditions::

  !-----------------------------------------------------------------------
  &namtsd        !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  !              !  file name                 ! frequency (hours) ! variable ! time interp.!  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                            !  (if <0  months)  !   name   !  (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_tem  = 'initcd_votemper',         -1        ,'votemper' ,    .false.    , .true. , 'yearly'   , ''       ,   ''    ,    ''
     sn_sal  = 'initcd_vosaline',         -1        ,'vosaline' ,    .false.    , .true. , 'yearly'   , ''       ,   ''    ,    ''
     !
     cn_dir        = '../../../../INPUTS/'     !  root directory for the location of the runoff files
     ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)

No lateral diffusion in momentum::

  ln_dynldf_lap =  .false.    !    laplacian operator
  ln_dynldf_blp =  .false.    !  bilaplacian operator

Don't need a boundaries mask file. But the rimwidth must match that used in the
 PyNEMO process::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
    ...
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
    ...
    nn_rimwidth   = 9                    !  width of the relaxation zone

Quadratic bottom friction - This means changing to nonlinear friction, with a coeff 2.5E-3, log layer and
bottom roughness as follows
::
  !-----------------------------------------------------------------------
  &nambfr        !   bottom friction
  !-----------------------------------------------------------------------
     nn_bfr      =    2      !  type of bottom friction :   = 0 : free slip,  = 1 : linear friction
                             !                              = 2 : nonlinear friction
     rn_bfri2    =    2.5e-3 !  bottom drag coefficient (non linear case)
     rn_bfeb2    =    0.0e0  !  bottom turbulent kinetic energy background  (m2/s2)
     ln_loglayer =    .true. !  loglayer bottom friction (only effect when nn_bfr = 2)
     rn_bfrz0    =    0.003  !  bottom roughness (only effect when ln_loglayer = .true.)
  /

I also changed the rn_bfeb2 = 2.5e-3 to zero. This would depend on the tke scheme.
I am not sure what it should be...




Rebuild the output and inspect
++++++++++++++++++++++++++++++

Rebuild the SSH files using old tools::

  export WDIR=/work/n01/n01/jelt/LBay/
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 Lbay_1h_20000101_20000130_SSH 5
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 Lbay_1d_20000101_20000130_Tides 5

Should remove individual processor files once the build is verified::

  rm Lbay_1h_20000101_20000130_SSH_00??.nc
  rm Lbay_1d_20000101_20000130_Tides_00??.nc

Inspect locally e.g.::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/Lbay_1d_20000101_20000130_Tides.nc .
  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/Lbay_1h_20000101_20000130_SSH.nc .

  ferret
  use Lbay_1d_20000101_20000110_Tides.nc
  plot /i=25/j=70 SOSSHEIG

Use some python to inspection of the domain_cfg.nc file or ssh, Tide output. See::

  cd /Users/jeff/GitLab/NEMO-RELOC/docs/source
  $ipython
  >> run quickplotNEMO

---













---

.. note::

  **TO DO** another time / for Solent config

  * Add in met forcing
  * Turn on T,S at boundares
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
