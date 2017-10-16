============================================================
Setting up a Liverpool Bay NEMO configuration, inc. livljob4
============================================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/LBay_livljobs4.html
*Actually I've taken this off readthedocs as it is internal. However writing in*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Originally tried building everything on ARCHER but ran into a java problem with
NRCT (nee PyNEMO) that they needed to fix. Also MDP can't use ARCHER.

Then got it working on a combination of Liverpool machines for setup with
 simulations on ARCHER.

Then tried to update this recipe with the pre-release of NEMOv4, which in
particular has a new method for handling the domain configuration. For this I use
the ORCHESTRA realease of the trunk (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. ARCHER: Run simulation


Issues that arose
=================

* PyNEMO doesn't yet deal with sigma parent grids.
* PyNEMO needs a gdept function. DOMAINcfg (will) only store e3t.

Other related recipes (all NEMOv3.6):
Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Follow PyNEMO recipe for Lighthouse Reef: ``http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html``
Follow PyNEMO recipe for LBay on ARCHER (not complete because PyNEMO java issue): ``http://nemo-reloc.readthedocs.io/en/latest/LBay.html``.
(Perhaps this other note is superceded by this present one, though the other
note includes information of using PyNEMO with different data sources).

.. note:

  It is very easy to break the code with bad edits to the iodef.xml file. Don't
  change the iodef.xml file at the same time as something else.

.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO bits are done on livljobs4.
(There was a problem with some java bits working on ARCHER)
Starting on ARCHER::

  ssh login.archer.ac.uk

  export CONFIG=LBay
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

  export JINPUTS=/work/n01/n01/jdha/2017/INPUTS/ODA/E-AFRICA
  export JEXP=/work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/EXP00/
  export OCDIR=/work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/


.. note:
 I will remove these links to James' files when I've figured out how to build my own

---

Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

Or just link XIOS executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---

Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

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

---

*(16 Oct 17)*

No gdept output which is needed for PyNEMO. Try outputting mesh files.
 Edit namelist_cfg::

  ln_e3_dep   = .false.    ! =T : e3=dk[depth] in discret sens.


Resubmit job

---



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
       <ns0:variable name="e3t" orgName="e3t_0" />
       <ns0:variable name="e3w" orgName="e3w_0" />
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
  ncrename -v deptht,gdept initcd_votemper.nc
  ncrename -v deptht,gdept initcd_vosaline.nc
  module unload nco

Copy the new files back onto ARCHER
::

  livljobs4$
  cd $INPUTS
  for file in LBay*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/. ; done
  for file in initcd_vo*nc; do scp $file jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/. ; done
  scp coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/.

8. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

*(21 Sept 2017 / 6 Oct 17)*

Open a terminal on **ARCHER**. Redefine PATHS. Reload modules::

  export CONFIG=LBay
  export WORK=/work/n01/n01
  export WDIR=$WORK/$USER/$CONFIG
  export INPUTS=$WORK/$USER/$CONFIG/INPUTS
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

OPA and XIOS are already built.


Link the boundary data to the EXP direcory::

 cd $EXP
 ln -s $INPUTS/coordinates.bdy.nc       $EXP/coordinates.bdy.nc
 ln -s $INPUTS/LBay_bdyT_y2000m01.nc    $EXP/LBay_bdyT_y2000m01.nc
 ln -s $INPUTS/LBay_bdyU_y2000m01.nc    $EXP/LBay_bdyU_y2000m01.nc
 ln -s $INPUTS/LBay_bdyV_y2000m01.nc    $EXP/LBay_bdyV_y2000m01.nc
 ln -s $INPUTS/LBay_bt_bdyT_y2000m01.nc $EXP/LBay_bt_bdyT_y2000m01.nc
 ln -s $INPUTS                          $EXP/bdydta

.. old:  and update the namelist_cfg for running, not mesh generation
 #sed -e 's/nn_msh      =    3/nn_msh      =    0/' namelist_cfg > tmp
 #sed -e 's/nn_itend    =      1/nn_itend    =       1440 /' tmp > namelist_cfg

There was a problem running with the namelist_cfg fresh from DOMAINcfg. First with
with namcfg group. So I cut it back to just reading the cfg file::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
     cn_domcfg = "domain_cfg"         ! domain configuration filename


Then there was an extra variable in namdom. Comment out ldbletanh (which was essential in DOMAINcfg)::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
  ...
  !    ldbletanh   =    .false.             !  Use/do not use double tanf function for vertical coordinates


Also noted that the sbc namelist variables have changed. Now use ``ln_blk`` and ``ln_COARE_3p5= .true.``
instead of ``ln_blk_core``

Also need to make sure the harmonic tidal boundary files are consistent with the
 harmonics expected e.g.::

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


Edit the output to have 1hrly SSH::

 vi file_def_nemo.xml
 ...
 <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
  <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
    <field field_ref="ssh"          name="zos"   />
  </file>
 </file_group>



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

Then submit::

 cd $EXP
 qsub -q short runscript

 4806706.sdb

---


**ERROR**::

  > Error [CAttributeMap::operator[](const StdString& key)] : In file '/work/n01/n01/jelt/xios-2.0_r1080/src/attribute_m
  ap.cpp', line 56 -> [ key = time_origin] key not found !
  > Error [CAttributeMap::operator[](const StdString& key)] : In file '/work/n01/n01/jelt/xios-2.0_r1080/src/attribute_m
  ap.cpp', line 56 -> [ key = time_origin] key not found !
  > Error [CAttributeMap::operator[](const StdString& key)] : In file '/work/n01/n01/jelt/xios-2.0_r1080/src/attribute_m
  ap.cpp', line 56 -> [ key = time_origin] key not found !
---

Looks like a XML problem. Copy working XML files from EAfrica
::
  export CONFIG=LBay
  export WORK=/work/n01/n01
  export WDIR=$WORK/$USER/$CONFIG
  export INPUTS=$WORK/$USER/$CONFIG/INPUTS
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00


  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

  export JINPUTS=/work/n01/n01/jdha/2017/INPUTS/ODA/E-AFRICA
  export JEXP=/work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/EXP00/

---

Copy working XML files from EAfrica::

  cd $EXP
  mv *.xml XML/.
  cp /work/n01/n01/jelt/ACCORD/trunk_NEMOGCM_r8395/CONFIG/ACCORD/EXP_EAFRICA/*xml .

This works. Highlights missing EOS choice in namelist_cfg. Add in::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
  ln_teos10   = .false.         !  = Use TEOS-10 equation of state
  ln_eos80    = .true.         !  = Use EOS80 equation of state
  ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)
                               !
  !                     ! S-EOS coefficients (ln_seos=T):
  !                             !  rd(T,S,Z)*rau0 = -a0*(1+.5*lambda*dT+mu*Z+nu*dS)*dT+b0*dS
  rn_a0       =  1.6550e-1      !  thermal expension coefficient
  rn_b0       =  7.6554e-1      !  saline  expension coefficient
  rn_lambda1  =  5.9520e-2      !  cabbeling coeff in T^2  (=0 for linear eos)
  rn_lambda2  =  7.4914e-4      !  cabbeling coeff in S^2  (=0 for linear eos)
  rn_mu1      =  1.4970e-4      !  thermobaric coeff. in T (=0 for linear eos)
  rn_mu2      =  1.1090e-5      !  thermobaric coeff. in S (=0 for linear eos)
  rn_nu       =  2.4341e-3      !  cabbeling coeff in T*S  (=0 for linear eos)

Odd conflict in notation between namelist_ref in my new build and in JAmes'
Copy James' here::

  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/EXP00/namelist_ref $EXP/namelist_ref


  vi namelist_ref
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .false.         !  = Use TEOS-10 equation of state
     ln_eos80    = .false.         !  = Use EOS80 equation of state
     ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)

(Other other format had integers to choose the scheme.)

Problem in namdom
::

  The IOIPSL calendar is "gregorian", i.e. leap year

  ===>>> : E R R O R
       ===========

  misspelled variable in namelist namdom in configuration namelist iostat =    19

  Namelist namdom : space & time domain
     linear free surface (=T)              ln_linssh  =  F
     suppression of closed seas (=0)       nn_closea  =            0
     create mesh/mask file(s)              nn_msh     =            0
          = 0   no file created
          = 1   mesh_mask
          = 2   mesh and mask
          = 3   mesh_hgr, msh_zgr and mask
     treshold to open the isf cavity       rn_isfhmin =
  1.00000000000000       (m)
     ocean time step                       rn_rdt     =
  60.0000000000000
     asselin time filter parameter         rn_atfp    =
  0.100000000000000
     online coarsening of dynamical fields ln_crs     =  F


Try the new format::

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



The vertical coordintes choice thing seems to have disappeared in the new build::

  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zps      = .false.   !  z-coordinate - partial steps   (T/F)
     ln_sco      = .true.    !  s- or hybrid z-s-coordinate    (T/F)
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
     ln_s_sh94   = .false.   !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
     ln_s_sf12   = .true.    !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
     ln_sigcrit  = .true.    !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                             !  stretching coefficients for all functions
     rn_hc       =   50.0    !  critical depth for transition to stretched coordinates

There are a few references to s-coordinates in the HPG and later diffusion variables.

Hmm. Try the other way around start with James' namelist_cfg_R12 and change bits to
match my LBay (old working) example

Update boundary mask file name::

  /
  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .true.              !  =T : read mask from file
      cn_mask_file   = 'bdydta/LBay_bdyT_y2000m01.nc'     !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !

Automatically set the processor decomposition (not sure if it is used)::

  &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)


Switch the vertical grid thing off. Comment it out as the default is .true. BTW
the run really doesn't work without this action::

  vi namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
  !   ln_e3_dep   = .false.    ! This will be obsolete soon. See namelist_ref

Change some lateral diffusion settings::

  !-----------------------------------------------------------------------
  &namtra_ldf    !   lateral diffusion scheme for tracers                 (default: NO diffusion)
  !-----------------------------------------------------------------------
     !                       !  Operator type:
     !                           !  no diffusion: set ln_traldf_lap=..._blp=F
     ln_traldf_lap   =  .true.  !    laplacian operator
     ln_traldf_blp   =  .false.  !  bilaplacian operator
     !
     !                       !  Direction of action:
     ln_traldf_lev   =  .false.  !  iso-level
     ln_traldf_hor   =  .false.  !  horizontal (geopotential)
     ln_traldf_iso   =  .true.  !  iso-neutral (standard operator)

Gets further. Now the ocean.output ends with::

  dia_25h_init : Output 25 hour mean diagnostics
  ~~~~~~~~~~~~
  Namelist nam_dia25h : set 25h outputs
  Switch for 25h diagnostics (T) or not (F)  ln_dia25h  =  F

  AAAAAAAA


  sbc_tide : Update of the components and (re)Init. the potential at kt=
            1
  ~~~~~~~~
  Q1    -12.3894431662406       0.908086877990601      -0.702799902800797
   6.495854101908828E-005
  O1    -12.3894431662406       0.908086877990601        1.76695796869312
   6.759774402887834E-005
  P1    0.000000000000000E+000   1.00000000000000      -0.189080734230733
   7.252294578606445E-005
  S1    0.000000000000000E+000   1.00000000000000        3.14377431515479
   7.272205216643040E-005
  K1   -0.138134509178184       0.943678543708499        6.47662936454030
   7.292115854679635E-005
  2N2   -12.5967638158724        1.02169282172100        3.30407159024559
   1.352404965560946E-004
  MU2   -12.5967638158724        1.02169282172100        10.1996260361573
   1.355937008184885E-004
  N2    -12.5967638158724        1.02169282172100        5.77382946173951
   1.378796995658846E-004
  NU2   -12.5967638158724        1.02169282172100        12.6693839076512
   1.382329038282786E-004
  M2    -12.5967638158724        1.02169282172100        8.24358733323342
   1.405189025756747E-004
  L2    -12.5967638158724        1.19824874449528        13.8549378583171
   1.431581055854647E-004
  T2    0.000000000000000E+000   1.00000000000000        6.32212956235245
   1.452450074605893E-004
  S2    0.000000000000000E+000   1.00000000000000        6.28754863030957
   1.454441043328608E-004
  K2   -0.264695210962853       0.854177079157964        16.0948513826704
   1.458423170935927E-004
  M4    -25.1935276317449        1.04385622195621        16.4871746664668
   2.810378051513493E-004

Try setting tides to false::

  !-----------------------------------------------------------------------
&nam_tide      !   tide parameters
!-----------------------------------------------------------------------
   ln_tide     = .false.
   ln_tide_pot = .false.    !  use tidal potential forcing

Caused problems with Flather bc etc.
Turned ln_tide = .true., keep tidal potential off. Simulation terminates with
same output (above), having listed the harmonic components. Hmmm

Tide forcing directory::

  !-----------------------------------------------------------------------
  &nambdy_tide     ! tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/LBay_bdytide_rotT_'         !  file name root of tidal forcing files

Turned ln_tide_pot = .true. (I think that the tidal boundary files are way more
likely to give trouble than tidal potential forcing)

Added some more constituents to the TPXO forcing list. Regenerated with PyNEMO.
Not sure about which constituents to include under key_tide in the OPA namelist_cfg
Not sure because some of the constituent names differ. Is it looking for TPXO files
of this name, or is it setting the internal harmonic frequencies. On reflection I guess
the harmonic analyisis is entirely separate and will only pick out freq requested
in the harm namelist.

Submit::

  qsub runscript  # changed to a 4 minute walltime request
  4834214.sdb

**SAME ERROR / NON ERROR AS ABOVE. What next. How to get past this point....?**

*Fixed lots of stuff. Mostly putting in missing stuff into namelist_cfg*

---

*10 Oct 2017*

Fail
::

 sbc_tide : Update of the components and (re)Init. the potential at kt=
           1
 ~~~~~~~~
 Q1    -12.3894431662406       0.908086877990601      -0.702799902800797
  6.495854101908828E-005
 O1    -12.3894431662406       0.908086877990601        1.76695796869312
  6.759774402887834E-005
 P1    0.000000000000000E+000   1.00000000000000      -0.189080734230733
  7.252294578606445E-005
 K1   -0.138134509178184       0.943678543708499        6.47662936454030
  7.292115854679635E-005
 N2    -12.5967638158724        1.02169282172100        5.77382946173951
  1.378796995658846E-004
 M2    -12.5967638158724        1.02169282172100        8.24358733323342
  1.405189025756747E-004
 S2    0.000000000000000E+000   1.00000000000000        6.28754863030957
  1.454441043328608E-004
 K2   -0.264695210962853       0.854177079157964        16.0948513826704
  1.458423170935927E-004
 M4    -25.1935276317449        1.04385622195621        16.4871746664668
  2.810378051513493E-004
                     iom_nf90_open ~~~ open existing file: ../../../../INPUTS/cu
 tdown_drowned_u10_DFS5.1.1_y2000.nc in READ mode
                    ---> ../../../../INPUTS/cutdown_drowned_u10_DFS5.1.1_y2000.n
 c OK
                     iom_nf90_open ~~~ open existing file: ../../../../INPUTS/cu
 tdown_drowned_u10_DFS5.1.1_y2000.nc in READ mode
                    ---> ../../../../INPUTS/cutdown_drowned_u10_DFS5.1.1_y2000.n
 c OK
                     iom_close ~~~ close file: ../../../../INPUTS/cutdown_drowne
 d_u10_DFS5.1.1_y2000.nc ok
                     iom_nf90_open ~~~ open existing file: ../../../../INPUTS/we
 ights_bicubic_atmos.nc in READ mode
                    ---> ../../../../INPUTS/weights_bicubic_atmos.nc OK
           read src01 (rec:      1) in ../../../../INPUTS/weights_bicubic_atmos.nc ok



Inspection of v3.6 ocean.output suggests there is a problem with
 weight_bicublic_atmos.nc. The output would have continued as::

  ls /work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LBay/EXP00/ocean.output
  ...

  iom_nf90_open ~~~ open existing file: ../../../../../INPUTS
  /cutdown_drowned_u10_DFS5.1.1_y2000.nc in READ mode
  ---> ../../../../../INPUTS/cutdown_drowned_u10_DFS5.1.1_y200
  0.nc OK
  iom_close ~~~ close file: ../../../../../INPUTS/cutdown_dro
  wned_u10_DFS5.1.1_y2000.nc ok
  iom_nf90_open ~~~ open existing file: ../../../../../INPUTS
  /weights_bicubic_atmos.nc in READ mode
  ---> ../../../../../INPUTS/weights_bicubic_atmos.nc OK
  read src01 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read src02 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read src03 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read src04 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read wgt01 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read wgt02 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok
  read wgt03 (rec:      1) in ../../../../../INPUTS/weights_bicubic_atmos.nc ok



Try switch to CORE v3.0 ln_COARE_3p0
No joy. Try and check what comes next. Is this reading in the right file or should it be looking for tide data?

.. note: There is a bug with the namelist implementation for COARE forcing.
ln_COARE_3p0= .true.   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
ln_COARE_3p5= .false.   ! "COARE 3.5" algorithm   (Edson et al. 2013)
--> both as true in ocean.output

namelist_cfg:   ln_COARE_3p0= .false.   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
namelist_cfg:   ln_COARE_3p5= .true.   ! "COARE 3.5" algorithm   (Edson et al. 2013)
--> both false in ocean.output

Typo in sbcblk.F90 line 251/ Logical flag pointing to wrong variable. See::

  WRITE(numout,*) '      "COARE 3.5" algorithm   (Edson et al. 2013)         ln_COARE_3p5 = ', ln_COARE_3p0


Try ln_NCAR instead...
*(10 Oct 2017)*

Change: nn_dyn2d_data = 2 —> 3. This just means that the 'LBay_bt_bdyT etc in
&nambdy_dta are read in::

  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      ...
      nn_dyn2d_dta   =  3                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing


**Does it work?**

Waiting for standard queue to proces, I note that

There are a number of differences in the tidal boundary conditions between James’ ORCHESTRA run and my LBay simulations.

Check nambdy::

  /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/EXP00/namelist_cfg
  /work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/namelist_cfg

**PENDING**

Need to keep a track of differences between namelist_cfg in EXP and DOMAINcfg


---

*(11 Oct 2017)*

Update the mask file::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'domain_cfg.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)


Make some changes in the open boundary files. Originally::

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'Tbdy',                   -1        , 'sossheig',    .false.   , .true. ,  'yearly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'Ubdy',                   -1        , 'vobtcrtx',    .false.   , .true. ,  'yearly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'Vbdy',                   -1        , 'vobtcrty',    .false.   , .true. ,  'yearly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'amm12_bdyU_u3d',         24        , 'vozocrtx',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'amm12_bdyV_u3d',         24        , 'vomecrty',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'amm12_bdyT_tra',         24        , 'votemper',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'amm12_bdyT_tra',         24        , 'vosaline',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''


Change to::

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
  bn_ssh      = 'LBay_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_u2d      = 'LBay_bdyU',  24        , 'vobtcrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_v2d      = 'LBay_bdyV',  24        , 'vobtcrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_u3d      = 'LBay_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_v3d      = 'LBay_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_tem      = 'LBay_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  bn_sal      = 'LBay_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''

Though I need to check that these contain the correct variables. In particular
 the 2D and 3D currents are in the same file?

Though I don’t use the 3D data - other than setting it to the initial state.


Observation.
Boundary 2D tides looks suspiciously like something I want to switch on. Currently::

!-----------------------------------------------------------------------
&nambdy_tide   !  tidal forcing at open boundaries
!-----------------------------------------------------------------------
   filtide      = 'bdydta/LBay_bdytide_rotT_'         !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
   ln_bdytide_conj  = .false.                   !


Things now look like they fail with the bulk forcing. Only read in one variable
from weights_bicubic_atmos.nc whereaas in the old code I read in ten or so...

::

  tail ocean.output
                       iom_nf90_open ~~~ open existing file: ../../../../INPUTS/we
  ights_bicubic_atmos.nc in READ mode
                   ---> ../../../../INPUTS/weights_bicubic_atmos.nc OK
          read src01 (rec:      1) in ../../../../INPUTS/weights_bicubic_atmos.nc ok

---

.. note:

  Not sure I have the enthusiasm to debug this COARE implentation. Try copying old
   SBC code from Maria and hope is compiles::

    cd /work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/MY_SRC
    cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/NEMO/OPA_SRC/SBC/sbcmod.F90 .
    cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/NEMO/OPA_SRC/SBC/sbcblk_core.F90 .

  Will need to make changes in namelist
  &namsbc_blk   --> &namsbc_core

  !ln_NCAR     = .false.   ! "NCAR"      algorithm   (Large and Yeager 2008)
  !ln_COARE_3p0= .true.   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
  !ln_COARE_3p5= .false.   ! "COARE 3.5" algorithm   (Edson et al. 2013)
  !ln_ECMWF    = .false.   ! "ECMWF"     algorithm   (IFS cycle 31)

  &namsbc
  ln_blk -->    ln_blk_core = .true.

  Compilation expects a bdy_par file. Bah!...


Instead copy the usrdef_sbc.F90 file to impose zero surface forcing::
  rm $EXP/../MY_SRC/*
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/MY_SRC/usrdef_sbc.F90 $EXP/../MY_SRC/.

Rebuild works. Edit namelist::

  vi namelist_cfg
  &namsbc
    ln_usr      = .true.
    ln_blk      = .false.    !  Bulk formulation                          (T => fill namsbc_blk )



  &nambdy
    ln_bdy         = .false.
    nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
Resubmit
----

Trying turning on tidal forcing at boundaries in namelist (Though James had this set false) ::

  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/LBay_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .true.                   !

This seems to do nothing. Keep it **FALSE**

Try::
   !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'domain_cfg.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state

Error. Seg fault.::

  tail ocean.output

  M2    -12.5967638158724        1.02169282172100        8.24358733323342
   1.405189025756747E-004
  S2    0.000000000000000E+000   1.00000000000000        6.28754863030957
   1.454441043328608E-004
  K2   -0.264695210962853       0.854177079157964        16.0948513826704
   1.458423170935927E-004
  M4    -25.1935276317449        1.04385622195621        16.4871746664668
   2.810378051513493E-004
   usr_sbc : WAD_TEST_CASES case: NO surface forcing
   ~~~~~~~~~~~   utau = vtau = taum = wndm = qns = qsr = emp = sfx = 0

Try::

 !-----------------------------------------------------------------------
 &nambdy        !  unstructured open boundaries
 !-----------------------------------------------------------------------
     ln_bdy         = .true.              !  Use unstructured open boundaries
     nb_bdy         = 1                    !  number of open boundary sets
     ln_coords_file = .true.               !  =T : read bdy coordinates from file
     cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
     ln_mask_file   = .false.              !  =T : read mask from file
     cn_mask_file   = 'domain_cfg.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
     cn_dyn2d       = 'flather'               !
     nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                           !  = 1, bdy data are read in 'bdydata   .nc' files
                                           !  = 2, use tidal harmonic forcing data from files


Same Error. Seg fault. Though ran for 40s::

  tail ocean.output

  M2    -12.5967638158724        1.02169282172100        8.24358733323342
 1.405189025756747E-004
  S2    0.000000000000000E+000   1.00000000000000        6.28754863030957
   1.454441043328608E-004
  K2   -0.264695210962853       0.854177079157964        16.0948513826704
   1.458423170935927E-004
  M4    -25.1935276317449        1.04385622195621        16.4871746664668
   2.810378051513493E-004
   usr_sbc : WAD_TEST_CASES case: NO surface forcing
   ~~~~~~~~~~~   utau = vtau = taum = wndm = qns = qsr = emp = sfx = 0


Try reducing the timestep for 60s to 10s::

  &namdom
    rn_rdt      =  10.     !  time step for the dynamics (and tracer if nn_acc=0)

Same problem. The timestep does not seem to have any effect.
Looking at an abort file that was generated when ``ln_bdy=.false.`` and the
velocities blew up showed that it was apparent that a coastal point was blowing up in velocity.
ocean.output abort message (ln_bdy = false)::

  ...
  1.458423170935927E-004
  M4    -25.1935276193127        1.04385624018260        16.4795866457279
  2.810378051513493E-004
  usr_sbc : WAD_TEST_CASES case: NO surface forcing
  ~~~~~~~~~~~   utau = vtau = taum = wndm = qns = qsr = emp = sfx = 0
           nit000-1 surface forcing fields set to nit000

  zdf_evd : Enhanced Vertical Diffusion (evd)
  ~~~~~~~


  zdf_mxl : mixed layer depth
  ~~~~~~~

  ssh_nxt : after sea surface height
  ~~~~~~~

  div_hor : horizontal velocity divergence
  ...




This test highlights the next output step that is missing when ln_bdy=T:
``           nit000-1 surface forcing fields set to nit000``. This must be in
an sbc module (``sbcmod.F90``). It is not clear why ln_bdy would break that...


* James has a bdy_mask.nc file. This looks like top_level in domain_cfg.nc

I can make a bdy_mask.nc file (I am not confident which modules to load to get
working on ARCHER so I'll do it on livljobs4)::

  livljobs$
  cd Desktop
  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/domain_cfg.nc  .

  module load nco/gcc/4.4.2.ncwa
  ncks -v top_level,nav_lat,nav_lon  domain_cfg.nc tmp.nc
  ncrename -v top_level,bdy_msk tmp.nc
  ncwa -a t tmp.nc bdy_mask.nc

  livmap$
  scp jelt@livljobs4:Desktop/bdy_mask.nc ~/Desktop/.
  ferret
  use bdy_mask.nc
  shade bdy_msk

However I notice that the (wet) western boundary is masked out by this mask. Is
that true in James' mask? Yes it appears so.

Copy new file into EXP directory (having copied it into $INPUTS)::

 cp $INPUTS/bdy_mask.nc $EXP/.

Then edit the namelist_cfg.nc to include this new file::

 !-----------------------------------------------------------------------
 &nambdy        !  unstructured open boundaries
 !-----------------------------------------------------------------------
     ln_bdy         = .true.              !  Use unstructured open boundaries
     nb_bdy         = 1                    !  number of open boundary sets
     ln_coords_file = .true.               !  =T : read bdy coordinates from file
     cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
     ln_mask_file   = .true.              !  =T : read mask from file
     cn_mask_file   = 'bdy_mask.nc'

This still dies in the same place::

  tail ocean.output
  ...
  M4    -25.1935276193127        1.04385624018260        16.4795866457279
   2.810378051513493E-004
   usr_sbc : WAD_TEST_CASES case: NO surface forcing
   ~~~~~~~~~~~   utau = vtau = taum = wndm = qns = qsr = emp = sfx = 0


Inspection of ``sbcmod.F90`` shows that thing might be different if it is a restart.
Copy a restart from old code base. Update namelist_cfg for a restart and resubmit

::
   &namrun
   nn_date0    =  20000106   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
   ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
   cn_ocerst_out   = "restart_oce_out"   !  suffix of ocean restart name (output)
   nrstdt = 0

It turns out that the restart file is missing e3t and has a wrong variable name for time.
Perhaps easier to create a new restart from a single timestep...?

no even with one timestep there is a seg fault problem.

---

change initialisation to false for T and S
ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)


Test. Turn off tides::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .false.
     ln_tide_pot = .false.    !  use tidal potential forcing

  &nambdy
       nn_dyn2d_dta   =  0  ! was 2

Didn't change the simulation crash point.

Change rimwidth. It is different in the boundary files to namelist_cfg::

  &nambdy
  nn_rimwidth   = 9                    !  width of the relaxation zone

(though the bdy fields, I think, are turned off).

Oh the rimwidth variable matching the boundary files seem to fix the problem. Next turn things back on.
::

  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.

  &nambdy        !  unstructured open boundaries
     nn_dyn2d_dta   =  1                   !  = 0, bdy data are equal to the initial state

This runs but blows up in momentum. Try nn_dyn2d_dta = 2 (like James)
Still blows up.

Try tiny time step  rn_rdt=5. No.

Try nn_dyn_dta = 3. Stabilize with boundary velocities? No

Update bottom friction to match old LBay run (rather than James' E-Africa run).
THis means changing to nonlinear friction, with a coeff 2.5E-3, log layer and
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
Though I don't have an idea of what it should be.

Turn off time splitting: ln_bt_nn_auto=.false.,
I think this has moved to a new namelist group &namdyn_spg
Instead turn it off there: ln_bt_auto    = .false.

This makes the run go further. 9 time steps...
Try and manually fix barotropic time step scale. nn_baro=1 (not 30)
This blew up in 5 steps. Revert to nn_baro=30

Switch bilaplacian to laplacian to see if it suddenly works

It works better - 19 steps.

Output every 3 steps and inspect. Didn't output.

Try and use ramp_tide to ramp up tides over 1 day.
(checked that rdttideramp corresponds to the number of ramping days.)

**It worked!!** 1 day in 40s. Extend to 5 days to run in 5 mins

Switch from Laplacian to Bilapacian diffusion for momentum.

No motion
try n_dyn2d_dta   =  2

Take that boundary mask off: ln_mask_file   = .false.
Didn't help.


Perhaps it is an issue with the XML files. Define a shortcut::

   EEXP=/work/n01/n01/jelt/ACCORD/trunk_NEMOGCM_r8395/CONFIG/ACCORD/EXP_EAFRICA

Hmm couldn't spot any differences wih the E-Africa XML files - they are identical

Try the tideramp = 0.5

No joy.
Try a restart.
Bad filename. Check and try again.

Take off the tideramp and it blows up (this is probably not to do with the restart)

There is a variable nb_harmo=0 in the ocean.output. It is not being defined properly.

Tried and changed the frequency of the boundary data: e..g bn_ssh      = 'LBay_bt_bdyT', -1   (was 24)

This blows up again. On the first time step.

It looks like the ramp is too severe. If it blows up in 19 steps = 19mins, with
a 1 day ramp. Maybe a 10 ramp would be better.

Do a cold start. With a 10 sim. (ramp can not exceed simulation time)
Signal dies out within a day. (Min SSS--> 100)

Switch off both laplacian and bilaplacian diffusion. Output info every hr::

  ln_dynldf_lap =  .false.    !    laplacian operator
  ln_dynldf_blp =  .false.    !  bilaplacian operator
  nn_write    =   24

Dies fast. Misinterpreted how the ramp worked. Never had a non-zero simulation with
the ramp. Make it smaller not larger: ``   rdttideramp =    1.     ``

SSH blows up::

  ==>> time-step=            1  ssh max:  1.311225785606566E-004
  ==>> time-step=           61  abs(U) max:   9.365770301892096E-002
  ==>> time-step=           61  SSS min:   36.8418750679304
  ==>> time-step=           61  ssh max:  6.350050976376775E-002
  ==>> time-step=          121  abs(U) max:   0.120795267289029
  ==>> time-step=          121  SSS min:   36.8414750346018
  ==>> time-step=          121  ssh max:  0.113158707022428

  ===>>> : E R R O R
         ===========

  stp_ctl : the ssh is larger than 10m
  =======
  kt=   154 max ssh:    Infinity, i j:    32   90


Reduce harmonic to just M2. Blows up in SSH in the same way.

Reduce the timestep rn_rdt=6.::

  ==>> time-step=         1321  ssh max:  0.111035814265442
  ==>> time-step=         1381  abs(U) max:   0.117910378601910
  ==>> time-step=         1381  SSS min:   36.8413812872391
  ==>> time-step=         1381  ssh max:  0.112485247560570
  ==>> time-step=         1441  abs(U) max:   0.116803563870784
  ==>> time-step=         1441  SSS min:   36.8413491215503
  ==>> time-step=         1441  ssh max:  0.113826565377745

  ===>>> : E R R O R
         ===========

  stpctl: the zonal velocity is larger than 20 m/s
  ======
  kt=  1489 max abs(U):   22.29    , i j k:    32  104   21

Ran for 2.5 hours instead of about 2hrs.

Need to look at some fields. Output initial conditions and evolving fields
Check CFL. Grid is about 1.3km x 0.7km
Therefore a 60s timestep 700m/60s > 10m/s >> u. Is OK
cp = sqrt(g.h) = sqrt(10*50) = 22 m/s.

Inspection of the domain_cfg.nc file shows that the e1x and e2x variables are wrong.
See: cd /Users/jeff/GitLab/NEMO-RELOC/docs/source
$ipython
>> run quickplotNEMO


**PENDING**

.. Ideas:
 * Try and restore things I changed --
 * Check tides are in
 * initial conditions in T and S
 * remove the bdy mask
 * surface forcing
 * rn_rdt=60
 * rn_hc=10000.0, --> 50






Rebuild the output and inspect
++++++++++++++++++++++++++++++

Rebuild the SSH files using old tools::

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

**Is there a semi-diurnal SSH signal?**

---













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
