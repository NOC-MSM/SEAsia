
In the following are specific instructions for a 3D baroclinic simulation of
NEMO with ERSEM (via the FABM coupler). The regional configuration is
demonstrated for South East Asia (75E to 135E and -20N to +25N). The model grid
has 1/12&deg; lat-lon resolution and 75 hybrid sigma-z-partial-step vertical levels.

Machines used: CentOS7 linux box, Cray XC30 HPC (ARCHER)

===========================================
Preparing for a new NEMO vp4 configuration
===========================================

Having cloned the NEMO-RELOC repository, the  entire build process can be run
with a single script ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh

However the process is taken through step by step here.


a. Clone the NEMO-RELOC repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step is set up a new configuration directory by git cloning the
repository (e.g. to your workspace)::

  cd /work/n01/n01/$USER
  git clone https://github.com/NOC-MSM/NEMO-RELOC


b. Set script paths, make paths and directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure that the first two lines  in ``SCRIPTS/make_paths.sh`` are appropriately
defined:

.. literalinclude:: ../SCRIPTS/make_paths.sh
  :start-at: export CONFIG
  :end-at: export WORK

Then ``SCRIPTS/make_paths.sh`` and ``SCRIPTS/make_directories.sh`` will create
the expected directory structure.

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: echo "Making Paths"
  :end-at: make_directories.sh

c. Build XIOS
^^^^^^^^^^^^^

Then we build XIOS for managing the input/output on distributed architecture.
Note that the modules need to be loaded appropriately for the architecture. The
settings here are for ARCHER:

.. literalinclude:: ../SCRIPTS/make_xios.sh

The script builds two versions of XIOS. XIOS2.5 is used for the NEMO simulations
but an earlier version XIOS1 (where performance is not an issue) is used for
some of the tooling required to build the configuration.

This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Installing XIOS
  :end-at: make_xios


d. Build NEMO
^^^^^^^^^^^^^

Then build the NEMO executable. This is the dynamical ocean core the NEMO
framework. Technically OPA is the ocean dynamical core and NEMO is the framework
it sits within (including, for example, modules for ice, biogeochemistry etc).

In this documentation there are examples for a NEMO physics only build (called
``SCRIPTS/make_nemo.sh``) and a NEMO-FABM-ERSEM build
(called ``SCRIPTS/make_nemo_fabm_ersem.sh``). Here we walk through the physics
only build:

.. literalinclude:: ../SCRIPTS/make_nemo.sh


This is executed by ``SCRIPTS/main_setup.sh`` (commenting out the appropriate
make NEMO script):

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Installing NEMO-FABM-ERSEM
  :end-at: ./make_nemo.sh



e. Build TOOLS
^^^^^^^^^^^^^^

The next stage is to build tools that will be used for a number of processes
during the build. Specifically tools to assist with domain configuration file
generation (NESTING, which can be used to interpolate existing grids, and
DOMAINcfg, which builds the domain configuration netCDF file). At the same time
we build tools to reconstruct NEMO output (REBUILD_NEMO) and tools for
generating grid appropriate weights files to use with external surface forcing
data (WEIGHTS):

.. literalinclude:: ../SCRIPTS/make_tools.sh

This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
 :start-at: various grid tools
 :end-at: ./make_tools.sh


======================================================
Build a domain configuration file for SE Asia NEMO vp4
======================================================

To build the netCDF domain configuration file (``domain_cfg.nc``) a
``coordinates.nc`` and a ``bathymetry.nc`` file must be first generated.

Here we show:

  * how to make a child configuration from an existing NEMO run.
  * how to generate a new regional configuration from GEBCO (global) bathymetry.


a. Generate coordinates file from parent configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
Here we use the tool ``agrif_create_coordinates.exe``, previously built, which
reads cutting indices, resolution scale factor on parent, and grid location from
 ``namelist.input`` and outputs a new files with new resolution grid elements.

The cutting indices are can be obtained by inspection of the parent coordinates
and should be chosen to avoid problems arising from forcing the child with the
new proposed boundaries. Potential issues include having tidal amphidromes near
the boundary (where slight errors in the boundary tidal forcing can magnify
within the domain), isolated sea mounts or bathymetry features that might be
present in both parent and child bathymetries. Thin sections of wetpoints (3 or
fewer boxes wide) that are trapped between the boundary and interior land points.

For the SE Asia domain we select a region spanning from the middle of the Indian
subcontinent in the west to West Papua in the east. The region spans from
Taiwan in the north to Australia in the south. Intersections with the land are
chosen to be near-orthogonal. (75E to 135E and -20N to +25N).

.. literalinclude:: ../SCRIPTS/make_coordinates_from_parent.sh

Editting the ``START_FILES/DOMAIN/namelist.input`` for the cutting coordinates
and setting the resolution scale factor to 1 to preserve the parent resolution:

.. literalinclude:: ../START_FILES/DOMAIN/namelist.input
  :start-at: imin = 50
  :end-at: rhot = 1


This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Creating coordinate file
  :end-at: make_coordinates_from_parent



bi. Generate bathymetry file from GEBCO bathymetry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download some `GEBCO <https://www.gebco.net/data_and_products/gridded_bathymetry_data/>`_
bathymetry data. Here for the region (75E,-21N,134E,25N). Save locally. the
script ``SCRIPTS/make_bathymetry_from_gebco.sh`` assumes the downloaded file is
called ``GRIDONE_2008_2D_74.0_-21.0_134.0_25.0.nc``


For a domain as large as SE Asia, the 30-minute GEBCO
data might be too large to process and additionally needs some spatial filtering.
BODC also host a 1-minute data set (2008) which is easier to handle.

The workflow with GEBCO data involves pre-processing the raw bathymetry data to 1) flatten
land elevations, and 2) make depths positive. The SCRIP tools are used to
interpolate the bathymetry onto the target coordinates. The SCRIP tools use a
namelist file (``namelist_reshape_bilin_gebco``) which might need to be editted (if your application is a modification
from this example). Specifically it is assumed that the name of pre-processes
GEBCO ``input_file``, the
coordinates file, and the latitude and lonigitude variable names are correctly
specified in ``namelist_reshape_bilin_gebco`` as:

.. literalinclude:: ../START_FILES/DOMAIN/namelist_reshape_bilin_gebco
  :start-at: &grid_inputs
  :end-at: nemo_lat =

and are in the ``$DOMAIN`` directory. It is also that the elevation variabe is
appropriately defined so that it can be interpolated:

.. literalinclude:: ../START_FILES/DOMAIN/namelist_reshape_bilin_gebco
  :start-at: &interp_inputs
  :end-at: input_name =


.. literalinclude:: ../SCRIPTS/make_bathymetry_from_gebco.sh


bii. Generate bathymetry file from parent configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If alternatively the bathymetry is created from a parent model then the pre-processing
steps are not required.

.. literalinclude:: ../SCRIPTS/make_bathymetry_from_parent.sh


This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Creating bathymetry
  :end-at: make_bathymetry



c. Generate a domain configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Having created the horizontal coordinates, grid spacings and the bathymetry for
that grid, this step creates the vertical grid and aggregates all the grid
configuration variables into a single domain configuration file called ``domain_cfg.nc``.

In previous version of NEMO the vertical grid discretisation was computed at
model run time. For code management purposes it is now abstracted as a pre-processing
step, potentially making the choice of grid more flexible.

As with other steps in the configuration build process, this is also managed
by running some code with settings that are specified in a namelist file.

To date we have namelists for 3 types of vertical coordinate setups that we have
experimented with in various configurations. These are the traditional `stretched
sigma <../START_FILES/DOMAIN/s-sig_DOMAINcfg_namelist_cfg>`_ and
`z-partial step <../START_FILES/DOMAIN/z-ps_DOMAINcfg_namelist_cfg>`_, and also
a `hybrid s-z coordinate <../START_FILES/DOMAIN/hyb-z-s_DOMAINcfg_namelist_cfg>`_ system that we will
proceed with in this worked example configuration. Examples of the other namelists
are given in the ``START_FILES/DOMAIN`` directory of the NEMO-RELOC repository.

DISCUSSION ABOUT HYBRID Z-SIGMA COORDINATES.
Specific namelist settings include:

.. literalinclude:: ../START_FILES/DOMAIN/hyb-z-s_DOMAINcfg_namelist_cfg
  :start-at: &namzgr
  :end-at: ln_s_melange

The domain configuration file is generated with:

.. literalinclude:: ../SCRIPTS/make_domain_cfg.sh


This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Creating domain
  :end-at: make_domain

This creates the desired ``domain_cfg.nc`` files which contains all the vertical and
horizontal grid information and spacings.

==============================
Generate initial conditions
==============================

For a new configuration you probably want to start with idealised, or homogenous
initial conditions. This is done with user defined initial conditions ``ln_usr=T``
with the expression being compiled into the executable. (In ``$CDIR/$CONFIG/MY_SRC``:
``usrdef_istate.F90``).

To use initial conditions from an existing T,S field you might need to do a bit
of interpolation. It is advisable to let NEMO do the heavy lifting for vertical
interpolation (rquiring some FORTRAN modifictions), though SOSIE tools can be user
to do simple horizontal interpolation. See notes in

This will generates files typically called ``initcd_votemper.nc``
and ``initcd_vosaline.nc``, and corresponding mask and depth variables files
 ``initcd_mask.nc`` and ``initcd_depth.nc``.

---


The regional configuration can be initialised from user defined 3D tempearture
and salinity fields. During the configuration development and robustness checking
phase these initial conditions might be idealised, such as contant temperature
and salinity or
constant stratification. During production runs the initial conditions might be
taken from a parent simulation such as the ORCA12 global simulation or obtained
from CMEMS catalogue. In this description the velocities are initialised from zero.

Idealised temperature and salinity initial conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To initialise with idealised, or user-defined, temperature and salinity the source
code is compiled with ``usrdef_istate.F90`` (in ``MY_SRC``), wherein the user
predefines the fields before compilation. E.g.

.. literalinclude:: ../MY_SRC/usrdef_istate.F90_constTS
  :start-at:  usr_def_istate : analytical definition of initial state
  :end-before:  END SUBROUTINE usr_def_istate

This is activated at runtime by a logical parameter ``ln_tsd_init = .false.``
(i.e. do not use initial temperature and salinity data) in the simulation namelist ``namelist_cfg``.


Initial conditions from CMEMS database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The initial conditions from CMEMS database are generated with:

.. literalinclude:: ../SCRIPTS/make_CMEMS_IC.sh


This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Creating IC
  :end-at: make_CMEMS_IC

This will generates files typically called ``initcd_votemper.nc``
and ``initcd_vosaline.nc``, and corresponding mask and depth variables files
 ``initcd_mask.nc`` and ``initcd_depth.nc``.


======================================================
Atmospheric Forcing Data
======================================================

If ERA5 forcing data is used, these must be first obtained and pre-processed to
generate a land-sea mask

.. literalinclude:: ../SCRIPTS/prep_ERA5_forcing.sh

After pre-processing weights files must be generated to map the variables to the
target grid. This workflow is the same for pre-processed ERA5 data or for the
DFS atmospheric data that forced the global parent model.

.. literalinclude:: ../SCRIPTS/make_ERA5_weights.sh

In this demonstration we use the DFS forcing, which forced the parent global model.

** I NEED TO MAKE SURE THE REPO NAMELISTS ARE FOR DFS AND NOT ERA5 OR CMEMS DATA **

The regional configuration can be forced with atmospheric fields from external
files, or with idealised forcings. For example, during the configuration development and robustness checking
phase it is convenient to be able to switch off atmospheric forcing (I.e. set
forcing to zero). As with user defined initial conditions, a user-defined
surface boundary forcing can be set in ``usrdef_sbc.F90`` (in ``MY_SRC``),
which is fixed in the compiled code. The user-defined surface forcing is
activated by the logical parameter ``ln_usr=T`` in the simulation namelist ``namelist_cfg``.


======================================================
Generate Other Forcing Data
======================================================

a. Generate river forcing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Matlab script to generate river forcing <Generate_river_forcing.rst>`_

Link the river forcing output into ``$INPUTS`` e.g.::

 scp river_test.nc $USER@login.archer.ac.uk:$INPUTS/$CONFIG_rivers.nc

Old method: run GCOMS ``runoff_interactive_iterative.m`` to generate a
``rivers_test.nc`` file for freshwater forcing.


b. ERSEM boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To be written






======================================================
Open Boundary Conditions
======================================================

In this configuration a number of modifications were made to the tidal forcing.
The modifications are implemented in the FORTRAN source code as additions to the
``MY_SRC`` directory. Modifications relative to the NEMO trunk include:

  * Use of FES2014 tides
  * Restartable harmonic analysis
  * Addition of long period tides

These modifications are detailed in `Tidal code modifications <tidal_code_modifications.rst>`_


Discussion about PyNEMO.

`install PyNEMO <SCRIPTS/install_pynemo.sh>`_ which was called during the
`make_tools.sh` process.

Generate Tidal Open Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate the tidal boundary conditions, access to FES2014 data is required.
However, at present, everything is hardcoded inside PyNEMO, even the directory this data are (i.e., in
PyNEMO that you downloaded above in a file that is called
fes\_extract\_HC.py). So if you need to modify these things you will
have to modify the code of PyNEMO and particularly the
fes\_extract\_HC.py.

The following works for livljobs servers and it requires that you have
access to the files in directory: /projectsa/NEMO/Forcing/FES2014/

Get the namelist for tides:
`namelist\_FES14.bdy <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/namelist_FES14.bdy>`__

define some paths (based on where you have install your pynemo change
these paths)

.. code:: bash

    export WORK=/work/annkat
    export PYTHON=$WORK/nrct/Python
    cd $WORK

load your virtual environment for the pynemo that you have created in
step 1

.. code:: bash

    module load anaconda/2.1.0  # Want python2
    source activate nrct_env
    export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk-1.7.0.241-2.6.20.0.el7_7.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

change branch to the tides branch of pynemo

.. code:: bash

    git checkout Generalise-tide-input

build my PyNEMO and ready to use:

.. code:: bash

    cd $PYTHON
    python setup.py build
    export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
    python setup.py install --prefix ~/.conda/envs/nrct_env

Run pyNEMO

.. code:: bash

    pynemo -s namelist_FES14.bdy > Tide_pynemo.txt 2>&1



Generate 3D Open Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

a. Open boundary conditions from CMEMS data
+++++++++++++++++++++++++++++++++++++++++++

Having already downloaded CMEMS data for the initial conditions the first step
in creating time varying open boundary conditions is to generate mesh and mask
files for these data. This is advisable over using the mask file from CMEMS
since the regridding process onto c-points, for the served data, leads to
possible problems with the velocity data. We therefore generate mesh and mask
files directly from the downloaded variable files.

**Generate CMEMS mesh and mask files**

Using
`generate\_CMEMS\_coordinates.sh <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/generate_CMEMS_coordinates.sh>`__

Go to the directory you have downloaded your data:

.. code:: bash

    cd CMEMS/2017

load nco module:

.. code:: bash

    module load nco/gcc/4.4.2.ncwa

create the coordinates using nco tool:

.. code:: bash

    ncks -v longitude,latitude,depth,time 2017/CMEMS_2017_01_01_download.nc  tmp1.nc
    ncrename -d time,t -d latitude,y -d longitude,x tmp1.nc
    ncap2 -O -s 'glamt[t,y,x]=longitude' tmp1.nc
    ncap2 -O -s 'glamu[t,y,x]=longitude' tmp1.nc
    ncap2 -O -s 'glamv[t,y,x]=longitude' tmp1.nc
    ncap2 -O -s 'gphit[t,y,x]=latitude' tmp1.nc
    ncap2 -O -s 'gphiu[t,y,x]=latitude' tmp1.nc
    ncap2 -O -s 'gphiv[t,y,x]=latitude' tmp1.nc
    mv tmp1.nc CMEMS_subdomain_coordinates.nc

create the masks:

.. code:: bash

    module load nco/gcc/4.4.2.ncwa
    ncks -v longitude,latitude,depth,time,so 2017/CMEMS_2017_01_01_download.nc tmp2.nc
    ncks -A -v uo,vo 2017/CMEMS_2017_01_01_download_UV.nc tmp2.nc
    ncatted -a _FillValue,,d,, tmp2.nc
    ncap2 -O -s 'where(so>0.) so=1' tmp2.nc tmp2.nc
    ncap2 -O -s 'where(so<=0.) so=0' tmp2.nc tmp2.nc
    ncap2 -O -s 'where(uo>-10.) uo=1' tmp2.nc tmp2.nc
    ncap2 -O -s 'where(uo<=-10.) uo=0' tmp2.nc tmp2.nc
    ncap2 -O -s 'where(vo>-10.) vo=1' tmp2.nc tmp2.nc
    ncap2 -O -s 'where(vo<=-10.) vo=0' tmp2.nc tmp2.nc
    ncrename -d time,t -d latitude,y -d longitude,x tmp2.nc
    ncrename -v so,tmask tmp2.nc
    ncrename -v uo,umask tmp2.nc
    ncrename -v vo,vmask tmp2.nc
    mv tmp2.nc CMEMS_subdomain_mask.nc

**Generate CMEMS OBC**

using
`Run\_Pynemo.sh <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/Run_Pynemo.sh>`__

Make sure that you have your domain file (domain\_cfg\_ORCA12.nc) and
your bathymetry of your domain file (bathy\_meter\_ORCA12.nc). You
should have generated these fields previously.

get the namelists that you will need for pynemo:
`namelist\_2017 <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/namelist_2017.bdy>`__,\ `inputs\_src\_zgr.ncml <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/inputs_src_zgr.ncml>`__,
`inputs\_dst.ncml <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/inputs_dst.ncml>`__.
Also get the namelist that defines the input in
`NCML <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/tree/master/FILES_START/OBC-TIDES/NCML>`__,
for example for 2017 get
`CMEMS\_2017.ncml <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/NCML/CMEMS_2017.ncml>`__.
**ATTENTION** modify them accordingly to reflect the name of your files
and directories and years/months you may want to use.

load anaconda and the virtual environment you have created associated
with installing PyNEMO in step 1:

.. code:: bash

    module load anaconda/2.1.0  # Want python2
    source activate nrct_env
    export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

Run your PyNEMO to generate your open boundary conditions for T,S,
U,V,SSH (note this is slow for a whole year it may take a couple of
hours):

.. code:: bash

    pynemo -s namelist_2017.bdy > T2017_pynemo.txt 2>&1



b. Open boundary conditions from parent ORCA12 data
+++++++++++++++++++++++++++++++++++++++++++

======================================================
Run Experiments
======================================================

a. Tide only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


directory ``EXP_tideonly``

Only tidal forcing. constant T and S
Include: key_harm_ana
EXEC: nemo_tideonyl_TSconst.exe


Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  cp $CONFIG/MY_SRC/usrdef_istate.F90_horizTS $CONFIG/MY_SRC/usrdef_istate.F90

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move executable to something permanent::

  cd $CDIR/$CONFIG/BLD/bin
  mv nemo.exe nemo_tideonly_TSconst.exe

  cd $CDIR/$CONFIG/EXP_tideonly
  ln -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tideonly_TSconst.exe opa

Edit namelist_cfg::

  nn_itend = 14400
  ln_rstart = .false.
  ...
  ln_tide     = .true.
  ln_tide_pot = .true.    !  use tidal potential forcing
  ...
  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state

Didn't bother with the tidal harmonics. It will run but I am spinning up.

Run for 4 hours. nt = 14400, dt =360, 60 days. Completed in 3hr 31.

Edit XML output to produce 5d output.
Resubmit::

  qsub runscript

*(23 Mar 2018)*
Turn on 19 harmonics using the POLCOMS harmonic analysis (Nico's instructions and edits)
Run for another 60 days with harmonic analysis restarting capbability.

Works. Did tidal analysis plots ::

  ~/GitLab/JMMP_tools
  python Tidal_analysis_amplitude.py --verbose
  python Tidal_analysis_plot.py --verbose


Need to continue run

::

  namelist_cfg

  ...
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  28801   !  first time step
     nn_itend    =  43200 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00028800_restart"   !  suffix of ocean restart name (input)

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                   !
                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.                   ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                     ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.true.                      ! Read haronic analyisis from a restart



   !-----------------------------------------------------------------------
   &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
   !-----------------------------------------------------------------------
       nit000_han = 28801         ! First time step used for harmonic analysis
       nitend_han = 43200 ! 1440 !      ! Last time step used for harmonic analysis



Submit::

  qsub runscript

Took 3:36 mins. The tidal analysis (https://www.evernote.com/shard/s652/nl/85824826/674115d9-29be-480a-ba71-6814de98df4b/) doesn't show significant improvement
::

  ~/GitLab/JMMP_tools
  python Tidal_analysis_amplitude.py --verbose
  python Tidal_analysis_plot.py --verbose

But 4 months of simulation might still be on the short side. Run for another two months.
::

  vi namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  43201   !  first time step
     nn_itend    =  57600 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00043200_restart"   !  suffix of ocean restart name (input)




   !-----------------------------------------------------------------------
   &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
   !-----------------------------------------------------------------------
       nit000_han = 43201         ! First time step used for harmonic analysis
       nitend_han = 57600 ! 1440 !      ! Last time step used for harmonic analysis



Submit::

  qsub runscript

Ran with minor change to statistcs.



b. No tidal forcing
^^^^^^^^^^^^^^^^^^^^^^^^^^

Need to turn off tidal forcing. See ``&nam_tide``::

  ln_tide     = .false.
  ln_tide_pot = .false.    !  use tidal potential forcing


Comment out tidal analysis line ::

  vi ../cpp_SEAsia.fcm
    bld::tool::fppkeys key_zdfgls        \
  //REMOVE//                 key_harm_ana       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Will have to recompile.


c. HPG test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Submit a 30 day simulation, from rest, with depth varying spatially homogeneous
temperature and salinity profiles, with no forcing, boundary conditions off:
``ln_bdy = F``


Edit runscript: 2hrs walltime. It took 1h 50mins

Edit namelist_cfg
360s, nt=7200 --> 30 days::

  ln_tide     = .false.
  ln_tide_pot = .false.    !  use tidal potential forcing

cd $EXP/../EXP_hpg_err

Scrape ``umax`` from ``solver.state`` using plot_solver_stat.py

Some along rim currents started but these are small compared to interior currents.
Restart for another 30 days.
After 30 days umax is still growing. Restart run and continue::

  mv solver.stat solver.stat_part1

Check progress with::

   hpg_error_plotNEMO.py
   plot_solver_stat.py

Edit namelist_cfg for restarting::

  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  7201   !  first time step
     nn_itend    =  14400 ! 30day=7200   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        cn_ocerst_in    = "SEAsia_00007200_restart"   !  suffix of ocean restart name (input)


Resubmit::

  qsub runscript

Ran in 1hr 46

Joing the solver.stat files together::

cp solver.stat solver.stat_part2
cp solver.stat_part1 solver.stat
cat solver.stat_part2 >> solver.stat

module load anaconda
python plot_solver_stat.py


..note::

  In the end I didn't do the restart for the last run I did.



d. Horizontally constant T and S initial condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Emulate James ORCHESTRA method (first moved usrdef_istate.F90 to usrdef_istate.F90_constTS for safe keeping)::

  cd $CDIR/$CONFIG/MY_SRC
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/usrdef_istate.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/lapack.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/splines.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/types.F90 .
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/utils.F90 .

Edit usrdef_istate.F90 to put in profile
Data from Hamburg WOCE Climatology Live Access Server at (-2N, 95E).
With constant variable padding below 4000m to make it up to 75 levels::

      zdep(:) = (/     0,    10,    20,    30,    40,    50,   75,    100,   125,   &
                &    150,   175,   200,   250,   300,   350,  400,    500,   600,   &
                &    700,   800,   900,  1000,  1100,  1200,  1300,  1400,  1500,   &
                &   1750,  2000,  2250,  2500,  2750,  3000,  3250,  3500,  3750,   &
                &   4000,  4100,  4200,  4300,  4400,  4500,  4550,  4600,  4700,   &
                &   4800,  4900,  5000,  5100,  5200,  5300,  5400,  5500,  5600,   &
                &   5700,  5800,  5900,  6000,  6100,  6200,  6300,  6400,  6500,   &
                &   6600,  6700,  6800,  6900,  7000,  7100,  7200,  7300,  7400,   &
                &   7500,  7600,  7700 /)

      zsal(:) = (/ 34.05, 34.05, 34.10, 34.13, 34.25, 34.42, 34.88, 35.08, 35.13,   &
                &  35.08, 35.07, 35.06, 35.06, 35.03, 35.01, 34.99, 34.96, 34.97,   &
                &  34.97, 34.95, 34.92, 34.91, 34.88, 34.87, 34.85, 34.83, 34.82,   &
                &  34.80, 34.77, 34.76, 34.75, 34.74, 34.73, 34.73, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72 /)

      ztmp(:) = (/ 28.87, 28.87, 28.87, 28.74, 28.33, 28.01, 25.21, 21.99, 18.51,   &
                &  15.55, 14.39, 13.43, 12.27, 11.48, 11.02, 10.51,  9.58,  8.95,   &
                &   8.35,  7.78,  7.16,  6.52,  5.88,  5.44,  5.02,  4.57,  4.14,   &
                &   3.34,  2.64,  2.31,  2.05,  1.86,  1.69,  1.58,  1.41,  1.23,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15  /)

I fixed the top two layers to be constant salinity and temperature to avoid possible cabling.
The raw data had inversions at 10m

NB avoding mistakes: namely a duplicate depth in the input depth funciton --> spline failure


e. Set boudaries to initial condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(nn_dyn2d_dta was 2)::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .false.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state



Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move the executable to a special name::

  mv $CONFIG/BLD/bin/nemo.exe $CONFIG/BLD/bin/nemo_notide_TSprofile.exe

Move to experiment dir and link executable file in::

  cd $EXP/../EXP-hpg_err
  ln -s $CDIR/$CONFIG/BLD/bin/nemo_notide_TSprofile.exe opa



James spotted that I didn't have lateral diffusion of momentum. Made some changes (following ORCHESTRA namelist_cfg)
Submitted run (EXP01) to test timestep. rn_rdt=360 ran 1304 in 20mins ==> 5.4 days

There an SSH instability in the NE corner when using::

 cn_dyn3d      =  'neumann'

Switch to::

 cn_dyn3d      =  'zerograd'

NB Tried masking, using mask_bdy.nc, but this didn't work.

Cold start run::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  01   !  first time step
     nn_itend    =  7200 ! 30day=7200   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)


Run::

  qsub runscript


Yes. That works.

Though the SSS min decreases by 0.017 over 30 days. A bit odd? Perhaps it is
the non conservative nature of the advective schemes...

NB The effect is coastal. There is no problem around the boundaries.






f. Internal tide with idealised stratification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Try with tides turned on.
(Recompile a tide and no tide version. Save in $CONFIG/BLD/bin as nemo_tide.exe
and nemo_notide.exe, then link as appropriate)::

  cd EXP00
  ls -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tide.exe opa

  <    ln_tide     = .true.
  <    ln_tide_pot = .true.    !  use tidal potential forcing
  ...
  <     nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state

  qsub -q short runscript





g. Tide with idealised stratification and rivers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

directory ``EXP_tide_TSprofile_river``

Only tidal forcing. profile of T and S
Include: key_harm_ana
EXEC: nemo_tide_TSprofile_nomet.exe

::

  cp EXP_tideonly/* EXP_tide_TSprofile_river/.
  ln -s $INPUTS bdydta

Recompile the code.
Resubmit with dt=60s and nt = 60 (ie, 1 hr)::

  cd $CDIR
  cp $CONFIG/MY_SRC/usrdef_istate.F90_horizTS $CONFIG/MY_SRC/usrdef_istate.F90

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move executable to something permanent::

  cd $CDIR/$CONFIG/BLD/bin
  mv nemo.exe nemo_tide_TSprofile_nomet.exe

  cd $CDIR/$CONFIG/EXP_tide_TSprofile_river
  ln -s /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/BLD/bin/nemo_tide_TSprofile_nomet.exe opa

Edit namelist_cfg::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =    "SEAsia"  !  experience name
     nn_it000    =  57601   !  first time step
     nn_itend    =  58560 ! 10day=14400   !  last  time step (std 5475)
     nn_date0    =  20000102   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
     nn_time0    =       0   !  initial time of day in hhmm
     nn_leapy    =       1   !  Leap year calendar (1) or not (0)
     ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
        !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
        !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
        cn_ocerst_in    = "SEAsia_00057600_restart"   !  suffix of ocean restart name (input)

  ...

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing
  ...

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_bdy         = .true.              !  Use unstructured open boundaries
      nb_bdy         = 1                    !  number of open boundary sets
      ln_coords_file = .true.               !  =T : read bdy coordinates from file
      cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      cn_dyn2d       = 'flather'               !
      nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state
      ...

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                   !
                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.                   ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                     ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.true.                      ! Read haronic analyisis from a restart

  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 57601         ! First time step used for harmonic analysis
      nitend_han = 58560 ! 1440 !      ! Last time step used for harmonic analysis




Finally turn rivers on::

  !-----------------------------------------------------------------------
  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
  ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)
  ...
  !-----------------------------------------------------------------------
  &namsbc_rnf    !   runoffs namelist surface boundary condition          (ln_rnf=T)
  !-----------------------------------------------------------------------
  !              !  file name           ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                      !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_rnf      = 'bdydta/SEAsia_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_cnf      = 'bdydta/SEAsia_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''
     sn_s_rnf    = 'runoffs'            ,        24         , 'rosaline',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_t_rnf    = 'runoffs'            ,        24         , 'rotemper',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_dep_rnf  = 'runoffs'            ,         0         , 'rodepth' ,   .false.    , .true. , 'yearly'  , ''       , ''       , ''

     cn_dir      = './'      !  root directory for the location of the runoff files
     ln_rnf_mouth= .false.    !  specific treatment at rivers mouths
     ...
     ln_rnf_depth_ini = .true. ! compute depth at initialisation from runoff file
      rn_rnf_max  = 0.6865698 ! 5.735e-4  !  max value of the runoff climatologie over global domain ( ln_rnf_depth_ini = .true )
      rn_dep_max  = 150.      !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )

Had to add dir path to filename otherwise file was not found...

Edit runscript::

  vi runscript
  #PBS -l walltime=00:30:00

Resubmit::

  qsub runscript

**PENDING**
*(27 Apr 2018)*

/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_tide_TSprofile_river

Run for 30 mins. nt = 960, dt =360, 4 days. Completed in  15 mins 06s.

CHECK OUTPUT - This works, but 4 days is not very exciting.


h. Full forcing: Initial conditions, tides, rivers, open boundaries + met
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

directory: EXP_fullforcing

In v4 SLP is needed to calculate air density. However in the old code it appears that
air density was set at a constant, in the bulk core formulations, set at 1.22
Plan:
* Replace SLP variable in the namelist with a dummy variable, like temperature, so that the code does not break
* In the fortran fix air density to 1.22



EXEC: ``nemo_FES14-tides_diaharm-fast.exe``

Note this is using old style of tidal implementation and can be updated to nemo_FES14-tides_diaharm-fast.exe later.
::
  mkdir EXP_fullforcing




i. Physics and biogeochemisty
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
