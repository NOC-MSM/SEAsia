
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

This is activated at runtime by a logical parameter ``ln_usr=T`` in the simulation
namelist ``namelist_cfg``.


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


c. Something about tidal forcing modification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Use of FES tides
* Restartable harmonics
* Long period tides (though not visible).

Update tides code with Nico's version.
++++++++++++++++++++++++++++++++++++++

Add the POLCOMS harmonic analysis to the executable (as in now in)
`<build_opa_orchestra.rst>`_
This requires some changes to the standard ``namelist_cfg``

Add the final (extra) three variables in your namelist_cfg / nambdy_tide ::

  vi $EXP/namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                    !
                                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.          ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                 ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.false.                    ! Read haronic analyisis from a restart
  /


Edit xml files to output harmonics as amplitudes and phases (e.g.)::

  vi file_def_nemo.xml
  ...
  <file_group id="tidal_harmonics" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1d files -->
    <file id="tidalanalysis.grid_T" name="harmonic_grid_T" description="ocean T grid variables"  enabled=".TRUE.">

      <field field_ref="O1amp"         name="O1amp"       operation="instant" enabled=".TRUE." />
      <field field_ref="O1phase"       name="O1phase"     operation="instant" enabled=".TRUE." />


  vi field_def_nemo-opa.xml
  ...
      <field_group id="Tides_T" grid_ref="grid_T_2D" operation="once" >
      <!-- tidal composante -->
      ...
      <field id="Q1amp"        long_name="Q1 Elevation harmonic Amplitude"                              unit="m"        />
      <field id="Q1phase"      long_name="Q1 Elevation harmonic Phase"                                  unit="degree"   />

*Recall there are elevation, u-vel and v-vel harmonics*. Also editted suffixes
 in velocity fields, adding ``_2D``.


* As before the constituents you want to analyse are set-up in ``nam_diaharm``
 namelist.

* The harmonic analysis is done at the end only as well as the restart dumping
so you can only restart from the last time step so make sure you output the full
 restart at the end. To restart, you just need to turn on the ``ln_harmana_read``
  and to map the files to something like ``restart_harm_ana_*``  as this bit as
   not been developed with a prefix to load the files. You can look at this
    python script if needed:
  ``/work/n01/n01/nibrun/RUNS/SWPacific/SIMU/01_harm_links.py``



======================================================
Generate Open Boundary Conditions
======================================================

Discussion about PyNEMO.

a. Open boundary conditions from CMEMS data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

b. Open boundary conditions from parent ORCA12 data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Copy files from EXP_fullocean to EXP_fullforcing **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_fullocean/* EXP_fullforcing/

Need to add in:

* met forcing

Get the executable in place (I think that the standard nemo.exe is OK)::

  cd EXP_fullforcing
  rm opa
  ln -s ../BLD/bin/nemo_FES14-tides_diaharm-fast.exe opa

Edit namelist_cfg::

  ln_usr      = .false.    !  user defined formulation                  (T => check usrdef_sbc)
  ln_flx      = .false.   !  flux formulation                          (T => fill namsbc_flx )
  ln_blk      = .true.    !  Bulk formulation                          (T => fill namsbc_blk )
  ...
  ln_NCAR     = .true.   ! "NCAR"      algorithm   (Large and Yeager 2008)


Add in a dummy variable for Sea Level Pressure. This is expected. If it was real
it would be used to calculate atmospheric density. Instead this is fixed to a constant = 1.22::

  sn_slp      = 'cutdown_drowned_t2_DFS5.1.1'        ,         3         , 't2',   .false.    , .false. , 'yearly'  , 'weights_bilinear_atmos.nc'  , ''       , '' ! DUMMY variable. Fill with T

Check the reference heights for the wind and T,q data::

  rn_zqt      = 2.       !  Air temperature and humidity reference height (m)
  rn_zu       = 10.       !  Wind vector reference height (m)


Edit ``MY_SRC/sbcblk.F90`` to make air desnity constnat::

  vi sbcblk.F90
  FUNCTION rho_air( ptak, pqa, pslp )
     ...
     !
     rho_air = 1.22 ! pslp / (  R_dry*ptak * ( 1._wp + rctv0*pqa )  ) ! jelt: 25 Oct 2018. Fix to make air density constant.
     !
  END FUNCTION rho_air


Recompile::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Update executables so that:: ``../EXP_fullforcing/opa -> ../BLD/bin/nemo_FES14-tides_diaharm-fast.exe``

submit on short queue::

  qsub runscript_short
  5695821.sdb

Rivers get really hot. Try turning them off...

ln_rnf      = .false.    !  runoffs                                   (T => fill namsbc_rnf)

The SST heats up quite a but - but I am not sure about when my initial condition was
taken. Temperatures are all sensible at least...

Try and fix the rivers. Turn them on again and address how the rivers enter the domain

Edit namelist_cfg::

  ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)
  ...
  sn_rnf      = 'bdydta/SEAsia_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
  sn_cnf      = 'bdydta/SEAsia_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''
  sn_s_rnf    = 'runoffs'            ,        24         , 'rosaline',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
  sn_t_rnf    = 'runoffs'            ,        24         , 'rotemper',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
  sn_dep_rnf  = 'runoffs'            ,         0         , 'rodepth' ,   .false.    , .true. , 'yearly'  , ''       , ''       , ''

  cn_dir      = './'      !  root directory for the location of the runoff files
  ...
  ln_rnf_depth_ini = .true. ! compute depth at initialisation from runoff file
   rn_rnf_max  = 0.6865698 ! 5.735e-4  !  max value of the runoff climatologie over global domain ( ln_rnf_depth_ini = .true )
   rn_dep_max  = 150.      !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )

Had to add dir name into filename as it couldn't otherwise find the file...!
rn_rnf_max value was calculated by loading in the river file and summing over
all space to get a sum per month and picking the max value
::
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset('bdydta/SEAsia_rivers.nc')
    np.sum(np.sum(ro,axis=2),axis=1)


qsub runscript_short

Ran for 1162 w/ dt=360. Blew up with vel
kt=  1162 max abs(U):   2.633    , i j k:   535  320    1

Resubmit with dt=100s. Runs for 20mins (1170 steps) without blow up.

Ran until 1221 for 20mins with 200s Steps

Resubmit to standard queue. (2hrs)::

 stpctl: the speed is larger than 20 m/s
 ======
 kt=  2532 max abs(U):   2.841    , i j k:   535  320    1


Ah. The error trapping assumes vel > 20, when it is triggering for vel^2 > 20.
Fix in ``MY_SRC/stpctl.F90`` and recompile::

    cd $CDIR
    ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Save in ``BLD/bin`` as ``nemo_FES14-tides_diaharm-fast.exe``

Resubmit
PENDING **5748307.sdb **

It works! Though the SST does heat up a bit (about 4 deg in 9 days)

----

*(1 March 2019)*
Edit namelist for queue
Submit 6026543.sdb
Ran to 20min completion. (2.8 days)

Submit for 2hrs: 6026571.sdb
**THIS WORKED. HAVE STABLE RESTARTS**

*(10 May 2019)*
Last run seemed to have crashed with a  core dump. Try the 20min short queue
again,
XIOS had been updated.
6195159.sdb

Ran for 20 mins (2.9 days)

Try from Jan 1960
Generate initial conditions from 1 5d output in Jan 1960
`<generate_initial_conditions.rst>`_

Generate open boundary conditions
`<generate_openboundaryconditions.rst>`_

Need to copy pynemo files across and submit::

livljobs4:
cd /scratch/jelt
rsync -uvt jelt@jasmin-xfer1.ceda.ac.uk:/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/coordinates.bdy.nc .
rsync -uvt jelt@jasmin-xfer1.ceda.ac.uk:/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/SEAsia_bdyT_y1960m01.nc .
rsync -uvt jelt@jasmin-xfer1.ceda.ac.uk:/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/SEAsia_bdyU_y1960m01.nc .
rsync -uvt jelt@jasmin-xfer1.ceda.ac.uk:/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/SEAsia_bdyV_y1960m01.nc .
rsync -uvt jelt@jasmin-xfer1.ceda.ac.uk:/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/SEAsia_bt_bdyT_y1960m01.nc .

rsync -uvt coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.
rsync -uvt SEAsia_bdyT_y1960m01.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.
rsync -uvt SEAsia_bdyU_y1960m01.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.
rsync -uvt SEAsia_bdyV_y1960m01.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.
rsync -uvt SEAsia_bt_bdyT_y1960m01.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.

rm SEAsia_bdyT_y1960m01.nc SEAsia_bdyU_y1960m01.nc SEAsia_bdyV_y1960m01.nc
rm SEAsia_bt_bdyT_y1960m01.nc coordinates.bdy.nc


Qsub 6195579.sdb

::

  stp_ctl : the ssh is larger than 10m
   =======
   kt=   535 max ssh:   16.00    , i j:   682  351

Stable with tides + 2D forcing (ec87ad0)

It looks like the surface heating is too strong.

Made some modifications to light penetration schemes (in line with N06
namelist_cfg `<http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N06/domain/namelist_cfg>`_)
git:commit f0fcc7f

*PENDING* 6197300.sdb
Does it improve the SST drift?
Also noted that N06 uses SST and SSS restoring (it looks to a LEVITUS climatology
monthly for sss and contant for sst...)
This has not been implemented.

Drift still occurs at the surface.
Try and turn on later baroclinic boundaries just to check it is not a vertical
discretisation issue.
kt=     7 max abs(vel)**2:   3099.    , i j k:   682  338   58

I think that the intial conditions and bcs are on different vertical grids.
Test with no mometum at the boundary just temperature.
It runs with initial condtions T/S + 3D vel

Try with vel=0 but T/S at bdy
 kt=     7 max abs(vel)**2:   3099.    , i j k:   682  338   58
It looks like a sea mount near the boundary. Smooth or something else (use ORCA bathy?)

Try running for 4 days (960 steps) with specified baroclinc T/S at initial values.
Does it iron out the sea mount issue for a restart with real T/S?
Looks OK. Restart with 3d vel turned on...
Blows up again...
Prbably need to sort out the bathmetry.
Having said all that it looks like specified open boundary conditions might help.
Though NaNs do emerge in solver.stat whilst remaining stable....

----

**Notes to add git to version control EXP FILES**

# Now change to CONFIG directory
cd NEMOGCM/CONFIG

# Checkout configuration directory structure
git init .
git remote add origin git@github.com:NOC-MSM/NEMO_cfgs.git
git config core.sparsecheckout true
echo "SEAsia/*" >> .git/info/sparse-checkout
git pull --depth=1 origin master

Populate a git ignore file::

  # ignores all *nc files recursive from the current folder
  **nc

  # Ignore all log and error files
  **/SEAsia.o*
  **/SEAsia.e*

  # Ingore forcing data
  **/bdydta/

Add MY_SRC
Add README

**Next Steps:**

* Also new plan to just use Global rivers rather than try and create a new (and
different) set.

* Use global met, not cutout DFS files. LINK IN DFS5.2 into INPUTS.

* Add important EXP files to git. Add repo to GitHub.



i. Physics and biogeochemisty
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




======================================================
OLD STUFF
======================================================






Initial conditions
++++++++++++++++++

directory: EXP_hpg_initcd

Switch in initial conditions from the existing (working) tide only homogeneous run.
::
  mkdir EXP_hpg_initcd

Copy files from ??? to EXP_hpg_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_tideonly/* EXP_hpg_initcd/


Try from rest.


Rerun the last successful tidal only simulation but with initial conditions.
Switch off bcs and tides run for 60 days::


  vi namelist_cfg

  ...
  cn_exp      =    "SEAsia"  !  experience name
  nn_it000    =  1   !  first time step
  nn_itend    =  7200 ! 30day=7200; dt=360  !  last  time step (std 5475)

  ln_restart = .false.

Turn off open boundaries::

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
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state


Turn off tides::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .false.
     ln_tide_pot = .false.    !  use tidal potential forcing

Make sure the harmonic restart is also off (not so important yet since ln_tide=F)::

  ln_harmana_read=.false.                      ! Read haronic analyisis from a restart
  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 1         ! First time step used for harmonic analysis
      nitend_han = 2400 ! 1440 !      ! Last time step used for harmonic analysis

Turn on intial conditions::

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


Add some daily output::

  vi file_def_nemo.xml

    <file_group id="1d" output_freq="1d"  output_level="10" enabled=".TRUE.">  <!-- 1d files -->
    <file id="file19" name_suffix="_SS" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"   />
     <field field_ref="sss"          name="sss"   />
    </file>
  </file_group>


Link in missing sybolic link files::

  rm opa
  ln -s ../BLD/bin/nemo_tide_nomet.exe opa

I also fudged the dates on the boundary conditions files::

  cd $INPUTS
  ln -s SEAsia_bdyT_y1979m11.nc SEAsia_bdyT_y2000m01.nc
  ln -s SEAsia_bdyU_y1979m11.nc SEAsia_bdyU_y2000m01.nc
  ln -s SEAsia_bdyV_y1979m11.nc SEAsia_bdyV_y2000m01.nc
  ln -s SEAsia_bt_bdyT_y1979m11.nc SEAsia_bt_bdyT_y2000m01.nc

Run on short queue::

  cd SEAsia/EXP_tide_initcd
  qsub runscript_short

Later in evening run on standard queue (2hrs)::

  cd SEAsia/EXP_tide_initcd
  qsub runscript


**THIS WORKS** 13 May 2018



Initial conditions + tides
--------------------------

directory: EXP_tide_initcd

EXEC: ``nemo_tide_nomet.exe``

Switch in initial conditions from the existing (working) tide only homogeneous run.
::
  mkdir EXP_tide_initcd

Copy files from ??? to EXP_tide_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_hpg_initcd/* EXP_tide_initcd/


Try from rest.::

  vi namelist_cfg

  ...
  cn_exp      =    "SEAsia"  !  experience name
  nn_it000    =  1   !  first time step
  nn_itend    =  7200 ! 30day=7200; dt=360  !  last  time step (std 5475)

  ln_restart = .false.

Turn off open boundaries::

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
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'zerograd'               !
      nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state


Turn on tides::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing

Make sure the harmonic restart is also off ::

  ln_harmana_read=.false.                      ! Read haronic analyisis from a restart
  ...
  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 1         ! First time step used for harmonic analysis
      nitend_han = 2400 ! 1440 !      ! Last time step used for harmonic analysis



cd /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_tide_initcd

**PENDING** 5324836.sdb
13 May 2018

completed 1hr 39




Try lateral boundary conditions T,S,u
++++++++++++++++++++++++++++++++++++++


**Next steps**

* lateral boundary conditions (T,S,u)
* tides + lateral boundary conditions


---

*(24 Apr 2018)*
Build the 3D boundary conditions from ORCA0083-N01
==================================================

**The following gets a bit messy towards the end of this section and needs tidying up. But it seems to work**

Use PyNEMO to generate 3D bcs on livljobs4::

  ssh livljobs4

  . ~/temporary_path_names_for_NEMO_build
  cd $INPUTS

Edit the namelist.bdy for 3D fields. There are a few key things to note:

* Turn off the tides ``ln_tide =.false.`` and change ``rimwidth``: 1 --> 9
(tides don't work with rimwidth != 1 and baroclinc bcs don't work with rimwidth=1)

* For the ORCA0083-N001 run there was a glitch with the model timestamp. It is
 out by 36 days and 16 hours. This can be adjusted with a new namelist variable
 ``nn_src_time_adj``.

* The parent files could not be accessed on the thredd server because it is broken.
 For now I've copied Nov 1978 into ``/projectsa/accord/ORCA0083``. This is reflected
 in local_inputs_src.ncml

* Make sure the namelist time variables (years and months) match the timestamps in the files you
 load otherwise not much happens...

* You might need to have the time data spanning a whole month, so that a month
 can be extracted.

* Finally, there may have been an issue with a pynemo datetime function that
 determined the datetime for end of the month in profile.py
I  made change to the master branch to use datetime and timedelta instead of
 netcdfdatime, or something. This only matters is pynemo breaks with an
associated datetime error. Then this is the fix.

::

  vi namelist.bdy

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
     sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
     sn_src_zgr = './mesh_zgr_src.nc'   !  parent
     sn_dst_hgr = './domain_cfg.nc'
     sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
     sn_src_msk = './mask_src.nc'       ! parent
     sn_bathy   = './bathy_meter.nc'

  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
     sn_src_dir = './local_inputs_src.ncml'       ! src_files/'
     sn_dst_dir = '/work/jelt/NEMO/SEAsia/INPUTS/'
     sn_fn      = 'SEAsia'                 ! prefix for output files
     nn_fv      = -1e20                     !  set fill value for output files
     nn_src_time_adj = 0                                    ! src time adjustment
     sn_dst_metainfo = 'metadata info: jelt'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_coords_file = .true.               !  =T : produce bdy coordinates files
      cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
      ln_tra         = .true.               !  boundary conditions for T and S
      ln_ice         = .false.               !  ice boundary condition
      nn_rimwidth    = 9                    !  width of the relaxation zone

  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
      ln_tide        = .false.               !  =T : produce bdy tidal conditions
      clname(1) =  'M2'
      clname(2) =  'S2'
      clname(3) =  'N2'
      clname(4) =  'K2'
      clname(5) =  'K1'
      clname(6) =  'O1'
      clname(7) =  'P1'
      clname(8) =  'Q1'
      clname(9) =  'M4'
      ln_trans       = .false.
      sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
      sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'

  !-----------------------------------------------------------------------
  !  Time information
  !-----------------------------------------------------------------------
      nn_year_000     = 1979        !  year start
      nn_year_end     = 1979        !  year end
      nn_month_000    = 11          !  month start (default = 1 is years>1)
      nn_month_end    = 11         !  month end (default = 12 is years>1)
      sn_dst_calendar = 'gregorian' !  output calendar format
      nn_base_year    = 1978        !  base year for time counter
      sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
      nn_src_time_adj    = -3254400  ! -3168000 - 86400 ! fix to align model time stamp
  !-----------------------------------------------------------------------
  !  Additional parameters
  !-----------------------------------------------------------------------
      nn_wei  = 1                   !  smoothing filter weights
      rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                    !  smoothing onto dst points. Need to
                                    !  make this a funct. of dlon
      sn_history  = 'bdy files produced by jelt from ORCA0083-N01'
                                    !  history for netcdf file
      ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
      nn_alpha    = 0               !  Euler rotation angle
      nn_beta     = 0               !  Euler rotation angle
      nn_gamma    = 0               !  Euler rotation angle
      rn_mask_max_depth = 300.0     !  Maximum depth to be ignored for the mask
      rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break


The ORCA0083-N01 parent files. (This could be streamlined with some wildcards)::

  vi local_inputs_src.ncml::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05U.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05V.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791111d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791116d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791121d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791126d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791201d05T.nc" />
          <ns0:netcdf location="file:/projectsa/accord/ORCA0083/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>


Run pynemo (on branch ORCA0083)::


  cd $INPUTS

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  pynemo -s namelist.bdy


*(25 Apr 2018)*

* ``nemo_bdy_extr_tm3.py`` Editted ind_NAN variable and caught exception for
empty index list. Resubmitted.

Outputs::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc

This took 2 hours! To produce a month of boundary conditions.


Copy files from SAN to ARCHER::

  livljobs4
  cd $INPUTS
  for file in  SEAsia_b*nc; do rsync -uvt $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done





Tide + 2D forcing simulation
++++++++++++++++++++++++++++

directory ``EXP_openbcs``

EXEC: ``nemo_tide_nomet.exe``


Switch in initial conditions from the existing (working) tide only homogeneous
 run::

  mkdir EXP_openbcs
  cd EXP_openbcs

Copy files from ??? to EXP_tide_initcd **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=$CONFIG*.nc EXP_tide_initcd/* EXP_openbcs/

Link in missing sybolic link files::

  rm opa
  ln -s ../BLD/bin/nemo_tide_nomet.exe opa

Try from rest.
First try tides only with ln_bdy=T. ``nn_dyn2d_dta   =  2``. Ran to 400 steps (40hrs) OK
Second try with 2D bcs not tides: ``nn_dyn2d_dta   =  1``.

===>>> : E R R O R
        ===========

 stpctl: the speed is larger than 20 m/s
 ======
kt=     1 max abs(U):  1.0000E+20, i j k:   288    2    1

          output of last fields in numwso

--- **PROBLEM WITH 2D BCS**

Third try with tides + tracers::

  cn_dyn2d       = 'flather'               !
  nn_dyn2d_dta   =  2                   !  = 0, bdy data are equal to the initial state
  cn_dyn3d      =  'zerograd'               !
  nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
  cn_tra        =  'frs'               !
  nn_tra_dta    =  1                    !  = 0, bdy data are equal to the initial state

Also blows immediately.

Revisit PyNEMO. Turn off tidal output and run.

----

trouble
+++++++

19 May 2018

Make a new mask file to only process the southern boundary.
Turn off vertical interpolation in PyNEMO so that it creates output on the parent
vertical grid.

livljobs4/6


The mask variable takes values (-1 mask, 1 wet/active domain, 0 land). Need to
only mask a single point around the edge since the rimwidth is considered to be
part of the active domain.

PyNEMO looks for the interface between -1 and 1 to generate boundary forcings. Get a
template from domain_cfg.nc and then modify as desired around the boundary.

Maked the entire boundary "land=0", except for the wet bits along the southern boundary
which are "mask=-1"::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython::

import netCDF4
import numpy as np
dset = netCDF4.Dataset('bdy_mask.nc','a')
[ny,nx] = np.shape(dset.variables['mask'][:])
for i in range(nx):
  if dset.variables['mask'][1,i] == 1:
    dset.variables['mask'][0,i] = -1
  else:
    dset.variables['mask'][0,i] = 0
dset.close()

Rebuld pynemo as necessary `<install_nrct.rst>`_ with any fixes.

Rebuild boundary conditions::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

----

**21 May 2018**

Get James' to generate open bcs and try them out.
Also created a new domain_cfg.nc file.
From /work/jdha/jelt copy to $INPUTS::

  cd /work/jdha/jelt/
  rsync -uvat  SEAsia_full_bdyT_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bt_bdyT_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdyU_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdyV_y1979m11.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  SEAsia_full_bdytide*nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  coordinates.bdy.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  rsync -uvat  domain_cfg_jdha.nc jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.

Also has a bdy_mask.nc and a domain_cfg.nc that were the same as mine.

Edit the namlist for accomdate new files::

  cd $EXP/../EXP_openbcs
  vi namelist_cfg

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
        !                   !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
        cn_domcfg = "domain_cfg_jdha"         ! domain configuration filename

  ...

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_full_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_full_bdyU',  24        , 'vobtcrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_full_bdyV',  24        , 'vobtcrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_full_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_full_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_full_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_full_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''


 !-----------------------------------------------------------------------
 &nambdy_tide   !  tidal forcing at open boundaries
 !-----------------------------------------------------------------------
    filtide      = 'bdydta/SEAsia_full_bdytide_'         !  file name root of tidal forcing files





Catch up with James' edits
++++++++++++++++++++++++++
*4 June 2018*

Add Git repo of MY_SRC and do a "sparse checkout" of the MY_SRC folder::

  mv $CDIR/$CONFIG/MY_SRC $CDIR/$CONFIG/MY_SRC/MY_SRC_safe
  cd $WORK/$USER/$CONFIG
  mkdir git_repo
  cd git_repo
  git init
  git remote add -f origin https://github.com/NOC-MSM/NEMO_cfgs.git  # Creates and empty repo without checking out the features
  git config core.sparseCheckout true

Specify the folders I want::

  echo $CONFIG/MY_SRC >> .git/info/sparse-checkout

Last but not least, update your empty repo with the state from the remote. (And
set the upstream remote master)::

  git pull origin master
  git push --set-upstream origin master

Symbolic link to where MY_SRC actually sits::

  ln -s $WORK/$USER/$CONFIG/git_repo/$CONFIG/MY_SRC $CDIR/$CONFIG/MY_SRC


.. note : the above git integration needs to be abstracted to the setup phase.

Rebuild opa (I removed the tidal analysis key from the cpp file, as we aren't
using it yet) At this time also removed Nico's harmonic analysis changes::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Link in executable::

  cd $EXP/../EXP_openbcs
  rm opa
  ln -s ../BLD/bin/nemo.exe opa

Edit namelist_cfg for boundary conditions. Use full velocities with FRS for tracers
and specified for 3d velcities::


  vi namelist_cfg
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
      nn_dyn2d_dta   =  3                  !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
                                            !  = 2, use tidal harmonic forcing data from files
                                            !  = 3, use external data AND tidal harmonic forcing
      cn_dyn3d      =  'specified'               !
      nn_dyn3d_dta  =  1                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_tra        =  'frs'               !
      nn_tra_dta    =  1                    !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      cn_ice_lim      =  'none'             !
      nn_ice_lim_dta  =  0                  !  = 0, bdy data are equal to the initial state
                                            !  = 1, bdy data are read in 'bdydata   .nc' files
      rn_ice_tem      = 270.                !  lim3 only: arbitrary temperature of incoming sea ice
      rn_ice_sal      = 10.                 !  lim3 only:      --   salinity           --
      rn_ice_age      = 30.                 !  lim3 only:      --   age                --

      ln_tra_dmp    =.false.                !  open boudaries conditions for tracers
      ln_dyn3d_dmp  =.true.                !  open boundary condition for baroclinic velocities
      rn_time_dmp   =  1.                   ! Damping time scale in days
      rn_time_dmp_out =  150.                 ! Outflow damping time scale
      nn_rimwidth   = 9                    !  width of the relaxation zone
      ln_vol        = .false.               !  total volume correction (see nn_volctl parameter)
      nn_volctl     = 1                     !  = 0, the total water flux across open boundaries is zero
      nb_jpk_bdy    =  75                   ! number of levels in the bdy data (set < 0 if consistent with planned run)
  /

Presently need to use James' forcing files::

  ...
  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_full_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_full2_bdyU',  24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_full2_bdyV',  24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_full2_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_full2_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_full_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_full_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
  ! for lim2
  !   bn_frld    = 'amm12_bdyT_ice',         24        , 'ileadfra',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_hicif   = 'amm12_bdyT_ice',         24        , 'iicethic',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_hsnif   = 'amm12_bdyT_ice',         24        , 'isnowthi',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  ! for lim3
  !   bn_a_i     = 'amm12_bdyT_ice',         24        , 'ileadfra',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_ht_i    = 'amm12_bdyT_ice',         24        , 'iicethic',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''
  !   bn_ht_s    = 'amm12_bdyT_ice',         24        , 'isnowthi',    .true.   , .false. ,  'daily'  ,    ''    ,   ''     ,     ''

     cn_dir      = 'bdydta/' !  root directory for the location of the bulk files
     ln_full_vel = .true.   !
  /

Submit for 3 hours (nt=7200; 30 days) and see what happens::

  qsub runscript
  5388562.sdb

*11 June 18* Got a lead from Andrew C about XIOS2. Sumbmitted and waiting
**PENDING**
It should output a sub domain...
Hmm it doesn't seem to output anything.
Perhaps start again with fresh XML files.


---

*19 June 2018*

Check that my PyNEMO implementation matches James' output.
::

  livljobs6:
  cd $INPUTS
  vi namelist.bdy
  ...
  ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
  ln_tra         = .true.               !  boundary conditions for T and S
  ...
  ln_tide = .false.


Rebuild PyNEMO for 3d fields (branch ORCA0083).
run time 1729.83568501 = 28 mins

Output::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc


The rimwidth doesn't get picked up for the u2d and u3d files. Though it does
 for the T,S. Try again using the master branch.
 Execution Time: 1601.68916202::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc

Not sure about the size of the velocity arrays. They don't use rimwidth != 1 ..


*(3 Oct 18)*
Rebuild PyNEMO for 3d fields (branch ORCA0083). And execute (namelist.bdy as above)
Build in environment ```nrct_env2``

Try building tides separately ``ln_tides=.False.``. Otherwise problem with the rimwidth variable.
Had to edit ``inputs_dst.ncml``  to remove e3t_0 and e3w_0 renaming.

(30mins)::

  pynemo -s namelist.bdy

Generates::

  coordinates.bdy.nc
  SEAsia_bdyT_y1979m11.nc
  SEAsia_bt_bdyT_y1979m11.nc
  SEAsia_bdyU_y1979m11.nc
  SEAsia_bdyV_y1979m11.nc



.. also had to regenerate the tidal files and copy.
 Caution the coordinates.bdy.nc for tides only should not be used as one needs
 to set rimwidth=1 to get it to work

Copy files from SAN to ARCHER::

  livljobs4
  cd $INPUTS
  for file in  SEAsia_b*nc; do rsync -uvt $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done
  rsync -uvt coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.

Also copy coordinates.bdy.nc to EXP directory

Move to ARCHER::

  cd $EXP/../EXP_openbcs

Remove links to *_full* in forcing files

vi namelist_cfg::

  !-----------------------------------------------------------------------
  &nambdy_dta    !  open boundaries - external data
  !-----------------------------------------------------------------------
  !              !  file name      ! frequency (hours) ! variable  ! time interp.!  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                 !  (if <0  months)  !   name    !  (logical)  !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
     bn_ssh      = 'SEAsia_bt_bdyT', 24      , 'sossheig',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u2d      = 'SEAsia_bdyU',  24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v2d      = 'SEAsia_bdyV',  24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_u3d      = 'SEAsia_bdyU'   24        , 'vozocrtx',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_v3d      = 'SEAsia_bdyV'   24        , 'vomecrty',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_tem      = 'SEAsia_bdyT'   24        , 'votemper',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''
     bn_sal      = 'SEAsia_bdyT'   24        , 'vosaline',    .true.   , .false. ,  'monthly'  ,    ''    ,   ''     ,     ''


     ...

  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files

This means that NEMO will take the 3d velocities and construct its own depth averaged
and perturbation components.

Restore missing namelist_cfg.nc file.
Rsync and make a symbolic link from INPUTS::

  livljobs6 INPUTS $ rsync -uvt domain_cfg.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/.
  ln -s /work/n01/n01/jelt/SEAsia/INPUTS/domain_cfg.nc domain_cfg.nc

PENDING

Submit (nt=240) and see what happens::

 qsub runscript_short
 5649089.sdb


**The above needs tidying up. But it seems to work**


The full ocean: Initial conditions, tides, rivers, open boundaries - no met
===========================================================================

directory: EXP_fullocean

EXEC: ``nemo.exe``

::
  mkdir EXP_fullocean

Copy files from EXP_tide_initcd to EXP_fullocean **ADDING IN -K keeps symlinks. I THINK**::

  rsync -aPvt --exclude=*restart*nc --exclude=*_?d_*grid_?.nc EXP_tide_initcd/* EXP_fullocean/

Need to add in:

* open boundary conditions
* rivers

Get the executable in place (I think that the standard nemo.exe is OK)::

  cd EXP_fullocean
  rm opa
  ln -s ../BLD/bin/nemo.exe opa


Add rivers
++++++++++

::

  cp ../EXP_openbcs/namelist_cfg namelist_cfg
  vi namelist_cfg

  ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)

  sn_rnf      = 'SEAsia_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
  sn_cnf      = 'SEAsia_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''

  cn_dir      = 'bdydta/'      !  root directory for the location of the runoff files
  ln_rnf_mouth= .false.    !  specific treatment at rivers mouths

Run on the standard queue::

  jpni        =  12    ! standardqueue:12 ! shortqueue:11      !  jpni   number of processors following i (set automatically if < 1)
  jpnj        =  8     !  jpnj   number of processors following j (set automatically if < 1)
  jpnij       =  92    ! standardqueue:92 ! shortqueue:88 !  jpnij  number of local domains (set automatically if < 1)

Submit::

  qsub runscript

NB I don't have Nico's tide mods in. They need to be restored.
30 days ran in 1h44.


Add in FES tides. Tidal Potential and harmonic analyisis
========================================================

Following `<FES2014_Nemo.rst>`_.
Changes to use FES2014 tides::

  cp /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_FES14.h90 $CDIR/$CONFIG/MY_SRC/.

Compile with flag *key_FES14_tides*

Long period tidal potential and variable Love number adjustments (latter including a new namelist variable)::

  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/sbctide.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_mod.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tideini.F90 $CDIR/$CONFIG/MY_SRC/.

For restarts, 3d and fast harmonic analysis::

  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm_fast.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step_oce.F90 $CDIR/$CONFIG/MY_SRC/.

with new compilation key *key_diaharm_fast*

The resulting compilation file is now::

  vi $CDIR/$CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                key_mpp_mpi       \
                key_iomput        \
                key_nosignedzero  \
                key_FES14_tides   \
                key_diaharm_fast


Namelist changes::

  !-----------------------------------------------------------------------
  &nam_diaharm_fast   !   Harmonic analysis of tidal constituents               ("key_diaharm_fast")
  !-----------------------------------------------------------------------
  ln_diaharm_store = .true.
  ln_diaharm_compute = .true.
  ln_diaharm_read_restart = .false.
  ln_ana_ssh   = .true.
  ln_ana_uvbar = .false.
  ln_ana_bfric = .false.
  ln_ana_rho  = .false.
  ln_ana_uv3d = .false.
  ln_ana_w3d  = .false.
  tname(1) = 'O1',
  tname(2) = 'M2',
  /
  ...
  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
     ln_tide     = .true.
     ln_tide_pot = .true.    !  use tidal potential forcing
     ln_tide_ramp= .true.   !
     rdttideramp =    1.     !
     dn_love_number = 0.69!   clname(1)   =  'M2'   !  name of constituent - all tidal components must be set in namelist_cfg
   ...

Changes to the XML file: ``file_def_nemo.xml`` and ``field_def_nemo-opa.xml`` which I already had

Move the working executable *BLD/bin/nemo.exe --> nemo_8Oct18.exe*. Build::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

Move the executable to something sensible::

  mv nemo.exe nemo_FES14-tides_diaharm-fast.exe


Generate FES tides boundary FILES
=================================

Rebuild PyNEMO for 3d fields (branch ORCA0083). And execute (namelist.bdy as above)
Build in environment ```nrct_env2`` (wish I'd kept it as nrct_env...)

On branch ``Generalise-tide-input`` . Then manually set FES as data source in ``Python/pynemo/tide/nemo_bdy_tide3.py``.
Compile pynemo etc as in `<install_nrct.rst>`_

In ``namelist.bdy`` set rimwidth=1 and turn off all the open bcs stuff::

  ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
  ln_tra         = .false.               !  boundary conditions for T and S

Just use::

   ln_tide = T

Try building tides separately ``ln_tides=.False.``. Otherwise problem with the rimwidth variable.
Had to edit ``inputs_dst.ncml``  to remove e3t_0 and e3w_0 renaming.

(37mins)::

  pynemo -s namelist.bdy

Generates::

  coordinates.bdy.nc
  SEAsia_bdytide*nc

Leave the coordinates.bdy.nc (want to use the rimwidth=9 one) and copy to ARCHER::

  livljobs4
  cd $INPUTS
  for file in  SEAsia_bdytide*nc; do rsync -uvt $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done



Submitted run to check if it works. PENDING.
Only runs for 15 mins.

Try running EXP_openbcs for 2hours... This will read FES tidal boundary files...

Link in the appropriate executable.
ln -s BLD/bin/nemo_8Oct18.exe $EXP/../EXP_openbcs/opa

PENDING. 5661999.sdb. (Only ran for 13 mins with Nico's tides in the executable).
COMPLETED



Also try EXP_fullocean without Nico's executable. (Namelist options might break it)
cd EXP_fullocean
ln -s BLD/bin/nemo_8Oct18.exe $EXP/../EXP_fullocean/opa

qsub runscript
Didn't work. Namelist issues.
Recompile and resubmit for 2hrs

Blows up after 4 days (15mins) with exceesive velocity through the sumba straits.
Tried rdttideramp = 5. (formerly 1)

**PENDING**




Backup to repo key files
========================

::

  cd ~/GitLab/NEMO-RELOC/docs/source
  # DOMANcfg namelist_cfg for domain_cfg.nc (for hybrid z-s coordinates)
  rsync -utv jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/namelist_cfg SEAsia_DOMAINcfg_namelist_cfg

  # EXP namelist_cfg
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP00/namelist_cfg SEAsia_EXP_namelist_cfg

  # PyNEMO namelist.bdy
  #rsync -utv jelt@livljobs4:/work/jelt/NEMO/SEAsia/INPUTS/namelist.bdy SEAsia_namelist.bdy

  # Python quick plot of SSH in the output.abort.nc file
  rsync -uvt jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP00/quickplotNEMO.py quickplotNEMO.py

Add PyNEMO namelist.bdy to git repo::

  cd $INPUTS/../..

  git init .
  git remote add origin git@github.com:NOC-MSM/NEMO_cfgs.git
  git config core.sparsecheckout true
  echo "SEAsia/INPUTS/*" >> .git/info/sparse-checkout

  git add SEAsia/INPUTS/namelist.bdy
  git commit -m 'add namelist.bdy'
---
