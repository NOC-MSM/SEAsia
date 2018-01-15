=============================================
Setting up a SW Pacific NEMO v4 configuration
=============================================

Machines: livljobs4, ARCHER

URL:: *to add*
*Actually I've taken this off readthedocs as it is internal. However writing in*
*ReStrucTed text means it is trivial to deploy anywhere with whatever style*

* Build notes with:: ~/GitLab/NEMO-RELOC/docs$ make html

Builds a SW Pacific regional tide-only model using GEBCO bathymetry, TPXO tidal
boundaries.

Build on a combination of livljobs4 and ARCHER.

Uses a pre-release of NEMO v4 (@r8395)

The summary procedure:
#. ARCHER: Get code. Build tools. Generate coordinates, bathymetry, domain_cfg.nc
#. ARCHER: Generate initial conditions and atmospheric forcings
#. LIVLJOBS4: Generate boundary conditions with NRCT/PyNEMO
#. ARCHER: Run simulation

It is a tide-only run.

Issues that arose
=================

* I had real trouble getting this to be stable is s-ccordinates, thought it worked
for the SEAsia doman. S-coordinates affects the domain_cfg.nc generation and PyNEMO
process. So I have made "protected" copies of the ``ARCHER:$TIDR/DOMAINcfg/SCO``
``LIVLJOBS:$INPUTS/SCO`` directories, and tried again with z-coords (following \
James' examples). When I got z-coords working I did not go back to see if s-coords
where also fixed.

* This configuration uses 31 levels. At one point I used 5 levels but it didn't
work for, perhaps, other reasons. I haven't revisited the number of levels.

* The rimwidth variable is not used for barotropic runs, though it must be consistent between
PyNEMO and NEMO namelists. The rimwdith determines the relaxation zone for U,V (+T,S?) fields.

* There was a mask fix missing from the NEMO trunk that meant the bdy mask was not imposed
back on the global mask. This led to near boundary instabilities.

* PyNEMO does not expect users to only want tidal output and error trapping is in
a nascent stage. With the current namelist.bdy (PyNEMO namelist) it crashes after
generating the tidal boundary conditions while generating the 3D boundary conditions
that I don't care about. I did not bother to resolve this.

* Trouble shooting issues are covered in `trouble_shooting.rst`_

.. note: PyNEMO is interchangabably called NRCT (NEMO Relocatable Configuration Tool)


----

Recipe Notes
============

In the following I build most stuff on ARCHER but the PyNEMO bits are done on livljobs4.
(There was a problem with some java bits working on ARCHER)
Starting on ARCHER::

  ssh login.archer.ac.uk

It is quite convenient to define a temporary file with all the path names in.
The first line defines the configuration name (assuming bash)::

cat > ~/temporary_path_names_for_NEMO_build << EOL
export CONFIG=LBay180
export WORK=/work/n01/n01
export WDIR=\$WORK/\$USER/\$CONFIG
export INPUTS=\$WDIR/INPUTS
export START_FILES=\$WDIR/START_FILES
export CDIR=\$WDIR/trunk_NEMOGCM_r8395/CONFIG
export TDIR=\$WDIR/trunk_NEMOGCM_r8395/TOOLS
export EXP=\$CDIR/\$CONFIG/EXP00

module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf-hdf5parallel cray-hdf5-parallel
EOL

Execute path settings::

  . ~/temporary_path_names_for_NEMO_build

---

Collect essential files
=======================

Note you might have to mkdir the odd directory or two...::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES

At this point you need to beg, steal or borrow to fill your $START_FILES. A good place
to look would be ``$WORK/jelt/$CONFIG/START_FILES``. E.g.::

  cp $WORK/jelt/LBay/START_FILES/dommsk.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/bdyini.F90 $START_FILES/.
  cp $WORK/jelt/LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WORK/jelt/LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.
  cp $WORK/jelt/SEAsia/START_FILES/usrdef_istate.F90 $START_FILES/.
  cp $WORK/jelt/SEAsia/START_FILES/usrdef_sbc.F90    $START_FILES/.

Checkout and build NEMO (ORCHESTRA) trunk @ r8395 `build_opa_orchestra.html`_.
Or just build (if it is already downloaded). Note here we use user defined
 functions for the initial state (constant T and S) and surface forcing (zero forcing)::

  cd $CDIR
  cp $START_FILES/usrdef_istate.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/usrdef_sbc.F90    $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dommsk.F90 $CDIR/$CONFIG/MY_SRC/.
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**Note** Make sure to copy across any inital state files into MY_SRC before building, this is important as the SWPacific requires a constant temp and salinity to be set at the start of the simulation. I used the generic recipe to build NEMO and missed these files.

---
**NOTE** You need to do this before making nemo
Checkout and build XIOS2 @ r1080 `build_XIOS2.html`_::

Or just link XIOS executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---

Note 
The NEMO executable is called ‘opa’ and found in your $EXP directory. 
This is automatically symbolically linked to 'nemo.exe' found in $CDIR/$CONFIG/BLD/bin/


Build TOOLS
===========

To generate domain coords and rebuild tools we first need
to compile some of the NEMO TOOLS.

.. note: These are compiled with XIOS2. However DOMAINcfg has to be compiled
  with XIOS1. There is a README in the $TDIR/DOMAINcfg on what to do.

First build DOMAINcfg (which is relatively new and in NEMOv4). Use my XIOS1 file
(see userid and path in variable ``%XIOS_HOME``). Copy from ARCH *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $CDIR/../ARCH/.
  cd $TDIR

  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO

For the generation of bathymetry and met forcing weights files we need to patch
the code (to allow direct passing of arguments. NB this code has not been
updated in 7 years.)::

  cd $TDIR/WEIGHTS/src
  patch -b < $START_FILES/scripinterp_mod.patch
  patch -b < $START_FILES/scripinterp.patch
  patch -b < $START_FILES/scrip.patch
  patch -b < $START_FILES/scripshape.patch
  patch -b < $START_FILES/scripgrid.patch

  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n WEIGHTS



.. note :
  Could skip 1 and 2 - COPY bathy_meter.nc and coordinates.nc from Tom with a
  resolution refinement of x3
    /work/thopri/NEMO/SWPacific_ver3.6/INPUTS/bathy_meter.nc
    /work/thopri/NEMO/SWPacific_ver3.6/INPUTS/coordinates.nc

    Otherwise I'd follow the following instuctions.



1. Generate new coordinates file
++++++++++++++++++++++++++++++++

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements. Originally we tried to interpolate to a finer grid with
``rho=rhot=3``. Now keeping resolution at R12 incase it sorts out a stability issue.

.. warning:
  Using the GRIDGEN/create_coordinates.exe tool runs into a problem for zoom factor
  >1, since the horizontal spacing metric e.g. e[12]t always match
  the parent grid. I think that this is a bug. The agrif version works.

Edit namelist::

  cd $TDIR/NESTING
  vi namelist.input

  &input_output
      iom_activated = true
  /
  &coarse_grid_files
      parent_coordinate_file = 'coordinates_ORCA_R12.nc'
  /
  &bathymetry
  /
  &nesting
      imin = 865
      imax = 1405
      jmin = 1116
      jmax = 1494
      rho  = 1
      rhot = 1
      bathy_update = false
  /
  &vertical_grid
  /
  &partial_cells
  /
  &nemo_coarse_grid
  /
  &forcing_files
  /
  &interp
  /
  &restart
  /
  &restart_trc
  /

Build and execute agrif version of create_coordinates.exe.
See `build_and_create_coordinates.rst`_

This creates a new coordinatesfile with contents, which is now copied to
INPUTS::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc

Now we need to generate a bathymetry on this new grid.



2. Generate bathymetry file
+++++++++++++++++++++++++++

thopri downloaded, by parts, and merged 1 minute GEBCO data for (-30N :0N , -170E : 145E ).
Method in ``/work/thopri/NEMO/SWPacific/START_FILES/gebco_lon_convertor.py``
Copy this to ARCHER::

  livljobs4$ rsync -utv /work/thopri/NEMO/SWPacific_ver3.6/START_FILES/GRIDONE_2D_140_-35.0_-165.0_5.0.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/SWPacific/INPUTS/.

Copy namelist for reshaping GEBCO data::

  rsync -utv $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.

Edit namelist to point to correct input file. Edit lat and lon variable names to
 make sure they match the nc file content (used e.g.
``ncdump -h GRIDONE_2D_140_-35.0_-165.0_5.0.nc`` to get input
variable names)::

  vi $INPUTS/namelist_reshape_bilin_gebco
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


Do some things to 1) flatten out land elevations, 2) make depths positive.
Have to swap around with the modules to get nco working *(James
noted a problem with the default nco module)*::

  cd $INPUTS

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5

  module load nco/4.5.0
  ncap2 -s 'where(elevation > 0) elevation=0' GRIDONE_2D_140_-35.0_-165.0_5.0.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc

Restore the original parallel modules::

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

Skip this for a tide-only run. Using user defined constant T and S.

For constant T and S use the user defined functions in ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``.


4. Generate a domain configuration file
=======================================

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg``
directory along with all the inputs files that would have previously been needed
get v3.6 running. The reason being that all the non-time stepping stuff, like
grid generating, has been abstracted from the core OPA code and is now done as
a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory::

    rsync -uvt $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
    rsync -uvt $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template namelist_cfg with only the essenetial domain building stuff.
Get the indices from ``ncdump -h coordinates.nc``.

When the increase in resolution was x3, the size of the coordinate were::

  jpidta      =    1624   !  1st lateral dimension ( >= jpi )
  jpjdta      =    1138   !  2nd    "         "    ( >= jpj )

When the increase in resolution is x1, the size of the coordinates were::

  jpidta      =    544   !  1st lateral dimension ( >= jpi )
  jpjdta      =    382   !  2nd    "         "    ( >= jpj )


Somewhat arbitrarily I am going to use **31** ``(jpkdta=31; rn_jpk=31)`` levels.
(Earlier I used 5 levels. This config gave me trouble. I've not revisited whether
the number of levels was a problem).


s-coordinates `SWPacific_DOMAINcfg_namelist_cfg`_
-------------

skip this and jump to the z-coordinates heading

.. note: At one point I tried to use s-coords but it wasn't working. Perhaps it
 would work now that the boundary problems are fixed. I did not investigate it
 though. However I leave these s-coordinate namelist options here because they
 will be useful when I need s-coordinates to work in v4.

(didn't work, though pehaps for reasons nothing to do with the coordinates...)

 (If I tried 5 levels then the function to compute the depth range breaks).
 Added in s-coordinate settings from AMM60. (There is a function that tapers dz
 near the equation. This is activated in ``domzgr.F90`` and not by a logical flag)::

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
     ln_e3_dep   = .false.   ! =T : e3=dk[depth] in discret sens.
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !
     cp_cfg      =  "orca"   !  name of the configuration
     jp_cfg      =      12   !  resolution of the configuration
     jpidta      =     544   !  1st lateral dimension ( >= jpi )
     jpjdta      =     382   !  2nd    "         "    ( >= jpj )
     jpkdta      =      31   !  number of levels      ( >= jpk )
     jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     382   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
    ln_zco      = .false.   !  z-coordinate - full    steps
    ln_zps      = .false.   !  z-coordinate - partial steps
    ln_sco      = .true.   !  s- or hybrid z-s-coordinate
    ln_isfcav   = .false.   !  ice shelf cavity
    ln_linssh   = .false.   !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
    ln_s_sh94   = .false.    !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
    ln_s_sf12   = .true.   !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
    ln_sigcrit  = .true.    !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                            !  stretching coefficients for all functions
    rn_jpk      =   31       ! Number of S levels
    !ln_eq_taper = .false.   !  Tapering of S coords near equator
    cn_coord_hgr = 'coordinates.nc'  ! File containing gphit (latitude) coordinate for use if ln_eq_taper=.true.
    rn_sbot_min =   10.0    !  minimum depth of s-bottom surface (>0) (m)
    rn_sbot_max = 7000.0    !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
    rn_hc       =   50.0    !  critical depth for transition to stretched coordinates
                         !!!!!!!  Envelop bathymetry
    rn_rmax     =    0.3    !  maximum cut-off r-value allowed (0<r_max<1)
                         !!!!!!!  SH94 stretching coefficients  (ln_s_sh94 = .true.)
    rn_theta    =    6.0    !  surface control parameter (0<=theta<=20)
    rn_bb       =    0.8    !  stretching with SH94 s-sigma
                         !!!!!!!  SF12 stretching coefficient  (ln_s_sf12 = .true.)
    rn_alpha    =    4.4    !  stretching with SF12 s-sigma
    rn_efold    =    0.0    !  efold length scale for transition to stretched coord
    rn_zs       =    1.0    !  depth of surface grid box
                            !  bottom cell depth (Zb) is a linear function of water depth Zb = H*a + b
    rn_zb_a     =    0.024  !  bathymetry scaling factor for calculating Zb
    rn_zb_b     =   -0.2    !  offset for calculating Zb
                         !!!!!!!! Other stretching (not SH94 or SF12) [also uses rn_theta above]
    rn_thetb    =    1.0    !  bottom control parameter  (0<=thetb<= 1)
  /


z-coordinates
-------------

::
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
     ln_e3_dep   = .false.   ! =T : e3=dk[depth] in discret sens.
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !
     cp_cfg      =  "orca"   !  name of the configuration
     jp_cfg      =      12   !  resolution of the configuration
     jpidta      =     544   !  1st lateral dimension ( >= jpi )
     jpjdta      =     382   !  2nd    "         "    ( >= jpj )
     jpkdta      =      31   !  number of levels      ( >= jpk )
     jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     382   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
    ln_zco      = .false.   !  z-coordinate - full    steps
    ln_zps      = .true.   !  z-coordinate - partial steps
    ln_sco      = .false.   !  s- or hybrid z-s-coordinate
    ln_isfcav   = .false.   !  ice shelf cavity
    ln_linssh   = .false.   !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     jphgr_msh   =       0               !  type of horizontal mesh
     ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
     ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
     ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
     ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
     ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
     ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
     ppsur       =   -4762.96143546300   !  ORCA r4, r2 and r05 coefficients
     ppa0        =     255.58049070440   ! (default coefficients)
     ppa1        =     245.58132232490   !
     ppkth       =      21.43336197938   !
     ppacr       =       3.0             !
     ppdzmin     =  999999.              !  Minimum vertical spacing
     pphmax      =  999999.              !  Maximum depth
     ldbletanh   =  .FALSE.              !  Use/do not use double tanf function for vertical coordinates
     ppa2        =  999999.              !  Double tanh function parameters
     ppkth2      =  999999.              !
     ppacr2      =  999999.              !
  /

----



.. note:

  No gdept output in the offical v4 release. Though it was acheived here setting
  ln_e3_dep = F. This is needed for PyNEMO, though could be constructed from e3[tw].

Build a script to run the executable::

  vi $TDIR/DOMAINcfg/rs

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

Put a copy of the namelist_cfg in $INPUTS for safe keeping::

  cp $TDIR/DOMAINcfg/namelist_cfg $INPUTS/namelist_cfg_generateDOMAINcfg

#Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild). Also fix the permission as it seems they are
 not readable by team members otherwise::

  chmod a+rx $TDIR/DOMAINcfg/domain_cfg.nc
  rsync -uvt $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  rsync -uvt $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.




5. Generate weights for atm forcing
+++++++++++++++++++++++++++++++++++

Skip this for a tide-only simulation

---



THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO
If all the files are ready to go jump straight to `7. Generate boundary conditions with PyNEMO: Run PyNEMO`_

----

Follow SEAsia notes `SEAsia_archer_livljobs4.rst`_ for PyNEMO bit. Changing
instances of `SEAsia` for `SWPacific`.
Note a lot of the files are actually used for a tide-only simulation but it is not
trivial to prevent PyNEMO looking for them. *Update 16Nov17: I think the following is now stand alone)*

----

Set up directory structure and files in livljobs4
=================================================

Define paths::

cat > ~/temporary_path_names_for_NEMO_build << EOL
export CONFIG=SWPacific
export WORK=/work
export WDIR=\$WORK/$USER/NEMO/\$CONFIG
export INPUTS=\$WDIR/INPUTS
export START_FILES=\$WDIR/START_FILES
#export CDIR=\$WDIR/trunk_NEMOGCM_r8395/CONFIG
#export TDIR=\$WDIR/trunk_NEMOGCM_r8395/TOOLS
#export EXP=\$CDIR/\$CONFIG/EXP00

EOL

Execute path settings::

  . ~/temporary_path_names_for_NEMO_build

Synchronise with key files from ARCHER::

  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/bathy_meter.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.nc $INPUTS/.
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/domain_cfg.nc  $INPUTS/.


Statement about external forcing and files
==========================================

Use ORCA 1/12 data via a thredds server. Also use ORCA 1/12 mesh and mask files.
(Originally tried to get them via a thredds server but gave up trying to figure
out how to write the NCML file).
 Copy necessary files into INPUTS. (Be careful of symbolic links in PyNEMO).::

   ls -lh $INPUTS/bathy_meter.nc
   ls -lh $INPUTS/coordinates.nc
   ls -lh $INPUTS/domain_cfg.nc

Need to generate 6 more files: A ``namelist.bdy`` (archived as `SWPacific_namelist.bdy`_)
which drives PyNEMO and which has two input files: ``inputs_src.ncml``
which points to the remote ORCA data source and ``inputs_dst.ncml`` which
remaps some variable names in the destination files. Then generate ncml files
to get the mesh and mask files (``mask_src.ncml  mesh_hgr_src.ncml  mesh_zgr_src.ncml``).
*NB In the end I copied the source mesh and mask files locally. So obly need to
generate 3 files*

Create ncml file for ORCA12 source data::

  vi  $INPUTS/inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791206d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>


Create NCML file for mapping variables in DESTINATION grid to what PyNEMO expects::

  vi inputs_dst.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf location="file:domain_cfg.nc">
      <ns0:variable name="mbathy" orgName="top_level" />
      <ns0:variable name="gdept" orgName="gdept_0" />
      <ns0:variable name="gdepw" orgName="gdepw_0" />
      <ns0:variable name="e3u" orgName="e3u_0" />
      <ns0:variable name="e3v" orgName="e3v_0" />
      <ns0:variable name="e3t" orgName="e3t_0" />
      <ns0:variable name="e3w" orgName="e3w_0" />
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>

Create a PyNEMO namelist.bdy file. Note this this is tides only and so z-coordinates
are OK as it is only 2D variables. I have archived a working version of this file
here `SWPacific_namelist.bdy`_::

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
     ln_zco      = .true.   !  z-coordinate - full    steps   (T/F)
     ln_zps      = .false.    !  z-coordinate - partial steps   (T/F)
     ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
     rn_hmin     =   -5     !  min depth of the ocean (>0) or
                             !  min number of ocean level (<0)

  !-----------------------------------------------------------------------
  !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
     rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
     rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
                             !  (= ocean depth) (>0) (m)
     ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
     rn_hc       =  50.0    !  critical depth with s-sigma

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
     sn_src_dir = './inputs_src.ncml'       ! src_files/'
     sn_dst_dir = '/work/jelt/NEMO/SWPacific/INPUTS/'
     sn_fn      = 'SWPacific'                 ! prefix for output files
     nn_fv      = -1e20                     !  set fill value for output files
     nn_src_time_adj = 0                                    ! src time adjustment
     sn_dst_metainfo = 'metadata info: jelt'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
  ln_coords_file = .true.               !  =T : produce bdy coordinates files
  cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
  ln_mask_file   = .true.              !  =T : read mask from file
  cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
  ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
  ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
  ln_tra         = .false.               !  boundary conditions for T and S
  ln_ice         = .false.               !  ice boundary condition
  nn_rimwidth    = 1                    !  width of the relaxation zone

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
  nn_year_000     = 1979        !  year start
  nn_year_end     = 1979        !  year end
  nn_month_000    = 11          !  month start (default = 1 is years>1)
  nn_month_end    = 11          !  month end (default = 12 is years>1)
  sn_dst_calendar = 'gregorian' !  output calendar format
  nn_base_year    = 1978        !  base year for time counter
  sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'

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
  rn_mask_max_depth = 7000.0    !  Maximum depth to be ignored for the mask
  rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break


Tried to put the mesh and mask files as remote access via ncml files. Failed, so
downloaded them instead. (Perhaps the mesh and mask variables could all go in
the inputs_src.ncml ?)...

.. note:

    Finally, create separate ncml files for each mesh mask file (because I don't know
    how to join them...) and they are on a thredds server::

      vi mesh_hgr_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mesh_hgr.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>


      vi mesh_zgr_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mesh_zgr.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>


      vi mask_src.ncml

      <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
      <ns0:aggregation type="union">
        <ns0:netcdf>
              <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/domain/mask.nc" />
        </ns0:netcdf>
      </ns0:aggregation>
      </ns0:netcdf>

Finally I am going to create a boundary mask file. This can also be done with
the PyNEMO GUI. The mask variable takes values (-1 mask, 1 wet, 0 land). Get a
template from domain_cfg.nc and then modify as desired around the boundary::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython::

  import netCDF4
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  dset.variables['mask'][0:4,:]  = -1     # Southern boundary
  dset.variables['mask'][-4:-1,:] = -1    # Northern boundary
  dset.variables['mask'][:,-4:-1] = -1    # Eastern boundary
  dset.variables['mask'][:,0] = -1        # Western boundary
  dset.close()

.. note : 8 Nov. Fiddling with the bdy_msk is probably not such a good idea unless
  there are clear issues with where the boundary lies E.g. island or other
  bathymetric or amphidromic features. I now think that the problems I had were
  dues to T,S boundary issues, which are now set and held constant.

Run PyNEMO / NRCT to generate boundary conditions
+++++++++++++++++++++++++++++++++++++++++++++++++

First install PyNEMO `install_nrct`_ if not already done so.

Generate the boundary conditions with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy

.. note : Can use the ``-g`` option if you want the GUI.

.. note : I have a PyNEMO mod to use FES tides instead of TPXO tides for these boundary
  forcing. It is currently a hardwire fix in ``tide/nemo_bdy_tide3.py``

---

(By the time it crashes on some 3D stuff I am not interested in) this has already
 generated what I need::

  ls -1 $INPUTS

  coordinates.bdy.nc
  SWPacific_bdytide_rotT_M2_grid_T.nc
  SWPacific_bdytide_rotT_K2_grid_T.nc
  SWPacific_bdytide_rotT_S2_grid_T.nc
  SWPacific_bdytide_rotT_M2_grid_U.nc
  SWPacific_bdytide_rotT_K2_grid_U.nc
  SWPacific_bdytide_rotT_S2_grid_U.nc
  SWPacific_bdytide_rotT_M2_grid_V.nc
  SWPacific_bdytide_rotT_K2_grid_V.nc
  SWPacific_bdytide_rotT_S2_grid_V.nc

.. note :
    with ``ln_tra = .false.`` it crashed earlier (and looks fixable) but still produces the
    above files. I assume that that is OK too but haven't checked.

Copy the new files back onto ARCHER::

  livljobs4$
  cd $INPUTS
  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.bdy.nc
  for file in $CONFIG*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done


8. Run the configuration ON ARCHER. Turn on the tides
+++++++++++++++++++++++++++++++++++++++++++++++++++++

Get set up::

  ssh archer
  . ~/temporary_path_names_for_NEMO_build

Get important files into EXP directory. Should already have ``domain_cfg.nc``::


  cd $EXP
  rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.nc $EXP/.
  rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.
  rsync -tuv $START_FILES/namelist_cfg $EXP/.

Create symbolic links from EXP directory::

  ln -s $INPUTS $EXP/bdydta

Edit the output to have 1hrly SSH::

 vi file_def_nemo.xml
 ...
 <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
  <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
    <field field_ref="ssh"          name="zos"   />
  </file>
 </file_group>

---

Create a short queue runscript. (Note: PBS -N jobname, PBS -m email)::

  vi runscript

  #!/bin/bash
  #PBS -N SWPacific
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
  ulimit -c unlimited
  ulimit -s unlimited

  rm -f core
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

  exit

Change the notification email to your own address::

  sed -i "s/xxx@noc/$USER@noc/g" runscript

---

Edit ``namelist_cfg`` to make sure it is OK

Synchronise the namelist_cfg with the GitLab repo. (This is also done at the end
of this repo)::

  e.g. rsync -uvt $USER@login.archer.ac.uk:/work/n01/n01/$USER/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/namelist_cfg SWPacific_EXP_namelist_cfg

---

Submit::

 cd $EXP
 qsub -q short runscript

 **WOO HOO THIS WORKS. THE FOLLOWING ARE LIVE NOTES THAT MAY OR MAY NOT BE INTERESTING**



Live notes (Probably best to skip this section)
================================================

.. note : 27 Oct. Switched to z-coords. Need to archive s-coord namelists

Iterative debugging process. Configuration blew up around the boundary
E.g.::

  stp_ctl : the ssh is larger than 10m
  =======
  kt=   410 max ssh:   10.08    , i j:   446   11

  stp_ctl : the ssh is larger than 10m
  =======
  kt=   101 max ssh:   10.32    , i j:   172    4

  stpctl: the zonal velocity is larger than 20 m/s
  ======
  kt=  6374 max abs(U):   21.05    , i j k:   542   57   30

  stpctl: the zonal velocity is larger than 20 m/s
  ======
  kt=  6209 max abs(U):   23.41    , i j k:   532  379   10

Edit the bdy_mask.nc file to blank out (-1) the edges with problems. (I blanked them
out preserving a rectangular grid. This may not be necessary). Rerun PyNEMO to
generate new boundary forcing.

.. note : 8 Nov. Fiddling with the bdy_msk is probably not such a good idea unless
  there are clear issues with where the boundary lies E.g. island or other
  bathymetric or amphidromic features. I now think that the problems I had were
  dues to T,S boundary issues, which are now set and held constant.



In the end I ran with ``rn_rdt=360`` (6 minutes). I also played with the number of
processors used to do the XIOS, but I am not sure that it did anything to speed up
the simulation. (Changed from ``-n 5 -N 5`` to ``-n 20 -N 20``).

I was able to run 20 days in 17 mins on the short queue.

Then optimising the MPI decomposition (93 processors - saving 3 from land). Competed
in 18 mins::

  aprun -b -n 20 -N 20 ./xios_server.exe : -n 93 -N 24 ./opa




MPP decomposition for land suppression
++++++++++++++++++++++++++++++++++++++

`MPP_decomp_lanf_suppression.rst`_

Before doing long runs it is important to optimise MPP decompositoin by invoking
 land supression to save redundant ocean processors.
Resulting decomposition::

   vi namelist_cfg
   ...
   !-----------------------------------------------------------------------
   &nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
   !-----------------------------------------------------------------------
      ...
      jpni        =  12       !  jpni   number of processors following i (set automatically if < 1)
      jpnj        =  8    !  jpnj   number of processors following j (set automatically if < 1)
      jpnij       =  1    !  jpnij  number of local domains (set automatically if < 1)

Inspect ``ocean_output`` to find ``jpnij``. In my simulation ``jpni=8, jpnj=12 --> jpnij = 93``
Update OCEANCORES in runscript (make sure the ``aprun`` statement is as expected too)::

  vi runscript
  ...
  OCEANCORES=93

And submit again.


Restart run
+++++++++++

Rebuild restart file `rebuild_and_inspect_NEMO_output.rst`_
Modify the namelist_cfg file for a restart (also doubled number of steps to 9600).
Using nn_restctl = 2 don't have to be too careful about dates in the namelist_cfg
as it comes from the restart file instead. (Though it seems to throw an error if
you get a mismatch between namelist and restart dates..)
::
  vi namelist_cfg
  ...
  nn_it000    =  14401   !  first time step
  nn_itend    =  19200 ! 10day=14400   !  last  time step (std 5475)
  nn_date0    =  20000302   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
  ln_rstart   = .true.   !  start from rest (F) or from a restart file (T)
  ...
   nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
   !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
  ...
  cn_ocerst_in    = "SWPacific_00014400_restart_tide"   !  suffix of ocean restart name (input)


  !-----------------------------------------------------------------------
  &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
  !-----------------------------------------------------------------------
      nit000_han = 4801         ! First time step used for harmonic analysis
      nitend_han = 14400 ! 7200      ! Last time step used for harmonic analysis
      nstep_han  = 1        ! Time step frequency for harmonic analysis
      tname(1)   = 'M2'      ! Name of tidal constituents
      ...


Double the wall time in ``runscript`` (40mins) and submit::

  qsub runscript

---


Backup to repo key files
========================

When I make changes I backup version of the key files in the repo to keep a track
of what has changed::

  cd ~/GitLab/NEMO-RELOC/docs/source
  # DOMANcfg namelist_cfg for domain_cfg.nc (for s-coordinates)
  rsync -utv $USER@login.archer.ac.uk:/work/n01/n01/$USER/SWPacific/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/namelist_cfg SWPacific_DOMAINcfg_namelist_cfg

  # EXP namelist_cfg (for s-coordinates)
  rsync -uvt $USER@login.archer.ac.uk:/work/n01/n01/$USER/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/namelist_cfg SWPacific_EXP_namelist_cfg

  # PyNEMO namelist.bdy (for s-coordinates)
  rsync -utv $USER@livljobs4:/work/$USER/NEMO/SWPacific/INPUTS/namelist.bdy SWPacific_namelist.bdy

  # Python quick plot of SSH in the output.abort.nc file
  rsync -uvt $USER@login.archer.ac.uk:/work/n01/n01/$USER/SWPacific/trunk_NEMOGCM_r8395/CONFIG/SWPacific/EXP00/quickplotNEMO.py quickplotNEMO.py
