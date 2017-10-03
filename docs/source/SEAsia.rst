=======================================
Setting up a SE Asia NEMO configuration
=======================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html

Issues that arose
=================

* Fix the sosie make.macro command

Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Follow PyNEMO recipe for LBay on ARCHER (not complete because PyNEMO java issue): ``http://nemo-reloc.readthedocs.io/en/latest/LBay.html``.
Follow PyNEMO recipe for LBay on ARCHER/Livljobs4: ``http://nemo-reloc.readthedocs.io/en/latest/LBay_livljobs4.html``.


Summary / Plan
==============

#. Update to ORCHESTRA code base

#. Extend domain in a new experiments: 1/12 degree




Recipe Notes
============

Define working directory, and other useful shortcuts::

  export CONFIG=ACCORD
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export CDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_SEAsia

#Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel


Build NEMO ORCHESTRA branch @ r8395::

  mkdir $WDIR
  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/NERC/dev_r6998_ORCHESTRA@8395

Use Dave's XIOS file (see ``%XIOS_HOME``)::

  cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Copy Maria's MY_SRC::

  rsync -uartv /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/CONFIG/TIDE/MY_SRC/ $CDIR/$CONFIG/MY_SRC

Copy Maria's cpp flags (without the lim flag)::

  vi $CDIR/$CONFIG/cpp_ACCORD.fcm

   bld::tool::fppkeys key_mpp_mpi          \
                      key_bdy              \
                      key_tide            \
                      key_zdftke           \
                      key_netcdf4          \
                      key_iomput           \
                      key_nosignedzero    \
                      key_trabbl        \
                      key_zdfddm        \
                      key_diaharm      \
                      key_xios2


Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**It compiles!!**

Make an EXP directory::

  mkdir $EXP


----


Build some NEMO tools (on ARCHER)
---------------------------------

For the old build ``dev_r4621_NOC4_BDY_VERT_INTERP/`` we applied patches.
For the new build ``dev_r6998_ORCHESTRA``, we do not. For some reason GRIDGEN doesn’t like INTEL::

  cd $TDIR
  ./maketools -n WEIGHTS -m XC_ARCHER_INTEL
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module swap PrgEnv-intel PrgEnv-cray
  module load cray-netcdf cray-hdf5
  ./maketools -n GRIDGEN -m XC_ARCHER

  module swap PrgEnv-cray PrgEnv-intel

----

*(3 Oct 2017)*


5. Generate mesh and mask files for open boundary conditions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Run the model to generate the mesh and mask files. Copy
 the input files from Liverpool to ARCHER::

  ssh livljobs4

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
  export START_FILES=$WDIR/START_FILES # generic stuff for making more stuff. Mostly code.
  export INPUTS=$WDIR/INPUTS         # config specific stuff that gets made and is for running NEMO

  # Copy into the $EXP directory
  for file in bathy_meter.nc coordinates.nc; do  scp $INPUTS/$file  jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/ACCORD/EXP_SEAsia/$file; done
  for file in runscript_archer namelist_cfg namelist_ref iodef.xml; do  scp $START_FILES/$file  jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/ACCORD/EXP_SEAsia/$file; done

Return to ARCHER. Make sure the executables are in the EXP dir.
Using XIOS from Dave::

  ln -s /work/n01/n01/munday/XIOS/bin/xios_server.exe $EXP/.
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa

Also link in the extra XML files::

  ln -s $EXP/../../SHARED/field_def.xml $EXP/.
  ln -s $EXP/../../SHARED/domain_def.xml $EXP/.



Edit the namelist files for this configuration::

  cd $EXP
  ncdump -h coordinates.nc
  x = 683 ;
  y = 553 ;

  vi namelist.cfg
  ...
  cn_exp      =   "SEAsia"  !  experience name
  ...
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     cp_cfg      =  "seasia"                !  name of the configuration
     jp_cfg      =     012               !  resolution of the configuration
     jpidta      =     683               !  1st lateral dimension ( >= jpi )
     jpjdta      =     553               !  2nd    "         "    ( >= jpj )
     jpkdta      =      51               !  number of levels      ( >= jpk )
     jpiglo      =     683               !  1st dimension of global domain --> i =jpidta
     jpjglo      =     553               !  2nd    -                  -    --> j  =jpjdta

**ACTION: There are further edits to be made for when the model is actually run**
**E.g. other filename instances of LBay**

**NOT SURE WHAT jp_cfg does. It might be passive? I changed it to 012, representing 1/12**


Edit the runscript to include modules and the Account name (n01-NOCL)::

  vi runscript_archer

  #!/bin/bash
  #PBS -N SEAsia
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel
  ...

Submit::

  qsub -q short runscript_archer


*IT WORKS IF I USE THE LBAY EXECUTABLE!!**
----

*(6 March 2017)*

If that works, we then need to rebuild the mesh and mask files in to single files for the next step::

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_zgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_hgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mask 96
  mv mesh_zgr.nc mesh_hgr.nc mask.nc $WDIR/INPUTS
  rm mesh_* mask_* LBay_0000*
  cd $WDIR/INPUTS


THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO


----







Old notes
---------






----

*(27 Sept 2017)*

Build the new SE Asia configuration at 1/12 degree, R12
=======================================================

Generate new coordinates file
=============================

Inspect TPXO harmonic amplitudes to find a good cut off location for boundaries:

cd /work/jelt/tpxo7.2
ferret
go  plot_SEAsia_harmonics.jnl

... note::
  ! plot_SEAsia_harmonics.jnl
  ! Plot tpxo harmonics for the SE Asia region.
  ! Want to build a NEMO config without significant amphidromes on the boundary

  use h_tpxo7.2.nc

  set win 1
  set viewport ul
  shade/k=1/j=300:700/i=250:500/levels=(0,1,0.1)/title="M2" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=2/j=300:700/i=250:500/levels=(0,1,0.1)/title="S2" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=3/j=300:700/i=250:500/levels=(0,1,0.1)/title="N2" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=4/j=300:700/i=250:500/levels=(0,1,0.1)/title="K2" HA, lon_z, lat_z; go fland

  set win 2
  set viewport ul
  shade/k=5/j=300:700/i=250:500/levels=(0,1,0.1)/title="K1" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=6/j=300:700/i=250:500/levels=(0,1,0.1)/title="O1" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=7/j=300:700/i=250:500/levels=(0,1,0.1)/title="P1" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=8/j=300:700/i=250:500/levels=(0,1,0.1)/title="Q1" HA, lon_z, lat_z; go fland


Conclusion. Plot the proposed domain::

  $livljobs2$ scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.

  ferret
  use coordinates_ORCA_R12.nc
  set win 1; shade/X=50:730/Y=1250:1800 E2T, nav_lon, nav_lat ; go fland
  set win 2; set viewport upper; shade/i=50:730/j=1250:1800 NAV_LAT
  set win 2; set viewport lower; shade/i=50:730/j=1250:1800 NAV_LON




---

----
