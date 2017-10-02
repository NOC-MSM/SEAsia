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

* start with LBay using ORCHESTRA code base
* Check James' /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/
configuration for suggested namelist options for a lighter tide only run


Or

* Replicate ODA_E-AFRICA config:
    1) copy namelists, execs, forcings. Test.
    2) generate execs. (James: trunk @ r7853 and XIOS2 @ r1080). Test.

* Extend domain in a new experiments: 1/12 degree
    3) generate coords, forcings

* Update to ORCHESTRA code base
    4) Later.


diff EA and LBay cpp keys.
NB Tide only wont require all the T,S forcing fields.


Recipe Notes
============

Define working directory, and other useful shortcuts::

  export CONFIG=ACCORD
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export CDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_EAFRICA

  export JINPUTS=/work/n01/n01/jdha/2017/INPUTS/ODA/E-AFRICA
  export JEXP=/work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/EXP00/

#Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel



.. Use Dave's XIOS executable:
  Make xios or copy from James:
        ln -s ??? xios_server.exe

      Build XIOS2 @ r1080::

        cd $WDIR
        svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0
        cd $WDIR/xios-2.0
        cp ../../LBay/xios-1.0/arch/arch-XC30_ARCHER.* ./arch

      Implement make command::

        ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par

      Link xios executable to the EXP directory::

        ln -s  $WDIR/xios-2.0/bin/xios_server.exe $EXP/xios_server.exe


Build NEMO ORCHESTRA branch @ r8395::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/NERC/dev_r6998_ORCHESTRA@8395

Use Dave's XIOS file (see ``%XIOS_HOME``)::

  cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Copy Maria's MY_SRC::

  cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/CONFIG/TIDE/MY_SRC/* $CDIR/$CONFIG/MY_SRC/.

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


Copy stuff from James' simulation::

  cp $JEXP/namelist_cfg_R12 $EXP/namelist_cfg   # copy namelist_cfg
  ln -s $JEXP/../../SHARED/namelist_ref $EXP/.

Edit namelist for self determining processors assignment::

  vi namelist_cfg
  ...
  jpni        =  -20       !  jpni   number of processors following i (set automatically if < 1)
  jpnj        =  -40    !  jpnj   number of processors following j (set automatically if < 1)
  jpnij       =  -550    !  jpnij  number of local domains (set automatically if < 1)


Link other setup and forcing files::

  ln -s $JINPUTS/R12/coordinates_E-AFRICA_R12.bdy.nc $EXP/coordinates.bdy.nc
  ln -s $JINPUTS/R12/bdy_mask_E-AFRICA_R12.nc $EXP/bdy_mask.nc
  ln -s $JINPUTS/R12/domain_cfg_R12.nc $EXP/domain_cfg.nc
  ln -s $JINPUTS/R12/TIDES $EXP/TIDES


Copy in ``iodef.xml`` file and dependencies::

  rm $EXP/*xml
  #ln -s $JEXP/context_nemo.xml $EXP/.
  #ln -s $JEXP/field_def_nemo-opa.xml $EXP/.
  #ln -s $JEXP/iodef.xml $EXP/.
  #ln -s $JEXP/../../AMM12/EXP00/file_def_nemo-opa.xml $EXP/.
  #ln -s $JEXP/../../AMM12/EXP00/domain_def_nemo.xml $EXP/.

  ln -s /work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LBay/EXP00/*xml $EXP/.

Edit the runscript::

  vi runscript

  #!/bin/bash
  # ---------------------------
  #===============================================================
  # CLUSTER BITS
  #===============================================================
  #PBS -N EA_R12
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL
  #PBS -j oe
  #PBS -r n

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
  # Change to the direcotry that the job was submitted from
  cd $PBS_O_WORKDIR


  # Set the number of threads to 1
  #   This prevents any system libraries from automatically
  #   using threading.
  export OMP_NUM_THREADS=1
  # Change to the directory that the job was submitted from
  ulimit -s unlimited
  ulimit -c unlimited

  export NEMOproc=96 #550
  export XIOSproc=1

  #===============================================================
  # LAUNCH JOB
  #===============================================================
  echo `date` : Launch Job
  aprun -b -n 5 -N 5 ./xios_server.exe : -n $NEMOproc -N 24 ./opa
  exit



Fix the links with the xios (from Dave) and opa exectutables::

  ln -s /work/n01/n01/munday/XIOS/bin/xios_server.exe $EXP/.
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa

Submit::

  cd $EXP
  #qsub -q short runscript
  qsub runscript

**PENDING. DOES IT WORK?**


.. Rebuild the SSH files::

   export WDIR=/work/n01/n01/jelt/LBay/
   export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

   $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 LBay_1h_20000102_20000106_grid_T 5


  Should remove individual processor files once the build is verified::

   rm LBay_1h_20000102_20000106_grid_?_*nc

  Inspect locally e.g.::

   scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/LBay/EXP00/LBay_1h_20000102_20000106_grid_T.nc .

   ferret
   use LBay_1h_20000102_20000106_grid_T.nc
   plot /i=25/j=70 SOSSHEIG




Nasty crashing. (Not good leads to follow). Not sure this is worth pursing.





Old notes
---------

Build NEMO trunk @ r7853::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@7853
  #cp xios-2.0/arch/arch-XC30_ARCHER.* NEMOGCM/ARCH


Copy compiler keys from James::

  cd $WDIR/NEMOGCM/CONFIG


  vi ACCORD/cpp_ACCORD.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_diaharm       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero


Copy James' entire WORK directory to MY_SRC::

  cp $JEXP/../WORK/* /work/n01/n01/jelt/ACCORD/NEMOGCM/CONFIG/ACCORD/MY_SRC/.

Edit XIOS_HOME in compiler options::

  vi $WDIR/NEMOGCM/ARCH/arch-XC_ARCHER_INTEL.fcm
  ...
  %XIOS_HOME           /work/n01/n01/jelt/ACCORD/xios-2.0

On first make only choose OPA_SRC::

  ./makenemo -n ACCORD -m XC_ARCHER_INTEL -j 10 clean
  ./makenemo -n ACCORD -m XC_ARCHER_INTEL -j 10

It might break if directory structure is built from makenemo. Then remove
``key_lim2`` from cpp*fcm file and remake.


**It does compile.**





----

Look at runscript. Add module load commands::

  vi rs_12
  ...
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel
  ...
  echo `date` : Launch Job
  touch stdouterr
  rm coordinates.bdy.nc
  rm bdy_mask.nc
  rm domain_cfg.nc
  rm TIDES
  ln -s $JINPUTS/R12/coordinates_E-AFRICA_R12.bdy.nc coordinates.bdy.nc
  ln -s $JINPUTS/R12/bdy_mask_E-AFRICA_R12.nc bdy_mask.nc
  ln -s $JINPUTS/R12/domain_cfg_R12.nc domain_cfg.nc
  #ln -s $JINPUTS/R24/TIDES TIDES
  ln -s $JINPUTS/R12/TIDES TIDES
  cp namelist_cfg_R12 namelist_cfg
  aprun -b -n $NEMOproc -N 24 ./opa   >&  stdouterr_nemo : -N 1 -n $XIOSproc ./xios_server.exe >&  stdouterr_xios

---

Submit run::

  cd $EXP
  qsub rs_R12


  4819100.sdb


**PENDING: 28 Sept 2017. DOES IT WORK?**






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
