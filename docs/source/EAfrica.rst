===========================================
Setting up a East Africa NEMO configuration
===========================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/EAfrica.html


Summary / Plan
==============

* Check James' /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/
configuration for suggested namelist options for a lighter tide only run.

* Try and replicate run either with his code or ORCHESTRA code. The idea being to
copy a config that is tide only.

**FAILED to replicate**. Diminishing returns. Abort thread.


Recipe Notes
============

Define working directory, and other useful shortcuts::

  export CONFIG=ACCORD
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_EAFRICA

  export JINPUTS=/work/n01/n01/jdha/2017/INPUTS/ODA/E-AFRICA
  export JEXP=/work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ODA_E-AFRICA/EXP00/

#Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

---

Build XIOS2 @ r1080::

  cd /work/n01/n01/$USER
  svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1080
  cd xios-2.0_r1080
  #cp ../../LBay/xios-1.0/arch/arch-XC30_ARCHER.* ./arch
  cp ../LBay/xios-2.0/arch/arch-XC30_ARCHER* arch/.

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par

Link the xios-2.0_r1080 to a generic XIOS directory name::

  ln -s  /work/n01/n01/$USER/xios-2.0_r108  /work/n01/n01/$USER/XIOS

Link xios executable to the EXP directory::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

---

Build NEMO trunk @ ORCHESTRA r8395::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

Use my XIOS file (see ``%XIOS_HOME``)::

  cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.
  vi $CDIR/../ARCH/arch-XC_ARCHER_INTEL.fcm
  ..
  %XIOS_HOME           /work/n01/n01/jelt/XIOS


Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Edit the CPP flags::

  vi ACCORD/cpp_ACCORD.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_diaharm       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**It compiles!!**


Copy INPUT stuff from James' simulation::

  mkdir $EXP

  cp $JEXP/namelist_cfg_R12 $EXP/namelist_cfg   # copy namelist_cfg
  ln -s $JEXP/../../SHARED/namelist_ref $EXP/.

Edit namelist for self determining processors assignment::

  vi $EXP/namelist_cfg
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
  ln -s $JEXP/context_nemo.xml $EXP/.
  ln -s $JEXP/field_def_nemo-opa.xml $EXP/.
  ln -s $JEXP/iodef.xml $EXP/.
  ln -s $JEXP/../../AMM12/EXP00/file_def_nemo-opa.xml $EXP/.
  ln -s $JEXP/../../AMM12/EXP00/domain_def_nemo.xml $EXP/.

  #ln -s /work/n01/n01/jelt/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LBay/EXP00/*xml $EXP/.

Edit/create the runscript::

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



Fix the links with the xios (if not already done) and opa exectutables::

  ln -s /work/n01/n01/jelt/XIOS/bin/xios_server.exe $EXP/.
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa

Submit::

  cd $EXP
  qsub -q short runscript
  #qsub runscript

*(4 Oct 2017)*
**It runs and outputs **

---

Rebuild the files and inspect locally
=====================================

Rebuild the SSH files (use an already compiled TOOLS)::

  export WDIR=/work/n01/n01/jelt/LBay/
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 EA_v3_1d_20010101_20010112_grid_T 5

Should remove individual processor files once the build is verified::

  rm EA_v3_1d_20010101_20010112_grid_T_*.nc

Inspect locally e.g.::

  scp jelt@login.archer.ac.uk:$EXP/EA_v3_1d_20010101_20010112_grid_T.nc .
  #scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/trunk_NEMOGCM_r8395/CONFIG/ACCORD/EXP_EAFRICA/EA_v3_1d_20010101_20010112_grid_T.nc .

  ferret
  use EA_v3_1d_20010101_20010112_grid_T.nc
  plot /i=25/j=70 SOSSHEIG
