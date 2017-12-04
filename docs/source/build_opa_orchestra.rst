Build NEMO (ORCHESTRA) trunk @ r8395
++++++++++++++++++++++++++++++++++++

You need to obtain a nemo account http://forge.ipsl.jussieu.fr/nemo/register
Suggest using the same unsernam as ARCHER account

::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

Use my XIOS file (see ``%XIOS_HOME``). Copy from *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure (only say YES to OPA_SRC, unless you have other plans)::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Edit the CPP flags::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_diaharm       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Add a fix to the mask variables, from the bdy mask variable::

  cp $START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dommsk.F90 $CDIR/$CONFIG/MY_SRC/.

Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**Note** Make sure that any inital state files that are required by the configuration are copied into the MY_SRC folder before building NEMO. E.g. for SWPacific, inital state files are required so that a constant temp and salinity are set at the start of the simulation.
