Build NEMO (ORCHESTRA) trunk @ r8395
++++++++++++++++++++++++++++++++++++
::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

Use my XIOS file (see ``%XIOS_HOME``). Copy from *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Edit the CPP flags::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_diaharm       \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10
