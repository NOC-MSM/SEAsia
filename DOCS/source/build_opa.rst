.. _build_nemo_label:

Build NEMO
==========

Build NEMOvp4 (ORCHESTRA branch) @trunk r8395

Physics Only:
-------------

You need to obtain a NEMO account http://forge.ipsl.jussieu.fr/nemo/register
Suggest using the same username as ARCHER account

::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

Use my XIOS file (see ``%XIOS_HOME``). Copy from *store*. Note that you are
compiling with a particular version of XIOS in mind. This must match the
xios_server.exe that is called at model run time::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure (only say YES to OPA_SRC, unless you have other plans)::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Edit the CPP flags (USE **key_diaharm_fast** instead of **key_harm_ana** for FES tides)::

  vi $CONFIG/cpp_$CONFIG.fcm

  bld::tool::fppkeys key_zdfgls        \
                   key_harm_ana       \   # key_diaharm_fast
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Add a fix to the mask variables, from the bdy mask variable::

  cp $START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dommsk.F90 $CDIR/$CONFIG/MY_SRC/.

Add a fix to permit vertical interpolation on-the-fly from initial conditions fields
onto child grid::

  cp $START_FILES/par_oce.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dtatsd.F90  $CDIR/$CONFIG/MY_SRC/.


---

Choose one or the other of the following for treatment of tides:

.. note : jelt: I think that the harmonic analysis instructions here are out of date.
 I think that Nico updated it, but I am still using his first version which I
 stored in START_FILES.

Add in POLCOMS harmonic analyisis (but old internal NEMO tides). This should speed things up...
::

  cp $START_FILES/diaharmana.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step_oce.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/bdytides.F90 $CDIR/$CONFIG/MY_SRC/.

Don't take ``sbctide.F90``, ``tide.h90``, ``tide_mod.F90``.

I editted the output to replace the ``*_x`` and ``*_y`` components::

  vi diaharmana.F90
  ...
  !      CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'x_new'//TRIM(suffix), cosamp2D(ih,:,:,jgrid) )
  !      CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'y_new'//TRIM(suffix), sinamp2D(ih,:,:,jgrid) )
        CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'x'//TRIM(suffix), cosamp2D(ih,:,:,jgrid) )
        CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'y'//TRIM(suffix), sinamp2D(ih,:,:,jgrid) )

---

Updated tides following `<FES2014_NEMO.rst>`_ ::

  cp $START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/diaharm_fast.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/sbctide.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step_oce.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tide_FES14.h90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tideini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tide_mod.F90 $CDIR/$CONFIG/MY_SRC/.


Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**Note** Make sure that any inital state files that are required by the
 configuration are copied into the MY_SRC folder before building NEMO. E.g. for
 SWPacific, inital state files are required so that a constant temp and salinity
  are set at the start of the simulation.


Physics and Biogeochemistry:
----------------------------

Get and compile NEMO along with ERSEM (BGC model) and FABM

This can be done using `make_nemo.sh <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/make_nemo.sh>`_ script::

  cd $NEMO


1. Get/download all necessary files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get NEMO::

  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

replace the TOP_SRC in NEMO::

  mv trunk_NEMOGCM_r8395/NEMO/TOP_SRC trunk_NEMOGCM_r8395/NEMO/TOP_SRC_old
  cp -r $GITCLONE/NEMO-FABM-ERSM/TOP_SRC_r8395_FABM trunk_NEMOGCM_r8395/NEMO/TOP_SRC

get ERSEM::

  cp -r $GITCLONE/NEMO-FABM-ERSM/ERSEM-master ./

get FABM::

  git clone https://github.com/fabm-model/fabm.git


2. Load modules
^^^^^^^^^^^^^^^
::

  module load cdt/15.11
  module unload PrgEnv-cray PrgEnv-gnu
  module load PrgEnv-intel/5.2.82
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.3.3.1
  module load cray-hdf5-parallel/1.8.14


3. compile FABM and copy
^^^^^^^^^^^^^^^^^^^^^^^^

Compile FABM and copy the associated libraries from your home directory to your
 fabm directory::

  module load cmake
  cd $FABM
  cmake $FABM/fabm/src -DFABM_HOST=nemo -DCMAKE_Fortran_COMPILER=ifort -DFABM_ERSEM_BASE=$WDIR/ERSEM-master -DFABM_EMBED_VERSION=ON
  make install
  cp -r /home/n01/n01/$USER/local/fabm/nemo/lib ./
  cp -r /home/n01/n01/$USER/local/fabm/nemo/include ./


4. compile NEMO and update your TOP source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get and set up your architecture files (**ATTENTION**: this file will need to be modified to have your correct paths)
::

  cp $NEMO/trunk_NEMOGCM_r8395/NEMO/TOP_SRC/arch-XC_ARCHER_INTEL_FABM.fcm $NEMO/trunk_NEMOGCM_r8395/ARCH/arch-XC_ARCHER_INTEL_FABM.fcm

make your configuration (**ATTENTION**: here we choose to include in our configuration only OPA and TOP (no ice etc.))
::

  cd $CDIR
  #make configuration first
  printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -m XC_ARCHER_INTEL_FABM -n $CONFIG -j 0

change your configuration keys to the ones we need/want and update the code source for your configuration
::

  cd $CDIR/$CONFIG
  cp $GITCLONE/NEMO-FABM-ERSM/cpp_SEAsia_FABM.fcm cpp_$CONFIG.fcm
  cp -r -f $GITCLONE/NEMO-FABM-ERSM/MY_SRC ./

You have to modify ``bldxag_cfg`` to include FABM to your compiler you can do this by simply updating this file rather than manually add them:
::

  cp $GITCLONE/NEMO-FABM-ERSM/bldxag_FABM.cfg $NEMO/trunk_NEMOGCM_r8395/TOOLS/COMPILE/bldxag.cfg

Now make the configuration with all the above updates included::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL_FABM -j 4 clean
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL_FABM -j 4


A successful compilation will generate a ``nemo.exe`` executable in
 ``$NEMO/trunk_NEMOGCM_r8395/$CONFIG/BLD/bin/``

Further information on the NEMO-ERSEM coupling system can be found here:
 `NEMO-ERSEM <https://github.com/NOC-MSM/NEMO_ERSEM>`_
