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
