Generate a domain configuration file
====================================

.. note :

  The plan is to write this generally, with differet cases. However I am starting
  with SEAsia hybrid s-z coordiates.

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg``
directory along with all the inputs files that would have previously been needed
get v3.6 running. The reason being that all the non-time stepping stuff, like
grid generating, has been abstracted from the core OPA code and is now done as
a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory::

    ln -s $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
    ln -s $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template namelist_cfg with only the essenetial domain building stuff.
Get the size of the new domain from ``ncdump -h bathy_meter.nc``::

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
     ln_e3_dep   = .false.    ! =T : e3=dk[depth] in discret sens.
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !
     cp_cfg      =  "orca"   !  name of the configuration
     jp_cfg      =       0   !  resolution of the configuration
     jpidta      =     684   !  1st lateral dimension ( >= jpi )
     jpjdta      =     554   !  2nd    "         "    ( >= jpj )
     jpkdta      =      51   !  number of levels      ( >= jpk )
     jpiglo      =     684   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     554   !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1   !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1   !  in data domain indices
     jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_sco      = .true.    !  s-coordinate
     ln_linssh   = .false.    !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
     jphgr_msh   =       0               !  type of horizontal mesh
  ...

.. note:

  No gdept output in the offical v4 release. Though it was acheived here setting
  ln_e3_dep = F. This is needed for PyNEMO, though could be constructed from e3[tw].


---

*8 Nov 17* Revisit this DOMAINcfg generation since I am not sure that the present
settings are the most suitable.

As a starting point iterate knowledge from SWPacific dev:


.. note: At one point I tried to use s-coords but it wasn't working. Perhaps it
 would work now that the boundary problems are fixed. I did not investigate it
 though. However I leave these s-coordinate namelist options here because they
 will be useful when I need s-coordinates to work in v4.

(didn't work, though pehaps for reasons nothing to do with the coordinates...)

 (If I tried 5 levels then the function to compute the depth range breaks).
 Added in s-coordinate settings from AMM60. (There is a function that tapers dz
 near the equation. This is activated in ``domzgr.F90`` and not by a logical flag).

How about::

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
     jpidta      =     684   !  1st lateral dimension ( >= jpi )
     jpjdta      =     554   !  2nd    "         "    ( >= jpj )
     jpkdta      =      31   !  number of levels      ( >= jpk )
     jpiglo      =     684   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     554   !  2nd    -                  -    --> j  =jpjdta
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
    rn_jpk      =   51       ! Number of S levels
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

This worked (once the levels were consistent). Copy locally and inspect.

.. note :
    Alternatively try and just copy the AMM60 namelist_cfg::

      cp /work/n01/n01/jelt/NEMO/NEMOGCM/CONFIG/AMM60smago/EXP_NSea/namelist_cfg .

    Then edit for the new grid size and removal all the processor decomposition
     stuff (we are just generating the mesh and mask info)::

       diff $TDIR/DOMAINcfg/namelist_cfg diff /work/n01/n01/jelt/NEMO/NEMOGCM/CONFIG/AMM60smago/EXP_NSea/namelist_cfg

      <    jpidta      =     684   !  1st lateral dimension ( >= jpi )
      <    jpjdta      =     554   !  2nd    "         "    ( >= jpj )
      ---
      >    jpidta      =     1120               !  1st lateral dimension ( >= jpi )
      >    jpjdta      =     1440               !  2nd    "         "    ( >= jpj )
      33,34c33,34
      <    jpiglo      =     684               !  1st dimension of global domain --> i =jpidta
      <    jpjglo      =     554              !  2nd    -                  -    --> j  =jpjdta
      ---
      >    jpiglo      =     1120               !  1st dimension of global domain --> i =jpidta
      >    jpjglo      =     1440              !  2nd    -                  -    --> j  =jpjdta
      499a500,506
      >    cn_mpi_send =  'I'      !  mpi send/recieve type   ='S', 'B', or 'I' for standard send,
      >                            !  buffer blocking send or immediate non-blocking sends, resp.
      >    nn_buffer   =   0       !  size in bytes of exported buffer ('B' case), 0 no exportation
      >    ln_nnogather=  .false.  !  activate code to avoid mpi_allgather use at the northfold
      >    jpni        =   40       !  jpni   number of processors following i (set automatically if < 1)
      >    jpnj        =   50     !  jpnj   number of processors following j (set automatically if < 1)
      >    jpnij       =   2000     !  jpnij  number of local domains (set automatically if < 1)

---

Build a script to run the executable::

  vi $TDIR/DOMAINcdf/rs

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

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

  cp $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  cp $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.

.. mote :  should check the difference between the homemade sco version the AMM60
  verison did:      ``diff namelist_cfg_sco_WIP namelist_cfg_AMM60``

.. note : alternativly should check the difference between the AMM60 and local
  output.namelist.dyn: ``diff output.namelist.dyn /work/n01/n01/jelt/NEMO/NEMOGCM/CONFIG/AMM60smago/EXP_NSea/output.namelist.dyn``
  I notice that rmax is different.

.. note : There are a lot of unknown parameters in these settings. And I don't
  want to find i made some naive error in six months. Looking at the domain there
  are some serious trenches near land. S-coords will not work well there. Conversely,
  ODA is all abount the near-coastal environment. There is a strong case for using
  hybrid s-z coordinates a la NNA...

  James noted that high resolution neat the bed caused significant difficulty in
  deep water stability. Whereas you want it on the shelf. Hence regular stretched
  s-coords wont really work.

  PyNEMO outputs boundary conditions on the parent z-grid. This can be interpolated
  at run-time to the child grid.
