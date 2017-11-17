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
For example:

   * For hybrid z-s coordinates: `SEAsia_DOMAINcfg_namelist_cfg`_
   * For z-partial step coordinates: `SWPacific_DOMAINcfg_namelist_cfg`_
   * For stretched sigma coordinates: `GThai_DOMAINcfg_namelist_cfg`_ (?Jason to edit coords etc)


For any vertical coordinate option, edit the size of the new domain. Get the
 size from ``ncdump -h bathy_meter.nc`` (E.g. for SE Asia 1/12 degree, with 75 levels)::

  cd $TDIR/DOMAINcfg
  vi namelist_cfg

  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ...
     jpidta      =     684   !  1st lateral dimension ( >= jpi )
     jpjdta      =     554   !  2nd    "         "    ( >= jpj )
     jpkdta      =      75   !  number of levels      ( >= jpk )
     jpiglo      =     684   !  1st dimension of global domain --> i =jpidta
     jpjglo      =     554   !  2nd    -                  -    --> j  =jpjdta
     ...


 .. note:

   No gdept output in the offical v4 release. Though it was acheived here setting
   ln_e3_dep = F. This is needed for PyNEMO, though could be constructed from e3[tw].

   In hybrid s-z example not outputting gdep[wt]. Will have to fix PyNEMO. Others
   to be updated in due coarse.

In the following special choices are made according to the choice of vertical coordinates

Hybrid z-s coordinates
======================

Vertical coordinates options (e.g. `SEAsia_DOMAINcfg_namelist_cfg`_)::

  vi namelist_cfg
  ...
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)      ("key_zco" may also be defined)
     ln_zps      = .false.    !  z-coordinate - partial steps
     ln_sco      = .true.   !  s- or hybrid z-s-coordinate    (T/F)
     ln_isfcav   = .false.   !  ice shelf cavity
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate                (default F)
  !-----------------------------------------------------------------------
     ln_s_sh94   = .true.    !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
     ln_s_sf12   = .false.   !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
     rn_sbot_min =   10.0    !  minimum depth of s-bottom surface (>0) (m)
     rn_sbot_max = 6000.0    !  maximum depth of s-bottom surface (= ocean depth) (>0) (m)
     rn_hc       =  39.0     !  critical depth for transition to stretched coordinates
                             !!!!!!!  Envelop bathymetry
     rn_rmax     =    0.05   !  maximum cut-off r-value allowed (0<r_max<1)
                             !!!!!!!  hybrid z-s parameters
     nn_sig_lev = 39
     ln_s_melange = .true.
   /
   !-----------------------------------------------------------------------
   &namdom        !   space and time domain (bathymetry, mesh, timestep)
   !-----------------------------------------------------------------------
      nn_closea   =    1      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
      nn_msh      =    0      !  create (=1) a mesh file or not (=0)
      rn_hmin     =   -10.    !  min depth of the ocean (>0) or min number of ocean level (<0)
      rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice
      rn_e3zps_min=   25.     !  partial step thickness is set larger than the minimum of
      rn_e3zps_rat=    0.2    !  rn_e3zps_min and rn_e3zps_rat*e3t, with 0<rn_e3zps_rat<1
      rn_rdt      =   300.    !  time step for the dynamics (and tracer if nacc=0)   ==> 5760
      jphgr_msh   =       0               !  type of horizontal mesh
                                          !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
                                          !  = 1 geographical mesh on the sphere with regular grid-spacing
                                          !  = 2 f-plane with regular grid-spacing
                                          !  = 3 beta-plane with regular grid-spacing
                                          !  = 4 Mercator grid with T/U point at the equator
      ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
      ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
      ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
      ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
      ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
      ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
      ppsur       =    -3958.951371276829 !  ORCA r4, r2 and r05 coefficients
      ppa0        =     103.9530096000000 ! (default coefficients)
      ppa1        =     2.415951269000000 !
      ppkth       =     15.35101370000000 !
      ppacr       =        7.0            !
      ppdzmin     =  999999.              !  Minimum vertical spacing
      pphmax      =  999999.              !  Maximum depth
      ldbletanh   =    .TRUE.             !  Use/do not use double tanf function for vertical coordinates
      ppa2        =      100.760928500000 !  Double tanh function parameters
      ppkth2      =       48.029893720000 !
      ppacr2      =       13.000000000000 !
      /

.. note: 17 Nov 17. Commented out the ztaper action near the equator in domzgr.f90

The new parameter sig_lev = 39 levels sets the merge level (counting from the
top down) at which point sigma transitions to z-partial step.
This should be set equal to rn_hc, which is the depth (m) for transition
 to from stretch (deeper) to pure (shallower) sigma coordinates.

Check number of levels ``jpkdta`` exceeds ``nn_sig_lev``.

The transition between zps and sco is here at 400m.

Therefore for bathymetry >= 39m the top layer with be 1m thick. (For water
 shallower than 39m it will be proportionately less.)

.. note: If the zps coords have a 0.5 metre top grid box then rn_hc has to come down to 19.5m.
 This makes the set up slightly inflexible, although it was never written to be a generic solution.
 *There is no reason to take a different underlying function that defines the vertical levels. At the moment it is the double tanh function from the ZPS code - it could quite easily be the sh94.*


The ``namdom`` are specified here to ensure the match the ORCA grid. Though are
probably not strickly necessary.

To get the new option to work I have copied a file into ``src``. This will
eventually be in the trunk but for now::

  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/TOOLS/DOMAINcfg/src/domzgr.f90.jelt $TDIR/DOMAIN_cfg/src/domzgr.f90

Recompile the tool e.g.::

  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg clean
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg



z-partial step coordinates
==========================

Vertical coordinates options (e.g. `SWPacific_DOMAINcfg_namelist_cfg`_)::

  vi namelist_cfg
  ...






Stetched sigma coordinates
==========================

.. note : At present these notes actually correspond to an early version of the
    SEAsia domain, though could be updated to the Gulf of Thailand.

Vertical coordinates options (e.g. `GThai_DOMAINcfg_namelist_cfg`_)::


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


*(8 Nov 2017)*
This worked (once the levels between EXP nemo and domain_cfg were consistent).

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


Build a script to run the executable
====================================

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
  #PBS -M xxx@noc.ac.uk
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

Change the notification email to your own address::

  sed -i "s/xxx@noc/$USER@noc/g" rs

Try running it::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Copy domain_cfg.nc to the EXP directory (also copy it to the INPUTS directory, which stores
 the bits and bobs for a rebuild)::

  rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  rsync -utv $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.

.. note :  James ran into memory issues with a big domain and needed to the run
 script to work in parallel.
 See /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/TOOLS/DOMAINcfg/rs
 and /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/TOOLS/DOMAINcfg/rs_rebuild


I checked the ``domain_cfg.nc`` output with a python script `inspect_domain_cfg.py`_::

  module load anaconda
  python inspect_domain_cfg.py

which produces a transect of the grid across a section.
