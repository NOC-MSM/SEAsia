Changes to the namelist_cfg when migrating from NEMO v3.6 to v4
===============================================================

In order to make a simulation with the new v4 namelist_cfg,  a number of
changes appear with its format. Some of these are noted below. But since they are
a cut down from a debug log they may not be exhaustive.

---

I might have added some tidal harmonic constituent names to the template::

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters (#ifdef key_tide)
  !-----------------------------------------------------------------------
  clname(1)    = 'Q1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(2)    = 'O1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(3)    = 'P1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(4)    = 'K1'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(5)    = 'N2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(6)   =  'M2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(7)   = 'S2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(8)   = 'K2'   !  name of constituent - all tidal components must be set in namelist_cfg
  clname(9)   = 'M4'   !  name of constituent - all tidal components must be set in namelist_cfg


Also noted that the sbc namelist variables have changed. Now use ``ln_blk`` and
 ``ln_COARE_3p5= .true.`` instead of ``ln_blk_core``. *(Except that I use*
 *ln_blk=F because I don't have a SLP variable for the COARE algorithm to work...)*

There was an extra variable in namdom. Comment out ldbletanh (which was essential in DOMAINcfg)::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
  ...
  !    ldbletanh   =    .false.             !  Use/do not use double tanf function for vertical coordinates



The configuration namelist variable is simplified with the new domain_cfg.nc file::

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.   !  (=T) read the domain configuration file
        !                   !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
        cn_domcfg = "domain_cfg"         ! domain configuration filename
        !
     ln_write_cfg= .false.   !  (=T) create the domain configuration file
        cn_domcfg_out = "domain_cfg_out" ! newly created domain configuration filename
        !
     ln_use_jattr = .false.  !  use (T) the file attribute: open_ocean_jstart, if present
     !                       !  in netcdf input files, as the start j-row for reading
  /


Equation of state needed fixing::

  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
     ln_teos10   = .false.         !  = Use TEOS-10 equation of state
     ln_eos80    = .true.         !  = Use EOS80 equation of state
     ln_seos     = .false.         !  = Use simplified equation of state (S-EOS)


``&namdom`` has a new format (I think)::

  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
   ln_linssh   = .false.   !  =T  linear free surface  ==>>  model level are fixed in time
   nn_closea   =    0      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
   !
   nn_msh      =    0      !  create (>0) a mesh file or not (=0)
   rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice
   !
   rn_rdt      =  60.     !  time step for the dynamics (and tracer if nn_acc=0)
   rn_atfp     =    0.1    !  asselin time filter parameter
   !
   ln_crs      = .false.   !  Logical switch for coarsening module


Note ``ln_linssh`` seems to have moved to ``&namdom``. (I think it was in ``&namzgr_sco``
 before)

I added in s-coordinate option and which seemed to have vanished from the default
namelist_cfg template::

  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate                                  (default: NO selection)
  !-----------------------------------------------------------------------
     ln_zco      = .false.   !  z-coordinate - full    steps
     ln_zps      = .false.   !  z-coordinate - partial steps
     ln_sco      = .true.   !  s- or hybrid z-s-coordinate
     ln_isfcav   = .false.   !  ice shelf cavity
  /
  !-----------------------------------------------------------------------
  &namzgr_sco
  !-----------------------------------------------------------------------
     ln_s_sf12=.true.,
     ln_s_sh94=.false.,
     ln_sigcrit=.true.,
     rn_alpha=4.4,
     rn_bb=0.8,
     rn_efold=0.0,
     rn_hc=10000.0,
     rn_rmax=1.0,
     rn_sbot_max=10000.0,
     rn_sbot_min=1.0e-5,
     rn_theta=6.0,
     rn_thetb=1.0,
     rn_zb_a=0.024,
     rn_zb_b=-0.2,
     rn_zs=1.0,
  /


Though the fields are set up (except for SLP field), no met forcing is used::

  !-----------------------------------------------------------------------
  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
     ln_usr      = .true.    !  user defined formulation                  (T => check usrdef_sbc)
     ln_flx      = .false.   !  flux formulation                          (T => fill namsbc_flx )
     ln_blk      = .false.    !  Bulk formulation                          (T => fill namsbc_blk )

Though the fields are ready to go I have not used inital T,S conditions::

  !-----------------------------------------------------------------------
  &namtsd        !   data : Temperature  & Salinity
  !-----------------------------------------------------------------------
  !              !  file name                 ! frequency (hours) ! variable ! time interp.!  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                            !  (if <0  months)  !   name   !  (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_tem  = 'initcd_votemper',         -1        ,'votemper' ,    .false.    , .true. , 'yearly'   , ''       ,   ''    ,    ''
     sn_sal  = 'initcd_vosaline',         -1        ,'vosaline' ,    .false.    , .true. , 'yearly'   , ''       ,   ''    ,    ''
     !
     cn_dir        = '../../../../INPUTS/'     !  root directory for the location of the runoff files
     ln_tsd_init   = .false.   !  Initialisation of ocean T & S with T &S input data (T) or not (F)
     ln_tsd_tradmp = .false.   !  damping of ocean T & S toward T &S input data (T) or not (F)

No lateral diffusion in momentum::

  ln_dynldf_lap =  .false.    !    laplacian operator
  ln_dynldf_blp =  .false.    !  bilaplacian operator

Don't need a boundaries mask file. But the rimwidth must match that used in the
 PyNEMO process::

  !-----------------------------------------------------------------------
  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
    ...
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = 'bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
    ...
    nn_rimwidth   = 9                    !  width of the relaxation zone

Quadratic bottom friction - This means changing to nonlinear friction, with a coeff 2.5E-3, log layer and
bottom roughness as follows
::
  !-----------------------------------------------------------------------
  &nambfr        !   bottom friction
  !-----------------------------------------------------------------------
     nn_bfr      =    2      !  type of bottom friction :   = 0 : free slip,  = 1 : linear friction
                             !                              = 2 : nonlinear friction
     rn_bfri2    =    2.5e-3 !  bottom drag coefficient (non linear case)
     rn_bfeb2    =    0.0e0  !  bottom turbulent kinetic energy background  (m2/s2)
     ln_loglayer =    .true. !  loglayer bottom friction (only effect when nn_bfr = 2)
     rn_bfrz0    =    0.003  !  bottom roughness (only effect when ln_loglayer = .true.)
  /

I also changed the rn_bfeb2 = 2.5e-3 to zero. This would depend on the tke scheme.
I am not sure what it should be...
