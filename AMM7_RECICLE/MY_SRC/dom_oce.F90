MODULE dom_oce
   !!======================================================================
   !!                       ***  MODULE dom_oce  ***
   !!       
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!======================================================================
   !! History :  1.0  ! 2005-10  (A. Beckmann, G. Madec)  reactivate s-coordinate 
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            4.0  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012     (S. Mocavero, I. Epicoco) Add arrays associated
   !!                             to the optimization of BDY communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Root    : dummy function used when lk_agrif=F
   !!   Agrif_CFixed  : dummy function used when lk_agrif=F
   !!   dom_oce_alloc : dynamical allocation of dom_oce arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PUBLIC             ! allows the acces to par_oce when dom_oce is used
   !                  ! exception to coding rules... to be suppressed ???

   PUBLIC dom_oce_alloc  ! Called from nemogcm.F90

   !!----------------------------------------------------------------------
   !! time & space domain namelist
   !! ----------------------------
   !                                    !!* Namelist namdom : time & space domain *
   INTEGER , PUBLIC ::   nn_bathy        !: = 0/1 ,compute/read the bathymetry file
   REAL(wp), PUBLIC ::   rn_bathy        !: depth of flat bottom (active if nn_bathy=0; if =0 depth=jpkm1)
   REAL(wp), PUBLIC ::   rn_hmin         !: minimum ocean depth (>0) or minimum number of ocean levels (<0)
   REAL(wp), PUBLIC ::   rn_e3zps_min    !: miminum thickness for partial steps (meters)
   REAL(wp), PUBLIC ::   rn_e3zps_rat    !: minimum thickness ration for partial steps
   INTEGER , PUBLIC ::   nn_msh          !: = 1 create a mesh-mask file
   INTEGER , PUBLIC ::   nn_acc          !: = 0/1 use of the acceleration of convergence technique
   REAL(wp), PUBLIC ::   rn_atfp         !: asselin time filter parameter
   REAL(wp), PUBLIC ::   rn_rdt          !: time step for the dynamics (and tracer if nacc=0)
   REAL(wp), PUBLIC ::   rn_rdtmin       !: minimum time step on tracers
   REAL(wp), PUBLIC ::   rn_rdtmax       !: maximum time step on tracers
   REAL(wp), PUBLIC ::   rn_rdth         !: depth variation of tracer step
   INTEGER , PUBLIC ::   nn_closea       !: =0 suppress closed sea/lake from the ORCA domain or not (=1)
   INTEGER , PUBLIC ::   nn_euler        !: =0 start with forward time step or not (=1)
   LOGICAL , PUBLIC ::   ln_crs          !: Apply grid coarsening to dynamical model output or online passive tracers

   !! Time splitting parameters
   !! =========================
   LOGICAL,  PUBLIC :: ln_bt_fw          !: Forward integration of barotropic sub-stepping
   LOGICAL,  PUBLIC :: ln_bt_av          !: Time averaging of barotropic variables
   LOGICAL,  PUBLIC :: ln_bt_nn_auto     !: Set number of barotropic iterations automatically
   INTEGER,  PUBLIC :: nn_bt_flt         !: Filter choice
   INTEGER,  PUBLIC :: nn_baro           !: Number of barotropic iterations during one baroclinic step (rdt)
   REAL(wp), PUBLIC :: rn_bt_cmax        !: Maximum allowed courant number (used if ln_bt_nn_auto=T)

   !! Horizontal grid parameters for domhgr
   !! =====================================
   INTEGER       ::   jphgr_msh        !: type of horizontal mesh
   !                                       !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
   !                                       !  = 1 geographical mesh on the sphere with regular grid-spacing
   !                                       !  = 2 f-plane with regular grid-spacing
   !                                       !  = 3 beta-plane with regular grid-spacing
   !                                       !  = 4 Mercator grid with T/U point at the equator

   REAL(wp)      ::   ppglam0              !: longitude of first raw and column T-point (jphgr_msh = 1)
   REAL(wp)      ::   ppgphi0              !: latitude  of first raw and column T-point (jphgr_msh = 1)
   !                                                        !  used for Coriolis & Beta parameters (jphgr_msh = 2 or 3)
   REAL(wp)      ::   ppe1_deg             !: zonal      grid-spacing (degrees)
   REAL(wp)      ::   ppe2_deg             !: meridional grid-spacing (degrees)
   REAL(wp)      ::   ppe1_m               !: zonal      grid-spacing (degrees)
   REAL(wp)      ::   ppe2_m               !: meridional grid-spacing (degrees)

   !! Vertical grid parameter for domzgr
   !! ==================================
   REAL(wp)      ::   ppsur                !: ORCA r4, r2 and r05 coefficients
   REAL(wp)      ::   ppa0                 !: (default coefficients)
   REAL(wp)      ::   ppa1                 !:
   REAL(wp)      ::   ppkth                !:
   REAL(wp)      ::   ppacr                !:
   !
   !  If both ppa0 ppa1 and ppsur are specified to 0, then
   !  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
   REAL(wp)      ::   ppdzmin              !: Minimum vertical spacing
   REAL(wp)      ::   pphmax               !: Maximum depth
   !
   LOGICAL       ::   ldbletanh            !: Use/do not use double tanf function for vertical coordinates
   REAL(wp)      ::   ppa2                 !: Double tanh function parameters
   REAL(wp)      ::   ppkth2               !:
   REAL(wp)      ::   ppacr2               !:

   !                                    !! old non-DOCTOR names still used in the model
   INTEGER , PUBLIC ::   ntopo           !: = 0/1 ,compute/read the bathymetry file
   REAL(wp), PUBLIC ::   e3zps_min       !: miminum thickness for partial steps (meters)
   REAL(wp), PUBLIC ::   e3zps_rat       !: minimum thickness ration for partial steps
   INTEGER , PUBLIC ::   nmsh            !: = 1 create a mesh-mask file
   INTEGER , PUBLIC ::   nacc            !: = 0/1 use of the acceleration of convergence technique
   REAL(wp), PUBLIC ::   atfp            !: asselin time filter parameter
   REAL(wp), PUBLIC ::   rdt             !: time step for the dynamics (and tracer if nacc=0)
   REAL(wp), PUBLIC ::   rdtmin          !: minimum time step on tracers
   REAL(wp), PUBLIC ::   rdtmax          !: maximum time step on tracers
   REAL(wp), PUBLIC ::   rdth            !: depth variation of tracer step

   !                                                  !!! associated variables
   INTEGER , PUBLIC                 ::   neuler        !: restart euler forward option (0=Euler)
   REAL(wp), PUBLIC                 ::   atfp1         !: asselin time filter coeff. (atfp1= 1-2*atfp)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   rdttra  !: vertical profile of tracer time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   r2dtra  !: = 2*rdttra except at nit000 (=rdttra) if neuler=0

   !                                         !!* Namelist namcla : cross land advection
   INTEGER, PUBLIC ::   nn_cla               !: =1 cross land advection for exchanges through some straits (ORCA2)

   !!----------------------------------------------------------------------
   !! space domain parameters
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   lzoom      =  .FALSE.   !: zoom flag
   LOGICAL, PUBLIC ::   lzoom_e    =  .FALSE.   !: East  zoom type flag
   LOGICAL, PUBLIC ::   lzoom_w    =  .FALSE.   !: West  zoom type flag
   LOGICAL, PUBLIC ::   lzoom_s    =  .FALSE.   !: South zoom type flag
   LOGICAL, PUBLIC ::   lzoom_n    =  .FALSE.   !: North zoom type flag

   !                                     !!! domain parameters linked to mpp
   INTEGER, PUBLIC ::   nperio            !: type of lateral boundary condition
   INTEGER, PUBLIC ::   nimpp, njmpp      !: i- & j-indexes for mpp-subdomain left bottom
   INTEGER, PUBLIC ::   nreci, nrecj      !: overlap region in i and j
   INTEGER, PUBLIC ::   nproc             !: number for local processor
   INTEGER, PUBLIC ::   narea             !: number for local area
   INTEGER, PUBLIC ::   nbondi, nbondj    !: mark of i- and j-direction local boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy(:)    !: mark i-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy(:)    !: mark j-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy_b(:)  !: mark i-direction of neighbours local boundaries for BDY open boundaries  
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy_b(:)  !: mark j-direction of neighbours local boundaries for BDY open boundaries  

   INTEGER, PUBLIC ::   npolj             !: north fold mark (0, 3 or 4)
   INTEGER, PUBLIC ::   nlci, nldi, nlei  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   nlcj, nldj, nlej  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   noea, nowe        !: index of the local neighboring processors in
   INTEGER, PUBLIC ::   noso, nono        !: east, west, south and north directions
   INTEGER, PUBLIC ::   npne, npnw        !: index of north east and north west processor
   INTEGER, PUBLIC ::   npse, npsw        !: index of south east and south west processor
   INTEGER, PUBLIC ::   nbne, nbnw        !: logical of north east & north west processor
   INTEGER, PUBLIC ::   nbse, nbsw        !: logical of south east & south west processor
   INTEGER, PUBLIC ::   nidom             !: ???

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig        !: local  ==> global domain i-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg        !: local  ==> global domain j-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mi0, mi1   !: global ==> local  domain i-index    !!bug ==> other solution?
   !                                                  ! (mi0=1 and mi1=0 if the global index is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mj0, mj1   !: global ==> local  domain j-index     !!bug ==> other solution?
   !                                                  ! (mi0=1 and mi1=0 if the global index is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nimppt, njmppt   !: i-, j-indexes for each processor
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ibonit, ibonjt   !: i-, j- processor neighbour existence
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nlcit , nlcjt    !: dimensions of every subdomain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nldit , nldjt    !: first, last indoor index for each i-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nleit , nlejt    !: first, last indoor index for each j-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nfiimpp, nfipproc, nfilcit

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  glamt, glamu   !: longitude of t-, u-, v- and f-points (degre)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  glamv, glamf   !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  gphit, gphiu   !: latitude  of t-, u-, v- and f-points (degre)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  gphiv, gphif   !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::  e1t, e2t, r1_e1t, r1_e2t   !: horizontal scale factors and inverse at t-point (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::  e1u, e2u, r1_e1u, r1_e2u   !: horizontal scale factors and inverse at u-point (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::  e1v, e2v, r1_e1v, r1_e2v   !: horizontal scale factors and inverse at v-point (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::  e1f, e2f, r1_e1f, r1_e2f   !: horizontal scale factors and inverse at f-point (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  e1e2t          !: surface at t-point (m2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ff             !: coriolis factor (2.*omega*sin(yphi) ) (s-1)

   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! ---------------------------------------------------------------------
   !                                 !!* Namelist namzgr : vertical coordinate *
   LOGICAL, PUBLIC ::   ln_zco        !: z-coordinate - full step
   LOGICAL, PUBLIC ::   ln_zps        !: z-coordinate - partial step
   LOGICAL, PUBLIC ::   ln_sco        !: s-coordinate or hybrid z-s coordinate
   LOGICAL, PUBLIC ::   ln_sco_interp        !: s-coordinate or hybrid z-s coordinate
   LOGICAL, PUBLIC ::   ln_isfcav     !: presence of ISF 

   !! All coordinates
   !! ---------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdep3w_0           !: depth of t-points (sum of e3w) (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdept_0, gdepw_0   !: analytical (time invariant) depth at t-w  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3v_0  , e3f_0     !: analytical (time invariant) vertical scale factors at  v-f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3t_0  , e3u_0     !:                                      t-u  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3vw_0             !: analytical (time invariant) vertical scale factors at  vw
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3w_0  , e3uw_0    !:                                      w-uw points (m)
#if defined key_vvl
   LOGICAL, PUBLIC, PARAMETER ::   lk_vvl = .TRUE.    !: variable grid flag

   !! All coordinates
   !! ---------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdep3w_n           !: now depth of T-points (sum of e3w) (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdept_n, gdepw_n   !: now depth at T-W  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdept_b, gdepw_b   !: before depth at T-W  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3t_n              !: now    vertical scale factors at  t       point  (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3u_n  , e3v_n     !:            -      -      -    -   u --v   points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3w_n  , e3f_n     !:            -      -      -    -   w --f   points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3uw_n , e3vw_n    !:            -      -      -    -   uw--vw  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3t_b              !: before     -      -      -    -   t       points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3w_b              !: before     -      -      -    -   t       points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3u_b  , e3v_b     !:   -        -      -      -    -   u --v   points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3uw_b , e3vw_b    !:   -        -      -      -    -   uw--vw  points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3t_a              !: after      -      -      -    -   t       point  (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   e3u_a  , e3v_a     !:   -        -      -      -    -   u --v   points (m)
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vvl = .FALSE.   !: fixed grid flag
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hur  , hvr     !: Now    inverse of u and v-points ocean depth (1/m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hu   , hv      !:        depth at u- and v-points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ht             !:        depth at t-points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ehur_a, ehvr_a !: After  inverse of u and v-points ocean depth (1/m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ehu_a , ehv_a  !:        depth at u- and v-points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ehur_b, ehvr_b !: Before inverse of u and v-points ocean depth (1/m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ehu_b , ehv_b  !:        depth at u- and v-points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ht_0           !: reference depth at t-       points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hu_0 , hv_0    !: reference depth at u- and v-points (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   re2u_e1u       !: scale factor coeffs at u points (e2u/e1u)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   re1v_e2v       !: scale factor coeffs at v points (e1v/e2v)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e12t , r1_e12t !: horizontal cell surface and inverse at t points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e12u , r1_e12u !: horizontal cell surface and inverse at u points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e12v , r1_e12v !: horizontal cell surface and inverse at v points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e12f , r1_e12f !: horizontal cell surface and inverse at f points

   INTEGER, PUBLIC ::   nla10              !: deepest    W level Above  ~10m (nlb10 - 1)
   INTEGER, PUBLIC ::   nlb10              !: shallowest W level Bellow ~10m (nla10 + 1) 

   !! z-coordinate with full steps (also used in the other cases as reference z-coordinate)
   !! =-----------------====------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gdept_1d, gdepw_1d !: reference depth of t- and w-points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_1d  , e3w_1d   !: reference vertical scale factors at T- and W-pts (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e3tp    , e3wp     !: ocean bottom level thickness at T and W points

   !! s-coordinate and hybrid z-s-coordinate
   !! =----------------======---------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gsigt, gsigw       !: model level depth coefficient at t-, w-levels (analytic)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gsi3w              !: model level depth coefficient at w-level (sum of gsigw)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   esigt, esigw       !: vertical scale factor coef. at t-, w-levels

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hbatv , hbatf      !: ocean depth at the vertical of  v--f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hbatt , hbatu      !:                                 t--u points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   scosrf, scobot     !: ocean surface and bottom topographies 
   !                                                                           !  (if deviating from coordinate surfaces in HYBRID)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hifv  , hiff       !: interface depth between stretching at v--f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hift  , hifu       !: and quasi-uniform spacing             t--u points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rx1                !: Maximum grid stiffness ratio

   !!----------------------------------------------------------------------
   !! masks, bathymetry
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbathy             !: number of ocean level (=0, 1, ... , jpk-1)
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbkt               !: vertical index of the bottom last T- ocean level
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbku, mbkv         !: vertical index of the bottom last U- and W- ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bathy                              !: ocean depth (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmask_i, umask_i, vmask_i, fmask_i !: interior domain T-point mask
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bmask                              !: land/ocean mask of barotropic stream function

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   misfdep                 !: top first ocean level                (ISF)
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mikt, miku, mikv, mikf  !: first wet T-, U-, V-, F- ocean level (ISF)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   risfdep                 !: Iceshelf draft                       (ISF)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssmask                   !: surface domain T-point mask 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: tmask, umask, vmask, fmask   !: land/ocean mask at T-, U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: wmask, wumask, wvmask        !: land/ocean mask at WT-, WU- and WV-pts

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   tpol, fpol          !: north fold mask (jperio= 3 or 4)

#if defined key_noslip_accurate
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  )  :: npcoa              !: ???
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nicoa, njcoa       !: ???
#endif

   !!----------------------------------------------------------------------
   !! calendar variables
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC ::   nyear         !: current year
   INTEGER , PUBLIC ::   nmonth        !: current month
   INTEGER , PUBLIC ::   nday          !: current day of the month
   INTEGER , PUBLIC ::   ndastp        !: time step date in yyyymmdd format
   INTEGER , PUBLIC ::   nday_year     !: current day counted from jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_year     !: current time step counted in second since 00h jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_month    !: current time step counted in second since 00h 1st day of the current month
   INTEGER , PUBLIC ::   nsec_week     !: current time step counted in second since 00h of last monday
   INTEGER , PUBLIC ::   nsec_day      !: current time step counted in second since 00h of the current day
   REAL(wp), PUBLIC ::   fjulday       !: current julian day 
   REAL(wp), PUBLIC ::   fjulstartyear !: first day of the current year in julian days
   REAL(wp), PUBLIC ::   adatrj        !: number of elapsed days since the begining of the whole simulation
   !                                   !: (cumulative duration of previous runs that may have used different time-step size)
   INTEGER , PUBLIC, DIMENSION(0: 2) ::   nyear_len     !: length in days of the previous/current/next year
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_len    !: length in days of the months of the current year
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_half   !: second since Jan 1st 0h of the current year and the half of the months
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_end    !: second since Jan 1st 0h of the current year and the end of the months
   INTEGER , PUBLIC                  ::   nsec1jan000   !: second since Jan 1st 0h of nit000 year and Jan 1st 0h the current year

   !!----------------------------------------------------------------------
   !! mpp reproducibility
   !!----------------------------------------------------------------------
#if defined key_mpp_rep
   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp_rep = .TRUE.    !: agrif flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp_rep = .FALSE.   !: agrif flag
#endif

   !!----------------------------------------------------------------------
   !! agrif domain
   !!----------------------------------------------------------------------
#if defined key_agrif
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .TRUE.    !: agrif flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .FALSE.   !: agrif flag
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_agrif
   !!----------------------------------------------------------------------
   !! NOT 'key_agrif'      dummy function                     No AGRIF zoom
   !!----------------------------------------------------------------------
   LOGICAL FUNCTION Agrif_Root()
      Agrif_Root = .TRUE.
   END FUNCTION Agrif_Root

   CHARACTER(len=3) FUNCTION Agrif_CFixed()
      Agrif_CFixed = '0' 
   END FUNCTION Agrif_CFixed
#endif

   INTEGER FUNCTION dom_oce_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(12) :: ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( rdttra(jpk), r2dtra(jpk), mig(jpi), mjg(jpj), nfiimpp(jpni,jpnj),  &
         &      nfipproc(jpni,jpnj), nfilcit(jpni,jpnj), STAT=ierr(1) )
         !
      ALLOCATE( nimppt(jpnij) , ibonit(jpnij) , nlcit(jpnij) , nlcjt(jpnij) ,     &
         &      njmppt(jpnij) , ibonjt(jpnij) , nldit(jpnij) , nldjt(jpnij) ,     &
         &                                      nleit(jpnij) , nlejt(jpnij) ,     &
         &      mi0(jpidta)   , mi1 (jpidta),  mj0(jpjdta)   , mj1 (jpjdta),      &
         &      tpol(jpiglo)  , fpol(jpiglo)                               , STAT=ierr(2) )
         !
      ALLOCATE( glamt(jpi,jpj) , gphit(jpi,jpj) , e1t(jpi,jpj) , e2t(jpi,jpj) , r1_e1t(jpi,jpj) , r1_e2t(jpi,jpj) ,   & 
         &      glamu(jpi,jpj) , gphiu(jpi,jpj) , e1u(jpi,jpj) , e2u(jpi,jpj) , r1_e1u(jpi,jpj) , r1_e2u(jpi,jpj) ,   &  
         &      glamv(jpi,jpj) , gphiv(jpi,jpj) , e1v(jpi,jpj) , e2v(jpi,jpj) , r1_e1v(jpi,jpj) , r1_e2v(jpi,jpj) ,   &  
         &      glamf(jpi,jpj) , gphif(jpi,jpj) , e1f(jpi,jpj) , e2f(jpi,jpj) , r1_e1f(jpi,jpj) , r1_e2f(jpi,jpj) ,   &
         &      e1e2t(jpi,jpj) , ff   (jpi,jpj) , STAT=ierr(3) )     
         !
      ALLOCATE( gdep3w_0(jpi,jpj,jpk) , e3v_0(jpi,jpj,jpk) , e3f_0 (jpi,jpj,jpk) ,                         &
         &      gdept_0 (jpi,jpj,jpk) , e3t_0(jpi,jpj,jpk) , e3u_0 (jpi,jpj,jpk) ,                         &
         &      gdepw_0 (jpi,jpj,jpk) , e3w_0(jpi,jpj,jpk) , e3vw_0(jpi,jpj,jpk) , e3uw_0(jpi,jpj,jpk) , STAT=ierr(4) )
         !
#if defined key_vvl
      ALLOCATE( gdep3w_n(jpi,jpj,jpk) , e3t_n (jpi,jpj,jpk) , e3u_n (jpi,jpj,jpk) ,                           &
         &      gdept_n (jpi,jpj,jpk) , e3v_n (jpi,jpj,jpk) , e3w_n (jpi,jpj,jpk) ,                           &
         &      gdepw_n (jpi,jpj,jpk) , e3f_n (jpi,jpj,jpk) , e3vw_n(jpi,jpj,jpk) , e3uw_n(jpi,jpj,jpk) ,     &
         &      e3t_b   (jpi,jpj,jpk) , e3u_b (jpi,jpj,jpk) , e3v_b (jpi,jpj,jpk) ,                           &
         &      e3uw_b  (jpi,jpj,jpk) , e3vw_b(jpi,jpj,jpk) ,                                                 &
         &      gdept_b (jpi,jpj,jpk) ,gdepw_b(jpi,jpj,jpk) , e3w_b (jpi,jpj,jpk) ,                           &
         &      e3t_a   (jpi,jpj,jpk) , e3u_a (jpi,jpj,jpk) , e3v_a (jpi,jpj,jpk) ,                           &
         &      ehu_a    (jpi,jpj)    , ehv_a  (jpi,jpj),                                                     &
         &      ehur_a   (jpi,jpj)    , ehvr_a (jpi,jpj),                                                     &
         &      ehu_b    (jpi,jpj)    , ehv_b  (jpi,jpj),                                                     &
         &      ehur_b   (jpi,jpj)    , ehvr_b (jpi,jpj),                                  STAT=ierr(5) )                          
#endif
         !
      ALLOCATE( hu      (jpi,jpj) , hur     (jpi,jpj) , hu_0(jpi,jpj) , ht_0  (jpi,jpj) ,     &
         &      hv      (jpi,jpj) , hvr     (jpi,jpj) , hv_0(jpi,jpj) , ht    (jpi,jpj) ,     &
         &      re2u_e1u(jpi,jpj) , re1v_e2v(jpi,jpj) ,                                       &
         &      e12t    (jpi,jpj) , r1_e12t (jpi,jpj) ,                                       &
         &      e12u    (jpi,jpj) , r1_e12u (jpi,jpj) ,                                       &
         &      e12v    (jpi,jpj) , r1_e12v (jpi,jpj) ,                                       &
         &      e12f    (jpi,jpj) , r1_e12f (jpi,jpj) ,                                   STAT=ierr(6)  )
         !
      ALLOCATE( gdept_1d(jpk) , gdepw_1d(jpk) ,                                     &
         &      e3t_1d  (jpk) , e3w_1d  (jpk) , e3tp (jpi,jpj), e3wp(jpi,jpj) ,     &
         &      gsigt   (jpk) , gsigw   (jpk) , gsi3w(jpk)    ,                     &
         &      esigt   (jpk) , esigw   (jpk)                                 , STAT=ierr(7) )
         !
      ALLOCATE( hbatv (jpi,jpj) , hbatf (jpi,jpj) ,     &
         &      hbatt (jpi,jpj) , hbatu (jpi,jpj) ,     &
         &      scosrf(jpi,jpj) , scobot(jpi,jpj) ,     &
         &      hifv  (jpi,jpj) , hiff  (jpi,jpj) ,     &
         &      hift  (jpi,jpj) , hifu  (jpi,jpj) , rx1 (jpi,jpj) , STAT=ierr(8) )

      ALLOCATE( mbathy(jpi,jpj) , bathy(jpi,jpj) ,                                      &
         &     tmask_i(jpi,jpj) , umask_i(jpi,jpj), vmask_i(jpi,jpj), fmask_i(jpi,jpj), &
         &     bmask(jpi,jpj)   ,                                                       &
         &     mbkt   (jpi,jpj) , mbku (jpi,jpj) , mbkv(jpi,jpj) , STAT=ierr(9) )

! (ISF) Allocation of basic array   
      ALLOCATE( misfdep(jpi,jpj) , risfdep(jpi,jpj),     &
         &     mikt(jpi,jpj), miku(jpi,jpj), mikv(jpi,jpj) ,           &
         &     mikf(jpi,jpj), ssmask(jpi,jpj), STAT=ierr(10) )

      ALLOCATE( tmask(jpi,jpj,jpk) , umask(jpi,jpj,jpk),     & 
         &      vmask(jpi,jpj,jpk) , fmask(jpi,jpj,jpk), STAT=ierr(11) )

      ALLOCATE( wmask(jpi,jpj,jpk) , wumask(jpi,jpj,jpk), wvmask(jpi,jpj,jpk) , STAT=ierr(12) )

#if defined key_noslip_accurate
      ALLOCATE( npcoa(4,jpk), nicoa(2*(jpi+jpj),4,jpk), njcoa(2*(jpi+jpj),4,jpk), STAT=ierr(12) )
#endif
      !
      dom_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION dom_oce_alloc

   !!======================================================================
END MODULE dom_oce

