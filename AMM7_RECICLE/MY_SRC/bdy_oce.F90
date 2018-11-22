MODULE bdy_oce
   !!======================================================================
   !!                       ***  MODULE bdy_oce   ***
   !! Unstructured Open Boundary Cond. :   define related variables
   !!======================================================================
   !! History :  1.0  !  2001-05  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version     
   !!            3.3  !  2010-09  (D. Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2012-01  (C. Rousset) add ice boundary conditions for lim3
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy'                      Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters
   USE bdy_par         ! Unstructured boundary parameters
   USE lib_mpp         ! distributed memory computing

   IMPLICIT NONE
   PUBLIC

   TYPE, PUBLIC ::   OBC_INDEX    !: Indices and weights which define the open boundary
      INTEGER,          DIMENSION(jpbgrd) ::  nblen
      INTEGER,          DIMENSION(jpbgrd) ::  nblenrim
      INTEGER, POINTER, DIMENSION(:,:)   ::  nbi
      INTEGER, POINTER, DIMENSION(:,:)   ::  nbj
      INTEGER, POINTER, DIMENSION(:,:)   ::  nbr
      INTEGER, POINTER, DIMENSION(:,:)   ::  nbmap
      REAL(wp)   , POINTER, DIMENSION(:,:)   ::  nbw
      REAL(wp)   , POINTER, DIMENSION(:,:)   ::  nbd
      REAL(wp)   , POINTER, DIMENSION(:,:)   ::  nbdout
      REAL(wp)   , POINTER, DIMENSION(:,:)   ::  flagu
      REAL(wp)   , POINTER, DIMENSION(:,:)   ::  flagv
   END TYPE OBC_INDEX

   !! Logicals in OBC_DATA structure are true if the chosen algorithm requires this
   !! field as external data. If true the data can come from external files
   !! or model initial conditions. If false then no "external" data array
   !! is required for this field. 

   TYPE, PUBLIC ::   OBC_DATA     !: Storage for external data
      INTEGER,       DIMENSION(2)     ::  nread
      LOGICAL                         ::  ll_ssh
      LOGICAL                         ::  ll_u2d
      LOGICAL                         ::  ll_v2d
      LOGICAL                         ::  ll_u3d
      LOGICAL                         ::  ll_v3d
      LOGICAL                         ::  ll_tem
      LOGICAL                         ::  ll_sal
      LOGICAL                           ::  ll_fvl
      REAL(wp), POINTER, DIMENSION(:)     ::  ssh
      REAL(wp), POINTER, DIMENSION(:)     ::  u2d
      REAL(wp), POINTER, DIMENSION(:)     ::  v2d
      REAL(wp), POINTER, DIMENSION(:,:)   ::  u3d
      REAL(wp), POINTER, DIMENSION(:,:)   ::  v3d
      REAL(wp), POINTER, DIMENSION(:,:)   ::  tem
      REAL(wp), POINTER, DIMENSION(:,:)   ::  sal
#if defined key_lim2
      LOGICAL                         ::  ll_frld
      LOGICAL                         ::  ll_hicif
      LOGICAL                         ::  ll_hsnif
      REAL(wp), POINTER, DIMENSION(:)     ::  frld
      REAL(wp), POINTER, DIMENSION(:)     ::  hicif
      REAL(wp), POINTER, DIMENSION(:)     ::  hsnif
#elif defined key_lim3
      LOGICAL                         ::  ll_a_i
      LOGICAL                         ::  ll_ht_i
      LOGICAL                         ::  ll_ht_s
      REAL, POINTER, DIMENSION(:,:)   ::  a_i   !: now ice leads fraction climatology
      REAL, POINTER, DIMENSION(:,:)   ::  ht_i  !: Now ice  thickness climatology
      REAL, POINTER, DIMENSION(:,:)   ::  ht_s  !: now snow thickness
#endif
#if defined key_top
      CHARACTER(LEN=20)                   :: cn_obc  !: type of boundary condition to apply
      REAL(wp)                            :: rn_fac  !: multiplicative scaling factor
      REAL(wp), POINTER, DIMENSION(:,:)   :: trc     !: now field of the tracer
      LOGICAL                             :: dmp     !: obc damping term
#endif

   END TYPE OBC_DATA

   !!----------------------------------------------------------------------
   !! Namelist variables
   !!----------------------------------------------------------------------
   CHARACTER(len=80), DIMENSION(jp_bdy) ::   cn_coords_file !: Name of bdy coordinates file
   CHARACTER(len=80)                    ::   cn_mask_file   !: Name of bdy mask file
   !
   LOGICAL, DIMENSION(jp_bdy) ::   ln_coords_file           !: =T read bdy coordinates from file; 
   !                                                        !: =F read bdy coordinates from namelist
   LOGICAL                    ::   ln_mask_file             !: =T read bdymask from file
   LOGICAL                    ::   ln_vol                   !: =T volume correction             
   !
   INTEGER                    ::   nb_bdy                   !: number of open boundary sets
   INTEGER, DIMENSION(jp_bdy) ::   nb_jpk_bdy               !: number of levels in the bdy data (set < 0 if consistent with planned run)
   INTEGER, DIMENSION(jp_bdy) ::   nn_rimwidth              !: boundary rim width for Flow Relaxation Scheme
   INTEGER                    ::   nn_volctl                !: = 0 the total volume will have the variability of the surface Flux E-P 
   !                                                        !  = 1 the volume will be constant during all the integration.
   CHARACTER(len=20), DIMENSION(jp_bdy) ::   cn_dyn2d       ! Choice of boundary condition for barotropic variables (U,V,SSH)
   INTEGER, DIMENSION(jp_bdy)           ::   nn_dyn2d_dta   !: = 0 use the initial state as bdy dta ; 
                                                            !: = 1 read it in a NetCDF file
                                                            !: = 2 read tidal harmonic forcing from a NetCDF file
                                                            !: = 3 read external data AND tidal harmonic forcing from NetCDF files
   CHARACTER(len=20), DIMENSION(jp_bdy) ::   cn_dyn3d       ! Choice of boundary condition for baroclinic velocities 
   INTEGER, DIMENSION(jp_bdy)           ::   nn_dyn3d_dta   !: = 0 use the initial state as bdy dta ; 
                                                            !: = 1 read it in a NetCDF file
   CHARACTER(len=20), DIMENSION(jp_bdy) ::   cn_tra         ! Choice of boundary condition for active tracers (T and S)
   INTEGER, DIMENSION(jp_bdy)           ::   nn_tra_dta     !: = 0 use the initial state as bdy dta ; 
                                                            !: = 1 read it in a NetCDF file
   LOGICAL, DIMENSION(jp_bdy) ::   ln_tra_dmp               !: =T Tracer damping
   LOGICAL, DIMENSION(jp_bdy) ::   ln_dyn3d_dmp             !: =T Baroclinic velocity damping
   REAL(wp),    DIMENSION(jp_bdy) ::   rn_time_dmp              !: Damping time scale in days
   REAL(wp),    DIMENSION(jp_bdy) ::   rn_time_dmp_out          !: Damping time scale in days at radiation outflow points

   CHARACTER(len=20), DIMENSION(jp_bdy) ::   cn_ice_lim       ! Choice of boundary condition for sea ice variables 
   INTEGER, DIMENSION(jp_bdy)           ::   nn_ice_lim_dta   !: = 0 use the initial state as bdy dta ; 
                                                              !: = 1 read it in a NetCDF file
   REAL(wp),    DIMENSION(jp_bdy) ::   rn_ice_tem             !: choice of the temperature of incoming sea ice
   REAL(wp),    DIMENSION(jp_bdy) ::   rn_ice_sal             !: choice of the salinity    of incoming sea ice
   REAL(wp),    DIMENSION(jp_bdy) ::   rn_ice_age             !: choice of the age         of incoming sea ice
   !
   
   !!----------------------------------------------------------------------
   !! Global variables
   !!----------------------------------------------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), TARGET ::   bdytmask   !: Mask defining computational domain at T-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), TARGET ::   bdyumask   !: Mask defining computational domain at U-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), TARGET ::   bdyvmask   !: Mask defining computational domain at V-points

   REAL(wp)                                    ::   bdysurftot !: Lateral surface of unstructured open boundary

   !!----------------------------------------------------------------------
   !! open boundary data variables
   !!----------------------------------------------------------------------

   INTEGER,  DIMENSION(jp_bdy)                     ::   nn_dta            !: =0 => *all* data is set to initial conditions
                                                                          !: =1 => some data to be read in from data files
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global        !: workspace for reading in global data arrays (unstr.  bdy)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global_z      !: workspace for reading in global depth arrays (unstr.  bdy)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global_dz     !: workspace for reading in global depth arrays (unstr.  bdy)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global2       !: workspace for reading in global data arrays (struct. bdy)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global2_z     !: workspace for reading in global depth arrays (struct. bdy)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   dta_global2_dz    !: workspace for reading in global depth arrays (struct. bdy)
!$AGRIF_DO_NOT_TREAT
   TYPE(OBC_INDEX), DIMENSION(jp_bdy), TARGET      ::   idx_bdy           !: bdy indices (local process)
   TYPE(OBC_DATA) , DIMENSION(jp_bdy), TARGET      ::   dta_bdy           !: bdy external data (local process)
!$AGRIF_END_DO_NOT_TREAT
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION bdy_oce_alloc()
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn, mpp_sum
      !
      INTEGER :: bdy_oce_alloc
      !!----------------------------------------------------------------------
      !
      ALLOCATE( bdytmask(jpi,jpj) , bdyumask(jpi,jpj), bdyvmask(jpi,jpj),     &  
         &      STAT=bdy_oce_alloc )
      !
      ! Initialize masks 
      bdytmask(:,:) = 1._wp
      bdyumask(:,:) = 1._wp
      bdyvmask(:,:) = 1._wp
      ! 
      IF( lk_mpp             )   CALL mpp_sum ( bdy_oce_alloc )
      IF( bdy_oce_alloc /= 0 )   CALL ctl_warn('bdy_oce_alloc: failed to allocate arrays.')
      !
   END FUNCTION bdy_oce_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                NO Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   LOGICAL ::   ln_tides = .false.  !: =T apply tidal harmonic forcing along open boundaries
#endif

   !!======================================================================
END MODULE bdy_oce

