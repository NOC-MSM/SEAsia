MODULE bdy_oce
   !!======================================================================
   !!                       ***  MODULE bdy_oce   ***
   !! Unstructured Open Boundary Cond. :   define related variables
   !!======================================================================
   !! History :  1.0  !  2001-05  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version     
   !!            3.3  !  2010-09  (D. Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2014-01  (C. Rousset) add ice boundary conditions for new model
   !!            4.0  !  2018     (C. Rousset) SI3 compatibility
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters

   IMPLICIT NONE
   PUBLIC

   INTEGER, PUBLIC, PARAMETER ::   jp_bdy  = 10       !: Maximum number of bdy sets
   INTEGER, PUBLIC, PARAMETER ::   jpbgrd  = 3        !: Number of horizontal grid types used  (T, U, V)

   TYPE, PUBLIC ::   OBC_INDEX    !: Indices and weights which define the open boundary
      INTEGER ,          DIMENSION(jpbgrd) ::  nblen
      INTEGER ,          DIMENSION(jpbgrd) ::  nblenrim
      INTEGER ,          DIMENSION(jpbgrd) ::  nblenrim0
      INTEGER , POINTER, DIMENSION(:,:)    ::  nbi
      INTEGER , POINTER, DIMENSION(:,:)    ::  nbj
      INTEGER , POINTER, DIMENSION(:,:)    ::  nbr
      INTEGER , POINTER, DIMENSION(:,:)    ::  nbmap
      INTEGER , POINTER, DIMENSION(:,:)    ::  ntreat
      REAL(wp), POINTER, DIMENSION(:,:)    ::  nbw
      REAL(wp), POINTER, DIMENSION(:,:)    ::  nbd
      REAL(wp), POINTER, DIMENSION(:,:)    ::  nbdout
      REAL(wp), POINTER, DIMENSION(:,:)    ::  flagu
      REAL(wp), POINTER, DIMENSION(:,:)    ::  flagv
   END TYPE OBC_INDEX

   !! Logicals in OBC_DATA structure are true if the chosen algorithm requires this
   !! field as external data. If true the data can come from external files
   !! or model initial conditions. If false then no "external" data array
   !! is required for this field. 

   TYPE, PUBLIC ::   OBC_DATA     !: Storage for external data
      INTEGER          , DIMENSION(2)   ::  nread
      LOGICAL                           ::  lneed_ssh
      !--- davbyr
      LOGICAL                           ::  lforced_ssh 
	  !--- END davbyr
      LOGICAL                           ::  lneed_dyn2d
      LOGICAL                           ::  lneed_dyn3d
      LOGICAL                           ::  lneed_tra
      LOGICAL                           ::  lneed_ice
      REAL(wp), POINTER, DIMENSION(:)   ::  ssh
      REAL(wp), POINTER, DIMENSION(:)   ::  u2d
      REAL(wp), POINTER, DIMENSION(:)   ::  v2d
      REAL(wp), POINTER, DIMENSION(:,:) ::  u3d
      REAL(wp), POINTER, DIMENSION(:,:) ::  v3d
      REAL(wp), POINTER, DIMENSION(:,:) ::  tem
      REAL(wp), POINTER, DIMENSION(:,:) ::  sal
      REAL(wp), POINTER, DIMENSION(:,:) ::  a_i    !: now ice leads fraction climatology
      REAL(wp), POINTER, DIMENSION(:,:) ::  h_i    !: Now ice  thickness climatology
      REAL(wp), POINTER, DIMENSION(:,:) ::  h_s    !: now snow thickness
      REAL(wp), POINTER, DIMENSION(:,:) ::  t_i    !: now ice  temperature
      REAL(wp), POINTER, DIMENSION(:,:) ::  t_s    !: now snow temperature
      REAL(wp), POINTER, DIMENSION(:,:) ::  tsu    !: now surf temperature
      REAL(wp), POINTER, DIMENSION(:,:) ::  s_i    !: now ice  salinity
      REAL(wp), POINTER, DIMENSION(:,:) ::  aip    !: now ice  pond concentration
      REAL(wp), POINTER, DIMENSION(:,:) ::  hip    !: now ice  pond depth
      REAL(wp), POINTER, DIMENSION(:,:) ::  hil    !: now ice  pond lid depth
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
   !                                                   !!** nambdy **
   LOGICAL, PUBLIC            ::   ln_bdy                   !: Unstructured Ocean Boundary Condition

   CHARACTER(len=80), DIMENSION(jp_bdy) ::   cn_coords_file !: Name of bdy coordinates file
   CHARACTER(len=80)                    ::   cn_mask_file   !: Name of bdy mask file
   !
   LOGICAL, DIMENSION(jp_bdy) ::   ln_coords_file           !: =T read bdy coordinates from file; 
   !                                                        !: =F read bdy coordinates from namelist
   LOGICAL                    ::   ln_mask_file             !: =T read bdymask from file
   LOGICAL                    ::   ln_vol                   !: =T volume correction             
   !
   INTEGER                    ::   nb_bdy                   !: number of open boundary sets
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
   LOGICAL , DIMENSION(jp_bdy) ::   ln_tra_dmp              !: =T Tracer damping
   LOGICAL , DIMENSION(jp_bdy) ::   ln_dyn3d_dmp            !: =T Baroclinic velocity damping
   REAL(wp), DIMENSION(jp_bdy) ::   rn_time_dmp             !: Damping time scale in days
   REAL(wp), DIMENSION(jp_bdy) ::   rn_time_dmp_out         !: Damping time scale in days at radiation outflow points

   CHARACTER(len=20), DIMENSION(jp_bdy) ::   cn_ice         ! Choice of boundary condition for sea ice variables 
   INTEGER , DIMENSION(jp_bdy)          ::   nn_ice_dta     !: = 0 use the initial state as bdy dta ; 
                                                            !: = 1 read it in a NetCDF file
   ! 
   !                                                   !!** nambdy_dta **
   REAL(wp), DIMENSION(jp_bdy) ::   rice_tem                !: temperature of incoming sea ice
   REAL(wp), DIMENSION(jp_bdy) ::   rice_sal                !: salinity    of incoming sea ice
   REAL(wp), DIMENSION(jp_bdy) ::   rice_age                !: age         of incoming sea ice
   REAL(wp), DIMENSION(jp_bdy) ::   rice_apnd               !: pond conc.  of incoming sea ice
   REAL(wp), DIMENSION(jp_bdy) ::   rice_hpnd               !: pond thick. of incoming sea ice
   REAL(wp), DIMENSION(jp_bdy) ::   rice_hlid               !: pond lid thick. of incoming sea ice
   
   !  davbyr 
   LOGICAL, DIMENSION(jp_bdy) ::   ln_ssh_bdy               !: =T USE SSH BDY - name list switch
   REAL(wp), DIMENSION(jp_bdy) ::  rn_ssh_shift             !: =F SHIFT SSH AT A BORDER BY rn_ssh_shift m_
   !  END davbyr
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
!$AGRIF_DO_NOT_TREAT
   TYPE(OBC_INDEX), DIMENSION(jp_bdy), TARGET      ::   idx_bdy           !: bdy indices (local process)
   TYPE(OBC_DATA) , DIMENSION(jp_bdy), TARGET      ::   dta_bdy           !: bdy external data (local process)
!$AGRIF_END_DO_NOT_TREAT
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lsend_bdy      !: mark needed communication for given boundary, grid and neighbour
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lrecv_bdy      !:  when searching in any direction
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lsend_bdyint   !: mark needed communication for given boundary, grid and neighbour
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lrecv_bdyint   !:  when searching towards the interior of the computational domain
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lsend_bdyext   !: mark needed communication for given boundary, grid and neighbour
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::   lrecv_bdyext   !:  when searching towards the exterior of the computational domain
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION bdy_oce_alloc()
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_stop, mpp_sum
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
      CALL mpp_sum ( 'bdy_oce', bdy_oce_alloc )
      IF( bdy_oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'bdy_oce_alloc: failed to allocate arrays.' )
      !
   END FUNCTION bdy_oce_alloc

   !!======================================================================
END MODULE bdy_oce

