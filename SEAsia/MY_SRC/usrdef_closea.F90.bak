MODULE usrdef_closea
   !!======================================================================
   !!                   ***  MODULE  usrdef_closea  ***
   !!
   !!                      ===  ORCA configuration  ===
   !!                         (2, 1 and 1/4 degrees)
   !!
   !! User define : specific treatments associated with closed seas
   !!======================================================================
   !! History :   8.2  !  2000-05  (O. Marti)  Original code
   !!   NEMO      1.0  !  2002-06  (E. Durand, G. Madec)  F90
   !!             3.0  !  2006-07  (G. Madec)  add clo_rnf, clo_ups, clo_bat
   !!             3.4  !  2014-12  (P.G. Fogli) sbc_clo bug fix & mpp reproducibility
   !!             4.0  !  2016-06  (G. Madec)  move to usrdef_closea, remove clo_ups
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_clo    : modification of the ocean domain for closed seas cases
   !!   sbc_clo    : Special handling of closed seas
   !!   clo_rnf    : set close sea outflows as river mouths (see sbcrnf)
   !!   clo_bat    : set to zero a field over closed sea (see domzrg)
   !!----------------------------------------------------------------------
   USE oce             ! dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! ocean surface boundary conditions
   !
   USE in_out_manager  ! I/O manager
   USE lib_fortran,    ONLY: glob_sum, DDPDD
   USE lbclnk          ! lateral boundary condition - MPP exchanges
   USE lib_mpp         ! MPP library
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dom_clo      ! called by domain module
   PUBLIC sbc_clo      ! called by step module
   PUBLIC clo_rnf      ! called by sbcrnf module
   PUBLIC clo_bat      ! called in domzgr module

   INTEGER, PUBLIC, PARAMETER          ::   jpncs   = 4      !: number of closed sea
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncstt            !: Type of closed sea
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsi1, ncsj1     !: south-west closed sea limits (i,j)
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsi2, ncsj2     !: north-east closed sea limits (i,j)
   INTEGER, PUBLIC, DIMENSION(jpncs)   ::   ncsnr            !: number of point where run-off pours
   INTEGER, PUBLIC, DIMENSION(jpncs,4) ::   ncsir, ncsjr     !: Location of runoff

   REAL(wp), DIMENSION (jpncs+1)       ::   surf             ! closed sea surface

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id: usrdef_closea.F90 7754 2017-03-03 11:51:06Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_clo( cd_cfg, kcfg )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_clo  ***
      !!        
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!                just the thermodynamic processes are applied.
      !!
      !! ** Action  :   ncsi1(), ncsj1() : south-west closed sea limits (i,j)
      !!                ncsi2(), ncsj2() : north-east Closed sea limits (i,j)
      !!                ncsir(), ncsjr() : Location of runoff
      !!                ncsnr            : number of point where run-off pours
      !!                ncstt            : Type of closed sea
      !!                                   =0 spread over the world ocean
      !!                                   =2 put at location runoff
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER         , INTENT(in   ) ::   kcfg     ! configuration identifier 
      !
      INTEGER ::   jc      ! dummy loop indices
      INTEGER ::   isrow   ! local index
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)'dom_clo : closed seas '
      IF(lwp) WRITE(numout,*)'~~~~~~~'
      !
      ! initial values
      ncsnr(:) = 1  ;  ncsi1(:) = 1  ;  ncsi2(:) = 1  ;  ncsir(:,:) = 1
      ncstt(:) = 0  ;  ncsj1(:) = 1  ;  ncsj2(:) = 1  ;  ncsjr(:,:) = 1
      !
      ! set the closed seas (in data domain indices)
      ! -------------------
      !
       ncsnr(1)   =   1  ;  ncstt(1)   =   0           ! spread over the globe
       ncsi1(1)   =   1200  ;  ncsj1(1)   =   2
       ncsi2(1)   =   1205  ;  ncsj2(1)   = 500
       ncsir(1,1) =   1  ;  ncsjr(1,1) =   1 

      ! convert the position in local domain indices
      ! --------------------------------------------
      DO jc = 1, jpncs
         ncsi1(jc)   = mi0( ncsi1(jc) )
         ncsj1(jc)   = mj0( ncsj1(jc) )
         !
         ncsi2(jc)   = mi1( ncsi2(jc) )   
         ncsj2(jc)   = mj1( ncsj2(jc) )  
      END DO
      !
   END SUBROUTINE dom_clo


   SUBROUTINE sbc_clo( kt, cd_cfg, kcfg )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_clo  ***
      !!                    
      !! ** Purpose :   Special handling of closed seas
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action  :   emp updated surface freshwater fluxes and associated heat content at kt
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kt       ! ocean model time step
      CHARACTER(len=*), INTENT(in   ) ::   cd_cfg   ! configuration name
      INTEGER         , INTENT(in   ) ::   kcfg     ! configuration identifier 
      !
      INTEGER             ::   ji, jj, jc, jn   ! dummy loop indices
      REAL(wp), PARAMETER ::   rsmall = 1.e-20_wp    ! Closed sea correction epsilon
      REAL(wp)            ::   zze2, ztmp, zcorr     ! 
      REAL(wp)            ::   zcoef, zcoef1         ! 
      COMPLEX(wp)         ::   ctmp 
      REAL(wp), DIMENSION(jpncs) ::   zfwf   ! 1D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_clo')
      !                                                   !------------------!
      IF( kt == nit000 ) THEN                             !  Initialisation  !
         !                                                !------------------!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'sbc_clo : closed seas '
         IF(lwp) WRITE(numout,*)'~~~~~~~'
 
      END IF

      IF( nn_timing == 1 )  CALL timing_stop('sbc_clo')
      !
   END SUBROUTINE sbc_clo


   SUBROUTINE clo_rnf( p_rnfmsk )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!                    
      !! ** Purpose :   allow the treatment of closed sea outflow grid-points
      !!                to be the same as river mouth grid-points
      !!
      !! ** Method  :   set to 1 the runoff mask (mskrnf, see sbcrnf module)
      !!                at the closed sea outflow grid-point.
      !!
      !! ** Action  :   update (p_)mskrnf (set 1 at closed sea outflow)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   p_rnfmsk   ! river runoff mask (rnfmsk array)
      !
      INTEGER  ::   jc, jn, ji, jj      ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jc = 1, jpncs
         IF( ncstt(jc) >= 1 ) THEN            ! runoff mask set to 1 at closed sea outflows
             DO jn = 1, 4
                DO jj =    mj0( ncsjr(jc,jn) ), mj1( ncsjr(jc,jn) )
                   DO ji = mi0( ncsir(jc,jn) ), mi1( ncsir(jc,jn) )
                      p_rnfmsk(ji,jj) = p_rnfmsk(ji,jj)
                   END DO
                END DO
            END DO 
         ENDIF 
      END DO 
      !
   END SUBROUTINE clo_rnf
   
      
   SUBROUTINE clo_bat( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE clo_bat  ***
      !!                    
      !! ** Purpose :   suppress closed sea from the domain
      !!
      !! ** Method  :   set first and last ocean level to 0 over the closed seas.
      !!
      !! ** Action  :   set pbat=0 and kbat=0 over closed seas
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(inout) ::   k_top, k_bot   ! ocean first and last level indices
      !
      INTEGER  ::   jc, ji, jj      ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jc = 1, jpncs
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               k_top(ji,jj) = 0   
               k_bot(ji,jj) = 0   
            END DO 
         END DO 
       END DO 
       !
   END SUBROUTINE clo_bat

   !!======================================================================
END MODULE usrdef_closea

