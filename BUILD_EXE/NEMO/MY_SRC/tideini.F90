MODULE tideini
   !!======================================================================
   !!                       ***  MODULE  tideini  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE daymod         ! calandar
   USE tide_mod       ! 
   !
   USE in_out_manager ! I/O units
   USE iom            ! xIOs server
   USE ioipsl         ! NetCDF IPSL library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PUBLIC

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   omega_tide   !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   v0tide       !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   utide        !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ftide        !:

   LOGICAL , PUBLIC ::   ln_tide         !:
   LOGICAL , PUBLIC ::   ln_tide_pot     !:
   LOGICAL , PUBLIC ::   ln_tide_ramp    !:
   INTEGER , PUBLIC ::   nb_harmo        !:
   INTEGER , PUBLIC ::   kt_tide         !:
   REAL(wp), PUBLIC ::   rdttideramp     !:
   ! NB - read love number from namelist
   REAL(wp), PUBLIC ::   dn_love_number  !:
   ! END NB
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntide   !:

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id: tideini.F90 7646 2017-02-06 09:25:03Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE tide_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init  ***
      !!----------------------------------------------------------------------      
      INTEGER  :: ji, jk
      CHARACTER(LEN=4), DIMENSION(jpmax_harmo) :: clname
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !
      ! NB - read love number from namelist
      !NAMELIST/nam_tide/ln_tide, ln_tide_pot, ln_tide_ramp, rdttideramp, clname
      NAMELIST/nam_tide/ln_tide, ln_tide_pot, ln_tide_ramp, rdttideramp, dn_love_number, clname
      ! END NB
      !!----------------------------------------------------------------------
      !
      ! Read Namelist nam_tide
      REWIND( numnam_ref )              ! Namelist nam_tide in reference namelist : Tides
      READ  ( numnam_ref, nam_tide, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_tide in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist nam_tide in configuration namelist : Tides
      READ  ( numnam_cfg, nam_tide, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_tide in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_tide )
      !
      IF (ln_tide) THEN
         IF (lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'tide_init : Initialization of the tidal components'
            WRITE(numout,*) '~~~~~~~~~ '
            WRITE(numout,*) '   Namelist nam_tide'
            WRITE(numout,*) '              Use tidal components : ln_tide      = ', ln_tide
            WRITE(numout,*) '      Apply astronomical potential : ln_tide_pot  = ', ln_tide_pot
!            WRITE(numout,*) '                                     nb_harmo     = ', nb_harmo
            WRITE(numout,*) '                                     ln_tide_ramp = ', ln_tide_ramp
! NB - Love number
            WRITE(numout,*) '                                     dn_love_number = ', dn_love_number
! End NB
         ENDIF
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tide_init : tidal components not used (ln_tide = F)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
         RETURN
      ENDIF
      !
      CALL tide_init_Wave
      !
      nb_harmo=0
      DO jk = 1, jpmax_harmo
         DO ji = 1,jpmax_harmo
            IF( TRIM(clname(jk)) == Wave(ji)%cname_tide )   nb_harmo = nb_harmo + 1
         END DO
      END DO
      IF (ln_tide .and.lwp) WRITE(numout,*) '                                     nb_harmo     = ', nb_harmo

      ! Ensure that tidal components have been set in namelist_cfg
      IF( nb_harmo == 0 )   CALL ctl_stop( 'tide_init : No tidal components set in nam_tide' )
      !
      IF( ln_tide_ramp.AND.((nitend-nit000+1)*rdt/rday < rdttideramp) )   &
         &   CALL ctl_stop('rdttideramp must be lower than run duration')
      IF( ln_tide_ramp.AND.(rdttideramp<0.) ) &
         &   CALL ctl_stop('rdttideramp must be positive')
      !
      ALLOCATE( ntide(nb_harmo) )
      DO jk = 1, nb_harmo
         DO ji = 1, jpmax_harmo
            IF( TRIM(clname(jk)) == Wave(ji)%cname_tide ) THEN
               ntide(jk) = ji
               EXIT
            ENDIF
         END DO
      END DO
      !
      ALLOCATE( omega_tide(nb_harmo), v0tide    (nb_harmo),   &
         &      utide     (nb_harmo), ftide     (nb_harmo)  )
      kt_tide = nit000
      !
   END SUBROUTINE tide_init
     
   !!======================================================================
END MODULE tideini
