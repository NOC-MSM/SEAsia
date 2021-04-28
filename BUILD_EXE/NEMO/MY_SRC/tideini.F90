MODULE tideini
   !!======================================================================
   !!                       ***  MODULE  tideini  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE daymod         ! calendar
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
   LOGICAL , PUBLIC ::   ln_read_load    !:
   LOGICAL , PUBLIC ::   ln_scal_load    !:
   LOGICAL , PUBLIC ::   ln_tide_ramp    !:
   INTEGER , PUBLIC ::   nb_harmo        !:
   INTEGER , PUBLIC ::   kt_tide         !:
   REAL(wp), PUBLIC ::   rdttideramp     !:
   REAL(wp), PUBLIC ::   rn_scal_load    !:
   ! davbyr - read love number from namelist
   REAL(wp), PUBLIC ::   dn_love_number  !:
   ! END davbyr
   CHARACTER(lc), PUBLIC ::   cn_tide_load   !: 

   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntide   !:

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
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
      NAMELIST/nam_tide/ln_tide, ln_tide_pot, ln_scal_load, ln_read_load, cn_tide_load, &
                  &     ln_tide_ramp, rn_scal_load, rdttideramp, dn_love_number, clname
      !!----------------------------------------------------------------------
      !
      ! Read Namelist nam_tide
      REWIND( numnam_ref )              ! Namelist nam_tide in reference namelist : Tides
      READ  ( numnam_ref, nam_tide, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_tide in reference namelist' )
      !
      REWIND( numnam_cfg )              ! Namelist nam_tide in configuration namelist : Tides
      READ  ( numnam_cfg, nam_tide, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_tide in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_tide )
      !
      IF( ln_tide ) THEN
         IF (lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'tide_init : Initialization of the tidal components'
            WRITE(numout,*) '~~~~~~~~~ '
            WRITE(numout,*) '   Namelist nam_tide'
            WRITE(numout,*) '      Use tidal components                       ln_tide      = ', ln_tide
            WRITE(numout,*) '         Apply astronomical potential            ln_tide_pot  = ', ln_tide_pot
            WRITE(numout,*) '         Use scalar approx. for load potential   ln_scal_load = ', ln_scal_load
            WRITE(numout,*) '         Read load potential from file           ln_read_load = ', ln_read_load
            WRITE(numout,*) '         Apply ramp on tides at startup          ln_tide_ramp = ', ln_tide_ramp
            WRITE(numout,*) '         Fraction of SSH used in scal. approx.   rn_scal_load = ', rn_scal_load
            WRITE(numout,*) '         Duration (days) of ramp                 rdttideramp  = ', rdttideramp
            ! davbyr - Love number (one line)
            WRITE(numout,*) '         Love Number                             dn_love_number = ', dn_love_number
         ENDIF
      ELSE
         rn_scal_load = 0._wp 

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
      !       
      ! Ensure that tidal components have been set in namelist_cfg
      IF( nb_harmo == 0 )   CALL ctl_stop( 'tide_init : No tidal components set in nam_tide' )
      !
      IF( ln_read_load.AND.(.NOT.ln_tide_pot) ) &
          &   CALL ctl_stop('ln_read_load requires ln_tide_pot')
      IF( ln_scal_load.AND.(.NOT.ln_tide_pot) ) &
          &   CALL ctl_stop('ln_scal_load requires ln_tide_pot')
      IF( ln_scal_load.AND.ln_read_load ) &
          &   CALL ctl_stop('Choose between ln_scal_load and ln_read_load')
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
      IF (.NOT.ln_scal_load ) rn_scal_load = 0._wp
      !
   END SUBROUTINE tide_init
     
   !!======================================================================
END MODULE tideini
