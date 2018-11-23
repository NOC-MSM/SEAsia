MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE usrdef_sbc  ***
   !! 
   !!                  ===  AMM7_SURGE configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!            4.0   ! 2017-12  (C. O'Neill)  add necessary options for surge work - either no fluxes 
   !!                                           (for tide-only run) or wind and pressure only
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in LOCK_EXCHANGE case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ocean fields
   USE fldread         ! read input fields
   USE phycst          ! physical constants
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 
   !
   USE in_out_manager  ! I/O manager
   USE iom
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo       ! work arrays
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce    ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by sbcice_lim.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by sbcice_lim.F90 for ice thermo
   PUBLIC   surge_oce    ! routine called by usrdef_sbc_oce (if required)
   
   
   INTEGER , PARAMETER ::   jpfld   = 2           ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_wndi = 1           ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2           ! index of 10m wind velocity (j-component) (m/s)    at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)

   REAL(wp), PARAMETER ::   rhoa =    1.22        ! air density

   !                                  !!* Namelist namsbc_usr
   REAL(wp) ::   rn_vfac     ! multiplication factor for ice/ocean velocity in the calculation of wind stress (clem)
   REAL(wp) ::   rn_charn_const 
   LOGICAL  ::   ln_use_sbc  ! Surface fluxes on or not 


   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usr_def_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   all 0 fields, for AMM7_SURGE case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - if tide-only case - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!              - if tide+surge case - read in wind and air pressure      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step

      INTEGER  ::   ierror   ! return error code
      INTEGER  ::   ifpr     ! dummy loop indice
      INTEGER  ::   ios      ! Local integer output status for namelist read
      !
      CHARACTER(len=100) ::  cn_dir   !   Root directory for location of flux files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj                        ! informations about the fields to be read

      NAMELIST/namsbc_usr/ ln_use_sbc, cn_dir , rn_vfac,  &
         &                   sn_wndi, sn_wndj, rn_charn_const
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         
         
         REWIND( numnam_cfg )              ! Namelist namsbc_usr in configuration namelist
         READ  ( numnam_cfg, namsbc_usr, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_surge in configuration namelist', lwp )

         IF(lwm) WRITE( numond, namsbc_usr )

         IF(ln_use_sbc) THEN
             IF(lwp) WRITE(numout,*)' usr_sbc : AMM7_SURGE tide + surge case: surface wind and pressure (assuming ln_dyn_apr=T) applied'
             IF(lwp) WRITE(numout,*)' ~~~~~~~~~~~ '

             !                                         ! store namelist information in an array
             slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
             !
             !
             ALLOCATE( sf(jpfld), STAT=ierror )         ! set sf structure
             IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_surge: unable to allocate sf structure' )
             DO ifpr= 1, jpfld
                ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
                IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
             END DO
             !                                         ! fill sf with slf_i and control print
             CALL fld_fill( sf, slf_i, cn_dir, 'sbc_surge', 'flux formulation for ocean surface boundary condition', 'namsbc_surge' )

         ELSE
             IF(lwp) WRITE(numout,*)' usr_sbc : AMM7_SURGE tide only case: NO surface forcing'
             IF(lwp) WRITE(numout,*)' ~~~~~~~~~~~   utau = vtau = taum = wndm = qns = qsr = emp = sfx = 0'

             utau(:,:) = 0._wp
             vtau(:,:) = 0._wp
             taum(:,:) = 0._wp
             wndm(:,:) = 0._wp
             !
             emp (:,:) = 0._wp
             sfx (:,:) = 0._wp
             qns (:,:) = 0._wp
             qsr (:,:) = 0._wp
             !         
             uwnd(:,:) = 0._wp
             vwnd(:,:) = 0._wp
         ENDIF

      ENDIF
      !
      IF(ln_use_sbc) THEN
          CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step

          !                                            ! compute the surface ocean fluxes using CORE bulk formulea
          IF( MOD( kt - 1, nn_fsbc ) == 0 )   CALL surge_oce( kt, sf, ssu_m, ssv_m, rn_charn_const )
         
      ENDIF
   END SUBROUTINE usrdef_sbc_oce


    
   SUBROUTINE surge_oce( kt, sf, pu, pv, rn_charn_const )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE surge_oce  ***
      !!
      !! ** Purpose :   provide the momentum fluxes at
      !!      the ocean surface at each time step
      !!
      !! ** Method  :   Charnock formulea for the ocean using atmospheric
      !!      fields read in sbc_read
      !! 
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!
      !!  ** Nota  :   sf has to be a dummy argument for AGRIF on NEC
      !!---------------------------------------------------------------------
      INTEGER  , INTENT(in   )                 ::   kt    ! time step index
      TYPE(fld), INTENT(inout), DIMENSION(:)   ::   sf    ! input data
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pu    ! surface current at U-point (i-component) [m/s]
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pv    ! surface current at V-point (j-component) [m/s]
      REAL(wp) , INTENT(in)                    ::   rn_charn_const! local variable
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp                ! local variable
      REAL(wp) ::   z_z0, z_Cd1          ! local variable
      REAL(wp) ::   zi                   ! local variable
      REAL(wp), DIMENSION(:,:), POINTER ::   zwnd_i, zwnd_j    ! wind speed components at T-point
      REAL(wp), DIMENSION(:,:), POINTER ::   Cd                ! transfer coefficient for momentum      (tau)
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('surge_oce')
      !
      CALL wrk_alloc( jpi,jpj, zwnd_i, zwnd_j )
      CALL wrk_alloc( jpi,jpj, Cd )
      !
      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
      zwnd_i(:,:) = 0.e0  
      zwnd_j(:,:) = 0.e0
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            uwnd(ji,jj) = (  sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pu(ji-1,jj  ) + pu(ji,jj) )  )
            vwnd(ji,jj) = (  sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pv(ji  ,jj-1) + pv(ji,jj) )  )
         END DO
      END DO
      zwnd_i(:,:) = uwnd(:,:)
      zwnd_j(:,:) = vwnd(:,:)

      CALL lbc_lnk( zwnd_i(:,:) , 'T', -1. )
      CALL lbc_lnk( zwnd_j(:,:) , 'T', -1. )
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
      wndm(:,:) = SQRT(  zwnd_i(:,:) * zwnd_i(:,:)   &
         &             + zwnd_j(:,:) * zwnd_j(:,:)  ) * tmask(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      I   Radiative FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !
      
      qsr(:,:)=0._wp

      ! ----------------------------------------------------------------------------- !
      !     II    Turbulent FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !
      Cd(:,:)=0.0001_wp
      DO jj = 1,jpj
         DO ji = 1,jpi
            z_Cd1=0._wp
            zi=1
            !Iterate
            DO WHILE((abs(Cd(ji,jj)-z_Cd1))>1E-6)
               z_Cd1=Cd(ji,jj)
               z_z0=rn_charn_const*z_Cd1*wndm(ji,jj)**2/grav
               Cd(ji,jj)=(0.41_wp/log(10._wp/z_z0))**2
               zi=zi+1
            ENDDO
         ENDDO
      ENDDO

      ! ... tau module, i and j component
      DO jj = 1, jpj
         DO ji = 1, jpi
            zztmp = rhoa * wndm(ji,jj) * Cd(ji,jj)
            taum  (ji,jj) = zztmp * wndm  (ji,jj)
            zwnd_i(ji,jj) = zztmp * zwnd_i(ji,jj)
            zwnd_j(ji,jj) = zztmp * zwnd_j(ji,jj)
         END DO
      END DO

      CALL iom_put( "taum_oce", taum )   ! output wind stress module
      CALL iom_put( "uwnd", uwnd )   ! output wind stress module
      CALL iom_put( "vwnd", vwnd )   ! output wind stress module

      ! ... utau, vtau at U- and V_points, resp.
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      !     Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
         END DO
      END DO
      CALL lbc_lnk( utau(:,:), 'U', -1. )
      CALL lbc_lnk( vtau(:,:), 'V', -1. )

    
      IF(ln_ctl) THEN
         CALL prt_ctl( tab2d_1=utau  , clinfo1=' surge_oce: utau   : ', mask1=umask,   &
            &          tab2d_2=vtau  , clinfo2=           ' vtau : '  , mask2=vmask )
         CALL prt_ctl( tab2d_1=wndm  , clinfo1=' surge_oce: wndm   : ')
      ENDIF
       
      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      emp (:,:) = 0._wp
      qns(:,:)  = 0._wp
      sfx(:,:) = 0._wp                          ! salt flux; zero unless ice is present (computed in limsbc(_2).F90)
      !
!       IF ( nn_ice == 0 ) THEN
!          CALL iom_put( "qns_oce" ,   qns  )                 ! output downward non solar heat over the ocean
!          CALL iom_put( "qsr_oce" ,   qsr  )                 ! output downward solar heat over the ocean
!          CALL iom_put( "qt_oce"  ,   qns+qsr )              ! output total downward heat over the ocean
!       ENDIF
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=utau , clinfo1=' surge_oce: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau , clinfo2=           ' vtau  : ' , mask2=vmask )
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, zwnd_i, zwnd_j )
      CALL wrk_dealloc( jpi,jpj, Cd )
      !
      IF( nn_timing == 1 )  CALL timing_stop('surge_oce')

      !
   END SUBROUTINE surge_oce


   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau

   SUBROUTINE usrdef_sbc_ice_flx( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
