MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of kbat and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.6  ! 2015-05  (P. Mathiot) ISF: add wmask,wumask and wvmask
   !!            4.0  ! 2016-06  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk       : compute land/ocean mask
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE usrdef_fmask   ! user defined fmask
   USE bdy_oce      
   USE in_out_manager ! I/O manager
   USE iom
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! Massively Parallel Processing library
   USE wrk_nemo       ! Memory allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk    ! routine called by inidom.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
   REAL(wp)        :: rn_shlat   ! type of lateral boundary condition on velocity
   LOGICAL, PUBLIC :: ln_vorlat  !  consistency of vorticity boundary condition 
   !                                            with analytical eqs.

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: dommsk.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_msk( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) points.
      !!
      !! ** Method  :   The ocean/land mask  at t-point is deduced from ko_top 
      !!      and ko_bot, the indices of the fist and last ocean t-levels which 
      !!      are either defined in usrdef_zgr or read in zgr_read.
      !!                The velocity masks (umask, vmask, wmask, wumask, wvmask) 
      !!      are deduced from a product of the two neighboring tmask.
      !!                The vorticity mask (fmask) is deduced from tmask taking
      !!      into account the choice of lateral boundary condition (rn_shlat) :
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !!      tmask_i : interior ocean mask at t-point, i.e. excluding duplicated
      !!                rows/lines due to cyclic or North Fold boundaries as well
      !!                as MPP halos.
      !!      tmask_h : halo mask at t-point, i.e. excluding duplicated rows/lines
      !!                due to cyclic or North Fold boundaries as well as MPP halos.
      !!
      !! ** Action :   tmask, umask, vmask, wmask, wumask, wvmask : land/ocean mask 
      !!                         at t-, u-, v- w, wu-, and wv-points (=0. or 1.)
      !!               fmask   : land/ocean mask at f-point (=0., or =1., or 
      !!                         =rn_shlat along lateral boundaries)
      !!               tmask_i : interior ocean mask 
      !!               tmask_h : halo mask
      !!               ssmask , ssumask, ssvmask, ssfmask : 2D ocean mask
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(in) ::   k_top, k_bot   ! first and last ocean level
      !
      INTEGER  ::   ji, jj, jk     ! dummy loop indices
      INTEGER  ::   iif, iil       ! local integers
      INTEGER  ::   ijf, ijl       !   -       -
      INTEGER  ::   iktop, ikbot   !   -       -
      INTEGER  ::   ios, inum
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwf   ! 2D workspace
      !!
      NAMELIST/namlbc/ rn_shlat, ln_vorlat
      NAMELIST/nambdy/ ln_bdy ,nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice_lim, nn_ice_lim_dta,                           &
         &             rn_ice_tem, rn_ice_sal, rn_ice_age,                 &
         &             ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_msk')
      !
      REWIND( numnam_ref )              ! Namelist namlbc in reference namelist : Lateral momentum boundary condition
      READ  ( numnam_ref, namlbc, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namlbc in configuration namelist : Lateral momentum boundary condition
      READ  ( numnam_cfg, namlbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namlbc )
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat 
      ENDIF

      IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  free-slip '
      ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  no-slip '
      ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  partial-slip '
      ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  strong-slip '
      ELSE
         WRITE(ctmp1,*) ' rn_shlat is negative = ', rn_shlat
         CALL ctl_stop( ctmp1 )
      ENDIF


      !  Ocean/land mask at t-point  (computed from ko_top and ko_bot)
      ! ----------------------------
      !
      tmask(:,:,:) = 0._wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            iktop = k_top(ji,jj)
            ikbot = k_bot(ji,jj)
            IF( iktop /= 0 ) THEN       ! water in the column
               tmask(ji,jj,iktop:ikbot  ) = 1._wp
            ENDIF
         END DO  
      END DO  
!SF  add here lbc_lnk: bug not still understood : cause now domain configuration is read !
!!gm I don't understand why...  
      CALL lbc_lnk( tmask  , 'T', 1._wp )      ! Lateral boundary conditions

     ! Mask corrections for bdy (read in mppini2)
      REWIND( numnam_ref )              ! Namelist nambdy in reference namelist :Unstructured open boundaries
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nambdy in configuration namelist :Unstructured open boundaries
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy in configuration namelist', lwp )
      ! ------------------------
      IF ( ln_bdy .AND. ln_mask_file ) THEN
         CALL iom_open( cn_mask_file, inum )
         CALL iom_get ( inum, jpdom_data, 'bdy_msk', bdytmask(:,:) )
         CALL iom_close( inum )
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  tmask(ji,jj,jk) = tmask(ji,jj,jk) * bdytmask(ji,jj)
               END DO
            END DO
         END DO
      ENDIF
         
      !  Ocean/land mask at u-, v-, and f-points   (computed from tmask)
      ! ----------------------------------------
      ! NB: at this point, fmask is designed for free slip lateral boundary condition
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
            END DO
            DO ji = 1, jpim1      ! NO vector opt.
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( umask  , 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( vmask  , 'V', 1._wp )
      CALL lbc_lnk( fmask  , 'F', 1._wp )

 
      ! Ocean/land mask at wu-, wv- and w points    (computed from tmask)
      !-----------------------------------------
      wmask (:,:,1) = tmask(:,:,1)     ! surface
      wumask(:,:,1) = umask(:,:,1)
      wvmask(:,:,1) = vmask(:,:,1)
      DO jk = 2, jpk                   ! interior values
         wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
         wumask(:,:,jk) = umask(:,:,jk) * umask(:,:,jk-1)   
         wvmask(:,:,jk) = vmask(:,:,jk) * vmask(:,:,jk-1)
      END DO


      ! Ocean/land column mask at t-, u-, and v-points   (i.e. at least 1 wet cell in the vertical)
      ! ----------------------------------------------
      ssmask (:,:) = MAXVAL( tmask(:,:,:), DIM=3 )
      ssumask(:,:) = MAXVAL( umask(:,:,:), DIM=3 )
      ssvmask(:,:) = MAXVAL( vmask(:,:,:), DIM=3 )


      ! Interior domain mask  (used for global sum)
      ! --------------------
      !
      iif = jpreci   ;   iil = nlci - jpreci + 1
      ijf = jprecj   ;   ijl = nlcj - jprecj + 1
      !
      !                          ! halo mask : 0 on the halo and 1 elsewhere
      tmask_h(:,:) = 1._wp                  
      tmask_h( 1 :iif,   :   ) = 0._wp      ! first columns
      tmask_h(iil:jpi,   :   ) = 0._wp      ! last  columns (including mpp extra columns)
      tmask_h(   :   , 1 :ijf) = 0._wp      ! first rows
      tmask_h(   :   ,ijl:jpj) = 0._wp      ! last  rows (including mpp extra rows)
      !
      !                          ! north fold mask
      tpol(1:jpiglo) = 1._wp 
      fpol(1:jpiglo) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0._wp
         fpol(     1    :jpiglo) = 0._wp
         IF( mjg(nlej) == jpjglo ) THEN                  ! only half of the nlcj-1 row for tmask_h
            DO ji = iif+1, iil-1
               tmask_h(ji,nlej-1) = tmask_h(ji,nlej-1) * tpol(mig(ji))
            END DO
         ENDIF
      ENDIF
      !
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0._wp
         fpol(jpiglo/2+1:jpiglo) = 0._wp
      ENDIF
      !
      !                          ! interior mask : 2D ocean mask x halo mask 
      tmask_i(:,:) = ssmask(:,:) * tmask_h(:,:)


      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------  
      IF( rn_shlat /= 0 ) THEN      ! Not free-slip lateral boundary condition
         !
         CALL wrk_alloc( jpi,jpj,   zwf )
         !
         DO jk = 1, jpk
            zwf(:,:) = fmask(:,:,jk)         
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fmask(ji,jj,jk) == 0._wp ) THEN
                     fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                        &                                           zwf(ji-1,jj), zwf(ji,jj-1)  )  )
                  ENDIF
               END DO
            END DO
            DO jj = 2, jpjm1
               IF( fmask(1,jj,jk) == 0._wp ) THEN
                  fmask(1  ,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
               ENDIF
               IF( fmask(jpi,jj,jk) == 0._wp ) THEN
                  fmask(jpi,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
               ENDIF
            END DO         
            DO ji = 2, jpim1
               IF( fmask(ji,1,jk) == 0._wp ) THEN
                  fmask(ji, 1 ,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
               ENDIF
               IF( fmask(ji,jpj,jk) == 0._wp ) THEN
                  fmask(ji,jpj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
               ENDIF
            END DO
         END DO
         !
         CALL wrk_dealloc( jpi,jpj,   zwf )
         !
         CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
         !
         ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 ) depending on ln_vorlat
         !
      ENDIF
      
      ! User defined alteration of fmask (use to reduce ocean transport in specified straits)
      ! -------------------------------- 
      !
      CALL usr_def_fmask( cn_cfg, nn_cfg, fmask )
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_msk')
      !
   END SUBROUTINE dom_msk
   
   !!======================================================================
END MODULE dommsk
