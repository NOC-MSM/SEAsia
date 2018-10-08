MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec) Merge of dtatem and dtasal & suppression of tb,tn/sb,sn 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_tem    : analytical profile for initial Temperature
   !!   istate_sal    : analytical profile for initial Salinity
   !!   istate_eel    : initial state setting of EEL R5 configuration
   !!   istate_gyre   : initial state setting of GYRE configuration
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE c1d             ! 1D vertical configuration
   USE daymod          ! calendar
   USE eosbn2          ! eq. of state, Brunt Vaisala frequency (eos     routine)
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE dtatsd          ! data temperature and salinity   (dta_tsd routine)
   USE dtauvd          ! data: U & V current             (dta_uvd routine)
   USE zpshde          ! partial step: hor. derivative (zps_hde routine)
   USE eosbn2          ! equation of state            (eos bn2 routine)
   USE domvvl          ! varying vertical mesh
   USE dynspg_oce      ! pressure gradient schemes
   USE dynspg_flt      ! filtered free surface
   USE sol_oce         ! ocean solver variables
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE restart         ! restart
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zuvd    ! U & V data workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('istate_init')
      !

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) 'istate_ini : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      
      CALL day_init                           ! model calendar (using both namelist and restart infos)

      CALL dta_tsd_init                       ! Initialisation of T & S input data
      IF( lk_c1d ) CALL dta_uvd_init          ! Initialization of U & V input data

      rhd  (:,:,:  ) = 0._wp   ;   rhop (:,:,:  ) = 0._wp      ! set one for all to 0 at level jpk
      rn2b (:,:,:  ) = 0._wp   ;   rn2  (:,:,:  ) = 0._wp      ! set one for all to 0 at levels 1 and jpk
      tsa  (:,:,:,:) = 0._wp                                   ! set one for all to 0 at level jpk
      rab_b(:,:,:,:) = 0._wp   ;   rab_n(:,:,:,:) = 0._wp      ! set one for all to 0 at level jpk

      IF( ln_rstart ) THEN                    ! Restart from a file
         !                                    ! -------------------
         CALL rst_read                           ! Read the restart file
      ELSE
         !                                    ! Start from rest
         !                                    ! ---------------
         numror = 0                              ! define numror = 0 -> no restart file to read
         neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
         !                                       ! Initialization of ocean to zero
         !   before fields      !       now fields     
         sshb (:,:)   = 0._wp   ;   sshn (:,:)   = 0._wp
         ub   (:,:,:) = 0._wp   ;   un   (:,:,:) = 0._wp
         vb   (:,:,:) = 0._wp   ;   vn   (:,:,:) = 0._wp  
         rotb (:,:,:) = 0._wp   ;   rotn (:,:,:) = 0._wp
         hdivb(:,:,:) = 0._wp   ;   hdivn(:,:,:) = 0._wp
         !
         IF( cp_cfg == 'eel' ) THEN
            CALL istate_eel                      ! EEL   configuration : start from pre-defined U,V T-S fields
         ELSEIF( cp_cfg == 'gyre' ) THEN         
            CALL istate_gyre                     ! GYRE  configuration : start from pre-defined T-S fields
         ELSE                                    ! Initial T-S, U-V fields read in files
            IF ( ln_tsd_init ) THEN              ! read 3D T and S data at nit000
               CALL dta_tsd( nit000, tsb )  
               tsn(:,:,:,:) = tsb(:,:,:,:)
               !
            ELSE                                 ! Initial T-S fields defined analytically
               CALL istate_t_s
            ENDIF
            IF ( ln_uvd_init .AND. lk_c1d ) THEN ! read 3D U and V data at nit000
               CALL wrk_alloc( jpi, jpj, jpk, 2, zuvd )
               CALL dta_uvd( nit000, zuvd )
               ub(:,:,:) = zuvd(:,:,:,1) ;  un(:,:,:) = ub(:,:,:)
               vb(:,:,:) = zuvd(:,:,:,2) ;  vn(:,:,:) = vb(:,:,:)
               CALL wrk_dealloc( jpi, jpj, jpk, 2, zuvd )
            ENDIF
         ENDIF
         !
         CALL eos( tsb, rhd, rhop, gdept_0(:,:,:) )        ! before potential and in situ densities
#if ! defined key_c1d
         IF( ln_zps .AND. .NOT. ln_isfcav)                                 &
            &            CALL zps_hde    ( nit000, jpts, tsb, gtsu, gtsv,  &    ! Partial steps: before horizontal gradient
            &                                            rhd, gru , grv    )  ! of t, s, rd at the last ocean level
         IF( ln_zps .AND.       ln_isfcav)                                 &
            &            CALL zps_hde_isf( nit000, jpts, tsb, gtsu, gtsv,  &    ! Partial steps for top cell (ISF)
            &                                            rhd, gru , grv , aru , arv , gzu , gzv , ge3ru , ge3rv ,   &
            &                                     gtui, gtvi, grui, grvi, arui, arvi, gzui, gzvi, ge3rui, ge3rvi    ) ! of t, s, rd at the last ocean level
#endif
         !   
         ! - ML - sshn could be modified by istate_eel, so that initialization of fse3t_b is done here
         IF( lk_vvl ) THEN
            DO jk = 1, jpk
               fse3t_b(:,:,jk) = fse3t_n(:,:,jk)
            ENDDO
         ENDIF
         ! 
      ENDIF
      !
      IF( lk_agrif ) THEN                  ! read free surface arrays in restart file
         IF( ln_rstart ) THEN
            IF( lk_dynspg_flt )  THEN      ! read or initialize the following fields
               !                           ! gcx, gcxb for agrif_opa_init
               IF( sol_oce_alloc()  > 0 )   CALL ctl_stop('agrif sol_oce_alloc: allocation of arrays failed')
               CALL flt_rst( nit000, 'READ' )
            ENDIF
         ENDIF                             ! explicit case not coded yet with AGRIF
      ENDIF
      !
      ! 
      ! Initialize "now" and "before" barotropic velocities:
      ! Do it whatever the free surface method, these arrays
      ! being eventually used
      !
      !
      un_b(:,:) = 0._wp ; vn_b(:,:) = 0._wp
      ub_b(:,:) = 0._wp ; vb_b(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               un_b(ji,jj) = un_b(ji,jj) + fse3u_n(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
               vn_b(ji,jj) = vn_b(ji,jj) + fse3v_n(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
               !
               ub_b(ji,jj) = ub_b(ji,jj) + fse3u_b(ji,jj,jk) * ub(ji,jj,jk) * umask(ji,jj,jk)
               vb_b(ji,jj) = vb_b(ji,jj) + fse3v_b(ji,jj,jk) * vb(ji,jj,jk) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      un_b(:,:) = un_b(:,:) * hur  (:,:)
      vn_b(:,:) = vn_b(:,:) * hvr  (:,:)
      !
      ub_b(:,:) = ub_b(:,:) * hur_b(:,:)
      vb_b(:,:) = vb_b(:,:) * hvr_b(:,:)
      !
      !
      IF( nn_timing == 1 )   CALL timing_stop('istate_init')
      !
   END SUBROUTINE istate_init


   SUBROUTINE istate_t_s
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE istate_t_s  ***
      !!   
      !! ** Purpose :   Intialization of the temperature field with an 
      !!      analytical profile or a file (i.e. in EEL configuration)
      !!
      !! ** Method  : - temperature: use Philander analytic profile
      !!              - salinity   : use to a constant value 35.5
      !!
      !! References :  Philander ???
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) ::   zsal = 35.50
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_t_s : Philander s initial temperature profile'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~   and constant salinity (',zsal,' psu)'
      !
      DO jk = 1, jpk
         tsn(:,:,jk,jp_tem) = (  ( ( 7.5 - 0. * ABS( gphit(:,:) )/30. ) * ( 1.-TANH((fsdept(:,:,jk)-80.)/30.) )   &
            &                + 10. * ( 5000. - fsdept(:,:,jk) ) /5000.)  ) * tmask(:,:,jk)
         tsb(:,:,jk,jp_tem) = tsn(:,:,jk,jp_tem)
      END DO
      tsn(:,:,:,jp_sal) = zsal * tmask(:,:,:)
      tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
      !
   END SUBROUTINE istate_t_s


   SUBROUTINE istate_eel
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_eel  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers for EEL R5
      !!      configuration (channel with or without a topographic bump)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity field
      !!              - set velocity field including horizontal divergence
      !!                and relative vorticity fields
      !!----------------------------------------------------------------------
      USE divcur     ! hor. divergence & rel. vorticity      (div_cur routine)
      USE iom
      !
      INTEGER  ::   inum              ! temporary logical unit
      INTEGER  ::   ji, jj, jk        ! dummy loop indices
      INTEGER  ::   ijloc
      REAL(wp) ::   zh1, zh2, zslope, zcst, zfcor   ! temporary scalars
      REAL(wp) ::   zt1  = 15._wp                   ! surface temperature value (EEL R5)
      REAL(wp) ::   zt2  =  5._wp                   ! bottom  temperature value (EEL R5)
      REAL(wp) ::   zsal = 35.0_wp                  ! constant salinity (EEL R2, R5 and R6)
      REAL(wp) ::   zueel = 0.1_wp                  ! constant uniform zonal velocity (EEL R5)
      REAL(wp), DIMENSION(jpiglo,jpjglo) ::   zssh  ! initial ssh over the global domain
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( jp_cfg ) 
         !                                              ! ====================
         CASE ( 5 )                                     ! EEL R5 configuration
            !                                           ! ====================
            !
            ! set temperature field with a linear profile
            ! -------------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: linear temperature profile'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            zh1 = gdept_1d(  1  )
            zh2 = gdept_1d(jpkm1)
            !
            zslope = ( zt1 - zt2 ) / ( zh1 - zh2 )
            zcst   = ( zt1 * ( zh1 - zh2) - ( zt1 - zt2 ) * zh1 ) / ( zh1 - zh2 )
            !
            DO jk = 1, jpk
               tsn(:,:,jk,jp_tem) = ( zt2 + zt1 * exp( - fsdept(:,:,jk) / 1000 ) ) * tmask(:,:,jk)
               tsb(:,:,jk,jp_tem) = tsn(:,:,jk,jp_tem)
            END DO
            !
            IF(lwp) CALL prizre( tsn(:,:,:,jp_tem), jpi   , jpj   , jpk   , jpj/2 ,   &
               &                             1     , jpi   , 5     , 1     , jpk   ,   &
               &                             1     , 1.    , numout                  )
            !
            ! set salinity field to a constant value
            ! --------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: constant salinity field, S = ', zsal
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            tsn(:,:,:,jp_sal) = zsal * tmask(:,:,:)
            tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
            !
            ! set the dynamics: U,V, hdiv, rot (and ssh if necessary)
            ! ----------------
            ! Start EEL5 configuration with barotropic geostrophic velocities 
            ! according the sshb and sshn SSH imposed.
            ! we assume a uniform grid (hence the use of e1t(1,1) for delta_y)
            ! we use the Coriolis frequency at mid-channel.   
            ub(:,:,:) = zueel * umask(:,:,:)
            un(:,:,:) = ub(:,:,:)
            ijloc = mj0(INT(jpjglo-1)/2)
            zfcor = ff(1,ijloc)
            !
            DO jj = 1, jpjglo
               zssh(:,jj) = - (FLOAT(jj)- FLOAT(jpjglo-1)/2.)*zueel*e1t(1,1)*zfcor/grav 
            END DO
            !
            IF(lwp) THEN
               WRITE(numout,*) ' Uniform zonal velocity for EEL R5:',zueel
               WRITE(numout,*) ' Geostrophic SSH profile as a function of y:'
               WRITE(numout,'(12(1x,f6.2))') zssh(1,:)
            ENDIF
            !
            DO jj = 1, nlcj
               DO ji = 1, nlci
                  sshb(ji,jj) = zssh( mig(ji) , mjg(jj) ) * tmask(ji,jj,1)
               END DO
            END DO
            sshb(nlci+1:jpi,      :   ) = 0.e0      ! set to zero extra mpp columns
            sshb(      :   ,nlcj+1:jpj) = 0.e0      ! set to zero extra mpp rows
            !
            sshn(:,:) = sshb(:,:)                   ! set now ssh to the before value
            !
            IF( nn_rstssh /= 0 ) THEN  
               nn_rstssh = 0                        ! hand-made initilization of ssh 
               CALL ctl_warn( 'istate_eel: force nn_rstssh = 0' )
            ENDIF
            !
            CALL div_cur( nit000 )                  ! horizontal divergence and relative vorticity (curl)
            ! N.B. the vertical velocity will be computed from the horizontal divergence field
            ! in istate by a call to wzv routine


            !                                     ! ==========================
         CASE ( 2 , 6 )                           ! EEL R2 or R6 configuration
            !                                     ! ==========================
            !
            ! set temperature field with a NetCDF file
            ! ----------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R2 or R6: read initial temperature in a NetCDF file'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            CALL iom_open ( 'eel.initemp', inum )
            CALL iom_get ( inum, jpdom_data, 'initemp', tsb(:,:,:,jp_tem) ) ! read before temprature (tb)
            CALL iom_close( inum )
            !
            tsn(:,:,:,jp_tem) = tsb(:,:,:,jp_tem)                            ! set nox temperature to tb
            !
            IF(lwp) CALL prizre( tsn(:,:,:,jp_tem), jpi   , jpj   , jpk   , jpj/2 ,   &
               &                            1     , jpi   , 5     , 1     , jpk   ,   &
               &                            1     , 1.    , numout                  )
            !
            ! set salinity field to a constant value
            ! --------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: constant salinity field, S = ', zsal
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            tsn(:,:,:,jp_sal) = zsal * tmask(:,:,:)
            tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
            !
            !                                    ! ===========================
         CASE DEFAULT                            ! NONE existing configuration
            !                                    ! ===========================
            WRITE(ctmp1,*) 'EEL with a ', jp_cfg,' km resolution is not coded'
            CALL ctl_stop( ctmp1 )
            !
      END SELECT
      !
   END SUBROUTINE istate_eel


   SUBROUTINE istate_gyre
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_gyre  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers for GYRE
      !!      configuration (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity field
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk  ! dummy loop indices
      INTEGER            ::   inum          ! temporary logical unit
      INTEGER, PARAMETER ::   ntsinit = 0   ! (0/1) (analytical/input data files) T&S initialization
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( ntsinit)
      !
      CASE ( 0 )                  ! analytical T/S profil deduced from LEVITUS
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'istate_gyre : initial analytical T and S profil deduced from LEVITUS '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         !
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  tsn(ji,jj,jk,jp_tem) = (  16. - 12. * TANH( (fsdept(ji,jj,jk) - 400) / 700 )         )   &
                       &           * (-TANH( (500-fsdept(ji,jj,jk)) / 150 ) + 1) / 2               &
                       &       + (      15. * ( 1. - TANH( (fsdept(ji,jj,jk)-50.) / 1500.) )       &
                       &                - 1.4 * TANH((fsdept(ji,jj,jk)-100.) / 100.)               &    
                       &                + 7.  * (1500. - fsdept(ji,jj,jk)) / 1500.             )   & 
                       &           * (-TANH( (fsdept(ji,jj,jk) - 500) / 150) + 1) / 2
                  tsn(ji,jj,jk,jp_tem) = tsn(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)
                  tsb(ji,jj,jk,jp_tem) = tsn(ji,jj,jk,jp_tem)

                  tsn(ji,jj,jk,jp_sal) =  (  36.25 - 1.13 * TANH( (fsdept(ji,jj,jk) - 305) / 460 )  )  &
                     &              * (-TANH((500 - fsdept(ji,jj,jk)) / 150) + 1) / 2          &
                     &          + (  35.55 + 1.25 * (5000. - fsdept(ji,jj,jk)) / 5000.         &
                     &                - 1.62 * TANH( (fsdept(ji,jj,jk) - 60.  ) / 650. )       &
                     &                + 0.2  * TANH( (fsdept(ji,jj,jk) - 35.  ) / 100. )       &
                     &                + 0.2  * TANH( (fsdept(ji,jj,jk) - 1000.) / 5000.)    )  &
                     &              * (-TANH((fsdept(ji,jj,jk) - 500) / 150) + 1) / 2 
                  tsn(ji,jj,jk,jp_sal) = tsn(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
                  tsb(ji,jj,jk,jp_sal) = tsn(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                  ! T/S data fields read in dta_tem.nc/data_sal.nc files
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'istate_gyre : initial T and S read from dta_tem.nc/data_sal.nc files'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*) '              NetCDF FORMAT'

         ! Read temperature field
         ! ----------------------
         CALL iom_open ( 'data_tem', inum )
         CALL iom_get ( inum, jpdom_data, 'votemper', tsn(:,:,:,jp_tem) ) 
         CALL iom_close( inum )

         tsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem) * tmask(:,:,:) 
         tsb(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)

         ! Read salinity field
         ! -------------------
         CALL iom_open ( 'data_sal', inum )
         CALL iom_get ( inum, jpdom_data, 'vosaline', tsn(:,:,:,jp_sal) ) 
         CALL iom_close( inum )

         tsn(:,:,:,jp_sal) = tsn(:,:,:,jp_sal) * tmask(:,:,:) 
         tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
         !
      END SELECT
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '              Initial temperature and salinity profiles:'
         WRITE(numout, "(9x,' level   gdept_1d   temperature   salinity   ')" )
         WRITE(numout, "(10x, i4, 3f10.2)" ) ( jk, gdept_1d(jk), tsn(2,2,jk,jp_tem), tsn(2,2,jk,jp_sal), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE istate_gyre


   SUBROUTINE istate_uvg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE istate_uvg  ***
      !!
      !! ** Purpose :   Compute the geostrophic velocities from (tn,sn) fields
      !!
      !! ** Method  :   Using the hydrostatic hypothesis the now hydrostatic 
      !!      pressure is computed by integrating the in-situ density from the
      !!      surface to the bottom.
      !!                 p=integral [ rau*g dz ]
      !!----------------------------------------------------------------------
      USE dynspg          ! surface pressure gradient             (dyn_spg routine)
      USE divcur          ! hor. divergence & rel. vorticity      (div_cur routine)
      USE lbclnk          ! ocean lateral boundary condition (or mpp link)
      !
      INTEGER ::   ji, jj, jk        ! dummy loop indices
      INTEGER ::   indic             ! ???
      REAL(wp) ::   zmsv, zphv, zmsu, zphu, zalfg     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprn
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi, jpj, jpk, zprn)
      !
      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'istate_uvg : Start from Geostrophy'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      ! Compute the now hydrostatic pressure
      ! ------------------------------------

      zalfg = 0.5 * grav * rau0
      
      zprn(:,:,1) = zalfg * fse3w(:,:,1) * ( 1 + rhd(:,:,1) )       ! Surface value

      DO jk = 2, jpkm1                                              ! Vertical integration from the surface
         zprn(:,:,jk) = zprn(:,:,jk-1)   &
            &         + zalfg * fse3w(:,:,jk) * ( 2. + rhd(:,:,jk) + rhd(:,:,jk-1) )
      END DO  

      ! Compute geostrophic balance
      ! ---------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vertor opt.
               zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
                               + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1.  )
               zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                    + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
               zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)

               zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
                               + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1.  )
               zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
                    + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
               zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)

               ! Compute the geostrophic velocities
               un(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
               vn(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
            END DO
         END DO
      END DO

      IF(lwp) WRITE(numout,*) '         we force to zero bottom velocity'

      ! Susbtract the bottom velocity (level jpk-1 for flat bottom case)
      ! to have a zero bottom velocity

      DO jk = 1, jpkm1
         un(:,:,jk) = ( un(:,:,jk) - un(:,:,jpkm1) ) * umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) - vn(:,:,jpkm1) ) * vmask(:,:,jk)
      END DO

      CALL lbc_lnk( un, 'U', -1. )
      CALL lbc_lnk( vn, 'V', -1. )
      
      ub(:,:,:) = un(:,:,:)
      vb(:,:,:) = vn(:,:,:)
      
      ! WARNING !!!!!
      ! after initializing u and v, we need to calculate the initial streamfunction bsf.
      ! Otherwise, only the trend will be computed and the model will blow up (inconsistency).
      ! to do that, we call dyn_spg with a special trick:
      ! we fill ua and va with the velocities divided by dt, and the streamfunction will be brought to the
      ! right value assuming the velocities have been set up in one time step.
      ! we then set bsfd to zero (first guess for next step is d(psi)/dt = 0.)
      !  sets up s false trend to calculate the barotropic streamfunction.

      ua(:,:,:) = ub(:,:,:) / rdt
      va(:,:,:) = vb(:,:,:) / rdt

      ! calls dyn_spg. we assume euler time step, starting from rest.
      indic = 0
      CALL dyn_spg( nit000, indic )       ! surface pressure gradient

      ! the new velocity is ua*rdt

      CALL lbc_lnk( ua, 'U', -1. )
      CALL lbc_lnk( va, 'V', -1. )

      ub(:,:,:) = ua(:,:,:) * rdt
      vb(:,:,:) = va(:,:,:) * rdt
      ua(:,:,:) = 0.e0
      va(:,:,:) = 0.e0
      un(:,:,:) = ub(:,:,:)
      vn(:,:,:) = vb(:,:,:)
       
      ! Compute the divergence and curl

      CALL div_cur( nit000 )            ! now horizontal divergence and curl

      hdivb(:,:,:) = hdivn(:,:,:)       ! set the before to the now value
      rotb (:,:,:) = rotn (:,:,:)       ! set the before to the now value
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zprn)
      !
   END SUBROUTINE istate_uvg

   !!=====================================================================
END MODULE istate
