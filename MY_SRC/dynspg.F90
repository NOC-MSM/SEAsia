MODULE dynspg
   !!======================================================================
   !!                       ***  MODULE  dynspg  ***
   !! Ocean dynamics:  surface pressure gradient control
   !!======================================================================
   !! History :  1.0  ! 2005-12  (C. Talandier, G. Madec, V. Garnier)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg     : update the dynamics trend with surface pressure gradient 
   !!   dyn_spg_init: initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE c1d            ! 1D vertical configuration
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition: ocean
   USE sbcapr         ! surface boundary condition: atmospheric pressure
   USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
   USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
   USE sbctide        ! 
   USE updtide        ! 
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE prtctl         ! Print control                     (prt_ctl routine)
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg        ! routine called by step module
   PUBLIC   dyn_spg_init   ! routine called by opa module

   INTEGER ::   nspg = 0   ! type of surface pressure gradient scheme defined from lk_dynspg_... 
!jth
   LOGICAL, PUBLIC ::  ln_ulimit
   REAL(wp), PUBLIC :: cn_ulimit,cnn_ulimit
!
   !                       ! Parameter to control the surface pressure gradient scheme
   INTEGER, PARAMETER ::   np_TS  = 1   ! split-explicit time stepping (Time-Splitting)
   INTEGER, PARAMETER ::   np_EXP = 0   !       explicit time stepping
   INTEGER, PARAMETER ::   np_NO  =-1   ! no surface pressure gradient, no scheme

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: dynspg.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_spg( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg  ***
      !!
      !! ** Purpose :   compute surface pressure gradient including the 
      !!              atmospheric pressure forcing (ln_apr_dyn=T).
      !!
      !! ** Method  :   Two schemes:
      !!              - explicit       : the spg is evaluated at now
      !!              - split-explicit : a time splitting technique is used
      !!
      !!              ln_apr_dyn=T : the atmospheric pressure forcing is applied 
      !!             as the gradient of the inverse barometer ssh:
      !!                apgu = - 1/rau0 di[apr] = 0.5*grav di[ssh_ib+ssh_ibb]
      !!                apgv = - 1/rau0 dj[apr] = 0.5*grav dj[ssh_ib+ssh_ibb]
      !!             Note that as all external forcing a time averaging over a two rdt
      !!             period is used to prevent the divergence of odd and even time step.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk                             ! dummy loop indices
      REAL(wp) ::   z2dt, zg_2, zintp, zgrau0r             ! temporary scalar
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zpice
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg')
      !
      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi,jpj,jpk,   ztrdu, ztrdv ) 
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF
      !
      IF(      ln_apr_dyn                                                &   ! atmos. pressure
         .OR.  ( .NOT.ln_dynspg_ts .AND. (ln_tide_pot .AND. ln_tide) )   &   ! tide potential (no time slitting)
         .OR.  nn_ice_embd == 2  ) THEN                                      ! embedded sea-ice
         !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               spgu(ji,jj) = 0._wp
               spgv(ji,jj) = 0._wp
            END DO
         END DO         
         !
         IF( ln_apr_dyn .AND. .NOT.ln_dynspg_ts ) THEN   !==  Atmospheric pressure gradient (added later in time-split case) ==!
            zg_2 = grav * 0.5
            DO jj = 2, jpjm1                          ! gradient of Patm using inverse barometer ssh
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + zg_2 * (  ssh_ib (ji+1,jj) - ssh_ib (ji,jj)    &
                     &                      + ssh_ibb(ji+1,jj) - ssh_ibb(ji,jj)  ) * r1_e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + zg_2 * (  ssh_ib (ji,jj+1) - ssh_ib (ji,jj)    &
                     &                      + ssh_ibb(ji,jj+1) - ssh_ibb(ji,jj)  ) * r1_e2v(ji,jj)
               END DO
            END DO
         ENDIF
         !
         !                                    !==  tide potential forcing term  ==!
         IF( .NOT.ln_dynspg_ts .AND. ( ln_tide_pot .AND. ln_tide )  ) THEN   ! N.B. added directly at sub-time-step in ts-case
            !
            CALL upd_tide( kt )                      ! update tide potential
            !
            DO jj = 2, jpjm1                         ! add tide potential forcing
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
               END DO 
            END DO
         ENDIF
         !
         IF( nn_ice_embd == 2 ) THEN          !== embedded sea ice: Pressure gradient due to snow-ice mass ==!
            CALL wrk_alloc( jpi,jpj,   zpice )
            !                                            
            zintp = REAL( MOD( kt-1, nn_fsbc ) ) / REAL( nn_fsbc )
            zgrau0r     = - grav * r1_rau0
            zpice(:,:) = (  zintp * snwice_mass(:,:) + ( 1.- zintp ) * snwice_mass_b(:,:)  ) * zgrau0r
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu(ji,jj) = spgu(ji,jj) + ( zpice(ji+1,jj) - zpice(ji,jj) ) * r1_e1u(ji,jj)
                  spgv(ji,jj) = spgv(ji,jj) + ( zpice(ji,jj+1) - zpice(ji,jj) ) * r1_e2v(ji,jj)
               END DO
            END DO
            !
            CALL wrk_dealloc( jpi,jpj,   zpice )         
         ENDIF
         !
         DO jk = 1, jpkm1                    !== Add all terms to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
                  va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
               END DO
            END DO
         END DO    
         !
!!gm add here a call to dyn_trd for ice pressure gradient, the surf pressure trends ????
         !    
      ENDIF
      !
      SELECT CASE ( nspg )                   !== surface pressure gradient computed and add to the general trend ==!
      CASE ( np_EXP )   ;   CALL dyn_spg_exp( kt )              ! explicit
      CASE ( np_TS  )   ;   CALL dyn_spg_ts ( kt )              ! time-splitting
      END SELECT
      !                    
      IF( l_trddyn )   THEN                  ! save the surface pressure gradient trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_spg, kt )
         CALL wrk_dealloc( jpi,jpj,jpk,   ztrdu, ztrdv ) 
      ENDIF
      !                                      ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg  - Ua: ', mask1=umask, &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg')
      !
   END SUBROUTINE dyn_spg


   SUBROUTINE dyn_spg_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_init  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options for 
      !!              surface pressure gradient schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! local integers
      !
      NAMELIST/namdyn_spg/ ln_dynspg_exp       , ln_dynspg_ts,   &
      &                    ln_bt_fw, ln_bt_av  , ln_bt_auto  ,   &
      &                    nn_baro , rn_bt_cmax, nn_bt_flt,ln_ulimit,cn_ulimit,cnn_ulimit
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_init')
      !
      REWIND( numnam_ref )              ! Namelist namdyn_spg in reference namelist : Free surface
      READ  ( numnam_ref, namdyn_spg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_spg in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist namdyn_spg in configuration namelist : Free surface
      READ  ( numnam_cfg, namdyn_spg, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_spg in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdyn_spg )
      !
      IF(lwp) THEN             ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_spg_init : choice of the surface pressure gradient scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '     Explicit free surface                  ln_dynspg_exp = ', ln_dynspg_exp
         WRITE(numout,*) '     Free surface with time splitting       ln_dynspg_ts  = ', ln_dynspg_ts

         write(numout,*) '     Limit velocities                       ln_ulimit     = ',ln_ulimit
         write(numout,*) '     Limit velocities      max inverse Courant number     = ',cn_ulimit
         write(numout,*) '     Limit velocities   multiplier for divergant flow     = ',cnn_ulimit

      ENDIF
      !                          ! Control of surface pressure gradient scheme options
                                     nspg =  np_NO    ;   ioptio = 0
      IF( ln_dynspg_exp ) THEN   ;   nspg =  np_EXP   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynspg_ts  ) THEN   ;   nspg =  np_TS    ;   ioptio = ioptio + 1   ;   ENDIF
      !
      IF( ioptio  > 1 )   CALL ctl_stop( 'Choose only one surface pressure gradient scheme' )
      IF( ioptio == 0 )   CALL ctl_warn( 'NO surface pressure gradient trend in momentum Eqs.' )
      IF( ln_dynspg_exp .AND. ln_isfcav )   &
           &   CALL ctl_stop( ' dynspg_exp not tested with ice shelf cavity ' )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         IF( nspg == np_EXP )   WRITE(numout,*) '      ===>>   explicit free surface'
         IF( nspg == np_TS  )   WRITE(numout,*) '      ===>>   free surface with time splitting scheme'
         IF( nspg == np_NO  )   WRITE(numout,*) '      ===>>   No surface surface pressure gradient trend in momentum Eqs.'
      ENDIF
      !
      IF( nspg == np_TS ) THEN   ! split-explicit scheme initialisation
         CALL dyn_spg_ts_init          ! do it first: set nn_baro used to allocate some arrays later on
         IF( dyn_spg_ts_alloc() /= 0  )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays' )
         IF( neuler/=0 .AND. ln_bt_fw )   CALL ts_rst( nit000, 'READ' )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_init')
      !
   END SUBROUTINE dyn_spg_init

  !!======================================================================
END MODULE dynspg
