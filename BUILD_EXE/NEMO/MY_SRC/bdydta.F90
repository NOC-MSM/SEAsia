MODULE bdydta
   !!======================================================================
   !!                       ***  MODULE bdydta  ***
   !! Open boundary data : read the data for the unstructured open boundaries.
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2012-01  (C. Rousset) add ice boundary conditions for sea ice
   !!            4.0  !  2018     (C. Rousset) SI3 compatibility
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    bdy_dta      : read external data along open boundaries from file
   !!    bdy_dta_init : initialise arrays etc for reading of external data
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbcapr         ! atmospheric pressure forcing
   USE sbctide        ! Tidal forcing or not
   USE bdy_oce        ! ocean open boundary conditions  
   USE bdytides       ! tidal forcing at boundaries
#if defined key_si3
   USE ice            ! sea-ice variables
   USE icevar         ! redistribute ice input into categories
#endif
   !
   USE lib_mpp, ONLY: ctl_stop, ctl_nam
   USE fldread        ! read input fields
   USE iom            ! IOM library
   USE in_out_manager ! I/O logical units
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta          ! routine called by step.F90 and dynspg_ts.F90
   PUBLIC   bdy_dta_init     ! routine called by nemogcm.F90

   INTEGER , PARAMETER ::   jpbdyfld  = 17    ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_bdyssh = 1     ! 
   INTEGER , PARAMETER ::   jp_bdyu2d = 2     ! 
   INTEGER , PARAMETER ::   jp_bdyv2d = 3     !
   INTEGER , PARAMETER ::   jp_bdyu3d = 4     !
   INTEGER , PARAMETER ::   jp_bdyv3d = 5     !
   INTEGER , PARAMETER ::   jp_bdytem = 6     ! 
   INTEGER , PARAMETER ::   jp_bdysal = 7     ! 
   INTEGER , PARAMETER ::   jp_bdya_i = 8     ! 
   INTEGER , PARAMETER ::   jp_bdyh_i = 9     ! 
   INTEGER , PARAMETER ::   jp_bdyh_s = 10    ! 
   INTEGER , PARAMETER ::   jp_bdyt_i = 11    ! 
   INTEGER , PARAMETER ::   jp_bdyt_s = 12    ! 
   INTEGER , PARAMETER ::   jp_bdytsu = 13    ! 
   INTEGER , PARAMETER ::   jp_bdys_i = 14    ! 
   INTEGER , PARAMETER ::   jp_bdyaip = 15    ! 
   INTEGER , PARAMETER ::   jp_bdyhip = 16    ! 
   INTEGER , PARAMETER ::   jp_bdyhil = 17    ! 
#if ! defined key_si3
   INTEGER , PARAMETER ::   jpl = 1
#endif

!$AGRIF_DO_NOT_TREAT
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:,:), TARGET ::   bf   ! structure of input fields (file informations, fields read)
!$AGRIF_END_DO_NOT_TREAT

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dta( kt, kit, kt_offset )
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta  ***
      !!                    
      !! ** Purpose :   Update external data for open boundary conditions
      !!
      !! ** Method  :   Use fldread.F90
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           ::   kt           ! ocean time-step index 
      INTEGER, INTENT(in), OPTIONAL ::   kit          ! subcycle time-step index (for timesplitting option)
      INTEGER, INTENT(in), OPTIONAL ::   kt_offset    ! time offset in units of timesteps. NB. if kit
      !                                               ! is present then units = subcycle timesteps.
      !                                               ! kt_offset = 0 => get data at "now" time level
      !                                               ! kt_offset = -1 => get data at "before" time level
      !                                               ! kt_offset = +1 => get data at "after" time level
      !                                               ! etc.
      !
      INTEGER ::  jbdy, jfld, jstart, jend, ib, jl    ! dummy loop indices
      INTEGER ::  ii, ij, ik, igrd, ipl               ! local integers
      TYPE(OBC_DATA)         , POINTER ::   dta_alias        ! short cut
      TYPE(FLD), DIMENSION(:), POINTER ::   bf_alias
      !!---------------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_dta')
      !
      ! Initialise data arrays once for all from initial conditions where required
      !---------------------------------------------------------------------------
      IF( kt == nit000 .AND. .NOT.PRESENT(kit) ) THEN

         ! Calculate depth-mean currents
         !-----------------------------

         DO jbdy = 1, nb_bdy
            !
            IF( nn_dyn2d_dta(jbdy) == 0 ) THEN 
               IF( dta_bdy(jbdy)%lneed_ssh ) THEN 
                  igrd = 1
                  DO ib = 1, idx_bdy(jbdy)%nblenrim(igrd)   ! ssh is allocated and used only on the rim
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%ssh(ib) = sshn(ii,ij) * tmask(ii,ij,1)         
                  END DO
               ENDIF
               IF( ASSOCIATED(dta_bdy(jbdy)%u2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
                  igrd = 2
                  DO ib = 1, SIZE(dta_bdy(jbdy)%u2d)      ! u2d is used either over the whole bdy or only on the rim
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%u2d(ib) = un_b(ii,ij) * umask(ii,ij,1)         
                  END DO
               ENDIF
               IF( ASSOCIATED(dta_bdy(jbdy)%v2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
                  igrd = 3
                  DO ib = 1, SIZE(dta_bdy(jbdy)%v2d)      ! v2d is used either over the whole bdy or only on the rim
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%v2d(ib) = vn_b(ii,ij) * vmask(ii,ij,1)         
                  END DO
               ENDIF
            ENDIF
            !
            IF( nn_dyn3d_dta(jbdy) == 0 ) THEN 
               IF( dta_bdy(jbdy)%lneed_dyn3d ) THEN 
                  igrd = 2 
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%u3d(ib,ik) =  ( un(ii,ij,ik) - un_b(ii,ij) ) * umask(ii,ij,ik)         
                     END DO
                  END DO
                  igrd = 3 
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%v3d(ib,ik) =  ( vn(ii,ij,ik) - vn_b(ii,ij) ) * vmask(ii,ij,ik)         
                     END DO
                  END DO
               ENDIF
            ENDIF

            IF( nn_tra_dta(jbdy) == 0 ) THEN 
               IF( dta_bdy(jbdy)%lneed_tra ) THEN
                  igrd = 1 
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%tem(ib,ik) = tsn(ii,ij,ik,jp_tem) * tmask(ii,ij,ik)         
                        dta_bdy(jbdy)%sal(ib,ik) = tsn(ii,ij,ik,jp_sal) * tmask(ii,ij,ik)         
                     END DO
                  END DO
               ENDIF
            ENDIF

#if defined key_si3
            IF( nn_ice_dta(jbdy) == 0 ) THEN    ! set ice to initial values
               IF( dta_bdy(jbdy)%lneed_ice ) THEN
                  igrd = 1   
                  DO jl = 1, jpl
                     DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%a_i(ib,jl) =  a_i (ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%h_i(ib,jl) =  h_i (ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%h_s(ib,jl) =  h_s (ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%t_i(ib,jl) =  SUM(t_i (ii,ij,:,jl)) * r1_nlay_i * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%t_s(ib,jl) =  SUM(t_s (ii,ij,:,jl)) * r1_nlay_s * tmask(ii,ij,1)
                        dta_bdy(jbdy)%tsu(ib,jl) =  t_su(ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%s_i(ib,jl) =  s_i (ii,ij,jl) * tmask(ii,ij,1)
                        ! melt ponds
                        dta_bdy(jbdy)%aip(ib,jl) =  a_ip(ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%hip(ib,jl) =  h_ip(ii,ij,jl) * tmask(ii,ij,1) 
                        dta_bdy(jbdy)%hil(ib,jl) =  h_il(ii,ij,jl) * tmask(ii,ij,1) 
                     END DO
                  END DO
               ENDIF
            ENDIF
#endif
         END DO ! jbdy
         !
      ENDIF ! kt == nit000

      ! update external data from files
      !--------------------------------

      DO jbdy = 1, nb_bdy

         dta_alias => dta_bdy(jbdy)
         bf_alias  => bf(:,jbdy)

         ! read/update all bdy data
         ! ------------------------
         CALL fld_read( kt, 1, bf_alias, kit = kit, kt_offset = kt_offset )

         ! apply some corrections in some specific cases...
         ! --------------------------------------------------
         !
         ! if runoff condition: change river flow we read (in m3/s) into barotropic velocity (m/s)
         IF( cn_tra(jbdy) == 'runoff' ) THEN   ! runoff
            !
            IF( ASSOCIATED(dta_bdy(jbdy)%u2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
               igrd = 2                         ! zonal flow (m3/s) to barotropic zonal velocity (m/s)
               DO ib = 1, SIZE(dta_alias%u2d)   ! u2d is used either over the whole bdy or only on the rim
                  ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                  ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                  dta_alias%u2d(ib) = dta_alias%u2d(ib) / ( e2u(ii,ij) * hu_0(ii,ij) )
               END DO
            ENDIF
            IF( ASSOCIATED(dta_bdy(jbdy)%v2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
               igrd = 3                         ! meridional flow (m3/s) to barotropic meridional velocity (m/s)
               DO ib = 1, SIZE(dta_alias%v2d)   ! v2d is used either over the whole bdy or only on the rim
                  ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                  ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                  dta_alias%v2d(ib) = dta_alias%v2d(ib) / ( e1v(ii,ij) * hv_0(ii,ij) )
               END DO
            ENDIF
         ENDIF

         ! tidal harmonic forcing ONLY: initialise arrays
         IF( nn_dyn2d_dta(jbdy) == 2 ) THEN   ! we did not read ssh, u/v2d 
            IF( ASSOCIATED(dta_alias%ssh) ) dta_alias%ssh(:) = 0._wp
            IF( ASSOCIATED(dta_alias%u2d) ) dta_alias%u2d(:) = 0._wp
            IF( ASSOCIATED(dta_alias%v2d) ) dta_alias%v2d(:) = 0._wp
         ENDIF

         ! If full velocities in boundary data, then split it into barotropic and baroclinic component
         IF( bf_alias(jp_bdyu3d)%ltotvel ) THEN     ! if we read 3D total velocity (can be true only if u3d was read)
            !
            igrd = 2                       ! zonal velocity
            DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
               ii   = idx_bdy(jbdy)%nbi(ib,igrd)
               ij   = idx_bdy(jbdy)%nbj(ib,igrd)
               dta_alias%u2d(ib) = 0._wp   ! compute barotrope zonal velocity and put it in u2d
               DO ik = 1, jpkm1
                  dta_alias%u2d(ib) = dta_alias%u2d(ib) + e3u_n(ii,ij,ik) * umask(ii,ij,ik) * dta_alias%u3d(ib,ik)
               END DO
               dta_alias%u2d(ib) =  dta_alias%u2d(ib) * r1_hu_n(ii,ij)
               DO ik = 1, jpkm1            ! compute barocline zonal velocity and put it in u3d
                  dta_alias%u3d(ib,ik) = dta_alias%u3d(ib,ik) - dta_alias%u2d(ib)
               END DO
            END DO
            igrd = 3                       ! meridional velocity
            DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
               ii   = idx_bdy(jbdy)%nbi(ib,igrd)
               ij   = idx_bdy(jbdy)%nbj(ib,igrd)
               dta_alias%v2d(ib) = 0._wp   ! compute barotrope meridional velocity and put it in v2d
               DO ik = 1, jpkm1
                  dta_alias%v2d(ib) = dta_alias%v2d(ib) + e3v_n(ii,ij,ik) * vmask(ii,ij,ik) * dta_alias%v3d(ib,ik)
               END DO
               dta_alias%v2d(ib) =  dta_alias%v2d(ib) * r1_hv_n(ii,ij)
               DO ik = 1, jpkm1            ! compute barocline meridional velocity and put it in v3d
                  dta_alias%v3d(ib,ik) = dta_alias%v3d(ib,ik) - dta_alias%v2d(ib)
               END DO
            END DO
         ENDIF   ! ltotvel

         ! update tidal harmonic forcing
         IF( PRESENT(kit) .AND. nn_dyn2d_dta(jbdy) .GE. 2 ) THEN
            CALL bdytide_update( kt = kt, idx = idx_bdy(jbdy), dta = dta_alias, td = tides(jbdy),   & 
               &                 kit = kit, kt_offset = kt_offset )
         ENDIF

         !  atm surface pressure : add inverted barometer effect to ssh if it was read
         IF ( ln_apr_obc .AND. TRIM(bf_alias(jp_bdyssh)%clrootname) /= 'NOT USED' ) THEN
            igrd = 1
            DO ib = 1, idx_bdy(jbdy)%nblenrim(igrd)   ! ssh is used only on the rim
               ii   = idx_bdy(jbdy)%nbi(ib,igrd)
               ij   = idx_bdy(jbdy)%nbj(ib,igrd)
               dta_alias%ssh(ib) = dta_alias%ssh(ib) + ssh_ib(ii,ij)
            END DO
         ENDIF

#if defined key_si3
         IF( dta_alias%lneed_ice .AND. idx_bdy(jbdy)%nblen(1) > 0 ) THEN
            ! fill temperature and salinity arrays
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyt_i)%fnow(:,1,:) = rice_tem (jbdy)
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyt_s)%fnow(:,1,:) = rice_tem (jbdy)
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' )   bf_alias(jp_bdytsu)%fnow(:,1,:) = rice_tem (jbdy)
            IF( TRIM(bf_alias(jp_bdys_i)%clrootname) == 'NOT USED' )   bf_alias(jp_bdys_i)%fnow(:,1,:) = rice_sal (jbdy)
            IF( TRIM(bf_alias(jp_bdyaip)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyaip)%fnow(:,1,:) = rice_apnd(jbdy) * & ! rice_apnd is the pond fraction
               &                                                                         bf_alias(jp_bdya_i)%fnow(:,1,:)     !   ( a_ip = rice_apnd * a_i )
            IF( TRIM(bf_alias(jp_bdyhip)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyhip)%fnow(:,1,:) = rice_hpnd(jbdy)
            IF( TRIM(bf_alias(jp_bdyhil)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyhil)%fnow(:,1,:) = rice_hlid(jbdy)

            ! if T_i is read and not T_su, set T_su = T_i
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdytsu)%fnow(:,1,:) = bf_alias(jp_bdyt_i)%fnow(:,1,:)
            ! if T_s is read and not T_su, set T_su = T_s
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdytsu)%fnow(:,1,:) = bf_alias(jp_bdyt_s)%fnow(:,1,:)
            ! if T_i is read and not T_s, set T_s = T_i
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_s)%fnow(:,1,:) = bf_alias(jp_bdyt_i)%fnow(:,1,:)
            ! if T_su is read and not T_s, set T_s = T_su
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_s)%fnow(:,1,:) = bf_alias(jp_bdytsu)%fnow(:,1,:)
            ! if T_su is read and not T_i, set T_i = (T_su + T_freeze)/2
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_i)%fnow(:,1,:) = 0.5_wp * ( bf_alias(jp_bdytsu)%fnow(:,1,:) + 271.15 )
            ! if T_s is read and not T_i, set T_i = (T_s + T_freeze)/2
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_i)%fnow(:,1,:) = 0.5_wp * ( bf_alias(jp_bdyt_s)%fnow(:,1,:) + 271.15 )

            ! make sure ponds = 0 if no ponds scheme
            IF ( .NOT.ln_pnd ) THEN
               bf_alias(jp_bdyaip)%fnow(:,1,:) = 0._wp
               bf_alias(jp_bdyhip)%fnow(:,1,:) = 0._wp
               bf_alias(jp_bdyhil)%fnow(:,1,:) = 0._wp
            ENDIF
            IF ( .NOT.ln_pnd_lids ) THEN
               bf_alias(jp_bdyhil)%fnow(:,1,:) = 0._wp
            ENDIF
            
            ! convert N-cat fields (input) into jpl-cat (output)
            ipl = SIZE(bf_alias(jp_bdya_i)%fnow, 3)            
            IF( ipl /= jpl ) THEN      ! ice: convert N-cat fields (input) into jpl-cat (output)
               CALL ice_var_itd( bf_alias(jp_bdyh_i)%fnow(:,1,:), bf_alias(jp_bdyh_s)%fnow(:,1,:), bf_alias(jp_bdya_i)%fnow(:,1,:), & ! in
                  &              dta_alias%h_i                  , dta_alias%h_s                  , dta_alias%a_i                  , & ! out
                  &              bf_alias(jp_bdyt_i)%fnow(:,1,:), bf_alias(jp_bdyt_s)%fnow(:,1,:), &                                  ! in (optional)
                  &              bf_alias(jp_bdytsu)%fnow(:,1,:), bf_alias(jp_bdys_i)%fnow(:,1,:), &                                  ! in     -
                  &              bf_alias(jp_bdyaip)%fnow(:,1,:), bf_alias(jp_bdyhip)%fnow(:,1,:), bf_alias(jp_bdyhil)%fnow(:,1,:), & ! in     -
                  &              dta_alias%t_i                  , dta_alias%t_s                  , &                                  ! out    -
                  &              dta_alias%tsu                  , dta_alias%s_i                  , &                                  ! out    -
                  &              dta_alias%aip                  , dta_alias%hip                  , dta_alias%hil )                    ! out    -
            ENDIF
         ENDIF
#endif
      END DO  ! jbdy

      IF ( ln_tide ) THEN
         IF (ln_dynspg_ts) THEN      ! Fill temporary arrays with slow-varying bdy data                           
            DO jbdy = 1, nb_bdy      ! Tidal component added in ts loop
               IF ( nn_dyn2d_dta(jbdy) .GE. 2 ) THEN
                  IF( ASSOCIATED(dta_bdy(jbdy)%ssh) ) dta_bdy_s(jbdy)%ssh(:) = dta_bdy(jbdy)%ssh(:)
                  IF( ASSOCIATED(dta_bdy(jbdy)%u2d) ) dta_bdy_s(jbdy)%u2d(:) = dta_bdy(jbdy)%u2d(:)
                  IF( ASSOCIATED(dta_bdy(jbdy)%v2d) ) dta_bdy_s(jbdy)%v2d(:) = dta_bdy(jbdy)%v2d(:)
               ENDIF
            END DO
         ELSE ! Add tides if not split-explicit free surface else this is done in ts loop
            !
            CALL bdy_dta_tides( kt=kt, kt_offset=kt_offset )
         ENDIF
      ENDIF
      
      ! davbyr - add a shift to the boundary + free elevation Enda, JT from NEMO RAN 3.6
      DO jbdy = 1, nb_bdy
         IF( dta_bdy(jbdy)%lneed_ssh ) THEN
            igrd  = 1
            DO ib = 1, idx_bdy(jbdy)%nblenrim(igrd)   ! ssh is used only on the rim
                ii = idx_bdy(jbdy)%nbi(ib,igrd)
                ij = idx_bdy(jbdy)%nbj(ib,igrd)
                dta_bdy(jbdy)%ssh(ib) = dta_bdy(jbdy)%ssh(ib) + rn_ssh_shift(jbdy) * tmask(ii,ij,1)
                IF( .NOT. dta_bdy(jbdy)%lforced_ssh ) dta_bdy(jbdy)%ssh(ib) = sshn(ii,ij) * tmask(ii,ij,1)
             END DO
         END IF
      END DO
      !--- END davbyr
      
      !
      IF( ln_timing )   CALL timing_stop('bdy_dta')
      !
   END SUBROUTINE bdy_dta
   

   SUBROUTINE bdy_dta_init
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta_init  ***
      !!                    
      !! ** Purpose :   Initialise arrays for reading of external data 
      !!                for open boundary conditions
      !!
      !! ** Method  :   
      !!                
      !!----------------------------------------------------------------------
      INTEGER ::   jbdy, jfld    ! Local integers
      INTEGER ::   ierror, ios     ! 
      !
      CHARACTER(len=3)                       ::   cl3           ! 
      CHARACTER(len=100)                     ::   cn_dir        ! Root directory for location of data files
      LOGICAL                                ::   ln_full_vel   ! =T => full velocities in 3D boundary data
      !                                                         ! =F => baroclinic velocities in 3D boundary data
      LOGICAL                                ::   ln_zinterp    ! =T => requires a vertical interpolation of the bdydta
      REAL(wp)                               ::   rn_ice_tem, rn_ice_sal, rn_ice_age, rn_ice_apnd, rn_ice_hpnd, rn_ice_hlid
      INTEGER                                ::   ipk,ipl       !
      INTEGER                                ::   idvar         ! variable ID
      INTEGER                                ::   indims        ! number of dimensions of the variable
      INTEGER                                ::   iszdim        ! number of dimensions of the variable
      INTEGER, DIMENSION(4)                  ::   i4dimsz       ! size of variable dimensions 
      INTEGER                                ::   igrd          ! index for grid type (1,2,3 = T,U,V)
      LOGICAL                                ::   lluld         ! is the variable using the unlimited dimension
      LOGICAL                                ::   llneed        !
      LOGICAL                                ::   llread        !
      LOGICAL                                ::   llfullbdy     !
      TYPE(FLD_N), DIMENSION(1), TARGET  ::   bn_tem, bn_sal, bn_u3d, bn_v3d   ! must be an array to be used with fld_fill
      TYPE(FLD_N), DIMENSION(1), TARGET  ::   bn_ssh, bn_u2d, bn_v2d           ! informations about the fields to be read
      TYPE(FLD_N), DIMENSION(1), TARGET  ::   bn_a_i, bn_h_i, bn_h_s, bn_t_i, bn_t_s, bn_tsu, bn_s_i, bn_aip, bn_hip, bn_hil       
      TYPE(FLD_N), DIMENSION(:), POINTER ::   bn_alias                        ! must be an array to be used with fld_fill
      TYPE(FLD  ), DIMENSION(:), POINTER ::   bf_alias
      !
      NAMELIST/nambdy_dta/ cn_dir, bn_tem, bn_sal, bn_u3d, bn_v3d, bn_ssh, bn_u2d, bn_v2d,                 &
                         & bn_a_i, bn_h_i, bn_h_s, bn_t_i, bn_t_s, bn_tsu, bn_s_i, bn_aip, bn_hip, bn_hil, &
                         & rn_ice_tem, rn_ice_sal, rn_ice_age, rn_ice_apnd, rn_ice_hpnd, rn_ice_hlid,      &
                         & ln_full_vel, ln_zinterp
      !!---------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'bdy_dta_ini : initialization of data at the open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) ''

      ALLOCATE( bf(jpbdyfld,nb_bdy), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate bf structure' )   ;   RETURN  
      ENDIF
      bf(:,:)%clrootname = 'NOT USED'   ! default definition used as a flag in fld_read to do nothing.
      bf(:,:)%lzint      = .FALSE.      ! default definition
      bf(:,:)%ltotvel    = .FALSE.      ! default definition
 
      ! Read namelists
      ! --------------
      REWIND(numnam_cfg)
      DO jbdy = 1, nb_bdy

         WRITE(ctmp1, '(a,i2)') 'BDY number ', jbdy
         WRITE(ctmp2, '(a,i2)') 'block nambdy_dta number ', jbdy

         ! There is only one nambdy_dta block in namelist_ref -> use it for each bdy so we do a rewind 
         REWIND(numnam_ref)
         READ  ( numnam_ref, nambdy_dta, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy_dta in reference namelist' )

         !   by-pass nambdy_dta reading if no input data used in this bdy   
         IF(       ( dta_bdy(jbdy)%lneed_dyn2d .AND. MOD(nn_dyn2d_dta(jbdy),2) == 1 )   &
            & .OR. ( dta_bdy(jbdy)%lneed_dyn3d .AND.     nn_dyn3d_dta(jbdy)    == 1 )   &
            & .OR. ( dta_bdy(jbdy)%lneed_tra   .AND.       nn_tra_dta(jbdy)    == 1 )   &
            & .OR. ( dta_bdy(jbdy)%lneed_ice   .AND.       nn_ice_dta(jbdy)    == 1 )   )   THEN
            ! WARNING: we don't do a rewind here, each bdy reads its own nambdy_dta block one after another
            READ  ( numnam_cfg, nambdy_dta, IOSTAT = ios, ERR = 902 )
902         IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy_dta in configuration namelist' )
            IF(lwm) WRITE( numond, nambdy_dta )           
         ENDIF

         ! get the number of ice categories in bdy data file (use a_i information to do this)
         ipl = jpl   ! default definition
         IF( dta_bdy(jbdy)%lneed_ice ) THEN    ! if we need ice bdy data
            IF( nn_ice_dta(jbdy) == 1 ) THEN   ! if we get ice bdy data from netcdf file
               CALL fld_fill(  bf(jp_bdya_i,jbdy:jbdy), bn_a_i, cn_dir, 'bdy_dta', 'a_i'//' '//ctmp1, ctmp2 )   ! use namelist info
               CALL fld_clopn( bf(jp_bdya_i,jbdy), nyear, nmonth, nday )   ! not a problem when we call it again after
               idvar = iom_varid( bf(jp_bdya_i,jbdy)%num, bf(jp_bdya_i,jbdy)%clvar, kndims=indims, kdimsz=i4dimsz, lduld=lluld )
               IF( indims == 4 .OR. ( indims == 3 .AND. .NOT. lluld ) ) THEN   ;   ipl = i4dimsz(3)   ! xylt or xyl
               ELSE                                                            ;   ipl = 1            ! xy or xyt
               ENDIF
               bf(jp_bdya_i,jbdy)%clrootname = 'NOT USED'   ! reset to default value as this subdomain may not need to read this bdy
            ENDIF
         ENDIF

#if defined key_si3
         IF( .NOT.ln_pnd ) THEN
            rn_ice_apnd = 0. ; rn_ice_hpnd = 0. ; rn_ice_hlid = 0.
            CALL ctl_warn( 'rn_ice_apnd & rn_ice_hpnd = 0 & rn_ice_hlid = 0 when no ponds' )
         ENDIF
         IF( .NOT.ln_pnd_lids ) THEN
            rn_ice_hlid = 0.
         ENDIF
#endif

         ! temp, salt, age and ponds of incoming ice
         rice_tem (jbdy) = rn_ice_tem
         rice_sal (jbdy) = rn_ice_sal
         rice_age (jbdy) = rn_ice_age
         rice_apnd(jbdy) = rn_ice_apnd
         rice_hpnd(jbdy) = rn_ice_hpnd
         rice_hlid(jbdy) = rn_ice_hlid

         
         DO jfld = 1, jpbdyfld

            ! =====================
            !          ssh 
            ! =====================
            IF( jfld == jp_bdyssh ) THEN
               cl3 = 'ssh'
               igrd = 1                                                    ! T point
               ipk = 1                                                     ! surface data
               llneed = dta_bdy(jbdy)%lneed_ssh                            ! dta_bdy(jbdy)%ssh will be needed
               llread = MOD(nn_dyn2d_dta(jbdy),2) == 1                     ! get data from NetCDF file
               bf_alias => bf(jp_bdyssh,jbdy:jbdy)                         ! alias for ssh structure of bdy number jbdy
               bn_alias => bn_ssh                                          ! alias for ssh structure of nambdy_dta 
               iszdim = idx_bdy(jbdy)%nblenrim(igrd)                       ! length of this bdy on this MPI processus : used only on the rim
            ENDIF
            ! =====================
            !         dyn2d
            ! =====================
            IF( jfld == jp_bdyu2d ) THEN
               cl3 = 'u2d'
               igrd = 2                                                    ! U point
               ipk = 1                                                     ! surface data
               llneed = dta_bdy(jbdy)%lneed_dyn2d                          ! dta_bdy(jbdy)%u2d will be needed
               llread = .NOT. ln_full_vel .AND. MOD(nn_dyn2d_dta(jbdy),2) == 1   ! don't get u2d from u3d and read NetCDF file
               bf_alias => bf(jp_bdyu2d,jbdy:jbdy)                         ! alias for u2d structure of bdy number jbdy
               bn_alias => bn_u2d                                          ! alias for u2d structure of nambdy_dta
               llfullbdy = ln_full_vel .OR. cn_dyn2d(jbdy) == 'frs'        ! need u2d over the whole bdy or only over the rim?
               IF( llfullbdy ) THEN  ;   iszdim = idx_bdy(jbdy)%nblen(igrd)
               ELSE                  ;   iszdim = idx_bdy(jbdy)%nblenrim(igrd)
               ENDIF
            ENDIF
            IF( jfld == jp_bdyv2d ) THEN
               cl3 = 'v2d'
               igrd = 3                                                    ! V point
               ipk = 1                                                     ! surface data
               llneed = dta_bdy(jbdy)%lneed_dyn2d                          ! dta_bdy(jbdy)%v2d will be needed
               llread = .NOT. ln_full_vel .AND. MOD(nn_dyn2d_dta(jbdy),2) == 1   ! don't get v2d from v3d and read NetCDF file
               bf_alias => bf(jp_bdyv2d,jbdy:jbdy)                         ! alias for v2d structure of bdy number jbdy
               bn_alias => bn_v2d                                          ! alias for v2d structure of nambdy_dta 
               llfullbdy = ln_full_vel .OR. cn_dyn2d(jbdy) == 'frs'        ! need v2d over the whole bdy or only over the rim?
               IF( llfullbdy ) THEN  ;   iszdim = idx_bdy(jbdy)%nblen(igrd)
               ELSE                  ;   iszdim = idx_bdy(jbdy)%nblenrim(igrd)
               ENDIF
            ENDIF
            ! =====================
            !         dyn3d
            ! =====================
            IF( jfld == jp_bdyu3d ) THEN
               cl3 = 'u3d'
               igrd = 2                                                    ! U point
               ipk = jpk                                                   ! 3d data
               llneed = dta_bdy(jbdy)%lneed_dyn3d .OR.               &     ! dta_bdy(jbdy)%u3d will be needed
                  &   ( dta_bdy(jbdy)%lneed_dyn2d .AND. ln_full_vel )      !   u3d needed to compute u2d
               llread = nn_dyn3d_dta(jbdy) == 1                            ! get data from NetCDF file
               bf_alias => bf(jp_bdyu3d,jbdy:jbdy)                         ! alias for u3d structure of bdy number jbdy
               bn_alias => bn_u3d                                          ! alias for u3d structure of nambdy_dta 
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
           ENDIF
            IF( jfld == jp_bdyv3d ) THEN
               cl3 = 'v3d'
               igrd = 3                                                    ! V point
               ipk = jpk                                                   ! 3d data
               llneed = dta_bdy(jbdy)%lneed_dyn3d .OR.               &     ! dta_bdy(jbdy)%v3d will be needed
                  &   ( dta_bdy(jbdy)%lneed_dyn2d .AND. ln_full_vel )      !   v3d needed to compute v2d
               llread = nn_dyn3d_dta(jbdy) == 1                            ! get data from NetCDF file
               bf_alias => bf(jp_bdyv3d,jbdy:jbdy)                         ! alias for v3d structure of bdy number jbdy
               bn_alias => bn_v3d                                          ! alias for v3d structure of nambdy_dta 
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
           ENDIF

            ! =====================
            !          tra
            ! =====================
            IF( jfld == jp_bdytem ) THEN
               cl3 = 'tem'
               igrd = 1                                                    ! T point
               ipk = jpk                                                   ! 3d data
               llneed = dta_bdy(jbdy)%lneed_tra                            ! dta_bdy(jbdy)%tem will be needed
               llread = nn_tra_dta(jbdy) == 1                              ! get data from NetCDF file
               bf_alias => bf(jp_bdytem,jbdy:jbdy)                         ! alias for ssh structure of bdy number jbdy
               bn_alias => bn_tem                                          ! alias for ssh structure of nambdy_dta 
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
            ENDIF
            IF( jfld == jp_bdysal ) THEN
               cl3 = 'sal'
               igrd = 1                                                    ! T point
               ipk = jpk                                                   ! 3d data
               llneed = dta_bdy(jbdy)%lneed_tra                            ! dta_bdy(jbdy)%sal will be needed
               llread = nn_tra_dta(jbdy) == 1                              ! get data from NetCDF file
               bf_alias => bf(jp_bdysal,jbdy:jbdy)                         ! alias for ssh structure of bdy number jbdy
               bn_alias => bn_sal                                          ! alias for ssh structure of nambdy_dta 
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
            ENDIF

            ! =====================
            !          ice
            ! =====================
            IF(  jfld == jp_bdya_i .OR. jfld == jp_bdyh_i .OR. jfld == jp_bdyh_s .OR. &
               & jfld == jp_bdyt_i .OR. jfld == jp_bdyt_s .OR. jfld == jp_bdytsu .OR. &
               & jfld == jp_bdys_i .OR. jfld == jp_bdyaip .OR. jfld == jp_bdyhip .OR. jfld == jp_bdyhil ) THEN
               igrd = 1                                                    ! T point
               ipk = ipl                                                   ! jpl-cat data
               llneed = dta_bdy(jbdy)%lneed_ice                            ! ice will be needed
               llread = nn_ice_dta(jbdy) == 1                              ! get data from NetCDF file
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
            ENDIF
            IF( jfld == jp_bdya_i ) THEN
               cl3 = 'a_i'
               bf_alias => bf(jp_bdya_i,jbdy:jbdy)                         ! alias for a_i structure of bdy number jbdy
               bn_alias => bn_a_i                                          ! alias for a_i structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyh_i ) THEN
               cl3 = 'h_i'
               bf_alias => bf(jp_bdyh_i,jbdy:jbdy)                         ! alias for h_i structure of bdy number jbdy
               bn_alias => bn_h_i                                          ! alias for h_i structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyh_s ) THEN
               cl3 = 'h_s'
               bf_alias => bf(jp_bdyh_s,jbdy:jbdy)                         ! alias for h_s structure of bdy number jbdy
               bn_alias => bn_h_s                                          ! alias for h_s structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyt_i ) THEN
               cl3 = 't_i'
               bf_alias => bf(jp_bdyt_i,jbdy:jbdy)                         ! alias for t_i structure of bdy number jbdy
               bn_alias => bn_t_i                                          ! alias for t_i structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyt_s ) THEN
               cl3 = 't_s'
               bf_alias => bf(jp_bdyt_s,jbdy:jbdy)                         ! alias for t_s structure of bdy number jbdy
               bn_alias => bn_t_s                                          ! alias for t_s structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdytsu ) THEN
               cl3 = 'tsu'
               bf_alias => bf(jp_bdytsu,jbdy:jbdy)                         ! alias for tsu structure of bdy number jbdy
               bn_alias => bn_tsu                                          ! alias for tsu structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdys_i ) THEN
               cl3 = 's_i'
               bf_alias => bf(jp_bdys_i,jbdy:jbdy)                         ! alias for s_i structure of bdy number jbdy
               bn_alias => bn_s_i                                          ! alias for s_i structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyaip ) THEN
               cl3 = 'aip'
               bf_alias => bf(jp_bdyaip,jbdy:jbdy)                         ! alias for aip structure of bdy number jbdy
               bn_alias => bn_aip                                          ! alias for aip structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyhip ) THEN
               cl3 = 'hip'
               bf_alias => bf(jp_bdyhip,jbdy:jbdy)                         ! alias for hip structure of bdy number jbdy
               bn_alias => bn_hip                                          ! alias for hip structure of nambdy_dta 
            ENDIF
            IF( jfld == jp_bdyhil ) THEN
               cl3 = 'hil'
               bf_alias => bf(jp_bdyhil,jbdy:jbdy)                         ! alias for hil structure of bdy number jbdy
               bn_alias => bn_hil                                          ! alias for hil structure of nambdy_dta 
            ENDIF

            IF( llneed .AND. iszdim > 0 ) THEN                             ! dta_bdy(jbdy)%xxx will be needed
               !                                                           !   -> must be associated with an allocated target
               ALLOCATE( bf_alias(1)%fnow( iszdim, 1, ipk ) )              ! allocate the target
               !
               IF( llread ) THEN                                           ! get data from NetCDF file
                  CALL fld_fill( bf_alias, bn_alias, cn_dir, 'bdy_dta', cl3//' '//ctmp1, ctmp2 )   ! use namelist info
                  IF( bf_alias(1)%ln_tint ) ALLOCATE( bf_alias(1)%fdta( iszdim, 1, ipk, 2 ) )
                  bf_alias(1)%imap    => idx_bdy(jbdy)%nbmap(1:iszdim,igrd)   ! associate the mapping used for this bdy
                  bf_alias(1)%igrd    = igrd                                  ! used only for vertical integration of 3D arrays
                  bf_alias(1)%ibdy    = jbdy                                  !  "    "    "     "          "      "  "    "    
                  bf_alias(1)%ltotvel = ln_full_vel                           ! T if u3d is full velocity
                  bf_alias(1)%lzint   = ln_zinterp                            ! T if it requires a vertical interpolation
               ENDIF

               ! associate the pointer and get rid of the dimensions with a size equal to 1
               IF( jfld == jp_bdyssh )        dta_bdy(jbdy)%ssh => bf_alias(1)%fnow(:,1,1)
               IF( jfld == jp_bdyu2d )        dta_bdy(jbdy)%u2d => bf_alias(1)%fnow(:,1,1)
               IF( jfld == jp_bdyv2d )        dta_bdy(jbdy)%v2d => bf_alias(1)%fnow(:,1,1)
               IF( jfld == jp_bdyu3d )        dta_bdy(jbdy)%u3d => bf_alias(1)%fnow(:,1,:)
               IF( jfld == jp_bdyv3d )        dta_bdy(jbdy)%v3d => bf_alias(1)%fnow(:,1,:)
               IF( jfld == jp_bdytem )        dta_bdy(jbdy)%tem => bf_alias(1)%fnow(:,1,:)
               IF( jfld == jp_bdysal )        dta_bdy(jbdy)%sal => bf_alias(1)%fnow(:,1,:)
               IF( jfld == jp_bdya_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%a_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%a_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyh_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%h_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%h_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyh_s ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%h_s => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%h_s(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyt_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%t_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%t_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyt_s ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%t_s => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%t_s(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdytsu ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%tsu => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%tsu(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdys_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%s_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%s_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyaip ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%aip => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%aip(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyhip ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%hip => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%hip(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyhil ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%hil => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%hil(iszdim,jpl) )
                  ENDIF
               ENDIF
            ENDIF

         END DO   ! jpbdyfld
         !
      END DO ! jbdy 
      !
   END SUBROUTINE bdy_dta_init
   
   !!==============================================================================
END MODULE bdydta
