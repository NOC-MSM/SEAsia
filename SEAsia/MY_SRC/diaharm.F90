MODULE diaharm 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.1  !  2007  (O. Le Galloudec, J. Chanut)  Original code
   !!----------------------------------------------------------------------
#if defined key_diaharm 
   !!----------------------------------------------------------------------
   !!   'key_diaharm'
   !!
   !!   NB: 2017-12 : add 3D harmonic analysis of velocities
   !!                 integration of Maria Luneva's development
   !!   'key_3Ddiaharm'
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE phycst
   USE daymod
   USE tide_mod
   USE sbctide         ! Tidal forcing or not
   !
# if defined key_3Ddiaharm
   USE zdf_oce
#endif
   !
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE ioipsl          ! NetCDF IPSL library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! preformance summary
   USE wrk_nemo        ! working arrays

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PUBLIC, PARAMETER :: lk_diaharm  = .TRUE.
   
   INTEGER, PARAMETER :: jpincomax    = 2.*jpmax_harmo
   INTEGER, PARAMETER :: jpdimsparse  = jpincomax*300*24

   !                         !!** namelist variables **
   INTEGER ::   nit000_han    ! First time step used for harmonic analysis
   INTEGER ::   nitend_han    ! Last time step used for harmonic analysis
   INTEGER ::   nstep_han     ! Time step frequency for harmonic analysis
   INTEGER ::   nb_ana        ! Number of harmonics to analyse


   INTEGER , ALLOCATABLE, DIMENSION(:)           ::   name
   REAL(wp), ALLOCATABLE, DIMENSION(:)           ::   ana_freq, ut   , vt   , ft
# if defined key_3Ddiaharm
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:,:)   ::   ana_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)     ::   out_eta , out_u, out_v , out_w , out_dzi
# else
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)     ::   ana_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       ::   out_eta , out_u, out_v
# endif

   INTEGER ::   ninco, nsparse
   INTEGER ,       DIMENSION(jpdimsparse)         ::   njsparse, nisparse
   INTEGER , SAVE, DIMENSION(jpincomax)           ::   ipos1
   REAL(wp),       DIMENSION(jpdimsparse)         ::   valuesparse
   REAL(wp),       DIMENSION(jpincomax)           ::   ztmp4 , ztmp7
   REAL(wp), SAVE, DIMENSION(jpincomax,jpincomax) ::   ztmp3 , zpilier
   REAL(wp), SAVE, DIMENSION(jpincomax)           ::   zpivot

   CHARACTER (LEN=4), DIMENSION(jpmax_harmo) ::   tname   ! Names of tidal constituents ('M2', 'K1',...)

   PUBLIC   dia_harm   ! routine called by step.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.5 , NEMO Consortium (2013)
   !! $Id: diaharm.F90 5585 2015-07-10 14:19:11Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_harm_init 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_harm_init  ***
      !!         
      !! ** Purpose :   Initialization of tidal harmonic analysis
      !!
      !! ** Method  :   Initialize frequency array and  nodal factor for nit000_han
      !!
      !!--------------------------------------------------------------------
      INTEGER :: jh, nhan, jl
      INTEGER ::   ios                 ! Local integer output status for namelist read

      NAMELIST/nam_diaharm/ nit000_han, nitend_han, nstep_han, tname
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_harm_init: Tidal harmonic analysis initialization'
# if defined key_3Ddiaharm
         WRITE(numout,*) '  - 3D harmonic analysis of currents actovated (key_3Ddiaharm)'
#endif
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( .NOT. ln_tide )   CALL ctl_stop( 'dia_harm_init : ln_tide must be true for harmonic analysis')
      !
      CALL tide_init_Wave
      !
      REWIND( numnam_ref )              ! Namelist nam_diaharm in reference namelist : Tidal harmonic analysis
      READ  ( numnam_ref, nam_diaharm, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diaharm in configuration namelist : Tidal harmonic analysis
      READ  ( numnam_cfg, nam_diaharm, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diaharm )
      !
      IF(lwp) THEN
         WRITE(numout,*) 'First time step used for analysis:  nit000_han= ', nit000_han
         WRITE(numout,*) 'Last  time step used for analysis:  nitend_han= ', nitend_han
         WRITE(numout,*) 'Time step frequency for harmonic analysis:  nstep_han= ', nstep_han
      ENDIF

      ! Basic checks on harmonic analysis time window:
      ! ----------------------------------------------
      IF( nit000 > nit000_han )   CALL ctl_stop( 'dia_harm_init : nit000_han must be greater than nit000',   &
         &                                       ' restart capability not implemented' )
      IF( nitend < nitend_han )   CALL ctl_stop( 'dia_harm_init : nitend_han must be lower than nitend',   &
         &                                       'restart capability not implemented' )

      IF( MOD( nitend_han-nit000_han+1 , nstep_han ) /= 0 )   &
         &                        CALL ctl_stop( 'dia_harm_init : analysis time span must be a multiple of nstep_han' )

      nb_ana = 0
      DO jh=1,jpmax_harmo
         DO jl=1,jpmax_harmo
            IF(TRIM(tname(jh)) == Wave(jl)%cname_tide) THEN
               nb_ana=nb_ana+1
            ENDIF
         END DO
      END DO
      !
      IF(lwp) THEN
         WRITE(numout,*) '        Namelist nam_diaharm'
         WRITE(numout,*) '        nb_ana    = ', nb_ana
         CALL flush(numout)
      ENDIF
      !
      IF (nb_ana > jpmax_harmo) THEN
        IF(lwp) WRITE(numout,*) ' E R R O R dia_harm_init : nb_ana must be lower than jpmax_harmo, stop'
        IF(lwp) WRITE(numout,*) ' jpmax_harmo= ', jpmax_harmo
        nstop = nstop + 1
      ENDIF

      ALLOCATE(name    (nb_ana))
      DO jh=1,nb_ana
       DO jl=1,jpmax_harmo
          IF (TRIM(tname(jh)) .eq. Wave(jl)%cname_tide) THEN
             name(jh) = jl
             EXIT
          END IF
       END DO
      END DO

      ! Initialize frequency array:
      ! ---------------------------
      ALLOCATE( ana_freq(nb_ana), ut(nb_ana), vt(nb_ana), ft(nb_ana) )

      CALL tide_harmo( ana_freq, vt, ut, ft, name, nb_ana )

      IF(lwp) WRITE(numout,*) 'Analysed frequency  : ',nb_ana ,'Frequency '

      DO jh = 1, nb_ana
        IF(lwp) WRITE(numout,*) '                    : ',tname(jh),' ',ana_freq(jh)
      END DO

      ! Initialize temporary arrays:
      ! ----------------------------
# if defined key_3Ddiaharm
      ALLOCATE( ana_temp( jpi, jpj, 2*nb_ana, 5, jpk ) )
      ana_temp(:,:,:,:,:) = 0._wp
# else
      ALLOCATE( ana_temp( jpi, jpj, 2*nb_ana, 3      ) )
      ana_temp(:,:,:,:  ) = 0._wp
#endif

   END SUBROUTINE dia_harm_init


   SUBROUTINE dia_harm ( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_harm  ***
      !!         
      !! ** Purpose :   Tidal harmonic analysis main routine
      !!
      !! ** Action  :   Sums ssh/u/v over time analysis [nit000_han,nitend_han]
      !!
      !!--------------------------------------------------------------------
      INTEGER, INTENT( IN ) :: kt
      !
      INTEGER  :: ji, jj, jh, jc, nhc
# if defined key_3Ddiaharm
      INTEGER  :: jk
# endif
      REAL(wp) :: ztime, ztemp
      !!--------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_harm')

      IF( kt == nit000 ) CALL dia_harm_init

      IF( kt >= nit000_han .AND. kt <= nitend_han .AND. MOD(kt,nstep_han) == 0 ) THEN

         ztime = (kt-nit000+1) * rdt 

         !IF(lwp) WRITE(numout,*) "ztime OLD", kt, ztime, sshn(25,25)
 
         nhc = 0
         DO jh = 1, nb_ana
            DO jc = 1, 2
               nhc = nhc+1
               ztemp =(     MOD(jc,2) * ft(jh) *COS(ana_freq(jh)*ztime + vt(jh) + ut(jh))  &
                  &    +(1.-MOD(jc,2))* ft(jh) *SIN(ana_freq(jh)*ztime + vt(jh) + ut(jh)))

! ssh, ub, vb are stored at the last level of 5d array
               DO jj = 1,jpj
                  DO ji = 1,jpi
                     ! Elevation and currents
# if defined key_3Ddiaharm
                     ana_temp(ji,jj,nhc,1,jpk) = ana_temp(ji,jj,nhc,1,jpk) + ztemp*sshn(ji,jj)*ssmask (ji,jj)        
                     ana_temp(ji,jj,nhc,2,jpk) = ana_temp(ji,jj,nhc,2,jpk) + ztemp*un_b(ji,jj)*ssumask(ji,jj)
                     ana_temp(ji,jj,nhc,3,jpk) = ana_temp(ji,jj,nhc,3,jpk) + ztemp*vn_b(ji,jj)*ssvmask(ji,jj)

                     ana_temp(ji,jj,nhc,5,jpk) = ana_temp(ji,jj,nhc,5,jpk)                               &
                   &                              + ztemp*bfrva(ji,jj)*vn(ji,jj,mbkv(ji,jj))*ssvmask(ji,jj)
                     ana_temp(ji,jj,nhc,4,jpk) = ana_temp(ji,jj,nhc,4,jpk)                               & 
                   &                              + ztemp*bfrua(ji,jj)*un(ji,jj,mbku(ji,jj))*ssumask(ji,jj)
# else
                      ana_temp(ji,jj,nhc,1) = ana_temp(ji,jj,nhc,1) + ztemp*sshn(ji,jj)*ssmask (ji,jj)        
                      ana_temp(ji,jj,nhc,2) = ana_temp(ji,jj,nhc,2) + ztemp*un_b(ji,jj)*ssumask(ji,jj)
                      ana_temp(ji,jj,nhc,3) = ana_temp(ji,jj,nhc,3) + ztemp*vn_b(ji,jj)*ssvmask(ji,jj)
# endif
                  END DO
               END DO
               !
# if defined key_3Ddiaharm
! 3d velocity and density:
             DO jk=1,jpk-1
               DO jj = 1,jpj
                  DO ji = 1,jpi
                     ! density and velocity
                     ana_temp(ji,jj,nhc,1,jk) = ana_temp(ji,jj,nhc,1,jk) + ztemp*rhd(ji,jj,jk)
                     ana_temp(ji,jj,nhc,2,jk) = ana_temp(ji,jj,nhc,2,jk) + ztemp*(un(ji,jj,jk)-un_b(ji,jj)) &
                &                                          *umask(ji,jj,jk)
                     ana_temp(ji,jj,nhc,3,jk) = ana_temp(ji,jj,nhc,3,jk) + ztemp*(vn(ji,jj,jk)-vn_b(ji,jj)) &
                &                                          *vmask(ji,jj,jk) 
                     ana_temp(ji,jj,nhc,4,jk) = ana_temp(ji,jj,nhc,4,jk) + ztemp*wn(ji,jj,jk)
 
                     ana_temp(ji,jj,nhc,5,jk) = ana_temp(ji,jj,nhc,5,jk) - 0.5*grav*ztemp*(rhd(ji,jj,jk)+rhd(ji,jj,jk+1) )/max(rn2(ji,jj,jk),1.e-8_wp)
!                     IF(jk<=mbathy(ji,jj) )      ana_temp(ji,jj,nhc,5,jk) = ana_temp(ji,jj,nhc,5,jk) -      &
!                &          0.5*grav*ztemp*(rhd(ji,jj,jk)+rhd(ji,jj,jk+1) )/max(rn2(ji,jj,jk),1.e-8_wp)
                  END DO
               END DO
             ENDDO
# endif

            END DO
         END DO
         !       
      END IF

      IF ( kt == nitend_han )   CALL dia_harm_end

      IF( nn_timing == 1 )   CALL timing_stop('dia_harm')
 
   END SUBROUTINE dia_harm


   SUBROUTINE dia_harm_end
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE diaharm_end  ***
      !!         
      !! ** Purpose :  Compute the Real and Imaginary part of tidal constituents
      !!
      !! ** Action  :  Decompose the signal on the harmonic constituents 
      !!
      !!--------------------------------------------------------------------
      INTEGER :: ji, jj, jh, jc, jn, nhan, jl 
# if defined key_3Ddiaharm
      INTEGER  :: jk
# endif
      INTEGER :: ksp, kun, keq
      REAL(wp) :: ztime, ztime_ini, ztime_end
      REAL(wp) :: X1,X2
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ana_amp
      !!--------------------------------------------------------------------
      CALL wrk_alloc( jpi , jpj , jpmax_harmo , 2 , ana_amp )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'anharmo_end: kt=nitend_han: Perform harmonic analysis'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      ztime_ini = nit000_han*rdt                 ! Initial time in seconds at the beginning of analysis
      ztime_end = nitend_han*rdt                 ! Final time in seconds at the end of analysis
      nhan = (nitend_han-nit000_han+1)/nstep_han ! Number of dumps used for analysis

# if defined key_3Ddiaharm
      ALLOCATE( out_eta(jpi,jpj,jpk,2*nb_ana),   &
         &      out_u  (jpi,jpj,jpk,2*nb_ana),   &
         &      out_v  (jpi,jpj,jpk,2*nb_ana),   &
         &      out_w  (jpi,jpj,jpk,2*nb_ana),   &
         &      out_dzi(jpi,jpj,jpk,2*nb_ana) )
# else
      ALLOCATE( out_eta(jpi,jpj,2*nb_ana),   &
         &      out_u  (jpi,jpj,2*nb_ana),   &
         &      out_v  (jpi,jpj,2*nb_ana)  )
# endif

      IF(lwp) WRITE(numout,*) 'ANA F OLD', ft 
      IF(lwp) WRITE(numout,*) 'ANA U OLD', ut
      IF(lwp) WRITE(numout,*) 'ANA V OLD', vt


      ninco = 2*nb_ana
      ksp = 0
      keq = 0        
      DO jn = 1, nhan
         ztime=( (nhan-jn)*ztime_ini + (jn-1)*ztime_end )/FLOAT(nhan-1)
         keq = keq + 1
         kun = 0
         DO jh = 1, nb_ana
            DO jc = 1, 2
               kun = kun + 1
               ksp = ksp + 1
               nisparse(ksp) = keq
               njsparse(ksp) = kun
               valuesparse(ksp) = (   MOD(jc,2) * ft(jh) * COS(ana_freq(jh)*ztime + vt(jh) + ut(jh))   &
                  &             + (1.-MOD(jc,2))* ft(jh) * SIN(ana_freq(jh)*ztime + vt(jh) + ut(jh)) )
            END DO
         END DO
      END DO

      nsparse = ksp

      ! Density and Elevation:
# if defined key_3Ddiaharm
    DO jk=1,jpk
# endif
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Fill input array
            kun = 0
            DO jh = 1, nb_ana
               DO jc = 1, 2
                  kun = kun + 1
# if defined key_3Ddiaharm
                  ztmp4(kun)=ana_temp(ji,jj,kun,1,jk)
# else
                  ztmp4(kun)=ana_temp(ji,jj,kun,1)
# endif
               END DO
            END DO

            CALL SUR_DETERMINE(jj)

            ! Fill output array
            DO jh = 1, nb_ana
               ana_amp(ji,jj,jh,1)=ztmp7((jh-1)*2+1)
               ana_amp(ji,jj,jh,2)=ztmp7((jh-1)*2+2)
            END DO
         END DO
      END DO


      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jh = 1, nb_ana 
               X1 = ana_amp(ji,jj,jh,1)
               X2 =-ana_amp(ji,jj,jh,2)
# if defined key_3Ddiaharm
               out_eta(ji,jj,jk,jh       ) = X1 * tmask_i(ji,jj)
               out_eta(ji,jj,jk,jh+nb_ana) = X2 * tmask_i(ji,jj)
# else
               out_eta(ji,jj   ,jh       ) = X1 * tmask_i(ji,jj)
               out_eta(ji,jj   ,jh+nb_ana) = X2 * tmask_i(ji,jj)
# endif
            END DO
         END DO
      END DO

      ! u-component of velocity
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Fill input array
            kun=0
            DO jh = 1,nb_ana
               DO jc = 1,2
                  kun = kun + 1
# if defined key_3Ddiaharm
                  ztmp4(kun)=ana_temp(ji,jj,kun,2,jk)
# else
                  ztmp4(kun)=ana_temp(ji,jj,kun,2)
# endif
               END DO
            END DO

            CALL SUR_DETERMINE(jj+1)

            ! Fill output array
            DO jh = 1, nb_ana
               ana_amp(ji,jj,jh,1) = ztmp7((jh-1)*2+1)
               ana_amp(ji,jj,jh,2) = ztmp7((jh-1)*2+2)
            END DO

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jh = 1, nb_ana 
               X1= ana_amp(ji,jj,jh,1)
               X2=-ana_amp(ji,jj,jh,2)
# if defined key_3Ddiaharm
               out_u(ji,jj,jk,       jh) = X1 * ssumask(ji,jj)
               out_u(ji,jj,jk,nb_ana+jh) = X2 * ssumask(ji,jj)
# else
               out_u(ji,jj,          jh) = X1 * ssumask(ji,jj)
               out_u(ji,jj,   nb_ana+jh) = X2 * ssumask(ji,jj)
# endif
            ENDDO
         ENDDO
      ENDDO

      ! v- velocity
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Fill input array
            kun=0
            DO jh = 1,nb_ana
               DO jc = 1,2
                  kun = kun + 1
# if defined key_3Ddiaharm
                  ztmp4(kun)=ana_temp(ji,jj,kun,3,jk)
# else
                  ztmp4(kun)=ana_temp(ji,jj,kun,3)
# endif
               END DO
            END DO

            CALL SUR_DETERMINE(jj+1)

            ! Fill output array
            DO jh = 1, nb_ana
               ana_amp(ji,jj,jh,1)=ztmp7((jh-1)*2+1)
               ana_amp(ji,jj,jh,2)=ztmp7((jh-1)*2+2)
            END DO

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jh = 1, nb_ana 
               X1=ana_amp(ji,jj,jh,1)
               X2=-ana_amp(ji,jj,jh,2)
# if defined key_3Ddiaharm
               out_v(ji,jj,jk,       jh)=X1 * ssvmask(ji,jj)
               out_v(ji,jj,jk,nb_ana+jh)=X2 * ssvmask(ji,jj)
# else
               out_v(ji,jj,          jh)=X1 * ssvmask(ji,jj)
               out_v(ji,jj,   nb_ana+jh)=X2 * ssvmask(ji,jj)
# endif
            END DO
         END DO
      END DO

# if defined key_3Ddiaharm
      ! w- velocity
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Fill input array
            kun=0
            DO jh = 1,nb_ana
               DO jc = 1,2
                  kun = kun + 1
                  ztmp4(kun)=ana_temp(ji,jj,kun,4,jk)
               END DO
            END DO

            CALL SUR_DETERMINE(jj+1)

            ! Fill output array
            DO jh = 1, nb_ana
               ana_amp(ji,jj,jh,1)=ztmp7((jh-1)*2+1)
               ana_amp(ji,jj,jh,2)=ztmp7((jh-1)*2+2)
            END DO

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jh = 1, nb_ana
               X1=ana_amp(ji,jj,jh,1)
               X2=-ana_amp(ji,jj,jh,2)
               out_w(ji,jj,jk,       jh)=X1 * tmask_i(ji,jj)
               out_w(ji,jj,jk,nb_ana+jh)=X2 * tmask_i(ji,jj)
            END DO
         END DO
      END DO

       ! dzi- isopycnal displacements
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Fill input array
            kun=0
            DO jh = 1,nb_ana
               DO jc = 1,2
                  kun = kun + 1
                  ztmp4(kun)=ana_temp(ji,jj,kun,5,jk)
               END DO
            END DO

            CALL SUR_DETERMINE(jj+1)

            ! Fill output array
            DO jh = 1, nb_ana
               ana_amp(ji,jj,jh,1)=ztmp7((jh-1)*2+1)
               ana_amp(ji,jj,jh,2)=ztmp7((jh-1)*2+2)
            END DO

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jh = 1, nb_ana
               X1=ana_amp(ji,jj,jh,1)
               X2=-ana_amp(ji,jj,jh,2)
               out_dzi(ji,jj,jk,       jh)=X1 * tmask_i(ji,jj)
               out_dzi(ji,jj,jk,nb_ana+jh)=X2 * tmask_i(ji,jj)
            END DO
         END DO
      END DO

   ENDDO ! jk
# endif

      CALL dia_wri_harm ! Write results in files
      CALL wrk_dealloc( jpi , jpj , jpmax_harmo , 2 , ana_amp )
      !
   END SUBROUTINE dia_harm_end


   SUBROUTINE dia_wri_harm
      !!--------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_harm  ***
      !!         
      !! ** Purpose : Write tidal harmonic analysis results in a netcdf file
      !!--------------------------------------------------------------------
      CHARACTER(LEN=lc) :: cltext
      CHARACTER(LEN=lc) ::   &
         cdfile_name_T   ,   & ! name of the file created (T-points)
         cdfile_name_U   ,   & ! name of the file created (U-points)
         cdfile_name_V         ! name of the file created (V-points)
      INTEGER  ::   jh

# if defined key_3Ddiaharm
      CHARACTER(LEN=lc) :: cdfile_name_W         ! name of the file created (W-points)
      INTEGER  :: jk
      REAL(WP), ALLOCATABLE, DIMENSION (:,:,:) :: z3real, z3im 
      REAL(WP), ALLOCATABLE, DIMENSION (:,:)   :: z2real, z2im      
# endif
!!----------------------------------------------------------------------

#if defined key_dimgout
      cdfile_name_T = TRIM(cexper)//'_Tidal_harmonics_gridT.dimgproc'
      cdfile_name_U = TRIM(cexper)//'_Tidal_harmonics_gridU.dimgproc'
      cdfile_name_V = TRIM(cexper)//'_Tidal_harmonics_gridV.dimgproc'
#   if defined key_3Ddiaharm
      cdfile_name_W = TRIM(cexper)//'_Tidal_harmonics_gridW.dimgproc'
#   endif
#endif

      IF(lwp) WRITE(numout,*) '  '
      IF(lwp) WRITE(numout,*) 'dia_wri_harm : Write harmonic analysis results'
#if defined key_dimgout
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~  Output files: ', TRIM(cdfile_name_T)
      IF(lwp) WRITE(numout,*) '                             ', TRIM(cdfile_name_U)
      IF(lwp) WRITE(numout,*) '                             ', TRIM(cdfile_name_V)
#   if defined key_3Ddiaharm
      IF(lwp) WRITE(numout,*) '                             ', TRIM(cdfile_name_W)
#   endif
#endif
      IF(lwp) WRITE(numout,*) '  '

# if defined key_3Ddiaharm
      ALLOCATE( z3real(jpi,jpj,jpk),z3im(jpi,jpj,jpk),z2real(jpi,jpj),z2im(jpi,jpj))
# endif

      ! A) density and elevation
      !/////////////
      !
#if defined key_dimgout
      cltext='density amplitude and phase; elevation is level=jpk '
      CALL dia_wri_dimg(TRIM(cdfile_name_T), TRIM(cltext), out_eta, 2*nb_ana, '2')
#else
#   if defined key_3Ddiaharm
      z3real(:,:,:) = 0._wp; z3im(:,:,:) = 0._wp
#   endif
      DO jh = 1, nb_ana
#   if defined key_3Ddiaharm
        DO jk=1,jpkm1
          z3real(:,:,jk)=out_eta(:,:,jk,jh)
          z3im  (:,:,jk)=out_eta(:,:,jk,jh+nb_ana)
        ENDDO
      z2real(:,:)=out_eta(:,:,jpk,jh); z2im(:,:)=out_eta(:,:,jpk,jh+nb_ana)
      CALL iom_put( TRIM(tname(jh))//'x_ro', z3real(:,:,:) )
      CALL iom_put( TRIM(tname(jh))//'y_ro', z3im  (:,:,:) )
      CALL iom_put( TRIM(tname(jh))//'x'   , z2real(:,:  ) )
      CALL iom_put( TRIM(tname(jh))//'y'   , z2im  (:,:  ) )
#   else 
      WRITE(numout,*) "OUTPUT ORI: ", TRIM(tname(jh))//'x', ' & ', TRIM(tname(jh))//'y', MAXVAL(out_eta(:,:,jh))
      CALL iom_put( TRIM(tname(jh))//'x', out_eta(:,:,jh) )
      CALL iom_put( TRIM(tname(jh))//'y', out_eta(:,:,nb_ana+jh) )
#   endif
      END DO
#endif

      ! B) u
      !/////////
      !
#if defined key_dimgout
      cltext='3d u amplitude and phase; ubar is the last level'
      CALL dia_wri_dimg(TRIM(cdfile_name_U), TRIM(cltext), out_u, 2*nb_ana, '2')
#else
#   if defined key_3Ddiaharm
      z3real(:,:,:) = 0._wp; z3im(:,:,:) = 0._wp
#   endif
      DO jh = 1, nb_ana
#   if defined key_3Ddiaharm
        DO jk=1,jpkm1
          z3real(:,:,jk)=out_u(:,:,jk,jh)
          z3im  (:,:,jk)=out_u(:,:,jk,jh+nb_ana)
        ENDDO
        z2real(:,:)=out_u(:,:,jpk,jh); z2im(:,:)=out_u(:,:,jpk,jh+nb_ana)
        CALL iom_put( TRIM(tname(jh))//'x_u3d', z3real(:,:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_u3d', z3im (:,:,:)  )
        CALL iom_put( TRIM(tname(jh))//'x_u2d', z2real(:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_u2d', z2im (:,:)  )
        z2real(:,:)=out_w(:,:,jpk,jh); z2im(:,:)=out_w(:,:,jpk,jh+nb_ana)
        CALL iom_put( TRIM(tname(jh))//'x_tabx', z2real(:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_tabx', z2im (:,:)  )
#   else
        CALL iom_put( TRIM(tname(jh))//'x_u2d', out_u(:,:,jh) )
        CALL iom_put( TRIM(tname(jh))//'y_u2d', out_u(:,:,nb_ana+jh) )
#   endif
      END DO
#endif

      ! C) v
      !/////////
      !
#if defined key_dimgout
      cltext='3d v amplitude and phase; vbar is the last level'
      CALL dia_wri_dimg(TRIM(cdfile_name_V), TRIM(cltext), out_v, 2*nb_ana, '2')
#else
#   if defined key_3Ddiaharm
      z3real(:,:,:) = 0._wp; z3im(:,:,:) = 0._wp
#   endif
      DO jh = 1, nb_ana
#   if defined key_3Ddiaharm
        DO jk=1,jpkm1
          z3real(:,:,jk)=out_v(:,:,jk,jh)
          z3im  (:,:,jk)=out_v(:,:,jk,jh+nb_ana)
        ENDDO
        z2real(:,:)=out_v(:,:,jpk,jh); z2im(:,:)=out_v(:,:,jpk,jh+nb_ana)
        CALL iom_put( TRIM(tname(jh))//'x_v3d', z3real(:,:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_v3d', z3im (:,:,:)  )
        CALL iom_put( TRIM(tname(jh))//'x_v2d'  , z2real(:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_v2d'  , z2im (:,:)  )
        z2real(:,:)=out_dzi(:,:,jpk,jh); z2im(:,:)=out_dzi(:,:,jpk,jh+nb_ana)
        CALL iom_put( TRIM(tname(jh))//'x_taby', z2real(:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_taby', z2im (:,:)  )
#   else
         CALL iom_put( TRIM(tname(jh))//'x_v2d', out_v(:,:,jh       ) )
         CALL iom_put( TRIM(tname(jh))//'y_v2d', out_v(:,:,jh+nb_ana) )
#   endif
       END DO

#endif
      ! D) w
# if defined key_3Ddiaharm
#   if defined key_dimgout
      cltext='3d w amplitude and phase; vort_baro is the last level'
      CALL dia_wri_dimg(TRIM(cdfile_name_W), TRIM(cltext), out_w, 2*nb_ana, '2')
#   else
      DO jh = 1, nb_ana
        DO jk=1,jpkm1
         z3real(:,:,jk)=out_w(:,:,jk,jh)
         z3im(:,:,jk)=out_w(:,:,jk,jh+nb_ana)
        ENDDO
        CALL iom_put( TRIM(tname(jh))//'x_w3d', z3real(:,:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_w3d', z3im(:,:,:) )
      END DO
#   endif

!       E) dzi + tau_bot
#   if defined key_dimgout
      cltext='dzi=g*ro/N2 amplitude and phase'
      CALL dia_wri_dimg(TRIM(cdfile_name_W), TRIM(cltext), out_w, 2*nb_ana, '2')
#   else
      DO jh = 1, nb_ana
        DO jk=1,jpkm1
         z3real(:,:,jk)=out_dzi(:,:,jk,jh)
         z3im(:,:,jk)=out_dzi(:,:,jk,jh+nb_ana)
        ENDDO
        CALL iom_put( TRIM(tname(jh))//'x_dzi', z3real(:,:,:) )
        CALL iom_put( TRIM(tname(jh))//'y_dzi', z3im(:,:,:) )
      END DO
#   endif
# endif 

      !
# if defined key_3Ddiaharm
   DEALLOCATE(z3real, z3im, z2real,z2im)
# endif

   END SUBROUTINE dia_wri_harm


   SUBROUTINE SUR_DETERMINE(init)
      !!---------------------------------------------------------------------------------
      !!                      *** ROUTINE SUR_DETERMINE ***
      !!    
      !!    
      !!       
      !!---------------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   init 
      !
      INTEGER                         :: ji_sd, jj_sd, ji1_sd, ji2_sd, jk1_sd, jk2_sd
      REAL(wp)                        :: zval1, zval2, zx1
      REAL(wp), POINTER, DIMENSION(:) :: ztmpx, zcol1, zcol2
      INTEGER , POINTER, DIMENSION(:) :: ipos2, ipivot
      !---------------------------------------------------------------------------------
      CALL wrk_alloc( jpincomax , ztmpx , zcol1 , zcol2 )
      CALL wrk_alloc( jpincomax , ipos2 , ipivot        )
            
      IF( init == 1 ) THEN
         IF( nsparse > jpdimsparse )   CALL ctl_stop( 'STOP', 'SUR_DETERMINE : nsparse .GT. jpdimsparse')
         IF( ninco   > jpincomax   )   CALL ctl_stop( 'STOP', 'SUR_DETERMINE : ninco .GT. jpincomax')
         !
         ztmp3(:,:) = 0._wp
         !
         DO jk1_sd = 1, nsparse
            DO jk2_sd = 1, nsparse
               nisparse(jk2_sd) = nisparse(jk2_sd)
               njsparse(jk2_sd) = njsparse(jk2_sd)
               IF( nisparse(jk2_sd) == nisparse(jk1_sd) ) THEN
                  ztmp3(njsparse(jk1_sd),njsparse(jk2_sd)) = ztmp3(njsparse(jk1_sd),njsparse(jk2_sd))  &
                     &                                     + valuesparse(jk1_sd)*valuesparse(jk2_sd)
               ENDIF
            END DO
         END DO
         !
         DO jj_sd = 1 ,ninco
            ipos1(jj_sd) = jj_sd
            ipos2(jj_sd) = jj_sd
         END DO
         !
         DO ji_sd = 1 , ninco
            !
            !find greatest non-zero pivot:
            zval1 = ABS(ztmp3(ji_sd,ji_sd))
            !
            ipivot(ji_sd) = ji_sd
            DO jj_sd = ji_sd, ninco
               zval2 = ABS(ztmp3(ji_sd,jj_sd))
               IF( zval2.GE.zval1 )THEN
                  ipivot(ji_sd) = jj_sd
                  zval1         = zval2
               ENDIF
            END DO
            !
            DO ji1_sd = 1, ninco
               zcol1(ji1_sd)               = ztmp3(ji1_sd,ji_sd)
               zcol2(ji1_sd)               = ztmp3(ji1_sd,ipivot(ji_sd))
               ztmp3(ji1_sd,ji_sd)         = zcol2(ji1_sd)
               ztmp3(ji1_sd,ipivot(ji_sd)) = zcol1(ji1_sd)
            END DO
            !
            ipos2(ji_sd)         = ipos1(ipivot(ji_sd))
            ipos2(ipivot(ji_sd)) = ipos1(ji_sd)
            ipos1(ji_sd)         = ipos2(ji_sd)
            ipos1(ipivot(ji_sd)) = ipos2(ipivot(ji_sd))
            zpivot(ji_sd)        = ztmp3(ji_sd,ji_sd)
            DO jj_sd = 1, ninco
               ztmp3(ji_sd,jj_sd) = ztmp3(ji_sd,jj_sd) / zpivot(ji_sd)
            END DO
            !
            DO ji2_sd = ji_sd+1, ninco
               zpilier(ji2_sd,ji_sd)=ztmp3(ji2_sd,ji_sd)
               DO jj_sd=1,ninco
                  ztmp3(ji2_sd,jj_sd)=  ztmp3(ji2_sd,jj_sd) - ztmp3(ji_sd,jj_sd) * zpilier(ji2_sd,ji_sd)
               END DO
            END DO
            !
         END DO
         !
      ENDIF ! End init==1

      DO ji_sd = 1, ninco
         ztmp4(ji_sd) = ztmp4(ji_sd) / zpivot(ji_sd)
         DO ji2_sd = ji_sd+1, ninco
            ztmp4(ji2_sd) = ztmp4(ji2_sd) - ztmp4(ji_sd) * zpilier(ji2_sd,ji_sd)
         END DO
      END DO

      !system solving: 
      ztmpx(ninco) = ztmp4(ninco) / ztmp3(ninco,ninco)
      ji_sd = ninco
      DO ji_sd = ninco-1, 1, -1
         zx1 = 0._wp
         DO jj_sd = ji_sd+1, ninco
            zx1 = zx1 + ztmpx(jj_sd) * ztmp3(ji_sd,jj_sd)
         END DO
         ztmpx(ji_sd) = ztmp4(ji_sd)-zx1
      END DO

      DO jj_sd =1, ninco
         ztmp7(ipos1(jj_sd))=ztmpx(jj_sd)
      END DO

      CALL wrk_dealloc( jpincomax , ztmpx , zcol1 , zcol2 )
      CALL wrk_dealloc( jpincomax , ipos2 , ipivot        )
      !
   END SUBROUTINE SUR_DETERMINE

#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diaharm = .FALSE.
CONTAINS
   SUBROUTINE dia_harm ( kt )     ! Empty routine
      INTEGER, INTENT( IN ) :: kt  
      WRITE(*,*) 'dia_harm: you should not have seen this print'
   END SUBROUTINE dia_harm
#endif

   !!======================================================================
END MODULE diaharm
