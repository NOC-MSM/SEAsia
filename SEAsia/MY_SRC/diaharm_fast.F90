MODULE diaharm_fast 
   !!======================================================================
   !!                       ***  MODULE  example  ***
   !! Ocean physics:  On line harmonic analyser
   !!                 
   !!=====================================================================

#if defined key_diaharm_fast 

   !!----------------------------------------------------------------------
   !!   'key_harm_ana'  :                Calculate harmonic analysis
   !!----------------------------------------------------------------------
   !!   harm_ana        :
   !!   harm_ana_init   :
   !!   NB: 2017-12 : add 3D harmonic analysis of velocities
   !!                 integration of Maria Luneva's development
   !!   'key_3Ddiaharm'
   !!----------------------------------------------------------------------

   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal bdy forcing
   USE daymod          ! calendar
   USE tideini
   USE restart
   USE ioipsl, ONLY : ju2ymds    ! for calendar
   !
   !
   USE timing          ! preformance summary
   USE zdf_oce

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC dia_harm_fast                                      ! routine called in step.F90 module
   LOGICAL, PUBLIC, PARAMETER :: lk_diaharm_fast  = .TRUE.   ! to be run or not
   LOGICAL, PUBLIC :: lk_diaharm_2D   ! = .TRUE.   ! to run 2d
   LOGICAL, PUBLIC :: lk_diaharm_3D   ! = .TRUE.   ! to run 3d

   !! * Module variables
   INTEGER, PARAMETER ::  nharm_max  = jpmax_harmo  ! max number of harmonics to be analysed 
   INTEGER, PARAMETER ::  nhm_max    = 2*nharm_max+1 
   INTEGER, PARAMETER ::  nvab       = 2 ! number of 3D variables
   INTEGER            ::  nharm
   INTEGER            ::  nhm 
   INTEGER ::                 & !!! ** toto namelist (namtoto) **
      nflag  =  1                ! default value of nflag 
   REAL(wp), DIMENSION(nharm_max) ::                & 
      om_tide                     ! tidal frequencies ( rads/sec)
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:)   ::                & 
      bzz,c,x    ! work arrays
   REAL(wp) :: cca,ssa,zm,bt,dd_cumul
!
   REAL(wp), PUBLIC ::   fjulday_startharm       !: Julian Day since start of harmonic analysis
   REAL(wp), PUBLIC, ALLOCATABLE,DIMENSION(:) :: anau, anav, anaf   ! nodel/phase corrections used by diaharmana
   REAL(WP), ALLOCATABLE,SAVE,DIMENSION(:,:)   :: cc,a
!
   INTEGER ::  nvar_2d, nvar_3d    !: number of 2d and 3d variables to analyse
   INTEGER, ALLOCATABLE,DIMENSION(:) :: m_posi_2d, m_posi_3d

!  Name of variables used in the restart
   CHARACTER( LEN = 10 ), DIMENSION(5), PARAMETER :: m_varName2d = (/'ssh','u2d','v2d','ubfr','vbfr'/)
   CHARACTER( LEN = 10 ), DIMENSION(4), PARAMETER :: m_varName3d = (/'rho','u3d','v3d','w3d'/)
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:  ) :: g_cosamp2D, g_sinamp2D, g_cumul_var2D
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:) :: g_cosamp3D, g_sinamp3D, g_cumul_var3D  
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:)       :: g_out2D,h_out2D  ! arrays for output
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:)     :: g_out3D,h_out3D  ! arrays for 3D output
!
!  NAMELIST
   LOGICAL, PUBLIC :: ln_diaharm_store           !: =T  Stores data for harmonic Analysis
   LOGICAL, PUBLIC :: ln_diaharm_compute         !: =T  Compute harmonic Analysis
   LOGICAL, PUBLIC :: ln_diaharm_read_restart   !: =T  Read restart from a previous run 
   LOGICAL, PUBLIC :: ln_ana_ssh, ln_ana_uvbar, ln_ana_bfric, ln_ana_rho, ln_ana_uv3d, ln_ana_w3d
   INTEGER ::   nb_ana        ! Number of harmonics to analyse
   CHARACTER (LEN=4), DIMENSION(jpmax_harmo) ::   tname   ! Names of tidal constituents ('M2', 'K1',...)
   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_all ! INDEX within the full set of constituents (tide.h90)
   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_sub ! INDEX within the subset of constituents pass in input

   !! * Substitutions

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! or LIM 2.0 , UCL-LOCEAN-IPSL (2005)
   !! or  TOP 1.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/module_example,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_harm_fast( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Harmonic analyser
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!        !  02-08  (Author names)  brief description of modifications
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * arguments
      INTEGER, INTENT( in  ) ::   &  
         kt                          ! describe it!!!

      !! * local declarations
      INTEGER  :: ji, jk, jj          ! dummy loop arguments
      INTEGER  :: jh, i1, i2, jgrid
      INTEGER  :: j2d, j3d
      REAL(WP) :: sec2start
      !!--------------------------------------------------------------------

      IF( nn_timing == 1 )   CALL timing_start( 'dia_harm_fast' )
      IF( kt == nit000   )   CALL harm_ana_init    ! Initialization (first time-step only)

     IF ( ln_diaharm_store .and. ( lk_diaharm_2D .or. lk_diaharm_3D) ) THEN

      ! this bit done every time step
      nhm=2*nb_ana+1
      c(1) = 1.0

      sec2start = nint( (fjulday-fjulday_startharm)*86400._wp ) 
      !IF(lwp) WRITE(numout,*) "ztime NEW", kt, sec2start, fjulday_startharm

      DO jh=1,nb_ana
         c(2*jh  ) = anaf(jh)*cos( sec2start*om_tide(jh) + anau(jh) + anav(jh) )
         c(2*jh+1) = anaf(jh)*sin( sec2start*om_tide(jh) + anau(jh) + anav(jh) )
      ENDDO 

      !IF(lwp) WRITE(numout,*) "c init", c, "c end", sec2start, om_tide(1), anau(1), anav(1),"end nodal"


      ! CUMULATE
      DO ji=1,jpi         ! loop lon
         DO jj=1,jpj      ! loop lat
            DO jh=1,nhm   ! loop harmonic

               DO j2d=1,nvar_2d
                  IF ( m_posi_2d(j2d) .eq. 1 ) dd_cumul = c(jh) * sshn(ji,jj) * ssmask (ji,jj)             ! analysis elevation
                  IF ( m_posi_2d(j2d) .eq. 2 ) dd_cumul = c(jh) * un_b(ji,jj) * ssumask(ji,jj)             ! analysis depth average velocities 
                  IF ( m_posi_2d(j2d) .eq. 3 ) dd_cumul = c(jh) * vn_b(ji,jj) * ssvmask(ji,jj)
                  IF ( m_posi_2d(j2d) .eq. 4 ) dd_cumul = c(jh) * bfrua(ji,jj) * un(ji,jj,mbku(ji,jj)) * ssumask(ji,jj) ! analysis bottom friction
                  IF ( m_posi_2d(j2d) .eq. 5 ) dd_cumul = c(jh) * bfrva(ji,jj) * vn(ji,jj,mbkv(ji,jj)) * ssvmask(ji,jj)
                  g_cumul_var2D(jh,ji,jj,j2d) = g_cumul_var2D(jh,ji,jj,j2d) + dd_cumul
               ENDDO

               DO j3d=1,nvar_3d
                  DO jk=1,jpkm1
                     IF ( m_posi_3d(j3d) .eq. 1 ) dd_cumul = c(jh) *  rhd(ji,jj,jk)               * tmask(ji,jj,jk)   
                     IF ( m_posi_3d(j3d) .eq. 2 ) dd_cumul = c(jh) * ( un(ji,jj,jk)-un_b(ji,jj) ) * umask(ji,jj,jk) 
                     IF ( m_posi_3d(j3d) .eq. 3 ) dd_cumul = c(jh) * ( vn(ji,jj,jk)-vn_b(ji,jj) ) * vmask(ji,jj,jk)
                     IF ( m_posi_3d(j3d) .eq. 4 ) dd_cumul = c(jh) *   wn(ji,jj,jk)               * wmask(ji,jj,jk)
                     g_cumul_var3D(jh,ji,jj,jk,j3d) = g_cumul_var3D(jh,ji,jj,jk,j3d) + dd_cumul
                  ENDDO
               ENDDO

            ENDDO     ! end loop harmonic
         ENDDO        ! end loop lat
      ENDDO           ! end loop lon

      ! Compute nodal factor cumulative cross-product
      DO i1=1,nhm
         DO i2=1,nhm
            cc(i1,i2)=cc(i1,i2)+c(i1)*c(i2)
         ENDDO
      ENDDO

      ! Output RESTART
      IF( kt == nitrst ) THEN
         CALL harm_rst_write(kt) ! Dump out data for a restarted run 
      ENDIF

      ! At End of run
      IF ( kt ==  nitend ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'harm_ana : harmonic analysis of tides at end of run'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'

         IF( ln_diaharm_compute ) THEN

             ! INITIALISE TABLE TO 0
             IF ( nvar_2d .gt. 0 ) THEN
                g_cosamp2D = 0.0_wp
                g_sinamp2D = 0.0_wp
             ENDIF
             IF ( nvar_3d .gt. 0 ) THEN
                g_cosamp3D = 0.0_wp
                g_sinamp3D = 0.0_wp
             ENDIF

             ! FIRST OUTPUT 2D VARIABLES
             DO jgrid=1,nvar_2d    ! loop number of 2d variables (ssh, U2d, V2d, UVfric) to analyse harmonically
                DO ji=1,jpi        ! loop lon
                   DO jj=1,jpj     ! loop lat
                      bt = 1.0_wp; bzz(:) = 0.0_wp
                      DO jh=1,nhm  ! loop harmonic
                         bzz(jh) = g_cumul_var2D(jh,ji,jj,jgrid)
                         bt = bt*bzz(jh)
                      ENDDO
                      ! Copy back original cumulated nodal factor
                      a(:,:) = cc(:,:)
!                     now do gaussian elimination of the system
!                     a * x = b
!                     the matrix x is (a0,a1,b1,a2,b2 ...)
!                     the matrix a and rhs b solved here for x
                      x=0.0_wp
                      IF(bt.ne.0.) THEN
                        CALL gelim( a, bzz, x, nhm )
!                       Backup output in variables
                        DO jh=1,nb_ana
                           g_cosamp2D(jh,ji,jj,jgrid) = x(jh*2  )
                           g_sinamp2D(jh,ji,jj,jgrid) = x(jh*2+1)
                        ENDDO
                        g_cosamp2D( 0,ji,jj,jgrid) = x(1)
                        g_sinamp2D( 0,ji,jj,jgrid) = 0.0_wp
                      ENDIF     ! bt.ne.0.
                   ENDDO        ! jj
                ENDDO           ! ji
             ENDDO              ! jgrid

             ! SECOND OUTPUT 3D VARIABLES
             DO jgrid=1,nvar_3d     ! loop number of 3d variables rho, U, V, W
                DO jk=1,jpkm1       ! loop over vertical level
                   DO ji=1,jpi      ! loop over lon
                      DO jj=1,jpj   ! loop over lat
                         bt = 1.0_wp; bzz(:) = 0.0_wp
                         DO jh=1,nhm
                            bzz(jh) = g_cumul_var3D(jh,ji,jj,jk,jgrid)
                            bt = bt*bzz(jh)
                         ENDDO
                         ! Copy back original cumulated nodal factor
                         a(:,:) = cc(:,:)                      
!                        now do gaussian elimination of the system
!                        a * x = b
!                        the matrix x is (a0,a1,b1,a2,b2 ...)
!                        the matrix a and rhs b solved here for x
                         x=0.0_wp
                         IF(bt.ne.0.) THEN
                           CALL gelim( a, bzz, x, nhm )
!                          Backup output in variables
                           DO jh=1,nb_ana
                              g_cosamp3D(jh,ji,jj,jk,jgrid) = x(jh*2  )
                              g_sinamp3D(jh,ji,jj,jk,jgrid) = x(jh*2+1)
                           ENDDO
                           g_cosamp3D   ( 0,ji,jj,jk,jgrid) = x(1)
                           g_sinamp3D   ( 0,ji,jj,jk,jgrid) = 0.0_wp
                        ENDIF     ! bt.ne.0.
                      ENDDO       ! jj
                   ENDDO          ! ji
                ENDDO             ! jk
             ENDDO                ! jgrid

             CALL harm_ana_out     ! output analysis (last time step)

         ELSE    ! ln_harmana_compute = False 
             IF(lwp) WRITE(numout,*) " Skipping Computing harmonics at last step"

         ENDIF   ! ln_harmana_compute 
      ENDIF      ! kt ==  nitend

     ENDIF

      IF( nn_timing == 1 )   CALL timing_stop( 'dia_harm_fast' )

   END SUBROUTINE dia_harm_fast 

   SUBROUTINE harm_ana_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::   ji, jk, jh  ! dummy loop indices
      INTEGER ::   ios                  ! Local integer output status for namelist read
      INTEGER ::   k2d, k3d             ! dummy number of analysis
      NAMELIST/nam_diaharm_fast/ ln_diaharm_store, ln_diaharm_compute, ln_diaharm_read_restart, ln_ana_ssh, ln_ana_uvbar, ln_ana_bfric, ln_ana_rho, ln_ana_uv3d, ln_ana_w3d, tname
      !!----------------------------------------------------------------------

      lk_diaharm_2D    = .TRUE.   ! to run 2d
      lk_diaharm_3D    = .TRUE.   ! to run 3d

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'harm_init : initialization of harmonic analysis of tides'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

      ! GET NAMELIST DETAILS
      REWIND( numnam_ref )              ! Namelist nam_diaharm_fast in reference namelist : Tidal harmonic analysis
      READ  ( numnam_ref, nam_diaharm_fast, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm_fast in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diaharm_fast in configuration namelist : Tidal harmonic analysis
      READ  ( numnam_cfg, nam_diaharm_fast, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm_fast in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diaharm_fast )

      ! GET NUMBER OF HARMONIC TO ANALYSE - from diaharm.F90
      nb_ana = 0
      DO jk=1,jpmax_harmo
         DO ji=1,nb_harmo
            IF(TRIM(tname(jk)) == Wave( ntide(ji) )%cname_tide ) THEN
               nb_ana=nb_ana+1
            ENDIF
         END DO
      END DO
      !
      IF(lwp) THEN
         WRITE(numout,*) '        Namelist nam_diaharm_fast'
         WRITE(numout,*) '        nb_ana    = ', nb_ana
         CALL flush(numout)
      ENDIF
      !
      IF (nb_ana > nharm_max) THEN
        IF(lwp) WRITE(numout,*) ' E R R O R harm_ana : nb_ana must be lower than nharm_max, stop'
        IF(lwp) WRITE(numout,*) ' nharm_max = ', nharm_max
        nstop = nstop + 1
      ENDIF

      ALLOCATE(ntide_all(nb_ana))
      ALLOCATE(ntide_sub(nb_ana))

      DO jk=1,nb_ana
       DO ji=1,nb_harmo
          IF (TRIM(tname(jk)) .eq. Wave( ntide(ji) )%cname_tide ) THEN
             ntide_sub(jk) = ji
             ntide_all(jk) = ntide(ji)
             EXIT
          END IF
       END DO
      END DO

      ! SEARCH HOW MANY VARIABLES 2D AND 3D TO PROCESS
      nvar_2d = 0; nvar_3d = 0
      IF ( ln_ana_ssh   ) nvar_2d = nvar_2d + 1       ! analysis elevation
      IF ( ln_ana_uvbar ) nvar_2d = nvar_2d + 2       ! analysis depth-averaged velocity
      IF ( ln_ana_bfric ) nvar_2d = nvar_2d + 2       ! analysis bottom friction 
            
      IF ( ln_ana_rho   ) nvar_3d = nvar_3d + 1       ! analysis density
      IF ( ln_ana_uv3d  ) nvar_3d = nvar_3d + 2       ! analysis 3D horizontal velocities
      IF ( ln_ana_w3d   ) nvar_3d = nvar_3d + 1       ! analysis 3D vertical velocity

      ! CHECK IF SOMETHING TO RUN
      IF ( nvar_2d .eq. 0 ) lk_diaharm_2D = .FALSE.   ! no 2d to run
      IF ( nvar_3d .eq. 0 ) lk_diaharm_3D = .FALSE.   ! no 3d to run
!      IF ( nvar_2d .gt. 0 .and. nvar_3d .gt. 0 ) lk_diaharm_fast = .FALSE.
!      IF ( .NOT. ln_diaharm_store ) lk_diaharm_fast = .FALSE.

      IF ( ln_diaharm_store .and. ( lk_diaharm_2D .or. lk_diaharm_3D) ) THEN

         ! DO ALLOCATIONS
         IF ( lk_diaharm_2D ) THEN
            ALLOCATE( g_cumul_var2D(nb_ana*2+1,jpi,jpj,    nvar_2d) )
            ALLOCATE( g_cosamp2D( 0:nb_ana*2+1,jpi,jpj,    nvar_2d) )
            ALLOCATE( g_sinamp2D( 0:nb_ana*2+1,jpi,jpj,    nvar_2d) )
            ALLOCATE( g_out2D (jpi,jpj) )
            ALLOCATE( h_out2D (jpi,jpj) )
            ALLOCATE( m_posi_2d( nvar_2d ) ); m_posi_2d(:)=0
         ENDIF
 
         IF ( lk_diaharm_3D ) THEN
            ALLOCATE( g_cumul_var3D(nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
            ALLOCATE( g_cosamp3D( 0:nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
            ALLOCATE( g_sinamp3D( 0:nb_ana*2+1,jpi,jpj,jpk,nvar_3d) )
            ALLOCATE( g_out3D (jpi,jpj,jpk) )
            ALLOCATE( h_out3D (jpi,jpj,jpk) )
            ALLOCATE( m_posi_3d( nvar_3d ) ); m_posi_3d(:)=0
         ENDIF

         ALLOCATE( cc(nb_ana*2+1,nb_ana*2+1) )
         ALLOCATE( a (nb_ana*2+1,nb_ana*2+1) )
         ALLOCATE( bzz(nb_ana*2+1) )
         ALLOCATE( x  (nb_ana*2+1) )
         ALLOCATE( c  (nb_ana*2+1) )
         ALLOCATE( anau(nb_ana) )
         ALLOCATE( anav(nb_ana) )
         ALLOCATE( anaf(nb_ana) )
         ! END ALLOCATE 

         ! STORE INDEX OF WHAT TO PRODUCE DEPENDING ON ACTIVATED LOGICAL
         ! MAKES THINGS EASIER AND FASTER LATER
         ! !!! UGLY !!!
         jh = 1; k2d = 0; 
         IF ( ln_ana_ssh   ) THEN
            k2d = k2d + 1; m_posi_2d(k2d) = jh
            IF(lwp) WRITE(numout,*) "   - ssh harmonic analysis activated (ln_ana_ssh)"
         ENDIF
         jh = jh + 1
         IF ( ln_ana_uvbar ) THEN
            k2d = k2d + 1; m_posi_2d(k2d) = jh
            jh  = jh  + 1 
            k2d = k2d + 1; m_posi_2d(k2d) = jh
            IF(lwp) WRITE(numout,*) "   - barotropic currents harmonic analysis activated (ln_ana_uvbar)"
         ELSE
            jh  = jh  + 1
         ENDIF
         jh = jh + 1
         IF ( ln_ana_bfric ) THEN
            k2d = k2d + 1; m_posi_2d(k2d) = jh
            jh  = jh  + 1; 
            k2d = k2d + 1; m_posi_2d(k2d) = jh
            IF(lwp) WRITE(numout,*) "   - bottom friction harmonic analysis activated (ln_ana_vbfr)"
         ELSE
            jh  = jh  + 1
         ENDIF

         ! and for 3D
         jh = 1; k3d = 0; 
         IF ( ln_ana_rho  ) THEN
            k3d = k3d + 1; m_posi_3d(k3d) = jh
            IF(lwp) WRITE(numout,*) "   - 3D density harmonic analysis activated (ln_ana_rho)"
         ENDIF
         jh = jh + 1
         IF ( ln_ana_uv3d )  THEN
            k3d = k3d + 1; m_posi_3d(k3d) = jh
            jh  = jh  + 1 
            k3d = k3d + 1; m_posi_3d(k3d) = jh
            IF(lwp) WRITE(numout,*) "   - 3D horizontal currents harmonic analysis activated (ln_ana_uv3d)"
         ELSE
            jh  = jh  + 1
         ENDIF
         jh = jh + 1
         IF ( ln_ana_w3d ) THEN
            k3d = k3d + 1; m_posi_3d(k3d) = jh
            IF(lwp) WRITE(numout,*) "   - 3D vertical currents harmonic analysis activated (ln_ana_w3d)"
         ENDIF

         ! SELECT AND STORE FREQUENCIES
         IF(lwp)    WRITE(numout,*) 'Analysed frequency  : ',nb_ana ,'Frequency '
         DO jh=1,nb_ana
            om_tide(jh) = omega_tide( ntide_sub(jh) ) 
            IF(lwp) WRITE(numout,*) '        - ',tname(jh),' ',om_tide(jh)
         ENDDO

         ! READ RESTART IF 
         IF ( ln_diaharm_read_restart ) THEN
            IF (lwp) WRITE(numout,*) "Reading previous harmonic data from previous run"
            ! Need to read in  bssh bz, cc anau anav and anaf 
            call harm_rst_read  ! This reads in from the previous day
                                ! Currrently the data in in assci format
         ELSE 

            IF (lwp) WRITE(numout,*) "Starting harmonic analysis from Fresh "
 
            IF ( lk_diaharm_2D ) g_cumul_var2D(:,:,:,:  ) = 0.0_wp
            IF ( lk_diaharm_3D ) g_cumul_var3D(:,:,:,:,:) = 0.0_wp
            cc           = 0.0_wp
            a    (:,:)   = 0.0_wp ! NB
            bzz  (:)     = 0.0_wp
            x    (:)     = 0.0_wp
            c    (:)     = 0.0_wp
            anau (:)     = 0.0_wp
            anav (:)     = 0.0_wp
            anaf (:)     = 0.0_wp

            DO jh = 1, nb_ana
               anau(jh) = utide ( ntide_sub(jh) )
               anav(jh) = v0tide( ntide_sub(jh) )
               anaf(jh) = ftide ( ntide_sub(jh) )
            END DO

            fjulday_startharm=fjulday !Set this at very start and store

            IF (lwp) THEN
               WRITE(numout,*) '--------------------------'
               WRITE(numout,*) '   - Output anaf for check'
               WRITE(numout,*) 'ANA F', anaf
               WRITE(numout,*) 'ANA U', anau
               WRITE(numout,*) 'ANA V', anav
               WRITE(numout,*) fjulday_startharm
               WRITE(numout,*) '--------------------------'
            ENDIF

         ENDIF

      ELSE

         IF (lwp) WRITE(numout,*) "No variable setup for harmonic analysis"

      ENDIF

   END SUBROUTINE harm_ana_init
!
   SUBROUTINE gelim (a,b,x,n)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Guassian elimination
      !!
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
        implicit none
!
        integer  :: n
        REAL(WP) :: b(nb_ana*2+1), a(nb_ana*2+1,nb_ana*2+1)
        REAL(WP) :: x(nb_ana*2+1)
        INTEGER  :: row,col,prow,pivrow,rrow
        REAL(WP) :: atemp
        REAL(WP) :: pivot
        REAL(WP) :: m

        do row=1,n-1
           pivrow=row
           pivot=a(row,n-row+1)
           do prow=row+1,n
              if (abs(a(prow,n-row+1)).gt.abs(pivot)  ) then
                 pivot=a(prow,n-row+1)
                 pivrow=prow
              endif
           enddo
!	swap row and prow
           if ( pivrow .ne. row ) then
              atemp=b(pivrow)
              b(pivrow)=b(row)
              b(row)=atemp
              do col=1,n
                 atemp=a(pivrow,col)
                 a(pivrow,col)=a(row,col)
                 a(row,col)=atemp
              enddo
           endif

           do rrow=row+1,n
              if (a(row,row).ne.0) then
   
                 m=-a(rrow,n-row+1)/a(row,n-row+1)
                 do col=1,n
                    a(rrow,col)=m*a(row,col)+a(rrow,col)
                 enddo
                 b(rrow)=m*b(row)+b(rrow)
              endif
           enddo
        enddo
!	back substitution now

        x(1)=b(n)/a(n,1)
        do row=n-1,1,-1
           x(n-row+1)=b(row)
           do col=1,(n-row)
              x(n-row+1)=(x(n-row+1)-a(row,col)*x(col)) 
           enddo

           x(n-row+1)=(x(n-row+1)/a(row,(n-row)+1))
        enddo

        return
   END SUBROUTINE gelim

   SUBROUTINE harm_ana_out
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
        USE dianam          ! build name of file (routine)
 
      !! * local declarations
      INTEGER :: ji, jj, jk, jgrid, jh    ! dummy loop indices
!      INTEGER :: nh_T
!      INTEGER :: nid_harm
!      CHARACTER (len=40) :: clhstnamt, clop1, clop2 ! temporary names 
!      CHARACTER (len=40) :: clhstnamu, clhstnamv    ! temporary names 
      CHARACTER (len=40) :: suffix
!      REAL(wp) :: zsto1, zsto2, zout, zmax, zjulian, zdt, zmdi  ! temporary scalars

      do jgrid=1,nvar_2d
          do jh=1,nb_ana
             h_out2D = 0.0
             g_out2D = 0.0
             do jj=1,nlcj
                do ji=1,nlci
                   cca=g_cosamp2D(jh,ji,jj,jgrid)
                   ssa=g_sinamp2D(jh,ji,jj,jgrid)
                   h_out2D(ji,jj)=sqrt(cca**2+ssa**2)
                   IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN 
                      g_out2D(ji,jj)= 0.0_wp
                   ELSE
                      g_out2D(ji,jj)=(180.0/rpi)*atan2(ssa,cca)       
                   ENDIF 
                   IF (h_out2D(ji,jj).ne.0) THEN
                       h_out2D(ji,jj)=h_out2D(ji,jj)/anaf(jh)
                   ENDIF
                   IF (g_out2D(ji,jj).ne.0) THEN  !Correct and take modulus
                       g_out2D(ji,jj) = g_out2D(ji,jj) + MOD( (anau(jh)+anav(jh))/rad , 360.0)
                       if (g_out2D(ji,jj).gt.360.0) then
                           g_out2D(ji,jj)=g_out2D(ji,jj)-360.0
                       else if (g_out2D(ji,jj).lt.0.0) then
                           g_out2D(ji,jj)=g_out2D(ji,jj)+360.0
                       endif
                   ENDIF
                enddo
             enddo
             !
             ! NETCDF OUTPUT
             suffix = TRIM( m_varName2d( m_posi_2d(jgrid) ) )
             CALL iom_put( TRIM(Wave(ntide_all(jh))%cname_tide)//'amp_'//TRIM(suffix), h_out2D(:,:) )
             CALL iom_put( TRIM(Wave(ntide_all(jh))%cname_tide)//'pha_'//TRIM(suffix), g_out2D(:,:) )

          enddo
      enddo
!
! DO THE SAME FOR 3D VARIABLES
!
      do jgrid=1,nvar_3d
          do jh=1,nb_ana
             h_out3D = 0.0
             g_out3D = 0.0
             DO jk=1,jpkm1
                do jj=1,nlcj
                   do ji=1,nlci
                      cca=g_cosamp3D(jh,ji,jj,jk,jgrid)
                      ssa=g_sinamp3D(jh,ji,jj,jk,jgrid)
                      h_out3D(ji,jj,jk)=sqrt(cca**2+ssa**2)
                      IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN
                         g_out3D(ji,jj,jk) = 0.0_wp
                      ELSE
                         g_out3D(ji,jj,jk) = (180.0/rpi)*atan2(ssa,cca)
                      ENDIF
                      IF (h_out3D(ji,jj,jk).ne.0) THEN
                          h_out3D(ji,jj,jk) = h_out3D(ji,jj,jk)/anaf(jh)
                      ENDIF
                      IF (g_out3D(ji,jj,jk).ne.0) THEN  !Correct and take modulus
                          g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk) + MOD( (anau(jh)+anav(jh))/rad , 360.0)
                          if      (g_out3D(ji,jj,jk).gt.360.0) then
                                   g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk)-360.0
                          else if (g_out3D(ji,jj,jk).lt.0.0) then
                                   g_out3D(ji,jj,jk) = g_out3D(ji,jj,jk)+360.0
                          endif
                      ENDIF
                   enddo    ! ji
                enddo       ! jj
             ENDDO          ! jk
             !
             ! NETCDF OUTPUT
             suffix = TRIM( m_varName3d( m_posi_3d(jgrid) ) )
             IF(lwp) WRITE(numout,*) "harm_ana_out", suffix
             CALL iom_put( TRIM(Wave(ntide_all(jh))%cname_tide)//'amp_'//TRIM(suffix), h_out3D(:,:,:) )
             CALL iom_put( TRIM(Wave(ntide_all(jh))%cname_tide)//'pha_'//TRIM(suffix), g_out3D(:,:,:) )
          enddo             ! jh 
      enddo                 ! jgrid
!
   END SUBROUTINE harm_ana_out
!
   SUBROUTINE harm_rst_write(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To write out cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   restart files will be dated by default
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      INTEGER             ::   jh, j2d, j3d
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name

      !restart file
      DO j2d=1,nvar_2d
         CALL iom_rstput( kt, nitrst, numrow, 'Mean_'//TRIM(m_varName2d( m_posi_2d(j2d) )), g_cumul_var2D( 1, :, :, j2d ) )
         DO jh=1,nb_ana
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_cos', g_cumul_var2D( jh*2  , :, :, j2d ) )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_sin', g_cumul_var2D( jh*2+1, :, :, j2d ) )
         ENDDO
      ENDDO

      DO j3d=1,nvar_3d
         CALL iom_rstput( kt, nitrst, numrow, 'Mean_'//TRIM(m_varName2d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         DO jh=1,nb_ana
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_cos', g_cumul_var3D( jh*2  , :, :, :, j3d ) )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_sin', g_cumul_var3D( jh*2+1, :, :, :, j3d ) )
         ENDDO
      ENDDO

      IF(lwp) THEN
        IF( kt > 999999999 ) THEN ; WRITE(clkt, *       ) kt
        ELSE                      ; WRITE(clkt, '(i8.8)') kt
        ENDIF
        clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
        clpath = TRIM(cn_ocerst_outdir)
        IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
        IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for writing: ',TRIM(clpath)//clname

        WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
        OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
        WRITE(66) cc
        WRITE(66) anau
        WRITE(66) anav
        WRITE(66) anaf
        WRITE(66) fjulday_startharm
        CLOSE(66)
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   harm_rst_write: DONE '
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_write

   SUBROUTINE harm_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To read in  cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name
      INTEGER             ::   jh, j2d, j3d

      IF( nit000 > 999999999 ) THEN ; WRITE(clkt, *       ) nit000-1
      ELSE                      ; WRITE(clkt, '(i8.8)') nit000-1
      ENDIF
      clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
      clpath = TRIM(cn_ocerst_outdir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'

      IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for reading: ',TRIM(clpath)//clname

      DO j2d=1,nvar_2d
         CALL iom_get( numror,jpdom_autoglo, 'Mean_'//TRIM(m_varName2d( m_posi_2d(j2d) )), g_cumul_var2D( 1, :, :, j2d ) )
         IF(lwp) WRITE(numout,*) "2D", j2d, m_posi_2d(j2d), m_varName2d( m_posi_2d(j2d) )
         DO jh=1,nb_ana
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_cos', g_cumul_var2D( jh*2  , :, :, j2d ) )
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName2d( m_posi_2d(j2d) ))//'_sin', g_cumul_var2D( jh*2+1, :, :, j2d ) )
         ENDDO
      ENDDO

      DO j3d=1,nvar_3d
         CALL iom_get( numror,jpdom_autoglo, 'Mean_'//TRIM(m_varName2d( m_posi_3d(j3d) )), g_cumul_var3D( 1, :, :, :, j3d ) )
         IF(lwp) WRITE(numout,*) "3D", j3d,  m_posi_3d(j3d), m_varName3d( m_posi_3d(j3d) )

         DO jh=1,nb_ana
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_cos', g_cumul_var3D( jh*2  , :, :, :, j3d ) )
            CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(jh))%cname_tide)//"_"//TRIM(m_varName3d( m_posi_3d(j3d) ))//'_sin', g_cumul_var3D( jh*2+1, :, :, :, j3d ) )
         ENDDO
      ENDDO

      WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
      OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
      READ(66) cc
      READ(66) anau
      READ(66) anav
      READ(66) anaf
      READ(66) fjulday_startharm
      CLOSE(66)

      IF(lwp) THEN
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   Checking anaf is correct'
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_read

   !!======================================================================
#else
!!---------------------------------------------------------------------------------
!!   Dummy module                                   NO harmonic Analysis
!!---------------------------------------------------------------------------------
        LOGICAL, PUBLIC, PARAMETER :: lk_diaharm_fast  = .FALSE.   ! to be run or not

        CONTAINS
           SUBROUTINE harm_rst_write(kt)     ! Dummy routine
           END SUBROUTINE harm_rst_write
           SUBROUTINE harm_rst_read    ! Dummy routine
           END SUBROUTINE harm_rst_read
           SUBROUTINE harm_ana_out      ! Dummy routine
           END SUBROUTINE harm_ana_out
           SUBROUTINE harm_ana_init
           END SUBROUTINE harm_ana_init
           SUBROUTINE harm_ana( kt )
!--- NB : end call not properly written
           END SUBROUTINE harm_ana
!           END SUBROUTINE harm_ana_init
!--- END NB
           SUBROUTINE gelim (a,b,x,n)
!--- NB : end call not properly written
           END SUBROUTINE gelim
!           END SUBROUTINE gelim (a,b,x,n)
!--- END NB           
#endif

END MODULE diaharm_fast 
