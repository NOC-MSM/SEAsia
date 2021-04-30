MODULE zdfgls
   !!======================================================================
   !!                       ***  MODULE  zdfgls  ***
   !! Ocean physics:  vertical mixing coefficient computed from the gls 
   !!                 turbulent closure parameterization
   !!======================================================================
   !! History :  3.0  !  2009-09  (G. Reffray)  Original code
   !!            3.3  !  2010-10  (C. Bricaud)  Add in the reference
   !!            4.0  !  2017-04  (G. Madec)  remove CPP keys & avm at t-point only 
   !!             -   !  2017-05  (G. Madec)  add top friction as boundary condition
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_gls       : update momentum and tracer Kz from a gls scheme
   !!   zdf_gls_init  : initialization, namelist read, and parameters control
   !!   gls_rst       : read/write gls restart in ocean restart file
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain
   USE domvvl         ! ocean space and time domain : variable volume layer
   USE zdfdrg  , ONLY : ln_drg_OFF            ! top/bottom free-slip flag
   USE zdfdrg  , ONLY : r_z0_top , r_z0_bot   ! top/bottom roughness
   USE zdfdrg  , ONLY : rCdU_top , rCdU_bot   ! top/bottom friction
   USE sbc_oce        ! surface boundary condition: ocean
   USE phycst         ! physical constants
   USE zdfmxl         ! mixed layer
   USE sbcwave , ONLY : hsw   ! significant wave height
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP manager
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_gls        ! called in zdfphy
   PUBLIC   zdf_gls_init   ! called in zdfphy
   PUBLIC   gls_rst        ! called in zdfphy

   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hmxl_n    !: now mixing length
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zwall   !: wall function
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ustar2_surf !: Squared surface velocity scale at T-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ustar2_top  !: Squared top     velocity scale at T-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ustar2_bot  !: Squared bottom  velocity scale at T-points

   !                              !! ** Namelist  namzdf_gls  **
   LOGICAL  ::   ln_length_lim     ! use limit on the dissipation rate under stable stratification (Galperin et al. 1988)
   LOGICAL  ::   ln_sigpsi         ! Activate Burchard (2003) modification for k-eps closure & wave breaking mixing
   INTEGER  ::   nn_bc_surf        ! surface boundary condition (=0/1)
   INTEGER  ::   nn_bc_bot         ! bottom boundary condition (=0/1)
   INTEGER  ::   nn_z0_met         ! Method for surface roughness computation
   INTEGER  ::   nn_z0_ice         ! Roughness accounting for sea ice
   INTEGER  ::   nn_stab_func      ! stability functions G88, KC or Canuto (=0/1/2)
   INTEGER  ::   nn_clos           ! closure 0/1/2/3 MY82/k-eps/k-w/gen
   REAL(wp) ::   rn_clim_galp      ! Holt 2008 value for k-eps: 0.267
   REAL(wp) ::   rn_epsmin         ! minimum value of dissipation (m2/s3)
   REAL(wp) ::   rn_emin           ! minimum value of TKE (m2/s2)
   REAL(wp) ::   rn_charn          ! Charnock constant for surface breaking waves mixing : 1400. (standard) or 2.e5 (Stacey value)
   REAL(wp) ::   rn_crban          ! Craig and Banner constant for surface breaking waves mixing
   REAL(wp) ::   rn_hsro           ! Minimum surface roughness
   REAL(wp) ::   rn_hsri           ! Ice ocean roughness
   REAL(wp) ::   rn_frac_hs        ! Fraction of wave height as surface roughness (if nn_z0_met > 1) 

   REAL(wp) ::   rcm_sf        =  0.73_wp     ! Shear free turbulence parameters
   REAL(wp) ::   ra_sf         = -2.0_wp      ! Must be negative -2 < ra_sf < -1 
   REAL(wp) ::   rl_sf         =  0.2_wp      ! 0 <rl_sf<vkarmn    
   REAL(wp) ::   rghmin        = -0.28_wp
   REAL(wp) ::   rgh0          =  0.0329_wp
   REAL(wp) ::   rghcri        =  0.03_wp
   REAL(wp) ::   ra1           =  0.92_wp
   REAL(wp) ::   ra2           =  0.74_wp
   REAL(wp) ::   rb1           = 16.60_wp
   REAL(wp) ::   rb2           = 10.10_wp         
   REAL(wp) ::   re2           =  1.33_wp         
   REAL(wp) ::   rl1           =  0.107_wp
   REAL(wp) ::   rl2           =  0.0032_wp
   REAL(wp) ::   rl3           =  0.0864_wp
   REAL(wp) ::   rl4           =  0.12_wp
   REAL(wp) ::   rl5           = 11.9_wp
   REAL(wp) ::   rl6           =  0.4_wp
   REAL(wp) ::   rl7           =  0.0_wp
   REAL(wp) ::   rl8           =  0.48_wp
   REAL(wp) ::   rm1           =  0.127_wp
   REAL(wp) ::   rm2           =  0.00336_wp
   REAL(wp) ::   rm3           =  0.0906_wp
   REAL(wp) ::   rm4           =  0.101_wp
   REAL(wp) ::   rm5           = 11.2_wp
   REAL(wp) ::   rm6           =  0.4_wp
   REAL(wp) ::   rm7           =  0.0_wp
   REAL(wp) ::   rm8           =  0.318_wp
   REAL(wp) ::   rtrans        =  0.1_wp
   REAL(wp) ::   rc02, rc02r, rc03, rc04                          ! coefficients deduced from above parameters
   REAL(wp) ::   rsbc_tke1, rsbc_tke2, rfact_tke                  !     -           -           -        -
   REAL(wp) ::   rsbc_psi1, rsbc_psi2, rfact_psi                  !     -           -           -        -
   REAL(wp) ::   rsbc_zs1, rsbc_zs2                               !     -           -           -        -
   REAL(wp) ::   rc0, rc2, rc3, rf6, rcff, rc_diff                !     -           -           -        -
   REAL(wp) ::   rs0, rs1, rs2, rs4, rs5, rs6                     !     -           -           -        -
   REAL(wp) ::   rd0, rd1, rd2, rd3, rd4, rd5                     !     -           -           -        -
   REAL(wp) ::   rsc_tke, rsc_psi, rpsi1, rpsi2, rpsi3, rsc_psi0  !     -           -           -        -
   REAL(wp) ::   rpsi3m, rpsi3p, rpp, rmm, rnn                    !     -           -           -        -
   !
   REAL(wp) ::   r2_3 = 2._wp/3._wp   ! constant=2/3

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_gls_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_gls_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( hmxl_n(jpi,jpj,jpk) , ustar2_surf(jpi,jpj) ,                     &
         &      zwall (jpi,jpj,jpk) , ustar2_top (jpi,jpj) , ustar2_bot(jpi,jpj) , STAT= zdf_gls_alloc )
         !
      CALL mpp_sum ( 'zdfgls', zdf_gls_alloc )
      IF( zdf_gls_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_gls_alloc: failed to allocate arrays' )
   END FUNCTION zdf_gls_alloc


   SUBROUTINE zdf_gls( kt, p_sh2, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_gls  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!              coefficients using the GLS turbulent closure scheme.
      !!----------------------------------------------------------------------
      USE zdf_oce , ONLY : en, avtb, avmb   ! ocean vertical physics
      !!
      INTEGER                   , INTENT(in   ) ::   kt             ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   p_sh2          ! shear production term
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   p_avm, p_avt   !  momentum and tracer Kz (w-points)
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop arguments
      INTEGER  ::   ibot, ibotm1  ! local integers
      INTEGER  ::   itop, itopp1  !   -       -
      REAL(wp) ::   zesh2, zsigpsi, zcoef, zex1 , zex2  ! local scalars
      REAL(wp) ::   ztx2, zty2, zup, zdown, zcof, zdir  !   -      - 
      REAL(wp) ::   zratio, zrn2, zflxb, sh     , z_en  !   -      -
      REAL(wp) ::   prod, buoy, diss, zdiss, sm         !   -      -
      REAL(wp) ::   gh, gm, shr, dif, zsqen, zavt, zavm !   -      -
      REAL(wp) ::   zmsku, zmskv                        !   -      -
      REAL(wp), DIMENSION(jpi,jpj)     ::   zdep
      REAL(wp), DIMENSION(jpi,jpj)     ::   zkar
      REAL(wp), DIMENSION(jpi,jpj)     ::   zflxs       ! Turbulence fluxed induced by internal waves 
      REAL(wp), DIMENSION(jpi,jpj)     ::   zhsro       ! Surface roughness (surface waves)
      REAL(wp), DIMENSION(jpi,jpj)     ::   zice_fra    ! Tapering of wave breaking under sea ice
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   eb          ! tke at time before
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   hmxl_b      ! mixing length at time before
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   eps         ! dissipation rate
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwall_psi   ! Wall function use in the wb case (ln_sigpsi)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   psi         ! psi at time now
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zd_lw, zd_up, zdiag   ! lower, upper  and diagonal of the matrix
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zstt, zstm  ! stability function on tracer and momentum
      !!--------------------------------------------------------------------
      !
      ! Preliminary computing

      ustar2_surf(:,:) = 0._wp   ;         psi(:,:,:) = 0._wp   
      ustar2_top (:,:) = 0._wp   ;   zwall_psi(:,:,:) = 0._wp
      ustar2_bot (:,:) = 0._wp

      SELECT CASE ( nn_z0_ice )
      CASE( 0 )   ;   zice_fra(:,:) = 0._wp
      CASE( 1 )   ;   zice_fra(:,:) =        TANH( fr_i(:,:) * 10._wp )
      CASE( 2 )   ;   zice_fra(:,:) =              fr_i(:,:)
      CASE( 3 )   ;   zice_fra(:,:) = MIN( 4._wp * fr_i(:,:) , 1._wp )
      END SELECT
      
      ! Compute surface, top and bottom friction at T-points
      DO jj = 2, jpjm1              !==  surface ocean friction
         DO ji = fs_2, fs_jpim1           ! vector opt.         
            ustar2_surf(ji,jj) = r1_rau0 * taum(ji,jj) * tmask(ji,jj,1)
         END DO
      END DO
      !   
!!gm Rq we may add here r_ke0(_top/_bot) ?  ==>> think about that...
      !    
      IF( .NOT.ln_drg_OFF ) THEN    !== top/bottom friction   (explicit before friction)
         DO jj = 2, jpjm1                      ! bottom friction
            DO ji = fs_2, fs_jpim1   ! vector opt.         
               zmsku = 0.5*( 2._wp - umask(ji-1,jj,mbkt(ji,jj)) * umask(ji,jj,mbkt(ji,jj)) )
               zmskv = 0.5*( 2._wp - vmask(ji,jj-1,mbkt(ji,jj)) * vmask(ji,jj,mbkt(ji,jj)) )     ! (CAUTION: CdU<0)
               ustar2_bot(ji,jj) = - rCdU_bot(ji,jj) * SQRT(  ( zmsku*( ub(ji,jj,mbkt(ji,jj))+ub(ji-1,jj,mbkt(ji,jj)) ) )**2  &
                  &                                         + ( zmskv*( vb(ji,jj,mbkt(ji,jj))+vb(ji,jj-1,mbkt(ji,jj)) ) )**2  )
            END DO
         END DO
         IF( ln_isfcav ) THEN       !top friction
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zmsku = ( 2._wp - umask(ji-1,jj,mikt(ji,jj)) * umask(ji,jj,mikt(ji,jj)) )
                  zmskv = ( 2._wp - vmask(ji,jj-1,mikt(ji,jj)) * vmask(ji,jj,mikt(ji,jj)) )     ! (CAUTION: CdU<0)
                  ustar2_top(ji,jj) = - rCdU_top(ji,jj) * SQRT(  ( zmsku*( ub(ji,jj,mikt(ji,jj))+ub(ji-1,jj,mikt(ji,jj)) ) )**2  &
                     &                                         + ( zmskv*( vb(ji,jj,mikt(ji,jj))+vb(ji,jj-1,mikt(ji,jj)) ) )**2  )
               END DO
            END DO
         ENDIF
      ENDIF
   
      SELECT CASE ( nn_z0_met )      !==  Set surface roughness length  ==!
      CASE ( 0 )                          ! Constant roughness          
         zhsro(:,:) = rn_hsro
      CASE ( 1 )             ! Standard Charnock formula
         zhsro(:,:) = MAX( rsbc_zs1 * ustar2_surf(:,:) , rn_hsro )
      CASE ( 2 )             ! Roughness formulae according to Rascle et al., Ocean Modelling (2008)
!!gm faster coding : the 2 comment lines should be used
!!gm         zcof = 2._wp * 0.6_wp / 28._wp
!!gm         zdep(:,:)  = 30._wp * TANH(  zcof/ SQRT( MAX(ustar2_surf(:,:),rsmall) )  )       ! Wave age (eq. 10)
         zdep (:,:) = 30.*TANH( 2.*0.3/(28.*SQRT(MAX(ustar2_surf(:,:),rsmall))) )         ! Wave age (eq. 10)
         zhsro(:,:) = MAX(rsbc_zs2 * ustar2_surf(:,:) * zdep(:,:)**1.5, rn_hsro)          ! zhsro = rn_frac_hs * Hsw (eq. 11)
      CASE ( 3 )             ! Roughness given by the wave model (coupled or read in file)
         zhsro(:,:) = MAX(rn_frac_hs * hsw(:,:), rn_hsro)   ! (rn_frac_hs=1.6 see Eq. (5) of Rascle et al. 2008 )
      END SELECT
      !
      ! adapt roughness where there is sea ice
      zhsro(:,:) = ( (1._wp-zice_fra(:,:)) * zhsro(:,:) + zice_fra(:,:) * rn_hsri )*tmask(:,:,1)  + (1._wp - tmask(:,:,1))*rn_hsro
      !
      DO jk = 2, jpkm1              !==  Compute dissipation rate  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               eps(ji,jj,jk)  = rc03 * en(ji,jj,jk) * SQRT( en(ji,jj,jk) ) / hmxl_n(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Save tke at before time step
      eb    (:,:,:) = en    (:,:,:)
      hmxl_b(:,:,:) = hmxl_n(:,:,:)

      IF( nn_clos == 0 ) THEN    ! Mellor-Yamada
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zup   = hmxl_n(ji,jj,jk) * gdepw_n(ji,jj,mbkt(ji,jj)+1)
                  zdown = vkarmn * gdepw_n(ji,jj,jk) * ( -gdepw_n(ji,jj,jk) + gdepw_n(ji,jj,mbkt(ji,jj)+1) )
                  zcoef = ( zup / MAX( zdown, rsmall ) )
                  zwall (ji,jj,jk) = ( 1._wp + re2 * zcoef*zcoef ) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      !!---------------------------------!!
      !!   Equation to prognostic k      !!
      !!---------------------------------!!
      !
      ! Now Turbulent kinetic energy (output in en)
      ! -------------------------------
      ! Resolution of a tridiagonal linear system by a "methode de chasse"
      ! computation from level 2 to jpkm1  (e(1) computed after and e(jpk)=0 ).
      ! The surface boundary condition are set after
      ! The bottom boundary condition are also set after. In standard e(bottom)=0.
      ! zdiag : diagonal zd_up : upper diagonal zd_lw : lower diagonal
      ! Warning : after this step, en : right hand side of the matrix

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               buoy = - p_avt(ji,jj,jk) * rn2(ji,jj,jk)     ! stratif. destruction
               !
               diss = eps(ji,jj,jk)                         ! dissipation
               !
               zdir = 0.5_wp + SIGN( 0.5_wp, p_sh2(ji,jj,jk) + buoy )   ! zdir =1(=0) if shear(ji,jj,jk)+buoy >0(<0)
               !
               zesh2 = zdir*(p_sh2(ji,jj,jk)+buoy)+(1._wp-zdir)*p_sh2(ji,jj,jk)          ! production term
               zdiss = zdir*(diss/en(ji,jj,jk))   +(1._wp-zdir)*(diss-buoy)/en(ji,jj,jk) ! dissipation term
!!gm better coding, identical results
!               zesh2 =   p_sh2(ji,jj,jk) + zdir*buoy               ! production term
!               zdiss = ( diss - (1._wp-zdir)*buoy ) / en(ji,jj,jk) ! dissipation term
!!gm
               !
               ! Compute a wall function from 1. to rsc_psi*zwall/rsc_psi0
               ! Note that as long that Dirichlet boundary conditions are NOT set at the first and last levels (GOTM style)
               ! there is no need to set a boundary condition for zwall_psi at the top and bottom boundaries.
               ! Otherwise, this should be rsc_psi/rsc_psi0
               IF( ln_sigpsi ) THEN
                  zsigpsi = MIN( 1._wp, zesh2 / eps(ji,jj,jk) )     ! 0. <= zsigpsi <= 1.
                  zwall_psi(ji,jj,jk) = rsc_psi /   & 
                     &     (  zsigpsi * rsc_psi + (1._wp-zsigpsi) * rsc_psi0 / MAX( zwall(ji,jj,jk), 1._wp )  )
               ELSE
                  zwall_psi(ji,jj,jk) = 1._wp
               ENDIF
               !
               ! building the matrix
               zcof = rfact_tke * tmask(ji,jj,jk)
               !                                        ! lower diagonal, in fact not used for jk = 2 (see surface conditions)
               zd_lw(ji,jj,jk) = zcof * ( p_avm(ji,jj,jk  ) + p_avm(ji,jj,jk-1) ) / ( e3t_n(ji,jj,jk-1) * e3w_n(ji,jj,jk) )
               !                                        ! upper diagonal, in fact not used for jk = ibotm1 (see bottom conditions)
               zd_up(ji,jj,jk) = zcof * ( p_avm(ji,jj,jk+1) + p_avm(ji,jj,jk  ) ) / ( e3t_n(ji,jj,jk  ) * e3w_n(ji,jj,jk) )
               !                                        ! diagonal
               zdiag(ji,jj,jk) = 1._wp - zd_lw(ji,jj,jk) - zd_up(ji,jj,jk)  + rdt * zdiss * wmask(ji,jj,jk) 
               !                                        ! right hand side in en
               en(ji,jj,jk) = en(ji,jj,jk) + rdt * zesh2 * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      zdiag(:,:,jpk) = 1._wp
      !
      ! Set surface condition on zwall_psi (1 at the bottom)
      zwall_psi(:,:, 1 ) = zwall_psi(:,:,2)
      zwall_psi(:,:,jpk) = 1._wp
      !
      ! Surface boundary condition on tke
      ! ---------------------------------
      !
      SELECT CASE ( nn_bc_surf )
      !
      CASE ( 0 )             ! Dirichlet boundary condition (set e at k=1 & 2) 
      ! First level
      en   (:,:,1) = MAX(  rn_emin , rc02r * ustar2_surf(:,:) * (1._wp + (1._wp-zice_fra(:,:))*rsbc_tke1)**r2_3  )
      zd_lw(:,:,1) = en(:,:,1)
      zd_up(:,:,1) = 0._wp
      zdiag(:,:,1) = 1._wp
      ! 
      ! One level below
      en   (:,:,2) =  MAX(  rc02r * ustar2_surf(:,:) * (  1._wp + (1._wp-zice_fra(:,:))*rsbc_tke1 * ((zhsro(:,:)+gdepw_n(:,:,2)) &
         &                 / zhsro(:,:) )**(1.5_wp*ra_sf)  )**(2._wp/3._wp) , rn_emin   )
      zd_lw(:,:,2) = 0._wp 
      zd_up(:,:,2) = 0._wp
      zdiag(:,:,2) = 1._wp
      !
      !
      CASE ( 1 )             ! Neumann boundary condition (set d(e)/dz)
      !
      ! Dirichlet conditions at k=1
      en   (:,:,1) = MAX(  rc02r * ustar2_surf(:,:) * (1._wp + (1._wp-zice_fra(:,:))*rsbc_tke1)**r2_3 , rn_emin  )
      zd_lw(:,:,1) = en(:,:,1)
      zd_up(:,:,1) = 0._wp
      zdiag(:,:,1) = 1._wp
      !
      ! at k=2, set de/dz=Fw
      !cbr
      zdiag(:,:,2) = zdiag(:,:,2) +  zd_lw(:,:,2) ! Remove zd_lw from zdiag
      zd_lw(:,:,2) = 0._wp
      zkar (:,:)   = (rl_sf + (vkarmn-rl_sf)*(1.-EXP(-rtrans*gdept_n(:,:,1)/zhsro(:,:)) ))
      zflxs(:,:)   = rsbc_tke2 * (1._wp-zice_fra(:,:)) * ustar2_surf(:,:)**1.5_wp * zkar(:,:) &
          &                    * (  ( zhsro(:,:)+gdept_n(:,:,1) ) / zhsro(:,:)  )**(1.5_wp*ra_sf)
!!gm why not   :                        * ( 1._wp + gdept_n(:,:,1) / zhsro(:,:) )**(1.5_wp*ra_sf)
      en(:,:,2) = en(:,:,2) + zflxs(:,:) / e3w_n(:,:,2)
      !
      !
      END SELECT

      ! Bottom boundary condition on tke
      ! --------------------------------
      !
      SELECT CASE ( nn_bc_bot )
      !
      CASE ( 0 )             ! Dirichlet 
         !                      ! en(ibot) = u*^2 / Co2 and hmxl_n(ibot) = rn_lmin
         !                      ! Balance between the production and the dissipation terms
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
!!gm This means that bottom and ocean w-level above have a specified "en" value.   Sure ????
!!   With thick deep ocean level thickness, this may be quite large, no ???
!!   in particular in ocean cavities where top stratification can be large...
               ibot   = mbkt(ji,jj) + 1      ! k   bottom level of w-point
               ibotm1 = mbkt(ji,jj)          ! k-1 bottom level of w-point but >=1
               !
               z_en =  MAX( rc02r * ustar2_bot(ji,jj), rn_emin )
               !
               ! Dirichlet condition applied at: 
               !     Bottom level (ibot)      &      Just above it (ibotm1)   
               zd_lw(ji,jj,ibot) = 0._wp   ;   zd_lw(ji,jj,ibotm1) = 0._wp
               zd_up(ji,jj,ibot) = 0._wp   ;   zd_up(ji,jj,ibotm1) = 0._wp
               zdiag(ji,jj,ibot) = 1._wp   ;   zdiag(ji,jj,ibotm1) = 1._wp
               en   (ji,jj,ibot) = z_en    ;   en   (ji,jj,ibotm1) = z_en
            END DO
         END DO
         !
         IF( ln_isfcav) THEN     ! top boundary   (ocean cavity)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  itop   = mikt(ji,jj)       ! k   top w-point
                  itopp1 = mikt(ji,jj) + 1   ! k+1 1st w-point below the top one
                  !                                                ! mask at the ocean surface points
                  z_en = MAX( rc02r * ustar2_top(ji,jj), rn_emin ) * ( 1._wp - tmask(ji,jj,1) )
                  !
 !!gm TO BE VERIFIED !!!
                  ! Dirichlet condition applied at: 
                  !     top level (itop)         &      Just below it (itopp1)   
                  zd_lw(ji,jj,itop) = 0._wp   ;   zd_lw(ji,jj,itopp1) = 0._wp
                  zd_up(ji,jj,itop) = 0._wp   ;   zd_up(ji,jj,itopp1) = 0._wp
                  zdiag(ji,jj,itop) = 1._wp   ;   zdiag(ji,jj,itopp1) = 1._wp
                  en   (ji,jj,itop) = z_en    ;   en   (ji,jj,itopp1) = z_en
               END DO
            END DO
         ENDIF
         !
      CASE ( 1 )             ! Neumman boundary condition
         !                      
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ibot   = mbkt(ji,jj) + 1      ! k   bottom level of w-point
               ibotm1 = mbkt(ji,jj)          ! k-1 bottom level of w-point but >=1
               !
               z_en =  MAX( rc02r * ustar2_bot(ji,jj), rn_emin )

!CEOD This is not set in default code .. bug.
               en(ji,jj,ibot) = MAX( rc02r * ustar2_bot(ji,jj), rn_emin )

               !
               ! Bottom level Dirichlet condition:
               !     Bottom level (ibot)      &      Just above it (ibotm1)   
               !         Dirichlet            !         Neumann
               zd_lw(ji,jj,ibot) = 0._wp   !   ! Remove zd_up from zdiag
               zdiag(ji,jj,ibot) = 1._wp   ;   zdiag(ji,jj,ibotm1) = zdiag(ji,jj,ibotm1) + zd_up(ji,jj,ibotm1)
               zd_up(ji,jj,ibot) = 0._wp   ;   zd_up(ji,jj,ibotm1) = 0._wp
            END DO
         END DO
         IF( ln_isfcav) THEN     ! top boundary   (ocean cavity)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  itop   = mikt(ji,jj)       ! k   top w-point
                  itopp1 = mikt(ji,jj) + 1   ! k+1 1st w-point below the top one
                  !                                                ! mask at the ocean surface points
                  z_en = MAX( rc02r * ustar2_top(ji,jj), rn_emin ) * ( 1._wp - tmask(ji,jj,1) )
                  !
                  ! Bottom level Dirichlet condition:
                  !     Bottom level (ibot)      &      Just above it (ibotm1)   
                  !         Dirichlet            !         Neumann
                  zd_lw(ji,jj,itop) = 0._wp   !   ! Remove zd_up from zdiag
                  zdiag(ji,jj,itop) = 1._wp   ;   zdiag(ji,jj,itopp1) = zdiag(ji,jj,itopp1) + zd_up(ji,jj,itopp1)
                  zd_up(ji,jj,itop) = 0._wp   ;   zd_up(ji,jj,itopp1) = 0._wp
               END DO
            END DO
         ENDIF
         !
      END SELECT

      ! Matrix inversion (en prescribed at surface and the bottom)
      ! ----------------------------------------------------------
      !
      DO jk = 2, jpkm1                             ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               zdiag(ji,jj,jk) = zdiag(ji,jj,jk) - zd_lw(ji,jj,jk) * zd_up(ji,jj,jk-1) / zdiag(ji,jj,jk-1)
            END DO
         END DO
      END DO
      DO jk = 2, jpkm1                             ! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               zd_lw(ji,jj,jk) = en(ji,jj,jk) - zd_lw(ji,jj,jk) / zdiag(ji,jj,jk-1) * zd_lw(ji,jj,jk-1)
            END DO
         END DO
      END DO
      DO jk = jpkm1, 2, -1                         ! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               en(ji,jj,jk) = ( zd_lw(ji,jj,jk) - zd_up(ji,jj,jk) * en(ji,jj,jk+1) ) / zdiag(ji,jj,jk)
            END DO
         END DO
      END DO
      !                                            ! set the minimum value of tke 
      en(:,:,:) = MAX( en(:,:,:), rn_emin )

      !!----------------------------------------!!
      !!   Solve prognostic equation for psi    !!
      !!----------------------------------------!!

      ! Set psi to previous time step value
      !
      SELECT CASE ( nn_clos )
      !
      CASE( 0 )               ! k-kl  (Mellor-Yamada)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  psi(ji,jj,jk)  = eb(ji,jj,jk) * hmxl_b(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      CASE( 1 )               ! k-eps
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  psi(ji,jj,jk)  = eps(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      CASE( 2 )               ! k-w
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  psi(ji,jj,jk)  = SQRT( eb(ji,jj,jk) ) / ( rc0 * hmxl_b(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
      CASE( 3 )               ! generic
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  psi(ji,jj,jk)  = rc02 * eb(ji,jj,jk) * hmxl_b(ji,jj,jk)**rnn 
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      ! Now gls (output in psi)
      ! -------------------------------
      ! Resolution of a tridiagonal linear system by a "methode de chasse"
      ! computation from level 2 to jpkm1  (e(1) already computed and e(jpk)=0 ).
      ! zdiag : diagonal zd_up : upper diagonal zd_lw : lower diagonal
      ! Warning : after this step, en : right hand side of the matrix

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               !
               ! psi / k
               zratio = psi(ji,jj,jk) / eb(ji,jj,jk) 
               !
               ! psi3+ : stable : B=-KhN²<0 => N²>0 if rn2>0 zdir = 1 (stable) otherwise zdir = 0 (unstable)
               zdir = 0.5_wp + SIGN( 0.5_wp, rn2(ji,jj,jk) )
               !
               rpsi3 = zdir * rpsi3m + ( 1._wp - zdir ) * rpsi3p
               !
               ! shear prod. - stratif. destruction
               prod = rpsi1 * zratio * p_sh2(ji,jj,jk)
               !
               ! stratif. destruction
               buoy = rpsi3 * zratio * (- p_avt(ji,jj,jk) * rn2(ji,jj,jk) )
               !
               ! shear prod. - stratif. destruction
               diss = rpsi2 * zratio * zwall(ji,jj,jk) * eps(ji,jj,jk)
               !
               zdir = 0.5_wp + SIGN( 0.5_wp, prod + buoy )     ! zdir =1(=0) if shear(ji,jj,jk)+buoy >0(<0)
               !
               zesh2 = zdir * ( prod + buoy )          + (1._wp - zdir ) * prod                        ! production term
               zdiss = zdir * ( diss / psi(ji,jj,jk) ) + (1._wp - zdir ) * (diss-buoy) / psi(ji,jj,jk) ! dissipation term
               !                                                        
               ! building the matrix
               zcof = rfact_psi * zwall_psi(ji,jj,jk) * tmask(ji,jj,jk)
               !                                               ! lower diagonal
               zd_lw(ji,jj,jk) = zcof * ( p_avm(ji,jj,jk  ) + p_avm(ji,jj,jk-1) ) / ( e3t_n(ji,jj,jk-1) * e3w_n(ji,jj,jk) )
               !                                               ! upper diagonal
               zd_up(ji,jj,jk) = zcof * ( p_avm(ji,jj,jk+1) + p_avm(ji,jj,jk  ) ) / ( e3t_n(ji,jj,jk  ) * e3w_n(ji,jj,jk) )
               !                                               ! diagonal
               zdiag(ji,jj,jk) = 1._wp - zd_lw(ji,jj,jk) - zd_up(ji,jj,jk) + rdt * zdiss * wmask(ji,jj,jk)
               !                                               ! right hand side in psi
               psi(ji,jj,jk) = psi(ji,jj,jk) + rdt * zesh2 * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      zdiag(:,:,jpk) = 1._wp

      ! Surface boundary condition on psi
      ! ---------------------------------
      !
      SELECT CASE ( nn_bc_surf )
      !
      CASE ( 0 )             ! Dirichlet boundary conditions
         !
         ! Surface value
         zdep    (:,:)   = zhsro(:,:) * rl_sf ! Cosmetic
         psi     (:,:,1) = rc0**rpp * en(:,:,1)**rmm * zdep(:,:)**rnn * tmask(:,:,1)
         zd_lw(:,:,1) = psi(:,:,1)
         zd_up(:,:,1) = 0._wp
         zdiag(:,:,1) = 1._wp
         !
         ! One level below
         zkar    (:,:)   = (rl_sf + (vkarmn-rl_sf)*(1._wp-EXP(-rtrans*gdepw_n(:,:,2)/zhsro(:,:) )))
         zdep    (:,:)   = (zhsro(:,:) + gdepw_n(:,:,2)) * zkar(:,:)
         psi     (:,:,2) = rc0**rpp * en(:,:,2)**rmm * zdep(:,:)**rnn * tmask(:,:,1)
         zd_lw(:,:,2) = 0._wp
         zd_up(:,:,2) = 0._wp
         zdiag(:,:,2) = 1._wp
         ! 
      CASE ( 1 )             ! Neumann boundary condition on d(psi)/dz
         !
         ! Surface value: Dirichlet
         zdep    (:,:)   = zhsro(:,:) * rl_sf
         psi     (:,:,1) = rc0**rpp * en(:,:,1)**rmm * zdep(:,:)**rnn * tmask(:,:,1)
         zd_lw(:,:,1) = psi(:,:,1)
         zd_up(:,:,1) = 0._wp
         zdiag(:,:,1) = 1._wp
         !
         ! Neumann condition at k=2
         zdiag(:,:,2) = zdiag(:,:,2) +  zd_lw(:,:,2) ! Remove zd_lw from zdiag
         zd_lw(:,:,2) = 0._wp
         !
         ! Set psi vertical flux at the surface:
         zkar (:,:)   = rl_sf + (vkarmn-rl_sf)*(1._wp-EXP(-rtrans*gdept_n(:,:,1)/zhsro(:,:) )) ! Lengh scale slope
         zdep (:,:)   = ((zhsro(:,:) + gdept_n(:,:,1)) / zhsro(:,:))**(rmm*ra_sf)
         zflxs(:,:)   = (rnn + (1._wp-zice_fra(:,:))*rsbc_tke1 * (rnn + rmm*ra_sf) * zdep(:,:)) &
            &           *(1._wp + (1._wp-zice_fra(:,:))*rsbc_tke1*zdep(:,:))**(2._wp*rmm/3._wp-1_wp)
         zdep (:,:)   = rsbc_psi1 * (zwall_psi(:,:,1)*p_avm(:,:,1)+zwall_psi(:,:,2)*p_avm(:,:,2)) * &
            &           ustar2_surf(:,:)**rmm * zkar(:,:)**rnn * (zhsro(:,:) + gdept_n(:,:,1))**(rnn-1.)
         zflxs(:,:)   = zdep(:,:) * zflxs(:,:)
         psi  (:,:,2) = psi(:,:,2) + zflxs(:,:) / e3w_n(:,:,2)
         !
      END SELECT

      ! Bottom boundary condition on psi
      ! --------------------------------
      !
!!gm should be done for ISF (top boundary cond.)
!!gm so, totally new staff needed      ===>>> think about that !
!
      SELECT CASE ( nn_bc_bot )     ! bottom boundary
      !
      CASE ( 0 )             ! Dirichlet 
         !                      ! en(ibot) = u*^2 / Co2 and hmxl_n(ibot) = vkarmn * r_z0_bot
         !                      ! Balance between the production and the dissipation terms
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ibot   = mbkt(ji,jj) + 1      ! k   bottom level of w-point
               ibotm1 = mbkt(ji,jj)          ! k-1 bottom level of w-point but >=1
               zdep(ji,jj) = vkarmn * r_z0_bot
               psi (ji,jj,ibot) = rc0**rpp * en(ji,jj,ibot)**rmm * zdep(ji,jj)**rnn
               zd_lw(ji,jj,ibot) = 0._wp
               zd_up(ji,jj,ibot) = 0._wp
               zdiag(ji,jj,ibot) = 1._wp
               !
               ! Just above last level, Dirichlet condition again (GOTM like)
               zdep(ji,jj) = vkarmn * ( r_z0_bot + e3t_n(ji,jj,ibotm1) )
               psi (ji,jj,ibotm1) = rc0**rpp * en(ji,jj,ibot  )**rmm * zdep(ji,jj)**rnn
               zd_lw(ji,jj,ibotm1) = 0._wp
               zd_up(ji,jj,ibotm1) = 0._wp
               zdiag(ji,jj,ibotm1) = 1._wp
            END DO
         END DO
         !
      CASE ( 1 )             ! Neumman boundary condition
         !                      
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ibot   = mbkt(ji,jj) + 1      ! k   bottom level of w-point
               ibotm1 = mbkt(ji,jj)          ! k-1 bottom level of w-point but >=1
               !
               ! Bottom level Dirichlet condition:
               zdep(ji,jj) = vkarmn * r_z0_bot
               psi (ji,jj,ibot) = rc0**rpp * en(ji,jj,ibot)**rmm * zdep(ji,jj)**rnn
               !
               zd_lw(ji,jj,ibot) = 0._wp
               zd_up(ji,jj,ibot) = 0._wp
               zdiag(ji,jj,ibot) = 1._wp
               !
               ! Just above last level: Neumann condition with flux injection
               zdiag(ji,jj,ibotm1) = zdiag(ji,jj,ibotm1) + zd_up(ji,jj,ibotm1) ! Remove zd_up from zdiag
               zd_up(ji,jj,ibotm1) = 0.
               !
               ! Set psi vertical flux at the bottom:
               zdep(ji,jj) = r_z0_bot + 0.5_wp*e3t_n(ji,jj,ibotm1)
               zflxb = rsbc_psi2 * ( p_avm(ji,jj,ibot) + p_avm(ji,jj,ibotm1) )   &
                  &  * (0.5_wp*(en(ji,jj,ibot)+en(ji,jj,ibotm1)))**rmm * zdep(ji,jj)**(rnn-1._wp)
               psi(ji,jj,ibotm1) = psi(ji,jj,ibotm1) + zflxb / e3w_n(ji,jj,ibotm1)
            END DO
         END DO
         !
      END SELECT

      ! Matrix inversion
      ! ----------------
      !
      DO jk = 2, jpkm1                             ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               zdiag(ji,jj,jk) = zdiag(ji,jj,jk) - zd_lw(ji,jj,jk) * zd_up(ji,jj,jk-1) / zdiag(ji,jj,jk-1)
            END DO
         END DO
      END DO
      DO jk = 2, jpkm1                             ! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               zd_lw(ji,jj,jk) = psi(ji,jj,jk) - zd_lw(ji,jj,jk) / zdiag(ji,jj,jk-1) * zd_lw(ji,jj,jk-1)
            END DO
         END DO
      END DO
      DO jk = jpkm1, 2, -1                         ! Third recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               psi(ji,jj,jk) = ( zd_lw(ji,jj,jk) - zd_up(ji,jj,jk) * psi(ji,jj,jk+1) ) / zdiag(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Set dissipation
      !----------------

      SELECT CASE ( nn_clos )
      !
      CASE( 0 )               ! k-kl  (Mellor-Yamada)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  eps(ji,jj,jk) = rc03 * en(ji,jj,jk) * en(ji,jj,jk) * SQRT( en(ji,jj,jk) ) / MAX( psi(ji,jj,jk), rn_epsmin)
               END DO
            END DO
         END DO
         !
      CASE( 1 )               ! k-eps
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  eps(ji,jj,jk) = psi(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      CASE( 2 )               ! k-w
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  eps(ji,jj,jk) = rc04 * en(ji,jj,jk) * psi(ji,jj,jk) 
               END DO
            END DO
         END DO
         !
      CASE( 3 )               ! generic
         zcoef = rc0**( 3._wp  + rpp/rnn )
         zex1  =      ( 1.5_wp + rmm/rnn )
         zex2  = -1._wp / rnn
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  eps(ji,jj,jk) = zcoef * en(ji,jj,jk)**zex1 * psi(ji,jj,jk)**zex2
               END DO
            END DO
         END DO
         !
      END SELECT

      ! Limit dissipation rate under stable stratification
      ! --------------------------------------------------
      DO jk = 1, jpkm1 ! Note that this set boundary conditions on hmxl_n at the same time
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               ! limitation
               eps   (ji,jj,jk)  = MAX( eps(ji,jj,jk), rn_epsmin )
               hmxl_n(ji,jj,jk)  = rc03 * en(ji,jj,jk) * SQRT( en(ji,jj,jk) ) / eps(ji,jj,jk)
               ! Galperin criterium (NOTE : Not required if the proper value of C3 in stable cases is calculated) 
               zrn2 = MAX( rn2(ji,jj,jk), rsmall )
               IF( ln_length_lim )   hmxl_n(ji,jj,jk) = MIN(  rn_clim_galp * SQRT( 2._wp * en(ji,jj,jk) / zrn2 ), hmxl_n(ji,jj,jk) )
            END DO
         END DO
      END DO 

      !
      ! Stability function and vertical viscosity and diffusivity
      ! ---------------------------------------------------------
      !
      SELECT CASE ( nn_stab_func )
      !
      CASE ( 0 , 1 )             ! Galperin or Kantha-Clayson stability functions
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! zcof =  l²/q²
                  zcof = hmxl_b(ji,jj,jk) * hmxl_b(ji,jj,jk) / ( 2._wp*eb(ji,jj,jk) )
                  ! Gh = -N²l²/q²
                  gh = - rn2(ji,jj,jk) * zcof
                  gh = MIN( gh, rgh0   )
                  gh = MAX( gh, rghmin )
                  ! Stability functions from Kantha and Clayson (if C2=C3=0 => Galperin)
                  sh = ra2*( 1._wp-6._wp*ra1/rb1 ) / ( 1.-3.*ra2*gh*(6.*ra1+rb2*( 1._wp-rc3 ) ) )
                  sm = ( rb1**(-1._wp/3._wp) + ( 18._wp*ra1*ra1 + 9._wp*ra1*ra2*(1._wp-rc2) )*sh*gh ) / (1._wp-9._wp*ra1*ra2*gh)
                  !
                  ! Store stability function in zstt and zstm
                  zstt(ji,jj,jk) = rc_diff * sh * tmask(ji,jj,jk)
                  zstm(ji,jj,jk) = rc_diff * sm * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      CASE ( 2, 3 )               ! Canuto stability functions
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! zcof =  l²/q²
                  zcof = hmxl_b(ji,jj,jk)*hmxl_b(ji,jj,jk) / ( 2._wp * eb(ji,jj,jk) )
                  ! Gh = -N²l²/q²
                  gh = - rn2(ji,jj,jk) * zcof
                  gh = MIN( gh, rgh0   )
                  gh = MAX( gh, rghmin )
                  gh = gh * rf6
                  ! Gm =  M²l²/q² Shear number
                  shr = p_sh2(ji,jj,jk) / MAX( p_avm(ji,jj,jk), rsmall )
                  gm = MAX( shr * zcof , 1.e-10 )
                  gm = gm * rf6
                  gm = MIN ( (rd0 - rd1*gh + rd3*gh*gh) / (rd2-rd4*gh) , gm )
                  ! Stability functions from Canuto
                  rcff = rd0 - rd1*gh +rd2*gm + rd3*gh*gh - rd4*gh*gm + rd5*gm*gm
                  sm = (rs0 - rs1*gh + rs2*gm) / rcff
                  sh = (rs4 - rs5*gh + rs6*gm) / rcff
                  !
                  ! Store stability function in zstt and zstm
                  zstt(ji,jj,jk) = rc_diff * sh * tmask(ji,jj,jk)
                  zstm(ji,jj,jk) = rc_diff * sm * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      END SELECT

      ! Boundary conditions on stability functions for momentum (Neumann):
      ! Lines below are useless if GOTM style Dirichlet conditions are used

      zstm(:,:,1) = zstm(:,:,2)

      ! default value, in case jpk > mbkt(ji,jj)+1. Not needed but avoid a bug when looking for undefined values (-fpe0)
      zstm(:,:,jpk) = 0.  
      DO jj = 2, jpjm1                ! update bottom with good values
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zstm(ji,jj,mbkt(ji,jj)+1) = zstm(ji,jj,mbkt(ji,jj))
         END DO
      END DO

      zstt(:,:,  1) = wmask(:,:,  1)  ! default value not needed but avoid a bug when looking for undefined values (-fpe0)
      zstt(:,:,jpk) = wmask(:,:,jpk)  ! default value not needed but avoid a bug when looking for undefined values (-fpe0)

!!gm should be done for ISF (top boundary cond.)
!!gm so, totally new staff needed!!gm

      ! Compute diffusivities/viscosities
      ! The computation below could be restrained to jk=2 to jpkm1 if GOTM style Dirichlet conditions are used
      !  -> yes BUT p_avm(:,:1) and p_avm(:,:jpk) are used when we compute zd_lw(:,:2) and zd_up(:,:jpkm1). These values are
      !     later overwritten by surface/bottom boundaries conditions, so we don't really care of p_avm(:,:1) and p_avm(:,:jpk)
      !     for zd_lw and zd_up but they have to be defined to avoid a bug when looking for undefined values (-fpe0)
      DO jk = 1, jpk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zsqen = SQRT( 2._wp * en(ji,jj,jk) ) * hmxl_n(ji,jj,jk)
               zavt  = zsqen * zstt(ji,jj,jk)
               zavm  = zsqen * zstm(ji,jj,jk)
               p_avt(ji,jj,jk) = MAX( zavt, avtb(jk) ) * wmask(ji,jj,jk) ! apply mask for zdfmxl routine
               p_avm(ji,jj,jk) = MAX( zavm, avmb(jk) )                   ! Note that avm is not masked at the surface and the bottom
            END DO
         END DO
      END DO
      p_avt(:,:,1) = 0._wp
      !
      IF(ln_ctl) THEN
         CALL prt_ctl( tab3d_1=en   , clinfo1=' gls  - e: ', tab3d_2=p_avt, clinfo2=' t: ', kdim=jpk)
         CALL prt_ctl( tab3d_1=p_avm, clinfo1=' gls  - m: ', kdim=jpk )
      ENDIF
      !
   END SUBROUTINE zdf_gls


   SUBROUTINE zdf_gls_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_gls_init  ***
      !!                     
      !! ** Purpose :   Initialization of the vertical eddy diffivity and 
      !!              viscosity computed using a GLS turbulent closure scheme
      !!
      !! ** Method  :   Read the namzdf_gls namelist and check the parameters
      !!
      !! ** input   :   Namlist namzdf_gls
      !!
      !! ** Action  :   Increase by 1 the nstop flag is setting problem encounter
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   jk    ! dummy loop indices
      INTEGER ::   ios   ! Local integer output status for namelist read
      REAL(wp)::   zcr   ! local scalar
      !!
      NAMELIST/namzdf_gls/rn_emin, rn_epsmin, ln_length_lim,       &
         &            rn_clim_galp, ln_sigpsi, rn_hsro, rn_hsri,   &
         &            rn_crban, rn_charn, rn_frac_hs,              &
         &            nn_bc_surf, nn_bc_bot, nn_z0_met, nn_z0_ice, &
         &            nn_stab_func, nn_clos
      !!----------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namzdf_gls in reference namelist : Vertical eddy diffivity and viscosity using gls turbulent closure scheme
      READ  ( numnam_ref, namzdf_gls, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namzdf_gls in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namzdf_gls in configuration namelist : Vertical eddy diffivity and viscosity using gls turbulent closure scheme
      READ  ( numnam_cfg, namzdf_gls, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namzdf_gls in configuration namelist' )
      IF(lwm) WRITE ( numond, namzdf_gls )

      IF(lwp) THEN                     !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_gls_init : GLS turbulent closure scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_gls : set gls mixing parameters'
         WRITE(numout,*) '      minimum value of en                           rn_emin        = ', rn_emin
         WRITE(numout,*) '      minimum value of eps                          rn_epsmin      = ', rn_epsmin
         WRITE(numout,*) '      Limit dissipation rate under stable stratif.  ln_length_lim  = ', ln_length_lim
         WRITE(numout,*) '      Galperin limit (Standard: 0.53, Holt: 0.26)   rn_clim_galp   = ', rn_clim_galp
         WRITE(numout,*) '      TKE Surface boundary condition                nn_bc_surf     = ', nn_bc_surf
         WRITE(numout,*) '      TKE Bottom boundary condition                 nn_bc_bot      = ', nn_bc_bot
         WRITE(numout,*) '      Modify psi Schmidt number (wb case)           ln_sigpsi      = ', ln_sigpsi
         WRITE(numout,*) '      Craig and Banner coefficient                  rn_crban       = ', rn_crban
         WRITE(numout,*) '      Charnock coefficient                          rn_charn       = ', rn_charn
         WRITE(numout,*) '      Surface roughness formula                     nn_z0_met      = ', nn_z0_met
         WRITE(numout,*) '      surface wave breaking under ice               nn_z0_ice      = ', nn_z0_ice
         SELECT CASE( nn_z0_ice )
         CASE( 0 )   ;   WRITE(numout,*) '   ==>>>   no impact of ice cover on surface wave breaking'
         CASE( 1 )   ;   WRITE(numout,*) '   ==>>>   roughness uses rn_hsri and is weigthed by 1-TANH( fr_i(:,:) * 10 )'
         CASE( 2 )   ;   WRITE(numout,*) '   ==>>>   roughness uses rn_hsri and is weighted by 1-fr_i(:,:)'
         CASE( 3 )   ;   WRITE(numout,*) '   ==>>>   roughness uses rn_hsri and is weighted by 1-MIN( 1, 4 * fr_i(:,:) )'
         CASE DEFAULT
            CALL ctl_stop( 'zdf_gls_init: wrong value for nn_z0_ice, should be 0,1,2, or 3')
         END SELECT
         WRITE(numout,*) '      Wave height frac. (used if nn_z0_met=2)       rn_frac_hs     = ', rn_frac_hs
         WRITE(numout,*) '      Stability functions                           nn_stab_func   = ', nn_stab_func
         WRITE(numout,*) '      Type of closure                               nn_clos        = ', nn_clos
         WRITE(numout,*) '      Surface roughness (m)                         rn_hsro        = ', rn_hsro
         WRITE(numout,*) '      Ice-ocean roughness (used if nn_z0_ice/=0)    rn_hsri        = ', rn_hsri
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namdrg_top/_bot:   used values:'
         WRITE(numout,*) '      top    ocean cavity roughness (m)             rn_z0(_top)   = ', r_z0_top
         WRITE(numout,*) '      Bottom seafloor     roughness (m)             rn_z0(_bot)   = ', r_z0_bot
         WRITE(numout,*)
      ENDIF

      !                                !* allocate GLS arrays
      IF( zdf_gls_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_gls_init : unable to allocate arrays' )

      !                                !* Check of some namelist values
      IF( nn_bc_surf < 0   .OR. nn_bc_surf   > 1 )              CALL ctl_stop( 'zdf_gls_init: bad flag: nn_bc_surf is 0 or 1' )
      IF( nn_bc_surf < 0   .OR. nn_bc_surf   > 1 )              CALL ctl_stop( 'zdf_gls_init: bad flag: nn_bc_surf is 0 or 1' )
      IF( nn_z0_met  < 0   .OR. nn_z0_met    > 3 )              CALL ctl_stop( 'zdf_gls_init: bad flag: nn_z0_met is 0, 1, 2 or 3' )
      IF( nn_z0_met == 3  .AND. .NOT. (ln_wave .AND. ln_sdw ) ) CALL ctl_stop( 'zdf_gls_init: nn_z0_met=3 requires ln_wave=T and ln_sdw=T' )
      IF( nn_stab_func < 0 .OR. nn_stab_func > 3 )              CALL ctl_stop( 'zdf_gls_init: bad flag: nn_stab_func is 0, 1, 2 and 3' )
      IF( nn_clos      < 0 .OR. nn_clos      > 3 )              CALL ctl_stop( 'zdf_gls_init: bad flag: nn_clos is 0, 1, 2 or 3' )

      SELECT CASE ( nn_clos )          !* set the parameters for the chosen closure
      !
      CASE( 0 )                              ! k-kl  (Mellor-Yamada)
         !
         IF(lwp) WRITE(numout,*) '   ==>>   k-kl closure chosen (i.e. closed to the classical Mellor-Yamada)'
         IF(lwp) WRITE(numout,*)
         rpp     = 0._wp
         rmm     = 1._wp
         rnn     = 1._wp
         rsc_tke = 1.96_wp
         rsc_psi = 1.96_wp
         rpsi1   = 0.9_wp
         rpsi3p  = 1._wp
         rpsi2   = 0.5_wp
         !
         SELECT CASE ( nn_stab_func )
         CASE( 0, 1 )   ;   rpsi3m = 2.53_wp       ! G88 or KC stability functions
         CASE( 2 )      ;   rpsi3m = 2.62_wp       ! Canuto A stability functions
         CASE( 3 )      ;   rpsi3m = 2.38          ! Canuto B stability functions (caution : constant not identified)
         END SELECT
         !
      CASE( 1 )                              ! k-eps
         !
         IF(lwp) WRITE(numout,*) '   ==>>   k-eps closure chosen'
         IF(lwp) WRITE(numout,*)
         rpp     =  3._wp
         rmm     =  1.5_wp
         rnn     = -1._wp
         rsc_tke =  1._wp
         rsc_psi =  1.2_wp  ! Schmidt number for psi
         rpsi1   =  1.44_wp
         rpsi3p  =  1._wp
         rpsi2   =  1.92_wp
         !
         SELECT CASE ( nn_stab_func )
         CASE( 0, 1 )   ;   rpsi3m = -0.52_wp      ! G88 or KC stability functions
         CASE( 2 )      ;   rpsi3m = -0.629_wp     ! Canuto A stability functions
         CASE( 3 )      ;   rpsi3m = -0.566        ! Canuto B stability functions
         END SELECT
         !
      CASE( 2 )                              ! k-omega
         !
         IF(lwp) WRITE(numout,*) '   ==>>   k-omega closure chosen'
         IF(lwp) WRITE(numout,*)
         rpp     = -1._wp
         rmm     =  0.5_wp
         rnn     = -1._wp
         rsc_tke =  2._wp
         rsc_psi =  2._wp
         rpsi1   =  0.555_wp
         rpsi3p  =  1._wp
         rpsi2   =  0.833_wp
         !
         SELECT CASE ( nn_stab_func )
         CASE( 0, 1 )   ;   rpsi3m = -0.58_wp       ! G88 or KC stability functions
         CASE( 2 )      ;   rpsi3m = -0.64_wp       ! Canuto A stability functions
         CASE( 3 )      ;   rpsi3m = -0.64_wp       ! Canuto B stability functions caution : constant not identified)
         END SELECT
         !
      CASE( 3 )                              ! generic
         !
         IF(lwp) WRITE(numout,*) '   ==>>   generic closure chosen'
         IF(lwp) WRITE(numout,*)
         rpp     = 2._wp
         rmm     = 1._wp
         rnn     = -0.67_wp
         rsc_tke = 0.8_wp
         rsc_psi = 1.07_wp
         rpsi1   = 1._wp
         rpsi3p  = 1._wp
         rpsi2   = 1.22_wp
         !
         SELECT CASE ( nn_stab_func )
         CASE( 0, 1 )   ;   rpsi3m = 0.1_wp         ! G88 or KC stability functions
         CASE( 2 )      ;   rpsi3m = 0.05_wp        ! Canuto A stability functions
         CASE( 3 )      ;   rpsi3m = 0.05_wp        ! Canuto B stability functions caution : constant not identified)
         END SELECT
         !
      END SELECT

      !
      SELECT CASE ( nn_stab_func )     !* set the parameters of the stability functions
      !
      CASE ( 0 )                             ! Galperin stability functions
         !
         IF(lwp) WRITE(numout,*) '   ==>>   Stability functions from Galperin'
         rc2     =  0._wp
         rc3     =  0._wp
         rc_diff =  1._wp
         rc0     =  0.5544_wp
         rcm_sf  =  0.9884_wp
         rghmin  = -0.28_wp
         rgh0    =  0.0233_wp
         rghcri  =  0.02_wp
         !
      CASE ( 1 )                             ! Kantha-Clayson stability functions
         !
         IF(lwp) WRITE(numout,*) '   ==>>   Stability functions from Kantha-Clayson'
         rc2     =  0.7_wp
         rc3     =  0.2_wp
         rc_diff =  1._wp
         rc0     =  0.5544_wp
         rcm_sf  =  0.9884_wp
         rghmin  = -0.28_wp
         rgh0    =  0.0233_wp
         rghcri  =  0.02_wp
         !
      CASE ( 2 )                             ! Canuto A stability functions
         !
         IF(lwp) WRITE(numout,*) '   ==>>   Stability functions from Canuto A'
         rs0 = 1.5_wp * rl1 * rl5*rl5
         rs1 = -rl4*(rl6+rl7) + 2._wp*rl4*rl5*(rl1-(1._wp/3._wp)*rl2-rl3) + 1.5_wp*rl1*rl5*rl8
         rs2 = -(3._wp/8._wp) * rl1*(rl6*rl6-rl7*rl7)
         rs4 = 2._wp * rl5
         rs5 = 2._wp * rl4
         rs6 = (2._wp/3._wp) * rl5 * ( 3._wp*rl3*rl3 - rl2*rl2 ) - 0.5_wp * rl5*rl1 * (3._wp*rl3-rl2)   &
            &                                                    + 0.75_wp * rl1 * ( rl6 - rl7 )
         rd0 = 3._wp * rl5*rl5
         rd1 = rl5 * ( 7._wp*rl4 + 3._wp*rl8 )
         rd2 = rl5*rl5 * ( 3._wp*rl3*rl3 - rl2*rl2 ) - 0.75_wp*(rl6*rl6 - rl7*rl7 )
         rd3 = rl4 * ( 4._wp*rl4 + 3._wp*rl8)
         rd4 = rl4 * ( rl2 * rl6 - 3._wp*rl3*rl7 - rl5*(rl2*rl2 - rl3*rl3 ) ) + rl5*rl8 * ( 3._wp*rl3*rl3 - rl2*rl2 )
         rd5 = 0.25_wp * ( rl2*rl2 - 3._wp *rl3*rl3 ) * ( rl6*rl6 - rl7*rl7 )
         rc0 = 0.5268_wp
         rf6 = 8._wp / (rc0**6._wp)
         rc_diff = SQRT(2._wp) / (rc0**3._wp)
         rcm_sf  =  0.7310_wp
         rghmin  = -0.28_wp
         rgh0    =  0.0329_wp
         rghcri  =  0.03_wp
         !
      CASE ( 3 )                             ! Canuto B stability functions
         !
         IF(lwp) WRITE(numout,*) '   ==>>   Stability functions from Canuto B'
         rs0 = 1.5_wp * rm1 * rm5*rm5
         rs1 = -rm4 * (rm6+rm7) + 2._wp * rm4*rm5*(rm1-(1._wp/3._wp)*rm2-rm3) + 1.5_wp * rm1*rm5*rm8
         rs2 = -(3._wp/8._wp) * rm1 * (rm6*rm6-rm7*rm7 )
         rs4 = 2._wp * rm5
         rs5 = 2._wp * rm4
         rs6 = (2._wp/3._wp) * rm5 * (3._wp*rm3*rm3-rm2*rm2) - 0.5_wp * rm5*rm1*(3._wp*rm3-rm2) + 0.75_wp * rm1*(rm6-rm7)
         rd0 = 3._wp * rm5*rm5
         rd1 = rm5 * (7._wp*rm4 + 3._wp*rm8)
         rd2 = rm5*rm5 * (3._wp*rm3*rm3 - rm2*rm2) - 0.75_wp * (rm6*rm6 - rm7*rm7)
         rd3 = rm4 * ( 4._wp*rm4 + 3._wp*rm8 )
         rd4 = rm4 * ( rm2*rm6 -3._wp*rm3*rm7 - rm5*(rm2*rm2 - rm3*rm3) ) + rm5 * rm8 * ( 3._wp*rm3*rm3 - rm2*rm2 )
         rd5 = 0.25_wp * ( rm2*rm2 - 3._wp*rm3*rm3 ) * ( rm6*rm6 - rm7*rm7 )
         rc0 = 0.5268_wp            !!       rc0 = 0.5540_wp (Warner ...) to verify !
         rf6 = 8._wp / ( rc0**6._wp )
         rc_diff = SQRT(2._wp)/(rc0**3.)
         rcm_sf  =  0.7470_wp
         rghmin  = -0.28_wp
         rgh0    =  0.0444_wp
         rghcri  =  0.0414_wp
         !
      END SELECT
    
      !                                !* Set Schmidt number for psi diffusion in the wave breaking case
      !                                     ! See Eq. (13) of Carniel et al, OM, 30, 225-239, 2009
      !                                     !  or Eq. (17) of Burchard, JPO, 31, 3133-3145, 2001
      IF( ln_sigpsi ) THEN
         ra_sf = -1.5 ! Set kinetic energy slope, then deduce rsc_psi and rl_sf 
         ! Verification: retrieve Burchard (2001) results by uncomenting the line below:
         ! Note that the results depend on the value of rn_cm_sf which is constant (=rc0) in his work
         ! ra_sf = -SQRT(2./3.*rc0**3./rn_cm_sf*rn_sc_tke)/vkarmn
         rsc_psi0 = rsc_tke/(24.*rpsi2)*(-1.+(4.*rnn + ra_sf*(1.+4.*rmm))**2./(ra_sf**2.))
      ELSE
         rsc_psi0 = rsc_psi
      ENDIF
 
      !                                !* Shear free turbulence parameters
      !
      ra_sf  = -4._wp*rnn*SQRT(rsc_tke) / ( (1._wp+4._wp*rmm)*SQRT(rsc_tke) &
               &                              - SQRT(rsc_tke + 24._wp*rsc_psi0*rpsi2 ) )

      IF ( rn_crban==0._wp ) THEN
         rl_sf = vkarmn
      ELSE
         rl_sf = rc0 * SQRT(rc0/rcm_sf) * SQRT( ( (1._wp + 4._wp*rmm + 8._wp*rmm**2_wp) * rsc_tke        &
            &                                            + 12._wp*rsc_psi0*rpsi2 - (1._wp + 4._wp*rmm)   &
            &                                                     *SQRT(rsc_tke*(rsc_tke                 &
            &                                                        + 24._wp*rsc_psi0*rpsi2)) )         &
            &                                              /(12._wp*rnn**2.)                             )
      ENDIF

      !
      IF(lwp) THEN                     !* Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Limit values :'
         WRITE(numout,*) '      Parameter  m = ', rmm
         WRITE(numout,*) '      Parameter  n = ', rnn
         WRITE(numout,*) '      Parameter  p = ', rpp
         WRITE(numout,*) '      rpsi1    = ', rpsi1
         WRITE(numout,*) '      rpsi2    = ', rpsi2
         WRITE(numout,*) '      rpsi3m   = ', rpsi3m
         WRITE(numout,*) '      rpsi3p   = ', rpsi3p
         WRITE(numout,*) '      rsc_tke  = ', rsc_tke
         WRITE(numout,*) '      rsc_psi  = ', rsc_psi
         WRITE(numout,*) '      rsc_psi0 = ', rsc_psi0
         WRITE(numout,*) '      rc0      = ', rc0
         WRITE(numout,*)
         WRITE(numout,*) '   Shear free turbulence parameters:'
         WRITE(numout,*) '      rcm_sf   = ', rcm_sf
         WRITE(numout,*) '      ra_sf    = ', ra_sf
         WRITE(numout,*) '      rl_sf    = ', rl_sf
      ENDIF

      !                                !* Constants initialization
      rc02  = rc0  * rc0   ;   rc02r = 1. / rc02
      rc03  = rc02 * rc0
      rc04  = rc03 * rc0
      rsbc_tke1 = -3._wp/2._wp*rn_crban*ra_sf*rl_sf                      ! Dirichlet + Wave breaking
      rsbc_tke2 = rdt * rn_crban / rl_sf                                 ! Neumann + Wave breaking 
      zcr = MAX(rsmall, rsbc_tke1**(1./(-ra_sf*3._wp/2._wp))-1._wp )
      rtrans = 0.2_wp / zcr                                              ! Ad. inverse transition length between log and wave layer 
      rsbc_zs1  = rn_charn/grav                                          ! Charnock formula for surface roughness
      rsbc_zs2  = rn_frac_hs / 0.85_wp / grav * 665._wp                  ! Rascle formula for surface roughness 
      rsbc_psi1 = -0.5_wp * rdt * rc0**(rpp-2._wp*rmm) / rsc_psi
      rsbc_psi2 = -0.5_wp * rdt * rc0**rpp * rnn * vkarmn**rnn / rsc_psi ! Neumann + NO Wave breaking 
      !
      rfact_tke = -0.5_wp / rsc_tke * rdt                                ! Cst used for the Diffusion term of tke
      rfact_psi = -0.5_wp / rsc_psi * rdt                                ! Cst used for the Diffusion term of tke
      !
      !                                !* Wall proximity function
!!gm tmask or wmask ????
      zwall(:,:,:) = 1._wp * tmask(:,:,:)

      !                                !* read or initialize all required files  
      CALL gls_rst( nit000, 'READ' )      ! (en, avt_k, avm_k, hmxl_n)
      !
      IF( lwxios ) THEN
         CALL iom_set_rstw_var_active('en')
         CALL iom_set_rstw_var_active('avt_k')
         CALL iom_set_rstw_var_active('avm_k')
         CALL iom_set_rstw_var_active('hmxl_n')
      ENDIF
      !
   END SUBROUTINE zdf_gls_init


   SUBROUTINE gls_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE gls_rst  ***
      !!                     
      !! ** Purpose :   Read or write TKE file (en) in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain TKE, en is either 
      !!                set to rn_emin or recomputed (nn_igls/=0)
      !!----------------------------------------------------------------------
      USE zdf_oce , ONLY : en, avt_k, avm_k   ! ocean vertical physics
      !!
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   jit, jk   ! dummy loop indices
      INTEGER ::   id1, id2, id3, id4
      INTEGER ::   ji, jj, ikbu, ikbv
      REAL(wp)::   cbx, cby
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            id1 = iom_varid( numror, 'en'    , ldstop = .FALSE. )
            id2 = iom_varid( numror, 'avt_k' , ldstop = .FALSE. )
            id3 = iom_varid( numror, 'avm_k' , ldstop = .FALSE. )
            id4 = iom_varid( numror, 'hmxl_n', ldstop = .FALSE. )
            !
            IF( MIN( id1, id2, id3, id4 ) > 0 ) THEN        ! all required arrays exist
               CALL iom_get( numror, jpdom_autoglo, 'en'    , en    , ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'avt_k' , avt_k , ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'avm_k' , avm_k , ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'hmxl_n', hmxl_n, ldxios = lrxios )
            ELSE                        
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>   previous run without GLS scheme, set en and hmxl_n to background values'
               en    (:,:,:) = rn_emin
               hmxl_n(:,:,:) = 0.05_wp
               ! avt_k, avm_k already set to the background value in zdf_phy_init
            ENDIF
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest, set en and hmxl_n by background values'
            en    (:,:,:) = rn_emin
            hmxl_n(:,:,:) = 0.05_wp
            ! avt_k, avm_k already set to the background value in zdf_phy_init
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- gls-rst ----'
         IF( lwxios ) CALL iom_swap(      cwxios_context         )
         CALL iom_rstput( kt, nitrst, numrow, 'en'    , en    , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'avt_k' , avt_k , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'avm_k' , avm_k , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'hmxl_n', hmxl_n, ldxios = lwxios )
         IF( lwxios ) CALL iom_swap(      cxios_context          )
         !
      ENDIF
      !
   END SUBROUTINE gls_rst

   !!======================================================================
END MODULE zdfgls

