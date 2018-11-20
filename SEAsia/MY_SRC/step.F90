MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             -   !  1991-11  (G. Madec)
   !!             -   !  1992-06  (M. Imbard)  add a first output record
   !!             -   !  1996-04  (G. Madec)  introduction of dynspg
   !!             -   !  1996-04  (M.A. Foujols)  introduction of passive tracer
   !!            8.0  !  1997-06  (G. Madec)  new architecture of call
   !!            8.2  !  1997-06  (G. Madec, M. Imbard, G. Roullet)  free surface
   !!             -   !  1999-02  (G. Madec, N. Grima)  hpg implicit
   !!             -   !  2000-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
   !!   NEMO     1.0  !  2002-06  (G. Madec)  free form, suppress macro-tasking
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-01  (C. Ethe) Add the KPP closure scheme
   !!             -   !  2005-11  (G. Madec)  Reorganisation of tra and dyn calls
   !!             -   !  2006-01  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   !  2006-07  (S. Masson)  restart using iom
   !!            3.2  !  2009-02  (G. Madec, R. Benshila)  reintroduicing z*-coordinate
   !!             -   !  2009-06  (S. Masson, G. Madec)  TKE restart compatible with key_cpl
   !!            3.3  !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal
   !!            3.6  !  2012-07  (J. Simeon, G. Madec. C. Ethe)  Online coarsening of outputs
   !!            3.6  !  2014-04  (F. Roquet, G. Madec) New equations of state
   !!            3.6  !  2014-10  (E. Clementi, P. Oddo) Add Qiao vertical mixing in case of waves
   !!            3.7  !  2014-10  (G. Madec)  LDF simplication 
   !!             -   !  2014-12  (G. Madec) remove KPP scheme
   !!             -   !  2015-11  (J. Chanut) free surface simplification
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp             : OPA system time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules
   !
   USE iom              ! xIOs server

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nemogcm.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: step.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
      !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
      !!              - Tme stepping  of TRC (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, hdiv,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji,jj,jk ! dummy loop indice
      INTEGER ::   indic    ! error indicator if < 0
      INTEGER ::   kcall    ! optional integer argument (dom_vvl_sf_nxt)
      !! ---------------------------------------------------------------------
#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar 
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             indic = 0                ! reset to no error condition
                             
      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init(      cxios_context          )  ! for model grid (including passible AGRIF zoom)
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update external forcing (tides, open boundaries, and surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_tide    )   CALL sbc_tide( kstp )                   ! update tide potential
      IF( ln_apr_dyn )   CALL sbc_apr ( kstp )                   ! atmospheric pressure (NB: call before bdy_dta which needs ssh_ib) 
      IF( ln_bdy     )   CALL bdy_dta ( kstp, time_offset=+1 )   ! update dynamic & tracer data at open boundaries
                         CALL sbc     ( kstp )                   ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update stochastic parameters and random T/S fluctuations
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_sto_eos ) CALL sto_par( kstp )          ! Stochastic parameters
      IF( ln_sto_eos ) CALL sto_pts( tsn  )          ! Random T/S fluctuations

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update                (ua, va, tsa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  THERMODYNAMICS
                         CALL eos_rab( tsb, rab_b )       ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( tsn, rab_n )       ! now    local thermal/haline expension ratio at T-points
                         CALL bn2    ( tsb, rab_b, rn2b ) ! before Brunt-Vaisala frequency
                         CALL bn2    ( tsn, rab_n, rn2  ) ! now    Brunt-Vaisala frequency

      !
      !  VERTICAL PHYSICS
                         CALL zdf_bfr( kstp )         ! bottom friction (if quadratic)
      !                                               ! Vertical eddy viscosity and diffusivity coefficients
      IF( lk_zdfric  )   CALL zdf_ric ( kstp )             ! Richardson number dependent Kz
      IF( lk_zdftke  )   CALL zdf_tke ( kstp )             ! TKE closure scheme for Kz
      IF( lk_zdfgls  )   CALL zdf_gls ( kstp )             ! GLS closure scheme for Kz
      IF( ln_zdfqiao )   CALL zdf_qiao( kstp )             ! Qiao vertical mixing 
      !
      IF( lk_zdfcst  ) THEN                                ! Constant Kz (reset avt, avm[uv] to the background value)
         avt (:,:,:) = rn_avt0 * wmask (:,:,:)
         avmu(:,:,:) = rn_avm0 * wumask(:,:,:)
         avmv(:,:,:) = rn_avm0 * wvmask(:,:,:)
      ENDIF

      IF( ln_rnf_mouth ) THEN                         ! increase diffusivity at rivers mouths
         DO jk = 2, nkrnf   ;   avt(:,:,jk) = avt(:,:,jk) + 2._wp * rn_avt_rnf * rnfmsk(:,:) * tmask(:,:,jk)   ;   END DO
      ENDIF
      IF( ln_zdfevd  )   CALL zdf_evd( kstp )         ! enhanced vertical eddy diffusivity

      IF( lk_zdftmx  )   CALL zdf_tmx( kstp )         ! tidal vertical mixing

      IF( lk_zdfddm  )   CALL zdf_ddm( kstp )         ! double diffusive mixing

                         CALL zdf_mxl( kstp )         ! mixed layer depth

                                                      ! write TKE or GLS information in the restart file
      IF( lrst_oce .AND. lk_zdftke )   CALL tke_rst( kstp, 'WRITE' )
      IF( lrst_oce .AND. lk_zdfgls )   CALL gls_rst( kstp, 'WRITE' )
      !
      !  LATERAL  PHYSICS
      !
      IF( l_ldfslp ) THEN                             ! slope of lateral mixing
                         CALL eos( tsb, rhd, gdept_0(:,:,:) )               ! before in situ density

         IF( ln_zps .AND. .NOT. ln_isfcav)                               &
            &            CALL zps_hde    ( kstp, jpts, tsb, gtsu, gtsv,  &  ! Partial steps: before horizontal gradient
            &                                          rhd, gru , grv    )  ! of t, s, rd at the last ocean level

         IF( ln_zps .AND.       ln_isfcav)                               &
            &            CALL zps_hde_isf( kstp, jpts, tsb, gtsu, gtsv, gtui, gtvi,  &  ! Partial steps for top cell (ISF)
            &                                          rhd, gru , grv , grui, grvi   )  ! of t, s, rd at the first ocean level
         IF( ln_traldf_triad ) THEN 
                         CALL ldf_slp_triad( kstp )                       ! before slope for triad operator
         ELSE     
                         CALL ldf_slp     ( kstp, rhd, rn2b )             ! before slope for standard operator
         ENDIF
      ENDIF
      !                                                                   ! eddy diffusivity coeff.
      IF( l_ldftra_time .OR. l_ldfeiv_time )   CALL ldf_tra( kstp )       !       and/or eiv coeff.
      IF( l_ldfdyn_time                    )   CALL ldf_dyn( kstp )       ! eddy viscosity coeff. 

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            CALL ssh_nxt       ( kstp )  ! after ssh (includes call to div_hor)
      IF(.NOT.ln_linssh )   CALL dom_vvl_sf_nxt( kstp )  ! after vertical scale factors 
                            CALL wzv           ( kstp )  ! now cross-level velocity 
                            CALL eos    ( tsn, rhd, rhop, gdept_n(:,:,:) )  ! now in situ density for hpg computation
                            
!!jc: fs simplification
!!jc: lines below are useless if ln_linssh=F. Keep them here (which maintains a bug if ln_linssh=T and ln_zps=T, cf ticket #1636) 
!!                                         but ensures reproductible results
!!                                         with previous versions using split-explicit free surface          
            IF( ln_zps .AND. .NOT. ln_isfcav )                               &
               &            CALL zps_hde    ( kstp, jpts, tsn, gtsu, gtsv,   &  ! Partial steps: before horizontal gradient
               &                                          rhd, gru , grv     )  ! of t, s, rd at the last ocean level
            IF( ln_zps .AND.       ln_isfcav )                                          &
               &            CALL zps_hde_isf( kstp, jpts, tsn, gtsu, gtsv, gtui, gtvi,  &  ! Partial steps for top cell (ISF)
               &                                          rhd, gru , grv , grui, grvi   )  ! of t, s, rd at the first ocean level
!!jc: fs simplification
                            
                         ua(:,:,:) = 0._wp            ! set dynamics trends to zero
                         va(:,:,:) = 0._wp

      IF(  lk_asminc .AND. ln_asmiau .AND. ln_dyninc )   &
               &         CALL dyn_asm_inc   ( kstp )  ! apply dynamics assimilation increment
      IF( ln_bdy     )   CALL bdy_dyn3d_dmp ( kstp )  ! bdy damping trends
#if defined key_agrif
      IF(.NOT. Agrif_Root())  & 
               &         CALL Agrif_Sponge_dyn        ! momentum sponge
#endif
                         CALL dyn_adv       ( kstp )  ! advection (vector or flux form)
                         CALL dyn_vor       ( kstp )  ! vorticity term including Coriolis
                         CALL dyn_ldf       ( kstp )  ! lateral mixing
                         CALL dyn_hpg       ( kstp )  ! horizontal gradient of Hydrostatic pressure
                         CALL dyn_spg       ( kstp )  ! surface pressure gradient

                                                      ! With split-explicit free surface, since now transports have been updated and ssha as well
      IF( ln_dynspg_ts ) THEN                         ! vertical scale factors and vertical velocity need to be updated
                            CALL div_hor    ( kstp )              ! Horizontal divergence  (2nd call in time-split case)
         IF(.NOT.ln_linssh) CALL dom_vvl_sf_nxt( kstp, kcall=2 )  ! after vertical scale factors (update depth average component)
                            CALL wzv        ( kstp )              ! now cross-level velocity 
      ENDIF

                         CALL dyn_bfr       ( kstp )  ! bottom friction
                         CALL dyn_zdf       ( kstp )  ! vertical diffusion

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! cool skin
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      IF ( ln_diurnal )  CALL stp_diurnal( kstp )
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs             (ua, va, tsa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_floats  )   CALL flo_stp( kstp )         ! drifting Floats
      IF( nn_diacfl == 1 )   CALL dia_cfl( kstp )         ! Courant number diagnostics
      IF( lk_diahth  )   CALL dia_hth( kstp )         ! Thermocline depth (20 degres isotherm depth)
      IF( lk_diadct  )   CALL dia_dct( kstp )         ! Transports
                         CALL dia_ar5( kstp )         ! ar5 diag
      IF( lk_diaharm )   CALL dia_harm( kstp )        ! Tidal harmonic analysis
      ! NB - new harmonic analysis 
      IF( lk_diaharm_fast )                           &
            &            CALL dia_harm_fast( kstp )   ! Tidal harmonic analysis - restart and faster version
      ! END NB
                         CALL dia_wri( kstp )         ! ocean model: outputs
      !
      IF( ln_crs     )   CALL crs_fld       ( kstp )  ! ocean model: online field coarsening & output
      
#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL trc_stp       ( kstp )  ! time-stepping
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers                              
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         tsa(:,:,:,:) = 0._wp         ! set tracer trends to zero

      IF(  lk_asminc .AND. ln_asmiau .AND. &
         & ln_trainc )   CALL tra_asm_inc   ( kstp )  ! apply tracer assimilation increment
                         CALL tra_sbc       ( kstp )  ! surface boundary condition
      IF( ln_traqsr  )   CALL tra_qsr       ( kstp )  ! penetrative solar radiation qsr
      IF( ln_trabbc  )   CALL tra_bbc       ( kstp )  ! bottom heat flux
      IF( lk_trabbl  )   CALL tra_bbl       ( kstp )  ! advective (and/or diffusive) bottom boundary layer scheme
      IF( ln_tradmp  )   CALL tra_dmp       ( kstp )  ! internal damping trends
      IF( ln_bdy     )   CALL bdy_tra_dmp   ( kstp )  ! bdy damping trends
#if defined key_agrif
      IF(.NOT. Agrif_Root())  & 
               &         CALL Agrif_Sponge_tra        ! tracers sponge
#endif
                         CALL tra_adv       ( kstp )  ! horizontal & vertical advection
                         CALL tra_ldf       ( kstp )  ! lateral mixing

!!gm : why CALL to dia_ptr has been moved here??? (use trends info?)
      IF( ln_diaptr  )   CALL dia_ptr                 ! Poleward adv/ldf TRansports diagnostics
!!gm
                         CALL tra_zdf       ( kstp )  ! vertical mixing and after tracer fields
      IF( ln_zdfnpc  )   CALL tra_npc       ( kstp )  ! update after fields by non-penetrative convection

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Set boundary conditions and Swap
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!jc1: For agrif, it would be much better to finalize tracers/momentum here (e.g. bdy conditions) and move the swap 
!!    (and time filtering) after Agrif update. Then restart would be done after and would contain updated fields. 
!!    If so: 
!!    (i) no need to call agrif update at initialization time
!!    (ii) no need to update "before" fields 
!!
!!    Apart from creating new tra_swp/dyn_swp routines, this however: 
!!    (i) makes boundary conditions at initialization time computed from updated fields which is not the case between 
!!    two restarts => restartability issue. One can circumvent this, maybe, by assuming "interface separation", 
!!    e.g. a shift of the feedback interface inside child domain. 
!!    (ii) requires that all restart outputs of updated variables by agrif (e.g. passive tracers/tke/barotropic arrays) are done at the same
!!    place.
!! 
!!jc2: dynnxt must be the latest call. e3t_b are indeed updated in that routine
                         CALL tra_nxt       ( kstp )  ! finalize (bcs) tracer fields at next time step and swap
                         CALL dyn_nxt       ( kstp )  ! finalize (bcs) velocities at next time step and swap
                         CALL ssh_swp       ( kstp )  ! swap of sea surface height
      IF(.NOT.ln_linssh) CALL dom_vvl_sf_swp( kstp )  ! swap of vertical scale factors
      !
      IF( ln_diahsb        )   CALL dia_hsb( kstp )         ! - ML - global conservation diagnostics

!!gm : This does not only concern the dynamics ==>>> add a new title
!!gm2: why ouput restart before AGRIF update?
!!
!!jc: That would be better, but see comment above
!!
      IF( lrst_oce         )   CALL rst_write    ( kstp )   ! write output ocean restart file
      IF( ln_sto_eos       )   CALL sto_rst_write( kstp )   ! write restart file for stochastic parameters

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
                         CALL Agrif_Integrate_ChildGrids( stp )  

      IF( Agrif_NbStepint() == 0 ) THEN               ! AGRIF Update 
!!jc in fact update is useless at last time step, but do it for global diagnostics
                         CALL Agrif_Update_Tra()      ! Update active tracers
                         CALL Agrif_Update_Dyn()      ! Update momentum
      ENDIF
#endif
      IF( ln_diaobs  )   CALL dia_obs( kstp )         ! obs-minus-model (assimilation) diagnostics (call after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL stp_ctl       ( kstp, indic )
      IF( indic < 0  ) THEN
                         CALL ctl_stop( 'step: indic < 0' )
                         CALL dia_wri_state( 'output.abort', kstp )
      ENDIF
!#if defined key_harm_ana
!--- NB Restart for the tidal harmonic analysis
!      IF( ln_harm_ana_store   )   CALL harm_ana( kstp )        ! Harmonic analysis of tides 
!--- END NB -----------------------------------
!# endif
      IF( kstp == nit000 ) THEN
                 CALL iom_close( numror )     ! close input  ocean restart file
         IF(lwm) CALL FLUSH    ( numond )     ! flush output namelist oce
         IF(lwm.AND.numoni /= -1 )   &
            &    CALL FLUSH    ( numoni )     ! flush output namelist ice (if exist)
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!gm why lk_oasis and not lk_cpl ????
      IF( lk_oasis   )   CALL sbc_cpl_snd( kstp )     ! coupled mode : field exchanges
      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 ) THEN 
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
         IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) ! 
      ENDIF
#endif
      !
      IF( nn_timing == 1 .AND.  kstp == nit000  )   CALL timing_reset
      !
   END SUBROUTINE stp
   
END MODULE step
