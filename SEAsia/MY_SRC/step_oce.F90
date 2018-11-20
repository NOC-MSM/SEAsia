MODULE step_oce
   !!======================================================================
   !!                       ***  MODULE step_oce  ***
   !! Ocean time-stepping : module used in both initialisation phase and time stepping
   !!======================================================================
   !! History :   3.3  !  2010-08  (C. Ethe)  Original code - reorganisation of the initial phase
   !!             3.7  !  2014-01  (G. Madec) LDF simplication 
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean space and time domain variables
   USE zdf_oce          ! ocean vertical physics variables

   USE daymod           ! calendar                         (day     routine)

   USE sbc_oce          ! surface boundary condition: ocean
   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE sbcrnf           ! surface boundary condition: runoff variables
   USE sbccpl           ! surface boundary condition: coupled formulation (call send at end of step)
   USE sbcapr           ! surface boundary condition: atmospheric pressure
   USE sbctide          ! Tide initialisation
   USE sbcwave          ! Wave intialisation

   USE traqsr           ! solar radiation penetration      (tra_qsr routine)
   USE trasbc           ! surface boundary condition       (tra_sbc routine)
   USE trabbc           ! bottom boundary condition        (tra_bbc routine)
   USE trabbl           ! bottom boundary layer            (tra_bbl routine)
   USE tradmp           ! internal damping                 (tra_dmp routine)
   USE traadv           ! advection scheme control     (tra_adv_ctl routine)
   USE traldf           ! lateral mixing                   (tra_ldf routine)
   USE trazdf           ! vertical mixing                  (tra_zdf routine)
   USE tranxt           ! time-stepping                    (tra_nxt routine)
   USE tranpc           ! non-penetrative convection       (tra_npc routine)

   USE eosbn2           ! equation of state                (eos_bn2 routine)

   USE divhor           ! horizontal divergence            (div_hor routine)
   USE dynadv           ! advection                        (dyn_adv routine)
   USE dynbfr           ! Bottom friction terms            (dyn_bfr routine)
   USE dynvor           ! vorticity term                   (dyn_vor routine)
   USE dynhpg           ! hydrostatic pressure grad.       (dyn_hpg routine)
   USE dynldf           ! lateral momentum diffusion       (dyn_ldf routine)
   USE dynzdf           ! vertical diffusion               (dyn_zdf routine)
   USE dynspg           ! surface pressure gradient        (dyn_spg routine)

   USE dynnxt           ! time-stepping                    (dyn_nxt routine)

   USE stopar           ! Stochastic parametrization       (sto_par routine)
   USE stopts 

   USE bdy_oce    , ONLY: ln_bdy
   USE bdydta           ! open boundary condition data     (bdy_dta routine)
   USE bdytra           ! bdy cond. for tracers            (bdy_tra routine)
   USE bdydyn3d         ! bdy cond. for baroclinic vel.  (bdy_dyn3d routine)

   USE sshwzv           ! vertical velocity and ssh        (ssh_nxt routine)
   !                                                       (ssh_swp routine)
   !                                                       (wzv     routine)
   USE domvvl           ! variable vertical scale factors  (dom_vvl_sf_nxt routine)
   !                                                       (dom_vvl_sf_swp routine)

   USE ldfslp           ! iso-neutral slopes               (ldf_slp routine)
   USE ldfdyn           ! lateral eddy viscosity coef.     (ldf_dyn routine)
   USE ldftra           ! lateral eddy diffusive coef.     (ldf_tra routine)

   USE zdftmx           ! tide-induced vertical mixing     (zdf_tmx routine)
   USE zdfbfr           ! bottom friction                  (zdf_bfr routine)
   USE zdftke           ! TKE vertical mixing              (zdf_tke routine)
   USE zdfgls           ! GLS vertical mixing              (zdf_gls routine)
   USE zdfddm           ! double diffusion mixing          (zdf_ddm routine)
   USE zdfevd           ! enhanced vertical diffusion      (zdf_evd routine)
   USE zdfric           ! Richardson vertical mixing       (zdf_ric routine)
   USE zdfmxl           ! Mixed-layer depth                (zdf_mxl routine)
   USE zdfqiao          !Qiao module wave induced mixing   (zdf_qiao routine)

   USE step_diu        ! Time stepping for diurnal sst
   USE diurnal_bulk    ! diurnal SST bulk routines  (diurnal_sst_takaya routine) 
   USE cool_skin       ! diurnal cool skin correction (diurnal_sst_coolskin routine)   
   USE sbc_oce         ! surface fluxes  
   
   USE zpshde           ! partial step: hor. derivative     (zps_hde routine)

   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE diaptr           ! poleward transports              (dia_ptr routine)
   USE diadct           ! sections transports              (dia_dct routine)
   USE diaar5           ! AR5 diagnosics                   (dia_ar5 routine)
   USE diahth           ! thermocline depth                (dia_hth routine)
   USE diahsb           ! heat, salt and volume budgets    (dia_hsb routine)
   USE diaharm
!--- NB for restart hamonic analysis
   USE diaharm_fast     ! harmonic analysis of tides (harm_ana routine) 
!--- END NB -----------------------------------
   USE diacfl
   USE flo_oce          ! floats variables
   USE floats           ! floats computation               (flo_stp routine)

   USE crsfld           ! Standard output on coarse grid   (crs_fld routine)

   USE asminc           ! assimilation increments      (tra_asm_inc routine)
   !                                                   (dyn_asm_inc routine)
   USE asmbkg
   USE stpctl           ! time stepping control            (stp_ctl routine)
   USE restart          ! ocean restart                    (rst_wri routine)
   USE prtctl           ! Print control                    (prt_ctl routine)

   USE diaobs           ! Observation operator

   USE in_out_manager   ! I/O manager
   USE iom              !
   USE lbclnk
   USE timing           ! Timing

#if defined key_iomput
   USE xios
#endif
#if defined key_agrif
   USE agrif_opa_sponge ! Momemtum and tracers sponges
   USE agrif_opa_update ! Update (2-way nesting)
#endif
#if defined key_top
   USE trcstp           ! passive tracer time-stepping      (trc_stp routine)
#endif
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: step_oce.F90 7646 2017-02-06 09:25:03Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE step_oce
