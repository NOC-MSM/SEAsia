MODULE ice
   !!======================================================================
   !!                        ***  MODULE ice  ***
   !! LIM-3 Sea Ice physics:  diagnostics variables of ice defined in memory
   !!=====================================================================
   !! History :  3.0  ! 2008-03  (M. Vancoppenolle) original code LIM-3
   !!            4.0  ! 2011-02  (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC    ice_alloc  !  Called in sbc_lim_init

   !!======================================================================
   !! LIM3 by the use of sweat, agile fingers and sometimes brain juice, 
   !!  was developed in Louvain-la-Neuve by : 
   !!    * Martin Vancoppenolle (UCL-ASTR, Belgium)
   !!    * Sylvain Bouillon (UCL-ASTR, Belgium)
   !!    * Miguel Angel Morales Maqueda (NOC-L, UK)
   !! 
   !! Based on extremely valuable earlier work by
   !!    * Thierry Fichefet
   !!    * Hugues Goosse
   !!
   !! The following persons also contributed to the code in various ways
   !!    * Gurvan Madec, Claude Talandier, Christian Ethe (LOCEAN, France)
   !!    * Xavier Fettweis (UCL-ASTR), Ralph Timmermann (AWI, Germany)
   !!    * Bill Lipscomb (LANL), Cecilia Bitz (UWa) 
   !!      and Elisabeth Hunke (LANL), USA.
   !! 
   !! For more info, the interested user is kindly invited to consult the following references
   !!    For model description and validation :
   !!    * Vancoppenolle et al., Ocean Modelling, 2008a.
   !!    * Vancoppenolle et al., Ocean Modelling, 2008b.
   !!    For a specific description of EVP :
   !!    * Bouillon et al., Ocean Modelling 2009.
   !!
   !!    Or the reference manual, that should be available by 2011
   !!======================================================================
   !!                                                                     |
   !!              I C E   S T A T E   V A R I A B L E S                  |
   !!                                                                     |
   !! Introduction :                                                      |
   !! --------------                                                      |
   !! Every ice-covered grid cell is characterized by a series of state   |
   !! variables. To account for unresolved spatial variability in ice     |
   !! thickness, the ice cover in divided in ice thickness categories.    |
   !!                                                                     |
   !! Sea ice state variables depend on the ice thickness category        |
   !!                                                                     |
   !! Those variables are divided into two groups                         |
   !! * Extensive (or global) variables.                                  |
   !!   These are the variables that are transported by all means         |
   !! * Intensive (or equivalent) variables.                              |
   !!   These are the variables that are either physically more           |
   !!   meaningful and/or used in ice thermodynamics                      |
   !!                                                                     |
   !! Routines in limvar.F90 perform conversions                          |
   !!  - lim_var_glo2eqv  : from global to equivalent variables           |
   !!  - lim_var_eqv2glo  : from equivalent to global variables           |
   !!                                                                     |
   !! For various purposes, the sea ice state variables have sometimes    |
   !! to be aggregated over all ice thickness categories. This operation  |
   !! is done in :                                                        |
   !!  - lim_var_agg                                                      |
   !!                                                                     |
   !! in icestp.F90, the routines that compute the changes in the ice     |
   !! state variables are called                                          |
   !! - lim_dyn : ice dynamics                                            |
   !! - lim_trp : ice transport                                           |
   !! - lim_itd_me : mechanical redistribution (ridging and rafting)      |
   !! - lim_thd : ice halo-thermodynamics                                 |
   !! - lim_itd_th : thermodynamic changes in ice thickness distribution  |
   !!                and creation of new ice                              |
   !!                                                                     |
   !! See the associated routines for more information                    |
   !!                                                                     |
   !! List of ice state variables :                                       |
   !! -----------------------------                                       |
   !!                                                                     |
   !!-------------|-------------|---------------------------------|-------|
   !!   name in   |   name in   |              meaning            | units |
   !! 2D routines | 1D routines |                                 |       |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Dynamical variables (prognostic)                    *** |
   !! ******************************************************************* |
   !!                                                                     |
   !! u_ice       |      -      |    Comp. U of the ice velocity  | m/s   |
   !! v_ice       |      -      |    Comp. V of the ice velocity  | m/s   |
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Category dependent state variables (prognostic)     *** |
   !! ******************************************************************* |
   !!                                                                     |
   !! ** Global variables                                                 |
   !!-------------|-------------|---------------------------------|-------|
   !! a_i         | a_i_1d      |    Ice concentration            |       |
   !! v_i         |      -      |    Ice volume per unit area     | m     |
   !! v_s         |      -      |    Snow volume per unit area    | m     |
   !! smv_i       |      -      |    Sea ice salt content         | ppt.m |
   !! oa_i        !      -      !    Sea ice areal age content    | day   |
   !! e_i         !      -      !    Ice enthalpy                 | J/m2  | 
   !!      -      ! q_i_1d      !    Ice enthalpy per unit vol.   | J/m3  | 
   !! e_s         !      -      !    Snow enthalpy                | J/m2  | 
   !!      -      ! q_s_1d      !    Snow enthalpy per unit vol.  | J/m3  | 
   !!                                                                     |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! ** Equivalent variables                                             |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! ht_i        | ht_i_1d     |    Ice thickness                | m     |
   !! ht_s        ! ht_s_1d     |    Snow depth                   | m     |
   !! sm_i        ! sm_i_1d     |    Sea ice bulk salinity        ! ppt   |
   !! s_i         ! s_i_1d      |    Sea ice salinity profile     ! ppt   |
   !! o_i         !      -      |    Sea ice Age                  ! days  |
   !! t_i         ! t_i_1d      |    Sea ice temperature          ! K     |
   !! t_s         ! t_s_1d      |    Snow temperature             ! K     |
   !! t_su        ! t_su_1d     |    Sea ice surface temperature  ! K     |
   !!                                                                     |
   !! notes: the ice model only sees a bulk (i.e., vertically averaged)   |
   !!        salinity, except in thermodynamic computations, for which    |
   !!        the salinity profile is computed as a function of bulk       |
   !!        salinity                                                     |
   !!                                                                     |
   !!        the sea ice surface temperature is not associated to any     |
   !!        heat content. Therefore, it is not a state variable and      |
   !!        does not have to be advected. Nevertheless, it has to be     |
   !!        computed to determine whether the ice is melting or not      |
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Category-summed state variables (diagnostic)        *** |
   !! ******************************************************************* |
   !! at_i        | at_i_1d     |    Total ice concentration      |       |
   !! vt_i        |      -      |    Total ice vol. per unit area | m     |
   !! vt_s        |      -      |    Total snow vol. per unit ar. | m     |
   !! smt_i       |      -      |    Mean sea ice salinity        | ppt   |
   !! tm_i        |      -      |    Mean sea ice temperature     | K     |
   !! et_i        !      -      !    Total ice enthalpy           | J/m2  | 
   !! et_s        !      -      !    Total snow enthalpy          | J/m2  | 
   !! bv_i        !      -      !    relative brine volume        | ???   | 
   !!=====================================================================

   LOGICAL, PUBLIC ::   con_i = .false.   ! switch for conservation test

   !!--------------------------------------------------------------------------
   !! * Share Module variables
   !!--------------------------------------------------------------------------
   !                                     !!** ice-generic parameters namelist (namicerun) **
   INTEGER           , PUBLIC ::   jpl             !: number of ice  categories 
   INTEGER           , PUBLIC ::   nlay_i          !: number of ice  layers 
   INTEGER           , PUBLIC ::   nlay_s          !: number of snow layers 
   REAL(wp)          , PUBLIC ::   rn_amax_n       !: maximum ice concentration Northern hemisphere
   REAL(wp)          , PUBLIC ::   rn_amax_s       !: maximum ice concentration Southern hemisphere
   ! DRM, 04/07/17 - increase length of cn_icerst_in, cn_icerst_out strings to prevent overruns.
   CHARACTER(len=256) , PUBLIC ::   cn_icerst_in    !: suffix of ice restart name (input)
   CHARACTER(len=256) , PUBLIC ::   cn_icerst_out   !: suffix of ice restart name (output)
   CHARACTER(len=256), PUBLIC ::   cn_icerst_indir !: ice restart input directory
   CHARACTER(len=256), PUBLIC ::   cn_icerst_outdir!: ice restart output directory
   LOGICAL           , PUBLIC ::   ln_limthd       !: flag for ice thermo (T) or not (F)
   LOGICAL           , PUBLIC ::   ln_limdyn       !: flag for ice dynamics (T) or not (F)
   INTEGER           , PUBLIC ::   nn_limdyn       !: flag for ice dynamics
   REAL(wp)          , PUBLIC ::   rn_uice         !: prescribed u-vel (case nn_limdyn=0)
   REAL(wp)          , PUBLIC ::   rn_vice         !: prescribed v-vel (case nn_limdyn=0)
   
   !                                     !!** ice-diagnostics namelist (namicediag) **
   LOGICAL , PUBLIC ::   ln_limdiachk     !: flag for ice diag (T) or not (F)
   LOGICAL , PUBLIC ::   ln_limdiahsb     !: flag for ice diag (T) or not (F)
   LOGICAL , PUBLIC ::   ln_limctl        !: flag for sea-ice points output (T) or not (F)
   INTEGER , PUBLIC ::   iiceprt          !: debug i-point
   INTEGER , PUBLIC ::   jiceprt          !: debug j-point

   !                                     !!** ice-init namelist (namiceini) **
                                          ! -- limistate -- !
   LOGICAL , PUBLIC ::   ln_limini        ! initialization or not
   LOGICAL , PUBLIC ::   ln_limini_file   ! Ice initialization state from 2D netcdf file
   REAL(wp), PUBLIC ::   rn_thres_sst     ! threshold water temperature for initial sea ice
   REAL(wp), PUBLIC ::   rn_hts_ini_n     ! initial snow thickness in the north
   REAL(wp), PUBLIC ::   rn_hts_ini_s     ! initial snow thickness in the south
   REAL(wp), PUBLIC ::   rn_hti_ini_n     ! initial ice thickness in the north
   REAL(wp), PUBLIC ::   rn_hti_ini_s     ! initial ice thickness in the south
   REAL(wp), PUBLIC ::   rn_ati_ini_n     ! initial leads area in the north
   REAL(wp), PUBLIC ::   rn_ati_ini_s     ! initial leads area in the south
   REAL(wp), PUBLIC ::   rn_smi_ini_n     ! initial salinity 
   REAL(wp), PUBLIC ::   rn_smi_ini_s     ! initial salinity
   REAL(wp), PUBLIC ::   rn_tmi_ini_n     ! initial temperature
   REAL(wp), PUBLIC ::   rn_tmi_ini_s     ! initial temperature
   
   !                                     !!** ice-thickness distribution namelist (namiceitd) **
   INTEGER , PUBLIC ::   nn_catbnd        !: categories distribution following: tanh function (1), or h^(-alpha) function (2)
   REAL(wp), PUBLIC ::   rn_himean        !: mean thickness of the domain (used to compute the distribution, nn_itdshp = 2 only)

   !                                     !!** ice-dynamics namelist (namicedyn) **
                                          ! -- limtrp & limadv -- !
   INTEGER , PUBLIC ::   nn_limadv        !: choose the advection scheme (-1=Prather ; 0=Ultimate-Macho)
   INTEGER , PUBLIC ::   nn_limadv_ord    !: choose the order of the advection scheme (if Ultimate-Macho)   
                                          ! -- limitd_me -- !
   INTEGER , PUBLIC ::   nn_icestr        !: ice strength parameterization (0=Hibler79 1=Rothrock75)
   REAL(wp), PUBLIC ::   rn_pe_rdg        !: ridging work divided by pot. energy change in ridging, nn_icestr = 1
   REAL(wp), PUBLIC ::   rn_pstar         !: determines ice strength, Hibler JPO79
   REAL(wp), PUBLIC ::   rn_crhg          !: determines changes in ice strength
   LOGICAL , PUBLIC ::   ln_icestr_bvf    !: use brine volume to diminish ice strength
                                          ! -- limdyn & limrhg -- !
   REAL(wp), PUBLIC ::   rn_cio           !: drag coefficient for oceanic stress
   REAL(wp), PUBLIC ::   rn_creepl        !: creep limit : has to be under 1.0e-9
   REAL(wp), PUBLIC ::   rn_ecc           !: eccentricity of the elliptical yield curve
   INTEGER , PUBLIC ::   nn_nevp          !: number of iterations for subcycling
   REAL(wp), PUBLIC ::   rn_relast        !: ratio => telast/rdt_ice (1/3 or 1/9 depending on nb of subcycling nevp) 
   LOGICAL , PUBLIC ::   ln_landfast      !: landfast ice parameterization (T or F) 
   REAL(wp), PUBLIC ::   rn_gamma         !: fraction of ocean depth that ice must reach to initiate landfast ice
   REAL(wp), PUBLIC ::   rn_icebfr        !: maximum bottom stress per unit area of contact (landfast ice) 
   REAL(wp), PUBLIC ::   rn_lfrelax       !: relaxation time scale (s-1) to reach static friction (landfast ice) 

   !                                     !!** ice-diffusion namelist (namicehdf) **
   INTEGER , PUBLIC ::   nn_ahi0          !: sea-ice hor. eddy diffusivity coeff. (3 ways of calculation)
   REAL(wp), PUBLIC ::   rn_ahi0_ref      !: sea-ice hor. eddy diffusivity coeff. (m2/s)

   !                                     !!** ice-thermodynamics namelist (namicethd) **
                                          ! -- limthd_dif -- !
   REAL(wp), PUBLIC ::   rn_kappa_i       !: coef. for the extinction of radiation Grenfell et al. (2006) [1/m]
   REAL(wp), PUBLIC ::   nn_conv_dif      !: maximal number of iterations for heat diffusion
   REAL(wp), PUBLIC ::   rn_terr_dif      !: maximal tolerated error (C) for heat diffusion
   INTEGER , PUBLIC ::   nn_ice_thcon     !: thermal conductivity: =0 Untersteiner (1964) ; =1 Pringle et al (2007)
   LOGICAL , PUBLIC ::   ln_it_qnsice     !: iterate surface flux with changing surface temperature or not (F)
   INTEGER , PUBLIC ::   nn_monocat       !: virtual ITD mono-category parameterizations (1) or not (0)
   REAL(wp), PUBLIC ::   rn_cdsn          !: thermal conductivity of the snow [W/m/K]
                                          ! -- limthd_dh -- !
   LOGICAL , PUBLIC ::   ln_limdH         !: activate ice thickness change from growing/melting (T) or not (F)
   REAL(wp), PUBLIC ::   rn_betas         !: coef. for partitioning of snowfall between leads and sea ice
                                          ! -- limthd_da -- !
   LOGICAL , PUBLIC ::   ln_limdA         !: activate lateral melting param. (T) or not (F)
   REAL(wp), PUBLIC ::   rn_beta          !: coef. beta for lateral melting param.
   REAL(wp), PUBLIC ::   rn_dmin          !: minimum floe diameter for lateral melting param.
                                          ! -- limthd_lac -- !
   LOGICAL , PUBLIC ::   ln_limdO         !: activate ice growth in open-water (T) or not (F)
   REAL(wp), PUBLIC ::   rn_hnewice       !: thickness for new ice formation (m)
   LOGICAL , PUBLIC ::   ln_frazil        !: use of frazil ice collection as function of wind (T) or not (F)
   REAL(wp), PUBLIC ::   rn_maxfrazb      !: maximum portion of frazil ice collecting at the ice bottom
   REAL(wp), PUBLIC ::   rn_vfrazb        !: threshold drift speed for collection of bottom frazil ice
   REAL(wp), PUBLIC ::   rn_Cfrazb        !: squeezing coefficient for collection of bottom frazil ice
                                          ! -- limitd_th -- !
   REAL(wp), PUBLIC ::   rn_himin         !: minimum ice thickness

   !                                     !!** ice-salinity namelist (namicesal) **
   LOGICAL , PUBLIC ::   ln_limdS         !: activate gravity drainage and flushing (T) or not (F)
   INTEGER , PUBLIC ::   nn_icesal        !: salinity configuration used in the model
   !                                      ! 1 - constant salinity in both space and time
   !                                      ! 2 - prognostic salinity (s(z,t))
   !                                      ! 3 - salinity profile, constant in time
   REAL(wp), PUBLIC ::   rn_icesal        !: bulk salinity (ppt) in case of constant salinity
   REAL(wp), PUBLIC ::   rn_sal_gd        !: restoring salinity for gravity drainage [PSU]
   REAL(wp), PUBLIC ::   rn_time_gd       !: restoring time constant for gravity drainage (= 20 days) [s]
   REAL(wp), PUBLIC ::   rn_sal_fl        !: restoring salinity for flushing [PSU]
   REAL(wp), PUBLIC ::   rn_time_fl       !: restoring time constant for gravity drainage (= 10 days) [s]
   REAL(wp), PUBLIC ::   rn_simax         !: maximum ice salinity [PSU]
   REAL(wp), PUBLIC ::   rn_simin         !: minimum ice salinity [PSU]

   !                                     !!** ice-mechanical redistribution namelist (namiceitdme)
   REAL(wp), PUBLIC ::   rn_cs            !: fraction of shearing energy contributing to ridging            
   INTEGER , PUBLIC ::   nn_partfun       !: participation function: =0 Thorndike et al. (1975), =1 Lipscomb et al. (2007)
   REAL(wp), PUBLIC ::   rn_gstar         !: fractional area of young ice contributing to ridging
   REAL(wp), PUBLIC ::   rn_astar         !: equivalent of G* for an exponential participation function
   LOGICAL , PUBLIC ::   ln_ridging       !: ridging of ice or not                        
   REAL(wp), PUBLIC ::   rn_hstar         !: thickness that determines the maximal thickness of ridged ice
   REAL(wp), PUBLIC ::   rn_por_rdg       !: initial porosity of ridges (0.3 regular value)
   REAL(wp), PUBLIC ::   rn_fsnowrdg      !: fractional snow loss to the ocean during ridging
   LOGICAL , PUBLIC ::   ln_rafting       !: rafting of ice or not                        
   REAL(wp), PUBLIC ::   rn_hraft         !: threshold thickness (m) for rafting / ridging 
   REAL(wp), PUBLIC ::   rn_craft         !: coefficient for smoothness of the hyperbolic tangent in rafting
   REAL(wp), PUBLIC ::   rn_fsnowrft      !: fractional snow loss to the ocean during ridging

   !                                     !!** some other parameters 
   INTEGER , PUBLIC ::   nstart           !: iteration number of the begining of the run 
   INTEGER , PUBLIC ::   nlast            !: iteration number of the end of the run 
   INTEGER , PUBLIC ::   nitrun           !: number of iteration
   INTEGER , PUBLIC ::   numit            !: iteration number
   REAL(wp), PUBLIC ::   rdt_ice          !: ice time step
   REAL(wp), PUBLIC ::   r1_rdtice        !: = 1. / rdt_ice
   REAL(wp), PUBLIC ::   r1_nlay_i        !: 1 / nlay_i
   REAL(wp), PUBLIC ::   r1_nlay_s        !: 1 / nlay_s 
   REAL(wp), PUBLIC ::   rswitch          !: switch for the presence of ice (1) or not (0)
   REAL(wp), PUBLIC, PARAMETER ::   epsi06   = 1.e-06_wp  !: small number 
   REAL(wp), PUBLIC, PARAMETER ::   epsi10   = 1.e-10_wp  !: small number 
   REAL(wp), PUBLIC, PARAMETER ::   epsi20   = 1.e-20_wp  !: small number 

   !                                     !!** define arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   u_oce, v_oce !: surface ocean velocity used in ice dynamics
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ahiu , ahiv !: hor. diffusivity coeff. at U- and V-points [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hicol       !: ice collection thickness accreted in leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   strength    !: ice strength
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   stress1_i, stress2_i, stress12_i   !: 1st, 2nd & diagonal stress tensor element
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   delta_i     !: ice rheology elta factor (Flato & Hibler 95) [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   divu_i      !: Divergence of the velocity field [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   shear_i     !: Shear of the velocity field [s-1]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   t_bo        !: Sea-Ice bottom temperature [Kelvin]     
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   frld        !: Leads fraction = 1 - ice fraction
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   pfrld       !: Leads fraction at previous time  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   phicif      !: Old ice thickness
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qlead       !: heat balance of the lead (or of the open ocean)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fhtur       !: net downward heat flux from the ice to the ocean
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fhld        !: heat flux from the lead used for bottom melting

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_snw     !: snow-ocean mass exchange   [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_spr     !: snow precipitation on ice  [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_sub     !: snow/ice sublimation       [kg.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_ice     !: ice-ocean mass exchange                   [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_sni     !: snow ice growth component of wfx_ice      [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_opw     !: lateral ice growth component of wfx_ice   [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_bog     !: bottom ice growth component of wfx_ice    [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_dyn     !: dynamical ice growth component of wfx_ice [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_bom     !: bottom melt component of wfx_ice          [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_sum     !: surface melt component of wfx_ice         [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_lam     !: lateral melt component of wfx_ice         [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_res     !: residual component of wfx_ice             [kg.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   afx_tot     !: ice concentration tendency (total)          [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   afx_thd     !: ice concentration tendency (thermodynamics) [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   afx_dyn     !: ice concentration tendency (dynamics)       [s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_bog     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_bom     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_lam     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_sum     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_sni     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_opw     !: salt flux due to ice growth/melt                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_bri     !: salt flux due to brine rejection                      [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_dyn     !: salt flux due to porous ridged ice formation          [PSU/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_res     !: residual salt flux due to correction of ice thickness [PSU/m2/s]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx_sub     !: salt flux due to ice sublimation

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_bog     !: total heat flux causing bottom ice growth        [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_bom     !: total heat flux causing bottom ice melt          [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_sum     !: total heat flux causing surface ice melt         [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_opw     !: total heat flux causing open water ice formation [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_dif     !: total heat flux causing Temp change in the ice   [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_snw     !: heat flux for snow melt                          [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_err     !: heat flux error after heat diffusion             [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_err_dif !: heat flux remaining due to change in non-solar flux [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_err_rem !: heat flux error after heat remapping             [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_in      !: heat flux available for thermo transformations   [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_out     !: heat flux remaining at the end of thermo transformations  [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wfx_err_sub !: mass flux error after sublimation [kg.m-2.s-1]
   
   ! heat flux associated with ice-atmosphere mass exchange
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_sub     !: heat flux for sublimation  [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_spr     !: heat flux of the snow precipitation  [W.m-2]

   ! heat flux associated with ice-ocean mass exchange
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_thd     !: ice-ocean heat flux from thermo processes (limthd_dh)  [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_dyn     !: ice-ocean heat flux from mecanical processes (limitd_me)  [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hfx_res     !: residual heat flux due to correction of ice thickness [W.m-2]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rn_amax_2d     !: maximum ice concentration 2d array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ftr_ice        !: transmitted solar radiation under ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   pahu3D, pahv3D !: ice hor. eddy diffusivity coef. at U- and V-points

   !!--------------------------------------------------------------------------
   !! * Ice global state variables
   !!--------------------------------------------------------------------------
   !! Variables defined for each ice category
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ht_i      !: Ice thickness (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i       !: Ice fractional areas (concentration)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   v_i       !: Ice volume per unit area (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   v_s       !: Snow volume per unit area(m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ht_s      !: Snow thickness (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   t_su      !: Sea-Ice Surface Temperature (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sm_i      !: Sea-Ice Bulk salinity (ppt)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   smv_i     !: Sea-Ice Bulk salinity times volume per area (ppt.m)
   !                                                                    !  this is an extensive variable that has to be transported
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   o_i       !: Sea-Ice Age (days)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   oa_i      !: Sea-Ice Age times ice area (days)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   bv_i      !: brine volume

   !! Variables summed over all categories, or associated to all the ice in a single grid cell
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   u_ice, v_ice !: components of the ice velocity (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vt_i , vt_s  !: ice and snow total volume per unit area (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   at_i         !: ice total fractional area (ice concentration)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ato_i        !: =1-at_i ; total open water fractional area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   et_i , et_s  !: ice and snow total heat content
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tm_i         !: mean ice temperature over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bvm_i        !: brine volume averaged over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   smt_i        !: mean sea ice salinity averaged over all categories [PSU]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tm_su        !: mean surface temperature over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htm_i        !: mean ice  thickness over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htm_s        !: mean snow thickness over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   om_i         !: mean ice age over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tau_icebfr   !: ice friction with bathy (landfast param activated)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   t_s      !: Snow temperatures [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_s      !: Snow ...      
      
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   t_i      !: ice temperatures          [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_i      !: ice thermal contents    [J/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   s_i      !: ice salinities          [PSU]

   !!--------------------------------------------------------------------------
   !! * Moments for advection
   !!--------------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sxopw, syopw, sxxopw, syyopw, sxyopw   !: open water in sea ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxice, syice, sxxice, syyice, sxyice   !: ice thickness 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxsn , sysn , sxxsn , syysn , sxysn    !: snow thickness
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxa  , sya  , sxxa  , syya  , sxya     !: lead fraction
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxc0 , syc0 , sxxc0 , syyc0 , sxyc0    !: snow thermal content
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxsal, sysal, sxxsal, syysal, sxysal   !: ice salinity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxage, syage, sxxage, syyage, sxyage   !: ice age
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxe  , sye  , sxxe  , syye  , sxye     !: ice layers heat content

   !!--------------------------------------------------------------------------
   !! * Old values of global variables
   !!--------------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_s_b, v_i_b               !: snow and ice volumes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_i_b, smv_i_b, oa_i_b     !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_s_b                      !: snow heat content
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_i_b                      !: ice temperatures
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   u_ice_b, v_ice_b           !: ice velocity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   at_i_b                     !: ice concentration (total)
            
   !!--------------------------------------------------------------------------
   !! * Ice thickness distribution variables
   !!--------------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   hi_max         !: Boundary of ice thickness categories in thickness space
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   hi_mean        !: Mean ice thickness in catgories 
   !
   !!--------------------------------------------------------------------------
   !! * Ice diagnostics
   !!--------------------------------------------------------------------------
   ! thd refers to changes induced by thermodynamics
   ! trp   ''         ''     ''       advection (transport of ice)
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_trp_vi   !: transport of ice volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_trp_vs   !: transport of snw volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_trp_ei   !: transport of ice enthalpy (W/m2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_trp_es   !: transport of snw enthalpy (W/m2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_trp_smv  !: transport of salt content
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_heat     !: snw/ice heat content variation   [W/m2] 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_smvi     !: ice salt content variation   [] 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_vice     !: ice volume variation   [m/s] 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_vsnw     !: snw volume variation   [m/s] 
   !
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2010)
   !! $Id: ice.F90 7813 2017-03-20 16:17:45Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION ice_alloc()
      !!-----------------------------------------------------------------
      !!               *** Routine ice_alloc ***
      !!-----------------------------------------------------------------
      INTEGER :: ice_alloc
      !
      INTEGER :: ierr(15), ii
      !!-----------------------------------------------------------------

      ierr(:) = 0

      ! What could be one huge allocate statement is broken-up to try to
      ! stay within Fortran's max-line length limit.
      ii = 1
      ALLOCATE( u_oce   (jpi,jpj) , v_oce    (jpi,jpj) ,                                             &
         &      ahiu    (jpi,jpj) , ahiv     (jpi,jpj) , hicol    (jpi,jpj) ,                        &
         &      strength(jpi,jpj) , stress1_i(jpi,jpj) , stress2_i(jpi,jpj) , stress12_i(jpi,jpj) ,  &
         &      delta_i (jpi,jpj) , divu_i   (jpi,jpj) , shear_i  (jpi,jpj) , STAT=ierr(ii) )

      ii = ii + 1
      ALLOCATE( t_bo   (jpi,jpj) , frld   (jpi,jpj) , pfrld  (jpi,jpj) , phicif (jpi,jpj) ,     &
         &      wfx_snw(jpi,jpj) , wfx_ice(jpi,jpj) , wfx_sub(jpi,jpj) , wfx_lam(jpi,jpj) ,     &
         &      wfx_bog(jpi,jpj) , wfx_dyn(jpi,jpj) , wfx_bom(jpi,jpj) , wfx_sum(jpi,jpj) ,     &
         &      wfx_res(jpi,jpj) , wfx_sni(jpi,jpj) , wfx_opw(jpi,jpj) , wfx_spr(jpi,jpj) ,     &
         &      afx_tot(jpi,jpj) , afx_thd(jpi,jpj),  afx_dyn(jpi,jpj) , rn_amax_2d(jpi,jpj),   &
         &      fhtur  (jpi,jpj) , qlead  (jpi,jpj) ,                                           &
         &      sfx_res(jpi,jpj) , sfx_bri(jpi,jpj) , sfx_dyn(jpi,jpj) , sfx_sub(jpi,jpj) , sfx_lam(jpi,jpj) ,  &
         &      sfx_bog(jpi,jpj) , sfx_bom(jpi,jpj) , sfx_sum(jpi,jpj) , sfx_sni(jpi,jpj) , sfx_opw(jpi,jpj) ,  &
         &      hfx_res(jpi,jpj) , hfx_snw(jpi,jpj) , hfx_sub(jpi,jpj) , hfx_err(jpi,jpj) ,     & 
         &      hfx_in (jpi,jpj) , hfx_out(jpi,jpj) , fhld   (jpi,jpj) ,                        &
         &      hfx_sum(jpi,jpj) , hfx_bom(jpi,jpj) , hfx_bog(jpi,jpj) , hfx_dif(jpi,jpj) ,     &
         &      hfx_opw(jpi,jpj) , hfx_thd(jpi,jpj) , hfx_dyn(jpi,jpj) , hfx_spr(jpi,jpj) ,     &
         &      hfx_err_dif(jpi,jpj) , hfx_err_rem(jpi,jpj) , wfx_err_sub(jpi,jpj)        ,  STAT=ierr(ii) )

      ! * Ice global state variables
      ii = ii + 1
      ALLOCATE( ftr_ice(jpi,jpj,jpl) , pahu3D(jpi,jpj,jpl+1) , pahv3D(jpi,jpj,jpl+1) , &
         &      ht_i   (jpi,jpj,jpl) , a_i   (jpi,jpj,jpl) , v_i   (jpi,jpj,jpl) ,     &
         &      v_s    (jpi,jpj,jpl) , ht_s  (jpi,jpj,jpl) , t_su  (jpi,jpj,jpl) ,     &
         &      sm_i   (jpi,jpj,jpl) , smv_i (jpi,jpj,jpl) , o_i   (jpi,jpj,jpl) ,     &
         &      oa_i   (jpi,jpj,jpl) , bv_i  (jpi,jpj,jpl) ,  STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( u_ice(jpi,jpj) , v_ice(jpi,jpj) ,                                       &
         &      vt_i (jpi,jpj) , vt_s (jpi,jpj) , at_i (jpi,jpj) , ato_i(jpi,jpj) ,     &
         &      et_i (jpi,jpj) , et_s (jpi,jpj) , tm_i (jpi,jpj) , bvm_i(jpi,jpj) ,     &
         &      smt_i(jpi,jpj) , tm_su(jpi,jpj) , htm_i(jpi,jpj) , htm_s(jpi,jpj) ,     &
         &      om_i (jpi,jpj) , tau_icebfr(jpi,jpj)                              , STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( t_s(jpi,jpj,nlay_s,jpl) , e_s(jpi,jpj,nlay_s,jpl) , STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( t_i(jpi,jpj,nlay_i,jpl) , e_i(jpi,jpj,nlay_i,jpl) , s_i(jpi,jpj,nlay_i,jpl) , STAT=ierr(ii) )

      ! * Moments for advection
      ii = ii + 1
      ALLOCATE( sxopw(jpi,jpj) , syopw(jpi,jpj) , sxxopw(jpi,jpj) , syyopw(jpi,jpj) , sxyopw(jpi,jpj) , STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( sxice(jpi,jpj,jpl) , syice(jpi,jpj,jpl) , sxxice(jpi,jpj,jpl) , syyice(jpi,jpj,jpl) , sxyice(jpi,jpj,jpl) ,   &
         &      sxsn (jpi,jpj,jpl) , sysn (jpi,jpj,jpl) , sxxsn (jpi,jpj,jpl) , syysn (jpi,jpj,jpl) , sxysn (jpi,jpj,jpl) ,   &
         &      STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( sxa  (jpi,jpj,jpl) , sya  (jpi,jpj,jpl) , sxxa  (jpi,jpj,jpl) , syya  (jpi,jpj,jpl) , sxya  (jpi,jpj,jpl) ,   &
         &      sxc0 (jpi,jpj,jpl) , syc0 (jpi,jpj,jpl) , sxxc0 (jpi,jpj,jpl) , syyc0 (jpi,jpj,jpl) , sxyc0 (jpi,jpj,jpl) ,   &
         &      sxsal(jpi,jpj,jpl) , sysal(jpi,jpj,jpl) , sxxsal(jpi,jpj,jpl) , syysal(jpi,jpj,jpl) , sxysal(jpi,jpj,jpl) ,   &
         &      sxage(jpi,jpj,jpl) , syage(jpi,jpj,jpl) , sxxage(jpi,jpj,jpl) , syyage(jpi,jpj,jpl) , sxyage(jpi,jpj,jpl) ,   &
         &      STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( sxe (jpi,jpj,nlay_i,jpl) , sye (jpi,jpj,nlay_i,jpl) , sxxe(jpi,jpj,nlay_i,jpl) ,     &
         &      syye(jpi,jpj,nlay_i,jpl) , sxye(jpi,jpj,nlay_i,jpl)                            , STAT=ierr(ii) )

      ! * Old values of global variables
      ii = ii + 1
      ALLOCATE( v_s_b  (jpi,jpj,jpl) , v_i_b  (jpi,jpj,jpl) , e_s_b(jpi,jpj,nlay_s,jpl) ,     &
         &      a_i_b  (jpi,jpj,jpl) , smv_i_b(jpi,jpj,jpl) , e_i_b(jpi,jpj,nlay_i,jpl) ,     &
         &      oa_i_b (jpi,jpj,jpl)                                                    , STAT=ierr(ii) )
      ii = ii + 1
      ALLOCATE( u_ice_b(jpi,jpj) , v_ice_b(jpi,jpj) , at_i_b(jpi,jpj) , STAT=ierr(ii) )
      
      ! * Ice thickness distribution variables
      ii = ii + 1
      ALLOCATE( hi_max(0:jpl), hi_mean(jpl),  STAT=ierr(ii) )

      ! * Ice diagnostics
      ii = ii + 1
      ALLOCATE( diag_trp_vi(jpi,jpj) , diag_trp_vs (jpi,jpj) , diag_trp_ei(jpi,jpj),   & 
         &      diag_trp_es(jpi,jpj) , diag_trp_smv(jpi,jpj) , diag_heat  (jpi,jpj),   &
         &      diag_smvi  (jpi,jpj) , diag_vice   (jpi,jpj) , diag_vsnw  (jpi,jpj), STAT=ierr(ii) )

      ice_alloc = MAXVAL( ierr(:) )
      IF( ice_alloc /= 0 )   CALL ctl_warn('ice_alloc: failed to allocate arrays.')
      !
   END FUNCTION ice_alloc

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module            NO LIM sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE ice
