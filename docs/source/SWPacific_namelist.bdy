!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  : namelist for BDY generation tool
!!            
!!             User inputs for generating open boundary conditions
!!             employed by the BDY module in NEMO. Boundary data
!!             can be set up for v3.2 NEMO and above.
!!            
!!             More info here.....
!!            
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!-----------------------------------------------------------------------
!   vertical coordinate
!-----------------------------------------------------------------------
   ln_zco      = .true.   !  z-coordinate - full    steps   (T/F)  
   ln_zps      = .false.    !  z-coordinate - partial steps   (T/F)
   ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
   rn_hmin     =   -5     !  min depth of the ocean (>0) or 
                           !  min number of ocean level (<0)

!-----------------------------------------------------------------------
!   s-coordinate or hybrid z-s-coordinate
!-----------------------------------------------------------------------
   rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
   rn_sbot_max = 7000.     !  maximum depth of s-bottom surface 
                           !  (= ocean depth) (>0) (m)
   ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
   rn_hc       =  50.0    !  critical depth with s-sigma

!-----------------------------------------------------------------------
!  grid information 
!-----------------------------------------------------------------------
   sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
   sn_src_zgr = './mesh_zgr_src.nc'   !  parent
   sn_dst_hgr = './domain_cfg.nc'
   sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
   sn_src_msk = './mask_src.nc'       ! parent
   sn_bathy   = './bathy_meter.nc'

!-----------------------------------------------------------------------
!  I/O 
!-----------------------------------------------------------------------
   sn_src_dir = './inputs_src.ncml'       ! src_files/'
   sn_dst_dir = '/work/jelt/NEMO/SWPacific/INPUTS/'
   sn_fn      = 'SWPacific'                 ! prefix for output files
   nn_fv      = -1e20                     !  set fill value for output files
   nn_src_time_adj = 0					  ! src time adjustment
   sn_dst_metainfo = 'metadata info: jelt'

!-----------------------------------------------------------------------
!  unstructured open boundaries                         
!-----------------------------------------------------------------------
    ln_coords_file = .true.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
    ln_mask_file   = .true.              !  =T : read mask from file
    cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
    ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
    ln_tra         = .false.               !  boundary conditions for T and S
    ln_ice         = .false.               !  ice boundary condition   
    nn_rimwidth    = 1                    !  width of the relaxation zone

!-----------------------------------------------------------------------
!  unstructured open boundaries tidal parameters                        
!-----------------------------------------------------------------------
    ln_tide        = .true.               !  =T : produce bdy tidal conditions
    clname(1)      = 'M2'                 ! constituent name
    clname(2)      = 'S2'         
    clname(3)      = 'K2' 
    ln_trans       = .false.
    sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
    sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'

!-----------------------------------------------------------------------
!  Time information
!-----------------------------------------------------------------------
    nn_year_000     = 1979        !  year start
    nn_year_end     = 1979        !  year end
    nn_month_000    = 11          !  month start (default = 1 is years>1)
    nn_month_end    = 11          !  month end (default = 12 is years>1)
    sn_dst_calendar = 'gregorian' !  output calendar format
    nn_base_year    = 1978        !  base year for time counter
    sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'

!-----------------------------------------------------------------------
!  Additional parameters
!-----------------------------------------------------------------------
    nn_wei  = 1                   !  smoothing filter weights 
    rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                  !  smoothing onto dst points. Need to 
                                  !  make this a funct. of dlon
    sn_history  = 'bdy files produced by jelt from ORCA0083-N01'
                                  !  history for netcdf file
    ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
    nn_alpha    = 0               !  Euler rotation angle
    nn_beta     = 0               !  Euler rotation angle
    nn_gamma    = 0               !  Euler rotation angle
    rn_mask_max_depth = 7000.0	  !  Maximum depth to be ignored for the mask
    rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break
