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
   ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)  
   ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
   ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
   rn_hmin     =   -10     !  min depth of the ocean (>0) or 
                           !  min number of ocean level (<0)

!-----------------------------------------------------------------------
!   s-coordinate or hybrid z-s-coordinate
!-----------------------------------------------------------------------
   rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
   rn_sbot_max = 6000.     !  maximum depth of s-bottom surface 
                           !  (= ocean depth) (>0) (m)
   ln_s_sigma  = .true.   !  hybrid s-sigma coordinates
   rn_hc       =  39.0    !  critical depth with s-sigma

!-----------------------------------------------------------------------
!  grid information 
!-----------------------------------------------------------------------
   sn_src_hgr = './mesh_hgr_src_N06.nc'   !  parent /grid/
   sn_src_zgr = './mesh_zgr_src_N06_renamed.nc'   !  parent
   sn_dst_hgr = './domain_cfg.nc'
   sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
   sn_src_msk = './mask_src_N06.nc'       ! parent
   sn_bathy   = './bathy_meter.nc'        ! child 

!-----------------------------------------------------------------------
!  I/O 
!-----------------------------------------------------------------------
   sn_src_dir = 'ORCA0083_N06_1960.ncml'       ! src_files/'
   sn_dst_dir = '/projectsa/accord/SEAsia/START_FILES/OPEN_BOUNDARIES/'
   sn_fn      = 'SEAsia'                 ! prefix for output files
   nn_fv      = -1e20                     !  set fill value for output files
   nn_src_time_adj = 0					  ! src time adjustment
   sn_dst_metainfo = 'metadata info: jelt'

!-----------------------------------------------------------------------
!  unstructured open boundaries                         
!-----------------------------------------------------------------------
    ln_coords_file = .true.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .false.               !  boundary conditions for barotropic fields
    ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
    ln_tra         = .false.               !  boundary conditions for T and S
    ln_ice         = .false.               !  ice boundary condition   
    nn_rimwidth    = 1                    !  width of the relaxation zone

!-----------------------------------------------------------------------
!  unstructured open boundaries tidal parameters                        
!-----------------------------------------------------------------------
    ln_tide        = .true.               !  =T : produce bdy tidal conditions
!FES
        clname(1) ='2N2'
        clname(2)='EPS2'
        clname(3)='J1'
        clname(4)='K1'
        clname(5)='K2'
        clname(6)='L2'
        clname(7)='LA2'
        clname(8)='M2'
        clname(9)='M3'
        clname(10)='M4'
        clname(11)='M6'
        clname(12)='M8'
        clname(13)='MF'
        clname(14)='MKS2'
        clname(15)='MM'
        clname(16)='MN4'
        clname(17)='MS4'
        clname(18)='MSF'
        clname(19)='MSQM'
        clname(20)='MTM'
        clname(21)='MU2'
        clname(22)='N2'
        clname(23)='N4'
        clname(24)='NU2'
        clname(25)='O1'
        clname(26)='P1'
        clname(27)='Q1'
        clname(28)='R2'
        clname(29)='S1'
        clname(30)='S2'
        clname(31)='S4'
        clname(32)='SA'
        clname(33)='SSA'
        clname(34)='T2'
    ln_trans       = .false.
    sn_tide_h     = '' # Hardwired with FES14 path
    sn_tide_u     = '' # Hardwired with FES14 path

!-----------------------------------------------------------------------
!  Time information
!-----------------------------------------------------------------------
    nn_year_000     = 1960        !  year start
    nn_year_end     = 1960          !  year end
    nn_month_000    = 1          !  month start (default = 1 is years>1)
    nn_month_end    = 12         !  month end (default = 12 is years>1)
    sn_dst_calendar = 'gregorian' !  output calendar format
    nn_base_year    = 1900        !  base year for time counter in output data
    sn_tide_grid    = '' # Hard wired for FES14 paths #/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
    nn_src_time_adj = 0  != -3168000 - 86400 ! N01: fix to align model time stamp
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
    rn_mask_max_depth = 300.0	  !  Maximum depth to be ignored for the mask
    rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break
