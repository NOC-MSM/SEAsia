Anatomy of the namelist.bdy
+++++++++++++++++++++++++++

This is a brief introduction to the content and use of the ``namelist.bdy`` file
which controls the Nemo Relocatable Configuration Tool [NRCT] nee PyNEMO.

Header::

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

First the vertical coordinates are defined. This is for the *TARGET* grid.
Only one boolean field can be true.**IS THIS TRUE? OR ln_zps x ln_sco=1 => hybrid z-s?**
The minimum depth of the ocean can be set
in terms of the number of model levels (useful for limiting the number of levels
lost to outcropping over sloping bathymtry, in z-coords) or as a depth.
 **CHECK: DOES rn_hmin WORK IN S-COORDS?**
::

  !-----------------------------------------------------------------------
  !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zco      = .true.   !  z-coordinate - full    steps   (T/F)
     ln_zps      = .false.    !  z-coordinate - partial steps   (T/F)
     ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
     rn_hmin     =   -5     !  min depth of the ocean (>0) or
                             !  min number of ocean level (<0)

The next section add details if either ``ln_sco = T`` or ``ln_sco=T & ln_zco = T``
**CORRECT?**
The max and min limits of the top of the deepest layer are set. e.g. the bottom
 layer is here 10m (below the surface) from outcropping and bathymetry below 7000m
 is filled in (if ``ln_s_sigma = F``) or treated a z-coordinates (if ``ln_s_sigma = T``).
 **IS THIS CORRECT?**
 The depth ``rn_hc`` controls .... **NO IDEA**
::

  !-----------------------------------------------------------------------
  !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
     rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
     rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
                             !  (= ocean depth) (>0) (m)
     ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
     rn_hc       =  50.0    !  critical depth with s-sigma


This section controls where the tool is to look for *source* (src) and
 *destination* (dst) datasets. A neat thing with NRCT is that the use of NetCDF Markup Language (ncml)
 files allows you to easily relabel variables in netcdf files.

The first two files are the horizontal and vertical mesh files from the parent grid.
The ``domain_cfg.nc`` file contains the bundled mesh and mask information
(NEMOv4 style) for the target grid.
The next file, ``inputs_dst.ncml`` ...
::

  !-----------------------------------------------------------------------
  !  grid information
  !-----------------------------------------------------------------------
     sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
     sn_src_zgr = './mesh_zgr_src.nc'   !  parent
     sn_dst_hgr = './domain_cfg.nc'
     sn_dst_zgr = './inputs_dst.ncml' ! rename output variables
     sn_src_msk = './mask_src.nc'       ! parent
     sn_bathy   = './bathy_meter.nc'

::

  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
     sn_src_dir = './inputs_src.ncml'       ! src_files/'
     sn_dst_dir = '/work/jelt/NEMO/SWPacific/INPUTS/'
     sn_fn      = 'SWPacific'                 ! prefix for output files
     nn_fv      = -1e20                     !  set fill value for output files
     nn_src_time_adj = 0					  ! src time adjustment
     sn_dst_metainfo = 'metadata info: jelt'

::

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_coords_file = .true.               !  =T : produce bdy coordinates files
      cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
      ln_mask_file   = .true.              !  =T : read mask from file
      cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .false.               !  boundary conditions for baroclinic velocities
      ln_tra         = .true.               !  boundary conditions for T and S
      ln_ice         = .false.               !  ice boundary condition
      nn_rimwidth    = 1                    !  width of the relaxation zone

::

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

Time information: The start and end years are self explanatory. **The comments
against the months don't make sense; are 12 months per year always computed if more
than one year is requested?**

``sn_dst_calendar`` refers to the calendar setting that NEMO uses**? Not sure what the options are**

``nn_base_year`` refers to the reference date in the NEMO simulations:
 typically ``time_counter:time_origin``

 ``sn_tide_grid`` is the grid information for the tidal forcing
 **WHY IS IT NOT WITH ``sn_tide_h`` etc?**
::

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


Additional parameters. It is a good idea to put as much information into the history
tag as possible **I DO NOT KNOW HOW ANY OF THE OTHER TERMS FUNCTION**
::

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
