MODULE par_oce
   !!======================================================================
   !!                        ***  par_oce  ***
   !! Ocean :   set the ocean parameters
   !!======================================================================
   !! History :  OPA  !  1991     (Imbard, Levy, Madec)  Original code
   !!   NEMO     1.0  !  2004-01  (G. Madec, J.-M. Molines)  Free form and module
   !!            3.3  !  2010-09  (C. Ethe) TRA-TRC merge: add jpts, jp_tem & jp_sal
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!   Domain decomposition
   !!----------------------------------------------------------------------
   !! if we dont use massively parallel computer (parameters jpni=jpnj=1) so jpiglo=jpi and jpjglo=jpj
   INTEGER, PUBLIC            ::   jpni         !: number of processors following i 
   INTEGER, PUBLIC            ::   jpnj         !: number of processors following j
   INTEGER, PUBLIC            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )
   INTEGER, PUBLIC, PARAMETER ::   jpr2di = 0   !: number of columns for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpr2dj = 0   !: number of rows    for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpreci = 1   !: number of columns for overlap 
   INTEGER, PUBLIC, PARAMETER ::   jprecj = 1   !: number of rows    for overlap 

   !!----------------------------------------------------------------------
   !!                   namcfg namelist parameters
   !!----------------------------------------------------------------------
   CHARACTER(lc) ::   cp_cfg           !: name of the configuration
   CHARACTER(lc) ::   cp_cfz           !: name of the zoom of configuration
   INTEGER       ::   jp_cfg           !: resolution of the configuration

   ! data size                                       !!! * size of all input files *
   INTEGER       ::   jpidta           !: 1st lateral dimension ( >= jpi )
   INTEGER       ::   jpjdta           !: 2nd    "         "    ( >= jpj )
   INTEGER       ::   jpkdta           !: number of levels      ( >= jpk )

   ! global or zoom domain size                      !!! * computational domain *
   INTEGER       ::   jpiglo           !: 1st dimension of global domain --> i
   INTEGER       ::   jpjglo           !: 2nd    -                  -    --> j

   ! zoom starting position 
   INTEGER       ::   jpizoom          !: left bottom (i,j) indices of the zoom
   INTEGER       ::   jpjzoom          !: in data domain indices

   ! Domain characteristics
   INTEGER       ::   jperio           !: lateral cond. type (between 0 and 6)
   !                                       !  = 0 closed                 ;   = 1 cyclic East-West
   !                                       !  = 2 equatorial symmetric   ;   = 3 North fold T-point pivot
   !                                       !  = 4 cyclic East-West AND North fold T-point pivot
   !                                       !  = 5 North fold F-point pivot
   !                                       !  = 6 cyclic East-West AND North fold F-point pivot

   ! Input file read offset
   LOGICAL       ::   ln_use_jattr     !: Use file global attribute: open_ocean_jstart to determine start j-row 
                                           ! when reading input from those netcdf files that have the 
                                           ! attribute defined. This is designed to enable input files associated 
                                           ! with the extended grids used in the under ice shelf configurations to 
                                           ! be used without redundant rows when the ice shelves are not in use.

   !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
   !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
   REAL(wp)      ::   pp_not_used       = 999999._wp   !: vertical grid parameter
   REAL(wp)      ::   pp_to_be_computed = 999999._wp   !:    -      -       -




   !!---------------------------------------------------------------------
   !! Active tracer parameters
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER ::   jpts   = 2    !: Number of active tracers (=2, i.e. T & S )
   INTEGER, PUBLIC, PARAMETER ::   jp_tem = 1    !: indice for temperature
   INTEGER, PUBLIC, PARAMETER ::   jp_sal = 2    !: indice for salinity
   INTEGER, PUBLIC, PARAMETER ::   jp_dep = 3    !: indice for depth
   INTEGER, PUBLIC, PARAMETER ::   jp_msk = 4    !: indice for depth


   !!---------------------------------------------------------------------
   !! Domain Matrix size  (if AGRIF, they are not all parameters)
   !!---------------------------------------------------------------------
#if defined key_agrif
   INTEGER, PUBLIC, PARAMETER ::   nbghostcells = 1                             !: number of ghost cells
   INTEGER, PUBLIC            ::   nbcellsx     = jpiglo - 2 - 2*nbghostcells   !: number of cells in i-direction
   INTEGER, PUBLIC            ::   nbcellsy     = jpjglo - 2 - 2*nbghostcells   !: number of cells in j-direction
   !
#endif
   INTEGER, PUBLIC  ::   jpi   ! = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   !: first  dimension
   INTEGER, PUBLIC  ::   jpj   ! = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   !: second dimension
   INTEGER, PUBLIC  ::   jpk   ! = jpkdta
   INTEGER, PUBLIC  ::   jpim1 ! = jpi-1                                            !: inner domain indices
   INTEGER, PUBLIC  ::   jpjm1 ! = jpj-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpkm1 ! = jpk-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpij  ! = jpi*jpj                                          !:  jpi x jpj

   !!---------------------------------------------------------------------
   !! Optimization/control flags
   !!---------------------------------------------------------------------
#if defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .TRUE.   !: flag to activate the all options
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_oce
