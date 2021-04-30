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
   !!                   namcfg namelist parameters
   !!----------------------------------------------------------------------
   LOGICAL       ::   ln_read_cfg      !: (=T) read the domain configuration file or (=F) not
   CHARACTER(lc) ::      cn_domcfg        !: filename the configuration file to be read
   LOGICAL       ::   ln_write_cfg     !: (=T) create the domain configuration file
   CHARACTER(lc) ::      cn_domcfg_out    !: filename the configuration file to be read
   !
   LOGICAL       ::   ln_use_jattr     !: input file read offset
   !                                   !  Use file global attribute: open_ocean_jstart to determine start j-row 
   !                                   !  when reading input from those netcdf files that have the 
   !                                   !  attribute defined. This is designed to enable input files associated 
   !                                   !  with the extended grids used in the under ice shelf configurations to 
   !                                   !  be used without redundant rows when the ice shelves are not in use.
   ! 

   !!---------------------------------------------------------------------
   !! Domain Matrix size 
   !!---------------------------------------------------------------------
   ! configuration name & resolution   (required only in ORCA family case)
   CHARACTER(lc) ::   cn_cfg           !: name of the configuration
   INTEGER       ::   nn_cfg           !: resolution of the configuration 

   ! global domain size               !!! * total computational domain *
   INTEGER       ::   jpiglo           !: 1st dimension of global domain --> i-direction
   INTEGER       ::   jpjglo           !: 2nd    -                  -    --> j-direction
   INTEGER       ::   jpkglo           !: 3nd    -                  -    --> k levels

   ! global domain size for AGRIF     !!! * total AGRIF computational domain *
   INTEGER, PUBLIC            ::   nbug_in_agrif_conv_do_not_remove_or_modify = 1 - 1
   INTEGER, PUBLIC, PARAMETER ::   nbghostcells = 3                             !: number of ghost cells
   INTEGER, PUBLIC            ::   nbcellsx   ! = jpiglo - 2 - 2*nbghostcells   !: number of cells in i-direction
   INTEGER, PUBLIC            ::   nbcellsy   ! = jpjglo - 2 - 2*nbghostcells   !: number of cells in j-direction

   ! local domain size                !!! * local computational domain *
   INTEGER, PUBLIC ::   jpi   !                                                    !: first  dimension
   INTEGER, PUBLIC ::   jpj   !                                                    !: second dimension
   INTEGER, PUBLIC ::   jpk   ! = jpkglo                                           !: third  dimension
   INTEGER, PUBLIC ::   jpim1 ! = jpi-1                                            !: inner domain indices
   INTEGER, PUBLIC ::   jpjm1 ! = jpj-1                                            !:   -     -      -
   INTEGER, PUBLIC ::   jpkm1 ! = jpk-1                                            !:   -     -      -
   INTEGER, PUBLIC ::   jpij  ! = jpi*jpj                                          !:  jpi x jpj
   INTEGER, PUBLIC ::   jpimax! = ( jpiglo-2*nn_hls + (jpni-1) ) / jpni + 2*nn_hls !: maximum jpi
   INTEGER, PUBLIC ::   jpjmax! = ( jpjglo-2*nn_hls + (jpnj-1) ) / jpnj + 2*nn_hls !: maximum jpj

   !!---------------------------------------------------------------------
   !! Active tracer parameters
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER ::   jpts   = 2    !: Number of active tracers (=2, i.e. T & S )
   INTEGER, PUBLIC, PARAMETER ::   jp_tem = 1    !: indice for temperature
   INTEGER, PUBLIC, PARAMETER ::   jp_sal = 2    !: indice for salinity
   !annkat
   INTEGER, PUBLIC, PARAMETER ::   jp_dep = 3    !: indice for depth
   INTEGER, PUBLIC, PARAMETER ::   jp_msk = 4    !: indice for mask
   !!----------------------------------------------------------------------
   !!   Domain decomposition
   !!----------------------------------------------------------------------
   !! if we dont use massively parallel computer (parameters jpni=jpnj=1) so jpiglo=jpi and jpjglo=jpj
   INTEGER, PUBLIC            ::   jpni         !: number of processors following i 
   INTEGER, PUBLIC            ::   jpnj         !: number of processors following j
   INTEGER, PUBLIC            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )
   INTEGER, PUBLIC, PARAMETER ::   jpr2di = 0   !: number of columns for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpr2dj = 0   !: number of rows    for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   nn_hls = 1   !: halo width (applies to both rows and columns)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE par_oce
