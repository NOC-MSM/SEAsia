MODULE iom_nf90
   !!=====================================================================
   !!                    ***  MODULE  iom_nf90 ***
   !! Input/Output manager :  Library to read input files with NF90 (only fliocom module)
   !!====================================================================
   !! History :  9.0  ! 05 12  (J. Belier) Original code
   !!            9.0  ! 06 02  (S. Masson) Adaptation to NEMO
   !!             "   ! 07 07  (D. Storkey) Changes to iom_nf90_gettime
   !!--------------------------------------------------------------------
   !!gm  caution add !DIR nec: improved performance to be checked as well as no result changes

   !!--------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_gettime    : read the time axis kvid in the file
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!--------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE netcdf          ! NetCDF library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC iom_nf90_open, iom_nf90_close, iom_nf90_varid, iom_nf90_get, iom_nf90_gettime, iom_nf90_rstput
   PUBLIC iom_nf90_getatt, iom_nf90_putatt

   INTERFACE iom_nf90_get
      MODULE PROCEDURE iom_nf90_g0d, iom_nf90_g123d
   END INTERFACE
   INTERFACE iom_nf90_getatt
      MODULE PROCEDURE iom_nf90_giatt, iom_nf90_gratt, iom_nf90_gcatt
   END INTERFACE
   INTERFACE iom_nf90_putatt
      MODULE PROCEDURE iom_nf90_piatt, iom_nf90_pratt, iom_nf90_pcatt
   END INTERFACE
   INTERFACE iom_nf90_rstput
      MODULE PROCEDURE iom_nf90_rp0123d
   END INTERFACE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE iom_nf90_open( cdname, kiomid, ldwrt, ldok, kdompar )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose : open an input file with NF90
      !!---------------------------------------------------------------------
      CHARACTER(len=*)       , INTENT(inout)           ::   cdname      ! File name
      INTEGER                , INTENT(  out)           ::   kiomid      ! nf90 identifier of the opened file
      LOGICAL                , INTENT(in   )           ::   ldwrt       ! read or write the file?
      LOGICAL                , INTENT(in   )           ::   ldok        ! check the existence 
      INTEGER, DIMENSION(2,5), INTENT(in   ), OPTIONAL ::   kdompar     ! domain parameters: 

      CHARACTER(LEN=256) ::   clinfo           ! info character
      CHARACTER(LEN=256) ::   cltmp            ! temporary character
      INTEGER            ::   iln              ! lengths of character
      INTEGER            ::   istop            ! temporary storage of nstop
      INTEGER            ::   if90id           ! nf90 identifier of the opened file
      INTEGER            ::   idmy             ! dummy variable
      INTEGER            ::   jl               ! loop variable
      INTEGER            ::   ichunk           ! temporary storage of nn_chunksz
      INTEGER            ::   imode            ! creation mode flag: NF90_CLOBBER or NF90_NOCLOBBER or NF90_HDF5
      INTEGER            ::   ihdf5            ! local variable for retrieval of value for NF90_HDF5
      LOGICAL            ::   llclobber        ! local definition of ln_clobber
      !---------------------------------------------------------------------

      clinfo = '                    iom_nf90_open ~~~  '
      istop = nstop   ! store the actual value of nstop
      IF( nn_chunksz > 0 ) THEN   ;   ichunk = nn_chunksz
      ELSE                        ;   ichunk = NF90_SIZEHINT_DEFAULT
      ENDIF
      !
      llclobber = ldwrt .AND. ln_clobber
      IF( ldok .AND. .NOT. llclobber ) THEN      ! Open existing file...
         !                 ! =============
         IF( ldwrt ) THEN  ! ... in write mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in WRITE mode'
            IF( snc4set%luse ) THEN
               CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_WRITE  , if90id ), clinfo)
            ELSE
               CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_WRITE  , if90id, chunksize = ichunk ), clinfo)
            ENDIF
            CALL iom_nf90_check(NF90_SET_FILL( if90id, NF90_NOFILL, idmy                          ), clinfo)
         ELSE              ! ... in read mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in READ mode'
            CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_NOWRITE, if90id, chunksize = ichunk ), clinfo)
         ENDIF
      ELSE                                       ! the file does not exist (or we overwrite it)
         !                 ! =============
         iln = INDEX( cdname, '.nc' )
         IF( ldwrt ) THEN  ! the file should be open in write mode so we create it...
            IF( jpnij > 1 ) THEN
               WRITE(cltmp,'(a,a,i4.4,a)') cdname(1:iln-1), '_', narea-1, '.nc'
               cdname = TRIM(cltmp)
            ENDIF
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' create new file: '//TRIM(cdname)//' in WRITE mode'

            IF( llclobber ) THEN   ;   imode = IOR( NF90_64BIT_OFFSET, NF90_CLOBBER   )
            ELSE                   ;   imode = IOR( NF90_64BIT_OFFSET, NF90_NOCLOBBER ) 
            ENDIF
            IF( snc4set%luse ) THEN
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' creating file: '//TRIM(cdname)//' in hdf5 (netcdf4) mode'
               CALL GET_NF90_SYMBOL("NF90_HDF5", ihdf5)
               IF( llclobber ) THEN   ;   imode = IOR(ihdf5, NF90_CLOBBER)
               ELSE                   ;   imode = IOR(ihdf5, NF90_NOCLOBBER)
               ENDIF
               CALL iom_nf90_check(NF90_CREATE( TRIM(cdname), imode, if90id ), clinfo)
            ELSE
               CALL iom_nf90_check(NF90_CREATE( TRIM(cdname), imode, if90id, chunksize = ichunk ), clinfo)
            ENDIF
            CALL iom_nf90_check(NF90_SET_FILL( if90id, NF90_NOFILL, idmy                     ), clinfo)
            ! define dimensions
            CALL iom_nf90_check(NF90_DEF_DIM( if90id, 'x', kdompar(1,1)  , idmy ), clinfo)
            CALL iom_nf90_check(NF90_DEF_DIM( if90id, 'y', kdompar(2,1)  , idmy ), clinfo)
            CALL iom_nf90_check(NF90_DEF_DIM( if90id, 'z', jpk           , idmy ), clinfo)
            CALL iom_nf90_check(NF90_DEF_DIM( if90id, 't', NF90_UNLIMITED, idmy ), clinfo)
            ! global attributes
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_number_total'   , jpnij              ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_number'         , narea-1            ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , (/1     , 2     /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_size_global'    , (/jpiglo, jpjglo/) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_size_local'     , kdompar(:,1)       ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_position_first' , kdompar(:,2)       ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_position_last'  , kdompar(:,3)       ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_halo_size_start', kdompar(:,4)       ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , kdompar(:,5)       ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_type'           , 'BOX'              ), clinfo)
         ELSE              ! the file should be open for read mode so it must exist...
            CALL ctl_stop( TRIM(clinfo), ' should be impossible case...' )
         ENDIF
      ENDIF
      ! start to fill file informations
      ! =============
      IF( istop == nstop ) THEN   ! no error within this routine
!does not work with some compilers         kiomid = MINLOC(iom_file(:)%nfid, dim = 1)
         kiomid = 0
         DO jl = jpmax_files, 1, -1
            IF( iom_file(jl)%nfid == 0 )   kiomid = jl
         ENDDO
         iom_file(kiomid)%name   = TRIM(cdname)
         iom_file(kiomid)%nfid   = if90id
         iom_file(kiomid)%iolib  = jpnf90
         iom_file(kiomid)%nvars  = 0
         iom_file(kiomid)%irec   = -1   ! useless for NetCDF files, used to know if the file is in define mode 
         CALL iom_nf90_check(NF90_Inquire(if90id, unlimitedDimId = iom_file(kiomid)%iduld), clinfo)
         IF ( iom_file(kiomid)%iduld .GE. 0 ) THEN
           CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, iom_file(kiomid)%iduld,   &
        &                                               name = iom_file(kiomid)%uldname), clinfo)
         ENDIF
         IF(lwp) WRITE(numout,*) '                   ---> '//TRIM(cdname)//' OK'
      ELSE
         kiomid = 0               ! return error flag
      ENDIF
      !
   END SUBROUTINE iom_nf90_open


   SUBROUTINE iom_nf90_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_nf90_close  ***
      !!
      !! ** Purpose : close an input file with NF90
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kiomid   ! iom identifier of the file to be closed
      CHARACTER(LEN=100)  ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      clinfo = '      iom_nf90_close    , file: '//TRIM(iom_file(kiomid)%name)
      CALL iom_nf90_check(NF90_CLOSE(iom_file(kiomid)%nfid), clinfo)
      !    
   END SUBROUTINE iom_nf90_close


   FUNCTION iom_nf90_varid ( kiomid, cdvar, kiv, kdimsz, kndims )  
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file with NF90
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER              , INTENT(in   )           ::   kiv   ! 
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of the dimensions
      INTEGER,               INTENT(  out), OPTIONAL ::   kndims   ! size of the dimensions
      !
      INTEGER                        ::   iom_nf90_varid   ! iom variable Id
      INTEGER                        ::   if90id           ! nf90 file identifier
      INTEGER                        ::   ji               ! dummy loop index
      INTEGER                        ::   ivarid           ! NetCDF  variable Id
      INTEGER                        ::   i_nvd            ! number of dimension of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   idimid           ! dimension ids of the variable
      LOGICAL                        ::   llok             ! ok  test
      CHARACTER(LEN=100)             ::   clinfo           ! info character
      !!-----------------------------------------------------------------------
      clinfo = '          iom_nf90_varid, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
      iom_nf90_varid = 0                    ! default definition
      IF( PRESENT(kdimsz) ) kdimsz(:) = 0   ! default definition
      if90id = iom_file(kiomid)%nfid        ! get back NetCDF file id
      !
      llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr   ! does the variable exist in the file
      IF( llok ) THEN
         iom_nf90_varid = kiv
         iom_file(kiomid)%nvars       = kiv
         iom_file(kiomid)%nvid(kiv)   = ivarid
         iom_file(kiomid)%cn_var(kiv) = TRIM(cdvar)
         CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, ndims = i_nvd), clinfo)   ! number of dimensions
         iom_file(kiomid)%ndims(kiv)  = i_nvd
         CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, dimids = idimid(1:i_nvd)), clinfo)   ! dimensions ids
         iom_file(kiomid)%luld(kiv) = .FALSE.   ! default value
         iom_file(kiomid)%dimsz(:,kiv) = 0      ! reset dimsz in case previously used
         DO ji = 1, i_nvd                       ! dimensions size
            CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, idimid(ji), len = iom_file(kiomid)%dimsz(ji,kiv)), clinfo)   
            IF( idimid(ji) == iom_file(kiomid)%iduld ) iom_file(kiomid)%luld(kiv) = .TRUE.   ! unlimited dimension? 
         END DO
         !---------- Deal with scale_factor and add_offset
         llok = NF90_Inquire_attribute(if90id, ivarid, 'scale_factor') == nf90_noerr
         IF( llok) THEN
            CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'scale_factor', iom_file(kiomid)%scf(kiv)), clinfo)
         ELSE
            iom_file(kiomid)%scf(kiv) = 1.
         END IF
         llok = NF90_Inquire_attribute(if90id, ivarid, 'add_offset') == nf90_noerr
         IF( llok ) THEN
            CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'add_offset', iom_file(kiomid)%ofs(kiv)), clinfo)
         ELSE
            iom_file(kiomid)%ofs(kiv) = 0.
         END IF
         ! return the simension size
         IF( PRESENT(kdimsz) ) THEN 
            IF( i_nvd == SIZE(kdimsz) ) THEN
               kdimsz(:) = iom_file(kiomid)%dimsz(1:i_nvd,kiv)
            ELSE
               WRITE(ctmp1,*) i_nvd, SIZE(kdimsz)
               CALL ctl_stop( TRIM(clinfo), 'error in kdimsz size'//TRIM(ctmp1) )
            ENDIF
         ENDIF
         IF( PRESENT(kndims) )  kndims = iom_file(kiomid)%ndims(kiv)
      ELSE  
         iom_nf90_varid = -1   !   variable not found, return error code: -1
      ENDIF
      !
   END FUNCTION iom_nf90_varid


   SUBROUTINE iom_nf90_g0d( kiomid, kvid, pvar, kstart )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_g0d  ***
      !!
      !! ** Purpose : read a scalar with NF90
      !!-----------------------------------------------------------------------
      INTEGER ,               INTENT(in   )            ::   kiomid   ! Identifier of the file
      INTEGER ,               INTENT(in   )            ::   kvid     ! variable id
      REAL(wp),               INTENT(  out)            ::   pvar     ! read field
      INTEGER , DIMENSION(1), INTENT(in   ), OPTIONAL  ::   kstart   ! start position of the reading in each axis
      !
      CHARACTER(LEN=100)      ::   clinfo   ! info character
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_g0d , file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), pvar, start = kstart), clinfo )
      ! 
   END SUBROUTINE iom_nf90_g0d


   SUBROUTINE iom_nf90_g123d( kiomid, kvid, knbdim, kstart, kcount, kx1, kx2, ky1, ky2,   &
         &                    pv_r1d, pv_r2d, pv_r3d )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_g123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable with NF90
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid    ! iom identifier of the file
      INTEGER                    , INTENT(in   )           ::   kvid      ! Name of the variable
      INTEGER                    , INTENT(in   )           ::   knbdim    ! number of dimensions of the variable
      INTEGER , DIMENSION(:)     , INTENT(in   )           ::   kstart    ! start position of the reading in each axis 
      INTEGER , DIMENSION(:)     , INTENT(in   )           ::   kcount    ! number of points to be read in each axis
      INTEGER ,                    INTENT(in   )           ::   kx1, kx2, ky1, ky2   ! subdomain indexes
      REAL(wp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pv_r1d    ! read field (1D case)
      REAL(wp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pv_r2d    ! read field (2D case)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pv_r3d    ! read field (3D case)
      !
      CHARACTER(LEN=100) ::   clinfo               ! info character
      INTEGER            ::   if90id               ! nf90 identifier of the opened file
      INTEGER            ::   ivid                 ! nf90 variable id
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_g123d , file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      if90id = iom_file(kiomid)%nfid         ! get back NetCDF file id
      ivid   = iom_file(kiomid)%nvid(kvid)   ! get back NetCDF var id
      !
      IF(     PRESENT(pv_r1d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pv_r1d(:                ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pv_r2d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pv_r2d(kx1:kx2,ky1:ky2  ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pv_r3d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pv_r3d(kx1:kx2,ky1:ky2,:), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ENDIF
      !
   END SUBROUTINE iom_nf90_g123d

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_nf90_getatt
   !!----------------------------------------------------------------------

   SUBROUTINE iom_nf90_giatt( kiomid, cdatt, pv_i0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_giatt  ***
      !!
      !! ** Purpose : read an integer attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      INTEGER         , INTENT(  out) ::   pv_i0d   ! read field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( llok ) THEN
            ! check the variable has the attribute required
            llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
         ELSE
            CALL ctl_warn('iom_nf90_getatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_getatt, file: '//TRIM(iom_file(kiomid)%name)//', giatt: '//TRIM(cdatt)
         CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values=pv_i0d), clinfo)
      ELSE
         CALL ctl_warn('iom_nf90_getatt: no attribute '//cdatt//' found')
         pv_i0d = -999
      ENDIF
      !
   END SUBROUTINE iom_nf90_giatt

   SUBROUTINE iom_nf90_gratt( kiomid, cdatt, pv_r0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_gratt  ***
      !!
      !! ** Purpose : read a real attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      REAL(wp)        , INTENT(  out) ::   pv_r0d   ! read field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( llok ) THEN
            ! check the variable has the attribute required
            llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
         ELSE
            CALL ctl_warn('iom_nf90_getatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_getatt, file: '//TRIM(iom_file(kiomid)%name)//', gratt: '//TRIM(cdatt)
         CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values=pv_r0d), clinfo)
      ELSE
         CALL ctl_warn('iom_nf90_getatt: no attribute '//cdatt//' found')
         pv_r0d = -999._wp
      ENDIF
      !
   END SUBROUTINE iom_nf90_gratt

   SUBROUTINE iom_nf90_gcatt( kiomid, cdatt, pv_c0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_gcatt  ***
      !!
      !! ** Purpose : read a character attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      CHARACTER(len=*), INTENT(  out) ::   pv_c0d   ! read field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( llok ) THEN
            ! check the variable has the attribute required
            llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
         ELSE
            CALL ctl_warn('iom_nf90_getatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_getatt, file: '//TRIM(iom_file(kiomid)%name)//', gcatt: '//TRIM(cdatt)
         CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values=pv_c0d), clinfo)
      ELSE
         CALL ctl_warn('iom_nf90_getatt: no attribute '//cdatt//' found')
         pv_c0d = '!'
      ENDIF
      !
   END SUBROUTINE iom_nf90_gcatt

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_nf90_putatt
   !!----------------------------------------------------------------------

   SUBROUTINE iom_nf90_piatt( kiomid, cdatt, pv_i0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_piatt  ***
      !!
      !! ** Purpose : write an integer attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      INTEGER         , INTENT(in   ) ::   pv_i0d   ! write field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      LOGICAL                         ::   lenddef  ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      lenddef = .false.
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( .NOT. llok ) THEN
            CALL ctl_warn('iom_nf90_putatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = .true.
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_putatt, file: '//TRIM(iom_file(kiomid)%name)//', piatt: '//TRIM(cdatt)
         IF( iom_file(kiomid)%irec /= -1 ) THEN   
            ! trick: irec used to know if the file is in define mode or not
            ! if it is not then temporarily put it into define mode
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)
            lenddef = .true.
         ENDIF
         !
         CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values=pv_i0d), clinfo)
         !
         IF( lenddef ) THEN   
            ! file was in data mode on entry; put it back in that mode
            CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)
         ENDIF
      ELSE
         CALL ctl_warn('iom_nf90_putatt: no attribute '//cdatt//' written')
      ENDIF
      !
   END SUBROUTINE iom_nf90_piatt

   SUBROUTINE iom_nf90_pratt( kiomid, cdatt, pv_r0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_pratt  ***
      !!
      !! ** Purpose : write a real attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      REAL(wp)        , INTENT(in   ) ::   pv_r0d   ! write field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      LOGICAL                         ::   lenddef  ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      lenddef = .false.
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( .NOT. llok ) THEN
            CALL ctl_warn('iom_nf90_putatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = .true.
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_putatt, file: '//TRIM(iom_file(kiomid)%name)//', pratt: '//TRIM(cdatt)
         IF( iom_file(kiomid)%irec /= -1 ) THEN   
            ! trick: irec used to know if the file is in define mode or not
            ! if it is not then temporarily put it into define mode
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)
            lenddef = .true.
         ENDIF
         !
         CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values=pv_r0d), clinfo)
         !
         IF( lenddef ) THEN   
            ! file was in data mode on entry; put it back in that mode
            CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)
         ENDIF
      ELSE
         CALL ctl_warn('iom_nf90_putatt: no attribute '//cdatt//' written')
      ENDIF
      !
   END SUBROUTINE iom_nf90_pratt

   SUBROUTINE iom_nf90_pcatt( kiomid, cdatt, pv_c0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_pcatt  ***
      !!
      !! ** Purpose : write a character attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      CHARACTER(len=*), INTENT(in   ) ::   pv_c0d   ! write field
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      LOGICAL                         ::   lenddef  ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      lenddef = .false.
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( .NOT. llok ) THEN
            CALL ctl_warn('iom_nf90_putatt: no variable '//cdvar//' found')
         ENDIF
      ELSE
         llok = .true.
         ivarid = NF90_GLOBAL
      ENDIF
!
      IF( llok) THEN
         clinfo = 'iom_nf90_putatt, file: '//TRIM(iom_file(kiomid)%name)//', pcatt: '//TRIM(cdatt)
         IF( iom_file(kiomid)%irec /= -1 ) THEN   
            ! trick: irec used to know if the file is in define mode or not
            ! if it is not then temporarily put it into define mode
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)
            lenddef = .true.
         ENDIF
         !
         CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values=pv_c0d), clinfo)
         !
         IF( lenddef ) THEN   
            ! file was in data mode on entry; put it back in that mode
            CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)
         ENDIF
      ELSE
         CALL ctl_warn('iom_nf90_putatt: no attribute '//cdatt//' written')
      ENDIF
      !
   END SUBROUTINE iom_nf90_pcatt


   SUBROUTINE iom_nf90_gettime( kiomid, kvid, ptime, cdunits, cdcalendar )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_gettime  ***
      !!
      !! ** Purpose : read the time axis kvid in the file with NF90
      !!--------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kiomid     ! file Identifier
      INTEGER                   , INTENT(in   ) ::   kvid       ! variable id
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   ptime      ! the time axis
      CHARACTER(len=*), OPTIONAL, INTENT(  out) ::   cdunits    ! units attribute
      CHARACTER(len=*), OPTIONAL, INTENT(  out) ::   cdcalendar ! calendar attribute
      !
      CHARACTER(LEN=100) ::   clinfo     ! info character
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_gettime, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), ptime(:),   &
            &                           start=(/ 1 /), count=(/ iom_file(kiomid)%dimsz(1, kvid) /)), clinfo)
      IF ( PRESENT(cdunits) ) THEN 
         CALL iom_nf90_check(NF90_GET_ATT(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), "units", &
            &                           values=cdunits), clinfo)
      ENDIF
      IF ( PRESENT(cdcalendar) ) THEN 
         CALL iom_nf90_check(NF90_GET_ATT(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), "calendar", &
            &                           values=cdcalendar), clinfo)
      ENDIF
      !
   END SUBROUTINE iom_nf90_gettime


   SUBROUTINE iom_nf90_rp0123d( kt, kwrite, kiomid, cdvar , kvid  , ktype,   &
         &                               pv_r0d, pv_r1d, pv_r2d, pv_r3d )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_nf90_rstput  ***
      !!
      !! ** Purpose : read the time axis cdvar in the file 
      !!--------------------------------------------------------------------
      INTEGER                     , INTENT(in)           ::   kt       ! ocean time-step
      INTEGER                     , INTENT(in)           ::   kwrite   ! writing time-step
      INTEGER                     , INTENT(in)           ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*)            , INTENT(in)           ::   cdvar    ! variable name
      INTEGER                     , INTENT(in)           ::   kvid     ! variable id
      INTEGER                     , INTENT(in), OPTIONAL ::   ktype    ! variable type (default R8)
      REAL(wp)                    , INTENT(in), OPTIONAL ::   pv_r0d   ! written Od field
      REAL(wp), DIMENSION(      :), INTENT(in), OPTIONAL ::   pv_r1d   ! written 1d field
      REAL(wp), DIMENSION(:, :   ), INTENT(in), OPTIONAL ::   pv_r2d   ! written 2d field
      REAL(wp), DIMENSION(:, :, :), INTENT(in), OPTIONAL ::   pv_r3d   ! written 3d field
      !
      INTEGER               :: idims                ! number of dimension
      INTEGER               :: idvar                ! variable id
      INTEGER               :: jd                   ! dimension loop counter   
      INTEGER               :: ix1, ix2, iy1, iy2   ! subdomain indexes   
      INTEGER, DIMENSION(4) :: idimsz               ! dimensions size  
      INTEGER, DIMENSION(4) :: idimid               ! dimensions id
      CHARACTER(LEN=256)    :: clinfo               ! info character
      CHARACTER(LEN= 12), DIMENSION(4) :: cltmp     ! temporary character
      INTEGER               :: if90id               ! nf90 file identifier
      INTEGER               :: idmy                 ! dummy variable
      INTEGER               :: itype                ! variable type
      INTEGER, DIMENSION(4) :: ichunksz             ! NetCDF4 chunk sizes. Will be computed using
                                                    ! nn_nchunks_[i,j,k,t] namelist parameters
      INTEGER               :: ichunkalg, ishuffle,&
                               ideflate, ideflate_level
                                                    ! NetCDF4 internally fixed parameters
      LOGICAL               :: lchunk               ! logical switch to activate chunking and compression
                                                    ! when appropriate (currently chunking is applied to 4d fields only)
      !---------------------------------------------------------------------
      !
      clinfo = '          iom_nf90_rp0123d, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
      if90id = iom_file(kiomid)%nfid
      !
      ! define dimension variables if it is not already done
      ! ==========================
      IF( iom_file(kiomid)%nvars == 0 ) THEN
         ! are we in define mode?
         IF( iom_file(kiomid)%irec /= -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = -1
         ENDIF
         ! define the dimension variables if it is not already done
         cltmp = (/ 'nav_lon     ', 'nav_lat     ', 'nav_lev     ', 'time_counter' /)
         CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cltmp(1)), NF90_FLOAT , (/ 1, 2 /), iom_file(kiomid)%nvid(1) ), clinfo)
         CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cltmp(2)), NF90_FLOAT , (/ 1, 2 /), iom_file(kiomid)%nvid(2) ), clinfo)
         CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cltmp(3)), NF90_FLOAT , (/ 3    /), iom_file(kiomid)%nvid(3) ), clinfo)
         CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cltmp(4)), NF90_DOUBLE, (/ 4    /), iom_file(kiomid)%nvid(4) ), clinfo)
         ! update informations structure related the dimension variable we just added...
         iom_file(kiomid)%nvars       = 4
         iom_file(kiomid)%luld(1:4)   = (/ .FALSE., .FALSE., .FALSE., .TRUE. /)
         iom_file(kiomid)%cn_var(1:4) = cltmp
         iom_file(kiomid)%ndims(1:4)  = (/ 2, 2, 1, 1 /)  
         ! trick: defined to 0 to say that dimension variables are defined but not yet written
         iom_file(kiomid)%dimsz(1, 1)  = 0   
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' define dimension variables done'
      ENDIF
      ! define the data if it is not already done
      ! ===============
      IF( kvid <= 0 ) THEN
         !
         ! NetCDF4 chunking and compression fixed settings
         ichunkalg = 0
         ishuffle = 1
         ideflate = 1
         ideflate_level = 1
         !
         idvar = iom_file(kiomid)%nvars + 1
         ! are we in define mode?
         IF( iom_file(kiomid)%irec /= -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = -1
         ENDIF
         ! variable definition
         IF(     PRESENT(pv_r0d) ) THEN   ;   idims = 0
         ELSEIF( PRESENT(pv_r1d) ) THEN   ;   idims = 2   ;   idimid(1:idims) = (/    3,4/)
         ELSEIF( PRESENT(pv_r2d) ) THEN   ;   idims = 3   ;   idimid(1:idims) = (/1,2  ,4/)
         ELSEIF( PRESENT(pv_r3d) ) THEN   ;   idims = 4   ;   idimid(1:idims) = (/1,2,3,4/)
         ENDIF
         IF( PRESENT(ktype) ) THEN   ! variable external type
            SELECT CASE (ktype)
            CASE (jp_r8)  ;   itype = NF90_DOUBLE
            CASE (jp_r4)  ;   itype = NF90_FLOAT
            CASE (jp_i4)  ;   itype = NF90_INT
            CASE (jp_i2)  ;   itype = NF90_SHORT
            CASE (jp_i1)  ;   itype = NF90_BYTE
            CASE DEFAULT   ;   CALL ctl_stop( TRIM(clinfo)//' unknown variable type' )
            END SELECT
         ELSE
            itype = NF90_DOUBLE
         ENDIF
         IF( PRESENT(pv_r0d) ) THEN
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cdvar), itype,                    &
                 &                            iom_file(kiomid)%nvid(idvar) ), clinfo)
         ELSE
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cdvar), itype, idimid(1:idims),   &
                 &                            iom_file(kiomid)%nvid(idvar) ), clinfo)
         ENDIF
         lchunk = .false.
         IF( snc4set%luse .AND. idims.eq.4 ) lchunk = .true.
         ! update informations structure related the new variable we want to add...
         iom_file(kiomid)%nvars         = idvar
         iom_file(kiomid)%cn_var(idvar) = TRIM(cdvar)
         iom_file(kiomid)%scf(idvar)    = 1.
         iom_file(kiomid)%ofs(idvar)    = 0.
         iom_file(kiomid)%ndims(idvar)  = idims
         IF( .NOT. PRESENT(pv_r0d) ) THEN   ;   iom_file(kiomid)%luld(idvar) = .TRUE.
         ELSE                               ;   iom_file(kiomid)%luld(idvar) = .FALSE.
         ENDIF
         DO jd = 1, idims
            CALL iom_nf90_check(NF90_INQUIRE_DIMENSION( if90id, idimid(jd), len = iom_file(kiomid)%dimsz(jd,idvar) ), clinfo)
            IF ( lchunk ) ichunksz(jd) = iom_file(kiomid)%dimsz(jd,idvar)
         END DO
         IF ( lchunk ) THEN
            ! Calculate chunk sizes by partitioning each dimension as requested in namnc4 namelist
            ! Disallow very small chunk sizes and prevent chunk sizes larger than each individual dimension
            ichunksz(1) = MIN( ichunksz(1),MAX( (ichunksz(1)-1)/snc4set%ni + 1 ,16 ) ) ! Suggested default nc4set%ni=4
            ichunksz(2) = MIN( ichunksz(2),MAX( (ichunksz(2)-1)/snc4set%nj + 1 ,16 ) ) ! Suggested default nc4set%nj=2
            ichunksz(3) = MIN( ichunksz(3),MAX( (ichunksz(3)-1)/snc4set%nk + 1 , 1 ) ) ! Suggested default nc4set%nk=6
            ichunksz(4) = 1                                                            ! Do not allow chunks to span the
                                                                                       ! unlimited dimension
            CALL iom_nf90_check(SET_NF90_DEF_VAR_CHUNKING(if90id, idvar, ichunkalg, ichunksz), clinfo)
            CALL iom_nf90_check(SET_NF90_DEF_VAR_DEFLATE(if90id, idvar, ishuffle, ideflate, ideflate_level), clinfo)
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' chunked ok. Chunks sizes: ', ichunksz
         ENDIF
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' defined ok'
      ELSE
         idvar = kvid
      ENDIF

      ! time step kwrite : write the variable
      IF( kt == kwrite ) THEN
         ! are we in write mode?
         IF( iom_file(kiomid)%irec == -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = 0
         ENDIF
         ! on what kind of domain must the data be written?
         IF( PRESENT(pv_r2d) .OR. PRESENT(pv_r3d) ) THEN
            idimsz(1:2) = iom_file(kiomid)%dimsz(1:2,idvar)
            IF(     idimsz(1) == (nlei - nldi + 1) .AND. idimsz(2) == (nlej - nldj + 1) ) THEN
               ix1 = nldi   ;   ix2 = nlei   ;   iy1 = nldj   ;   iy2 = nlej
            ELSEIF( idimsz(1) == nlci              .AND. idimsz(2) == nlcj              ) THEN
               ix1 = 1      ;   ix2 = nlci   ;   iy1 = 1      ;   iy2 = nlcj
            ELSEIF( idimsz(1) == jpi               .AND. idimsz(2) == jpj               ) THEN
               ix1 = 1      ;   ix2 = jpi    ;   iy1 = 1      ;   iy2 = jpj
            ELSE 
               CALL ctl_stop( 'iom_nf90_rp0123d: should have been an impossible case...' )
            ENDIF

            ! write dimension variables if it is not already done
            ! =============
            ! trick: is defined to 0 => dimension variable are defined but not yet written
            IF( iom_file(kiomid)%dimsz(1, 1) == 0 ) THEN
               CALL iom_nf90_check(NF90_INQ_VARID( if90id, 'nav_lon'     , idmy ), clinfo)
               CALL iom_nf90_check(NF90_PUT_VAR( if90id, idmy, glamt(ix1:ix2, iy1:iy2) ), clinfo)
               CALL iom_nf90_check(NF90_INQ_VARID( if90id, 'nav_lat'     , idmy ), clinfo)
               CALL iom_nf90_check(NF90_PUT_VAR( if90id, idmy, gphit(ix1:ix2, iy1:iy2) ), clinfo)
               CALL iom_nf90_check(NF90_INQ_VARID( if90id, 'nav_lev'     , idmy ), clinfo)
               CALL iom_nf90_check(NF90_PUT_VAR( if90id, idmy, gdept_1d                ), clinfo)
               ! +++ WRONG VALUE: to be improved but not really useful...
               CALL iom_nf90_check(NF90_INQ_VARID( if90id, 'time_counter', idmy ), clinfo)
               CALL iom_nf90_check(NF90_PUT_VAR( if90id, idmy, kt                      ), clinfo)   
               ! update the values of the variables dimensions size
               CALL iom_nf90_check(NF90_INQUIRE_DIMENSION( if90id, 1, len = iom_file(kiomid)%dimsz(1,1) ), clinfo)
               CALL iom_nf90_check(NF90_INQUIRE_DIMENSION( if90id, 2, len = iom_file(kiomid)%dimsz(2,1) ), clinfo)
               iom_file(kiomid)%dimsz(1:2, 2) = iom_file(kiomid)%dimsz(1:2, 1)
               CALL iom_nf90_check(NF90_INQUIRE_DIMENSION( if90id, 3, len = iom_file(kiomid)%dimsz(1,3) ), clinfo)
               iom_file(kiomid)%dimsz(1  , 4) = 1   ! unlimited dimension
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' write dimension variables done'
            ENDIF
         ENDIF

         ! write the data
         ! =============
         IF(     PRESENT(pv_r0d) ) THEN
            CALL iom_nf90_check(NF90_PUT_VAR( if90id, idvar, pv_r0d                      ), clinfo)
         ELSEIF( PRESENT(pv_r1d) ) THEN
            CALL iom_nf90_check(NF90_PUT_VAR( if90id, idvar, pv_r1d(                  :) ), clinfo)
         ELSEIF( PRESENT(pv_r2d) ) THEN
            CALL iom_nf90_check(NF90_PUT_VAR( if90id, idvar, pv_r2d(ix1:ix2, iy1:iy2   ) ), clinfo)
         ELSEIF( PRESENT(pv_r3d) ) THEN
            CALL iom_nf90_check(NF90_PUT_VAR( if90id, idvar, pv_r3d(ix1:ix2, iy1:iy2, :) ), clinfo)
         ENDIF
         ! add 1 to the size of the temporal dimension (not really useful...)
         IF( iom_file(kiomid)%luld(idvar) )   iom_file(kiomid)%dimsz(iom_file(kiomid)%ndims(idvar), idvar)    &
               &                            = iom_file(kiomid)%dimsz(iom_file(kiomid)%ndims(idvar), idvar) + 1
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' written ok'
      ENDIF
      !     
   END SUBROUTINE iom_nf90_rp0123d


   SUBROUTINE iom_nf90_check( kstatus, cdinfo )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_nf90_check  ***
      !!
      !! ** Purpose :   check nf90 errors
      !!--------------------------------------------------------------------
      INTEGER,          INTENT(in) :: kstatus
      CHARACTER(LEN=*), INTENT(in) :: cdinfo
      !---------------------------------------------------------------------
      IF(kstatus /= nf90_noerr)   CALL ctl_stop( 'iom_nf90_check : '//TRIM(nf90_strerror(kstatus)), TRIM(cdinfo) )
   END SUBROUTINE iom_nf90_check

   !!======================================================================
END MODULE iom_nf90
