.. contents:: Table of Contents

*****
Add new SBC forcing to NEMO 
*****

For any reason, one might want to add extra forcing like wave conditions here. It's fairly easy through the 
NEMO **Input Data generic interface**. 

This has been done in NEMO 4, revision 9395 and code can be found under 
**/work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC**
on ARCHER.

Create a module for initialisation and read of the variables
========================

You can use the **OPA/SBC/sbcwave.F90** as an example or the following ::
 
   MODULE sbcNOCwave
     USE iom            ! I/O manager library
     USE in_out_manager ! I/O manager
     USE lib_mpp        ! distribued memory computing library
     USE fldread         ! read input fields

     IMPLICIT NONE
     PRIVATE

     PUBLIC   sbc_noc_wave        ! routine called in sbcmod
     PUBLIC   sbc_noc_wave_init   ! routine called in sbcmod
     LOGICAL, PUBLIC ::   ln_hs

     TYPE(FLD)       , ALLOCATABLE, DIMENSION(:)   :: sf_hs     ! structure of input fields (file informations, fields read) Drag Coefficient
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: hs_wave

     CONTAINS

        SUBROUTINE sbc_noc_wave( kt )
           USE wrk_nemo
           INTEGER, INTENT( in  ) ::  kt       ! ocean time step
           IF ( ln_hs ) THEN
              CALL fld_read( kt, nn_fsbc, sf_hs )       !* read hs from external forcing
              hs_wave(:,:) = sf_hs(1)%fnow(:,:,1) * tmask (:,:,1)
              CALL iom_put( 'Hs_ERA5', hs_wave(:,:) )
           ENDIF
        END SUBROUTINE sbc_noc_wave

        SUBROUTINE sbc_noc_wave_init
           INTEGER                ::  ierror      ! return error code
           INTEGER                ::  ios         ! Local integer output status for namelist read
           CHARACTER(len=100)     ::  cn_dir      ! Root directory for location of drag coefficient files
           TYPE(FLD_N)            ::  sn_hs       ! informations about the fields to be read
           !!---------------------------------------------------------------------
           NAMELIST/namsbc_noc_wave/  cn_dir, ln_hs, sn_hs
           !!---------------------------------------------------------------------
           REWIND( numnam_ref )              ! Namelist namsbc_wave in reference namelist : File for wave model
           READ  ( numnam_ref, namsbc_noc_wave, IOSTAT = ios, ERR = 901 )
      901  IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_noc_wave in reference namelist', lwp )
           REWIND( numnam_cfg )              ! Namelist namsbc_wave in configuration namelist : File for wave model
           READ  ( numnam_cfg, namsbc_noc_wave, IOSTAT = ios, ERR = 902 )
      902  IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_noc_wave in configuration namelist', lwp )
           IF(lwm) WRITE ( numond, namsbc_noc_wave )

           IF ( ln_hs ) THEN
              ALLOCATE( sf_hs(1), STAT=ierror )           !* allocate and fill sf_wave with sn_cdg
              IF( ierror > 0 )    CALL ctl_stop( 'STOP', 'sbc_noc_wave_init: unable to allocate sf_wave structure for hs' )
              ALLOCATE( sf_hs(1)%fnow(jpi,jpj,1)   )
              IF( sn_hs%ln_tint ) ALLOCATE( sf_hs(1)%fdta(jpi,jpj,1,2) )
              CALL fld_fill( sf_hs, (/ sn_hs /), cn_dir, 'sbc_noc_wave_init', 'NOC Wave module ', 'namsbc_noc_wave' )
              ALLOCATE( hs_wave(jpi,jpj) )
              hs_wave(:,:) = 0.0
           ENDIF
           IF(lwp) WRITE(numout,*) "sbc_noc_wave_init"
         END SUBROUTINE sbc_noc_wave_init

   END MODULE sbcNOCwave

Call the function in sbc module
========================

In **OPA/SBC/sbcmod.F90**, the new module needs to be called in the header 

::

   USE sbcNOCwave

Defined a new logical in the namelist (**sbc_init** routine)

::

      NAMELIST/namsbc/ nn_fsbc  ,                                                    &
         &             ln_usr   , ln_flx   , ln_blk       ,                          &
         &             ln_cpl   , ln_mixcpl, nn_components, nn_limflx,               &
         &             nn_ice   , nn_ice_embd,                                       &
         &             ln_traqsr, ln_dm2dc ,                                         &
         &             ln_rnf   , nn_fwb   , ln_ssr   , ln_isf    , ln_apr_dyn ,     &
         &             ln_wave  ,                                                    &
      !--- NB ---
         &             ln_noc_wave  ,                                                &
      !--- END NB ---
         &             ln_cdgw  , ln_sdw   , ln_tauoc  , ln_stcor   ,                &
         &             nn_lsm

You can also add some comments in the section for the *ocean.output*. Then load your module

::

      IF( ln_wave     )   CALL sbc_wave_init              ! surface wave initialisation
      !
      !--- NB ---
      IF( ln_noc_wave )   CALL sbc_noc_wave_init      ! surface wave initialisation
      !--- END NB ---
      END SUBROUTINE sbc_init

Update the field at each time step

::

      IF( ln_wave     )   CALL sbc_wave( kt )            ! surface waves
      !--- NB ---
      IF( ln_noc_wave )   CALL sbc_noc_wave( kt )        ! surface waves
      !--- END NB ---
      
Finally, declare the logical you have created in **OPA/SBC/sbc_oce.F90**

::

   LOGICAL , PUBLIC ::   ln_wave        !: true if some coupling with wave model
   !--- NB
   LOGICAL , PUBLIC ::   ln_noc_wave    !: true if some coupling with wave model
   !--- END NB


Others / to run
===============

To output the parameters, you load, you can add in the module you created something like

::

           IF ( ln_hs ) THEN
              CALL fld_read( kt, nn_fsbc, sf_hs )       !* read hs from external forcing
              hs_wave(:,:) = sf_hs(1)%fnow(:,:,1) * tmask (:,:,1)
              CALL iom_put( 'Hs_ERA5', hs_wave(:,:) )
           ENDIF

and ifor example, you will need to update the **field_def_nemo-opa.xml**

::

      <!-- WAVES -->
      <field_group id="Waves_grid_T" grid_ref="grid_T_2D" operation="once" >
         <field id="Hs_ERA5"    long_name="Hs from ERA5"  unit="m"      />
      </field_group>

and in **file_def_nemo.xml**

::

      <file_group id="1ts" output_freq="1ts"  output_level="10" enabled=".TRUE."> <!-- 1 time step files -->
        <file id="Waves_grid_T" name="Waves_grid_T" description="ocean T grid variables"  enabled=".TRUE.">
          <field field_ref="Hs_ERA5"   name="Hs_ERA5"    operation="instant" enabled=".TRUE." />
        </file>
      </file_group>

The namelist will need to specify a few things

::

  !-----------------------------------------------------------------------
  &namsbc_noc_wave   ! External fields from wave model
  !-----------------------------------------------------------------------
  !              !  file name  ! frequency (hours) ! variable     ! time interp. !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !             !  (if <0  months)  !   name       !   (logical)  !  (T/F)  ! 'monthly' ! filename ! pairing  ! filename      !
     sn_hs       =  'SWPacific_hs0',    1          ,   'swh'      ,     .true.   , .false. , 'yearly'  ,  'weights_ERA5_SWpacific_bicubic.nc' , ''       , ''
  !
     cn_dir  = './'  !  root directory for the location of drag coefficient files
     ln_hs   = .true.
  /

and in **&namsbc**

::

     ln_noc_wave = .true. 






