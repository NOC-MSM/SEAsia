MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE splines
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      REAL(wp), DIMENSION(75)  ::   zdep, zsal, ztmp
      INTEGER :: ji, jj, jk  ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : analytical definition of initial state '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp
      !

      zdep(:) = (/     0.,    10.,    20.,    30.,    40.,    50.,    75.,   100.,   125.,   &
                &    150.,   175.,   200.,   250.,   300.,   350.,   400.,   500.,   600.,   &
                &    700.,   800.,   900.,  1000.,  1100.,  1200.,  1300.,  1400.,  1500.,   &
                &   1750.,  2000.,  2250.,  2500.,  2750.,  3000.,  3250.,  3500.,  3750.,   &
                &   4000.,  4100.,  4200.,  4300.,  4400.,  4500.,  4550.,  4600.,  4700.,   &
                &   4800.,  4900.,  5000.,  5100.,  5200.,  5300.,  5400.,  5500.,  5600.,   &
                &   5700.,  5800.,  5900.,  6000.,  6100.,  6200.,  6300.,  6400.,  6500.,   &
                &   6600.,  6700.,  6800.,  6900.,  7000.,  7100.,  7200.,  7300.,  7400.,   &
                &   7500.,  7600.,  7700. /)

      zsal(:) = (/ 34.05, 34.05, 34.10, 34.13, 34.25, 34.42, 34.88, 35.08, 35.13,   &
                &  35.08, 35.07, 35.06, 35.06, 35.03, 35.01, 34.99, 34.96, 34.97,   &
                &  34.97, 34.95, 34.92, 34.91, 34.88, 34.87, 34.85, 34.83, 34.82,   &
                &  34.80, 34.77, 34.76, 34.75, 34.74, 34.73, 34.73, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72,   &
                &  34.72, 34.72, 34.72 /)

      ztmp(:) = (/ 28.87, 28.87, 28.87, 28.74, 28.33, 28.01, 25.21, 21.99, 18.51,   &
                &  15.55, 14.39, 13.43, 12.27, 11.48, 11.02, 10.51,  9.58,  8.95,   &
                &   8.35,  7.78,  7.16,  6.52,  5.88,  5.44,  5.02,  4.57,  4.14,   &
                &   3.34,  2.64,  2.31,  2.05,  1.86,  1.69,  1.58,  1.41,  1.23,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,  1.15,   &
                &   1.15,  1.15,  1.15  /)

      !DO jk = 1, jpk             ! constant T & S 
      !   DO jj = 1, jpj
      !      DO ji = 1, jpi
      !         pts(ji,jj,jk,jp_tem) = 1.15_wp * ptmask(ji,jj,jk)
      !         pts(ji,jj,jk,jp_sal) = 34.72_wp * ptmask(ji,jj,jk)
      !      END DO
      !   END DO
      !END DO

      !
    !  DO jk = 1, jpk             ! horizontally uniform T & S profiles
         DO jj = 1, jpj
            DO ji = 1, jpi
               pts(ji,jj,:,jp_tem) = spline3(zdep,ztmp,pdept(ji,jj,:)) * ptmask(ji,jj,:)
               pts(ji,jj,:,jp_sal) = spline3(zdep,zsal,pdept(ji,jj,:)) * ptmask(ji,jj,:) 
            END DO
         END DO
    !  END DO
      !   
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
