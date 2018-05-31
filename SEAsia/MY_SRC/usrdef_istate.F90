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
      zdep(:) = (/   0.00,    1.56,    2.67,    3.86,    5.14,    6.54,    8.09,    9.82, &
                &   11.77,   13.99,   16.53,   19.43,   22.76,   26.56,   30.87,   35.74, &
                &   41.18,   47.21,   53.85,   61.11,   69.02,   77.61,   86.93,   97.04, &
                &  108.03,  120.00,  133.08,  147.41,  163.16,  180.55,  199.79,  221.14, &
                &  244.89,  271.36,  300.89,  333.86,  370.69,  411.79,  457.63,  508.64, &
                &  565.29,  628.03,  697.26,  773.37,  856.68,  947.45, 1045.85, 1151.99, &
                & 1265.86, 1387.38, 1516.36, 1652.57, 1795.67, 1945.30, 2101.03, 2262.42, &
                & 2429.03, 2600.38, 2776.04, 2955.57, 3138.56, 3324.64, 3513.45, 3704.66, &
                & 3897.98, 4093.16, 4289.95, 4488.15, 4687.58, 4888.07, 5089.48, 5291.68, &
                & 5494.58, 5698.06, 6400.00 /)
      zsal(:) = (/  33.90,   33.90,   33.90,   33.90,   33.90,   33.91,   33.91,   33.91, &
                &   33.91,   33.92,   33.92,   33.93,   33.93,   33.94,   33.94,   33.95, &
                &   33.95,   33.96,   33.97,   33.98,   33.99,   34.01,   34.03,   34.07, &
                &   34.10,   34.14,   34.18,   34.21,   34.24,   34.27,   34.30,   34.32, &
                &   34.35,   34.37,   34.39,   34.41,   34.44,   34.46,   34.48,   34.50, &
                &   34.53,   34.55,   34.57,   34.59,   34.61,   34.63,   34.64,   34.66, &
                &   34.67,   34.68,   34.69,   34.70,   34.71,   34.71,   34.71,   34.71, &
                &   34.71,   34.70,   34.70,   34.70,   34.70,   34.69,   34.69,   34.69, &
                &   34.69,   34.68,   34.68,   34.68,   34.68,   34.68,   34.67,   34.67, &
                &   34.67,   34.65,   34.65 /)
      ztmp(:) = (/   0.80,    0.80,    0.80,    0.80,    0.80,    0.80,    0.80,    0.79, &
                &    0.79,    0.79,    0.78,    0.78,    0.78,    0.77,    0.76,    0.75, &
                &    0.73,    0.71,    0.68,    0.66,    0.63,    0.63,    0.70,    0.82, &
                &    0.95,    1.10,    1.24,    1.37,    1.50,    1.61,    1.70,    1.77, &
                &    1.82,    1.86,    1.90,    1.94,    1.97,    2.00,    2.02,    2.04, &
                &    2.05,    2.03,    1.98,    1.92,    1.84,    1.75,    1.66,    1.57, &
                &    1.48,    1.38,    1.29,    1.19,    1.10,    1.00,    0.90,    0.80, &
                &    0.69,    0.59,    0.50,    0.40,    0.31,    0.22,    0.13,    0.04, &
                &   -0.05,   -0.13,   -0.20,   -0.27,   -0.32,   -0.36,   -0.41,   -0.48, &
                &   -0.42,   -0.64,   -0.64 /) 
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
