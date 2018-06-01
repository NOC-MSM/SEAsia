.. contents:: Table of Contents

*****
Change to tides in NEMO
*****

The documentation and final changes are mainly for NEMO 4 (trunk revision 9395 used for ORCHESTRA).
This describes two significant changes:

* Internal treatment of tides
 - Increase the number of harmonic constituents processed by to 34 harmonic (to accomodate FES2014 tides)
 - Switch to Doodson number notation from Schureman.
 - Add in long period tides to the potential forcing.
 - Add in variable Love number.

* Harmonic analysis of tides
 - Switch to POLCOMS method
 - Accomodate restarts in harmonic analysis (namelist contolled)
 - Accomodate 3D harmonic analysis (namelist contolled)


Where to access FES 2014
========================

**FES 2014** (Finit-Element Solution) global tide database can be accessed, after registration, via the **AVISO** portal :

   https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html

The FES2014 tides database includes 3 components:  tide elevations (amplitude and phase), tide currents (u and v) and tide
loading on a 1/16°x1/16° grid including including 34 tidal components : K1, M2, M4 , N2 , O1 , P1 , Q1 , S1, S2 , K2 , 2N2,
E2, J1, L2, La2, M3, M6, M8, Mf, MKS2, Mm, MN4, MS4, MSf, MSqm, Mtm, Mu2, N4, Nu2, R2, S4, Sa, Ssa, T2.


Code changes for FES 2014
=========================

In order to use the FES data (or any other dataset), one needs to provide the information of the new constituents.
The file **NEMOGCM/NEMO/OPA_SRC/SBC/tide.h90** needs to be  updated to contain information for the full set of tidal
constituents used by NEMO. In the new package the **tide.h90** is named **tide_FES14.h90** and is activated through the
cpp flag **key_FES14_tides**.

This table mainly provides the equilibrium tide coefficients (*equitide*), the *Schureman*
coefficients (similar to *Doodson* numbers but different) and the formulation needed. As an example, the Schureman vs
Doodson numbers for M2 ::
   M2 : Doodson   : 2  0 0 0 0 0
        Schureman : 2 -2 2 0 0 0

The following document (see page 189 for some proposed values) seems to match implemented code :
   https://docs.lib.noaa.gov/rescue/cgs_specpubs/QB275U35no981924.pdf

The equitide parameter seems to be the equilibrium tide amplitude corrected with the C_n^m coefficient (see PDF p66). As
I couldn't figure out where the equitide component were coming from, I recomputed all of them using the same method
for consistency. So all equitide were computed using the last epocs from Cartwright and Tayer (1971) multiplied by the
corresponding coefficient of their Table 2 and using the equation 12. The equitide is only used for the tidal potential.
And the tidal potential only account for lond-period, diurnal and semi-diurnal effects.

As an example in their Table 4c (p66), M2 (Doodson 200000) has an amplitude of around 0.63186 m
Table 2, give us a correction of m = 2, n = 2 (semi-diurnal) ::
   0.63186 * 3 * sqrt( 5 / 96 / pi ) = 0.24407

Very close to the one defined originally in NEMO (which was **0.242297**). So to correct (to match what
is implemented in **NEMOGCM/NEMO/OPA_SRC/SBC/sbctide.F90** - take care CT71 uses co-latitude) ::
   long wave : Amplitude from CT71 * [ -1   * sqrt( 5 /  4 / pi ) ]
   diurnal   : Amplitude from CT71 * [ -3/2 * sqrt( 5 / 24 / pi ) ]
   semi-diur : Amplitude from CT71 * [  3   * sqrt( 5 / 96 / pi ) ]

ATTENTION: convention seems to be to have a positive coefficient and a 180 shift to represent negative
value. To be confirmed though.

**nutide** is a flag integer [0 - longwave; 1 - diurnal; 2 - semi-diurnal; o/w
 - no tide potential] used to trigger the formulation of the tide potential.
 See **NEMOGCM/NEMO/OPA_SRC/SBC/sbctide.F90** in function tide_init_potential.


Tidal potential due to long-period and love number
==================================================

The original code of NEMO was not accounting for tidal potential due to long-period. This has been added in
**NEMOGCM/NEMO/OPA_SRC/SBC/sbctide.F90**::
     IF    ( Wave(ntide(jk))%nutide == 1 )  THEN  ;  zcs = zcons * SIN( 2._wp*zlat )
     ELSEIF( Wave(ntide(jk))%nutide == 2 )  THEN  ;  zcs = zcons * COS( zlat )**2
     !--- NB 11/2017
     ! Add tide potential for long period tides
     ELSEIF( Wave(ntide(jk))%nutide == 0 )  THEN  ;  zcs = zcons * (0.5_wp-1.5_wp*SIN(zlat)**2._wp)
     !--- END NB
     ELSE                                         ;  zcs = 0._wp
     ENDIF

Further as the maximum number of constituents is hard coded, other routines were changed to extend this number:
    * NEMOGCM/NEMO/OPA_SRC/SBC/tide_mod.F90
    * NEMOGCM/NEMO/OPA_SRC/BDY/bdytides.F90

The **Love number** was fixed to 0.7. It is now set-up in the namelist **namelist nam_tide** under **dn_love_number**.


Restart, 3D and fast harmonic analysis
======================================

Two versions of tidal harmonic analysis co-exist in NEMO, the original one :
    * **NEMOGCM/NEMO/OPA_SRC/SBC/diaharm.F90**

and a second one (developed by Enda from POLCOMS) which allows to output the harmonic analysis to a restart file :
    * **NEMOGCM/NEMO/OPA_SRC/SBC/diaharm_fast.F90**

Here we modified the second one in order to be able to select the constituents we want as the version
analysis the full set of input constituents leading to *crazy* results if the period simulated is not
appropriate (long enough). However we want to keep the possibility to have more constituents and long
period constituents in our input harmonics for more accurate simulations and only do a tidal analysis
of a few selcted constituents (in addition it fastens the computation).

So now the code allow you to select which constituents you want to output. They need to be included in
the input one. (we can imagine in the future extending this to any constituents could be outputed).

To activate this harmonic analysis, you need to compile with the cpp flag **key_diaharm_fast** (This replaces either
**key_diaharm** or **key_diaharm_ana**). Then a namelist need
to be added::
     !-----------------------------------------------------------------------
     &nam_diaharm_fast   !   Harmonic analysis of tidal constituents               ("key_diaharm_fast")
     !-----------------------------------------------------------------------
     ln_diaharm_store = .true.
     ln_diaharm_compute = .true.
     ln_diaharm_read_restart = .true.
     ln_ana_ssh   = .true.
     ln_ana_uvbar = .false.
     ln_ana_bfric = .true.
     ln_ana_rho  = .false.
     ln_ana_uv3d = .true.
     ln_ana_w3d  = .false.
     tname(1) = 'O1',
     tname(2) = 'M2',
     /

The code is able to analyse ssh (**ln_ana_ssh**), barotropic currents (**ln_ana_uvbar**), bottom friction (**ln_ana_bfric**),
and 3D components: density (**ln_ana_rho**), horizontal velocities (**ln_ana_uv3d**) and vertical velocity (**ln_ana_w3d**).
To activate it and write restart, **ln_diaharm_store** needs to be on. To output the results, **ln_diaharm_compute** and to
restart from a previous run with harmonic analysis, **ln_diaharm_read_restart**.

The code is been modified to allow flexibility and to be able to add new analysis in a straightforward way.

Miscellaneous / Important
=========================

You need to update your xml files to include each components !!!


Examples
========

You can find the source code in ARCHER under **/work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC** using :
  - bdyini.F90
  - diaharm.F90
  - diaharm_fast.F90
  - sbctide.F90
  - step.F90
  - step_oce.F90
  - tide_FES14.h90
  - tideini.F90
  - tide_mod.F90

Compilation flags
  - key_FES14_tides
  - key_diaharm_fast

Namelists and xml file examples can be found in **/work/n01/n01/nibrun/RUNS/SWPacific/SIMU** with :
  - namelist_cfg
  - namelist_ref
  - file_def_nemo.xml
  - field_def_nemo-opa.xml
