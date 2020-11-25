

Update tides code with Nico's version.
++++++++++++++++++++++++++++++++++++++

Add the POLCOMS harmonic analysis to the executable (as in now in)
`<build_opa_orchestra.rst>`_
This requires some changes to the standard ``namelist_cfg``

Add the final (extra) three variables in your namelist_cfg / nambdy_tide ::

  vi $EXP/namelist_cfg
  ...
  !-----------------------------------------------------------------------
  &nambdy_tide   !  tidal forcing at open boundaries
  !-----------------------------------------------------------------------
     filtide      = 'bdydta/SEAsia_bdytide_rotT_'         !  file name root of tidal forcing files
     ln_bdytide_2ddta = .false.                   !
     ln_bdytide_conj  = .false.                    !
                                                                  ! Harmonic analysis with restart from polcom
     ln_harm_ana_compute=.true.          ! Compute the harmonic analysis at the last time step
     ln_harm_ana_store=.true.                 ! Store the harmonic analysis at the last time step for restart
     ln_harmana_read=.false.                    ! Read haronic analyisis from a restart
  /


Edit xml files to output harmonics as amplitudes and phases (e.g.)::

  vi file_def_nemo.xml
  ...
  <file_group id="tidal_harmonics" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1d files -->
    <file id="tidalanalysis.grid_T" name="harmonic_grid_T" description="ocean T grid variables"  enabled=".TRUE.">

      <field field_ref="O1amp"         name="O1amp"       operation="instant" enabled=".TRUE." />
      <field field_ref="O1phase"       name="O1phase"     operation="instant" enabled=".TRUE." />


  vi field_def_nemo-opa.xml
  ...
      <field_group id="Tides_T" grid_ref="grid_T_2D" operation="once" >
      <!-- tidal composante -->
      ...
      <field id="Q1amp"        long_name="Q1 Elevation harmonic Amplitude"                              unit="m"        />
      <field id="Q1phase"      long_name="Q1 Elevation harmonic Phase"                                  unit="degree"   />

*Recall there are elevation, u-vel and v-vel harmonics*. Also editted suffixes
 in velocity fields, adding ``_2D``.


* As before the constituents you want to analyse are set-up in ``nam_diaharm``
 namelist.

* The harmonic analysis is done at the end only as well as the restart dumping
so you can only restart from the last time step so make sure you output the full
 restart at the end. To restart, you just need to turn on the ``ln_harmana_read``
  and to map the files to something like ``restart_harm_ana_*``  as this bit as
   not been developed with a prefix to load the files. You can look at this
    python script if needed:
  ``/work/n01/n01/nibrun/RUNS/SWPacific/SIMU/01_harm_links.py``
