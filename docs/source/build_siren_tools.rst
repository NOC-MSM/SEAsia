Build and run SIREN tools
+++++++++++++++++++++++++

Attempt to build and run SIREN tools. I want to see if they are better, in particular
for generating a domain_cfg.nc file in s-coordinates.

From the SIREN docs::

  Functions/Subroutines
  • program create_meshmask
  • subroutine create__mask (td_nam, jpi, jpj, jpk, ld_domcfg)

  This subroutine compute land/ocean mask arrays at tracer points, horizontal
   velocity points (u & v), vorticity points (f) and barotropic stream function
    points (b).
  • type(tatt) function, dimension(ip_maxatt) create__gloatt (cd_bathy, cd_coord, cd_isfdep, td_namh, td_namz) this function create array of global attributes.

  12.4.1 Detailed Description
  This program creates the NetCDF file(s) which contain(s) all the ocean domain
   informations. It allows to create the domain_cfg.nc file needed to run NEMO,
    or the mesh_mask file(s).

  12.4.2 method
  Bathymetry (and optionally ice shelf draft) is read on input file.
  Horizontal grid-point position and scale factors, and the coriolis factor are
   read in coordinates file or computed. Vertical coordinate is defined, and
   the bathymetry recomputed to fit the vertical grid.
  Finally the masks from the bathymetry are computed.
  All the arrays generated, are writen in one to three file(s) depending on output option.
  Note the file contain depends on the vertical coordinate used (z-coord, partial steps, s-coord)

  12.4.3 how to to create domain_cfg or meshmask file:
  ./SIREN/bin/create_meshmask create_meshmask.nam

Note To create the domain_cfg file, you should put in_msh=0

Open terminal on ARCHER and load paths::

  ssh archer

  . ~/temporary_path_names_for_NEMO_build
  echo $CONFIG

First SIREN. Use my XIOS1 file
(see userid and path in variable ``%XIOS_HOME``). Copy from ARCH *store*::

  cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $CDIR/../ARCH/.
  cd $TDIR

  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n SIREN


---

Create a namelist (Zps)

::

  vi ../SIREN/create_meshmask.nam

  &namlog /
  &namcfg
     cn_varcfg = "../SIREN/cfg/variable.cfg"
     cn_dimcfg = "../SIREN/cfg/dimension.cfg"
  /
  &namin
     cn_bathy = "../DOMAINcfg/bathy_meter.nc"
     cn_coord = "../DOMAINcfg/coordinates.nc"
     in_perio = 0
  /
  &namhgr
    in_mshhgr = 0
  /
  &namzgr
    ln_zps   = .TRUE.
    in_nlevel= 10
  /
  &namdmin
    dn_hmin=10.
  /
  &namzco
  !dn_ppsur =-3958.951371276829
  !dn_ppa0 =103.9530096
  !dn_ppa1    =2.415951269
  dn_ppkth   =15.3510137
  dn_ppacr   =7.000
  dn_ppdzmin = 0.1
  dn_pphmax  = 10.
  ln_dbletanh= .TRUE.
  dn_ppa2    = 100.760928500000
  dn_ppkth2  =  48.029893720000
  dn_ppacr2  = 13.000000000000
  /
  &namzps
    dn_e3zps_min = 0.5
    dn_e3zps_rat = 0.2
  /
  &namsco /
  &namlbc /
  &namwd /
  &namgrd /
  &namout
    in_msh = 0
  /


.. note: The above numbers are not understood. Do not assume they are tuned for this config...

Try it out::

  cd $TDIR/SIREN
  ./create_meshmask.exe create_meshmask.nam

This seems to work OK.
 Try and switch to s-coordinates.

::

  vi ../SIREN/create_domaincfg.nam
  ...

**DOESNT WORK AT THE MO**

Copied variables form AMM60 namelist_cfg and namelist_ref

Created a domain_cfg.nc that looks OK (though the top layer is not at constant depth...)
