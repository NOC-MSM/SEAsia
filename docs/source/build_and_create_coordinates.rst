Generate new coordinates file
=============================

Generate a ``coordinates.nc`` file from a parent NEMO grid at some resolution.
**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and
parent grid location from ``namelist.input`` and outputs a new files with new
resolution grid elements.

.. note: This uses NEMO v4. I noticed that the namelist was different for v3.6

First we need to figure out the indices for the new domain, from the parent grid.
Move parent grid into INPUTS::

  cp $START_FILES/coordinates_ORCA_R12.nc $INPUTS/.

Inspect this parent coordinates file to define the boundary indices for the new config.

Note, I used FERRET locally::

  $livljobs2$ scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.
  ferret etc
  shade/i=3385:3392/j=2251:2266 NAV_LAT
  shade/i=3385:3392/j=2251:2266 NAV_LON


---

There was a problem with the GRIDGEN create_coordinates.exe tool. So here we
use the AGRIF version. The problem was that the spacing variables were not being
scaled by the zoom factor.

Move to the TOOL directory e.g.::

  export TDIR=/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/TOOLS
  cd $TDIR

Copy in the right ARCH file::

  cp /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL_NOXIOS.fcm ../ARCH/.

Build the NESTING tool::

  ./maketools -n NESTING -m XC_ARCHER_INTEL_NOXIOS -j 6

This makes a number of executables in NESTING::

  cd NESTING

Link in parent coordinates file (from somewhere)::

  ln -s $INPUTS/coordinates_ORCA_R12.nc $TDIR/NESTING/.

Write a namelist file::

  vi namelist.input

  &input_output
      iom_activated = true
  /
  &coarse_grid_files
      parent_coordinate_file = 'coordinates_ORCA_R12.nc'
  /
  &bathymetry
  /
  &nesting
      imin = 3385
      imax = 3392
      jmin = 2251
      jmax = 2266
      rho  = 7
      rhot = 7
      bathy_update = false
  /
  &vertical_grid
  /
  &partial_cells
  /
  &nemo_coarse_grid
  /
  &forcing_files
  /
  &interp
  /
  &restart
  /
  &restart_trc
  /

---

Execute tool::

  ./agrif_create_coordinates.exe

This creates a coordinate file::

 ncdump -h 1_coordinates_ORCA_R12.nc

 netcdf \1_coordinates_ORCA_R12 {
  dimensions:
  	x = 53 ;
  	y = 109 ;
  variables:
  	float nav_lon(y, x) ;
  		nav_lon:units = "degrees_east" ;
  		nav_lon:valid_min = -3.540191f ;
  		nav_lon:valid_max = -2.716762f ;
  		nav_lon:long_name = "Longitude" ;
  	float nav_lat(y, x) ;
  		nav_lat:units = "degrees_north" ;
  		nav_lat:valid_min = 53.11438f ;
  		nav_lat:valid_max = 53.96194f ;
  		nav_lat:long_name = "Latitude" ;
  	double glamt(y, x) ;
  		glamt:missing_value = 1.e+20f ;
  	double glamu(y, x) ;
  		glamu:missing_value = 1.e+20f ;
  	double glamv(y, x) ;
  		glamv:missing_value = 1.e+20f ;
  	double glamf(y, x) ;
  		glamf:missing_value = 1.e+20f ;
  	double gphit(y, x) ;
  		gphit:missing_value = 1.e+20f ;
  	double gphiu(y, x) ;
  		gphiu:missing_value = 1.e+20f ;
  	double gphiv(y, x) ;
  		gphiv:missing_value = 1.e+20f ;
  	double gphif(y, x) ;
  		gphif:missing_value = 1.e+20f ;
  	double e1t(y, x) ;
  		e1t:missing_value = 1.e+20f ;
  	double e1u(y, x) ;
  		e1u:missing_value = 1.e+20f ;
  	double e1v(y, x) ;
  		e1v:missing_value = 1.e+20f ;
  	double e1f(y, x) ;
  		e1f:missing_value = 1.e+20f ;
  	double e2t(y, x) ;
  		e2t:missing_value = 1.e+20f ;
  	double e2u(y, x) ;
  		e2u:missing_value = 1.e+20f ;
  	double e2v(y, x) ;
  		e2v:missing_value = 1.e+20f ;
  	double e2f(y, x) ;
  		e2f:missing_value = 1.e+20f ;
  }

Copy it to the $INPUTS directory::

  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc
