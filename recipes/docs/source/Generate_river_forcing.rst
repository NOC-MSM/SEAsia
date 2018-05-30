Create river forcing
====================

A matlab (GCOMS?) script for generating a netCDF file of river forcing exists.
For now I have put it in the PyNEMO repo.

To run it needs a bathy_meter.nc and coordinates.nc file, as well as a river database.

Check out PyNEMO / nrct if have already::

  cd $WORK/$USER
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd

  cd nrct/rivers

The idea is that the source data sits in ``DATA`` and the scripts sit in ``mfiles``.

Edit the run script::

  vi mfiles/runoff_interactive_iterative.m

Change the directory settings and the resolution of the configuration.
Open matlab and run this script.

Output::

  ``river_test.nc``

Copy to $INPUTS ()::

  scp river_test.nc $USER@login.archer.ac.uk:$INPUTS/$CONFIG_rivers.nc
  # scp river_test.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/SEAsia_rivers.nc
