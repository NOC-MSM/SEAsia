#!/bin/bash

#:'
#
#*************************************
#adjust_bathymetry_SEAsia.sh
#*************************************
#
# In the following the bathmetry file is adjusted to set a minimum depth and remove any unexpected negative values.
# You can either use python (as shown here), or ``nco`` commands (commented out below) which are not available
# on ARCHER2, but are commonly available.
#::

  cd $DOMAIN

  # activate the python environment and set python path
  conda activate nemo-reloc
  PYTHONPATH=/work/n01/n01/$USER/miniconda3/envs/pynemo3/lib/python3.7/site-packages:$PYTHONPATH
  
  # process bathymetry
  python adjust_bathymetry_SEAsia.py
  cd $WDIR/SCRIPTS






#  # load nco modules. Modules do not currently exist on ARCHER2 so do elsewhere.
#  module load nco
#
#  cd $DOMAIN
#
#  # Remove weirdness with negative bathymetry and make minimum bathymetry
#  # equal to 10 m (resolve any possible wet-drying problems)
#  ncap2 -s 'where(Bathymetry < 0) Bathymetry=0' bathy_meter.nc tmp1.nc
#  ncap2 -s 'where(Bathymetry < 10 && Bathymetry > 0) Bathymetry=10' tmp1.nc -O bathy_meter.nc
#  rm tmp1.nc
#
#  # Copy it if you want for safe keeping
#  cp bathy_meter.nc bathy_meter_SEAsia.nc
#
#  # Fix bathymetry to deal with instabilities (opening some straights that
#  # have only 2 grid points)
#  ncap2 -s 'Bathymetry(140,464)=200' bathy_meter_SEAsia.nc bathy_meter_SEAsia.nc -O
#  ncap2 -s 'Bathymetry(141,464)=200' bathy_meter_SEAsia.nc bathy_meter_SEAsia.nc -O
#  ncap2 -s 'Bathymetry(145,563)=400' bathy_meter_SEAsia.nc bathy_meter_SEAsia.nc -O
#  ncap2 -s 'Bathymetry(145,564)=400' bathy_meter_SEasia.nc bathy_meter_SEAsia.nc -O
#  ncap2 -s 'Bathymetry(140,467)=80'  bathy_meter_SEAsia.nc bathy_meter_SEAsia.nc -O
#
#  cd $WDIR/SCRIPTS
