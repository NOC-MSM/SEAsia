# adjust_bathymetry_SEAsia.py
#
# In the following the bathmetry file is adjusted to set a minimum depth and remove any unexpected negative values.
# Execute in a fresh shell as the python requirements might otherwise conflict with previously used, or yet to be used modules
#

import numpy as np
import netCDF4

dset = netCDF4.Dataset('bathy_meter.nc','a')
bathy = dset.variables['Bathymetry'][:]
# flatten the land
bathy[bathy < 0 ] = 0
# set a minimum depth of 10m
bathy[(bathy < 10) & (bathy > 0)] = 10
dset.variables['Bathymetry'] = bathy

# modify finescale bathymetry appropriate to the configuration resolution
dset.variables['Bathymetry'][140,464] = 200 
dset.variables['Bathymetry'][141,464] = 200 
dset.variables['Bathymetry'][145,563] = 400 
dset.variables['Bathymetry'][145,564] = 400 
dset.variables['Bathymetry'][140,467] = 80 
dset.close()
exit()

