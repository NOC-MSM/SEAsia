
from   netCDF4 import Dataset
import numpy as np

coorfile = 'ERA5_MSL_y2013.nc'      ## One ERA forcing file generated previously
maskfile = 'my_era5_LSM.nc'         ## Land Sea Mask from ncks
outfile  = 'ERA5_LSM.nc'            ## Output file

#----------------------------------------------------------------------------

## READ SRC BATHYMETRY
#nc_c  = Dataset( coorfile, 'r' )
#lon_src = nc_c.variables[ 'lon' ][:]
#lat_src = nc_c.variables[ 'lat' ][:]
#nc_c.close()
#print(coorfile, "loaded", lon_src.shape)

## READ SRC BATHYMETRY
nc_src  = Dataset( maskfile, 'r' )
msk_src = nc_src.variables[ 'lsm' ][0,::-1] ## lat to be reverse as it was done in the generation of the forcing files
lat_src = nc_src.variables['latitude'][::-1]
lon_src = nc_src.variables['longitude'][:]
print(maskfile, "loaded", msk_src.shape)
#msk_src[(msk_src==0.)] = -1
#msk_src[(msk_src<1)] = -1
seas = msk_src <  0.5
land = msk_src >= 0.5
msk_src[seas] = -1
msk_src[land] =  1
nlat = msk_src.shape[0] 
nlon = msk_src.shape[1] 

## NETCDF OUTPUT
ncout = Dataset( outfile, 'w', format='NETCDF3_CLASSIC' )
ncout.createDimension( 'nlat', nlat )
ncout.createDimension( 'nlon', nlon )
lon = ncout.createVariable( 'lon', 'f4', ('nlat', 'nlon',), zlib='True' )
lat = ncout.createVariable( 'lat', 'f4', ('nlat', 'nlon',), zlib='True' )
lon[:] = np.broadcast_to(lon_src.reshape(1,nlon), (nlat,nlon))
lat[:] = np.broadcast_to(lat_src.reshape(nlat,1), (nlat,nlon))
bout    = ncout.createVariable( "LSM", 'f4', ('nlat','nlon',), zlib='True', fill_value=-999. )
bout[:] = msk_src
ncout.close()
