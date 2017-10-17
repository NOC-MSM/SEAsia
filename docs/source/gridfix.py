# fix the grid spacing valaues
##############################
#
# jp 13 Oct 2017
#
# Quick script to read in netCDF coordinates.nc files
#  Compute the e1* and e2* variables and overwrite


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt  # plotting

ra = 6371229 # Radius of the Earth

dirname = '/Users/jeff/Desktop/temp'
filename = dirname + '/' + 'coordinates.nc'
f = Dataset(filename)

#print f.variables


## Load in data
###############
e1t = f.variables['e1t'][:] # (1, 1, y, x)
e2t = f.variables['e2t'][:] # (1, 1, y, x)
glamt = f.variables['glamt'][:] # (1, 1, y, x)
gphit = f.variables['gphit'][:] # (1, 1, y, x)

#load lat and lon
nav_lat = f.variables['nav_lat'][:] # (y,x)
nav_lon = f.variables['nav_lon'][:] # (y,x)

# Show data shape
print 'e1t', np.shape(e1t)
print 'e2t', np.shape(e2t)
print 'gphit', np.shape(gphit)
print 'glamt', np.shape(glamt)


## Do calculations
##################

# Define dx, dy
dx = (glamt - np.roll( glamt, 1, axis=3))*2*np.pi/180.*ra
dy = (gphit - np.roll( gphit, 1, axis=2))*2*np.pi/180.*ra

# Plot data
############
fig = plt.figure()
plt.rcParams['figure.figsize'] = (10.0, 5.0)

ax = fig.add_subplot(321)
plt.pcolormesh(nav_lon, nav_lat, glamt.squeeze() ); plt.colorbar
plt.title('glamt')
#plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()

ax = fig.add_subplot(322)
plt.pcolormesh(nav_lon, nav_lat, gphit.squeeze() ); plt.colorbar
plt.title('gphit')
#plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()

##
ax = fig.add_subplot(323)
plt.pcolormesh(nav_lon, nav_lat, dx.squeeze() ); plt.colorbar
plt.title('dx')
#plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()

ax = fig.add_subplot(324)
plt.pcolormesh(nav_lon, nav_lat, dy.squeeze() ); plt.colorbar
plt.title('dy')
#plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()


##
ax = fig.add_subplot(325)
plt.pcolormesh(nav_lon, nav_lat, e1t.squeeze() ); plt.colorbar
plt.title('e1t')
plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()

ax = fig.add_subplot(326)
plt.pcolormesh(nav_lon, nav_lat, e2t.squeeze() ); plt.colorbar
plt.title('e2t')
plt.xlabel('long'); plt.ylabel('lat')
plt.colorbar()

plt.show()
