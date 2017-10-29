# quick plot NEMO output
#########################
#
# jp 12 Oct 2017
#
# Quick script to read in nerCDF output from v4 NEMO and see if it is realistic


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt  # plotting
import matplotlib.cm as cm  # colormaps
import sys # Exit command
#####%matplotlib inline

# Set path
dirname = ''#/Users/jeff/Desktop/'
filename = dirname + 'output.abort.nc'

## Load file and variables
f = Dataset(filename)

#load lat and lon
nav_lat = f.variables['nav_lat'][:] # (y,x)
nav_lon = f.variables['nav_lon'][:] # (y,x)
zos = f.variables['sossheig'][:].squeeze() # (time_counter, y, x)
lim = np.max(np.abs(zos[:])) # Find extrema

# Plot data
cmap = cm.Spectral
fig = plt.figure()
plt.rcParams['figure.figsize'] = (10.0, 10.0)

ax = fig.add_subplot(211)
plt.pcolormesh( zos, cmap=cmap )
plt.clim([-lim,lim])
plt.colorbar()
ax = fig.add_subplot(212)
plt.pcolormesh( zos, cmap=cmap )
plt.xlim([159,199])
plt.ylim([51,91])
plt.clim([-lim,lim])
plt.colorbar()
plt.title('SSH')
#plt.xlabel('long'); plt.ylabel('lat')
plt.show()


sys.exit("End script")

###############################################################################
###############################################################################

#var = 'ssh'
#var = 'tide'
var = 'grid'

ra = 6371229 # Radius of the Earth

dirname = '/Users/jeff/Desktop/'
if var == 'ssh':
    filename = dirname + 'Lbay_1d_20000101_20000105_grid_T.nc'
elif var == 'tide':
    filename = dirname + 'Lbay_1d_20000101_20000105_Tides.nc'
elif var == 'grid':
    filename = dirname + 'domain_cfg.nc'

f = Dataset(filename)



## Load in data

if var == 'ssh':
    zos = f.variables['zos'][:] # (time_counter, y, x)
    sos = f.variables['sos'][:] # (time_counter, deptht, y, x)

    # load in time
    time_counter = f.variables['time_counter'][:] # vector
    [nt,ny,nx] = np.shape(zos)

elif var == 'tide':
    m2x = f.variables['M2x'][:] # (y, x)
    m2y = f.variables['M2y'][:] # (y, x)
    e3t = f.variables['e3t'][:] # (t,z,y,x)

elif var == 'grid':
    e1t = f.variables['e1t'][:].squeeze() # (t,y,x)
    e2t = f.variables['e2t'][:].squeeze() # (t,y,x)
    e3t = f.variables['e3t_0'][:] # (t,z,y,x)
    depth = np.sum(e3t,1).squeeze()


#load lat and lon
nav_lat = f.variables['nav_lat'][:] # (y,x)
nav_lon = f.variables['nav_lon'][:] # (y,x)


# Plot data
fig = plt.figure()
plt.rcParams['figure.figsize'] = (10.0, 5.0)

if var == 'ssh':
    ax = fig.add_subplot(121)
    plt.pcolormesh(nav_lon, nav_lat, zos[0,:,:] ); plt.colorbar
    plt.title('SSH')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

    ax = fig.add_subplot(122)
    plt.pcolormesh(nav_lon, nav_lat, sos[0,:,:] ); plt.colorbar
    plt.title('SOS')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

elif var == 'tide':
    ax = fig.add_subplot(131)
    plt.pcolormesh(nav_lon, nav_lat, m2x[:,:] ); plt.colorbar
    plt.title('M2x')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

    ax = fig.add_subplot(132)
    plt.pcolormesh(nav_lon, nav_lat, m2y[:,:] ); plt.colorbar
    plt.title('M2y')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

    ax = fig.add_subplot(132)
    plt.pcolormesh(nav_lon, nav_lat, np.sum(e3t[0,:,:,:],0) ); plt.colorbar
    plt.title('sum e3t')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

elif var == 'grid':
    ax = fig.add_subplot(131)
    plt.pcolormesh(nav_lon, nav_lat, depth ); plt.colorbar
    plt.title('depth')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

    ax = fig.add_subplot(132)
    plt.pcolormesh(nav_lon, nav_lat, e1t ); plt.colorbar
    plt.title('e1t')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()

    ax = fig.add_subplot(133)
    plt.pcolormesh(nav_lon, nav_lat, e2t ); plt.colorbar
    plt.title('e2t')
    plt.xlabel('long'); plt.ylabel('lat')
    plt.colorbar()
plt.show()
