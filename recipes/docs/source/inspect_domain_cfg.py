"""
inspect_domain_cfg.py

** Summary: **
Check the contents of NEMO domain_cfg.nc

Aim:
* Does it look right?

** Notebook: **
* source: domain_cfg.ipynb

** Author: ** JP 17 Nov 2017

** Changelog: **
* 17 Nov 2017: start.


"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt  # plotting
import matplotlib.colors as mc   # fancy symetric colours on log scale
import matplotlib.cm as cm   # colormap functionality
import seaborn as sns # cyclic colormap for harmonic phase
import numpy.ma as ma

#%matplotlib inline

##############################################################################
# Set up paths
#
# Check host name and username.
# Depending on machine and modify path tree and flag to load partial data set
import socket
hostname = socket.gethostname()

import getpass
username = getpass.getuser()

if 'livmap' in hostname.lower() and username in ['jeff','jelt']:
    rootdir = '/Users/jeff/Desktop'
    dirname = '/'

elif 'livljobs' in hostname.lower() and username in ['jeff','jelt']:
    rootdir = '/scratch/jelt/tmp'
    dirname = '/'

elif 'archer' in hostname.lower() or 'eslogin' in hostname.lower():
    rootdir = ''
    dirname = '/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg/'

else:
    rootdir = ''
    dirname = '/'
    print 'Note ready for this option {} on machine {}'.format(username, hostname)


##############################################################################

def findJI(lat, lon, lat_grid, lon_grid):
    """
    Simple routine to find the nearest J,I coordinates for given lat lon
    Usage: [J,I] = findJI(49, -12, nav_lat_grid_T, nav_lon_grid_T)
    """
    dist2 = np.square(lat_grid - lat) + np.square(lon_grid - lon)
    [J,I] = np.unravel_index( dist2.argmin(), dist2.shape  )
    return [J,I]


def e3_to_depth(pe3t, pe3w, jpk):
    '''
    funtion e3_to_depth
    Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
    Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp.
    Action  :   pe3t, pe3w : scale factor of t- and w-point (m)
    Useage: [gdept, gdepw] = e3_to_depth(e3t, e3w, nz)
    '''

    pdepw      = np.zeros_like(pe3w)
    pdepw[0,:] = 0.
    pdept      = np.zeros_like(pe3t)
    pdept[0,:] = 0.5 * pe3w[0,:]

    for jk in np.arange(1,jpk,1):
        pdepw[jk,:] = pdepw[jk-1,:] + pe3t[jk-1,:]
        pdept[jk,:] = pdept[jk-1,:] + pe3w[jk  ,:]

    return pdept, pdepw


###################
###################

# Set parameters
nn_sig_lev = 39
print 'nn_sig_lev = {}. This must be hardwired. Obtained from namelist_cfg'.format(nn_sig_lev)


## Load in data
print 'load in data'
f = Dataset(rootdir + dirname + 'domain_cfg.nc')

e3t = f.variables['e3t_0'][:].squeeze() #(time=1, depth, y, x)
e3w = f.variables['e3w_0'][:].squeeze() #(time=1, depth, y, x)
bottom_level = f.variables['bottom_level'][:].squeeze() #(time=1, y, x)
top_level = f.variables['top_level'][:].squeeze() #(time=1, y, x)

#glamt = f.variables['glamt'][:].squeeze() #(time=1, depth, y, x)

#load lat and lon and time
nav_lat = f.variables['nav_lat'][:] # (y,x)
nav_lon = f.variables['nav_lon'][:] # (y,x)

[nz, ny, nx] = np.shape(e3t)

# Define depth coordinates
[gdept, gdepw] = e3_to_depth(e3t, e3w, nz)

# Contruct bathymetry
H = np.reshape([gdept[bottom_level[J,I],J,I] for J in range(ny) for I in range(nx) ],(ny,nx))

# Construct level number in 3D
nlev_3d = np.tile( range(nz),(nx,ny,1)).swapaxes(0,2) # level number

## Get e3t and level number onto (Z,L) sections
var1 = e3t[:,:,:]
var2 = nlev_3d # level number




## Choose transects

lab1 = 'onshelf section'
X1 = [100, 105]
Y1= [-7.5, -5]

lab1 = 'I.Oc to Viet'
X1 = [92, 105]
Y1 = [5.5, 10]

lab1 = 'Malay to S.China Sea'
X1 = [95,115]
Y1 = [-5,10]




## Define section A:
[lon1, lon2] = X1
[lat1, lat2] = Y1

[J1,I1] = findJI(lat1, lon1, nav_lat, nav_lon)
[J2,I2] = findJI(lat2, lon2, nav_lat, nav_lon)

npts = max( np.abs(J2-J1), np.abs(I2-I1))
JJ = [int(jj) for jj in np.rint(np.linspace(J1,J2,num=npts))] # force to nearest integers
II = [int(ii) for ii in np.rint(np.linspace(I1,I2,num=npts))] # force to nearest integers


# Build the section coordinates and variable between the two poits
bath_sec1 = np.transpose([H[JJ[i],II[i]] for i in range(npts)])
dep_sec1 = np.transpose([gdept[:,JJ[i],II[i]] for i in range(npts)])
lat_sec1 = np.tile([nav_lat[JJ[i],II[i]] for i in range(npts)], (nz,1))
lon_sec1 = np.tile([nav_lon[JJ[i],II[i]] for i in range(npts)], (nz,1))
var1_sec1 = np.transpose([var1[:,JJ[i],II[i]]  for i in range(npts) ])
var2_sec1 = np.transpose([var2[:,JJ[i],II[i]]  for i in range(npts) ])




# Plot sections
################
print 'plot sections'

#cmap = cm.jet
cmap = cm.Spectral_r
cmap.set_bad('white',1.)
cmap.set_over('red',1.)

thing = { 'log10(e3t)'   : np.log10(var1_sec1) }
thing = { 'level number' : var2_sec1 }

fig = plt.figure()
plt.rcParams['figure.figsize'] = (15.0, 15.0)

ax = fig.add_subplot(111)
plt.contour(lon_sec1, np.log10(dep_sec1), var2_sec1, levels=[0,nn_sig_lev])

plt.pcolormesh(lon_sec1, np.log10(dep_sec1), np.squeeze(thing.values()), cmap=cmap)
plt.plot(lon_sec1[0,:], np.log10(bath_sec1), label='bathymetry' )
plt.plot(lon_sec1[0,:], np.log10(np.sum(var1_sec1,axis=0)), label='sum e3t')
plt.ylim(np.log10([7000,10]))
#plt.xlim([ min(X1), max(X1)])
plt.title(thing.keys())
plt.ylabel('log10 [depth (m)]')
plt.xlabel('latitude')
#plt.clim([-6,-3])
plt.colorbar()
plt.legend( loc='lower right' )


## Save output
print 'save figure'
fname = rootdir + dirname + lab1.replace(" ", "") + '.png'
plt.savefig(fname)
