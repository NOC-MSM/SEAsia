#%load EA_SSH_anim.py
#
# SEAsia_SSH_anim.py
#
# jp 22 Oct 2017
"""
** Summary: **
Build 2D animation of South East Asia SSH

Load in some existing output SEAsia_1h_20000101_20000130_SSH.nc and plot it and make animation

** origin: ** cotidalchart_NEMO.ipynb

** Author: ** JP 31 Jan 2017

** Changelog: **
* 31 Jan 2017: Start with AMM60_MASSMO_fronts.ipynb
* 24 Mar 2017: Copied from AMM60_SSH_anim.py
* 22 Oct 2017: Modified for SE Asia

** Issues: **
* The longitude label is missing from the final plots
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt  # plotting
import sys
sys.path.append('../jcomp_tools_dev/') # Add the directory with the amm60_data_tools.py file to path. Not used on the NOCL SAN
from AMM60_tools import NEMO_fancy_datestr # NEMO date strings
from mpl_toolkits.mplot3d import Axes3D # 3d plotting
from itertools import product, combinations # plotting 3d frame
import os #, sys  # removing files and building animation
import datetime
##%matplotlib inline

## Check host name and locate file path
import socket
if str('livmaf') in socket.gethostname():
    dirname = '/Volumes/archer/jelt/NEMO/NEMOGCM_jdha/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/XIOS_AMM60_nemo_harmIT2/EXP_harmIT2/OUTPUT/'
#    dirname = '/Users/jeff/Desktop/OneWeekExpiry/tmp/'
elif str('livljobs') in socket.gethostname():
    dirname = '/scratch/jelt/tmp/'
elif str('eslogin') in socket.gethostname():
    dirname = '/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_openbcs/'
else:
    print 'There is no working ELSE option'

# Fix for EAfrica
filename = 'EA_v2_1ts_00010106_00010106_grid_T.nc'
variable = 'sosheig'

# Fix for SEAsia
dirname = '/work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_openbcs/'
filename = 'SEAsia_1h_19791101_19791130_grid_T.nc'
variable = 'zos'
xlim = [75,135]
ylim = [-20.,20.]
levs = np.arange(-2,2+0.1,0.1)

ofile = 'FIGURES/SEAsia_SSH.gif'



f = Dataset(dirname+filename)

## Load in data
var = f.variables[variable][:] # (t, y, x)

# load in time
time_counter = f.variables['time_counter'][:] # vector
time_origin = f.variables['time_counter'].time_origin
time_calendar = f.variables['time_counter'].calendar
time_units = f.variables['time_counter'].units

#load lat and lon
nav_lat = f.variables['nav_lat'][:] # (y,x)
nav_lon = f.variables['nav_lon'][:] # (y,x)


# Process the time data
################################
# Note that the Error flag doesn't work and I haven't actually checked it. What happens with leap years etc...
[time_str, time_datetime, flag_err] = NEMO_fancy_datestr( time_counter, time_origin ) # internal tide data



def make_gif(files,output,delay=100, repeat=True,**kwargs):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    """

    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s'%(delay,loop," ".join(files),output))

def sshplot(X,Y,var,count):
    xlim = [75,135]
    ylim = [-20.,20.]
    #levs = np.arange(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,1)[:-1]
    dc = .2
    levs = np.arange(-1,1+dc,dc)
    #levs = [-1.5,-1,-0.5,0,0.5,1,1.5]

#    levs = [-6,-5,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6]

    cset = ax.contourf(X, Y, var[count,:,:], levels=levs, cmap='Spectral')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('\nlongitude (deg)')
    ax.set_ylabel('\nlatitude (deg)')

    ax.hold(True)
    cs = ax.contour(X, Y, var[count,:,:], levels=levs, colors='k')


    fig.colorbar(cset) # Add colorbar
    adjustFigAspect(fig,aspect=(xlim[1]-xlim[0])/(ylim[1]-ylim[0]))

def sliceplotshade(X,Y,var,count,xlim,ylim,levs):

    cset = ax.pcolormesh(X, Y, var[count,:,:], cmap='Spectral')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('\nlongitude (deg)')
    ax.set_ylabel('\nlatitude (deg)')
    cset.set_clim([levs[0],levs[-1]])
    ax.hold(True)
    #cs = ax.contour(X, Y, var[count,:,:], levels=levs, colors='k')

    fig.colorbar(cset) # Add colorbar
    adjustFigAspect(fig,aspect=(xlim[1]-xlim[0])/(ylim[1]-ylim[0]))


def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)


## Now do the main routine stuff
if __name__ == '__main__':

    count = 1

    [nt,ny,nx] = np.shape(var)

    X_arr = nav_lon
    Y_arr = nav_lat
    var_arr = var[:,:,:] # Note this also has time dimension

    var_arr[var==0] = np.nan
    files = []

    print 'Hardwired timestep in variable dat'
    ## Time timeseries of frames
    for count in range(nt-25,nt):
        print count

        # Plot SST map
        plt.close('all')
        fig = plt.figure(figsize=(10,10))
        ax = fig.gca()
        #sshplot(X_arr,Y_arr,var_arr,count)
        sliceplotshade(X_arr,Y_arr,var,count,xlim,ylim,levs)
        dat = datetime.datetime.strftime(time_datetime[count], '%d %b %Y: %H:%M')
        #dat = 'hrs: '+str(count)
        #dat = 'hour: '+str(count/10)
        plt.title('SEast Asia: Sea Surface Height (m) '+str(dat))
        ax.text(80, 18, str(dat), fontsize=10)
        plt.axis("on")
        #plt.show() # Can not save to file with this command present
        fname = filename.replace('.nc','_'+str(count).zfill(4)+'.png')

        plt.savefig(fname, dpi=100)
        files.append(fname)


    # Make the animated gif and clean up the files
    make_gif(files,ofile,delay=20)

    for f in files:
        os.remove(f)
