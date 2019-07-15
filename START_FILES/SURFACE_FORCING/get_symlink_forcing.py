
# get_symlink_forcing.py
"""
Generate Symbolic links to the DFS5.2 forcing fields.

Note that NEMO can use the global forcing files and internally cut out a
subdomain using predefined weights files.

Here symbolic links are generated to put the global forcing files into the
configurations forcing directory path.

JP: July 2019

"""

import numpy as np
import calendar # testing for leap year
import os


DST_DIR = '/work/n01/n01/jelt/SEAsia/INPUTS/SBC/'
SRC_DIR = '/work/n01/n01/acc/DFS/DFS5.2/'

START_YEAR = 1960
END_YEAR = 2015


for year in np.arange(START_YEAR, 1978+1):
    print year
    #print 'Processing year: {}'.format(str(year))

    for var in ['t2', 'u10', 'v10', 'q2']:
        file = 'drowned_'+var+'_DFS5.2_y'+str(year)+'.nc'
	#os.system('rm '+DST_DIR+file)
        os.system('ln -s '+SRC_DIR+str(year)+'/'+file+' '+DST_DIR+file)

    for var in ['precip', 'radlw', 'radsw', 'snow']:
        ofile = 'drowned_'+var+'_DFS5.2_y'+str(year)+'.nc'
        if calendar.isleap(year):
            file = 'drowned_'+var+'_DFS5.2_y1958-1978.leapYears_filtered.nc'
        else:
            file = 'drowned_'+var+'_DFS5.2_y1958-1978_filtered.nc'
	#os.system('rm '+DST_DIR+ofile)
        os.system('ln -s '+SRC_DIR+'1958-1978/'+file+' '+DST_DIR+ofile)


for year in np.arange(1979, END_YEAR+1):
    #print('Processing year: {}'.format(year))
    print year

    for var in ['t2', 'u10', 'v10', 'q2']+['precip', 'radlw', 'radsw', 'snow']:

        file = 'drowned_'+var+'_DFS5.2_y'+str(year)+'.nc'
        os.system('ln -s '+SRC_DIR+str(year)+'/'+file+' '+DST_DIR+file)
