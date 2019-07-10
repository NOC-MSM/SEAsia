
# Generate namelist files for annual open boundary conditions.
"""

ln -s
Edit the pynemo namelist file to use these grid files::

  sn_src_hgr = './mesh_hgr_src_N06.nc'   !  parent /grid/
  sn_src_zgr = './mesh_zgr_src_N06_renamed.nc'   !  parent
  sn_src_msk = './mask_src_N06.nc'       ! parent


Edit the IO files to use SSHFS links::

  sn_src_dir = 'ORCA0083_N06_1960.ncml'       ! src_files/'
  sn_dst_dir = '/projectsa/accord/SEAsia/START_FILES/OPEN_BOUNDARIES/'

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

