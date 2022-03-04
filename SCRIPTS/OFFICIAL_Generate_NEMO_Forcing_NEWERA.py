#
#====================== DOCSTRING ============================
"""
Generate ERA5 atmospheric forcing for NEMO
So far to prduce a year by a year - need to be automated
--------------------------------------------------------------
"""
__author__      = "Nicolas Bruneau"
__copyright__   = "Copyright 2018, NOC"
__email__       = "nibrun@noc.ac.uk"
__date__        = "2018-05"

#====================== USR PARAMS ===========================

Year_init    = 2013                                             ## First year to process          
Year_end     = 2014                                             ## Last one [included]
East         =   0                                              ## East Border
West         =  -10                                             ## West Border
North        =   55                                             ## North Border
South        =   45                                             ## South Border
path_ERA5    = '/projectsa/NEMO/Forcing/ERA5/SURFACE_FORCING/'  ## ROOT PATH OD ERA5 DATA
path_EXTRACT = './BUILD_CFG/SURFACE_FORCING/EXTRACTION/'        ## WHERE TO EXTRACT YOUR REGION
path_FORCING = './BUILD_CFG/SURFACE_FORCING/'                   ## NEMO FORCING TARGET DIR
clean        = False                                            ## Clean extraction (longest bit)
sph_ON       = True                                             ## Compute specific humidity or not

#================== NEMO DOCUMENTATION =======================

"""
See the manual in section SBC for more details on the way data
are used in NEMO
The time variable from the netcdf is not used
"""

#====================== LOAD MODULES =========================

import os, sys, glob
#import subprocess
import numpy as np
import datetime
from   netCDF4 import Dataset, MFDataset
import netcdftime
#import matplotlib.pyplot as plt
#from   matplotlib.mlab import griddata
#import scipy.spatial.qhull as qhull

#====================== VARIABLE DEF =========================

var_path = { "10m_u_component_of_wind" : "u10", \
             "10m_v_component_of_wind" : "v10", \
             "2m_temperature"          : "t2m", \
             "mean_sea_level_pressure" : "msl", \
             "mean_snowfall_rate"      : "msr" , \
             "mean_surface_downward_long_wave_radiation_flux"  : "msdwlwrf", \
             "mean_surface_downward_short_wave_radiation_flux" : "msdwswrf", \
             "mean_total_precipitation_rate" : "mtpr" }

if sph_ON :
   var_path[ "surface_pressure"  ] = 'sp'
   var_path[ "2m_dewpoint_temperature" ] = 'd2m'

#===================== INTERNAL FCTNS ========================

help_tosec = np.vectorize( lambda x : x.total_seconds() )

def Read_NetCDF_Concatenate( fname, KeyVar ) :
    'Read NetCDF file'
    lfiles = sorted( glob.glob( fname ) )
    for iF, myfile in enumerate(lfiles) :
        nc = Dataset( myfile, 'r' )
        ## Get time using the time variable
        Time_Var = nc.variables[ 'time']
        dt = Time_Var[:][1] - Time_Var[:][0]
        Time_H   = np.arange( Time_Var[:][0], Time_Var[:][0]+dt*Time_Var[:].size, dt )
        try :    Tref = netcdftime.utime( Time_Var.units, calendar = Time_Var.calendar )
        except : Tref = netcdftime.utime( Time_Var.units, calendar = "gregorian" )
        Time = Tref.num2date( Time_H )
        ## Get Coordinates
        if iF == 0:
           try :
             Lon = nc.variables[ 'longitude' ][:]
             Lat = nc.variables[ 'latitude' ][:]
             LON, LAT = np.meshgrid( Lon, Lat )
           except :
             LON = nc.variables[ 'lon' ][:]
             LAT = nc.variables[ 'lat' ][:]
        ## Get Variable
        dum = nc.variables[ KeyVar ]
        Var = dum[:]; ind = ( Var == dum._FillValue ); Var[ind] = np.nan
        ## save
        if iF == 0 : out = Var; tout = Time
        else       : out = np.concatenate( [out,Var], axis=0 ); tout = np.concatenate( [tout,Time], axis=0 )
    print( tout[0], tout[-1], tout.shape, out.shape, LON.shape )
    try    : return tout, LON, LAT, out, dum.units, dum.long_name
    except : return tout, LON, LAT, out, dum.units, dum.standard_name


def Read_NetCDF( fname, KeyVar ) : 
    'Read NetCDF file'
    if "*" in fname : nc = MFDataset( fname, 'r' )
    else            : nc =   Dataset( fname, 'r' )
    ## Get time using the time variable
    Time_Var = nc.variables[ 'time']
    dt = Time_Var[:][1] - Time_Var[:][0]
    Time_H   = np.arange( Time_Var[:][0], Time_Var[:][0]+dt*Time_Var[:].size, dt ) 
    try :    Tref = netcdftime.utime( Time_Var.units, calendar = Time_Var.calendar )
    except : Tref = netcdftime.utime( Time_Var.units, calendar = "gregorian" )
    Time = Tref.num2date( Time_H )
    print( "====================++++")
    ## Get Coordinates
    try :
       Lon = nc.variables[ 'longitude' ][:]
       Lat = nc.variables[ 'latitude' ][:]
       LON, LAT = np.meshgrid( Lon, Lat )
    except :
       LON = nc.variables[ 'lon' ][:]
       LAT = nc.variables[ 'lat' ][:]

    ## Get Variable
    dum = nc.variables[ KeyVar ]
    Var = dum[:]
    ind = ( Var == dum._FillValue )
    Var[ind] = np.nan 
    #Var = Var / dum.scale_factor + dum.add_offset
    ind = (np.isnan(Var))
    Var[ind] = -9999999

    print( Time[0], Time[-1], Var.shape, Time.shape, np.sum(ind))
    try    : return Time, LON, LAT, Var, dum.units, dum.long_name
    except : return Time, LON, LAT, Var, dum.units, dum.standard_name

#=================== MANIPULATE NetCDF =======================

def compute_scale_and_offset( Var, n ):
    'http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/'
    Vmin = np.nanmin( Var )
    Vmax = np.nanmax( Var )
    print( "scaleoffset", Vmin, Vmax)
    # stretch/compress data to the available packed range
    scale_factor = (Vmax - Vmin) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add_offset = Vmin + 2 ** (n - 1) * scale_factor
    return scale_factor, add_offset

def Add_Variable( nc, vName, vDim, vVal, long_name=None, units=None, standard_name=None, fill_value=None) :
    "Add a variable with its attributes in a netcdf file"
    if vName not in ['time','lon','lat',] : fprec = 'i'
    else : fprec = 'f8'

    if fill_value != None : nc.createVariable( vName, fprec, vDim, fill_value=fill_value, zlib=True, complevel=5 )
    else                  : nc.createVariable( vName, fprec, vDim, zlib=True, complevel=5 )
    if long_name     != None : nc.variables[ vName ].long_name     = long_name
    if units         != None : nc.variables[ vName ].units         = units
    if standard_name != None : nc.variables[ vName ].standard_name = standard_name
    if vName not in ['time','lon','lat',] :
       sc, off = compute_scale_and_offset( vVal, 16 )
       nc.variables[ vName ].scale_factor = sc
       nc.variables[ vName ].add_offset   = off

    nc.variables[ vName ][:] = vVal   # np.floor((vVal-off)/sc)

def Create_Dimensions( nc, lon_name, nLon, lat_name, nLat ) :
    "Create NetCDF dimensions time, nx, ny"
    nc.createDimension( lon_name , nLon )
    nc.createDimension( lat_name , nLat )
    nc.createDimension( 'time'   , None )

def Create_NetCDF_core( nc, tDim, tRef, tVal, sDim, sVal_lon, sVal_lat ) :
    "Create Lon, Lat and Time variables"
    # WRITE TIME INFO
    tUnit = "days since {0} UTC".format( tRef.strftime( "%Y-%m-%d %H:%M:%S" ) ); tCal  = "standard"
    Tref  = netcdftime.utime( tUnit, calendar = tCal )
    Add_Variable( nc, 'time', ('time'), Tref.date2num( tVal ), \
                  long_name = "time since {0}".format(tUnit) , \
                  units = tUnit )
    nc.variables['time'].calendar  = tCal
    #nc.variables['time'].base_date = np.array( [ tRef.year, tRef.month, tRef.day, tRef.hour ] )

    # WRITE LON INFO
    Add_Variable( nc, 'lon', sDim, sVal_lon, long_name = 'Longitude', \
                  units = 'degree_east', standard_name = 'longitude'  )

    # WRITE L INFOT
    Add_Variable( nc, 'lat', sDim, sVal_lat, long_name = 'Latitude', \
                  units = 'degree_north', standard_name = 'latitude' )

def Create_Attributes( nc ) :
    "Add some info - I do it at the end as I had issue with not properly readable netcdf if not"
    nc.Description = 'ERA5 Atmospheric conditions for AMM15 NEMO3.6'
    nc.Author      = 'Prepare_ERA5.py'
    nc.Created     = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.Conventions = "CF-1.0"
    nc.close()

#======================= EXTRACTION ==========================

def Extract( fin, fout, clean=True ) :
    if clean : os.system( "rm {0}".format( fout ) )
    if not os.path.exists( fout ) :
       command = "ncks -d latitude,{0},{1} -d longitude,{2},{3} {4} {5}".format( np.float(South), np.float(North),\
                                                                                 np.float(West), np.float(East), fin, fout )
       print( command )
       os.system( command )
      
def datetime_range(start, end, delta):
    current = [start, ]
    while current[-1] < end:
        current.append( current[-1]+delta )
    return np.array(current)

#======================= CORE PROGR ==========================

## load NCO
os.system( "module load nco/gcc/4.4.2.ncwa" )
os.system( "mkdir {1} {0}".format( path_EXTRACT, path_FORCING ) )
if West < 0 : West = 360.+West
if East < 0 : East = 360.+East

## Loop over each variable
for dirVar, nameVar in var_path.items() :

    print( "================== {0} - {1} ==================".format( dirVar, nameVar ))

    ##---------- EXTRACT ALL DATA FOR DOMAIN ----------------
    for iY in range( Year_init, Year_end+1 ) :
        ## Files
        finput  = "{0}/{1}/{2}_{1}.nc".format( path_ERA5, dirVar, iY )
        foutput = "{2}/{0}_{1}.nc".format( nameVar, iY, path_EXTRACT )
        ## Extract the subdomain
        Extract( finput, foutput, clean=clean ) 

    ##---------- LOAD FULLL TIME SERIES IN MEMORY -----------
    Time, Lon, Lat, dum, Units, Name = Read_NetCDF_Concatenate( "{1}/{0}_*.nc".format( nameVar, path_EXTRACT ), nameVar )
    print("Time" , Time)

    dt  = Time[1] - Time[0]   ## assume to be constant in time
    dt2 = datetime.timedelta( seconds=dt.total_seconds() / 2. )
    print( "dt", dt, dt2 )

    ##---------- SOME PREPROCESSING -------------------------
    ## Add time step for last hour - copy the last input
    dumA  = np.concatenate( [  dum,  dum[-1][np.newaxis,...] ], axis = 0 )
    TimeA = np.array( Time.tolist() + [Time[-1],]  )

    print( "Time" , Time )
    print( "TimeA", TimeA )
    ## instantaneous field every hour. we center it in mid-time step (00:30) as it
    ## is what NEMO assumes according to documentation
    dumC  = ( dumA[0:-1] + dumA[1::] ) / 2.0   
    TimeC =  TimeA[0:-1] + dt2           ## shift half time step positively due to averaging
    suffix = ''

    print( "TimeC", TimeC )

    ##---------- OUTPUT A FILE PER YEAR ---------------------
    for iY in range( Year_init, Year_end+1 ) :

        print( datetime.datetime( iY  ,1,1 ), datetime.datetime( iY+1,1,1 ) )
        indT = ( np.array(TimeC) >= datetime.datetime( iY  ,1,1,0,0,0 ) ) \
             * ( np.array(TimeC) <  datetime.datetime( iY+1,1,1,0,0,0 ) )
        print( "indT",np.sum(indT))
 
        if nameVar in [ "d2m", "sp" ] :
               Fout = "{2}/forSPH_ERA5_{0}_y{1}.nc".format( nameVar.upper(), iY, path_FORCING )
        else : Fout = "{2}/ERA5_{0}_y{1}.nc".format( nameVar.upper(), iY, path_FORCING )
        nc = Dataset( Fout, 'w', format='NETCDF4_CLASSIC')
        Create_Dimensions ( nc, 'nLon', Lon.shape[1], 'nLat' , Lat.shape[0] )
        Create_NetCDF_core( nc, ('time'), TimeC[indT][0], TimeC[indT], ('nLat', 'nLon'), Lon[::-1,:], Lat[::-1,:] )
        Add_Variable( nc, nameVar.upper(), ( 'time', 'nLat', 'nLon'), dumC[indT,::-1,:], units=Units+suffix, standard_name=Name, fill_value=-999999 )
        Create_Attributes( nc )


##---------- PROCESS SPECIFIC HUMIDITY ----------------------     
## Compute Specific Humidity according to ECMWF documentation

if sph_ON : 

   for iY in range( Year_init, Year_end+1 ) :
       Time, Lon, Lat, d2m, dUnits, dName = Read_NetCDF( "{1}/forSPH_ERA5_D2M_y{0}.nc".format( iY, path_FORCING ), 'D2M' )
       Time, Lon, Lat, sp , dUnits, dName = Read_NetCDF( "{1}/forSPH_ERA5_SP_y{0}.nc" .format( iY, path_FORCING ), 'SP'  )
       esat = 611.21 * np.exp( 17.502 * (d2m-273.16) / (d2m-32.19) )
       dyrvap = 287.0597 / 461.5250
       dVar = dyrvap * esat / ( sp - (1-dyrvap) * esat)
       Units = "1"; Name = "Specific Humidity"

       indT = ( Time >= datetime.datetime( iY  ,1,1 ) ) \
            * ( Time <  datetime.datetime( iY+1,1,1 ) )

       Fout = "{1}/ERA5_SPH_y{0}.nc".format( iY, path_FORCING )
       nc = Dataset( Fout, 'w', format='NETCDF4_CLASSIC')
       Create_Dimensions ( nc, 'nLon', Lon.shape[1], 'nLat' , Lat.shape[0] )
       Create_NetCDF_core( nc, ('time'), Time[indT][0], Time[indT], ('nLat', 'nLon'), Lon[:,:], Lat[:,:] )
       Add_Variable( nc, "SPH", ( 'time', 'nLat', 'nLon'), dVar[indT,:,:], units=Units, standard_name=Name, fill_value=-999999 )
       Create_Attributes( nc )


