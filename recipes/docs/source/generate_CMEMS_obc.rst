Generating open boundary condition files from COPERNICUS data
=============================================================

* parent data source: COPERNICUS
* generation tool: PyNEMO
* parent data on JASMIN.
* Generate multiple years in yearly blocks.

* Parent data: GLOBAL_ANALYSIS_FORECAST_PHY_001_024 (daily variables)
* Access script::

  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2016-01-01 12:00:00" --date-max "2016-01-02 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable thetao --variable so --variable zos --variable uo --variable vo --out-dir <OUTPUT_DIRECTORY> --out-name <OUTPUT_FILENAME> --user <USERNAME> --pwd <PASSWORD>


For the purposes of progress I am going to do the open bcs extraction using
pynemo install on JASMIN ``<install_nrct.rst>_``.

Workflow:
=========

* build parent grid files
* build NCML and namelist.bdy files for each year
* Generate monthly boundary files



Building parent grid files
--------------------------

PyNEMO expects information about the parent horizontal and vertical grids. These
are specified in the namelist.bdy file::

 sn_src_hgr = './mesh_hgr_src.nc'   !  parent /grid/
 sn_src_zgr = './mesh_zgr_src.nc'   !  parent

This would normally all be obtrained from the domain_cfg.nc file, if it was
available.

ftp://nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics

get the statics files::

  wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_coordinates.nc
  wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_mask_bathy.nc
  #wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_mdt.nc

Copy them into $INPUTS directory.
jasmin-xfer1.ceda.ac.uk
cd  /gws/nopw/j04/campus/pseudoDropBox/BoBEAS/INPUTS
wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024/2016/01/mercatorpsy4v3r1_gl12_mean_20160101_R20160113.nc


Try and create a post-hoc bathymetry file for child model (ideally this would be
the source of bathymetry in the domain_cfg.nc file, but it is hard to recover
what Ash did so I will reverse engineer it). First create a variable of the
appropriate name, dimension and type::

  ncks -v bottom_level,nav_lon,nav_lat domain_cfg.nc tmp.nc
  ncrename -h -v bottom_level,Bathymetry tmp.nc tmp2.nc
  ncwa -a t tmp2.nc tmp3.nc
  ncap2 -s 'Bathymetry=float(Bathymetry)' tmp3.nc fake_bathy.nc

In python. Reconstruct bathymetry by summing e3w::

  import netCDF4
  import numpy as np
  dset = netCDF4.Dataset('domain_cfg.nc','r')
  dout = netCDF4.Dataset('fake_bathy.nc','a')
  [ny,nx] = np.shape(dset.variables['nav_lat'][:])

  #e3w = dset.variables['e3w_0'][:].squeeze() # z,y,x
  e3t = dset.variables['e3t_0'][:].squeeze() # z,y,x
  bathymetry = np.zeros((ny,nx))

  bottom_level = dset.variables['bottom_level'][:].squeeze() # y,x
  for i in range(nx):
    for j in range(ny):
      bathymetry[j,i] = np.sum(e3t[0:bottom_level[j,i],j,i],0)
      #bathymetry[j,i] = np.sum(e3w[0:bottom_level[j,i],j,i],0)

  dout.variables['Bathymetry'][:,:] = np.squeeze(bathymetry)

  dset.close()
  dout.close()


When PyNEMO runs it, by default, creates a mask with a 1pt border around the
domain, unless otherwise stated. At present invoking this default action throws
up an error to do with QWidget, which I image is a java issue. We will circumvent
this issue by manually making a 1pt mask.

NB The mask variable takes values (-1 mask, 1 wet/active domain, 0 land). Need to
only mask a single point around the edge since the rimwidth is considered to be
part of the active domain.

PyNEMO looks for the interface between -1 and 1 to generate boundary forcings. Get a
template from domain_cfg.nc and then modify as desired around the boundary.

Maked the entire boundary "land=0", except for the wet bits along the southern boundary
which are "mask=-1"::

  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

Then in python::

  import netCDF4
  import numpy as np
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  [ny,nx] = np.shape(dset.variables['mask'][:])
  for i in range(nx):
    if dset.variables['mask'][1,i] == 1:
      dset.variables['mask'][0,i] = -1
    else:
      dset.variables['mask'][0,i] = 0

  for j in range(ny):
    if dset.variables['mask'][j,1] == 1:
      dset.variables['mask'][j,0] = -1
    else:
      dset.variables['mask'][j,0] = 0

  dset.close()

The mask is activated in the namelist.bdy file with ``ln_mask_file =.true.``


Pynemo expects the parent grid coordinates to come from model output and therefore
have a singleton time dimension...

Preprocess the parent coordinates file so that it has a time dimension [t,y,x].
For now, this is easier than fixing pynemo to accept generalised dimensional
 files ...::

  livljobs6:
  module load nco/gcc/4.4.2.ncwa
  ncap2 -s 'defdim("time",1);time[time]=74875.0;t@long_name="Dummy time"' -O GLO-MFC_001_024_coordinates.nc tmp.nc

Then add the time dimension to latitude and longitude variables, recreating the glamt and gphit::

  ncap2 -O -s 'gphit[time,latitude,longitude]=latitude' tmp.nc tmp2.nc

  ncap2 -O -s 'glamt[time,latitude,longitude]=longitude' tmp2.nc tmp3.nc
  ncrename -d time,t -d latitude,y -d longitude,x tmp3.nc tmp4.nc

  ncap2 -O -s 'glamu[t,y,x]=longitude' tmp4.nc tmp5.nc
  ncap2 -O -s 'glamv[t,y,x]=longitude' tmp5.nc tmp6.nc

  ncap2 -O -s 'gphiu[t,y,x]=latitude' tmp6.nc tmp7.nc
  ncap2 -O -s 'gphiv[t,y,x]=latitude' tmp7.nc GLO-MFC_001_024_coordinates_v2.nc
  rm tmp*.nc


Manipulate mask variables to recreate expected template::


  livljobs6:
  module load nco/gcc/4.4.2.ncwa
  rm tmp*.nc
  ncap2 -s 'defdim("time",1);time[time]=74875.0;t@long_name="Dummy time"' -O GLO-MFC_001_024_mask_bathy.nc tmp.nc

Then add the time dimension to latitude and longitude variables, recreating the glamt and gphit::

  ncap2 -O -s 'tmask[time,depth,latitude,longitude]=mask' tmp.nc tmp2.nc
  ncap2 -O -s 'umask[time,depth,latitude,longitude]=mask' tmp2.nc tmp3.nc
  ncap2 -O -s 'vmask[time,depth,latitude,longitude]=mask' tmp3.nc tmp4.nc
  ncrename -d time,t -d latitude,y -d longitude,x tmp4.nc GLO-MFC_001_024_mask_bathy_v2.nc
  rm tmp*.nc

Redefining mask to tmask with an NCML file didnt work




Parent data files
=================

Try downloading (locally e.g. Liverpool) monthly files::

  cd /projectsa/accord/BoBEAS/INPUTS
  wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-monthly/2016/mercatorpsy4v3r1_gl12_mean_201603.nc

Try downloading jan, feb, march 2016 and run pynemo for Feb only. Does it work?




Run PyNEMO (on JASMIN)
======================

Generate pynemo namelist files in INPUTS directory where the necessary pynemo
input files sit (you have to put them there)
::

  livljobs6 ~ $
  exec ssh-agent $SHELL
  ssh-add ~/.ssh/id_rsa_jasmin
   Enter passphrase for /login/jelt/.ssh/id_rsa_jasmin:
  Identity added: /login/jelt/.ssh/id_rsa_jasmin (/login/jelt/.ssh/id_rsa_jasmin)
  ssh -A jelt@jasmin-login1.ceda.ac.uk
  ssh -A jelt@jasmin-sci1.ceda.ac.uk


Execute pynemo (jasmin-sci1.ceda.ac.uk). Example::

    cd /gws/nopw/j04/campus/pseudoDropBox/BoBEAS/INPUTS/
    export PATH=/home/users/jelt/anaconda/bin:$PATH
    source activate nrct_env # If required

    export PYTHONPATH=$HOME/anaconda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH

    module load java/1.8.0
    export LD_LIBRARY=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.201.b09-2.el6_10.x86_64/jre/lib/amd64/server/:$LD_LIBRARY_PATH

    pynemo -s namelist_2016.bdy


Run PyNEMO on livljobs6
=======================

JASMIN didn't have enough memory to complete the dev cycle / debugging so I have
rsynced the ``BOBEAS/INPUTS`` directory to Liverpool:
 ``/projectsa/accord/BOBEAS/INPUTS/``::

  cd $INPUTS

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist_2016.bdy


  ---

Notes for 2019: 20 April – 10 May ish

Get the state variable data::

  livljobs4 $
  cd /projectsa/accord/BoBEAS/INPUTS/
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-03-31 12:00:00" --date-max "2019-05-10 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable thetao --out-name CMEMS_2019-03-31_2019-05-10_download_T.nc --user jpolton --pwd JeffPCMEMS2018
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-03-31 12:00:00" --date-max "2019-05-10 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable so --out-name CMEMS_2019-03-31_2019-05-10_download_S.nc --user jpolton --pwd JeffPCMEMS2018
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-03-31 12:00:00" --date-max "2019-05-10 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --out-name CMEMS_2019-03-31_2019-05-10_download_U.nc --user jpolton --pwd JeffPCMEMS2018
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-03-31 12:00:00" --date-max "2019-05-10 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable vo --out-name CMEMS_2019-03-31_2019-05-10_download_V.nc --user jpolton --pwd JeffPCMEMS2018
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 60 --longitude-max 110 --latitude-min 0 --latitude-max 30 --date-min "2019-03-31 12:00:00" --date-max "2019-05-10 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable zos --out-name CMEMS_2019-03-31_2019-05-10_download_Z.nc --user jpolton --pwd JeffPCMEMS2018

Get the static files::

    wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_coordinates.nc
    wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_mask_bathy.nc
    #wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024-statics/GLO-MFC_001_024_mdt.nc

Cut them down::
  module load nco/gcc/4.4.2.ncwa
  ncks -d longitude,60.,110. -d latitude,0.,30. GLO-MFC_001_024_coordinates.nc  cut_coordinates.nc
  ncks -d longitude,60.,110. -d latitude,0.,30. GLO-MFC_001_024_mask_bathy.nc  cut_mask_bathy.nc
  #ncks -d longitude,60.,110. -d latitude,0.,30. GLO-MFC_001_024_mdt.nc  cut_mdt.nc # Don't use

Run pynemo::

  cd $INPUTS

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist_apr2019.bdy

This completes but something fishy has happened as the extreme values are way too large.
There must be a mess up with enveloping bathymetry and fill values polluting the
interpolated values.

Implement a threshold value on the velocities as a temporary fix::

    ncap2 -O -s 'where(vomecrty<-0.5) vomecrty=-0.5' BoBEAS_bdyV_y2019m04.nc BoBEAS_bdyV_y2019m04.nc
    ncap2 -O -s 'where(vomecrty>0.5)  vomecrty=0.5'  BoBEAS_bdyV_y2019m04.nc BoBEAS_bdyV_y2019m04.nc
    ncap2 -O -s 'where(vozocrtx<-0.5) vozocrtx=-0.5' BoBEAS_bdyU_y2019m04.nc BoBEAS_bdyU_y2019m04.nc
    ncap2 -O -s 'where(vozocrtx>0.5)  vozocrtx=0.5'  BoBEAS_bdyU_y2019m04.nc BoBEAS_bdyU_y2019m04.nc
