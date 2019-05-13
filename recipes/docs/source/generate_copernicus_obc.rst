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

  GLO_MFC_001_024_coordinates.nc
  GLO_MFC_001_024_mask_bathy.nc
  GLO_MFC_001_024_mdt.nc

Copy them into $INPUTS directory.
jasmin-xfer1.ceda.ac.uk
cd  /gws/nopw/j04/campus/pseudoDropBox/BoBEAS/INPUTS
wget ftp://jpolton:JeffCMEMS2018@nrt.cmems-du.eu/Core/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/global-analysis-forecast-phy-001-024/2016/01/mercatorpsy4v3r1_gl12_mean_20160101_R20160113.nc



The following started as an implementation of James' notes: https://github.com/jdha/PyNEMO/wiki/Accessing-data:-Hints-and-Tips
But then he produced NCML files that took away all the pain. So jump to the next
section

Building NCML files
===================

On the target machine, where the files are to be generated, create a directory
to store the symbolically linked inputs and outputs. The inputs should be stored
in yearly directories for easier processing, with capping files from the last and first
outputs from the surrounding years to enable complete interpolation for the year.
::
