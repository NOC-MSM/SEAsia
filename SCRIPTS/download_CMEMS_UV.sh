#!/bin/bash

:'

***********************
download_CMEMS_UV.sh
***********************

Files, alternatively, could be manually downloaded.
This script will install the `motuclient` python package.
You will need to register with COPERNICUS to obtain an ID and password.

Substitute ``username`` and ``password`` in the motuclient calls below.
This example downloads a whole year of daily eastward and northward
velocities for the subdomain.
'

#::

  # Define some paths
  export CONFIG=SEAsia_CMEMS
  export WORK=/projectsa/accord
  export CMEMS=$WORK/$CONFIG/CMEMS

  # Create python environment
  module load anaconda
  conda create --name motu_env python=3.7 anaconda
  source activate motu_env

  # Install motuclient package
  python -m pip install motuclient

  # Download a whole year of daily files for the subdomain.
  #  username and password must be editted for your copernicus username and password
  year="2017"
  for mon in $(seq -f "%02g" 1 12)
  do
  if [ "$mon" = "01" -o "$mon" = "03" -o "$mon" = "05" -o "$mon" = "07" -o "$mon" = "08" -o "$mon" = "10" -o "$mon" = "12" ]
  then
  echo $mon
   for day in $(seq -f "%02g" 1 31)
    do
    echo $day
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 55 --longitude-max 140 --latitude-min -25 --latitude-max 30 --date-min "${year}-${mon}-${day} 12:00:00" --date-max "${year}-$mon-$day 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --variable vo --out-name "CMEMS_${year}_${mon}_${day}_download_UV.nc" --user username --pwd password
  done
  fi
  if [ "$mon" = "04" -o "$mon" = "06" -o "$mon" = "09" -o "$mon" = "11" ]
  then
  echo $mon
   for day in $(seq -f "%02g" 1 30)
    do
    echo $day
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 55 --longitude-max 140 --latitude-min -25 --latitude-max 30 --date-min "${year}-${mon}-${day} 12:00:00" --date-max "${year}-$mon-$day 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --variable vo --out-name "CMEMS_${year}_${mon}_${day}_download_UV.nc" --user username --pwd password
  done
  fi
  if [ "$mon" = "02" ]
  then
  echo $mon
   for day in $(seq -f "%02g" 1 28)
    do
    echo $day
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 55 --longitude-max 140 --latitude-min -25 --latitude-max 30 --date-min "${year}-${mon}-${day} 12:00:00" --date-max "${year}-$mon-$day 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --variable vo --out-name "CMEMS_${year}_${mon}_${day}_download_UV.nc" --user username --pwd password
  done
  fi
  done

  # Finally download the last day from the previous year and the first day of the
  #  next year to facilitate time interpolating, if required
  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 55 --longitude-max 140 --latitude-min -25 --latitude-max 30 --date-min "2018-01-01 12:00:00" --date-max "2018-01-01 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --variable vo --out-name "CMEMS_2018_01_01_download_UV.nc" --user username --pwd password

  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 55 --longitude-max 140 --latitude-min -25 --latitude-max 30 --date-min "2016-12-31 12:00:00" --date-max "2016-12-31 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable uo --variable vo --out-name "CMEMS_2016_12_31_download_UV.nc" --user username --pwd password
