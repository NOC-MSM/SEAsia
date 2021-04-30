#!/bin/bash

:'

***********************
make_CMEMS_IC.sh
***********************

Generate Initial conditions for T and S from the Copernicus dataset
(GLOBAL_ANALYSIS_FORECAST_PHY_001_024)

First download data manaually or use the scripts `download_CMEMS_TSZ.sh <../SCRIPTS/download_CMEMS_TSZ.sh>`_
and `download_CMEMS_UV.sh <../SCRIPTS/download_CMEMS_UV.sh>`_ (This process is
not part of the main script). The data is expected in ``$ICS``

The following script generates 3D temperature and salinity initial conditions.
The procedure is to flood fill the land in the parent data and then interpolate
onto the child grid.
'

#::

  ## Flood original fields using SOSIE to avoid any problems with mask and interpolation
  # move to the location of the download data
  cd $ICS
  #load modules
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  # Download temperature and salinity from GLOBAL_ANALYSIS_FORECAST_PHY_001_024
  #  from CMEMS catalogue for 01/01/2017 (the initialisation for our experiment).

  # Alternatively. this file exists on ARCHER at : /work/n01/n01/annkat/EXTRA_TOOLS/CMEMS_2017_01_01_download.nc
  # cp /work/n01/n01/annkat/EXTRA_TOOLS/CMEMS_2017_01_01_download.nc $ICS/.

  #####################################################################
  ###############Flooding/fixing original CMEMS files
  ####################################################################
  #ATTENTION some of these files have land as 10^19 (not nan) and weird fillvalues
  ncatted -a _FillValue,,m,f,0 CMEMS_2017_01_01_download.nc
  ncap2 -O -s 'where(so >100.) so=0' CMEMS_2017_01_01_download.nc CMEMS_2017_01_01_download.nc
  ncap2 -O -s 'where(thetao >100.) thetao=0' CMEMS_2017_01_01_download.nc CMEMS_2017_01_01_download.nc
  ncap2 -O -s 'where(so <-100.) so=0' CMEMS_2017_01_01_download.nc CMEMS_2017_01_01_download.nc
  ncap2 -O -s 'where(thetao <-100.) thetao=0' CMEMS_2017_01_01_download.nc CMEMS_2017_01_01_download.nc

  #flood fields before anything
  #copy namelists
  cp $GITCLONE/IC/initcd_votemper_orig.namelist $ICS/initcd_votemper_orig.namelist
  cp $GITCLONE/IC/initcd_vosaline_orig.namelist $ICS/initcd_vosaline_orig.namelist
  #change namelist (you may want to ensure manually namelist correct)
  #sed -i "s/l_reg_trg  = F/l_reg_trg  = T/g" initcd_votemper_orig.namelist
  #sed -i "s/l_reg_trg  = F/l_reg_trg  = T/g" initcd_vosaline_orig.namelist

  #create mask
  ncks -d time,0,0,1 -v so CMEMS_2017_01_01_download.nc CMEMS_mask.nc
  ncatted -a _FillValue,,d,, CMEMS_mask.nc
  ncap2 -O -s 'where(so<25.) so=0' CMEMS_mask.nc CMEMS_mask.nc
  ncap2 -O -s 'where(so >=25.) so=1' CMEMS_mask.nc CMEMS_mask.nc
  ncrename -v so,mask CMEMS_mask.nc

  # SOSIE tool land flooding
  cp $SCRIPTS/job_sosie_orig_CMEMS.sh $ICS/job_sosie_orig_CMEMS.sh
  #ATTENTION modify job for sossie paths etc.
  qsub -q serial job_sosie_orig_CMEMS.sh
  cp job_sosie_orig_CMEMS.sh job_sosie_orig_CMEMS_sal.sh
  sed -i "s/initcd_votemper_orig.namelist/initcd_vosaline_orig.namelist/g" job_sosie_orig_CMEMS_sal.sh
  qsub -q serial job_sosie_orig_CMEMS_sal.sh

  #ATTENTION!!!!: THE SUBMITTED JOBS WILL TAKE SOME TIME (10 MIN) WAIT FOR THEM TO
  #FINISH BEFORE YOU PROCEED
  ###################################################################################

  ##################################################################################
  #### interpolate CMEMS in model grid
  ##################################################################################

  #load modules
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

  #copy the namelists
  cp $GITCLONE/IC/namelist_reshape_bilin_initcd_votemper $ICS/.
  cp $GITCLONE/IC/namelist_reshape_bilin_initcd_vosaline $ICS/.
  #change the namelist as you want
  #sed -i 's/ORCA0083-N06_19600105d05T_SEAsia_grid_T.nc/thetao_CMEMS_2017_01_01_download.nc-SEAsia_2017.nc/g' namelist_reshape_bilin_initcd_votemper
  #sed -i 's/ORCA0083-N06_19600105d05T_SEAsia_grid_T.nc/so_CMEMS_2017_01_01_download.nc-SEAsia_2017.nc/g' namelist_reshape_bilin_initcd_vosaline

  #interpolation: creates weights and remaps
  ln -s $DOMAIN/domain_cfg_ORCA12.nc $ICS/.
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

  #interpolation: mask
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  ncks -d time_counter,0,0,1 -v vosaline initcd_vosaline.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline<25.) vosaline=0' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline >=25.) vosaline=1' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncrename -v vosaline,mask sosie_initcd_mask.nc

  #interpolation: flood and vertical interpolate using sosie tool
  cp $GITCLONE/IC/initcd_votemper.namelist $ICS/initcd_votemper.namelist
  cp $GITCLONE/IC/initcd_vosaline.namelist $ICS/initcd_vosaline.namelist

  #change the namelist (manually too if you have to check/ensure changes)
  #sed -i 's/idrown   = 5,5,.false./idrown   = 300,0,.false./g' initcd_vosaline.namelist
  #sed -i 's/idrown   = 5,5,.false./idrown   = 300,0,.false./g' initcd_votemper.namelist

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

  #Sossie to flood and vertical interpolate to have the same
  #vertical levels with the model
  cp $GITCLONE/IC/job_sosie_CMEMS.sh $ICS/job_sosie_CMEMS.sh
  #ATTENTION modify job for sossie paths etc.
  qsub -q serial job_sosie_CMEMS.sh
  /home/n01/n01/annkat/local/bin/sosie3.x -f initcd_vosaline.namelist
  #/home/n01/n01/annkat/local/bin/sosie3.x -f initcd_votemper.namelist

  #prepare for z interpolation on the fly that ensures stability
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  #make the gdept as a 4D
  ncap2 -O -s 'gdept_3D[z,y,x]=nav_lev' vosaline_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-SEAsia_2017.nc tmp.nc
  ncap2 -O -s 'gdept_4D[time_counter,z,y,x]=gdept_3D' tmp.nc initcd_depth.nc
  rm tmp.nc

  #create mask that corresponds to the new fields
  ncks -d time_counter,0,0,1 -v vosaline vosaline_CMEMS-GLOBAL_ANALYSIS_FORECAST_PHY_001_024-SEAsia_2017.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline<25.) vosaline=0' initcd_mask.nc initcd_mask.nc
  ncap2 -O -s 'where(vosaline >=25.) vosaline=1' initcd_mask.nc initcd_mask.nc
  ncrename -v vosaline,mask initcd_mask.nc
  ############################################################################################

  cd $WORK
