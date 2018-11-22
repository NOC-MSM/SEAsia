
.. _AMM7_surge:

*****************************************
Setting up AMM7 surge model
*****************************************

**NOTE** This receipe is a copy from :ref:`EAfrica_Surge` and `Solent_surge`. It
is a bit of a hack job to initially get some data size estimates for a proposal.


1a) Get and build Met Surge Config
=================================

Login to Archer ::

  ssh -l $USER login.archer.ac.uk

Create file 'temporary_path_names_for_NEMO_build_surge' and add the following ::

  export CONFIG=AMM7_surge
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/CONFIG
  export TDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

execute e.g. ::

  . ~/temporary_path_names_for_NEMO_build_surge

Make the paths that would normally be created in the buid method.
 (Here we grab prebuilt files)::

    mkdir /work/n01/n01/jelt/AMM7_surge
    mkdir /work/n01/n01/jelt/AMM7_surge/INPUTS
    mkdir /work/n01/n01/jelt/AMM7_surge/dev_r8814_surge_modelling_Nemo4/
    mkdir /work/n01/n01/jelt/AMM7_surge/dev_r8814_surge_modelling_Nemo4/CONFIG
    mkdir /work/n01/n01/jelt/AMM7_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_surge
    mkdir /work/n01/n01/jelt/AMM7_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_surge/EXP0

    ln -s $INPUTS $EXP/bdydta

Copy XIOS executable from working domain::

  ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe

Copy compiled surge code from MASSMO experiment::

  # A version that will run with Nico tide modifications
  ln -s /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/MASSMO5_surge/BLD/bin/nemo.exe $EXP/opa
  # A version that will run with the _Old Tides_
  #ln -s /work/n01/n01/jelt/Solent_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/Solent_surge/BLD/bin/nemo.exe $EXP/opa

Copy across some EXP files (added a PBS -q short to runscript)::

  cp /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_SURGE/EXP00/* $EXP/.
  cp /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/MASSMO5_surge/EXP00/runscript $EXP/.



Copy across domain_cfg.nc
MONSooN::

  cd /projects/jcomp/fred/SURGE/AMM7_INPUTS
  scp amm7_surge_domain_cfg.nc jelt@login.archer.ac.uk:$EXP/.


Copy the tides from a AMM7 run (Note these are TPXO tides)::

  cp /work/n01/n01/nibrun/RUNS/AMM7/TEST/bdy/amm7_bdytide*nc $INPUTS/.
  cp /work/n01/n01/nibrun/RUNS/AMM7/TEST/coordinates.bdy.nc .

Get running with tides only. Then add met from Nico's wiki and :ref:`EAfrica_Surge`

Something has broken. It worked but do not now. I think I've reverted all minor mods...
It can only be a namlist_cfg change...


Generate tidal boundary conditions
==================================


livljobs4:
PyNEMO (at the time of writing was not designed to output only 2D tidal forcing,
so some of the error checking for 3D boundary conditions is not needed but has
to be satisfied. So, get all the necessary files onto this machine.
This contains the grid::

  cd $INPUTS
  rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/domain_cfg.nc .

This is used to ...::

  #rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.nc .

This contains the bathymetry, though it will be extracted from domain_cfg.nc instead::

  #rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/bathy_meter.nc .

Need to generate 3 more files: A ``namelist.bdy`` which drives PyNEMO and which
has two input files: ``inputs_src.ncml`` which points to the data source and
``inputs_dst.ncml`` which remaps some variable names in the destination files::

  cp ../../Solent/INPUTS/namelist.bdy .


First install PyNEMO `install_nrct`_ if not already done so. Use branch ``Generalise-tide-input``::

  cd /work/$USER/nrct
  git checkout Generalise-tide-input

Copy across some parent mesh files and a mask file (even though they are not
used. This is because this old version of PyNEMO didn't anticipate tide-only usage)::

  cp ../../SEAsia/INPUTS/mesh_?gr_src.nc $INPUTS/.
  cp ../../SEAsia/INPUTS/mask_src.nc $INPUTS/.
  cp ../../SEAsia/INPUTS/inputs_dst.ncml $INPUTS/.
  cp ../../SEAsia/INPUTS/cut_inputs_src.ncml $INPUTS/.


If I don't make a boundary mask then it doesn't work... This can also be done with
the PyNEMO GUI. The mask variable takes values (-1 mask, 1 wet, 0 land). Get a
template from domain_cfg.nc and then modify as desired around the boundary.

For this domain there was an issue with the top right corner being too near the amphidrome
(I think) so I chopped it out here::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython::

  import netCDF4, numpy
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  dset.variables['mask'][0,:]  = -1     # Southern boundary
  dset.variables['mask'][-1,:] = -1    # Northern boundary
  dset.variables['mask'][:,-1] = -1    # Eastern boundary
  dset.variables['mask'][:,0] = -1        # Western boundary
  dset.close()

.. delete
  ny,nx = numpy.shape(dset.variables['mask'][:])

  [x1,y1] = [500, ny]
  [x2,y2] = [nx, 300]
  for i in range(x1,nx):
    for j in range(y2,ny):
      if j*(x2-x1) + i*(y1-y2) -y1*x2+y2*x1 > 0:
        dset.variables['mask'][j,i] = -1
  dset.close()


Make a bathymetry file from envolope bathymetry variable ``hbatt``
 (I think this is OK to do..)::


  module load nco/gcc/4.4.2.ncwa
  rm -f hbatt.nc tmp1.nc tmp2.nc
  ncks -v hbatt, nav_lat, nav_lon domain_cfg.nc tmp1.nc
  ncrename -h -v hbatt,Bathymetry tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc hbatt.nc


FES2014 tidal data is used as the tidal data source. This is clumsily set in
``nemo_bdy_tide3.py`` before pynemo is built, though the following namelist.bdy has redundant
references to TPXO.

Get the INPUTS/namelist.bdy file, either as a checkout::

  cd $INPUTS/../..

  git init .
  git remote add origin git@github.com:NOC-MSM/NEMO_cfgs.git
  git config core.sparsecheckout true
  echo "AMM7_surge/INPUTS/*" >> .git/info/sparse-checkout
  git pull --depth=1 origin master


Generate the boundary conditions with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH

  pynemo -s namelist.bdy


This creates::

  coordinates.bdy.nc
  AMM7_surge_bdytide_rotT_NU2_grid_T.nc
  AMM7_surge_bdytide_rotT_O1_grid_T.nc
  AMM7_surge_bdytide_rotT_P1_grid_T.nc
  AMM7_surge_bdytide_rotT_Q1_grid_T.nc
  AMM7_surge_bdytide_rotT_MTM_grid_T.nc
  AMM7_surge_bdytide_rotT_MU2_grid_T.nc
  AMM7_surge_bdytide_rotT_N2_grid_T.nc
  AMM7_surge_bdytide_rotT_N4_grid_T.nc
  AMM7_surge_bdytide_rotT_R2_grid_T.nc
  AMM7_surge_bdytide_rotT_S1_grid_T.nc
  AMM7_surge_bdytide_rotT_2N2_grid_T.nc
  AMM7_surge_bdytide_rotT_J1_grid_T.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_T.nc
  AMM7_surge_bdytide_rotT_K2_grid_T.nc
  AMM7_surge_bdytide_rotT_K1_grid_T.nc
  AMM7_surge_bdytide_rotT_LA2_grid_T.nc
  AMM7_surge_bdytide_rotT_L2_grid_T.nc
  AMM7_surge_bdytide_rotT_M3_grid_T.nc
  AMM7_surge_bdytide_rotT_M2_grid_T.nc
  AMM7_surge_bdytide_rotT_M6_grid_T.nc
  AMM7_surge_bdytide_rotT_M4_grid_T.nc
  AMM7_surge_bdytide_rotT_MF_grid_T.nc
  AMM7_surge_bdytide_rotT_M8_grid_T.nc
  AMM7_surge_bdytide_rotT_MM_grid_T.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_T.nc
  AMM7_surge_bdytide_rotT_MS4_grid_T.nc
  AMM7_surge_bdytide_rotT_MN4_grid_T.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_T.nc
  AMM7_surge_bdytide_rotT_MSF_grid_T.nc
  AMM7_surge_bdytide_rotT_S4_grid_T.nc
  AMM7_surge_bdytide_rotT_S2_grid_T.nc
  AMM7_surge_bdytide_rotT_T2_grid_T.nc
  AMM7_surge_bdytide_rotT_SSA_grid_T.nc
  AMM7_surge_bdytide_rotT_SA_grid_T.nc
  AMM7_surge_bdytide_rotT_NU2_grid_U.nc
  AMM7_surge_bdytide_rotT_O1_grid_U.nc
  AMM7_surge_bdytide_rotT_P1_grid_U.nc
  AMM7_surge_bdytide_rotT_Q1_grid_U.nc
  AMM7_surge_bdytide_rotT_MTM_grid_U.nc
  AMM7_surge_bdytide_rotT_MU2_grid_U.nc
  AMM7_surge_bdytide_rotT_N2_grid_U.nc
  AMM7_surge_bdytide_rotT_N4_grid_U.nc
  AMM7_surge_bdytide_rotT_R2_grid_U.nc
  AMM7_surge_bdytide_rotT_S1_grid_U.nc
  AMM7_surge_bdytide_rotT_2N2_grid_U.nc
  AMM7_surge_bdytide_rotT_J1_grid_U.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_U.nc
  AMM7_surge_bdytide_rotT_K2_grid_U.nc
  AMM7_surge_bdytide_rotT_K1_grid_U.nc
  AMM7_surge_bdytide_rotT_LA2_grid_U.nc
  AMM7_surge_bdytide_rotT_L2_grid_U.nc
  AMM7_surge_bdytide_rotT_M3_grid_U.nc
  AMM7_surge_bdytide_rotT_M2_grid_U.nc
  AMM7_surge_bdytide_rotT_M6_grid_U.nc
  AMM7_surge_bdytide_rotT_M4_grid_U.nc
  AMM7_surge_bdytide_rotT_MF_grid_U.nc
  AMM7_surge_bdytide_rotT_M8_grid_U.nc
  AMM7_surge_bdytide_rotT_MM_grid_U.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_U.nc
  AMM7_surge_bdytide_rotT_MS4_grid_U.nc
  AMM7_surge_bdytide_rotT_MN4_grid_U.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_U.nc
  AMM7_surge_bdytide_rotT_MSF_grid_U.nc
  AMM7_surge_bdytide_rotT_S4_grid_U.nc
  AMM7_surge_bdytide_rotT_S2_grid_U.nc
  AMM7_surge_bdytide_rotT_T2_grid_U.nc
  AMM7_surge_bdytide_rotT_SSA_grid_U.nc
  AMM7_surge_bdytide_rotT_SA_grid_U.nc
  AMM7_surge_bdytide_rotT_NU2_grid_V.nc
  AMM7_surge_bdytide_rotT_O1_grid_V.nc
  AMM7_surge_bdytide_rotT_P1_grid_V.nc
  AMM7_surge_bdytide_rotT_Q1_grid_V.nc
  AMM7_surge_bdytide_rotT_MTM_grid_V.nc
  AMM7_surge_bdytide_rotT_MU2_grid_V.nc
  AMM7_surge_bdytide_rotT_N2_grid_V.nc
  AMM7_surge_bdytide_rotT_N4_grid_V.nc
  AMM7_surge_bdytide_rotT_R2_grid_V.nc
  AMM7_surge_bdytide_rotT_S1_grid_V.nc
  AMM7_surge_bdytide_rotT_2N2_grid_V.nc
  AMM7_surge_bdytide_rotT_J1_grid_V.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_V.nc
  AMM7_surge_bdytide_rotT_K2_grid_V.nc
  AMM7_surge_bdytide_rotT_K1_grid_V.nc
  AMM7_surge_bdytide_rotT_LA2_grid_V.nc
  AMM7_surge_bdytide_rotT_L2_grid_V.nc
  AMM7_surge_bdytide_rotT_M3_grid_V.nc
  AMM7_surge_bdytide_rotT_M2_grid_V.nc
  AMM7_surge_bdytide_rotT_M6_grid_V.nc
  AMM7_surge_bdytide_rotT_M4_grid_V.nc
  AMM7_surge_bdytide_rotT_MF_grid_V.nc
  AMM7_surge_bdytide_rotT_M8_grid_V.nc
  AMM7_surge_bdytide_rotT_MM_grid_V.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_V.nc
  AMM7_surge_bdytide_rotT_MS4_grid_V.nc
  AMM7_surge_bdytide_rotT_MN4_grid_V.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_V.nc
  AMM7_surge_bdytide_rotT_MSF_grid_V.nc
  AMM7_surge_bdytide_rotT_S4_grid_V.nc
  AMM7_surge_bdytide_rotT_S2_grid_V.nc
  AMM7_surge_bdytide_rotT_T2_grid_V.nc
  AMM7_surge_bdytide_rotT_SSA_grid_V.nc
  AMM7_surge_bdytide_rotT_SA_grid_V.nc





Copy the new files back onto ARCHER::

  livljobs4$
  cd $INPUTS
  #rsync -utv namelist.bdy $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/namelist.bdy
  rsync -utv coordinates.bdy.nc $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/coordinates.bdy.nc
  #for file in $CONFIG*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done
  for file in AMM7_surge_bdytide*nc; do rsync -utv $file $USER@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/$file ; done

git commit namelist.bdy::

  cd $INPUTS
  git add namelist.bdy
  git commit -m 'Add namelist.bdy'


Updates to the namelist_cfg to reflect the new files.
ARCHER:

Does not timestep::

  tail ocean.output

  ...
  dia_25h_init : Output 25 hour mean diagnostics
  ~~~~~~~~~~~~
  Namelist nam_dia25h : set 25h outputs
  Switch for 25h diagnostics (T) or not (F)  ln_dia25h  =  F

  AAAAAAAA
