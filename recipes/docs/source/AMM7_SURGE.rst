
.. _AMM7_SURGE:

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

  export CONFIG=AMM7_SURGE
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/CONFIG
  export TDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/TOOLS
  export EXP=$CDIR/$CONFIG/EXP00

execute e.g. ::

  . ~/temporary_path_names_for_NEMO_build_surge

Get the code::

  mkdir $WDIR
  mkdir $INPUTS
  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM dev_r8814_surge_modelling_Nemo4

Or make the paths that would normally be created in the buid method.
 (Here we grab prebuilt files)::

    #mkdir $CDIR/$CONFIG/MY_SRC/
    #mkdir /work/n01/n01/jelt/AMM7_SURGE/dev_r8814_surge_modelling_Nemo4/
    #mkdir /work/n01/n01/jelt/AMM7_SURGE/dev_r8814_surge_modelling_Nemo4/CONFIG
    #mkdir /work/n01/n01/jelt/AMM7_SURGE/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_surge
    #mkdir /work/n01/n01/jelt/AMM7_SURGE/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_surge/EXP00

    ln -s $INPUTS $EXP/bdydta

Put files into MY_SRC / get (git) files from MY_SRC::

  rsync -uvt $WORK/jelt/LBay/START_FILES/dommsk.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt $WORK/jelt/LBay/START_FILES/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/MY_SRC/stpctl.F90 $CDIR/$CONFIG/MY_SRC/.
  #I.e.   rsync -uvt $WORK/$USER/*/git_repo/SEAsia/MY_SRC/stpctl.F90 $CDIR/$CONFIG/MY_SRC/.

Files for improved tides and (last 2 items for) harmonic analysis::

  cd /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC
  #cp bdyini.F90 $CDIR/$CONFIG/MY_SRC/. # Already have this file
  rsync -uvt step.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt step_oce.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt tideini.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt tide_FES14.h90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt sbctide.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt tide_mod.F90 $CDIR/$CONFIG/MY_SRC/.
  #cp bdytides.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt diaharm.F90 $CDIR/$CONFIG/MY_SRC/.
  rsync -uvt diaharm_fast.F90 $CDIR/$CONFIG/MY_SRC/.

Put *xml file into EXP. If git cloning it is already there::

  rsync -uvt /work/n01/n01/nibrun/RUNS/SWPacific/SIMU/field_def_nemo-opa.xml $EXP/.


Add a couple of extra lines into the field_def files. This is a glitch in the surge code,
because it doesn't expect you to not care about the winds::

  vi $EXP/field_def_nemo-opa.xml
  line 338
  <field id="wspd"         long_name="wind speed module"                     standard_name="wind_speed"                                                           unit="m/s"                            />
  <field id="uwnd"         long_name="u component of wind"       unit="m/s"         />
  <field id="vwnd"         long_name="v component of wind"       unit="m/s"        />


Load modules ::

  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

1b) Build XIOS2 @ r1242
======================

Follow instructions at :ref:`build_XIOS2`
(Note the final instruction to link the xios_server.exe may not work if the file structure has not been set
up, leave it, we do it here anyway)


Copy XIOS executable from working domain::

  ln -s /work/n01/n01/$USER/xios-2.0_r1242/bin/xios_server.exe $EXP/xios_server.exe
  #ln -s  /work/n01/n01/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe


1c) Build NEMO
=============

Alreadty got NEMO branch ::

    #cd $WDIR
    #svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM dev_r8814_surge_modelling_Nemo4



Copy files required to build nemo.exe . Or get it from git repo. Or get it here.
(Use the FES ready version).
The compile flags::

  vi $CDIR/$CONFIG/cpp_AMM7_SURGE.fcm
  bld::tool::fppkeys  key_nosignedzero key_diainstant key_mpp_mpi key_iomput  \
                      key_diaharm_fast key_FES14_tides

The compiler options (get from git repo download, or find elsewhere)::

  cp $CDIR/$CONFIG/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.
  #cp $WORK/$USER/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.


Make NEMO ::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10


#Copy compiled surge code from MASSMO experiment::

  # A version that will run with Nico tide modifications
  #ln -s /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/MASSMO5_surge/BLD/bin/nemo.exe $EXP/opa
  # A version that will run with the _Old Tides_
  #ln -s /work/n01/n01/jelt/Solent_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/Solent_surge/BLD/bin/nemo.exe $EXP/opa

Copy across some EXP files (added a PBS -q short to runscript)::

  #cp /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_SURGE/EXP00/* $EXP/.
  cp /work/n01/n01/jelt/MASSMO5_surge/dev_r8814_surge_modelling_Nemo4/CONFIG/MASSMO5_surge/EXP00/runscript $EXP/.



Copy across domain_cfg.nc
MONSooN::

  cd /projects/jcomp/fred/SURGE/AMM7_INPUTS
  scp amm7_surge_domain_cfg.nc jelt@login.archer.ac.uk:$INPUTS/domain_cfg.nc



Generate tidal boundary conditions
==================================


livljobs4:
PyNEMO (at the time of writing was not designed to output only 2D tidal forcing,
so some of the error checking for 3D boundary conditions is not needed but has
to be satisfied. So, get all the necessary files onto this machine.
This contains the grid::

  ls $INPUTS/domain_cfg.nc
  #rsync -uvrt jelt@login.archer.ac.uk:/work/n01/n01/$USER/$CONFIG/INPUTS/domain_cfg.nc .

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

It works!!


Backup to git repo
==================

Not sure about the following. It is definiately a good idea not to try and write
to master, especially if --force is used on a spare checkout; you might
overwrite the the whole repo with null ...
Perhaps the safest way is to create new CONFIG direction from the recipe repo
(from a linux box), push changes and then pull it down as a sparse repo to a
 new branch on ARCHER.
::

  cd $CDIR

  git init .
  git remote add -f origin git@github.com:NOC-MSM/NEMO_cfgs.git
  git config core.sparsecheckout true


Specify the folders I want::

  echo $CONFIG/* >> .git/info/sparse-checkout

Last but not least, update your empty repo with the state from the remote. (And
set the upstream remote master)::

git checkout -b amm7
cd AMM7_SURGE/EXP00
git add namelist_* *xml runscript
git commit -m 'add runtime files'

git push --set-upstream origin amm7

To merge amm7 with master do something like::

  git checkout master
  git checkout other_branch .
  (with period)

This will make full copy of other_branch to current (master): Then make regular commit::

  git add --all
  git commit -m "* copy other_branch to master working tree"

---
