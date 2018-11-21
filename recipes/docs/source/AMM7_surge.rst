
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



Copy across domain_cdf.nc
MONSooN::

  cd /projects/jcomp/fred/SURGE/AMM7_INPUTS
  scp amm7_surge_domain_cfg.nc jelt@login.archer.ac.uk:$EXP/.


Copy the tides from a AMM7 run (Note these are TPXO tides)::

  cp /work/n01/n01/nibrun/RUNS/AMM7/TEST/bdy/amm7_bdytide*nc $INPUTS/.
  cp /work/n01/n01/nibrun/RUNS/AMM7/TEST/coordinates.bdy.nc .

Get running with tides only. Then add met from Nico's wiki and :ref:`EAfrica_Surge`
