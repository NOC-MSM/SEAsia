
#!/bin/bash

#:'
#
#***********************
#run_EXP_unforcedStrat.sh
#***********************
#'

# Run the experiment from rest with horizontally contant T,S initial condition
# with no met forcing

#::

export CONFIG=NEMOhorizTS
export EXP=$WDIR/RUN_DIRECTORIES/EXP_unforcedStrat

# Choose an appropriate directory for your EXP installation
if [ ! -d "$EXP" ]; then
  mkdir $EXP
  mkdir $EXP/RESTART
fi


cp $NEMO/cfgs/SHARED/*namelist* $EXP/.
cp $NEMO/cfgs/SHARED/*.xml $EXP/.

# Copy in NEMO/XIOS executables
ln -s $NEMO/cfgs/$CONFIG/BLD/bin/nemo.exe $EXP/nemo.exe
ln -s $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe


# namelist_cfg
# nambdy: freeze the boundary conditions. Set to initial state
# ln_usr = true. User defined initial state and surface forcing. Here we use
# a stratification profile that is horizontally contant, and no wind.
# These are compiled into the executable. (In
#  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``).

# Submit job
cd $EXP
sbatch submit.slurm

## Check on queue
# squeue -u $USER
