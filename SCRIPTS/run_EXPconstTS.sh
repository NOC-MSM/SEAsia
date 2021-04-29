
#!/bin/bash

#:'
#
#***********************
#run_EXPconstTS.sh
#***********************
#'

# Run the experiment with contant T,S initial condition with tides-only.

#::

export CONFIG=NEMOconstTS
export EXP=$WDIR/RUN_DIRECTORIES/EXPconstTS

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


## Need namelist_cfg, tscript.slurm

# namelist_cfg
# nambdy: Except for tides, freeze the boundary conditions. Set to initial state
# ln_usr = true. User defined initial state and surface forcing. Here we use
# or homogenous initial conditions.
# with the expression being compiled into the executable. (In
#  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``).

# Submit job
sbatch submit.slurm
