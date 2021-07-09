
#!/bin/bash

#:'
#
#***********************
#run_EXP_unforced.sh
#***********************
#'

# Run the experiment from rest with contant T,S initial condition
# with no met forcing

#::

export CONFIG=NEMOconstTS
export EXP=$WDIR/RUN_DIRECTORIES/EXP_unforced

# Choose an appropriate directory for your EXP installation
if [ ! -d "$EXP/RESTART" ]; then
  #mkdir $EXP
  mkdir $EXP/RESTART
fi

rsync -av --ignore-existing $NEMO/cfgs/SHARED/*namelist* $EXP/. # only get the files not already in the repo.
rsync -av --ignore-existing $NEMO/cfgs/SHARED/*.xml $EXP/. 

# Copy in NEMO/XIOS executables
ln -s $NEMO/cfgs/$CONFIG/BLD/bin/nemo.exe $EXP/nemo.exe
ln -s $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe

# Link in domain_cfg file
ln -s $DOMAIN/domain_cfg_SEVERN.nc $EXP/domain_cfg.nc

# Link in boundary files (just coordinates.bdy.nc)
ln -s $WDIR/INPUTS/OBC/coordinates.bdy.nc $EXP/.

# namelist_cfg
# nambdy: freeze the boundary conditions. Set to initial state
# ln_usr = true. User defined initial state and surface forcing. Here we use
# a contant density, and no wind.
# These are compiled into the executable. (In
#  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``).

# Submit job
cd $EXP
sbatch submit.slurm

## Check on queue
# squeue -u $USER
