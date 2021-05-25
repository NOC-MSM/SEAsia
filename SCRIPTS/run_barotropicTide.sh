
#!/bin/bash

#:'
#
#***********************
#run_EXP_barotropicTide.sh
#***********************
#'

# Run the experiment with contant T,S initial condition with tides-only.

#::

export CONFIG=NEMOconstTS
export EXP=$WDIR/RUN_DIRECTORIES/EXP_barotropicTide

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

# Link in tidal bondary forcing
#ln -s /work/n01/n01/annkat/SEAsia_HadGEM_R12/TIDES $EXP/.
ln -s $WDIR/INPUTS/TIDES $EXP/.

# Link in boundary files (just coordinates.bdy.nc)
ln -s $WDIR/INPUTS/OBC/coordinates.bdy.nc $EXP/.

# namelist_cfg
# nambdy: Except for tides, freeze the boundary conditions. Set to initial state
# ln_usr = true. User defined initial state and surface forcing. Here we use
# homogenous initial conditions and no met forcing.
# with the expression being compiled into the executable. (In
#  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``).

# Submit job
sbatch submit.slurm

## Check on queue
# squeue -u $USER
