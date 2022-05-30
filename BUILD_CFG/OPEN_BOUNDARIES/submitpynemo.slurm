#!/bin/bash
#SBATCH --job-name=pynemo
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --account=n01-ACCORD
#SBATCH --partition=serial
#SBATCH --qos=serial

source ~/.bashrc
conda activate pynemo
cd /work/n01/n01/$USER/SEVERN-SWOT/BUILD_CFG/OPEN_BOUNDARIES/
pynemo -s namelist_AMM15Severn.bdy
