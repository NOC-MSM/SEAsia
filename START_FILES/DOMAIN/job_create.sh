#!/bin/bash
#PBS -N domain_cfg
#PBS -l walltime=00:20:00
#PBS -l select=1
#PBS -j oe
#PBS -A n01-ACCORD
# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m bea
#PBS -M annkat@noc.ac.uk
#! -----------------------------------------------------------------------------

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Set the number of threads to 1
#   This prevents any system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
# Change to the directory that the job was submitted from
ulimit -s unlimited

#===============================================================
# LAUNCH JOB
#===============================================================
echo `date` : Launch Job
aprun -n 1 -N 1 ./make_domain_cfg.exe

exit

