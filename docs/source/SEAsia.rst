=======================================
Setting up a SE Asia NEMO configuration
=======================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html

Issues that arose
=================

* Helpful to add in extra module command for first step
* Fix the sosie make.macro command
* Queue accound is now n01-NOCL in runscript

Following the recipe
====================

Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Define working directory::

  export WDIR=/work/n01/n01/jelt/lighthousereef/

Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

**James: You should add the above 3 lines as an e.g. in the docs**. The *swap* one threw me for a while.

Follow recipe. Step 1::

  cd $WDIR
  mkdir INPUTS
  cd INPUTS
  wget ftp.nerc-liv.ac.uk:/pub/general/jdha/inputs.tar.gz
  tar xvfz inputs.tar.gz
  rm inputs.tar.gz
  cd ../
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP@5709
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@629
  cd xios-1.0
  cp $WDIR/INPUTS/arch-XC30_ARCHER.* ./arch

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par

The comile works. But there are lots of warnings of three types::

  WARNING: /fs4/n01/n01/jelt/lighthousereef/xios-1.0/bld.cfg: LINE 40:
           ${PWD}: cyclic dependency, variable not expanded, ${PWD}: cyclic dependency, variable not expanded, ${PWD}: cyclic dependency, variable not expanded

  WARNING: CFLAGS__date.flags: duplicated targets for building:

  /home/n01/n01/jelt/work/lighthousereef/xios-1.0/Makefile:1418: warning: ignoring old commands for target `FPPKEYS.flags'


Step 2, as far as the ``makenemo`` call::

  export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS
  cd $CDIR/../NEMO/OPA_SRC/SBC
  patch -b < $WDIR/INPUTS/fldread.patch
  cd ../DOM
  patch -b < $WDIR/INPUTS/dommsk.patch
  cd ../BDY
  patch -b < $WDIR/INPUTS/bdyini.patch
  cd $CDIR
  rm $CDIR/../NEMO/OPA_SRC/TRD/trdmod.F90
  cp $WDIR/INPUTS/arch-* ../ARCH
  ./makenemo -n LH_REEF -m XC_ARCHER_INTEL -j 10

Build fails (as described) so remove ``key_lim2`` (as described) and re-issue the make command::

  ./makenemo -n LH_REEF -m XC_ARCHER_INTEL -j 10

  cp $WDIR/INPUTS/cpp_LH_REEF.fcm ./LH_REEF
  cp $WDIR/INPUTS/dtatsd.F90 LH_REEF/MY_SRC/

Step 3 (seemed to work OK)::

  cd $WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/WEIGHTS/src
  patch -b < $WDIR/INPUTS/scripinterp_mod.patch
  patch -b < $WDIR/INPUTS/scripinterp.patch
  patch -b < $WDIR/INPUTS/scrip.patch
  patch -b < $WDIR/INPUTS/scripshape.patch
  patch -b < $WDIR/INPUTS/scripgrid.patch
  cd ../../
  ./maketools -n WEIGHTS -m XC_ARCHER_INTEL
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module swap PrgEnv-intel PrgEnv-cray
  module load cray-netcdf cray-hdf5
  ./maketools -n GRIDGEN -m XC_ARCHER
  module swap PrgEnv-cray PrgEnv-intel
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

Step 4::

  cd $TDIR/GRIDGEN
  cp $WDIR/INPUTS/namelist_R12 ./
  ln -s namelist_R12 namelist.input
  ./create_coordinates.exe
  cp 1_coordinates_ORCA_R12.nc $WDIR/INPUTS/coordinates.nc

Step 5::

  cd $WDIR/INPUTS
  module load nco/4.5.0
  ncap2 -s 'where(topo > 0) topo=0' gebco_1_cutdown.nc tmp.nc
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
  rm tmp.nc
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco


Step 6 didn't find the ``make.macro`` file. Copy it from ``jdha`` and edit the path::

  cp /home/n01/n01/jdha/sosie/make.macro /home/n01/n01/jelt/sosie/.

  vi /home/n01/n01/jelt/sosie/make.macro
  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/jelt/local

Proceed with Step 6::

  cd ~
  mkdir local
  svn co svn://svn.code.sf.net/p/sosie/code/trunk sosie
  cd sosie

  FIX (copied from jdha instead): cp $WDIR/INPUTS/make.macro ./

  make
  make install
  export PATH=~/local/bin:$PATH
  cd $WDIR/INPUTS
  sosie.x -f initcd_votemper.namelist
  sosie.x -f initcd_vosaline.namelist
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper
  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Step 7, weight files for atmospheric forcing::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos
  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos

Step 8 (mesh and mask for open boundary condition generation)::

  cd $CDIR
  cp $WDIR/INPUTS/cpp_LH_REEF.fcm LH_REEF/
  ln -s $WDIR/INPUTS/bathy_meter.nc $CDIR/LH_REEF/EXP00/bathy_meter.nc
  ln -s $WDIR/INPUTS/coordinates.nc $CDIR/LH_REEF/EXP00/coordinates.nc
  cp $WDIR/INPUTS/runscript $CDIR/LH_REEF/EXP00
  cp $WDIR/INPUTS/namelist_cfg $CDIR/LH_REEF/EXP00/namelist_cfg
  cp $WDIR/INPUTS/namelist_ref $CDIR/LH_REEF/EXP00/namelist_ref
  ./makenemo clean
  ./makenemo -n LH_REEF -m XC_ARCHER_INTEL -j 10
  cd LH_REEF/EXP00
  ln -s $WDIR/xios-1.0/bin/xios_server.exe xios_server.exe

*(16 Jan 2017)* Edit the runscript to include modules and the Account name (n01-NOCL)::

  vi runscript

  #!/bin/bash
  #PBS -N LH_REEF
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel
  ...

Submit::

  qsub -q short runscript


*(17 Jan 17)*

Moved module load to .bashrc::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel







Doesn't work. No output. I've also tried a fresh rebuild of everything::

  execve error: No such file or directory
  aprun: Apid 24880812: Commands are not supported in MPMD mode
  aprun: Apid 24880812: Exiting due to errors. Application aborted

It looks, to me, like the runscript is missing some flags for something similar.

Just incase tried reloading modules and resubmitting (to standard queue after 8pm)::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  qsub runscript
  4195460.sdb

Still no joy :-(

---



Tried using James' xios executable::

  cd ~/work/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/EXP00
  ln -s /work/n01/n01/jdha/ST/xios-1.0/bin/xios_server.exe xios_server.exe

*(16 Feb 2017)*::

  cd /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/EXP00
  ln -s /work/n01/n01/jdha/TEST2/xios-1.0/bin/xios_server.exe xios_server.exe
  ln -s /work/n01/n01/jelt/lighthousereef/INPUTS/bathy_meter.nc bathy_meter.nc
  ln -s /work/n01/n01/jelt/lighthousereef/INPUTS/coordinates.nc coordinates.nc

  ln -s /work/n01/n01/jdha/TEST2/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/BLD/bin/nemo.exe opa

Spotted symlink issue in WDIR definition in ARCH file. Fix::

  cd /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG/LH_REEF/WORK>
  vi ../../../ARCH/arch-XC_ARCHER_INTEL.fcm
  ...
  %XIOS_HOME           /work/n01/n01/jdha/TEST2/xios-1.0

Recomile::

  cd /work/n01/n01/jelt/lighthousereef/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

  ./makenemo clean
  ./makenemo -n LH_REEF -m XC_ARCHER_INTEL -j 10

  cd LH_REEF/EXP00
  qsub -q short runscript



---

Get the BDY stuff together::

  cd LH_REEF/EXP00

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_zgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_hgr 96
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mask 96
  mv mesh_zgr.nc mesh_hgr.nc mask.nc $WDIR/INPUTS
  rm mesh_* mask_* LH_REEF_0000*
  cd $WDIR/INPUTS

install PyNEMO (**Note need to use https://ccpforge.cse.rl.ac.uk**)::

  cd ~
  module load anaconda
  conda create --name pynemo_env python scipy numpy matplotlib basemap netcdf4
  source activate pynemo_env
  conda install -c https://conda.anaconda.org/srikanthnagella seawater
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

New *(18 Feb 2017)*::
  conda install mock nose
::

  export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  svn checkout https://ccpforge.cse.rl.ac.uk/svn/pynemo
  cd pynemo/trunk/Python
  python setup.py build

Insert change to PYTHONPATH::

  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/pynemo/lib/python2.7/site-packages/:$PYTHONPATH

Proceed::

  python setup.py install --prefix ~/.conda/envs/pynemo
  cd $WDIR/INPUTS

Startup the PyNEMO and generate boundary conditions::

  ssh -Y espp1
  module load anaconda
  source activate pynemo_env
  export WDIR=/work/n01/n01/jelt/lighthousereef/
  cd $WDIR/INPUTS

Call to `pynemo_ncml_generator` doesn't work::

  (pynemo_env) jelt@esPP001:/work/n01/n01/jelt/lighthousereef/INPUTS> pynemo_ncml_generator

Error::

 Traceback (most recent call last):
   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/bin/pynemo_ncml_generator", line 11, in <module>
     load_entry_point('pynemo==0.2', 'console_scripts', 'pynemo_ncml_generator')()
   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 565, in load_entry_point

   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2598, in load_entry_point

   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2258, in load

   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2264, in resolve

   File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_ncml_generator.py", line 10, in <module>
     from PyQt4.QtGui import *
 ImportError: /usr/lib64/libxcb-xlib.so.0: undefined symbol: _xcb_unlock_io


Maybe I'm not supposed to be able to execute `pynemo_ncml_generator`; the comments suggest it is not necessary

Pressing on::

  (pynemo_env) jelt@esPP001:/work/n01/n01/jelt/lighthousereef/INPUTS> export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
  (pynemo_env) jelt@esPP001:/work/n01/n01/jelt/lighthousereef/INPUTS> export PYTHONPATH=~/.conda/envs/pynemo_env/lib/python2.7/site-packages:$PYTHONPATH
  (pynemo_env) jelt@esPP001:/work/n01/n01/jelt/lighthousereef/INPUTS> pynemo -g -s namelist.bdy

  Traceback (most recent call last):
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/bin/pynemo", line 11, in <module>
      load_entry_point('pynemo==0.2', 'console_scripts', 'pynemo')()
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 565, in load_entry_point

    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2598, in load_entry_point

    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2258, in load

    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/setuptools-27.2.0-py2.7.egg/pkg_resources/__init__.py", line 2264, in resolve

    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_exe.py", line 8, in <module>
      import profile
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 22, in <module>
      import numpy as np
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/__init__.py", line 170, in <module>
      from . import add_newdocs
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/add_newdocs.py", line 13, in <module>
      from numpy.lib import add_newdoc
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/lib/__init__.py", line 18, in <module>
      from .polynomial import *
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/lib/polynomial.py", line 19, in <module>
      from numpy.linalg import eigvals, lstsq, inv
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/linalg/__init__.py", line 51, in <module>
      from .linalg import *
    File "/home/n01/n01/jelt/.conda/envs/pynemo_env/lib/python2.7/site-packages/numpy/linalg/linalg.py", line 29, in <module>
      from numpy.linalg import lapack_lite, _umath_linalg
  ImportError: libgfortran.so.1: cannot open shared object file: No such file or directory

Hmm things are not working as I'd hoped

---

`Next steps: install pyNEMO <install_pynemo.html>`_
