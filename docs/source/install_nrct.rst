Install NEMO Relocatable Configuration Tool (NRCT)
==================================================

**LIVLJOBS4**::

  ssh -Y livljobs4

  export CONFIG=LBay
  export WORK=/work
  export WDIR=$WORK/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  #export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
  #export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
  #export EXP=$CDIR/$CONFIG/EXP00

  cd $WORK/$USER
  mkdir $WDIR
  module load anaconda/2.1.0  # Want python2
  conda create --name nrct_env scipy=0.16.0 numpy matplotlib=1.5.1 basemap netcdf4 libgfortran=1.0.0
  source activate nrct_env
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 # Note had to add https path
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

Find java object by doing a which java and then following the trail
find  /usr/lib/jvm/jre-1.7.0-openjdk.x86_64/ -name libjvm.so -print
::

  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS # Didn't need this on ARCHER...
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd nrct/Python
  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $INPUTS

---

On **ARCHER**::

#From the archer login node::

  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct

#Compute a compute node::

  ssh -Y espp1
  cd ~
  module load anaconda/4.3.1-python2
  conda create --name nrct_env scipy=0.16.0 numpy matplotlib=1.5.1 basemap netcdf4 libgfortran=1.0.0
  source activate nrct_env
  conda install -c conda-forge seawater=3.3.4
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

#This sort of helped::

  conda install libxcb

# Change the path where java is sought: $whereis java
#Find the shared object libjvm.so in
#/usr/lib64/jvm/jre-1.7.0-ibm/bin/j9vm
#and link to it::

  export LD_LIBRARY_PATH=/usr/lib64/jvm/jre-1.7.0-ibm/bin/j9vm:$LD_LIBRARY_PATH

  cd nrct/Python
  python setup.py build
  export PYTHONPATH=/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env


Try it out::

  pynemo -g -s /work/n01/n01/jelt/LBay/INPUTS/namelist.bdy

----

**ERROR**
.. warning:

  Didn't find a proxy environment variable
  Traceback (most recent call last):
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/bin/pynemo", line 11, in <module>
    load_entry_point('pynemo==0.2', 'console_scripts', 'pynemo')()
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 570, in load_entry_point
    return get_distribution(dist).load_entry_point(group, name)
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 2687, in load_entry_point
    return ep.load()
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 2341, in load
    return self.resolve()
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 2347, in resolve
    module = __import__(self.module_name, fromlist=['__name__'], level=0)
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_exe.py", line 8, in <module>
    import profile
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 44, in <module>
    from pynemo import pynemo_settings_editor
  File "/home/n01/n01/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_settings_editor.py", line 8, in <module>
    from PyQt4 import QtGui
  ImportError: /usr/lib64/libxcb-xlib.so.0: undefined symbol: _xcb_unlock_io
