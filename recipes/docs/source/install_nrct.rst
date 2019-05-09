Install NEMO Relocatable Configuration Tool (NRCT)
==================================================

**LIVLJOBS4**::

  ssh -Y livljobs4

  export CONFIG=SWPacific
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

Find java object by doing a ``which java`` and then following the trail
find  /usr/lib/jvm/jre-1.7.0-openjdk.x86_64/ -name libjvm.so -print
::

  cd $WORK/$USER
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS # Didn't need this on ARCHER...
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd $WORK/$USER/nrct/Python
  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $INPUTS

.. note : 6 Nov. Following git clone you might want to do:
    git fetch
    git checkout Generalise-tide-input
  to get the FES-tides enabled branch

  You have to manually set the TPXO or FES data source in Python/pynemo/tide/nemo_bdy_tide3.py

Just do it::

  pynemo -s namelist.bdy
---


**JASMIN**
*(8 May 2019)*
::


  livljobs6 ~ $ exec ssh-agent $SHELL
  livljobs6 ~ $ ssh-add ~/.ssh/id_rsa_jasmin
  Enter passphrase for /login/jelt/.ssh/id_rsa_jasmin:
  Identity added: /login/jelt/.ssh/id_rsa_jasmin (/login/jelt/.ssh/id_rsa_jasmin)
  ssh -A jelt@jasmin-login1.ceda.ac.uk
  ssh -A jelt@jasmin-sci1.ceda.ac.uk

  cd /home/users/$USER
  unset SSH_ASKPASS # Didn't need this on ARCHER...
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd


Now there are a couple of files that need updating because they are not in the
repo yet. From livljobs6 (having launched the ssh-agent)::

  livljobs6>
  scp /work/jdha/anaconda/envs/pynemo_jelt/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py jelt@jasmin-login1.ceda.ac.uk:nrct/Python/pynemo/.
  scp /work/jdha/anaconda/envs/pynemo_jelt/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/nemo_bdy_extr_tm3.py jelt@jasmin-login1.ceda.ac.uk:nrct/Python/pynemo/.

Try installing anaconda (into a gws as it is mega)::

  cd /gws/nopw/j04/campus/downloads/
  wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh
  bash Anaconda-2.3.0-Linux-x86_64.sh
  On option: Install to /gws/nopw/j04/campus/anaconda

  export PATH=/home/users/jelt/anaconda/bin:$PATH
  cd /home/users/$USER
  conda create --name nrct_env scipy=0.16.0 numpy
  source activate nrct_env
  conda install matplotlib
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius
  conda install libgfortran=1.0.0
  conda install basemap netcdf4

Check the library list with ``pip freeze``::

  basemap==1.0.7
  certifi==2016.2.28
  Cython==0.26
  jnius==1.1.dev0
  lxml==3.8.0
  matplotlib==1.4.3
  netCDF4==1.1.9
  numpy==1.9.2
  pyparsing==2.0.3
  python-dateutil==2.6.1
  pytz==2017.2
  requests==2.14.2
  scipy==0.16.0
  six==1.10.0
  thredds-crawler==1.0.0



  cd
  cd nrct/Python
  git checkout ORCA0083 # Get on the appropriate branch
  #export LD_LIBRARY=/opt/local/lib/gcc48/gcj-4.8.2-14/:$LD_LIBRARY_PATH # I had this but think it is wrong / did nothing
  python setup.py build
  export PYTHONPATH=$HOME/anaconda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix $HOME/anaconda/envs/nrct_env


That should complete the build. Now copy all the necessary input files to a target space
E.g. (with having setup sshfs)::

  livljobs6 INPUTS> scp namelist.bdy jelt@jasmin-login1.ceda.ac.uk:tmp/.

Then move to gws using compute server (jasmin-sci1.ceda.ac.uk)::

  mv ~/tmp/*  /gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/.


Execute pynemo (jasmin-sci1.ceda.ac.uk). Exmaple::

  cd /gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/
  export PATH=/home/users/jelt/anaconda/bin:$PATH
  source activate nrct_env # If required

  export PYTHONPATH=$HOME/anaconda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH

  module load java/1.8.0
  export LD_LIBRARY=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.201.b09-2.el6_10.x86_64/jre/lib/amd64/server/:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy # e.g. tested on thredds_namelist_jan14.bdy

Outputs:
  ---

In progress builds on other machines
====================================


  On **MacOSX**. *(26 Oct 2017)* (conda 4.3.30, python2.7) NB couldn't find ``libgfortran=1.0.0``. I used tcsh so you need to
  switch to bash to get ``source`` working.

  ::
    git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct
    conda create --name nrct_env scipy=0.16.0 numpy matplotlib=1.5.1 basemap netcdf4 libgfortran=3.0.0

    bash
    source $HOME/anaconda/envs/nrct_env/bin/activate $HOME/anaconda/envs/nrct_env
    conda install -c conda-forge seawater=3.3.4
    conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
    conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

  I am not convinced that the following worked properly; not sure that I ironed out
  all the problems with the java virtual machine (which is also trouble for ARCHER).
  Anyway, this is what I did::

    cd nrct/Python
    export LD_LIBRARY=/opt/local/lib/gcc48/gcj-4.8.2-14/:$LD_LIBRARY_PATH
    python setup.py build
    export PYTHONPATH=$HOME/anaconda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
    python setup.py install --prefix ~/anaconda/envs/nrct_env

    pynemo -s data/namelist.bdy

  It looked like it worked though I didn't try it with actualy bathymetry and coords files

  ---

  On **ARCHER**

  #From the archer login node::

    module load git/1.8.5.2
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
  ::

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



    ---
