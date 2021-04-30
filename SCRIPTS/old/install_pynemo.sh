#!/bin/bash

:'

***********************
install_pynemo.sh
***********************

The open boundary conditions are generated using
`PyNEMO <https://github.com/NOC-MSM/PyNEMO>`__. This can be done in
**Jasmin** or in **livljobs** (NOC Liverpool server). To my knowledge
(katavouta,annkat) PyNEMO has not been used successfully on archer (when
I tried the problem was something with the libraries/java path I think).

Here the following instructions are for generating open boundary
conditions in **livljobs7**. ## **1. Install PyNEMO** (using
`Install\_NRCT.sh <https://github.com/NOC-MSM/SEAsia_ERSEM_R12/blob/master/FILES_START/OBC-TIDES/Install_NRCT.sh>`__)

define some paths:

.. code:: bash

    export CONFIG=SEAsia_R12
    export WORK=/work/$USER
    export PYTHON=$WORK/nrct/Python
    cd $WORK

load anaconda (you want python 2 in this case):

.. code:: bash

    module load anaconda/2.1.0

create your virtual environment that includes the packages that you want
(you can name it as you want but here I name it nrct-env):

.. code:: bash

    conda create --name nrct_env scipy=0.16.0 numpy=1.9.2 matplotlib=1.4.3 basemap=1.0.7 netcdf4=1.1.9 libgfortran=1.0.0

activate your virtual environment:

.. code:: bash

    source activate nrct_env

install the python packages that you need:

.. code:: bash

    conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 # Note had to add https path
    conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
    conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

find and set the java path; you may have to search for this path but in
my case (and if you are using livljobs7 probably that is your case too)
it was:

.. code:: bash

     export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk-1.7.0.241-2.6.20.0.el7_7.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

remove the need for asking password (I am not sure we need that in
livljobs but include it just in case):

.. code:: bash

    unset SSH_ASKPASS

Obtain PyNEMO from the git repository (in this case I use a previous
version of pyNEMO, in the future we will probably want to migrate in the
current updated version in
`PyNEMO <https://github.com/NOC-MSM/PyNEMO>`__):

.. code:: bash

    git clone https://bitbucket.org/jdha/nrct.git

PyNEMO has different branches, first let us set he one we want for
T,S,U,V,SSH boundaries:

.. code:: bash

    cd $WORK/nrct
    git checkout ORCA0083

build my PyNEMO and ready to use:

.. code:: bash

    cd $PYTHON
    python setup.py build
    export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
    python setup.py install --prefix ~/.conda/envs/nrct_env

'
