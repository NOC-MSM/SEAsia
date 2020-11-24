#!/bin/bash

:'

*******************************
prep_ERA5_forcing.sh
*******************************

a. Extract ERA5 data
^^^^^^^^^^^^^^^^^^^^

The ERA5 data must first be downloaded from the BADC repository.
They can then be manipulated for use in forcing a regional configuration.

On a machine with adequate space, in ``SCRIPTS/OFFICIAL_GENERATE_NEMO_Forcing_NEWERA.py``,
 change the date range of the data and also the lat and lon coordinates of the
 region you want to extract. *IMPORTANT:* Make sure the extracted region is
 larger than your model domain size.

Also need to update ``path_EXTRACT`` and ``path_FORCING``, (but ``path_ERA5``
need not change if running on a NOC livljobs server).

Load the appropriate modules (I.e. nco, anaconda).
Create a python environment and run the script:

'

#::

    module load nco/gcc/4.4.2 # At NOC this is an updated module and is subsequently switched back to nco/gcc/4.4.2.ncwa in the main code_

    module load anaconda/3-5.1.0
    conda create -n ERA5_env numpy
    conda activate ERA5_env
    python OFFICIAL_Generate_NEMO_Forcing_NEWERA.py



:'


Depending on the size of the region and the data range of the data, this might
take a while.


b. Extract the Land-Sea Mask
^^^^^^^^^^^^^^^^^^^^

For better representation of the land/sea boundaries, extract a land-sea mask
from the ERA5 data. This can be done using one of the available ERA5 files.

First, use the ncks command to cut out the region you want (should be same as
that stated in ``OFFICIAL_Generate_NEMO_Forcing_NEWERA.py``):
'

#::

    export ERA5_dir = /projectsa/NEMO/Forcing/ERA5 # Or whatever you local path is
    cd FORCING
    ncks -d latitude,38.,68. -d longitude,332.,19. $ERA5_dir/era5_atmos_landseamask.nc ./my_era5_LSM.nc

:'

Next, use the ``create_LSM.py`` script to create a draft mask file (``ERA5_LSM_draft.nc``):
'

#::

    python create_LSM.py

:'
Finally, run the Matlab script ``draft_LSM_nc.m`` to create the
``ERA5_LSM.nc`` file.

.. note:: The mask for SEAsia had an anomalous point in the middle of the land
mass. This is removed in the matlab script. This part of the script will need
to be commented out, or edited, for different domains

'
