=======================================
Setting up a SE Asia NEMO configuration
=======================================

URL:: http://nemo-reloc.readthedocs.io/en/latest/SEAsia.html

Issues that arose
=================

* Fix the sosie make.macro command

Follow PyNEMO recipe for Lighthouse Reef: ``http://pynemo.readthedocs.io/en/latest/examples.html#example-2-lighthouse-reef``
Follow PyNEMO recipe for LBay on ARCHER (not complete because PyNEMO java issue): ``http://nemo-reloc.readthedocs.io/en/latest/LBay.html``.
Follow PyNEMO recipe for LBay on ARCHER/Livljobs4: ``http://nemo-reloc.readthedocs.io/en/latest/LBay_livljobs4.html``.


Summary / Plan
==============

#. Update to ORCHESTRA code base

#. Extend domain in a new experiments: 1/12 degree




Recipe Notes
============

Define working directory, and other useful shortcuts::

  export CONFIG=ACCORD
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG
  export TDIR=$WDIR/dev_r6998_ORCHESTRA/NEMOGCM/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_SEAsia

#Load modules::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

Make the directories::

  mkdir $WDIR
  mkdir $INPUTS

Build NEMO ORCHESTRA branch @ r8395::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/NERC/dev_r6998_ORCHESTRA@8395

Use Dave's XIOS file (see ``%XIOS_HOME``)::

  cp /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.

Make a new config directory structure::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean

Copy Maria's MY_SRC::

  rsync -uartv /work/n01/n01/mane1/ORCHESTRA/NEMOGCM/CONFIG/TIDE/MY_SRC/ $CDIR/$CONFIG/MY_SRC

Copy Maria's cpp flags (without the lim flag)::

  vi $CDIR/$CONFIG/cpp_ACCORD.fcm

   bld::tool::fppkeys key_mpp_mpi          \
                      key_bdy              \
                      key_tide            \
                      key_zdftke           \
                      key_netcdf4          \
                      key_iomput           \
                      key_nosignedzero    \
                      key_trabbl        \
                      key_zdfddm        \
                      key_diaharm      \
                      key_xios2


Build opa::

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

**It compiles!!**

Make an EXP directory::

  mkdir $EXP


----


Build some NEMO tools (on ARCHER)
---------------------------------

For the old build ``dev_r4621_NOC4_BDY_VERT_INTERP/`` we applied patches.
For the new build ``dev_r6998_ORCHESTRA``, we do not. For some reason GRIDGEN doesn’t like INTEL::

  cd $TDIR
  ./maketools -n WEIGHTS -m XC_ARCHER_INTEL
  ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module swap PrgEnv-intel PrgEnv-cray
  module load cray-netcdf cray-hdf5
  ./maketools -n GRIDGEN -m XC_ARCHER

  module swap PrgEnv-cray PrgEnv-intel

----

*(3 Oct 2017)*


5. Generate mesh and mask files for open boundary conditions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Run the model to generate the mesh and mask files. Copy
 the input files from Liverpool to ARCHER::

  ssh livljobs4

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
  export START_FILES=$WDIR/START_FILES # generic stuff for making more stuff. Mostly code.
  export INPUTS=$WDIR/INPUTS         # config specific stuff that gets made and is for running NEMO

  # Copy into the $EXP directory
  for file in bathy_meter.nc coordinates.nc; do  scp $INPUTS/$file  jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/ACCORD/EXP_SEAsia/$file; done
  for file in runscript_archer namelist_cfg namelist_ref iodef.xml; do  scp $START_FILES/$file  jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/ACCORD/EXP_SEAsia/$file; done

---

Return to **ARCHER**. Make sure the executables are in the EXP dir.
Using XIOS from Dave::

  ln -s /work/n01/n01/munday/XIOS/bin/xios_server.exe $EXP/.
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa

Also link in the extra XML files::

  ln -s $EXP/../../SHARED/field_def.xml $EXP/.
  ln -s $EXP/../../SHARED/domain_def.xml $EXP/.



Edit the namelist files for this configuration::

  cd $EXP
  ncdump -h coordinates.nc
  x = 683 ;
  y = 553 ;

  vi namelist.cfg
  ...
  cn_exp      =   "SEAsia"  !  experience name
  ...
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     cp_cfg      =  "seasia"                !  name of the configuration
     jp_cfg      =     012               !  resolution of the configuration
     jpidta      =     683               !  1st lateral dimension ( >= jpi )
     jpjdta      =     553               !  2nd    "         "    ( >= jpj )
     jpkdta      =      51               !  number of levels      ( >= jpk )
     jpiglo      =     683               !  1st dimension of global domain --> i =jpidta
     jpjglo      =     553               !  2nd    -                  -    --> j  =jpjdta

**ACTION: There are further edits to be made for when the model is actually run**
**E.g. other filename instances of LBay**

**NOT SURE WHAT jp_cfg does. It might be passive? I changed it to 012, representing 1/12**


Edit the runscript to include modules and the Account name (n01-NOCL)::

  vi runscript_archer

  #!/bin/bash
  #PBS -N SEAsia
  #PBS -l select=5
  #PBS -l walltime=00:20:00
  #PBS -A n01-NOCL

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel
  ...

Submit::

  qsub -q short runscript_archer


**IT WORKED WHEN I USED YESTERDAY'S LBAY EXECUTABLE. HOWEVER I LOST IT WHEN**
**I RECOMPILED THE CODE TO FIND OUT WHAT WAS DIFFERENT...**

----

*(6 March 2017)*

If that works, we then need to rebuild the mesh and mask files in to single files for the next step::

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_mask 96
  #mv mesh_mask.nc $WDIR/INPUTS
  #rm mesh_* mask_* LBay_0000*
  #cd $INPUTS


THIS IS WHERE START WITH LIVLJOBS4 to create boundary files with PyNEMO


----

On **livljobs4** copy the mesh_mask file from **ARCHER**::

  ssh livljobs4

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
  export START_FILES=$WDIR/START_FILES # generic stuff for making more stuff. Mostly code.
  export INPUTS=$WDIR/INPUTS         # config specific stuff that gets made and is for running NEMO

  # Copy from the $EXP directory
  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/ACCORD/dev_r6998_ORCHESTRA/NEMOGCM/CONFIG/ACCORD/EXP_SEAsia/mesh_mask.nc $INPUTS/mesh_mask.nc


6. Generate boundary conditions with PyNEMO: Create netcdf abstraction wrapper
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this section there are two stages.
* generate a ncml file which describes the files needed to create boundary conditions
* generate a namelist.bdy file which controls the actual boundary condition generation.

For each parent data set a new pair of (``*.ncml``, ``namelist.bdy``) are needed.
Here I attempt to use parent data from:
* AMM60 local data (doesn't yet work because of the sigma levels)
* thredds server (as in the LH_REEF example)
* NNA local data (easiest ?)

First install PyNEMO if not already done so. Full description::

  ssh -Y livljobs4
  cd /work/$USER
  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
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
  export PYTHONPATH=/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $WDIR/INPUTS


I suggest managing the namelist.bdy file after the ``ncml`` file is generated.
A fresh ``ncml`` file can be generated automatically or an existing one can be
edited.


6a. Generate ncml files: thredds_inputs_src.ncml
++++++++++++++++++++++++++++++++++++++++++++++++

**Untested**
In the pynemo_ncml_generator if using the thredds server use:
Source directory: ``http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data``

*(16 March 2017)*
Created a thredds_inputs_src.ncml file to access ORCA12 data from the
thredds server. Note that the pynemo_ncml_generator populates this file with available
files according to the input regular expressions::

  cp $START_FILES/thredds_inputs_src.ncml $INPUTS/.
  cd $INPUTS
  vi thredds_inputs_src.ncml

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791206d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791201d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791126d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791121d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791116d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791111d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791106d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791206d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791201d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791126d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791121d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791116d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791111d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791106d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791206d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791201d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791126d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791121d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791116d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791111d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791106d05U.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791101d05U.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791206d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791201d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791126d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791121d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791116d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791111d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791106d05V.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791101d05V.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791206d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791201d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791126d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791121d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791116d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791111d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791106d05T.nc" />
          <ns0:netcdf location="http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/data/ORCA025-N206_19791101d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>



6b. Generate the namelist.bdy file for PyNEMO
+++++++++++++++++++++++++++++++++++++++++++++


Copy the PyNEMO template namelist.bdy from the START_FILES dir::

  cd $INPUTS
  cp $START_FILES/namelist.bdy $INPUTS/.

Edit namelist.bdy to for the configuration name and ``ncml`` file name. **Note
need the slash following OUTPUT**::

  vi namelist.bdy
  sn_src_dir = './thredds_inputs_src.ncml'       ! src_files/'
  sn_dst_dir = '/work/jelt/NEMO/SEAsia/INPUTS/'
  sn_fn      = 'SEAsia'                 ! prefix for output files
  ...
  cn_mask_file   = './mesh_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)

Now edit the pynemo namelist file. Add location of grid information. Lots of
the separate mesh and mask files are combined into the new mesh_mask.nc output.
 Note use ncml to convert to variables without *_0.

 Make sure the timestamps correspond to the input data. (Not sure this is
  important yet *3-Oct-17*)
Turn off as many things as possible to help it along.
Turned off ``ln_mask_file``. James said it was for outputting a new mask file
but it might have given me trouble.

NB I have a namelist.bdy file for each ncml configuration
* namelist.bdy_AMM60 (should use for LBay and Solent)
* namelist.bdy_thredds (Used here. Uses global 1/12 degree data)
* namelist.bdy_NNA (used for LBay)



7. Generate boundary conditions with PyNEMO: Run PyNEMO
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Using livljobs4

*(3 Oct 2017)*
::

  export CONFIG=SEAsia
  export WDIR=/work/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES

  cp $START_FILES/namelist.bdy $INPUTS/namelist.bdy
  cp $START_FILES/thredds_inputs_src.ncml $INPUTS/thredds_inputs_src.ncml
  cp $START_FILES/inputs_dst.ncml $INPUTS/inputs_dst.ncml

.. Need to grab some INPUT files from ARCHER because I am not building them with livljobs4::

  rsync -uartv jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/ $WDIR/INPUTS
  cd $WDIR/INPUTS

.. Make sure the NNA data is available::

  mkdir $WDIR/INPUTS/NNA
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mesh_hgr.nc $WDIR/INPUTS/NNA/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/mask.nc $WDIR/INPUTS/NNA/.
  for file in NNA_*200001*nc ; do scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/NNA/$file $WDIR/INPUTS/NNA/. ; done

.. Make sure the destination meshes exist. These were generated by running the new config for a timestep::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/mesh_zgr.nc $WDIR/INPUTS/.
  scp jelt@login.archer.ac.uk:/work/n01/n01/jdha/LBay/INPUTS/mesh_hgr.nc $WDIR/INPUTS/.



cd $INPUTS
Edit the namelist.bdy
Sort out the dates



Generate the boundary conditions again, with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -g -s namelist.bdy

Sort of worked. Problem with access to the ORCA12 data::

  INFO:pynemo.reader.ncml:(1, 1, 553, 683)
  INFO:pynemo.reader.ncml:(1, 1, 553, 683)
  INFO:pynemo.reader.ncml:(1, 1, 553, 683)
  INFO:pynemo.reader.ncml:(1, 1, 553, 683)
  INFO:pynemo.profile:4.68
  INFO:pynemo.profile:['wv', 'wu', 'wt', 'u', 't', 'v']
  INFO:pynemo.profile:wv  is (51, 1427) and a nan  max: [ 18.71657181  18.68452454]
  INFO:pynemo.profile:wu  is (51, 1426) and a nan  max: [ 18.79060745  18.71673203]
  INFO:pynemo.profile:wt  is (51, 12704) and a nan  max: [ 18.7412262   18.69191742]
  INFO:pynemo.profile:u  is (51, 1426) and a nan  max: [ 21.30413628  21.21548653]
  INFO:pynemo.profile:t  is (51, 12704) and a nan  max: [ 21.24488068  21.185709  ]
  INFO:pynemo.profile:v  is (51, 1427) and a nan  max: [ 21.21529484  21.17683792]
  INFO:pynemo.profile:horizontal grid info
  INFO:pynemo.profile:0.02
  INFO:pynemo.profile:read and assign netcdf data, lon lat
  INFO:pynemo.profile:0.02
  INFO:pynemo.profile:range is 12704
  INFO:pynemo.profile:range is 1426
  INFO:pynemo.profile:range is 1427
  INFO:pynemo.profile:(12704,)
  ./thredds_inputs_src.ncml
  Oct 03, 2017 4:31:36 PM org.apache.http.impl.client.DefaultRequestDirector tryConnect
  INFO: I/O exception (java.net.NoRouteToHostException) caught when connecting to the target host: No route to host (Host unreachable)
  Oct 03, 2017 4:31:36 PM org.apache.http.impl.client.DefaultRequestDirector tryConnect
  INFO: Retrying connect
  Traceback (most recent call last):
    File "/login/jelt/.conda/envs/nrct_env/bin/pynemo", line 11, in <module>
      load_entry_point('pynemo==0.2', 'console_scripts', 'pynemo')()
    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/pynemo_exe.py", line 44, in main
      profile.process_bdy(setup_file, mask_gui)
    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/profile.py", line 271, in process_bdy
      reader = factory.GetReader(settings['src_dir'],acc)
    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/reader/factory.py", line 27, in GetReader
      return NcMLReader(uri,t_adjust)
    File "/login/jelt/.conda/envs/nrct_env/lib/python2.7/site-packages/pynemo-0.2-py2.7.egg/pynemo/reader/ncml.py", line 56, in __init__
      self.dataset = NetcdfDataset.openFile(self.uri, None)
    File "jnius_export_class.pxi", line 893, in jnius.JavaMultipleMethod.__call__ (jnius/jnius.c:23945)
    File "jnius_export_class.pxi", line 624, in jnius.JavaMethod.__call__ (jnius/jnius.c:20680)
    File "jnius_export_class.pxi", line 790, in jnius.JavaMethod.call_staticmethod (jnius/jnius.c:22522)
    File "jnius_utils.pxi", line 65, in jnius.check_exception (jnius/jnius.c:3815)
  jnius.JavaException: JVM exception occurred: org.apache.http.conn.ConnectTimeoutException: Connect to esurgeod.noc.soton.ac.uk:8080 timed out
  Exception AttributeError: "'Reader' object has no attribute 'dataset'" in <bound method Reader.__del__ of <pynemo.reader.ncml.Reader object at 0x7fe5e80a0190>> ignored

  (nrct_env)livljobs4 INPUTS $








Old notes
---------



----

*(27 Sept 2017)*

Build the new SE Asia configuration at 1/12 degree, R12
=======================================================

Generate new coordinates file
=============================

Inspect TPXO harmonic amplitudes to find a good cut off location for boundaries:

cd /work/jelt/tpxo7.2
ferret
go  plot_SEAsia_harmonics.jnl

... note::
  ! plot_SEAsia_harmonics.jnl
  ! Plot tpxo harmonics for the SE Asia region.
  ! Want to build a NEMO config without significant amphidromes on the boundary

  use h_tpxo7.2.nc

  set win 1
  set viewport ul
  shade/k=1/j=300:700/i=250:500/levels=(0,1,0.1)/title="M2" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=2/j=300:700/i=250:500/levels=(0,1,0.1)/title="S2" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=3/j=300:700/i=250:500/levels=(0,1,0.1)/title="N2" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=4/j=300:700/i=250:500/levels=(0,1,0.1)/title="K2" HA, lon_z, lat_z; go fland

  set win 2
  set viewport ul
  shade/k=5/j=300:700/i=250:500/levels=(0,1,0.1)/title="K1" HA, lon_z, lat_z; go fland
  set viewport ur
  shade/k=6/j=300:700/i=250:500/levels=(0,1,0.1)/title="O1" HA, lon_z, lat_z; go fland
  set viewport ll
  shade/k=7/j=300:700/i=250:500/levels=(0,1,0.1)/title="P1" HA, lon_z, lat_z; go fland
  set viewport lr
  shade/k=8/j=300:700/i=250:500/levels=(0,1,0.1)/title="Q1" HA, lon_z, lat_z; go fland


Conclusion. Plot the proposed domain::

  $livljobs2$ scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/INPUTS/coordinates_ORCA_R12.nc ~/Desktop/.

  ferret
  use coordinates_ORCA_R12.nc
  set win 1; shade/X=50:730/Y=1250:1800 E2T, nav_lon, nav_lat ; go fland
  set win 2; set viewport upper; shade/i=50:730/j=1250:1800 NAV_LAT
  set win 2; set viewport lower; shade/i=50:730/j=1250:1800 NAV_LON




---

----
