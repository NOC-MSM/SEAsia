# Running SWPacifc model on ARCHER

Copy Input and Start files from Jelt workspace on ARCHER to my own to allow me to build my NEMO config. Check SWPacific_archer_livljob4 recipe to find out what he did to make these.

    rsync --progress -r /work/n01/n01/jelt/SWPacific/INPUTS /work/n01/n01/thopri/SWPacific
    rsync --progress -r /work/n01/n01/jelt/SWPacific/START_FILES /work/n01/n01/thopir/SWPacific

Set environmental varaibles::

    export CONFIG=SWPacific
    export WORK=/work/n01/n01/thopri
    export WDIR=$WORK/$CONFIG
    export INPUTS=$WDIR/INPUTS
    export START_FILES=$WDIR/START_FILES
    export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
    export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
    export EXP=$CDIR/$CONFIG/EXP00
    
Load Modules::
    
    module swap PrgEnv-cray PrgEnv-intel
    module load cray-netcdf-hdf5parallel cray-hdf5-parallel

Follow build XIOS_2.rst and build_opa_orchestra.rst recipe on gitlab server to build XIOS2.0 and NEMO ver4 make sure to change $USER to jelt when copying ARCH files and fixes. Once built the relevent start files that Jeff has created can be copied over and a runscript created to run on ARCHER.

Copy over input and start files into NEMO config::

    cd $EXP
    rsync -tuv $INPUTS/bathy_meter.nc $EXP/.
    rsync -tuv $INPUTS/coordinates.nc $EXP/.
    rsync -tuv $INPUTS/coordinates.bdy.nc $EXP/.
    rsync -tuv $START_FILES/namelist_cfg $EXP/.
    rsync -tuv $INPUTS/domain_cfg.nc $EXP/.

Create symbolic links from EXP directory::

    ln -s $INPUTS $EXP/bdydta

Edit the output to have 1hrly SSH::

    vi file_def_nemo.xml
    ...
    <!-- 1h files -->
---

Create a short queue runscript. (Note: PBS -N jobname, PBS -m email)::
Note this is an updated script from the live notes at the bottom of Jeff's recipe.

    vi runscript

    #!/bin/bash
    #PBS -N SWPacific
    #PBS -l select=5
    #PBS -l walltime=00:20:00
    #PBS -A n01-NOCL
    # mail alert at (b)eginning, (e)nd and (a)bortion of execution
    #PBS -m bea
    #PBS -M thopri@noc.ac.uk

    module swap PrgEnv-cray PrgEnv-intel
    module load cray-netcdf-hdf5parallel
    module load cray-hdf5-parallel

    export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
    #  echo $(readlink -f $PBS_O_WORKDIR)
    # export OMP_NUM_THREADS=1

    cd $PBS_O_WORKDIR
    #
    echo " ";
    OCEANCORES=93
    XIOCORES=1
    ulimit -c unlimited
    ulimit -s unlimited

    rm -f core

    #aprun -n $OCEANCORES -N 24 ./opa
    aprun -b -n 20 -N 20 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa
    #aprun -b -n $XIOCORES -N 1 ./xios_server.exe : -n $OCEANCORES -N 24 ./opa

    exit

Edit ``namelist_cfg`` so that a fresh start is made::

    &namrun        !   parameters of the run
    !-----------------------------------------------------------------------
    cn_exp      =    "SWPacific"  !  experience name
    nn_it000    =  1   !  first time step
    nn_itend    =  4800 ! 10day=14400   !  last  time step (std 5475)
    nn_date0    =  20000101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
    nn_time0    =       0   !  initial time of day in hhmm
    nn_leapy    =       1   !  Leap year calendar (1) or not (0)
    ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)
    nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
    nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
    !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist

i.e. ln_rstart = .false. rather than .true. As Jeff has set it to restart. After first run I will also have to do this.

The tide harmonic calculation also needs to be changed::

    !-----------------------------------------------------------------------
    &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
    !-----------------------------------------------------------------------
    nit000_han = 1         ! First time step used for harmonic analysis
    nitend_han = 4800 ! 7200      ! Last time step used for harmonic analysis
    nstep_han  = 1        ! Time step frequency for harmonic analysis
    tname(1)   = 'M2'      ! Name of tidal constituents
    tname(2)   = 'S2'      ! Name of tidal constituents
    tname(3)   = 'K1'      ! Name of tidal constituents
    tname(4)   = 'O1'      ! Name of tidal constituents

Submit::

    cd $EXP
    qsub -q short runscript

Currently getting error in running on ARCHER::

    Error [CField* CField::getDirectFieldReference(void)] : In file '/work/n01/n01/thopri/xios-2.0_r1080/src/node/field.cpp', line 1452 -> field_ref="S1x" refers to an unknown field id.
    
Will investigate tomorrow. **TOMORROW** Have updated timesteps above to see if that fixes the issue. Sadly not :(

Yay! Jeff has helped so much today. The error above is due to the S1 harmonic being present in the file_def_nemo.xml file (I.e. its expecting to output the harmonic) however it was not listed in the namelist_cfg file so wasn't being computed.

I also had to reduce the time steps from 19200 to 4800 which at 360 seconds per time step equals 20 days. I think the 19200 went on too long for the short queue? Possibly related to a smaller time step earlier iteration of the model. (Not 100% sure as will need to increase it eventually). I also changed the start start to 01/01/2000 rather than the 02/03/2000. So far model run was a success. (Not checked output yet)

Next steps 16/11/2017
--------------------------

1. Add more harmonics to both file_def_nemo.xml and namelist_cfg (in terms of output just m2 and s2 for now until accuracy improves.) **DONE** 11 HC's added all the ones common to both NEMO and TOPEX
2. Check model output against M2 and S2 from FES and reality to see how model is doing. **DONE** its awful! on plus side additonal HC code in the script is working. WIll check FES output with truth next week.
3. Find out list of HC that are able to be computed within NEMO, compare with FES list **Written**
4. Increase length of simulation (ultimatly up to one year) **Need to get simulation running right**
5. Experiment with shorter time steps

Added more harmonics to namelist_cfg file::

    -----------------------------------------------------------------------
    &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
    !-----------------------------------------------------------------------
    nit000_han = 1         ! First time step used for harmonic analysis
    nitend_han = 4800 ! 7200      ! Last time step used for harmonic analysis
    nstep_han  = 1        ! Time step frequency for harmonic analysis
    tname(1)   = 'M2'      ! Name of tidal constituents
    tname(2)   = 'S2'      ! Name of tidal constituents
    tname(3)   = 'K1'      ! Name of tidal constituents
    tname(4)   = 'O1'      ! Name of tidal constituents
    tname(5)   = 'N2'      ! Name of tidal constituents
    tname(6)   = 'K2'      ! Name of tidal constituents
    tname(7)   = 'P1'     ! Name of tidal constituents
    tname(8)   = 'Q1'      ! Name of tidal constituents
    tname(9)   = 'M4'      ! Name of tidal constituents
    tname(10)  = 'Mm'      ! Name of tidal constituents
    tname(11)  = 'Mf'      ! Name of tidal constituents
    
These harmonics are the ones that are common to both NEMO and TXPO boundary forcing. Ideally will add some more particulary shallow water ones, Ms4 etc.

To output more harmonics the file_def_nemo file also needs to ammended currently only outputs M2 harmonics but will add S2. (More will be added later)::

    <file id="file20" name_suffix="_Tides" description="tidal harmonics" >
    <field field_ref="M2x"          name="M2x"      long_name="M2 Elevation harmonic real part"                       />
    <field field_ref="M2y"          name="M2y"      long_name="M2 Elevation harmonic imaginary part"                  />
    <field field_ref="M2x_u"        name="M2x_u"    long_name="M2 current barotrope along i-axis harmonic real part "       />
    <field field_ref="M2y_u"        name="M2y_u"    long_name="M2 current barotrope along i-axis harmonic imaginary part "  />
    <field field_ref="M2x_v"        name="M2x_v"    long_name="M2 current barotrope along j-axis harmonic real part "       />
    <field field_ref="M2y_v"        name="M2y_v"    long_name="M2 current barotrope along j-axis harmonic imaginary part "  />
    <field field_ref="S2x"          name="S2x"      long_name="S2 Elevation harmonic real part"                       />
    <field field_ref="S2y"          name="S2y"      long_name="S2 Elevation harmonic imaginary part"                  />
    <field field_ref="S2x_u"        name="S2x_u"    long_name="S2 current barotrope along i-axis harmonic real part "       />
    <field field_ref="S2y_u"        name="S2y_u"    long_name="S2 current barotrope along i-axis harmonic imaginary part "  />
    <field field_ref="S2x_v"        name="S2x_v"    long_name="S2 current barotrope along j-axis harmonic real part "       />
    <field field_ref="S2y_v"        name="S2y_v"    long_name="S2 currnet barotrope along j-axis haronic imaginary part "   />
    </file>

The additional lines with S2 in add the S2 harmonic real and imaginary parts for both elevation and u and v components of the current. This will need to be repeated for all harmonics. I also had to extend the wall time of the runscript to 25 mins as 20 mins was just sliightly too short. This means I have been running on the main queue of ARCHER. Currently takes just under 23 mins. Change line in runscript file from::

    #PBS -l walltime=00:20:00

too::

    #PBS -l walltime=00:25:00
    
Next steps are too lengthen the runtime to a couple of months rather than 20 days. Currently takes 23 mins for 20 days so an hour should give approx 60 days. Will check with Jeff appropiate run times.

Ok after chatting to Jeff its best to stay on short queue until model is ready for full run, i.e. all harmonics are being run and outputted. So I have halved the time step number to 2400 so the model runs for 10 days and just over 10 mins. Will add extra harmonics to output files and try restart model run next week. Will increase to 14 days so that full neap spring cycle is simulated. This equates to 3360 time steps. Also set wall time to 20 mins so model will run in short queue.

For Next Week: 20/11/17
-----------------------------

1. Output all harmonics that are being computed
2. Get model to use restart files rather than start from rest. **Not required for tides as restart doesn't help with harmonics**
3. Add extra HC's into NEMO, particulary shallow water ones Ms4 etc
4. QC model output check QA script throughly!
5. Longer run of at least 60 days (maybe 90?)

Monday:

Have added all the harmonics I can to the model, there are 19 so far and 11 are currently being outputted. Will update this to the full 17 this afternoon. Currently Jeff implies that a restart is not needed for tides as harmonics do not benefit from a restart. However the output from the model is not good so I am investigating to try and resolve as we need to match or ideally exceed FES.

Error in running the model, didn't realise I also needed to also edit the field_def_nemo-opa xml file as awell as the file_def_nemo xml file. Some additional harmonics were already in here (initial 11) which had confused me as the initial number of extra harmonics worked. Will add the fields I need to this file, the x and y fields for tide height and currents. Example of the field def xml file below::

    <field_group id="Tides_T" grid_ref="grid_T_2D" operation="once" >
    <!-- tidal composante -->
    <field id="M2x"          long_name="M2 Elevation harmonic real part "                             unit="m"        />
    <field id="M2y"          long_name="M2 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="S2x"          long_name="S2 Elevation harmonic real part "                             unit="m"        />
    <field id="S2y"          long_name="S2 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="N2x"          long_name="N2 Elevation harmonic real part "                             unit="m"        />
    <field id="N2y"          long_name="N2 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="K1x"          long_name="K1 Elevation harmonic real part "                             unit="m"        />
    <field id="K1y"          long_name="K1 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="O1x"          long_name="O1 Elevation harmonic real part "                             unit="m"        />
    <field id="O1y"          long_name="O1 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="Q1x"          long_name="Q1 Elevation harmonic real part "                             unit="m"        />
    <field id="Q1y"          long_name="Q1 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="M4x"          long_name="M4 Elevation harmonic real part "                             unit="m"        />
    <field id="M4y"          long_name="M4 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="K2x"          long_name="K2 Elevation harmonic real part "                             unit="m"        />
    <field id="K2y"          long_name="K2 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="P1x"          long_name="P1 Elevation harmonic real part "                             unit="m"        />
    <field id="P1y"          long_name="P1 Elevation harmonic imaginary part"                         unit="m"        />
    <field id="Mfx"          long_name="Mf Elevation harmonic real part "                             unit="m"        />
    <field id="Mfy"          long_name="Mf Elevation harmonic imaginary part"                         unit="m"        />
    <field id="Mmx"          long_name="Mm Elevation harmonic real part "                             unit="m"        />
    <field id="Mmy"          long_name="Mm Elevation harmonic imaginary part"                         unit="m"        />
    </field_group>

This section is also repeated for both U and V componets of the tidal currents. Once added for the extra HC's the model should run with the full 19 harmonics that are availble in NEMO.

Tomorrow I plan to investigate adding more harmonics to the NEMO code itself. In the config/WORK/tide.h90 file extra harmonics can be added using a alternative to doodson numbers. Nico has written about this so will read up. Would be useful to add the two shallow water harmonics that are in TPIO but not in NEMO. (Ms4 etc).

Must also make sure to start the harmonic analysis a few weeks into the model (one month) so that the model has settled down before analysing.

21/11/17
----------

Today I added the extra fields into the field_def_nemo-opa xml file so that the model now runs with all 19 harmonics for 10 days. The upshot of this is that the M2 and S2 amplitudes and phases are extreme (25,000 for S2 amp) suggesting that the model is either running wrong or that the harmonic analysis is struggling with having ten days with harmonics of 2 weeks or over periods. To check the model was working ok reran the model and changed the file_def_nemo xml file to output hourly sea surface heights::

    <file_group id="1h" output_freq="1h"  output_level="10" enabled=".TRUE."> <!-- 1h files -->
     <file id="file19" name_suffix="_SSH" description="ocean T grid variables" >
     <field field_ref="ssh"          name="zos"   />
     </file>
    </file_group>
    
This will output a netcdf of sea surface height at hourly intervals and can simply be commeted out by using ! in front of the three lines between the file group tags.

Upon checking this netcdf file, the SSH was found to have a ditribution and values that would be expected. (You can see the tide moving around the domain when advancing through the arrays in the nc file).

As the model seems to be working a longer simulation was set up to last for 5 weeks. One week spin up with 4 weeks dedicated to the harmonic analysis. This results in 8400 time steps of 360 seconds. It was estimated that as 10 days takes 15 mins, 35 days would take around 52 mins. Therefore the wall time in the runscript was changed to an hour (01:00:00) and the namelist_cfg file was edited with the new number of time steps::

    cn_exp      =    "SWPacific"  !  experience name
    nn_it000    =  1   !  first time step
    nn_itend    =  8400 ! 10day=14400   !  last  time step (std 5475)
    
    nn_stock    =   8400 ! 14400   !  frequency of creation of a restart file (modulo referenced to 1)
    nn_stocklist = 0,0,0,0,0,0,0,0,0,0 ! List of timesteps when a restart file is to be written
    nn_write    =  240 ! 14400   !  frequency of write in the output file   (modulo referenced to nn_it000)

nn_itend defines the length of the simulation so was set to 8400 time steps (8400 x 360 s = 35 days) The restart file was set to be generated at the end of the simulation, (nn_stock) and the output into the ocean.output file was set to daily updates (nn_write).

The harmonic analysis was also set within namelist_cfg::

    !-----------------------------------------------------------------------
    &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
    !-----------------------------------------------------------------------
    nit000_han = 2400         ! First time step used for harmonic analysis
    nitend_han = 8400 ! 7200      ! Last time step used for harmonic analysis
    nstep_han  = 1        ! Time step frequency for harmonic analysis
    tname(1)   = 'M2'      ! Name of tidal constituents
    tname(2)   = 'S2'      ! Name of tidal constituents
    tname(3)   = 'K1'      ! Name of tidal constituents
    tname(4)   = 'O1'      ! Name of tidal constituents
    tname(5)   = 'N2'      ! Name of tidal constituents
    tname(6)   = 'K2'      ! Name of tidal constituents
    tname(7)   = 'P1'      ! Name of tidal constituents
    tname(8)   = 'Q1'      ! Name of tidal constituents
    tname(9)   = 'M4'      ! Name of tidal constituents
    tname(10)  = 'Mm'      ! Name of tidal constituents
    tname(11)  = 'Mf'      ! Name of tidal constituents
    tname(12)  = '2N2'     ! Name of tidal constituents
    tname(13)  = 'MU2'     ! Name of tidal constituents
    tname(14)  = 'NU2'     ! Name of tidal constituents
    tname(15)  = 'L2'      ! Name of tidal constituents
    tname(16)  = 'T2'      ! Name of tidal constituents
    tname(17)  = 'S1'      ! Name of tidal constituents
    tname(18)  = 'Msqm'    ! Name of tidal constituents
    tname(19)  = 'Mtm'     ! Name of tidal consituents

The nitend_han is set to the last time step used for analysis (35 days) and nit000_han is set to define when to start harmonic analysis (7 days in). The model was then run and found to crash as time step 4801 where the zonal velocity exceeded 20 m/s :(  The case of this will be investigated tomorrow. I think maybe reducing the time step will help. Useful doc here::

    http://www.hector.ac.uk/cse/distributedcse/reports/nemo/nemo_notes/node46.html

=====================
Some Useful Commands
=====================

**GREP::**

    grep "what your looking for" filename
    e.g. grep "ssh max" ocean.output
    
This command finds all lines with the string and prints them to the terminal. In the example above the max ssh for each written output period is printed. (In this case every 240 time steps or daily.)::

    ==>> time-step=            1  ssh max:  0.477583651140116
    ==>> time-step=          241  ssh max:   1.93981543025109
    ==>> time-step=          481  ssh max:   2.97875556152862
    ==>> time-step=          721  ssh max:   3.41631958976275
    ==>> time-step=          961  ssh max:   3.63379636023758
    ==>> time-step=         1201  ssh max:   3.65780427222438
    ==>> time-step=         1441  ssh max:   3.51990549861162
    ==>> time-step=         1681  ssh max:   3.19143228284551
    ==>> time-step=         1921  ssh max:   2.60460558366317
    ==>> time-step=         2161  ssh max:   1.93733415603733
    ==>> time-step=         2401  ssh max:   1.50624956040312
    ==>> time-step=         2641  ssh max:   1.38879686447409
    ==>> time-step=         2881  ssh max:   2.29851817169739
    ==>> time-step=         3121  ssh max:   2.30086809134897
    ==>> time-step=         3361  ssh max:   1.74965464073241
    ==>> time-step=         3601  ssh max:   1.20885310471017
    ==>> time-step=         3841  ssh max:   1.20840138522246
    ==>> time-step=         4081  ssh max:   2.98698803564672
    ==>> time-step=         4321  ssh max:   4.15497319068015
    ==>> time-step=         4561  ssh max:   4.56194437857435
    ==>> time-step=         4801  ssh max:   4.42258997404563
    
For this monthly simulation the SSH does not seem to be excessive.

To view nc files on ARCHER::

    module load ncview
    ncview filename
    
**Note** Need to make sure the -Y flag is used when sshing in. (Currently this doesn't work for me so I scp and use panoply)

22/11/17
=======

ARCHER down for maintenance today so no NEMO modelling today!

**Some thoughts**

-- Need to reduce time step of model (to 180 seconds?) to see if that makes it more stable.
-- Investigate adding extra harmonics particularly shallow water ones.

23/11/17
=======

Running model for longer results in a blowup where the zonal velocity exceeds a threshold on 20 m/s. This occurs at around 3.8 weeks in. This is same issues found on the 21/11/17. Some research suggested that reducing the time step, currently 360 seconds could help. Both 240 and 180 seconds per time step have been used with no change in the outcome.

Currently trying to setup a restart so that model can be started just before blowup allowing output to be inspected to find issue. To start the following lines within namelist_cfg need to be changed::

    !-----------------------------------------------------------------------
    &namrun        !   parameters of the run
    !-----------------------------------------------------------------------
    cn_exp      =    "SWPacific"  !  experience name
    nn_it000    =  1   !  first time step
    nn_itend    =  16800 ! 10day=14400   !  last  time step (std 5475)
    nn_date0    =  20000101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
    nn_time0    =       0   !  initial time of day in hhmm
    nn_leapy    =       1   !  Leap year calendar (1) or not (0)
    ln_rstart   = .false.    !  start from rest (F) or from a restart file (T)
    nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
    nn_rstctl   =    2            !  restart control ==> activated only if ln_rstart=T
    !                             !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
    !                             !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
    !                             !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
    cn_ocerst_in    = "SWPacific_00012677_restart_tide"   !  suffix of ocean restart name (input)
    cn_ocerst_indir = "."         !  directory from which to read input ocean restarts
    cn_ocerst_out   = "restart_tide"   !  suffix of ocean restart name (output)
    cn_ocerst_outdir= "."         !  directory in which to write output ocean restarts
    ln_iscpl    = .true.   !  cavity evolution forcing or coupling to ice sheet model
    nn_istate   =       0   !  output the initial state (1) or not (0)
    ln_rst_list = .false.   !  output restarts at list of times using nn_stocklist (T) or at set frequency with nn_stock (F)
    nn_stock    =   12677 ! 14400   !  frequency of creation of a restart file (modulo referenced to 1)
    nn_stocklist = 0,0,0,0,0,0,0,0,0,0 ! List of timesteps when a restart file is to be written
    nn_write    =  480  ! 14400   !  frequency of write in the output file   (modulo referenced to nn_it000)
    ln_mskland  = .false.   !  mask land points in NetCDF outputs (costly: + ~15%)
    ln_cfmeta   = .false.   !  output additional data to netCDF files required for compliance with the CF metadata standard
    ln_clobber  = .true.    !  clobber (overwrite) an existing file
    nn_chunksz  =       0   !  chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)

ln_rstart needs setting to true. The cn_ocerst_in needs changing to the restart file name (in this case 12677 time steps into the simulation) This also assumes that the restart file has been created in the proper time period by running the model with nn_stock set to output an restart file. (in this case at time step 12677)

Once a set of restarts have been created, they can be combined into one file (a restart file is produced for all cores) this requires the rebuild nemo tools to be compiled and used. Please follow the guide "rebuild and inspect NEMO output" recipe. The tool is compiled using the following command::

    cd $TDIR
    ./maketools -m XC_ARCHER_INTEL -n REBUILD_NEMO
    
then to rebuild::
    
    cd $EXP
    export nproc=93
    $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 $CONFIG_000012677_restart_tide $nproc

The 93 indivdual files can now be deleted.

Finally to run the restart file the starting timesteps for the simulation and harmonic analysis need to be changed in the namelist_cfg as before::

    !-----------------------------------------------------------------------
    &namrun        !   parameters of the run
    !-----------------------------------------------------------------------
    cn_exp      =    "SWPacific"  !  experience name
    nn_it000    =  12678   !  first time step
    nn_itend    =  16800 ! 10day=14400   !  last  time step (std 5475)
    nn_date0    =  20000101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
    nn_time0    =       0   !  initial time of day in hhmm
    nn_leapy    =       1   !  Leap year calendar (1) or not (0)
    ln_rstart   = .true.    !  start from rest (F) or from a restart file (T)

As the restart file in this case occurs at time step 12677, the starting time step is this plus one (12678). The start date should not need changing (not 100% sure!)
The strating point of the harmonic analysis also needs changing::

    !-----------------------------------------------------------------------
    &nam_diaharm   !   Harmonic analysis of tidal constituents               ("key_diaharm")
    !-----------------------------------------------------------------------
    nit000_han = 12678         ! First time step used for harmonic analysis
    nitend_han = 16800 ! 7200      ! Last time step used for harmonic analysis
    nstep_han  = 1        ! Time step frequency for harmonic analysis

As the model is a restart and not starting from rest there is no need to stagger the start of the harmonic analysis so it can start at the same time as the simulation restart timestep plus one (12678).

Currently the model has errors when trying to run from restart. It produced an E R R O R saying tmask, vmask and umask variable among others were missing. This was identified to be potentially due to ln_mask_file being set to true::

    !-----------------------------------------------------------------------
    &nambdy        !  unstructured open boundaries
    !-----------------------------------------------------------------------
    ln_bdy         = .true.              !  Use unstructured open boundaries
    nb_bdy         = 1                    !  number of open boundary sets
    ln_coords_file = .true.               !  =T : read bdy coordinates from file
    cn_coords_file = 'coordinates.bdy.nc' !  bdy coordinates files
    ln_mask_file   = .true.              !  =T : read mask from file
    cn_mask_file   = 'bdydta/SWPacific_bdytide_rotT_M2_grid_T.nc'    !  name of mask file (if ln_mask_file=.TRUE.)
    cn_dyn2d       = 'flather'

However, setting this to false results in the model crashing with ssh in excess of 10 m. (Jeff.......!!!)

















