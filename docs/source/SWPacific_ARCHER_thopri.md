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

1. Add more harmonics to both file_def_nemo.xml and namelist_cfg (in terms of output just m2 and s2 for now until accuracy improves.)
2. Check model output against M2 and S2 from FES and reality to see how model is doing.
3. Find out list of HC that are able to be computed within NEMO, compare with FES list
4. Increase length of simulation (ultimatly up to one year)
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

The additional lines with S2 in add the S2 harmonic real and imaginary parts for both elevation and u and v components of the current. This will need to be repeated for all harmonics. I also had to extend the wall time of the runscript as 20 mins was just sliightly too short. This means I have been running on the main queue of ARCHER. Currently takes just under 23 mins.





















