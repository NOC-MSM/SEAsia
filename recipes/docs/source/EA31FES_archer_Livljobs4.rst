==========================================
Setting up EAfrica NEMO v4 configuration
==========================================

Machines: livljobs4, ARCHER

Notes 

The instructions presented here describe the steps I did to produce a regional model with NEMO v4 (@r8395) during my internship 
(July - October 2018). I use the case study of the East Africa coast (for the ACCORD project). I followed/adapted the method developed 
by Jeff for the SEAsia configuration.

I have a working version located here : /work/n01/n01/valegu/EA31FES/trunk_NEMOGCM_r8395/CONFIG/EA31FES/EXP_FullOcean

The outputs of the run are here : /work/n01/n01/valegu/EA31FES/trunk_NEMOGCM_r8395/CONFIG/EA31FES/EXP_FullOcean/OUTPUTS/FullOcean_rivers

Plan 

1. Download the code 

1.a) Creation of a temporary file with all the path names

1.b) Collect the essential files

1.c) Build XIOS @ r1080

1.d) Build NEMO (ORCHESTRA) trunk @r8395

1.e) Build tools

2. Generate new coordinates file (coordinates.nc)

2.a) Choose a domain

2.b) Create coordinates.nc

3. Generate bathymetry file

3.a) Parent bathymetry
 
3.b) Generate bathy_meter.nc

4. Generate domain file (domain_cfg.nc)

4.a) Configure namelist_cfg

4.b) Build a script to run the executable

5. Generate initial conditions

5.a) Building T,S field initial condition from existing fields

5.b) Rough cut some initial conditions from parent (global) dataset

5.c) Use SOSIE tools to flood fill the parent initial conditions

5.d) Use SCRIP tools to remap to the new grid

5.e) Create a land mask

5.f) Interpolate in z on the fly

6. Generation of the atmopsheric forcings

6.a) Where to find data 

6.b) Implementation of the atmospheric forcings

7. Implementation of the river forcings

7.a) Turn on the rivers in NEMO (ln_rnf = true)

8. Generate boundary conditions for tides (FES)

8.a) Install Pynemo (NRCT)

8.b) Generate the ncml file that points to the external data

8.c) Generate the namelist.bdy file for Pynemo (NRCT)

9. Generation of the 3D boundary conditions (LIVLJOBS4)
  
9.a) Parent dataset

9.b) Setting up Pynemo for ORCA0083 N06

9.c) Namelist.bdy set up

9.d) James notes

---

For any questions or any feedback, please feel free to contact me at (valegu@noc.ac.uk)


Recipe Notes
============

The following instructions are specific for a run using the FES tides and 31 levels with pure sigma coordinates.

1. Download the code

Uses a prerelease of NEMO v4 (@r8395). 
  
1.a) Creation of a temporary file with all the path names

Starting on ARCHER::

  ssh valegu@login.archer.ac.uk
  cat > ~/temporary_path_names_for_NEMO_EA31FES_build << EOL    
    export CONFIG=EA31FES
    export WORK=/work/n01/n01
    export WDIR=$WORK/\$USER/\$CONFIG
    export INPUTS=$WDIR/INPUTS
    export START_FILES=$WDIR/START_FILES
    export CDIR=$WDIR/trunk_NEMOGCM_r8395/CONFIG
    export TDIR=$WDIR/trunk_NEMOGCM_r8395/TOOLS
    export EXP=$CDIR/\$CONFIG/EXP00

  module swap PrgEnv-cray PrgEnv-intel
  module unload cray-netcdf 
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel
  EOL
  
  . ~/temporary_path_names_for_NEMO_EA31FES_build

I had to unload cray-netcdf (Jeff's instruction didn't need this step) before loading hdf5parallel and hdf5-parallel. 

---

1.b) Collect the essential files

Prior to get the code, we have to create the directories and grab some important files::

  mkdir $WDIR
  mkdir $INPUTS
  mkdir $START_FILES

  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/bdyini.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/diaharm_fast.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/sbctide.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/step_oce.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_FES14.h90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tideini.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/tide_mod.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/dommsk.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/usrdef_istate.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/nibrun/NEMO/NEMO_trunk_9395/NEMOGCM/CONFIG/SWPacific/MY_SRC/usrdef_sbc.F90 $START_FILES/.
  
As I am using FES tides, I went to have a look at Nico's implementation in the following file (GitHub): "FES2014_NEMO.rst". 
So, all the above files were taken from his repository. 

I took par_oce.F90, dtatsd.F90, coordinated_ORCA_R12.nc and namelist_reshape_bilin_gebco from the same place that Jeff's did::
  
  rsync -uvt /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/par_oce.F90 $START_FILES/.
  rsync -uvt /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/dtatsd.F90 $START_FILES/.

  cp $WORK/jelt/LBay/START_FILES/coordinates_ORCA_R12.nc $START_FILES/.
  cp $WORK/jelt/LBay/INPUTS/namelist_reshape_bilin_gebco $START_FILES/.

---

1.c) Build XIOS @ r1080

Gaby told me it was too long to try to build it myself and it's preferably to take it from someone else. 
This is why I took it from Jeff's repositories::

  cd /work/n01/n01/jelt
  scp -rp xios-2.0_r1080 /work/n01/n01/valegu/

I will copy xios_server.exe into my EXP directory only later, once the code will be downloaded.   

---

1.d) Build NEMO (ORCHESTRA) trunk @r8395

If you don't have a NEMO account, then you need to register here : http://forge.ipsl.jussieu.fr/nemo/register
Now it's time to download the code::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395
  cp $WORK/jelt/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/. 
  
I used the arch file from Jeff. Make sure that the line %XIOS_HOME point to Jeff's directory where XIOS is located. 
I added extra information on the FLAGS (as suggested by Sarah), which is important in the development phase of a new configuration, 
because it will help detect more easily where the error comes from (if there is one) ::

  vi arch-XC_ARCHER_INTEL.fcm
  
  # USER_INC    complete list of include files
  # USER_LIB    complete list of libraries to pass to the linker
  # CC          C compiler used to compile conv for AGRIF
  # CFLAGS      compiler flags used with CC
  #
  # Note that:
  #  - unix variables "$..." are accpeted and will be evaluated before calling fcm.
  #  - fcm variables are starting with a % (and not a $)
  #
  %NCDF_HOME           $NETCDF_DIR
  %HDF5_HOME           $HDF5_DIR
  %XIOS_HOME           /work/n01/n01/jelt/XIOS
  #OASIS_HOME

  %NCDF_INC            -I%NCDF_HOME/include -I%HDF5_HOME/include
  %NCDF_LIB            -L%HDF5_HOME/lib -L%NCDF_HOME/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios
  #OASIS_INC           -I%OASIS_HOME/build/lib/mct -I%OASIS_HOME/build/lib/psmile.MPI1
  #OASIS_LIB           -L%OASIS_HOME/lib -lpsmile.MPI1 -lmct -lmpeu -lscrip

  %CPP                 cpp
  %FC                  ftn
  %FCFLAGS             -integer-size 32 -real-size 64 -g -O3 -fp-model source -zero -fpp -warn all -traceback
  %FFLAGS              -integer-size 32 -real-size 64 -g -O3 -fp-model source -zero -fpp -warn all -funroll-loops -safe-cray-ptr -free -xHost -fp-model source -traceback
  %LD                  CC -Wl,"--allow-multiple-definition"
  %FPPFLAGS            -P -C -traditional
  %LDFLAGS
  %AR                  ar
  %ARFLAGS             -r
  %MK                  gmake
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB
  #USER_INC            %XIOS_INC %OASIS_INC %NCDF_INC
  #USER_LIB            %XIOS_LIB %OASIS_LIB %NCDF_LIB

  %CC                  cc
  %CFLAGS              -O0

Make a new config directory structure (only say YES to OPA_SRC, unless you have other plans) ::
  
  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 clean
  
Edit the CPP flags (USE key_diaharm_fast instead of key_harm_ana for FES tides) ::  
  
  vi $CONFIG/cpp_$CONFIG.fcm
  bld::tool::fppkeys key_zdfgls        \
                   key_FES14_tides   \
                   key_diaharm_fast  \
                   key_mpp_mpi       \
                   key_iomput        \
                   key_nosignedzero

Then here put all the important files inside MY_SRC. This has to be done before compiling ::

  cp $START_FILES/dommsk.F90  $CDIR/$CONFIG/MY_SRC/.

  cp $START_FILES/bdyini.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tideini.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tide_mod.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/tide_FES14.h90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step_oce.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/step.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/sbctide.F90  $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/diaharm_fast.F90    $CDIR/$CONFIG/MY_SRC/.

  cp $START_FILES/diaharm.F90    $CDIR/$CONFIG/MY_SRC/.   
  
Do I need this file (diaharm.F90)? I don't know... but I still put it. 

If you don't want to use idealized and constant initial conditions, you have to change the parameter ln_usr and put the value F 
inside the namelist_cfg inside the working directory $EXP. It doesn't matter if the files (usrdef_istate.F90 and usrdef_sbc.F90) 
are still inside MY_SRC when putting ln_usr = F ::

  cp $START_FILES/usrdef_istate.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/usrdef_sbc.F90 $CDIR/$CONFIG/MY_SRC/.

Add a fix to permit vertical interpolation on-the-fly from initial conditions fields onto child grid ::
  
  cp $START_FILES/par_oce.F90 $CDIR/$CONFIG/MY_SRC/.
  cp $START_FILES/dtatsd.F90  $CDIR/$CONFIG/MY_SRC/.

Copy the xios executable into $EXP directory. I didn't do the symbolic links (ln -s), just cp ::
  
  cd /work/n01/n01/valegu/xios-2.0_r1080/bin/
  cp xios_server.exe $EXP/xios_server.exe

Build opa ::
  
  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10 

After compiling, you should now have the executable inside your BLD/bin folder !

---

1.e) Build tools

To generate domain_cfg and rebuild tools we first need to compile some of the NEMO TOOLS.
Note : DOMAINcfg has to be compiled with XIOS1. There is a README in the $TDIR/DOMAINcfg on what to do.
First build DOMAINcfg (which is relatively new in NEMOv4). Use my XIOS1 file (see userid and path in variable ``%XIOS_HOME``). 
Copy from ARCH *store*::

  cp $WORK/jelt/ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm $CDIR/../ARCH/.
  vi ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm

  # compiler options for Archer CRAY XC-30 (using intel compiler)
  #
  # NCDF_HOME   root directory containing lib and include subdirectories for netcdf4
  # HDF5_HOME   root directory containing lib and include subdirectories for HDF5
  # XIOS_HOME   root directory containing lib for XIOS
  # OASIS_HOME  root directory containing lib for OASIS
  #
  # NCDF_INC    netcdf4 include file
  # NCDF_LIB    netcdf4 library
  # XIOS_INC    xios include file    (taken into accound only if key_iomput is activated)
  # XIOS_LIB    xios library         (taken into accound only if key_iomput is activated)
  # OASIS_INC   oasis include file   (taken into accound only if key_oasis3 is activated)
  # OASIS_LIB   oasis library        (taken into accound only if key_oasis3 is activated)
  #
  # FC          Fortran compiler command
  # FCFLAGS     Fortran compiler flags
  # FFLAGS      Fortran 77 compiler flags
  # LD          linker
  # LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries
  # FPPFLAGS    pre-processing flags
  # AR          assembler
  # ARFLAGS     assembler flags
  # MK          make
  # USER_INC    complete list of include files
  # USER_LIB    complete list of libraries to pass to the linker
  # CC          C compiler used to compile conv for AGRIF
  # CFLAGS      compiler flags used with CC
  #
  # Note that:
  #  - unix variables "$..." are accpeted and will be evaluated before calling fcm.
  #  - fcm variables are starting with a % (and not a $)
  #
  %NCDF_HOME           $NETCDF_DIR
  %HDF5_HOME           $HDF5_DIR
  %XIOS_HOME           /work/n01/n01/jelt/xios-1.0_r703
  #OASIS_HOME

  %NCDF_INC            -I%NCDF_HOME/include -I%HDF5_HOME/include
  %NCDF_LIB            -L%HDF5_HOME/lib -L%NCDF_HOME/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios
  #OASIS_INC           -I%OASIS_HOME/build/lib/mct -I%OASIS_HOME/build/lib/psmile.MPI1
  #OASIS_LIB           -L%OASIS_HOME/lib -lpsmile.MPI1 -lmct -lmpeu -lscrip

  %CPP                 cpp
  %FC                  ftn
  %FCFLAGS             -integer-size 32 -real-size 64 -g -O3 -fp-model source -zero -fpp -warn all
  %FFLAGS              -integer-size 32 -real-size 64 -g -O3 -fp-model source -zero -fpp -warn all
  %LD                  CC -Wl,"--allow-multiple-definition"
  %FPPFLAGS            -P -C -traditional
  %LDFLAGS
  %AR                  ar
  %ARFLAGS             -r
  %MK                  gmake
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB
  #USER_INC            %XIOS_INC %OASIS_INC %NCDF_INC
  #USER_LIB            %XIOS_LIB %OASIS_LIB %NCDF_LIB

  %CC                  cc
  %CFLAGS              -O0

Make sure that the XIOS_HOME is updated with Jeff's path.

Then build some tools ::

  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n DOMAINcfg
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n REBUILD_NEMO


For the generation of bathymetry and met forcing weights files we need to patch the code (to allow direct passing of arguments. 
NB this code has not been updated in 7 years.).
Before, the patching, we need to get the files, I took them from Jeff's repositories ::

  cd /work/n01/n01/jelt/SEAsia/START_FILES/
  cp scripinterp_mod.patch $START_FILES/.
  cp scripinterp.patch $START_FILES/.
  cp scrip.patch $START_FILES/.
  cp scripshape.patch $START_FILES/.
  cp scripgrid.patch $START_FILES/.

  cd $TDIR/WEIGHTS/src
  patch -b < $START_FILES/scripinterp_mod.patch
  patch -b < $START_FILES/scripinterp.patch
  patch -b < $START_FILES/scrip.patch
  patch -b < $START_FILES/scripshape.patch
  patch -b < $START_FILES/scripgrid.patch

  cd $TDIR
  ./maketools -m XC_ARCHER_INTEL_XIOS1 -n WEIGHTS
  
---

2. Generate new coordinates file (coordinates.nc)

In this section, we aim at generating the file named "coordinates.nc" from a parent NEMO grid at some resolution 
(here we take the parent file named "coordinates_ORCA_R12.nc" which was collected during step 1.b))

**Plan:** Use tool ``agrif_create_coordinates.exe`` which reads cutting indices and parent grid location from ``namelist.input``
and outputs a new files with new resolution grid elements. 

2.a) Choose a domain

First we need to figure out the indices for the new domain from the parent grid. Move parent grid into INPUTS ::

  cp $START_FILES/coordinates_ORCA_R12.nc $INPUTS/.
  
Inspect this parent coordinates file to define the boundary indices for the new configuration. You can do it in the way you feel more 
comfortable. Jeff did it on Ferret, I used Matlab. To define the boundaries, it is worth looking at the tides. Jeff looked at the TPXO 
harmonic amplitudes to find a good cut off location for boundaries. As I was reproducing some of Gaby's work, I didn't do it and suppose 
that taking Gaby's previous box was good enough. I just extended a little bit south to take into account a big river.  

---

2.b) Create coordinates.nc

Once you know at which indexes you want to cut the file coordinates_ORCA_R12.nc, you have to edit the namelist. I changed imin, imax, 
jmin, jmax so that it delimits the East African area of interest.  Rho & rhot are a scaling factor. If I want a 12th degree resolution,
let the value of 1. If I want 36th degree of resolution, I have to put 3. If I want 60th degree resolution, I have to put 5 ::

  cd $TDIR/NESTING
  vi namelist.input

  &input_output
      iom_activated = true
  /
  &coarse_grid_files
      parent_coordinate_file = 'coordinates_ORCA_R12.nc'
  /
  &bathymetry
  /
  &nesting
      imin = 3907
      imax = 3962
      jmin = 1356
      jmax = 1478
      rho  = 5
      rhot = 5
      bathy_update = false
  /
  &vertical_grid
  /
  &partial_cells
  /
  &nemo_coarse_grid
  /
  &forcing_files
  /
  &interp
  /
  &restart
  /
  &restart_trc
  /

Move to the TOOL directory ::

  cd $TDIR
  
Copy in the right ARCH file ::

  cp /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL_NOXIOS.fcm ../ARCH/.
  
Build the NESTING tool ::

  ./maketools -n NESTING -m XC_ARCHER_INTEL_NOXIOS -j 6
  
This makes a number of executables in NESTING. (I did this step on a previous configuration when trying to do the hybrid coordinates, 
so I didn't do it again on the configuration EA31FES. This is why I don't have the executables. I just directly copied the coordinates.nc. 
Go to NESTING ::

  cd NESTING

Link in parent coordinates file ::

  ln -s $START_FILES/coordinates_ORCA_R12.nc $TDIR/NESTING/.
  
Execute tool ::
 
  ./agrif_create_coordinates.exe 
  
This creates a coordinate file. Copy it to the $INPUTS directory ::
 
  cp 1_coordinates_ORCA_R12.nc $INPUTS/coordinates.nc
   
   
---

3. Generate bathymetry file

Once we have the coordinates.nc file, we need to generate the associated bathymetry.

3.a) Parent bathymetry

Go to the following link : https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data
You need to create an account and put a request for the access of the data you want to download. You get the response very quickly. 
We want the dataset to be spatially larger than the desired domain. 

I selected the 30 arc second (and not the one minute resolution grid), and entered the following box coordinates : 33.0000,-19.2000,59.0000,9.6000.

The file I obtain is called GEBCO_2014_2D_33.0_-19.2_59.0_9.6.nc. Put your file on archer. 
Then copy it into your INPUTS directory ::

  scp -rp GEBCO_2014_2D_33.0_-19.2_59.0_9.6.nc /work/n01/n01/$USER/EA31FES/INPUTS/

  ---
  
3.b) Generate bathy_meter.nc
  
Copy over namelist for reshaping bathymetry ::
 
  cp $START_FILES/namelist_reshape_bilin_gebco $INPUTS/.
 
Edit namelist to point to correct input file. Edit lat and lon variable names to make sure they match the nc file content
(GEBCO_2014_2D_33.0_-19.2_59.0_9.6.nc for me). You can use "ncdump -h" to get the input variable names ::

  ncdump -h GEBCO_2014_2D_33.0_-19.2_59.0_9.6.nc 

  vi $INPUTS/namelist_reshape_bilin_gebco
    ...
  /
  &grid_inputs
      input_file = 'gebco_in.nc'
      nemo_file = 'coordinates.nc'
      datagrid_file = 'remap_data_grid_gebco.nc'
      nemogrid_file = 'remap_nemo_grid_gebco.nc'
      method = 'regular'
      input_lon = 'lon'
      input_lat = 'lat'
      nemo_lon = 'glamt'
      nemo_lat = 'gphit'
      nemo_mask = 'none'
      nemo_mask_value =  0
      input_mask = 'none'
      input_mask_value = 0
  /
  ...

  &interp_inputs
      input_file = "gebco_in.nc"
      interp_file = "data_nemo_bilin_gebco.nc"
      input_name = "elevation"
      input_start = 1,1
      input_stride = 1,1
      input_stop = 0,0
  /

  &interp_outputs
      output_file = "bathy_meter.nc"
      output_mode = "create"
      output_dims = 'x', 'y'
      output_scaling = "topo|1.0"
      output_name = 'Bathymetry'
      output_lon = 'nav_lon'
      output_lat = 'nav_lat'
  /
  ...

Now we need to flatten out the land elevations and make the depths positive ::
  
  cd $INPUTS
  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0 
  ncap2 -s 'where(elevation > 0) elevation=0' GEBCO_2014_2D_33.0_-19.2_59.0_9.6.nc tmp.nc  
  ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc 
  rm tmp.nc 
  
Restore original modules ::
 
  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

Execute first scrip thing ::
  
$TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco

Output files ::

  remap_nemo_grid_gebco.nc
  remap_data_grid_gebco.nc

Execute second scrip thing ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco

Output files ::

  data_nemo_bilin_gebco.nc
  
 Execute third scrip thing ::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco

Output files ::

  bathy_meter.nc

Use ncview to make sure that it has been well created. 
  
---

4. Generate domain file (domain_cfg.nc)

4.a) Configure namelist_cfg 

The general idea is that you have to copy the ``namelist_cfg`` file into the ``DOMAINcfg`` directory along with all the inputs files 
that would have previously been needed get v3.6 running. The reason being that all the non-time stepping stuff, like grid generating, 
has been abstracted from the core OPA code and is now done as a pre-processing step, and output into an important file ``domain_cfg.nc``.

Copy essential files into DOMAINcfg directory ::

 ln -s $INPUTS/coordinates.nc $TDIR/DOMAINcfg/.
 ln -s $INPUTS/bathy_meter.nc $TDIR/DOMAINcfg/.

Edit the template ``namelist_cfg`` with only the essential domain building stuff. 
Get the size of the new domain from ::

  ncdump -h bathy_meter.nc
  
It took me time to figure out how to set up the namelist. My namelist is located here : /work/n01/n01/valegu/EA31FES/trunk_NEMOGCM_r8395/TOOLS/DOMAINcfg
Here, this is specific for sigma coordinates ::

  cd $TDIR/DOMAINcdf/
  vi namelist_cfg
  
  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
    nn_no       =       0   !  job number (no more used...)
    cn_exp      =  "domaincfg"  !  experience name
    nn_it000    =       1   !  first time step
    nn_itend    =      75   !  last  time step (std 5475)
  /
  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
    !
    ln_e3_dep   = .false.   ! =T : e3=dk[depth] in discret sens.
    !                       !      ===>>> will become the only possibility in v4.0
    !                       ! =F : e3 analytical derivative of depth function
    !                       !      only there for backward compatibility test with v3.6
    !                       !
    cp_cfg      =  "orca"   !  name of the configuration
    jp_cfg      =      60   !  resolution of the configuration
    jpidta      =     279   !  1st lateral dimension ( >= jpi )
    jpjdta      =     614   !  2nd    "         "    ( >= jpj )
    jpkdta      =      31    !  number of levels      ( >= jpk )
    jpiglo      =     279   !  1st dimension of global domain --> i =jpidta
    jpjglo      =     614   !  2nd    -                  -    --> j  =jpjdta
    jpizoom     =       1   !  left bottom (i,j) indices of the zoom
    jpjzoom     =       1   !  in data domain indices
    jperio      =       0   !  lateral cond. type (between 0 and 6)
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
    ln_zco      = .false.   !  z-coordinate - full    steps
    ln_zps      = .false.   !  z-coordinate - partial steps
    ln_sco      = .true.    !  s- or hybrid z-s-coordinate
    ln_isfcav   = .false.   !  ice shelf cavity
    ln_linssh   = .false.   !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
    ln_s_sh94   = .false.  !  Song & Haidvogel 1994 hybrid S-sigma   (T)|
    ln_s_sf12   = .true.   !  Siddorn & Furner 2012 hybrid S-z-sigma (T)| if both are false the NEMO tanh stretching is applied
    ln_sigcrit  = .true.   !  use sigma coordinates below critical depth (T) or Z coordinates (F) for Siddorn & Furner stretch
                           !  stretching coefficients for all functions
    rn_sbot_min =   10.0     !  minimum depth of s-bottom surface (>0) (m )               *** VAL : changed from 6 to 10 as in james namelist ***
    rn_sbot_max =   6000.0   !  maximum depth of s-bottom surface (= ocean depth) (>0) (m) *** VAL : changed from 7000 to 6000 as in james namelist ***
    rn_hc       =   50.0     !  critical depth for transition to stretched coordinates    *** VAL : Not sure of this value, I used the same as in jeff and james namelists ***
           !!!!!!!  Envelop bathymetry
    rn_rmax     =   0.05     !  maximum cut-off r-value allowed (0<r_max<1) *** VAL : changed from 0.05 to 0.3 as in james namelist ***
           !!!!!!!  SH94 stretching coefficients  (ln_s_sh94 = .true.)
    rn_theta    =   6.0     !  surface control parameter (0<=theta<=20) *** VAL : changed from 20 to 6 as in james namelist ***
    rn_bb       =   0.8     !  stretching with SH94 s-sigma
           !!!!!!!  SF12 stretching coefficient  (ln_s_sf12 = .true.) *** VAL : all this block below was added following Sarah idea ***
    rn_alpha    =    4.4    !  stretching with SF12 s-sigma
    rn_efold    =    0.0    !  efold length scale for transition to stretched coord
    rn_zs       =    1.0    !  depth of surface grid box
                            !  bottom cell depth (Zb) is a linear function of water depth Zb = H*a + b
    rn_zb_a     =    0.024  !  bathymetry scaling factor for calculating Zb
    rn_zb_b     =   -0.2    !  offset for calculating Zb
                         !!!!!!!! Other stretching (not SH94 or SF12) [also uses rn_theta above]
    rn_thetb = 1.0          ! bottom control parameter (0<=thetb<= 1)
  /
  !-----------------------------------------------------------------------
  &namdom        !   space and time domain (bathymetry, mesh, timestep)
  !-----------------------------------------------------------------------
    nn_closea   =    1      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)         *** VAL : added line as in james ***
    nn_msh      =    0      !  create (=1) a mesh file or not (=0)                           *** VAL : added line as in james ***
    rn_hmin     =   -10.    !  min depth of the ocean (>0) or min number of ocean level (<0) *** VAL : added line as in james ***
    rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice    *** VAL : added line as in james ***
    rn_e3zps_min=   25.     !  partial step thickness is set larger than the minimum of      *** VAL : added line as in james ***
    rn_e3zps_rat=    0.2    !  rn_e3zps_min and rn_e3zps_rat*e3t, with 0<rn_e3zps_rat<1      *** VAL : added line as in james ***
    rn_rdt      =   300.    !  time step for the dynamics (and tracer if nn_acc=0)           *** VAL : changed from 360 to 300 as in james ***
    jphgr_msh = 0           ! type of horizontal mesh                                        *** VAL : added line from jeff SWPacific
                                       !  type of horizontal mesh
                                       !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
                                       !  = 1 geographical mesh on the sphere with regular grid-spacing
                                       !  = 2 f-plane with regular grid-spacing
                                       !  = 3 beta-plane with regular grid-spacing
                                       !  = 4 Mercator grid with T/U point at the equator
    ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
    ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
    ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
    ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
    ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
    ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
    ppsur       =  999999.              !  ORCA r4, r2 and r05 coefficients   *** VAL : changed from -4762.96143546300 in jeffSWPacific namelist to 999999. ***
    ppa0        =  999999.              ! (default coefficients)              *** VAL : changed from  255.58049070440 in jeffSWPacific namelist to 999999. ***
    ppa1        =  999999.              !                                     *** VAL : changed from  245.58132232490 in jeffSWPacific namelist to 999999. ***
    ppkth       =  21.43336197938       !                                     *** VAL : use  21.43336197938 as in jeff SWPacific namelist ***
    ppacr       =  3.0                  !                                     *** VAL : use  3.0 as in jeff SWPacific namelist ***
    ppdzmin     =  5.                   !  Minimum vertical spacing           *** VAL : Changed from  999999. in jeffSWPacific to 5. ***
    pphmax      =  4000.                !  Maximum depth                      *** VAL : Changed from 999999. in jeffSWPacific to 4000.
    ldbletanh   =  .false.              !  Use/do not use double tanf function for vertical coordinates
    ppa2        =  999999.              !  Double tanh function parameters
    ppkth2      =  999999.              !
    ppacr2      =  999999.
  /
  !-----------------------------------------------------------------------
  &nameos        !   ocean physical parameters
  !-----------------------------------------------------------------------
    ln_teos10   = .true.         !  = Use TEOS-10 equation of state
  /

4.b) Build a script to run the executable

The script to run is named rs ::

  vi $TDIR/DOMAINcdf/rs

  #!/bin/bash
  #PBS -N domain_cfg
  #PBS -l walltime=00:20:00
  #PBS -l select=1
  #PBS -j oe
  #PBS -A n01-ACCORD
  # mail alert at (b)eginning, (e)nd and (a)bortion of execution
  #PBS -m bea
  #PBS -M valegu@noc.ac.uk
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
    aprun -n 1 -N 1 ./make_domain_cfg.exe >&  stdouterr_cfg

    exit

Then try running it ::

  cd $TDIR/DOMAINcfg
  qsub -q short rs

Copy domain_cfg.nc to the EXP directory and into the INPUTS directory. Jeff did an rsync -utc, I just did scp -rp ::

  scp -rp  $TDIR/DOMAINcfg/domain_cfg.nc $EXP/.
  scp -rp  $TDIR/DOMAINcfg/domain_cfg.nc $INPUTS/.

You should check your domain_cfg.nc (Matlab, ncview,...).

---

5. Generate initial conditions

For a new configuration you probably want to start with idealised, or homogenous initial conditions. This is done with user defined 
initial conditions ``ln_usr=T`` with the expression being compiled into the executable. (In ``$CDIR/$CONFIG/MY_SRC``:  ``usrdef_sbc.F90``
and ``usrdef_istate.F90``).

To use initial conditions from an existing T,S field you might need to do a bit of interpolation. It is advisable to let NEMO do the 
heavy lifting for vertical interpolation (requiring some FORTRAN modifications), though SOSIE tools can be user to do simple horizontal
interpolation.

5.a) Building T,S field initial condition from existing fields

My parent file is ORCA0083-N06_20130105d05T.nc, obtained from Gaby and located here : /projectsa/NEMO/gmaya/ORCA12/2013.

---

5.b) Rough cut some initial conditions from parent (global) dataset

Make cut down parent file using ORCA0083-N06. Copy parent file to ARCHER INPUTS.

Livljobs4 ::

  scp /projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130105d05T.nc valegu@login.archer.ac.uk:/work/n01/n01/valegu/EA31FES/INPUTS/.

Archer ::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0
  cd $WDIR/INPUTS
  ncks -d x,3906,4000 -d y,1353,1478 ORCA0083-N06_20130105d05T.nc $WDIR/INPUTS/cut_down_20130105d05_EA31FES_grid_T.nc

Restore the parallel modules ::

  module unload nco cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

---

5.c) Use SOSIE tools to flood fill the parent initial conditions

Interpolating the T,S on z-levels onto hybrid levels can create water where there was previously only land. Convert all the land in
the parent initial conditions to water by flooding the domain. This can be done with the SOSIE tool.
So you need to build the SOSIE tool. Jeff zipped the sosie directory and copied it in my directories. It is named sosie.tar, located on archer
here : /home/n01/n01/valegu/

Otherwise, you can get the sosie folder like that ::
 
  cd ~
  git clone https://github.com/brodeau/sosie.git

Go on your sosie folder and edit the path of the ``make.macro`` file ::

  vi /home/n01/n01/valegu/sosie/make.macro

  # Makefile for SOSIE with Intel Ifort on Linux
  # ============================================

  # Fortran compiler:
  FC = ftn

  # Root directory for Netcdf:
  #NETCDF_DIR = /opt/netcdf_intel11

  # Linking argument: usually -lnetcdf or -lnetcdff (or both):
  L_NCDF = -lnetcdf -lnetcdff

  # Fortran compilation flags:
  # -- Production
  FF = -O3 -i4 -xHost -module mod/
  # -- Debugging
  #FF = -O0 -i4 -traceback -CB -module mod/

  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/valegu/local

Install. This might be best done in a clean terminal ::

  cd ~
  mkdir local
  cd sosie

  make clean
  make
  make install
  export PATH=~/local/bin:$PATH
  cd $WDIR/INPUTS

Obtain the fields to interpolate. E.g interpolate ORCA data. Get the namelists ::
 
  cd /work/n01/n01/jelt/SEAsia/START_FILES
  cp initcd_votemper.namelist /work/n01/n01/valegu/EA31FES/INPUTS/
  cp initcd_vosaline.namelist /work/n01/n01/valegu/EA31FES/INPUTS/
  cd /work/n01/n01/valegu/EA31FES/INPUTS/

Edit namelists to the variables you want ::

  vi initcd_vosaline.namelist 
  ...
  &ninput
  ivect     = 0
  lregin    = F
  cf_in     = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_in     = 'so'
  cv_t_in   = 'time_counter'
  jt1       = 0
  jt2       = 0
  jplev     = 0
  cf_x_in   = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_lon_in = 'nav_lon'
  cv_lat_in = 'nav_lat'
  cf_lsm_in = 'missing_value'
  cv_lsm_in = ''
  ldrown    = T
  ewper     = -1
  vmax      =  1.E6
  vmin      = -1.E6
  ismooth   = 0
  /
  ...
  
  &n3d
  cf_z_in  = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_z_in  = 'deptht'
  cf_z_out = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_z_out = 'deptht'
  cv_z_out_name = 'gdept'
  ctype_z_in = 'z'
  ctype_z_out = 'z'
  ...
  
  &nhtarget
  lregout    = F
  cf_x_out   = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_lon_out = 'nav_lon'
  cv_lat_out = 'nav_lat'
  cf_lsm_out = ''
  cv_lsm_out = ''
  lmout      = F
  !rmaskvalue = -9999
  lct        = F
  t0         = 0.
  t_stp      = 0.
  ewper_out  = -1
  ...
  
  &noutput
  cmethod  = 'bilin'
  cv_t_out = 'time_counter'
  cv_out   = 'vosaline'
  cu_out   = 'psu'
  cln_out  = 'Salinity'
  cd_out   = '.'
  !!
  csource  = 'ORCA0083-N06'
  ctarget  = 'EAfrica'
  cextra   = '2013'

Similarly for initcd_votemper.namelist ::

  ...
  &ninput
  ivect     = 0
  lregin    = F
  cf_in     = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_in     = 'thetao'
  cv_t_in   = 'time_counter'
  jt1       = 0
  jt2       = 0
  jplev     = 0
  cf_x_in   = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_lon_in = 'nav_lon'
  cv_lat_in = 'nav_lat'
  cf_lsm_in = 'missing_value'
  cv_lsm_in = ''
  ldrown    = T
  ewper     = -1
  vmax      =  1.E6
  vmin      = -1.E6
  ismooth   = 0
  ...
  
  &n3d
  cf_z_in  = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_z_in  = 'deptht'
  cf_z_out = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_z_out = 'deptht'
  cv_z_out_name = 'gdept'
  ctype_z_in = 'z'
  ctype_z_out = 'z'
  ...
  
  &nhtarget
  lregout    = F
  cf_x_out   = 'cut_down_20130105d05T_EA31FES_grid_T.nc'
  cv_lon_out = 'nav_lon'
  cv_lat_out = 'nav_lat'
  cf_lsm_out = ''
  cv_lsm_out = ''
  lmout      = F
  !rmaskvalue = -9999
  lct        = F
  t0         = 0.
  t_stp      = 0.
  ewper_out  = -1
  ...
  
  &noutput
  cmethod  = 'bilin'
  cv_t_out = 'time_counter'
  cv_out   = 'votemper'
  cu_out   = 'C'
  cln_out  = 'Temperature'
  cd_out   = '.'
  !!
  csource  = 'ORCA0083-N06'
  ctarget  = 'EAfrica'
  cextra   = '2013'
  /

Then, you can use a PBS submission script. But I didn't used one, I just run from the command line ::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel

  cd /home/n01/n01/valegu/sosie
  make clean
  make
  make install

  cd /work/n01/n01/valegu/EA31FES/INPUTS
  /home/n01/n01/valegu/local/bin/sosie.x -f initcd_votemper.namelist
  /home/n01/n01/valegu/local/bin/sosie.x -f initcd_vosaline.namelist

Whether as a serial job or from the command line, the temperature process creates ::
  
  sosie_mapping_ORCA0083-N06-EAfrica.nc
  votemper_ORCA0083-N06-EAfrica_2013.nc

And the salinity process creates ::

  vosaline_ORCA0083-N06-EAfrica_2013.nc

Check these fields are OK.

---

5.d) Use SCRIP tools to remap to the new grid

The scrip tools are build in TDIR during the step explained in 1.e). Now do interpolation onto child lateral grid. First copy 
the namelists ::

  cd work/n01/n01/jelt/SEAsia/INPUTS
  cp namelist_reshape_bilin_initcd_votemper /work/n01/n01/valegu/EA31FES/INPUTS/
  cp namelist_reshape_bilin_initcd_vosaline /work/n01/n01/valegu/EA31FES/INPUTS/

Edit the input files ::

  vi $INPUTS/namelist_reshape_bilin_initcd_votemper
  ...
  &grid_inputs
      input_file = 'votemper_ORCA0083-N06-EAfrica_2013.nc'
      nemo_file = 'coordinates.nc'
      datagrid_file = 'remap_data_grid_R12.nc'
      nemogrid_file = 'remap_nemo_grid_R12.nc'
      method = 'regular'
      input_lon = 'nav_lon'
      input_lat = 'nav_lat'
      nemo_lon = 'glamt'
      nemo_lat = 'gphit'
      nemo_mask = 'none'
      nemo_mask_value =  0
      input_mask = 'none'
      input_mask_value = 0
  ...
  &interp_inputs
      input_file = "votemper_ORCA0083-N06-EAfrica_2013.nc"
      interp_file = "data_nemo_bilin_R12.nc"
      input_name = "votemper"
      input_start = 1,1,1,1
      input_stride = 1,1,1,1
      input_stop = 0,0,0,0
      input_vars = "gdept","time_counter"

Similarly for the *vosaline.nc file ::

  &grid_inputs
      input_file = 'vosaline_ORCA0083-N06-EAfrica_2013.nc'
      nemo_file = 'coordinates.nc'
      datagrid_file = 'remap_data_grid_R12.nc'
      nemogrid_file = 'remap_nemo_grid_R12.nc'
      method = 'regular'
      input_lon = 'nav_lon'
      input_lat = 'nav_lat'
      nemo_lon = 'glamt'
      nemo_lat = 'gphit'
      nemo_mask = 'none'
      nemo_mask_value =  0
      input_mask = 'none'
      input_mask_value = 0
  /
  ...

  &interp_inputs
      input_file = "vosaline_ORCA0083-N06-EAfrica_2013.nc"
      interp_file = "data_nemo_bilin_R12.nc"
      input_name = "vosaline"
      input_start = 1,1,1,1
      input_stride = 1,1,1,1
      input_stop = 0,0,0,0
      input_vars = "gdept","time_counter"
  /
  
Produce the remap files ::

  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates remap_nemo_grid_R12.nc and remap_data_grid_R12.nc. Then ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates data_nemo_bilin_R12.nc. Then ::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates initcd_votemper.nc. Then ::

  $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates initcd_vosaline.nc.

---

5.e) Create a land mask

You can create a land mask to tell SOSIE what needs flooding. Use the salinity field to do this since we know the salinity field is zero on land. 
Then you can specify the name of the file inside the parameter called  ``cf_lsm_in`` inside the file ``initcd_vosaline.namelist`` and 
``initcd_votemper.namelist``. But, it didn't have to use it as I put inside ``cf_lsm_in`` missing value.

If you want to create the mask, do this ::
 
 module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  ncks -d time_counter,0,0,1 -v vosaline initcd_vosaline.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline <=30.) vosaline=0' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncap2 -O -s 'where(vosaline >0.) vosaline=1' sosie_initcd_mask.nc sosie_initcd_mask.nc
  ncrename -v vosaline,mask sosie_initcd_mask.nc

Restore modules ::

  module unload nco/4.5.0
  module unload cray-netcdf cray-hdf5
  module load cray-netcdf-hdf5parallel cray-hdf5-parallel

This has created a file initcd_mask with a variable mask.

---

5.f) Interpolate in z on the fly

For vertical interpolation we let NEMO do the heavy lifting. This requires some changes to the FORTRAN using par_oce.F90 and dtatsd.F90 
in MY_SRC. To interpolate the initial conditions on-the-fly need to pass information to NEMO about the parent vertical grid and parent
mask file. Appropriate variables are created in external files that are read into the namelist.

These mask and depth variables need to be 4D variables, where length(t)=1. They can be created with NCO tools by manipulating a parent 
initial condition file. On archer, load the appropriate modules ::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

Jeff's configuration use the same number of levels as the parent file (75 levels). My configuration use only 31 levels and the 
parent file has 75 levels. I tried to create initcd_depth.nc with 75 levels and then interpolate the tracers (temperature
and salinity) on the fly by NEMO using the logical switch to do vertical interpolation ln_tsd_interp=T (namelist_cfg).
However, it didn't work (NEMO kept saying 75 is different from 31) so I did the interpolation of initcd_votemper.nc and initcd_vosaline
(both having 75 levels) into the depth of my domain_cfg.nc configuration (31 levels) using Matlab.

I obtained initcd_votemper_interp.nc and initcd_vosaline_interp.nc. 

If the depth (gdept) variable is 1D and the file has dimensions [time,z,y,x] then first we make it 3D and call it something like 
gdept_3D ::

  cd $INPUTS/
  ncap2 -O -s gdept_3D[z,y,x]=gdept initcd_votemper_interp.nc tmp.nc

Then add a time dimension ::
  
  ncap2 -O -s gdept_3D[time_counter,z,y,x]=gdept_3D tmp.nc initcd_depth_3.nc
  rm tmp.nc

The file initcd_depth_3.nc contains the depth values of my configuration (31 levels). It worked when the file provided inside namelist_cfg
(``sn_dep``) had the same number of levels as the domain_cfg.nc. Maybe re-visit this part later on if I have time, to try the
interpolation on the fly directly through NEMO (not doing it myself through Matlab). 

---

6. Generation of the atmopsheric forcings

6.a) Where to find data 

I first went to have a look inside the folder (ARCHER) : /work/n01/n01/acc/ORCA0083/NEMOGCM/CONFIG/R12_ORCA/EXP00/FORCING

Jeff took his atmospheric data from there (DFS5 data). I needed datasets for 2013 but there was only from 1979 to 2011. 
Also Gaby suggested high frequency data (hourly if possible), so I used Era5.
Another place you can have a look for atmospheric data for Era5 is here (LIVLJOBS4): /projectsa/NEMO/Forcing/ERA5/INST/

I extracted Era5 datasets, the steps are explained here : https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5+data+via+the+ECMWF+Web+API

Take care, the specific humidity has to be computed. I found on Era5 documentation the following : 
https://confluence.ecmwf.int/display/CKB/ERA+datasets%3A+near-surface+humidity

So I computed the specific humidity at 2m from the 2m dew point temperature and the 2m surface pressure and created the file 
EAfrica_q_calculated.nc.

---

6.b) Implementation of the atmospheric forcings

Once you have your data, you need to put them on your INPUT directory in ARCHER.
I didn't have to cut my atmospheric files using ncea because this was already done during the extraction through the Python script.
Obtain namelist files and data file ::

  cp /work/n01/n01/jelt/SEAsia/START_FILES/namelist_reshape_bilin_atmos $INPUTS/.
  cp /work/n01/n01/jelt/SEAsia/START_FILES/namelist_reshape_bicubic_atmos $INPUTS/.

Edit namelist to reflect source filenames ::

  vi $WDIR/INPUTS/namelist_reshape_bilin_atmos
  ...
  &grid_inputs
      input_file = 'cutdown_EAfrica_mtpr_y2013.nc'
      nemo_file = 'coordinates.nc'
      datagrid_file = 'remap_data_grid_atmos.nc'
      nemogrid_file = 'remap_nemo_grid_atmos.nc'
      method = 'regular'
      input_lon = 'longitude'
      input_lat = 'latitude'
      nemo_lon = 'glamt'
      nemo_lat = 'gphit'
      nemo_mask = 'none'
      nemo_mask_value =  0
      input_mask = 'none'
      input_mask_value = 0

  vi $WDIR/INPUTS/namelist_reshape_bicubic_atmos
    ...
    &grid_inputs
      input_file = 'cutdown_EAfrica_mtpr_y2013.nc'
      nemo_file = 'coordinates.nc'
      datagrid_file = 'remap_data_grid_atmos.nc'
      nemogrid_file = 'remap_nemo_grid_atmos.nc'
      method = 'regular'
      input_lon = 'longitude'
      input_lat = 'latitude'
      nemo_lon = 'glamt'
      nemo_lat = 'gphit'
      nemo_mask = 'none'
      nemo_mask_value =  0
      input_mask = 'none'
      input_mask_value = 0

Setup weights files for the atmospheric forcing. Use the pre-compiled tools ::

  cd $INPUTS
  $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos

Generate  remap files ``remap_nemo_grid_atmos.nc`` and ``remap_data_grid_atmos.nc``. Then ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos

Generates ``data_nemo_bilin_atmos.nc``. Then ::

  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos

Generates ``weights_bilinear_atmos.nc``. Then ::

  $TDIR/WEIGHTS/scrip.exe namelist_reshape_bicubic_atmos
  
Generates ``data_nemo_bicubic_atmos.nc``. Then ::

  $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bicubic_atmos
  
Generates ``weights_bicubic_atmos.nc``. 

---

7. Implementation of the river forcings

Sarah confirmed that the river runoff has to be in kg/m2/s. She provided a sample code to use to create the river forcing file. 
I used Matlab to generate the river forcings netcdf file for my 6 rivers along the EAfrica Coast ::

  Sample code:
  % to convert from flow in m3/s to kg/m2/s
  % density = 1000 kg/m3
  % rate = flow * density / area of  grid box  in kg/m2/s
  %for NEMO
  dx=ncread('coordinates_CS15r.nc','e1t');
  dy=ncread('coordinates_CS15r.nc','e2t');
  nemo_area=dx.*dy;

7.a) Turn on the rivers in NEMO (ln_rnf = true)

I changed the method to implement the rivers (because initially I had a huge increase of temperature at the river mouth).
The main change is ln_rnf_depth_ini (put to true instead of initially false). Important documentation that made me consider it : 
https://www.nemo-ocean.eu/doc/node53.html

Setup to run with rivers is as follow ::

  vi namelist_cfg
  &namsbc_rnf    !   runoffs namelist surface boundary condition          (ln_rnf=T)
  !-----------------------------------------------------------------------
  !              !  file name           ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
  !              !                      !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
     sn_rnf      = 'EAfrica_rivers',        -1         , 'rorunoff',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_cnf      = 'EAfrica_rivers',         0         , 'socoefr0',   .false.    , .true. , 'yearly'  , ''       , ''       , ''
     sn_s_rnf    = 'runoffs'            ,        24         , 'rosaline',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_t_rnf    = 'runoffs'            ,        24         , 'rotemper',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
     sn_dep_rnf  = 'runoffs'            ,         0         , 'rodepth' ,   .false.    , .true. , 'yearly'  , ''       , ''       , ''

     cn_dir      = 'bdydta/'      !  root directory for the location of the runoff files
     ln_rnf_mouth= .false.        !  specific treatment at rivers mouths ***val : put false, initially true ***
        rn_hrnf     =  15.e0      !  depth over which enhanced vertical mixing is used    (ln_rnf_mouth=T)
        rn_avt_rnf  =  1.e-3      !  value of the additional vertical mixing coef. [m2/s] (ln_rnf_mouth=T)
     rn_rfact    =   1.e0         !  multiplicative factor for runoff
     ln_rnf_depth= .false.        !  read in depth information for runoff
     ln_rnf_tem  = .false.        !  read in temperature information for runoff
     ln_rnf_sal  = .false.        !  read in salinity information for runoff
     ln_rnf_depth_ini = .true.    !  compute depth at initialisation from runoff file ***val : I need to activate that so that the rivers are not boiling anymore***
        rn_rnf_max  = 1.6224      !  5.735e-4    !  max value of the runoff climatologie over global domain ( ln_rnf_depth_ini = .true )
        rn_dep_max  = 50.         !  150.        !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )
        nn_rnf_depth_file = 0     !  create (=1) a runoff depth file or not (=0)
  /

When I turned ``ln_rnf_depth_ini`` to true, I had to put my EAfrica_rivers.nc file inside my running directory. Putting it inside
/bdydta was not enough.

---

8. Generate boundary conditions for tides (FES)

8.a) Install Pynemo (NRCT)

LIVLJOBS4

I created a new environment (named "new_env") that specifies the version of numpy and matlplotlib.
Set up pynemo with the Generalise-tide-input branch. OnGitHub, Jeff and James have a development branch called Generalise-tide-input 
that I have to use in order to create the tides boundaries ::

  ssh -X livljobs4

  export CONFIG=EA31FES_tide
  export WORK=/work
  export WDIR=$WORK/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES

  cd $WORK/$USER
  mkdir NEMO
  cd NEMO
  mkdir EA31FES_tide

  cd /work/valegu/NEMO/EA31FES_tide/INPUTS

  module load anaconda/2.1.0  # Want python2

  conda create --name new_env scipy=0.16.0 numpy=1.9.2 matplotlib=1.4.3 basemap netcdf4 libgfortran=1.0.0 # Another way that works (Jeff did it) is to not specify the version number : conda create --name new_env scipy=0.16.0 numpy matplotlib  basemap netcdf4 libgfortran=1.0.0

  source activate new_env
 
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius
  
  cd $WORK/$USER
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd $WORK/$USER/nrct/Python
  git checkout Generalise-tide-input # Go inside the Generalise-tide-input branch
  git status # make sure we are in the right branch
  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/new_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/new_env


Just do it ::

  pynemo -s namelist.bdy # before doing this command, you need to have you namelist.bdy set up and the associated ncml file created
(see section 8b and 8c)

If James or Jeff do any modifications on the python codes inside the branch of interest on GitLab, you can update your nrct folder by 
doing the following command ::

  git pull 
  (then put your password)
  python setup.py build clean
  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/new_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/new_env

Once your namelist.bdy file is correctly set up, just do ::

  module load anaconda/2.1.0  # Want python2
  source activate new_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH

  pynemo -s namelist.bdy # It should work now

Important note : 
Depending on which tides you want (TPXO or FES), you need to make a change on the nemo_bdy_tide3.py file ::

  cd /work/valegu/nrct/Python/pynemo/tide
  vi nemo_bdy_tide3.py
  ...
  tide_src = 'FES' # Implemented alternative tidal source dataset.
  #tide_src = 'TPXO' # Implemented alternative tidal source dataset.
  print 'tide_src: ',tide_src

Other : You Can use the ``-g`` option if you want the GUI.

---

8.b) Generate the ncml file that points to the external data 
  
You can provide your grid information through ncml files for sn_src_hgr, sn_src_zgr and sn_src_msk. Because my data where not on a 
thredds server, I just put the name of the netcdf file for sn_src_hgr, sn_src_zgr and sn_src_msk. Those files comes from Gaby, located
here : /projectsa/NEMO/gmaya/ORCA12/
The only ncml file I used is inputs_dst.ncml ::
 
 cd $INPUTS
  vi inputs_dst.ncml
  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
    <ns0:aggregation type="union">
      <ns0:netcdf location="file:domain_cfg.nc">
      <ns0:variable name="mbathy" orgName="bottom_level" />
      <ns0:variable name="e3u" orgName="e3u_0" />
      <ns0:variable name="e3v" orgName="e3v_0" />
      </ns0:netcdf>
    </ns0:aggregation>
  </ns0:netcdf>

---

8.c) Generate the namelist.bdy file for Pynemo (NRCT)

Set up your namelist_cfg ::

  vi namelist.bdy
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  : namelist for BDY generation tool
  !!
  !!             User inputs for generating open boundary conditions
  !!             employed by the BDY module in NEMO. Boundary data
  !!             can be set up for v3.2 NEMO and above.
  !!
  !!             More info here.....
  !!
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !-----------------------------------------------------------------------
  !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)
     ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
     ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
     rn_hmin     =   -10     !  min depth of the ocean (>0) or
                           !  min number of ocean level (<0)

  !-----------------------------------------------------------------------
  !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
     rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
     rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
                             !  (= ocean depth) (>0) (m)
     ln_s_sigma  = .true.    !  hybrid s-sigma coordinates
     rn_hc       =  150.0    !  critical depth with s-sigma

  !-----------------------------------------------------------------------
  !  grid information
  !-----------------------------------------------------------------------
     sn_src_hgr = './mesh_hgr.nc'           !  parent /grid/ *** val : obtained from Gaby (Jasmine repository), ORCA0083-N06 ***
     sn_src_zgr = './mesh_zgr_renamed.nc'   !  parent        *** val : obtained from Gaby (Jasmine repository), ORCA0083-N06 and renamed the variables inside so that it has the same as the one in ORCA0083-N01 ***
     sn_dst_hgr = './domain_cfg.nc'
     sn_dst_zgr = './inputs_dst.ncml'       ! rename output variables
     sn_src_msk = './mask.nc'               ! parent         *** val : obtained from Gaby (Jasmine repository), ORCA0083-N06 ***
     sn_bathy   = './bathy_meter.nc'

  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
     sn_src_dir = 'local_inputs_feb.ncml'       ! src_files/'  *** val : specify here how to do the extraction (local files, local files automatized or via a thredds server) ***
     sn_dst_dir = '/work/valegu/NEMO/EA31FES_tides/INPUTS/'
     sn_fn      = 'EA31FES'                     !  prefix for output files
     nn_fv      = -1e20                         !  set fill value for output files
     nn_src_time_adj = 0                        !  src time adjustment
     sn_dst_metainfo = 'metadata info: valegu'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_coords_file = .true.               !  =T : produce bdy coordinates files
      cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = './bdy_mask.nc'      !  name of mask file (if ln_mask_file=.TRUE.)
      ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
      ln_tra         = .false.              !  boundary conditions for T and S
      ln_ice         = .false.              !  ice boundary condition

      nn_rimwidth    = 1                    !  width of the relaxation zone

  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
      ln_tide        = .true.               !  =T : produce bdy tidal conditions
      clname(1) =  'M2'
      clname(2) =  'S2'
      clname(3) =  'N2'
      clname(4) =  'K2'
      clname(5) =  'K1'
      clname(6) =  'O1'
      clname(7) =  'P1'
      clname(8) =  'Q1'
      clname(9) =  'M4'
      ln_trans       = .false.
      sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
      sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'

  !-----------------------------------------------------------------------
  !  Time information 
  !-----------------------------------------------------------------------
      nn_year_000     = 2013        !  year start
      nn_year_end     = 2013        !  year end
      nn_month_000    = 01          !  month start (default = 1 is years>1)
      nn_month_end    = 12          !  month end (default = 12 is years>1)
      sn_dst_calendar = 'gregorian' !  output calendar format
      nn_base_year    = 1900        !  base year for time counter  *** val : this value is obtained via ncdump -h of your parent dataset ***
      sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
      nn_src_time_adj    = 0  !-3168000 !- 86400 ! fix to align model time stamp *** val : James confirmed that I should put the value 0 here ***
  !-----------------------------------------------------------------------
  !  Additional parameters
  !-----------------------------------------------------------------------
      nn_wei  = 1                   !  smoothing filter weights
      rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                    !  smoothing onto dst points. Need to
                                    !  make this a funct. of dlon
      sn_history  = 'bdy files produced by valegu from ORCA0083-N06'
                                    !  history for netcdf file
      ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
      nn_alpha    = 0               !  Euler rotation angle
      nn_beta     = 0               !  Euler rotation angle
      nn_gamma    = 0               !  Euler rotation angle
      rn_mask_max_depth = 300.0     !  Maximum depth to be ignored for the mask
      rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break

---

9. Generation of the 3D boundary conditions (LIVLJOBS4)
  
9.a) Parent dataset
  
Jeff is using ORCA0083 version N01, whereas I am using ORCA0083 version N06. The two datasets are not entirely identical and 
have different variables name. The ORCA0083 parent file datasets can be found here : http://gws-access.ceda.ac.uk/public/nemo/runs/
Notice that ORCA0083-N01 range from Jan-1978 to Dec-2010 and ORCA0083-N06 range from Jan-1958 to Dec-2012.
I am doing my runs from Jan-2013 to Dec-2015. My parent files have been obtained from Gaby as she has access to the repository on 
JASMIN. She copied them into a folder on LIVLJOBS4, so I can access them. The path on LIVLOJOBS4 is : /projectsa/NEMO/gmaya/ORCA12/2013. 

If you already have an account on JASMIN, you can find the data here : 
/group_workspaces/jasmin2/nemo/vol1/ORCA0083-N006/means

I have adapted the method so that it target specifically the version N06 of ORCA0083. If you use N01 datasets, 
you don't need to do all the variables renaming. 

Because of the different names inside the files, it can cause an issue later on, when trying to generate the 3D boundary conditions. 
This is why I renamed them so that it is similar as what is observed inside N01 ::

  module load nco/gcc/4.4.2.ncwa
  ncrename -v e3t_0,e3t -v e3u_0,e3u -v e3v_0,e3_v -v e3w_0,e3w -v gdept_0,gdept -v gdept_1d,gdept_0 -v gdepw_1d,gdepw_0 -v e3t_1d,e3t_0 -v e3w_1d,e3w_0 mesh_zgr.nc mesh_zgr_renamed.nc

Rename as well T,S,U,V and SSH. But you don't have to do it manually for each file... This is done through the ncml file automatically (see 9.3)  

---

9.b) Setting up Pynemo for ORCA0083 N06

LIVLJOBS4

I created a new environment (named "new_env") that specifies the version of numpy and matlplotlib.

Set up pynemo with the ORCA0083 branch. OnGitHub, Jeff and James have a development branch called ORCA0083 that I have to use in 
order to create the 3D boundaries :: 

  ssh -X livljobs4
  cd /work/valegu/NEMO/EA31FES/INPUTS
  module load anaconda/2.1.0  # Want python2
  conda create --name new_env scipy=0.16.0 numpy=1.9.2 matplotlib=1.4.3 basemap netcdf4 libgfortran=1.0.0 # Another way that works (Jeff did it) is to not specify the version number : conda create --name new_env scipy=0.16.0 numpy matplotlib basemap netcdf4 libgfortran=1.0.0
  source activate new_env
 
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius
  
  cd $WORK/$USER
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd $WORK/$USER/nrct/Python
  git checkout ORCA0083 # Go inside the ORCA0083 branch
  git status # make sure we are in the right branch
  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/new_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/new_env

Once the branch is well set up, do only the following steps (If it doesn't work, try to open another session on Livljobs4 and try
ONLY the following again) ::

  module load anaconda/2.1.0  # Want python2
  source activate new_env
  cd $INPUTS # For me it's /work/valegu/NEMO/EA31FES_3Dboun/INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  pynemo -s namelist.bdy # Now it should work. 

---

9.c) Namelist.bdy set up

Notice that there is a section about the vertical coordinates inside the namelist.bdy and that ln_zps is put to true. However, I am
doing s-coordinates, and ln_sco is false. I tried to put ln_sco to true, but it didn't work. It could only work for ln_zps = true.
After talking with Jeff, he mentioned that it doesn't matter here.

You have to change the grid information section. 
Mesh_hgr.nc, mesh_zgr.nc and mask.nc are obtained from Gaby through the JASMIN repository (ORCAOO83-N06 in my case) : 
/projectsa/NEMO/gmaya/ORCA12/. As a remark, those files regarding the domain, should be identical to the one you can obtain
online : http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N06/domain/
 
How to target the parent data? 

There is three ways this can be done : 

1. Use local files (maybe not the most efficient way as it suppose you have to download the data first) and specify manually each
of the files you want to use. Below is an example for a month of boundary extraction. Notice that if the 5days mean dataset doesn't
start on the day 1, you have to put the previous dataset from the month before. Do the same for the end of the month.
(/work/valegu/NEMO/EA31FES_3Dboun/INPUTS/local_inputs_may.ncml) ::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130430d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130505d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130510d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130515d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130520d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130525d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130530d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130604d05T.nc" />
          <ns0:variable name="votemper" orgName="thetao" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130430d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130505d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130510d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130515d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130520d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130525d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130530d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130604d05T.nc" />
          <ns0:variable name="vosaline" orgName="so" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130430d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130505d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130510d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130515d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130520d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130525d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130530d05U.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130604d05U.nc" />
          <ns0:variable name="vozocrtx" orgName="uo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130430d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130505d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130510d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130515d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130520d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130525d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130530d05V.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130604d05V.nc" />
          <ns0:variable name="vomecrty" orgName="vo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130430d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130505d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130510d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130515d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130520d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130525d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130530d05T.nc" />
          <ns0:netcdf location="file:/projectsa/NEMO/gmaya/ORCA12/2013/ORCA0083-N06_20130604d05T.nc" />
          <ns0:variable name="sossheig" orgName="zos" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>


2. Use local files but do an automatized extraction (avoid having to write down all the files names you want)
(/work/valegu/NEMO/EA31FES_3Dboun/INPUTS/automatized_inputs.ncml) ::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/2013/" suffix="T.nc" />
          <ns0:variable name="votemper" orgName="thetao" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/2013/" suffix="T.nc" />
          <ns0:variable name="vosaline" orgName="so" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/2013/" suffix="U.nc" />
          <ns0:variable name="vozocrtx" orgName="uo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/2013/" suffix="V.nc" />
          <ns0:variable name="vomecrty" orgName="vo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/2013/" suffix="T.nc" />
          <ns0:variable name="sossheig" orgName="zos" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>


3. Use a thredds server (This is what Jeff is using. I can't use that as my data for 2013 were not available on a web server during the
period of my internship). A template of the syntax is as below ::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791106d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791101d05T.nc" />
          <ns0:netcdf location="http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N01/means/1979/ORCA0083-N01_19791027d05T.nc" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
    ...


My working namelist.bdy is located here : /work/valegu/NEMO/EA31FES_3Dboun/INPUTS/namelist.bdy ::

  vi namelist.bdy
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !! NEMO/OPA  : namelist for BDY generation tool
  !!
  !!             User inputs for generating open boundary conditions
  !!             employed by the BDY module in NEMO. Boundary data
  !!             can be set up for v3.2 NEMO and above.
  !!
  !!             More info here.....
  !!
  !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !-----------------------------------------------------------------------
  !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)
     ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
     ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)   *** val : couldn't work with ln_sco set to true ***
     rn_hmin     =   -10     !  min depth of the ocean (>0) or
                           !  min number of ocean level (<0)

  !-----------------------------------------------------------------------
  !   s-coordinate or hybrid z-s-coordinate
  !-----------------------------------------------------------------------
     rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
     rn_sbot_max = 7000.     !  maximum depth of s-bottom surface
                             !  (= ocean depth) (>0) (m)
     ln_s_sigma  = .true.    !  hybrid s-sigma coordinates
     rn_hc       =  150.0    !  critical depth with s-sigma

  !-----------------------------------------------------------------------
  !  grid information
  !-----------------------------------------------------------------------
     sn_src_hgr = './mesh_hgr.nc'           !  parent /grid/ *** val : obtained from Gaby (Jasmine repository), ORCA0083-N06 ***
     sn_src_zgr = './mesh_zgr_renamed.nc'   !  parent        *** val : obtained from Gaby (Jasmine repository), ORCA0083-N06 and renamed the variables inside so that it has the same as the one in ORCA0083-N01 ***
     sn_dst_hgr = './domain_cfg.nc'
     sn_dst_zgr = './inputs_dst.ncml'       ! rename output variables
     sn_src_msk = './mask.nc'               ! parent
     sn_bathy   = './bathy_meter.nc'

  !-----------------------------------------------------------------------
  !  I/O
  !-----------------------------------------------------------------------
     sn_src_dir = 'local_inputs_feb.ncml'       ! src_files/'  *** val : specify here how to do the extraction (local files, local files automatized or via a thredds server) ***
     sn_dst_dir = '/work/valegu/NEMO/EA31FES_3Dboun/INPUTS/'
     sn_fn      = 'EAfrica3D'                   !  prefix for output files
     nn_fv      = -1e20                         !  set fill value for output files
     nn_src_time_adj = 0                        !  src time adjustment
     sn_dst_metainfo = 'metadata info: valegu'

  !-----------------------------------------------------------------------
  !  unstructured open boundaries
  !-----------------------------------------------------------------------
      ln_coords_file = .true.               !  =T : produce bdy coordinates files
      cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
      ln_mask_file   = .false.              !  =T : read mask from file
      cn_mask_file   = './bdy_mask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
      ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
      ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
      ln_tra         = .true.               !  boundary conditions for T and S
      ln_ice         = .false.               !  ice boundary condition
      nn_rimwidth    = 9                    !  width of the relaxation zone

  !-----------------------------------------------------------------------
  !  unstructured open boundaries tidal parameters
  !-----------------------------------------------------------------------
      ln_tide        = .false.               !  =T : produce bdy tidal conditions
      clname(1) =  'M2'
      clname(2) =  'S2'
      clname(3) =  'N2'
      clname(4) =  'K2'
      clname(5) =  'K1'
      clname(6) =  'O1'
      clname(7) =  'P1'
      clname(8) =  'Q1'
      clname(9) =  'M4'
      ln_trans       = .false.
      sn_tide_h     = '/work/jelt/tpxo7.2/h_tpxo7.2.nc'
      sn_tide_u     = '/work/jelt/tpxo7.2/u_tpxo7.2.nc'

  !-----------------------------------------------------------------------
  !  Time information
  !-----------------------------------------------------------------------
      nn_year_000     = 2013        !  year start
      nn_year_end     = 2013        !  year end
      nn_month_000    = 02          !  month start (default = 1 is years>1)
      nn_month_end    = 02          !  month end (default = 12 is years>1)
      sn_dst_calendar = 'gregorian' !  output calendar format
      nn_base_year    = 1900        !  base year for time counter  *** val : this value is obtained via ncdump -h of your parent dataset ***
      sn_tide_grid    = '/work/jelt/tpxo7.2/grid_tpxo7.2.nc'
      nn_src_time_adj    = 0  !-3168000 !- 86400 ! fix to align model time stamp *** val : James confirmed that I should put the value 0 here ***
  !-----------------------------------------------------------------------
  !  Additional parameters
  !-----------------------------------------------------------------------
      nn_wei  = 1                   !  smoothing filter weights
      rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                    !  smoothing onto dst points. Need to
                                    !  make this a funct. of dlon
      sn_history  = 'bdy files produced by jelt from ORCA0083-N01'
                                    !  history for netcdf file
      ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
      nn_alpha    = 0               !  Euler rotation angle
      nn_beta     = 0               !  Euler rotation angle
      nn_gamma    = 0               !  Euler rotation angle
      rn_mask_max_depth = 300.0     !  Maximum depth to be ignored for the mask
      rn_mask_shelfbreak_dist = 60    !  Distance from the shelf break
                                                                                    
Once you have your namelist.bdy, then do ::

  module load anaconda/2.1.0  # Want python2
  source activate new_env
  cd /work/valegu/NEMO/EA31FES_3Dboun/INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH 
  pynemo -s namelist.bdy

Important remarks:

I had to extract month by month the 3D boundary conditions. Changing nn_month_000 = 01 and nn_month_end = 12 didn't work for me. 

Each month extraction took around 2h !! 

Once you have created the 3D boundary conditions, copy them into ARCHER : 
eg ::
  rsync -utv EAfrica__bdyT_y2013m05.nc valegu@login.archer.ac.uk:/work/n01/n01/valegu/EA31FES/INPUTS/

Important note : My January 3D dataset starts on the 3rd of January (size of 29). However, it looks like NEMO needs a file with 
the length of 31. 

I tried to extract again the 3D files by adding the datasets for the previous year 2012 inside the .ncml file 
(ORCA0083-N06_20121230d05T.nc, ORCA0083-N06_20121230d05U.nc and ORCA0083-N06_20121230d05V.nc). But it didn't work.

So what I did, is that I appended at the beginning of my 3D files, 2 times the date 3rd of January to obtain a file with the size of 31 ::

  module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
  module load cray-netcdf cray-hdf5
  module load nco/4.5.0

  ncks -d time_counter,0,0 EAfrica3D_bdyT_y2013m01.nc temp.nc # this only extract the variables along the time dimension and for only the first
  time period (3rd of January 2013 here). 
  ncrcat temp.nc temp.nc EAfrica3D_bdyT_y2013m01.nc Merged.nc

Now the file Merged.nc has the size of 31 (with the three first days being identical). I then renamed it ::

  cp Merged.nc EAfrica3D_bdyT_y2013m01.nc

Then I did the same for bdyU, bdyV and bt_bdyT.

---

9.d) James notes

I asked James some help to understand two issues : 

1 -> Why my initial conditions start on the 3rd of January (instead of the 1st). James will have a look at it. It is due to a problem
of variable names for the year 2012. I need to look at it later when time.

2 -> Why I couldn't compute in one go all the 12 months for the year 2013.

Below, a syntax provided by James (24/10/2018) for the automatic extraction that even take into account subdirectories! I didn't try 
it as I already had my 3D boundary files (obtained months by months...) ::

  <ns0:netcdf xmlns:ns0="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2" title="NEMO aggregation">
  <ns0:aggregation type="union">
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="temperature" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/" regExp=".*T\.nc$" subdirs="true"/>
          <ns0:variable name="votemper" orgName="thetao" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="salinity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/" regExp=".*T\.nc$" subdirs="true"/>
          <ns0:variable name="vosaline" orgName="so" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="zonal_velocity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/" regExp=".*U\.nc$" subdirs="true"/>
          <ns0:variable name="vozocrtx" orgName="uo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="meridian_velocity" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/" regExp=".*V\.nc$" subdirs="true"/>
          <ns0:variable name="vomecrty" orgName="vo" />
      </ns0:aggregation>
    </ns0:netcdf>
    <ns0:netcdf>
      <ns0:aggregation dimName="time_counter" name="sea_surface_height" type="joinExisting">
          <ns0:scan location="/projectsa/NEMO/gmaya/ORCA12/" regExp=".*T\.nc$" subdirs="true"/>
          <ns0:variable name="sossheig" orgName="zos" />
      </ns0:aggregation>
    </ns0:netcdf>
  </ns0:aggregation>
  </ns0:netcdf>

If it was taking 2h for ONE month of 3D boundaries, it was because Pynemo was trying to read the entire global dataset.  
James managed to fix it. I tried it and now a month of 3D boundary takes 86sec instead of the 2h as previously.
