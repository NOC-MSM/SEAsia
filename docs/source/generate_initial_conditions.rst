Generate Initial conditions
+++++++++++++++++++++++++++
+++++++++++++++++++++++++++

For a new configuration you probably want to start with idealised, or homogenous
initial conditions. This is done with user defined initial conditions ``ln_usr``
with the expression being compiled into the executable

To use initial conditions from an existing T,S field you might need to do a bit
of interpolation. It is advisable to let NEMO do the heavy lifting for vertical
interpolation (rquiring some FORTRAN modifictions), though SOSIE tools can be user
to do simple horizontal interpolation.


User defined initial initial conditions
=======================================

For constant T and S use the user defined functions in ``$CDIR/$CONFIG/MY_SRC``:
  ``usrdef_sbc.F90``  and ``usrdef_istate.F90``. Compile and save executable with
  telegraphic names that point to compile options. e.g.::

    nemo_notide_TSprofile.exe
    nemo_tideonly_TSconst.exe


Building T,S field initial conditions from existing fields
==========================================================

Second time around we build 3D initial conditions
*(27 Apr 2018)*

*Since my parent and child are on the same grid I'm not sure I need the SOSIE
step. Try skipping it for now.*


Use SOSIE tools to interpolate onto a finer grid
------------------------------------------------

Copy ``make.macro`` file and edit the path if necessary::
**FIX** to the notes (copied from jdha instead): ``cp $WDIR/INPUTS/make.macro ./``::

  cp /home/n01/n01/jdha/sosie/make.macro /home/n01/n01/jelt/sosie/.

  vi /home/n01/n01/jelt/sosie/make.macro
  # Directory to install binaries:
  INSTALL_DIR = /home/n01/n01/jelt/local

Proceed with Step 6 (of Lighhouse Reef Readthedocs). This is best done in a clean terminal
::

  cd ~
  mkdir local
  svn co svn://svn.code.sf.net/p/sosie/code/trunk sosie
  cd sosie

  make
  make install
  export PATH=~/local/bin:$PATH
  cd $WDIR/INPUTS


Obtain the fields to interpolate. E.g interpolate AMM60 or ORCA
data. Get the namelists::

  cp $START_FILES/initcd_votemper.namelist $INPUTS/.
  cp $START_FILES/initcd_vosaline.namelist $INPUTS/.

The sosie routine is VERY slow. (2.5 hrs)
*UPDATE 30 Apr 18: possible because I was doing unnecessary 3D interpolation*.
So, make a cut down parent file using ORCA0083-N01.
Cut down based on coordintaes from create coordinates namelist. (Add a buffer as
I'm not sure how the sosie extraction works)::

    module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
    module load cray-netcdf cray-hdf5
    module load nco/4.5.0
    cd $WDIR/INPUTS

    ncks -d x,45,735 -d y,1245,1810 ORCA0083-N01_19791101d05T.nc $WDIR/INPUTS/cut_down_19791101d05_SEAsia_grid_T.nc

Average over time and restore the parallel modules (Not necessary for this data with 1 time point)::

    #ncwa -a time_counter $WDIR/INPUTS/cut_down_20131013_LBay_grid_T.nc  $WDIR/INPUTS/cut_down_201310_LBay_grid_T.nc

    module unload nco cray-netcdf cray-hdf5
    module load cray-netcdf-hdf5parallel cray-hdf5-parallel


----



Edit namelists to the variables you want.
It is advisable to let NEMO do the vertical interpolation so only use SOSIE tools for 2D
interpolation.

In this exmaple I am now pointlessly interpolating something onto itself. This is a template and would
be useful when I have a child grid at a finer resolution than the parent.

Set jplev = 1
To interpolate only on the first level (and output level..). **NOT SURE HOW YOU FIX IT TO DO 2D INTERPOLATION**
::

  vi initcd_votemper.namelist
  cf_in     = 'cut_down_19791101d05_SEAsia_grid_T.nc'
  cv_in     = 'votemper'
  cf_x_in   = 'cut_down_19791101d05_SEAsia_grid_T.nc'
  cv_out   = 'votemper'
  csource  = 'ORCA0083-N01'
  ctarget  = 'SEAsia'

  vi initcd_vosaline.namelist
  ...
  cv_out   = 'vosaline'
  ...

Copy parent file to ARCHER INPUTS (need to generalise / improve)::

  livljobs4
  scp /projectsa/accord/ORCA0083/ORCA0083-N01_19791101d05T.nc jelt@login.archer.ac.uk:/work/n01/n01/jelt/SEAsia/INPUTS/.


Do stuff (on ARCHER). I think the intention was for SOSIE to flood fill the land::

  cd $INPUTS
  sosie.x -f initcd_votemper.namelist
  sosie.x -f initcd_vosaline.namelist


I had trouble getting ARCHER to run this. (Though with the cut down parent the
 commandline is fine, though it runs out of walltime in Interactive Mode)
Initially, running in the commandline the job failed with insufficient memory,
 because I didn't cut the ORCA data down first.
 In the end I submitted it as a pair of serial jobs. **IT TOOK 4hrs 25m TO DO 3D**::

  vi $INPUTS/sosie_initcd_T

  #!/bin/bash
  #PBS -N init_T
  #PBS -l select=serial=true:ncpus=1
  #PBS -l walltime=06:00:00
  #PBS -o init_T.log
  #PBS -e init_T.err
  #PBS -A n01-ACCORD
  ###################################################

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel


  cd /home/n01/n01/jelt/sosie
  make clean
  make
  make install

  #set up paths
  cd /work/n01/n01/jelt/SEAsia/INPUTS

  /home/n01/n01/jelt/local/bin/sosie.x -f initcd_votemper.namelist
  #/home/n01/n01/jelt/local/bin/sosie.x -f initcd_vosaline.namelist


  # qsub -q serial <filename>
  ###################################################


Similarly for ``sosie_initcd_S``. Then::

  qsub -q serial sosie_initcd_T
  qsub -q serial sosie_initcd_S

3 hours not enough - resubmit with 6 hrs! *It took 4h 25min*

Whether as a serial job or from the commandline, the temperature process creates::

  sosie_mapping_ORCA0083-N01-SEAsia.nc
  votemper_ORCA0083-N01-SEAsia_1978.nc4

And the salinity process creates::

  vosaline_ORCA0083-N01-SEAsia_1978.nc4


Use SCRIP tools to remap to the new grid
----------------------------------------

Now do interpolation onto child grid.  The ``scrip`` tools are build in ``TDIR``
e.g. in `Build Tools<SEAsia_archer_livljobs4.rst>`_
::

  export OLD_TDIR=$WORK/$USER/LBay/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/

First copy the namelists::

  cp $START_FILES/namelist_reshape_bilin_initcd_votemper $INPUTS/.
  cp $START_FILES/namelist_reshape_bilin_initcd_vosaline $INPUTS/.

Edit the input files::

  vi $INPUTS/namelist_reshape_bilin_initcd_votemper
  &grid_inputs
    input_file = 'votemper_ORCA0083-N01-SEAsia_1978.nc4'
  ...
    input_name = "votemper"

  &interp_inputs
    input_file = "votemper_ORCA0083-N01-SEAsia_1978.nc4"
  ...

Similarly for the *vosaline.nc file::

  vi $INPUTS/namelist_reshape_bilin_initcd_vosaline
  &grid_inputs
    input_file = 'vosaline_ORCA0083-N01-SEAsia_1978.nc4'
    ...
    input_name = "vosaline"
  ...

  &interp_inputs
    input_file = "vosaline_ORCA0083-N01-SEAsia_1978.nc4"
  ...


.. Note: Alternatively since child and parent are on the same horizontal grid,
  skip the SOSIE step and do the SCRIP interpolation only. To do this will need
  ammend the input_file and input_vars parameters. E.g.

  &interp_inputs
      input_file = "cut_down_19791101d05_SEAsia_grid_T.nc4"
      interp_file = "data_nemo_bilin_R12.nc"
      input_name = "vosaline"
      input_start = 1,1,1,1
      input_stride = 1,1,1,1
      input_stop = 0,0,0,0
      input_vars = "deptht", "time_counter"



Produce the remap files::

  $OLD_TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper

Creates ``remap_nemo_grid_R12.nc`` and ``remap_data_grid_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper

Creates ``data_nemo_bilin_R12.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper

Creates ``initcd_votemper.nc``. Then::

  $OLD_TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Creates ``initcd_vosaline.nc``.

---

Interpolate in z on the fly
===========================


For vertical interpolation we let NEMO do the heavy lifting. This requires some changes
to the FORTRAN.

Maria's code worked for v3.6: ``cd /work/n01/n01/mane1/ARC25v3.6/OPA_SRC``

Make some changes ``MY_SRC/dtatsd.F90`` including


line 25::

    USE iom

dta_tsd_init
line 46::

  #if defined key_gen_IC
     REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   gdept_init, gdept_dta, sal_dta,temp_dta
     REAL(wp), ALLOCATABLE, DIMENSION(:)     ::   gdept_init_1d
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   ssh_dta
     INTEGER                                 ::   jpk_init , inum_dta
     LOGICAL                                 ::   ln_tsd3  !( T if depth is 3d, else 1d)
     INTEGER ::   id ,linum   ! local integers
     INTEGER                                 ::   ddims(4),dimsd(3)
  #endif

line 107 change::

        ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
  IF( sn_tem%ln_tint )   ALLOCATE( sf_tsd(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
        ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
  IF( sn_sal%ln_tint )   ALLOCATE( sf_tsd(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )

Into::

  #if defined key_gen_IC
         CALL iom_open ( 'bathy_meter', inum_dta )
         !! get dimensions
         id = iom_varid( inum_dta, 'gdept_glo', dimsd )
         jpk_init = dimsd(3)
         IF(lwp) WRITE(numout,*) 'Dimensions of ICs: ', dimsd, jpk_init
                                ALLOCATE( temp_dta(jpidta,jpjdta,jpk_init)                , STAT=ierr0 )
                                ALLOCATE( sal_dta(jpidta,jpjdta,jpk_init)                 , STAT=ierr1 )
                                ALLOCATE( ssh_dta(jpidta,jpjdta         )                 , STAT=ierr2 )
                                ALLOCATE( gdept_dta (jpidta,jpjdta,jpk_init),               STAT=ierr3 )
       !
                                ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk_init)   , STAT=ierr4 )
                                ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk_init)   , STAT=ierr5 )
                                ALLOCATE( gdept_init         (jpi,jpj,jpk_init),    STAT=ierr6 )

         CALL iom_close( inum_dta )   ! Close the input file
  #else
                                ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
        IF( sn_tem%ln_tint )   ALLOCATE( sf_tsd(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
        IF( sn_sal%ln_tint )   ALLOCATE( sf_tsd(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
  #endif


line 120:
Desparate times, hardcoded bathymetry variable::

  vi  dtatsd.F90
  line 120:
    id = iom_varid( inum_dta, 'Bathymetry', dimsd )

This is a hack and probably ought to be some whizzy function from ``domain_cfg.nc``

Hardwired initial condition file::

  line 190:
  CALL iom_open ('bdydta/cut_down_19791101d05_SEAsia_grid_T.nc', sf_tsd(jp_tem)%num )






Compile with ``key_gen_IC``

Edit cpp_SEAsia.fcm::

  bld::tool::fppkeys key_zdfgls        \
                key_harm_ana      \
                key_gen_IC        \
                key_mpp_mpi       \
                key_iomput        \
                key_nosignedzero




Compile and debug on short queue::


  cd $CDIR

  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10

  mv SEAsia/BLD/bin/nemo.exe SEAsia/BLD/bin/nemo_tide_genIC_nomet.exe

  cd SEAsia/EXP_tide_initcd

  rm opa
  ln -s ../BLD/bin/nemo_tide_genIC_nomet.exe opa

  qsub runscript_short



**PENDING**
(3 May 2018)*
DID IT CRASH?
cd /work/n01/n01/jelt/SEAsia/trunk_NEMOGCM_r8395/CONFIG/SEAsia/EXP_tide_initcd

Crash report::


  less ocean.output

  ...


-----

scratch
=======

domain_cfg.nc
x = 684 ; 76.9166: 133.833
y = 554 ;  -20.0761: 24.7641
z = 75 ;
