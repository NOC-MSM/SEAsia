.. contents:: Table of Contents

*****
ARCHER: compilation and useful tools
*****



Compile NEMO on ARCHER
======================

You need to compile both XIOS and NEMO (with the same compilers). For NEMO 3.6, you need XIOS 1.0 (rev 703 max). With NEMO 4.0, 
you need XIOS 2.0 (trunk). NEMO provides architecture files already for ARCHER.

Cray compiler
-------------

Load the set of modules ::

    module load cray-libsci/16.11.1
    module load cray-mpich/7.5.5
    module load cray-netcdf-hdf5parallel/4.4.1.1
    module load cray-hdf5-parallel/1.10.0.1
    module load PrgEnv-cray/5.2.82

And you can use of the ARCHER NEMO template for compiling (**ARCH/arch-XC_ARCHER.fcm**); only need to update the path to xios ::

  %NCDF_HOME           $NETCDF_DIR
  %HDF5_HOME           $HDF5_DIR
  %XIOS_HOME           /work/n01/n01/nibrun/NEMO/xios-1.0_r703/

  %NCDF_INC            -I%NCDF_HOME/include -I%HDF5_HOME/include
  %NCDF_LIB            -L%HDF5_HOME/lib -L%NCDF_HOME/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %CPP                 cpp
  %FC                  ftn
  %FCFLAGS             -em -s integer32 -s real64 -O0 -e0 -eZ
  %FFLAGS              -em -s integer32 -s real64 -O0 -e0 -eZ

For XIOS, you can use **arch-XC30_Cray** file to compile.

Intel compiler
-------------

Similarly ::

  module load cray-libsci/16.11.1
  module load cray-mpich/7.5.5
  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel

and for the compilation ::

  %NCDF_HOME           $NETCDF_DIR
  %HDF5_HOME           $HDF5_DIR
  %XIOS_HOME           /work/n01/n01/jelt/xios-1.0_r703/

  %NCDF_INC            -I%NCDF_HOME/include -I%HDF5_HOME/include
  %NCDF_LIB            -L%HDF5_HOME/lib -L%NCDF_HOME/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %CPP                 cpp
  %FC                  ftn
  %FCFLAGS             -integer-size 32 -real-size 64 -O3 -fp-model source -zero -fpp -warn all
  %FFLAGS              -integer-size 32 -real-size 64 -O3 -fp-model source -zero -fpp -warn all
  %LD                  CC -Wl,"--allow-multiple-definition"
  %FPPFLAGS            -P -C -traditional
  %LDFLAGS
  %AR                  ar
  %ARFLAGS             -r
  %MK                  gmake
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB

  %CC                  cc
  %CFLAGS              -O0

General comments
----------------

Results cannot be seen on the fly during the run. From a quick set of tests, killing the job is likely to lead to no output file access. 
If the job stops before the end, in my case, I couldn't access the data with cray compiled executable but could with intel...

Compiling XIOS with Cray was slow...

Running ncview
=========================

You need to slightly update / change your environement. If parrallel NetCDF and h5 loaded, they need to be switched to the non-parallel version.
In addition, python through anaconda should not be loaded as it breaks links required by ncview ::

   module swap cray-netcdf-hdf5parallel/4.4.1.1 cray-netcdf
   module swap cray-hdf5-parallel/1.10.0.1 cray-hdf5
   module unload anaconda
   module load ncview



