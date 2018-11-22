Build XIOS2 @ r1242
+++++++++++++++++++

Note when NEMO (nemo.exe / opa) is compiled it is done with reference to a particular version of
XIOS. So on NEMO run time the version of XIOS that built xios_server.exe must be compatible with the
version of XIOS that built nemo.exe / opa.

Modules::

  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

Download XIOS2 and prep::

  cd $WORK/$USER
  svn co -r1242 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1242
  cd xios-2.0_r1242
  cp $WORK/$USER/ARCH/arch-XC30_ARCHER* arch/.

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 8

Link the xios-2.0_r1242 to a generic XIOS directory name::

  ln -s  $WORK/$USER/xios-2.0_r1242  $WORK/$USER/XIOS

Link xios executable to the EXP directory::

  ln -s  $WORK/$USER/xios-2.0_r1242/bin/xios_server.exe $EXP/xios_server.exe
