Build XIOS2 @ r1080
+++++++++++++++++++

::

  cd $WORK/$USER
  svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1080
  cd xios-2.0_r1080
  cp $WORK/$USER/ARCH/arch-XC30_ARCHER* arch/.

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par

Link the xios-2.0_r1080 to a generic XIOS directory name::

  ln -s  $WORK/$USER/xios-2.0_r1080  $WORK/$USER/XIOS

Link xios executable to the EXP directory::

  ln -s  $WORK/$USER/xios-2.0_r1080/bin/xios_server.exe $EXP/xios_server.exe
