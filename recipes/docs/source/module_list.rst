.. _module_list:

**************************
Versions of ARCHER modules
**************************

Sometimes archer updates modules and then things break...
Here is a list of modules that have worked in the past.

A number of these recipes call for::

  module swap PrgEnv-cray PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  module load cray-hdf5-parallel

As of 22 Nov 2018 my currently loaded modules versions are::

  Currently Loaded Modulefiles:
    1) modules/3.2.10.6                      13) leave_time/1.3.0                      25) dmapp/7.0.1-1.0502.11080.8.76.ari
    2) eswrap/1.3.3-1.020200.1280.0          14) quickstart/1.0                        26) gni-headers/4.0-1.0502.10859.7.8.ari
    3) switch/1.0-1.0502.60522.1.61.ari      15) ack/2.14                              27) xpmem/0.1-2.0502.64982.5.3.ari
    4) intel/17.0.0.098                      16) xalt/0.6.0                            28) dvs/2.5_0.9.0-1.0502.2188.1.116.ari
    5) craype-network-aries                  17) openssl/1.1.0g_build1                 29) alps/5.2.4-2.0502.9774.31.11.ari
    6) craype-ivybridge                      18) curl/7.58.0_build1                    30) rca/1.0.0-2.0502.60530.1.62.ari
    7) craype/2.5.10                         19) git/2.16.2_build1                     31) atp/2.1.0
    8) pbs/12.2.401.141761                   20) epcc-tools/8.0                        32) PrgEnv-intel/5.2.82
    9) cray-mpich/7.5.5                      21) cray-libsci/16.11.1                   33) cray-hdf5-parallel/1.10.0.1
   10) packages-archer                       22) udreg/2.3.2-1.0502.10518.2.17.ari     34) cray-netcdf-hdf5parallel/4.4.1.1
   11) bolt/0.6                              23) ugni/6.0-1.0502.10863.8.29.ari        35) anaconda/python2
   12) nano/2.2.6                            24) pmi/5.0.13


Note (32), (33) ,(34)::
   PrgEnv-intel/5.2.82
   cray-hdf5-parallel/1.10.0.1
   cray-netcdf-hdf5parallel/4.4.1.1
