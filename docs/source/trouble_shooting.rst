.. contents:: Table of Contents

*****
Trouble Shooting
*****

Things to try if your new configuration isn't working

Mysterious fail. Either no clue of XIOS / XML / netcdf hint
===========================================================

The XML control of I/O, and in particular formatting and content of the ``iodef.xml`` files (and its kin) are **extremely** sensitive to errors.

Several times an inconsistency in the xml files, or a typo, can lead to a NEMO failure with little or no debugging info.

E.g. top line of ``iodef.xml``. Switching from ::

    <file_definition type="one_file" name="@expname@_@freq@_@startdate@_@enddate@" sync_freq="10d" min_digits="4">

to ::

    <file_definition type="multiple_file" name="@expname@_@freq@_@startdate@_@enddate@" sync_freq="10d" min_digits="4">

Can mean the difference between a nice weekend or a bad one. It seems that if you decide to use multi XIOS cores, you need 
to use the second option (**multiple_file**).


If a run suddenly stops without any errors in your NEMO **ocean.output** or in the cluster error/output files, it could 
come from the memory of the XIOS server (you probably get in the error file something like "OOM killer terminated this process."). 
A way to pass over this is to use more nodes, and using a single core per node to access the full memory of the node. 
An example for Archer, using 1 core of 12 nodes: ::

   aprun -b -n 12 -N 1 ./xios_server.exe : -n $OCEANCORES -N 24 ./nemo.exe


My configuration blows up
=========================

This can happen for a host of reasons...

#. Probaby the first thing to check is whether it is behaving without all the
add-on forcings. Is it stable if you turn everything off? Is it stable with
 clamped initial condition temperature and salinity, and boundary velocities
It is a good idea to try and 20 day 'initial condition' run to make sure there
 are no surprises.

Compile with ``usrdef_istate.F90`` and  ``usrdef_sbc.F90`` in ``MY_SRC``. Then
 edit the ``namelist_cfg`` to turn off forcings (checking there are no inconsistencies)::

  &namsbc        !   Surface Boundary Condition (surface module)
  !-----------------------------------------------------------------------
                     ! Type of air-sea fluxes
     ln_usr      = .true.    !  user defined formulation                  (T => check usrdef_sbc)

  ...

  !-----------------------------------------------------------------------
  &nam_tide      !   tide parameters
  !-----------------------------------------------------------------------
    ln_tide     = .true.
    ln_tide_pot = .false.    !  use tidal potential forcing

  ...

  &nambdy        !  unstructured open boundaries
  !-----------------------------------------------------------------------
     ln_bdy         = .true.              !  Use unstructured open boundaries
     nn_dyn2d_dta   =  0                   !  = 0, bdy data are equal to the initial state
     nn_dyn3d_dta  =  0                    !  = 0, bdy data are equal to the initial state
     nn_tra_dta    =  0                    !  = 0, bdy data are equal to the initial state


#. Tidal models can blow up at, or near, the boundaries. This can happen for a number
of reasons: Does the domain have amphidromes near the boundary? This can result
large spatial gradients in velocity that can cause problems.


#. If the configuration blows up at, or near, the boundary modify the boundary mask to blank it out.
This means that new boundary files will need to be generated for the new effective boundary location

**Update**. Modification of the mask is not such good advice. If there is a problem
with an island or some significant bathymetry feature appearing on the child grid
 but not in the parent, then this maybe a good idea.
Note the PyNEMO mask takes three value. One each for mask(-1), wet points (1) and dry points (0).
E.g. Start with with the land mask from ``domain_cfg.nc`` and introduce boundary masking. First
create a mask file from a template. (Using **livljobs4**)::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

In ipython manually edit the mask locations::

  import netCDF4
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  dset.variables['mask'][0:4,:]  = -1
  dset.variables['mask'][-1,:] = -1
  dset.variables['mask'][:,-1] = -1
  dset.variables['mask'][:,0] = -1
  dset.close()

Then ``bdy_mask.nc`` can be specified in the PyNEMO ``namelist.bdy``. The PyNEMO
 generated files contain the bdy_msk variable, for use in the NEMO ``namelist_cfg``

Run PyNEMO again. Run NEMO again.

---

If the bdy_msk does not appear to be functional. Perhaps missing updates to the
OPA source::

  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/bdyini.F90 $CDIR/$CONFIG/MY_SRC/.
  cp /work/n01/n01/jdha/2017/nemo/trunk/NEMOGCM/CONFIG/ORCHESTRA/MY_SRC/dommsk.F90 $CDIR/$CONFIG/MY_SRC/dommsk.F90

---

#. If the model is blowing up at the boundary and the water is deep. Check the time step. Deepwater waves are fast.

---

#. If the model is blowing up at the boundary and the water is shallow. Have the tidal transports be mapped from parent to child grid correctly?

---


Memory error for combining outputs (ARCHER)
===========================================

If your configuration becomes massive, combining the output might bring memory issues on **ARCHER** login nodes.
A solution is to submit on the post-processing node. However post-processing nodes and computing nodes have different
architecture and you need to recompile your tools for it. Basically on those node the compiler shortcuts (ftn, CC, ...) 
are not recognized so you need to alter them depending on the compiler. for example with intel, **ftn** becomes **ifort**.


More details can be find on the ARCHER documentation :
   http://www.archer.ac.uk/documentation/user-guide/development.php#sec-4.7


Pynemo on livl server
=====================

I (Nico) could not manage to install pynemo locally on my work computer (worked fine at home). I change my anaconda set-up to install the 
environments in the **/work** instead of **/login** through the **.condarc**. Finally the only way I manage to install it was to reverse 
my environment to the **/login** default one. It seems weird to not work in the other way and it's not very class was of sorting this but 
at least it works. I guess it's only due to *java virtual machine* path not properly path trhough but still... the error I got was the 
following ::
  
   File "/work/nibrun/nico-conda/nrct_env/lib/python2.7/site-packages/jnius/reflect.py", line 162, in autoclass
       c = find_javaclass(clsname)
   File "jnius_export_func.pxi", line 23, in jnius.find_javaclass (jnius/jnius.c:12356)
   jnius.JavaException: Class not found 'ucar/nc2/dataset/NetcdfDataset'











