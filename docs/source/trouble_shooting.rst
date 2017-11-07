Trouble Shooting
++++++++++++++++
++++++++++++++++

Things to try if your new configuration isn't working

Mysterious fail. Either no clue of XIOS / XML / netcdf hint
===========================================================

The XML control of I/O, and in particular formatting and content of the
 ``iodef.xml`` files (and its kin) are **extremely** sensitive to errors.
  Several times an inconsistency in the xml files, or a typo, can lead to a NEMO failure with little or no debuggin info


E.g. top line of iodef.xml. Switching from ::

    <file_definition type="one_file" name="@expname@_@freq@_@startdate@_@enddate@" sync_freq="10d" min_digits="4">

to::
      <file_definition type="multiple_file" name="@expname@_@freq@_@startdate@_@enddate@" sync_freq="10d" min_digits="4">

Can mean the difference between a nice weekend or a bad one.


My configuration blows up
+++++++++++++++++++++++++

If the configuration blows up at, or near, the boundary modify the boundary mask to blank it out.
This means that new boundary files will need to be generated for the new effective boundary location

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
  dset.variables['mask'][:,-4:-1] = -1
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

If the model is blowing up at the boundary and the water is deep. Check the time step. Deepwater waves are fast.

---

If the model is blowing up at the boundary and the water is shallow. Have the tidal transports be mapped from parent to child grid correctly?
