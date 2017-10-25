Rebuild the output and inspect
++++++++++++++++++++++++++++++

Rebuild the SSH files using old tools::

  export WDIR=/work/n01/n01/jelt/LBay/
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 Lbay_1h_20000101_20000130_SSH 5
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 Lbay_1d_20000101_20000130_Tides 5

Should remove individual processor files once the build is verified::

  rm Lbay_1h_20000101_20000130_SSH_00??.nc
  rm Lbay_1d_20000101_20000130_Tides_00??.nc

Inspect locally e.g.::

  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/Lbay_1d_20000101_20000130_Tides.nc .
  scp jelt@login.archer.ac.uk:/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/CONFIG/LBay/EXP00/Lbay_1h_20000101_20000130_SSH.nc .

  ferret
  use Lbay_1d_20000101_20000110_Tides.nc
  plot /i=25/j=70 SOSSHEIG

Use some python to inspection of the domain_cfg.nc file or ssh, Tide output. See::

  cd /Users/jeff/GitLab/NEMO-RELOC/docs/source
  $ipython
  >> run quickplotNEMO.py

.. note : it may be better to use ncview to inspect i.e.:

    module load ncview
    ncview file.nc

---
