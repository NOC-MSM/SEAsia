Rebuild the output and inspect
++++++++++++++++++++++++++++++

Rebuild the NEMO output. If you have followed the path definition thing::

  . ~/temporary_path_names_for_NEMO_build

Otherwise do something like::

  export WDIR=/work/n01/n01/jelt/LBay/
  export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

---

Rebuild the **OUTPUT**::

  export filehead=Lbay_1h_20000101_20000130
  export nproc=5

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 $filehead_SSH $nproc
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 $filehead_Tides $nproc

Should remove individual processor files once the build is verified::

  rm $filehead_SSH_00??.nc
  rm $filehead_Tides_00??.nc

---

Rebuild the **ABORT** files::

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 output.abort 93

#Should remove individual processor files once the build is verified::

  rm output.abort_00??.nc

---

Rebuild the **RESTART** files (note 3 processors were recovered / not used because of land)::

  $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 SWPacific_00004800_restart_tide 93

#Should remove individual processor files once the build is verified::

  rm SWPacific_00004800_restart_tide_00??.nc

---

Inspect on ARCHER::

  cd $EXP
  module load anaconda
  python quickplotNEMO.py


---

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
