Generate open boundary conditions using PyNEMO
==============================================

* parent data on JASMIN.
* Generate multiple years in yearly blocks.

* Parent data: ORCA0083-N06
* Data web viewable at: http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N06/means/1960/
* Data files readable at: /gws/nopw/j04/nemo_vol1/ORCA0083-N006/means/1960

For the purposes of progress I am going to do the open bcs extraction using
pynemo install on JASMIN ``<install_nrct.rst>_``. This might be better done
mounting the JASMIN fileserver over SSHFS https://github.com/jdha/PyNEMO/wiki/SSHFS

This is an implementation of James' notes: https://github.com/jdha/PyNEMO/wiki/Accessing-data:-Hints-and-Tips
