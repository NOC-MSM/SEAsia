Generate open boundary conditions using PyNEMO
==============================================

* parent data on JASMIN.
* Generate multiple years in yearly blocks.

* Parent data: ORCA0083-N06
* Data web viewable at: http://gws-access.ceda.ac.uk/public/nemo/runs/ORCA0083-N06/means/1960/
* Data files readable at e.g.: /gws/nopw/j04/nemo_vol1/ORCA0083-N006/means/1960

For the purposes of progress I am going to do the open bcs extraction using
pynemo install on JASMIN ``<install_nrct.rst>_``. This might be better done
mounting the JASMIN fileserver over SSHFS https://github.com/jdha/PyNEMO/wiki/SSHFS

The following started as an implementation of James' notes: https://github.com/jdha/PyNEMO/wiki/Accessing-data:-Hints-and-Tips
But then he produced NCML files that took away all the pain. So jump to the next
section

Building NCML files
===================

On the target machine, where the files are to be generated, create a directory
to store the symbolically linked inputs and outputs. The inputs should be stored
in yearly directories for easier processing, with capping files from the last and first
outputs from the surrounding years to enable complete interpolation for the year.
::

  #!/bin/bash
  export PATH_TO_LOCAL_DATA=/gws/nopw/j04/campus/pseudoDropBox/SEAsia/ORCA0083-N006/
  export PATH_TO_REMOTE_DATA=/gws/nopw/j04/nemo_vol1/ORCA0083-N006/means/

  export YEAR_START=1960
  export YEAR_END=1970

  mkdir $PATH_TO_LOCAL_DATA

  # loop over years and symbolically link in all the files from the remote source.
  for y in {YEAR_START..YEAR_END}
    do
    mkdir $PATH_TO_LOCAL_DATA/$y
    cd $PATH_TO_REMOTE_DATA/$y
    for fn in *d05T.nc
      do
      ln -s $PATH_TO_REMOTE_DATA/$y/$fn $PATH_TO_LOCAL_DATA/$y/$fn
    done
  done


Then for each year add the end cap files from the neighbouring years to allow PyNEMO to
interpolate to the start and end of the year.
::

  #!/bin/bash

  for y in {YEAR_START..YEAR_END}
    do
    ym1=`expr $y - 1`
    yp1=`expr $y + 1`
    # Copy in Jan file from next year
    cd $PATH_TO_LOCAL_DATA
    ln -s $yp1/ORCA0083-N006_$yp1\0105d05T.nc $y/ORCA0083-N006_$yp1\0105d05T.nc

    # Copy in Dec file from prev year, avoiding issues with explicitly declaring
    #  Dec 31 or Dec 30 file due to the leap year
    cd $PATH_TO_LOCAL_DATA/$ym1
    for fn in ORCA0083-N006_$ym1\123?d05T.nc
      do
      ln -s $PATH_TO_LOCAL_DATA/$ym1/$fn $PATH_TO_LOCAL_DATA/$y/$fn
    done
  done

*INCOMPLETE because James did it*


Using James' NCML files
=======================

Now James has already set up NCML files for 1970 - 2015
/gws/nopw/j04/nemo_vol5/jdha/ORCA0083-N006/NCML
E.g. ORCA0083_N06_2001.ncml


Up to end of 2012 time counter seconds from 1950
From 2013 time counter from 1900. Need to accomodate this.



Generate pynemo namelist files in INPUTS directory where the necessary pynemo
input files sit (you have to put them there)
::

  livljobs6 ~ $
  exec ssh-agent $SHELL
  ssh-add ~/.ssh/id_rsa_jasmin
   Enter passphrase for /login/jelt/.ssh/id_rsa_jasmin:
  Identity added: /login/jelt/.ssh/id_rsa_jasmin (/login/jelt/.ssh/id_rsa_jasmin)
  ssh -A jelt@jasmin-login1.ceda.ac.uk
  ssh -A jelt@jasmin-sci1.ceda.ac.uk

  cd /gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS

Then run the scipt to generate namelist files::

  #!/bin/bash

  export PATH_PYNEMO_OUTPUTS=/gws/nopw/j04/campus/pseudoDropBox/SEAsia/INPUTS/
  export PATH_TO_NCML=/gws/nopw/j04/nemo_vol5/jdha/ORCA0083-N006/NCML/
  export YEAR_START=1960
  export YEAR_END=2015

  cd $PATH_PYNEMO_OUTPUTS

  for y in $(seq $YEAR_START $YEAR_END)
    do
    # Write the namelist file for the whole year
    echo year: $y
    sed "s/YEAR_START/$y/g" namelist.bdy_TEMPLATE > tmp1.bdy
    sed "s/YEAR_END/$y/g"   tmp1.bdy > tmp2.bdy
    sed "s?YYYY_NCML_FILE?$PATH_TO_NCML\ORCA0083_N06_$y.ncml?" tmp2.bdy > tmp3.bdy
    sed "s?DST_DIR?$PATH_PYNEMO_OUTPUTS?" tmp3.bdy > namelist_$y.bdy

    rm tmp?.bdy

    # Submit PyNEMO
    #pynemo -s namelist_$y.bdy

  done
