============================
Caribbean setup instructions
============================

Ash Brereton, June 2018
Machines: ARCHER
Based on: `<SEAsia_archer_livljobs4.rst>`_ and `<EAfricaSurge_jobs4_archer.rst>`_

Because I’m lazy and don’t get much joy from complicated compilation procedures,
I’ve wrote a script which does most of the nitty gritty for you – so you can
make a cup of tea while it’s working hard. However, you’re not excluded from
 doing any work. So here are your preliminary jobs (on archer).

---

1)	Go into archer and go to your working directory::

  cd /work/n01/n01/$USER

2)	Make a directory for your project and go into it. Call your project
 something other than ``MY_CONFIG`` if you don’t want this as your project name::

  export CONFIG=MY_CONFIG
  mkdir $CONFIG
  cd $CONFIG

3)	In this directory, copy in all the *.sh scripts::

  cp /work/n01/n01/ashbre/CARIB/START_FILES_SH/* .

4)	Go into make_paths.sh and change ``export CONFIG=   `` to your project name e.g.::

  sed -i  '/export CONFIG/c\export CONFIG=MY_CONFIG'  make_paths.sh

5)	Copy across my ``START_FILES`` directory, then go into START_FILES::

  cp -r /work/n01/n01/ashbre/CARIB/START_FILES_COPY START_FILES
  cd  START_FILES

6)	Here, you’ll find a global co-ordinate file coordinates_ORCA_R12.nc,
interrogate this file to see where your subdomain fits. My subdomain is between
 95W (i=2306) and 55W (i=2786) 8N (j=1579) and 25N (j=1846).


7)	 Edit ``namelist.input``, changing the following lines based on the
co-ordinate indexes you found::

    vi  namelist.input
    ...
    imin = 2245
    imax = 2785
    jmin = 1555
    jmax = 1893

8)	 Go into ``namelist_cfg`` and edit the co-ordinate indexes::

  vi  namelist_cfg
  ...
  jpidta      =     544   !  1st lateral dimension ( >= jpi )
  jpjdta      =     342   !  2nd    "         "    ( >= jpj )
  jpkdta      =      75   !  number of levels      ( >= jpk )
  jpiglo      =     544   !  1st dimension of global domain --> i =jpidta
  jpjglo      =     342   !  2nd    -                  -    --> j  =jpjdta

Note, for reasons I don’t understand::

  jpidta=(imax-imin)+4        (your own imin and imax from instruction 6)
  jpjdta=(jmax-jmin)+4
  jpiglo=jpidta
  jpjglo=jpjdta

9)	Go the GEBCO and download bathymetry data a little bigger than your domain.
Head to https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/
Then put this file into the ``START_FILES`` directory. It should have the
filename ``GRIDONE_2D...``


10)	Everything should be ready to go. You should have signed up for a NEMO
account already, if not, sign up at http://forge.ipsl.jussieu.fr/nemo/wiki/Users
as you’ll need to enter your password for first time use. If you have, then run
 the script::

  cd /work/n01/n01/$USER/$CONFIG/
  ./main.sh

11)	This should put your bathymetry files, coordinate files etc into the INPUT
directory. Note this it will take some time to install xios (30 mins) and
 compile nemo (15 mins) etc.

12)	 If you want to run my Caribbean instructions, copy and paste the following
 blindly::

  cd /work/n01/n01/$USER
  export CONFIG=ASH_CARIB
  mkdir $CONFIG
  cd $CONFIG
  cp /work/n01/n01/ashbre/CARIB/START_FILES_SH/* .
  sed -i  '/export CONFIG/c\export CONFIG=ASH_CARIB'  make_paths.sh
  cp -r /work/n01/n01/ashbre/CARIB/START_FILES_COPY START_FILES
  cd  START_FILES
  cd /work/n01/n01/$USER/$CONFIG
  ./main.sh
