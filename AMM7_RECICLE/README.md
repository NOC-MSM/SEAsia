# AMM7_RECICLE

The following code was used in this configuration:

git@gitlab.ecosystem-modelling.pml.ac.uk:momm/NEMO-shelf.git

The initial conditions and boundary data can be downloaded from JASMIN:

http://gws-access.ceda.ac.uk/public/recicle/config

NB This recipe has be written with the ARCHER HPC INTEL environment in mind.

```
# Change to some working directory of choice
cd $WORK_DIR

# Clone SSB code 
git clone git@gitlab.ecosystem-modelling.pml.ac.uk:momm/NEMO-shelf.git
cd NEMO-shelf
git checkout e3a60987d075d05eaacb9ba6b7b8116b2f99a0a6

# Now change to CONFIG directory
cd NEMOGCM/CONFIG

# Checkout configuration directory structure
git init .
git remote add origin git@github.com:NOC-MSM/NEMO_cfgs.git
git config core.sparsecheckout true
echo "AMM7_RECICLE/*" >> .git/info/sparse-checkout
git pull --depth=1 origin master

# Add it to the configuration list
echo "AMM7_RECICLE OPA_SRC" >> cfg.txt
```

At this point you can checkout and compile XIOS or use a version you already have. If you're starting from scratch:

```
# Choose an appropriate directory for your XIOS installation
cd $TO_XIOS_DIR
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@703
cd xios-1.0
mv $WORK_DIR/NEMO-shelf/NEMOGCM/CONFIG/AMM7_RECICLE/arch_xios/* ./arch
rm -rf $WORK_DIR/NEMO-shelf/NEMOGCM/CONFIG/AMM7_RECICLE/arch_xios
./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 4
```

You can fold the ```make_xios``` command into a serial job. NB ```$NETCDF_DIR``` and ```$HDF5_DIR``` must be part of your environment. This should be the case if you've used ```modules``` to setup the netcdf and hdf5 e.g. 

```
module swap PrgEnv-cray PrgEnv-intel
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
```

Next, compile the NEMO code itself. First we copy the arch files into the appropriate directory.

```
cd $WORK_DIR/NEMO-shelf/NEMOGCM/CONFIG/AMM7_RECICLE/
mv ARCH/* ../../ARCH
rm -rf ARCH
```

NB you will either have to edit ```../../ARCH/arch-XC_ARCHER_INTEL_XIOS1.fcm``` replacing ```$XIOS_DIR``` with the expanded ```$TO_XIOS_DIR/xios-1.0``` or define ```$XIOS_DIR``` in your environment.

```
cd ../
./makenemo -n AMM7_RECICLE -m XC_ARCHER_INTEL_XIOS1 -j 4
```

That should be enough to produce a valid executable. Now to copy the forcing data from JASMIN. 

```
cd AMM7_RECICLE/ENSEMBLE_INPUTS
wget -r -np -nH --cut-dirs=3 -erobots=off --reject="index.html*" http://gws-access.ceda.ac.uk/public/recicle/config/
```

And finally link the XIOS binary to the configuration directory.

```
cd ../ENSEMBLE_CONTROL
rm xios_server.exe
ln -s $TO_XIOS_DIR/xios-1.0/bin/xios_server.exe xios_server.exe
```

Edit and run the ```ensemble.pbs``` script in ```../ENSEMBLE_CONTROL``` accordingly.
