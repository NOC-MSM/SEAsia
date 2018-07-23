# AMM7_RECICLE

The following code was used in this configuration:

git@gitlab.ecosystem-modelling.pml.ac.uk:momm/NEMO-shelf.git

The initial conditions and boundary data can be downloaded from JASMIN:

http://gws-access.ceda.ac.uk/public/recicle/config

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
echo "AMM7_RECICLE OPA_SRC" >> cfg.txtgit init AMM7_RECICLE
```

At this point you can checkout and compile XIOS or use a version you already have. If you're starting from scratch:

```
# Choose an appropriate directory for your XIOS installation
cd $XIOS_DIR
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@703
cd xios-1.0
cp $WORK_DIR/NEMO-shelf/NEMOGCM/CONFIG/AMM7_RECICLE/arch_xios/* ./arch
./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 4
```

You can fold the ```make_xios``` command into a serial job.


