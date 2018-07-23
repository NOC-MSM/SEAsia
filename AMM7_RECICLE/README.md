# AMM7_RECICLE

The following code was used in this configuration:

git clone git@gitlab.ecosystem-modelling.pml.ac.uk:momm/NEMO-shelf.git
cd NEMO-shelf
git checkout e3a60987d075d05eaacb9ba6b7b8116b2f99a0a6

The initial conditions and boundary data can be downloaded from JASMIN:

http://

```
git init AMM7_RECICLE
cd AMM7_RECICLE
git remote add origin git@github.com:NOC-MSM/NEMO_cfgs.git
git config core.sparsecheckout true
echo "AMM&_RECICLE/*" >> .git/info/sparse-checkout
git pull --depth=1 origin master
```
