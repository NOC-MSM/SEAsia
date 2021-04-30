********
# SEAsia
********

This model configuration has been developed through the ACCORD (Addressing Challenges of Coastal Communities through Ocean Research for Developing Economies) Project, funded by [Natural Environment Research Council, under a National Capability Official Development Assistance](http://gotw.nerc.ac.uk/list_full.asp?pcode=NE%2FR000123%2F1)

*************************************************
## NEMO regional configuration of South East Asia
*************************************************

### Model Summary

A specific region of focus includes exploring South East Asia (75E to 135E and -20N to +25N)

The model grid has 1/12&deg; lat-lon resolution and 75 hybrid sigma-z-partial-step vertical levels.

![SE Asia bathymetry](https://github.com/NOC-MSM/SEAsia/wiki/FIGURES/ACCORD_SEAsia_bathy.png)

### Model Setup

The following code was used in this configuration:

svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.6

svn checkout -r 1964 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5

The initial conditions and boundary data can be downloaded from JASMIN:

http://  ...

### Experiment Summary

EXP_barotropicTide
==================
Only tidal forcing. constant T and S


EXP_unforcedStrat
=================
No forcing. T(z),S(z) profiles. Clamped T(z),S(z) boundaries. Start from rest.


EXP_biogeochem
==============
...

EXP_fullforcing
===============
...

### Repository structure

The repository is structure as follows: **NEEDS TO BE UPDATED**
<pre>
MYCONFIG
|____ARCH
| |____NEMO
| | |___arch-XC_ARCHER.fcm
| |___XIOS
|   |____arch-XC_ARCHER.env
|   |____arch-XC_ARCHER.fcm
|   |____arch-XC_ARCHER.path
|
|____cpp_MYCONFIG.fcm
|
|____EXP00
| |____1_namelist_cfg
| |____1_namelist_ice_cfg
| |____1_namelist_ice_ref
| |____1_namelist_ref
| |____context_nemo.xml
| |____domain_def_nemo.xml
| |____field_def_nemo-lim.xml
| |____field_def_nemo-opa.xml
| |____field_def_nemo-pisces.xml
| |____file_def_nemo.xml
| |____iodef.xml
| |____namelist_cfg
| |____namelist_ice_cfg
| |____namelist_ice_ref
| |____namelist_pisces_cfg
| |____namelist_pisces_ref
| |____namelist_ref
| |____namelist_top_cfg
| |____namelist_top_ref
| |____runscript
|
|____MY_SRC
| |____*.F90
|
|____SCRIPTS
|
|____INPUTS (a place where forcing files are put)
|
|____DOCS
|
|____START_FILES
|
|____TOOLS
|____LICENSE
|____README.md
</pre>
