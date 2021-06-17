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


...

### Repository structure

