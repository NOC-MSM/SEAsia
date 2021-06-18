******************
# Relocatable NEMO
******************

An example configuration of SE Asia, demonstrating how to setup new regional domains in the NEMO framework.
This model configuration has been developed through the ACCORD (Addressing Challenges of Coastal Communities through Ocean Research for Developing Economies) Project, funded by [Natural Environment Research Council, under a National Capability Official Development Assistance](http://gotw.nerc.ac.uk/list_full.asp?pcode=NE%2FR000123%2F1).

*************************************************
## NEMO regional configuration of South East Asia
*************************************************

### Model Summary

A specific region of focus includes exploring South East Asia (75E to 135E and -20N to +25N)

The model grid has 1/12&deg; lat-lon resolution and 75 hybrid sigma-z-partial-step vertical levels. Featuring:

* FES2014 tides
* Boundary conditions from ... (in prog.)
* Freshwater forcing (in prog.)
* ERA5 wind and sea level pressure (in prog.)

![SE Asia bathymetry](https://github.com/NOC-MSM/SEAsia/wiki/FIGURES/ACCORD_SEAsia_bathy.png)

### Model Setup


The following process is followed to build and get started with this configuration

``git clone https://github.com/NOC-MSM/NEMO-RELOC.git``

Then follow descritptions in: https://github.com/NOC-MSM/NEMO-RELOC/wiki

The example is based on NEMO v4.0.6 and XIOS v2.5:



### Experiment Summary

* ``EXP_barotropicTide``
Only tidal forcing. Constant T and S

* ``EXP_unforced``
No forcing. Stratification varies only with depth. Start from rest.


...

### Repository structure

* ...
