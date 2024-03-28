[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6483231.svg)](https://zenodo.org/record/6483231)

## NEMO regional configuration of South East Asia

An example configuration of SE Asia, demonstrating how to setup new regional domains in the NEMO framework.
This model configuration has been developed through the ACCORD (Addressing Challenges of Coastal Communities through Ocean Research for Developing Economies) Project, funded by [Natural Environment Research Council, under a National Capability Official Development Assistance](http://gotw.nerc.ac.uk/list_full.asp?pcode=NE%2FR000123%2F1).


### Model Summary

A specific region of focus includes exploring South East Asia (75E to 135E and -20N to +25N)

The model grid has 1/12&deg; lat-lon resolution and 75 hybrid sigma-z-partial-step vertical levels. Featuring:

* FES2014 tides
* Boundary conditions
* Freshwater forcing 
* ERA5 wind and sea level pressure

![SE Asia bathymetry](https://github.com/NOC-MSM/SEAsia/wiki/FIGURES/ACCORD_SEAsia_bathy.png)

### Model Setup


The following process is followed to build and get started with this configuration

``git clone https://github.com/NOC-MSM/SEAsia.git``

Then follow descritptions in: https://github.com/NOC-MSM/SEAsia/wiki

The example is based on NEMO v4.0.6 and XIOS v2.5

### Domain file

The domain configuration file can be [downloaded](https://gws-access.jasmin.ac.uk/public/jmmp/SEAsia_R12/)
