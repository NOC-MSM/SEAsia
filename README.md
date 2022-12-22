[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7473198.svg)](https://zenodo.org/record/7473198))

## NEMO regional configuration of the Severn Estuary


This model configuration has been developed in order to test the NEMO model in a tidal coastal area at 500 m resolution 
(in preparation for the UK wide NEMO 500 m model development).


### Model Summary

Severn Estuary in the UK  (-5E to -2E and 50.1N to 51.8N)

The model grid target resolution 500m; with 31 sigma vertical levels.
Featuring:

* FES2014 tides
* ERA5 atmospheric forcing
* CMEMS open boundaries & Initial conditions
* Freshwater forcing 
* Wetting and Drying 

![Severn bathymetry](https://github.com/JMMP-Group/SEVERN-SWOT/wiki/FIGURES/SEVERN-SWOT_bathy.png)

### Model Setup

The following process is followed to build and get started with this configuration

``git clone https://github.com/JMMP-Group/SEVERN-SWOT.git``

Then follow descriptions in: https://github.com/JMMP-Group/SEVERN-SWOT/wiki


### Experiment Summary
In the master branch:
* ``EXP_unforced_constTS``
Constant initial T and S. No tides. No met.

* ``EXP_barotropicTide``
Constant initial T and S. Only tidal forcing. No met.

* ``EXP_Tide_ERA5``
Constant initial T and S. Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.).

* ``EXP_Tide_ERA5_BDY``
Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.), 3D boundary conditions and initial conditions.

In the WAD branch:
* ``EXP_barotropicTide_WAD``
Constant initial T and S. Only tidal forcing. No met. With Wetting And Drying.

* ``EXP_Tide_ERA5_WAD``
Constant initial T and S. Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.). With Wetting And Drying.

* ``EXP_Tide_ERA5_River_IC_BDY_WAD``
Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.), 3D boundary conditions and initial conditions.
River forcing. With Wetting And Drying.

