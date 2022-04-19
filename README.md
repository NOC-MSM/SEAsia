****************************
# Severn regional NEMO model
****************************

This model configuration has been developed in order to ...

*************************************************
## NEMO regional configuration of Severn Estuary
*************************************************

### Model Summary

Severn Estaury in the UK  (-5E to -2E and 50.1N to 51.8N)

The model grid target resolution 500m; with 31 sigma vertical levels, with a baroptropic ocean. Featuring:

* FES2014 tides
* Wave coupling (in prog.)
* Freshwater forcing (in prog.)
* ERA5 wind and sea level pressure (in prog.)
* Wetting and Drying (in prog.)

![Severn bathymetry](https://github.com/JMMP-Group/SEVERN-SWOT/wiki/FIGURES/severn.png)

### Model Setup

The following process is followed to build and get started with this configuration

``git clone https://github.com/JMMP-Group/SEVERN-SWOT.git``

Then follow descritptions in: https://github.com/JMMP-Group/SEVERN-SWOT/wiki


### Experiment Summary

* ``EXP_unforced_constTS``
Constant T and S. No tides. No met.

* ``EXP_unforced_horizTS``
Constant T, S stratification. No tides. No met.

* ``EXP_barotropicTide``
Only tidal forcing. constant T and S

* ``EXP_Tide_ERA5``
Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.), constant T,S.

* ``EXP_Tide_ERA5_BDY``
Tidal forcing. ERA5 forcing (sea-level pressure, 10m winds, fluxes etc.), 3D boundary conditions and initial conditions.

### Repository structure

* ...
