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

![tidal SSH(m)](https://www.dropbox.com/s/25yb47pjibtfkpz/SEAsia_domain.png?dl=0)

### Model Setup

The following code was used in this configuration:

svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM@8395 trunk_NEMOGCM_r8395

svn co -r1080 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1080

The initial conditions and boundary data can be downloaded from JASMIN:

http://  ...



Repository structure

The repository is structure as follows: **UPDATE**
<pre>
MYCONFIG
|____ARCH
| |____arch-XC_ARCHER.fcm
|____arch_xios
| |____arch-XC_ARCHER.env
| |____arch-XC_ARCHER.fcm
| |____arch-XC_ARCHER.path
|____cpp_MYCONFIG.fcm
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
|____MY_SRC
| |____*.F90
|____INPUTS
| |____namelist.bdy
|____README.md
</pre>
