# NEMO_cfgs
NEMO configurations

Each configuration directory should be laid out in the following manner:

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
|____README.md
</pre>
