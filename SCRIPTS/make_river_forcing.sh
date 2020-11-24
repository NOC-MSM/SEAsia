#!/bin/bash

:'

*******************************
make_river_forcing.sh
*******************************

Generate the river data

River data (discharge + nutrients) was generated using the GlobalNEWS2 model.
Instructions for how to generate the rivers can be found on the
`SANH Wiki page <https://github.com/NOC-MSM/SANH/wiki/6.-Get-rivers>`_.
'
#::

  cd $RIVER

  # copy the rivers but see description in the wiki on how they have
  #been generated
  cp /work/n01/n01/jenjar93/accord/SEAsia_rivertest_v3.nc $RIVER

  cd $WORK
