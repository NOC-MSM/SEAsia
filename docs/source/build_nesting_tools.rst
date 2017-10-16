Build AGRIF tools
=================

There was a problem with the GRIDGEN coordinate.nc output file in that the spacing
variables e.g. e1t where not scaled with the zoom from the parent grid.

James' found a fix using the AGRIF tools.

Move to the TOOL directory e.g.::

  export TDIR=/work/n01/n01/jelt/LBay/trunk_NEMOGCM_r8395/TOOLS
  cd $TDIR

Copy in the right ARCH file::

  cp /work/n01/n01/jelt/ARCH/arch-XC_ARCHER_INTEL_NOXIOS.fcm ../ARCH/.

Build the NESTING tool::

  ./maketools -n NESTING -m XC_ARCHER_INTEL_NOXIOS -j 6

This makes a number of executables in NESTING::

  cd NESTING

Link in parent coordinates file (from somewhere)::

  ln -s $INPUTS/coordinates_ORCA_R12.nc $TDIR/NESTING/.
