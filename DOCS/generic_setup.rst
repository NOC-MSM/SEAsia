Machines used: CentOS7 linux box, Cray XC30 HPC (ARCHER)

Having cloned the NEMO-RELOC repository, the  entire build process can be run
with a single script ``SCRIPTS/main_setup.sh``. However the process is taken
through step by step in the following. These first steps are generic and
 only need be done once for any number of configurations on the same architecture
and code base. The second group of steps are configuration specific steps.


TEST :cite:`Gil:02`

===========================================
Preparing for a new NEMO vp4 configuration
===========================================


.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: ****************
  :end-at: ./make_nemo.sh


a. Clone the NEMO-RELOC repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step is set up a new configuration directory by git cloning the
repository (e.g. to your workspace)::

  cd /work/n01/n01/$USER
  git clone https://github.com/NOC-MSM/NEMO-RELOC


b. Set script paths, make paths and directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure that the first two lines  in ``SCRIPTS/make_paths_setup.sh`` are appropriately
defined:

.. literalinclude:: ../SCRIPTS/make_paths_setup.sh
  :start-at: export CONFIG
  :end-at: export WORK

Then ``SCRIPTS/make_paths_setup.sh`` and ``SCRIPTS/make_directories_setup.sh`` will create
the expected directory structure.

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: echo "Making Paths"
  :end-at: make_directories_setup.sh

c. Build XIOS
^^^^^^^^^^^^^

Then we build XIOS for managing the input/output on distributed architecture.
Note that the modules need to be loaded appropriately for the architecture. The
settings here are for ARCHER:

.. literalinclude:: ../SCRIPTS/make_xios.sh

The script builds two versions of XIOS. XIOS2.5 is used for the NEMO simulations
but an earlier version XIOS1 (where performance is not an issue) is used for
some of the tooling required to build the configuration.

This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
  :start-at: Installing XIOS
  :end-at: make_xios


d. Build TOOLS
^^^^^^^^^^^^^^

The next stage is to build tools that will be used for a number of processes
during the build. Specifically tools to assist with domain configuration file
generation (NESTING, which can be used to interpolate existing grids, and
DOMAINcfg, which builds the domain configuration netCDF file). At the same time
we build tools to reconstruct NEMO output (REBUILD_NEMO) and tools for
generating grid appropriate weights files to use with external surface forcing
data (WEIGHTS):

.. literalinclude:: ../SCRIPTS/make_tools.sh

This is executed by ``SCRIPTS/main_setup.sh``:

.. literalinclude:: ../SCRIPTS/main_setup.sh
 :start-at: various grid tools
 :end-at: ./make_tools.sh


e. Build NEMO
^^^^^^^^^^^^^

Then build the NEMO executable. This is the dynamical ocean core the NEMO
framework. Technically OPA is the ocean dynamical core and NEMO is the framework
it sits within (including, for example, modules for ice, biogeochemistry etc).

In this documentation there are examples for a NEMO physics only build (called
``SCRIPTS/make_nemo.sh``) and a NEMO-FABM-ERSEM build
(called ``SCRIPTS/make_nemo_fabm_ersem.sh``). Here we walk through the physics
only build:

.. literalinclude:: ../SCRIPTS/make_nemo.sh


This is executed by ``SCRIPTS/main_setup.sh`` (commenting out the appropriate
make NEMO script):

.. literalinclude:: ../SCRIPTS/main_setup.sh
 :start-at: Installing NEMO-FABM-ERSEM
 :end-at: ./make_nemo.sh
