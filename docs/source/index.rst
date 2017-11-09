.. NEMO-RELOC documentation master file, created by
   sphinx-quickstart on Wed Oct  5 09:17:15 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NEMO-RELOC's documentation!
======================================

There are two README.md files is this project. Their contents are loaded here and should describe the project (rather than have intro text duplicated in this index file). *It would have been better if the README.md files loaded their contents here without needing to click on a link. Like <INCLUDE=""> in HTML-speak. I don't know how to do this in Markdown)* **EDIT** I think a simple pre-"make html" command would be needed to paste in the contents of the README files.


`Include the root README.md file <../../../README.md>`_

`Include the tools README.md file <../../../tools/README.md>`_

Note however that the markdown in the README.md files is not converted into html with the ``make html`` rendering command.


Complete Recipes:
=================

When a new note file is created manually add it to this list to aid finding it.

.. toctree::
   :maxdepth: 1

   template

   LBay_archer_livljobs4

   LBay_180m

   SEAsia_archer_livljobs4

   SWPacific_archer_livljobs4

   EAfrica


Modules:
========
.. toctree::
  :maxdepth: 1

  install_nrct

  build_and_create_coordinates

  build_opa_orchestra

  build_siren_tools

  build_XIOS2

  MPP_decomp_land_suppression

  rebuild_and_inspect_NEMO_output


Tools:
======

.. toctree::
  :maxdepth: 1

  quickplotNEMO.py

  SEAsia_SSH_anim.py


WIP:
====

.. toctree::
  :maxdepth: 1

  Solent

  todo




Trouble Shooting:
=================
.. toctree::
  :maxdepth: 2

  trouble_shooting

  updating_namelist_cfg_from_v3.6_to_v4


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
