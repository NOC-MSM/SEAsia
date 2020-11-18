===========================
Compiling the Documentation
===========================

How to edit the docs:

``cd /Users/jeff/GitHub/NEMO_docs/docs``

* Edit the ``*rst`` files:

Really good ``*.rst`` cheatsheet ``https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst``

Automatically generating the HTML
=================================

* Preview the ``*.rst`` files:
 ``$ make html``
* Commit change to github and sync NEMO-docs:
* Check files in readthedocs: http://nemo-reloc.readthedocs.io/en/latest/?# *(Not 100% sure about this link as I don't have t'internet)*

Automatically generating the PDF
================================

Having install `rst2pdf` (https://github.com/rst2pdf/rst2pdf) (in a conda environment)::

  cd recipes/docs
  sphinx-build -b pdf source build
