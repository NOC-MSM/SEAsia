===========================
Compiling the Documentation
===========================

How to edit the docs:

``cd DOCS``

* Edit the ``*rst`` files:

Really good ``*.rst`` cheatsheet ``https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst``

Automatically generating the HTML
=================================

* Preview the ``*.rst`` files:
 ``$ make html``
* Commit change to github and sync NEMO-docs:
* View at ``build/html/index.html`` in your favourite browser.

Automatically generating the PDF
================================

Having installed `rst2pdf` (https://github.com/rst2pdf/rst2pdf) (in a conda
environment)::

  cd DOCS
  sphinx-build -b pdf source build

View pdf in ``build``.

Automatically generate latex files
==================================

(Only tried within the above conda environment, having installed `rst2pdf`)::

  sphinx-build -b latex DOCS build/latex
