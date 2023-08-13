.. NIRSimple documentation master file

Welcome to NIRSimple's documentation!
=====================================

`NIRSimple` is a Python 3 package for fNIRS with a focus on signal processing. It
gives control on preprocessing and processing, enabling to adjust parameters
and use different methods from the scientific literature.

This package handles data as numpy arrays. It can be used with
`MNE <https://mne.tools/stable/index.html>`_ pipelines by creating MNE-specific
objects (Raw, Epochs, ...) from numpy arrays.

Features:
  * conversion from light intensity to optical density changes
  * conversion from optical density changes to hemoglobin concentration changes
    with the modified Beer-Lambert law (options for different extinction
    coefficient tables)
  * signal correction with correlation based signal improvement (CBSI)


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   modules
   example


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
