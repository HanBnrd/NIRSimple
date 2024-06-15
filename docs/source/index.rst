.. NIRSimple documentation master file

NIRSimple
=========

`fNIRS signal processing simplified`

`NIRSimple` is a Python 3 package for fNIRS with a focus on signal processing. It
gives control on preprocessing and processing, enabling to adjust parameters
and use different methods from the scientific literature.

This package handles data as numpy arrays. It can be used with
`MNE <https://mne.tools/stable/index.html>`_ pipelines by creating MNE-specific
objects (Raw, Epochs, ...) from numpy arrays.

Features:
  * conversion from light intensity to optical density changes
  * differential pathlength factor (DPF) from wavelength and age
  * conversion from optical density changes to hemoglobin concentration changes
    with the modified Beer-Lambert law (options for different extinction
    coefficient tables)
  * signal correction with correlation based signal improvement (CBSI)

.. role::  raw-html(raw)
    :format: html

:raw-html:`&rarr;` `Source code on GitHub <https://github.com/HanBnrd/NIRSimple>`_

.. image:: https://img.shields.io/badge/doi-10.3389%2Ffnrgo.2023.994969-blue
  :target: https://doi.org/10.3389/fnrgo.2023.994969

.. image:: https://img.shields.io/badge/license-MIT-lightgrey
  :target: https://github.com/HanBnrd/NIRSimple/blob/master/LICENSE

.. image:: https://img.shields.io/github/actions/workflow/status/HanBnrd/NIRSimple/sphinx.yml?label=pipeline
  :target: https://github.com/HanBnrd/NIRSimple

.. image:: https://img.shields.io/pypi/v/nirsimple
  :target: https://pypi.org/project/nirsimple/

.. toctree::
   :maxdepth: 1
   :caption: Contents

   install
   modules
   example


Acknowledgements
----------------

Until there is an article specifically on `NIRSimple`, please cite `this article <https://doi.org/10.3389/fnrgo.2023.994969>`_.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
