SpectralRadex
==============

SpectralRadex is a python module made up of two parts: a RADEX wrapper and a spectral model.

A number of libraries exist for the first purpose. However, SpectralRadex uses F2PY to compile a version of RADEX written in modern Fortran. As a result, running a RADEX model creates no subprocesses, no text files, and can be easily parallelized. Check our tutorials section for an example of running a grid of RADEX models quickly and entirely within Python.


For the second purpose, we use the RADEX calculated line opacities and excitation temperatures to calculate the brightness temperature as a function of frequency. This allows observed spectra to be modelled in python in a non-LTE fashion. See our spectral modelling tutorial for more.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user/install
   user/formalism
   user/referencing
   user/trouble

.. toctree::
	:maxdepth: 1
	:caption: API

	api/spectralradex
	api/radex


.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/radex
   tutorials/spectralmodelling

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
