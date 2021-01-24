.. _install:

Installation
============
Pypi
----------------
We recommend the simple approach of using pypi:

``pip install spectralradex``

Note, you may receive an error of the kind 


``could not find a version that satisfies the requirement spectralradex``


which is caused by a bug in how the requirements for spectralradex is configured. Pip will be in the process of installing a library (pandas/numpy) when this occurs and it can be solved by installing that library through pip before installing SpectralRadex.


Manual Install
----------------
However, if you wish to install manually, clone the repo and from the main directory run the following

``python3 setup.py install``

optionally, specify a path for the installation using

``python3 setup.py install --prefix=/path/to/my/install``

making sure that the install path is part of your PYTHONPATH environmental variable. You'll need to add the .egg directory eg. ```path/to/my/install/lib/python3.7/site-packages/spectralradex-0.0.2-py3.7-linux-x86_64.egg```