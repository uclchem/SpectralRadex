.. _trouble:

Trouble Shooting
=================


Malloc() Error
--------------
SpectralRadex can return results for up to 500 transitions. This number is hard coded because Fortran cannot use variable sized arrays as part of the python interface and so we had to choose a number which trades off a reasonably high maximum with the fact a massive array would take up a lot of memory without being needed in 99% of cases. However, if you run a species such as CH3OH over a very large frequency range, you can have more transitions than this. This will result in an error

.. code:: shell

    malloc(): corrupted top size
    Abort (core dumped)

which can be resolved by setting fmin and fmax such that there are fewer than 500 transitions in the range of interest. If you require more than 500 transitions, please contact us via github or email.

pip cannot find version
-----------------------
**Could not find a version that satisfies the requirement spectralradex**

If you get an error like this when trying to install spectralradex, you may need to install it from source. This can be done by running

.. code:: shell
    git clone https://github.com/uclchem/SpectralRadex.git
    pip install ./SpectralRadex


This happens because we use Github Actions to pre-build the library for various python versions and OS combinations. Not every combination is possible and so if your combination doesn't exist, you need to build it from source.

Mac Issues
-----------
**Import Error... library not loaded**


Recent updates to Mac OS have resulted in many Mac user's python distributions expecting the standard Fortran libraries to be in one place when they're actually in another. The resulting error message looks like

.. parsed-literal::
    Exception has occurred: ImportErrordlopen(/usr/local/lib/python3.9/site-packages/radexwrap.cpython-39-darwin.so, 2): Library not loaded: /usr/local/opt/gcc/lib/gcc/10/libgfortran.5.dylib   Referenced from: /usr/local/lib/python3.9/site-packages/radexwrap.cpython-39-darwin.so   Reason: image not found


In this case, SpectralRadex wants the libgfortran.5.dylib library and can't find it. You can solve this with locate

.. code:: shell

    locate libgfortran.5.dylib

which will tell you the actual location of the required library and then you can create a symbolic link to the expected location.

.. code:: python

    ln /actual/path/to/libgfortran.5.dylib /the/path/python/expected/libgfortran.5.dylib
