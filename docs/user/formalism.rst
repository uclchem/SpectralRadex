.. _formalism:

Formalism
=============

RADEX
---------
RADEX is a non-LTE radiative transfer solver that calculates the intensities of molecular lines assuming an homogeneous medium with a simple geometry. For a full description of the code, please see  `their release paper <http://dx.doi.org/10.1051/0004-6361:20066820>`_. The RADEX code has been modified to meet modern Fortran specifications but is otherwise unchanged in SpectralRadex.

Spectral Modelling
-------------------

In order to calculate the emission from a molecular transition as a function of frequency, we need the excitation temperature and the optical depth as a function of velocity. This allows us to calculate the brightness temperature as a function of velocity:

.. math::
	T_B = [J_{\nu}(T_{ex})-J_{\nu}(T_{BG})](1-\exp(-\tau_v))

Where  :math:`T_{ex}` is the excitation temperature and  :math:`T_{BG}` is the background temperature, likely 2.73 K. In LTE, the optical depth at line centre can be calculated from the column density and Boltzmann distribution whilst the excitation temperature is assumed to be the LTE temperature. We can then calculate  :math:`\tau_v` assuming a gaussian line profile:

.. math::
	\tau_v = \tau_0 e^{\left(-4ln(2)\frac{(v-v_0)^2}{\Delta v^2}\right)}

However, using RADEX, we can do better than to assume LTE. For a given set of physical parameters RADEX will provide the optical depth at line centre for every transition and the excitation temperature that gives the correct brightness temperature at line centre.
    
Thus, rather than using our gas kinetic temperature and an LTE derived :math:`\tau_0`, we can take the values for each line from an appropriate RADEX output. In the high density limit, this tends to the LTE solution but at lower densities it can deviate significantly.

In SpectralRadex, we do this for each transition in a collisional datafile between a minimum and maximum frequency set by the user. :math:`T_B` is calculated as a function of frequency for each line and then combined to give the overall spectrum of the molecule. 


Finally, we need to consider what to do with overlapping lines. We follow `Hsieh et al 2015 <https://iopscience.iop.org/article/10.1088/0004-637X/802/2/126>`_ and use an opacity weighted radiation temperature:

.. math::
	T_B = \left(\frac{\Sigma_i J{\nu}(T^i_{ex})\tau^i_v}{\Sigma_i \tau^i_v}-J_{\nu}(T_{BG})\right)(1-\exp(-\tau_v))

We can multiply :math:`T_B` by the filling factor to get the main beam temperature.
