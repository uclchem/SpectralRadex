Radex solves the radiative transfer problem without assuming LTE and then provides various line properties. Two of those properties are key:
    1. The excitation temperature of each line
    2. The optical depth at line centre
    
The first is the excitation temperature of an LTE model that would give the same amount of emission. Which means we can use it and the optical depth to get the brightness temperature as a function of velocity.

The brightness temperature is:
$$T_B = [J_{\nu}(T_{ex})-J_{\nu}(T_{BG})](1-\exp(-\tau_v))$$

Where $T_{ex}$ is the excitation temperature and $T_{BG}$ is the background temperature, likely 2.73 K. So RADEX provides the value of $T_{ex}$ such that this equation gives the same value as the RADEX model. It also provides us with the optical depth at line centre, $\tau_0$ so we simply need to calculate $tau_v$ assuming a gaussian line profile:

$$\tau_v = \tau_0 e^{\left(-4ln(2)\frac{(v-v_0)^2}{\Delta v^2}\right)}$$

Finally, for CCH in particular, we need to consider what to do with overlapping lines as the hyperfine lines commonly overlap. We follow [Hsieh et al 2015](https://iopscience.iop.org/article/10.1088/0004-637X/802/2/126) and use an opacity weighted radiation temperature:


$$T_B = \left(\frac{\Sigma_i J{\nu}(T^i_{ex})\tau^i_v}{\Sigma_i \tau^i_v}-J_{\nu}(T_{BG})\right)(1-\exp(-\tau_v))$$

We can multiply $T_B$ by the filling factor to get the main beam temperature.