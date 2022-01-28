from spectralradex import radex
from multiprocessing import Pool

import time

# Single run using just the basic run() method from Spectral Radex

params = radex.get_default_parameters()
# params["molfile"] = "co.dat"
# output = radex.run(params)
# print(output)
# params = radex.get_default_parameters()

# #try to exceed temperature limit
# params["tkin"]=1e5
# output=radex.run(params)
# print(output)


# # Try with no collisional partners
# params["tkin"]=20.0
# params["h2"]=0.0
# params["e-"]=100
# output=radex.run(params)
# print(output)

# # check run_grid() method and how to use or not use multiprocessing with it.
params["tkin"]=[10,50]
res=radex.run_grid(params,
             target_value="FLUX (K*km/s)", pool=Pool(4))
print(res)
try:
    res=radex.run_grid(params,
                 target_value=1, pool=Pool(4))
except:
    pass
res=radex.run_grid(params, target_value="FLUXy (K*km/s)", pool=Pool(4))
print(res)