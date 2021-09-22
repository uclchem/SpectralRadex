from spectralradex import radex
import numpy as np
from time import perf_counter


params=radex.get_default_parameters()
def new():
    start=perf_counter()
    radex.new_run("co.dat",30.0,cdmol=1e16,nh=0.0,nh2=1e5,op_ratio=3.0,ne=0.0,nhe=0.0,nhx=0.0,
                    linewidth=1.0,fmin=0.0,fmax=500.0,tbg=2.73,geometry=1)
    stop=perf_counter()
    return stop-start

def old():
    start=perf_counter()
    radex.run(params)
    stop=perf_counter()
    return stop-start

new_times=[new() for it in range(100)]
old_times=[old() for it in range(100)]

print(f"new runs in {np.mean(new_times)}")
print(f"old runs in {np.mean(old_times)}")