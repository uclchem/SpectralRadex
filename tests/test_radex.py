from spectralradex import radex

params=radex.get_default_parameters()
params["molfile"]="catom.dat"
radex.run(params,"test.csv")