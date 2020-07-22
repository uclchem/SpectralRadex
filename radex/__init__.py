from .radexwrap import *
from pandas import DataFrame

def run(parameters,output_file=None):
	columns=['E_UP (K)','freq','WAVEL (um)','T_ex','tau',
			'T_R (K)','POP UP','POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']
	  
	nlines,qup,qlow,output=radex(parameters)
	output=DataFrame(columns=columns,data=output[:,:nlines].T)
	output["QN Upper"]=qup.reshape(-1, 6).view('S6')[:nlines]
	output["QN Lower"]=qlow.reshape(-1, 6).view('S6')[:nlines]
	output["QN Upper"]=output["QN Upper"].map(lambda x: str(x, 'utf-8')).str.strip()
	output["QN Lower"]=output["QN Lower"].map(lambda x: str(x, 'utf-8')).str.strip()
	if output_file is not None:
		output.to_csv(output_file,index=False)
	return output


def get_default_parameters():
	parameters={
		"molfile":"co.dat",
		"tkin":30.0,
		"tbg":2.73,
		"cdmol":1.0e13,
		"h2":1.0e5,
		"h":0.0,
		"e-":0.0,
		"p-h2":0.0,
		"o-h2":0.0,
		"h+":0.0,
		"linewidth":1.0,
		"fmin":0.0,
		"fmax":3.0e7
	}
	return parameters