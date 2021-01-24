from . import radex
from pandas import DataFrame,read_csv
import numpy as np
import os 

package_directory = os.path.dirname(os.path.abspath(__file__))


light_speed=2.99792e5
planck=6.62607e-34
boltzman=1.38e-23

#c^3/8pi and geometric factor for gaussian integral
lte_constants=(light_speed*light_speed*light_speed)/(8.515736*np.pi) 


#calculate the brightness temperature for each observed frequency for a given set of parameters
def calculate_spectrum(obs_freqs,v0,radex_params):
	#user supplies the observed frequency so doppler shift to emitted
	#tau dist makes this unnecessary
	emit_freqs=obs_freqs*(1.0+v0/light_speed)

	#we'll return a dataframe of Frequency, Intensity
	new_df=DataFrame({"Frequency":obs_freqs})
	new_df["Intensity"]=0.0
	new_df=new_df.sort_values("Frequency",ascending=False)
	#solve the radex model and get all line properties
	tau_0_df=get_radex_taus(radex_params)

	delta_v=radex_params["linewidth"]

	#now loop through line and build up the tau weighted radiation temperature average
	for i,line in tau_0_df.iterrows():
		#get the relative velocity of all the emitting frequencies
		velocities=((line["freq"]/new_df["Frequency"])-1.0)*light_speed

		if velocities[1]-velocities[0]>delta_v:
			print("Velocity bins larger than linewidth")
		taus=get_tau_dist(v0,delta_v,line["tau"],velocities)

		#store tau weighted radiation temp
		new_df[f"{line.freq:.3f}"]=rad_temp(line["T_ex"],emit_freqs)*taus
		#and add tau to running total
		new_df["Intensity"]+=taus


	#sum our tau weighted temperatures and divide by sum of taus
	line_cols=[x for x in new_df if x not in ["Intensity","Frequency"]]
	new_df["temp"]=new_df[line_cols].sum(axis=1)/new_df["Intensity"]
	#now get brightness temperature as a function of frequency
	new_df["Intensity"]=(new_df["temp"]-rad_temp(2.73,emit_freqs))*(1.0-np.exp(-new_df["Intensity"]))
	new_df["Intensity"]=new_df["Intensity"].fillna(0.0)
	return new_df[["Frequency","Intensity"]].sort_values("Frequency")

#This runs radex to get the excitation temperature and optical depth for every line
def get_radex_taus(params):
	columns=['E_UP (K)','freq','WAVEL (um)','T_ex','tau','T_R (K)','POP UP',
			'POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']
	  
	output=radex.run(params)
	idx=(output["freq"]>0.0) & (output["tau"]>0)
	return output.loc[idx,["freq","tau","T_ex"]]

#calculate the optical depth at all our observed frequencies for a given line
#based on velocity relative to line centre
def get_tau_dist(v0,delta_v,tau_0,velocities):
	taus=np.exp(-4.0*np.log(2.0)*(((velocities-v0)**2.0)/(delta_v*delta_v)))
	taus=taus*tau_0
	return taus

#Calculate raditaion temperature for a given excitation temperature
def rad_temp(t,frequencies):
	hvk=(planck*frequencies*1.0e9)/boltzman
	hvk=hvk/(np.exp(hvk/t)-1)
	return hvk




def get_lte_taus(N,T,delta_v):
	tau_df=DataFrame(line_df["Frequency"])
	line_df["tau"]=0.0
	line_df["T_ex"]=T
	
	line_df["tau"]=line_df["Frequency"].values*1.0e9
	line_df["tau"]=line_df["tau"].pow(3)
	
	#then divide column density in m-2 by partition function
	N_Q=(1e4*(N))/part_func(T)
	
	
	line_df["tau"]=line_df["Aij"]*line_df["gu"]/(line_df["tau"]*delta_v*1.0e3)

	line_df["boltz"]=np.exp(-line_df["E_L"].values/T)
	line_df["boltz"]=line_df["boltz"]-np.exp(-line_df["E_U"].values/T)

	line_df["tau"]=constants*line_df["tau"]*N_Q*line_df["boltz"]
	return line_df[["Frequency","tau","T_ex"]].rename({"Frequency":"freq"},axis=1)


def spectral_error(x,obs_freqs,obs_intensity,noise,ncomponents=1):
	#straight away we can get first component
	intensity=calculate_spectrum(*x[:6],obs_freqs)
	#then if ncomponents > 1 we can add additional
	for i in range(1,ncomponents):
		y=x[6*i:6*(i+1)]
		fit=calculate_spectrum(*y,obs_freqs)
		intensity=intensity+fit

	indx=np.where(intensity>0.0)[0]
	error=obs_intensity[indx]-intensity[indx]
	error=(error*error)/(noise*noise)
	chi=np.sum(error)/len(obs_intensity[indx])
	return chi
	
def noise_from_spectrum(intensities):
	noise=np.median(intensities)
	ints=intensities[np.where(intensities<noise)[0]]
	noise=np.mean((noise-ints)**2.0)
	noise=np.sqrt(noise)
	return noise