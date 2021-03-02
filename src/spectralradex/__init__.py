from . import radex
from pandas import DataFrame,read_csv
import numpy as np
import os 

package_directory = os.path.dirname(os.path.abspath(__file__))


light_speed=2.99792e5
planck=6.62607e-34
boltzman=1.38e-23

#calculate the optical depth at all our observed frequencies for a given line
#based on velocity relative to line centre
def maxwellian_distribution(v0,delta_v,tau_0,velocities):
    """
    Returns the optical depth as a function of velocity, assuming gaussian line profiles and given an optical depth a line centre

    :param v0: Peak velocity of the emission
    :type v0: float

    :param delta_v: FWHM of the peaks, taken from linewidth parameter of RADEX when called via :func:`model_spectrum`
    :type v0: float

    :param tau_0: The optical depth at line centre. Taken from RADEX when called via :func:`model_spectrum`
    :type tau_0:float

    :param velocities: An iterable containing the velocity values at which to calculate tau

    :return: An array with the tau value at each velocity in velocities
    :rtype: ndarray,float
    """
    taus=np.exp(-4.0*np.log(2.0)*(((velocities-v0)**2.0)/(delta_v*delta_v)))
    taus=taus*tau_0
    return taus

#This runs radex to get the excitation temperature and optical depth for every line
def get_radex_taus(params):
    columns=['E_UP (K)','freq','WAVEL (um)','T_ex','tau','T_R (K)','POP UP',
            'POP LOW', 'FLUX (K*km/s)', 'FLUX (erg/cm2/s)']
      
    output=radex.run(params)
    idx=(output["freq"]>0.0)
    return output.loc[idx,["freq","tau","T_ex"]]

def model_spectrum(obs_freqs,v0,radex_params,tau_profile=maxwellian_distribution):
    """
    Calculates the brightness temperature as a function of frequency for given input frequencies, :math:`V_{LSR}`
    velocity and RADEX parameters.

    :param obs_freqs: An array of frequency values in GHz at which the brightness temperature should be calculated.
    :type obs_freqs: iterable, float

    :param v0: The :math:`V_{LSR}` velocity of the emitting object to be modelled in km/s
    :type v0: float

    :param radex_params: A dictionary containing the inputs for the RADEX model. See :func:`radex.get_default_parameters` for a list of possible parameters. Note this includes the linewidth in km/s that will be used to set the shape of the emission lines.
    :type radex_params: dict

    :param tau_profile: A function with the same arguments as :func:`maxwellian_distribution` that returns the optical depth as a function of velocity. If not set, spectralradex will assume gaussian line profiles centred on :param:`v0` and a FWHM taken from the RADEX parameters.
    :type tau_profile: function, optional
    """

    #make sure input frquencies are a sorted array
    obs_freqs=np.asarray(obs_freqs)
    obs_freqs.sort()
    
    #solve the radex model and get all line properties
    tau_0_df=get_radex_taus(radex_params)
    tau_0_df=tau_0_df[tau_0_df["freq"]<1.1*obs_freqs.max()]
    delta_v=radex_params["linewidth"]

    #obtain tau as a function of frequency and tau*radiation temperature for each line
    vfunc=lambda x: get_tau_distribution(x,v0,delta_v,obs_freqs,tau_profile)
    solve=np.apply_along_axis(vfunc,axis=1,arr=tau_0_df.values)

    #sum to get total optical depth
    taus=[x[0] for x in solve]
    taus=np.sum(taus,axis=0)  
    
    #sum radiation temperatures that have been multiplied by line taus and divide by total tau
    #this is tau weighted radiation temperature - helps for overlapping lines
    rad_weights=[x[1] for x in solve]
    rad_weights=np.sum(rad_weights,axis=0)/taus

    #finally, calculated observed brightness temperature
    taus=(rad_weights-rad_temp(2.73,obs_freqs))*(1.0-np.exp(-taus))
    #we'll return a dataframe of Frequency, Intensity
    new_df=DataFrame({"Frequency":obs_freqs,"Intensity":taus}).fillna(0.0)
    return new_df
    
def get_tau_distribution(x,v0,delta_v,frequencies,tau_profile):
    """
    Internal function meant to turn radex output into frequency
    dependent optical depth and radiation temperature
    
    """
    #unpack tau_df row
    line_freq,tau_0,t_ex=x

    #get the relative velocity of all the emitting frequencies
    velocities=((line_freq/frequencies)-1.0)*light_speed
    
    #warn user if they're frequencies are too far apart
    if velocities[1]-velocities[0]>delta_v:
        print("Velocity bins larger than linewidth")
    
    #calculate optical depth as function of v
    taus=tau_profile(v0,delta_v,tau_0,velocities)
    #calculate tau weighted radiation temperature
    rad_weights=rad_temp(t_ex,frequencies*(1.0+v0/light_speed))*taus
    return taus,rad_weights

def rad_temp(t,frequencies):
    """
    Calculate radiation temperature for a given excitation temperature
    :meta private:

    """
    hvk=(planck*frequencies*1.0e9)/boltzman
    hvk=hvk/(np.exp(hvk/t)-1)
    return hvk

def chi_squared(obs_freqs,obs_intensity,noise,v0,params):
    intensity=model_spectrum(obs_freqs,v0,params)

    indx=np.where(intensity>0.0)[0]
    error=obs_intensity[indx]-intensity[indx]
    error=(error*error)/(noise*noise)
    chi=np.sum(error)/len(obs_intensity[indx])
    return chi
    
def noise_from_spectrum(intensities):
    """
    Estimate the rms noise level from a spectrum by assuming it is a gaussian noise distribution plus positive signal.
        If this is true, the median should be the peak of the noise distribution and values smaller are just noise. Thus, the mean
        square difference between the median and smaller values is the square of the noise rms.
    :param intensities: An array of the intensity values representing a spectrm
    :type intensities: float, iterable

    :return: The rms noise value
    :rtype: float
    """
    noise=np.median(intensities)
    ints=intensities[np.where(intensities<noise)[0]]
    noise=np.mean((noise-ints)**2.0)
    noise=np.sqrt(noise)
    return noise