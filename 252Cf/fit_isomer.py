import numpy as np
import ROOT 
from mama_to_python import read_mama_1D, read_mama_2D
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time


file = ROOT.TFile.Open("252Cf_22nov2021.root"," READ ")
h = file.Get('time_isomer_doublegate_134Te')

x_bins = h.GetNbinsX()

x_array = np.zeros(x_bins)
y_array = np.zeros(x_bins)
    
for i in range(x_bins):
    x_array[i] = h.GetBinCenter(i+1)
    y_array[i] = h.GetBinContent(i+1)

#M_gate_double = y_array
#x_gate_double = x_array

M_gate_double, C_gate_double, x_gate_double = read_mama_1D("252Cf_time_isomer_gate_double_134Te_18may2021.m")
#M_spec, C_spec, x_spec = read_mama_1D("252Cf_time_isomer_gate_double_20may2021.m")


#For 134Te
tau = 164.1/np.log(2)

def expgaussian_plus_gauss_plus_smeared_exp(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
    """
    exponentially modified Gaussian (prompt production and isomeric decay) 
    + Gaussian component (prompt production and prompt decay)
    + Smeared exp decay from second fragment
    """
    dx = mean-x
    
    return amplitude_conv*np.exp(dx/tau)*erfc(dx/(np.sqrt(2)*sigma))+amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) + gaussian_filter1d(np.piecewise(x, [x < mean + 6, x >= mean + 6], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)]),sigma)


def expgaussian_plus_gauss_plus_exp(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
    """
    exponentially modified Gaussian (prompt production and isomeric decay) 
    + Gaussian component (prompt production and prompt decay) 
    + Exp decay from second fragment
    """
    dx = mean-x
    
    return amplitude_conv*np.exp(dx/tau)*erfc(dx/(np.sqrt(2)*sigma))+amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) + np.piecewise(x, [x < mean + 6, x >= mean + 6], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)])

def expgaussian_plus_gauss(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
    """exponentially modified Gaussian (prompt production and isomeric decay) 
    + Gaussian component (prompt production and prompt decay)
    """
    dx = mean-x
    return amplitude_conv*np.exp(dx/tau)*erfc(dx/(np.sqrt(2)*sigma))+amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))

def expgaussian(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
    """exponentially modified Gaussian (prompt production and isomeric decay)"""
    dx = mean-x
    return amplitude_conv*np.exp(dx/tau)*erfc(dx/(np.sqrt(2)*sigma))

def gauss(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    dx = mean-x
    return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))

def exp(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
	""" Exponential decay """
	dx = mean-x
	return np.piecewise(x, [x < mean + 6, x >= mean + 6], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)])

def smeared_exp(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
	""" Exponential decay """
	dx = mean-x
	return gaussian_filter1d(np.piecewise(x, [x < mean + 6, x >= mean + 6], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)]),sigma)

def sum_smeared_exp_gauss(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
	""" Sum of exponential decay and gaussian """
	dx = mean-x
	return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)]),sigma)

def sum_two_smeared_exp_gauss(x, amplitude_conv=1, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp=1.0, amplitude_exp2=1.0):
	""" Sum of exponential decay and gaussian """
	dx = mean-x
	return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)]),sigma) + gaussian_filter1d(np.piecewise(x, [x < mean+6, x >= mean+6], [lambda x:0, lambda x:amplitude_exp2*np.exp((mean-x)/tau)]),sigma)


func = sum_smeared_exp_gauss
print(str(func))


##########################
## 		Fit curve 		## 
##########################

#Remember: will need some help to find correct fit parameters
#amplitude_conv, mean, sigma, amplitude_gauss, amplitude_exp, amplitude_exp2
P, cov = curve_fit(func, x_gate_double, M_gate_double, bounds=([0,970,5,20,10,0],[100,1100,20,100,50,50]))

print("\n")
print("amplitude_conv: %.4f" % P[0])
print("mean: %.4f" % P[1])
print("sigma: %.4f" % P[2])
print("amplitude_gauss: %.4f" % P[3])
print("amplitude_exp: %.4f" % P[4])
print("amplitude_exp2: %.4f" % P[5])

##########################
## 		Chi squared 	## 
##########################

"""
Remember: higher p-value, better fit; 
H0: the two distr are the same
Ha: the two distr are significantly different 
degrees of freedom = (#rows-1)*(#colums-1) = (len(M_gate_double[1990:2300])-1)*1
"""

#chi2, p_value = chisquare(f_obs=M_gate_double[1990:2300], f_exp=func(x_gate_double, P[0], P[1], P[2], P[3], P[4], P[5])[1990:2300])
#print("\n")
#print("Chi^2 from %.d ns to %.d ns is %.3f, p-value is %.6f, df is %.d" % (x_gate_double[1990], x_gate_double[2300], chi2, p_value, len(M_gate_double[1990:2300])-1 ))

##########################
## 		Find IYR 		## 
##########################

x_arr = np.linspace(0,3000,3000)

area_tot = np.trapz(func(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), x_arr)
area_exp = np.trapz(smeared_exp(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), x_arr)
IYR = (area_exp/2.0)/area_tot
print("\n")
print("total area of fit: %.6f" % area_tot)
print("exp area of fit: %.6f" % area_exp)
print("IYR: %.6f " % IYR)
print("\n")

##########################
## 			Plot 		## 
##########################

plt.plot(x_gate_double, M_gate_double, label="time_isomer_gate_double")
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.title(str(func))
plt.plot(x_arr, func(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="fit, total")
#plt.plot(x_arr, expgaussian(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="convolution")
plt.plot(x_arr, gauss(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="gaussian")
plt.plot(x_arr, smeared_exp(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="smeared exp from in-flight fragment")
#plt.plot(x_arr, exp(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="exp from in-flight fragment")
plt.axis([800,2000,-10,150])
plt.legend(fontsize=14)
plt.grid()
plt.show()



