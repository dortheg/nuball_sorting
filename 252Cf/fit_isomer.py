import numpy as np
import ROOT 
from mama_to_python import read_mama_1D, read_mama_2D
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time


##########################
##     Read in data     ## 
##########################

file = ROOT.TFile.Open("252Cf_23nov2021.root"," READ ")

#####  134Te #####

#Doublegate
hist_doublegate_134Te = file.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_134Te.GetNbinsX()

x_doublegate_134Te = np.zeros(x_bins)
y_doublegate_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_134Te[i] = hist_doublegate_134Te.GetBinCenter(i+1)
    y_doublegate_134Te[i] = hist_doublegate_134Te.GetBinContent(i+1)

#Doublegate all
hist_doublegate_all_134Te = file.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_134Te.GetNbinsX()

x_doublegate_all_134Te = np.zeros(x_bins)
y_doublegate_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinCenter(i+1)
    y_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinContent(i+1)

#Doublegate bg
hist_doublegate_bg_134Te = file.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_134Te.GetNbinsX()

x_doublegate_bg_134Te = np.zeros(x_bins)
y_doublegate_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinCenter(i+1)
    y_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinContent(i+1)

##############

M_gate_double, C_gate_double, x_gate_double = read_mama_1D("252Cf_time_isomer_gate_double_134Te_18may2021.m")


###################################
##      Define fitting func      ## 
###################################

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
#P, cov = curve_fit(func, x_gate_double, M_gate_double, bounds=([0,970,5,20,10,0],[100,1100,20,100,50,50]))

P, cov = curve_fit(func, x_doublegate_134Te, y_doublegate_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50]))
print("\n")
print(" *****  True spectrum fit ***** \n")
print("amplitude_conv: %.4f" % P[0])
print("mean: %.4f" % P[1])
print("sigma: %.4f" % P[2])
print("amplitude_gauss: %.4f" % P[3])
print("amplitude_exp: %.4f" % P[4])
print("amplitude_exp2: %.4f" % P[5])

# P_all, cov_all = curve_fit(func, x_doublegate_all_134Te, y_doublegate_all_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50]))
# print("\n")
# print(" *****  All spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_all[0])
# print("mean: %.4f" % P_all[1])
# print("sigma: %.4f" % P_all[2])
# print("amplitude_gauss: %.4f" % P_all[3])
# print("amplitude_exp: %.4f" % P_all[4])
# print("amplitude_exp2: %.4f" % P_all[5])


# P_bg, cov_bg = curve_fit(func, x_doublegate_bg_134Te, y_doublegate_bg_134Te, bounds=([0,950,0,0,0,0],[1000,1100,40,40,40,50]))
# print("\n")
# print(" *****  BG spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_bg[0])
# print("mean: %.4f" % P_bg[1])
# print("sigma: %.4f" % P_bg[2])
# print("amplitude_gauss: %.4f" % P_bg[3])
# print("amplitude_exp: %.4f" % P_bg[4])
# print("amplitude_exp2: %.4f" % P_bg[5])

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

x_arr = np.linspace(900,4000,4000)

area_tot = np.trapz(func(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), x_arr)
area_gauss = np.trapz(gauss(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), x_arr)
area_exp = np.trapz(smeared_exp(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), x_arr)

IYR = (area_exp)/(2*area_gauss + area_exp)
print("\n")
print("*************** \n")
print("total area of true fit: %.6f" % area_tot)
print("gauss area of true fit: %.6f" % area_gauss)
print("exp area of true fit: %.6f" % area_exp)
print("IYR: %.6f " % IYR)
print("\n")

diff = area_tot - area_gauss - area_exp
print("Difference between area_tot and area_gauss+area_exp: %.4f" % diff)
print("\n")

##########################
## 			Plot 		## 
##########################

plt.plot(x_doublegate_134Te, y_doublegate_134Te, label="doublegate_134Te", color="blue")
#plt.plot(x_doublegate_all_134Te, y_doublegate_all_134Te, label="doublegate_all_134Te", color="black")
#plt.plot(x_doublegate_bg_134Te, y_doublegate_bg_134Te, label="doublegate_bg_134Te", color="pink")

plt.plot(x_arr, func(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp(x_arr, P[0], P[1], P[2], P[3], P[4], P[5]), label="true smeared exp from in-flight fragment", color="red")

# plt.plot(x_arr, func(x_arr, P_all[0], P_all[1], P_all[2], P_all[3], P_all[4], P_all[5]), label="all fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_all[0], P_all[1], P_all[2], P_all[3], P_all[4], P_all[5]), label="all gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_all[0], P_all[1], P_all[2], P_all[3], P_all[4], P_all[5]), label="all smeared exp from in-flight fragment", color="red")

# plt.plot(x_arr, func(x_arr, P_bg[0], P_bg[1], P_bg[2], P_bg[3], P_bg[4], P_bg[5]), label="bg fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_bg[0], P_bg[1], P_bg[2], P_bg[3], P_bg[4], P_bg[5]), label="bg gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_bg[0], P_bg[1], P_bg[2], P_bg[3], P_bg[4], P_bg[5]), label="bg smeared exp from in-flight fragment", color="red")


plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.title(str(func))
plt.axis([800,2000,-10,250])
plt.legend(fontsize=14)
plt.grid()
plt.show()



