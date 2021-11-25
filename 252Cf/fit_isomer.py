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

file = ROOT.TFile.Open("252Cf_24nov2021.root"," READ ")

#####################
###     134Te      ##
#####################

#Isomer_1 gate (297keV) true
hist_isomer_1_gate_134Te = file.Get('time_isomer_1_gate_134Te')
x_bins = hist_isomer_1_gate_134Te.GetNbinsX()

x_isomer_1_gate_134Te = np.zeros(x_bins)
y_isomer_1_gate_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_134Te[i] = hist_isomer_1_gate_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_134Te[i] = hist_isomer_1_gate_134Te.GetBinContent(i+1)

#Isomer_1 gate (297keV) all
hist_isomer_1_gate_all_134Te = file.Get('time_isomer_1_gate_all_134Te')
x_bins = hist_isomer_1_gate_all_134Te.GetNbinsX()

x_isomer_1_gate_all_134Te = np.zeros(x_bins)
y_isomer_1_gate_all_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_all_134Te[i] = hist_isomer_1_gate_all_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_all_134Te[i] = hist_isomer_1_gate_all_134Te.GetBinContent(i+1)

#Isomer_1 gate (297keV) bg
hist_isomer_1_gate_bg_134Te = file.Get('time_isomer_1_gate_bg_134Te')
x_bins = hist_isomer_1_gate_bg_134Te.GetNbinsX()

x_isomer_1_gate_bg_134Te = np.zeros(x_bins)
y_isomer_1_gate_bg_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_bg_134Te[i] = hist_isomer_1_gate_bg_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_bg_134Te[i] = hist_isomer_1_gate_bg_134Te.GetBinContent(i+1)


#######################

#Doublegate true
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

############################


###################################
##      Define fitting func      ## 
###################################

#For 134Te
tau = 164.1/np.log(2)

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
	return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp*np.exp((mean-x)/tau)]),sigma)

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
##     IYR functions    ## 
##########################

def IYR(prompt, delayed):
    return (delayed)/(2*prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed)
    return np.sqrt( ((2*prompt/(2*prompt+delayed)**2)*sigma_delayed)**2 + ((-2*delayed/(2*prompt+delayed)**2)*sigma_prompt)**2 )

##########################
## 		Fit curve 		## 
##########################

#Remember: will need some help to find correct fit parameters
#amplitude_conv, mean, sigma, amplitude_gauss, amplitude_exp, amplitude_exp2

######################
#134Te isomer_1 gate fits
######################

P_isomer_1, cov_isomer_1 = curve_fit(func, x_isomer_1_gate_134Te, y_isomer_1_gate_134Te, bounds=([0,950,0,20,10,0],[1000,1020,20,10000,10000,50]))
# print("\n")
# print(" *****  True spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_isomer_1[0])
# print("mean: %.4f" % P_isomer_1[1])
# print("sigma: %.4f" % P_isomer_1[2])
# print("amplitude_gauss: %.4f" % P_isomer_1[3])
# print("amplitude_exp: %.4f" % P_isomer_1[4])
# print("amplitude_exp2: %.4f" % P_isomer_1[5])

P_isomer_1_all, cov_isomer_1_all = curve_fit(func, x_isomer_1_gate_all_134Te, y_isomer_1_gate_all_134Te, bounds=([0,990,0,140000,0,0],[1000,1000,30,170000,18000,50]))
# print("\n")
# print(" *****  All spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_isomer_1_all[0])
# print("mean: %.4f" % P_isomer_1_all[1])
# print("sigma: %.4f" % P_isomer_1_all[2])
# print("amplitude_gauss: %.4f" % P_isomer_1_all[3])
# print("amplitude_exp: %.4f" % P_isomer_1_all[4])
# print("amplitude_exp2: %.4f" % P_isomer_1_all[5])

P_isomer_1_bg, cov_isomer_1_bg = curve_fit(func, x_isomer_1_gate_bg_134Te, y_isomer_1_gate_bg_134Te, bounds=([0,990,0,140000,0,0],[1000,1000,30,170000,18000,50]))
print("\n")
print(" *****  BG spectrum fit ***** \n")
print("amplitude_conv: %.4f" % P_isomer_1_bg[0])
print("mean: %.4f" % P_isomer_1_bg[1])
print("sigma: %.4f" % P_isomer_1_bg[2])
print("amplitude_gauss: %.4f" % P_isomer_1_bg[3])
print("amplitude_exp: %.4f" % P_isomer_1_bg[4])
print("amplitude_exp2: %.4f" % P_isomer_1_bg[5])

######################
#134Te doublegate fits
######################
P_double, cov_double = curve_fit(func, x_doublegate_134Te, y_doublegate_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50]))
# print("\n")
# print(" *****  True spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_double[0])
# print("mean: %.4f" % P_double[1])
# print("sigma: %.4f" % P_double[2])
# print("amplitude_gauss: %.4f" % P_double[3])
# print("amplitude_exp: %.4f" % P_double[4])
# print("amplitude_exp2: %.4f" % P_double[5])

P_double_all, cov_double_all = curve_fit(func, x_doublegate_all_134Te, y_doublegate_all_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50]))
# print("\n")
# print(" *****  All spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_double_all[0])
# print("mean: %.4f" % P_double_all[1])
# print("sigma: %.4f" % P_double_all[2])
# print("amplitude_gauss: %.4f" % P_double_all[3])
# print("amplitude_exp: %.4f" % P_double_all[4])
# print("amplitude_exp2: %.4f" % P_double_all[5])


P_double_bg, cov_double_bg = curve_fit(func, x_doublegate_bg_134Te, y_doublegate_bg_134Te, bounds=([0,950,0,0,0,0],[1000,1100,40,40,40,50]))
# print("\n")
# print(" *****  BG spectrum fit ***** \n")
# print("amplitude_conv: %.4f" % P_double_bg[0])
# print("mean: %.4f" % P_double_bg[1])
# print("sigma: %.4f" % P_double_bg[2])
# print("amplitude_gauss: %.4f" % P_double_bg[3])
# print("amplitude_exp: %.4f" % P_double_bg[4])
# print("amplitude_exp2: %.4f" % P_double_bg[5])



######################################
## 		 Find IYR + uncertainty 	## 
######################################

x_arr = np.linspace(0,3000,3000)

#Doublegated
area_double_true = np.trapz(func(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)
area_double_true_prompt = np.trapz(gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)
area_double_true_delayed = np.trapz(smeared_exp(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)

area_double_all = np.trapz(func(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)
area_double_all_prompt = np.trapz(gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)
area_double_all_delayed = np.trapz(smeared_exp(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)

area_double_bg = np.trapz(func(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)
area_double_bg_prompt = np.trapz(gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)
area_double_bg_delayed = np.trapz(smeared_exp(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)

IYR_double = IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed)
sigma_IYR_double = sigma_IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed, all_prompt=area_double_all_prompt, all_delayed=area_double_all_delayed, bg_prompt=area_double_bg_prompt, bg_delayed=area_double_bg_delayed) #Make this into a function


#Isomer_1_gate
area_isomer_1_gate_true = np.trapz(func(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), x_arr)
area_isomer_1_gate_true_prompt = np.trapz(gauss(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), x_arr)
area_isomer_1_gate_true_delayed = np.trapz(smeared_exp(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), x_arr)


IYR_isomer_1_gate = IYR(prompt=area_isomer_1_gate_true_prompt, delayed=area_isomer_1_gate_true_delayed)




print("\n")
print(" ***** Isomeric Yield Ratios ****")
print("IYR_double:        %.4f +/- %.4f " % (IYR_double, sigma_IYR_double) )
print("IYR_isomer_1 (297 keV):  %.4f" % IYR_isomer_1_gate )
print("\n")


##########################
## 			Plot 		## 
##########################

#Doublegated
#plt.errorbar(x_doublegate_134Te, y_doublegate_134Te, yerr=np.sqrt(abs(y_doublegate_134Te)),fmt=".", label="doublegate_134Te", color="royalblue")
#plt.plot(x_doublegate_134Te, y_doublegate_134Te, label="doublegate_134Te", color="royalblue")
#plt.plot(x_doublegate_all_134Te, y_doublegate_all_134Te, label="doublegate_all_134Te", color="black")
#plt.plot(x_doublegate_bg_134Te, y_doublegate_bg_134Te, label="doublegate_bg_134Te", color="pink")

# plt.plot(x_arr, func(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="all fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="all gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="all smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="bg fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="bg gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="bg smeared exp", color="red")

# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.title(str(func))
# plt.axis([800,2000,-10,250])
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


#Isomer 1 gated
#plt.plot(x_isomer_1_gate_134Te, y_isomer_1_gate_134Te, label="isomer_1_gate_134Te", color="royalblue")
#plt.plot(x_isomer_1_gate_all_134Te, y_isomer_1_gate_all_134Te, label="isomer_1_gate_all_134Te", color="black")
#plt.plot(x_isomer_1_gate_bg_134Te, y_isomer_1_gate_bg_134Te, label="isomer_1_gate_bg_134Te", color="pink")

# plt.plot(x_arr, func(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4], P_isomer_1[5]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4], P_isomer_1_all[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4], P_isomer_1_all[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4], P_isomer_1_all[5]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4], P_isomer_1_bg[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4], P_isomer_1_bg[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4], P_isomer_1_bg[5]), label="true smeared exp", color="red")

# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()

