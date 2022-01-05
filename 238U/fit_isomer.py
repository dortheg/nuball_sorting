import numpy as np
import ROOT 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time


##########################
##     Read in data     ## 
##########################

file = ROOT.TFile.Open("CubeSort.root"," READ ")

#Define lower and upper fit limit
x_lower = 100
x_upper = 300

#######################

#Doublegate true
hist_doublegate_134Te = file.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_134Te.GetNbinsX()

x_doublegate_134Te = np.zeros(x_bins)
y_doublegate_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_134Te[i] = hist_doublegate_134Te.GetBinCenter(i+1)
    y_doublegate_134Te[i] = hist_doublegate_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_134Te = 2*x_doublegate_134Te

# x_doublegate_134Te = x_doublegate_134Te[x_lower:x_upper]
# y_doublegate_134Te = y_doublegate_134Te[x_lower:x_upper]

# #Doublegate all
# hist_doublegate_all_134Te = file.Get('time_isomer_doublegate_all_134Te')
# x_bins = hist_doublegate_all_134Te.GetNbinsX()

# x_doublegate_all_134Te = np.zeros(x_bins)
# y_doublegate_all_134Te = np.zeros(x_bins)
    
# for i in range(x_bins):
#     x_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinCenter(i+1)
#     y_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinContent(i+1)


# #Doublegate bg
# hist_doublegate_bg_134Te = file.Get('time_isomer_doublegate_bg_134Te')
# x_bins = hist_doublegate_bg_134Te.GetNbinsX()

# x_doublegate_bg_134Te = np.zeros(x_bins)
# y_doublegate_bg_134Te = np.zeros(x_bins)
    
# for i in range(x_bins):
#     x_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinCenter(i+1)
#     y_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinContent(i+1)


# #Doublegate bg ridge
# hist_doublegate_bg_ridge_134Te = file.Get('time_isomer_doublegate_bg_ridge_134Te')
# x_bins = hist_doublegate_bg_ridge_134Te.GetNbinsX()

# x_doublegate_bg_ridge_134Te = np.zeros(x_bins)
# y_doublegate_bg_ridge_134Te = np.zeros(x_bins)
    
# for i in range(x_bins):
#     x_doublegate_bg_ridge_134Te[i] = hist_doublegate_bg_ridge_134Te.GetBinCenter(i+1)
#     y_doublegate_bg_ridge_134Te[i] = hist_doublegate_bg_ridge_134Te.GetBinContent(i+1)

# #Doublegate bg random
# hist_doublegate_bg_random_134Te = file.Get('time_isomer_doublegate_bg_random_134Te')
# x_bins = hist_doublegate_bg_random_134Te.GetNbinsX()

# x_doublegate_bg_random_134Te = np.zeros(x_bins)
# y_doublegate_bg_random_134Te = np.zeros(x_bins)
    
# for i in range(x_bins):
#     x_doublegate_bg_random_134Te[i] = hist_doublegate_bg_random_134Te.GetBinCenter(i+1)
#     y_doublegate_bg_random_134Te[i] = hist_doublegate_bg_random_134Te.GetBinContent(i+1)


###################################
##            Lifetimes          ##
###################################

#For 134Te
tau_134Te = 164.1/np.log(2)
sigma_tau_134Te = 0.9/np.log(2)

##########################
##     IYR functions    ## 
##########################

def IYR(prompt, delayed):
    return (delayed)/(2*prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.05*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.05*bg_delayed)**2)
    return np.sqrt( ((2*prompt/(2*prompt+delayed)**2)*sigma_delayed)**2 + ((-2*delayed/(2*prompt+delayed)**2)*sigma_prompt)**2 )

###################################
##      Define fitting func      ## 
###################################

def gauss(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))

def smeared_exp_decay(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma)

def const_bg(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
	return const_bg + x*0

def sum_smeared_exp_gauss_const_bg(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Sum of exponential decay and gaussian, with a constant background """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) 
        + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma) 
        + const_bg)


####################################################
## 		             Fit data 		              ## 
####################################################

#True
mean_lower = 300
mean_upper = 400
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = 0#tau_134Te
tau_decay_upper = 400#tau_134Te+0.0001


P_double, cov_double = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_134Te, y_doublegate_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
print("\n")
print(" ***** 134Te:  Doublegate true spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
print("mean:                     %.2f         [%.d,%.d]" % (P_double[0], mean_lower, mean_upper))
print("sigma:                    %.2f         [%.d,%.d]" % (P_double[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f         [%.d,%.d]" % (P_double[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double[5]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")



##########################
## 			Plot 		## 
##########################

x_arr = np.linspace(0,700,1000)

#True spectrum
plt.plot(x_doublegate_134Te, y_doublegate_134Te, label="doublegate_134Te", color="royalblue")
#plt.errorbar(x_doublegate_1_134Te, y_doublegate_1_134Te, yerr=sigma_doublegate_remove_0(y_doublegate_1_all_134Te, y_doublegate_1_bg_ridge_134Te, y_doublegate_1_bg_random_134Te), label="doublegate_1_134Te", color="royalblue")
plt.plot(x_arr, sum_smeared_exp_gauss_const_bg(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true smeared exp decay", color="red")
plt.plot(x_arr, const_bg(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_134Te[x_lower],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_134Te[x_upper],0,6000, color="black")

plt.title("134Te: Doublegate true spectrum fit")
#plt.axis([800,2000,0,500])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()
