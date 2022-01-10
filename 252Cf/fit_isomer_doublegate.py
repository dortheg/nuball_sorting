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

file = ROOT.TFile.Open("252Cf_8des2021_A.root"," READ ") #2ns bins
#ile = ROOT.TFile.Open("252Cf_8des2021.root"," READ ") #4ns bins

#######################

#Define lower and upper fit limit
x_lower = (1000 + 980)//2
x_upper = (1000 + 1400)//2


#######################

#Doublegated with 297 time increment

#Doublegate true
hist_doublegate_1_134Te = file.Get('time_isomer_doublegate_1_134Te')
x_bins = hist_doublegate_1_134Te.GetNbinsX()

x_doublegate_1_134Te = np.zeros(x_bins)
y_doublegate_1_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_134Te[i] = hist_doublegate_1_134Te.GetBinCenter(i+1)
    y_doublegate_1_134Te[i] = hist_doublegate_1_134Te.GetBinContent(i+1)


# x_doublegate_1_134Te = x_doublegate_1_134Te[x_lower:x_upper]
# y_doublegate_1_134Te = y_doublegate_1_134Te[x_lower:x_upper]

#Doublegate all
hist_doublegate_1_all_134Te = file.Get('time_isomer_doublegate_1_all_134Te')
x_bins = hist_doublegate_1_all_134Te.GetNbinsX()

x_doublegate_1_all_134Te = np.zeros(x_bins)
y_doublegate_1_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_all_134Te[i] = hist_doublegate_1_all_134Te.GetBinCenter(i+1)
    y_doublegate_1_all_134Te[i] = hist_doublegate_1_all_134Te.GetBinContent(i+1)

# x_doublegate_1_all_134Te = x_doublegate_1_all_134Te[x_lower:x_upper]
# y_doublegate_1_all_134Te = y_doublegate_1_all_134Te[x_lower:x_upper]

#Doublegate bg
hist_doublegate_1_bg_134Te = file.Get('time_isomer_doublegate_1_bg_134Te')
x_bins = hist_doublegate_1_bg_134Te.GetNbinsX()

x_doublegate_1_bg_134Te = np.zeros(x_bins)
y_doublegate_1_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_134Te[i] = hist_doublegate_1_bg_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_134Te[i] = hist_doublegate_1_bg_134Te.GetBinContent(i+1)

# x_doublegate_1_bg_134Te = x_doublegate_1_bg_134Te[x_lower:x_upper]
# y_doublegate_1_bg_134Te = y_doublegate_1_bg_134Te[x_lower:x_upper]

#Doublegate bg ridge
hist_doublegate_1_bg_ridge_134Te = file.Get('time_isomer_doublegate_1_bg_ridge_134Te')
x_bins = hist_doublegate_1_bg_ridge_134Te.GetNbinsX()

x_doublegate_1_bg_ridge_134Te = np.zeros(x_bins)
y_doublegate_1_bg_ridge_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_ridge_134Te[i] = hist_doublegate_1_bg_ridge_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_ridge_134Te[i] = hist_doublegate_1_bg_ridge_134Te.GetBinContent(i+1)

# x_doublegate_1_bg_ridge_134Te = x_doublegate_1_bg_ridge_134Te[x_lower:x_upper]
# y_doublegate_1_bg_ridge_134Te = y_doublegate_1_bg_ridge_134Te[x_lower:x_upper]

#Doublegate bg random
hist_doublegate_1_bg_random_134Te = file.Get('time_isomer_doublegate_1_bg_random_134Te')
x_bins = hist_doublegate_1_bg_random_134Te.GetNbinsX()

x_doublegate_1_bg_random_134Te = np.zeros(x_bins)
y_doublegate_1_bg_random_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_random_134Te[i] = hist_doublegate_1_bg_random_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_random_134Te[i] = hist_doublegate_1_bg_random_134Te.GetBinContent(i+1)

# x_doublegate_1_bg_random_134Te = x_doublegate_1_bg_random_134Te[x_lower:x_upper]
# y_doublegate_1_bg_random_134Te = y_doublegate_1_bg_random_134Te[x_lower:x_upper]

#1279 spectra

#Doublegate true
hist_doublegate_2_134Te = file.Get('time_isomer_doublegate_2_134Te')
x_bins = hist_doublegate_2_134Te.GetNbinsX()

x_doublegate_2_134Te = np.zeros(x_bins)
y_doublegate_2_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_134Te[i] = hist_doublegate_2_134Te.GetBinCenter(i+1)
    y_doublegate_2_134Te[i] = hist_doublegate_2_134Te.GetBinContent(i+1)

#Doublegate all
hist_doublegate_2_all_134Te = file.Get('time_isomer_doublegate_2_all_134Te')
x_bins = hist_doublegate_2_all_134Te.GetNbinsX()

x_doublegate_2_all_134Te = np.zeros(x_bins)
y_doublegate_2_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_all_134Te[i] = hist_doublegate_2_all_134Te.GetBinCenter(i+1)
    y_doublegate_2_all_134Te[i] = hist_doublegate_2_all_134Te.GetBinContent(i+1)

#Doublegate bg
hist_doublegate_2_bg_134Te = file.Get('time_isomer_doublegate_2_bg_134Te')
x_bins = hist_doublegate_2_bg_134Te.GetNbinsX()

x_doublegate_2_bg_134Te = np.zeros(x_bins)
y_doublegate_2_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_bg_134Te[i] = hist_doublegate_2_bg_134Te.GetBinCenter(i+1)
    y_doublegate_2_bg_134Te[i] = hist_doublegate_2_bg_134Te.GetBinContent(i+1)

#Doublegate bg ridge
hist_doublegate_2_bg_ridge_134Te = file.Get('time_isomer_doublegate_2_bg_ridge_134Te')
x_bins = hist_doublegate_2_bg_ridge_134Te.GetNbinsX()

x_doublegate_2_bg_ridge_134Te = np.zeros(x_bins)
y_doublegate_2_bg_ridge_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_bg_ridge_134Te[i] = hist_doublegate_2_bg_ridge_134Te.GetBinCenter(i+1)
    y_doublegate_2_bg_ridge_134Te[i] = hist_doublegate_2_bg_ridge_134Te.GetBinContent(i+1)

#Doublegate bg random
hist_doublegate_2_bg_random_134Te = file.Get('time_isomer_doublegate_2_bg_random_134Te')
x_bins = hist_doublegate_2_bg_random_134Te.GetNbinsX()

x_doublegate_2_bg_random_134Te = np.zeros(x_bins)
y_doublegate_2_bg_random_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_bg_random_134Te[i] = hist_doublegate_2_bg_random_134Te.GetBinCenter(i+1)
    y_doublegate_2_bg_random_134Te[i] = hist_doublegate_2_bg_random_134Te.GetBinContent(i+1)


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
    "Only get prompt contribution from one fragment, as the other is in-flight"
    return (delayed)/(2*prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.05*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.05*bg_delayed)**2)
    return np.sqrt( ((2*prompt/(2*prompt+delayed)**2)*sigma_delayed)**2 + ((-2*delayed/(2*prompt+delayed)**2)*sigma_prompt)**2 )

def sigma_data_singlegate(data_all, data_bg):
    return np.sqrt(data_all + data_bg + (0.05*data_bg)**2)
    #return np.sqrt(data_all + data_bg)

def sigma_singlegate_remove_0(data_all, data_bg):
    sigma = sigma_data_singlegate(data_all, data_bg)
    for i in range(len(sigma)): 
        if sigma[i] == 0:
            sigma[i] = 1
    return sigma

def sigma_data_doublegate(data_all, data_bg_ridge, data_bg_random):
    return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random + (0.03*data_bg_ridge)**2 + (0.03*data_bg_random)**2)
    #return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random)

def sigma_doublegate_remove_0(data_all, data_bg_ridge, data_bg_random):
    sigma = sigma_data_doublegate(data_all, data_bg_ridge, data_bg_random)
    for i in range(len(sigma)):
        if sigma[i] == 0:
            sigma[i] = 1
    return sigma

def sigma_fill_0(data_incoming):
    sigma = np.ones(len(data_incoming))
    for i in range(len(data_incoming)):
        if data_incoming[i] > 0:
            sigma[i] = np.sqrt(data_incoming[i])
    return sigma

###################################
##      Define fitting func      ## 
###################################

def gauss(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))


def exp_Ge(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
	""" Exponential decay """
	return np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)])


def exp_decay(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)])


def smeared_exp_Ge(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
	""" Exponential decay """
	return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)]),sigma)


def smeared_exp_decay(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma)


def sum_smeared_exp_gauss_Ge(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Sum of exponential decay and gaussian """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) 
        + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)]),sigma) )


def sum_two_smeared_exp_gauss(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Sum of exponential decay and gaussian """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)]),sigma)
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma) )

####################################################
## 		             Fit data 		              ## 
####################################################


######################
### Doublegate_1
######################

#True
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 500
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 500
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te#0
tau_decay_upper = tau_134Te+0.0001#1000

# mean_lower = 1000.00
# mean_upper = 1000.00+0.000001
# sigma_lower = 11
# sigma_upper = 11+0.000001
# amplitude_gauss_lower = 100
# amplitude_gauss_upper = 100+0.0001
# amplitude_exp_Ge_lower = 10
# amplitude_exp_Ge_upper = 10+0.000001
# tau_Ge_lower = 8.65
# tau_Ge_upper = 8.65+0.00001
# amplitude_exp_decay_lower = 65.92
# amplitude_exp_decay_upper = 65.92+0.00001
# tau_decay_lower = tau_134Te#0
# tau_decay_upper = tau_134Te+0.0001#1000

P_double_1, cov_double_1 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_1_134Te, y_doublegate_1_134Te, sigma=sigma_doublegate_remove_0(y_doublegate_1_all_134Te, y_doublegate_1_bg_ridge_134Te, y_doublegate_1_bg_random_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper])) #30nov
#P_double_1, cov_double_1 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_1_134Te, y_doublegate_1_134Te, bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper])) #30nov
print("\n")
print(" ***** 134Te:  Doublegate_1 true spectrum fit ***** ")
print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
print("mean:                     %.2f         [%.d,%.d]" % (P_double_1[0], mean_lower, mean_upper))
print("sigma:                    %.2f         [%.d,%.d]" % (P_double_1[1], sigma_lower, sigma_upper))
print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_1[2], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_1[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_1[4], tau_Ge_lower, tau_Ge_upper))
print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_1[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_1[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")

#all
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 2000
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = 0
tau_decay_upper = 1000

P_double_1_all, cov_double_1_all = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_1_all_134Te, y_doublegate_1_all_134Te, sigma=sigma_fill_0(y_doublegate_1_all_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("\n")
# print(" ***** 134Te:  Doublegate_1 all spectrum fit ***** ")
# print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_1_all[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_1_all[1], sigma_lower, sigma_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_1_all[2], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_1_all[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
# print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_1_all[4], tau_Ge_lower, tau_Ge_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_1_all[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_1_all[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

#bg
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 1000
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = 0
tau_decay_upper = 1000

P_double_1_bg, cov_double_1_bg = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_1_bg_134Te, y_doublegate_1_bg_134Te, sigma=sigma_fill_0(y_doublegate_1_bg_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("\n")
# print(" ***** 134Te:  Doublegate_1 bg spectrum fit ***** ")
# print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_1_bg[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_1_bg[1], sigma_lower, sigma_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_1_bg[2], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_1_bg[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
# print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_1_bg[4], tau_Ge_lower, tau_Ge_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_1_bg[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_1_bg[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


######################
### Doublegate_2
######################

#True
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 3000
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te#0
tau_decay_upper = tau_134Te+0.0001#1000

P_double_2, cov_double_2 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_2_134Te, y_doublegate_2_134Te, sigma=sigma_doublegate_remove_0(y_doublegate_2_all_134Te, y_doublegate_2_bg_ridge_134Te, y_doublegate_2_bg_random_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#P_double_2, cov_double_2 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_2_134Te, y_doublegate_2_134Te, bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
print("\n")
print(" ***** 134Te:  Doublegate_2 true spectrum fit ***** ")
print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
print("mean:                     %.2f         [%.d,%.d]" % (P_double_2[0], mean_lower, mean_upper))
print("sigma:                    %.2f         [%.d,%.d]" % (P_double_2[1], sigma_lower, sigma_upper))
print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_2[2], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_2[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_2[4], tau_Ge_lower, tau_Ge_upper))
print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_2[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_2[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")

#all
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 3000
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = 0
tau_decay_upper = 1000

P_double_2_all, cov_double_2_all = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_2_all_134Te, y_doublegate_2_all_134Te, sigma=sigma_fill_0(y_doublegate_2_all_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("\n")
# print(" ***** 134Te:  Doublegate_2 all spectrum fit ***** ")
# print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_2_all[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_2_all[1], sigma_lower, sigma_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_2_all[2], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_2_all[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
# print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_2_all[4], tau_Ge_lower, tau_Ge_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_2_all[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_2_all[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

#bg
mean_lower = 950
mean_upper = 1100
sigma_lower = 0
sigma_upper = 40
amplitude_gauss_lower = 0
amplitude_gauss_upper = 3000
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = 0
tau_decay_upper = 1000

P_double_2_bg, cov_double_2_bg = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_2_bg_134Te, y_doublegate_2_bg_134Te, sigma=sigma_fill_0(y_doublegate_2_bg_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("\n")
# print(" ***** 134Te:  Doublegate_2 bg spectrum fit ***** ")
# print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_2_bg[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_2_bg[1], sigma_lower, sigma_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_2_bg[2], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_2_bg[3], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
# print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_2_bg[4], tau_Ge_lower, tau_Ge_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_2_bg[5], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_2_bg[6]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


######################################
## 		 Find IYR + uncertainty 	## 
######################################

x_arr = np.linspace(-1000,3000,4000)


############
### Doublegate_1
############

area_double_1_true = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)
area_double_1_true_prompt = np.trapz(gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6])+smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)
area_double_1_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)

# print("Area double true tot: %.3f" % area_double_1_true)
# print("Area double true prompt: %.3f" % area_double_1_true_prompt)
# print("Area double true delayed: %.3f" % area_double_1_true_delayed)

area_double_1_all = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), x_arr)
area_double_1_all_prompt = np.trapz(gauss(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6])+smeared_exp_Ge(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), x_arr)
area_double_1_all_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), x_arr)

area_double_1_bg = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), x_arr)
area_double_1_bg_prompt = np.trapz(gauss(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6])+smeared_exp_Ge(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), x_arr)
area_double_1_bg_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), x_arr)

IYR_double_1 = IYR(prompt=area_double_1_true_prompt, delayed=area_double_1_true_delayed)
sigma_IYR_double_1 = sigma_IYR(prompt=area_double_1_true_prompt, delayed=area_double_1_true_delayed, all_prompt=area_double_1_all_prompt, all_delayed=area_double_1_all_delayed, bg_prompt=area_double_1_bg_prompt, bg_delayed=area_double_1_bg_delayed)

############
### Doublegate_2
############

area_double_2_true = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)
area_double_2_true_prompt = np.trapz(gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6])+smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)
area_double_2_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)

# print("Area double true tot: %.3f" % area_double_1_true)
# print("Area double true prompt: %.3f" % area_double_1_true_prompt)
# print("Area double true delayed: %.3f" % area_double_1_true_delayed)

# print("\n")

area_double_2_all = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), x_arr)
area_double_2_all_prompt = np.trapz(gauss(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6])+smeared_exp_Ge(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), x_arr)
area_double_2_all_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), x_arr)

area_double_2_bg = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), x_arr)
area_double_2_bg_prompt = np.trapz(gauss(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6])+smeared_exp_Ge(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), x_arr)
area_double_2_bg_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), x_arr)


IYR_double_2 = IYR(prompt=area_double_2_true_prompt, delayed=area_double_2_true_delayed)
sigma_IYR_double_2 = sigma_IYR(prompt=area_double_2_true_prompt, delayed=area_double_2_true_delayed, all_prompt=area_double_2_all_prompt, all_delayed=area_double_2_all_delayed, bg_prompt=area_double_2_bg_prompt, bg_delayed=area_double_2_bg_prompt)


##################################################

print(" ***** Isomeric Yield Ratios ****")
print("IYR_double_1:               %.3f +/- %.3f" % (IYR_double_1, sigma_IYR_double_1))
print("IYR_double_2:               %.3f +/- %.3f" % (IYR_double_2, sigma_IYR_double_2))
print("\n")

##########################
## 			Plot 		## 
##########################



#################
###  Doublegate_1
################

#True spectrum
plt.plot(x_doublegate_1_134Te, y_doublegate_1_134Te, label="doublegate_1_134Te", color="royalblue")
#plt.errorbar(x_doublegate_1_134Te, y_doublegate_1_134Te, yerr=sigma_doublegate_remove_0(y_doublegate_1_all_134Te, y_doublegate_1_bg_ridge_134Te, y_doublegate_1_bg_random_134Te), label="doublegate_1_134Te", color="royalblue")
plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true smeared exp Ge", color="lime")
plt.plot(x_arr, gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6])+smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="sum prompt", color="darkolivegreen")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true smeared exp decay", color="red")
plt.title("134Te: Doublegate_1 true spectrum fit")
plt.axis([800,2000,0,500])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()

# #All spectrum
# plt.errorbar(x_doublegate_1_all_134Te, y_doublegate_1_all_134Te, yerr=sigma_fill_0(y_doublegate_1_all_134Te), label="doublegate_1_all_134Te", color="black")
# plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), label="true fit, total", color="orange")
# #plt.plot(x_arr, gauss(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), label="true gaussian", color="green")
# #plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), label="true smeared exp Ge", color="lime")
# plt.plot(x_arr, gauss(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6])+smeared_exp_Ge(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), label="sum prompt", color="darkolivegreen")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_1_all[0], P_double_1_all[1], P_double_1_all[2], P_double_1_all[3], P_double_1_all[4], P_double_1_all[5], P_double_1_all[6]), label="true smeared exp decay", color="red")
# plt.title("134Te: Doublegate_1 all spectrum fit")
# plt.axis([800,2000,0,250])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

# #BG spectrum
# plt.errorbar(x_doublegate_1_bg_134Te, y_doublegate_1_bg_134Te, yerr=sigma_fill_0(y_doublegate_1_bg_134Te), label="doublegate_1_bg_134Te", color="pink")
# plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), label="true fit, total", color="orange")
# #plt.plot(x_arr, gauss(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), label="true gaussian", color="green")
# #plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), label="true smeared exp Ge", color="lime")
# plt.plot(x_arr, gauss(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6])+smeared_exp_Ge(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), label="sum prompt", color="darkolivegreen")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_1_bg[0], P_double_1_bg[1], P_double_1_bg[2], P_double_1_bg[3], P_double_1_bg[4], P_double_1_bg[5], P_double_1_bg[6]), label="true smeared exp decay", color="red")
# plt.title("134Te: Doublegate_1 BG spectrum fit")
# plt.axis([800,2000,0,250])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


#################
###  Doublegate_2
################

#true
#plt.plot(x_doublegate_2_134Te, y_doublegate_2_134Te, label="doublegate_2_134Te", color="royalblue")
plt.errorbar(x_doublegate_2_134Te, y_doublegate_2_134Te, yerr=sigma_doublegate_remove_0(y_doublegate_2_all_134Te, y_doublegate_2_bg_ridge_134Te, y_doublegate_2_bg_random_134Te), label="doublegate_2_134Te", color="aqua")
plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true fit, total", color="k")
#plt.plot(x_arr, gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true gaussian", color="green")
#plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true smeared exp Ge", color="lime")
plt.plot(x_arr, gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6])+smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="sum prompt", color="darkolivegreen")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true smeared exp decay", color="red")
plt.title("134Te: Doublegate_2 true spectrum fit")
plt.axis([800,2000,0,500])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()

#all
# plt.errorbar(x_doublegate_2_all_134Te, y_doublegate_2_all_134Te, yerr=sigma_fill_0(y_doublegate_2_all_134Te), label="doublegate_2_all_134Te", color="black")
# plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), label="true fit, total", color="orange")
# #plt.plot(x_arr, gauss(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), label="true gaussian", color="green")
# #plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), label="true smeared exp Ge", color="lime")
# plt.plot(x_arr, gauss(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6])+smeared_exp_Ge(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), label="sum prompt", color="darkolivegreen")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_2_all[0], P_double_2_all[1], P_double_2_all[2], P_double_2_all[3], P_double_2_all[4], P_double_2_all[5], P_double_2_all[6]), label="true smeared exp decay", color="red")
# plt.title("134Te: Doublegate_2 all spectrum fit")
# plt.axis([800,2000,0,500])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


# #bg
# plt.errorbar(x_doublegate_2_bg_134Te, y_doublegate_2_bg_134Te, yerr=sigma_fill_0(y_doublegate_2_bg_134Te), label="doublegate_2_bg_134Te", color="pink")
# plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), label="true fit, total", color="orange")
# #plt.plot(x_arr, gauss(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), label="true gaussian", color="green")
# #plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), label="true smeared exp Ge", color="lime")
# plt.plot(x_arr, gauss(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6])+smeared_exp_Ge(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), label="sum prompt", color="darkolivegreen")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_2_bg[0], P_double_2_bg[1], P_double_2_bg[2], P_double_2_bg[3], P_double_2_bg[4], P_double_2_bg[5], P_double_2_bg[6]), label="true smeared exp decay", color="red")
# plt.title("134Te: Doublegate_2 all spectrum fit")
# plt.axis([800,2000,0,500])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

