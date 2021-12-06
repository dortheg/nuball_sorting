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

file = ROOT.TFile.Open("252Cf_6des2021.root"," READ ")
#file = ROOT.TFile.Open("FineSort.root"," READ ")

#######################

#Isomer_2 gate (1279 keV) true
hist_isomer_2_gate_134Te = file.Get('time_isomer_2_gate_134Te')
x_bins = hist_isomer_2_gate_134Te.GetNbinsX()

x_isomer_2_gate_134Te = np.zeros(x_bins)
y_isomer_2_gate_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_2_gate_134Te[i] = hist_isomer_2_gate_134Te.GetBinCenter(i+1)
    y_isomer_2_gate_134Te[i] = hist_isomer_2_gate_134Te.GetBinContent(i+1)

#Isomer_2 gate (1279 keV) all
hist_isomer_2_gate_all_134Te = file.Get('time_isomer_2_gate_all_134Te')
x_bins = hist_isomer_2_gate_all_134Te.GetNbinsX()

x_isomer_2_gate_all_134Te = np.zeros(x_bins)
y_isomer_2_gate_all_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_2_gate_all_134Te[i] = hist_isomer_2_gate_all_134Te.GetBinCenter(i+1)
    y_isomer_2_gate_all_134Te[i] = hist_isomer_2_gate_all_134Te.GetBinContent(i+1)

#Isomer_2 gate (1279 keV) bg
hist_isomer_2_gate_bg_134Te = file.Get('time_isomer_2_gate_bg_134Te')
x_bins = hist_isomer_2_gate_bg_134Te.GetNbinsX()

x_isomer_2_gate_bg_134Te = np.zeros(x_bins)
y_isomer_2_gate_bg_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_2_gate_bg_134Te[i] = hist_isomer_2_gate_bg_134Te.GetBinCenter(i+1)
    y_isomer_2_gate_bg_134Te[i] = hist_isomer_2_gate_bg_134Te.GetBinContent(i+1)

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

#Doublegate all
hist_doublegate_1_all_134Te = file.Get('time_isomer_doublegate_1_all_134Te')
x_bins = hist_doublegate_1_all_134Te.GetNbinsX()

x_doublegate_1_all_134Te = np.zeros(x_bins)
y_doublegate_1_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_all_134Te[i] = hist_doublegate_1_all_134Te.GetBinCenter(i+1)
    y_doublegate_1_all_134Te[i] = hist_doublegate_1_all_134Te.GetBinContent(i+1)

#Doublegate bg
hist_doublegate_1_bg_134Te = file.Get('time_isomer_doublegate_1_bg_134Te')
x_bins = hist_doublegate_1_bg_134Te.GetNbinsX()

x_doublegate_1_bg_134Te = np.zeros(x_bins)
y_doublegate_1_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_134Te[i] = hist_doublegate_1_bg_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_134Te[i] = hist_doublegate_1_bg_134Te.GetBinContent(i+1)

#Doublegate bg ridge
hist_doublegate_1_bg_ridge_134Te = file.Get('time_isomer_doublegate_1_bg_ridge_134Te')
x_bins = hist_doublegate_1_bg_ridge_134Te.GetNbinsX()

x_doublegate_1_bg_ridge_134Te = np.zeros(x_bins)
y_doublegate_1_bg_ridge_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_ridge_134Te[i] = hist_doublegate_1_bg_ridge_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_ridge_134Te[i] = hist_doublegate_1_bg_ridge_134Te.GetBinContent(i+1)

#Doublegate bg random
hist_doublegate_1_bg_random_134Te = file.Get('time_isomer_doublegate_1_bg_random_134Te')
x_bins = hist_doublegate_1_bg_random_134Te.GetNbinsX()

x_doublegate_1_bg_random_134Te = np.zeros(x_bins)
y_doublegate_1_bg_random_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_bg_random_134Te[i] = hist_doublegate_1_bg_random_134Te.GetBinCenter(i+1)
    y_doublegate_1_bg_random_134Te[i] = hist_doublegate_1_bg_random_134Te.GetBinContent(i+1)


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
    return (delayed)/(2*prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.05*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.05*bg_delayed)**2)
    return np.sqrt( ((2*prompt/(2*prompt+delayed)**2)*sigma_delayed)**2 + ((-2*delayed/(2*prompt+delayed)**2)*sigma_prompt)**2 )

def sigma_data(data_all, data_bg):
    return np.sqrt(data_all + data_bg + (0.05*data_bg)**2)

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


def smeared_exp_Ge(x, mean=0, sigma=1.0, amplitude_gauss=1.0, rel_amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
	""" Exponential decay """
	return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:rel_amplitude_exp_Ge*amplitude_gauss*np.exp((mean-x)/tau_Ge)]),sigma)


def smeared_exp_decay(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma)


def sum_smeared_exp_gauss_Ge(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0):
    """ Sum of exponential decay and gaussian """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) 
        + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)]),sigma) )


def sum_two_smeared_exp_gauss(x, mean=0, sigma=1.0, amplitude_gauss=1.0, rel_amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Sum of exponential decay and gaussian """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:rel_amplitude_exp_Ge*amplitude_gauss*np.exp((mean-x)/tau_Ge)]),sigma)
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma) )

####################################################
## 		             Fit data 		              ## 
####################################################


######################
#134Te isomer_2 gate fits
######################
# P_isomer_2_bg_gaussexp, cov_isomer_2_bg_gaussexp = curve_fit(sum_smeared_exp_gauss_Ge, x_isomer_2_gate_bg_134Te, y_isomer_2_gate_bg_134Te, sigma=sigma_fill_0(y_isomer_2_gate_bg_134Te), bounds=([960,0,0,0,0],[1100,10,100000,100000,50]))
# P_isomer_2_bg_gaussexp_unc = np.sqrt(np.diag(cov_isomer_2_bg_gaussexp))
# print("\n")
# print("  ***** 134Te: Isomer_2 BG spectrum fit ***** \n")
# print("          -- GAUSS + SMEARED EXP FIT --   ")
# print("mean:              %.4f +/- %.4f" % (P_isomer_2_bg_gaussexp[0], P_isomer_2_bg_gaussexp_unc[0] ))
# print("sigma:             %.4f +/- %.4f" % (P_isomer_2_bg_gaussexp[1], P_isomer_2_bg_gaussexp_unc[1]))
# print("amplitude_gauss:   %.4f +/- %.4f" % (P_isomer_2_bg_gaussexp[2], P_isomer_2_bg_gaussexp_unc[2]))
# print("amplitude_exp_Ge:  %.4f +/- %.4f" % (P_isomer_2_bg_gaussexp[3], P_isomer_2_bg_gaussexp_unc[3]))
# print("tau_Ge:            %.4f +/- %.4f" % (P_isomer_2_bg_gaussexp[4], P_isomer_2_bg_gaussexp_unc[4]))

# rel_amplitude_exp_Ge_unc = np.sqrt( (1*P_isomer_2_bg_gaussexp_unc[3]/P_isomer_2_bg_gaussexp[2])**2 + (P_isomer_2_bg_gaussexp[3]*P_isomer_2_bg_gaussexp_unc[2]/(P_isomer_2_bg_gaussexp[2])**2)**2   )
#                                                                                 #mean=0, sigma=1.0, amplitude_gauss=1.0, rel_amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0
# mean_lower = P_isomer_2_bg_gaussexp[0]
# mean_upper = P_isomer_2_bg_gaussexp[0]+0.0000001
# sigma_lower = P_isomer_2_bg_gaussexp[1]
# sigma_upper = P_isomer_2_bg_gaussexp[1]+0.0000001
# amplitude_gauss_lower = 2000
# amplitude_gauss_upper = 10000
# rel_amplitude_exp_Ge_lower = P_isomer_2_bg_gaussexp[3]/P_isomer_2_bg_gaussexp[2]
# rel_amplitude_exp_Ge_upper = P_isomer_2_bg_gaussexp[3]/P_isomer_2_bg_gaussexp[2]+0.0000001
# tau_Ge_lower = P_isomer_2_bg_gaussexp[4]
# tau_Ge_upper = P_isomer_2_bg_gaussexp[4]+0.0000001
# amplitude_exp_decay_lower = 0
# amplitude_exp_decay_upper = 10000
# tau_decay_lower = tau_134Te
# tau_decay_upper = tau_134Te+0.0000001

# P_isomer_2, cov_isomer_2 = curve_fit(sum_two_smeared_exp_gauss, x_isomer_2_gate_134Te, y_isomer_2_gate_134Te, sigma=sigma_fill_0(y_isomer_2_gate_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,rel_amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,rel_amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper])) #30nov
# print("\n")
# print(" ***** 134Te: Isomer_2 true spectrum fit ***** ")
# print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
# print("mean:                   %.4f" % P_isomer_2[0])
# print("sigma:                  %.4f" % P_isomer_2[1])
# print("amplitude_gauss:        %.4f" % P_isomer_2[2])
# print("rel_amplitude_exp_Ge:   %.4f" % P_isomer_2[3])
# print("tau_Ge:                 %.4f" % P_isomer_2[4])
# print("amplitude_exp_decay:    %.4f" % P_isomer_2[5])
# print("tau_decay:              %.4f,  in half_life:  %.4f" % (P_isomer_2[6],P_isomer_2[6]*np.log(2)))
# print("\n")

######################
#134Te doublegate fits
######################
#mean, sigma, amplitude_gauss, amplitude_exp_Ge, tau_Ge, amplitude_exp_decay, tau_decay

### 297 keV incremented doublegate
P_double_1_bg_random_gaussexp, cov_double_1_bg_random_gaussexp = curve_fit(sum_smeared_exp_gauss_Ge, x_doublegate_1_bg_random_134Te, y_doublegate_1_bg_random_134Te, sigma=sigma_fill_0(y_doublegate_1_bg_random_134Te), bounds=([960,0,0,0,0],[1100,40,1000,1000,50]))
P_double_1_bg_random_gaussexp_unc = np.sqrt(np.diag(cov_double_1_bg_random_gaussexp))
print("\n")
print(" ***** 134Te: Doublegated_1 BG random spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP FIT --   ")
print("mean:            %.4f +/- %.4f" % (P_double_1_bg_random_gaussexp[0], P_double_1_bg_random_gaussexp_unc[0] ))
print("sigma:           %.4f +/- %.4f" % (P_double_1_bg_random_gaussexp[1], P_double_1_bg_random_gaussexp_unc[1]))
print("amplitude_gauss:     %.4f +/- %.4f" % (P_double_1_bg_random_gaussexp[2], P_double_1_bg_random_gaussexp_unc[2]))
print("amplitude_exp_Ge:    %.4f +/- %.4f" % (P_double_1_bg_random_gaussexp[3], P_double_1_bg_random_gaussexp_unc[3]))
print("tau_Ge:      %.4f +/- %.4f" % (P_double_1_bg_random_gaussexp[4], P_double_1_bg_random_gaussexp_unc[4]))


rel_amplitude_exp_Ge_unc = np.sqrt( (1*P_double_1_bg_random_gaussexp_unc[3]/P_double_1_bg_random_gaussexp[2])**2 + (P_double_1_bg_random_gaussexp[3]*P_double_1_bg_random_gaussexp_unc[2]/(P_double_1_bg_random_gaussexp[2])**2)**2   )
                                                                                #mean=0, sigma=1.0, amplitude_gauss=1.0, rel_amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0
mean_lower = P_double_1_bg_random_gaussexp[0]
mean_upper = P_double_1_bg_random_gaussexp[0]+0.0000001
sigma_lower = P_double_1_bg_random_gaussexp[1]
sigma_upper = P_double_1_bg_random_gaussexp[1]+0.0000001
amplitude_gauss_lower = 0
amplitude_gauss_upper = 1000
rel_amplitude_exp_Ge_lower = P_double_1_bg_random_gaussexp[3]/P_double_1_bg_random_gaussexp[2]
rel_amplitude_exp_Ge_upper = P_double_1_bg_random_gaussexp[3]/P_double_1_bg_random_gaussexp[2]+0.0000001
tau_Ge_lower = P_double_1_bg_random_gaussexp[4]
tau_Ge_upper = P_double_1_bg_random_gaussexp[4]+0.0000001
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0000001

P_double_1, cov_double_1 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_1_134Te, y_doublegate_1_134Te, sigma=sigma_fill_0(y_doublegate_1_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,rel_amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,rel_amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper])) #30nov
print("\n")
print(" ***** 134Te:  Doublegate_1 true spectrum fit ***** ")
print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
print("mean:                   %.4f" % P_double_1[0])
print("sigma:                  %.4f" % P_double_1[1])
print("amplitude_gauss:        %.4f" % P_double_1[2])
print("rel_amplitude_exp_Ge:   %.4f" % P_double_1[3])
print("tau_Ge:                 %.4f" % P_double_1[4])
print("amplitude_exp_decay:    %.4f" % P_double_1[5])
print("tau_decay:              %.4f,  in half_life:  %.4f" % (P_double_1[6],P_double_1[6]*np.log(2)))
print("\n")



### 1279 keV incremented doublegate

P_double_2_bg_random_gaussexp, cov_double_2_bg_random_gaussexp = curve_fit(sum_smeared_exp_gauss_Ge, x_doublegate_2_bg_random_134Te, y_doublegate_2_bg_random_134Te, sigma=sigma_fill_0(y_doublegate_2_bg_random_134Te), bounds=([960,0,0,0,0],[1100,40,1000,1000,50]))
P_double_2_bg_random_gaussexp_unc = np.sqrt(np.diag(cov_double_2_bg_random_gaussexp))
print("\n")
print(" ***** 134Te: Doublegated_2 BG random spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP FIT --   ")
print("mean: 			%.4f +/- %.4f" % (P_double_2_bg_random_gaussexp[0], P_double_2_bg_random_gaussexp_unc[0] ))
print("sigma:			%.4f +/- %.4f" % (P_double_2_bg_random_gaussexp[1], P_double_2_bg_random_gaussexp_unc[1]))
print("amplitude_gauss: 	%.4f +/- %.4f" % (P_double_2_bg_random_gaussexp[2], P_double_2_bg_random_gaussexp_unc[2]))
print("amplitude_exp_Ge: 	%.4f +/- %.4f" % (P_double_2_bg_random_gaussexp[3], P_double_2_bg_random_gaussexp_unc[3]))
print("tau_Ge: 		%.4f +/- %.4f" % (P_double_2_bg_random_gaussexp[4], P_double_2_bg_random_gaussexp_unc[4]))



rel_amplitude_exp_Ge_unc = np.sqrt( (1*P_double_2_bg_random_gaussexp_unc[3]/P_double_2_bg_random_gaussexp[2])**2 + (P_double_2_bg_random_gaussexp[3]*P_double_2_bg_random_gaussexp_unc[2]/(P_double_2_bg_random_gaussexp[2])**2)**2   )
                                                                                #mean=0, sigma=1.0, amplitude_gauss=1.0, rel_amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0
mean_lower = P_double_2_bg_random_gaussexp[0]
mean_upper = P_double_2_bg_random_gaussexp[0]+0.0000001
sigma_lower = P_double_2_bg_random_gaussexp[1]
sigma_upper = P_double_2_bg_random_gaussexp[1]+0.0000001
amplitude_gauss_lower = 0
amplitude_gauss_upper = 1000
rel_amplitude_exp_Ge_lower = P_double_2_bg_random_gaussexp[3]/P_double_2_bg_random_gaussexp[2]
rel_amplitude_exp_Ge_upper = P_double_2_bg_random_gaussexp[3]/P_double_2_bg_random_gaussexp[2]+0.0000001
tau_Ge_lower = P_double_2_bg_random_gaussexp[4]
tau_Ge_upper = P_double_2_bg_random_gaussexp[4]+0.0000001
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0000001

P_double_2, cov_double_2 = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_2_134Te, y_doublegate_2_134Te, sigma=sigma_fill_0(y_doublegate_2_134Te), bounds=([mean_lower,sigma_lower,amplitude_gauss_lower,rel_amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,amplitude_gauss_upper,rel_amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper])) #30nov
print("\n")
print(" ***** 134Te:  Doublegate_2 true spectrum fit ***** ")
print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
print("mean:                   %.4f" % P_double_2[0])
print("sigma:                  %.4f" % P_double_2[1])
print("amplitude_gauss:        %.4f" % P_double_2[2])
print("rel_amplitude_exp_Ge:   %.4f" % P_double_2[3])
print("tau_Ge:                 %.4f" % P_double_2[4])
print("amplitude_exp_decay:    %.4f" % P_double_2[5])
print("tau_decay:              %.4f,  in half_life:  %.4f" % (P_double_2[6],P_double_2[6]*np.log(2)))
print("\n")


######################################
## 		 Find IYR + uncertainty 	## 
######################################

x_arr = np.linspace(-1000,3000,3000)

area_double_1_true = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)
area_double_1_true_prompt = np.trapz(gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6])+smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)
area_double_1_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), x_arr)

# print("Area double true tot: %.3f" % area_double_1_true)
# print("Area double true prompt: %.3f" % area_double_1_true_prompt)
# print("Area double true delayed: %.3f" % area_double_1_true_delayed)

area_double_2_true = np.trapz(sum_two_smeared_exp_gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)
area_double_2_true_prompt = np.trapz(gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6])+smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)
area_double_2_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), x_arr)

# print("Area double true tot: %.3f" % area_double_1_true)
# print("Area double true prompt: %.3f" % area_double_1_true_prompt)
# print("Area double true delayed: %.3f" % area_double_1_true_delayed)

print("\n")
IYR_double_1 = IYR(prompt=area_double_1_true_prompt, delayed=area_double_1_true_delayed)
IYR_double_2 = IYR(prompt=area_double_2_true_prompt, delayed=area_double_2_true_delayed)

print(" ***** Isomeric Yield Ratio ****")
print("IYR_double_1:               %.4f" % IYR_double_1 )
print("IYR_double_2:               %.4f" % IYR_double_2 )
print("\n")

##########################
## 			Plot 		## 
##########################

# #1279Gate BG-fit
# plt.errorbar(x_isomer_2_gate_bg_134Te, y_isomer_2_gate_bg_134Te, yerr=sigma_data(y_isomer_2_gate_all_134Te, y_isomer_2_gate_bg_134Te), label="isomer_2_gate_bg_134Te", color="orange")
# plt.plot(x_arr, sum_smeared_exp_gauss_Ge(x_arr, P_isomer_2_bg_gaussexp[0], P_isomer_2_bg_gaussexp[1], P_isomer_2_bg_gaussexp[2], P_isomer_2_bg_gaussexp[3], P_isomer_2_bg_gaussexp[4]), label="gauss + smeared exp", color="red")
# plt.plot(x_arr, gauss(x_arr, P_isomer_2_bg_gaussexp[0], P_isomer_2_bg_gaussexp[1], P_isomer_2_bg_gaussexp[2], P_isomer_2_bg_gaussexp[3], P_isomer_2_bg_gaussexp[4]), label="gauss", color="green")
# plt.plot(x_arr, smeared_exp_Ge(x_arr, P_isomer_2_bg_gaussexp[0], P_isomer_2_bg_gaussexp[1], P_isomer_2_bg_gaussexp[2], P_isomer_2_bg_gaussexp[3]/P_isomer_2_bg_gaussexp[2], P_isomer_2_bg_gaussexp[4]), label="smeared Ge exp", color="lime")

# plt.title("134Te: Isomer_2 (1279keV) BG-fit")
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.axis([800,2000,0,30000])
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

# #1279Gate true spectrum
# plt.errorbar(x_isomer_2_gate_134Te, y_isomer_2_gate_134Te, yerr=sigma_fill_0(y_isomer_2_gate_134Te), label="isomer_2_134Te", color="royalblue")
# plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6]), label="true fit, total", color="orange")
# #plt.plot(x_arr, gauss(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6]), label="true gaussian", color="green")
# #plt.plot(x_arr, smeared_exp_Ge(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6]), label="true smeared exp Ge", color="lime")
# plt.plot(x_arr, gauss(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6])+smeared_exp_Ge(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6]), label="sum prompt", color="darkolivegreen")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4], P_isomer_2[5], P_isomer_2[6]), label="true smeared exp decay", color="red")

# plt.title("134Te: Isomer_2 (1279keV) true fit")
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.axis([800,2000,0,10000])
# plt.grid()
# plt.show()


#Doublegate

#297
#Random BG-fit
plt.errorbar(x_doublegate_1_bg_random_134Te, y_doublegate_1_bg_random_134Te, yerr=np.sqrt(y_doublegate_1_bg_random_134Te), label="doublegate_1_bg_random_134Te", color="orange")
plt.plot(x_arr, sum_smeared_exp_gauss_Ge(x_arr, P_double_1_bg_random_gaussexp[0], P_double_1_bg_random_gaussexp[1], P_double_1_bg_random_gaussexp[2], P_double_1_bg_random_gaussexp[3], P_double_1_bg_random_gaussexp[4]), label="gauss + smeared exp", color="red")
plt.plot(x_arr, gauss(x_arr, P_double_1_bg_random_gaussexp[0], P_double_1_bg_random_gaussexp[1], P_double_1_bg_random_gaussexp[2], P_double_1_bg_random_gaussexp[3], P_double_1_bg_random_gaussexp[4]), label="gauss", color="green")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_1_bg_random_gaussexp[0], P_double_1_bg_random_gaussexp[1], P_double_1_bg_random_gaussexp[2], P_double_1_bg_random_gaussexp[3]/P_double_1_bg_random_gaussexp[2], P_double_1_bg_random_gaussexp[4]), label="smeared Ge exp", color="lime")

plt.title("134Te: Doublegate_1 random BG-fit")
plt.axis([800,2000,0,250])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()

#Fit of true spectrum
#plt.plot(x_doublegate_1_134Te, y_doublegate_1_134Te, label="doublegate_1_134Te", color="royalblue")
plt.errorbar(x_doublegate_1_134Te, y_doublegate_1_134Te, yerr=sigma_fill_0(y_doublegate_1_134Te), label="doublegate_1_134Te", color="royalblue")
plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true smeared exp Ge", color="lime")
plt.plot(x_arr, gauss(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6])+smeared_exp_Ge(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="sum prompt", color="darkolivegreen")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_1[0], P_double_1[1], P_double_1[2], P_double_1[3], P_double_1[4], P_double_1[5], P_double_1[6]), label="true smeared exp decay", color="red")

plt.title("134Te: Doublegate_1 true spectrum fit")
plt.axis([800,2000,0,250])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()


#1279
#Random BG-fit
plt.errorbar(x_doublegate_2_bg_random_134Te, y_doublegate_2_bg_random_134Te, yerr=np.sqrt(y_doublegate_2_bg_random_134Te), label="doublegate_2_bg_random_134Te", color="orange")
plt.plot(x_arr, sum_smeared_exp_gauss_Ge(x_arr, P_double_2_bg_random_gaussexp[0], P_double_2_bg_random_gaussexp[1], P_double_2_bg_random_gaussexp[2], P_double_2_bg_random_gaussexp[3], P_double_2_bg_random_gaussexp[4]), label="gauss + smeared exp", color="red")
plt.plot(x_arr, gauss(x_arr, P_double_2_bg_random_gaussexp[0], P_double_2_bg_random_gaussexp[1], P_double_2_bg_random_gaussexp[2], P_double_2_bg_random_gaussexp[3], P_double_2_bg_random_gaussexp[4]), label="gauss", color="green")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2_bg_random_gaussexp[0], P_double_2_bg_random_gaussexp[1], P_double_2_bg_random_gaussexp[2], P_double_2_bg_random_gaussexp[3]/P_double_2_bg_random_gaussexp[2], P_double_2_bg_random_gaussexp[4]), label="smeared Ge exp", color="lime")

plt.title("134Te: Doublegate_2 random BG-fit")
plt.axis([800,2000,0,250])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()

#Fit of true spectrum
#plt.plot(x_doublegate_2_134Te, y_doublegate_2_134Te, label="doublegate_2_134Te", color="royalblue")
plt.errorbar(x_doublegate_2_134Te, y_doublegate_2_134Te, yerr=sigma_fill_0(y_doublegate_2_134Te), label="doublegate_2_134Te", color="royalblue")
plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true smeared exp Ge", color="lime")
plt.plot(x_arr, gauss(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6])+smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="sum prompt", color="darkolivegreen")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6]), label="true smeared exp decay", color="red")

plt.title("134Te: Doublegate_2 true spectrum fit")
plt.axis([800,2000,0,250])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()



