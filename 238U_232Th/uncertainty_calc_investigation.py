import numpy as np
import ROOT 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time
from prettytable import PrettyTable
import scipy.stats as stats

BOOTSTRAP = False


####################################################
##     				Read in data                  ## START_DATA
####################################################

file_238U_lowE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.bin_final.root"," READ ")
rebin_value = 1

# file_238U_lowE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.binrebin_2.root"," READ ")
# rebin_value = 2

#file_238U_lowE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.binrebin_5.root"," READ ")
#rebin_value = 5

print("\n")
print("238U lowE_highT file: ", file_238U_lowE_highT)



print("\n")
print("FIT CHOICES")

#Define lower and upper fit limit
x_lower = 330#330
x_upper = 640

bin_lower = x_lower//(2*rebin_value)
bin_upper = x_upper//(2*rebin_value)


################   238U lowE_highT -  134Te   #################
#Doublegate true
hist_doublegate_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_lowE_highT_134Te[i] = hist_doublegate_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_238U_lowE_highT_134Te[i] = hist_doublegate_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_lowE_highT_134Te = 2*x_doublegate_238U_lowE_highT_134Te

#Rebin
print(len(y_doublegate_238U_lowE_highT_134Te))

x_doublegate_238U_lowE_highT_134Te_long = x_doublegate_238U_lowE_highT_134Te
y_doublegate_238U_lowE_highT_134Te_long = y_doublegate_238U_lowE_highT_134Te

x_doublegate_238U_lowE_highT_134Te = x_doublegate_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_238U_lowE_highT_134Te = y_doublegate_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_238U_lowE_highT_134Te[0],x_doublegate_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT all -  134Te   #################

#doublegate_all
hist_doublegate_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_238U_lowE_highT_134Te[i] = hist_doublegate_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_all_238U_lowE_highT_134Te[i] = hist_doublegate_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_238U_lowE_highT_134Te = 2*x_doublegate_all_238U_lowE_highT_134Te

x_doublegate_all_238U_lowE_highT_134Te_long = x_doublegate_all_238U_lowE_highT_134Te
y_doublegate_all_238U_lowE_highT_134Te_long = y_doublegate_all_238U_lowE_highT_134Te

x_doublegate_all_238U_lowE_highT_134Te = x_doublegate_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_all_238U_lowE_highT_134Te = y_doublegate_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_all_238U_lowE_highT_134Te[0],x_doublegate_all_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT bg -  134Te   #################

#doublegate_bg
hist_doublegate_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_238U_lowE_highT_134Te[i] = hist_doublegate_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_238U_lowE_highT_134Te[i] = hist_doublegate_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_238U_lowE_highT_134Te = 2*x_doublegate_bg_238U_lowE_highT_134Te

x_doublegate_bg_238U_lowE_highT_134Te_long = x_doublegate_bg_238U_lowE_highT_134Te
y_doublegate_bg_238U_lowE_highT_134Te_long = y_doublegate_bg_238U_lowE_highT_134Te

x_doublegate_bg_238U_lowE_highT_134Te = x_doublegate_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_238U_lowE_highT_134Te = y_doublegate_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_238U_lowE_highT_134Te[0],x_doublegate_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te   #################

#doublegate_bg_ridge
hist_doublegate_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_ridge_134Te')
x_bins = hist_doublegate_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_bg_ridge_238U_lowE_highT_134Te

x_doublegate_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_bg_ridge_238U_lowE_highT_134Te
y_doublegate_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_bg_ridge_238U_lowE_highT_134Te

x_doublegate_bg_ridge_238U_lowE_highT_134Te = x_doublegate_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_ridge_238U_lowE_highT_134Te = y_doublegate_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_bg_ridge_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_random -  134Te   #################

#doublegate_bg_random
hist_doublegate_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_random_134Te')
x_bins = hist_doublegate_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_bg_random_238U_lowE_highT_134Te

x_doublegate_bg_random_238U_lowE_highT_134Te_long = x_doublegate_bg_random_238U_lowE_highT_134Te
y_doublegate_bg_random_238U_lowE_highT_134Te_long = y_doublegate_bg_random_238U_lowE_highT_134Te

x_doublegate_bg_random_238U_lowE_highT_134Te = x_doublegate_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_random_238U_lowE_highT_134Te = y_doublegate_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_random_238U_lowE_highT_134Te[0],x_doublegate_bg_random_238U_lowE_highT_134Te[-1]))



####################################################
##         			 Lifetimes                    ##
####################################################

#134Te
tau_134Te = 164.1/np.log(2)
sigma_tau_134Te = 0.9/np.log(2)




####################################################
##    			   IYR functions                  ## 
####################################################

def IYR(prompt, delayed):
	"Not same as for 252Cf, due to both fragments being stopped"
	return (delayed)/(prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.03*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.03*bg_delayed)**2)
    return np.sqrt( ((prompt/(prompt+delayed)**2)*sigma_delayed)**2 + ((delayed/(prompt+delayed)**2)*sigma_prompt)**2 )

# def sigma_data_doublegate(data_all, data_bg_ridge, data_bg_random):
#     return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random + (0.03*data_bg_ridge)**2 + (0.03*data_bg_random)**2)
#     #return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random)

def sigma_data_doublegate(data):
    unc = np.zeros(len(data))

    for i in range(len(data)):
        if data[i] <= 0:
            unc[i] = 1
        else:
            unc[i] = np.sqrt(data[i])
    return unc

def sigma_data_doublegate_all_bg(data_all, data_bg_ridge, data_bg_random):
    unc = np.zeros(len(data_all))

    for i in range(len(data_all)):
        if data_all[i] <= 0 and data_bg_ridge[i]<= 0 and data_bg_random[i]<= 0:
            unc[i] = 1 #Should it be 1?
        else:
            unc[i] = np.sqrt(data_all[i]+data_bg_ridge[i]+data_bg_random[i])
    return unc


def reduced_chisquare_func(f_obs, f_exp, N):
    dof = N-1-5

    warning = 0
    for i in range(len(f_exp)):
        if f_exp[i] < 15:
            warning = 1

    if warning==1:
        print("Warning: in chisquared-calc, the expected is < 15")

    return np.sum((f_obs-f_exp)**2/(f_exp))/dof


####################################################
##     			 Define fitting func              ## 
####################################################

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
## 		             Fit data          		      ##  START_FIT 
####################################################


################   238U lowE_highT -  134Te   #################

mean_lower = 300
mean_upper = 400
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 30000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 10000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_238U_lowE_highT_134Te, cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_highT_134Te, y_doublegate_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te Using uncertainty-weighted fit")

P_double_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_238U_lowE_highT_134Te))

print("\n")
print(" ***** 238U lowE_highT - 134Te:  Doublegate true spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[0], P_double_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[1], P_double_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[2], P_double_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[3], P_double_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[4], P_double_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")

####################################################
## 		     		 Find IYR 	                  ## 
####################################################


################   238U lowE_highT -  134Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_134Te = np.linspace(0,100*round(tau_134Te),10000)

area_double_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)

print("Area all: %.2f \n" % area_double_true_238U_lowE_highT_134Te)
print("Area prompt: %.2f \n" % area_double_true_prompt_238U_lowE_highT_134Te)
print("Area delayed: %.2f \n" % area_double_true_delayed_238U_lowE_highT_134Te)


IYR_double_238U_lowE_highT_134Te = IYR(prompt=area_double_true_prompt_238U_lowE_highT_134Te, delayed=area_double_true_delayed_238U_lowE_highT_134Te)


####################################################
##                BOOTSTRAPPING                   ## START_BOOTSTRAP
####################################################

if BOOTSTRAP==True:

    N_BOOTSTRAP = 1000

    print("\n Starting BOOTSTRAPPING... Number of iterations: %.d" % N_BOOTSTRAP)

    ### Find uncertainty on data set ###

    #238U lowE_highT - 134Te
    unc_y_doublegate_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te_long)
    unc_y_doublegate_238U_lowE_highT_134Te = sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te)

    ### Arrays to store resampled data ###

    #238U lowE_highT - 134Te
    resampled_y_doublegate_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_238U_lowE_highT_134Te_long))
    resampled_IYR_array_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    for n in range(N_BOOTSTRAP):


################   238U lowE_highT -  134Te   #################

        y_doublegate_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_all_238U_lowE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-sigma spans a BG-variation of 5%
        y_doublegate_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_all_238U_lowE_highT_134Te_long-y_doublegate_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)
        #y_doublegate_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_all_238U_lowE_highT_134Te_long-y_doublegate_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.05)


        for i in range (len(y_doublegate_all_238U_lowE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_238U_lowE_highT_134Te_long[i] = y_doublegate_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
        sigma_lower = 0
        sigma_upper = 40
        const_bg_lower = 0
        const_bg_upper = 1000
        amplitude_gauss_lower = 0
        amplitude_gauss_upper = 30000
        amplitude_exp_decay_lower = 0
        amplitude_exp_decay_upper = 10000
        tau_decay_lower = tau_134Te
        tau_decay_upper = tau_134Te+0.0001

        #Define lower and upper fit limit
        #x_lower = 330 + int(np.random.normal(0,5))
        x_lower = 330 + int(np.random.normal(0,5))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        resampled_P_double_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

        resampled_area_double_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
        resampled_area_double_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
        resampled_area_double_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)

        resampled_IYR_array_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_true_delayed_238U_lowE_highT_134Te)

        if(n%100==0):
            print("Now finished iteration %.d" % n)

    ################   END OF BOOTSTRAP-SAMPLING   #################

    #238U lowE_highT - 134Te
    sigma_bootstrap_IYR_238U_lowE_highT_134Te = np.std(resampled_IYR_array_238U_lowE_highT_134Te)

    # Plot distribution of calculated IYRs
    IYR_min = 0.5
    IYR_max = 1.0 + 2*0.001
    IYR_binwith = 0.001
    N_IYR_bins = int((IYR_max-IYR_min)//IYR_binwith)

    IYR_array = np.linspace(IYR_min,IYR_max,N_IYR_bins)
    IYR_histogram = np.zeros(len(IYR_array))

    for i in range(N_BOOTSTRAP):
    	IYR_bin = int(resampled_IYR_array_238U_lowE_highT_134Te[i]//IYR_binwith - IYR_min//IYR_binwith)
    	IYR_histogram[IYR_bin] += 1

    plt.plot(IYR_array,IYR_histogram, label="Calculated IYR distribution")
    plt.axvline(x=IYR_double_238U_lowE_highT_134Te, ymin=0, ymax=max(IYR_histogram), label="IYR", color="red")
    plt.axvline(x=IYR_double_238U_lowE_highT_134Te-sigma_bootstrap_IYR_238U_lowE_highT_134Te, ymin=0, ymax=max(IYR_histogram), label="IYR-sigma", color="black")
    plt.axvline(x=IYR_double_238U_lowE_highT_134Te+sigma_bootstrap_IYR_238U_lowE_highT_134Te, ymin=0, ymax=max(IYR_histogram), label="IYR+sigma", color="black")
    plt.xlabel("Calculated IYR", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.legend(fontsize=14)
    plt.grid()
    plt.show()

    # Calculate standard deviation by hand
    std_value = np.linspace(0,0.5,10000)
    std_found = False
    std_handcalc = 0

    #print(len(np.where((resampled_IYR_array_238U_lowE_highT_134Te >= 0.75) & (resampled_IYR_array_238U_lowE_highT_134Te < 0.77))[0]))


    for i in range(len(std_value)):
    	within_std = np.where((resampled_IYR_array_238U_lowE_highT_134Te >= IYR_double_238U_lowE_highT_134Te-std_value[i]) & (resampled_IYR_array_238U_lowE_highT_134Te < IYR_double_238U_lowE_highT_134Te+std_value[i]))[0]

    	percentage_covered = len(within_std)/len(resampled_IYR_array_238U_lowE_highT_134Te)

    	if 0.675 <= percentage_covered < 0.685:
    		std_found = True
    		std_handcalc = std_value[i]
    		print("By-hand standard deviation is %.3f" % std_value[i])
    		break



   #  plt.plot(IYR_array,IYR_histogram, label="Calculated IYR distribution")
   #  #plt.axvline(x=IYR_double_238U_lowE_highT_134Te, xmin=0, xmax=N_BOOTSTRAP, label="IYR")
   #  #plt.axvline(x=IYR_double_238U_lowE_highT_134Te-sigma_bootstrap_IYR_238U_lowE_highT_134Te, xmin=0, xmax=N_BOOTSTRAP, label="IYR-sigma")
   #  #plt.axvline(x=IYR_double_238U_lowE_highT_134Te+sigma_bootstrap_IYR_238U_lowE_highT_134Te, xmin=0, xmax=N_BOOTSTRAP, label="IYR+sigma")
  	# plt.xlabel("Calculated IYR", fontsize=14)
   #  plt.ylabel("Counts", fontsize=14)
   #  plt.legend(fontsize=14)
   #  plt.grid()
   #  plt.show()


if BOOTSTRAP==False:
    #238U lowE_highT - 134Te
    sigma_bootstrap_IYR_238U_lowE_highT_134Te = 0
    std_handcalc = 0


####################################################
###     Calculate chi-squared values for fits    ###
####################################################

chisquare_double_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_238U_lowE_highT_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), N=len(y_doublegate_238U_lowE_highT_134Te))


####################################################
###        		Print table of IYRs              ###
####################################################

print("\n")
print(" **********************************************************************************")
print("                              Isomeric Yield Ratios                                ")

t = PrettyTable(['System', 'Nucleus', 'Gate', 'Energies', 'Rebin value', 'IYR', 'unc_bootstrap np.std', 'unc_bootstrap hand-calc', 'red. chi^2'])

t.add_row(['238U - lowE highT', '134Te', 'Double', '1279-297', rebin_value , round(IYR_double_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_238U_lowE_highT_134Te,3), round(std_handcalc,3) , round(chisquare_double_238U_lowE_highT_134Te,3)])
print(t)

####################################################
## 						Plot 		              ## START_PLOT
####################################################

x_array_plot = np.linspace(0,1000,10000)


################   238U lowE_highT -  134Te   #################

#plt.plot(x_doublegate_238U_lowE_highT_134Te_long, y_doublegate_238U_lowE_highT_134Te_long, label="doublegate_238U_lowE_highT_134Te", color="royalblue")

plt.errorbar(x_doublegate_238U_lowE_highT_134Te_long, y_doublegate_238U_lowE_highT_134Te_long, yerr=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te_long), label="doublegate_238U_lowE_highT_134Te", color="royalblue")


plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_238U_lowE_highT_134Te[-1],0,6000, color="black")
#plt.yscale("log")
plt.title("238U lowE_highT - 134Te: Doublegate true spectrum fit")
#plt.axis([0,700,1,10**(4)])
#plt.axis([0,650,-50,4*10**(3)])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
#plt.legend(fontsize=12)
plt.grid()
plt.show()
