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

file_232Th_highT = ROOT.TFile.Open("Sorted_files/232Thcube_hit4_highE_highT_28feb2022.bin.root"," READ ")
file_232Th_lowT = ROOT.TFile.Open("Sorted_files/232Thcube_hit4_highE_lowT_28feb2022.bin.root"," READ ")

#Define lower and upper fit limit
# x_lower = 330
# x_upper = 640

x_lower = 330
x_upper = 640


#############################################################
### 				   	Read in data		    		  ###
#############################################################


################   232Th highT -  134Te 1279-297  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_232Th_highT_134Te = file_232Th_highT.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_232Th_highT_134Te.GetNbinsX()

x_doublegate_232Th_highT_134Te = np.zeros(x_bins)
y_doublegate_232Th_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_232Th_highT_134Te[i] = hist_doublegate_232Th_highT_134Te.GetBinCenter(i+1)
    y_doublegate_232Th_highT_134Te[i] = hist_doublegate_232Th_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_232Th_highT_134Te = 2*x_doublegate_232Th_highT_134Te

x_doublegate_232Th_highT_134Te_long = x_doublegate_232Th_highT_134Te
y_doublegate_232Th_highT_134Te_long = y_doublegate_232Th_highT_134Te

x_doublegate_232Th_highT_134Te = x_doublegate_232Th_highT_134Te[bin_lower:bin_upper]
y_doublegate_232Th_highT_134Te = y_doublegate_232Th_highT_134Te[bin_lower:bin_upper]



################   232Th highT -  134Te 1279-3n  #################

#doublegate_2 true
hist_doublegate_2_232Th_highT_134Te = file_232Th_highT.Get('time_isomer_doublegate_2_134Te')
x_bins = hist_doublegate_2_232Th_highT_134Te.GetNbinsX()

x_doublegate_2_232Th_highT_134Te = np.zeros(x_bins)
y_doublegate_2_232Th_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_2_232Th_highT_134Te[i] = hist_doublegate_2_232Th_highT_134Te.GetBinCenter(i+1)
    y_doublegate_2_232Th_highT_134Te[i] = hist_doublegate_2_232Th_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_2_232Th_highT_134Te = 2*x_doublegate_2_232Th_highT_134Te

x_doublegate_2_232Th_highT_134Te_long = x_doublegate_2_232Th_highT_134Te
y_doublegate_2_232Th_highT_134Te_long = y_doublegate_2_232Th_highT_134Te

x_doublegate_2_232Th_highT_134Te = x_doublegate_2_232Th_highT_134Te[bin_lower:bin_upper]
y_doublegate_2_232Th_highT_134Te = y_doublegate_2_232Th_highT_134Te[bin_lower:bin_upper]

#print("* 232Th_highT - 134Te 2 fit range %d - %.d" % (x_doublegate_2_232Th_highT_134Te[0],x_doublegate_2_232Th_highT_134Te[-1]))



################   232Th lowT -  134Te 297-3n  #################

#doublegate_2 true
hist_doublegate_1_232Th_lowT_134Te = file_232Th_lowT.Get('time_isomer_doublegate_1_134Te')
x_bins = hist_doublegate_1_232Th_lowT_134Te.GetNbinsX()

x_doublegate_1_232Th_lowT_134Te = np.zeros(x_bins)
y_doublegate_1_232Th_lowT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1_232Th_lowT_134Te[i] = hist_doublegate_1_232Th_lowT_134Te.GetBinCenter(i+1)
    y_doublegate_1_232Th_lowT_134Te[i] = hist_doublegate_1_232Th_lowT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1_232Th_lowT_134Te = 2*x_doublegate_1_232Th_lowT_134Te

x_doublegate_1_232Th_lowT_134Te_long = x_doublegate_1_232Th_lowT_134Te
y_doublegate_1_232Th_lowT_134Te_long = y_doublegate_1_232Th_lowT_134Te

x_doublegate_1_232Th_lowT_134Te = x_doublegate_1_232Th_lowT_134Te[bin_lower:bin_upper]
y_doublegate_1_232Th_lowT_134Te = y_doublegate_1_232Th_lowT_134Te[bin_lower:bin_upper]

#print("* 232Th_lowT - 134Te 2 fit range %d - %.d" % (x_doublegate_1_232Th_lowT_134Te[0],x_doublegate_1_232Th_lowT_134Te[-1]))

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


#############################################################
### 				     Fit spectra				      ###
#############################################################



################   232Th highT -  134Te 1279-297  #################

mean_lower = 0
mean_upper = 700
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 5000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_232Th_highT_134Te, cov_double_232Th_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_232Th_highT_134Te, y_doublegate_232Th_highT_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

P_double_unc_232Th_highT_134Te = np.sqrt(np.diag(cov_double_232Th_highT_134Te))

# print("\n")
# print(" ***** 232Th_highT - 134Te 1279-297 ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_232Th_highT_134Te[0], P_double_unc_232Th_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_highT_134Te[1], P_double_unc_232Th_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_highT_134Te[2], P_double_unc_232Th_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_highT_134Te[3], P_double_unc_232Th_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_highT_134Te[4], P_double_unc_232Th_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_highT_134Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")



################   232Th highT -  134Te 1279-3n   #################

mean_lower = 0
mean_upper = 700
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 5000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_2_232Th_highT_134Te, cov_double_2_232Th_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_2_232Th_highT_134Te, y_doublegate_2_232Th_highT_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

P_double_2_unc_232Th_highT_134Te = np.sqrt(np.diag(cov_double_2_232Th_highT_134Te))

# print("\n")
# print(" ***** 232Th_highT - 134Te 1279-3n  ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_2_232Th_highT_134Te[0], P_double_2_unc_232Th_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_2_232Th_highT_134Te[1], P_double_2_unc_232Th_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_2_232Th_highT_134Te[2], P_double_2_unc_232Th_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_2_232Th_highT_134Te[3], P_double_2_unc_232Th_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_2_232Th_highT_134Te[4], P_double_2_unc_232Th_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_2_232Th_highT_134Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   232Th lowT -  134Te 297-3n   #################

mean_lower = 0
mean_upper = 700
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 5000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_1_232Th_lowT_134Te, cov_double_1_232Th_lowT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1_232Th_lowT_134Te, y_doublegate_1_232Th_lowT_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

P_double_1_unc_232Th_lowT_134Te = np.sqrt(np.diag(cov_double_1_232Th_lowT_134Te))

# print("\n")
# print(" ***** 232Th_lowT - 134Te 297-3n ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[0], P_double_1_unc_232Th_lowT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[1], P_double_1_unc_232Th_lowT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[2], P_double_1_unc_232Th_lowT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[3], P_double_1_unc_232Th_lowT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[4], P_double_1_unc_232Th_lowT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_1_232Th_lowT_134Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")



####################################################
## 		     		 Find IYR 	                  ## 
####################################################

x_arr_134Te = np.linspace(0,40*round(tau_134Te),5000)


################   232Th highT -  134Te 1279-297   #################

area_double_true_232Th_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_highT_134Te[0], P_double_232Th_highT_134Te[1], P_double_232Th_highT_134Te[2], P_double_232Th_highT_134Te[3], P_double_232Th_highT_134Te[4], P_double_232Th_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_232Th_highT_134Te[0], P_double_232Th_highT_134Te[1], P_double_232Th_highT_134Te[2], P_double_232Th_highT_134Te[3], P_double_232Th_highT_134Te[4], P_double_232Th_highT_134Te[5]), x_arr_134Te)
area_double_true_prompt_232Th_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_highT_134Te[0], P_double_232Th_highT_134Te[1], P_double_232Th_highT_134Te[2], P_double_232Th_highT_134Te[3], P_double_232Th_highT_134Te[4], P_double_232Th_highT_134Te[5]), x_arr_134Te)
area_double_true_delayed_232Th_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_232Th_highT_134Te[0], P_double_232Th_highT_134Te[1], P_double_232Th_highT_134Te[2], P_double_232Th_highT_134Te[3], P_double_232Th_highT_134Te[4], P_double_232Th_highT_134Te[5]), x_arr_134Te)

IYR_double_232Th_highT_134Te = IYR(prompt=area_double_true_prompt_232Th_highT_134Te, delayed=area_double_true_delayed_232Th_highT_134Te)

print("\n")
print("* For 232Th highT - 134Te 1279-297: *")
print("		Area total: %.d" % area_double_true_232Th_highT_134Te)
print("		Area prompt: %.d" % area_double_true_prompt_232Th_highT_134Te)
print("		Area delayed: %.d" % area_double_true_delayed_232Th_highT_134Te)
print("		IYR 232Th highT -  134Te 1279-297: %.2f" % IYR_double_232Th_highT_134Te)


################   232Th highT -  134Te 1279-3n  #################

area_double_2_true_232Th_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_2_232Th_highT_134Te[0], P_double_2_232Th_highT_134Te[1], P_double_2_232Th_highT_134Te[2], P_double_2_232Th_highT_134Te[3], P_double_2_232Th_highT_134Te[4], P_double_2_232Th_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_2_232Th_highT_134Te[0], P_double_2_232Th_highT_134Te[1], P_double_2_232Th_highT_134Te[2], P_double_2_232Th_highT_134Te[3], P_double_2_232Th_highT_134Te[4], P_double_2_232Th_highT_134Te[5]), x_arr_134Te)
area_double_2_true_prompt_232Th_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_2_232Th_highT_134Te[0], P_double_2_232Th_highT_134Te[1], P_double_2_232Th_highT_134Te[2], P_double_2_232Th_highT_134Te[3], P_double_2_232Th_highT_134Te[4], P_double_2_232Th_highT_134Te[5]), x_arr_134Te)
area_double_2_true_delayed_232Th_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_2_232Th_highT_134Te[0], P_double_2_232Th_highT_134Te[1], P_double_2_232Th_highT_134Te[2], P_double_2_232Th_highT_134Te[3], P_double_2_232Th_highT_134Te[4], P_double_2_232Th_highT_134Te[5]), x_arr_134Te)

IYR_double_2_232Th_highT_134Te = IYR(prompt=area_double_2_true_prompt_232Th_highT_134Te, delayed=area_double_2_true_delayed_232Th_highT_134Te)

print("\n")
print("* For 232Th highT - 134Te 1279-3n: *")
print("		Area total: %.d" % area_double_2_true_232Th_highT_134Te)
print("		Area prompt: %.d" % area_double_2_true_prompt_232Th_highT_134Te)
print("		Area delayed: %.d" % area_double_2_true_delayed_232Th_highT_134Te)
print("		IYR 232Th highT -  134Te 1279-3n: %.2f" % IYR_double_2_232Th_highT_134Te)


################   232Th lowT -  134Te 297-3n  #################

area_double_1_true_232Th_lowT_134Te = np.trapz(gauss(x_arr_134Te, P_double_1_232Th_lowT_134Te[0], P_double_1_232Th_lowT_134Te[1], P_double_1_232Th_lowT_134Te[2], P_double_1_232Th_lowT_134Te[3], P_double_1_232Th_lowT_134Te[4], P_double_1_232Th_lowT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_1_232Th_lowT_134Te[0], P_double_1_232Th_lowT_134Te[1], P_double_1_232Th_lowT_134Te[2], P_double_1_232Th_lowT_134Te[3], P_double_1_232Th_lowT_134Te[4], P_double_1_232Th_lowT_134Te[5]), x_arr_134Te)
area_double_1_true_prompt_232Th_lowT_134Te = np.trapz(gauss(x_arr_134Te, P_double_1_232Th_lowT_134Te[0], P_double_1_232Th_lowT_134Te[1], P_double_1_232Th_lowT_134Te[2], P_double_1_232Th_lowT_134Te[3], P_double_1_232Th_lowT_134Te[4], P_double_1_232Th_lowT_134Te[5]), x_arr_134Te)
area_double_1_true_delayed_232Th_lowT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_1_232Th_lowT_134Te[0], P_double_1_232Th_lowT_134Te[1], P_double_1_232Th_lowT_134Te[2], P_double_1_232Th_lowT_134Te[3], P_double_1_232Th_lowT_134Te[4], P_double_1_232Th_lowT_134Te[5]), x_arr_134Te)

IYR_double_1_232Th_lowT_134Te = IYR(prompt=area_double_1_true_prompt_232Th_lowT_134Te, delayed=area_double_1_true_delayed_232Th_lowT_134Te)

print("\n")
print("* For 232Th lowT - 134Te 1279-3n: *")
print("		Area total: %.d" % area_double_1_true_232Th_lowT_134Te)
print("		Area prompt: %.d" % area_double_1_true_prompt_232Th_lowT_134Te)
print("		Area delayed: %.d" % area_double_1_true_delayed_232Th_lowT_134Te)
print("		IYR 232Th lowT -  134Te 297-3n: %.2f" % IYR_double_1_232Th_lowT_134Te)


#############################################################
### 				Calculate efficiency				  ###
#############################################################

#A: Time-gateon delayed portion, 200:, no constant BG subtraction
counts_highT_1279_3n_delayed_A = np.trapz(y_doublegate_2_232Th_highT_134Te_long[200:], x_doublegate_2_232Th_highT_134Te_long[200:])
counts_lowT_297_3n_delayed_A = np.trapz(y_doublegate_1_232Th_lowT_134Te_long[200:], x_doublegate_1_232Th_lowT_134Te_long[200:])
eff_A = counts_highT_1279_3n_delayed_A/counts_lowT_297_3n_delayed_A

#B: Time-gateon delayed portion, 250:, no constant BG subtraction
counts_highT_1279_3n_delayed_B = np.trapz(y_doublegate_2_232Th_highT_134Te_long[250:], x_doublegate_2_232Th_highT_134Te_long[250:])
counts_lowT_297_3n_delayed_B = np.trapz(y_doublegate_1_232Th_lowT_134Te_long[250:], x_doublegate_1_232Th_lowT_134Te_long[250:])
eff_B = counts_highT_1279_3n_delayed_B/counts_lowT_297_3n_delayed_B



#C: Use fitted delayed spectrum
eff_C = area_double_2_true_delayed_232Th_highT_134Te/area_double_1_true_delayed_232Th_lowT_134Te


#Print table
print("\n")
t = PrettyTable(["Symbol", "Method", "Efficiency"])
t.add_row(['A','Time-gate bin 200:, no const BG-sub', round(eff_A,2)])
t.add_row(['B','Time-gate bin 250:, no const BG-sub', round(eff_B,2)])
t.add_row(['C','Fitted delayed', round(eff_C,2)])
print(t)



#############################################################
### 				Calculate 2+ feeding				  ###
#############################################################

eff_used = 0.64

print("\n")
print("Efficiency used: %.2f" % eff_used)

print("\n")
print("Integrated number of counts in highT 1279-3n: %.d "  % area_double_2_true_232Th_highT_134Te)
print("Integrated number of counts in lowT 297-3n: %.d " % area_double_1_true_232Th_lowT_134Te)

area_double_2_true_232Th_highT_134Te_corrected = area_double_2_true_232Th_highT_134Te*(1.0/eff_used)
print("After eff-correction, Integrated number of counts in highT 1279-3n: %.d " % area_double_2_true_232Th_highT_134Te_corrected)

print("\n")
counts_2plusfeed = area_double_2_true_232Th_highT_134Te_corrected-area_double_1_true_232Th_lowT_134Te
print("Counts from 2+ side feeding: %.d " % counts_2plusfeed)

#Correct the 1279-297 IYR with the 2+-feeding
IYR_double_232Th_highT_134Te_corrected = IYR(prompt=area_double_true_prompt_232Th_highT_134Te + counts_2plusfeed, delayed=area_double_true_delayed_232Th_highT_134Te)
print("2+feeding corrected IYR 232Th highT -  134Te 1279-297: %.2f" % IYR_double_232Th_highT_134Te_corrected)
print("\n")


#############################################################
### 						Plot				          ###
#############################################################


################   232Th highT -  134Te 1279-297   #################

# plt.plot(x_doublegate_232Th_highT_134Te_long, y_doublegate_232Th_highT_134Te_long)
# plt.plot(x_arr_134Te, sum_smeared_exp_gauss_const_bg(x_arr_134Te, P_double_232Th_highT_134Te[0], P_double_232Th_highT_134Te[1], P_double_232Th_highT_134Te[2], P_double_232Th_highT_134Te[3], P_double_232Th_highT_134Te[4], P_double_232Th_highT_134Te[5]))
# plt.title("232Th highT -  134Te 1279-297")
# plt.grid()
# plt.axis([0,700,0,10000])
# plt.show()


################   232Th highT -  134Te 1279-3n  #################

# plt.plot(x_doublegate_2_232Th_highT_134Te_long, y_doublegate_2_232Th_highT_134Te_long)
# plt.plot(x_arr_134Te, sum_smeared_exp_gauss_const_bg(x_arr_134Te, P_double_2_232Th_highT_134Te[0], P_double_2_232Th_highT_134Te[1], P_double_2_232Th_highT_134Te[2], P_double_2_232Th_highT_134Te[3], P_double_2_232Th_highT_134Te[4], P_double_2_232Th_highT_134Te[5]))
# plt.title("232Th highT -  134Te 1279-3n")
# plt.grid()
# plt.axis([0,700,0,10000])
# plt.show()


################   232Th lowT -  134Te 297-3n  #################

plt.plot(x_doublegate_1_232Th_lowT_134Te_long, y_doublegate_1_232Th_lowT_134Te_long)
plt.plot(x_arr_134Te, sum_smeared_exp_gauss_const_bg(x_arr_134Te, P_double_1_232Th_lowT_134Te[0], P_double_1_232Th_lowT_134Te[1], P_double_1_232Th_lowT_134Te[2], P_double_1_232Th_lowT_134Te[3], P_double_1_232Th_lowT_134Te[4], P_double_1_232Th_lowT_134Te[5]))
plt.title("232Th lowT -  134Te 1279-3n")
plt.grid()
plt.axis([0,700,0,10000])
plt.show()


