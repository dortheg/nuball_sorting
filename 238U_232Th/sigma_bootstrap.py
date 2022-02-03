import numpy as np
import ROOT 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time


####################################################
##     				Read in data                  ## 
####################################################

file_238U_lowE = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_2ns_lowE_12jan2022_03feb22.root"," READ ")

################   238U lowE -  134Te   #################

#Define lower and upper fit limit
x_lower = 330
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2


#Doublegate true
hist_doublegate_238U_lowE_134Te = file_238U_lowE.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_238U_lowE_134Te.GetNbinsX()

x_doublegate_238U_lowE_134Te = np.zeros(x_bins)
y_doublegate_238U_lowE_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_lowE_134Te[i] = hist_doublegate_238U_lowE_134Te.GetBinCenter(i+1)
    y_doublegate_238U_lowE_134Te[i] = hist_doublegate_238U_lowE_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_lowE_134Te = 2*x_doublegate_238U_lowE_134Te

x_doublegate_238U_lowE_134Te_long = x_doublegate_238U_lowE_134Te
y_doublegate_238U_lowE_134Te_long = y_doublegate_238U_lowE_134Te

x_doublegate_238U_lowE_134Te = x_doublegate_238U_lowE_134Te[bin_lower:bin_upper]
y_doublegate_238U_lowE_134Te = y_doublegate_238U_lowE_134Te[bin_lower:bin_upper]


#doublegate_all
hist_doublegate_all_238U_lowE_134Te = file_238U_lowE.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_238U_lowE_134Te.GetNbinsX()

x_doublegate_all_238U_lowE_134Te = np.zeros(x_bins)
y_doublegate_all_238U_lowE_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_238U_lowE_134Te[i] = hist_doublegate_all_238U_lowE_134Te.GetBinCenter(i+1)
    y_doublegate_all_238U_lowE_134Te[i] = hist_doublegate_all_238U_lowE_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_238U_lowE_134Te = 2*x_doublegate_all_238U_lowE_134Te

x_doublegate_all_238U_lowE_134Te_long = x_doublegate_all_238U_lowE_134Te
y_doublegate_all_238U_lowE_134Te_long = y_doublegate_all_238U_lowE_134Te

x_doublegate_all_238U_lowE_134Te = x_doublegate_all_238U_lowE_134Te[bin_lower:bin_upper]
y_doublegate_all_238U_lowE_134Te = y_doublegate_all_238U_lowE_134Te[bin_lower:bin_upper]


#doublegate_bg_ridge
hist_doublegate_bg_238U_lowE_134Te = file_238U_lowE.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_238U_lowE_134Te.GetNbinsX()

x_doublegate_bg_238U_lowE_134Te = np.zeros(x_bins)
y_doublegate_bg_238U_lowE_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_238U_lowE_134Te[i] = hist_doublegate_bg_238U_lowE_134Te.GetBinCenter(i+1)
    y_doublegate_bg_238U_lowE_134Te[i] = hist_doublegate_bg_238U_lowE_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_238U_lowE_134Te = 2*x_doublegate_bg_238U_lowE_134Te

x_doublegate_bg_238U_lowE_134Te_long = x_doublegate_bg_238U_lowE_134Te
y_doublegate_bg_238U_lowE_134Te_long = y_doublegate_bg_238U_lowE_134Te

x_doublegate_bg_238U_lowE_134Te = x_doublegate_bg_238U_lowE_134Te[bin_lower:bin_upper]
y_doublegate_bg_238U_lowE_134Te = y_doublegate_bg_238U_lowE_134Te[bin_lower:bin_upper]


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
## 		        Fit data first time		          ## 
####################################################


################   238U lowE -  134Te   #################

absolute_sigma_value = True

print("\n * Absolute sigma is %.d \n" % absolute_sigma_value)

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

P_double_238U_lowE_134Te, cov_double_238U_lowE_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_134Te, y_doublegate_238U_lowE_134Te, sigma=sigma_data_doublegate(y_doublegate_238U_lowE_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = absolute_sigma_value)
#print("* 238U lowE - 134Te Using uncertainty-weighted fit")
# P_double_238U_lowE_134Te, cov_double_238U_lowE_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_134Te, y_doublegate_238U_lowE_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_unc_238U_lowE_134Te = np.sqrt(np.diag(cov_double_238U_lowE_134Te))

# print("\n")
# print(" ***** 238U lowE - 134Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_134Te[0], P_double_unc_238U_lowE_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_134Te[1], P_double_unc_238U_lowE_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_134Te[2], P_double_unc_238U_lowE_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_134Te[3], P_double_unc_238U_lowE_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_134Te[4], P_double_unc_238U_lowE_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

####################################################
## 		     		 Find IYR 	                  ## 
####################################################


################   238U lowE -  134Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_134Te = np.linspace(0,40*round(tau_134Te),5000)

area_double_true_238U_lowE_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)
area_double_true_prompt_238U_lowE_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)
area_double_true_delayed_238U_lowE_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)

IYR_double_238U_lowE_134Te = IYR(prompt=area_double_true_prompt_238U_lowE_134Te, delayed=area_double_true_delayed_238U_lowE_134Te)


####################################################
###       Calculate uncertainty on IYR by MC     ###
####################################################

N = 500

#238U lowE - 134Te
P_double_new_238U_lowE_134Te = np.zeros(len(P_double_238U_lowE_134Te))
IYR_array_238U_lowE_134Te = np.zeros(N)


#Do N iterations
for n in range(N):

    for i in range(len(P_double_238U_lowE_134Te)):

    	#238U lowE - 134Te
        P_double_new_238U_lowE_134Te[i] = P_double_238U_lowE_134Te[i] + np.random.normal(0, P_double_unc_238U_lowE_134Te[i])


    #238U lowE - 134Te
    area_double_true_new_238U_lowE_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_238U_lowE_134Te[0], P_double_new_238U_lowE_134Te[1], P_double_new_238U_lowE_134Te[2], P_double_new_238U_lowE_134Te[3], P_double_new_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_new_238U_lowE_134Te[0], P_double_new_238U_lowE_134Te[1], P_double_new_238U_lowE_134Te[2], P_double_new_238U_lowE_134Te[3], P_double_new_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)
    area_double_true_prompt_new_238U_lowE_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_238U_lowE_134Te[0], P_double_new_238U_lowE_134Te[1], P_double_new_238U_lowE_134Te[2], P_double_new_238U_lowE_134Te[3], P_double_new_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)
    area_double_true_delayed_new_238U_lowE_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_new_238U_lowE_134Te[0], P_double_new_238U_lowE_134Te[1], P_double_new_238U_lowE_134Te[2], P_double_new_238U_lowE_134Te[3], P_double_new_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), x_arr_134Te)
    IYR_array_238U_lowE_134Te[n] = IYR(prompt=area_double_true_prompt_new_238U_lowE_134Te, delayed=area_double_true_delayed_new_238U_lowE_134Te)

#238U lowE - 134Te
sigma_IYR_238U_lowE_134Te = np.std(IYR_array_238U_lowE_134Te) #(np.max(IYR_array_238U_lowE_134Te) - np.min(IYR_array_238U_lowE_134Te))/2.0


####################################################
## 						Plot 		              ## 
####################################################

x_array_plot = np.linspace(0,1000,10000)


################   238U lowE -  134Te   #################

# plt.plot(x_doublegate_238U_lowE_134Te_long, y_doublegate_238U_lowE_134Te_long, label="doublegate_238U_lowE_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_134Te[0], P_double_238U_lowE_134Te[1], P_double_238U_lowE_134Te[2], P_double_238U_lowE_134Te[3], P_double_238U_lowE_134Te[4], P_double_238U_lowE_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_lowE_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_lowE_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE - 134Te: Doublegate true spectrum fit")
# #plt.axis([0,700,1,10**(4)])
# plt.axis([0,700,0,4*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()



####################################################
## 				  BOOTSTRAPPING 		          ## 
####################################################


#First calculate uncertainty in each data point
def func_sigma_data(data_all, data_bg, BG):
    return np.sqrt(data_all + data_bg + (BG*data_bg)**2)

BG_percent = 0.05
sigma_data = func_sigma_data(x_doublegate_all_238U_lowE_134Te, x_doublegate_bg_238U_lowE_134Te, BG_percent)
print(" \n * Uncertainty on BG is %.2f \n" % BG_percent)

N = 500
print(" \n * Number of resamples is %.d \n" % N)

IYR_resampled = np.zeros(N)
data = y_doublegate_238U_lowE_134Te

for n in range(N):
	data_resampled = np.zeros(len(x_doublegate_238U_lowE_134Te))

	#Resample the data
	for i in range(len(data_resampled)):
		data_resampled[i] = data[i] + np.random.normal(0, sigma_data[i])

	#Fit the resampled data
	P_resample, cov_resample = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_134Te, data_resampled, sigma=sigma_data, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = absolute_sigma_value)

	#Calculate IYR of resampled data
	area_resampled_true = np.trapz(gauss(x_arr_134Te, P_resample[0], P_resample[1], P_resample[2], P_resample[3], P_resample[4], P_resample[5]) + smeared_exp_decay(x_arr_134Te, P_resample[0], P_resample[1], P_resample[2], P_resample[3], P_resample[4], P_resample[5]), x_arr_134Te)
	area_resampled_prompt = np.trapz(gauss(x_arr_134Te, P_resample[0], P_resample[1], P_resample[2], P_resample[3], P_resample[4], P_resample[5]), x_arr_134Te)
	area_resampled_delayed = np.trapz(smeared_exp_decay(x_arr_134Te, P_resample[0], P_resample[1], P_resample[2], P_resample[3], P_resample[4], P_resample[5]), x_arr_134Te)

	IYR_resampled[n] = IYR(prompt=area_resampled_prompt, delayed=area_resampled_delayed)


sigma_bootstrap = np.std(IYR_resampled)


# plt.errorbar(x_doublegate_238U_lowE_134Te, y_doublegate_238U_lowE_134Te, yerr=sigma_data, label="doublegate_238U_lowE_134Te", color="royalblue")
# plt.plot(x_doublegate_238U_lowE_134Te, data_resampled, label="resampled", color="yellowgreen")
# plt.axis([300,700,0,4*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()



####################################################
###        		Print table of IYRs              ###
####################################################

from prettytable import PrettyTable
t = PrettyTable(['System', 'Nucleus', 'IYR', 'sigma_curvefit', 'sigma_bootstrap'])
t.add_row(['238U - lowE', '134Te', round(IYR_double_238U_lowE_134Te,3), round(sigma_IYR_238U_lowE_134Te,3), round(sigma_bootstrap,3)])
print(t)


print("\n")


