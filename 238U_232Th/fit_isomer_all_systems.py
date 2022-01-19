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

#file_238U = ROOT.TFile.Open("Sorted_files/238Ucube_hit3_2ns_lowE_12jan2022.root"," READ ")
file_238U = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_2ns_lowE_12jan2022.root"," READ ")
#file_238U = ROOT.TFile.Open("Sorted_files/238Ucube_hit3_2ns_lowE_12jan2022.root"," READ ")

file_232Th = ROOT.TFile.Open("Sorted_files/232Thcube_hit4_2ns_17jan2022.root"," READ ")

print("\n")
print("238U file: ", file_238U)
print("232Th file: ", file_232Th)


print("\n")
print("FIT CHOICES")



################   238U -  134Te   #################

#Define lower and upper fit limit
x_lower = 330
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_238U_134Te = file_238U.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_238U_134Te.GetNbinsX()

x_doublegate_238U_134Te = np.zeros(x_bins)
y_doublegate_238U_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_134Te[i] = hist_doublegate_238U_134Te.GetBinCenter(i+1)
    y_doublegate_238U_134Te[i] = hist_doublegate_238U_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_134Te = 2*x_doublegate_238U_134Te

x_doublegate_238U_134Te_long = x_doublegate_238U_134Te
y_doublegate_238U_134Te_long = y_doublegate_238U_134Te

x_doublegate_238U_134Te = x_doublegate_238U_134Te[bin_lower:bin_upper]
y_doublegate_238U_134Te = y_doublegate_238U_134Te[bin_lower:bin_upper]

print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_238U_134Te[0],x_doublegate_238U_134Te[-1]))

################   238U -  134Te 1n  #################

#doublegate_1n true
hist_doublegate_1n_238U_134Te = file_238U.Get('time_isomer_doublegate_1n_134Te')
x_bins = hist_doublegate_1n_238U_134Te.GetNbinsX()

x_doublegate_1n_238U_134Te = np.zeros(x_bins)
y_doublegate_1n_238U_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_238U_134Te[i] = hist_doublegate_1n_238U_134Te.GetBinCenter(i+1)
    y_doublegate_1n_238U_134Te[i] = hist_doublegate_1n_238U_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_238U_134Te = 2*x_doublegate_1n_238U_134Te

x_doublegate_1n_238U_134Te_long = x_doublegate_1n_238U_134Te
y_doublegate_1n_238U_134Te_long = y_doublegate_1n_238U_134Te

x_doublegate_1n_238U_134Te = x_doublegate_1n_238U_134Te[bin_lower:bin_upper]
y_doublegate_1n_238U_134Te = y_doublegate_1n_238U_134Te[bin_lower:bin_upper]

print("* 238U - 134Te 1n fit range %d - %.d" % (x_doublegate_1n_238U_134Te[0],x_doublegate_1n_238U_134Te[-1]))



################   238U -  134Te 3n  #################

#doublegate_3n true
hist_doublegate_3n_238U_134Te = file_238U.Get('time_isomer_doublegate_3n_134Te')
x_bins = hist_doublegate_3n_238U_134Te.GetNbinsX()

x_doublegate_3n_238U_134Te = np.zeros(x_bins)
y_doublegate_3n_238U_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_238U_134Te[i] = hist_doublegate_3n_238U_134Te.GetBinCenter(i+1)
    y_doublegate_3n_238U_134Te[i] = hist_doublegate_3n_238U_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_238U_134Te = 2*x_doublegate_3n_238U_134Te

x_doublegate_3n_238U_134Te_long = x_doublegate_3n_238U_134Te
y_doublegate_3n_238U_134Te_long = y_doublegate_3n_238U_134Te

x_doublegate_3n_238U_134Te = x_doublegate_3n_238U_134Te[bin_lower:bin_upper]
y_doublegate_3n_238U_134Te = y_doublegate_3n_238U_134Te[bin_lower:bin_upper]

print("* 238U - 134Te 3n fit range %d - %.d" % (x_doublegate_3n_238U_134Te[0],x_doublegate_3n_238U_134Te[-1]))


################   238U -  134Te 5n  #################

#doublegate_5n true
hist_doublegate_5n_238U_134Te = file_238U.Get('time_isomer_doublegate_5n_134Te')
x_bins = hist_doublegate_5n_238U_134Te.GetNbinsX()

x_doublegate_5n_238U_134Te = np.zeros(x_bins)
y_doublegate_5n_238U_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_238U_134Te[i] = hist_doublegate_5n_238U_134Te.GetBinCenter(i+1)
    y_doublegate_5n_238U_134Te[i] = hist_doublegate_5n_238U_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_238U_134Te = 2*x_doublegate_5n_238U_134Te

x_doublegate_5n_238U_134Te_long = x_doublegate_5n_238U_134Te
y_doublegate_5n_238U_134Te_long = y_doublegate_5n_238U_134Te

x_doublegate_5n_238U_134Te = x_doublegate_5n_238U_134Te[bin_lower:bin_upper]
y_doublegate_5n_238U_134Te = y_doublegate_5n_238U_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 5n fit range %d - %.d" % (x_doublegate_5n_238U_134Te[0],x_doublegate_5n_238U_134Te[-1]))



################   232Th -  134Te   #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_232Th_134Te = file_232Th.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_232Th_134Te.GetNbinsX()

x_doublegate_232Th_134Te = np.zeros(x_bins)
y_doublegate_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_232Th_134Te[i] = hist_doublegate_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_232Th_134Te[i] = hist_doublegate_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_232Th_134Te = 2*x_doublegate_232Th_134Te

x_doublegate_232Th_134Te_long = x_doublegate_232Th_134Te
y_doublegate_232Th_134Te_long = y_doublegate_232Th_134Te

x_doublegate_232Th_134Te = x_doublegate_232Th_134Te[bin_lower:bin_upper]
y_doublegate_232Th_134Te = y_doublegate_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_232Th_134Te[0],x_doublegate_232Th_134Te[-1]))


################   238U -  135Te   #################

#Define lower and upper fit limit
x_lower = 320
x_upper = 600

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_238U_135Te = file_238U.Get('time_isomer_doublegate_135Te')
x_bins = hist_doublegate_238U_135Te.GetNbinsX()

x_doublegate_238U_135Te = np.zeros(x_bins)
y_doublegate_238U_135Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_135Te[i] = hist_doublegate_238U_135Te.GetBinCenter(i+1)
    y_doublegate_238U_135Te[i] = hist_doublegate_238U_135Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_135Te = 2*x_doublegate_238U_135Te

x_doublegate_238U_135Te_long = x_doublegate_238U_135Te
y_doublegate_238U_135Te_long = y_doublegate_238U_135Te

x_doublegate_238U_135Te = x_doublegate_238U_135Te[bin_lower:bin_upper]
y_doublegate_238U_135Te = y_doublegate_238U_135Te[bin_lower:bin_upper]

#print("* 135Te fit range %d - %.d" % (x_doublegate_238U_135Te[0],x_doublegate_238U_135Te[-1]))



####################################################
##         			 Lifetimes                    ##
####################################################

#134Te
tau_134Te = 164.1/np.log(2)
sigma_tau_134Te = 0.9/np.log(2)

#135Te
tau_135Te = 511.0/np.log(2)
sigma_tau_135Te = 20.0/np.log(2)



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
## 		             Fit data 		              ## 
####################################################


################   238U -  134Te   #################

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

P_double_238U_134Te, cov_double_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_134Te, y_doublegate_238U_134Te, sigma=sigma_data_doublegate(y_doublegate_238U_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
print("\n* 238U - 134Te Using uncertainty-weighted fit")
# P_double_238U_134Te, cov_double_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_134Te, y_doublegate_238U_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_unc_238U_134Te = np.sqrt(np.diag(cov_double_238U_134Te))

print("\n")
print(" ***** 238U - 134Te:  Doublegate true spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_134Te[0], P_double_unc_238U_134Te[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_134Te[1], P_double_unc_238U_134Te[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_134Te[2], P_double_unc_238U_134Te[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_134Te[3], P_double_unc_238U_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_134Te[4], P_double_unc_238U_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")

################   238U -  134Te 1n  #################

mean_lower = 341
mean_upper = 342
sigma_lower = 4
sigma_upper = 10
const_bg_lower = 0
const_bg_upper = 10
amplitude_gauss_lower = 550
amplitude_gauss_upper = 1000
amplitude_exp_decay_lower = 40
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_1n_238U_134Te, cov_double_1n_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_238U_134Te, y_doublegate_1n_238U_134Te, sigma=sigma_data_doublegate(y_doublegate_1n_238U_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#print("\n* 238U - 134Te 1n Using uncertainty-weighted fit")
# P_double_1n_238U_134Te, cov_double_1n_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_238U_134Te, y_doublegate_1n_238U_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_1n_unc_238U_134Te = np.sqrt(np.diag(cov_double_1n_238U_134Te))

print("\n")
print(" ***** 238U - 134Te 1n:  doublegate true spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_1n_238U_134Te[0], P_double_1n_unc_238U_134Te[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_134Te[1], P_double_1n_unc_238U_134Te[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_134Te[2], P_double_1n_unc_238U_134Te[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_134Te[3], P_double_1n_unc_238U_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_134Te[4], P_double_1n_unc_238U_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")


################   238U -  134Te 3n  #################

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

P_double_3n_238U_134Te, cov_double_3n_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_238U_134Te, y_doublegate_3n_238U_134Te, sigma=sigma_data_doublegate(y_doublegate_3n_238U_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
print("\n* 238U - 134Te 3n Using uncertainty-weighted fit")
# P_double_3n_238U_134Te, cov_double_3n_238U_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_238U_134Te, y_doublegate_3n_238U_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_3n_unc_238U_134Te = np.sqrt(np.diag(cov_double_3n_238U_134Te))

print("\n")
print(" ***** 238U - 134Te 3n:  doublegate true spectrum fit ***** ")
print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_238U_134Te[0], P_double_3n_unc_238U_134Te[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_134Te[1], P_double_3n_unc_238U_134Te[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_134Te[2], P_double_3n_unc_238U_134Te[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_134Te[3], P_double_3n_unc_238U_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_134Te[4], P_double_3n_unc_238U_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")


################   232Th -  134Te   #################

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

P_double_232Th_134Te, cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_232Th_134Te, y_doublegate_232Th_134Te, sigma=sigma_data_doublegate(y_doublegate_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#print("\n * 232Th - 134Te Using uncertainty-weighted fit")
#P_double_232Th_134Te, cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_232Th_134Te, y_doublegate_232Th_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_unc_232Th_134Te = np.sqrt(np.diag(cov_double_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_232Th_134Te[0], P_double_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_134Te[1], P_double_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_134Te[2], P_double_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_134Te[3], P_double_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_134Te[4], P_double_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_232Th_134Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U -  135Te   #################

mean_lower = 0
mean_upper = 700
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 20
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 5000
tau_decay_lower = tau_135Te
tau_decay_upper = tau_135Te+0.0001

#P_double_238U_135Te, cov_double = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_135Te, y_doublegate_135Te, sigma=sigma_data_doublegate(y_doublegate_all_135Te, y_doublegate_bg_ridge_135Te, y_doublegate_bg_random_135Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#print("\n Using uncertainty-weighted fit")
P_double_238U_135Te, cov_double_238U_135Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_135Te, y_doublegate_238U_135Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#print("* Not using uncertainty-weighted fit")

P_double_unc_238U_135Te = np.sqrt(np.diag(cov_double_238U_135Te))

# print("\n")
# print(" ***** 238U - 135Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_135Te[0], P_double_unc_238U_135Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_135Te[1], P_double_unc_238U_135Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_135Te[2], P_double_unc_238U_135Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_135Te[3], P_double_unc_238U_135Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_135Te[4], P_double_unc_238U_135Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_135Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


####################################################
## 		     		 Find IYR 	                  ## 
####################################################


################   238U -  134Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_134Te = np.linspace(0,40*round(tau_134Te),5000)

area_double_true_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)
area_double_true_prompt_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)
area_double_true_delayed_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)

IYR_double_238U_134Te = IYR(prompt=area_double_true_prompt_238U_134Te, delayed=area_double_true_delayed_238U_134Te)


################   238U -  134Te 1n  #################

area_double_1n_true_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)
area_double_1n_true_prompt_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)
area_double_1n_true_delayed_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)

IYR_double_1n_238U_134Te = IYR(prompt=area_double_1n_true_prompt_238U_134Te, delayed=area_double_1n_true_delayed_238U_134Te)



################   238U -  134Te 3n  #################

area_double_3n_true_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)
area_double_3n_true_prompt_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)
area_double_3n_true_delayed_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)

IYR_double_3n_238U_134Te = IYR(prompt=area_double_3n_true_prompt_238U_134Te, delayed=area_double_3n_true_delayed_238U_134Te)



################   232Th -  134Te   #################

area_double_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
area_double_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
area_double_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)

IYR_double_232Th_134Te = IYR(prompt=area_double_true_prompt_232Th_134Te, delayed=area_double_true_delayed_232Th_134Te)




################   238U -  135Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_135Te = np.linspace(0,40*round(tau_135Te),40*round(tau_135Te))

area_double_true_238U_135Te = np.trapz(gauss(x_arr_135Te, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]) + smeared_exp_decay(x_arr_135Te, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)
area_double_true_prompt_238U_135Te = np.trapz(gauss(x_arr_135Te, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)
area_double_true_delayed_238U_135Te = np.trapz(smeared_exp_decay(x_arr_135Te, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)

IYR_double_238U_135Te = IYR(prompt=area_double_true_prompt_238U_135Te, delayed=area_double_true_delayed_238U_135Te)



####################################################
###       Calculate uncertainty on IYR by MC     ###
####################################################

N = 500

#238U - 134Te
P_double_new_238U_134Te = np.zeros(len(P_double_238U_134Te))
IYR_array_238U_134Te = np.zeros(N)

#238U - 134Te 1n
P_double_1n_new_238U_134Te = np.zeros(len(P_double_1n_238U_134Te))
IYR_array_1n_238U_134Te = np.zeros(N)

#238U - 134Te 3n
P_double_3n_new_238U_134Te = np.zeros(len(P_double_3n_238U_134Te))
IYR_array_3n_238U_134Te = np.zeros(N)

#232Th - 134Te
P_double_new_232Th_134Te = np.zeros(len(P_double_232Th_134Te))
IYR_array_232Th_134Te = np.zeros(N)

#238U - 135Te
P_double_new_238U_135Te = np.zeros(len(P_double_238U_135Te))
IYR_array_238U_135Te = np.zeros(N)

#Do N iterations
for n in range(N):

    for i in range(len(P_double_238U_134Te)):

    	#238U - 134Te
        P_double_new_238U_134Te[i] = P_double_238U_134Te[i] + P_double_unc_238U_134Te[i]*np.random.uniform(-1, 1)

        #238U - 134Te 1n
        P_double_1n_new_238U_134Te[i] = P_double_1n_238U_134Te[i] + P_double_1n_unc_238U_134Te[i]*np.random.uniform(-1, 1)

        #238U - 134Te 3n
        P_double_3n_new_238U_134Te[i] = P_double_3n_238U_134Te[i] + P_double_3n_unc_238U_134Te[i]*np.random.uniform(-1, 1)

        #232Th - 134Te
        P_double_new_232Th_134Te[i] = P_double_232Th_134Te[i] + P_double_unc_232Th_134Te[i]*np.random.uniform(-1, 1)


        #238U - 135Te
        P_double_new_238U_135Te[i] = P_double_238U_135Te[i] + P_double_unc_238U_135Te[i]*np.random.uniform(-1, 1)


    #238U - 134Te
    area_double_true_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_238U_134Te[0], P_double_new_238U_134Te[1], P_double_new_238U_134Te[2], P_double_new_238U_134Te[3], P_double_new_238U_134Te[4], P_double_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_new_238U_134Te[0], P_double_new_238U_134Te[1], P_double_new_238U_134Te[2], P_double_new_238U_134Te[3], P_double_new_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)
    area_double_true_prompt_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_238U_134Te[0], P_double_new_238U_134Te[1], P_double_new_238U_134Te[2], P_double_new_238U_134Te[3], P_double_new_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)
    area_double_true_delayed_new_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_new_238U_134Te[0], P_double_new_238U_134Te[1], P_double_new_238U_134Te[2], P_double_new_238U_134Te[3], P_double_new_238U_134Te[4], P_double_238U_134Te[5]), x_arr_134Te)
    IYR_array_238U_134Te[n] = IYR(prompt=area_double_true_prompt_new_238U_134Te, delayed=area_double_true_delayed_new_238U_134Te)

    #238U - 134Te 1n
    area_double_1n_true_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_new_238U_134Te[0], P_double_1n_new_238U_134Te[1], P_double_1n_new_238U_134Te[2], P_double_1n_new_238U_134Te[3], P_double_1n_new_238U_134Te[4], P_double_1n_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_1n_new_238U_134Te[0], P_double_1n_new_238U_134Te[1], P_double_1n_new_238U_134Te[2], P_double_1n_new_238U_134Te[3], P_double_1n_new_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)
    area_double_1n_true_prompt_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_new_238U_134Te[0], P_double_1n_new_238U_134Te[1], P_double_1n_new_238U_134Te[2], P_double_1n_new_238U_134Te[3], P_double_1n_new_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)
    area_double_1n_true_delayed_new_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_1n_new_238U_134Te[0], P_double_1n_new_238U_134Te[1], P_double_1n_new_238U_134Te[2], P_double_1n_new_238U_134Te[3], P_double_1n_new_238U_134Te[4], P_double_1n_238U_134Te[5]), x_arr_134Te)
    IYR_array_1n_238U_134Te[n] = IYR(prompt=area_double_1n_true_prompt_new_238U_134Te, delayed=area_double_1n_true_delayed_new_238U_134Te)

    #238U - 134Te 3n
    area_double_3n_true_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_new_238U_134Te[0], P_double_3n_new_238U_134Te[1], P_double_3n_new_238U_134Te[2], P_double_3n_new_238U_134Te[3], P_double_3n_new_238U_134Te[4], P_double_3n_238U_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_new_238U_134Te[0], P_double_3n_new_238U_134Te[1], P_double_3n_new_238U_134Te[2], P_double_3n_new_238U_134Te[3], P_double_3n_new_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)
    area_double_3n_true_prompt_new_238U_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_new_238U_134Te[0], P_double_3n_new_238U_134Te[1], P_double_3n_new_238U_134Te[2], P_double_3n_new_238U_134Te[3], P_double_3n_new_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)
    area_double_3n_true_delayed_new_238U_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_new_238U_134Te[0], P_double_3n_new_238U_134Te[1], P_double_3n_new_238U_134Te[2], P_double_3n_new_238U_134Te[3], P_double_3n_new_238U_134Te[4], P_double_3n_238U_134Te[5]), x_arr_134Te)
    IYR_array_3n_238U_134Te[n] = IYR(prompt=area_double_3n_true_prompt_new_238U_134Te, delayed=area_double_3n_true_delayed_new_238U_134Te)

    #232Th - 134Te
    area_double_true_new_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_232Th_134Te[0], P_double_new_232Th_134Te[1], P_double_new_232Th_134Te[2], P_double_new_232Th_134Te[3], P_double_new_232Th_134Te[4], P_double_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_new_232Th_134Te[0], P_double_new_232Th_134Te[1], P_double_new_232Th_134Te[2], P_double_new_232Th_134Te[3], P_double_new_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
    area_double_true_prompt_new_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_new_232Th_134Te[0], P_double_new_232Th_134Te[1], P_double_new_232Th_134Te[2], P_double_new_232Th_134Te[3], P_double_new_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
    area_double_true_delayed_new_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_new_232Th_134Te[0], P_double_new_232Th_134Te[1], P_double_new_232Th_134Te[2], P_double_new_232Th_134Te[3], P_double_new_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
    IYR_array_232Th_134Te[n] = IYR(prompt=area_double_true_prompt_new_232Th_134Te, delayed=area_double_true_delayed_new_232Th_134Te)

    #238U - 135Te
    area_double_true_new_238U_135Te = np.trapz(gauss(x_arr_135Te, P_double_new_238U_135Te[0], P_double_new_238U_135Te[1], P_double_new_238U_135Te[2], P_double_new_238U_135Te[3], P_double_new_238U_135Te[4], P_double_238U_135Te[5]) + smeared_exp_decay(x_arr_135Te, P_double_new_238U_135Te[0], P_double_new_238U_135Te[1], P_double_new_238U_135Te[2], P_double_new_238U_135Te[3], P_double_new_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)
    area_double_true_prompt_new_238U_135Te = np.trapz(gauss(x_arr_135Te, P_double_new_238U_135Te[0], P_double_new_238U_135Te[1], P_double_new_238U_135Te[2], P_double_new_238U_135Te[3], P_double_new_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)
    area_double_true_delayed_new_238U_135Te = np.trapz(smeared_exp_decay(x_arr_135Te, P_double_new_238U_135Te[0], P_double_new_238U_135Te[1], P_double_new_238U_135Te[2], P_double_new_238U_135Te[3], P_double_new_238U_135Te[4], P_double_238U_135Te[5]), x_arr_135Te)
    IYR_array_238U_135Te[n] = IYR(prompt=area_double_true_prompt_new_238U_135Te, delayed=area_double_true_delayed_new_238U_135Te)


#238U - 134Te
sigma_IYR_array_238U_134Te = (np.max(IYR_array_238U_134Te) - np.min(IYR_array_238U_134Te))/2.0

#238U - 134Te 1n
sigma_IYR_array_1n_238U_134Te = (np.max(IYR_array_1n_238U_134Te) - np.min(IYR_array_1n_238U_134Te))/2.0

#238U - 134Te 3n
sigma_IYR_array_3n_238U_134Te = (np.max(IYR_array_3n_238U_134Te) - np.min(IYR_array_3n_238U_134Te))/2.0

#232Th - 134Te
sigma_IYR_array_232Th_134Te = (np.max(IYR_array_232Th_134Te) - np.min(IYR_array_232Th_134Te))/2.0

#238U - 135Te
sigma_IYR_array_238U_135Te = (np.max(IYR_array_238U_135Te) - np.min(IYR_array_238U_135Te))/2.0



####################################################
###        		Print table of IYRs              ###
####################################################

print("\n")
print(" ***** Isomeric Yield Ratios ****")
print("IYR 238U - 134Te:               %.2f +/- %.2f" % (np.average(IYR_array_238U_134Te), sigma_IYR_array_238U_134Te))
print("IYR 238U - 134Te 1n:            %.2f +/- %.2f" % (np.average(IYR_array_1n_238U_134Te), sigma_IYR_array_1n_238U_134Te))
print("IYR 238U - 134Te 3n:            %.2f +/- %.2f" % (np.average(IYR_array_3n_238U_134Te), sigma_IYR_array_3n_238U_134Te))
#print("IYR 232Th - 134Te:               %.2f +/- %.2f" % (np.average(IYR_array_232Th_134Te), sigma_IYR_array_232Th_134Te))

#print("IYR 238U - 135Te:               %.3f +/- %.3f" % (np.average(IYR_array_238U_135Te), sigma_IYR_array_238U_135Te))

print("\n")




####################################################
## 						Plot 		              ## 
####################################################

x_array_plot = np.linspace(0,1000,10000)


################   238U -  134Te   #################

plt.plot(x_doublegate_238U_134Te_long, y_doublegate_238U_134Te_long, label="doublegate_238U_134Te", color="royalblue")

plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), label="true fit, total", color="orange")
plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), label="true gaussian", color="green")
plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), label="true smeared exp decay", color="red")
plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_134Te[0], P_double_238U_134Te[1], P_double_238U_134Te[2], P_double_238U_134Te[3], P_double_238U_134Te[4], P_double_238U_134Te[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_238U_134Te[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_238U_134Te[-1],0,6000, color="black")
#plt.yscale("log")
plt.title("238U - 134Te: Doublegate true spectrum fit")
#plt.axis([0,700,1,10**(4)])
plt.axis([0,600,0,4*10**(3)])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=14)
plt.grid()
plt.show()

################   238U -  134Te  1n  #################

plt.plot(x_doublegate_1n_238U_134Te_long, y_doublegate_1n_238U_134Te_long, label="doublegate_1n_238U_134Te", color="royalblue")

plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), label="true fit, total", color="orange")
plt.plot(x_array_plot, gauss(x_array_plot, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), label="true gaussian", color="green")
plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), label="true smeared exp decay", color="red")
plt.plot(x_array_plot, const_bg(x_array_plot, P_double_1n_238U_134Te[0], P_double_1n_238U_134Te[1], P_double_1n_238U_134Te[2], P_double_1n_238U_134Te[3], P_double_1n_238U_134Te[4], P_double_1n_238U_134Te[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_1n_238U_134Te[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_1n_238U_134Te[-1],0,6000, color="black")
#plt.yscale("log")
plt.title("238U - 134Te 1n: Doublegate true spectrum fit")
plt.axis([0,700,0,2000])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()

################   238U -  134Te  3n  #################

plt.plot(x_doublegate_3n_238U_134Te_long, y_doublegate_3n_238U_134Te_long, label="doublegate_3n_238U_134Te", color="royalblue")

plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), label="true fit, total", color="orange")
plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), label="true gaussian", color="green")
plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), label="true smeared exp decay", color="red")
plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_238U_134Te[0], P_double_3n_238U_134Te[1], P_double_3n_238U_134Te[2], P_double_3n_238U_134Te[3], P_double_3n_238U_134Te[4], P_double_3n_238U_134Te[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_3n_238U_134Te[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_3n_238U_134Te[-1],0,6000, color="black")
#plt.yscale("log")
plt.title("238U - 134Te 3n: Doublegate true spectrum fit")
plt.axis([0,700,0,5*10**(3)])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()


################   232Th -  134Te   #################

# plt.plot(x_doublegate_232Th_134Te_long, y_doublegate_232Th_134Te_long, label="doublegate_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th - 134Te: Doublegate true spectrum fit")
# plt.axis([0,700,1,10**(4)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()



################   238U -  135Te   #################

# plt.plot(x_doublegate_238U_135Te_long, y_doublegate_238U_135Te_long, label="doublegate_238U_135Te", color="royalblue")
# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_135Te[0], P_double_238U_135Te[1], P_double_238U_135Te[2], P_double_238U_135Te[3], P_double_238U_135Te[4], P_double_238U_135Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_135Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_135Te[-1],0,6000, color="black")
# plt.yscale("log")
# plt.title("238U - 135Te: Doublegate true spectrum fit")
# plt.axis([0,700,1,10**(4)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()



# plt.plot(x_doublegate_238U_134Te_long, 6.8*y_doublegate_238U_134Te_long, label="doublegate_238U_134Te", color="royalblue")
# plt.plot(x_doublegate_232Th_134Te_long, y_doublegate_232Th_134Te_long, label="doublegate_232Th_134Te", color="orchid")

# plt.title("134Te: 232Th vs 238U")
# #plt.axis([0,700,1,10**(4)])
# plt.axis([0,700,0,10**(5)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()





