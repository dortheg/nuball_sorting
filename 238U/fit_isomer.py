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

x_doublegate_134Te_long = x_doublegate_134Te
y_doublegate_134Te_long = y_doublegate_134Te

x_doublegate_134Te = x_doublegate_134Te[x_lower:x_upper]
y_doublegate_134Te = y_doublegate_134Te[x_lower:x_upper]

#######################

#Doublegate all
hist_doublegate_all_134Te = file.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_134Te.GetNbinsX()

x_doublegate_all_134Te = np.zeros(x_bins)
y_doublegate_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinCenter(i+1)
    y_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_134Te = 2*x_doublegate_all_134Te

x_doublegate_all_134Te_long = x_doublegate_all_134Te
y_doublegate_all_134Te_long = y_doublegate_all_134Te

x_doublegate_all_134Te = x_doublegate_all_134Te[x_lower:x_upper]
y_doublegate_all_134Te = y_doublegate_all_134Te[x_lower:x_upper]

#######################

#Doublegate bg
hist_doublegate_bg_134Te = file.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_134Te.GetNbinsX()

x_doublegate_bg_134Te = np.zeros(x_bins)
y_doublegate_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinCenter(i+1)
    y_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinContent(i+1)


#Make up for 2ns bins
x_doublegate_bg_134Te = 2*x_doublegate_bg_134Te

x_doublegate_bg_134Te_long = x_doublegate_bg_134Te
y_doublegate_bg_134Te_long = y_doublegate_bg_134Te

x_doublegate_bg_134Te = x_doublegate_bg_134Te[x_lower:x_upper]
y_doublegate_bg_134Te = y_doublegate_bg_134Te[x_lower:x_upper]


#######################

#Doublegate bg ridge
hist_doublegate_bg_ridge_134Te = file.Get('time_isomer_doublegate_bg_ridge_134Te')
x_bins = hist_doublegate_bg_ridge_134Te.GetNbinsX()

x_doublegate_bg_ridge_134Te = np.zeros(x_bins)
y_doublegate_bg_ridge_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_ridge_134Te[i] = hist_doublegate_bg_ridge_134Te.GetBinCenter(i+1)
    y_doublegate_bg_ridge_134Te[i] = hist_doublegate_bg_ridge_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_ridge_134Te = 2*x_doublegate_bg_ridge_134Te

x_doublegate_bg_ridge_134Te_long = x_doublegate_bg_ridge_134Te
y_doublegate_bg_ridge_134Te_long = y_doublegate_bg_ridge_134Te

x_doublegate_bg_ridge_134Te = x_doublegate_bg_ridge_134Te[x_lower:x_upper]
y_doublegate_bg_ridge_134Te = y_doublegate_bg_ridge_134Te[x_lower:x_upper]

#######################

#Doublegate bg random
hist_doublegate_bg_random_134Te = file.Get('time_isomer_doublegate_bg_random_134Te')
x_bins = hist_doublegate_bg_random_134Te.GetNbinsX()

x_doublegate_bg_random_134Te = np.zeros(x_bins)
y_doublegate_bg_random_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_random_134Te[i] = hist_doublegate_bg_random_134Te.GetBinCenter(i+1)
    y_doublegate_bg_random_134Te[i] = hist_doublegate_bg_random_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_random_134Te = 2*x_doublegate_bg_random_134Te

x_doublegate_bg_random_134Te_long = x_doublegate_bg_random_134Te
y_doublegate_bg_random_134Te_long = y_doublegate_bg_random_134Te

x_doublegate_bg_random_134Te = x_doublegate_bg_random_134Te[x_lower:x_upper]
y_doublegate_bg_random_134Te = y_doublegate_bg_random_134Te[x_lower:x_upper]


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
	"Not same as for 252Cf, due to both fragments being stopped"
	return (delayed)/(prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.03*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.03*bg_delayed)**2)
    return np.sqrt( ((prompt/(prompt+delayed)**2)*sigma_delayed)**2 + ((delayed/(prompt+delayed)**2)*sigma_prompt)**2 )

def sigma_data_doublegate(data_all, data_bg_ridge, data_bg_random):
    return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random + (0.03*data_bg_ridge)**2 + (0.03*data_bg_random)**2)
    #return np.sqrt(data_all + 0.25*data_bg_ridge + 0.125*data_bg_random)


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

#P_double, cov_double = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_134Te, y_doublegate_134Te, sigma=sigma_data_doublegate(y_doublegate_all_134Te, y_doublegate_bg_ridge_134Te, y_doublegate_bg_random_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
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


#All
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

P_double_all, cov_double_all = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_all_134Te, y_doublegate_all_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))

# print("\n")
# print(" ***** 134Te:  Doublegate all spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_all[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_all[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f         [%.d,%.d]" % (P_double_all[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_all[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_all[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_all[5]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


#bg
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

P_double_bg, cov_double_bg = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_bg_134Te, y_doublegate_bg_134Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))

# print("\n")
# print(" ***** 134Te:  Doublegate bg spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f         [%.d,%.d]" % (P_double_bg[0], mean_lower, mean_upper))
# print("sigma:                    %.2f         [%.d,%.d]" % (P_double_bg[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f         [%.d,%.d]" % (P_double_bg[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f         [%.d,%.d]" % (P_double_bg[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_bg[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_bg[5]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


######################################
## 		 Find IYR + uncertainty 	## 
######################################

x_arr = np.linspace(0,1500,4000)

area_double_true = np.trapz(gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]) + smeared_exp_decay(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)
area_double_true_prompt = np.trapz(gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)
area_double_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), x_arr)

# print("Area double true tot: %.3f" % area_double_true)
# print("Area double true prompt: %.3f" % area_double_true_prompt)
# print("Area double true delayed: %.3f" % area_double_true_delayed)

area_double_all = np.trapz(gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]) + smeared_exp_decay(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)
area_double_all_prompt = np.trapz(gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)
area_double_all_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), x_arr)

area_double_bg = np.trapz(gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]) + smeared_exp_decay(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)
area_double_bg_prompt = np.trapz(gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)
area_double_bg_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), x_arr)

IYR_double = IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed)
sigma_IYR_double = sigma_IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed, all_prompt=area_double_all_prompt, all_delayed=area_double_all_delayed, bg_prompt=area_double_bg_prompt, bg_delayed=area_double_bg_delayed)

print("\n")
print(" ***** Isomeric Yield Ratio ****")
print("IYR_double:               %.3f +/- %.3f" % (IYR_double, sigma_IYR_double))
print("\n")

##########################
## 			Plot 		## 
##########################

#True spectrum
#plt.plot(x_doublegate_134Te, y_doublegate_134Te, label="doublegate_134Te", color="royalblue")
#plt.plot(x_doublegate_134Te_long, y_doublegate_134Te_long, label="doublegate_134Te", color="royalblue")
plt.errorbar(x_doublegate_134Te_long, y_doublegate_134Te_long, yerr=sigma_data_doublegate(y_doublegate_all_134Te_long, y_doublegate_bg_ridge_134Te_long, y_doublegate_bg_random_134Te_long), label="doublegate_134Te", color="royalblue")

plt.plot(x_arr, sum_smeared_exp_gauss_const_bg(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true fit, total", color="orange")
plt.plot(x_arr, gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true gaussian", color="green")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="true smeared exp decay", color="red")
plt.plot(x_arr, const_bg(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5]), label="constant BG", color="hotpink")

plt.vlines(x_doublegate_134Te[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_134Te[-1],0,6000, color="black")

plt.title("134Te: Doublegate true spectrum fit")
#plt.axis([800,2000,0,500])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()


# #All spectrum
# plt.plot(x_doublegate_all_134Te_long, y_doublegate_all_134Te_long, label="doublegate_all_134Te", color="black")

# plt.plot(x_arr, sum_smeared_exp_gauss_const_bg(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="true smeared exp decay", color="red")
# plt.plot(x_arr, const_bg(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4], P_double_all[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_all_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_all_134Te[-1],0,6000, color="black")

# plt.title("134Te: Doublegate all spectrum fit")
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


# #BG spectrum
# plt.plot(x_doublegate_bg_134Te_long, y_doublegate_bg_134Te_long, label="doublegate_bg_134Te", color="lightpink")

# plt.plot(x_arr, sum_smeared_exp_gauss_const_bg(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="true smeared exp decay", color="red")
# plt.plot(x_arr, const_bg(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4], P_double_bg[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_bg_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_bg_134Te[-1],0,6000, color="black")

# plt.title("134Te: Doublegate bg spectrum fit")
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()



