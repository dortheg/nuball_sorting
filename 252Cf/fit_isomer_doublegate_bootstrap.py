import numpy as np
import ROOT 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time

BOOTSTRAP = False

##########################
##     Read in data     ## 
##########################

file = ROOT.TFile.Open("252Cf_11jan2022.root"," READ ")


####################################################
###           1279-time gated spectra
####################################################

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

def sigma_data_doublegate_all_bg(data_all, data_bg):
    unc = np.zeros(len(data_all))

    for i in range(len(data_all)):
        if data_all[i] <= 0 and data_bg[i]<= 0:
            unc[i] = 1
        else:
            unc[i] = np.sqrt(data_all[i]+data_bg[i])
    return unc

###################################
##      Define fitting func      ## 
###################################


def sum_two_smeared_exp_two_gauss_const_bg(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    return amplitude_gauss1*np.exp(-(x-mean1)**2/(2*sigma1**2)) \
    + amplitude_gauss2*np.exp(-(x-mean2)**2/(2*sigma2**2)) \
    + gaussian_filter1d(np.piecewise(x, [x < mean1, x >= mean1], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean1-x)/tau_Ge)]),sigma1) \
    + gaussian_filter1d(np.piecewise(x, [x < mean1, x >= mean1], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean1-x)/tau_decay)]),sigma1)  \
    + const_bg


def gauss_1(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            """Gaussian component (prompt production and prompt decay)"""
            return amplitude_gauss1*np.exp(-(x-mean1)**2/(2*sigma1**2))

def gauss_2(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            """Gaussian component (prompt production and prompt decay)"""
            return amplitude_gauss2*np.exp(-(x-mean2)**2/(2*sigma2**2))

def const_bg(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            return const_bg + x*0

def smeared_exp_Ge(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean1, x >= mean1], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean1-x)/tau_Ge)]),sigma1)


def smeared_exp_decay(x, mean1=0, mean2=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean1, x >= mean1], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean1-x)/tau_decay)]),sigma1)


####################################################
## 		             Fit data 		              ## 
####################################################

######################
### Doublegate_2
######################

#True
mean1_lower = 950
mean1_upper = 1100
mean2_lower = 950
mean2_upper = 1100
sigma1_lower = 0
sigma1_upper = 40
amplitude_gauss1_lower = 0
amplitude_gauss1_upper = 3000
sigma2_lower = 0
sigma2_upper = 40
amplitude_gauss2_lower = 0
amplitude_gauss2_upper = 40
const_bg_lower = 0
const_bg_upper = 0+0.001
amplitude_exp_Ge_lower = 0
amplitude_exp_Ge_upper = 1000
tau_Ge_lower = 0
tau_Ge_upper = 40
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 1000
tau_decay_lower = tau_134Te#0
tau_decay_upper = tau_134Te+0.0001#1000

unc_y_doublegate_2_134Te = sigma_data_doublegate_all_bg(data_all=y_doublegate_2_all_134Te, data_bg=y_doublegate_2_bg_134Te)

P_double_2, cov_double_2 = curve_fit(sum_two_smeared_exp_two_gauss_const_bg, x_doublegate_2_134Te, y_doublegate_2_134Te, sigma=unc_y_doublegate_2_134Te, bounds=([mean1_lower,mean2_lower,sigma1_lower,amplitude_gauss1_lower,sigma2_lower,amplitude_gauss2_lower,const_bg_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean1_upper,mean1_upper,sigma1_upper,amplitude_gauss1_upper,sigma2_upper,amplitude_gauss2_upper,const_bg_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]))

print("\n")
print(" ***** 134Te:  Doublegate_2 true spectrum fit ***** ")
print("          -- GAUSS + TWO SMEARED EXP FIT --   ")
print("mean1:                     %.2f         [%.d,%.d]" % (P_double_2[0], mean1_lower, mean1_upper))
print("sigma1:                    %.2f         [%.d,%.d]" % (P_double_2[2], sigma1_lower, sigma1_upper))
print("amplitude_gauss1:          %.2f         [%.d,%.d]" % (P_double_2[3], amplitude_gauss1_lower, amplitude_gauss1_upper))
print("sigma2:                    %.2f         [%.d,%.d]" % (P_double_2[4], sigma2_lower, sigma2_upper))
print("amplitude_gauss2:          %.2f         [%.d,%.d]" % (P_double_2[5], amplitude_gauss2_lower, amplitude_gauss2_upper))
print("const_bg:                   %.2f         [%.d,%.d]" % (P_double_2[6], const_bg_lower, const_bg_upper))
print("amplitude_exp_Ge:         %.2f         [%.d,%.d]" % (P_double_2[7], amplitude_exp_Ge_lower, amplitude_exp_Ge_upper))
print("tau_Ge:                   %.2f         [%.d,%.d]" % (P_double_2[8], tau_Ge_lower, tau_Ge_upper))
print("amplitude_exp_decay:      %.2f         [%.d,%.d]" % (P_double_2[9], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
print("tau_decay, in half_life:  %.2f         [%.d,%.d]" % (P_double_2[10]*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
print("\n")


####################################################
## 		     		 Find IYR 	                  ## 
####################################################

x_arr = np.linspace(-1000,10000,10000)

############
### Doublegate_2
############

area_double_2_true_prompt = np.trapz( \
    gauss_1(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]) \
    +gauss_2(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]) \
    +smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), x_arr)
area_double_2_true_delayed = np.trapz(smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), x_arr)

IYR_double_2 = IYR(prompt=area_double_2_true_prompt, delayed=area_double_2_true_delayed)



####################################################
##                BOOTSTRAPPING                   ## START_BOOTSTRAP
####################################################

if BOOTSTRAP==True:

    N_BOOTSTRAP = 5000
    print("\n Starting BOOTSTRAPPING... Number of iterations: %.d" % N_BOOTSTRAP)

    resampled_y_doublegate_2_134Te = np.zeros(len(y_doublegate_2_134Te))
    resampled_IYR_array_double_2 = np.zeros(N_BOOTSTRAP)

    for n in range(N_BOOTSTRAP):

        y_doublegate_2_134Te_bgvaried = np.zeros(len(y_doublegate_2_134Te))

        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_2_134Te_bgvaried = y_doublegate_2_all_134Te-y_doublegate_2_bg_134Te*np.random.normal(1,0.025)
        
        for i in range (len(y_doublegate_2_134Te_bgvaried)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_2_134Te[i] = y_doublegate_2_134Te_bgvaried[i] + np.random.normal(0, unc_y_doublegate_2_134Te[i]) 


        #True
        mean_lower = 950
        mean_upper = 1100
        sigma1_lower = 0
        sigma1_upper = 40
        amplitude_gauss1_lower = 0
        amplitude_gauss1_upper = 3000
        sigma2_lower = 0
        sigma2_upper = 40
        amplitude_gauss2_lower = 0
        amplitude_gauss2_upper = 40
        const_bg_lower = 0
        const_bg_upper = 0+0.001
        amplitude_exp_Ge_lower = 0
        amplitude_exp_Ge_upper = 1000
        tau_Ge_lower = 0
        tau_Ge_upper = 40
        amplitude_exp_decay_lower = 0
        amplitude_exp_decay_upper = 1000
        tau_decay_lower = tau_134Te#0
        tau_decay_upper = tau_134Te+0.0001#1000

        try:
            resampled_P_double_2, resampled_cov_double_2 = curve_fit(sum_two_smeared_exp_two_gauss_const_bg, x_doublegate_2_134Te, resampled_y_doublegate_2_134Te, sigma=unc_y_doublegate_2_134Te, bounds=([mean_lower,sigma1_lower,amplitude_gauss1_lower,sigma2_lower,amplitude_gauss2_lower,const_bg_lower,amplitude_exp_Ge_lower,tau_Ge_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma1_upper,amplitude_gauss1_upper,sigma2_upper,amplitude_gauss2_upper,const_bg_upper,amplitude_exp_Ge_upper,tau_Ge_upper,amplitude_exp_decay_upper,tau_decay_upper]), maxfev=10000)

            resampled_area_double_2_true_prompt = np.trapz( \
                gauss_1(x_arr, resampled_P_double_2[0], resampled_P_double_2[1], resampled_P_double_2[2], resampled_P_double_2[3], resampled_P_double_2[4], resampled_P_double_2[5], resampled_P_double_2[6], resampled_P_double_2[7], resampled_P_double_2[8], resampled_P_double_2[9]) \
                +gauss_2(x_arr, resampled_P_double_2[0], resampled_P_double_2[1], resampled_P_double_2[2], resampled_P_double_2[3], resampled_P_double_2[4], resampled_P_double_2[5], resampled_P_double_2[6], resampled_P_double_2[7], resampled_P_double_2[8], resampled_P_double_2[9]) \
                +smeared_exp_Ge(x_arr, resampled_P_double_2[0], resampled_P_double_2[1], resampled_P_double_2[2], resampled_P_double_2[3], resampled_P_double_2[4], resampled_P_double_2[5], resampled_P_double_2[6], resampled_P_double_2[7], resampled_P_double_2[8], resampled_P_double_2[9]), x_arr)
            resampled_area_double_2_true_delayed = np.trapz(smeared_exp_decay(x_arr, resampled_P_double_2[0], resampled_P_double_2[1], resampled_P_double_2[2], resampled_P_double_2[3], resampled_P_double_2[4], resampled_P_double_2[5], resampled_P_double_2[6], resampled_P_double_2[7], resampled_P_double_2[8], resampled_P_double_2[9]), x_arr)

            resampled_IYR_array_double_2[n] = IYR(prompt=resampled_area_double_2_true_prompt, delayed=resampled_area_double_2_true_delayed)

        except RuntimeError:
            resampled_IYR_array_double_2[n] = IYR_double_2

        if(n%100==0):
            print("Now finished iteration %.d" % n)

sigma_bootstrap_IYR_double_2 = 0
sigma_IYR_handcalc = 0 

if BOOTSTRAP==True:

    sigma_bootstrap_IYR_double_2 = np.std(resampled_IYR_array_double_2) #Change this 


    # Plot distribution of calculated IYRs
    IYR_binwith = 0.001
    IYR_min = 0.0
    IYR_max = 1.0 + 2*IYR_binwith
    N_IYR_bins = int((IYR_max-IYR_min)//IYR_binwith)

    IYR_array = np.linspace(IYR_min,IYR_max,N_IYR_bins)
    IYR_histogram = np.zeros(len(IYR_array))

    for i in range(N_BOOTSTRAP):
        IYR_bin = int(resampled_IYR_array_double_2[i]//IYR_binwith - IYR_min//IYR_binwith)
        IYR_histogram[IYR_bin] += 1

        # Calculate standard deviation by hand
        std_value = np.linspace(0,0.5,10000)
        std_found = False
        sigma_IYR_handcalc = 0

        for i in range(len(std_value)):
            within_std = np.where((resampled_IYR_array_double_2 >= IYR_double_2-std_value[i]) & (resampled_IYR_array_double_2 < IYR_double_2+std_value[i]))[0]
            percentage_covered = len(within_std)/len(resampled_IYR_array_double_2)

            if 0.675 <= percentage_covered < 0.685:
                std_found = True
                sigma_IYR_handcalc = std_value[i]
                #print("By-hand standard deviation is %.3f" % std_value[i])
                break

    plt.plot(IYR_array,IYR_histogram, label="Calculated IYR distribution")
    plt.axvline(x=IYR_double_2, ymin=0, ymax=max(IYR_histogram), label="IYR center", color="red")
    plt.axvline(x=IYR_double_2-sigma_bootstrap_IYR_double_2, ymin=0, ymax=max(IYR_histogram), label="IYR-sigma_np", color="black")
    plt.axvline(x=IYR_double_2+sigma_bootstrap_IYR_double_2, ymin=0, ymax=max(IYR_histogram), label="IYR+sigma_np", color="black")
    plt.axvline(x=IYR_double_2-sigma_IYR_handcalc, ymin=0, ymax=max(IYR_histogram), label="IYR-sigma_handcalc", color="grey")
    plt.axvline(x=IYR_double_2+sigma_IYR_handcalc, ymin=0, ymax=max(IYR_histogram), label="IYR+sigma_handcalc", color="grey")
    plt.xlabel("IYR")
    plt.ylabel("Counts")
    plt.legend()
    plt.title(" Bootstrapping IYR, N= %d " % N_BOOTSTRAP)
    plt.grid()
    plt.show()


print(" ***** Isomeric Yield Ratios ****")
print("IYR_double_2:               %.3f +/- %.3f" % (IYR_double_2, sigma_bootstrap_IYR_double_2))
print("IYR_double_2:               %.3f +/- %.3f" % (IYR_double_2, sigma_IYR_handcalc))


#################
###  Doublegate_2
################

#true
#plt.plot(x_doublegate_2_134Te, y_doublegate_2_134Te, label="doublegate_2_134Te", color="royalblue")
#plt.plot(x_doublegate_2_134Te, resampled_y_doublegate_2_134Te, label="resampled doublegate_2_134Te", color="purple")
#plt.errorbar(x_doublegate_2_134Te, y_doublegate_2_134Te, yerr=unc_y_doublegate_2_134Te, label="doublegate_2_134Te", color="aqua")
plt.errorbar(x_doublegate_2_134Te, y_doublegate_2_134Te, yerr=unc_y_doublegate_2_134Te, label="doublegate_2_134Te", color="aqua", linestyle="None")

plt.plot(x_arr, sum_two_smeared_exp_two_gauss_const_bg(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), label="true fit, total", color="k")
plt.plot(x_arr, gauss_1(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), label="gauss_1")
plt.plot(x_arr, gauss_2(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), label="gauss_2")
plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), label="smeared exp Ge")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double_2[0], P_double_2[1], P_double_2[2], P_double_2[3], P_double_2[4], P_double_2[5], P_double_2[6], P_double_2[7], P_double_2[8], P_double_2[9], P_double_2[10]), label="smeared exp decay")


plt.title("134Te: Doublegate_2 true spectrum fit")
plt.axis([800,2000,0,500])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=10)
plt.grid()
plt.show()
