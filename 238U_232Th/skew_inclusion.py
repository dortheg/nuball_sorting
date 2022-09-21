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

####################################################
##     				Read in data                  ## START_DATA
####################################################

file_238U_lowE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.bin_final.root"," READ ")

#Define lower and upper fit limit
x_lower = 275
x_upper = 400

bin_lower = x_lower//2
bin_upper = x_upper//2


################   238U lowE_highT -  140Xe   #################
#Doublegate true
hist_doublegate_238U_lowE_highT_140Xe = file_238U_lowE_highT.Get('time_isomer_doublegate_140Xe')
x_bins = hist_doublegate_238U_lowE_highT_140Xe.GetNbinsX()

x_doublegate_238U_lowE_highT_140Xe = np.zeros(x_bins)
y_doublegate_238U_lowE_highT_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_lowE_highT_140Xe[i] = hist_doublegate_238U_lowE_highT_140Xe.GetBinCenter(i+1)
    y_doublegate_238U_lowE_highT_140Xe[i] = hist_doublegate_238U_lowE_highT_140Xe.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_lowE_highT_140Xe = 2*x_doublegate_238U_lowE_highT_140Xe

x_doublegate_238U_lowE_highT_140Xe_long = x_doublegate_238U_lowE_highT_140Xe
y_doublegate_238U_lowE_highT_140Xe_long = y_doublegate_238U_lowE_highT_140Xe

x_doublegate_238U_lowE_highT_140Xe = x_doublegate_238U_lowE_highT_140Xe[bin_lower:bin_upper]
y_doublegate_238U_lowE_highT_140Xe = y_doublegate_238U_lowE_highT_140Xe[bin_lower:bin_upper]

#print("* 238U - 140Xe fit range %d - %.d" % (x_doublegate_238U_lowE_highT_140Xe[0],x_doublegate_238U_lowE_highT_140Xe[-1]))


################   238U lowE_highT all -  140Xe   #################

#doublegate_all
hist_doublegate_all_238U_lowE_highT_140Xe = file_238U_lowE_highT.Get('time_isomer_doublegate_all_140Xe')
x_bins = hist_doublegate_all_238U_lowE_highT_140Xe.GetNbinsX()

x_doublegate_all_238U_lowE_highT_140Xe = np.zeros(x_bins)
y_doublegate_all_238U_lowE_highT_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_238U_lowE_highT_140Xe[i] = hist_doublegate_all_238U_lowE_highT_140Xe.GetBinCenter(i+1)
    y_doublegate_all_238U_lowE_highT_140Xe[i] = hist_doublegate_all_238U_lowE_highT_140Xe.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_238U_lowE_highT_140Xe = 2*x_doublegate_all_238U_lowE_highT_140Xe

x_doublegate_all_238U_lowE_highT_140Xe_long = x_doublegate_all_238U_lowE_highT_140Xe
y_doublegate_all_238U_lowE_highT_140Xe_long = y_doublegate_all_238U_lowE_highT_140Xe

x_doublegate_all_238U_lowE_highT_140Xe = x_doublegate_all_238U_lowE_highT_140Xe[bin_lower:bin_upper]
y_doublegate_all_238U_lowE_highT_140Xe = y_doublegate_all_238U_lowE_highT_140Xe[bin_lower:bin_upper]

#print("* 238U - 140Xe fit range %d - %.d" % (x_doublegate_all_238U_lowE_highT_140Xe[0],x_doublegate_all_238U_lowE_highT_140Xe[-1]))



################   238U lowE_highT bg -  140Xe   #################

#doublegate_bg
hist_doublegate_bg_238U_lowE_highT_140Xe = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_140Xe')
x_bins = hist_doublegate_bg_238U_lowE_highT_140Xe.GetNbinsX()

x_doublegate_bg_238U_lowE_highT_140Xe = np.zeros(x_bins)
y_doublegate_bg_238U_lowE_highT_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_238U_lowE_highT_140Xe.GetBinCenter(i+1)
    y_doublegate_bg_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_238U_lowE_highT_140Xe.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_238U_lowE_highT_140Xe = 2*x_doublegate_bg_238U_lowE_highT_140Xe

x_doublegate_bg_238U_lowE_highT_140Xe_long = x_doublegate_bg_238U_lowE_highT_140Xe
y_doublegate_bg_238U_lowE_highT_140Xe_long = y_doublegate_bg_238U_lowE_highT_140Xe

x_doublegate_bg_238U_lowE_highT_140Xe = x_doublegate_bg_238U_lowE_highT_140Xe[bin_lower:bin_upper]
y_doublegate_bg_238U_lowE_highT_140Xe = y_doublegate_bg_238U_lowE_highT_140Xe[bin_lower:bin_upper]

#print("* 238U - 140Xe fit range %d - %.d" % (x_doublegate_bg_238U_lowE_highT_140Xe[0],x_doublegate_bg_238U_lowE_highT_140Xe[-1]))


################   238U lowE_highT bg_ridge -  140Xe   #################

#doublegate_bg_ridge
hist_doublegate_bg_ridge_238U_lowE_highT_140Xe = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_ridge_140Xe')
x_bins = hist_doublegate_bg_ridge_238U_lowE_highT_140Xe.GetNbinsX()

x_doublegate_bg_ridge_238U_lowE_highT_140Xe = np.zeros(x_bins)
y_doublegate_bg_ridge_238U_lowE_highT_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_ridge_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_ridge_238U_lowE_highT_140Xe.GetBinCenter(i+1)
    y_doublegate_bg_ridge_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_ridge_238U_lowE_highT_140Xe.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_ridge_238U_lowE_highT_140Xe = 2*x_doublegate_bg_ridge_238U_lowE_highT_140Xe

x_doublegate_bg_ridge_238U_lowE_highT_140Xe_long = x_doublegate_bg_ridge_238U_lowE_highT_140Xe
y_doublegate_bg_ridge_238U_lowE_highT_140Xe_long = y_doublegate_bg_ridge_238U_lowE_highT_140Xe

x_doublegate_bg_ridge_238U_lowE_highT_140Xe = x_doublegate_bg_ridge_238U_lowE_highT_140Xe[bin_lower:bin_upper]
y_doublegate_bg_ridge_238U_lowE_highT_140Xe = y_doublegate_bg_ridge_238U_lowE_highT_140Xe[bin_lower:bin_upper]

#print("* 238U - 140Xe fit range %d - %.d" % (x_doublegate_bg_ridge_238U_lowE_highT_140Xe[0],x_doublegate_bg_ridge_238U_lowE_highT_140Xe[-1]))


################   238U lowE_highT bg_random -  140Xe   #################

#doublegate_bg_random
hist_doublegate_bg_random_238U_lowE_highT_140Xe = file_238U_lowE_highT.Get('time_isomer_doublegate_bg_random_140Xe')
x_bins = hist_doublegate_bg_random_238U_lowE_highT_140Xe.GetNbinsX()

x_doublegate_bg_random_238U_lowE_highT_140Xe = np.zeros(x_bins)
y_doublegate_bg_random_238U_lowE_highT_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_random_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_random_238U_lowE_highT_140Xe.GetBinCenter(i+1)
    y_doublegate_bg_random_238U_lowE_highT_140Xe[i] = hist_doublegate_bg_random_238U_lowE_highT_140Xe.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_random_238U_lowE_highT_140Xe = 2*x_doublegate_bg_random_238U_lowE_highT_140Xe

x_doublegate_bg_random_238U_lowE_highT_140Xe_long = x_doublegate_bg_random_238U_lowE_highT_140Xe
y_doublegate_bg_random_238U_lowE_highT_140Xe_long = y_doublegate_bg_random_238U_lowE_highT_140Xe

x_doublegate_bg_random_238U_lowE_highT_140Xe = x_doublegate_bg_random_238U_lowE_highT_140Xe[bin_lower:bin_upper]
y_doublegate_bg_random_238U_lowE_highT_140Xe = y_doublegate_bg_random_238U_lowE_highT_140Xe[bin_lower:bin_upper]

#print("* 238U - 140Xe fit range %d - %.d" % (x_doublegate_bg_random_238U_lowE_highT_140Xe[0],x_doublegate_bg_random_238U_lowE_highT_140Xe[-1]))



####################################################
##    			     Functions                    ## 
####################################################


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


def gauss(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    return amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))

def const_bg(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0):
	return const_bg + x*0


def sum_gauss_const_bg(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0):
    """ Sum of exponential decay and gaussian, with a constant background """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) 
        + const_bg)

def pdf_skewnormal(x):
    return 1.0/np.sqrt(2*np.pi) * np.exp(-x**2/2)

def cdf_skewnormal(x):
    return (1 + erf(x/np.sqrt(2))) / 2

def skewnormal(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, alpha = 1.0):
    t = (x-mean)/sigma
    return (amplitude_gauss*2.0/sigma * pdf_skewnormal(t) * cdf_skewnormal(alpha*t))

def sum_skewnormal_const_bg(x, mean=0, sigma=1.0, const_bg=1.0, amplitude_gauss=1.0, alpha = 1.0):
	t = (x-mean)/sigma
	return (amplitude_gauss*2.0/sigma * pdf_skewnormal(t) * cdf_skewnormal(alpha*t)) + const_bg



####################################################
## 		             Fit data          		      ##  START_FIT 
####################################################


################   238U lowE_highT -  140Xe - normal  #################

mean_lower = 300
mean_upper = 400
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 30000

P_double_238U_lowE_highT_140Xe, cov_double_238U_lowE_highT_140Xe = curve_fit(sum_gauss_const_bg, x_doublegate_238U_lowE_highT_140Xe, y_doublegate_238U_lowE_highT_140Xe, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_140Xe, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_140Xe, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_140Xe), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 140Xe Using uncertainty-weighted fit")

P_double_unc_238U_lowE_highT_140Xe = np.sqrt(np.diag(cov_double_238U_lowE_highT_140Xe))

print("\n")
print(" ***** 238U lowE_highT - 140Xe:  Doublegate true spectrum fit: ***** ")
print("                  -- GAUSS + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_highT_140Xe[0], P_double_unc_238U_lowE_highT_140Xe[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_140Xe[1], P_double_unc_238U_lowE_highT_140Xe[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_140Xe[2], P_double_unc_238U_lowE_highT_140Xe[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_140Xe[3], P_double_unc_238U_lowE_highT_140Xe[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("\n")



################   238U lowE_highT -  140Xe - skew  #################

mean_lower = 300
mean_upper = 400
sigma_lower = 0
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 1000
amplitude_gauss_lower = 0
amplitude_gauss_upper = 300000
alpha_lower = -100
alpha_upper = 100

P_double_skew_238U_lowE_highT_140Xe, cov_double_skew_238U_lowE_highT_140Xe = curve_fit(sum_skewnormal_const_bg, x_doublegate_238U_lowE_highT_140Xe, y_doublegate_238U_lowE_highT_140Xe, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_140Xe, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_140Xe, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_140Xe), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower, alpha_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper, alpha_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 140Xe Using uncertainty-weighted fit")

P_double_skew_unc_238U_lowE_highT_140Xe = np.sqrt(np.diag(cov_double_skew_238U_lowE_highT_140Xe))

print("\n")
print(" ***** 238U lowE_highT - 140Xe:  Doublegate true spectrum fit: ***** ")
print("                  -- SKEW GAUSS + CONST_BG FIT --   ")
print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_skew_238U_lowE_highT_140Xe[0], P_double_skew_unc_238U_lowE_highT_140Xe[0], mean_lower, mean_upper))
print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_skew_238U_lowE_highT_140Xe[1], P_double_skew_unc_238U_lowE_highT_140Xe[1], sigma_lower, sigma_upper))
print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_skew_238U_lowE_highT_140Xe[2], P_double_skew_unc_238U_lowE_highT_140Xe[2], const_bg_lower, const_bg_upper))
print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_skew_238U_lowE_highT_140Xe[3], P_double_skew_unc_238U_lowE_highT_140Xe[3], amplitude_gauss_lower, amplitude_gauss_upper))
print("alpha:          			 %.2f +/- %.2f         [%.d,%.d]" % (P_double_skew_238U_lowE_highT_140Xe[4], P_double_skew_unc_238U_lowE_highT_140Xe[4], alpha_lower, alpha_upper))

print("\n")





####################################################
###     Calculate chi-squared values for fits    ###
####################################################

chisquare_double_238U_lowE_highT_140Xe = reduced_chisquare_func(f_obs=y_doublegate_238U_lowE_highT_140Xe, f_exp=sum_gauss_const_bg(x_doublegate_238U_lowE_highT_140Xe, P_double_238U_lowE_highT_140Xe[0], P_double_238U_lowE_highT_140Xe[1], P_double_238U_lowE_highT_140Xe[2], P_double_238U_lowE_highT_140Xe[3]), N=len(y_doublegate_238U_lowE_highT_140Xe))

print(" Red chi**2: %.2f \n" % chisquare_double_238U_lowE_highT_140Xe)


####################################################
## 						Plot 		              ## START_PLOT
####################################################

x_array_plot = np.linspace(0,1000,10000)

# plt.plot(x_array_plot, skewnormal(x_array_plot,mean=100,sigma=5,alpha=2))
# plt.grid()
# plt.show()


# ################   238U lowE_highT -  140Xe   #################

plt.plot(x_doublegate_238U_lowE_highT_140Xe_long, y_doublegate_238U_lowE_highT_140Xe_long, label="doublegate_238U_lowE_highT_140Xe", color="royalblue")
#plt.errorbar(x_doublegate_238U_lowE_highT_140Xe_long, y_doublegate_238U_lowE_highT_140Xe_long, yerr=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_140Xe_long, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_140Xe_long, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_140Xe_long), label="doublegate_238U_lowE_highT_140Xe", color="royalblue")


plt.plot(x_array_plot, sum_gauss_const_bg(x_array_plot, P_double_238U_lowE_highT_140Xe[0], P_double_238U_lowE_highT_140Xe[1], P_double_238U_lowE_highT_140Xe[2], P_double_238U_lowE_highT_140Xe[3]), label="sum gauss + const bg", color="black")
#plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_highT_140Xe[0], P_double_238U_lowE_highT_140Xe[1], P_double_238U_lowE_highT_140Xe[2], P_double_238U_lowE_highT_140Xe[3]), label="true gaussian", color="grey")
#plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_highT_140Xe[0], P_double_238U_lowE_highT_140Xe[1], P_double_238U_lowE_highT_140Xe[2], P_double_238U_lowE_highT_140Xe[3]), label="constant BG", color="grey")

plt.plot(x_array_plot, sum_skewnormal_const_bg(x_array_plot, P_double_skew_238U_lowE_highT_140Xe[0], P_double_skew_238U_lowE_highT_140Xe[1], P_double_skew_238U_lowE_highT_140Xe[2], P_double_skew_238U_lowE_highT_140Xe[3], P_double_skew_238U_lowE_highT_140Xe[4]), label="sum skew gauss + const bg", color="red")
#plt.plot(x_array_plot, skewnormal(x_array_plot, P_double_skew_238U_lowE_highT_140Xe[0], P_double_skew_238U_lowE_highT_140Xe[1], P_double_skew_238U_lowE_highT_140Xe[2], P_double_skew_238U_lowE_highT_140Xe[3], P_double_skew_238U_lowE_highT_140Xe[4]), label="skew gaussian", color="orange")
#plt.plot(x_array_plot, const_bg(x_array_plot, P_double_skew_238U_lowE_highT_140Xe[0], P_double_skew_238U_lowE_highT_140Xe[1], P_double_skew_238U_lowE_highT_140Xe[2], P_double_skew_238U_lowE_highT_140Xe[3]), label="constant BG", color="orange")

plt.vlines(x_doublegate_238U_lowE_highT_140Xe[0],0,6000, label="fit range", color="black")
plt.vlines(x_doublegate_238U_lowE_highT_140Xe[-1],0,6000, color="black")
#plt.yscale("log")
plt.title("238U lowE_highT - 140Xe: Doublegate true spectrum fit")
#plt.axis([0,700,1,10**(4)])
#plt.axis([0,650,-50,4*10**(3)])
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.legend(fontsize=12)
plt.grid()
plt.show()



