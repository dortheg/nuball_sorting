import numpy as np
import ROOT 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d
import time


##########################
##     Read in data     ## 
##########################

#file = ROOT.TFile.Open("252Cf_24nov2021.root"," READ ")
#file = ROOT.TFile.Open("252Cf_25nov2021.root"," READ ")
file = ROOT.TFile.Open("252Cf_30nov2021.root"," READ ")

#####################
###     134Te      ##
#####################

#Define lower and upper fit limit
x_lower = (1000 + 950)//2
x_upper = (1000 + 1600)//2

#Isomer_1 gate (297keV) true
hist_isomer_1_gate_134Te = file.Get('time_isomer_1_gate_134Te')
x_bins = hist_isomer_1_gate_134Te.GetNbinsX()

x_isomer_1_gate_134Te = np.zeros(x_bins)
y_isomer_1_gate_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_134Te[i] = hist_isomer_1_gate_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_134Te[i] = hist_isomer_1_gate_134Te.GetBinContent(i+1)

#x_isomer_1_gate_134Te = x_isomer_1_gate_134Te[x_lower:x_upper] #Slice to fit limits
#y_isomer_1_gate_134Te = y_isomer_1_gate_134Te[x_lower:x_upper]


#Isomer_1 gate (297keV) all
hist_isomer_1_gate_all_134Te = file.Get('time_isomer_1_gate_all_134Te')
x_bins = hist_isomer_1_gate_all_134Te.GetNbinsX()

x_isomer_1_gate_all_134Te = np.zeros(x_bins)
y_isomer_1_gate_all_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_all_134Te[i] = hist_isomer_1_gate_all_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_all_134Te[i] = hist_isomer_1_gate_all_134Te.GetBinContent(i+1)

#x_isomer_1_gate_all_134Te = x_isomer_1_gate_all_134Te[x_lower:x_upper]
#y_isomer_1_gate_all_134Te = y_isomer_1_gate_all_134Te[x_lower:x_upper]


#Isomer_1 gate (297keV) bg
hist_isomer_1_gate_bg_134Te = file.Get('time_isomer_1_gate_bg_134Te')
x_bins = hist_isomer_1_gate_bg_134Te.GetNbinsX()

x_isomer_1_gate_bg_134Te = np.zeros(x_bins)
y_isomer_1_gate_bg_134Te = np.zeros(x_bins)

for i in range(x_bins):
    x_isomer_1_gate_bg_134Te[i] = hist_isomer_1_gate_bg_134Te.GetBinCenter(i+1)
    y_isomer_1_gate_bg_134Te[i] = hist_isomer_1_gate_bg_134Te.GetBinContent(i+1)

#x_isomer_1_gate_bg_134Te = x_isomer_1_gate_bg_134Te[x_lower:x_upper]
#y_isomer_1_gate_bg_134Te = y_isomer_1_gate_bg_134Te[x_lower:x_upper]

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

#Doublegate true
hist_doublegate_134Te = file.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_134Te.GetNbinsX()

x_doublegate_134Te = np.zeros(x_bins)
y_doublegate_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_134Te[i] = hist_doublegate_134Te.GetBinCenter(i+1)
    y_doublegate_134Te[i] = hist_doublegate_134Te.GetBinContent(i+1)

# x_doublegate_134Te = x_doublegate_134Te[x_lower:x_upper]
# y_doublegate_134Te = y_doublegate_134Te[x_lower:x_upper]


#Doublegate all
hist_doublegate_all_134Te = file.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_134Te.GetNbinsX()

x_doublegate_all_134Te = np.zeros(x_bins)
y_doublegate_all_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinCenter(i+1)
    y_doublegate_all_134Te[i] = hist_doublegate_all_134Te.GetBinContent(i+1)

# x_doublegate_all_134Te = x_doublegate_all_134Te[x_lower:x_upper]
# y_doublegate_all_134Te = y_doublegate_all_134Te[x_lower:x_upper]


#Doublegate bg
hist_doublegate_bg_134Te = file.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_134Te.GetNbinsX()

x_doublegate_bg_134Te = np.zeros(x_bins)
y_doublegate_bg_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinCenter(i+1)
    y_doublegate_bg_134Te[i] = hist_doublegate_bg_134Te.GetBinContent(i+1)

# x_doublegate_bg_134Te = x_doublegate_bg_134Te[x_lower:x_upper]
# y_doublegate_bg_134Te = y_doublegate_bg_134Te[x_lower:x_upper]


#####################
###     140Xe      ##
#####################

#Doublegate true
hist_doublegate_140Xe = file.Get('time_gamma_doublegate_140Xe')
x_bins = hist_doublegate_140Xe.GetNbinsX()

x_doublegate_140Xe = np.zeros(x_bins)
y_doublegate_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_140Xe[i] = hist_doublegate_140Xe.GetBinCenter(i+1)
    y_doublegate_140Xe[i] = hist_doublegate_140Xe.GetBinContent(i+1)

#Doublegate all
hist_doublegate_all_140Xe = file.Get('time_gamma_doublegate_all_140Xe')
x_bins = hist_doublegate_all_140Xe.GetNbinsX()

x_doublegate_all_140Xe = np.zeros(x_bins)
y_doublegate_all_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_140Xe[i] = hist_doublegate_all_140Xe.GetBinCenter(i+1)
    y_doublegate_all_140Xe[i] = hist_doublegate_all_140Xe.GetBinContent(i+1)

#Doublegate bg
hist_doublegate_bg_140Xe = file.Get('time_gamma_doublegate_bg_140Xe')
x_bins = hist_doublegate_bg_140Xe.GetNbinsX()

x_doublegate_bg_140Xe = np.zeros(x_bins)
y_doublegate_bg_140Xe = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_140Xe[i] = hist_doublegate_bg_140Xe.GetBinCenter(i+1)
    y_doublegate_bg_140Xe[i] = hist_doublegate_bg_140Xe.GetBinContent(i+1)



###################################
##      Define fitting func      ## 
###################################

#For 134Te
tau_134Te = 164.1/np.log(2)
sigma_tau_134Te = 0.9/np.log(2)


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


def sum_smeared_exp_gauss(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
	""" Sum of exponential decay and gaussian """
	return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2)) 
        + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma) )


def sum_two_smeared_exp_gauss(x, mean=0, sigma=1.0, amplitude_gauss=1.0, amplitude_exp_Ge=1.0, tau_Ge=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Sum of exponential decay and gaussian """
    return ( amplitude_gauss*np.exp(-(x-mean)**2/(2*sigma**2))
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_Ge*np.exp((mean-x)/tau_Ge)]),sigma)
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma) )


func = sum_smeared_exp_gauss
#print(str(func))

##########################
##     IYR functions    ## 
##########################

def IYR(prompt, delayed):
    return (delayed)/(2*prompt + delayed)

def sigma_IYR(prompt, delayed, all_prompt, all_delayed, bg_prompt, bg_delayed):
    sigma_prompt = np.sqrt(all_prompt + bg_prompt + (0.05*bg_prompt)**2)
    sigma_delayed = np.sqrt(all_delayed + bg_delayed + (0.05*bg_delayed)**2)
    #sigma_prompt = np.sqrt(all_prompt + bg_prompt)
    #sigma_delayed = np.sqrt(all_delayed + bg_delayed)
    return np.sqrt( ((2*prompt/(2*prompt+delayed)**2)*sigma_delayed)**2 + ((-2*delayed/(2*prompt+delayed)**2)*sigma_prompt)**2 )

def sigma_data(data_all, data_bg):
    return np.sqrt(data_all + data_bg + (0.05*data_bg)**2)


####################################################
## 		             Fit data 		              ## 
####################################################

#Remember: will need some help to find correct fit parameters

##################################################################

#########################
#134Te isomer_1 gate fits
#########################
#P_isomer_1, cov_isomer_1 = curve_fit(func, x_isomer_1_gate_134Te, y_isomer_1_gate_134Te, bounds=([950,0,20,10,tau_134Te-2*sigma_tau_134Te],[1020,20,10000,10000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_1 true spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_1[0])
# print("sigma: %.4f" % P_isomer_1[1])
# print("amplitude_gauss: %.4f" % P_isomer_1[2])
# print("amplitude_exp: %.4f" % P_isomer_1[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_1[4],P_isomer_1[4]*np.log(2)))


#P_isomer_1_all, cov_isomer_1_all = curve_fit(func, x_isomer_1_gate_all_134Te, y_isomer_1_gate_all_134Te, bounds=([990,0,140000,0,0],[1000,30,170000,18000,50]))

#P_isomer_1_all, cov_isomer_1_all = curve_fit(func, x_isomer_1_gate_all_134Te, y_isomer_1_gate_all_134Te, bounds=([990,10,120000,0,tau_134Te-2*sigma_tau_134Te],[1000,20,170000,30000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_1 all spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_1_all[0])
# print("sigma: %.4f" % P_isomer_1_all[1])
# print("amplitude_gauss: %.4f" % P_isomer_1_all[2])
# print("amplitude_exp: %.4f" % P_isomer_1_all[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_1_all[4],P_isomer_1_all[4]*np.log(2)))

#P_isomer_1_bg, cov_isomer_1_bg = curve_fit(func, x_isomer_1_gate_bg_134Te, y_isomer_1_gate_bg_134Te, bounds=([990,0,140000,0,tau_134Te-2*sigma_tau_134Te],[1000,30,170000,18000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_1 BG spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_1_bg[0])
# print("sigma: %.4f" % P_isomer_1_bg[1])
# print("amplitude_gauss: %.4f" % P_isomer_1_bg[2])
# print("amplitude_exp: %.4f" % P_isomer_1_bg[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_1_bg[4],P_isomer_1_bg[4]*np.log(2)))

######################
#134Te isomer_2 gate fits
######################
#P_isomer_2, cov_isomer_2 = curve_fit(func, x_isomer_2_gate_134Te, y_isomer_2_gate_134Te, bounds=([950,5,20,10,tau_134Te-2*sigma_tau_134Te],[1020,20,4000,4000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_2 true spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_2[0])
# print("sigma: %.4f" % P_isomer_2[1])
# print("amplitude_gauss: %.4f" % P_isomer_2[2])
# print("amplitude_exp: %.4f" % P_isomer_2[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_2[4],P_isomer_2[4]*np.log(2)))

#P_isomer_2_all, cov_isomer_2_all = curve_fit(func, x_isomer_2_gate_all_134Te, y_isomer_2_gate_all_134Te, bounds=([990,0,0,0,tau_134Te-2*sigma_tau_134Te],[1000,30,170000,10000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_2 all spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_2_all[0])
# print("sigma: %.4f" % P_isomer_2_all[1])
# print("amplitude_gauss: %.4f" % P_isomer_2_all[2])
# print("amplitude_exp: %.4f" % P_isomer_2_all[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_2_all[4],P_isomer_2_all[4]*np.log(2)))

#P_isomer_2_bg, cov_isomer_2_bg = curve_fit(func, x_isomer_2_gate_bg_134Te, y_isomer_2_gate_bg_134Te, bounds=([995,0,0,0,tau_134Te-2*sigma_tau_134Te],[1000,30,170000,10000,tau_134Te+2*sigma_tau_134Te]))
# print("\n")
# print(" *****  Isomer_2 BG spectrum fit ***** \n")
# print("mean: %.4f" % P_isomer_2_bg[0])
# print("sigma: %.4f" % P_isomer_2_bg[1])
# print("amplitude_gauss: %.4f" % P_isomer_2_bg[2])
# print("amplitude_exp: %.4f" % P_isomer_2_bg[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_isomer_2_bg[4],P_isomer_2_bg[4]*np.log(2)))

######################
#134Te doublegate fits
######################
#mean, sigma, amplitude_gauss, amplitude_exp_Ge, tau_Ge, amplitude_exp_decay, tau_decay
#P_double, cov_double = curve_fit(func, x_doublegate_134Te, y_doublegate_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50])) #24nov
#P_double, cov_double = curve_fit(func, x_doublegate_134Te, y_doublegate_134Te, bounds=([950,0,0,0,tau_134Te-2*sigma_tau_134Te],[1100,40,3000,1000,tau_134Te+2*sigma_tau_134Te])) #25nov
P_double, cov_double = curve_fit(sum_two_smeared_exp_gauss, x_doublegate_134Te, y_doublegate_134Te, bounds=([950,0,0,0,0,0,0],[1100,40,3000,1000,40,1000,1000])) #30nov
print("\n")
print(" ***** 134Te:  Doublegated true spectrum fit ***** \n")
print("mean: %.4f" % P_double[0])
print("sigma: %.4f" % P_double[1])
print("amplitude_gauss: %.4f" % P_double[2])
print("amplitude_exp_Ge: %.4f" % P_double[3])
print("tau_Ge: %.4f" % P_double[4])
print("amplitude_exp_decay: %.4f" % P_double[5])
print("tau_decay: %.4f,  in half_life:  %.4f" % (P_double[6],P_double[6]*np.log(2)))

#P_double_all, cov_double_all = curve_fit(func, x_doublegate_all_134Te, y_doublegate_all_134Te, bounds=([0,950,0,20,10,0],[1000,1100,40,300,200,50])) #24nov
#P_double_all, cov_double_all = curve_fit(func, x_doublegate_all_134Te, y_doublegate_all_134Te, bounds=([950,0,20,10,tau_134Te-2*sigma_tau_134Te],[1100,40,300,200,tau_134Te+2*sigma_tau_134Te])) #25nov
# print("\n")
# print(" *****  Doublegated all spectrum fit ***** \n")
# print("mean: %.4f" % P_double_all[0])
# print("sigma: %.4f" % P_double_all[1])
# print("amplitude_gauss: %.4f" % P_double_all[2])
# print("amplitude_exp: %.4f" % P_double_all[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_double_all[4],P_double_all[4]*np.log(2)))


#P_double_bg, cov_double_bg = curve_fit(func, x_doublegate_bg_134Te, y_doublegate_bg_134Te, bounds=([0,950,0,0,0,0],[1000,1100,40,40,40,50])) #24.nov
#P_double_bg, cov_double_bg = curve_fit(func, x_doublegate_bg_134Te, y_doublegate_bg_134Te, bounds=([950,0,0,0,tau_134Te-2*sigma_tau_134Te],[1100,20,40,40,tau_134Te+2*sigma_tau_134Te])) #25.nov
# print("\n")
# print(" *****  Doublegated BG spectrum fit ***** \n")
# print("mean: %.4f" % P_double_bg[0])
# print("sigma: %.4f" % P_double_bg[1])
# print("amplitude_gauss: %.4f" % P_double_bg[2])
# print("amplitude_exp: %.4f" % P_double_bg[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_double_bg[4],P_double_bg[4]*np.log(2)))

########################################################################################

######################
#140Xe doublegate fits
######################
#mean, sigma, amplitude_gauss, amplitude_exp, tau
#P_double_140Xe, cov_double_140Xe = curve_fit(func, x_doublegate_140Xe, y_doublegate_140Xe, bounds=([950,0,0,0,0],[1100,40,3000,1000,1000]))
# print("\n")
# print(" *****  140Xe: Doublegated true spectrum fit ***** \n")
# print("mean: %.4f" % P_double_140Xe[0])
# print("sigma: %.4f" % P_double_140Xe[1])
# print("amplitude_gauss: %.4f" % P_double_140Xe[2])
# print("amplitude_exp: %.4f" % P_double_140Xe[3])
# print("tau: %.4f,  in half_life:  %.4f" % (P_double_140Xe[4],P_double_140Xe[4]*np.log(2)))




######################################
## 		 Find IYR + uncertainty 	## 
######################################

x_arr = np.linspace(-1000,3000,3000)

#Isomer_1_gate
# area_isomer_1_gate_true = np.trapz(func(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4]), x_arr)
# area_isomer_1_gate_true_prompt = np.trapz(gauss(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4]), x_arr)
# area_isomer_1_gate_true_delayed = np.trapz(smeared_exp(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4]), x_arr)

# area_isomer_1_gate_all = np.trapz(func(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), x_arr)
# area_isomer_1_gate_all_prompt = np.trapz(gauss(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), x_arr)
# area_isomer_1_gate_all_delayed = np.trapz(smeared_exp(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), x_arr)

# area_isomer_1_gate_bg = np.trapz(func(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), x_arr)
# area_isomer_1_gate_bg_prompt = np.trapz(gauss(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), x_arr)
# area_isomer_1_gate_bg_delayed = np.trapz(smeared_exp(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), x_arr)

# IYR_isomer_1_gate = IYR(prompt=area_isomer_1_gate_true_prompt, delayed=area_isomer_1_gate_true_delayed)
# sigma_IYR_isomer_1_gate = sigma_IYR(prompt=area_isomer_1_gate_true_prompt, delayed=area_isomer_1_gate_true_delayed, all_prompt=area_isomer_1_gate_all_prompt, all_delayed=area_isomer_1_gate_all_delayed, bg_prompt=area_isomer_1_gate_bg_prompt, bg_delayed=area_isomer_1_gate_bg_delayed)

# print("\n")
# print("Area isomer_1 true tot: %.3f" % area_isomer_1_gate_true)
# print("Area isomer_1 true prompt: %.3f" % area_isomer_1_gate_true_prompt)
# print("Area isomer_1 true delayed: %.3f" % area_isomer_1_gate_true_delayed)
# area_isomer_1_data_true = np.trapz(y_isomer_1_gate_134Te, x_isomer_1_gate_134Te)
# print("Area isomer_1 true tot data: %.3f" % area_isomer_1_data_true)
# rel_fit_isomer_1 = abs(area_isomer_1_data_true-area_isomer_1_gate_true)/(area_isomer_1_data_true)*100
# print("Percent difference between data and fit: %.3f percent" % rel_fit_isomer_1)

#Isomer_2_gate
# area_isomer_2_gate_true = np.trapz(func(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), x_arr)
# area_isomer_2_gate_true_prompt = np.trapz(gauss(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), x_arr)
# area_isomer_2_gate_true_delayed = np.trapz(smeared_exp(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), x_arr)

# area_isomer_2_gate_all = np.trapz(func(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), x_arr)
# area_isomer_2_gate_all_prompt = np.trapz(gauss(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), x_arr)
# area_isomer_2_gate_all_delayed = np.trapz(smeared_exp(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), x_arr)

# area_isomer_2_gate_bg = np.trapz(func(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), x_arr)
# area_isomer_2_gate_bg_prompt = np.trapz(gauss(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), x_arr)
# area_isomer_2_gate_bg_delayed = np.trapz(smeared_exp(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), x_arr)

# IYR_isomer_2_gate = IYR(prompt=area_isomer_2_gate_true_prompt, delayed=area_isomer_2_gate_true_delayed)
# sigma_IYR_isomer_2_gate = sigma_IYR(prompt=area_isomer_2_gate_true_prompt, delayed=area_isomer_2_gate_true_delayed, all_prompt=area_isomer_2_gate_all_prompt, all_delayed=area_isomer_2_gate_all_delayed, bg_prompt=area_isomer_2_gate_bg_prompt, bg_delayed=area_isomer_2_gate_bg_delayed)

# print("\n")
# print("Area isomer_2 true tot: %.3f" % area_isomer_2_gate_true)
# print("Area isomer_2 true prompt: %.3f" % area_isomer_2_gate_true_prompt)
# print("Area isomer_2 true delayed: %.3f" % area_isomer_2_gate_true_delayed)
# area_isomer_2_data_true = np.trapz(y_isomer_2_gate_134Te, x_isomer_2_gate_134Te)
# print("Area isomer_2 true tot data: %.3f" % area_isomer_2_data_true)
# rel_fit_isomer_2 = abs(area_isomer_2_data_true-area_isomer_2_gate_true)/(area_isomer_2_data_true)*100
# print("Percent difference between data and fit: %.3f percent" % rel_fit_isomer_2)



#Doublegated
# area_double_true = np.trapz(func(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4]), x_arr)
# area_double_true_prompt = np.trapz(gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4]), x_arr)
# area_double_true_delayed = np.trapz(smeared_exp(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4]), x_arr)

# area_double_all = np.trapz(func(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), x_arr)
# area_double_all_prompt = np.trapz(gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), x_arr)
# area_double_all_delayed = np.trapz(smeared_exp(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), x_arr)

# area_double_bg = np.trapz(func(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), x_arr)
# area_double_bg_prompt = np.trapz(gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), x_arr)
# area_double_bg_delayed = np.trapz(smeared_exp(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), x_arr)

# IYR_double = IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed)
# sigma_IYR_double = sigma_IYR(prompt=area_double_true_prompt, delayed=area_double_true_delayed, all_prompt=area_double_all_prompt, all_delayed=area_double_all_delayed, bg_prompt=area_double_bg_prompt, bg_delayed=area_double_bg_delayed) #Make this into a function

#print("\n")
# print("Area double true tot: %.3f" % area_double_true)
# print("Area double true prompt: %.3f" % area_double_true_prompt)
# print("Area double true delayed: %.3f" % area_double_true_delayed)
#area_double_data_true = np.trapz(y_doublegate_134Te, x_doublegate_134Te)
#print("Area double true tot data: %.3f" % area_double_data_true)
#area_double_true = np.trapz(func(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4]), x_arr)
#rel_fit_double = abs(area_double_data_true-area_double_true)/(area_double_data_true)*100
#print("Percent difference between data and fit: %.3f percent" % rel_fit_double)


# print("\n")
# print(" ***** Isomeric Yield Ratios ****")
# print("IYR_isomer_1 (297 keV):   %.4f +/- %.4f" % (IYR_isomer_1_gate, sigma_IYR_isomer_1_gate) )
# print("IYR_isomer_2 (1279 keV):  %.4f +/- %.4f" % (IYR_isomer_2_gate, sigma_IYR_isomer_2_gate) )
# print("IYR_double:               %.4f +/- %.4f " % (IYR_double, sigma_IYR_double) )
# print("\n")


##########################
## 			Plot 		## 
##########################

#134Te isomer 1 gated, fit
#plt.plot(x_isomer_1_gate_134Te, y_isomer_1_gate_134Te, label="isomer_1_gate_134Te", color="royalblue")
#plt.plot(x_isomer_1_gate_all_134Te, y_isomer_1_gate_all_134Te, label="isomer_1_gate_all_134Te", color="black")
#plt.plot(x_isomer_1_gate_bg_134Te, y_isomer_1_gate_bg_134Te, label="isomer_1_gate_bg_134Te", color="pink")

# plt.plot(x_arr, func(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3]), P_isomer_1[4], label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1[0], P_isomer_1[1], P_isomer_1[2], P_isomer_1[3], P_isomer_1[4]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1_all[0], P_isomer_1_all[1], P_isomer_1_all[2], P_isomer_1_all[3], P_isomer_1_all[4]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_1_bg[0], P_isomer_1_bg[1], P_isomer_1_bg[2], P_isomer_1_bg[3], P_isomer_1_bg[4]), label="true smeared exp", color="red")

# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()

#134Te isomer 2 gated, fit
#plt.plot(x_isomer_2_gate_134Te, y_isomer_2_gate_134Te, label="isomer_2_gate_134Te", color="royalblue")
#plt.plot(x_isomer_2_gate_all_134Te, y_isomer_2_gate_all_134Te, label="isomer_2_gate_all_134Te", color="black")
#plt.plot(x_isomer_2_gate_bg_134Te, y_isomer_2_gate_bg_134Te, label="isomer_2_gate_bg_134Te", color="pink")

# plt.plot(x_arr, func(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_2[0], P_isomer_2[1], P_isomer_2[2], P_isomer_2[3], P_isomer_2[4]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_2_all[0], P_isomer_2_all[1], P_isomer_2_all[2], P_isomer_2_all[3], P_isomer_2_all[4]), label="true smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_isomer_2_bg[0], P_isomer_2_bg[1], P_isomer_2_bg[2], P_isomer_2_bg[3], P_isomer_2_bg[4]), label="true smeared exp", color="red")

# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


#134Te doublegated, fit
plt.errorbar(x_doublegate_134Te, y_doublegate_134Te, yerr=sigma_data(y_doublegate_all_134Te, y_doublegate_bg_134Te), fmt=".", label="doublegate_134Te", color="royalblue")
#plt.plot(x_doublegate_134Te, y_doublegate_134Te, label="doublegate_134Te", color="royalblue")
#plt.plot(x_doublegate_all_134Te, y_doublegate_all_134Te, label="doublegate_all_134Te", color="black")
#plt.plot(x_doublegate_bg_134Te, y_doublegate_bg_134Te, label="doublegate_bg_134Te", color="pink")

plt.plot(x_arr, sum_two_smeared_exp_gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6]), label="true fit, total", color="orange")
#plt.plot(x_arr, gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6]), label="true gaussian", color="green")
#plt.plot(x_arr, smeared_exp_Ge(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6]), label="true smeared exp Ge", color="lime")
plt.plot(x_arr, smeared_exp_decay(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6]), label="true smeared exp decay", color="red")
plt.plot(x_arr, gauss(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6])+smeared_exp_Ge(x_arr, P_double[0], P_double[1], P_double[2], P_double[3], P_double[4], P_double[5], P_double[6]), label="sum prompt", color="darkolivegreen")

# plt.plot(x_arr, func(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), label="all fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), label="all gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double_all[0], P_double_all[1], P_double_all[2], P_double_all[3], P_double_all[4]), label="all smeared exp", color="red")

# plt.plot(x_arr, func(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), label="bg fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), label="bg gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double_bg[0], P_double_bg[1], P_double_bg[2], P_double_bg[3], P_double_bg[4]), label="bg smeared exp", color="red")

plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.title(str(func))
plt.axis([800,2000,-10,250])
plt.legend(fontsize=14)
plt.grid()
plt.show()


#140Xe doublegated
#plt.plot(x_doublegate_140Xe, y_doublegate_140Xe, label="doublegate_140Xe", color="royalblue")
#plt.plot(x_doublegate_all_140Xe, y_doublegate_all_140Xe, label="doublegate_all_140Xe", color="black")
#plt.plot(x_doublegate_bg_140Xe, y_doublegate_bg_140Xe, label="doublegate_bg_140Xe", color="pink")

# plt.plot(x_arr, func(x_arr, P_double_140Xe[0], P_double_140Xe[1], P_double_140Xe[2], P_double_140Xe[3], P_double_140Xe[4]), label="true fit, total", color="orange")
# plt.plot(x_arr, gauss(x_arr, P_double_140Xe[0], P_double_140Xe[1], P_double_140Xe[2], P_double_140Xe[3], P_double_140Xe[4]), label="true gaussian", color="green")
# plt.plot(x_arr, smeared_exp(x_arr, P_double_140Xe[0], P_double_140Xe[1], P_double_140Xe[2], P_double_140Xe[3], P_double_140Xe[4]), label="true smeared exp", color="red")

# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.title(str(func))
# plt.axis([800,1600,-10,1500])
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


#Plot the IYR
# plt.errorbar([1], [IYR_isomer_1_gate], yerr=[sigma_IYR_isomer_1_gate], fmt="ro", label="isomer_1_gate")
# plt.errorbar([2], [IYR_isomer_2_gate], yerr=[sigma_IYR_isomer_2_gate], fmt="bo", label="isomer_2_gate")
# plt.errorbar([3], [IYR_double], yerr=[sigma_IYR_double], fmt="go", label="doublegate")
# plt.axis([0,4,0,1.2])
# plt.xlabel("Case", fontsize=14)
# plt.ylabel("IYR (uncorrected)", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


