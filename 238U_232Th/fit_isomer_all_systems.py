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

BOOTSTRAP = False 

####################################################
##     				Read in data                  ## START_DATA
####################################################

file_238U_lowE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.bin_final.root"," READ ")
file_238U_lowE_lowT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_2ns_lowE_lowT_7feb2022.bin.root"," READ ")

file_238U_highE_highT = ROOT.TFile.Open("Sorted_files/238Ucube_hit4_2ns_highE_highT_19jan2022.bin.root"," READ ")

file_232Th = ROOT.TFile.Open("Sorted_files/232Thcube_hit4_lowE_highT_28feb2022.bin_final.root"," READ ")

print("\n")
print("238U lowE_highT file: ", file_238U_lowE_highT)
print("238U highE file: ", file_238U_highE_highT)
print("232Th file: ", file_232Th)


print("\n")
print("FIT CHOICES")

#Define lower and upper fit limit
x_lower = 330
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2


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

x_doublegate_238U_lowE_highT_134Te_long = x_doublegate_238U_lowE_highT_134Te
y_doublegate_238U_lowE_highT_134Te_long = y_doublegate_238U_lowE_highT_134Te

x_doublegate_238U_lowE_highT_134Te = x_doublegate_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_238U_lowE_highT_134Te = y_doublegate_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_238U_lowE_highT_134Te[0],x_doublegate_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT all -  134Te   #################

#doublegate_all true
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

#doublegate_bg true
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



################   238U lowE_highT -  134Te 1n  #################

#doublegate_1n true
hist_doublegate_1n_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_1n_134Te')
x_bins = hist_doublegate_1n_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_1n_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_1n_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_238U_lowE_highT_134Te[i] = hist_doublegate_1n_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_1n_238U_lowE_highT_134Te[i] = hist_doublegate_1n_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_238U_lowE_highT_134Te = 2*x_doublegate_1n_238U_lowE_highT_134Te

x_doublegate_1n_238U_lowE_highT_134Te_long = x_doublegate_1n_238U_lowE_highT_134Te
y_doublegate_1n_238U_lowE_highT_134Te_long = y_doublegate_1n_238U_lowE_highT_134Te

x_doublegate_1n_238U_lowE_highT_134Te = x_doublegate_1n_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_1n_238U_lowE_highT_134Te = y_doublegate_1n_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 1n fit range %d - %.d" % (x_doublegate_1n_238U_lowE_highT_134Te[0],x_doublegate_1n_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT all -  134Te 1n  #################

#doublegate_1n_all true
hist_doublegate_1n_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_1n_all_134Te')
x_bins = hist_doublegate_1n_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_1n_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_1n_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_all_238U_lowE_highT_134Te[i] = hist_doublegate_1n_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_1n_all_238U_lowE_highT_134Te[i] = hist_doublegate_1n_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_all_238U_lowE_highT_134Te = 2*x_doublegate_1n_all_238U_lowE_highT_134Te

x_doublegate_1n_all_238U_lowE_highT_134Te_long = x_doublegate_1n_all_238U_lowE_highT_134Te
y_doublegate_1n_all_238U_lowE_highT_134Te_long = y_doublegate_1n_all_238U_lowE_highT_134Te

x_doublegate_1n_all_238U_lowE_highT_134Te = x_doublegate_1n_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_1n_all_238U_lowE_highT_134Te = y_doublegate_1n_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 1n_all fit range %d - %.d" % (x_doublegate_1n_all_238U_lowE_highT_134Te[0],x_doublegate_1n_all_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg -  134Te 1n  #################

#doublegate_1n_bg true
hist_doublegate_1n_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_1n_bg_134Te')
x_bins = hist_doublegate_1n_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_1n_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_238U_lowE_highT_134Te = 2*x_doublegate_1n_bg_238U_lowE_highT_134Te

x_doublegate_1n_bg_238U_lowE_highT_134Te_long = x_doublegate_1n_bg_238U_lowE_highT_134Te
y_doublegate_1n_bg_238U_lowE_highT_134Te_long = y_doublegate_1n_bg_238U_lowE_highT_134Te

x_doublegate_1n_bg_238U_lowE_highT_134Te = x_doublegate_1n_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_238U_lowE_highT_134Te = y_doublegate_1n_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 1n_bg fit range %d - %.d" % (x_doublegate_1n_bg_238U_lowE_highT_134Te[0],x_doublegate_1n_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te 1n  #################

#doublegate_1n_bg_ridge true
hist_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_1n_bg_ridge_134Te')
x_bins = hist_doublegate_1n_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te

x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te
y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te

x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te = y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 1n_bg_ridge fit range %d - %.d" % (x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_1n_bg_ridge_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_random -  134Te 1n  #################

#doublegate_1n_bg_random true
hist_doublegate_1n_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_1n_bg_random_134Te')
x_bins = hist_doublegate_1n_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_1n_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_1n_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_1n_bg_random_238U_lowE_highT_134Te

x_doublegate_1n_bg_random_238U_lowE_highT_134Te_long = x_doublegate_1n_bg_random_238U_lowE_highT_134Te
y_doublegate_1n_bg_random_238U_lowE_highT_134Te_long = y_doublegate_1n_bg_random_238U_lowE_highT_134Te

x_doublegate_1n_bg_random_238U_lowE_highT_134Te = x_doublegate_1n_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_random_238U_lowE_highT_134Te = y_doublegate_1n_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 1n_bg_random fit range %d - %.d" % (x_doublegate_1n_bg_random_238U_lowE_highT_134Te[0],x_doublegate_1n_bg_random_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT -  134Te 3n  #################

#doublegate_3n true
hist_doublegate_3n_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_134Te')
x_bins = hist_doublegate_3n_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_238U_lowE_highT_134Te[i] = hist_doublegate_3n_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_238U_lowE_highT_134Te[i] = hist_doublegate_3n_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_238U_lowE_highT_134Te = 2*x_doublegate_3n_238U_lowE_highT_134Te

x_doublegate_3n_238U_lowE_highT_134Te_long = x_doublegate_3n_238U_lowE_highT_134Te
y_doublegate_3n_238U_lowE_highT_134Te_long = y_doublegate_3n_238U_lowE_highT_134Te

x_doublegate_3n_238U_lowE_highT_134Te = x_doublegate_3n_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_238U_lowE_highT_134Te = y_doublegate_3n_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n fit range %d - %.d" % (x_doublegate_3n_238U_lowE_highT_134Te[0],x_doublegate_3n_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT all -  134Te 3n  #################

#doublegate_3n_all true
hist_doublegate_3n_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_all_134Te')
x_bins = hist_doublegate_3n_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_all_238U_lowE_highT_134Te = 2*x_doublegate_3n_all_238U_lowE_highT_134Te

x_doublegate_3n_all_238U_lowE_highT_134Te_long = x_doublegate_3n_all_238U_lowE_highT_134Te
y_doublegate_3n_all_238U_lowE_highT_134Te_long = y_doublegate_3n_all_238U_lowE_highT_134Te

x_doublegate_3n_all_238U_lowE_highT_134Te = x_doublegate_3n_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_all_238U_lowE_highT_134Te = y_doublegate_3n_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_all fit range %d - %.d" % (x_doublegate_3n_all_238U_lowE_highT_134Te[0],x_doublegate_3n_all_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT bg -  134Te 3n  #################

#doublegate_3n_bg true
hist_doublegate_3n_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_bg_134Te')
x_bins = hist_doublegate_3n_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_238U_lowE_highT_134Te = 2*x_doublegate_3n_bg_238U_lowE_highT_134Te

x_doublegate_3n_bg_238U_lowE_highT_134Te_long = x_doublegate_3n_bg_238U_lowE_highT_134Te
y_doublegate_3n_bg_238U_lowE_highT_134Te_long = y_doublegate_3n_bg_238U_lowE_highT_134Te

x_doublegate_3n_bg_238U_lowE_highT_134Te = x_doublegate_3n_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_238U_lowE_highT_134Te = y_doublegate_3n_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_bg fit range %d - %.d" % (x_doublegate_3n_bg_238U_lowE_highT_134Te[0],x_doublegate_3n_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te 3n  #################

#doublegate_3n_bg_ridge true
hist_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_bg_ridge_134Te')
x_bins = hist_doublegate_3n_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te
y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te = y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_bg_ridge fit range %d - %.d" % (x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_3n_bg_ridge_238U_lowE_highT_134Te[-1]))

################   238U lowE_highT bg_random -  134Te 3n  #################

#doublegate_3n_bg_random true
hist_doublegate_3n_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_bg_random_134Te')
x_bins = hist_doublegate_3n_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_3n_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_bg_random_238U_lowE_highT_134Te_long = x_doublegate_3n_bg_random_238U_lowE_highT_134Te
y_doublegate_3n_bg_random_238U_lowE_highT_134Te_long = y_doublegate_3n_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_bg_random_238U_lowE_highT_134Te = x_doublegate_3n_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_random_238U_lowE_highT_134Te = y_doublegate_3n_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_bg_random fit range %d - %.d" % (x_doublegate_3n_bg_random_238U_lowE_highT_134Te[0],x_doublegate_3n_bg_random_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT -  134Te 3n_2plus  #################

#doublegate_3n_2plus true
hist_doublegate_3n_2plus_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_2plus_134Te')
x_bins = hist_doublegate_3n_2plus_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_2plus_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_238U_lowE_highT_134Te = 2*x_doublegate_3n_2plus_238U_lowE_highT_134Te

x_doublegate_3n_2plus_238U_lowE_highT_134Te_long = x_doublegate_3n_2plus_238U_lowE_highT_134Te
y_doublegate_3n_2plus_238U_lowE_highT_134Te_long = y_doublegate_3n_2plus_238U_lowE_highT_134Te

x_doublegate_3n_2plus_238U_lowE_highT_134Te = x_doublegate_3n_2plus_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_238U_lowE_highT_134Te = y_doublegate_3n_2plus_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_2plus fit range %d - %.d" % (x_doublegate_3n_2plus_238U_lowE_highT_134Te[0],x_doublegate_3n_2plus_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT all -  134Te 3n_2plus  #################

#doublegate_3n_2plus_all true
hist_doublegate_3n_2plus_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_2plus_all_134Te')
x_bins = hist_doublegate_3n_2plus_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_2plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_all_238U_lowE_highT_134Te = 2*x_doublegate_3n_2plus_all_238U_lowE_highT_134Te

x_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long = x_doublegate_3n_2plus_all_238U_lowE_highT_134Te
y_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long = y_doublegate_3n_2plus_all_238U_lowE_highT_134Te

x_doublegate_3n_2plus_all_238U_lowE_highT_134Te = x_doublegate_3n_2plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_all_238U_lowE_highT_134Te = y_doublegate_3n_2plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_2plus_all fit range %d - %.d" % (x_doublegate_3n_2plus_all_238U_lowE_highT_134Te[0],x_doublegate_3n_2plus_all_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT bg -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg true
hist_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_2plus_bg_134Te')
x_bins = hist_doublegate_3n_2plus_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = 2*x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te_long = x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te
y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te_long = y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te = y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_2plus_bg fit range %d - %.d" % (x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[0],x_doublegate_3n_2plus_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg_ridge true
hist_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_2plus_bg_ridge_134Te')
x_bins = hist_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te
y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te = y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_2plus_bg_ridge fit range %d - %.d" % (x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_random -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg_random true
hist_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_2plus_bg_random_134Te')
x_bins = hist_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te_long = x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te
y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te_long = y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te = y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_2plus_bg_random fit range %d - %.d" % (x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[0],x_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT -  134Te 3n_4plus  #################

#doublegate_3n_4plus true
hist_doublegate_3n_4plus_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_4plus_134Te')
x_bins = hist_doublegate_3n_4plus_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_4plus_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_238U_lowE_highT_134Te = 2*x_doublegate_3n_4plus_238U_lowE_highT_134Te

x_doublegate_3n_4plus_238U_lowE_highT_134Te_long = x_doublegate_3n_4plus_238U_lowE_highT_134Te
y_doublegate_3n_4plus_238U_lowE_highT_134Te_long = y_doublegate_3n_4plus_238U_lowE_highT_134Te

x_doublegate_3n_4plus_238U_lowE_highT_134Te = x_doublegate_3n_4plus_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_238U_lowE_highT_134Te = y_doublegate_3n_4plus_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_4plus fit range %d - %.d" % (x_doublegate_3n_4plus_238U_lowE_highT_134Te[0],x_doublegate_3n_4plus_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT all -  134Te 3n_4plus  #################

#doublegate_3n_4plus_all true
hist_doublegate_3n_4plus_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_4plus_all_134Te')
x_bins = hist_doublegate_3n_4plus_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_4plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_all_238U_lowE_highT_134Te = 2*x_doublegate_3n_4plus_all_238U_lowE_highT_134Te

x_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long = x_doublegate_3n_4plus_all_238U_lowE_highT_134Te
y_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long = y_doublegate_3n_4plus_all_238U_lowE_highT_134Te

x_doublegate_3n_4plus_all_238U_lowE_highT_134Te = x_doublegate_3n_4plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_all_238U_lowE_highT_134Te = y_doublegate_3n_4plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_4plus_all fit range %d - %.d" % (x_doublegate_3n_4plus_all_238U_lowE_highT_134Te[0],x_doublegate_3n_4plus_all_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT bg -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg true
hist_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_4plus_bg_134Te')
x_bins = hist_doublegate_3n_4plus_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = 2*x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te_long = x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te
y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te_long = y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te = y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_4plus_bg fit range %d - %.d" % (x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[0],x_doublegate_3n_4plus_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg_ridge true
hist_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_4plus_bg_ridge_134Te')
x_bins = hist_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te
y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te = y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_4plus_bg_ridge fit range %d - %.d" % (x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_random -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg_random true
hist_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_4plus_bg_random_134Te')
x_bins = hist_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te_long = x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te
y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te_long = y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te = y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_4plus_bg_random fit range %d - %.d" % (x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[0],x_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT -  134Te 3n_6plus  #################

#doublegate_3n_6plus true
hist_doublegate_3n_6plus_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_6plus_134Te')
x_bins = hist_doublegate_3n_6plus_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_6plus_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_6plus_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_6plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_6plus_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_6plus_238U_lowE_highT_134Te = 2*x_doublegate_3n_6plus_238U_lowE_highT_134Te

x_doublegate_3n_6plus_238U_lowE_highT_134Te_long = x_doublegate_3n_6plus_238U_lowE_highT_134Te
y_doublegate_3n_6plus_238U_lowE_highT_134Te_long = y_doublegate_3n_6plus_238U_lowE_highT_134Te

x_doublegate_3n_6plus_238U_lowE_highT_134Te = x_doublegate_3n_6plus_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_6plus_238U_lowE_highT_134Te = y_doublegate_3n_6plus_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_6plus fit range %d - %.d" % (x_doublegate_3n_6plus_238U_lowE_highT_134Te[0],x_doublegate_3n_6plus_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT all -  134Te 3n_6plus  #################

#doublegate_3n_6plus_all true
hist_doublegate_3n_6plus_all_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_6plus_all_134Te')
x_bins = hist_doublegate_3n_6plus_all_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_6plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_6plus_all_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_6plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_all_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_6plus_all_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_all_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_6plus_all_238U_lowE_highT_134Te = 2*x_doublegate_3n_6plus_all_238U_lowE_highT_134Te

x_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long = x_doublegate_3n_6plus_all_238U_lowE_highT_134Te
y_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long = y_doublegate_3n_6plus_all_238U_lowE_highT_134Te

x_doublegate_3n_6plus_all_238U_lowE_highT_134Te = x_doublegate_3n_6plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_6plus_all_238U_lowE_highT_134Te = y_doublegate_3n_6plus_all_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_6plus_all fit range %d - %.d" % (x_doublegate_3n_6plus_all_238U_lowE_highT_134Te[0],x_doublegate_3n_6plus_all_238U_lowE_highT_134Te[-1]))



################   238U lowE_highT bg -  134Te 3n_6plus  #################

#doublegate_3n_6plus_bg true
hist_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_6plus_bg_134Te')
x_bins = hist_doublegate_3n_6plus_bg_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = 2*x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te_long = x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te
y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te_long = y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te = y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_6plus_bg fit range %d - %.d" % (x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[0],x_doublegate_3n_6plus_bg_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_ridge -  134Te 3n_6plus  #################

#doublegate_3n_6plus_bg_ridge true
hist_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_6plus_bg_ridge_134Te')
x_bins = hist_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = 2*x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te_long = x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te
y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te_long = y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te = y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_6plus_bg_ridge fit range %d - %.d" % (x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[0],x_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te[-1]))


################   238U lowE_highT bg_random -  134Te 3n_6plus  #################

#doublegate_3n_6plus_bg_random true
hist_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = file_238U_lowE_highT.Get('time_isomer_doublegate_3n_6plus_bg_random_134Te')
x_bins = hist_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te.GetNbinsX()

x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[i] = hist_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = 2*x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te_long = x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te
y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te_long = y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te

x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]
y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te = y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te 3n_6plus_bg_random fit range %d - %.d" % (x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[0],x_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te[-1]))



################   238U lowE_lowT -  134Te   #################

#Define lowE_lowT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_238U_lowE_lowT_134Te = file_238U_lowE_lowT.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_238U_lowE_lowT_134Te.GetNbinsX()

x_doublegate_238U_lowE_lowT_134Te = np.zeros(x_bins)
y_doublegate_238U_lowE_lowT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_lowE_lowT_134Te[i] = hist_doublegate_238U_lowE_lowT_134Te.GetBinCenter(i+1)
    y_doublegate_238U_lowE_lowT_134Te[i] = hist_doublegate_238U_lowE_lowT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_lowE_lowT_134Te = 2*x_doublegate_238U_lowE_lowT_134Te

x_doublegate_238U_lowE_lowT_134Te_long = x_doublegate_238U_lowE_lowT_134Te
y_doublegate_238U_lowE_lowT_134Te_long = y_doublegate_238U_lowE_lowT_134Te

x_doublegate_238U_lowE_lowT_134Te = x_doublegate_238U_lowE_lowT_134Te[bin_lower:bin_upper]
y_doublegate_238U_lowE_lowT_134Te = y_doublegate_238U_lowE_lowT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_238U_lowE_lowT_134Te[0],x_doublegate_238U_lowE_lowT_134Te[-1]))




################   238U highE -  134Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_238U_highE_highT_134Te = file_238U_highE_highT.Get('time_isomer_doublegate_134Te')
x_bins = hist_doublegate_238U_highE_highT_134Te.GetNbinsX()

x_doublegate_238U_highE_highT_134Te = np.zeros(x_bins)
y_doublegate_238U_highE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_highE_highT_134Te[i] = hist_doublegate_238U_highE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_238U_highE_highT_134Te[i] = hist_doublegate_238U_highE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_highE_highT_134Te = 2*x_doublegate_238U_highE_highT_134Te

x_doublegate_238U_highE_highT_134Te_long = x_doublegate_238U_highE_highT_134Te
y_doublegate_238U_highE_highT_134Te_long = y_doublegate_238U_highE_highT_134Te

x_doublegate_238U_highE_highT_134Te = x_doublegate_238U_highE_highT_134Te[bin_lower:bin_upper]
y_doublegate_238U_highE_highT_134Te = y_doublegate_238U_highE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_238U_highE_highT_134Te[0],x_doublegate_238U_highE_highT_134Te[-1]))


################   238U highE all -  134Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_all_238U_highE_highT_134Te = file_238U_highE_highT.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_238U_highE_highT_134Te.GetNbinsX()

x_doublegate_all_238U_highE_highT_134Te = np.zeros(x_bins)
y_doublegate_all_238U_highE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_238U_highE_highT_134Te[i] = hist_doublegate_all_238U_highE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_all_238U_highE_highT_134Te[i] = hist_doublegate_all_238U_highE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_238U_highE_highT_134Te = 2*x_doublegate_all_238U_highE_highT_134Te

x_doublegate_all_238U_highE_highT_134Te_long = x_doublegate_all_238U_highE_highT_134Te
y_doublegate_all_238U_highE_highT_134Te_long = y_doublegate_all_238U_highE_highT_134Te

x_doublegate_all_238U_highE_highT_134Te = x_doublegate_all_238U_highE_highT_134Te[bin_lower:bin_upper]
y_doublegate_all_238U_highE_highT_134Te = y_doublegate_all_238U_highE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_all_238U_highE_highT_134Te[0],x_doublegate_all_238U_highE_highT_134Te[-1]))


################   238U highE_highT bg -  134Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_bg_238U_highE_highT_134Te = file_238U_highE_highT.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_238U_highE_highT_134Te.GetNbinsX()

x_doublegate_bg_238U_highE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_238U_highE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_238U_highE_highT_134Te[i] = hist_doublegate_bg_238U_highE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_238U_highE_highT_134Te[i] = hist_doublegate_bg_238U_highE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_238U_highE_highT_134Te = 2*x_doublegate_bg_238U_highE_highT_134Te

x_doublegate_bg_238U_highE_highT_134Te_long = x_doublegate_bg_238U_highE_highT_134Te
y_doublegate_bg_238U_highE_highT_134Te_long = y_doublegate_bg_238U_highE_highT_134Te

x_doublegate_bg_238U_highE_highT_134Te = x_doublegate_bg_238U_highE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_238U_highE_highT_134Te = y_doublegate_bg_238U_highE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_238U_highE_highT_134Te[0],x_doublegate_bg_238U_highE_highT_134Te[-1]))


################   238U highE_highT bg_ridge -  134Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_bg_ridge_238U_highE_highT_134Te = file_238U_highE_highT.Get('time_isomer_doublegate_bg_ridge_134Te')
x_bins = hist_doublegate_bg_ridge_238U_highE_highT_134Te.GetNbinsX()

x_doublegate_bg_ridge_238U_highE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_ridge_238U_highE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_ridge_238U_highE_highT_134Te[i] = hist_doublegate_bg_ridge_238U_highE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_ridge_238U_highE_highT_134Te[i] = hist_doublegate_bg_ridge_238U_highE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_ridge_238U_highE_highT_134Te = 2*x_doublegate_bg_ridge_238U_highE_highT_134Te

x_doublegate_bg_ridge_238U_highE_highT_134Te_long = x_doublegate_bg_ridge_238U_highE_highT_134Te
y_doublegate_bg_ridge_238U_highE_highT_134Te_long = y_doublegate_bg_ridge_238U_highE_highT_134Te

x_doublegate_bg_ridge_238U_highE_highT_134Te = x_doublegate_bg_ridge_238U_highE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_ridge_238U_highE_highT_134Te = y_doublegate_bg_ridge_238U_highE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_ridge_238U_highE_highT_134Te[0],x_doublegate_bg_ridge_238U_highE_highT_134Te[-1]))


################   238U highE_highT bg_random -  134Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 640

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_bg_random_238U_highE_highT_134Te = file_238U_highE_highT.Get('time_isomer_doublegate_bg_random_134Te')
x_bins = hist_doublegate_bg_random_238U_highE_highT_134Te.GetNbinsX()

x_doublegate_bg_random_238U_highE_highT_134Te = np.zeros(x_bins)
y_doublegate_bg_random_238U_highE_highT_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_random_238U_highE_highT_134Te[i] = hist_doublegate_bg_random_238U_highE_highT_134Te.GetBinCenter(i+1)
    y_doublegate_bg_random_238U_highE_highT_134Te[i] = hist_doublegate_bg_random_238U_highE_highT_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_random_238U_highE_highT_134Te = 2*x_doublegate_bg_random_238U_highE_highT_134Te

x_doublegate_bg_random_238U_highE_highT_134Te_long = x_doublegate_bg_random_238U_highE_highT_134Te
y_doublegate_bg_random_238U_highE_highT_134Te_long = y_doublegate_bg_random_238U_highE_highT_134Te

x_doublegate_bg_random_238U_highE_highT_134Te = x_doublegate_bg_random_238U_highE_highT_134Te[bin_lower:bin_upper]
y_doublegate_bg_random_238U_highE_highT_134Te = y_doublegate_bg_random_238U_highE_highT_134Te[bin_lower:bin_upper]

#print("* 238U - 134Te fit range %d - %.d" % (x_doublegate_bg_random_238U_highE_highT_134Te[0],x_doublegate_bg_random_238U_highE_highT_134Te[-1]))



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


################   232Th all -  134Te   #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_all true
hist_doublegate_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_all_134Te')
x_bins = hist_doublegate_all_232Th_134Te.GetNbinsX()

x_doublegate_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_all_232Th_134Te[i] = hist_doublegate_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_all_232Th_134Te[i] = hist_doublegate_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_all_232Th_134Te = 2*x_doublegate_all_232Th_134Te

x_doublegate_all_232Th_134Te_long = x_doublegate_all_232Th_134Te
y_doublegate_all_232Th_134Te_long = y_doublegate_all_232Th_134Te

x_doublegate_all_232Th_134Te = x_doublegate_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_all_232Th_134Te = y_doublegate_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_all_232Th_134Te[0],x_doublegate_all_232Th_134Te[-1]))



################   232Th bg -  134Te   #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_bg true
hist_doublegate_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_bg_134Te')
x_bins = hist_doublegate_bg_232Th_134Te.GetNbinsX()

x_doublegate_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_232Th_134Te[i] = hist_doublegate_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_bg_232Th_134Te[i] = hist_doublegate_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_232Th_134Te = 2*x_doublegate_bg_232Th_134Te

x_doublegate_bg_232Th_134Te_long = x_doublegate_bg_232Th_134Te
y_doublegate_bg_232Th_134Te_long = y_doublegate_bg_232Th_134Te

x_doublegate_bg_232Th_134Te = x_doublegate_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_bg_232Th_134Te = y_doublegate_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_bg_232Th_134Te[0],x_doublegate_bg_232Th_134Te[-1]))


################   232Th bg_ridge -  134Te   #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_bg_ridge true
hist_doublegate_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_bg_ridge_134Te')
x_bins = hist_doublegate_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_ridge_232Th_134Te[i] = hist_doublegate_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_bg_ridge_232Th_134Te[i] = hist_doublegate_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_ridge_232Th_134Te = 2*x_doublegate_bg_ridge_232Th_134Te

x_doublegate_bg_ridge_232Th_134Te_long = x_doublegate_bg_ridge_232Th_134Te
y_doublegate_bg_ridge_232Th_134Te_long = y_doublegate_bg_ridge_232Th_134Te

x_doublegate_bg_ridge_232Th_134Te = x_doublegate_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_bg_ridge_232Th_134Te = y_doublegate_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_bg_ridge_232Th_134Te[0],x_doublegate_bg_ridge_232Th_134Te[-1]))


################   232Th bg_random -  134Te   #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_bg_random true
hist_doublegate_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_bg_random_134Te')
x_bins = hist_doublegate_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_bg_random_232Th_134Te[i] = hist_doublegate_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_bg_random_232Th_134Te[i] = hist_doublegate_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_bg_random_232Th_134Te = 2*x_doublegate_bg_random_232Th_134Te

x_doublegate_bg_random_232Th_134Te_long = x_doublegate_bg_random_232Th_134Te
y_doublegate_bg_random_232Th_134Te_long = y_doublegate_bg_random_232Th_134Te

x_doublegate_bg_random_232Th_134Te = x_doublegate_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_bg_random_232Th_134Te = y_doublegate_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_bg_random_232Th_134Te[0],x_doublegate_bg_random_232Th_134Te[-1]))



################   232Th -  134Te 1n  #################

#doublegate_1n true
hist_doublegate_1n_232Th_134Te = file_232Th.Get('time_isomer_doublegate_1n_134Te')
x_bins = hist_doublegate_1n_232Th_134Te.GetNbinsX()

x_doublegate_1n_232Th_134Te = np.zeros(x_bins)
y_doublegate_1n_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_232Th_134Te[i] = hist_doublegate_1n_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_1n_232Th_134Te[i] = hist_doublegate_1n_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_232Th_134Te = 2*x_doublegate_1n_232Th_134Te

x_doublegate_1n_232Th_134Te_long = x_doublegate_1n_232Th_134Te
y_doublegate_1n_232Th_134Te_long = y_doublegate_1n_232Th_134Te

x_doublegate_1n_232Th_134Te = x_doublegate_1n_232Th_134Te[bin_lower:bin_upper]
y_doublegate_1n_232Th_134Te = y_doublegate_1n_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 1n fit range %d - %.d" % (x_doublegate_1n_232Th_134Te[0],x_doublegate_1n_232Th_134Te[-1]))


################   232Th all -  134Te 1n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_1n_all true
hist_doublegate_1n_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_1n_all_134Te')
x_bins = hist_doublegate_1n_all_232Th_134Te.GetNbinsX()

x_doublegate_1n_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_1n_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_all_232Th_134Te[i] = hist_doublegate_1n_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_1n_all_232Th_134Te[i] = hist_doublegate_1n_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_all_232Th_134Te = 2*x_doublegate_1n_all_232Th_134Te

x_doublegate_1n_all_232Th_134Te_long = x_doublegate_1n_all_232Th_134Te
y_doublegate_1n_all_232Th_134Te_long = y_doublegate_1n_all_232Th_134Te

x_doublegate_1n_all_232Th_134Te = x_doublegate_1n_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_1n_all_232Th_134Te = y_doublegate_1n_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_1n_all_232Th_134Te[0],x_doublegate_1n_all_232Th_134Te[-1]))



################   232Th bg -  134Te 1n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_1n_bg true
hist_doublegate_1n_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_1n_bg_134Te')
x_bins = hist_doublegate_1n_bg_232Th_134Te.GetNbinsX()

x_doublegate_1n_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_232Th_134Te[i] = hist_doublegate_1n_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_232Th_134Te[i] = hist_doublegate_1n_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_232Th_134Te = 2*x_doublegate_1n_bg_232Th_134Te

x_doublegate_1n_bg_232Th_134Te_long = x_doublegate_1n_bg_232Th_134Te
y_doublegate_1n_bg_232Th_134Te_long = y_doublegate_1n_bg_232Th_134Te

x_doublegate_1n_bg_232Th_134Te = x_doublegate_1n_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_232Th_134Te = y_doublegate_1n_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_1n_bg_232Th_134Te[0],x_doublegate_1n_bg_232Th_134Te[-1]))

################   232Th bg_ridge -  134Te 1n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_1n_bg_ridge true
hist_doublegate_1n_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_1n_bg_ridge_134Te')
x_bins = hist_doublegate_1n_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_1n_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_ridge_232Th_134Te[i] = hist_doublegate_1n_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_ridge_232Th_134Te[i] = hist_doublegate_1n_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_ridge_232Th_134Te = 2*x_doublegate_1n_bg_ridge_232Th_134Te

x_doublegate_1n_bg_ridge_232Th_134Te_long = x_doublegate_1n_bg_ridge_232Th_134Te
y_doublegate_1n_bg_ridge_232Th_134Te_long = y_doublegate_1n_bg_ridge_232Th_134Te

x_doublegate_1n_bg_ridge_232Th_134Te = x_doublegate_1n_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_ridge_232Th_134Te = y_doublegate_1n_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_1n_bg_ridge_232Th_134Te[0],x_doublegate_1n_bg_ridge_232Th_134Te[-1]))


################   232Th bg_random -  134Te 1n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_1n_bg_random true
hist_doublegate_1n_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_1n_bg_random_134Te')
x_bins = hist_doublegate_1n_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_1n_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_1n_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_1n_bg_random_232Th_134Te[i] = hist_doublegate_1n_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_1n_bg_random_232Th_134Te[i] = hist_doublegate_1n_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_1n_bg_random_232Th_134Te = 2*x_doublegate_1n_bg_random_232Th_134Te

x_doublegate_1n_bg_random_232Th_134Te_long = x_doublegate_1n_bg_random_232Th_134Te
y_doublegate_1n_bg_random_232Th_134Te_long = y_doublegate_1n_bg_random_232Th_134Te

x_doublegate_1n_bg_random_232Th_134Te = x_doublegate_1n_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_1n_bg_random_232Th_134Te = y_doublegate_1n_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_1n_bg_random_232Th_134Te[0],x_doublegate_1n_bg_random_232Th_134Te[-1]))



################   232Th -  134Te 3n  #################

#doublegate_3n true
hist_doublegate_3n_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_134Te')
x_bins = hist_doublegate_3n_232Th_134Te.GetNbinsX()

x_doublegate_3n_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_232Th_134Te[i] = hist_doublegate_3n_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_232Th_134Te[i] = hist_doublegate_3n_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_232Th_134Te = 2*x_doublegate_3n_232Th_134Te

x_doublegate_3n_232Th_134Te_long = x_doublegate_3n_232Th_134Te
y_doublegate_3n_232Th_134Te_long = y_doublegate_3n_232Th_134Te

x_doublegate_3n_232Th_134Te = x_doublegate_3n_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_232Th_134Te = y_doublegate_3n_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n fit range %d - %.d" % (x_doublegate_3n_232Th_134Te[0],x_doublegate_3n_232Th_134Te[-1]))


################   232Th all -  134Te 3n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_3n_all true
hist_doublegate_3n_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_all_134Te')
x_bins = hist_doublegate_3n_all_232Th_134Te.GetNbinsX()

x_doublegate_3n_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_all_232Th_134Te[i] = hist_doublegate_3n_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_all_232Th_134Te[i] = hist_doublegate_3n_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_all_232Th_134Te = 2*x_doublegate_3n_all_232Th_134Te

x_doublegate_3n_all_232Th_134Te_long = x_doublegate_3n_all_232Th_134Te
y_doublegate_3n_all_232Th_134Te_long = y_doublegate_3n_all_232Th_134Te

x_doublegate_3n_all_232Th_134Te = x_doublegate_3n_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_all_232Th_134Te = y_doublegate_3n_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_3n_all_232Th_134Te[0],x_doublegate_3n_all_232Th_134Te[-1]))



################   232Th bg -  134Te 3n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_3n_bg true
hist_doublegate_3n_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_bg_134Te')
x_bins = hist_doublegate_3n_bg_232Th_134Te.GetNbinsX()

x_doublegate_3n_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_232Th_134Te[i] = hist_doublegate_3n_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_232Th_134Te[i] = hist_doublegate_3n_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_232Th_134Te = 2*x_doublegate_3n_bg_232Th_134Te

x_doublegate_3n_bg_232Th_134Te_long = x_doublegate_3n_bg_232Th_134Te
y_doublegate_3n_bg_232Th_134Te_long = y_doublegate_3n_bg_232Th_134Te

x_doublegate_3n_bg_232Th_134Te = x_doublegate_3n_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_232Th_134Te = y_doublegate_3n_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_3n_bg_232Th_134Te[0],x_doublegate_3n_bg_232Th_134Te[-1]))


################   232Th bg_ridge -  134Te 3n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_3n_bg_ridge true
hist_doublegate_3n_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_bg_ridge_134Te')
x_bins = hist_doublegate_3n_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_3n_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_ridge_232Th_134Te = 2*x_doublegate_3n_bg_ridge_232Th_134Te

x_doublegate_3n_bg_ridge_232Th_134Te_long = x_doublegate_3n_bg_ridge_232Th_134Te
y_doublegate_3n_bg_ridge_232Th_134Te_long = y_doublegate_3n_bg_ridge_232Th_134Te

x_doublegate_3n_bg_ridge_232Th_134Te = x_doublegate_3n_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_ridge_232Th_134Te = y_doublegate_3n_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_3n_bg_ridge_232Th_134Te[0],x_doublegate_3n_bg_ridge_232Th_134Te[-1]))


################   232Th bg_random -  134Te 3n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_3n_bg_random true
hist_doublegate_3n_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_bg_random_134Te')
x_bins = hist_doublegate_3n_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_3n_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_bg_random_232Th_134Te[i] = hist_doublegate_3n_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_bg_random_232Th_134Te[i] = hist_doublegate_3n_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_bg_random_232Th_134Te = 2*x_doublegate_3n_bg_random_232Th_134Te

x_doublegate_3n_bg_random_232Th_134Te_long = x_doublegate_3n_bg_random_232Th_134Te
y_doublegate_3n_bg_random_232Th_134Te_long = y_doublegate_3n_bg_random_232Th_134Te

x_doublegate_3n_bg_random_232Th_134Te = x_doublegate_3n_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_bg_random_232Th_134Te = y_doublegate_3n_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_3n_bg_random_232Th_134Te[0],x_doublegate_3n_bg_random_232Th_134Te[-1]))



################   232Th -  134Te 4n  #################

#doublegate_4n true
hist_doublegate_4n_232Th_134Te = file_232Th.Get('time_isomer_doublegate_4n_134Te')
x_bins = hist_doublegate_4n_232Th_134Te.GetNbinsX()

x_doublegate_4n_232Th_134Te = np.zeros(x_bins)
y_doublegate_4n_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_4n_232Th_134Te[i] = hist_doublegate_4n_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_4n_232Th_134Te[i] = hist_doublegate_4n_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_4n_232Th_134Te = 2*x_doublegate_4n_232Th_134Te

x_doublegate_4n_232Th_134Te_long = x_doublegate_4n_232Th_134Te
y_doublegate_4n_232Th_134Te_long = y_doublegate_4n_232Th_134Te

x_doublegate_4n_232Th_134Te = x_doublegate_4n_232Th_134Te[bin_lower:bin_upper]
y_doublegate_4n_232Th_134Te = y_doublegate_4n_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 4n fit range %d - %.d" % (x_doublegate_4n_232Th_134Te[0],x_doublegate_4n_232Th_134Te[-1]))


################   232Th all -  134Te 4n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_4n_all true
hist_doublegate_4n_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_4n_all_134Te')
x_bins = hist_doublegate_4n_all_232Th_134Te.GetNbinsX()

x_doublegate_4n_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_4n_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_4n_all_232Th_134Te[i] = hist_doublegate_4n_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_4n_all_232Th_134Te[i] = hist_doublegate_4n_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_4n_all_232Th_134Te = 2*x_doublegate_4n_all_232Th_134Te

x_doublegate_4n_all_232Th_134Te_long = x_doublegate_4n_all_232Th_134Te
y_doublegate_4n_all_232Th_134Te_long = y_doublegate_4n_all_232Th_134Te

x_doublegate_4n_all_232Th_134Te = x_doublegate_4n_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_4n_all_232Th_134Te = y_doublegate_4n_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_4n_all_232Th_134Te[0],x_doublegate_4n_all_232Th_134Te[-1]))



################   232Th bg -  134Te 4n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_4n_bg true
hist_doublegate_4n_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_4n_bg_134Te')
x_bins = hist_doublegate_4n_bg_232Th_134Te.GetNbinsX()

x_doublegate_4n_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_4n_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_4n_bg_232Th_134Te[i] = hist_doublegate_4n_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_4n_bg_232Th_134Te[i] = hist_doublegate_4n_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_4n_bg_232Th_134Te = 2*x_doublegate_4n_bg_232Th_134Te

x_doublegate_4n_bg_232Th_134Te_long = x_doublegate_4n_bg_232Th_134Te
y_doublegate_4n_bg_232Th_134Te_long = y_doublegate_4n_bg_232Th_134Te

x_doublegate_4n_bg_232Th_134Te = x_doublegate_4n_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_4n_bg_232Th_134Te = y_doublegate_4n_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_4n_bg_232Th_134Te[0],x_doublegate_4n_bg_232Th_134Te[-1]))


################   232Th bg_ridge -  134Te 4n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_4n_bg_ridge true
hist_doublegate_4n_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_4n_bg_ridge_134Te')
x_bins = hist_doublegate_4n_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_4n_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_4n_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_4n_bg_ridge_232Th_134Te[i] = hist_doublegate_4n_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_4n_bg_ridge_232Th_134Te[i] = hist_doublegate_4n_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_4n_bg_ridge_232Th_134Te = 2*x_doublegate_4n_bg_ridge_232Th_134Te

x_doublegate_4n_bg_ridge_232Th_134Te_long = x_doublegate_4n_bg_ridge_232Th_134Te
y_doublegate_4n_bg_ridge_232Th_134Te_long = y_doublegate_4n_bg_ridge_232Th_134Te

x_doublegate_4n_bg_ridge_232Th_134Te = x_doublegate_4n_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_4n_bg_ridge_232Th_134Te = y_doublegate_4n_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_4n_bg_ridge_232Th_134Te[0],x_doublegate_4n_bg_ridge_232Th_134Te[-1]))

################   232Th bg_random -  134Te 4n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_4n_bg_random true
hist_doublegate_4n_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_4n_bg_random_134Te')
x_bins = hist_doublegate_4n_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_4n_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_4n_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_4n_bg_random_232Th_134Te[i] = hist_doublegate_4n_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_4n_bg_random_232Th_134Te[i] = hist_doublegate_4n_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_4n_bg_random_232Th_134Te = 2*x_doublegate_4n_bg_random_232Th_134Te

x_doublegate_4n_bg_random_232Th_134Te_long = x_doublegate_4n_bg_random_232Th_134Te
y_doublegate_4n_bg_random_232Th_134Te_long = y_doublegate_4n_bg_random_232Th_134Te

x_doublegate_4n_bg_random_232Th_134Te = x_doublegate_4n_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_4n_bg_random_232Th_134Te = y_doublegate_4n_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_4n_bg_random_232Th_134Te[0],x_doublegate_4n_bg_random_232Th_134Te[-1]))



################   232Th -  134Te 5n  #################

#doublegate_5n true
hist_doublegate_5n_232Th_134Te = file_232Th.Get('time_isomer_doublegate_5n_134Te')
x_bins = hist_doublegate_5n_232Th_134Te.GetNbinsX()

x_doublegate_5n_232Th_134Te = np.zeros(x_bins)
y_doublegate_5n_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_232Th_134Te[i] = hist_doublegate_5n_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_5n_232Th_134Te[i] = hist_doublegate_5n_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_232Th_134Te = 2*x_doublegate_5n_232Th_134Te

x_doublegate_5n_232Th_134Te_long = x_doublegate_5n_232Th_134Te
y_doublegate_5n_232Th_134Te_long = y_doublegate_5n_232Th_134Te

x_doublegate_5n_232Th_134Te = x_doublegate_5n_232Th_134Te[bin_lower:bin_upper]
y_doublegate_5n_232Th_134Te = y_doublegate_5n_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 5n fit range %d - %.d" % (x_doublegate_5n_232Th_134Te[0],x_doublegate_5n_232Th_134Te[-1]))


################   232Th all -  134Te 5n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_5n_all true
hist_doublegate_5n_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_5n_all_134Te')
x_bins = hist_doublegate_5n_all_232Th_134Te.GetNbinsX()

x_doublegate_5n_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_5n_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_all_232Th_134Te[i] = hist_doublegate_5n_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_5n_all_232Th_134Te[i] = hist_doublegate_5n_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_all_232Th_134Te = 2*x_doublegate_5n_all_232Th_134Te

x_doublegate_5n_all_232Th_134Te_long = x_doublegate_5n_all_232Th_134Te
y_doublegate_5n_all_232Th_134Te_long = y_doublegate_5n_all_232Th_134Te

x_doublegate_5n_all_232Th_134Te = x_doublegate_5n_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_5n_all_232Th_134Te = y_doublegate_5n_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_5n_all_232Th_134Te[0],x_doublegate_5n_all_232Th_134Te[-1]))



################   232Th bg -  134Te 5n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_5n_bg true
hist_doublegate_5n_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_5n_bg_134Te')
x_bins = hist_doublegate_5n_bg_232Th_134Te.GetNbinsX()

x_doublegate_5n_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_5n_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_bg_232Th_134Te[i] = hist_doublegate_5n_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_5n_bg_232Th_134Te[i] = hist_doublegate_5n_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_bg_232Th_134Te = 2*x_doublegate_5n_bg_232Th_134Te

x_doublegate_5n_bg_232Th_134Te_long = x_doublegate_5n_bg_232Th_134Te
y_doublegate_5n_bg_232Th_134Te_long = y_doublegate_5n_bg_232Th_134Te

x_doublegate_5n_bg_232Th_134Te = x_doublegate_5n_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_5n_bg_232Th_134Te = y_doublegate_5n_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_5n_bg_232Th_134Te[0],x_doublegate_5n_bg_232Th_134Te[-1]))


################   232Th bg_ridge -  134Te 5n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_5n_bg_ridge true
hist_doublegate_5n_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_5n_bg_ridge_134Te')
x_bins = hist_doublegate_5n_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_5n_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_5n_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_bg_ridge_232Th_134Te[i] = hist_doublegate_5n_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_5n_bg_ridge_232Th_134Te[i] = hist_doublegate_5n_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_bg_ridge_232Th_134Te = 2*x_doublegate_5n_bg_ridge_232Th_134Te

x_doublegate_5n_bg_ridge_232Th_134Te_long = x_doublegate_5n_bg_ridge_232Th_134Te
y_doublegate_5n_bg_ridge_232Th_134Te_long = y_doublegate_5n_bg_ridge_232Th_134Te

x_doublegate_5n_bg_ridge_232Th_134Te = x_doublegate_5n_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_5n_bg_ridge_232Th_134Te = y_doublegate_5n_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_5n_bg_ridge_232Th_134Te[0],x_doublegate_5n_bg_ridge_232Th_134Te[-1]))

################   232Th bg_random -  134Te 5n  #################

bin_lower = x_lower//2
bin_upper = x_upper//2

#doublegate_5n_bg_random true
hist_doublegate_5n_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_5n_bg_random_134Te')
x_bins = hist_doublegate_5n_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_5n_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_5n_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_5n_bg_random_232Th_134Te[i] = hist_doublegate_5n_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_5n_bg_random_232Th_134Te[i] = hist_doublegate_5n_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_5n_bg_random_232Th_134Te = 2*x_doublegate_5n_bg_random_232Th_134Te

x_doublegate_5n_bg_random_232Th_134Te_long = x_doublegate_5n_bg_random_232Th_134Te
y_doublegate_5n_bg_random_232Th_134Te_long = y_doublegate_5n_bg_random_232Th_134Te

x_doublegate_5n_bg_random_232Th_134Te = x_doublegate_5n_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_5n_bg_random_232Th_134Te = y_doublegate_5n_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te fit range %d - %.d" % (x_doublegate_5n_bg_random_232Th_134Te[0],x_doublegate_5n_bg_random_232Th_134Te[-1]))


################   232Th -  134Te 3n_2plus  #################

#doublegate_3n_2plus true
hist_doublegate_3n_2plus_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_2plus_134Te')
x_bins = hist_doublegate_3n_2plus_232Th_134Te.GetNbinsX()

x_doublegate_3n_2plus_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_232Th_134Te[i] = hist_doublegate_3n_2plus_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_232Th_134Te[i] = hist_doublegate_3n_2plus_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_232Th_134Te = 2*x_doublegate_3n_2plus_232Th_134Te

x_doublegate_3n_2plus_232Th_134Te_long = x_doublegate_3n_2plus_232Th_134Te
y_doublegate_3n_2plus_232Th_134Te_long = y_doublegate_3n_2plus_232Th_134Te

x_doublegate_3n_2plus_232Th_134Te = x_doublegate_3n_2plus_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_232Th_134Te = y_doublegate_3n_2plus_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_2plus fit range %d - %.d" % (x_doublegate_3n_2plus_232Th_134Te[0],x_doublegate_3n_2plus_232Th_134Te[-1]))


################   232Th all -  134Te 3n_2plus  #################

#doublegate_3n_2plus_all true
hist_doublegate_3n_2plus_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_2plus_all_134Te')
x_bins = hist_doublegate_3n_2plus_all_232Th_134Te.GetNbinsX()

x_doublegate_3n_2plus_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_all_232Th_134Te[i] = hist_doublegate_3n_2plus_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_all_232Th_134Te[i] = hist_doublegate_3n_2plus_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_all_232Th_134Te = 2*x_doublegate_3n_2plus_all_232Th_134Te

x_doublegate_3n_2plus_all_232Th_134Te_long = x_doublegate_3n_2plus_all_232Th_134Te
y_doublegate_3n_2plus_all_232Th_134Te_long = y_doublegate_3n_2plus_all_232Th_134Te

x_doublegate_3n_2plus_all_232Th_134Te = x_doublegate_3n_2plus_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_all_232Th_134Te = y_doublegate_3n_2plus_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_2plus_all fit range %d - %.d" % (x_doublegate_3n_2plus_all_232Th_134Te[0],x_doublegate_3n_2plus_all_232Th_134Te[-1]))



################   232Th bg -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg true
hist_doublegate_3n_2plus_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_2plus_bg_134Te')
x_bins = hist_doublegate_3n_2plus_bg_232Th_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_232Th_134Te = 2*x_doublegate_3n_2plus_bg_232Th_134Te

x_doublegate_3n_2plus_bg_232Th_134Te_long = x_doublegate_3n_2plus_bg_232Th_134Te
y_doublegate_3n_2plus_bg_232Th_134Te_long = y_doublegate_3n_2plus_bg_232Th_134Te

x_doublegate_3n_2plus_bg_232Th_134Te = x_doublegate_3n_2plus_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_232Th_134Te = y_doublegate_3n_2plus_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_2plus_bg fit range %d - %.d" % (x_doublegate_3n_2plus_bg_232Th_134Te[0],x_doublegate_3n_2plus_bg_232Th_134Te[-1]))


################   232Th bg_ridge -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg_ridge true
hist_doublegate_3n_2plus_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_2plus_bg_ridge_134Te')
x_bins = hist_doublegate_3n_2plus_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_ridge_232Th_134Te = 2*x_doublegate_3n_2plus_bg_ridge_232Th_134Te

x_doublegate_3n_2plus_bg_ridge_232Th_134Te_long = x_doublegate_3n_2plus_bg_ridge_232Th_134Te
y_doublegate_3n_2plus_bg_ridge_232Th_134Te_long = y_doublegate_3n_2plus_bg_ridge_232Th_134Te

x_doublegate_3n_2plus_bg_ridge_232Th_134Te = x_doublegate_3n_2plus_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_ridge_232Th_134Te = y_doublegate_3n_2plus_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_2plus_bg_ridge fit range %d - %.d" % (x_doublegate_3n_2plus_bg_ridge_232Th_134Te[0],x_doublegate_3n_2plus_bg_ridge_232Th_134Te[-1]))


################   232Th bg_random -  134Te 3n_2plus  #################

#doublegate_3n_2plus_bg_random true
hist_doublegate_3n_2plus_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_2plus_bg_random_134Te')
x_bins = hist_doublegate_3n_2plus_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_3n_2plus_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_2plus_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_2plus_bg_random_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_2plus_bg_random_232Th_134Te[i] = hist_doublegate_3n_2plus_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_2plus_bg_random_232Th_134Te = 2*x_doublegate_3n_2plus_bg_random_232Th_134Te

x_doublegate_3n_2plus_bg_random_232Th_134Te_long = x_doublegate_3n_2plus_bg_random_232Th_134Te
y_doublegate_3n_2plus_bg_random_232Th_134Te_long = y_doublegate_3n_2plus_bg_random_232Th_134Te

x_doublegate_3n_2plus_bg_random_232Th_134Te = x_doublegate_3n_2plus_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_2plus_bg_random_232Th_134Te = y_doublegate_3n_2plus_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_2plus_bg_random fit range %d - %.d" % (x_doublegate_3n_2plus_bg_random_232Th_134Te[0],x_doublegate_3n_2plus_bg_random_232Th_134Te[-1]))



################   232Th -  134Te 3n_4plus  #################

#doublegate_3n_4plus true
hist_doublegate_3n_4plus_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_4plus_134Te')
x_bins = hist_doublegate_3n_4plus_232Th_134Te.GetNbinsX()

x_doublegate_3n_4plus_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_232Th_134Te[i] = hist_doublegate_3n_4plus_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_232Th_134Te[i] = hist_doublegate_3n_4plus_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_232Th_134Te = 2*x_doublegate_3n_4plus_232Th_134Te

x_doublegate_3n_4plus_232Th_134Te_long = x_doublegate_3n_4plus_232Th_134Te
y_doublegate_3n_4plus_232Th_134Te_long = y_doublegate_3n_4plus_232Th_134Te

x_doublegate_3n_4plus_232Th_134Te = x_doublegate_3n_4plus_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_232Th_134Te = y_doublegate_3n_4plus_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_4plus fit range %d - %.d" % (x_doublegate_3n_4plus_232Th_134Te[0],x_doublegate_3n_4plus_232Th_134Te[-1]))


################   232Th lowE_highT all -  134Te 3n_4plus  #################

#doublegate_3n_4plus_all true
hist_doublegate_3n_4plus_all_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_4plus_all_134Te')
x_bins = hist_doublegate_3n_4plus_all_232Th_134Te.GetNbinsX()

x_doublegate_3n_4plus_all_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_all_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_all_232Th_134Te[i] = hist_doublegate_3n_4plus_all_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_all_232Th_134Te[i] = hist_doublegate_3n_4plus_all_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_all_232Th_134Te = 2*x_doublegate_3n_4plus_all_232Th_134Te

x_doublegate_3n_4plus_all_232Th_134Te_long = x_doublegate_3n_4plus_all_232Th_134Te
y_doublegate_3n_4plus_all_232Th_134Te_long = y_doublegate_3n_4plus_all_232Th_134Te

x_doublegate_3n_4plus_all_232Th_134Te = x_doublegate_3n_4plus_all_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_all_232Th_134Te = y_doublegate_3n_4plus_all_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_4plus_all fit range %d - %.d" % (x_doublegate_3n_4plus_all_232Th_134Te[0],x_doublegate_3n_4plus_all_232Th_134Te[-1]))



################   232Th lowE_highT bg -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg true
hist_doublegate_3n_4plus_bg_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_4plus_bg_134Te')
x_bins = hist_doublegate_3n_4plus_bg_232Th_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_232Th_134Te = 2*x_doublegate_3n_4plus_bg_232Th_134Te

x_doublegate_3n_4plus_bg_232Th_134Te_long = x_doublegate_3n_4plus_bg_232Th_134Te
y_doublegate_3n_4plus_bg_232Th_134Te_long = y_doublegate_3n_4plus_bg_232Th_134Te

x_doublegate_3n_4plus_bg_232Th_134Te = x_doublegate_3n_4plus_bg_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_232Th_134Te = y_doublegate_3n_4plus_bg_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_4plus_bg fit range %d - %.d" % (x_doublegate_3n_4plus_bg_232Th_134Te[0],x_doublegate_3n_4plus_bg_232Th_134Te[-1]))


################   232Th lowE_highT bg_ridge -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg_ridge true
hist_doublegate_3n_4plus_bg_ridge_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_4plus_bg_ridge_134Te')
x_bins = hist_doublegate_3n_4plus_bg_ridge_232Th_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_ridge_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_ridge_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_ridge_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_ridge_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_ridge_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_ridge_232Th_134Te = 2*x_doublegate_3n_4plus_bg_ridge_232Th_134Te

x_doublegate_3n_4plus_bg_ridge_232Th_134Te_long = x_doublegate_3n_4plus_bg_ridge_232Th_134Te
y_doublegate_3n_4plus_bg_ridge_232Th_134Te_long = y_doublegate_3n_4plus_bg_ridge_232Th_134Te

x_doublegate_3n_4plus_bg_ridge_232Th_134Te = x_doublegate_3n_4plus_bg_ridge_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_ridge_232Th_134Te = y_doublegate_3n_4plus_bg_ridge_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_4plus_bg_ridge fit range %d - %.d" % (x_doublegate_3n_4plus_bg_ridge_232Th_134Te[0],x_doublegate_3n_4plus_bg_ridge_232Th_134Te[-1]))


################   232Th lowE_highT bg_random -  134Te 3n_4plus  #################

#doublegate_3n_4plus_bg_random true
hist_doublegate_3n_4plus_bg_random_232Th_134Te = file_232Th.Get('time_isomer_doublegate_3n_4plus_bg_random_134Te')
x_bins = hist_doublegate_3n_4plus_bg_random_232Th_134Te.GetNbinsX()

x_doublegate_3n_4plus_bg_random_232Th_134Te = np.zeros(x_bins)
y_doublegate_3n_4plus_bg_random_232Th_134Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_3n_4plus_bg_random_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_random_232Th_134Te.GetBinCenter(i+1)
    y_doublegate_3n_4plus_bg_random_232Th_134Te[i] = hist_doublegate_3n_4plus_bg_random_232Th_134Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_3n_4plus_bg_random_232Th_134Te = 2*x_doublegate_3n_4plus_bg_random_232Th_134Te

x_doublegate_3n_4plus_bg_random_232Th_134Te_long = x_doublegate_3n_4plus_bg_random_232Th_134Te
y_doublegate_3n_4plus_bg_random_232Th_134Te_long = y_doublegate_3n_4plus_bg_random_232Th_134Te

x_doublegate_3n_4plus_bg_random_232Th_134Te = x_doublegate_3n_4plus_bg_random_232Th_134Te[bin_lower:bin_upper]
y_doublegate_3n_4plus_bg_random_232Th_134Te = y_doublegate_3n_4plus_bg_random_232Th_134Te[bin_lower:bin_upper]

#print("* 232Th - 134Te 3n_4plus_bg_random fit range %d - %.d" % (x_doublegate_3n_4plus_bg_random_232Th_134Te[0],x_doublegate_3n_4plus_bg_random_232Th_134Te[-1]))




################   238U lowE_highT -  135Te   #################

#Define lowE_highT and upper fit limit
x_lower = 320
x_upper = 600

bin_lower = x_lower//2
bin_upper = x_upper//2

#Doublegate true
hist_doublegate_238U_lowE_highT_135Te = file_238U_lowE_highT.Get('time_isomer_doublegate_135Te')
x_bins = hist_doublegate_238U_lowE_highT_135Te.GetNbinsX()

x_doublegate_238U_lowE_highT_135Te = np.zeros(x_bins)
y_doublegate_238U_lowE_highT_135Te = np.zeros(x_bins)
    
for i in range(x_bins):
    x_doublegate_238U_lowE_highT_135Te[i] = hist_doublegate_238U_lowE_highT_135Te.GetBinCenter(i+1)
    y_doublegate_238U_lowE_highT_135Te[i] = hist_doublegate_238U_lowE_highT_135Te.GetBinContent(i+1)

#Make up for 2ns bins
x_doublegate_238U_lowE_highT_135Te = 2*x_doublegate_238U_lowE_highT_135Te

x_doublegate_238U_lowE_highT_135Te_long = x_doublegate_238U_lowE_highT_135Te
y_doublegate_238U_lowE_highT_135Te_long = y_doublegate_238U_lowE_highT_135Te

x_doublegate_238U_lowE_highT_135Te = x_doublegate_238U_lowE_highT_135Te[bin_lower:bin_upper]
y_doublegate_238U_lowE_highT_135Te = y_doublegate_238U_lowE_highT_135Te[bin_lower:bin_upper]

#print("* 135Te fit range %d - %.d" % (x_doublegate_238U_lowE_highT_135Te[0],x_doublegate_238U_lowE_highT_135Te[-1]))



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



def sum_smeared_exp_two_gauss_const_bg(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    return amplitude_gauss1*np.exp(-(x-mean)**2/(2*sigma1**2)) \
    + amplitude_gauss2*np.exp(-(x-mean)**2/(2*sigma2**2)) \
    + gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma1)  \
    + const_bg

def gauss_1(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    return amplitude_gauss1*np.exp(-(x-mean)**2/(2*sigma1**2))

def gauss_2(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """Gaussian component (prompt production and prompt decay)"""
    return amplitude_gauss2*np.exp(-(x-mean)**2/(2*sigma2**2))

def smeared_exp_decay_(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    """ Exponential decay """
    return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma1)

def const_bg_(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
    return const_bg + x*0

####################################################
## 		             Fit data          		      ##  START_FIT 
####################################################


################   238U lowE_highT -  134Te   #################

mean_lower = 300
mean_upper = 400
sigma1_lower = 0
sigma1_upper = 40
amplitude_gauss1_lower = 0
amplitude_gauss1_upper = 10000
sigma2_lower = 0
sigma2_upper = 40
amplitude_gauss2_lower = 0
amplitude_gauss2_upper = 10000
const_bg_lower = 0
const_bg_upper = 1000
amplitude_exp_decay_lower = 0
amplitude_exp_decay_upper = 5000
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

#mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):

#P_double_238U_lowE_highT_134Te, cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_highT_134Te, y_doublegate_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
P_double_238U_lowE_highT_134Te, cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_two_gauss_const_bg, x_doublegate_238U_lowE_highT_134Te, y_doublegate_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower, sigma1_lower, amplitude_gauss1_lower, sigma2_lower, amplitude_gauss2_lower, const_bg_lower, amplitude_exp_decay_lower, tau_decay_lower],[mean_upper, sigma1_upper, amplitude_gauss1_upper, sigma2_upper, amplitude_gauss2_upper, const_bg_upper, amplitude_exp_decay_upper, tau_decay_upper]), absolute_sigma = False)

#print("* 238U lowE_highT - 134Te Using uncertainty-weighted fit")

P_double_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[0], P_double_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma1:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[1], P_double_unc_238U_lowE_highT_134Te[1], sigma1_lower, sigma1_upper))
# print("amplitude_gauss1:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[2], P_double_unc_238U_lowE_highT_134Te[2], amplitude_gauss1_lower, amplitude_gauss1_upper))
# print("sigma2:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[3], P_double_unc_238U_lowE_highT_134Te[3], sigma2_lower, sigma2_upper))
# print("amplitude_gauss2:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[4], P_double_unc_238U_lowE_highT_134Te[4], amplitude_gauss2_lower, amplitude_gauss2_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[5], P_double_unc_238U_lowE_highT_134Te[5], const_bg_lower, const_bg_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[6], P_double_unc_238U_lowE_highT_134Te[6], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_134Te[7]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

# ################   238U lowE_highT -  134Te 1n  #################

mean_lower = 300
mean_upper = 400
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

P_double_1n_238U_lowE_highT_134Te, cov_double_1n_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_238U_lowE_highT_134Te, y_doublegate_1n_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_1n_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_1n_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te 1n Using uncertainty-weighted fit")

P_double_1n_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_1n_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te 1n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U lowE_highT -  134Te 3n  #################

mean_lower = 300
mean_upper = 400
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

P_double_3n_238U_lowE_highT_134Te, cov_double_3n_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_238U_lowE_highT_134Te, y_doublegate_3n_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_3n_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te 3n Using uncertainty-weighted fit")

P_double_3n_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_3n_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te 3n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U lowE_highT -  134Te 3n_2plus  #################

mean_lower = 300
mean_upper = 400
sigma_lower = 3
sigma_upper = 40
const_bg_lower = 0
const_bg_upper = 100
amplitude_gauss_lower = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lower = 30
amplitude_exp_decay_upper = 100
tau_decay_lower = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_3n_2plus_238U_lowE_highT_134Te, cov_double_3n_2plus_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_2plus_238U_lowE_highT_134Te, y_doublegate_3n_2plus_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_2plus_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_3n_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te 3n_2plus Using uncertainty-weighted fit")

P_double_3n_2plus_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_3n_2plus_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te 3n_2plus:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U lowE_highT -  134Te 3n_4plus  #################

mean_lower = 300
mean_upper = 400
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

P_double_3n_4plus_238U_lowE_highT_134Te, cov_double_3n_4plus_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_4plus_238U_lowE_highT_134Te, y_doublegate_3n_4plus_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_4plus_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te 3n_4plus Using uncertainty-weighted fit")

P_double_3n_4plus_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_3n_4plus_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te 3n_4plus:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U lowE_highT -  134Te 3n_6plus  #################

mean_lower = 300
mean_upper = 400
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

P_double_3n_6plus_238U_lowE_highT_134Te, cov_double_3n_6plus_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_6plus_238U_lowE_highT_134Te, y_doublegate_3n_6plus_238U_lowE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_6plus_all_238U_lowE_highT_134Te, data_bg_ridge=y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te, data_bg_random=y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_highT - 134Te 3n_6plus Using uncertainty-weighted fit")

P_double_3n_6plus_unc_238U_lowE_highT_134Te = np.sqrt(np.diag(cov_double_3n_6plus_238U_lowE_highT_134Te))

# print("\n")
# print(" ***** 238U lowE_highT - 134Te 3n_6plus:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_unc_238U_lowE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_unc_238U_lowE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_unc_238U_lowE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_unc_238U_lowE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_unc_238U_lowE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_6plus_238U_lowE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   238U lowE_lowT -  134Te   #################

mean_lowE_lowT = 0
mean_upper = 700
sigma_lowE_lowT = 0
sigma_upper = 40
const_bg_lowE_lowT = 0
const_bg_upper = 1000
amplitude_gauss_lowE_lowT = 0
amplitude_gauss_upper = 10000
amplitude_exp_decay_lowE_lowT = 0
amplitude_exp_decay_upper = 5000
tau_decay_lowE_lowT = tau_134Te
tau_decay_upper = tau_134Te+0.0001

P_double_238U_lowE_lowT_134Te, cov_double_238U_lowE_lowT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_lowT_134Te, y_doublegate_238U_lowE_lowT_134Te, sigma=sigma_data_doublegate(y_doublegate_238U_lowE_lowT_134Te), bounds=([mean_lowE_lowT,sigma_lowE_lowT,const_bg_lowE_lowT,amplitude_gauss_lowE_lowT,amplitude_exp_decay_lowE_lowT,tau_decay_lowE_lowT],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U lowE_lowT - 134Te Using uncertainty-weighted fit")
# P_double_238U_lowE_lowT_134Te, cov_double_238U_lowE_lowT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_lowT_134Te, y_doublegate_238U_lowE_lowT_134Te, bounds=([mean_lowE_lowT,sigma_lowE_lowT,const_bg_lowE_lowT,amplitude_gauss_lowE_lowT,amplitude_exp_decay_lowE_lowT,tau_decay_lowE_lowT],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
# print("* Not using uncertainty-weighted fit")

P_double_unc_238U_lowE_lowT_134Te = np.sqrt(np.diag(cov_double_238U_lowE_lowT_134Te))

# print("\n")
# print(" ***** 238U lowE_lowT - 134Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[0], P_double_unc_238U_lowE_lowT_134Te[0], mean_lowE_lowT, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[1], P_double_unc_238U_lowE_lowT_134Te[1], sigma_lowE_lowT, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[2], P_double_unc_238U_lowE_lowT_134Te[2], const_bg_lowE_lowT, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[3], P_double_unc_238U_lowE_lowT_134Te[3], amplitude_gauss_lowE_lowT, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[4], P_double_unc_238U_lowE_lowT_134Te[4], amplitude_exp_decay_lowE_lowT, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_lowT_134Te[5]*np.log(2), 0, tau_decay_lowE_lowT*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")



################   238U highE_highT -  134Te   #################

mean_lower = 300
mean_upper = 400
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

P_double_238U_highE_highT_134Te, cov_double_238U_highE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_highE_highT_134Te, y_doublegate_238U_highE_highT_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_highE_highT_134Te, data_bg_ridge=y_doublegate_bg_ridge_238U_highE_highT_134Te, data_bg_random=y_doublegate_bg_random_238U_highE_highT_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 238U highE_highT - 134Te Using uncertainty-weighted fit")

P_double_unc_238U_highE_highT_134Te = np.sqrt(np.diag(cov_double_238U_highE_highT_134Te))

# print("\n")
# print(" ***** 238U highE_highT - 134Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_highE_highT_134Te[0], P_double_unc_238U_highE_highT_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_highE_highT_134Te[1], P_double_unc_238U_highE_highT_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_highE_highT_134Te[2], P_double_unc_238U_highE_highT_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_highE_highT_134Te[3], P_double_unc_238U_highE_highT_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_highE_highT_134Te[4], P_double_unc_238U_highE_highT_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_highE_highT_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   232Th -  134Te   #################

mean_lower = 300
mean_upper = 400
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

P_double_232Th_134Te, cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_232Th_134Te, y_doublegate_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_232Th_134Te, data_bg_ridge=y_doublegate_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th - 134Te Using uncertainty-weighted fit")

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


################   232Th -  134Te 1n  #################

# mean_lower = 341
# mean_upper = 342
# sigma_lower = 4
# sigma_upper = 6
# const_bg_lower = 0
# const_bg_upper = 10
# amplitude_gauss_lower = 550
# amplitude_gauss_upper = 10000
# amplitude_exp_decay_lower = 150
# amplitude_exp_decay_upper = 1000
# tau_decay_lower = tau_134Te
# tau_decay_upper = tau_134Te+0.0001

mean_lower = 300
mean_upper = 400
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

P_double_1n_232Th_134Te, cov_double_1n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_232Th_134Te, y_doublegate_1n_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_1n_all_232Th_134Te, data_bg_ridge=y_doublegate_1n_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_1n_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th - 134Te 1n Using uncertainty-weighted fit")

P_double_1n_unc_232Th_134Te = np.sqrt(np.diag(cov_double_1n_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 1n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_1n_232Th_134Te[0], P_double_1n_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_232Th_134Te[1], P_double_1n_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_232Th_134Te[2], P_double_1n_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_232Th_134Te[3], P_double_1n_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_232Th_134Te[4], P_double_1n_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_1n_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

################   232Th -  134Te 3n  #################

mean_lower = 300
mean_upper = 400
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

P_double_3n_232Th_134Te, cov_double_3n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_232Th_134Te, y_doublegate_3n_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_all_232Th_134Te, data_bg_ridge=y_doublegate_3n_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_3n_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th - 134Te 3n Using uncertainty-weighted fit")

P_double_3n_unc_232Th_134Te = np.sqrt(np.diag(cov_double_3n_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 3n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_232Th_134Te[0], P_double_3n_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_232Th_134Te[1], P_double_3n_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_232Th_134Te[2], P_double_3n_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_232Th_134Te[3], P_double_3n_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_232Th_134Te[4], P_double_3n_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")

################   232Th -  134Te 4n  #################

# mean_lower = 341
# mean_upper = 342
# sigma_lower = 4
# sigma_upper = 6
# const_bg_lower = 0
# const_bg_upper = 10
# amplitude_gauss_lower = 1500
# amplitude_gauss_upper = 10000
# amplitude_exp_decay_lower = 100
# amplitude_exp_decay_upper = 1000
# tau_decay_lower = tau_134Te
# tau_decay_upper = tau_134Te+0.0001

mean_lower = 300
mean_upper = 400
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

P_double_4n_232Th_134Te, cov_double_4n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_4n_232Th_134Te, y_doublegate_4n_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_4n_all_232Th_134Te, data_bg_ridge=y_doublegate_4n_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_4n_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th - 134Te 4n Using uncertainty-weighted fit")

P_double_4n_unc_232Th_134Te = np.sqrt(np.diag(cov_double_4n_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 4n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_4n_232Th_134Te[0], P_double_4n_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_4n_232Th_134Te[1], P_double_4n_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_4n_232Th_134Te[2], P_double_4n_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_4n_232Th_134Te[3], P_double_4n_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_4n_232Th_134Te[4], P_double_4n_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_4n_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   232Th -  134Te 5n  #################

mean_lower = 300
mean_upper = 400
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

P_double_5n_232Th_134Te, cov_double_5n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_5n_232Th_134Te, y_doublegate_5n_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_5n_all_232Th_134Te, data_bg_ridge=y_doublegate_5n_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_5n_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th - 134Te 5n Using uncertainty-weighted fit")

P_double_5n_unc_232Th_134Te = np.sqrt(np.diag(cov_double_5n_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 5n:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_5n_232Th_134Te[0], P_double_5n_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_5n_232Th_134Te[1], P_double_5n_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_5n_232Th_134Te[2], P_double_5n_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_5n_232Th_134Te[3], P_double_5n_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_5n_232Th_134Te[4], P_double_5n_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_5n_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")



################   232Th -  134Te 3n_2plus  #################

mean_lower = 300
mean_upper = 400
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

P_double_3n_2plus_232Th_134Te, cov_double_3n_2plus_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_2plus_232Th_134Te, y_doublegate_3n_2plus_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_2plus_all_232Th_134Te, data_bg_ridge=y_doublegate_3n_2plus_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_3n_2plus_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th lowE_highT - 134Te 3n_2plus Using uncertainty-weighted fit")


P_double_3n_2plus_unc_232Th_134Te = np.sqrt(np.diag(cov_double_3n_2plus_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 3n_2plus:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_2plus_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")


################   232Th -  134Te 3n_4plus  #################

# mean_lower = 300
# mean_upper = 400
# sigma_lower = 3
# sigma_upper = 40
# const_bg_lower = 0
# const_bg_upper = 100
# amplitude_gauss_lower = 0
# amplitude_gauss_upper = 4000
# amplitude_exp_decay_lower = 100
# amplitude_exp_decay_upper = 1000
# tau_decay_lower = tau_134Te
# tau_decay_upper = tau_134Te+0.0001

mean_lower = 300
mean_upper = 400
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

P_double_3n_4plus_232Th_134Te, cov_double_3n_4plus_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_4plus_232Th_134Te, y_doublegate_3n_4plus_232Th_134Te, sigma=sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_4plus_all_232Th_134Te, data_bg_ridge=y_doublegate_3n_4plus_bg_ridge_232Th_134Te, data_bg_random=y_doublegate_3n_4plus_bg_random_232Th_134Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("* 232Th lowE_highT - 134Te 3n_4plus Using uncertainty-weighted fit")

P_double_3n_4plus_unc_232Th_134Te = np.sqrt(np.diag(cov_double_3n_4plus_232Th_134Te))

# print("\n")
# print(" ***** 232Th - 134Te 3n_4plus:  doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_unc_232Th_134Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_unc_232Th_134Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_unc_232Th_134Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_unc_232Th_134Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_unc_232Th_134Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_3n_4plus_232Th_134Te[5]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")




################   238U -  135Te   #################

mean_lower = 300
mean_upper = 400
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

#P_double_238U_lowE_highT_135Te, cov_double = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_135Te, y_doublegate_135Te, sigma=sigma_data_doublegate(y_doublegate_all_135Te, y_doublegate_bg_ridge_135Te, y_doublegate_bg_random_135Te), bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)
#print("\n Using uncertainty-weighted fit")
P_double_238U_lowE_highT_135Te, cov_double_238U_lowE_highT_135Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_highT_135Te, y_doublegate_238U_lowE_highT_135Te, bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]))
#print("* Not using uncertainty-weighted fit")

P_double_unc_238U_lowE_highT_135Te = np.sqrt(np.diag(cov_double_238U_lowE_highT_135Te))

# print("\n")
# print(" ***** 238U - 135Te:  Doublegate true spectrum fit ***** ")
# print("          -- GAUSS + SMEARED EXP + CONST_BG FIT --   ")
# print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[0], P_double_unc_238U_lowE_highT_135Te[0], mean_lower, mean_upper))
# print("sigma:                    %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[1], P_double_unc_238U_lowE_highT_135Te[1], sigma_lower, sigma_upper))
# print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[2], P_double_unc_238U_lowE_highT_135Te[2], const_bg_lower, const_bg_upper))
# print("amplitude_gauss:          %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[3], P_double_unc_238U_lowE_highT_135Te[3], amplitude_gauss_lower, amplitude_gauss_upper))
# print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[4], P_double_unc_238U_lowE_highT_135Te[4], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
# print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (P_double_238U_lowE_highT_135Te[5]*np.log(2), sigma_tau_134Te*np.log(2), tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
# print("\n")



####################################################
## 		     		 Find IYR 	                  ## 
####################################################


################   238U lowE_highT -  134Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_134Te = np.linspace(0,100*round(tau_134Te),10000)

area_double_true_238U_lowE_highT_134Te = np.trapz(gauss_1(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7])+gauss_2(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]) + smeared_exp_decay_(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), x_arr_134Te)
area_double_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss_1(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7])+gauss_2(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), x_arr_134Te)
area_double_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay_(x_arr_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), x_arr_134Te)

IYR_double_238U_lowE_highT_134Te = IYR(prompt=area_double_true_prompt_238U_lowE_highT_134Te, delayed=area_double_true_delayed_238U_lowE_highT_134Te)

print("Area all: %.2f " % area_double_true_238U_lowE_highT_134Te)
print("Area prompt: %.2f " % area_double_true_prompt_238U_lowE_highT_134Te)
print("Area delayed: %.2f " % area_double_true_delayed_238U_lowE_highT_134Te)
print("IYR: %.3f" % IYR_double_238U_lowE_highT_134Te)



################   238U lowE_highT -  134Te 1n  #################

area_double_1n_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_1n_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_1n_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)

IYR_double_1n_238U_lowE_highT_134Te = IYR(prompt=area_double_1n_true_prompt_238U_lowE_highT_134Te, delayed=area_double_1n_true_delayed_238U_lowE_highT_134Te)



################   238U lowE_highT -  134Te 3n  #################

area_double_3n_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)

IYR_double_3n_238U_lowE_highT_134Te = IYR(prompt=area_double_3n_true_prompt_238U_lowE_highT_134Te, delayed=area_double_3n_true_delayed_238U_lowE_highT_134Te)



################   238U lowE_highT -  134Te 3n_2plus  #################

area_double_3n_2plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_2plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_2plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

IYR_double_3n_2plus_238U_lowE_highT_134Te = IYR(prompt=area_double_3n_2plus_true_prompt_238U_lowE_highT_134Te, delayed=area_double_3n_2plus_true_delayed_238U_lowE_highT_134Te)


################   238U lowE_highT -  134Te 3n_4plus  #################

area_double_3n_4plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_4plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_4plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

IYR_double_3n_4plus_238U_lowE_highT_134Te = IYR(prompt=area_double_3n_4plus_true_prompt_238U_lowE_highT_134Te, delayed=area_double_3n_4plus_true_delayed_238U_lowE_highT_134Te)


################   238U lowE_highT -  134Te 3n_6plus  #################

area_double_3n_6plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_6plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
area_double_3n_6plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

IYR_double_3n_6plus_238U_lowE_highT_134Te = IYR(prompt=area_double_3n_6plus_true_prompt_238U_lowE_highT_134Te, delayed=area_double_3n_6plus_true_delayed_238U_lowE_highT_134Te)


################   238U lowE_lowT -  134Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_134Te = np.linspace(0,40*round(tau_134Te),5000)

area_double_true_238U_lowE_lowT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), x_arr_134Te)
area_double_true_prompt_238U_lowE_lowT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), x_arr_134Te)
area_double_true_delayed_238U_lowE_lowT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), x_arr_134Te)

IYR_double_238U_lowE_lowT_134Te = IYR(prompt=area_double_true_prompt_238U_lowE_lowT_134Te, delayed=area_double_true_delayed_238U_lowE_lowT_134Te)



################   238U highE -  134Te   #################

area_double_true_238U_highE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), x_arr_134Te)
area_double_true_prompt_238U_highE_highT_134Te = np.trapz(gauss(x_arr_134Te, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), x_arr_134Te)
area_double_true_delayed_238U_highE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), x_arr_134Te)

IYR_double_238U_highE_highT_134Te = IYR(prompt=area_double_true_prompt_238U_highE_highT_134Te, delayed=area_double_true_delayed_238U_highE_highT_134Te)



################   232Th -  134Te   #################

area_double_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
area_double_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)
area_double_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_232Th_134Te[0], P_double_232Th_134Te[1], P_double_232Th_134Te[2], P_double_232Th_134Te[3], P_double_232Th_134Te[4], P_double_232Th_134Te[5]), x_arr_134Te)

IYR_double_232Th_134Te = IYR(prompt=area_double_true_prompt_232Th_134Te, delayed=area_double_true_delayed_232Th_134Te)


################   232Th -  134Te 1n  #################

area_double_1n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), x_arr_134Te)
area_double_1n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), x_arr_134Te)
area_double_1n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), x_arr_134Te)

IYR_double_1n_232Th_134Te = IYR(prompt=area_double_1n_true_prompt_232Th_134Te, delayed=area_double_1n_true_delayed_232Th_134Te)


################   232Th -  134Te 3n  #################

area_double_3n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), x_arr_134Te)
area_double_3n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), x_arr_134Te)
area_double_3n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), x_arr_134Te)

IYR_double_3n_232Th_134Te = IYR(prompt=area_double_3n_true_prompt_232Th_134Te, delayed=area_double_3n_true_delayed_232Th_134Te)


################   232Th -  134Te 4n  #################

area_double_4n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), x_arr_134Te)
area_double_4n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), x_arr_134Te)
area_double_4n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), x_arr_134Te)

IYR_double_4n_232Th_134Te = IYR(prompt=area_double_4n_true_prompt_232Th_134Te, delayed=area_double_4n_true_delayed_232Th_134Te)


################   232Th -  134Te 5n  #################

area_double_5n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), x_arr_134Te)
area_double_5n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), x_arr_134Te)
area_double_5n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), x_arr_134Te)

IYR_double_5n_232Th_134Te = IYR(prompt=area_double_5n_true_prompt_232Th_134Te, delayed=area_double_5n_true_delayed_232Th_134Te)


################   232Th -  134Te 3n_2plus  #################

area_double_3n_2plus_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)
area_double_3n_2plus_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)
area_double_3n_2plus_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)

IYR_double_3n_2plus_232Th_134Te = IYR(prompt=area_double_3n_2plus_true_prompt_232Th_134Te, delayed=area_double_3n_2plus_true_delayed_232Th_134Te)


################   232Th -  134Te 3n_4plus  #################

area_double_3n_4plus_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)
area_double_3n_4plus_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)
area_double_3n_4plus_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)

IYR_double_3n_4plus_232Th_134Te = IYR(prompt=area_double_3n_4plus_true_prompt_232Th_134Te, delayed=area_double_3n_4plus_true_delayed_232Th_134Te)


################   238U -  135Te   #################

#NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
x_arr_135Te = np.linspace(0,40*round(tau_135Te),40*round(tau_135Te))

area_double_true_238U_lowE_highT_135Te = np.trapz(gauss(x_arr_135Te, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]) + smeared_exp_decay(x_arr_135Te, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), x_arr_135Te)
area_double_true_prompt_238U_lowE_highT_135Te = np.trapz(gauss(x_arr_135Te, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), x_arr_135Te)
area_double_true_delayed_238U_lowE_highT_135Te = np.trapz(smeared_exp_decay(x_arr_135Te, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), x_arr_135Te)

IYR_double_238U_lowE_highT_135Te = IYR(prompt=area_double_true_prompt_238U_lowE_highT_135Te, delayed=area_double_true_delayed_238U_lowE_highT_135Te)



####################################################
##                BOOTSTRAPPING                   ## START_BOOTSTRAP
####################################################

if BOOTSTRAP==True:

    N_BOOTSTRAP = 1000

    print("\n Starting BOOTSTRAPPING... Number of iterations: %.d" % N_BOOTSTRAP)

    ### Find uncertainty on data set ###

    #238U lowE_highT - 134Te
    unc_y_doublegate_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te_long)

    #238U lowE_highT - 134Te 1n
    unc_y_doublegate_1n_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_1n_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_1n_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_1n_bg_random_238U_lowE_highT_134Te_long)

    #238U lowE_highT - 134Te 3n
    unc_y_doublegate_3n_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_3n_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_3n_bg_random_238U_lowE_highT_134Te_long)

    #238U lowE_highT - 134Te 3n_2plus
    unc_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_3n_2plus_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_3n_2plus_bg_random_238U_lowE_highT_134Te_long)

    #238U lowE_highT - 134Te 3n_4plus
    unc_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_3n_4plus_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_3n_4plus_bg_random_238U_lowE_highT_134Te_long)

    #238U lowE_highT - 134Te 3n_6plus
    unc_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_3n_6plus_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_3n_6plus_bg_random_238U_lowE_highT_134Te_long)

    #238U highE_highT - 134Te
    unc_y_doublegate_238U_highE_highT_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_highE_highT_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_238U_highE_highT_134Te_long, data_bg_random=y_doublegate_bg_random_238U_highE_highT_134Te_long)


    #232Th - 134Te
    unc_y_doublegate_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_all_232Th_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_bg_random_232Th_134Te_long)

    #232Th - 134Te 1n
    unc_y_doublegate_1n_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_1n_all_232Th_134Te_long, data_bg_ridge=y_doublegate_1n_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_1n_bg_random_232Th_134Te_long)

    #232Th - 134Te 3n
    unc_y_doublegate_3n_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_all_232Th_134Te_long, data_bg_ridge=y_doublegate_3n_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_3n_bg_random_232Th_134Te_long)

    #232Th - 134Te 4n
    unc_y_doublegate_4n_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_4n_all_232Th_134Te_long, data_bg_ridge=y_doublegate_4n_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_4n_bg_random_232Th_134Te_long)

    #232Th - 134Te 5n
    unc_y_doublegate_5n_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_5n_all_232Th_134Te_long, data_bg_ridge=y_doublegate_5n_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_5n_bg_random_232Th_134Te_long)

    #232Th - 134Te 3n_2plus
    unc_y_doublegate_3n_2plus_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_2plus_all_232Th_134Te_long, data_bg_ridge=y_doublegate_3n_2plus_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_3n_2plus_bg_random_232Th_134Te_long)

    #232Th - 134Te 3n_4plus
    unc_y_doublegate_3n_4plus_232Th_134Te_long = sigma_data_doublegate_all_bg(data_all=y_doublegate_3n_4plus_all_232Th_134Te_long, data_bg_ridge=y_doublegate_3n_4plus_bg_ridge_232Th_134Te_long, data_bg_random=y_doublegate_3n_4plus_bg_random_232Th_134Te_long)


    ### Arrays to store resampled data ###

    #238U lowE_highT - 134Te
    resampled_y_doublegate_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_238U_lowE_highT_134Te_long))
    resampled_IYR_array_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U lowE_highT - 134Te 1n
    resampled_y_doublegate_1n_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_1n_238U_lowE_highT_134Te_long))
    resampled_IYR_array_1n_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U lowE_highT - 134Te 3n
    resampled_y_doublegate_3n_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_3n_238U_lowE_highT_134Te_long))
    resampled_IYR_array_3n_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U lowE_highT - 134Te 3n_2plus
    resampled_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_3n_2plus_238U_lowE_highT_134Te_long))
    resampled_IYR_array_3n_2plus_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U lowE_highT - 134Te 3n_4plus
    resampled_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_3n_4plus_238U_lowE_highT_134Te_long))
    resampled_IYR_array_3n_4plus_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U lowE_highT - 134Te 3n_6plus
    resampled_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long = np.zeros(len(y_doublegate_3n_6plus_238U_lowE_highT_134Te_long))
    resampled_IYR_array_3n_6plus_238U_lowE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #238U highE_highT - 134Te
    resampled_y_doublegate_238U_highE_highT_134Te_long = np.zeros(len(y_doublegate_238U_highE_highT_134Te_long))
    resampled_IYR_array_238U_highE_highT_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te
    resampled_y_doublegate_232Th_134Te_long = np.zeros(len(y_doublegate_232Th_134Te_long))
    resampled_IYR_array_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 1n
    resampled_y_doublegate_1n_232Th_134Te_long = np.zeros(len(y_doublegate_1n_232Th_134Te_long))
    resampled_IYR_array_1n_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 3n
    resampled_y_doublegate_3n_232Th_134Te_long = np.zeros(len(y_doublegate_3n_232Th_134Te_long))
    resampled_IYR_array_3n_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 4n
    resampled_y_doublegate_4n_232Th_134Te_long = np.zeros(len(y_doublegate_4n_232Th_134Te_long))
    resampled_IYR_array_4n_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 5n
    resampled_y_doublegate_5n_232Th_134Te_long = np.zeros(len(y_doublegate_5n_232Th_134Te_long))
    resampled_IYR_array_5n_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 3n_2plus
    resampled_y_doublegate_3n_2plus_232Th_134Te_long = np.zeros(len(y_doublegate_3n_2plus_232Th_134Te_long))
    resampled_IYR_array_3n_2plus_232Th_134Te = np.zeros(N_BOOTSTRAP)

    #232Th - 134Te 3n_4plus
    resampled_y_doublegate_3n_4plus_232Th_134Te_long = np.zeros(len(y_doublegate_3n_4plus_232Th_134Te_long))
    resampled_IYR_array_3n_4plus_232Th_134Te = np.zeros(N_BOOTSTRAP)



    for n in range(N_BOOTSTRAP):


################   238U lowE_highT -  134Te   #################

        y_doublegate_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_all_238U_lowE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_all_238U_lowE_highT_134Te_long-y_doublegate_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)    

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
        amplitude_gauss_upper = 10000
        amplitude_exp_decay_lower = 0
        amplitude_exp_decay_upper = 5000
        tau_decay_lower = tau_134Te
        tau_decay_upper = tau_134Te+0.0001

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,10))

        bin_lower = x_lower//2
        bin_upper = x_upper//2

        try:        

            resampled_P_double_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_lowE_highT_134Te[0], resampled_P_double_238U_lowE_highT_134Te[1], resampled_P_double_238U_lowE_highT_134Te[2], resampled_P_double_238U_lowE_highT_134Te[3], resampled_P_double_238U_lowE_highT_134Te[4], resampled_P_double_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_238U_lowE_highT_134Te[n] = IYR_double_238U_lowE_highT_134Te  

#################   238U lowE_highT 1n -  134Te   #################

        y_doublegate_1n_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_1n_all_238U_lowE_highT_134Te_long))
       
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_1n_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_1n_all_238U_lowE_highT_134Te_long-y_doublegate_1n_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)

        for i in range (len(y_doublegate_1n_all_238U_lowE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_1n_238U_lowE_highT_134Te_long[i] = y_doublegate_1n_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_1n_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        


        try:
            resampled_P_double_1n_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_1n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_1n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_1n_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_1n_238U_lowE_highT_134Te[0], resampled_P_double_1n_238U_lowE_highT_134Te[1], resampled_P_double_1n_238U_lowE_highT_134Te[2], resampled_P_double_1n_238U_lowE_highT_134Te[3], resampled_P_double_1n_238U_lowE_highT_134Te[4], resampled_P_double_1n_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_1n_238U_lowE_highT_134Te[0], resampled_P_double_1n_238U_lowE_highT_134Te[1], resampled_P_double_1n_238U_lowE_highT_134Te[2], resampled_P_double_1n_238U_lowE_highT_134Te[3], resampled_P_double_1n_238U_lowE_highT_134Te[4], resampled_P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_1n_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_1n_238U_lowE_highT_134Te[0], resampled_P_double_1n_238U_lowE_highT_134Te[1], resampled_P_double_1n_238U_lowE_highT_134Te[2], resampled_P_double_1n_238U_lowE_highT_134Te[3], resampled_P_double_1n_238U_lowE_highT_134Te[4], resampled_P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_1n_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_1n_238U_lowE_highT_134Te[0], resampled_P_double_1n_238U_lowE_highT_134Te[1], resampled_P_double_1n_238U_lowE_highT_134Te[2], resampled_P_double_1n_238U_lowE_highT_134Te[3], resampled_P_double_1n_238U_lowE_highT_134Te[4], resampled_P_double_1n_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_1n_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_1n_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_1n_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_1n_238U_lowE_highT_134Te[n] = IYR_double_1n_238U_lowE_highT_134Te
            print("Error - curve_fit failed")


#################   238U lowE_highT 3n -  134Te   #################

        y_doublegate_3n_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_all_238U_lowE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_3n_all_238U_lowE_highT_134Te_long-y_doublegate_3n_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)

        for i in range (len(y_doublegate_3n_all_238U_lowE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_238U_lowE_highT_134Te_long[i] = y_doublegate_3n_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2 

        try: 
            resampled_P_double_3n_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_3n_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_238U_lowE_highT_134Te[0], resampled_P_double_3n_238U_lowE_highT_134Te[1], resampled_P_double_3n_238U_lowE_highT_134Te[2], resampled_P_double_3n_238U_lowE_highT_134Te[3], resampled_P_double_3n_238U_lowE_highT_134Te[4], resampled_P_double_3n_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_238U_lowE_highT_134Te[0], resampled_P_double_3n_238U_lowE_highT_134Te[1], resampled_P_double_3n_238U_lowE_highT_134Te[2], resampled_P_double_3n_238U_lowE_highT_134Te[3], resampled_P_double_3n_238U_lowE_highT_134Te[4], resampled_P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_238U_lowE_highT_134Te[0], resampled_P_double_3n_238U_lowE_highT_134Te[1], resampled_P_double_3n_238U_lowE_highT_134Te[2], resampled_P_double_3n_238U_lowE_highT_134Te[3], resampled_P_double_3n_238U_lowE_highT_134Te[4], resampled_P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_238U_lowE_highT_134Te[0], resampled_P_double_3n_238U_lowE_highT_134Te[1], resampled_P_double_3n_238U_lowE_highT_134Te[2], resampled_P_double_3n_238U_lowE_highT_134Te[3], resampled_P_double_3n_238U_lowE_highT_134Te[4], resampled_P_double_3n_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_3n_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_3n_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_238U_lowE_highT_134Te[n] = IYR_double_3n_238U_lowE_highT_134Te
            print("Error - curve_fit failed")


#################   238U lowE_highT 3n_2plus -  134Te   #################

        y_doublegate_3n_2plus_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long))

        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_2plus_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long-y_doublegate_3n_2plus_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)
                    

        for i in range (len(y_doublegate_3n_2plus_all_238U_lowE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long[i] = y_doublegate_3n_2plus_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try: 
            resampled_P_double_3n_2plus_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_2plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_2plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_3n_2plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_2plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_2plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_2plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_2plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_2plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_2plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_2plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_2plus_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_3n_2plus_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_3n_2plus_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_2plus_238U_lowE_highT_134Te[n] = IYR_double_3n_2plus_238U_lowE_highT_134Te
            print("Error - curve_fit failed")

#################   238U lowE_highT 3n_4plus -  134Te   #################

        y_doublegate_3n_4plus_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_4plus_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long-y_doublegate_3n_4plus_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)   

        for i in range (len(y_doublegate_3n_4plus_all_238U_lowE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long[i] = y_doublegate_3n_4plus_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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


        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_3n_4plus_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_4plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_4plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_3n_4plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_4plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_4plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_4plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_4plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_4plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_4plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_4plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_4plus_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_3n_4plus_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_3n_4plus_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_4plus_238U_lowE_highT_134Te[n] = IYR_double_3n_4plus_238U_lowE_highT_134Te
            print("Error - curve_fit failed")

#################   238U lowE_highT 3n_6plus -  134Te   #################

        y_doublegate_3n_6plus_238U_lowE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_6plus_238U_lowE_highT_134Te_long_bgvaried = y_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long-y_doublegate_3n_6plus_bg_238U_lowE_highT_134Te_long*np.random.normal(1,0.025)

        for i in range (len(y_doublegate_3n_6plus_all_238U_lowE_highT_134Te_long)):
             #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long[i] = y_doublegate_3n_6plus_238U_lowE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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


        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_3n_6plus_238U_lowE_highT_134Te, resampled_cov_double_238U_lowE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_6plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_6plus_238U_lowE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_3n_6plus_true_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_6plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_6plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_6plus_true_prompt_238U_lowE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_6plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_6plus_true_delayed_238U_lowE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_6plus_238U_lowE_highT_134Te[0], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[1], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[2], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[3], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[4], resampled_P_double_3n_6plus_238U_lowE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_6plus_238U_lowE_highT_134Te[n] = IYR(prompt=resampled_area_double_3n_6plus_true_prompt_238U_lowE_highT_134Te, delayed=resampled_area_double_3n_6plus_true_delayed_238U_lowE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_6plus_238U_lowE_highT_134Te[n] = IYR_double_3n_6plus_238U_lowE_highT_134Te
            print("Error - curve_fit failed")


################   238U highE_highT -  134Te   #################

        y_doublegate_238U_highE_highT_134Te_long_bgvaried = np.zeros(len(y_doublegate_all_238U_highE_highT_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_238U_highE_highT_134Te_long_bgvaried = y_doublegate_all_238U_highE_highT_134Te_long-y_doublegate_bg_238U_highE_highT_134Te_long*np.random.normal(1,0.025)  

        for i in range (len(y_doublegate_all_238U_highE_highT_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_238U_highE_highT_134Te_long[i] = y_doublegate_238U_highE_highT_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_238U_highE_highT_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_238U_highE_highT_134Te, resampled_cov_double_238U_highE_highT_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_238U_highE_highT_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_238U_highE_highT_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_238U_highE_highT_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_true_238U_highE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_highE_highT_134Te[0], resampled_P_double_238U_highE_highT_134Te[1], resampled_P_double_238U_highE_highT_134Te[2], resampled_P_double_238U_highE_highT_134Te[3], resampled_P_double_238U_highE_highT_134Te[4], resampled_P_double_238U_highE_highT_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_highE_highT_134Te[0], resampled_P_double_238U_highE_highT_134Te[1], resampled_P_double_238U_highE_highT_134Te[2], resampled_P_double_238U_highE_highT_134Te[3], resampled_P_double_238U_highE_highT_134Te[4], resampled_P_double_238U_highE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_true_prompt_238U_highE_highT_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_238U_highE_highT_134Te[0], resampled_P_double_238U_highE_highT_134Te[1], resampled_P_double_238U_highE_highT_134Te[2], resampled_P_double_238U_highE_highT_134Te[3], resampled_P_double_238U_highE_highT_134Te[4], resampled_P_double_238U_highE_highT_134Te[5]), x_arr_134Te)
            resampled_area_double_true_delayed_238U_highE_highT_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_238U_highE_highT_134Te[0], resampled_P_double_238U_highE_highT_134Te[1], resampled_P_double_238U_highE_highT_134Te[2], resampled_P_double_238U_highE_highT_134Te[3], resampled_P_double_238U_highE_highT_134Te[4], resampled_P_double_238U_highE_highT_134Te[5]), x_arr_134Te)

            resampled_IYR_array_238U_highE_highT_134Te[n] = IYR(prompt=resampled_area_double_true_prompt_238U_highE_highT_134Te, delayed=resampled_area_double_true_delayed_238U_highE_highT_134Te)

        except RuntimeError:
            resampled_IYR_array_238U_highE_highT_134Te[n] = IYR_double_238U_highE_highT_134Te
            print("Error - curve_fit failed")

################   232Th -  134Te   #################

        y_doublegate_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_all_232Th_134Te_long))

        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_232Th_134Te_long_bgvaried = y_doublegate_all_232Th_134Te_long-y_doublegate_bg_232Th_134Te_long*np.random.normal(1,0.025)
                  
        for i in range (len(y_doublegate_all_232Th_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_232Th_134Te_long[i] = y_doublegate_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_232Th_134Te, resampled_cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_232Th_134Te[0], resampled_P_double_232Th_134Te[1], resampled_P_double_232Th_134Te[2], resampled_P_double_232Th_134Te[3], resampled_P_double_232Th_134Te[4], resampled_P_double_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_232Th_134Te[0], resampled_P_double_232Th_134Te[1], resampled_P_double_232Th_134Te[2], resampled_P_double_232Th_134Te[3], resampled_P_double_232Th_134Te[4], resampled_P_double_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_232Th_134Te[0], resampled_P_double_232Th_134Te[1], resampled_P_double_232Th_134Te[2], resampled_P_double_232Th_134Te[3], resampled_P_double_232Th_134Te[4], resampled_P_double_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_232Th_134Te[0], resampled_P_double_232Th_134Te[1], resampled_P_double_232Th_134Te[2], resampled_P_double_232Th_134Te[3], resampled_P_double_232Th_134Te[4], resampled_P_double_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_232Th_134Te[n] = IYR(prompt=resampled_area_double_true_prompt_232Th_134Te, delayed=resampled_area_double_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_232Th_134Te[n] = IYR_double_232Th_134Te
            print("Error - curve_fit failed")        

################   232Th -  134Te 1n  #################

        y_doublegate_1n_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_1n_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_1n_232Th_134Te_long_bgvaried = y_doublegate_1n_all_232Th_134Te_long-y_doublegate_1n_bg_232Th_134Te_long*np.random.normal(1,0.025)
           
        for i in range (len(y_doublegate_1n_all_232Th_134Te_long)):
             #Vary value of each bin within uncertainty
            resampled_y_doublegate_1n_232Th_134Te_long[i] = y_doublegate_1n_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_1n_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_1n_232Th_134Te, resampled_cov_double_1n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_1n_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_1n_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_1n_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_1n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_1n_232Th_134Te[0], resampled_P_double_1n_232Th_134Te[1], resampled_P_double_1n_232Th_134Te[2], resampled_P_double_1n_232Th_134Te[3], resampled_P_double_1n_232Th_134Te[4], resampled_P_double_1n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_1n_232Th_134Te[0], resampled_P_double_1n_232Th_134Te[1], resampled_P_double_1n_232Th_134Te[2], resampled_P_double_1n_232Th_134Te[3], resampled_P_double_1n_232Th_134Te[4], resampled_P_double_1n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_1n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_1n_232Th_134Te[0], resampled_P_double_1n_232Th_134Te[1], resampled_P_double_1n_232Th_134Te[2], resampled_P_double_1n_232Th_134Te[3], resampled_P_double_1n_232Th_134Te[4], resampled_P_double_1n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_1n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_1n_232Th_134Te[0], resampled_P_double_1n_232Th_134Te[1], resampled_P_double_1n_232Th_134Te[2], resampled_P_double_1n_232Th_134Te[3], resampled_P_double_1n_232Th_134Te[4], resampled_P_double_1n_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_1n_232Th_134Te[n] = IYR(prompt=resampled_area_double_1n_true_prompt_232Th_134Te, delayed=resampled_area_double_1n_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_1n_232Th_134Te[n] = IYR_double_1n_232Th_134Te
            print("Error - curve_fit failed")

################   232Th -  134Te 3n  #################

        y_doublegate_3n_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_232Th_134Te_long_bgvaried = y_doublegate_3n_all_232Th_134Te_long-y_doublegate_3n_bg_232Th_134Te_long*np.random.normal(1,0.025)
           
        for i in range (len(y_doublegate_3n_all_232Th_134Te_long)):
             #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_232Th_134Te_long[i] = y_doublegate_3n_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_3n_232Th_134Te, resampled_cov_double_3n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_3n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_232Th_134Te[0], resampled_P_double_3n_232Th_134Te[1], resampled_P_double_3n_232Th_134Te[2], resampled_P_double_3n_232Th_134Te[3], resampled_P_double_3n_232Th_134Te[4], resampled_P_double_3n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_232Th_134Te[0], resampled_P_double_3n_232Th_134Te[1], resampled_P_double_3n_232Th_134Te[2], resampled_P_double_3n_232Th_134Te[3], resampled_P_double_3n_232Th_134Te[4], resampled_P_double_3n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_232Th_134Te[0], resampled_P_double_3n_232Th_134Te[1], resampled_P_double_3n_232Th_134Te[2], resampled_P_double_3n_232Th_134Te[3], resampled_P_double_3n_232Th_134Te[4], resampled_P_double_3n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_232Th_134Te[0], resampled_P_double_3n_232Th_134Te[1], resampled_P_double_3n_232Th_134Te[2], resampled_P_double_3n_232Th_134Te[3], resampled_P_double_3n_232Th_134Te[4], resampled_P_double_3n_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_232Th_134Te[n] = IYR(prompt=resampled_area_double_3n_true_prompt_232Th_134Te, delayed=resampled_area_double_3n_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_232Th_134Te[n] = IYR_double_3n_232Th_134Te
            print("Error - curve_fit failed")


################   232Th -  134Te 4n  #################

        y_doublegate_4n_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_4n_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_4n_232Th_134Te_long_bgvaried = y_doublegate_4n_all_232Th_134Te_long-y_doublegate_4n_bg_232Th_134Te_long*np.random.normal(1,0.025)

        for i in range (len(y_doublegate_4n_all_232Th_134Te_long)):
             #Vary value of each bin within uncertainty
            resampled_y_doublegate_4n_232Th_134Te_long[i] = y_doublegate_4n_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_4n_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_4n_232Th_134Te, resampled_cov_double_4n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_4n_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_4n_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_4n_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_4n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_4n_232Th_134Te[0], resampled_P_double_4n_232Th_134Te[1], resampled_P_double_4n_232Th_134Te[2], resampled_P_double_4n_232Th_134Te[3], resampled_P_double_4n_232Th_134Te[4], resampled_P_double_4n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_4n_232Th_134Te[0], resampled_P_double_4n_232Th_134Te[1], resampled_P_double_4n_232Th_134Te[2], resampled_P_double_4n_232Th_134Te[3], resampled_P_double_4n_232Th_134Te[4], resampled_P_double_4n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_4n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_4n_232Th_134Te[0], resampled_P_double_4n_232Th_134Te[1], resampled_P_double_4n_232Th_134Te[2], resampled_P_double_4n_232Th_134Te[3], resampled_P_double_4n_232Th_134Te[4], resampled_P_double_4n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_4n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_4n_232Th_134Te[0], resampled_P_double_4n_232Th_134Te[1], resampled_P_double_4n_232Th_134Te[2], resampled_P_double_4n_232Th_134Te[3], resampled_P_double_4n_232Th_134Te[4], resampled_P_double_4n_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_4n_232Th_134Te[n] = IYR(prompt=resampled_area_double_4n_true_prompt_232Th_134Te, delayed=resampled_area_double_4n_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_4n_232Th_134Te[n] = IYR_double_4n_232Th_134Te
            print("Error - curve_fit failed")


################   232Th -  134Te 5n  #################

        y_doublegate_5n_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_5n_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_5n_232Th_134Te_long_bgvaried = y_doublegate_5n_all_232Th_134Te_long-y_doublegate_5n_bg_232Th_134Te_long*np.random.normal(1,0.025)    

        for i in range (len(y_doublegate_5n_all_232Th_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_5n_232Th_134Te_long[i] = y_doublegate_5n_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_5n_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_5n_232Th_134Te, resampled_cov_double_5n_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_5n_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_5n_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_5n_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_5n_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_5n_232Th_134Te[0], resampled_P_double_5n_232Th_134Te[1], resampled_P_double_5n_232Th_134Te[2], resampled_P_double_5n_232Th_134Te[3], resampled_P_double_5n_232Th_134Te[4], resampled_P_double_5n_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_5n_232Th_134Te[0], resampled_P_double_5n_232Th_134Te[1], resampled_P_double_5n_232Th_134Te[2], resampled_P_double_5n_232Th_134Te[3], resampled_P_double_5n_232Th_134Te[4], resampled_P_double_5n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_5n_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_5n_232Th_134Te[0], resampled_P_double_5n_232Th_134Te[1], resampled_P_double_5n_232Th_134Te[2], resampled_P_double_5n_232Th_134Te[3], resampled_P_double_5n_232Th_134Te[4], resampled_P_double_5n_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_5n_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_5n_232Th_134Te[0], resampled_P_double_5n_232Th_134Te[1], resampled_P_double_5n_232Th_134Te[2], resampled_P_double_5n_232Th_134Te[3], resampled_P_double_5n_232Th_134Te[4], resampled_P_double_5n_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_5n_232Th_134Te[n] = IYR(prompt=resampled_area_double_5n_true_prompt_232Th_134Te, delayed=resampled_area_double_5n_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_5n_232Th_134Te[n] = IYR_double_5n_232Th_134Te
            print("Error - curve_fit failed")


#################   232Th 3n_2plus -  134Te   #################

        y_doublegate_3n_2plus_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_2plus_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_2plus_232Th_134Te_long_bgvaried = y_doublegate_3n_2plus_all_232Th_134Te_long-y_doublegate_3n_2plus_bg_232Th_134Te_long*np.random.normal(1,0.025)
            

        for i in range (len(y_doublegate_3n_2plus_all_232Th_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_2plus_232Th_134Te_long[i] = y_doublegate_3n_2plus_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_2plus_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_3n_2plus_232Th_134Te, resampled_cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_2plus_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_2plus_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_2plus_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False, maxfev=10000)

            resampled_area_double_3n_2plus_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_2plus_232Th_134Te[0], resampled_P_double_3n_2plus_232Th_134Te[1], resampled_P_double_3n_2plus_232Th_134Te[2], resampled_P_double_3n_2plus_232Th_134Te[3], resampled_P_double_3n_2plus_232Th_134Te[4], resampled_P_double_3n_2plus_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_2plus_232Th_134Te[0], resampled_P_double_3n_2plus_232Th_134Te[1], resampled_P_double_3n_2plus_232Th_134Te[2], resampled_P_double_3n_2plus_232Th_134Te[3], resampled_P_double_3n_2plus_232Th_134Te[4], resampled_P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_2plus_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_2plus_232Th_134Te[0], resampled_P_double_3n_2plus_232Th_134Te[1], resampled_P_double_3n_2plus_232Th_134Te[2], resampled_P_double_3n_2plus_232Th_134Te[3], resampled_P_double_3n_2plus_232Th_134Te[4], resampled_P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_2plus_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_2plus_232Th_134Te[0], resampled_P_double_3n_2plus_232Th_134Te[1], resampled_P_double_3n_2plus_232Th_134Te[2], resampled_P_double_3n_2plus_232Th_134Te[3], resampled_P_double_3n_2plus_232Th_134Te[4], resampled_P_double_3n_2plus_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_2plus_232Th_134Te[n] = IYR(prompt=resampled_area_double_3n_2plus_true_prompt_232Th_134Te, delayed=resampled_area_double_3n_2plus_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_2plus_232Th_134Te[n] = IYR_double_3n_2plus_232Th_134Te
            print("Error - curve_fit failed")

#################   232Th 3n_4plus -  134Te   #################

        y_doublegate_3n_4plus_232Th_134Te_long_bgvaried = np.zeros(len(y_doublegate_3n_4plus_all_232Th_134Te_long))
        
        #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
        y_doublegate_3n_4plus_232Th_134Te_long_bgvaried = y_doublegate_3n_4plus_all_232Th_134Te_long-y_doublegate_3n_4plus_bg_232Th_134Te_long*np.random.normal(1,0.025)

        for i in range (len(y_doublegate_3n_4plus_all_232Th_134Te_long)):
            #Vary value of each bin within uncertainty
            resampled_y_doublegate_3n_4plus_232Th_134Te_long[i] = y_doublegate_3n_4plus_232Th_134Te_long_bgvaried[i] + np.random.normal(0, unc_y_doublegate_3n_4plus_232Th_134Te_long[i]) 


        #Fit resampled data
        mean_lower = 300
        mean_upper = 400
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

        #Define lower and upper fit limit
        x_lower = 320 + int(np.random.normal(0,10))
        x_upper = 640 + int(np.random.normal(0,5))

        bin_lower = x_lower//2
        bin_upper = x_upper//2        

        try:
            resampled_P_double_3n_4plus_232Th_134Te, resampled_cov_double_232Th_134Te = curve_fit(sum_smeared_exp_gauss_const_bg, x_doublegate_3n_4plus_232Th_134Te_long[bin_lower:bin_upper], resampled_y_doublegate_3n_4plus_232Th_134Te_long[bin_lower:bin_upper], sigma=unc_y_doublegate_3n_4plus_232Th_134Te_long[bin_lower:bin_upper], bounds=([mean_lower,sigma_lower,const_bg_lower,amplitude_gauss_lower,amplitude_exp_decay_lower,tau_decay_lower],[mean_upper,sigma_upper,const_bg_upper,amplitude_gauss_upper,amplitude_exp_decay_upper,tau_decay_upper]), absolute_sigma = False)

            resampled_area_double_3n_4plus_true_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_4plus_232Th_134Te[0], resampled_P_double_3n_4plus_232Th_134Te[1], resampled_P_double_3n_4plus_232Th_134Te[2], resampled_P_double_3n_4plus_232Th_134Te[3], resampled_P_double_3n_4plus_232Th_134Te[4], resampled_P_double_3n_4plus_232Th_134Te[5]) + smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_4plus_232Th_134Te[0], resampled_P_double_3n_4plus_232Th_134Te[1], resampled_P_double_3n_4plus_232Th_134Te[2], resampled_P_double_3n_4plus_232Th_134Te[3], resampled_P_double_3n_4plus_232Th_134Te[4], resampled_P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_4plus_true_prompt_232Th_134Te = np.trapz(gauss(x_arr_134Te, resampled_P_double_3n_4plus_232Th_134Te[0], resampled_P_double_3n_4plus_232Th_134Te[1], resampled_P_double_3n_4plus_232Th_134Te[2], resampled_P_double_3n_4plus_232Th_134Te[3], resampled_P_double_3n_4plus_232Th_134Te[4], resampled_P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)
            resampled_area_double_3n_4plus_true_delayed_232Th_134Te = np.trapz(smeared_exp_decay(x_arr_134Te, resampled_P_double_3n_4plus_232Th_134Te[0], resampled_P_double_3n_4plus_232Th_134Te[1], resampled_P_double_3n_4plus_232Th_134Te[2], resampled_P_double_3n_4plus_232Th_134Te[3], resampled_P_double_3n_4plus_232Th_134Te[4], resampled_P_double_3n_4plus_232Th_134Te[5]), x_arr_134Te)

            resampled_IYR_array_3n_4plus_232Th_134Te[n] = IYR(prompt=resampled_area_double_3n_4plus_true_prompt_232Th_134Te, delayed=resampled_area_double_3n_4plus_true_delayed_232Th_134Te)

        except RuntimeError:
            resampled_IYR_array_3n_4plus_232Th_134Te[n] = IYR_double_3n_4plus_232Th_134Te
            print("Error - curve_fit failed")

        if(n%100==0):
            print("Now finished iteration %.d" % n)

    ################   END OF BOOTSTRAP-SAMPLING   #################

    #238U lowE_highT - 134Te
    sigma_bootstrap_IYR_238U_lowE_highT_134Te = np.std(resampled_IYR_array_238U_lowE_highT_134Te)
    #238U lowE_highT - 134Te 1n
    sigma_bootstrap_IYR_1n_238U_lowE_highT_134Te = np.std(resampled_IYR_array_1n_238U_lowE_highT_134Te)
    #238U lowE_highT - 134Te 3n
    sigma_bootstrap_IYR_3n_238U_lowE_highT_134Te = np.std(resampled_IYR_array_3n_238U_lowE_highT_134Te)
    #238U lowE_highT - 134Te 3n_2plus
    sigma_bootstrap_IYR_3n_2plus_238U_lowE_highT_134Te = np.std(resampled_IYR_array_3n_2plus_238U_lowE_highT_134Te)
    #238U lowE_highT - 134Te 3n_4plus
    sigma_bootstrap_IYR_3n_4plus_238U_lowE_highT_134Te = np.std(resampled_IYR_array_3n_4plus_238U_lowE_highT_134Te)
    #238U lowE_highT - 134Te 3n_6plus
    sigma_bootstrap_IYR_3n_6plus_238U_lowE_highT_134Te = np.std(resampled_IYR_array_3n_6plus_238U_lowE_highT_134Te)
    #238U highE_highT - 134Te
    sigma_bootstrap_IYR_238U_highE_highT_134Te = np.std(resampled_IYR_array_238U_highE_highT_134Te)
    #232Th - 134Te
    sigma_bootstrap_IYR_232Th_134Te = np.std(resampled_IYR_array_232Th_134Te)
    #232Th - 134Te 1n
    sigma_bootstrap_IYR_1n_232Th_134Te = np.std(resampled_IYR_array_1n_232Th_134Te)
    #232Th - 134Te 3n
    sigma_bootstrap_IYR_3n_232Th_134Te = np.std(resampled_IYR_array_3n_232Th_134Te)
    #232Th - 134Te 4n
    sigma_bootstrap_IYR_4n_232Th_134Te = np.std(resampled_IYR_array_4n_232Th_134Te)
    #232Th - 134Te 5n
    sigma_bootstrap_IYR_5n_232Th_134Te = np.std(resampled_IYR_array_5n_232Th_134Te)
    #232Th - 134Te 3n_2plus
    sigma_bootstrap_IYR_3n_2plus_232Th_134Te = np.std(resampled_IYR_array_3n_2plus_232Th_134Te)
    #232Th - 134Te 3n_4plus
    sigma_bootstrap_IYR_3n_4plus_232Th_134Te = np.std(resampled_IYR_array_3n_4plus_232Th_134Te)    


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


if BOOTSTRAP==False:
    #238U lowE_highT - 134Te
    sigma_bootstrap_IYR_238U_lowE_highT_134Te = 0
    #238U lowE_highT - 134Te 1n
    sigma_bootstrap_IYR_1n_238U_lowE_highT_134Te = 0
    #238U lowE_highT - 134Te 3n
    sigma_bootstrap_IYR_3n_238U_lowE_highT_134Te = 0
    #238U lowE_highT - 134Te 3n_2plus
    sigma_bootstrap_IYR_3n_2plus_238U_lowE_highT_134Te = 0
    #238U lowE_highT - 134Te 3n_4plus
    sigma_bootstrap_IYR_3n_4plus_238U_lowE_highT_134Te = 0
    #238U lowE_highT - 134Te 3n_6plus
    sigma_bootstrap_IYR_3n_6plus_238U_lowE_highT_134Te = 0
    #238U highE_highT - 134Te
    sigma_bootstrap_IYR_238U_highE_highT_134Te = 0
    #232Th - 134Te
    sigma_bootstrap_IYR_232Th_134Te = 0
    #232Th - 134Te 1n
    sigma_bootstrap_IYR_1n_232Th_134Te = 0
    #232Th - 134Te 3n
    sigma_bootstrap_IYR_3n_232Th_134Te = 0
    #232Th - 134Te 4n
    sigma_bootstrap_IYR_4n_232Th_134Te = 0
    #232Th - 134Te 5n
    sigma_bootstrap_IYR_5n_232Th_134Te = 0
    #232Th - 134Te 3n_2plus
    sigma_bootstrap_IYR_3n_2plus_232Th_134Te = 0
    #232Th - 134Te 3n_4plus
    sigma_bootstrap_IYR_3n_4plus_232Th_134Te = 0  

####################################################
###     Calculate chi-squared values for fits    ###
####################################################

chisquare_double_238U_lowE_highT_134Te = 0 #reduced_chisquare_func(f_obs=y_doublegate_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_238U_lowE_highT_134Te, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), N=len(y_doublegate_238U_lowE_highT_134Te))
# chisquare_double_1n_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_1n_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_1n_238U_lowE_highT_134Te, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), N=len(y_doublegate_1n_238U_lowE_highT_134Te))
# chisquare_double_3n_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_3n_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_3n_238U_lowE_highT_134Te, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), N=len(y_doublegate_3n_238U_lowE_highT_134Te))

# chisquare_double_3n_2plus_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_3n_2plus_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_3n_2plus_238U_lowE_highT_134Te, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), N=len(y_doublegate_3n_2plus_238U_lowE_highT_134Te))
# chisquare_double_3n_4plus_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_3n_4plus_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_3n_4plus_238U_lowE_highT_134Te, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), N=len(y_doublegate_3n_4plus_238U_lowE_highT_134Te))
# chisquare_double_3n_6plus_238U_lowE_highT_134Te = reduced_chisquare_func(f_obs=y_doublegate_3n_6plus_238U_lowE_highT_134Te, f_exp=sum_smeared_exp_gauss_const_bg(x_doublegate_3n_6plus_238U_lowE_highT_134Te, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), N=len(y_doublegate_3n_6plus_238U_lowE_highT_134Te))


####################################################
###        		Print table of IYRs              ###
####################################################

print("\n")
print(" **********************************************************************************")
print("                              Isomeric Yield Ratios                                ")

t = PrettyTable(['System', 'Nucleus', 'Gate','Energies', 'IYR', 'unc_bootstrap', 'red. chi^2'])

t.add_row(['238U - lowE highT', '134Te', 'Double', '1279-297', round(IYR_double_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_238U_lowE_highT_134Te,3), round(chisquare_double_238U_lowE_highT_134Te,3)])
# t.add_row(['238U - lowE highT', '134Te', '1n', '1279-312', round(IYR_double_1n_238U_lowE_highT_134Te,3), round(sigma_IYR_1n_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_1n_238U_lowE_highT_134Te,3), round(chisquare_double_1n_238U_lowE_highT_134Te,3)])
# t.add_row(['238U - lowE highT', '134Te', '3n','1279-326', round(IYR_double_3n_238U_lowE_highT_134Te,3), round(sigma_IYR_3n_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_3n_238U_lowE_highT_134Te,3), round(chisquare_double_3n_238U_lowE_highT_134Te,3)])

# t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

# t.add_row(['238U - lowE highT', '134Te', '3n - 2plus', '1279-151', round(IYR_double_3n_2plus_238U_lowE_highT_134Te,3), round(sigma_IYR_3n_2plus_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_3n_2plus_238U_lowE_highT_134Te,3), round(chisquare_double_3n_2plus_238U_lowE_highT_134Te,3)])
# t.add_row(['238U - lowE highT', '134Te', '3n - 4plus', '1279-326', round(IYR_double_3n_4plus_238U_lowE_highT_134Te,3), round(sigma_IYR_3n_4plus_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_3n_4plus_238U_lowE_highT_134Te,3), round(chisquare_double_3n_4plus_238U_lowE_highT_134Te,3)])
# t.add_row(['238U - lowE highT', '134Te', '3n - 6plus', '1279-487', round(IYR_double_3n_6plus_238U_lowE_highT_134Te,3), round(sigma_IYR_3n_6plus_238U_lowE_highT_134Te,3), round(sigma_bootstrap_IYR_3n_6plus_238U_lowE_highT_134Te,3), round(chisquare_double_3n_6plus_238U_lowE_highT_134Te,3)])

# t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

# #t.add_row(['238U - lowE lowT', '134Te', 'Double', '1279-297', round(IYR_double_238U_lowE_lowT_134Te,3), round(sigma_IYR_238U_lowE_lowT_134Te,3), 0, 0])

# t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

# t.add_row(['238U - highE highT', '134Te', 'Double', '1279-297', round(IYR_double_238U_highE_highT_134Te,3), round(sigma_IYR_238U_highE_highT_134Te,3), round(sigma_bootstrap_IYR_238U_lowE_highT_134Te,3), 0])

# t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

# t.add_row(['232Th - lowE highT', '134Te', 'Double', '1279-297', round(IYR_double_232Th_134Te,3), round(sigma_IYR_232Th_134Te,3), round(sigma_bootstrap_IYR_232Th_134Te,3), 0])
# t.add_row(['232Th - lowE highT', '134Te', '1n', '1279-433', round(IYR_double_1n_232Th_134Te,3), round(sigma_IYR_1n_232Th_134Te,3), round(sigma_bootstrap_IYR_1n_232Th_134Te,3), 0])
# t.add_row(['232Th - lowE highT', '134Te', '3n', '1279-815', round(IYR_double_3n_232Th_134Te,3), round(sigma_IYR_3n_232Th_134Te,3), round(sigma_bootstrap_IYR_3n_232Th_134Te,3), 0])
# t.add_row(['232Th - lowE highT', '134Te', '4n', '1279-352', round(IYR_double_4n_232Th_134Te,3), round(sigma_IYR_4n_232Th_134Te,3), round(sigma_bootstrap_IYR_4n_232Th_134Te,3), 0])
# t.add_row(['232Th - lowE highT', '134Te', '5n', '1279-837', round(IYR_double_5n_232Th_134Te,3), round(sigma_IYR_5n_232Th_134Te,3), round(sigma_bootstrap_IYR_5n_232Th_134Te,3), 0])

# t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

# t.add_row(['232Th - lowE highT', '134Te', '3n - 2plus', '1279-815', round(IYR_double_3n_2plus_232Th_134Te,3), round(sigma_IYR_3n_2plus_232Th_134Te,3), round(sigma_bootstrap_IYR_3n_2plus_232Th_134Te,3), 0])
# t.add_row(['232Th - lowE highT', '134Te', '3n - 4plus', '1279-978', round(IYR_double_3n_4plus_232Th_134Te,3), round(sigma_IYR_3n_4plus_232Th_134Te,3), round(sigma_bootstrap_IYR_3n_4plus_232Th_134Te,3), 0])


print(t)



####################################################
## 						Plot 		              ## START_PLOT
####################################################

x_array_plot = np.linspace(0,1000,10000)


################   238U lowE_highT -  134Te   #################

# #plt.plot(x_doublegate_238U_lowE_highT_134Te_long, y_doublegate_238U_lowE_highT_134Te_long, label="doublegate_238U_lowE_highT_134Te", color="royalblue")
# plt.errorbar(x_doublegate_238U_lowE_highT_134Te_long, y_doublegate_238U_lowE_highT_134Te_long, yerr=sigma_data_doublegate_all_bg(data_all=y_doublegate_all_238U_lowE_highT_134Te_long, data_bg_ridge=y_doublegate_bg_ridge_238U_lowE_highT_134Te_long, data_bg_random=y_doublegate_bg_random_238U_lowE_highT_134Te_long), label="doublegate_238U_lowE_highT_134Te", color="royalblue")


# plt.plot(x_array_plot, sum_smeared_exp_two_gauss_const_bg(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss_1(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="gauss1", color="lightgreen")
# plt.plot(x_array_plot, gauss_2(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="gauss1", color="yellowgreen")
# #plt.plot(x_array_plot, gauss_1(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7])+gauss_2(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="two gauss", color="green")

# plt.plot(x_array_plot, smeared_exp_decay_(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="exp decay", color="red")
# plt.plot(x_array_plot, const_bg_(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5], P_double_238U_lowE_highT_134Te[6], P_double_238U_lowE_highT_134Te[7]), label="const bg", color="pink")


# #plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# #plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# #plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# #plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_highT_134Te[0], P_double_238U_lowE_highT_134Te[1], P_double_238U_lowE_highT_134Te[2], P_double_238U_lowE_highT_134Te[3], P_double_238U_lowE_highT_134Te[4], P_double_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te: Doublegate true spectrum fit")
# #plt.axis([0,700,1,10**(4)])
# plt.axis([0,700,-50,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=12)
# plt.grid()
# plt.show()

################   238U lowE_highT -  134Te  1n  #################

# plt.plot(x_doublegate_1n_238U_lowE_highT_134Te_long, y_doublegate_1n_238U_lowE_highT_134Te_long, label="doublegate_1n_238U_lowE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_1n_238U_lowE_highT_134Te[0], P_double_1n_238U_lowE_highT_134Te[1], P_double_1n_238U_lowE_highT_134Te[2], P_double_1n_238U_lowE_highT_134Te[3], P_double_1n_238U_lowE_highT_134Te[4], P_double_1n_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_1n_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_1n_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te 1n: Doublegate true spectrum fit")
# plt.axis([0,700,0,2000])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   238U lowE_highT -  134Te  3n  #################

# plt.plot(x_doublegate_3n_238U_lowE_highT_134Te_long, y_doublegate_3n_238U_lowE_highT_134Te_long, label="doublegate_3n_238U_lowE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_238U_lowE_highT_134Te[0], P_double_3n_238U_lowE_highT_134Te[1], P_double_3n_238U_lowE_highT_134Te[2], P_double_3n_238U_lowE_highT_134Te[3], P_double_3n_238U_lowE_highT_134Te[4], P_double_3n_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te 3n: Doublegate true spectrum fit")
# plt.axis([0,700,0,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   238U lowE_highT -  134Te  3n_2plus  #################

# plt.plot(x_doublegate_3n_2plus_238U_lowE_highT_134Te_long, y_doublegate_3n_2plus_238U_lowE_highT_134Te_long, label="doublegate_3n_2plus_238U_lowE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_2plus_238U_lowE_highT_134Te[0], P_double_3n_2plus_238U_lowE_highT_134Te[1], P_double_3n_2plus_238U_lowE_highT_134Te[2], P_double_3n_2plus_238U_lowE_highT_134Te[3], P_double_3n_2plus_238U_lowE_highT_134Te[4], P_double_3n_2plus_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_2plus_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_2plus_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te 3n_2plus: Doublegate true spectrum fit")
# plt.axis([0,700,0,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   238U lowE_highT -  134Te  3n_4plus  #################

# plt.plot(x_doublegate_3n_4plus_238U_lowE_highT_134Te_long, y_doublegate_3n_4plus_238U_lowE_highT_134Te_long, label="doublegate_3n_4plus_238U_lowE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_4plus_238U_lowE_highT_134Te[0], P_double_3n_4plus_238U_lowE_highT_134Te[1], P_double_3n_4plus_238U_lowE_highT_134Te[2], P_double_3n_4plus_238U_lowE_highT_134Te[3], P_double_3n_4plus_238U_lowE_highT_134Te[4], P_double_3n_4plus_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_4plus_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_4plus_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te 3n_4plus: Doublegate true spectrum fit")
# plt.axis([0,700,-50,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   238U lowE_highT -  134Te  3n_6plus  #################

# plt.plot(x_doublegate_3n_6plus_238U_lowE_highT_134Te_long, y_doublegate_3n_6plus_238U_lowE_highT_134Te_long, label="doublegate_3n_6plus_238U_lowE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_6plus_238U_lowE_highT_134Te[0], P_double_3n_6plus_238U_lowE_highT_134Te[1], P_double_3n_6plus_238U_lowE_highT_134Te[2], P_double_3n_6plus_238U_lowE_highT_134Te[3], P_double_3n_6plus_238U_lowE_highT_134Te[4], P_double_3n_6plus_238U_lowE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_6plus_238U_lowE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_6plus_238U_lowE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_highT - 134Te 3n_6plus: Doublegate true spectrum fit")
# plt.axis([0,700,0,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   238U lowE_lowT -  134Te   #################

# plt.plot(x_doublegate_238U_lowE_lowT_134Te_long, y_doublegate_238U_lowE_lowT_134Te_long, label="doublegate_238U_lowE_lowT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_lowT_134Te[0], P_double_238U_lowE_lowT_134Te[1], P_double_238U_lowE_lowT_134Te[2], P_double_238U_lowE_lowT_134Te[3], P_double_238U_lowE_lowT_134Te[4], P_double_238U_lowE_lowT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_lowE_lowT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_lowE_lowT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U lowE_lowT - 134Te: Doublegate true spectrum fit")
# #plt.axis([0,700,1,10**(4)])
# plt.axis([200,700,0,4*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


################   238U highE_highT -  134Te   #################

# plt.plot(x_doublegate_238U_highE_highT_134Te_long, y_doublegate_238U_highE_highT_134Te_long, label="doublegate_238U_highE_highT_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_highE_highT_134Te[0], P_double_238U_highE_highT_134Te[1], P_double_238U_highE_highT_134Te[2], P_double_238U_highE_highT_134Te[3], P_double_238U_highE_highT_134Te[4], P_double_238U_highE_highT_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_highE_highT_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_highE_highT_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("238U highE_highT - 134Te: Doublegate true spectrum fit")
# #plt.axis([0,700,1,10**(4)])
# plt.axis([0,600,0,4*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=14)
# plt.grid()
# plt.show()


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
# plt.axis([0,700,1,2*10**(4)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   232Th -  134Te  1n  #################

# plt.plot(x_doublegate_1n_232Th_134Te_long, y_doublegate_1n_232Th_134Te_long, label="doublegate_1n_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_1n_232Th_134Te[0], P_double_1n_232Th_134Te[1], P_double_1n_232Th_134Te[2], P_double_1n_232Th_134Te[3], P_double_1n_232Th_134Te[4], P_double_1n_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_1n_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_1n_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th - 134Te 1n: Doublegate true spectrum fit")
# plt.axis([0,700,0,2000])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   232Th -  134Te  3n  #################

# plt.plot(x_doublegate_3n_232Th_134Te_long, y_doublegate_3n_232Th_134Te_long, label="doublegate_3n_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_232Th_134Te[0], P_double_3n_232Th_134Te[1], P_double_3n_232Th_134Te[2], P_double_3n_232Th_134Te[3], P_double_3n_232Th_134Te[4], P_double_3n_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th - 134Te 3n: Doublegate true spectrum fit")
# plt.axis([0,700,0,10000])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   232Th -  134Te  4n  #################

# plt.plot(x_doublegate_4n_232Th_134Te_long, y_doublegate_4n_232Th_134Te_long, label="doublegate_4n_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_4n_232Th_134Te[0], P_double_4n_232Th_134Te[1], P_double_4n_232Th_134Te[2], P_double_4n_232Th_134Te[3], P_double_4n_232Th_134Te[4], P_double_4n_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_4n_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_4n_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th - 134Te 4n: Doublegate true spectrum fit")
# plt.axis([0,700,0,10000])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   232Th -  134Te  5n  #################

# plt.plot(x_doublegate_5n_232Th_134Te_long, y_doublegate_5n_232Th_134Te_long, label="doublegate_5n_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_5n_232Th_134Te[0], P_double_5n_232Th_134Te[1], P_double_5n_232Th_134Te[2], P_double_5n_232Th_134Te[3], P_double_5n_232Th_134Te[4], P_double_5n_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_5n_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_5n_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th - 134Te 5n: Doublegate true spectrum fit")
# plt.axis([0,700,0,10000])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()

################   232Th -  134Te  3n_2plus  #################

# plt.plot(x_doublegate_3n_2plus_232Th_134Te_long, y_doublegate_3n_2plus_232Th_134Te_long, label="doublegate_3n_2plus_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_2plus_232Th_134Te[0], P_double_3n_2plus_232Th_134Te[1], P_double_3n_2plus_232Th_134Te[2], P_double_3n_2plus_232Th_134Te[3], P_double_3n_2plus_232Th_134Te[4], P_double_3n_2plus_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_2plus_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_2plus_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th lowE_highT - 134Te 3n_2plus: Doublegate true spectrum fit")
# plt.axis([0,700,0,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   232Th -  134Te  3n_4plus  #################

# plt.plot(x_doublegate_3n_4plus_232Th_134Te_long, y_doublegate_3n_4plus_232Th_134Te_long, label="doublegate_3n_4plus_232Th_134Te", color="royalblue")

# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_3n_4plus_232Th_134Te[0], P_double_3n_4plus_232Th_134Te[1], P_double_3n_4plus_232Th_134Te[2], P_double_3n_4plus_232Th_134Te[3], P_double_3n_4plus_232Th_134Te[4], P_double_3n_4plus_232Th_134Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_3n_4plus_232Th_134Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_3n_4plus_232Th_134Te[-1],0,6000, color="black")
# #plt.yscale("log")
# plt.title("232Th lowE_highT - 134Te 3n_4plus: Doublegate true spectrum fit")
# plt.axis([0,700,0,5*10**(3)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()


################   238U -  135Te   #################

# plt.plot(x_doublegate_238U_lowE_highT_135Te_long, y_doublegate_238U_lowE_highT_135Te_long, label="doublegate_238U_lowE_highT_135Te", color="royalblue")
# plt.plot(x_array_plot, sum_smeared_exp_gauss_const_bg(x_array_plot, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), label="true fit, total", color="orange")
# plt.plot(x_array_plot, gauss(x_array_plot, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), label="true gaussian", color="green")
# plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), label="true smeared exp decay", color="red")
# plt.plot(x_array_plot, const_bg(x_array_plot, P_double_238U_lowE_highT_135Te[0], P_double_238U_lowE_highT_135Te[1], P_double_238U_lowE_highT_135Te[2], P_double_238U_lowE_highT_135Te[3], P_double_238U_lowE_highT_135Te[4], P_double_238U_lowE_highT_135Te[5]), label="constant BG", color="hotpink")

# plt.vlines(x_doublegate_238U_lowE_highT_135Te[0],0,6000, label="fit range", color="black")
# plt.vlines(x_doublegate_238U_lowE_highT_135Te[-1],0,6000, color="black")
# plt.yscale("log")
# plt.title("238U - 135Te: Doublegate true spectrum fit")
# plt.axis([0,700,1,10**(4)])
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Counts", fontsize=14)
# plt.legend(fontsize=10)
# plt.grid()
# plt.show()








