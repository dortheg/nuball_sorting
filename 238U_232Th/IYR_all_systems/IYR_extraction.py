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

from IYR_extraction_class import IYR_extraction_class

####################################################
##     				Files                  	      ## START_DATA
####################################################

file_238U_lowE= ROOT.TFile.Open("../Sorted_files/238Ucube_hit4_lowE_highT_8feb2022.bin_final.root"," READ ")
file_238U_highE = ROOT.TFile.Open("../Sorted_files/238Ucube_hit4_2ns_highE_highT_19jan2022.bin.root"," READ ")
file_232Th = ROOT.TFile.Open("../Sorted_files/232Thcube_hit4_lowE_highT_28feb2022.bin_final.root"," READ ")

#Paramterer choices for analysis
x_lower = 320
N_BOOTSTRAP = 1000


####################################################
##     				Analyze                       ## 
####################################################


double_238U_lowE_134Te = IYR_extraction_class(name="doublegate", gatetype="double", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=False, plot=False, calc_redchi2=True)
#double_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)

double_1n_238U_lowE_134Te = IYR_extraction_class(name="doublegate_1n", gatetype="partner", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_1n_238U_lowE_134Te.plot_spec()
#double_1n_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=False, plot=False)
#double_1n_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)

double_3n_238U_lowE_134Te = IYR_extraction_class(name="doublegate_3n", gatetype="partner", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_3n_238U_lowE_134Te.plot_spec()
#double_3n_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=False, plot=False)
#double_3n_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)

double_3n_2plus_238U_lowE_134Te = IYR_extraction_class(name="doublegate_3n_2plus", gatetype="partner", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_3n_2plus_238U_lowE_134Te.plot_spec()
#double_3n_2plus_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=False, plot=False)
#double_3n_2plus_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)

double_3n_4plus_238U_lowE_134Te = IYR_extraction_class(name="doublegate_3n_4plus", gatetype="partner", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_3n_4plus_238U_lowE_134Te.plot_spec()
#double_3n_4plus_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=False, plot=False)
#double_3n_4plus_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)

double_3n_6plus_238U_lowE_134Te = IYR_extraction_class(name="doublegate_3n_6plus", gatetype="partner", CN="238U", file=file_238U_lowE, x_lower=x_lower)
#double_3n_6plus_238U_lowE_134Te.plot_spec()
double_3n_6plus_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=True, plot=True)
#double_3n_6plus_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=N_BOOTSTRAP, plot=False)


####################################################
###        		Print table of IYRs              ###
####################################################

print("\n")
print(" **********************************************************************************")
print("                              Isomeric Yield Ratios                                ")

t = PrettyTable(['System', 'Nucleus', 'Gate','Energies', 'IYR', 'unc_numpy', 'unc_handcalc', 'red_chi^2'])

t.add_row(['238U - lowE', '134Te', 'Double', '1279-297', round(double_238U_lowE_134Te.IYR,3), round(double_238U_lowE_134Te.sigma_IYR,3), round(double_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_238U_lowE_134Te.red_chi2,3)])
t.add_row(['238U - lowE', '134Te', '1n', '1279-312', round(double_1n_238U_lowE_134Te.IYR,3), round(double_1n_238U_lowE_134Te.sigma_IYR,3), round(double_1n_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_1n_238U_lowE_134Te.red_chi2,3)])
t.add_row(['238U - lowE', '134Te', '3n', '1279-326', round(double_3n_238U_lowE_134Te.IYR,3), round(double_3n_238U_lowE_134Te.sigma_IYR,3), round(double_3n_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_3n_238U_lowE_134Te.red_chi2,3)])

t.add_row([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])

t.add_row(['238U - lowE', '134Te', '3n_2plus', '1279-151', round(double_3n_2plus_238U_lowE_134Te.IYR,3), round(double_3n_2plus_238U_lowE_134Te.sigma_IYR,3), round(double_3n_2plus_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_3n_2plus_238U_lowE_134Te.red_chi2,3)])
t.add_row(['238U - lowE', '134Te', '3n_4plus', '1279-326', round(double_3n_4plus_238U_lowE_134Te.IYR,3), round(double_3n_4plus_238U_lowE_134Te.sigma_IYR,3), round(double_3n_4plus_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_3n_4plus_238U_lowE_134Te.red_chi2,3)])
t.add_row(['238U - lowE', '134Te', '3n_6plus', '1279-487', round(double_3n_6plus_238U_lowE_134Te.IYR,3), round(double_3n_6plus_238U_lowE_134Te.sigma_IYR,3), round(double_3n_6plus_238U_lowE_134Te.sigma_IYR_handcalc,3), round(double_3n_6plus_238U_lowE_134Te.red_chi2,3)])


print(t)


