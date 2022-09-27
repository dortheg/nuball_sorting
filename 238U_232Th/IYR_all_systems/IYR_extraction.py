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



double_238U_lowE_134Te = IYR_extraction_class(name="doublegate", CN="238U", file=file_238U_lowE, x_lower=330)
double_238U_lowE_134Te.fit_spec_and_IYR_calc(print_fitparam=True, plot=False)
#double_238U_lowE_134Te.IYR_uncertainty(N_BOOTSTRAP=1000, plot=True)

# double_3n_238U_lowE_134Te = IYR_extraction_class(name="doublegate_3n", CN="238U", file=file_238U_lowE)
# double_3n_238U_lowE_134Te.load_spec_all()
# double_3n_238U_lowE_134Te.plot_spec()


