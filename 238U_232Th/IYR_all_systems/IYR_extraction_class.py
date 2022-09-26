import ROOT
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
from scipy.stats import chisquare
from scipy.ndimage import gaussian_filter1d

class IYR_extraction_class:

    def __init__(self, name, CN, file, x_lower):
        self.name = name
        self.CN = CN   
        self.file = file
        self.x_lower = x_lower
        self.x_upper = 640

        #Lifetimes we need
        self.tau_134Te = 164.1/np.log(2)
        self.sigma_tau_134Te = 0.9/np.log(2)

        def load_spec(name_spec):


            bin_lower = self.x_lower//2
            bin_upper = self.x_upper//2

            spec_name = "time_isomer_" + name_spec + "_134Te"

            histogram = self.file.Get(spec_name)
            x_bins = histogram.GetNbinsX()

            x_spec = np.zeros(x_bins)
            y_spec = np.zeros(x_bins)
    
            for i in range(x_bins):
                x_spec[i] = histogram.GetBinCenter(i+1)
                y_spec[i] = histogram.GetBinContent(i+1)

            #Make up for 2ns bins
            x_spec = 2*x_spec

            x_spec_long = x_spec
            y_spec_long = y_spec

            x_spec = x_spec[bin_lower:bin_upper]
            y_spec = y_spec[bin_lower:bin_upper]

            return x_spec, y_spec, x_spec_long, y_spec_long

        self.x_spec, self.y_spec, self.x_spec_long, self.y_spec_long = load_spec(name_spec=self.name)
        self.x_spec_all, self.y_spec_all, self.x_spec_all_long, self.y_spec_all_long = load_spec(name_spec=self.name+"_all")
        self.x_spec_bg, self.y_spec_bg, self.x_spec_bg_long, self.y_spec_bg_long = load_spec(name_spec=self.name+"_bg")
        self.x_spec_bg_ridge, self.y_spec_bg_ridge, self.x_spec_bg_ridge_long, self.y_spec_bg_ridge_long = load_spec(name_spec=self.name+"_bg_ridge")
        self.x_spec_bg_random, self.y_spec_bg_random, self.x_spec_bg_random_long, self.y_spec_bg_random_long = load_spec(name_spec=self.name+"_bg_random")
        
        def yerr(spec_all, spec_bg_ridge, spec_bg_random):
            unc = np.zeros(len(spec_all))

            for i in range(len(spec_all)):
                if spec_all[i] <= 0 and spec_bg_ridge[i]<= 0 and spec_bg_random[i]<= 0:
                    unc[i] = 1 #Should it be 1?
                else:
                    unc[i] = np.sqrt(spec_all[i]+spec_bg_ridge[i]+spec_bg_random[i])
            return unc

        self.y_spec_unc = yerr(self.x_spec_all, self.x_spec_bg_ridge, self.x_spec_bg_random)
        self.y_spec_long_unc = yerr(self.x_spec_all_long, self.x_spec_bg_ridge_long, self.x_spec_bg_random_long)


###########################################################################################


    def plot_spec(self):

        plt.errorbar(self.x_spec_long, self.y_spec_long, yerr=self.y_spec_long_unc, label="true")
        #plt.plot(self.x_spec_all_long, self.y_spec_all_long, label="all")
        #plt.plot(self.x_spec_bg_long, self.y_spec_bg_long, label="bg")
        #plt.plot(self.x_spec_bg_ridge_long, self.y_spec_bg_ridge_long - self.y_spec_bg_random_long, label="bg_ridge - bg_random")
        plt.legend()
        plt.xlabel("Time [ns]")
        plt.xlabel("Counts")
        plt.grid()
        plt.show()


###########################################################################################


    def fit_spec_and_IYR_calc(self, plot=False, show_fitparam=False):

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

        def smeared_exp_decay(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            """ Exponential decay """
            return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma1)

        def const_bg(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            return const_bg + x*0

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
        tau_decay_lower = self.tau_134Te
        tau_decay_upper = self.tau_134Te+0.0001

        self.P, self.cov = curve_fit(sum_smeared_exp_two_gauss_const_bg, self.x_spec, self.y_spec, sigma=self.y_spec_unc, bounds=([mean_lower, sigma1_lower, amplitude_gauss1_lower, sigma2_lower, amplitude_gauss2_lower, const_bg_lower, amplitude_exp_decay_lower, tau_decay_lower],[mean_upper, sigma1_upper, amplitude_gauss1_upper, sigma2_upper, amplitude_gauss2_upper, const_bg_upper, amplitude_exp_decay_upper, tau_decay_upper]), absolute_sigma = False)
        self.P_unc = np.sqrt(np.diag(self.cov))

        if show_fitparam==True:
            print("\n")
            print(" *****       %s    %s    134Te ***** " % (self.CN, self.name) )
            print("          -- 2 GAUSS + SMEARED EXP + CONST_BG FIT --   ")
            print("mean:                     %.2f +/- %.2f          [%.d,%.d]" % (self.P[0], self.P_unc[0], mean_lower, mean_upper))
            print("sigma1:                    %.2f +/- %.2f         [%.d,%.d]" % (self.P[1], self.P_unc[1], sigma1_lower, sigma1_upper))
            print("amplitude_gauss1:          %.2f +/- %.2f         [%.d,%.d]" % (self.P[2], self.P_unc[2], amplitude_gauss1_lower, amplitude_gauss1_upper))
            print("sigma2:                    %.2f +/- %.2f         [%.d,%.d]" % (self.P[3], self.P_unc[3], sigma2_lower, sigma2_upper))
            print("amplitude_gauss2:          %.2f +/- %.2f         [%.d,%.d]" % (self.P[4], self.P_unc[4], amplitude_gauss2_lower, amplitude_gauss2_upper))
            print("const_bg:                 %.2f +/- %.2f         [%.d,%.d]" % (self.P[5], self.P_unc[5], const_bg_lower, const_bg_upper))
            print("amplitude_exp_decay:      %.2f +/- %.2f         [%.d,%.d]" % (self.P[6], self.P_unc[6], amplitude_exp_decay_lower, amplitude_exp_decay_upper))
            print("tau_decay, in half_life:  %.2f +/- %.2f         [%.d,%.d]" % (self.P[7]*np.log(2), 0, tau_decay_lower*np.log(2), tau_decay_upper*np.log(2)))
            print("\n")

        if plot == True:

            x_array_plot = np.linspace(0,1000,10000)

            plt.errorbar(self.x_spec_long, self.y_spec_long, yerr=self.y_spec_long_unc, label="true")
            plt.plot(x_array_plot, sum_smeared_exp_two_gauss_const_bg(x_array_plot, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), label="true fit, total", color="orange")
            plt.plot(x_array_plot, gauss_1(x_array_plot, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), label="gauss1", color="lightgreen")
            plt.plot(x_array_plot, gauss_2(x_array_plot, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), label="gauss1", color="yellowgreen")

            plt.plot(x_array_plot, smeared_exp_decay(x_array_plot, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), label="exp decay", color="red")
            plt.plot(x_array_plot, const_bg(x_array_plot, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), label="const bg", color="pink")

            plt.vlines(self.x_spec[0],0,6000, label="fit range", color="black")
            plt.vlines(self.x_spec[-1],0,6000, color="black")
            
            plt.title("%s %s  134Te" % (self.CN, self.name) )  
            plt.axis([0,700,-50,5*10**(3)])
            plt.xlabel("Time [ns]", fontsize=14)
            plt.ylabel("Counts", fontsize=14)
            plt.legend(fontsize=12)
            plt.grid()
            plt.show()


        ####################################################
        ##                   Find IYR                     ## 
        ####################################################

        #NB: Upper integration limit should be high enough that no changes are observed in the IYR when increasing the range
        self.x_arr_134Te = np.linspace(0,100*round(self.tau_134Te),10000)

        area_all = np.trapz(gauss_1(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7])+gauss_2(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]) + smeared_exp_decay(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), self.x_arr_134Te)
        area_prompt = np.trapz(gauss_1(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7])+gauss_2(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), self.x_arr_134Te)
        area_delayed = np.trapz(smeared_exp_decay(self.x_arr_134Te, self.P[0], self.P[1], self.P[2], self.P[3], self.P[4], self.P[5], self.P[6], self.P[7]), self.x_arr_134Te)

        self.IYR = area_delayed/(area_prompt+area_delayed)

###########################################################################################


    def IYR_uncertainty(self, N_BOOTSTRAP=1000, plot=False):

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

        def smeared_exp_decay(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            """ Exponential decay """
            return gaussian_filter1d(np.piecewise(x, [x < mean, x >= mean], [lambda x:0, lambda x:amplitude_exp_decay*np.exp((mean-x)/tau_decay)]),sigma1)

        def const_bg(x, mean=0, sigma1=1.0, amplitude_gauss1=1.0, sigma2=1.0, amplitude_gauss2=1.0, const_bg=1.0, amplitude_exp_decay=1.0, tau_decay=1.0):
            return const_bg + x*0

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
        tau_decay_lower = self.tau_134Te
        tau_decay_upper = self.tau_134Te+0.0001

        print("\n Starting BOOTSTRAPPING... Number of iterations: %.d" % N_BOOTSTRAP)

        resampled_IYR_array = np.zeros(N_BOOTSTRAP)
        resampled_y_spec_long = np.zeros(len(self.y_spec_long))

        for n in range(N_BOOTSTRAP):
            #Vary BG withing a normal distribution of sigma = 0.025, then +-2sigma spans a BG-variation of 5%
            y_spec_long_bgvaried = self.y_spec_all_long-self.y_spec_bg_long*np.random.normal(1,0.025)    

            for i in range (len(self.y_spec_all_long)):
                #Vary value of each bin within uncertainty
                resampled_y_spec_long[i] = y_spec_long_bgvaried[i] + np.random.normal(0, self.y_spec_long_unc[i]) 

            x_lower_new = self.x_lower + int(np.random.normal(0,5))
            x_upper_new = self.x_upper + int(np.random.normal(0,5))

            bin_lower = x_lower_new//2
            bin_upper = x_upper_new//2

            try:        
                resampled_P, resampled_cov = curve_fit(sum_smeared_exp_two_gauss_const_bg, self.x_spec_long[bin_lower:bin_upper], resampled_y_spec_long[bin_lower:bin_upper], sigma=self.y_spec_long_unc[bin_lower:bin_upper], bounds=([mean_lower, sigma1_lower, amplitude_gauss1_lower, sigma2_lower, amplitude_gauss2_lower, const_bg_lower, amplitude_exp_decay_lower, tau_decay_lower],[mean_upper, sigma1_upper, amplitude_gauss1_upper, sigma2_upper, amplitude_gauss2_upper, const_bg_upper, amplitude_exp_decay_upper, tau_decay_upper]), absolute_sigma = False, maxfev=10000)

                resampled_area_all = np.trapz(gauss_1(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7])+gauss_2(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7]) + smeared_exp_decay(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7]), self.x_arr_134Te)
                resampled_area_prompt = np.trapz(gauss_1(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7])+gauss_2(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7]), self.x_arr_134Te)
                resampled_area_delayed = np.trapz(smeared_exp_decay(self.x_arr_134Te, resampled_P[0], resampled_P[1], resampled_P[2], resampled_P[3], resampled_P[4], resampled_P[5], resampled_P[6], resampled_P[7]), self.x_arr_134Te)

                resampled_IYR_array[n] = resampled_area_delayed /(resampled_area_delayed + resampled_area_prompt)

            except RuntimeError:
                resampled_IYR_array[n] = self.IYR 


            if(n%100==0):
                print("Now finished iteration %.d" % n)

        self.sigma_IYR = np.std(resampled_IYR_array)

        # Plot distribution of calculated IYRs
        IYR_binwith = 0.01
        IYR_min = 0.5
        IYR_max = 1.0 + 2*IYR_binwith
        N_IYR_bins = int((IYR_max-IYR_min)//IYR_binwith)

        IYR_array = np.linspace(IYR_min,IYR_max,N_IYR_bins)
        IYR_histogram = np.zeros(len(IYR_array))

        for i in range(N_BOOTSTRAP):
            IYR_bin = int(resampled_IYR_array[i]//IYR_binwith - IYR_min//IYR_binwith)
            IYR_histogram[IYR_bin] += 1

        if plot == True:    
            plt.plot(IYR_array,IYR_histogram, label="Calculated IYR distribution")
            plt.axvline(x=self.IYR, ymin=0, ymax=max(IYR_histogram), label="IYR", color="red")
            plt.axvline(x=self.IYR-self.sigma_IYR, ymin=0, ymax=max(IYR_histogram), label="IYR-sigma", color="black")
            plt.axvline(x=self.IYR+self.sigma_IYR, ymin=0, ymax=max(IYR_histogram), label="IYR+sigma", color="black")
            plt.xlabel("Calculated IYR", fontsize=14)
            plt.ylabel("Counts", fontsize=14)
            plt.legend(fontsize=14)
            plt.grid()
            plt.show()

    #     # Calculate standard deviation by hand
    #     std_value = np.linspace(0,0.5,10000)
    #     std_found = False
    #     std_handcalc = 0

    # for i in range(len(std_value)):
    #     within_std = np.where((resampled_IYR_array >= IYR_double_238U_lowE_highT_134Te-std_value[i]) & (resampled_IYR_array_238U_lowE_highT_134Te < IYR_double_238U_lowE_highT_134Te+std_value[i]))[0]
    #     percentage_covered = len(within_std)/len(resampled_IYR_array_238U_lowE_highT_134Te)

    #     if 0.675 <= percentage_covered < 0.685:
    #         std_found = True
    #         std_handcalc = std_value[i]
    #         print("By-hand standard deviation is %.3f" % std_value[i])
    #         break



