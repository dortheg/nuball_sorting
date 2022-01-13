import numpy as np 
import matplotlib.pyplot as plt 


def gatewidth(energy):
	#Jon's parameterisation of the energy-dependence of the FWHM
	NFWHM=1.0
	A=1.059
	B=2.814
	gatewidth = (0.5+(np.sqrt(A*A+(B*B*energy/1000.0))/1.0))*(NFWHM/2.0)
	return gatewidth 

def round_gatewidth(energy):
	#Jon's parameterisation of the energy-dependence of the FWHM
	NFWHM=1.0
	A=1.059
	B=2.814
	gatewidth = (0.5+(np.sqrt(A*A+(B*B*energy/1000.0))/1.0))*(NFWHM/2.0)
	gatewidth_round = np.round(gatewidth)
	return gatewidth_round 

energy_array = np.linspace(0,2000,5000)

plt.plot(energy_array, gatewidth(energy_array), label="function value")
plt.plot(energy_array, round_gatewidth(energy_array), label="gatewidth value")
plt.xlabel("Energy [keV]", fontsize=12)
plt.xlabel("Gatewidth, gate is Energy +/- Gatewidth [keV]", fontsize=12)
plt.legend(fontsize=12)
plt.grid()
plt.plot()
plt.show()