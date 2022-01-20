#include <sstream>
#include "TImage.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TROOT.h"
#include "TText.h"
#include "TColor.h"
#include "TAttImage.h"
#include "TEnv.h"
#include "TH2D.h"
#include "TF2.h"
#include "TH3I.h"
#include "TColor.h"
#include "TLine.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include <fstream>
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "CubeDDT.hxx"
#include "CubeDDT.cxx"

using namespace std;

#include "WriteRWSpec.cxx"

int TBINS=350;
int EBINS=2048;

short lookup_134Te[2048][2048] = {{0}};
short lookup_1n_134Te[2048][2048] = {{0}};
short lookup_2n_134Te[2048][2048] = {{0}};
short lookup_3n_134Te[2048][2048] = {{0}};
short lookup_4n_134Te[2048][2048] = {{0}};
short lookup_5n_134Te[2048][2048] = {{0}};
short lookup_6n_134Te[2048][2048] = {{0}};
short lookup_135Te[2048][2048] = {{0}};
short lookup_1_93Rb[2048][2048] = {{0}};
short lookup_2_93Rb[2048][2048] = {{0}};
short lookup_94Rb[2048][2048] = {{0}};
short lookup_140Xe[2048][2048] = {{0}};
short lookup_138Xe[2048][2048] = {{0}};
short lookup_92Sr[2048][2048] = {{0}};
short lookup_94Sr[2048][2048] = {{0}};


void fill_lookuptable(int energy_1, int energy_2, short lookup[2048][2048]);
void fill_spectra(int lookup_value, int cube_value, int k, TH1D *true_spec, TH1D *all_spec, TH1D *bg_spec, TH1D *bg_ridge_spec, TH1D *bg_random_spec);


int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_CubeSort.cxx"
	CubeDDT *Cube1=new CubeDDT("",EBINS,TBINS,2);

	int A = 252;

	if(A==238){
		//238U
		//Cube1->Read("Data_cubes/238Ucube_hit3_2ns_lowE_4jan2022.bin");
		Cube1->Read("Data_cubes/238Ucube_hit4_2ns_lowE_12jan2022.bin");
		//Cube1->Read("Data_cubes/238Ucube_hit5_2ns_lowE_12jan2022.bin");
	}

	else if(A==232){
		//232Th
		Cube1->Read("Data_cubes/232Thcube_hit4_2ns_17jan2022.bin");
	}

	else if(A==252){
		//232Th
		Cube1->Read("Data_cubes/252Cfcube_hit4_2ns_19jan2022.bin");
	}

 	//////////////////////////////////////////
	/// 		  Isomer properties	       ///
	//////////////////////////////////////////


	//134Te
 	int gamma_energy_1_134Te = 297; //279.0
	int gamma_energy_2_134Te = 1279; //1279.01

	int gamma_energy_1n = 0;
	int gamma_energy_2n = 0;
	int gamma_energy_3n = 0;
	int gamma_energy_4n = 0;
	int gamma_energy_5n = 0;
	int gamma_energy_6n = 0;

	if(A==238){
		//104Zr, 1n-partner with 134Te for 238U
		gamma_energy_1n = 312; //312.2

		//102Zr, 3n-partner with 134Te
		gamma_energy_3n = 326; //326.48

		//101Zr, 4n-partner with 134Te
		gamma_energy_4n = 217; //216.68

		//100Zr, 5n-partner with 134Te
		gamma_energy_5n = 212; //212.61
	}

	if(A==232){
		//98Sr
		gamma_energy_1n = 433; //144.70

		//96Sr
		gamma_energy_3n = 815; //815.0

		//95Sr
		gamma_energy_4n = 352; //352.01

		//94Sr
		gamma_energy_5n = 837; //836.9
	}

	if(A==252){
		//117Pd
		gamma_energy_1n = 438; //438.0, cannot go lower due to isomer

		//116Pd
		gamma_energy_2n = 340; //340.3

		//115Pd
		gamma_energy_3n = 128; //127.8

		//114Pd
		gamma_energy_4n = 333; //332.6

		//113Pd
		gamma_energy_5n = 128; //127.8

		//122Pd
		gamma_energy_6n = 128; //127.8
	}

	//135Te
 	int gamma_energy_1_135Te = 325; //325.0
	int gamma_energy_2_135Te = 1180; //1179.9

	//93Rb
	int gamma_energy_1_93Rb = 373; //372.5
	int gamma_energy_2_93Rb = 746; //746.4, can use this with pretty much everything 

	int gamma_energy_3_93Rb = 552; //551.8
	int gamma_energy_4_93Rb = 733; //733.4

	//94Rb
	int gamma_energy_1_94Rb = 217; //217.2
	int gamma_energy_2_94Rb = 339; //339.2

	//140Xe, no isomer
	int gamma_energy_1_140Xe = 458;
	int gamma_energy_2_140Xe = 377;

	//138Xe, no isomer
	int gamma_energy_1_138Xe = 589;
	int gamma_energy_2_138Xe = 484;

	//92Sr, no isomer
	int gamma_energy_1_92Sr = 815;
	int gamma_energy_2_92Sr = 570;

	//94Sr, no isomer
	int gamma_energy_1_94Sr = 1309;
	int gamma_energy_2_94Sr = 837;


	//////////////////////////////////////////
	/// 		   Lookup table		  	   ///
	//////////////////////////////////////////

	//Implement FWHM-dependent gates later

	fill_lookuptable(gamma_energy_1_134Te, gamma_energy_2_134Te, lookup_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_1n, lookup_1n_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_2n, lookup_2n_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_3n, lookup_3n_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_4n, lookup_4n_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_5n, lookup_5n_134Te);

	fill_lookuptable(gamma_energy_2_134Te, gamma_energy_6n, lookup_6n_134Te);

	fill_lookuptable(gamma_energy_1_135Te, gamma_energy_2_135Te, lookup_135Te);

	fill_lookuptable(gamma_energy_1_93Rb, gamma_energy_2_93Rb, lookup_1_93Rb);

	fill_lookuptable(gamma_energy_3_93Rb, gamma_energy_4_93Rb, lookup_2_93Rb);

	fill_lookuptable(gamma_energy_1_94Rb, gamma_energy_2_94Rb, lookup_94Rb);

	fill_lookuptable(gamma_energy_1_140Xe, gamma_energy_2_140Xe, lookup_140Xe);

	fill_lookuptable(gamma_energy_1_138Xe, gamma_energy_2_138Xe, lookup_138Xe);

	fill_lookuptable(gamma_energy_1_92Sr, gamma_energy_2_92Sr, lookup_92Sr);

	fill_lookuptable(gamma_energy_1_94Sr, gamma_energy_2_94Sr, lookup_94Sr);


	//Print loopup-table to make sure all elements have right value
/*	cout << "94Rb" << " " << endl;
	for(int j=gamma_energy_2_94Rb-10; j<=gamma_energy_2_94Rb+10; j++){
		for(int i=gamma_energy_1_94Rb-10; i<=gamma_energy_1_94Rb+10; i++){	
			cout << lookup_94Rb[i][j] << " ";
		}
		cout << "\n";
	}*/


	//////////////////////////////////////////
	/// 		  	 Sorting	     	   ///
	//////////////////////////////////////////

	for (int k=0; k < TBINS; k++){
		for (int i=0; i < EBINS; i++){
			for (int j=0; j < EBINS; j++){

				int cube_value = Cube1->Get(i,j,k);

				fill_spectra(lookup_134Te[i][j], cube_value, k, time_isomer_doublegate_134Te, time_isomer_doublegate_all_134Te, time_isomer_doublegate_bg_134Te, time_isomer_doublegate_bg_ridge_134Te, time_isomer_doublegate_bg_random_134Te);
				
				fill_spectra(lookup_1n_134Te[i][j], cube_value, k, time_isomer_doublegate_1n_134Te, time_isomer_doublegate_1n_all_134Te, time_isomer_doublegate_1n_bg_134Te, time_isomer_doublegate_1n_bg_ridge_134Te, time_isomer_doublegate_1n_bg_random_134Te);

				fill_spectra(lookup_2n_134Te[i][j], cube_value, k, time_isomer_doublegate_2n_134Te, time_isomer_doublegate_2n_all_134Te, time_isomer_doublegate_2n_bg_134Te, time_isomer_doublegate_2n_bg_ridge_134Te, time_isomer_doublegate_2n_bg_random_134Te);

				fill_spectra(lookup_3n_134Te[i][j], cube_value, k, time_isomer_doublegate_3n_134Te, time_isomer_doublegate_3n_all_134Te, time_isomer_doublegate_3n_bg_134Te, time_isomer_doublegate_3n_bg_ridge_134Te, time_isomer_doublegate_3n_bg_random_134Te);

				fill_spectra(lookup_4n_134Te[i][j], cube_value, k, time_isomer_doublegate_4n_134Te, time_isomer_doublegate_4n_all_134Te, time_isomer_doublegate_4n_bg_134Te, time_isomer_doublegate_4n_bg_ridge_134Te, time_isomer_doublegate_4n_bg_random_134Te);

				fill_spectra(lookup_5n_134Te[i][j], cube_value, k, time_isomer_doublegate_5n_134Te, time_isomer_doublegate_5n_all_134Te, time_isomer_doublegate_5n_bg_134Te, time_isomer_doublegate_5n_bg_ridge_134Te, time_isomer_doublegate_5n_bg_random_134Te);

				fill_spectra(lookup_6n_134Te[i][j], cube_value, k, time_isomer_doublegate_6n_134Te, time_isomer_doublegate_6n_all_134Te, time_isomer_doublegate_6n_bg_134Te, time_isomer_doublegate_6n_bg_ridge_134Te, time_isomer_doublegate_6n_bg_random_134Te);

				fill_spectra(lookup_135Te[i][j], cube_value, k, time_isomer_doublegate_135Te, time_isomer_doublegate_all_135Te, time_isomer_doublegate_bg_135Te, time_isomer_doublegate_bg_ridge_135Te, time_isomer_doublegate_bg_random_135Te);

				fill_spectra(lookup_1_93Rb[i][j], cube_value, k, time_isomer_doublegate_1_93Rb, time_isomer_doublegate_1_all_93Rb, time_isomer_doublegate_1_bg_93Rb, time_isomer_doublegate_1_bg_ridge_93Rb, time_isomer_doublegate_1_bg_random_93Rb);

				fill_spectra(lookup_2_93Rb[i][j], cube_value, k, time_isomer_doublegate_2_93Rb, time_isomer_doublegate_2_all_93Rb, time_isomer_doublegate_2_bg_93Rb, time_isomer_doublegate_2_bg_ridge_93Rb, time_isomer_doublegate_2_bg_random_93Rb);

				fill_spectra(lookup_94Rb[i][j], cube_value, k, time_isomer_doublegate_94Rb, time_isomer_doublegate_all_94Rb, time_isomer_doublegate_bg_94Rb, time_isomer_doublegate_bg_ridge_94Rb, time_isomer_doublegate_bg_random_94Rb);

				fill_spectra(lookup_140Xe[i][j], cube_value, k, time_isomer_doublegate_140Xe, time_isomer_doublegate_all_140Xe, time_isomer_doublegate_bg_140Xe, time_isomer_doublegate_bg_ridge_140Xe, time_isomer_doublegate_bg_random_140Xe);
			
				fill_spectra(lookup_138Xe[i][j], cube_value, k, time_isomer_doublegate_138Xe, time_isomer_doublegate_all_138Xe, time_isomer_doublegate_bg_138Xe, time_isomer_doublegate_bg_ridge_138Xe, time_isomer_doublegate_bg_random_138Xe);	

				fill_spectra(lookup_92Sr[i][j], cube_value, k, time_isomer_doublegate_92Sr, time_isomer_doublegate_all_92Sr, time_isomer_doublegate_bg_92Sr, time_isomer_doublegate_bg_ridge_92Sr, time_isomer_doublegate_bg_random_92Sr);	

				fill_spectra(lookup_94Sr[i][j], cube_value, k, time_isomer_doublegate_94Sr, time_isomer_doublegate_all_94Sr, time_isomer_doublegate_bg_94Sr, time_isomer_doublegate_bg_ridge_94Sr, time_isomer_doublegate_bg_random_94Sr);	
			
			}
		}
	}

	//////////////////////////////////////////
	/// 	  Write spectra to file	       ///
	//////////////////////////////////////////

	TFile *outputspectrafile = new TFile("CubeSort.root","RECREATE");

	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();
	time_isomer_doublegate_bg_134Te->Write();
	time_isomer_doublegate_bg_ridge_134Te->Write();
	time_isomer_doublegate_bg_random_134Te->Write();

	time_isomer_doublegate_1n_134Te->Write();
	time_isomer_doublegate_1n_all_134Te->Write();
	time_isomer_doublegate_1n_bg_134Te->Write();
	time_isomer_doublegate_1n_bg_ridge_134Te->Write();
	time_isomer_doublegate_1n_bg_random_134Te->Write();

	time_isomer_doublegate_2n_134Te->Write();
	time_isomer_doublegate_2n_all_134Te->Write();
	time_isomer_doublegate_2n_bg_134Te->Write();
	time_isomer_doublegate_2n_bg_ridge_134Te->Write();
	time_isomer_doublegate_2n_bg_random_134Te->Write();

	time_isomer_doublegate_3n_134Te->Write();
	time_isomer_doublegate_3n_all_134Te->Write();
	time_isomer_doublegate_3n_bg_134Te->Write();
	time_isomer_doublegate_3n_bg_ridge_134Te->Write();
	time_isomer_doublegate_3n_bg_random_134Te->Write();

	time_isomer_doublegate_4n_134Te->Write();
	time_isomer_doublegate_4n_all_134Te->Write();
	time_isomer_doublegate_4n_bg_134Te->Write();
	time_isomer_doublegate_4n_bg_ridge_134Te->Write();
	time_isomer_doublegate_4n_bg_random_134Te->Write();

	time_isomer_doublegate_5n_134Te->Write();
	time_isomer_doublegate_5n_all_134Te->Write();
	time_isomer_doublegate_5n_bg_134Te->Write();
	time_isomer_doublegate_5n_bg_ridge_134Te->Write();
	time_isomer_doublegate_5n_bg_random_134Te->Write();

	time_isomer_doublegate_6n_134Te->Write();
	time_isomer_doublegate_6n_all_134Te->Write();
	time_isomer_doublegate_6n_bg_134Te->Write();
	time_isomer_doublegate_6n_bg_ridge_134Te->Write();
	time_isomer_doublegate_6n_bg_random_134Te->Write();

	time_isomer_doublegate_135Te->Write();
	time_isomer_doublegate_all_135Te->Write();
	time_isomer_doublegate_bg_135Te->Write();
	time_isomer_doublegate_bg_ridge_135Te->Write();
	time_isomer_doublegate_bg_random_135Te->Write();

	time_isomer_doublegate_1_93Rb->Write();
	time_isomer_doublegate_1_all_93Rb->Write();
	time_isomer_doublegate_1_bg_93Rb->Write();
	time_isomer_doublegate_1_bg_ridge_93Rb->Write();
	time_isomer_doublegate_1_bg_random_93Rb->Write();

	time_isomer_doublegate_2_93Rb->Write();
	time_isomer_doublegate_2_all_93Rb->Write();
	time_isomer_doublegate_2_bg_93Rb->Write();
	time_isomer_doublegate_2_bg_ridge_93Rb->Write();
	time_isomer_doublegate_2_bg_random_93Rb->Write();

	time_isomer_doublegate_94Rb->Write();
	time_isomer_doublegate_all_94Rb->Write();
	time_isomer_doublegate_bg_94Rb->Write();
	time_isomer_doublegate_bg_ridge_94Rb->Write();
	time_isomer_doublegate_bg_random_94Rb->Write();

	time_isomer_doublegate_140Xe->Write();
	time_isomer_doublegate_all_140Xe->Write();
	time_isomer_doublegate_bg_140Xe->Write();
	time_isomer_doublegate_bg_ridge_140Xe->Write();
	time_isomer_doublegate_bg_random_140Xe->Write();

	time_isomer_doublegate_138Xe->Write();
	time_isomer_doublegate_all_138Xe->Write();
	time_isomer_doublegate_bg_138Xe->Write();
	time_isomer_doublegate_bg_ridge_138Xe->Write();
	time_isomer_doublegate_bg_random_138Xe->Write();

	time_isomer_doublegate_92Sr->Write();
	time_isomer_doublegate_all_92Sr->Write();
	time_isomer_doublegate_bg_92Sr->Write();
	time_isomer_doublegate_bg_ridge_92Sr->Write();
	time_isomer_doublegate_bg_random_92Sr->Write();

	time_isomer_doublegate_94Sr->Write();
	time_isomer_doublegate_all_94Sr->Write();
	time_isomer_doublegate_bg_94Sr->Write();
	time_isomer_doublegate_bg_ridge_94Sr->Write();
	time_isomer_doublegate_bg_random_94Sr->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();
}







void fill_lookuptable(int energy_1, int energy_2, short lookup[2048][2048]){
	
	//a b c
	//d e f
	//g h i

	//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
	//bg_ridge = (b+d+f+h)/2
	//bg_random = (a+c+g+i)/4

	//Number of FWHM in the energy gate
	double NFWHM=1.0;

	//Jon's parameterisation of the energy-dependence of the FWHM
	double A=1.059;
	double B=2.814;
	double gatewidth_1_float = (0.5+(sqrt(A*A+(B*B*energy_1/1000.0))/1.0))*(NFWHM/2.0); //1/2 FWHM.
	double gatewidth_2_float = (0.5+(sqrt(A*A+(B*B*energy_2/1000.0))/1.0))*(NFWHM/2.0); //1/2 FWHM.

	int gatewidth_1 = round(gatewidth_1_float);
	int gatewidth_2 = round(gatewidth_2_float);

/*	std::cout << "Number of channels in true gate for energy_1: " << energy_1 << " keV " << gatewidth_1 << std::endl;
	std::cout << "Number of channels in true gate for energy_2: " << energy_2 << " keV " << gatewidth_2 << std::endl;
*/
	gatewidth_1 = 2; // energy +/- gatewidth gives total width of 5
	gatewidth_2 = 2;

	//fill whole square with random bg-nr
	for(int i=energy_1-3*gatewidth_1-1; i<=energy_1+3*gatewidth_1+1; i++){
		for(int j=energy_2-3*gatewidth_2-1; j<=energy_2+3*gatewidth_2+1; j++){
			lookup[i][j] = 3;
		}
	}

	//peak gate; e
	for(int i=energy_1-gatewidth_1; i<=energy_1+gatewidth_1; i++){
		for(int j=energy_2-gatewidth_2; j<=energy_2+gatewidth_2; j++){
			lookup[i][j] = 2;
		}
	}

	//bg-ridge gate, energy_1 i bg_lower; d
	for(int i=energy_1-3*gatewidth_1-1; i<=energy_1-gatewidth_1-1; i++){
		for(int j=energy_2-gatewidth_2; j<=energy_2+gatewidth_2; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_1 i bg_upper; f
	for(int i=energy_1+gatewidth_1+1; i<=energy_1+3*gatewidth_1+1; i++){
		for(int j=energy_2-gatewidth_2; j<=energy_2+gatewidth_2; j++){
			lookup[i][j] = 1;
		}
	} 

	//bg-ridge gate, energy_2 i bg_lower; h
	for(int i=energy_1-gatewidth_1; i<=energy_1+gatewidth_1; i++){
		for(int j=energy_2-3*gatewidth_2-1; j<=energy_2-gatewidth_2-1; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_2 i bg_upper; b
	for(int i=energy_1-gatewidth_1; i<=energy_1+gatewidth_1; i++){
		for(int j=energy_2+gatewidth_2+1; j<=energy_2+3*gatewidth_2+1; j++){
			lookup[i][j] = 1;
		}
	}
}






void fill_spectra(int lookup_value, int cube_value, int k, TH1D *true_spec, TH1D *all_spec, TH1D *bg_spec, TH1D *bg_ridge_spec, TH1D *bg_random_spec){

	//a b c
	//d e f
	//g h i

	//2: true = e - (b+d+f+h)/2 + (a+c+g+i)/4
	//1: bg_ridge = (b+d+f+h)/2
	//3: bg_random = (a+c+g+i)/4


	//e
	if(lookup_value==2){
		true_spec->Fill(k, cube_value);
		all_spec->Fill(k, cube_value);
	}

	//b, d, f, h, only subtract half
	if(lookup_value==1){
		true_spec->Fill(k, -0.5*cube_value);
		bg_spec->Fill(k, 0.5*cube_value);
		bg_ridge_spec->Fill(k, 0.5*cube_value);
	}

	//a, c, g, i, add one quarter
	if(lookup_value==3){
		true_spec->Fill(k, 0.25*cube_value);
		bg_spec->Fill(k, -0.25*cube_value);
		bg_random_spec->Fill(k, 0.25*cube_value);
	}
}


//g++ -g -o CubeSort CubeSort.cxx ` root-config --cflags` `root-config --glibs` -lSpectrum
