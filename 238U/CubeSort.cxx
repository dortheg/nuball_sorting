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
short lookup_140Xe[2048][2048] = {{0}};
short lookup_138Xe[2048][2048] = {{0}};
short lookup_92Sr[2048][2048] = {{0}};
short lookup_94Sr[2048][2048] = {{0}};

//Number of FWHM in the energy gate
double NFWHM=1.0;


void fill_lookuptable(int energy_1, int energy_2, short lookup[2048][2048]);
void fill_spectra(int lookup_value, int cube_value, int k, TH1D *true_spec, TH1D *all_spec, TH1D *bg_spec, TH1D *bg_ridge_spec, TH1D *bg_random_spec);


int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_CubeSort.cxx"

	//Create & read in cubes
	CubeDDT *Cube1=new CubeDDT("",EBINS,TBINS,2);
	Cube1->Read("U238cube_n3_2ns_4jan2021.bin");


 	//////////////////////////////////////////
	/// 		  Isomer properties	       ///
	//////////////////////////////////////////

	//134Te
 	int gamma_energy_1_134Te = 297;
	int gamma_energy_2_134Te = 1279;
	int gamma_energy_3_134Te = 115;

	//140Xe
	int gamma_energy_1_140Xe = 458;
	int gamma_energy_2_140Xe = 377;

	//138Xe
	int gamma_energy_1_138Xe = 589;
	int gamma_energy_2_138Xe = 484;

	//92Sr
	int gamma_energy_1_92Sr = 815;
	int gamma_energy_2_92Sr = 570;

	//94Sr
	int gamma_energy_1_94Sr = 1309;
	int gamma_energy_2_94Sr = 837;


	//////////////////////////////////////////
	/// 		   Lookup table		  	   ///
	//////////////////////////////////////////

	//Implement FWHM-dependent gates later

	fill_lookuptable(gamma_energy_1_134Te, gamma_energy_2_134Te, lookup_134Te);

	fill_lookuptable(gamma_energy_1_140Xe, gamma_energy_2_140Xe, lookup_140Xe);

	fill_lookuptable(gamma_energy_1_138Xe, gamma_energy_2_138Xe, lookup_138Xe);

	fill_lookuptable(gamma_energy_1_92Sr, gamma_energy_2_92Sr, lookup_92Sr);

	fill_lookuptable(gamma_energy_1_94Sr, gamma_energy_2_94Sr, lookup_94Sr);


/*	//Print loopup-table to make sure all elements have right value
	for(int j=gamma_energy_2_140Xe-5; j<=gamma_energy_2_140Xe+6; j++){
		for(int i=gamma_energy_1_140Xe-5; i<=gamma_energy_1_140Xe+6; i++){
			//cout << lookup_134Te[i][j] << " ";
			cout << lookup_140Xe[i][j] << " ";
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

	//fill whole square with random bg-nr
	for(int i=energy_1-4; i<=energy_1+5; i++){
		for(int j=energy_2-4; j<=energy_2+5; j++){
			lookup[i][j] = 3;
		}
	}

	//peak gate; e
	for(int i=energy_1-2; i<=energy_1+2; i++){
		for(int j=energy_2-2; j<=energy_2+2; j++){
			lookup[i][j] = 2;
		}
	}

	//bg-ridge gate, energy_1 i bg_lower; d
	for(int i=energy_1-4; i<=energy_1-3; i++){
		for(int j=energy_2-2; j<=energy_2+2; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_1 i bg_upper; f
	for(int i=energy_1+3; i<=energy_1+5; i++){
		for(int j=energy_2-2; j<=energy_2+2; j++){
			lookup[i][j] = 1;
		}
	} 

	//bg-ridge gate, energy_2 i bg_lower; h
	for(int i=energy_1-2; i<=energy_1+2; i++){
		for(int j=energy_2-4; j<=energy_2-3; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_2 i bg_upper; b
	for(int i=energy_1-2; i<=energy_1+2; i++){
		for(int j=energy_2+3; j<=energy_2+5; j++){
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
