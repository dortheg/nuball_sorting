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
short lookup_isomer_134Te[2048][2048] = {{0}};

int coincidences[20][2];

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

 	int gamma_energy_1_134Te = 297;
	int gamma_energy_2_134Te = 1279;
	int gamma_energy_3_134Te = 115;



	//////////////////////////////////////////
	/// 		   Lookup table		  	   ///
	//////////////////////////////////////////

	//Implement FWHM-dependent gates later

	fill_lookuptable(gamma_energy_1_134Te, gamma_energy_2_134Te, lookup_134Te);
	//fill_lookuptable(gamma_energy_2_134Te, gamma_energy_3_134Te, lookup_134Te);


/*	//Print loopup-table to make sure all elements have right value
	for(int j=gamma_energy_2_134Te-5; j<=gamma_energy_2_134Te+6; j++){
		for(int i=gamma_energy_1_134Te-5; i<=gamma_energy_1_134Te+6; i++){
			//cout << lookup[i][j] << " ";
			cout << lookup_134Te[i][j] << " ";
		}
		cout << "\n";
	}*/


	//////////////////////////////////////////
	/// 		  	 Sorting	     	   ///
	//////////////////////////////////////////

	for (int k=0; k < TBINS; k++){
		for (int i=0; i < EBINS; i++){
			for (int j=0; j < EBINS; j++){

				//a b c
				//d e f
				//g h i

				//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
				//bg_ridge = (b+d+f+h)/2
				//bg_random = (a+c+g+i)/4

				//Also move this stuff to function

				//e
				if(lookup_134Te[i][j]==2){
					time_isomer_doublegate_134Te->Fill(k, Cube1->Get(i,j,k));
					time_isomer_doublegate_all_134Te->Fill(k, Cube1->Get(i,j,k));
				}

				//b, d, f, h
				if(lookup_134Te[i][j]==1){
					time_isomer_doublegate_134Te->Fill(k, -0.5*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_134Te->Fill(k, 0.5*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_ridge_134Te->Fill(k, 0.5*Cube1->Get(i,j,k));
				}

				//a, c, g, i
				if(lookup_134Te[i][j]==3){
					time_isomer_doublegate_134Te->Fill(k, 0.25*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_134Te->Fill(k, -0.25*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_random_134Te->Fill(k, 0.25*Cube1->Get(i,j,k));
				}

				//int cube_value = Cube1->Get(i,j,k);

				//fill_spectra(lookup_134Te[i][j], cube_value, k, time_isomer_doublegate_134Te, time_isomer_doublegate_all_134Te, time_isomer_doublegate_bg_134Te, time_isomer_doublegate_bg_ridge_134Te, time_isomer_doublegate_bg_random_134Te);
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

	//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
	//bg_ridge = (b+d+f+h)/2
	//bg_random = (a+c+g+i)/4


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
