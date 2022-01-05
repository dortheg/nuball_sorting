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

short lookup[2048][2048] = {0};
int coincidences[20][2];
int lookup1D[2048] = {0};

//Number of FWHM in the energy gate
double NFWHM=1.0;


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

	//Implement FWHM-dependent gates later


	//////////////////////////////////////////
	/// 		   Lookup table		  	   ///
	//////////////////////////////////////////

	//a b c
	//d e f
	//g h i

	//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
	//bg_ridge = (b+d+f+h)/2
	//bg_random = (a+c+g+i)/4

	//fill whole square with random bg-nr
	for(int i=gamma_energy_1_134Te-4; i<=gamma_energy_1_134Te+5; i++){
		for(int j=gamma_energy_2_134Te-4; j<=gamma_energy_2_134Te+5; j++){
			lookup[i][j] = 3;
		}
	}

	//peak gate; e
	for(int i=gamma_energy_1_134Te-2; i<=gamma_energy_1_134Te+2; i++){
		for(int j=gamma_energy_2_134Te-2; j<=gamma_energy_2_134Te+2; j++){
			lookup[i][j] = 2;
		}
	}

	//bg-ridge gate, energy_1 i bg_lower; d
	for(int i=gamma_energy_1_134Te-4; i<=gamma_energy_1_134Te-3; i++){
		for(int j=gamma_energy_2_134Te-2; j<=gamma_energy_2_134Te+2; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_1 i bg_upper; f
	for(int i=gamma_energy_1_134Te+3; i<=gamma_energy_1_134Te+5; i++){
		for(int j=gamma_energy_2_134Te-2; j<=gamma_energy_2_134Te+2; j++){
			lookup[i][j] = 1;
		}
	} 

	//bg-ridge gate, energy_2 i bg_lower; h
	for(int i=gamma_energy_1_134Te-2; i<=gamma_energy_1_134Te+2; i++){
		for(int j=gamma_energy_2_134Te-4; j<=gamma_energy_2_134Te-3; j++){
			lookup[i][j] = 1;
		}
	}

	//bg-ridge gate, energy_2 i bg_upper; b
	for(int i=gamma_energy_1_134Te-2; i<=gamma_energy_1_134Te+2; i++){
		for(int j=gamma_energy_2_134Te+3; j<=gamma_energy_2_134Te+5; j++){
			lookup[i][j] = 1;
		}
	}

/*	//Print loopup-table to make sure all elements have right value
	for(int j=gamma_energy_2_134Te-5; j<=gamma_energy_2_134Te+6; j++){
		for(int i=gamma_energy_1_134Te-5; i<=gamma_energy_1_134Te+6; i++){
			cout << lookup[i][j] << " ";
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

				//e
				if(lookup[i][j]==2){
					time_isomer_doublegate_134Te->Fill(k, Cube1->Get(i,j,k));
					time_isomer_doublegate_all_134Te->Fill(k, Cube1->Get(i,j,k));
				}

				//b, d, f, h
				if(lookup[i][j]==1){
					time_isomer_doublegate_134Te->Fill(k, -0.5*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_134Te->Fill(k, 0.5*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_ridge_134Te->Fill(k, 0.5*Cube1->Get(i,j,k));
				}

				//a, c, g, i
				if(lookup[i][j]==3){
					time_isomer_doublegate_134Te->Fill(k, 0.25*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_134Te->Fill(k, -0.25*Cube1->Get(i,j,k));
					time_isomer_doublegate_bg_random_134Te->Fill(k, 0.25*Cube1->Get(i,j,k));
				}
			}
		}
	}

	//////////////////////////////////////////
	/// 	  Write spectra to file	       ///
	//////////////////////////////////////////

	TFile *outputspectrafile = new TFile("CubeSort.root","RECREATE");
	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();


}

//g++ -g -o CubeSort CubeSort.cxx ` root-config --cflags` `root-config --glibs` -lSpectrum
