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
//#include "Spec.hxx"
//#include "Spec.cxx"
//#include "RWMat.hxx"
//#include "RWMat.cxx"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "CubeDDT.hxx"
#include "CubeDDT.cxx"
//#include "kr94isob.hxx"
//#include "sr95iso.hxx"
//#include "ba144.hxx"

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

/*	//Create & read in cubes
	CubeDDT *Cube1=new CubeDDT("",EBINS,TBINS,10);
	CubeDDT *Cube2=new CubeDDT("",EBINS,TBINS,10);
	Cube1->Read("U238cube_n3_4jan2021.bin");
	Cube2->Read("U238cube_n3_4jan2021-sub.bin");*/

	// Create coincidence matrix
	TH2F* ggm = new TH2F("ggm","ggm",EBINS,0,EBINS,EBINS,0,EBINS);

	//Create time spectrum
	TH1D **tspecs; tspecs = new TH1D*[EBINS];
	for (int i=1; i <= TBINS; i++){
 		TString name=Form("t%d",i);
 		tspecs[i]=new TH1D(name.Data(),name.Data(),EBINS,0,EBINS);
 	}

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

	

}

//g++ -g -o CubeSort CubeSort.cxx ` root-config --cflags` `root-config --glibs` -lSpectrum
