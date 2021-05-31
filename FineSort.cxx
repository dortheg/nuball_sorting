#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include "TImage.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TROOT.h"
#include "TColor.h"
#include "TAttImage.h"
#include "TEnv.h"
#include "TH2D.h"
#include "TH3S.h"
#include "TF2.h"
#include "TColor.h"
#include "TLine.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "Spec.hxx"
#include "Spec.cxx"
#include "RWMat.hxx"
#include "RWMat.cxx"
#include "WriteRWSpec.cxx"
#include "fwhmIC.hxx"
#include "TGraphErrors.h"
#include "CorBGO.hxx"
#include "idmap.hxx"
#include "TimeOffsets.hxx"
#include "CalBGO.hxx"
#include "EuPeaksBis.hxx" //Ge calibration
#include "NewLaBr3Calb.hxx" //LaBr3 calibration
#include "Names.cxx"
#include "Walks2.cxx"
#include "Egate.cxx"
#include "CubeDDT.cxx"
#include "CalCorrection.hxx"
#include <iostream>
#include <iomanip>

using namespace std;

unsigned int bsize=4000000000;
UShort_t *TheEvents_1 = new UShort_t[bsize];
UInt_t *tpointer_array_1 = new UInt_t[2];

void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_bg(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_all(double energy_isomer, double energy, double time, TH1D *h1);

void check_energy_doublegate_true(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_all(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);

int test_both_isomer(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2);


int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_FineSort.cxx"

	//Read TheEvents_1 from file
	ifstream infile_edata("edata.txt", ios::in | ios::binary);
	infile_edata.read((char *) TheEvents_1, bsize*sizeof(UShort_t));
	infile_edata.close();

	//Read tpointer_array from file
	ifstream tpointer_infile("tpointer.txt", ios::out | ios::binary);
	tpointer_infile.read((char *) tpointer_array_1, 2*sizeof(UInt_t));
	tpointer_infile.close();
	std::cout << "tpointer_1: " << tpointer_array_1[0] << std::endl;


	//////////////////////////////////////////
	/// 		  Isomer properties	       ///
	//////////////////////////////////////////

	//134Te
	double isomer_energy_1_134Te = 297.0;
	double isomer_energy_2_134Te = 1279.0;

	//////////////////////////////////////////
	/// 		  Create spectra	       ///
	//////////////////////////////////////////

	int aval_1, mult_1;
	double tot_mult_1 = 0;

	int i_count = 0;
	while(i_count<tpointer_array_1[0]){

		aval_1=TheEvents_1[i_count++];
		mult_1=TheEvents_1[i_count++];

		tot_mult_1 += mult_1;

		double energy_1[mult_1];
		double time_1[mult_1];

		for (int k=0; k < mult_1; k++){
			energy_1[k] = TheEvents_1[i_count++];
			time_1[k] = TheEvents_1[i_count++];

			single_gamma->Fill(energy_1[k]);

			//134Te
			check_energy_singlegate_true(isomer_energy_1_134Te, energy_1[k], time_1[k], time_isomer_gate_134Te);
			check_energy_singlegate_bg(isomer_energy_1_134Te, energy_1[k], time_1[k], time_isomer_gate_bg_134Te);
			check_energy_singlegate_all(isomer_energy_1_134Te, energy_1[k], time_1[k], time_isomer_gate_all_134Te);
		}

		int test_both_isomer_value = 0;
		//Loop over all pairs of gammas
		if(mult_1>1){
			for(int m=0; m < mult_1-1; m++){
				for(int n=m+1; n < mult_1; n++ ){
					double_gamma->Fill(energy_1[m], energy_1[n]);

					check_energy_doublegate_true(isomer_energy_1_134Te, isomer_energy_2_134Te, energy_1[m], energy_1[n], time_1[m], time_1[n], time_isomer_doublegate_134Te);
					check_energy_doublegate_bg(isomer_energy_1_134Te, isomer_energy_2_134Te, energy_1[m], energy_1[n], time_1[m], time_1[n], time_isomer_doublegate_bg_134Te);
					check_energy_doublegate_all(isomer_energy_1_134Te, isomer_energy_2_134Te, energy_1[m], energy_1[n], time_1[m], time_1[n], time_isomer_doublegate_all_134Te);
				
					test_both_isomer_value += test_both_isomer(isomer_energy_1_134Te, isomer_energy_2_134Te, energy_1[m], energy_1[n]);

				}
			}
		}

		if(test_both_isomer_value>0){
			for(int j=0;j<mult_1;j++){
				all_gammas_doublegated_134Te->Fill(energy_1[j]);
			}
		}
	}

	std::cout << "tot mult after writing: " << tot_mult_1 << std::endl;
	std::cout << "\n" << std::endl;


	//////////////////////////////////////////
	/// 	      Write spectra	           ///
	//////////////////////////////////////////

	string OutputSpectraFile="FineSort.root";
	TFile *outputspectrafile = new TFile(OutputSpectraFile.c_str(),"RECREATE");	
	
	single_gamma->Write();
	double_gamma->Write();

	time_isomer_gate_134Te->Write();
	time_isomer_gate_bg_134Te->Write();
	time_isomer_gate_all_134Te->Write(); 

	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_bg_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();


}

void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1){
	//BG, smaller than energy_isomer, subtract
	if( (energy>=energy_isomer-5.0) && (energy<energy_isomer-2.5) ){
		h1->Fill(time,-1);
	}	
	//True, inside energy gate, add
	if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		h1->Fill(time);
	}
	//BG, larger than energy_isomer, subtract
	if( (energy>=energy_isomer+2.5) && (energy<energy_isomer+5.0) ){
		h1->Fill(time,-1);
	}
}

void check_energy_singlegate_bg(double energy_isomer, double energy, double time, TH1D *h1){
	//BG, smaller than energy_isomer, add
	if( (energy>=energy_isomer-5.0) && (energy<energy_isomer-2.5) ){
		h1->Fill(time);
	}	
	//BG, larger than energy_isomer, add
	if( (energy>=energy_isomer+2.5) && (energy<energy_isomer+5.0) ){
		h1->Fill(time);
	}
}

void check_energy_singlegate_all(double energy_isomer, double energy, double time, TH1D *h1){
	//True, inside energy gate, add
	if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		h1->Fill(time);
	}
}


void check_energy_doublegate_true(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1){
	
	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////
	
	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1);
	}
	//BG 1: energy_1 in bg below isomer_A, energy_2 in isomer_B 
	if( (energy_1>=energy_isomer_A-5.0) && (energy_1<energy_isomer_A-2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 2: energy_1 in bg above isomer_A, energy_2 in isomer_B 
	if( (energy_1>=energy_isomer_A+2.5) && (energy_1<energy_isomer_A+5.0) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 3: energy_1 in isomer_A and energy_2 in bg below isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-5.0) && (energy_2<energy_isomer_B-2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 4: energy_1 in isomer_A and energy_2 in bg above isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B+2.5) && (energy_2<energy_isomer_B+5.0)){
		h1->Fill(time_1, -0.25);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_B and energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1);
	}
	//BG 1: energy_1 in bg below isomer_B, energy_2 in isomer_A 
	if( (energy_1>=energy_isomer_B-5.0) && (energy_1<energy_isomer_B-2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 2: energy_1 in bg above isomer_B, energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B+2.5) && (energy_1<energy_isomer_B+5.0) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 3: energy_1 in isomer_B and energy_2 in bg below isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-5.0) && (energy_2<energy_isomer_A-2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 4: energy_1 in isomer_B and energy_2 in bg above isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A+2.5) && (energy_2<energy_isomer_A+5.0)){
		h1->Fill(time_1, -0.25);
	}
}

void check_energy_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1){

	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////

	//BG 1: energy_1 in bg below isomer_A, energy_2 in isomer_B 
	if( (energy_1>=energy_isomer_A-5.0) && (energy_1<energy_isomer_A-2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 2: energy_1 in bg above isomer_A, energy_2 in isomer_B 
	if( (energy_1>=energy_isomer_A+2.5) && (energy_1<energy_isomer_A+5.0) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 3: energy_1 in isomer_A and energy_2 in bg below isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-5.0) && (energy_2<energy_isomer_B-2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 4: energy_1 in isomer_A and energy_2 in bg above isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B+2.5) && (energy_2<energy_isomer_B+5.0)){
		h1->Fill(time_1, 0.25);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//BG 1: energy_1 in bg below isomer_B, energy_2 in isomer_A 
	if( (energy_1>=energy_isomer_B-5.0) && (energy_1<energy_isomer_B-2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 2: energy_1 in bg above isomer_B, energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B+2.5) && (energy_1<energy_isomer_B+5.0) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 3: energy_1 in isomer_B and energy_2 in bg below isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-5.0) && (energy_2<energy_isomer_A-2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 4: energy_1 in isomer_B and energy_2 in bg above isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A+2.5) && (energy_2<energy_isomer_A+5.0)){
		h1->Fill(time_1, 0.25);
	}
}


void check_energy_doublegate_all(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1){
	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////
	
	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_B and energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1);
	}
}


int test_both_isomer(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2){
	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		return 1;
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_B and energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		return 1;
	}

	else{
		return 0;
	}
}


//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`