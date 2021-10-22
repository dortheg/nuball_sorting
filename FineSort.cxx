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
#include <chrono>

using namespace std;
using namespace std::chrono;

unsigned int bsize=4000000000;
UShort_t *TheEvents = new UShort_t[bsize];
UInt_t *tpointer_array = new UInt_t[2];

void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_bg(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_all(double energy_isomer, double energy, double time, TH1D *h1);

void check_energy_doublegate_true(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_all(int energy_isomer_A, int energy_isomer_B, int energy_1, int energy_2, int time_1, int time_2, TH1D *h1);


int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_FineSort.cxx"

	auto start_program = high_resolution_clock::now();

	// //Read TheEvents from file
	string OutputDirectory="/Applications/nuball_sorting/IYR_Data/";
	std::cout << "Reading data from file " << std::endl;
	ifstream infile_edata("/Applications/nuball_sorting/IYR_Data/edata_allfiles_10sep2021.txt", ios::in | ios::binary);
	infile_edata.read((char *) TheEvents, bsize*sizeof(UShort_t));
	infile_edata.close();

	//Read tpointer_array from file
	ifstream tpointer_infile("/Applications/nuball_sorting/IYR_Data/tpointer_allfiles_10sep2021.txt", ios::out | ios::binary);
	tpointer_infile.read((char *) tpointer_array, 2*sizeof(UInt_t));
	tpointer_infile.close();
	std::cout << "tpointer: " << tpointer_array[0] << std::endl;
	std::cout << "Sorting data " << std::endl;


	//////////////////////////////////////////
	/// 		  Isomer properties	       ///
	//////////////////////////////////////////

	//134Te
	int isomer_energy_1_134Te = 297;
	int isomer_energy_2_134Te = 1279;

	//////////////////////////////////////////
	/// 		   Lookup tables		   ///
	//////////////////////////////////////////

	//134Te, isomer_1
	int lookup_134Te_isomer_1[65535] = {0};
	//peak gate
	lookup_134Te_isomer_1[isomer_energy_1_134Te] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-1] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-2] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+1] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+2] = 2;
	//bg gate
	lookup_134Te_isomer_1[isomer_energy_1_134Te-3] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-4] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-5] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+3] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+4] = 1;

	//134Te, isomer_2
	int lookup_134Te_isomer_2[65535] = {0};
	//peak gate
	lookup_134Te_isomer_2[isomer_energy_2_134Te] = 2;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-1] = 2;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-2] = 2;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+1] = 2;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+2] = 2;
	//bg gate
	lookup_134Te_isomer_2[isomer_energy_2_134Te-3] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-4] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-5] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+3] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+4] = 1;


	//////////////////////////////////////////
	/// 		  	 Sorting	     	   ///
	//////////////////////////////////////////

	double tot_mult = 0;
	UShort_t avalmult;
	unsigned char new_aval, mult;
	UShort_t aval;

	auto start_sorting = high_resolution_clock::now();

	//dummy multiplicity data
	//single 1 true: 2
	//single 1 all: 2
	//single 1 bg: 1

	//sigle 2 true: 2
	//single 2 all: 2
	//single 2 bg: 0
	//int dummy_data[23] = {3, 150, 1000, 1100, 1000, 297, 1000, 2, 298, 1000, 1280, 1000, 5, 600, 1000, 400, 1000, 294, 1000, 100, 1000, 1279, 1000};

	int i_count = 0;
	while(i_count<tpointer_array[0]){
	//while(i_count<23){

		//Decompress aval and mult
		avalmult = TheEvents[i_count++];
		new_aval = avalmult & 0xFF;
		mult = avalmult >> 8;
		tot_mult += mult;
		aval = new_aval*5; //Because divided aval by 5

		//mult = dummy_data[i_count++];

		int energy[mult];
		int time[mult];

		//Loop over single gammas
		for (int k=0; k < mult; k++){
			energy[k] = TheEvents[i_count++];
			time[k] = TheEvents[i_count++];

			//energy[k] = dummy_data[i_count++];
			//time[k] = dummy_data[i_count++];

			single_gamma->Fill(energy[k]);

			//134Te, iftest isomer 1
			check_energy_singlegate_true(isomer_energy_1_134Te, energy[k], time[k], time_isomer_1_iftest_gate_134Te);
			check_energy_singlegate_bg(isomer_energy_1_134Te, energy[k], time[k], time_isomer_1_iftest_gate_bg_134Te);
			check_energy_singlegate_all(isomer_energy_1_134Te, energy[k], time[k], time_isomer_1_iftest_gate_all_134Te);
			
			//134Te, iftest isomer 2
			check_energy_singlegate_true(isomer_energy_2_134Te, energy[k], time[k], time_isomer_2_iftest_gate_134Te);
			check_energy_singlegate_bg(isomer_energy_2_134Te, energy[k], time[k], time_isomer_2_iftest_gate_bg_134Te);
			check_energy_singlegate_all(isomer_energy_2_134Te, energy[k], time[k], time_isomer_2_iftest_gate_all_134Te);
			

			//134Te, isomer_1
			if(lookup_134Te_isomer_1[energy[k]]==2){
				time_isomer_1_gate_134Te->Fill(time[k]);
				time_isomer_1_gate_all_134Te->Fill(time[k]);
			}
			else if(lookup_134Te_isomer_1[energy[k]]==1){
				time_isomer_1_gate_134Te->Fill(time[k],-1);
				time_isomer_1_gate_bg_134Te->Fill(time[k]);
			}

			//134Te, isomer_2
			if(lookup_134Te_isomer_2[energy[k]]==2){
				time_isomer_2_gate_134Te->Fill(time[k]);
				time_isomer_2_gate_all_134Te->Fill(time[k]);
			}
			else if(lookup_134Te_isomer_2[energy[k]]==1){
				time_isomer_2_gate_134Te->Fill(time[k],-1);
				time_isomer_2_gate_bg_134Te->Fill(time[k]);
			}
		}

		//Loop over all pairs of gammas
		if(mult>1){
			for(int m=0; m < mult; m++){
				for(int n=0; n < mult; n++){
					if(m!=n){

						double_gamma->Fill(energy[m], energy[n]);

						//123Te doublegate
						if(lookup_134Te_isomer_1[energy[m]]==2 && lookup_134Te_isomer_2[energy[n]]==2){
							time_isomer_doublegate_134Te->Fill(time[m]);
							time_isomer_doublegate_all_134Te->Fill(time[m]);
						}

						else if(lookup_134Te_isomer_1[energy[m]]==2 && lookup_134Te_isomer_2[energy[n]]==1){
							time_isomer_doublegate_134Te->Fill(time[m], -0.25);
							time_isomer_doublegate_bg_134Te->Fill(time[m]);
						}

						else if(lookup_134Te_isomer_1[energy[m]]==1 && lookup_134Te_isomer_2[energy[n]]==2){
							time_isomer_doublegate_134Te->Fill(time[m], -0.25);
							time_isomer_doublegate_bg_134Te->Fill(time[m]);
						}
					}
				}
			}
		}

		if(mult>1){
			for(int m=0; m < mult-1; m++){
				for(int n=m+1; n < mult; n++ ){
					//double_gamma->Fill(energy[m], energy[n]);

					check_energy_doublegate_true(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_iftest_doublegate_134Te);
					check_energy_doublegate_bg(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_iftest_doublegate_bg_134Te);
					check_energy_doublegate_all(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_iftest_doublegate_all_134Te);
				}
			}
		}
	}

	auto stop_sorting = high_resolution_clock::now();

	auto duration_sorting = duration_cast<milliseconds>(stop_sorting - start_sorting);
	cout << "Runtime sorting: " << duration_sorting.count() << endl;

	//std::cout << "tot mult after writing: " << tot_mult << std::endl;
	//std::cout << "\n" << std::endl;


	//////////////////////////////////////////
	/// 	      Write spectra	           ///
	//////////////////////////////////////////

	string OutputSpectraFile=OutputDirectory+"FineSort.root";
	TFile *outputspectrafile = new TFile(OutputSpectraFile.c_str(),"RECREATE");	
	
	single_gamma->Write();
	double_gamma->Write();

	time_isomer_1_gate_134Te->Write();
	time_isomer_1_gate_bg_134Te->Write();
	time_isomer_1_gate_all_134Te->Write(); 

	time_isomer_1_iftest_gate_134Te->Write();
	time_isomer_1_iftest_gate_bg_134Te->Write();
	time_isomer_1_iftest_gate_all_134Te->Write(); 

	time_isomer_2_gate_134Te->Write();
	time_isomer_2_gate_bg_134Te->Write();
	time_isomer_2_gate_all_134Te->Write();

	time_isomer_2_iftest_gate_134Te->Write();
	time_isomer_2_iftest_gate_bg_134Te->Write();
	time_isomer_2_iftest_gate_all_134Te->Write();

	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_bg_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();

	time_isomer_iftest_doublegate_134Te->Write();
	time_isomer_iftest_doublegate_bg_134Te->Write();
	time_isomer_iftest_doublegate_all_134Te->Write();

	aval_prompt_134Te->Write();
	aval_delayed_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();

	auto stop_program = high_resolution_clock::now();

	auto duration_program = duration_cast<milliseconds>(stop_program - start_program);
	cout << "Runtime program: " << duration_program.count() << endl;

}

void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1){
	//"""Single gated isomer time spectrum """

	//BG, smaller than energy_isomer, subtract
	if( (energy>=energy_isomer-5.0) && (energy<energy_isomer-2.5) ){
		h1->Fill(time,-1);
	}	
	//True, inside energy gate, add
	else if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		h1->Fill(time);
	}
	//BG, larger than energy_isomer, subtract
	else if( (energy>=energy_isomer+2.5) && (energy<energy_isomer+5.0) ){
		h1->Fill(time,-1);
	}
}

void check_energy_singlegate_bg(double energy_isomer, double energy, double time, TH1D *h1){
	//"""Single gated isomer time spectrum, background """

	//BG, smaller than energy_isomer, add
	if( (energy>=energy_isomer-5.0) && (energy<energy_isomer-2.5) ){
		h1->Fill(time);
	}	
	//BG, larger than energy_isomer, add
	else if( (energy>=energy_isomer+2.5) && (energy<energy_isomer+5.0) ){
		h1->Fill(time);
	}
}

void check_energy_singlegate_all(double energy_isomer, double energy, double time, TH1D *h1){
	//"""Single gated isomer time spectrum, before bg-sub """

	//True, inside energy gate, add
	if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		h1->Fill(time);
	}

}


void check_energy_doublegate_true(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1){
	//"""Create doublegated time spectrum, bg-subtracted"""

	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////
	
	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1);
	}
	//BG 1: energy_1 in bg below isomer_A, energy_2 in isomer_B 
	else if( (energy_1>=energy_isomer_A-5.0) && (energy_1<energy_isomer_A-2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 2: energy_1 in bg above isomer_A, energy_2 in isomer_B 
	else if( (energy_1>=energy_isomer_A+2.5) && (energy_1<energy_isomer_A+5.0) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 3: energy_1 in isomer_A and energy_2 in bg below isomer_B
	else if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-5.0) && (energy_2<energy_isomer_B-2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 4: energy_1 in isomer_A and energy_2 in bg above isomer_B
	else if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B+2.5) && (energy_2<energy_isomer_B+5.0)){
		h1->Fill(time_1, -0.25);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_B and energy_2 in isomer_A
	else if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1);
	}
	//BG 1: energy_1 in bg below isomer_B, energy_2 in isomer_A 
	else if( (energy_1>=energy_isomer_B-5.0) && (energy_1<energy_isomer_B-2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 2: energy_1 in bg above isomer_B, energy_2 in isomer_A
	else if( (energy_1>=energy_isomer_B+2.5) && (energy_1<energy_isomer_B+5.0) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 3: energy_1 in isomer_B and energy_2 in bg below isomer_A
	else if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-5.0) && (energy_2<energy_isomer_A-2.5)){
		h1->Fill(time_1, -0.25);
	}
	//BG 4: energy_1 in isomer_B and energy_2 in bg above isomer_A
	else if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A+2.5) && (energy_2<energy_isomer_A+5.0)){
		h1->Fill(time_1, -0.25);
	}
}

void check_energy_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1){
	//"""Doublegated time spectrum BG subtraction"""

	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////

	//BG 1: energy_1 in bg below isomer_A, energy_2 in isomer_B 
	if( (energy_1>=energy_isomer_A-5.0) && (energy_1<energy_isomer_A-2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 2: energy_1 in bg above isomer_A, energy_2 in isomer_B 
	else if( (energy_1>=energy_isomer_A+2.5) && (energy_1<energy_isomer_A+5.0) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 3: energy_1 in isomer_A and energy_2 in bg below isomer_B
	else if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-5.0) && (energy_2<energy_isomer_B-2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 4: energy_1 in isomer_A and energy_2 in bg above isomer_B
	else if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B+2.5) && (energy_2<energy_isomer_B+5.0)){
		h1->Fill(time_1, 0.25);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//BG 1: energy_1 in bg below isomer_B, energy_2 in isomer_A 
	else if( (energy_1>=energy_isomer_B-5.0) && (energy_1<energy_isomer_B-2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 2: energy_1 in bg above isomer_B, energy_2 in isomer_A
	else if( (energy_1>=energy_isomer_B+2.5) && (energy_1<energy_isomer_B+5.0) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 3: energy_1 in isomer_B and energy_2 in bg below isomer_A
	else if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-5.0) && (energy_2<energy_isomer_A-2.5)){
		h1->Fill(time_1, 0.25);
	}
	//BG 4: energy_1 in isomer_B and energy_2 in bg above isomer_A
	else if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A+2.5) && (energy_2<energy_isomer_A+5.0)){
		h1->Fill(time_1, 0.25);
	}
}


void check_energy_doublegate_all(int energy_isomer_A, int energy_isomer_B, int energy_1, int energy_2, int time_1, int time_2, TH1D *h1){
	//"""Create doublegated time spectrum, before bg-sub"""

	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////
	
	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2) && (energy_1<=energy_isomer_A+2) && (energy_2>=energy_isomer_B-2) && (energy_2<=energy_isomer_B+2)){
		h1->Fill(time_1);
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////

	//True: energy_1 in isomer_B and energy_2 in isomer_A
	else if( (energy_1>=energy_isomer_B-2) && (energy_1<=energy_isomer_B+2) && (energy_2>=energy_isomer_A-2) && (energy_2<=energy_isomer_A+2)){
		h1->Fill(time_1);
	}
}





//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`

