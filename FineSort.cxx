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
	lookup_134Te_isomer_1[isomer_energy_1_134Te-6] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-7] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+3] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+4] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+5] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+6] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+7] = 1;



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
	lookup_134Te_isomer_2[isomer_energy_2_134Te-6] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-7] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+3] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+4] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+5] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+6] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+7] = 1;

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
	
	//int N_dummy = 40;
	// int dummy_data[40] = {3, 150, 1000, 1100, 1000, 297, 1000, 
	// 	2, 297, 1000, 1279, 1000, 
	// 	5, 600, 1000, 400, 1000, 297, 1000, 100, 1000, 1279, 1000, 
	// 	2, 297, 1000, 1279, 1000,
	// 	2, 294, 1000, 1279, 1000,
	// 	3, 293, 1000, 1276, 1000, 450, 1000};

	int i_count = 0;
	while(i_count<tpointer_array[0]){
	//while(i_count<N_dummy){

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

			// energy[k] = dummy_data[i_count++];
			// time[k] = dummy_data[i_count++];

			single_gamma->Fill(energy[k]);

			//134Te, lookup isomer_1
			if(lookup_134Te_isomer_1[energy[k]]==2){
				time_isomer_1_gate_134Te->Fill(time[k]);
				time_isomer_1_gate_all_134Te->Fill(time[k]);
			}
			else if(lookup_134Te_isomer_1[energy[k]]==1){
				time_isomer_1_gate_134Te->Fill(time[k],-0.5);
				time_isomer_1_gate_bg_134Te->Fill(time[k],0.5);
			}

			//134Te, lookup isomer_2
			if(lookup_134Te_isomer_2[energy[k]]==2){
				time_isomer_2_gate_134Te->Fill(time[k]);
				time_isomer_2_gate_all_134Te->Fill(time[k]);
			}
			else if(lookup_134Te_isomer_2[energy[k]]==1){
				time_isomer_2_gate_134Te->Fill(time[k],-0.5);
				time_isomer_2_gate_bg_134Te->Fill(time[k],0.5);
			}
		}

		//Lookup doublegate
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
							time_isomer_doublegate_134Te->Fill(time[m], -0.125);
							time_isomer_doublegate_bg_134Te->Fill(time[m], 0.125);
						}

						else if(lookup_134Te_isomer_1[energy[m]]==1 && lookup_134Te_isomer_2[energy[n]]==2){
							time_isomer_doublegate_134Te->Fill(time[m], -0.125);
							time_isomer_doublegate_bg_134Te->Fill(time[m], 0.125);
						}
					}
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
	time_isomer_1_gate_bg_new_134Te->Write();
	time_isomer_1_gate_all_134Te->Write(); 

	time_isomer_2_gate_134Te->Write();
	time_isomer_2_gate_bg_134Te->Write();
	time_isomer_2_gate_bg_new_134Te->Write();
	time_isomer_2_gate_all_134Te->Write();

	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_bg_134Te->Write();
	time_isomer_doublegate_bg_new_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();


	aval_prompt_134Te->Write();
	aval_delayed_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();

	auto stop_program = high_resolution_clock::now();

	auto duration_program = duration_cast<milliseconds>(stop_program - start_program);
	cout << "Runtime program: " << duration_program.count() << endl;

}



//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`

