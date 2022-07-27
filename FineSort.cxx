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

void fill_lookuptable(int energy, int lookup[65535]);
void fill_spectra_doublegate(int lookup_value_1, int lookup_value_2, int time, TH1D *true_spec, TH1D *all_spec, TH1D *bg_spec, TH1D *bg_ridge_spec, TH1D *bg_random_spec);

int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_FineSort.cxx"

	auto start_program = high_resolution_clock::now();

	// //Read TheEvents from file
	string OutputDirectory="/Applications/nuball_sorting/IYR_Data/";
	std::cout << "Reading data from file " << std::endl;
	//ifstream infile_edata("/Applications/nuball_sorting/IYR_Data/edata_allfiles_10sep2021.txt", ios::in | ios::binary);
	ifstream infile_edata("/Applications/nuball_sorting/IYR_Data/edata_11jan2022.txt", ios::in | ios::binary);
	infile_edata.read((char *) TheEvents, bsize*sizeof(UShort_t));
	infile_edata.close();

	//Read tpointer_array from file
	//ifstream tpointer_infile("/Applications/nuball_sorting/IYR_Data/tpointer_allfiles_10sep2021.txt", ios::out | ios::binary);
	ifstream tpointer_infile("/Applications/nuball_sorting/IYR_Data/tpointer_11jan2022.txt", ios::out | ios::binary);
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

	//140Xe
	int gamma_energy_1_140Xe = 458;
	int gamma_energy_2_140Xe = 377;

	//132Sn
	int gamma_energy_132Sn = 4042;

	//////////////////////////////////////////
	/// 		   Lookup tables		   ///
	//////////////////////////////////////////

	int bg_param = 1;

	//134Te, isomer_1
	int lookup_134Te_isomer_1[65535] = {0};
	fill_lookuptable(isomer_energy_1_134Te, lookup_134Te_isomer_1);

	for(int i=isomer_energy_1_134Te-10; i<=isomer_energy_1_134Te+10; i++){
		cout << lookup_134Te_isomer_1[i] << endl;
	}

	//134Te, isomer_2
	int lookup_134Te_isomer_2[65535] = {0};
	fill_lookuptable(isomer_energy_2_134Te, lookup_134Te_isomer_2);

	//140Xe, gamma_1
	int lookup_140Xe_gamma_1[65535] = {0};
	fill_lookuptable(gamma_energy_1_140Xe, lookup_140Xe_gamma_1);

	//140Xe, gamma_2
	int lookup_140Xe_gamma_2[65535] = {0};
	fill_lookuptable(gamma_energy_2_140Xe, lookup_140Xe_gamma_2);

	//132Sn
	//Needs to be fixed the width of the gates
	int lookup_132Sn_gamma[65535] = {0};
	fill_lookuptable(gamma_energy_132Sn, lookup_132Sn_gamma);


	//////////////////////////////////////////
	/// 		  	 Sorting	     	   ///
	//////////////////////////////////////////

	double tot_mult = 0;
	UShort_t hitmult;
	unsigned char hit, mult;

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

		//Decompress hit and mult
		hitmult = TheEvents[i_count++];
		hit = hitmult & 0xFF;
		mult = hitmult >> 8;
		tot_mult += mult;

		//mult = dummy_data[i_count++];

		int energy[mult];
		int time[mult];

		mult_distr->Fill(mult);
		hit_distr->Fill(hit);

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
				time_isomer_1_gate_134Te->Fill(time[k],-bg_param);
				time_isomer_1_gate_bg_134Te->Fill(time[k],bg_param);
			}

			//134Te, lookup isomer_2
			if(lookup_134Te_isomer_2[energy[k]]==2){
				time_isomer_2_gate_134Te->Fill(time[k]);
				time_isomer_2_gate_all_134Te->Fill(time[k]);
			}
			else if(lookup_134Te_isomer_2[energy[k]]==1){
				time_isomer_2_gate_134Te->Fill(time[k],-bg_param);
				time_isomer_2_gate_bg_134Te->Fill(time[k],bg_param);
			}
		}

		//Lookup doublegate
		//Loop over all pairs of gammas
		int dt;
		if(mult>1){
			for(int m=0; m < mult; m++){
				for(int n=0; n < mult; n++){
					if(m!=n){
						double_gamma->Fill(energy[m], energy[n]);

						int dt = abs(time[m]-time[n]);

						if(lookup_134Te_isomer_1[energy[m]]==2 && lookup_134Te_isomer_2[energy[n]]==2){
							hit_doublegate->Fill(hit);
						}

						///////////////
						//   134Te   //
						///////////////

						//Fill with time_m
						fill_spectra_doublegate(lookup_134Te_isomer_1[energy[m]], lookup_134Te_isomer_2[energy[n]], time[m], time_isomer_doublegate_1_134Te, time_isomer_doublegate_1_all_134Te, time_isomer_doublegate_1_bg_134Te, time_isomer_doublegate_1_bg_ridge_134Te, time_isomer_doublegate_1_bg_random_134Te);
						//Fill with time_n
						fill_spectra_doublegate(lookup_134Te_isomer_1[energy[m]], lookup_134Te_isomer_2[energy[n]], time[n], time_isomer_doublegate_2_134Te, time_isomer_doublegate_2_all_134Te, time_isomer_doublegate_2_bg_134Te, time_isomer_doublegate_2_bg_ridge_134Te, time_isomer_doublegate_2_bg_random_134Te);
						
						if(hit>=3){
							fill_spectra_doublegate(lookup_134Te_isomer_1[energy[m]], lookup_134Te_isomer_2[energy[n]], time[n], time_isomer_doublegate_2_hit3_134Te, time_isomer_doublegate_2_all_hit3_134Te, time_isomer_doublegate_2_bg_hit3_134Te, time_isomer_doublegate_2_bg_ridge_hit3_134Te, time_isomer_doublegate_2_bg_random_hit3_134Te);	
						}

						if(hit>=4){
							fill_spectra_doublegate(lookup_134Te_isomer_1[energy[m]], lookup_134Te_isomer_2[energy[n]], time[n], time_isomer_doublegate_2_hit4_134Te, time_isomer_doublegate_2_all_hit4_134Te, time_isomer_doublegate_2_bg_hit4_134Te, time_isomer_doublegate_2_bg_ridge_hit4_134Te, time_isomer_doublegate_2_bg_random_hit4_134Te);	
						}

						///////////////
						//   140Xe   //
						///////////////

						fill_spectra_doublegate(lookup_140Xe_gamma_1[energy[m]], lookup_140Xe_gamma_2[energy[n]], time[m], time_gamma_doublegate_1_140Xe, time_gamma_doublegate_1_all_140Xe, time_gamma_doublegate_1_bg_140Xe, time_gamma_doublegate_1_bg_ridge_140Xe, time_gamma_doublegate_1_bg_random_140Xe);

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
	time_isomer_1_gate_all_134Te->Write(); 

	time_isomer_2_gate_134Te->Write();
	time_isomer_2_gate_bg_134Te->Write();
	time_isomer_2_gate_all_134Te->Write();

	time_isomer_doublegate_1_134Te->Write();
	time_isomer_doublegate_1_all_134Te->Write();
	time_isomer_doublegate_1_bg_134Te->Write();
	time_isomer_doublegate_1_bg_ridge_134Te->Write();
	time_isomer_doublegate_1_bg_random_134Te->Write();

	time_isomer_doublegate_1_dt70_134Te->Write();
	time_isomer_doublegate_1_all_dt70_134Te->Write();
	time_isomer_doublegate_1_bg_dt70_134Te->Write();

	time_isomer_doublegate_2_134Te->Write();
	time_isomer_doublegate_2_all_134Te->Write();
	time_isomer_doublegate_2_bg_134Te->Write();
	time_isomer_doublegate_2_bg_ridge_134Te->Write();
	time_isomer_doublegate_2_bg_random_134Te->Write();

	time_isomer_doublegate_2_dt70_134Te->Write();
	time_isomer_doublegate_2_all_dt70_134Te->Write();
	time_isomer_doublegate_2_bg_dt70_134Te->Write();

	time_isomer_doublegate_2_hit3_134Te->Write();
	time_isomer_doublegate_2_all_hit3_134Te->Write();
	time_isomer_doublegate_2_bg_hit3_134Te->Write();
	time_isomer_doublegate_2_bg_ridge_hit3_134Te->Write();
	time_isomer_doublegate_2_bg_random_hit3_134Te->Write();

	time_isomer_doublegate_2_hit4_134Te->Write();
	time_isomer_doublegate_2_all_hit4_134Te->Write();
	time_isomer_doublegate_2_bg_hit4_134Te->Write();
	time_isomer_doublegate_2_bg_ridge_hit4_134Te->Write();
	time_isomer_doublegate_2_bg_random_hit4_134Te->Write();

	time_gamma_doublegate_1_140Xe->Write();
	time_gamma_doublegate_1_all_140Xe->Write();
	time_gamma_doublegate_1_bg_140Xe->Write();
	time_gamma_doublegate_1_bg_ridge_140Xe->Write();
	time_gamma_doublegate_1_bg_random_140Xe->Write();

	time_gamma_gate_132Sn->Write();
	time_gamma_gate_all_132Sn->Write();
	time_gamma_gate_bg_132Sn->Write();

	aval_prompt_134Te->Write();
	aval_delayed_134Te->Write();

	mult_distr->Write();
	hit_distr->Write();
	hit_doublegate->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();

	auto stop_program = high_resolution_clock::now();

	auto duration_program = duration_cast<milliseconds>(stop_program - start_program);
	cout << "Runtime program: " << duration_program.count() << endl;

}




void fill_lookuptable(int energy, int lookup[65535]){
	///1 for bg, 2 for peak

	//fill whole area with 1 first
	for(int i=energy-7; i<=energy+7; i++){
		lookup[i] = 1;
	}

	//peak gate; e
	for(int i=energy-2; i<=energy+2; i++){
		lookup[i] = 2;
	}
}





void fill_spectra_doublegate(int lookup_value_1, int lookup_value_2, int time, TH1D *true_spec, TH1D *all_spec, TH1D *bg_spec, TH1D *bg_ridge_spec, TH1D *bg_random_spec){
	//a b c
	//d e f
	//g h i

	//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
	//bg_ridge = (b+d+f+h)/2
	//bg_random = (a+c+g+i)/4

	//NB: WITH ASYMMETRIC BG GATES, THE BG SUBTRACTION IS WRONG!

	//e
	if(lookup_value_1==2 && lookup_value_2==2){
		//Increment spec with time of one gamma-ray
		true_spec->Fill(time);
		all_spec->Fill(time);
	}

	//h & b
	else if(lookup_value_1==2 && lookup_value_2==1){
		true_spec->Fill(time, -0.5);
		bg_spec->Fill(time, 0.5);
		bg_ridge_spec->Fill(time, 0.5);
	}

	//d & f 
	else if(lookup_value_1==1 && lookup_value_2==2){
		true_spec->Fill(time, -0.5);
		bg_spec->Fill(time, 0.5);
		bg_ridge_spec->Fill(time, 0.5);
	}

	//a & c & g & i
	else if(lookup_value_1==1 && lookup_value_2==1){
		true_spec->Fill(time, 0.25);
		bg_spec->Fill(time, -0.25);
		bg_random_spec->Fill(time, 0.25);	
	}

}


//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`
