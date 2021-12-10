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
	//peak gate
	lookup_134Te_isomer_1[isomer_energy_1_134Te] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-1] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-2] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+1] = 2;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+2] = 2;
	//bg gate
	lookup_134Te_isomer_1[isomer_energy_1_134Te-3] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-4] = 1;
/*	lookup_134Te_isomer_1[isomer_energy_1_134Te-5] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-6] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te-7] = 1;*/
	lookup_134Te_isomer_1[isomer_energy_1_134Te+3] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+4] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+5] = 1;
/*	lookup_134Te_isomer_1[isomer_energy_1_134Te+6] = 1;
	lookup_134Te_isomer_1[isomer_energy_1_134Te+7] = 1;*/

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
/*	lookup_134Te_isomer_2[isomer_energy_2_134Te-5] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-6] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te-7] = 1;*/
	lookup_134Te_isomer_2[isomer_energy_2_134Te+3] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+4] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+5] = 1;
/*	lookup_134Te_isomer_2[isomer_energy_2_134Te+6] = 1;
	lookup_134Te_isomer_2[isomer_energy_2_134Te+7] = 1;*/

	//140Xe, gamma_1
	int lookup_140Xe_gamma_1[65535] = {0};
	//peak gate
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe] = 2;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-1] = 2;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-2] = 2;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+1] = 2;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+2] = 2;
	//bg gate
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-3] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-4] = 1;
/*	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-5] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-6] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe-7] = 1;*/
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+3] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+4] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+5] = 1;
/*	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+6] = 1;
	lookup_140Xe_gamma_1[gamma_energy_1_140Xe+7] = 1;*/

	//140Xe, gamma_2
	int lookup_140Xe_gamma_2[65535] = {0};
	//peak gate
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe] = 2;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-1] = 2;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-2] = 2;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+1] = 2;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+2] = 2;
	//bg gate
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-3] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-4] = 1;
	/*lookup_140Xe_gamma_2[gamma_energy_2_140Xe-5] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-6] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe-7] = 1;*/
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+3] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+4] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+5] = 1;
/*	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+6] = 1;
	lookup_140Xe_gamma_2[gamma_energy_2_140Xe+7] = 1;*/

	//132Sn
	int lookup_132Sn_gamma[65535] = {0};
	//peak gate
	lookup_132Sn_gamma[gamma_energy_132Sn-7] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-6] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-5] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-4] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-3] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-1] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn-2] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+1] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+2] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+3] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+4] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+5] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+6] = 2;
	lookup_132Sn_gamma[gamma_energy_132Sn+7] = 2;
	//bg gate
	lookup_132Sn_gamma[gamma_energy_132Sn-8] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-9] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-10] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-11] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-12] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-13] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn-14] = 1;

	lookup_132Sn_gamma[gamma_energy_132Sn+8] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+9] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+10] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+11] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+12] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+13] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+14] = 1;
	lookup_132Sn_gamma[gamma_energy_132Sn+15] = 1;


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

		mult_distr->Fill(mult);

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

			//132Sn
			if(lookup_132Sn_gamma[energy[k]]==2){
				time_gamma_gate_132Sn->Fill(time[k]);
				time_gamma_gate_all_132Sn->Fill(time[k]);
			}
			else if(lookup_132Sn_gamma[energy[k]]==1){
				time_gamma_gate_132Sn->Fill(time[k],-bg_param);
				time_gamma_gate_bg_132Sn->Fill(time[k],bg_param);
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
						/////////////////////
						// 134Te doublegate
						/////////////////////

						//a b c
						//d e f
						//g h i

						//true = e - (b+d+f+h)/2 + (a+c+g+i)/4
						//bg_ridge = (b+d+f+h)/2
						//bg_random = (a+c+g+i)/4

						//e
						if(lookup_134Te_isomer_1[energy[m]]==2 && lookup_134Te_isomer_2[energy[n]]==2){
							//Increment time[m] since m is 297keV line
							time_isomer_doublegate_1_134Te->Fill(time[m]);
							time_isomer_doublegate_1_all_134Te->Fill(time[m]);

							//Increment time[n] since n is 1279keV line
							time_isomer_doublegate_2_134Te->Fill(time[n]);
							time_isomer_doublegate_2_all_134Te->Fill(time[n]);

							if(dt<=70){

								time_isomer_doublegate_1_dt70_134Te->Fill(time[m]);
								time_isomer_doublegate_1_all_dt70_134Te->Fill(time[m]);

								time_isomer_doublegate_2_dt70_134Te->Fill(time[n]);
								time_isomer_doublegate_2_all_dt70_134Te->Fill(time[n]);
							}
						}
						//b & h, so bg_ridge
						else if(lookup_134Te_isomer_1[energy[m]]==2 && lookup_134Te_isomer_2[energy[n]]==1){
							time_isomer_doublegate_1_134Te->Fill(time[m], -bg_param*0.5);
							time_isomer_doublegate_1_bg_134Te->Fill(time[m], bg_param*0.5);
							time_isomer_doublegate_1_bg_ridge_134Te->Fill(time[m], bg_param);

							time_isomer_doublegate_2_134Te->Fill(time[n], -bg_param*0.5);
							time_isomer_doublegate_2_bg_134Te->Fill(time[n], bg_param*0.5);
							time_isomer_doublegate_2_bg_ridge_134Te->Fill(time[n], bg_param);

							if(dt<=70){
								time_isomer_doublegate_1_dt70_134Te->Fill(time[m], -bg_param*0.5);
								time_isomer_doublegate_1_bg_dt70_134Te->Fill(time[m], bg_param*0.5);

								time_isomer_doublegate_2_dt70_134Te->Fill(time[n], -bg_param*0.5);
								time_isomer_doublegate_2_bg_dt70_134Te->Fill(time[n], bg_param*0.5);
							}
						}

						// d & f, so bg_ridge
						else if(lookup_134Te_isomer_1[energy[m]]==1 && lookup_134Te_isomer_2[energy[n]]==2){
							time_isomer_doublegate_1_134Te->Fill(time[m], -bg_param*0.5);
							time_isomer_doublegate_1_bg_134Te->Fill(time[m], bg_param*0.5);
							time_isomer_doublegate_1_bg_ridge_134Te->Fill(time[m], bg_param);

							time_isomer_doublegate_2_134Te->Fill(time[n], -bg_param*0.5);
							time_isomer_doublegate_2_bg_134Te->Fill(time[n], bg_param*0.5);
							time_isomer_doublegate_2_bg_ridge_134Te->Fill(time[n], bg_param);

							if(dt<=70){
								time_isomer_doublegate_1_dt70_134Te->Fill(time[m], -bg_param*0.5);
								time_isomer_doublegate_1_bg_dt70_134Te->Fill(time[m], bg_param*0.5);

								time_isomer_doublegate_2_dt70_134Te->Fill(time[n], -bg_param*0.5);
								time_isomer_doublegate_2_bg_dt70_134Te->Fill(time[n], bg_param*0.5);
							}
						}

						//a + c + g + i, so bg_random
						else if(lookup_134Te_isomer_1[energy[m]]==1 && lookup_134Te_isomer_2[energy[n]]==1){
							time_isomer_doublegate_1_134Te->Fill(time[m], bg_param*0.25);
							time_isomer_doublegate_1_bg_134Te->Fill(time[m], -bg_param*0.25);
							time_isomer_doublegate_1_bg_random_134Te->Fill(time[m], bg_param);

							time_isomer_doublegate_2_134Te->Fill(time[n], bg_param*0.25);
							time_isomer_doublegate_2_bg_134Te->Fill(time[n], -bg_param*0.25);
							time_isomer_doublegate_2_bg_random_134Te->Fill(time[n], bg_param);

							if(dt<=70){
								time_isomer_doublegate_1_dt70_134Te->Fill(time[m], bg_param*0.25);
								time_isomer_doublegate_1_bg_dt70_134Te->Fill(time[m], -bg_param*0.25);

								time_isomer_doublegate_2_dt70_134Te->Fill(time[n], bg_param*0.25);
								time_isomer_doublegate_2_bg_dt70_134Te->Fill(time[n], -bg_param*0.25);
							}
						}

						/////////////////////
						// 140Xe doublegate
						/////////////////////

						//e
						if(lookup_140Xe_gamma_1[energy[m]]==2 && lookup_140Xe_gamma_2[energy[n]]==2){
							time_gamma_doublegate_140Xe->Fill(time[m]);
							time_gamma_doublegate_all_140Xe->Fill(time[m]);
						}
						//b & h, so bg_ridge
						else if(lookup_140Xe_gamma_1[energy[m]]==2 && lookup_140Xe_gamma_2[energy[n]]==1){
							time_gamma_doublegate_140Xe->Fill(time[m], -bg_param*0.5);
							time_gamma_doublegate_bg_140Xe->Fill(time[m], bg_param*0.5);
							time_gamma_doublegate_bg_ridge_140Xe->Fill(time[m], bg_param);
						}
						// d & f, so bg_ridge
						else if(lookup_140Xe_gamma_1[energy[m]]==1 && lookup_140Xe_gamma_2[energy[n]]==2){
							time_gamma_doublegate_140Xe->Fill(time[m], -bg_param*0.5);
							time_gamma_doublegate_bg_140Xe->Fill(time[m], bg_param*0.5);
							time_gamma_doublegate_bg_ridge_140Xe->Fill(time[m], bg_param);
						}
						//a + c + g + i, so bg_random
						else if(lookup_140Xe_gamma_1[energy[m]]==1 && lookup_140Xe_gamma_2[energy[n]]==1){
							time_gamma_doublegate_140Xe->Fill(time[m], bg_param*0.25);
							time_gamma_doublegate_bg_140Xe->Fill(time[m], -bg_param*0.25);
							time_gamma_doublegate_bg_random_140Xe->Fill(time[m], bg_param);
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

	mult_distr->Write();

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

	time_isomer_doublegate_mult3_134Te->Write();
	time_isomer_doublegate_bg_mult3_134Te->Write();
	time_isomer_doublegate_all_mult3_134Te->Write();

	time_gamma_doublegate_140Xe->Write();
	time_gamma_doublegate_all_140Xe->Write();
	time_gamma_doublegate_bg_140Xe->Write();
	time_gamma_doublegate_bg_ridge_140Xe->Write();
	time_gamma_doublegate_bg_random_140Xe->Write();

	time_gamma_gate_132Sn->Write();
	time_gamma_gate_all_132Sn->Write();
	time_gamma_gate_bg_132Sn->Write();

	aval_prompt_134Te->Write();
	aval_delayed_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();

	auto stop_program = high_resolution_clock::now();

	auto duration_program = duration_cast<milliseconds>(stop_program - start_program);
	cout << "Runtime program: " << duration_program.count() << endl;

}



//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`

