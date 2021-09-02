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
UShort_t *TheEvents = new UShort_t[bsize];
UInt_t *tpointer_array = new UInt_t[2];

void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_bg(double energy_isomer, double energy, double time, TH1D *h1);
void check_energy_singlegate_all(double energy_isomer, double energy, double time, TH1D *h1);

void aval_prompt_delayed(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time, double aval, TH1D *h1, TH1D *h2);

void check_energy_doublegate_true(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);
void check_energy_doublegate_all(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time_1, double time_2, TH1D *h1);

void coincident_gammas_doublegate(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1);
void coincident_gammas_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1);
void coincident_gammas_doublegate_all(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1);

int main(int argc, char **argv){

	//////////////////////////////////////////
	/// 		   Import files			   ///
	//////////////////////////////////////////

	#include "SpecDefs_FineSort.cxx"


	std::cout << "Reading data from file " << std::endl;
	//Read TheEvents from file
	ifstream infile_edata("/Volumes/240Pu_d_p/nuball/252Cf/Sorted/edata_31may2021.txt", ios::in | ios::binary);
	infile_edata.read((char *) TheEvents, bsize*sizeof(UShort_t));
	infile_edata.close();

	//Read tpointer_array from file
	ifstream tpointer_infile("/Volumes/240Pu_d_p/nuball/252Cf/Sorted/tpointer_31may2021.txt", ios::out | ios::binary);
	tpointer_infile.read((char *) tpointer_array, 2*sizeof(UInt_t));
	tpointer_infile.close();
	std::cout << "tpointer: " << tpointer_array[0] << std::endl;

	std::cout << "Sorting data " << std::endl;


	//////////////////////////////////////////
	/// 		  Isomer properties	       ///
	//////////////////////////////////////////

	//134Te
	double isomer_energy_1_134Te = 297.0;
	double isomer_energy_2_134Te = 1279.0;

	//////////////////////////////////////////
	/// 		  Create spectra	       ///
	//////////////////////////////////////////

	int aval, mult;
	double tot_mult = 0;

	int i_count = 0;
	while(i_count<tpointer_array[0]){

		aval=TheEvents[i_count++];
		mult=TheEvents[i_count++];

		tot_mult += mult;

		double energy[mult];
		double time[mult];

		//Loop over single gammas
		for (int k=0; k < mult; k++){
			energy[k] = TheEvents[i_count++];
			time[k] = TheEvents[i_count++];

			single_gamma->Fill(energy[k]);

			//134Te
			check_energy_singlegate_true(isomer_energy_1_134Te, energy[k], time[k], time_isomer_gate_134Te);
			check_energy_singlegate_bg(isomer_energy_1_134Te, energy[k], time[k], time_isomer_gate_bg_134Te);
			check_energy_singlegate_all(isomer_energy_1_134Te, energy[k], time[k], time_isomer_gate_all_134Te);

			check_energy_singlegate_true(isomer_energy_2_134Te, energy[k], time[k], time_isomer_gate_2_134Te);
		}

		//Loop over all pairs of gammas
		if(mult>1){
			for(int m=0; m < mult-1; m++){
				for(int n=m+1; n < mult; n++ ){
					double_gamma->Fill(energy[m], energy[n]);

					check_energy_doublegate_true(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_doublegate_134Te);
					check_energy_doublegate_bg(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_doublegate_bg_134Te);
					check_energy_doublegate_all(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], time[n], time_isomer_doublegate_all_134Te);
				
					aval_prompt_delayed(isomer_energy_1_134Te, isomer_energy_2_134Te, energy[m], energy[n], time[m], aval, aval_prompt_134Te, aval_delayed_134Te);
				}
			}
		}

		//coincident_gammas_doublegate(isomer_energy_1_134Te, isomer_energy_2_134Te, energy, mult, coincident_gammas_doublegated_134Te);
		//coincident_gammas_doublegate_bg(isomer_energy_1_134Te, isomer_energy_2_134Te, energy, mult, coincident_gammas_doublegated_bg_134Te);
		coincident_gammas_doublegate_all(isomer_energy_1_134Te, isomer_energy_2_134Te, energy, mult, coincident_gammas_doublegated_all_134Te);
	}


	std::cout << "tot mult after writing: " << tot_mult << std::endl;
	std::cout << "\n" << std::endl;


	//////////////////////////////////////////
	/// 	      Write spectra	           ///
	//////////////////////////////////////////

	string OutputSpectraFile="FineSort.root";
	TFile *outputspectrafile = new TFile(OutputSpectraFile.c_str(),"RECREATE");	
	
	single_gamma->Write();
	double_gamma->Write();

	time_isomer_gate_134Te->Write();
	time_isomer_gate_2_134Te->Write();
	time_isomer_gate_bg_134Te->Write();
	time_isomer_gate_all_134Te->Write(); 

	time_isomer_doublegate_134Te->Write();
	time_isomer_doublegate_bg_134Te->Write();
	time_isomer_doublegate_all_134Te->Write();

	coincident_gammas_doublegated_134Te->Write();
	coincident_gammas_doublegated_bg_134Te->Write();
	coincident_gammas_doublegated_all_134Te->Write();

	aval_prompt_134Te->Write();
	aval_delayed_134Te->Write();

	outputspectrafile->cd();
	outputspectrafile->Close();

}

/////////////////////////////////////////////////////////////////////////////////////////////
///									Function definitions 								  ///
/////////////////////////////////////////////////////////////////////////////////////////////


void check_energy_singlegate_true(double energy_isomer, double energy, double time, TH1D *h1){
	//"""Single gated isomer time spectrum """

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
	//"""Single gated isomer time spectrum, background """

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
	//"""Single gated isomer time spectrum, before bg-sub """

	//True, inside energy gate, add
	if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		h1->Fill(time);
	}

}

void aval_prompt_delayed(double energy_isomer_A, double energy_isomer_B, double energy_1, double energy_2, double time, double aval, TH1D *h1, TH1D *h2){

	//////////////////////////////////////////////////////////////
	//Check first if isomer_A: energy_1 and isomer_B: energy_2
	//////////////////////////////////////////////////////////////
	//True: energy_1 in isomer_A and energy_2 in isomer_B
	if( (energy_1>=energy_isomer_A-2.5) && (energy_1<energy_isomer_A+2.5) && (energy_2>=energy_isomer_B-2.5) && (energy_2<energy_isomer_B+2.5)){
		//Take bigger peaks, 500ns 950-1040 is prompt
		if(time>=970 && time<1000){
			h1->Fill(aval);
		}
		//Delayed
		if(time>=1020 && time<1050){
			h2->Fill(aval);
		}
	}

	//////////////////////////////////////////////////////////////
	//Then check opposite: isomer_A: energy_2 and isomer_B: energy_1
	//////////////////////////////////////////////////////////////
	//True: energy_1 in isomer_B and energy_2 in isomer_A
	if( (energy_1>=energy_isomer_B-2.5) && (energy_1<energy_isomer_B+2.5) && (energy_2>=energy_isomer_A-2.5) && (energy_2<energy_isomer_A+2.5)){
		if(time>=970 && time<1000){
			h1->Fill(aval);
		}
		//Delayed
		if(time>=1020 && time<1050){
			h2->Fill(aval);
		}
	}
/*	if( (energy>=energy_isomer-2.5) && (energy<energy_isomer+2.5) ){
		//Prompt
		if(time>=970 && time<1005){
			h1->Fill(aval);
		}
		//Delayed
		if(time>=1005 && time<1040){
			h2->Fill(aval);
		}
	}*/
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
	//"""Doublegated time spectrum BG subtraction"""

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
	//"""Create doublegated time spectrum, before bg-sub"""

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

void coincident_gammas_doublegate(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1){
	//"""Doublegated energy spectrum, true """

	bool isomer_A = false;
	bool isomer_B = false;

	bool isomer_A_lowbg = false;
	bool isomer_B_lowbg = false;

	bool isomer_A_highbg = false;
	bool isomer_B_highbg = false;

	for(int i=0; i<mult;i++){

		//If find isomer_A
		if((energy[i]>=energy_isomer_A-2.5) && (energy[i]<energy_isomer_A+2.5)){
			isomer_A = true;
		}
		//If find isomer_B
		if((energy[i]>=energy_isomer_B-2.5) && (energy[i]<energy_isomer_B+2.5)){
			isomer_B = true;
		}

		//Low bg isomer_A
		if((energy[i]>=energy_isomer_A-5.0) && (energy[i]<energy_isomer_A-2.5)){
			isomer_A_lowbg = true;
		}

		//Low bg isomer_B
		if((energy[i]>=energy_isomer_B-5.0) && (energy[i]<energy_isomer_B-2.5)){
			isomer_B_lowbg = true;
		}

		//High bg isomer_A
		if((energy[i]>=energy_isomer_A+2.5) && (energy[i]<energy_isomer_A+5.0)){
			isomer_A_highbg = true;
		}

		//High bg isomer_B
		if((energy[i]>=energy_isomer_B+2.5) && (energy[i]<energy_isomer_B+5.0)){
			isomer_B_highbg = true;
		}
	}

	int j;
	if(isomer_A==true && isomer_B==true){
		for(int j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
	if(isomer_A==true && isomer_B_lowbg==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j], -0.25);
		}
	}
	if(isomer_A==true && isomer_B_highbg==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j], -0.25);
		}
	}
	if(isomer_A==isomer_A_lowbg && isomer_B==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j], -0.25);
		}
	}
	if(isomer_A_highbg==true && isomer_B==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j], -0.25);
		}
	}
}


void coincident_gammas_doublegate_bg(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1){
	//"""Doublegated energy spectrum, the subtracted bg """

	bool isomer_A = false;
	bool isomer_B = false;

	bool isomer_A_lowbg = false;
	bool isomer_B_lowbg = false;

	bool isomer_A_highbg = false;
	bool isomer_B_highbg = false;

	for(int i=0; i<mult;i++){

		//If find isomer_A
		if((energy[i]>=energy_isomer_A-2.5) && (energy[i]<energy_isomer_A+2.5)){
			isomer_A = true;
		}
		//If find isomer_B
		if((energy[i]>=energy_isomer_B-2.5) && (energy[i]<energy_isomer_B+2.5)){
			isomer_B = true;
		}

		//Low bg isomer_A
		if((energy[i]>=energy_isomer_A-5.0) && (energy[i]<energy_isomer_A-2.5)){
			isomer_A_lowbg = true;
		}

		//Low bg isomer_B
		if((energy[i]>=energy_isomer_B-5.0) && (energy[i]<energy_isomer_B-2.5)){
			isomer_B_lowbg = true;
		}

		//High bg isomer_A
		if((energy[i]>=energy_isomer_A+2.5) && (energy[i]<energy_isomer_A+5.0)){
			isomer_A_highbg = true;
		}

		//High bg isomer_B
		if((energy[i]>=energy_isomer_B+2.5) && (energy[i]<energy_isomer_B+5.0)){
			isomer_B_highbg = true;
		}

	}

	int j;
	if(isomer_A==true && isomer_B_lowbg==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
	if(isomer_A==true && isomer_B_highbg==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
	if(isomer_A==isomer_A_lowbg && isomer_B==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
	if(isomer_A_highbg==true && isomer_B==true){
		for(j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
}


void coincident_gammas_doublegate_all(double energy_isomer_A, double energy_isomer_B, double energy[], int mult, TH1D *h1){
	//"""Doublegated energy spectrum, before bg-sub """

	bool isomer_A = false;
	bool isomer_B = false;

	bool isomer_A_lowbg = false;
	bool isomer_B_lowbg = false;

	bool isomer_A_highbg = false;
	bool isomer_B_highbg = false;

	for(int i=0; i<mult;i++){

		//If find isomer_A
		if((energy[i]>=energy_isomer_A-2.5) && (energy[i]<energy_isomer_A+2.5)){
			isomer_A = true;
		}
		//If find isomer_B
		if((energy[i]>=energy_isomer_B-2.5) && (energy[i]<energy_isomer_B+2.5)){
			isomer_B = true;
		}
	}

	if(isomer_A==true && isomer_B==true){
		for(int j=0; j<mult; j++){
			h1->Fill(energy[j]);
		}
	}
}

//g++ -g -o FineSort FineSort.cxx ` root-config --cflags` `root-config --glibs`

