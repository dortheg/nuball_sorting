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
#include<iostream>
#include <iomanip>
#include <string> 

inline bool FileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

using namespace std;

double TimePeaks[256][200];
double APeaks[3][200];
double tshifts[256][200]; //detector number, file number

// THE INPUT FILE
string FileList="runlist_allfiles.dat"; //The list of Event-built files to process
//string OutputDirectory="/Applications/nuball_sorting/252Cf/";
string OutputDirectory="/Users/dortheagjestvang/Documents/Author_articles/IYR/Data/";
string InputDirectory="/Volumes/240Pu_d_p/nuball_data/252Cf/252Cf/";


//Things needed for data storage
unsigned int bsize=4000000000;
UShort_t *TheEvents=new UShort_t[bsize]; //4 billion word buffer size
UShort_t *TheEvents_1 = new UShort_t[bsize];
UInt_t *tpointer_array = new UInt_t[2];
UInt_t *tpointer_array_1 = new UInt_t[2];
UInt_t tpointer=0;
double tot_mult = 0;
double tot_mult_char = 0;
double tot_mult_char_out = 0;


int NCHANS=4096;
int NDETS=255;
int TOFF=1311;
bool CAL_CORRECT=false;

// MAKE LOOKUP TABLES OF ALL THESE
int *ring=new int[NDETS];
int *alveole=new int[NDETS];
int *module=new int[NDETS];
bool *isthere=new bool[NDETS];
bool *isge=new bool[NDETS];
bool *isbgo=new bool[NDETS];
bool *isphase1ge=new bool[NDETS];
bool *isphase1bgo=new bool[NDETS];
bool *isphase1=new bool[NDETS];
bool *isclover=new bool[NDETS];
bool *iscloverge=new bool[NDETS];
bool *iscloverbgo=new bool[NDETS];
bool *islabr3=new bool[NDETS];
bool *ismadrid=new bool[NDETS];
bool *iseden=new bool[NDETS];
bool *ispulse=new bool[NDETS];
bool *bigwalk=new bool[NDETS];
bool *isbad=new bool[NDETS];
bool *isused=new bool[NDETS];
int *isring1=new int[NDETS];
int *isring2=new int[NDETS];
int *isring3=new int[NDETS];
int *isring4=new int[NDETS];

///////////////////////////////////////////////////////
/// 		           Start main 	 		        ///
///////////////////////////////////////////////////////

int main(int argc, char **argv)
{

//THE OUTPUT FILE
//ofstream outfile("/Applications/nuball/252Cf/252Cf_outputfile_IYR.txt");

double alti[256];
double elti[256];
double multi[356];

//CubeDDT *TheCube=new CubeDDT();

int ND=256;
int pcounter=0;

//read in the time shifts
ifstream infile("tshifts_new.dat");
int i,j;
double val;
while(infile >> i >> j >> val)
{
tshifts[i][j]=1000.0-val;
if (val > 0) 
	{
	tshifts[i][j]=1000.0-val;
	}
}
for (int i=0; i < 256; i++)
	{
	for (int j=0; j < 100; j++)
		{
		//cout << i << " " << j << " " << tshifts[i][j] << endl;
		}
}
//exit(1);

///////////////////////////////////////////////////////
/// 		         Calibration	 		        ///
///////////////////////////////////////////////////////

double *Gains=new double[256];
double *Offsets=new double[256];
double *Squared=new double[256];

double *Gains2=new double[256]; //Ge Calibration correction
double *Offsets2=new double[256];
double *Squared2=new double[256];

double *Gains3=new double[256]; //BGO Calibration correction
double *Offsets3=new double[256];
double *Squared3=new double[256];

//Ge calibration
int id,id2;
for (int i=0; i < 256; i++)
{
	id=int(EuPeaks[i][0]); //Calibration is for each faster channel. 
	if ((id > 0) && (EuPeaks[i][1] > 0))
	{
	Gains[id]=(1408.0110-121.7830)/(EuPeaks[i][2]-EuPeaks[i][1]);
	Offsets[id]=121.7830-EuPeaks[i][1]*Gains[id];
	//cout <<id << " " << Gains[id] << " " << Offsets[id] << endl;
	if ((id==61) || (id==144))
		{
		Gains[id]=(1408.0110-344.2785)/(EuPeaks[i][2]-EuPeaks[i][1]);
		Offsets[id]=344.2785-EuPeaks[i][1]*Gains[id];
		}
	if (id==180) 
		{
		Gains[id]=(1279.1-511.0)/(EuPeaks[i][2]-EuPeaks[i][1]);
		Offsets[id]=511.0-EuPeaks[i][1]*Gains[id];
		}
	}
	else {Gains[id]=1.0; Offsets[id]=0;}
	//if (i==180) {Gains[id]=0.1; Offsets[id]=0;}
	
	id2=CalCorrection[i][0];
	if (CalCorrection[i][0] > 0) 
		{
		Offsets2[id2]=CalCorrection[i][1];
		Gains2[id2]=CalCorrection[i][2];
		Squared[id2]=CalCorrection[i][3];
		}
	else {Gains2[id2]=0; Offsets2[id2]=0; Squared[id2]=0;}

}
// Calibrate the LaBr3 from my peak positions
for (int i=0; i <= 20; i++)
{
	int id=int(NewLaBr3Cal[i][0]);
	if (id > 0)
	{
	//Gains[id]=(1408.0110-121.7830)/((EuPeaksL[i][2]-EuPeaksL[i][1])*500.0);
	//Offsets[id]=121.7830-(EuPeaksL[i][1]*500.0)*Gains[id];
	Offsets[id]=NewLaBr3Cal[i][1];
	Gains[id]=NewLaBr3Cal[i][2];
	Squared[id]=NewLaBr3Cal[i][3];
	//cout <<id << " " << Gains[id] << " " << Offsets[id] << endl;
	}	
	//else {Gains[id]=1.0; Offsets[id]=0;}

}

srand48(time(0));

// Calibrate the BGO's from my 1172 peak positions
// Ring 4 are given 50 since were not there for 60Co.
// Ring 4 needs a correction so that mean spectral values correspond to clover BGO mean values.
for (int i=0; i < 60; i++)
{
if (CalBGO[i][0] > 0)
	{
	int theid=CalBGO[i][0];
	double val=CalBGO[i][1];
	Gains[theid]=(1172/val)/500.0; //channel# is divded by 500 to produce spectrum. val normally 50.
	Offsets[theid]=0;
	}
//NEW CORRECTION
if (CorBGO[i][0] > 0)
	{
	int theid=CorBGO[i][0];
	double offset=CorBGO[i][1];
	double gain=CorBGO[i][2];
	Gains3[theid]=gain;
	Offsets3[theid]=offset;
	}
}

///////////////////////////////////////////////////////
/// 		       File handling 	 		        ///
///////////////////////////////////////////////////////

#include "SpecDefs.cxx"
#include "Lookup.cxx"
//Plus two gamma-gamma matrices (prompt and delayed)
TH2F* matb = new TH2F("matb","matb",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr
TH2F* matp = new TH2F("matp","matp",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr
TH2F* mati = new TH2F("mati","mati",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr
TH2F* matl = new TH2F("matl","matl",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr

TH2F* lala = new TH2F("lala","lala",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr
TH2F* gela = new TH2F("gela","gela",4096,0,4096,4096,0,4096); //xchans,0,xr,ychans,0,yr

long EventCount=0;
long GoodEventCount=0;

////////////////////////////////////////////////// FILE HANDLING ////////////

// Make file names. Find run number
ofstream logfile("log.dat");
ifstream iflist2(FileList.c_str());
string line;
int countfiles=0;
while(iflist2 >> line){countfiles++;}
iflist2.close();

cout << "Processing " << countfiles << " files " <<endl;
ifstream iflist(FileList.c_str());

// Read through again and process each file in turn
int nfiles=0;
int RunNumber;
int OldRunNumber;
int FileNumber;
bool LastFile=false, RunChange=false; //Flags


while(iflist >> line)
{
nfiles++;
if (nfiles==countfiles) {cout << "Reading LAST file " << endl; LastFile=true;}
RunNumber=atoi(line.substr(4,3).c_str());
if (nfiles==1) {OldRunNumber=RunNumber;}
if (RunNumber != OldRunNumber) {RunChange=true; OldRunNumber=RunNumber;}

int fn1,fn2,fn3;

//D: not used for 252Cf (check before remove)
fn1=atoi(line.substr(line.size()-6,1).c_str()); //first digit
fn2=atoi(line.substr(line.size()-7,1).c_str()); //second digit
fn3=atoi(line.substr(line.size()-8,1).c_str()); //third digit
if ((line.substr(line.size()-7,1) == "_")) {FileNumber=fn1;}
if ((line.substr(line.size()-8,1) == "_")) {FileNumber=fn1+(10*fn2);}
if ((line.substr(line.size()-9,1) == "_")) {FileNumber=fn1+(10*fn2)+(100*fn3);}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START TO PROCESS THIS FILE !!!!!!!!!!!!!!!!
TStopwatch timer;
timer.Reset();
timer.Start();

string InputFile=InputDirectory+line;
cout << "Reading " << InputFile << endl;
string OutputSpectraFile=OutputDirectory+"spectra.root"; //D: output

//////////////////////////////////////// THE INPUT/OUTPUT TREES ////////////////////////////////
ULong64_t tab_tm;
UChar_t tab_label;
UInt_t tab_nrj;
UInt_t tab_nrj2;
Bool_t tab_pileup;

//DEFINE INPUT TREE
if (!FileExists(InputFile.c_str())) {cout << "File doesn't exit. Skipping to next"<<endl; continue;} // skip this file

TFile *fin = new TFile(InputFile.c_str(),"READ");

//DEFINE OUTPUT TREE
// New tree in which to write the sorted data
UInt_t size = 0; //D: due to mac memory storage, might have to set everything to 0 before sorting
UInt_t energies[512]; 
Double_t timestamps[512];
UInt_t ids[512]; //ids array
Bool_t pups[512]; //pileups array

UInt_t sz;
UInt_t VALUE[300];
Double_t ts[300];

UInt_t mult; // Multiplicity of the event

//D: where data format is described
TTree* inputtree = (TTree*)fin->Get("DataTree");;
inputtree->SetBranchAddress("size",&mult);            // Multiplicity of event, all events
inputtree->SetBranchAddress("ts",timestamps);         // Time
inputtree->SetBranchAddress("ID",ids);                // Detector number
inputtree->SetBranchAddress("VALUE",energies);        // Energy of hit

int entries=inputtree->GetEntries();
cout << "File " << nfiles << " contains " << entries << " events" << endl;
logfile << "File " << nfiles << " " << InputFile << " contains " << entries << " events" << endl;

int TriggerCount=0;

//D: global event param, adding everything that happened in event (group of hits) together
double TotalEnergy=0;
int LaMult=0; // LaBr3 Mult
int DLaMult=0;
int CMult=0; //Prompt clean Ge mult; D: clean means no surrounding BGO within a given time spectrum for a Ge
int TotalMult=0; //Total Event Multiplicity, of gammas, preprompt, prompt, delayed
int PromptMult=0; 
int DelayedMult=0;

//PROMPT: Used for Compton suppression and addback in clovers
int* spat=new int[35]; //BGO hit pattern in the event, D: 0 or 1 if not hit or hit
int* gpat=new int[35]; //Ge clover mult pattern in the event
int* ppat=new int[35]; //Ge clover pileup pattern in the event
double* Esum=new double[35]; //Ge module energies
double* Tsum=new double[35]; //Variable to make an "average" clover timestamp

//Make an array of energies deposited in each module
double* Msum=new double[35]; 
int modmult=0;
unsigned short modenergy[512]; 

//Isomer calorimetry
double* DelayedModuleSum=new double[35]; //Ge module energies
int DelayedModuleMult=0;
double DelayedModuleEnergy=0;
double LADsum=0;
int LADmult=0;

//DELAYED: Used for Compton suppression and addback in clovers
int* dspat=new int[35]; //BGO hit pattern in the event
int* dgpat=new int[35]; //Ge clover mult pattern in the event
int* dppat=new int[35]; //Ge clover pileup pattern in the event
double* dEsum=new double[35]; //Ge module energies
double* dTsum=new double[35]; //Variable to make an "average" clover timestamp


//ALL: No prompt/delayed time constraint: Used for Compton suppression and addback in clovers
int* aspat=new int[35]; //BGO hit pattern in the event
int* agpat=new int[35]; //Ge clover mult pattern in the event
int* appat=new int[35]; //Ge clover pileup pattern in the event
double* aEsum=new double[35]; //Ge module energies
double* aTsum=new double[35]; //Variable to make an "average" clover timestamp


//Arrays in which to extract the information
//D: arrays we'd like
double* penergy=new double[512]; //Array of prompt Ge module energies
double* ptime=new double[512]; //Array of prompt Ge module energies
double* pid=new double[512]; //Array of prompt Ge module energies
double* lenergy=new double[512]; //Array of delayed LaBr3 channels (for detection of neutrons)

double cnrj; //D: calibrated energy
int detid;
bool pileup;
double timestamp;
double aval,cval; //anode value, cathode value
unsigned int BGOEnergy=0, GeEnergy=0, LaEnergy=0;


///////////////////////////////////////////////////////////////////
/// 		    Looping over all entries in file 	 	        ///
///////////////////////////////////////////////////////////////////

//D: for all entries
for(ULong64_t i = 0; i < entries; i++){
	unsigned int NewTotMult=0; LADmult=0; LADsum=0; DelayedModuleMult=0; DelayedModuleEnergy=0;
	unsigned int Etot=0;
	modmult=0;
	inputtree->GetEntry(i);
	//if (i==1000000) {break; fin->Close(); continue;}	
	TriggerCount++;
	

	//Clear the arrays for storing event prompt and delayed info
	for(int k=1; k <= 34; k++) {
		spat[k]=0; gpat[k]=0; ppat[k]=0; Esum[k]=0; Tsum[k]=0; Msum[k]=0;//prompt 
		dspat[k]=0; dgpat[k]=0; dppat[k]=0; dEsum[k]=0; dTsum[k]=0; //delayed
		aspat[k]=0; agpat[k]=0; appat[k]=0; aEsum[k]=0; aTsum[k]=0; //all; no time constraints

		DelayedModuleSum[k]=0; 
	}
	
	for(int k=0; k <= 34; k++) {penergy[k]=0; pid[k]=0; ptime[k]=0; lenergy[k]=0; modenergy[k]=0;}
	
	LaMult=0; // LaBr3 Mult
	DLaMult=0;
	CMult=0; //Prompt clean Ge mult
	TotalMult=0;
	PromptMult=0;
	DelayedMult=0;
	BGOEnergy=0; LaEnergy=0; GeEnergy=0;
	Double_t tzero,tcat,tano;

	if (size > 250) {cout << "VERY LARGE EVENT SIZE!! EXITING " << size << endl; exit(1);}
        
        for (int j=0; j < mult; j++) //D: for all hits in all detectors
		{
		// Ionisation chamber
		if (ids[j] == 1) 
			{
			double rnum=drand48()-0.5;
			//ANODE fired. fasterchan 195			
			tano=timestamps[j];
			aval=(energies[j]+rnum)/25.0;
			if (nfiles >= 45) {aval*=1.038;} //Correct for the new gain
			energies[j]=aval;			
			//if (mult >= 15) {anode->Fill(aval);}
			anode->Fill(aval);
			}
		if (ids[j] == 2) 
			{
			double rnum=drand48()-0.5;
			//CATHODE fired. fasterchan 191
			tcat=timestamps[j];
			tzero=tcat;
			cval=(energies[j]+rnum)/10000.0;
			energies[j]=cval;
			//if (mult > 15) {cathode->Fill(cval);}
			cathode->Fill(cval);
			}
		}
	
	///////////////////////////////////////////////////////
	/// 		   Gate on ionization chamber 	 		///
	///////////////////////////////////////////////////////

	if (aval >= 900) {continue;} //question: is this condition right?
	//if (aval < 700) {continue;} //condition on the anode
	GoodEventCount++;
        
        //D: For all hits in event
        for (int j=0; j < mult; j++){
			
			// Singles spectra
			detid=idmap_old[ids[j]]; // mapping 20101 to a faster channel stored in idmap.hxx
			double nrj=energies[j];  //D: nrj er energy detected, not calibrated
			double rnum=drand48()-0.5; //D: random nr
			double rnum2=drand48()-0.5; //D:random nr
			timestamp=timestamps[j];
			//pileup=pups[j];
		
			if (isbad[detid]) {continue;}
			//cnrj=(nrj+rnum)*Gains[detid]+Offsets[detid];
		
			if ((Gains2[detid] > 0.5) && (CAL_CORRECT)){ //D:check that calib coeff is not 0	
				cnrj=(nrj+rnum)*Gains[detid]+Offsets[detid]; //D: Calibrating the energy, correcting with random nr
				double newnrj=(cnrj*cnrj*Squared2[detid])+(cnrj*Gains2[detid])+Offsets2[detid]; //D: correction to calib 
				cnrj=newnrj;
			}

			else {cnrj=(nrj+rnum)*Gains[detid]+Offsets[detid];} //D: basic calib
			//D: These BGO works as veto, but not for energy
			//if ((detid==61) && (cnrj < 200)) {cnrj=-1;} //bad BGO detector
			//if ((detid==144) && (cnrj < 130)) {cnrj=-1;} //bad BGO detector
			//if (detid==129) {cnrj=-1;} //bad BGO detector

			if (islabr3[detid]) {
				double nnrj=rnum+(nrj/500.0);
				cnrj=(nnrj*nnrj*Squared[detid])+(nnrj*Gains[detid])+Offsets[detid];
			}
		
		double tdif=tshifts[detid][nfiles]+TOFF+rnum2+(timestamp-tzero)-toffset[detid]; //D: corrected time, give relative time to IC(anode)
		
		if (isge[detid]) {tdif-=walkge(cnrj);}
		if (isbgo[detid]){ //Walk and energy corrections
			tdif-=walkbgo(cnrj);
			if (Gains3[detid] < 0) {
				double cnrj3=((Gains3[detid]*cnrj)+Offsets3[detid]-0.156374)/-0.00377435;
				cnrj=cnrj3;
			}
		}
		
		tspecs[detid]->Fill(tdif);

		if ((isge[detid]) && (tdif > 1000)){
			dirty_ge->Fill(cnrj); 
			getime->Fill(tdif); 
			ET_ge->Fill(tdif,cnrj);
		}

		if (isbgo[detid]) {
			//bgo->Fill(cnrj); 
			//bgotime->Fill(tdif); 
			//ETbgo->Fill(tdif,cnrj);
			if ((tdif > 950) && (tdif < 1050)) {
				especs[detid]->Fill(cnrj);
				bgo->Fill(cnrj); 
				bgotime->Fill(tdif); 
				ETbgo->Fill(tdif,cnrj);
			}
		}

		if (islabr3[detid]) {
			ETla->Fill(tdif,cnrj);
			
			//prompt neutrons in the LaBr3
			if ((tdif >=1008) && (tdif <=1060)) {
				labr3->Fill(cnrj); latime->Fill(tdif);
				lenergy[DLaMult]=cnrj; DLaMult++;
			}
			
			//prompt
			if ((tdif >=990) && (tdif <1008)&& (cnrj > 0)) {LaMult++; LaEnergy+=cnrj;} //Prompt La Mult
			//delayed
			if ((tdif > 1010) && (tdif < 4000)) {LADsum+=cnrj; LADmult++;}
		}


		///////////////////////////////////////////////////////////////
		/// 		  Split events in prompt, delayed, all 	 		///
		///////////////////////////////////////////////////////////////

		// Prompt Ge gammas: Do Compton supression, addback and pileup reject D: question: Compton supression further down?
		if ((isge[detid]) && (tdif >= 945) && (tdif < 1105)&& (cnrj > 0)){
			int mnum=module[detid];
			gpat[mnum]++; Esum[mnum]+=cnrj; Tsum[mnum]+=tdif; //Add Ge's to the module sume, D: here's the addback, so where the clovers of the Ge are summed?
			Msum[mnum]+=cnrj; //D: question; should Msum also be += cnrj?
			GeEnergy+=cnrj;
		}

		// Prompt BGO Energy
		if ((isbgo[detid]) && (tdif >= 945)&& (tdif < 1105) && (cnrj > 0)) {
			int mnum=module[detid]; spat[mnum]++; 
			Msum[mnum]+=(cnrj*10); //Add BGO's to the module sum
			BGOEnergy+=(cnrj*10);
		}
			
		//Delayed Ge
		if ((isge[detid]) && (tdif >= 1050) && (tdif < 1800)){
			int mnum=module[detid];
			dgpat[mnum]++; dEsum[mnum]+=cnrj; dTsum[mnum]+=tdif;
		}
		
		//Delayed BGO
		if ((isbgo[detid]) && (tdif >= 1050) && (tdif < 1800)) {int mnum=module[detid]; dspat[mnum]++;}
		
		//All times, Ge
		if ((isge[detid]) && (cnrj > 0)){
			int mnum=module[detid];
			agpat[mnum]++; aEsum[mnum]+=cnrj; aTsum[mnum]+=tdif; //Add Ge's to the module sume, D: here's the addback, so where the clovers of the Ge are summed? 
		}

		//All times, BGO
		if ((isbgo[detid])) {int mnum=module[detid]; aspat[mnum]++;}
		
	}// End of j loop, reading the event


	///////////////////////////////////////////////////////////////////////////////////////
	/// 		  Loop through Ge modules k;  find dirty/clean and prompt/delayed 	 	///
	///////////////////////////////////////////////////////////////////////////////////////

	modmult=0;

	//134Te
	int isomer_1_134Te = 0;
	int isomer_2_134Te = 0;
	double isomer_energy_1_134Te = 297.0;
	double isomer_energy_2_134Te = 1279.0;
	double time_diff;
	UShort_t clean_multiplicity = 0;

	for(int k = 1; k <= 34; k++){ //loop through the Ge modules in event
		
		// PROMPT
		if (spat[k] || gpat[k]) {NewTotMult++;}
		if (Msum[k] > 5) {modenergy[modmult]=Msum[k]; Etot+=Msum[k]; modmult++;} //20 bins. Max energy 10 MeV.
		

		//Clean prompt Ge only, D: Esum[k] > 5, think this is to avoid noise
		//D: spat[k] == 0; no BGO in module, aka no Compton out of detector
		if ((Esum[k] > 5) && (spat[k]==0) && (gpat[k] >= 1)){
			penergy[CMult]=Esum[k]; //D: question: don't think CMult is updated anywhere?
                	if (gpat[k] > 1) {ptime[CMult]=Tsum[k]/double(gpat[k]);} //Clover avg time
                	else {ptime[CMult]=Tsum[k];} //Normal case. 1 clover element, or phase1.
			pid[CMult]=k; //The module number
 	 		CMult++;
		}

		//D: Clean delayed Ge only, dEsum[k] > 5, think this is to avoid noise
		if ((dEsum[k] > 5) && (dspat[k]==0) && (dgpat[k] >= 1)){
			d_clean_ge->Fill(dEsum[k]);
		}

		double time_diff;
		if ((aEsum[k] > 5) && aspat[k]==0 && (agpat[k] >= 1)){
			clean_ge->Fill(aEsum[k]);

			clean_multiplicity +=1;

			//134Te single energy gate
			if( (aEsum[k]>=isomer_energy_1_134Te-2.5) && (aEsum[k]<isomer_energy_1_134Te+2.5) ){
				if (agpat[k] > 1) {time_diff = aTsum[k]/double(agpat[k]);} //Clover average time
                else {time_diff = aTsum[k];}
                time_isomer_gate_134Te->Fill(time_diff);
                time_isomer_gate_all_134Te->Fill(time_diff);
                isomer_1_134Te = 1;
			}
			//134Te BG-gate
			if( ((aEsum[k]>=isomer_energy_1_134Te-5) && (aEsum[k]<isomer_energy_1_134Te-2.5)) || ((aEsum[k]>=isomer_energy_1_134Te+2.5) && (aEsum[k]<isomer_energy_1_134Te+5)) ) {
				if (agpat[k] > 1) {time_diff = aTsum[k]/double(agpat[k]);} //Clover average time
                else {time_diff = aTsum[k];}
                time_isomer_gate_134Te->Fill(time_diff,-1);
                time_isomer_gate_bg_134Te->Fill(time_diff);
			}
		}

		//134Te create double energy gate
		if ((aEsum[k] > 5) && aspat[k]==0 && (agpat[k] >= 1) && (aEsum[k]>=isomer_energy_2_134Te-2.5) && (aEsum[k]<isomer_energy_2_134Te+2.5)){
			isomer_2_134Te = 1;
		}
	}//end for k loop


	///////////////////////////////////////////////////////
	/// 			  Output file for easy processing 	 			///
	///////////////////////////////////////////////////////

	//Jon's way of writing to file

	//Compress aval and mult into one UShort
	UShort_t new_aval = aval/5;
	unsigned char aval_char, mult_char;

	if (new_aval < 254){aval_char = new_aval;}
	else {aval_char = 254;}

	if (clean_multiplicity < 254){mult_char = clean_multiplicity;}
	else {mult_char = 254;}

	UShort_t avalmult = (mult_char << 8) | aval_char;
	TheEvents[tpointer++] = avalmult;

	//Checking that I get same back fro bitshifting
	//unsigned char aval_char_out, mult_char_out;
	//aval_char_out = avalmult & 0xFF;
	//mult_char_out = avalmult >> 8;

	//Before compressing aval & clean_mult
	//TheEvents[tpointer++]=aval; //anode pulse height; NB I changed order
	//TheEvents[tpointer++]=clean_multiplicity; //Clean Ge multiplicity
	
	tot_mult += clean_multiplicity;
	//tot_mult_char += mult_char;
	//tot_mult_char_out += mult_char_out;

	for (int k=0; k <=34; k++){ //For all Ge modules
		if ((aEsum[k] > 5) && aspat[k]==0 && (agpat[k] >= 1)){ //aspat[k]==0 no BGO fired, agpat: ge-mult larger than 1
			
			TheEvents[tpointer++]=aEsum[k];

			if (agpat[k] > 1) {time_diff = aTsum[k]/double(agpat[k]);} //Clover average time
     	else {time_diff = aTsum[k];}
		
			TheEvents[tpointer++]=time_diff;
		}
	}

	////////////////////////////////////////////////
  	
	NewTotMult+=LaMult;


	//Etot=(BGOEnergy+GeEnergy+LaEnergy)/2000.0;
	Etot+=LaEnergy;
	Etot/=1000;

	//isomer calorimetry
	DelayedModuleMult+=LADmult;
	DelayedModuleEnergy+=LADsum;
	

	///////////////////////////////////////////////
	/// 		  Incrementing spectra  	 	///
	///////////////////////////////////////////////

	//Prompt Ge gamma rays are now stored in the arrays penergy[],ptime[],pid[]
	dmult->Fill(DelayedMult); //How many LONG delayed gammas were there, including BGO?

	if (DLaMult >=2) {
		for (int k=0; k < DLaMult-1; k++) {
			for (int l=k+1; l < DLaMult; l++) {
				lala->Fill(lenergy[k],lenergy[l]);
				lala->Fill(lenergy[l],lenergy[k]);
			}
		}
	}
	
	//Asymmetric Ge-LaBr3 matrix
	if ((DLaMult >=1) && (CMult >=1)){
		for (int k=0; k < DLaMult; k++) {
			for (int l=0; l < CMult; l++) {
				gela->Fill(penergy[l],lenergy[k]);
			}
		}
	}

	if (NewTotMult == 15) {anode15->Fill(aval);}
	if (NewTotMult == 0) {anode0->Fill(aval);}
	if (NewTotMult == 1) {anode1->Fill(aval);}
	if (NewTotMult == 2) {anode2->Fill(aval);}
	if (NewTotMult < 256) {multi[NewTotMult]++;}
	if (Etot < 256) {elti[Etot]++;}
	
	int k=0;
	
	if  ((CMult >=2) && (NewTotMult >= 3)){
		for (int k=0; k < CMult-1; k++){
			for (int l=k+1; l < CMult; l++){	
				//circular time gate
				double tdist=sqrt(((ptime[k]-1025)*(ptime[k]-1025))+((ptime[l]-1025)*(ptime[l]-1025)));
				int deltat=abs(ptime[k]-ptime[l]);
			
				if ((tdist <= 80) && (penergy[k] < 4096) && (penergy[l] < 4096)){
				
					if (aval > 592){ //120 channels
						matp->Fill(penergy[k],penergy[l]);
						matp->Fill(penergy[l],penergy[k]);
					}

					if (aval <= 592){ //230 channels
						matb->Fill(penergy[k],penergy[l]);
						matb->Fill(penergy[l],penergy[k]);
					}
				}
			}
		}
	}
} //End of loop to read the Tree entries

//Write out all the spectra
inputtree->Delete();
fin->Close();
timer.Stop();
Double_t rtime2 = timer.RealTime();
Double_t ctime2 = timer.CpuTime();
cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<endl;
EventCount+=TriggerCount;
cout << "Total Events Read = " << EventCount << endl;
cout << "Total Good Events Processed = " << GoodEventCount << endl;
//RunChange=false; LastFile=false;

if (LastFile || (nfiles > 0)) //Write out the spectra
	{
	cout << "Writing output files " << endl;
	TFile *outputspectrafile = new TFile(OutputSpectraFile.c_str(),"RECREATE");	
	string outfname="bgofits"+nfiles;
ofstream bgofits(outfname);

//TheCube->Write("ddpmult.bin");

//WriteRadware(dirtyge);
//WriteRadware(anode300);
//WriteRadware(anode350);
//WriteRadware(anode400);
//WriteRadware(anode450);
//WriteRadware(anode500);
//WriteRadware(anode550);
//WriteRadware(anode600);
//WriteRadware(anode650);
//WriteRadware(anode700);
//WriteRadware(anode750);
//WriteRadware(anode800);

HK->Write(); //D: not filled, sum energy vs multiplicity
dirty_ge->Write(); //All ge, not Compton-suppressed
clean_ge->Write(); //Ge no time gate, clean (aka Compton-suppressed)
clean_ge_doublegate->Write(); //Ge no time gate, must have seen two isomer lines in event
d_clean_ge->Write(); //Ge delayed, clean


//134Te
time_isomer_gate_134Te->Write(); //Time spectrum gated on isomer energy, bg-subtracted 
time_isomer_gate_bg_134Te->Write(); //BG spectrum for time_isomer_gate
time_isomer_gate_all_134Te->Write(); 

time_isomer_gate_double_134Te->Write(); //Time spectrum gated on isomer energy when required second also seen, bg-subtracted 
time_isomer_gate_double_bg_134Te->Write(); //BG spectrum for time_isomer_gate
time_isomer_gate_double_all_134Te->Write();



getime->Write();
bgotime->Write();
latime->Write();
dge->Write(); //D: not filled, delayed ge 
dmult->Write(); //D: something weird here, only mult 1
totmult->Write(); //D: not filled
cmult->Write(); //D: not filled
sume->Write(); //D: not filled, sum energy, total event energy
bgo->Write();
labr3->Write();
ET_ge->Write();
ETla->Write();
ETbgo->Write();
lala->Write();

//RWMat *lalar=new RWMat(lala);
//lalar->Write();
//delete lalar;

//RWMat *gelar=new RWMat(gela);
//gelar->Write();
//delete gelar;

//RWMat *matbr=new RWMat(matb);
//matbr->Write();
//delete matbr;

//RWMat *matpr=new RWMat(matp);
//matpr->Write();
//delete matpr;

//RWMat *matir=new RWMat(mati);
//matir->Write();
//delete matir;

//RWMat *matlr=new RWMat(matl);
//matlr->Write();
//delete matlr;

WriteRadware(labr3);



int nge=0;
int nbgo=0;
int nla=0;
double avg_slope=0;
double avg_int=0;
double npar=0;
	for (int i=1; i < NDETS; i++)
		{
		if (especs[i]->Integral() > 10) 
			{
			
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Do the new BGO calibration calculation here
int xc2=1400, xc1=117;
double tot=especs[i]->Integral(xc1,xc2); //cout << "Tot = " << tot << " ";
Double_t p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,p7=0;
for (int j=xc2; j > xc1; j--)
	{
	double tot2=especs[i]->Integral(j,xc2);
	if ((tot2 > tot/10000.0) && (p1==0)) {p1=j;} //The top 0.1% point
	if ((tot2 > tot/3162.27766) && (p2==0)) {p2=j;} //The top 1% point
	if ((tot2 > tot/1000.0)&& (p3==0)) {p3=j;} //The top 10% point
	if ((tot2 > tot/316.227766)&& (p4==0)) {p4=j;} //The top 50% point
	if ((tot2 > tot/100.0) && (p5==0)) {p5=j;} //The top 0.1% point
	if ((tot2 > tot/31.6227766) && (p6==0)) {p6=j;} //The top 1% point
	if ((tot2 > tot/10.0)&& (p7==0)) {p7=j;} //The top 10% point
	}
//cout << "BGO " << i << " " << p1 << " " << p2 << " " << p3 << " " << p4 <<endl;
//cout << "-4 " << p1 << endl;
//cout << "-3.5 " <<p2 << endl;
//cout << "-3 " <<p3 << endl;
//cout << "-2.5 " <<p4 << endl;
//cout << "-2 " << p5 << endl;
//cout << "-1.5 " <<p6 << endl;
//cout << "-1 " <<p7 << endl;
Int_t n = 7;
Double_t *x = new Double_t[n];
Double_t *y = new Double_t[n];
Double_t *e = new Double_t[n];

x[0]=-4; x[1]=-3.5; x[2]=-3; x[3]=-2.5; x[4]=-2.0; x[5]=-1.5; x[6]=-1; //log10(cumulative counts integral > 2 MeV)
y[0]=p1; y[1]=p2; y[2]=p3; y[3]=p4; y[4]=p5; y[5]=p6; y[6]=p7; //Channel number of BGO detector
bgofits << y[0] << " " << x[0] << endl;
bgofits << y[1] << " " << x[1] << endl;
bgofits << y[2] << " " << x[2] << endl;
bgofits << y[3] << " " << x[3] << endl;
bgofits << y[4] << " " << x[4] << endl;
bgofits << y[5] << " " << x[5] << endl;
bgofits << y[6] << " " << x[6] << endl;

//TGraphErrors* gr = new TGraphErrors(n,x,y,0,0);
//TF1 *mf=new TF1("mf","pol1"); //Range and number of fit parameters
//gr->Fit(mf,"B"); //Do the fit
//double par0=mf->GetParameter(0); //Get the parameters
//double par1=mf->GetParameter(1); //Get the parameters
//bgofits << i << ", " << par0 << ", " << par1 << endl;
bgofits << "&" << endl;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			especs[i]->Write(); tspecs[i]->Write(); 
			if (isge[i]) {nge++;}
			if (isbgo[i]) {nbgo++;}
			if (islabr3[i]) {nla++;}
			}
		}
cout << nge << "/106 Ge (" << 100*nge/106.0 <<"%) crystals working" << endl;
cout << nbgo << "/58 BGO (" << 100*nbgo/58.0 <<"%) crystals working" << endl;
cout << nla << "/20 LaBr3 (" << 100*nla/20.0 <<"%) crystals working" << endl;
	outputspectrafile->cd();
	outputspectrafile->Close();
	RunChange=false; //reset the flag

	//for (int i=1; i < NDETS; i++)
		//{
		//especs[i]->Reset();
		//tspecs[i]->Reset(); 
		//}
	//anode->Reset();
	//cathode->Reset();
	//bgofits << avg_slope/npar << " " << avg_int/npar <<endl;
	bgofits.close();
	
	}
	

}//end of loop for getting the next file to process


///////////////////////////////////////////////////
/// 			Write output file 				///
///////////////////////////////////////////////////

tpointer_array[0] = tpointer;
tpointer_array[1] = tpointer;

//Write tpointer_array value to file
//ofstream tpointer_outfile("/Volumes/240Pu_d_p/nuball_data/252Cf/Sorted/tpointer.txt", ios::out | ios::binary);
ofstream tpointer_outfile("/Users/dortheagjestvang/Documents/Author_articles/IYR/Data/tpointer.txt", ios::out | ios::binary);
tpointer_outfile.write((char *) tpointer_array, 2*sizeof(UInt_t));
tpointer_outfile.close();

std::cout << "\n" << std::endl;
std::cout << "tot mult before writing: " << tot_mult << std::endl;
//std::cout << "tot mult char: " << tot_mult_char << std::endl;
//std::cout << "tot mult char out: " << tot_mult_char_out << std::endl;
std::cout << "tpointer: " << tpointer << std::endl;

//Write TheEvents to file
//ofstream ofile_edata("/Volumes/240Pu_d_p/nuball_data/252Cf/Sorted/edata.txt", ios::out | ios::binary);
ofstream ofile_edata("/Users/dortheagjestvang/Documents/Author_articles/IYR/Data/edata.txt", ios::out | ios::binary);
ofile_edata.write((char *) TheEvents, bsize*sizeof(UShort_t)); //TheEvents array have all info stored; two bytes per entry. 
ofile_edata.close();


//Read TheEvents_1 from file
//ifstream infile_edata("/Volumes/240Pu_d_p/nuball_data/252Cf/Sorted/edata.txt", ios::in | ios::binary);
ifstream infile_edata("/Users/dortheagjestvang/Documents/Author_articles/IYR/Data/edata.txt", ios::in | ios::binary);
infile_edata.read((char *) TheEvents_1, bsize*sizeof(UShort_t));
infile_edata.close();

//Read tpointer_array from file
//ifstream tpointer_infile("/Volumes/240Pu_d_p/nuball_data/252Cf/Sorted/tpointer.txt", ios::out | ios::binary);
ifstream tpointer_infile("/Users/dortheagjestvang/Documents/Author_articles/IYR/Data/tpointer.txt", ios::out | ios::binary);
tpointer_infile.read((char *) tpointer_array_1, 2*sizeof(UInt_t));
tpointer_infile.close();
std::cout << "tpointer_1: " << tpointer_array_1[0] << std::endl;

//Interpret the data in TheEvents_1
int aval_1;
UShort_t avalmult_1;
unsigned char new_aval_1, mult_1;
double tot_mult_1 = 0;

int i_count = 0;
while(i_count<tpointer_array_1[0]){

	//aval_1=TheEvents_1[i_count++];
	//mult_1=TheEvents_1[i_count++];
	avalmult_1 = TheEvents_1[i_count++];

	new_aval_1 = avalmult_1 & 0xFF;
	mult_1 = avalmult_1 >> 8;


	tot_mult_1 += mult_1;

	double energy_1[mult_1];
	double time_1[mult_1];

	for (int k=0; k < mult_1; k++){
		energy_1[k] = TheEvents_1[i_count++];
		time_1[k] = TheEvents_1[i_count++];
	}
}

std::cout << "tot mult after writing: " << tot_mult_1 << std::endl;
std::cout << "\n" << std::endl;


///////////////////////////////////
///   Bitshifting to compress   ///
///////////////////////////////////


// unsigned short MyShort;
// unsigned char Char1 = 127; // lower byte
// unsigned char Char2 = 128; // upper byte

// // merge two char into short
// MyShort = (Char2 << 8) | Char1;
// unsigned char Char1_out, Char2_out;

// // Split short into two char
// Char1_out = MyShort & 0xFF;
// Char2_out = MyShort >> 8;
// std::cout << "Char1_out, Char2_out: " << +Char1_out << " " << +Char2_out <<  std::endl;

}
/*


g++ -g -o DortheaSort01 DortheaSort01.cxx ` root-config --cflags` `root-config --glibs`


*/


