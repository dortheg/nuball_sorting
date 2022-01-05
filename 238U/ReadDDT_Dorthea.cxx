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

short lookup[2048][2048];
int coincs[20][2];
int lookup1D[2048];

double NFWHM=1.0;


int main(int argc, char **argv)
{
if (argc <= 2) {cout << "type ddtcube name, .iso file " << endl; exit(1);}
double eff1,eff2;
string fname1=argv[1];
string fname2=fname1.substr(0,fname1.length()-4);
fname2=fname2+"-sub.bin";
string froot2=fname1.substr(2,fname1.length()-6);
string halflife;
string fname3;
if (argc == 3) {fname3=argv[2];}
if (argc > 3) {cout << "too many parameters" << endl; exit(1);}

CubeDDT *Cube1=new CubeDDT("",EBINS,TBINS,10);
CubeDDT *Cube2=new CubeDDT("",EBINS,TBINS,10);
Cube1->Read(fname1);
Cube2->Read(fname2);

// Gamma gamm matrix
TH2F* ggm = new TH2F("ggm","ggm",EBINS,0,EBINS,EBINS,0,EBINS); //gated matrix


//time spectra
TH1D **tspecs; tspecs = new TH1D*[EBINS];
for (int i=1; i <= TBINS; i++)
 {
 TString name=Form("t%d",i);
 tspecs[i]=new TH1D(name.Data(),name.Data(),EBINS,0,EBINS);
 }

//if (argc==2)
//{
//cout << "Enter isomer coincidences file (e.g. kr94.iso) ";
//cin >> fname3;
//}

string froot=fname3.substr(0,fname3.length()-4);
ifstream infile(fname3.c_str());
if (!infile) {cout << "File doesn't exist " << fname3 << endl; exit(1);}

int e1,e2;
int k=0;
while(infile >> e1 >> e2)
{
cout << e1 << " " << e2 << endl;
if (e1 > 0) {coincs[k][0]=e1; coincs[k][1]=e2; k++; }
k++;
}
infile.close();

double count1=0,count2=0,count3=0,count4=0;
for (int k=0; k < 20; k++)
	{
	int e1=coincs[k][0];
	int e2=coincs[k][1];

	if (e1 > 0)
		{
		double A=1.059;
		double B=2.814;
		double PW1=(0.5+(sqrt(A*A+(B*B*e1/1000.0))/1.0))*(NFWHM/2.0); //1/2 FWHM.
		double PW2=(0.5+(sqrt(A*A+(B*B*e2/1000.0))/1.0))*(NFWHM/2.0); //1/2 FWHM.
		int low1=int(e1-PW1);
		int high1=int(e1+PW1);
		int low2=int(e2-PW2);
		int high2=int(e2+PW2);
		
		int blow1=int(e1-(PW1*3.0));
		int bhigh1=int(e1+(PW1*3.0));
		int blow2=int(e2-(PW2*3.0));
		int bhigh2=int(e2+(PW2*3.0));
		
		for (int i=blow1; i <=bhigh1; i++)
			{
			for (int j=blow2; j <=bhigh2; j++)
				{
				lookup1D[j]=2; //background regions either side of the 2nd peak in the list
				//Can add the two gates together if uncomment this next line
				//lookup1D[j]=1;
				lookup[i][j]=2;
				lookup[j][i]=2;
				cout << k << " " << i << "  " << j << " 2" <<endl;
				}
			}
		
		for (int i=low1; i <=high1; i++)
			{
			lookup1D[i]=1; //1D gate on first energy in the list
			for (int j=low2; j <=high2; j++)
				{
				//Can add the two gates together if uncomment this next line
				//lookup1D[j]=1;
				lookup[i][j]=1;
				lookup[j][i]=1;
				cout << k << " " << i << "  " << j << " 1" <<endl;
				}
			}
		}
	}
//cout << "gate width 1 = " << count1 << " gate width 2 " << count2 <<endl;

int bcount=0;
int pcount=0;

// Classic radware BK sub for each time bin
TH1D *px = new TH1D("px","px",EBINS,0,EBINS);  
TH1D *bx = new TH1D("bx","bx",EBINS,0,EBINS);  


string fname4=froot+froot2+".dat";
string fname5=froot+froot2+"2.dat";
string fname6=froot+froot2+"3.dat";
string fname7=froot+froot2+"4.dat";


ofstream outfile(fname4.c_str());
ofstream outfile2(fname5.c_str());
ofstream outfile3(fname6.c_str());
ofstream outfile4(fname7.c_str());

double peak_integral=0;
double bk_integral=0;
double bk_integral2=0;
double back_integral=0;
double tot_integral=0;
double tot_integral2=0;

double beta=0;

for (int k=0; k < TBINS; k++)
	{
double pcounts=0;
double raw_counts=0;
double pcounts_corr=0;
double raw_counts_corr=0;
//Note that Cube1 is the raw cube and Cube2 is the bk subtracted cube
	for (int i=0; i < EBINS; i++)
		{
		for (int j=0; j < EBINS; j++)
			{
			if (lookup[i][j]==1) //peak region
				{
				//bk sub cube
				float pval=Cube2->Get(i,j,k);
				pcounts+=pval;	
				
				//raw cube
				float raw_val=Cube1->Get(i,j,k);
				raw_counts+=raw_val;
				}
			if (lookup[i][j]==2) //bk regions
				{
				//bk sub cube
				float pval=Cube2->Get(i,j,k);
				pcounts_corr+=pval;	

				//raw cube
				float raw_val=Cube1->Get(i,j,k);
				raw_counts_corr+=raw_val;
				}
			}
		}
//NOTE: VARIABLE NAMES ARE CONFUSING and NEEDS CHANGING
	//raw_counts = peak integral in raw cube
	//pcounts = peak integral in bk sub cube
	double pcor=(pcounts_corr/8.0); //Integral bk regions bk sub cube (should be almost zero)
	double bcor=(raw_counts_corr/8.0); //Integral bk regions raw cube
	cout << k << " " << pcounts << " " << (raw_counts-bcor) << " " << endl;//pcounts from bk sub cube. (peak - bk) from raw cube.
// OUTPUT TO FILE HERE
	outfile << k << " " << pcounts << " " << sqrt(raw_counts) << endl;//bk sub cube peak counts
	pcounts=pcounts-pcor;
	outfile2 << k << " " << (pcounts-pcor) << " " << sqrt(raw_counts) << endl;//corrected peak counts from 8 adjacent zones
	outfile3 << k << " " << raw_counts << " " << sqrt(raw_counts) << endl;//raw cube (peak+bk) counts. TOTAL
	//cout << "t=" << k*10 << " old pcounts = " << pcounts << " " << pcor << " correction (%) " << 100.0*(pcounts-pcor)/pcounts << endl;
	//should be the same
	outfile4 << k << " "<< (raw_counts-bcor) << " " << sqrt(raw_counts) << endl;//Only Raw cube peak estimate
	}

cout << "Writing output data file " << fname4 << endl;
outfile.close();
//RWMat *ggmr=new RWMat(ggm);
//ggmr->Write();
//delete ggmr;

}


/*


g++ -g -o ReadDDT_Dorthea ReadDDT_Dorthea.cxx ` root-config --cflags` `root-config --glibs` -lSpectrum


*/


