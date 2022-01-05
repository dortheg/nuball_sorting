#include <sstream>
#include "TImage.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TROOT.h"
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
#include "Spec.hxx"
#include "Spec.cxx"
#include "RWMat.hxx"
#include "RWMat.cxx"
#include "WriteRWSpec.cxx"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum2.h"
#include "TSpectrum2Painter.h"
#include "TSpectrum2Transform.h"
#include "TSpectrum3.h"
#include "TSpectrumFit.h"
#include "TSpectrum.h"
#include "TSpectrumTransform.h"
#include "TGraphErrors.h"

bool GateOnDelayed=true;
using namespace std;

//Returns the background subtracted spectrum for TH1D
void backsub(TH1D* TheSpectrum, int nsmooth)
{
int j;
int nbins=TheSpectrum->GetNbinsX();
Double_t *source=new Double_t[nbins];
TSpectrum *s = new TSpectrum();
for (j=0;j<nbins;j++) source[j]=TheSpectrum->GetBinContent(j+1);
s->Background(source,nbins,nsmooth,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE,TSpectrum::kBackSmoothing3,kFALSE);
for (j=0;j<nbins;j++) {TheSpectrum->SetBinContent(j+1, TheSpectrum->GetBinContent(j+1)-source[j]);}
s->Delete();
delete [] source;
}

double FWHM(double energy)
{
double val=sqrt(1.059*1.059+(2.814*2.814*energy/1000.0));
return val;
}

double sigma(double energy)
{
double val=sqrt(1.059*1.059+(2.814*2.814*energy/1000.0));
double val2=val/2.355; //convert to sigma
return val2;
}

//Returns the background subtracted spectrum for RW Spec
void backsubS(Spec* TheSpec, int nsmooth)
{
int j;
int nbins=4096;
Double_t *source=new Double_t[nbins];
TSpectrum *s = new TSpectrum();
for (j=0;j<nbins;j++) source[j]=TheSpec->GetChan(j);
s->Background(source,nbins,nsmooth,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE,TSpectrum::kBackSmoothing3,kFALSE);
for (j=0;j<nbins;j++) {TheSpec->SetChan(j, TheSpec->GetChan(j)-source[j]);}
s->Delete();
delete [] source;
}

int main(int argc, char **argv)
{
if (argc <= 1) {cout << "type fname " << endl; exit(1);}

string fname1=argv[1];

string froot=fname1.substr(0,fname1.length()-4);

//int thresh=70;
int thresh=5;
Spec *eff=new Spec(4096);
eff->Read("eff-nuball.spe");

Spec *TheSpec=new Spec(4096);
Spec *TheSpecRaw=new Spec(4096);
Spec *TheSpecBkg=new Spec(4096);
Spec *Result=new Spec(4096);
TheSpec->Read("xe140e.spe");
TheSpecRaw->Read("xe140u.spe");
for (int i=508; i <= 514; i++) {TheSpec->SetChan(i,0);} //We're not going to be fitting 511's

TH1D *spec = TheSpec->Convert2Root(); //Keep the name and keep the size

//TH1D *specb = new TH1D("specb","specb",4096,0,4096); 
//TH1D *specraw = new TH1D("specraw","specraw",4096,0,4096); 

for (int i=0; i < 4096; i++) 
{
spec->SetBinContent(i+1,TheSpec->GetChan(i));
specb->SetBinContent(i+1,TheSpec->GetChan(i));
}
for (int i=0; i < 4096; i++) {specraw->SetBinContent(i+1,TheSpecRaw->GetChan(i));}


backsub(spec,8);
for (int i=0; i < 4096; i++) {TheSpecBkg->SetChan(i,specb->GetBinContent(i+1));}
TheSpecBkg->Write("bkg.spe");

Result=TheSpec->AddSpecs(TheSpecBkg,-1.0);
Result->Write("result.spe");
double widths[20];
int n=17;
double peaks[20]={376.805, 457.743, 582.514, 566.551, 607.047, 679.229, 728.451, 551.391, 510.516, 563.801, 619.757, 768.05, 752.689, 659.3, 585, 681, 731};
// contaminants at the end

for (int i=0; i < n; i++) {widths[i]=sigma(peaks[i]); cout << sigma(peaks[i]) << endl;}

TApplication theApp("App", &argc, argv);
TCanvas* c1 = new TCanvas("c1","Spectrum",200,10,700,500);
spec->GetXaxis()->SetRange(300,850);
spec->Draw();

//Vector of the individual fits
vector <TF1 *> g;
for (int i=0; i < n; i++)
{
string gs="g"+i;
g.push_back(new TF1(gs.c_str(),"gaus",peaks[i]-5,peaks[i]+5));   
g[i]->SetParameter(0,100);
g[i]->FixParameter(1,peaks[i]+0.5);
g[i]->FixParameter(2,widths[i]);
if (i >= 14) {g[i]->SetLineColor(8);}
spec->Fit(g[i],"R+B");
}

Double_t par[42];

stringstream ss;
for (int i=0; i < n; i++) 
{
ss << "gaus(" << i*3 <<")+";
}
string gstring=ss.str();
gstring=gstring.substr(0,gstring.size()-1); //lop off the last "+"

//cout << gstring << endl;

//TF1 *total = new TF1("total","gaus(0)+gaus(3)+gaus(6)",300,700);
TF1 *total = new TF1("total",gstring.c_str(),300,850);
int peaknum=2;
int areanum=peaknum*3;
int posnum=peaknum*3;
int widnum=peaknum*3;

for (int i=0; i < n; i++) 
{
int area=i*3;
int pos=i*3+1;
int wid=i*3+2;
total->SetParameter(area,100);
total->FixParameter(pos,peaks[i]+0.5);
total->FixParameter(wid,widths[i]);
}

total->SetLineColor(8);
for (int i=0; i < n; i++) 
{
double fitval=g[i]->GetParameter(0);
double intval=fitval*1000.0/eff->GetChan(int(peaks[i])); //Real intensity after effcor
double lo=peaks[i]+0.5-(widths[i]*2.0);
double hi=peaks[i]+0.5+(widths[i]*2.0);
specraw->GetXaxis()->SetRange(lo,hi);
double stat_error=specraw->Integral(lo,hi);

cout << peaks[i] << " " << intval << " " << (sqrt(stat_error)/fitval)*intval << endl;

}

//g[0]->GetParameters(&par[0]);
//g[1]->GetParameters(&par[3]);
//g[2]->GetParameters(&par[6]);
   
//total->SetParameters(par);

//spec->Fit(total,"R+B");

theApp.Run();

}

/*


g++ -g -o SpecFit2 SpecFit2.cxx  ` root-config --cflags` `root-config --glibs` -lSpectrum


*/


