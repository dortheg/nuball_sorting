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
Spec *TheSpec=new Spec(4096);
Spec *TheSpecRaw=new Spec(4096);
Spec *TheSpecBkg=new Spec(4096);
Spec *Result=new Spec(4096);
TheSpec->Read("xe140e.spe");
TheSpecRaw->Read("xe140u.spe");
for (int i=508; i <= 514; i++) {TheSpec->SetChan(i,0);} //We're not going to be fitting 511's

TH1D *spec = new TH1D("spec","spec",4096,0,4096); 
TH1D *specb = new TH1D("specb","specb",4096,0,4096); 
TH1D *specraw = new TH1D("spec","spec",4096,0,4096); 
for (int i=0; i < 4096; i++) 
{
spec->SetBinContent(i+1,TheSpec->GetChan(i));
specb->SetBinContent(i+1,TheSpec->GetChan(i));
}
//for (int i=0; i < 4096; i++) {specraw->SetBinContent(i+1,TheSpecRaw->GetChan(i));}


backsub(spec,8);
for (int i=0; i < 4096; i++) {TheSpecBkg->SetChan(i,specb->GetBinContent(i+1));}
TheSpecBkg->Write("bkg.spe");

Result=TheSpec->AddSpecs(TheSpecBkg,-1.0);
Result->Write("result.spe");
double widths[20];
double peaks[20]={376.805, 457.743, 582.514, 566.551, 607.047, 679.229, 728.451, 551.391, 510.516, 563.801, 619.757, 768.05, 752.689, 659.3};
peaks[1]=585;
for (int i=0; i < 14; i++) {widths[i]=sigma(peaks[i]); cout << sigma(peaks[i]) << endl;}

TApplication theApp("App", &argc, argv);
TCanvas* c1 = new TCanvas("c1","Spectrum",200,10,700,500);
spec->GetXaxis()->SetRange(300,850);
spec->Draw();
vector <TF1 *> TheFits;

   Double_t par[42];
   TF1 *g1 = new TF1("g1","gaus",peaks[0]-5,peaks[0]+5);
   g1->FixParameter(1,peaks[0]+0.5);
   g1->FixParameter(2,widths[0]);

   TF1 *g2 = new TF1("g2","gaus",peaks[1]-5,peaks[1]+5);
   g2->FixParameter(1,peaks[1]+0.5);
   g2->FixParameter(2,widths[1]);
   
   TF1 *g3 = new TF1("g3","gaus",peaks[2]-5,peaks[2]+5);
   g3->FixParameter(1,peaks[2]+0.5);
   g3->FixParameter(2,widths[2]);
   
   TF1 *total = new TF1("total","gaus(0)+gaus(3)+gaus(6)",300,700);
   //total->FixParameter(4,peaks[0]+0.5);
   total->FixParameter(5,widths[0]);
   
   total->FixParameter(7,peaks[2]+0.5);
   total->FixParameter(8,widths[2]);
  
   total->SetLineColor(8);
   //spec->Fit(g1,"RB");
   //spec->Fit(g2,"R+B");
   //spec->Fit(g3,"R+B");
   //spec->Fit(g4,"R+B");
   //spec->Fit(g5,"R+B");
   //spec->Fit(g6,"R+B");
   //spec->Fit(g7,"R+B");
   //spec->Fit(g8,"R+B");
   //spec->Fit(g9,"R+B");
   //spec->Fit(g10,"R+B");
   //spec->Fit(g11,"R+B");
   //spec->Fit(g12,"R+B");
   //spec->Fit(g13,"R+B");
   //spec->Fit(g14,"R+B");
   
   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);
   g3->GetParameters(&par[6]);
   
   total->SetParameters(par);
   spec->Fit(total,"R+B");

theApp.Run();

}

/*


g++ -g -o SpecFit SpecFit.cxx  ` root-config --cflags` `root-config --glibs` -lSpectrum


*/


