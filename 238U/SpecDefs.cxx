TH1D *cleange = new TH1D("cleange","cleange",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *pcleange = new TH1D("pcleange","pcleange",8192,0,8192); //Prompt clean Ge (suppressed anti-coincidences)
TH1D *dcleange = new TH1D("dcleange","dcleange",8192,0,8192); //Delayed clean Ge (suppressed anti-coincidences)
TH1D *lcleange = new TH1D("lcleange","lcleange",8192,0,8192); //LongDelayed clean Ge (suppressed anti-coincidences)

TH1D *nlab = new TH1D("nlab","nlab",8192,0,8192);


TH1D *ringw = new TH1D("ringw","ringw",50,0,50);
TH1D *back = new TH1D("back","back",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *middle = new TH1D("middle","middle",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *forward = new TH1D("forward","forward",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *evt = new TH1D("evt","evt",100,0,100);
TH2F *select=new TH2F("select","select",100,0,100,50,0,50);

TH1D *rw3 = new TH1D("rw3","rw3",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw4 = new TH1D("rw4","rw4",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw5 = new TH1D("rw5","rw5",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw6 = new TH1D("rw6","rw6",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw7 = new TH1D("rw7","rw7",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw8 = new TH1D("rw8","rw8",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw9 = new TH1D("rw9","rw9",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw10 = new TH1D("rw10","rw10",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw11 = new TH1D("rw11","rw11",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw12 = new TH1D("rw12","rw12",8192,0,8192); //All clean Ge (suppressed anti-coincidences)
TH1D *rw13 = new TH1D("rw13","rw13",8192,0,8192); //All clean Ge (suppressed anti-coincidences)

TH1D *pileupge = new TH1D("pileupge","pileupge",8192,0,8192); //All dirty ge but with pileup bit set

TH1D *dirtyge = new TH1D("dirtyge","dirtyge",8192,0,8192); //All dirty Ge 
TH1D *pdirtyge = new TH1D("pdirtyge","pdirtyge",8192,0,8192); //All dirty Ge 
TH1D *ddirtyge = new TH1D("ddirtyge","ddirtyge",8192,0,8192); //All dirty Ge 
TH1D *ldirtyge = new TH1D("ldirtyge","ldirtyge",8192,0,8192); //All dirty Ge 

TH1D *bgo = new TH1D("bgo","bgo",8192,0,8129); //All bgo 
TH1D *anti = new TH1D("anti","anti",8192,0,8129); //Ge in coincidence with BGO 
TH1D *dirtygeb = new TH1D("dirtygeb","dirtygeb",8192,0,8192); //All dirty Ge 
TH1D *labr3b = new TH1D("labr3b","labr3b",8192,0,8192); //All labr3 

TH1D *labr3 = new TH1D("labr3","labr3",8192,0,8192); //All labr3 
TH1D *plabr3 = new TH1D("plabr3","plabr3",8192,0,8192); //All labr3 
TH1D *nlabr3 = new TH1D("nlabr3","nlabr3",8192,0,8192); //All labr3 
TH1D *prompt = new TH1D("prompt","prompt",8192,0,8192);
TH1D *delayed = new TH1D("delayed","delayed",8192,0,8192);

TH1D *dmm = new TH1D("dmm","dmm",50,0,50); //Delayed module multiplicity
TH1D *dme = new TH1D("dme","dme",100,0,100); //Delayed module sum energy 100 keV per channel
TH1D *dmm1 = new TH1D("dmm1","dmm1",4096,0,4096); //Delayed module multiplicity
TH1D *dmm2 = new TH1D("dmm2","dmm2",4096,0,4096); //Delayed module multiplicity
TH1D *dmm3 = new TH1D("dmm3","dmm3",4096,0,4096); //Delayed module multiplicity
TH1D *dmm4 = new TH1D("dmm4","dmm4",4096,0,4096); //Delayed module multiplicity
TH1D *dmm5 = new TH1D("dmm5","dmm5",4096,0,4096); //Delayed module multiplicity
TH1D *dmm6 = new TH1D("dmm6","dmm6",4096,0,4096); //Delayed module multiplicity
TH1D *dme1 = new TH1D("dme1","dme1",4096,0,4096); //Delayed module sum energy 100 keV per channel
TH1D *dme2 = new TH1D("dme2","dme2",4096,0,4096); //Delayed module sum energy 100 keV per channel
TH1D *dme3 = new TH1D("dme3","dme3",4096,0,4096); //Delayed module sum energy 100 keV per channel
TH1D *dme4 = new TH1D("dme4","dme4",4096,0,4096); //Delayed module sum energy 100 keV per channel
TH1D *dme5 = new TH1D("dme5","dme5",4096,0,4096); //Delayed module sum energy 100 keV per channel


TH1D *preprompt = new TH1D("preprompt","preprompt",4000,0,4000); //300ns window before prompt for 1/60th of the data
//or 250 files (minus the run225-250) of 60 seconds. This will give delayed 691 keV per second in nuball. or the random
// rate of 691 keV's to subtract off. X prompt events with 1200ns window gives Y seconds of data corresponding to Z randoms.

TH1D *r2 = new TH1D("r2","r2",4000,0,4000); //All labr3 
TH1D *r3 = new TH1D("r3","r3",4000,0,4000); //All labr3 
TH1D *r4 = new TH1D("r4","r4",4000,0,4000); //All labr3 
TH1D *r2r2 = new TH1D("r2r2","r2r2",4000,0,4000); //All labr3 
TH1D *orings = new TH1D("orings","orings",4000,0,4000); //All labr3 

TH1D *tr3 = new TH1D("tr3","tr3",4096,0,4096);  
TH1D *tr4 = new TH1D("tr4","tr4",4096,0,4096);  
TH1D *tr5 = new TH1D("tr5","tr5",4096,0,4096);  
TH1D *tr6 = new TH1D("tr6","tr6",4096,0,4096); 

TH1D *ce148 = new TH1D("ce148","ce148",4096,0,4096); 
TH1D *ce148b = new TH1D("ce148b","ce148b",4096,0,4096); 

TH1D *dmult = new TH1D("dmult","dmult",50,0,50); //All Ge times
TH1D *dsume = new TH1D("dsume","dsume",512,0,512); //Sum energy 100 keV per bin
TH1D *dge = new TH1D("dge","dge",4000,0,4000); //All delayed Ge 
TH1D *ldge = new TH1D("ldge","ldge",4000,0,4000); //Long delayed Ge 

TH1D *cmult = new TH1D("cmult","cmult",50,0,50); //Clean Ge Mult 
TH1D *totmult = new TH1D("totmult","totmult",50,0,50); //Event Mult 
TH1D *sume = new TH1D("sume","sume",512,0,512); //Event Sum Energy (100kev/bin) 
TH1D *pmult = new TH1D("pmult","pmult",50,0,50); //prompt mult. plus or minus 1,2,3 fwhm 
TH1D *tmult = new TH1D("tmult","tmult",50,0,50); //total multiplicity in 400ns window
TH1D *psume = new TH1D("psume","psume",512,0,512); //prompt sume energy 100 keV per bin 
TH1D *lmult = new TH1D("lmult","lmult",50,0,50); //prompt LaBr3 mult.


TH1D *latime = new TH1D("latime","latime",2000,0,2000); //All Labr3 times
TH1D *bgotime = new TH1D("bgotime","bgotime",2000,0,2000); //All BGO times
TH1D *getime = new TH1D("getime","getime",2000,0,2000); //All Ge times

TH2F *gewalk=new TH2F("gewalk","gewalk",400,0,400,150,0,150);
TH2F *bgowalk=new TH2F("bgowalk","bgowalk",400,0,400,150,0,150);
TH2F *lawalk=new TH2F("lawalk","lawalk",400,0,400,150,0,150);

TH2F *HK=new TH2F("HK","HK",25,0,25,250,0,250);
TH2F *HK2=new TH2F("HK2","HK2",25,0,25,250,0,250);
TH2F *HK3=new TH2F("HK3","HK3",25,0,25,250,0,250);
TH2F *dHK=new TH2F("dHK","dHK",25,0,25,250,0,250);
TH2F *ETge=new TH2F("ETge","ETge",1200,0,1200,2000,0,2000);
TH2F *ETla=new TH2F("ETla","ETla",1200,0,1200,2000,0,2000);
TH2F *ETbgo=new TH2F("ETbgo","ETbgo",1200,0,1200,2000,0,2000);
TH2F *ETgea=new TH2F("ETgea","ETgea",1200,0,1200,2000,0,2000);

TH1D *TSdif = new TH1D("TSdif","TSdif",5000,0,5000); //count the number of pulses since the previous trigger

//walk spectra
TH2F **walk ; walk = new TH2F*[NCHANS];
for (int i=1; i <= NCHANS; i++)
 {
 TString name=Form("w%d",i);
 walk[i]=new TH2F(name.Data(),name.Data(),400,0,400,150,0,150);
 }
// RAW ENERGY AND TIME SINGLES SPECTRA
TH1D **especs ; especs = new TH1D*[NCHANS];
TString name,mst;
for (int i=1; i <= NDETS; i++)
 {
 name=Form("e%d",i);//Spectrum i should be fasterchan here for coherent numbering
 especs[i]=new TH1D(name.Data(),name.Data(),NCHANS,0,NCHANS);
 }
//time spectra
TH1D **tspecs ; tspecs = new TH1D*[NCHANS];
for (int i=1; i <= NDETS; i++)
 {
 name=Form("t%d",i);
 tspecs[i]=new TH1D(name.Data(),name.Data(),400,0,400);
 }
 
//prompt spectra
TH1D **pspecs ; pspecs = new TH1D*[NCHANS];
for (int i=1; i <= NDETS; i++)
 {
 name=Form("p%d",i);
 pspecs[i]=new TH1D(name.Data(),name.Data(),NCHANS,0,NCHANS);
 }
 
//labr3 neutron spectra
TH1D **nspecs ; nspecs = new TH1D*[NCHANS];
for (int i=1; i <= NDETS; i++)
 {
 name=Form("n%d",i);
 nspecs[i]=new TH1D(name.Data(),name.Data(),NCHANS,0,NCHANS);
 }
 
