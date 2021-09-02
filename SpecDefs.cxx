
//All clean Ge (suppressed anti-coincidences)
TH1D *clean_ge = new TH1D("clean_ge","clean_ge",4000,0,4000);
clean_ge->GetXaxis()->SetTitle("Energy [keV]");
clean_ge->GetYaxis()->SetTitle("Counts");

//All clean Ge, double isomer gate
TH1D *clean_ge_doublegate = new TH1D("clean_ge_doublegate","clean_ge_doublegate",4000,0,4000);
clean_ge_doublegate->GetXaxis()->SetTitle("Energy [keV]");
clean_ge_doublegate->GetYaxis()->SetTitle("Counts");

//Delayed clean Ge (so no Compton)
TH1D *d_clean_ge = new TH1D("d_clean_ge","d_clean_ge",4000,0,4000); //All delayed clean Ge (suppressed Compton)
d_clean_ge->GetXaxis()->SetTitle("Energy [keV]");
d_clean_ge->GetYaxis()->SetTitle("Counts");


///////////////////////////////////////////////////////////////////////////////////////
/// 							time_isomer_gate 								   ////
///////////////////////////////////////////////////////////////////////////////////////
//Gate on isomer energy, plot time 

//134Te
TH1D *time_isomer_gate_134Te = new TH1D("time_isomer_gate_134Te","time_isomer_gate_134Te",2000,-1000,3000);
time_isomer_gate_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_gate_bg_134Te = new TH1D("time_isomer_gate_bg_134Te","time_isomer_gate_bg_134Te",2000,-1000,3000);
time_isomer_gate_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_gate_all_134Te = new TH1D("time_isomer_gate_all_134Te","time_isomer_gate_all_134Te",2000,-1000,3000);
time_isomer_gate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_all_134Te->GetYaxis()->SetTitle("Counts");


///////////////////////////////////////////////////////////////////////////////////////
/// 							time_isomer_gate_double 						   ////
///////////////////////////////////////////////////////////////////////////////////////
//Gate on isomer energy when required that second gamma also seen, plot time

//134Te
TH1D *time_isomer_gate_double_134Te = new TH1D("time_isomer_gate_double_134Te","time_isomer_gate_double_134Te",2000,-1000,3000);
time_isomer_gate_double_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_double_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_gate_double_bg_134Te = new TH1D("time_isomer_gate_double_bg_134Te","time_isomer_gate_double_bg_134Te",2000,-1000,3000);
time_isomer_gate_double_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_double_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_gate_double_all_134Te = new TH1D("time_isomer_gate_double_all_134Te","time_isomer_gate_all_134Te",2000,-1000,3000);
time_isomer_gate_double_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_gate_double_all_134Te->GetYaxis()->SetTitle("Counts");


///////////////////////////////////////////////////////////////////////////////////////

//All dirty ge but with pileup bit set
TH1D *pileupge = new TH1D("pileupge","pileupge",4000,0,4000); 

//All dirty Ge
TH1D *dirty_ge = new TH1D("dirty_ge","dirty_ge",8192,0,8192);
dirty_ge->GetXaxis()->SetTitle("Energy [keV]");
dirty_ge->GetYaxis()->SetTitle("Counts"); 

TH1D *bgo = new TH1D("bgo","bgo",4000,0,4000); //All bgo 
TH1D *dirtygeb = new TH1D("dirtygeb","dirtygeb",8192,0,8192); //All dirty Ge 
TH1D *labr3b = new TH1D("labr3b","labr3b",8192,0,8192); //All labr3 

TH1D *anode = new TH1D("anode","anode",8192,0,8192); //IC anode pulse height
TH1D *cathode = new TH1D("cathode","cathode",8192,0,8192); //IC anode pulse height
TH1D *slice = new TH1D("slice","slice",2000,0,2000); //slice of the matrix

TH1D *anode0 = new TH1D("anode0","anode0",8192,0,8192); //IC anode pulse height
TH1D *anode1 = new TH1D("anode1","anode1",8192,0,8192); //IC anode pulse height
TH1D *anode2 = new TH1D("anode2","anode2",8192,0,8192); //IC anode pulse height
TH1D *anode15 = new TH1D("anode15","anode15",8192,0,8192); //IC anode pulse height
TH1D *anode300 = new TH1D("anode300","anode300",8192,0,8192); //IC anode pulse height
TH1D *anode350 = new TH1D("anode350","anode350",8192,0,8192); //IC anode pulse height
TH1D *anode400 = new TH1D("anode400","anode400",8192,0,8192); //IC anode pulse height
TH1D *anode450 = new TH1D("anode450","anode450",8192,0,8192); //IC anode pulse height
TH1D *anode500 = new TH1D("anode500","anode500",8192,0,8192); //IC anode pulse height
TH1D *anode550 = new TH1D("anode550","anode550",8192,0,8192); //IC anode pulse height
TH1D *anode600 = new TH1D("anode600","anode600",8192,0,8192); //IC anode pulse height
TH1D *anode650 = new TH1D("anode650","anode650",8192,0,8192); //IC anode pulse height
TH1D *anode700 = new TH1D("anode700","anode700",8192,0,8192); //IC anode pulse height
TH1D *anode750 = new TH1D("anode750","anode750",8192,0,8192); //IC anode pulse height
TH1D *anode800 = new TH1D("anode800","anode800",8192,0,8192); //IC anode pulse height


TH1D *labr3 = new TH1D("labr3","labr3",8192,0,8192); //All labr3 
TH1D *plabr3 = new TH1D("plabr3","plabr3",8192,0,8192); //All labr3 
TH1D *nlabr3 = new TH1D("nlabr3","nlabr3",8192,0,8192); //All labr3 
TH1D *prompt = new TH1D("prompt","prompt",8192,0,8192);
TH1D *delayed = new TH1D("delayed","delayed",8192,0,8192);

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

TH1D *cmult = new TH1D("cmult","cmult",50,0,50); //Clean Ge Mult 
TH1D *totmult = new TH1D("totmult","totmult",50,0,50); //Event Mult 
TH1D *sume = new TH1D("sume","sume",512,0,512); //Event Sum Energy (100kev/bin) 


TH1D *latime = new TH1D("latime","latime",2000,0,2000); //All Labr3 times
TH1D *bgotime = new TH1D("bgotime","bgotime",2000,0,2000); //All BGO times
TH1D *getime = new TH1D("getime","getime",2000,0,2000); //All Ge times

TH2F *gewalk=new TH2F("gewalk","gewalk",400,0,400,150,0,150);
TH2F *bgowalk=new TH2F("bgowalk","bgowalk",400,0,400,150,0,150);
TH2F *lawalk=new TH2F("lawalk","lawalk",400,0,400,150,0,150);

TH2F *HK=new TH2F("HK","HK",25,0,25,250,0,250);

//Energy-Time Germanium
TH2F *ET_ge=new TH2F("ET_ge","ET_ge",5000,0,5000,2000,0,2000);
ET_ge->GetXaxis()->SetTitle("Time [ns]");
ET_ge->GetYaxis()->SetTitle("Energy [keV]");
ET_ge->SetTitle("Energy-Time germanium");

TH2F *ETla=new TH2F("ETla","ETla",2000,0,2000,2000,0,2000);
TH2F *ETbgo=new TH2F("ETbgo","ETbgo",2000,0,2000,2000,0,2000);

TH1D *TSdif = new TH1D("TSdif","TSdif",5000,0,5000); //count the number of pulses since the previous trigger

//walk spectra
TH2F **walk ; walk = new TH2F*[NCHANS];
for (int i=1; i <= NCHANS; i++)
 {
 TString name=Form("w%d",i);
 //walk[i]=new TH2F(name.Data(),name.Data(),400,0,400,150,0,150); //This line causes bus error, do I need it?
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
 tspecs[i]=new TH1D(name.Data(),name.Data(),4000,0,4000);
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
TH1D *gmult0 = new TH1D("gmult0","gmult0",4096,0,4096);
TH1D *gmult1 = new TH1D("gmult1","gmult1",4096,0,4096);
TH1D *gmult2 = new TH1D("gmult2","gmult2",4096,0,4096);
TH1D *gmult3 = new TH1D("gmult3","gmult3",4096,0,4096);
TH1D *gmult4 = new TH1D("gmult4","gmult4",4096,0,4096);
TH1D *gmult5 = new TH1D("gmult5","gmult5",4096,0,4096);
TH1D *gmult6 = new TH1D("gmult6","gmult6",4096,0,4096);
TH1D *gmult7 = new TH1D("gmult7","gmult7",4096,0,4096);
TH1D *gmult8 = new TH1D("gmult8","gmult8",4096,0,4096);
TH1D *gmult9 = new TH1D("gmult9","gmult9",4096,0,4096);
TH1D *gmult10 = new TH1D("gmult10","gmult10",4096,0,4096);
TH1D *gmult11 = new TH1D("gmult11","gmult11",4096,0,4096);
TH1D *gmult12 = new TH1D("gmult12","gmult12",4096,0,4096);
TH1D *gmult13 = new TH1D("gmult13","gmult13",4096,0,4096);
TH1D *gmult14 = new TH1D("gmult14","gmult14",4096,0,4096);
TH1D *gmult15 = new TH1D("gmult15","gmult15",4096,0,4096);
TH1D *gmult16 = new TH1D("gmult16","gmult16",4096,0,4096);
TH1D *gmult17 = new TH1D("gmult17","gmult17",4096,0,4096);
TH1D *gmult18 = new TH1D("gmult18","gmult18",4096,0,4096);
TH1D *gmult19 = new TH1D("gmult19","gmult19",4096,0,4096);
TH1D *gmult20 = new TH1D("gmult20","gmult20",4096,0,4096);

TH1D *etot0 = new TH1D("etot0","etot0",4096,0,4096);
TH1D *etot1 = new TH1D("etot1","etot1",4096,0,4096);
TH1D *etot2 = new TH1D("etot2","etot2",4096,0,4096);
TH1D *etot3 = new TH1D("etot3","etot3",4096,0,4096);
TH1D *etot4 = new TH1D("etot4","etot4",4096,0,4096);
TH1D *etot5 = new TH1D("etot5","etot5",4096,0,4096);
TH1D *etot6 = new TH1D("etot6","etot6",4096,0,4096);
TH1D *etot7 = new TH1D("etot7","etot7",4096,0,4096);
TH1D *etot8 = new TH1D("etot8","etot8",4096,0,4096);
TH1D *etot9 = new TH1D("etot9","etot9",4096,0,4096);
TH1D *etot10 = new TH1D("etot10","etot10",4096,0,4096);
TH1D *etot11 = new TH1D("etot11","etot11",4096,0,4096);
TH1D *etot12 = new TH1D("etot12","etot12",4096,0,4096);
TH1D *etot13 = new TH1D("etot13","etot13",4096,0,4096);
TH1D *etot14 = new TH1D("etot14","etot14",4096,0,4096);
TH1D *etot15 = new TH1D("etot15","etot15",4096,0,4096);
TH1D *etot16 = new TH1D("etot16","etot16",4096,0,4096);
TH1D *etot17 = new TH1D("etot17","etot17",4096,0,4096);
TH1D *etot18 = new TH1D("etot18","etot18",4096,0,4096);
TH1D *etot19 = new TH1D("etot19","etot19",4096,0,4096);
TH1D *etot20 = new TH1D("etot20","etot20",4096,0,4096);

TH1D *mode0 = new TH1D("mode0","mode0",4096,0,4096);
TH1D *mode1 = new TH1D("mode1","mode1",4096,0,4096);
TH1D *mode2 = new TH1D("mode2","mode2",4096,0,4096);
TH1D *mode3 = new TH1D("mode3","mode3",4096,0,4096);
TH1D *mode4 = new TH1D("mode4","mode4",4096,0,4096);
TH1D *mode5 = new TH1D("mode5","mode5",4096,0,4096);
TH1D *mode6 = new TH1D("mode6","mode6",4096,0,4096);
TH1D *mode7 = new TH1D("mode7","mode7",4096,0,4096);
TH1D *mode8 = new TH1D("mode8","mode8",4096,0,4096);
TH1D *mode9 = new TH1D("mode9","mode9",4096,0,4096);
TH1D *mode10 = new TH1D("mode10","mode10",4096,0,4096);
TH1D *mode11 = new TH1D("mode11","mode11",4096,0,4096);
TH1D *mode12 = new TH1D("mode12","mode12",4096,0,4096);
TH1D *mode13 = new TH1D("mode13","mode13",4096,0,4096);
TH1D *mode14 = new TH1D("mode14","mode14",4096,0,4096);
TH1D *mode15 = new TH1D("mode15","mode15",4096,0,4096);
TH1D *mode16 = new TH1D("mode16","mode16",4096,0,4096);
TH1D *mode17 = new TH1D("mode17","mode17",4096,0,4096);
TH1D *mode18 = new TH1D("mode18","mode18",4096,0,4096);
TH1D *mode19 = new TH1D("mode19","mode19",4096,0,4096);
TH1D *mode20 = new TH1D("mode20","mode20",4096,0,4096);

