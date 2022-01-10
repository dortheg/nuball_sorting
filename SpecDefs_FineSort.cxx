//All single, clean ge-gammas
TH1D *single_gamma = new TH1D("single_gamma","single_gamma",5000,0,5000);
single_gamma->GetXaxis()->SetTitle("Energy [keV]");
single_gamma->GetYaxis()->SetTitle("Counts");

//All double, clean ge-gammas
TH2F *double_gamma = new TH2F("double_gamma","double_gamma",4000,0,4000, 4000, 0, 4000);
double_gamma->GetXaxis()->SetTitle("Energy [keV]");
double_gamma->GetYaxis()->SetTitle("Energy [keV]");
double_gamma->SetOption("colz");

//Multiplicity distribution
TH1D *mult_distr = new TH1D("mult_distr","mult_distr",20,0,20);
mult_distr->GetXaxis()->SetTitle("Multiplicity");
mult_distr->GetYaxis()->SetTitle("Counts");

//Distribution of hits in array
TH1D *hit_distr = new TH1D("hit_distr","hit_distr",20,0,20);
hit_distr->GetXaxis()->SetTitle("Hits");
hit_distr->GetYaxis()->SetTitle("Counts");


//////////////////////////////////////////////////////////////
//							 134Te   						//
//////////////////////////////////////////////////////////////

//Gate on isomer_1
TH1D *time_isomer_1_gate_134Te = new TH1D("time_isomer_1_gate_134Te","time_isomer_1_gate_134Te",2000,-1000,3000);
time_isomer_1_gate_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_1_gate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_1_gate_bg_134Te = new TH1D("time_isomer_1_gate_bg_134Te","time_isomer_1_gate_bg_134Te",2000,-1000,3000);
time_isomer_1_gate_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_1_gate_bg_134Te->GetYaxis()->SetTitle("Counts");


TH1D *time_isomer_1_gate_all_134Te = new TH1D("time_isomer_1_gate_all_134Te","time_isomer_1_gate_all_134Te",2000,-1000,3000);
time_isomer_1_gate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_1_gate_all_134Te->GetYaxis()->SetTitle("Counts");

//Gate on isomer_2
TH1D *time_isomer_2_gate_134Te = new TH1D("time_isomer_2_gate_134Te","time_isomer_2_gate_134Te",2000,-1000,3000);
time_isomer_2_gate_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_2_gate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_2_gate_bg_134Te = new TH1D("time_isomer_2_gate_bg_134Te","time_isomer_2_gate_bg_134Te",2000,-1000,3000);
time_isomer_2_gate_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_2_gate_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_2_gate_all_134Te = new TH1D("time_isomer_2_gate_all_134Te","time_isomer_2_gate_all_134Te",2000,-1000,3000);
time_isomer_2_gate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_2_gate_all_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegates, gate on 297
TH1D *time_isomer_doublegate_1_134Te = new TH1D("time_isomer_doublegate_1_134Te","time_isomer_doublegate_1_134Te",1000,-1000,3000);
time_isomer_doublegate_1_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_all_134Te = new TH1D("time_isomer_doublegate_1_all_134Te","time_isomer_doublegate_1_all_134Te",1000,-1000,3000);
time_isomer_doublegate_1_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_134Te = new TH1D("time_isomer_doublegate_1_bg_134Te","time_isomer_doublegate_1_bg_134Te",1000,-1000,3000);
time_isomer_doublegate_1_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_ridge_134Te = new TH1D("time_isomer_doublegate_1_bg_ridge_134Te","time_isomer_doublegate_1_bg_ridge_134Te",1000,-1000,3000);
time_isomer_doublegate_1_bg_ridge_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_random_134Te = new TH1D("time_isomer_doublegate_1_bg_random_134Te","time_isomer_doublegate_1_bg_random_134Te",1000,-1000,3000);
time_isomer_doublegate_1_bg_random_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_bg_random_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_dt70_134Te = new TH1D("time_isomer_doublegate_1_dt70_134Te","time_isomer_doublegate_1_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_1_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_dt70_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_all_dt70_134Te = new TH1D("time_isomer_doublegate_1_all_dt70_134Te","time_isomer_doublegate_1_all_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_1_all_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_all_dt70_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_dt70_134Te = new TH1D("time_isomer_doublegate_1_bg_dt70_134Te","time_isomer_doublegate_1_bg_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_1_bg_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_1_bg_dt70_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegates, gate on 1279
TH1D *time_isomer_doublegate_2_134Te = new TH1D("time_isomer_doublegate_2_134Te","time_isomer_doublegate_2_134Te",1000,-1000,3000);
time_isomer_doublegate_2_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_all_134Te = new TH1D("time_isomer_doublegate_2_all_134Te","time_isomer_doublegate_2_all_134Te",1000,-1000,3000);
time_isomer_doublegate_2_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_134Te = new TH1D("time_isomer_doublegate_2_bg_134Te","time_isomer_doublegate_2_bg_134Te",1000,-1000,3000);
time_isomer_doublegate_2_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_ridge_134Te = new TH1D("time_isomer_doublegate_2_bg_ridge_134Te","time_isomer_doublegate_2_bg_ridge_134Te",1000,-1000,3000);
time_isomer_doublegate_2_bg_ridge_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_random_134Te = new TH1D("time_isomer_doublegate_2_bg_random_134Te","time_isomer_doublegate_2_bg_random_134Te",1000,-1000,3000);
time_isomer_doublegate_2_bg_random_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_random_134Te->GetYaxis()->SetTitle("Counts");


TH1D *time_isomer_doublegate_2_dt70_134Te = new TH1D("time_isomer_doublegate_2_dt70_134Te","time_isomer_doublegate_2_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_2_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_dt70_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_all_dt70_134Te = new TH1D("time_isomer_doublegate_2_all_dt70_134Te","time_isomer_doublegate_2_all_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_2_all_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_all_dt70_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_dt70_134Te = new TH1D("time_isomer_doublegate_2_bg_dt70_134Te","time_isomer_doublegate_2_bg_dt70_134Te",1000,-1000,3000);
time_isomer_doublegate_2_bg_dt70_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_dt70_134Te->GetYaxis()->SetTitle("Counts");



//Multiplicity-restricted doublegated spectra
//hit3
TH1D *time_isomer_doublegate_2_hit3_134Te = new TH1D("time_isomer_doublegate_2_hit3_134Te","time_isomer_doublegate_2_hit3_134Te",2000,-1000,3000);
time_isomer_doublegate_2_hit3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_hit3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_all_hit3_134Te = new TH1D("time_isomer_doublegate_2_all_hit3_134Te","time_isomer_doublegate_2_all_hit3_134Te",2000,-1000,3000);
time_isomer_doublegate_2_all_hit3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_all_hit3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_hit3_134Te = new TH1D("time_isomer_doublegate_2_bg_hit3_134Te","time_isomer_doublegate_2_bg_hit3_134Te",2000,-1000,3000);
time_isomer_doublegate_2_bg_hit3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_hit3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_ridge_hit3_134Te = new TH1D("time_isomer_doublegate_2_bg_ridge_hit3_134Te","time_isomer_doublegate_2_bg_ridge_hit3_134Te",2000,-1000,3000);
time_isomer_doublegate_2_bg_ridge_hit3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_ridge_hit3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_random_hit3_134Te = new TH1D("time_isomer_doublegate_2_bg_random_hit3_134Te","time_isomer_doublegate_2_bg_random_hit3_134Te",2000,-1000,3000);
time_isomer_doublegate_2_bg_random_hit3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_2_bg_random_hit3_134Te->GetYaxis()->SetTitle("Counts");


//Energy spectrum of all gammas that arrive in events with both 134Te lines
TH1D *coincident_gammas_doublegated_134Te = new TH1D("coincident_gammas_doublegated_134Te","coincident_gammas_doublegated_134Te",4000,0,4000);
coincident_gammas_doublegated_134Te->GetXaxis()->SetTitle("Energy [keV]");
coincident_gammas_doublegated_134Te->GetYaxis()->SetTitle("Counts");

TH1D *coincident_gammas_doublegated_all_134Te = new TH1D("coincident_gammas_doublegated_all_134Te","coincident_gammas_doublegated_all_134Te",4000,0,4000);
coincident_gammas_doublegated_all_134Te->GetXaxis()->SetTitle("Energy [keV]");
coincident_gammas_doublegated_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *coincident_gammas_doublegated_bg_134Te = new TH1D("coincident_gammas_doublegated_bg_134Te","coincident_gammas_doublegated_bg_134Te",4000,0,4000);
coincident_gammas_doublegated_bg_134Te->GetXaxis()->SetTitle("Energy [keV]");
coincident_gammas_doublegated_bg_134Te->GetYaxis()->SetTitle("Counts");

//////////////////////////////////////////////////////////////
//							 140Xe   						//
//////////////////////////////////////////////////////////////

TH1D *time_gamma_doublegate_1_140Xe = new TH1D("time_gamma_doublegate_1_140Xe","time_gamma_doublegate_1_140Xe",2000,-1000,3000);
time_gamma_doublegate_1_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_1_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_1_all_140Xe = new TH1D("time_gamma_doublegate_1_all_140Xe","time_gamma_doublegate_1_all_140Xe",2000,-1000,3000);
time_gamma_doublegate_1_all_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_1_all_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_1_bg_140Xe = new TH1D("time_gamma_doublegate_1_bg_140Xe","time_gamma_doublegate_1_bg_140Xe",2000,-1000,3000);
time_gamma_doublegate_1_bg_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_1_bg_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_1_bg_ridge_140Xe = new TH1D("time_gamma_doublegate_1_bg_ridge_140Xe","time_gamma_doublegate_1_bg_ridge_140Xe",2000,-1000,3000);
time_gamma_doublegate_1_bg_ridge_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_1_bg_ridge_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_1_bg_random_140Xe = new TH1D("time_gamma_doublegate_1_bg_random_140Xe","time_gamma_doublegate_1_bg_random_140Xe",2000,-1000,3000);
time_gamma_doublegate_1_bg_random_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_1_bg_random_140Xe->GetYaxis()->SetTitle("Counts");



TH1D *time_gamma_doublegate_2_140Xe = new TH1D("time_gamma_doublegate_2_140Xe","time_gamma_doublegate_2_140Xe",2000,-1000,3000);
time_gamma_doublegate_2_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_2_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_2_all_140Xe = new TH1D("time_gamma_doublegate_2_all_140Xe","time_gamma_doublegate_2_all_140Xe",2000,-1000,3000);
time_gamma_doublegate_2_all_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_2_all_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_2_bg_140Xe = new TH1D("time_gamma_doublegate_2_bg_140Xe","time_gamma_doublegate_2_bg_140Xe",2000,-1000,3000);
time_gamma_doublegate_2_bg_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_2_bg_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_2_bg_ridge_140Xe = new TH1D("time_gamma_doublegate_2_bg_ridge_140Xe","time_gamma_doublegate_2_bg_ridge_140Xe",2000,-1000,3000);
time_gamma_doublegate_2_bg_ridge_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_2_bg_ridge_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_2_bg_random_140Xe = new TH1D("time_gamma_doublegate_2_bg_random_140Xe","time_gamma_doublegate_2_bg_random_140Xe",2000,-1000,3000);
time_gamma_doublegate_2_bg_random_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_2_bg_random_140Xe->GetYaxis()->SetTitle("Counts");

//////////////////////////////////////////////////////////////
//							 132Sn   						//
//////////////////////////////////////////////////////////////

TH1D *time_gamma_gate_132Sn = new TH1D("time_gamma_gate_132Sn","time_gamma_gate_132Sn",400,-1000,3000);
time_gamma_gate_132Sn->GetXaxis()->SetTitle("Time [ns]");
time_gamma_gate_132Sn->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_gate_bg_132Sn = new TH1D("time_gamma_gate_bg_132Sn","time_gamma_gate_bg_132Sn",400,-1000,3000);
time_gamma_gate_bg_132Sn->GetXaxis()->SetTitle("Time [ns]");
time_gamma_gate_bg_132Sn->GetYaxis()->SetTitle("Counts");


TH1D *time_gamma_gate_all_132Sn = new TH1D("time_gamma_gate_all_132Sn","time_gamma_gate_all_132Sn",400,-1000,3000);
time_gamma_gate_all_132Sn->GetXaxis()->SetTitle("Time [ns]");
time_gamma_gate_all_132Sn->GetYaxis()->SetTitle("Counts");


//////////////////////////////////////////////////////////////////////////




//Anode values for promt and delayed peaks
TH1D *aval_prompt_134Te = new TH1D("aval_prompt_134Te","aval_prompt_134Te",1000,0,1000);
aval_prompt_134Te->GetXaxis()->SetTitle("Energy [keV]");
aval_prompt_134Te->GetYaxis()->SetTitle("Counts");

TH1D *aval_delayed_134Te = new TH1D("aval_delayed_134Te","aval_delayed_134Te",1000,0,1000);
aval_delayed_134Te->GetXaxis()->SetTitle("Energy [keV]");
aval_delayed_134Te->GetYaxis()->SetTitle("Counts");

