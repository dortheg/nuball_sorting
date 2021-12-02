//All single, clean ge-gammas
TH1D *single_gamma = new TH1D("single_gamma","single_gamma",4000,0,4000);
single_gamma->GetXaxis()->SetTitle("Energy [keV]");
single_gamma->GetYaxis()->SetTitle("Counts");

//All double, clean ge-gammas
TH2F *double_gamma = new TH2F("double_gamma","double_gamma",4000,0,4000, 4000, 0, 4000);
double_gamma->GetXaxis()->SetTitle("Energy [keV]");
double_gamma->GetYaxis()->SetTitle("Energy [keV]");
double_gamma->SetOption("colz");

//Multiplicity distribution
TH1D *mult_distr = new TH1D("mult_distr","mult_distr",15,0,15);
mult_distr->GetXaxis()->SetTitle("Multiplicity");
mult_distr->GetYaxis()->SetTitle("Counts");


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

TH1D *time_isomer_1_gate_bg_new_134Te = new TH1D("time_isomer_1_gate_bg_new_134Te","time_isomer_1_gate_bg_new_134Te",2000,-1000,3000);
time_isomer_1_gate_bg_new_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_1_gate_bg_new_134Te->GetYaxis()->SetTitle("Counts");

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

TH1D *time_isomer_2_gate_bg_new_134Te = new TH1D("time_isomer_2_gate_bg_new_134Te","time_isomer_2_gate_bg_new_134Te",2000,-1000,3000);
time_isomer_2_gate_bg_new_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_2_gate_bg_new_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_2_gate_all_134Te = new TH1D("time_isomer_2_gate_all_134Te","time_isomer_2_gate_all_134Te",2000,-1000,3000);
time_isomer_2_gate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_2_gate_all_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegates
TH1D *time_isomer_doublegate_134Te = new TH1D("time_isomer_doublegate_134Te","time_isomer_doublegate_134Te",2000,-1000,3000);
time_isomer_doublegate_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_134Te = new TH1D("time_isomer_doublegate_all_134Te","time_isomer_doublegate_all_134Te",2000,-1000,3000);
time_isomer_doublegate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_134Te = new TH1D("time_isomer_doublegate_bg_134Te","time_isomer_doublegate_bg_134Te",2000,-1000,3000);
time_isomer_doublegate_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_134Te = new TH1D("time_isomer_doublegate_bg_ridge_134Te","time_isomer_doublegate_bg_ridge_134Te",2000,-1000,3000);
time_isomer_doublegate_bg_ridge_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_134Te = new TH1D("time_isomer_doublegate_bg_random_134Te","time_isomer_doublegate_bg_random_134Te",2000,-1000,3000);
time_isomer_doublegate_bg_random_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//Multiplicity-restricted doublegated spectra
//mult3
TH1D *time_isomer_doublegate_mult3_134Te = new TH1D("time_isomer_doublegate_mult3_134Te","time_isomer_doublegate_mult3_134Te",2000,-1000,3000);
time_isomer_doublegate_mult3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_mult3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_mult3_134Te = new TH1D("time_isomer_doublegate_all_mult3_134Te","time_isomer_doublegate_all_mult3_134Te",2000,-1000,3000);
time_isomer_doublegate_all_mult3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_all_mult3_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_mult3_134Te = new TH1D("time_isomer_doublegate_bg_mult3_134Te","time_isomer_doublegate_bg_mult3_134Te",2000,-1000,3000);
time_isomer_doublegate_bg_mult3_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_bg_mult3_134Te->GetYaxis()->SetTitle("Counts");

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

TH1D *time_gamma_doublegate_140Xe = new TH1D("time_gamma_doublegate_140Xe","time_gamma_doublegate_140Xe",2000,-1000,3000);
time_gamma_doublegate_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_all_140Xe = new TH1D("time_gamma_doublegate_all_140Xe","time_gamma_doublegate_all_140Xe",2000,-1000,3000);
time_gamma_doublegate_all_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_all_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_bg_140Xe = new TH1D("time_gamma_doublegate_bg_140Xe","time_gamma_doublegate_bg_140Xe",2000,-1000,3000);
time_gamma_doublegate_bg_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_bg_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_bg_ridge_140Xe = new TH1D("time_gamma_doublegate_bg_ridge_140Xe","time_gamma_doublegate_bg_ridge_140Xe",2000,-1000,3000);
time_gamma_doublegate_bg_ridge_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_bg_ridge_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_gamma_doublegate_bg_random_140Xe = new TH1D("time_gamma_doublegate_bg_random_140Xe","time_gamma_doublegate_bg_random_140Xe",2000,-1000,3000);
time_gamma_doublegate_bg_random_140Xe->GetXaxis()->SetTitle("Time [ns]");
time_gamma_doublegate_bg_random_140Xe->GetYaxis()->SetTitle("Counts");

//////////////////////////////////////////////////////////////////////////

//Anode values for promt and delayed peaks
TH1D *aval_prompt_134Te = new TH1D("aval_prompt_134Te","aval_prompt_134Te",1000,0,1000);
aval_prompt_134Te->GetXaxis()->SetTitle("Energy [keV]");
aval_prompt_134Te->GetYaxis()->SetTitle("Counts");

TH1D *aval_delayed_134Te = new TH1D("aval_delayed_134Te","aval_delayed_134Te",1000,0,1000);
aval_delayed_134Te->GetXaxis()->SetTitle("Energy [keV]");
aval_delayed_134Te->GetYaxis()->SetTitle("Counts");

