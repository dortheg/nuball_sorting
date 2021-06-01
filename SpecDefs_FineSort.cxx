//All single, clean ge-gammas
TH1D *single_gamma = new TH1D("single_gamma","single_gamma",4000,0,4000);
single_gamma->GetXaxis()->SetTitle("Energy [keV]");
single_gamma->GetYaxis()->SetTitle("Counts");

//All double, clean ge-gammas
TH2F *double_gamma = new TH2F("double_gamma","double_gamma",4000,0,4000, 4000, 0, 4000);
double_gamma->GetXaxis()->SetTitle("Energy [keV]");
double_gamma->GetYaxis()->SetTitle("Energy [keV]");
double_gamma->SetOption("colz");

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

TH1D *time_isomer_doublegate_134Te = new TH1D("time_isomer_doublegate_134Te","time_isomer_doublegate_134Te",2000,-1000,3000);
time_isomer_doublegate_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_134Te = new TH1D("time_isomer_doublegate_bg_134Te","time_isomer_doublegate_bg_134Te",2000,-1000,3000);
time_isomer_doublegate_bg_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_134Te = new TH1D("time_isomer_doublegate_all_134Te","time_isomer_doublegate_all_134Te",2000,-1000,3000);
time_isomer_doublegate_all_134Te->GetXaxis()->SetTitle("Time [ns]");
time_isomer_doublegate_all_134Te->GetYaxis()->SetTitle("Counts");

//Energy spectrum in coincidence with both 134Te lines
TH1D *coincident_gammas_doublegated_all_134Te = new TH1D("coincident_gammas_doublegated_all_134Te","coincident_gammas_doublegated_all_134Te",4000,0,4000);
coincident_gammas_doublegated_all_134Te->GetXaxis()->SetTitle("Energy [keV]");
coincident_gammas_doublegated_all_134Te->GetYaxis()->SetTitle("Counts");

