//Singlegate 3n partner
TH1D *singlegate_3n = new TH1D("singlegate_3n","singlegate_3n",2048,0,2048);
singlegate_3n->GetXaxis()->SetTitle("Energy [keV]");
singlegate_3n->GetYaxis()->SetTitle("Counts");

//Singlegate 3n partner all
TH1D *singlegate_3n_all = new TH1D("singlegate_3n_all","singlegate_3n_all",2048,0,2048);
singlegate_3n_all->GetXaxis()->SetTitle("Energy [keV]");
singlegate_3n_all->GetYaxis()->SetTitle("Counts");

//Singlegate 3n partner bg
TH1D *singlegate_3n_bg = new TH1D("singlegate_3n_bg","singlegate_3n_bg",2048,0,2048);
singlegate_3n_bg->GetXaxis()->SetTitle("Energy [keV]");
singlegate_3n_bg->GetYaxis()->SetTitle("Counts");


//134Te doublegate
TH1D *time_isomer_doublegate_134Te = new TH1D("time_isomer_doublegate_134Te","time_isomer_doublegate_134Te",350,0,350);
time_isomer_doublegate_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_134Te = new TH1D("time_isomer_doublegate_all_134Te","time_isomer_doublegate_all_134Te",350,0,350);
time_isomer_doublegate_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_134Te = new TH1D("time_isomer_doublegate_bg_134Te","time_isomer_doublegate_bg_134Te",350,0,350);
time_isomer_doublegate_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_134Te = new TH1D("time_isomer_doublegate_bg_ridge_134Te","time_isomer_doublegate_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_134Te = new TH1D("time_isomer_doublegate_bg_random_134Te","time_isomer_doublegate_bg_random_134Te",350,0,350);
time_isomer_doublegate_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_isomer_3n
TH1D *time_isomer_doublegate_isomer_3n_134Te = new TH1D("time_isomer_doublegate_isomer_3n_134Te","time_isomer_doublegate_isomer_3n_134Te",350,0,350);
time_isomer_doublegate_isomer_3n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_isomer_3n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_isomer_3n_all_134Te = new TH1D("time_isomer_doublegate_isomer_3n_all_134Te","time_isomer_doublegate_isomer_3n_all_134Te",350,0,350);
time_isomer_doublegate_isomer_3n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_isomer_3n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_isomer_3n_bg_134Te = new TH1D("time_isomer_doublegate_isomer_3n_bg_134Te","time_isomer_doublegate_isomer_3n_bg_134Te",350,0,350);
time_isomer_doublegate_isomer_3n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_isomer_3n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_isomer_3n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_isomer_3n_bg_ridge_134Te","time_isomer_doublegate_isomer_3n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_isomer_3n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_isomer_3n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_isomer_3n_bg_random_134Te = new TH1D("time_isomer_doublegate_isomer_3n_bg_random_134Te","time_isomer_doublegate_isomer_3n_bg_random_134Te",350,0,350);
time_isomer_doublegate_isomer_3n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_isomer_3n_bg_random_134Te->GetYaxis()->SetTitle("Counts");



//134Te doublegate_1
TH1D *time_isomer_doublegate_1_134Te = new TH1D("time_isomer_doublegate_1_134Te","time_isomer_doublegate_1_134Te",350,0,350);
time_isomer_doublegate_1_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_all_134Te = new TH1D("time_isomer_doublegate_1_all_134Te","time_isomer_doublegate_1_all_134Te",350,0,350);
time_isomer_doublegate_1_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_134Te = new TH1D("time_isomer_doublegate_1_bg_134Te","time_isomer_doublegate_1_bg_134Te",350,0,350);
time_isomer_doublegate_1_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_ridge_134Te = new TH1D("time_isomer_doublegate_1_bg_ridge_134Te","time_isomer_doublegate_1_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_1_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_random_134Te = new TH1D("time_isomer_doublegate_1_bg_random_134Te","time_isomer_doublegate_1_bg_random_134Te",350,0,350);
time_isomer_doublegate_1_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_2
TH1D *time_isomer_doublegate_2_134Te = new TH1D("time_isomer_doublegate_2_134Te","time_isomer_doublegate_2_134Te",350,0,350);
time_isomer_doublegate_2_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_all_134Te = new TH1D("time_isomer_doublegate_2_all_134Te","time_isomer_doublegate_2_all_134Te",350,0,350);
time_isomer_doublegate_2_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_134Te = new TH1D("time_isomer_doublegate_2_bg_134Te","time_isomer_doublegate_2_bg_134Te",350,0,350);
time_isomer_doublegate_2_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_ridge_134Te = new TH1D("time_isomer_doublegate_2_bg_ridge_134Te","time_isomer_doublegate_2_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_2_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_random_134Te = new TH1D("time_isomer_doublegate_2_bg_random_134Te","time_isomer_doublegate_2_bg_random_134Te",350,0,350);
time_isomer_doublegate_2_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_1n
TH1D *time_isomer_doublegate_1n_134Te = new TH1D("time_isomer_doublegate_1n_134Te","time_isomer_doublegate_1n_134Te",350,0,350);
time_isomer_doublegate_1n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1n_all_134Te = new TH1D("time_isomer_doublegate_1n_all_134Te","time_isomer_doublegate_1n_all_134Te",350,0,350);
time_isomer_doublegate_1n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1n_bg_134Te = new TH1D("time_isomer_doublegate_1n_bg_134Te","time_isomer_doublegate_1n_bg_134Te",350,0,350);
time_isomer_doublegate_1n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_1n_bg_ridge_134Te","time_isomer_doublegate_1n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_1n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1n_bg_random_134Te = new TH1D("time_isomer_doublegate_1n_bg_random_134Te","time_isomer_doublegate_1n_bg_random_134Te",350,0,350);
time_isomer_doublegate_1n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1n_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_2n
TH1D *time_isomer_doublegate_2n_134Te = new TH1D("time_isomer_doublegate_2n_134Te","time_isomer_doublegate_2n_134Te",350,0,350);
time_isomer_doublegate_2n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2n_all_134Te = new TH1D("time_isomer_doublegate_2n_all_134Te","time_isomer_doublegate_2n_all_134Te",350,0,350);
time_isomer_doublegate_2n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2n_bg_134Te = new TH1D("time_isomer_doublegate_2n_bg_134Te","time_isomer_doublegate_2n_bg_134Te",350,0,350);
time_isomer_doublegate_2n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_2n_bg_ridge_134Te","time_isomer_doublegate_2n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_2n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2n_bg_random_134Te = new TH1D("time_isomer_doublegate_2n_bg_random_134Te","time_isomer_doublegate_2n_bg_random_134Te",350,0,350);
time_isomer_doublegate_2n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2n_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_3n
TH1D *time_isomer_doublegate_3n_134Te = new TH1D("time_isomer_doublegate_3n_134Te","time_isomer_doublegate_3n_134Te",350,0,350);
time_isomer_doublegate_3n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_all_134Te = new TH1D("time_isomer_doublegate_3n_all_134Te","time_isomer_doublegate_3n_all_134Te",350,0,350);
time_isomer_doublegate_3n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_bg_134Te = new TH1D("time_isomer_doublegate_3n_bg_134Te","time_isomer_doublegate_3n_bg_134Te",350,0,350);
time_isomer_doublegate_3n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_3n_bg_ridge_134Te","time_isomer_doublegate_3n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_3n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_bg_random_134Te = new TH1D("time_isomer_doublegate_3n_bg_random_134Te","time_isomer_doublegate_3n_bg_random_134Te",350,0,350);
time_isomer_doublegate_3n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_3n_2plus
TH1D *time_isomer_doublegate_3n_2plus_134Te = new TH1D("time_isomer_doublegate_3n_2plus_134Te","time_isomer_doublegate_3n_2plus_134Te",350,0,350);
time_isomer_doublegate_3n_2plus_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_2plus_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_2plus_all_134Te = new TH1D("time_isomer_doublegate_3n_2plus_all_134Te","time_isomer_doublegate_3n_2plus_all_134Te",350,0,350);
time_isomer_doublegate_3n_2plus_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_2plus_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_2plus_bg_134Te = new TH1D("time_isomer_doublegate_3n_2plus_bg_134Te","time_isomer_doublegate_3n_2plus_bg_134Te",350,0,350);
time_isomer_doublegate_3n_2plus_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_2plus_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_2plus_bg_ridge_134Te = new TH1D("time_isomer_doublegate_3n_2plus_bg_ridge_134Te","time_isomer_doublegate_3n_2plus_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_3n_2plus_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_2plus_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_2plus_bg_random_134Te = new TH1D("time_isomer_doublegate_3n_2plus_bg_random_134Te","time_isomer_doublegate_3n_2plus_bg_random_134Te",350,0,350);
time_isomer_doublegate_3n_2plus_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_2plus_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_3n_4plus
TH1D *time_isomer_doublegate_3n_4plus_134Te = new TH1D("time_isomer_doublegate_3n_4plus_134Te","time_isomer_doublegate_3n_4plus_134Te",350,0,350);
time_isomer_doublegate_3n_4plus_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_4plus_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_4plus_all_134Te = new TH1D("time_isomer_doublegate_3n_4plus_all_134Te","time_isomer_doublegate_3n_4plus_all_134Te",350,0,350);
time_isomer_doublegate_3n_4plus_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_4plus_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_4plus_bg_134Te = new TH1D("time_isomer_doublegate_3n_4plus_bg_134Te","time_isomer_doublegate_3n_4plus_bg_134Te",350,0,350);
time_isomer_doublegate_3n_4plus_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_4plus_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_4plus_bg_ridge_134Te = new TH1D("time_isomer_doublegate_3n_4plus_bg_ridge_134Te","time_isomer_doublegate_3n_4plus_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_3n_4plus_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_4plus_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_4plus_bg_random_134Te = new TH1D("time_isomer_doublegate_3n_4plus_bg_random_134Te","time_isomer_doublegate_3n_4plus_bg_random_134Te",350,0,350);
time_isomer_doublegate_3n_4plus_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_4plus_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_3n_6plus
TH1D *time_isomer_doublegate_3n_6plus_134Te = new TH1D("time_isomer_doublegate_3n_6plus_134Te","time_isomer_doublegate_3n_6plus_134Te",350,0,350);
time_isomer_doublegate_3n_6plus_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_6plus_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_6plus_all_134Te = new TH1D("time_isomer_doublegate_3n_6plus_all_134Te","time_isomer_doublegate_3n_6plus_all_134Te",350,0,350);
time_isomer_doublegate_3n_6plus_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_6plus_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_6plus_bg_134Te = new TH1D("time_isomer_doublegate_3n_6plus_bg_134Te","time_isomer_doublegate_3n_6plus_bg_134Te",350,0,350);
time_isomer_doublegate_3n_6plus_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_6plus_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_6plus_bg_ridge_134Te = new TH1D("time_isomer_doublegate_3n_6plus_bg_ridge_134Te","time_isomer_doublegate_3n_6plus_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_3n_6plus_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_6plus_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_6plus_bg_random_134Te = new TH1D("time_isomer_doublegate_3n_6plus_bg_random_134Te","time_isomer_doublegate_3n_6plus_bg_random_134Te",350,0,350);
time_isomer_doublegate_3n_6plus_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_6plus_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_3n_8plus
TH1D *time_isomer_doublegate_3n_8plus_134Te = new TH1D("time_isomer_doublegate_3n_8plus_134Te","time_isomer_doublegate_3n_8plus_134Te",350,0,350);
time_isomer_doublegate_3n_8plus_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_8plus_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_8plus_all_134Te = new TH1D("time_isomer_doublegate_3n_8plus_all_134Te","time_isomer_doublegate_3n_8plus_all_134Te",350,0,350);
time_isomer_doublegate_3n_8plus_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_8plus_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_8plus_bg_134Te = new TH1D("time_isomer_doublegate_3n_8plus_bg_134Te","time_isomer_doublegate_3n_8plus_bg_134Te",350,0,350);
time_isomer_doublegate_3n_8plus_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_8plus_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_8plus_bg_ridge_134Te = new TH1D("time_isomer_doublegate_3n_8plus_bg_ridge_134Te","time_isomer_doublegate_3n_8plus_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_3n_8plus_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_8plus_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_3n_8plus_bg_random_134Te = new TH1D("time_isomer_doublegate_3n_8plus_bg_random_134Te","time_isomer_doublegate_3n_8plus_bg_random_134Te",350,0,350);
time_isomer_doublegate_3n_8plus_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_3n_8plus_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_4n
TH1D *time_isomer_doublegate_4n_134Te = new TH1D("time_isomer_doublegate_4n_134Te","time_isomer_doublegate_4n_134Te",350,0,350);
time_isomer_doublegate_4n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_4n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_4n_all_134Te = new TH1D("time_isomer_doublegate_4n_all_134Te","time_isomer_doublegate_4n_all_134Te",350,0,350);
time_isomer_doublegate_4n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_4n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_4n_bg_134Te = new TH1D("time_isomer_doublegate_4n_bg_134Te","time_isomer_doublegate_4n_bg_134Te",350,0,350);
time_isomer_doublegate_4n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_4n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_4n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_4n_bg_ridge_134Te","time_isomer_doublegate_4n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_4n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_4n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_4n_bg_random_134Te = new TH1D("time_isomer_doublegate_4n_bg_random_134Te","time_isomer_doublegate_4n_bg_random_134Te",350,0,350);
time_isomer_doublegate_4n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_4n_bg_random_134Te->GetYaxis()->SetTitle("Counts");


//134Te doublegate_5n
TH1D *time_isomer_doublegate_5n_134Te = new TH1D("time_isomer_doublegate_5n_134Te","time_isomer_doublegate_5n_134Te",350,0,350);
time_isomer_doublegate_5n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_5n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_5n_all_134Te = new TH1D("time_isomer_doublegate_5n_all_134Te","time_isomer_doublegate_5n_all_134Te",350,0,350);
time_isomer_doublegate_5n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_5n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_5n_bg_134Te = new TH1D("time_isomer_doublegate_5n_bg_134Te","time_isomer_doublegate_5n_bg_134Te",350,0,350);
time_isomer_doublegate_5n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_5n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_5n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_5n_bg_ridge_134Te","time_isomer_doublegate_5n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_5n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_5n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_5n_bg_random_134Te = new TH1D("time_isomer_doublegate_5n_bg_random_134Te","time_isomer_doublegate_5n_bg_random_134Te",350,0,350);
time_isomer_doublegate_5n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_5n_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_6n
TH1D *time_isomer_doublegate_6n_134Te = new TH1D("time_isomer_doublegate_6n_134Te","time_isomer_doublegate_6n_134Te",350,0,350);
time_isomer_doublegate_6n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_6n_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_6n_all_134Te = new TH1D("time_isomer_doublegate_6n_all_134Te","time_isomer_doublegate_6n_all_134Te",350,0,350);
time_isomer_doublegate_6n_all_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_6n_all_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_6n_bg_134Te = new TH1D("time_isomer_doublegate_6n_bg_134Te","time_isomer_doublegate_6n_bg_134Te",350,0,350);
time_isomer_doublegate_6n_bg_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_6n_bg_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_6n_bg_ridge_134Te = new TH1D("time_isomer_doublegate_6n_bg_ridge_134Te","time_isomer_doublegate_6n_bg_ridge_134Te",350,0,350);
time_isomer_doublegate_6n_bg_ridge_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_6n_bg_ridge_134Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_6n_bg_random_134Te = new TH1D("time_isomer_doublegate_6n_bg_random_134Te","time_isomer_doublegate_6n_bg_random_134Te",350,0,350);
time_isomer_doublegate_6n_bg_random_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_6n_bg_random_134Te->GetYaxis()->SetTitle("Counts");

//134Te doublegate_any_n
TH1D *time_isomer_doublegate_any_n_134Te = new TH1D("time_isomer_doublegate_any_n_134Te","time_isomer_doublegate_any_n_134Te",350,0,350);
time_isomer_doublegate_any_n_134Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_any_n_134Te->GetYaxis()->SetTitle("Counts");


//135Te doublegate
TH1D *time_isomer_doublegate_135Te = new TH1D("time_isomer_doublegate_135Te","time_isomer_doublegate_135Te",350,0,350);
time_isomer_doublegate_135Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_135Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_135Te = new TH1D("time_isomer_doublegate_all_135Te","time_isomer_doublegate_all_135Te",350,0,350);
time_isomer_doublegate_all_135Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_135Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_135Te = new TH1D("time_isomer_doublegate_bg_135Te","time_isomer_doublegate_bg_135Te",350,0,350);
time_isomer_doublegate_bg_135Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_135Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_135Te = new TH1D("time_isomer_doublegate_bg_ridge_135Te","time_isomer_doublegate_bg_ridge_135Te",350,0,350);
time_isomer_doublegate_bg_ridge_135Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_135Te->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_135Te = new TH1D("time_isomer_doublegate_bg_random_135Te","time_isomer_doublegate_bg_random_135Te",350,0,350);
time_isomer_doublegate_bg_random_135Te->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_135Te->GetYaxis()->SetTitle("Counts");


//93Rb doublegate_1
TH1D *time_isomer_doublegate_1_93Rb = new TH1D("time_isomer_doublegate_1_93Rb","time_isomer_doublegate_1_93Rb",350,0,350);
time_isomer_doublegate_1_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_all_93Rb = new TH1D("time_isomer_doublegate_1_all_93Rb","time_isomer_doublegate_1_all_93Rb",350,0,350);
time_isomer_doublegate_1_all_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_all_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_93Rb = new TH1D("time_isomer_doublegate_1_bg_93Rb","time_isomer_doublegate_1_bg_93Rb",350,0,350);
time_isomer_doublegate_1_bg_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_ridge_93Rb = new TH1D("time_isomer_doublegate_1_bg_ridge_93Rb","time_isomer_doublegate_1_bg_ridge_93Rb",350,0,350);
time_isomer_doublegate_1_bg_ridge_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_ridge_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_1_bg_random_93Rb = new TH1D("time_isomer_doublegate_1_bg_random_93Rb","time_isomer_doublegate_1_bg_random_93Rb",350,0,350);
time_isomer_doublegate_1_bg_random_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_1_bg_random_93Rb->GetYaxis()->SetTitle("Counts");

//93Rb doublegate_2
TH1D *time_isomer_doublegate_2_93Rb = new TH1D("time_isomer_doublegate_2_93Rb","time_isomer_doublegate_2_93Rb",350,0,350);
time_isomer_doublegate_2_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_all_93Rb = new TH1D("time_isomer_doublegate_2_all_93Rb","time_isomer_doublegate_2_all_93Rb",350,0,350);
time_isomer_doublegate_2_all_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_all_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_93Rb = new TH1D("time_isomer_doublegate_2_bg_93Rb","time_isomer_doublegate_2_bg_93Rb",350,0,350);
time_isomer_doublegate_2_bg_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_ridge_93Rb = new TH1D("time_isomer_doublegate_2_bg_ridge_93Rb","time_isomer_doublegate_2_bg_ridge_93Rb",350,0,350);
time_isomer_doublegate_2_bg_ridge_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_ridge_93Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_2_bg_random_93Rb = new TH1D("time_isomer_doublegate_2_bg_random_93Rb","time_isomer_doublegate_2_bg_random_93Rb",350,0,350);
time_isomer_doublegate_2_bg_random_93Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_2_bg_random_93Rb->GetYaxis()->SetTitle("Counts");

//94Rb doublegate
TH1D *time_isomer_doublegate_94Rb = new TH1D("time_isomer_doublegate_94Rb","time_isomer_doublegate_94Rb",350,0,350);
time_isomer_doublegate_94Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_94Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_94Rb = new TH1D("time_isomer_doublegate_all_94Rb","time_isomer_doublegate_all_94Rb",350,0,350);
time_isomer_doublegate_all_94Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_94Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_94Rb = new TH1D("time_isomer_doublegate_bg_94Rb","time_isomer_doublegate_bg_94Rb",350,0,350);
time_isomer_doublegate_bg_94Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_94Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_94Rb = new TH1D("time_isomer_doublegate_bg_ridge_94Rb","time_isomer_doublegate_bg_ridge_94Rb",350,0,350);
time_isomer_doublegate_bg_ridge_94Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_94Rb->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_94Rb = new TH1D("time_isomer_doublegate_bg_random_94Rb","time_isomer_doublegate_bg_random_94Rb",350,0,350);
time_isomer_doublegate_bg_random_94Rb->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_94Rb->GetYaxis()->SetTitle("Counts");

//140Xe doublegate
TH1D *time_isomer_doublegate_140Xe = new TH1D("time_isomer_doublegate_140Xe","time_isomer_doublegate_140Xe",350,0,350);
time_isomer_doublegate_140Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_140Xe = new TH1D("time_isomer_doublegate_all_140Xe","time_isomer_doublegate_all_140Xe",350,0,350);
time_isomer_doublegate_all_140Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_140Xe = new TH1D("time_isomer_doublegate_bg_140Xe","time_isomer_doublegate_bg_140Xe",350,0,350);
time_isomer_doublegate_bg_140Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_140Xe = new TH1D("time_isomer_doublegate_bg_ridge_140Xe","time_isomer_doublegate_bg_ridge_140Xe",350,0,350);
time_isomer_doublegate_bg_ridge_140Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_140Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_140Xe = new TH1D("time_isomer_doublegate_bg_random_140Xe","time_isomer_doublegate_bg_random_140Xe",350,0,350);
time_isomer_doublegate_bg_random_140Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_140Xe->GetYaxis()->SetTitle("Counts");


//138Xe doublegate
TH1D *time_isomer_doublegate_138Xe = new TH1D("time_isomer_doublegate_138Xe","time_isomer_doublegate_138Xe",350,0,350);
time_isomer_doublegate_138Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_138Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_138Xe = new TH1D("time_isomer_doublegate_all_138Xe","time_isomer_doublegate_all_138Xe",350,0,350);
time_isomer_doublegate_all_138Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_138Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_138Xe = new TH1D("time_isomer_doublegate_bg_138Xe","time_isomer_doublegate_bg_138Xe",350,0,350);
time_isomer_doublegate_bg_138Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_138Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_138Xe = new TH1D("time_isomer_doublegate_bg_ridge_138Xe","time_isomer_doublegate_bg_ridge_138Xe",350,0,350);
time_isomer_doublegate_bg_ridge_138Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_138Xe->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_138Xe = new TH1D("time_isomer_doublegate_bg_random_138Xe","time_isomer_doublegate_bg_random_138Xe",350,0,350);
time_isomer_doublegate_bg_random_138Xe->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_138Xe->GetYaxis()->SetTitle("Counts");


//92Sr doublegate
TH1D *time_isomer_doublegate_92Sr = new TH1D("time_isomer_doublegate_92Sr","time_isomer_doublegate_92Sr",350,0,350);
time_isomer_doublegate_92Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_92Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_92Sr = new TH1D("time_isomer_doublegate_all_92Sr","time_isomer_doublegate_all_92Sr",350,0,350);
time_isomer_doublegate_all_92Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_92Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_92Sr = new TH1D("time_isomer_doublegate_bg_92Sr","time_isomer_doublegate_bg_92Sr",350,0,350);
time_isomer_doublegate_bg_92Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_92Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_92Sr = new TH1D("time_isomer_doublegate_bg_ridge_92Sr","time_isomer_doublegate_bg_ridge_92Sr",350,0,350);
time_isomer_doublegate_bg_ridge_92Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_92Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_92Sr = new TH1D("time_isomer_doublegate_bg_random_92Sr","time_isomer_doublegate_bg_random_92Sr",350,0,350);
time_isomer_doublegate_bg_random_92Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_92Sr->GetYaxis()->SetTitle("Counts");



//94Sr doublegate
TH1D *time_isomer_doublegate_94Sr = new TH1D("time_isomer_doublegate_94Sr","time_isomer_doublegate_94Sr",350,0,350);
time_isomer_doublegate_94Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_94Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_all_94Sr = new TH1D("time_isomer_doublegate_all_94Sr","time_isomer_doublegate_all_94Sr",350,0,350);
time_isomer_doublegate_all_94Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_all_94Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_94Sr = new TH1D("time_isomer_doublegate_bg_94Sr","time_isomer_doublegate_bg_94Sr",350,0,350);
time_isomer_doublegate_bg_94Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_94Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_ridge_94Sr = new TH1D("time_isomer_doublegate_bg_ridge_94Sr","time_isomer_doublegate_bg_ridge_94Sr",350,0,350);
time_isomer_doublegate_bg_ridge_94Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_ridge_94Sr->GetYaxis()->SetTitle("Counts");

TH1D *time_isomer_doublegate_bg_random_94Sr = new TH1D("time_isomer_doublegate_bg_random_94Sr","time_isomer_doublegate_bg_random_94Sr",350,0,350);
time_isomer_doublegate_bg_random_94Sr->GetXaxis()->SetTitle("Time [2 ns]");
time_isomer_doublegate_bg_random_94Sr->GetYaxis()->SetTitle("Counts");

//All double, clean ge-gammas
TH2F *double_gamma = new TH2F("double_gamma","double_gamma",2048,0,2048, 2048, 0, 2048);
double_gamma->GetXaxis()->SetTitle("Energy [keV]");
double_gamma->GetYaxis()->SetTitle("Energy [keV]");
double_gamma->SetOption("colz");

