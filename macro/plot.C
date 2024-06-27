void plot()
{
  TFile* file = _file0;
  TTree* tree = (TTree*)file->Get("ana_B4");
  int EventNo;
  double E[1000];
  double Esum;
  tree->SetBranchAddress("EventNo", &EventNo);
  tree->SetBranchAddress("E", E);
  tree->SetBranchAddress("Esum", &Esum);
  TH1* h = new TH1D("h", ";Energy [MeV]", 80, 0, 20);
  TH1* h_count = new TH1D("h_count", ";#Hit Boxes", 10, 0, 10);
  TH1* h_sum = new TH1D("h_sum", ";Energy sum [MeV]", 80, 0, 20);
  double count2 = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    h->Fill(E[332]);
    h_sum->Fill(Esum);
    int count = 0;
    
    for (int n = 0; n < 1000; n++){
      if(E[n] > 0) {
	count++;
      }
    }
    h_count->Fill(count);
    count2 += count;
  }


  double avg = 0;
  avg = count2 / 1000000;
  printf("反応の平均値は%f\n", avg);

  
  TCanvas* c_count = new TCanvas("c_count", "c_count", 640, 640);
  h_count->Draw();
  TCanvas* c_sum = new TCanvas("c_sum", "c_sum", 640, 640);
  h_sum->Draw();
  TCanvas* c1 = new TCanvas("c1", "c1", 640, 640);
  h->Draw();
}
