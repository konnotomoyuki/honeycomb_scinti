void draw()
{
  TFile* file = _file0;
  TTree* tree = (TTree*)file->Get("B4");
  int EventNo, Id;
  double Edep;
  tree->SetBranchAddress("EventNo", &EventNo);
  tree->SetBranchAddress("Id", &Id);
  tree->SetBranchAddress("Edep", &Edep);
  TH1* h = new TH1D("h", ";Energy [MeV]", 80, 0, 4);
  TH1* h_sum = new TH1D("h_sum", ";Energy sum [MeV]", 80, 0, 4);
  TFile* o_file = new TFile(Form("ana_%s", file->GetName()), "recreate");
  TTree* o_tree = new TTree("ana_B4", "ana_B4");
  int o_EventNo = -1;
  double E[1000];
  double Esum = 0;
  o_tree->Branch("EventNo", &o_EventNo, "EventNo/I");
  o_tree->Branch("E", E, "E[1000]/D");
  o_tree->Branch("Esum", &Esum, "Esum/D");
  memset(E, 0, sizeof(double) * 1000);
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (Id == 2 && Edep > 0) {
      h->Fill(Edep);
    }
    if (o_EventNo >= 0 && o_EventNo != EventNo) {
      h_sum->Fill(Esum);
      o_tree->Fill();
      memset(E, 0, sizeof(double) * 1000);
      Esum = 0;
    }
    E[Id] = Edep;
    Esum += Edep;
    o_EventNo = EventNo;
  }
  h_sum->Draw();
  o_file->cd();
  o_tree->Write();
  o_file->Close();
}
