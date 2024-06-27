void ana(double x0=155, double y0=5, double z0=5)
{
  TFile* file = _file0;
  TTree* tree = (TTree*)file->Get("B4");
  int EventNo, Id;
  double Edep;
  double x, y, z;
  tree->SetBranchAddress("EventNo", &EventNo);
  tree->SetBranchAddress("Id", &Id);
  tree->SetBranchAddress("Edep", &Edep);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  TH1* h_Esum = new TH1D("h_Esum", ";Energy [MeV]", 100, 0, 3);
  TH1* h_Emax = new TH1D("h_Emax", ";Energy [MeV]", 100, 0, 1);
  TH1* h_Ef = new TH1D("h_Ef", ";Energy [MeV]", 100, 0, 1);
  TH1* h_Eb = new TH1D("h_Eb", ";Energy [MeV]", 100, 0, 1);
  TH1* h_cos_theta = new TH1D("h_cos_theta", ";cos theta", 100, -1, 1);
  TH2* h_Ef_Eb = new TH2D("h_Ef_Eb", ";Energy [MeV];Energy [MeV]", 100, 1, 0, 100, 0, 1);
  TH2* h_Esum_Eb = new TH2D("h_Esum_Eb", ";Energy [MeV];Energy [MeV]", 100, 1, 0, 100, 0, 1);
  int o_EventNo = -1;
  double Esum = 0, Emax = 0;
  struct vec4 {
    TVector3 v;
    double E;
  };
  std::vector<vec4> vs;
  vec4 vmax;
  double Ef = 0;
  double Eb = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (o_EventNo >= 0 && o_EventNo != EventNo) {
      h_Esum->Fill(Esum);
      h_Emax->Fill(Emax);
      TVector3 vf, vb;
      for (auto v : vs) {
	if (vmax.v.Dot(v.v) > 0) {
	  vf += v.E * v.v.Unit(); 
	  Ef += v.E;
	} else {
	  vb += v.E * v.v.Unit(); 
	  Eb += v.E;
	}
      }
      double cos_theta = 2;
      if (Ef > 0 && Eb > 0) {
	cos_theta = vf.Unit().Dot(vb.Unit());
	h_cos_theta->Fill(cos_theta);
      }
      if (Ef > 0.1) {
	if (cos_theta < -0.6) {
	  h_Ef->Fill(Ef);
	  if (Eb > 0.1)
	    h_Eb->Fill(Eb);
	}
      }
      h_Ef_Eb->Fill(Ef, Eb);
      h_Esum_Eb->Fill(Esum, Eb);
      Ef = Eb = 0;
      Esum = 0;
      Emax = 0;
      vs = std::vector<vec4>();
    }
    vec4 v;
    v.v = TVector3(x - x0, y - y0, z - z0);
    v.E = Edep;
    if ((abs(x - x0) < 2 && fabs(y - y0) < 2 && fabs(z - z0) < 2)) {
      Esum += Edep;
    } else {
      vs.push_back(v);
      if (Id > 0) {
	if (Emax < Edep) {
	  Emax = Edep;
	  vmax = v;
	}
      }
    }
     o_EventNo = EventNo;
  }
  TCanvas* c_Esum = new TCanvas("c_Esum", "c_Esum", 640, 640);
  h_Esum->Draw();
  TCanvas* c_Emax = new TCanvas("c_Emax", "c_Emax", 640, 640);
  h_Emax->Draw();
  TCanvas* c_Ef = new TCanvas("c_Ef", "c_Ef", 640, 640);
  h_Ef->Draw();
  TCanvas* c_Eb = new TCanvas("c_Eb", "c_Eb", 640, 640);
  h_Eb->Draw();
  TCanvas* c_cos_theta = new TCanvas("c_cos_theta", "c_cos_theta", 640, 640);
  h_cos_theta->Draw();
  TCanvas* c_Ef_Eb = new TCanvas("c_Ef_Eb", "c_Ef_Eb", 640, 640);
  h_Ef_Eb->Draw("colz");
  TCanvas* c_Esum_Eb = new TCanvas("c_Esum_Eb", "c_Esum_Eb", 640, 640);
  h_Esum_Eb->Draw("colz");
}
