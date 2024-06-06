void
hitmap(std::string recodata_infilename)
{

  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto fin = TFile::Open(recodata_infilename.c_str());
  auto tin = (TTree *)fin->Get("recodata");
  auto nev = tin->GetEntries();
  tin->SetBranchAddress("n", &n);
  tin->SetBranchAddress("x", &x);
  tin->SetBranchAddress("y", &y);
  tin->SetBranchAddress("t", &t);

  auto hXY = new TH2F("hMap", ";x (mm);y (mm)", 396, -99, 99, 396, -99, 99);
  hXY->SetStats(false);
  auto hT = new TH1F("hT", ";t (ns);", 50, -78.125, 78.125);  

  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);
    for (int i = 0 ; i < n; ++i) {
      hXY->Fill(gRandom->Uniform(x[i] - 1.5, x[i] + 1.5), gRandom->Uniform(y[i] - 1.5, y[i] + 1.5));
      hT->Fill(t[i]);
    } 
  }

  auto cXY = new TCanvas("cXY", "cXY", 800, 800);
  cXY->SetMargin(0.15, 0.15, 0.15, 0.15);
  cXY->SetLogz();
  hXY->Draw("colz");
  cXY->SaveAs("hitmap.png");
}
