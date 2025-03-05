void
filter_rollover(std::string recodata_infilename, std::string rolloverdata_infilename ,
		std::string recodata_outfilename, int min_required_good_lanes = 272)
{

  /** input recodata **/
  auto fin = TFile::Open(recodata_infilename.c_str());
  auto tin = (TTree *)fin->Get("recodata");
  auto nev = tin->GetEntries();
  std::cout << " input " << tin->GetEntries() << " events " << std::endl;
  unsigned short spill;
  tin->SetBranchAddress("spill", &spill);

  /** output recodata **/
  auto fout = TFile::Open(recodata_outfilename.c_str(), "RECREATE");
  auto tout = tin->CloneTree(0);
  
  /** input rolloverdata **/
  auto hRollover = new TH1F("hRollover", "", 100, 0, 100);
  auto hLanes = new TH1F("hLanes", "", 300, 0, 300);
  auto frin = TFile::Open(rolloverdata_infilename.c_str());
  auto lrin = frin->GetListOfKeys();
  for (auto i = 0; i < lrin->GetEntries(); ++i) {
    auto key = (TKey *)lrin->At(i);
    auto grin = (TGraph *)frin->Get(key->GetName());
    for (int ii = 0; ii < grin->GetN(); ++ii) {
      if (grin->GetY()[ii] != 5858) continue;
      hRollover->Fill(grin->GetX()[ii]);
    }
  }
  for (int ispill = 0; ispill < hRollover->GetNbinsX(); ++ispill) {
    if (hRollover->GetBinError(ispill + 1) <= 0) continue;
    hLanes->Fill(hRollover->GetBinContent(ispill + 1));
  }
  int required_good_lanes = hLanes->GetBinLowEdge(hLanes->GetMaximumBin());
  if (required_good_lanes < min_required_good_lanes) required_good_lanes = min_required_good_lanes;
  int expected_spills = hLanes->GetBinContent(hLanes->GetMaximumBin());
  std::cout << required_good_lanes << " required good lanes, will select " << expected_spills << " / " << hLanes->GetEntries() << " spills " << std::endl;
  
  /** loop over events **/
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** selection based on number of good lanes in spill **/
    auto good_lanes = hRollover->GetBinContent(spill + 1);
    if (good_lanes != required_good_lanes) continue;

    tout->Fill();

  }

  std::cout << " output " << tout->GetEntries() << " events " << std::endl;
  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();
  
}

