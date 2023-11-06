std::vector<int> devices_indices = {192, 193, 194, 195, 196, 197, 198, 207};

TF1 *f_fit = new TF1("f_fit", "[0]*(1./(exp(([1]-x)/[2])+1))*(1./(exp((x-[3])/[4])+1))", 10., 200.);

void finecalib_step0(std::string input_filename, std::string output_filename = "finecalib_step0.root", int minimum_entries = 1000, bool use_fit = true)
{
  //  Open file with TDC fine distributions
  TFile *input_file = new TFile(input_filename.c_str());
  TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
  //  Loop over available TDC distributions
  for (auto current_device_index : devices_indices)
  {
    TH2F *current_calib_histo = (TH2F *)(input_file->Get(Form("hFine_%i", current_device_index)));
    if (!current_calib_histo)
      continue;
    cout << "[INFO] Starting device " << current_device_index << endl;
    TH1F *hIF = new TH1F(Form("hIF_%i", current_device_index), "hIF", 768, 0, 768);
    TH1F *hCUT = new TH1F(Form("hCUT_%i", current_device_index), "hCUT", 768, 0, 768);
    //  Take 2D TDC distribution for fit
    for (auto iBin = 1; iBin <= current_calib_histo->GetNbinsX(); iBin++)
    {
      //  Take slice to act on single TDC
      auto current_histo = current_calib_histo->ProjectionY("tmp", iBin, iBin);
      //  Book results variable
      auto IF = 0.;
      auto CUT = 0.;
      auto OFFSET = 0.;
      //  If no entries, skip fit
      if (current_histo->GetEntries() < minimum_entries)
      {
        //  Assign dummy values in the calibration object
        hIF->SetBinContent(iBin, IF);
        hCUT->SetBinContent(iBin, CUT);
        continue;
      }
      //  Setup fit function
      auto height_guess = current_histo->GetBinContent(current_histo->GetMaximumBin());
      auto bin_step = 3;
      auto critical_value = 0.25 * height_guess;
      auto min_guess = 0.;
      for (int jBin = 1; jBin <= current_histo->GetNbinsX() - bin_step; jBin += bin_step)
      {
        double y1 = current_histo->GetBinContent(jBin);
        double y2 = current_histo->GetBinContent(jBin + bin_step);
        if (fabs(y1 - y2) > critical_value)
        {
          min_guess = (jBin + bin_step * 0.5);
          break;
        }
      }
      auto max_guess = 0.;
      for (int jBin = current_histo->GetNbinsX(); jBin >= bin_step; jBin -= bin_step)
      {
        double y1 = current_histo->GetBinContent(jBin);
        double y2 = current_histo->GetBinContent(jBin - bin_step);
        if (fabs(y1 - y2) > critical_value)
        {
          max_guess = (jBin - bin_step * 0.5);
          break;
        }
      }
      f_fit->SetParameter(0, height_guess);
      f_fit->SetParameter(1, min_guess);
      f_fit->SetParLimits(1, min_guess * 0.8, min_guess * 1.2);
      f_fit->SetParameter(2, 0.5);
      f_fit->SetParLimits(2, 0.1, 1.);
      f_fit->SetParameter(3, max_guess);
      f_fit->SetParLimits(3, max_guess * 0.8, max_guess * 1.2);
      f_fit->SetParameter(4, 0.5);
      f_fit->SetParLimits(4, 0.1, 1.);
      //  Fit Fine distribution
      if (use_fit)
        current_histo->Fit(f_fit, "IMRESQ", "", 20, 120);
      //  Recover MIN and MAX from fit
      auto minimum = f_fit->GetParameter(1);
      auto maximum = f_fit->GetParameter(3);
      //  Calculate IF and CUT
      IF = maximum - minimum;
      CUT = 0.5 * (minimum + maximum);
      //  Assign the values in the calibration object
      hIF->SetBinContent(iBin, IF);
      hCUT->SetBinContent(iBin, CUT);
    }
    outfile->cd();
    hIF->Write();
    hCUT->Write();
  }
  outfile->Close();
}