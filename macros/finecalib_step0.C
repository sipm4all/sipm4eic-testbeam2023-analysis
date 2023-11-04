std::vector<int> devices_indices = {192, 193, 194, 195, 196, 197, 198, 207};

TF1 *f_fit = new TF1("f_fit", "[0]*(1./(exp((x-[1])/[2])+1))*(1./(exp(([3]-x)/[4])+1))", 10., 200.);

void finecalib_step0(std::string input_filename, std::string output_filename = "finetune_calib0.root", int reference_bin_for_distribution_heigth = 70, int minimum_entries = 1)
{
  // Output
  std::map<std::string, TH1F *> _TH1F;

  //  Open file with TDC fine distributions
  TFile *input_file = new TFile(input_filename.c_str());
  //  Loop over available TDC distributions
  for (auto current_device_index : devices_indices)
  {
    TH2F *current_calib_histo = (TH2F *)(input_file->Get(Form("hFine_%i", current_device_index)));
    if (!current_calib_histo)
      continue;
    cout << "[INFO] Starting device " << current_device_index << endl;
    _TH1F[Form("hIF_%i", current_device_index)] = new TH1F(Form("hIF_%i", current_device_index), "hIF", 768, 0, 768);
    _TH1F[Form("hCUT_%i", current_device_index)] = new TH1F(Form("hCUT_%i", current_device_index), "hCUT", 768, 0, 768);
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
        _TH1F[Form("hIF_%i", current_device_index)]->SetBinContent(iBin, IF);
        _TH1F[Form("hCUT_%i", current_device_index)]->SetBinContent(iBin, CUT);
        continue;
      }
      //  Setup fit function
      f_fit->SetParameter(0, current_histo->GetBinContent(reference_bin_for_distribution_heigth));
      f_fit->SetParameter(1, 100);
      f_fit->SetParLimits(1, 80, 120);
      f_fit->SetParameter(2, 0.5);
      f_fit->SetParLimits(2, 0.1, 1.);
      f_fit->SetParameter(3, 35);
      f_fit->SetParLimits(3, 15, 55);
      f_fit->SetParameter(4, 0.5);
      f_fit->SetParLimits(4, 0.1, 1.);
      //  Fit Fine distribution
      current_histo->Fit(f_fit, "IMRESQ", "", 20, 120);
      //  Recover MIN and MAX from fit
      auto maximum = f_fit->GetParameter(1);
      auto minimum = f_fit->GetParameter(3);
      //  Calculate IF and CUT
      IF = maximum - minimum;
      CUT = 0.5 * (minimum + maximum);
      //  Assign the values in the calibration object
      _TH1F[Form("hIF_%i", current_device_index)]->SetBinContent(iBin, IF);
      _TH1F[Form("hCUT_%i", current_device_index)]->SetBinContent(iBin, CUT);
    }
  }
  TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
  for (auto [key, object] : _TH1F)
  {
    object->Write();
  }
  outfile->Close();
}