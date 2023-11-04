std::map<int, std::string> device_int2string = {
    {192, "kc705-192"},
    {193, "kc705-193"},
    {194, "kc705-194"},
    {195, "kc705-195"},
    {196, "kc705-196"},
    {197, "kc705-197"},
    {198, "kc705-198"},
    {207, "kc705-207"}};

std::map<std::string, int> device_string2int = {
    {"kc705-192", 192},
    {"kc705-193", 193},
    {"kc705-194", 194},
    {"kc705-195", 195},
    {"kc705-196", 196},
    {"kc705-197", 197},
    {"kc705-198", 198},
    {"kc705-207", 207}};

TF1 *fine_analysis_fit_function = new TF1("fine_analysis_fit_function", "[0]*(1./(exp((x-[1])/[2])+1))*(1./(exp(([3]-x)/[4])+1))", 10., 200.);

void finetune_macro(std::string filename)
{
  // Output
  std::map<std::string, TH1F *> _TH1F;

  //  Open file with TDC fine distributions
  TFile *input_file = new TFile(filename.c_str());
  std::map<std::string, TH2F *> input_TDC_dist;
  for (auto [current_device_index, current_device] : device_int2string)
  {
    //  Take 2D TDC distribution for fit
    TH2F *current_calib_histo = (TH2F *)(input_file->Get(Form("hFine_%i", current_device_index)));
    if (!current_calib_histo)
      continue;
    input_TDC_dist[current_device] = (TH2F *)current_calib_histo->Clone();
  }
  //  Loop over available TDC distributions
  for (auto [current_device, current_calib_histo] : input_TDC_dist)
  {
    cout << "[INFO] Starting device " << current_device.c_str() << endl;
    _TH1F[Form("hIF_%i", device_string2int[current_device])] = new TH1F(Form("hIF_%i", device_string2int[current_device]), "hIF", 768, 0, 768);
    _TH1F[Form("hCUT_%i", device_string2int[current_device])] = new TH1F(Form("hCUT_%i", device_string2int[current_device]), "hCUT", 768, 0, 768);
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
      if (current_histo->GetEntries() < 1)
      {
        //  Assign dummy values in the calibration object
        _TH1F[Form("hIF_%i", device_string2int[current_device])]->SetBinContent(iBin, IF);
        _TH1F[Form("hCUT_%i", device_string2int[current_device])]->SetBinContent(iBin, CUT);
        continue;
      }
      //  Setup fit function
      fine_analysis_fit_function->SetParameter(0, current_histo->GetBinContent(70));
      fine_analysis_fit_function->SetParameter(1, 100);
      fine_analysis_fit_function->SetParLimits(1, 80, 120);
      fine_analysis_fit_function->SetParameter(2, 0.5);
      fine_analysis_fit_function->SetParLimits(2, 0.1, 1.);
      fine_analysis_fit_function->SetParameter(3, 35);
      fine_analysis_fit_function->SetParLimits(3, 15, 55);
      fine_analysis_fit_function->SetParameter(4, 0.5);
      fine_analysis_fit_function->SetParLimits(4, 0.1, 1.);
      //  Fit Fine distribution
      current_histo->Fit(fine_analysis_fit_function, "IMRESQ", "", 20, 120);
      //  Recover MIN and MAX from fit
      auto maximum = fine_analysis_fit_function->GetParameter(1);
      auto minimum = fine_analysis_fit_function->GetParameter(3);
      //  Calculate IF and CUT
      IF = maximum - minimum;
      CUT = 0.5 * (minimum + maximum);
      //  Assign the values in the calibration object
      _TH1F[Form("hIF_%i", device_string2int[current_device])]->SetBinContent(iBin, IF);
      _TH1F[Form("hCUT_%i", device_string2int[current_device])]->SetBinContent(iBin, CUT);
    }
  }
  TFile *outfile = new TFile("finetune_calib0.root", "RECREATE");
  for (auto [key, object] : _TH1F)
  {
    object->Write();
  }
  outfile->Close();
}