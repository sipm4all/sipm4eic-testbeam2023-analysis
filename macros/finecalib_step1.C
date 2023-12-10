
#include "../lib/lightio.h"

TF2 *f2 = new TF2("f", "[3]*exp(-0.5*(((y-(x/[0]+[2]-(x>=[1])))/([4]))**2))", 30., 140., -6, 6);

std::tuple<std::pair<float, float>, std::pair<float, float>, std::pair<float, float>> finecalib_step1(std::string input_filename, std::string input_fine_calibration, int device, int cindex)
{
  TFile *input_file = new TFile(input_filename.c_str());
  auto h2 = (TH2F *)(input_file->Get("hDeltaT"));
  auto h1 = h2->ProfileX("h1");

  sipm4eic::lightdata::load_fine_calibration(input_fine_calibration);

  f2->SetParameter(0, sipm4eic::lightdata::fine_if[device - 192][cindex]);
  f2->SetParameter(1, sipm4eic::lightdata::fine_cut[device - 192][cindex]);
  f2->SetParameter(2, sipm4eic::lightdata::fine_off[device - 192][cindex]);
  f2->SetParameter(3, 1);
  f2->SetParameter(4, 0.1);

  h2->Fit(f2,"N0","N0");

  return {{f2->GetParameter(0), f2->GetParError(0)},{f2->GetParameter(1), f2->GetParError(1)},{f2->GetParameter(2), f2->GetParError(2)}};
}