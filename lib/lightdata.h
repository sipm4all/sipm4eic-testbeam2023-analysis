#pragma once

namespace sipm4eic {

class lightdata {

 public:

  unsigned char device = 0;
  unsigned char index = 0;
  unsigned char coarse = 0;
  unsigned char fine = 0;
  unsigned char tdc = 0;

  lightdata() = default;
  
  lightdata(unsigned char _device, unsigned char _index, unsigned char _coarse, unsigned char _fine, unsigned char _tdc) :
    device(_device),
    index(_index),
    coarse(_coarse),
    fine(_fine),
    tdc(_tdc) { };

  int chip() const { return index / 32; };
  int eoch() const { return index % 64; };
  int cindex() const { return tdc + 4 * index; };
  float time() const;

  /** calibration **/

  static float fine_iif[16][768];
  static float fine_cut[16][768];
  static float fine_off[16][768];
  static bool load_fine_calibration(std::string filename);

};

float lightdata::fine_iif[16][768] = {0.};
float lightdata::fine_cut[16][768] = {0.};
float lightdata::fine_off[16][768] = {0.};

float
lightdata::time() const
{
  auto ci = cindex();
  auto di = device - 192;
  float corr = (float)fine * fine_iif[di][ci] + fine_off[di][ci];
  if (fine >= std::round(fine_cut[di][ci])) corr -= 1.;
  auto time = (float)coarse - corr;
  return time;
}

bool
lightdata::load_fine_calibration(std::string filename)
{
  std::cout << " --- loading fine calibration: " << filename << std::endl;
  auto fin = TFile::Open(filename.c_str());
  for (int i = 0; i < 16; ++i) {
    auto device = 192 + i;
    auto hFine_iif = (TH1 *)fin->Get(Form("hIIF_%d", device));
    auto hFine_cut = (TH1 *)fin->Get(Form("hCUT_%d", device));
    auto hFine_off = (TH1 *)fin->Get(Form("hOFF_%d", device));
    for (int j = 0; j < 768; ++j) {
      fine_iif[i][j] = hFine_iif ? hFine_iif->GetBinContent(j + 1)  : 0.;
      fine_cut[i][j] = hFine_cut ? hFine_cut->GetBinContent(j + 1) : 0.;
      fine_off[i][j] = hFine_off ? hFine_off->GetBinContent(j + 1) : -fine_cut[i][j] * fine_iif[i][j];
    }
  }
  fin->Close();
  return true;
}

} /** namespace sipm4eic **/
