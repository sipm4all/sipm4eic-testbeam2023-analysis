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

  bool operator<(const lightdata &rhs) const { return time() < rhs.time(); };

  
  int chip() const { return index / 32; };
  int eoch() const { return index % 64; };
  int cindex() const { return tdc + 4 * index; };
  float time() const;
  bool near_cut(float dist = 2) const { return fabs(fine - fine_cut[device - 192][cindex()]) < dist; };

  /** calibration **/

  static constexpr float coarse_to_ns = 3.125;
  
  static float fine_iif[16][768];
  static float fine_cut[16][768];
  static float fine_off[16][768];
  static bool load_fine_calibration(std::string filename);
  static bool write_fine_calibration(std::string filename);

};

float lightdata::fine_iif[16][768] = {0.};
float lightdata::fine_cut[16][768] = {0.};
float lightdata::fine_off[16][768] = {0.};

float
lightdata::time() const
{
  auto ci = cindex();
  auto di = device - 192;
  if (di < 0 || di > 15) return (float)coarse;
  float corr = (float)fine * fine_iif[di][ci] + fine_off[di][ci];
  if (fine >= std::round(fine_cut[di][ci])) corr -= 1.;
  auto time = (float)coarse - corr;
  return time;
}

bool
lightdata::load_fine_calibration(std::string filename)
{
  if (filename.empty()) {
    std::cout << " --- fine calibration not requested " << std::endl;
    return true;
  }
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

bool
lightdata::write_fine_calibration(std::string filename)
{
  std::cout << " --- writing fine calibration: " << filename << std::endl;
  auto fin = TFile::Open(filename.c_str(), "RECREATE");
  for (int i = 0; i < 16; ++i) {
    auto device = 192 + i;
    auto hFine_iif = new TH1F(Form("hIIF_%d", device), "hIIF", 768, 0, 768);
    auto hFine_cut = new TH1F(Form("hCUT_%d", device), "hCUT", 768, 0, 768);
    auto hFine_off = new TH1F(Form("hOFF_%d", device), "hOFF", 768, 0, 768);
    for (int j = 0; j < 768; ++j) {
      hFine_iif->SetBinContent(j + 1, fine_iif[i][j]);
      hFine_cut->SetBinContent(j + 1, fine_cut[i][j]);
      hFine_off->SetBinContent(j + 1, fine_off[i][j]);
    }
    hFine_iif->Write();
    hFine_cut->Write();
    hFine_off->Write();
  }
  fin->Close();
  return true;
}

} /** namespace sipm4eic **/
