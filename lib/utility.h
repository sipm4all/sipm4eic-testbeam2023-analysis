#pragma once

//  Data structures
//  === recodata
//  === === Structure
struct recodata
{
  unsigned short n = 0;
  float x[1024];
  float y[1024];
  float t[1024];
};
//  === === Load data from tree
void load_data(TTree *reco_data_tree, recodata &target_data_struct)
{
  reco_data_tree->SetBranchAddress("n", &target_data_struct.n);
  reco_data_tree->SetBranchAddress("x", &target_data_struct.x);
  reco_data_tree->SetBranchAddress("y", &target_data_struct.y);
  reco_data_tree->SetBranchAddress("t", &target_data_struct.t);
}
//  === === Load data from file
TTree *load_data(std::string filename, recodata &target_data_struct)
{
  auto input_file = TFile::Open(filename.c_str());
  auto input_tree = (TTree *)input_file->Get("recodata");
  load_data(input_tree, target_data_struct);
  return input_tree;
}
//  === Devices
std::map<std::string, int> devices_name = {
    {"kc705-192", 192},
    {"kc705-193", 193},
    {"kc705-194", 194},
    {"kc705-195", 195},
    {"kc705-196", 196},
    {"kc705-197", 197},
    {"kc705-198", 198},
    {"kc705-199", 199},
    {"kc705-200", 200},
    {"kc705-201", 201},
    {"kc705-202", 202},
    {"kc705-207", 207}};
std::map<int, int> devices_enum = {
    {192, 0},
    {193, 1},
    {194, 2},
    {195, 3},
    {196, 4},
    {197, 5},
    {198, 6},
    {199, 7},
    {200, 8},
    {201, 9},
    {202, 10},
    {207, 11}};

// General Info
const float kSiPM_x_dimension = 1.5;
const float kSiPM_y_dimension = 1.5;

//  General utilities
//  === Cartesian to Polar
inline std::array<float, 2> cartesian_to_polar(std::array<float, 2> target, std::array<float, 2> center_shift = {0., 0.})
{
  float x_new_coordinate = target[0] - center_shift[0];
  float y_new_coordinate = target[1] - center_shift[1];
  float r_coordinate = TMath::Sqrt(x_new_coordinate * x_new_coordinate + y_new_coordinate * y_new_coordinate);
  float phi_coordinate = TMath::ATan2(y_new_coordinate, x_new_coordinate);
  return {r_coordinate, phi_coordinate};
}
inline std::array<float, 2> polar_to_cartesian(std::array<float, 2> target, std::array<float, 2> center_shift = {0., 0.})
{
  float r_new_coordinate = target[0] - center_shift[0];
  float phi_new_coordinate = target[1] - center_shift[1];
  float x_coordinate = target[0] * TMath::Cos(target[1]);
  float y_coordinate = target[0] * TMath::Sin(target[1]);
  return {x_coordinate, y_coordinate};
}
//  === Filler functions
bool is_within_SiPM(std::array<float, 2> cartesian_coordinates, std::array<float, 2> SiPM_center)
{
  bool is_within_SiPM = false;
  auto x_distance = fabs(cartesian_coordinates[0] - SiPM_center[0]);
  auto y_distance = fabs(cartesian_coordinates[1] - SiPM_center[1]);
  if ((x_distance) <= kSiPM_x_dimension && fabs(y_distance) <= kSiPM_y_dimension)
    is_within_SiPM = true;
  return is_within_SiPM;
}
template <bool rnd_smearing = true>
void fill_persistance(TH2F *hTarget, recodata reco_data)
{
  for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
    hTarget->Fill(rnd_smearing ? gRandom->Uniform(reco_data.x[iPnt] - 1.5, reco_data.x[iPnt] + 1.5) : reco_data.x[iPnt], rnd_smearing ? gRandom->Uniform(reco_data.y[iPnt] - 1.5, reco_data.y[iPnt] + 1.5) : reco_data.y[iPnt]);
}
template <bool allow_overlap = false>
void fill_with_SiPM_coverage(TH2F *target, std::array<float, 2> SiPM_center)
{
  float x_center_bin = target->GetXaxis()->FindBin(SiPM_center[0]);
  float y_center_bin = target->GetYaxis()->FindBin(SiPM_center[1]);
  float x_center_bin_center = target->GetXaxis()->GetBinCenter(x_center_bin);
  float y_center_bin_center = target->GetYaxis()->GetBinCenter(y_center_bin);
  float current_center_bin_content = target->GetBinContent(target->FindBin(x_center_bin, y_center_bin));

  //  Did we set this SiPM already?
  //  * No overlap is permitted
  if ((current_center_bin_content > 0) && !allow_overlap)
    return;

  //  Start with know center bin
  for (auto xBin = 0; xBin < target->GetNbinsX(); xBin++)
  {
    //  Check at least one fill has been done
    bool at_least_one_fill = false;

    float current_bin_center_x_high = (float)(target->GetXaxis()->GetBinCenter(x_center_bin + xBin));
    float current_bin_center_x_low_ = (float)(target->GetXaxis()->GetBinCenter(x_center_bin - xBin));
    for (auto yBin = 0; yBin < target->GetNbinsY(); yBin++)
    {
      float current_bin_center_y_high = (float)(target->GetYaxis()->GetBinCenter(y_center_bin + yBin));
      float current_bin_center_y_low_ = (float)(target->GetYaxis()->GetBinCenter(y_center_bin - yBin));

      //  Set bin contents
      if (is_within_SiPM({current_bin_center_x_high, current_bin_center_y_high}, SiPM_center))
      {
        target->SetBinContent(x_center_bin + xBin, y_center_bin + yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_high, current_bin_center_y_low_}, SiPM_center))
      {
        target->SetBinContent(x_center_bin + xBin, y_center_bin - yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_low_, current_bin_center_y_high}, SiPM_center))
      {
        target->SetBinContent(x_center_bin - xBin, y_center_bin + yBin, 1);
        at_least_one_fill = true;
      }
      if (is_within_SiPM({current_bin_center_x_low_, current_bin_center_y_low_}, SiPM_center))
      {
        target->SetBinContent(x_center_bin - xBin, y_center_bin - yBin, 1);
        at_least_one_fill = true;
      }

      //  return if we stopped filling
      if (!at_least_one_fill)
        break;
    }
    if (!is_within_SiPM({current_bin_center_x_high, y_center_bin_center}, SiPM_center) && !is_within_SiPM({current_bin_center_x_low_, y_center_bin_center}, SiPM_center))
      return;
  }
}
//  === Graphical helpers
TCanvas *get_std_canvas()
{
  TCanvas *result = new TCanvas("", "", 800, 800);
  gStyle->SetOptStat(1110);
  result->SetRightMargin(0.15);
  result->SetTopMargin(0.15);
  result->SetLeftMargin(0.15);
  result->SetBottomMargin(0.15);
  return result;
}