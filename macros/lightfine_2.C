#include "../lib/lightio.h"
//
//  Global variables
auto n4_fired_required = 32;
auto n5_fired_required = 32;
//
//  Fill fine tune calibration of device and cindex
void lightfine_2(std::string input_fine_data, std::string input_fine_calibration, int device, int cindex, std::string output_file)
{
  //  Output  ---
  TH2F *hN = new TH2F("hN", "hN", 33, -0.5, 32.5, 33, -0.5, 32.5);
  TH1F *hDeltaTref = new TH1F("hDeltaT", "hDeltaT", 1000, -8., 8.);
  TH2F *hDeltaT = new TH2F("hDeltaT", "hDeltaT", 256, 0, 256, 1000, -8., 8.);
  //  Read data
  sipm4eic::lightio io;
  io.read_from_tree(input_fine_data);
  //  Read calibration
  sipm4eic::lightdata::load_fine_calibration(input_fine_calibration);
  //  Loop on data
  while (io.next_spill())
  {
    while (io.next_frame())
    {
      //  Recover hits vectors
      auto timing_vector = io.get_timing_vector();
      auto cherenkov_vector = io.get_cherenkov_vector();
      float Tdut = -1.;
      float Tdutfine = -1.;
      //  Calculate reference time
      std::map<int, std::map<int, int>> actual_participants;
      std::map<int, std::map<int, sipm4eic::lightdata>> channel_map;
      for (auto &hit : timing_vector)
      {
        auto current_device = hit.device;
        auto current_chip = hit.chip();
        auto current_eoch = hit.eoch();
        auto current_cindex = hit.cindex();
        auto target_eoch = cindex / 4 - 32 * current_chip;
        actual_participants[current_chip][current_eoch] = 1;
        if ((device == (int)(current_device)) && ((int)(current_cindex) == cindex))
        {
          if ((Tdut < hit.coarse) && (Tdut > 0))
            continue;
          Tdut = hit.coarse;
          Tdutfine = hit.fine;
          continue;
        }
        if ((device == (int)(current_device)) && (target_eoch == current_eoch))
          continue;
        if ((channel_map[current_chip].count(current_eoch) == 0) || hit.time() < channel_map[current_chip][current_eoch].time())
          channel_map[current_chip][current_eoch] = hit;
      }
      //  Check required hits
      auto N4 = channel_map[4].size();
      auto N5 = channel_map[5].size();
      auto N4act = actual_participants[4].size();
      auto N5act = actual_participants[5].size();
      hN->Fill(N4act, N5act);
      if ((N4act < n4_fired_required) || (N5act < n5_fired_required))
        continue;
      auto T4 = 0.;
      auto T5 = 0.;
      for (auto [chip, channels] : channel_map)
      {
        for (auto [channel, hit] : channels)
        {
          if (chip == 4)
            T4 += hit.time();
          if (chip == 5)
            T5 += hit.time();
        }
      }
      T4 /= N4;
      T5 /= N5;
      auto Tref = 0.5 * (T4 + T5);
      hDeltaTref->Fill(T4 - T5);
      for (auto &hit : cherenkov_vector)
      {
        if ((device == (int)(hit.device)) && (hit.cindex() == cindex))
        {
          Tdut = hit.coarse;
          Tdutfine = hit.fine;
          break;
        }
      }
      if (Tdut > 0)
        hDeltaT->Fill(Tdutfine, Tdut - Tref);
    }
  }
  //  Save to file
  auto fout = TFile::Open(output_file.c_str(), "RECREATE");
  hN->Write();
  hDeltaT->Write();
  hDeltaTref->Write();
  fout->Close();
}
