#include "../lib/lightio.h"
#include "../lib/data.h"

void timing_res(std::string input_fine_data, std::string input_fine_calibration = "", std::string output_file = "output.root")
{
  //  Output
  std::map<std::string, TH1F *> _TH1F;
  _TH1F["hTiming_delta_references"] = new TH1F("hTiming_delta_references", "hTiming_delta_references", 100, -10., 10.);
  std::map<std::string, TH2F *> _TH2F;
  _TH2F["hTiming_hits_per_chip"] = new TH2F("hTiming_hits_per_chip", "hTiming_hits_per_chip", 33, -0.5, 32.5, 33, -0.5, 32.5);
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
      //  Calculate reference time
      std::map<int, float> average_time;
      std::map<int, std::map<int, float>> channel_map;
      for (auto &hit : timing_vector)
      {
        auto current_device = hit.device;
        auto current_chip = hit.chip();
        auto current_eoch = hit.eoch();
        auto current_cindex = hit.cindex();
        if ((channel_map[current_chip].count(current_eoch) == 0) || hit.time() < channel_map[current_chip][current_eoch])
          channel_map[current_chip][current_eoch] = hit.time();
      }
      //  Check required hits
      auto hits_on_chip4 = channel_map[4].size();
      auto hits_on_chip5 = channel_map[5].size();
      _TH2F["hTiming_hits_per_chip"]->Fill(hits_on_chip4, hits_on_chip5);
      //  Must be all firing
      if ((hits_on_chip4 < 32) || (hits_on_chip5 < 32))
        continue;
      //  Measure average times per chip
      for (auto [current_chip, all_channels] : channel_map)
      {
        for (auto [current_channel, current_hit_time] : all_channels)
        {
          average_time[current_chip] += current_hit_time;
        }
        average_time[current_chip] /= all_channels.size();
      }
      //  Check coincidence window
      _TH1F["hTiming_delta_references"]->Fill((average_time[4] - average_time[5]) * sipm4eic::data::coarse_to_ns);
      if (fabs(average_time[4] - average_time[5]) > 5 / (sipm4eic::data::coarse_to_ns))
        continue;
      //  Measure reference time
      auto reference_time = 0.5 * (average_time[4] + average_time[5]);
      //  Loop on single sensors
      std::map<int, std::map<int, std::map<int, float>>> channel_map_cherenkov;
      for (auto &hit : cherenkov_vector)
      {
        auto current_device = (int)(hit.device);
        auto current_chip = hit.chip();
        auto current_eoch = hit.eoch();
        auto current_cindex = hit.cindex();
        if ((channel_map_cherenkov[current_device][current_chip].count(current_eoch) == 0) || hit.time() < channel_map_cherenkov[current_device][current_chip][current_eoch])
          channel_map_cherenkov[current_device][current_chip][current_eoch] = hit.time();
      }
      for (auto [current_device, all_chips] : channel_map_cherenkov)
      {
        for (auto [current_chip, all_channels] : all_chips)
        {
          for (auto [current_channel, current_hit_time] : all_channels)
          {
            std::string current_histo_name = Form("hDelta_hit_minus_timing_device_%d_chip_%d_channel_%d", current_device, current_chip, current_channel);
            if (!_TH1F[current_histo_name.c_str()])
              _TH1F[current_histo_name.c_str()] = new TH1F(current_histo_name.c_str(), current_histo_name.c_str(), 500, -10*sipm4eic::data::coarse_to_ns, 10*sipm4eic::data::coarse_to_ns);
            _TH1F[current_histo_name.c_str()]->Fill((current_hit_time - reference_time)*sipm4eic::data::coarse_to_ns);
          }
        }
      }
    }
  }
  //
  //  Save to file
  auto fout = TFile::Open(output_file.c_str(), "RECREATE");
  for (auto [current_histo_name, current_histo] : _TH1F)
    current_histo->Write();
  for (auto [current_histo_name, current_histo] : _TH2F)
    current_histo->Write();
  fout->Close();
}
