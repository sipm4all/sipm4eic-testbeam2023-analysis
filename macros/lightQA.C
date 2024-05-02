#include "../lib/lightio.h"
#include "../lib/data.h"

#define TESTBEAM2023

int TRIGGER0_device = 192;
int TRIGGER0_offset = 112;

#ifdef TESTBEAM2023
int TIMING1_device = 207, TIMING1_chip = 4;
int TIMING2_device = 207, TIMING2_chip = 5;
#else
int TIMING1_device = 200, TIMING1_chip = 5;
int TIMING2_device = 201, TIMING2_chip = 5;
#endif

int TRACKING1_device = 200, TRACKING1_chip = 4;
int TRACKING2_device = 201, TRACKING2_chip = 4;

const float min_tspill = -0.1; // [s]
const float max_tspill = 0.6;  // [s]
const int max_nspill = 100;

const float min_tdelta = -25; // [clock]
const float max_tdelta = 25;  // [clock]

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

void lightQA(std::string input_file = "lightdata.root", std::string output_file = "out.root", std::string save_dir = "./images/")
{
  //  Define output objects
  //  === Trigger
  auto hTriggerHitsTimeInSpill = new TH1F("hTriggerHitsTimeInSpill", "TRIGGER readout;trigger time (s);entries", 1000 * (max_tspill - min_tspill), min_tspill, max_tspill);
  auto hTriggerChannelInFrameIntegrated = new TH1F("hTriggerChannelInFrame", "TRIGGER readout occupancy;number of triggers;frame", 10, 0, 10);

  //  === Tracking
  auto hTrackingHitsTimeInSpill = new TH1F("hTrackingHitsTimeInSpill", "TRACKING readout;trigger time (s);entries", 1000 * (max_tspill - min_tspill), min_tspill, max_tspill);
  auto hTrackingChannelInFrameIntegrated = new TH1F("hTrackingChannelInFrame", "TRACKING readout occupancy;number of triggers;frame", 10, 0, 10);

  //  === Timing
  auto hTimingHitsTimeInSpill = new TH1F("hTimingHitsTimeInSpill", "TIMING readout;timing time (s);entries", 1000 * (max_tspill - min_tspill), min_tspill, max_tspill);
  auto hTimingChannelInFrameIntegrated = new TH1F("hTimingChannelInFrame", "TIMING readout occupancy;number of channels;frame", 80, 0, 80);
  auto hTimingChannelMap = new TH2F("hTimingChannelMap", "TIMING readout occupancy;number of channels (TIMING 1);number of channels (TIMING 2)", 40, 0, 40, 40, 0, 40);
  auto hTimingTimeResolution = new TH1F("hTimingTimeResolution", "TIMING time coincidences;TIMING 1 - TIMING 2 (clock cycles);", max_tdelta - min_tdelta, min_tdelta, max_tdelta);

  //  === Cherenkov
  auto hCherenkovHitsTimeInSpill = new TH1F("hCherenkovHitsTimeInSpill", "CHERENKOV readout;cherenkov time (s);entries", 1000 * (max_tspill - min_tspill), min_tspill, max_tspill);
  auto hCherenkovChannelInFrameIntegrated = new TH1F("hCherenkovChannelInFrame", "CHERENKOV readout occupancy;number of channels;frame", 200, 0, 200);

  //  === General
  std::map<std::array<int, 2>, TH1F *> hGenericCoincidenceMapwTrigger;
  for (auto [device_id, enumerator] : devices_enum)
  {
    for (auto iChip = 0; iChip < 6; iChip++)
    {
      hGenericCoincidenceMapwTrigger[{enumerator, iChip}] = new TH1F(Form("hGenericCoincidenceMap_d%i_c%i", device_id, iChip), Form("Time coincidences (kc705-%i, chip-%i);hit - trigger time (clock cycles);entries", device_id, iChip), 256 * 2, -256, 256);
    }
  }

  cout << "[INFO] Starting the light QA: reading " << input_file << endl;
  sipm4eic::lightio *io = new sipm4eic::lightio();
  io->read_from_tree(input_file);
  cout << "[INFO] Finished reading from tree in file " << input_file << endl;

  //  Loop on spills
  auto full_frames = 0;
  auto current_spill = 0;
  while (io->next_spill() && (current_spill < 3))
  {
    current_spill++;
    cout << "[INFO] Start spill: " << current_spill << endl;

    //  Loop on frames
    auto current_frame = 0;
    while (io->next_frame())
    {
      //  === Frame info
      current_frame++;
      auto frame_id = io->frame[current_frame];

      // === Trigger
      auto trigger0_vector = io->get_trigger0_vector();
      hTriggerChannelInFrameIntegrated->Fill(trigger0_vector.size(), current_spill);
      auto trigger_coarse = 0.;
      if (trigger0_vector.size() != 0)
      {
        trigger_coarse = trigger0_vector[0].coarse;
        hTriggerHitsTimeInSpill->Fill((trigger_coarse + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9);
      }

      // === Tracking
      auto tracking_vector = io->get_tracking_vector();
      std::map<int, std::vector<float>> tracking_channels_times;
      for (auto &tracking : tracking_vector)
      {
        auto tracking_coarse = tracking.coarse;
        auto tracking_device = devices_enum[tracking.device];
        auto tracking_chip = tracking.chip();
        auto tracking_channel = tracking.eoch();
        tracking_channels_times[tracking_channel].push_back(tracking_coarse);
        hTrackingHitsTimeInSpill->Fill((trigger_coarse + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9);
      }
      hTrackingChannelInFrameIntegrated->Fill(tracking_channels_times.size());

      // === Timing
      auto timing_vector = io->get_timing_vector();
      std::map<int, std::vector<float>> timing_first_channels_times, timing_second_channels_times;
      auto contributors_first_chip = 0;
      auto contributors_second_chip = 0;
      auto time_first_chip = 0.;
      auto time_second_chip = 0.;
      for (auto &timing : timing_vector)
      {
        auto timing_coarse = timing.coarse;
        auto timing_device = devices_enum[timing.device];
        auto timing_chip = timing.chip();
        auto timing_channel = timing.eoch();
        if (trigger0_vector.size() != 0)
          hGenericCoincidenceMapwTrigger[{timing_device, timing_chip}]->Fill(timing_coarse - trigger_coarse);
        hTimingHitsTimeInSpill->Fill((timing_coarse + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9);
        if (timing.device == TIMING1_device && timing.chip() == TIMING1_chip)
          timing_first_channels_times[timing_channel].push_back(timing_coarse);
        else if (timing.device == TIMING2_device && timing.chip() == TIMING2_chip)
          timing_second_channels_times[timing_channel].push_back(timing_coarse);
      }
      for (auto [timing_channel, timing_time] : timing_first_channels_times)
      {
        time_first_chip += timing_time[0];
        contributors_first_chip++;
      }
      for (auto [timing_channel, timing_time] : timing_second_channels_times)
      {
        time_second_chip += timing_time[0];
        contributors_second_chip++;
      }
      if (contributors_second_chip != 0 && contributors_first_chip != 0) // otherwise we risk division by zero
        hTimingTimeResolution->Fill((time_first_chip / contributors_first_chip) - (time_second_chip / contributors_second_chip));

      hTimingChannelInFrameIntegrated->Fill(timing_first_channels_times.size() + timing_second_channels_times.size(), current_spill);
      hTimingChannelMap->Fill(contributors_first_chip, contributors_second_chip);

      // === Cherenkov
      auto cherenkov_vector = io->get_cherenkov_vector();
      std::map<int, std::vector<float>> cherenkov_channels_times;
      for (auto &cherenkov : cherenkov_vector)
      {
        auto cherenkov_coarse = cherenkov.coarse;
        auto cherenkov_device = devices_enum[cherenkov.device];
        auto cherenkov_chip = cherenkov.chip();
        auto cherenkov_channel = cherenkov.eoch();
        cherenkov_channels_times[cherenkov_channel].push_back(cherenkov_coarse);
        if (trigger0_vector.size() != 0)
          hGenericCoincidenceMapwTrigger[{cherenkov_device, cherenkov_chip}]->Fill(cherenkov_coarse - trigger_coarse);
        hCherenkovHitsTimeInSpill->Fill((cherenkov_coarse + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9);
      }
      hCherenkovChannelInFrameIntegrated->Fill(cherenkov_channels_times.size(), current_spill);
    }
  }

  //  === Graphics
  gROOT->SetBatch();
  system(Form("mkdir -p %s/", save_dir.c_str()));

  gStyle->SetPalette(kInvertedDarkBodyRadiator);

  //  === === Trigger
  auto current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTriggerHitsTimeInSpill->Draw();
  current_canvas->SaveAs(Form("%s/hTriggerHitsTimeInSpill.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTriggerChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTriggerChannelInFrameIntegrated.png", save_dir.c_str()));

  //  === === Tracking
  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTrackingHitsTimeInSpill->Draw();
  current_canvas->SaveAs(Form("%s/hTrackingHitsTimeInSpill.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTrackingChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTrackingChannelInFrameIntegrated.png", save_dir.c_str()));

  //  === === Timing
  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTimingHitsTimeInSpill->Draw();
  current_canvas->SaveAs(Form("%s/hTimingHitsTimeInSpill.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTimingChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTimingChannelInFrameIntegrated.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogz();
  hTimingChannelMap->Draw();
  current_canvas->SaveAs(Form("%s/hTimingChannelMap.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hTimingTimeResolution->Draw();
  current_canvas->SaveAs(Form("%s/hTimingTimeResolution.png", save_dir.c_str()));

  //  === === Cherenkov
  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hCherenkovHitsTimeInSpill->Draw();
  current_canvas->SaveAs(Form("%s/hCherenkovHitsTimeInSpill.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  hCherenkovChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hCherenkovChannelInFrameIntegrated.png", save_dir.c_str()));

  for (auto [iCoordinate, object] : hGenericCoincidenceMapwTrigger)
  {
    auto iDevice = iCoordinate[0];
    auto iChip = iCoordinate[1];
    current_canvas = get_std_canvas();
    gPad->SetLogy();
    object->Draw();
    current_canvas->SaveAs(Form("%s/%s.png", save_dir.c_str(), object->GetName()));
  }

  TFile *out = new TFile(output_file.c_str(), "RECREATE");
  hTriggerHitsTimeInSpill->Write();
  hTriggerChannelInFrameIntegrated->Write();
  hTimingHitsTimeInSpill->Write();
  hCherenkovHitsTimeInSpill->Write();
  for (auto [iCoordinate, object] : hGenericCoincidenceMapwTrigger)
  {
    object->Write();
  }
  out->Close();

  gROOT->SetBatch(kFALSE);
}
