#include "../lib/lightio.h"
#include "../lib/data.h"

const float min_tspill = -0.1; // [s]
const float max_tspill = 0.6; // [s]
const int max_nspill = 100;

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
  auto hTriggerHitsTimeInSpill = new TH2F("hTriggerHitsTimeInSpill", "Trigger readout;trigger time (s);entries", 1000 * (max_tspill - min_tspill), min_tspill, max_tspill, max_nspill, 0, max_nspill);
  auto hTriggerChannelInFrame = new TH2F("hTriggerChannelInFrame", "Trigger readout occupancy;number of triggers;spill;number of frames", 10, 0, 10, max_nspill, 0, max_nspill);
  //  === Timing
  auto hTimingChannelInFrame = new TH2F("hTimingChannelInFrame", "Timing readout occupancy;number of channels;spill;number of frames", 80, 0, 80, max_nspill, 0, max_nspill);
  auto hTimingChannelMap = new TH2F("hTimingChannelMap", "Timing readout occupancy;number of channels (TIME-1);number of channels (TIME-2)", 40, 0, 40, 40, 0, 40);
  auto hTimingTimeResolution = new TH1F("hTimingTimeResolution", ";time chip o - time chip 1 (clock)", 400, -50, 50);
  //  === Cherenkov
  auto hCherenkovChannelInFrame = new TH2F("hCherenkovChannelInFrame", "Cherenkov readout occupancy;number of channels;spill;number of frames", 200, 0, 200, max_nspill, 0, max_nspill);
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
  while (io->next_spill())
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
      hTriggerChannelInFrame->Fill(trigger0_vector.size(), current_spill);
      auto trigger_time = 0.;
      if (trigger0_vector.size() != 0)
      {
        trigger_time = trigger0_vector[0].coarse;
        hTriggerHitsTimeInSpill->Fill((trigger_time + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9, current_spill);
      }

      // === Timing
      auto timing_vector = io->get_timing_vector();
      std::map<int, std::vector<float>> timing_channels_times;
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
        timing_channels_times[timing_channel].push_back(timing_coarse);
        hGenericCoincidenceMapwTrigger[{timing_device, timing_chip}]->Fill(timing_coarse - trigger_time);
      }
      for (auto [timing_channel, timing_time] : timing_channels_times)
      {
        if (timing_channel > 31)
        {
          time_first_chip += timing_time[0];
          contributors_first_chip++;
        }
        else
        {
          time_second_chip += timing_time[0];
          contributors_second_chip++;
        }
      }
      hTimingTimeResolution->Fill((time_second_chip / contributors_second_chip) - (time_first_chip / contributors_first_chip));
      hTimingChannelInFrame->Fill(timing_channels_times.size(), current_spill);
      hTimingChannelMap->Fill(contributors_first_chip, contributors_second_chip);

      // === Cherenkov
      auto cherenkov_vector = io->get_cherenkov_vector();
      std::map<int, std::vector<float>> cherenkov_channels_times;
      auto average_cherenkov_time = 0.;
      for (auto &cherenkov : cherenkov_vector)
      {
        auto cherenkov_coarse = cherenkov.coarse;
        auto cherenkov_device = devices_enum[cherenkov.device];
        auto cherenkov_chip = cherenkov.chip();
        auto cherenkov_channel = cherenkov.eoch();
        cherenkov_channels_times[cherenkov_channel].push_back(cherenkov_coarse);
        hGenericCoincidenceMapwTrigger[{cherenkov_device, cherenkov_chip}]->Fill(cherenkov_coarse - trigger_time);
      }
      hCherenkovChannelInFrame->Fill(cherenkov_channels_times.size(), current_spill);
    }
  }

  gROOT->SetBatch();
  system(Form("mkdir -p %s/", save_dir.c_str()));

  gStyle->SetPalette(kInvertedDarkBodyRadiator);

  auto current_canvas = get_std_canvas();
  gPad->SetLogz();
  hTimingChannelMap->Draw();
  current_canvas->SaveAs(Form("%s/hTimingChannelMap.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);
  hTimingTimeResolution->Fit("gaus", "IMQ", "", -10, 10);
  hTimingTimeResolution->Draw();
  current_canvas->SaveAs(Form("%s/hTimingTimeResolution.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTriggerChannelTimeInSpillIntegrated = hTriggerHitsTimeInSpill->ProjectionX("hTriggerChannelTimeInSpillIntegrated", 0, 99);
  hTriggerChannelTimeInSpillIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTriggerChannelTimeInSpillIntegrated.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTriggerChannelInFrameIntegrated = hTriggerChannelInFrame->ProjectionX("hTriggerChannelInFrameIntegrated", 0, 99);
  hTriggerChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTriggerChannelInFrameIntegrated.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTimingChannelInFrameIntegrated = hTimingChannelInFrame->ProjectionX("hTimingChannelInFrameIntegrated", 0, 99);
  hTimingChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs(Form("%s/hTimingChannelInFrameIntegrated.png", save_dir.c_str()));

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hCherenkovChannelInFrameIntegrated = hCherenkovChannelInFrame->ProjectionX("hCherenkovChannelInFrameIntegrated", 0, 99);
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
  hTriggerChannelInFrame->Write();
  hTimingChannelInFrame->Write();
  hCherenkovChannelInFrame->Write();
  hTimingChannelMap->Write();
  hTriggerChannelTimeInSpillIntegrated->Write();
  hTriggerChannelInFrameIntegrated->Write();
  hTimingChannelInFrameIntegrated->Write();
  hCherenkovChannelInFrameIntegrated->Write();
  hTimingTimeResolution->Write();
  for (auto [iCoordinate, object] : hGenericCoincidenceMapwTrigger)
  {
    object->Write();
  }
  out->Close();

  gROOT->SetBatch(kFALSE);
}
